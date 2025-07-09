#ifndef FOF_H
#define FOF_H

#include "photoNs.h"

#define MyFloat float

#include <vector>
#include "MTkdtree.h"

using namespace std;

#define NTAB_MF 100
#define NDIM 3


extern FILE *fhabnd;
extern float FOF_WIDTH_BOUND;

typedef struct {
	int npart;
//	long long id;
//	double pot;
	float center[NDIM];
} fof_Group;

class Group {
public:
	int npart;
	float center[NDIM];
	vector<int> pidx;

	int nsubhalo;
	vector<int> subidx;
#ifdef HALOPROPERTY
//	float R200;
	float Mcrit200;
	float Mmean200;
//	float Mtophat;
#endif

	Group() {
		npart = 0;
		pidx.resize(0);
		nsubhalo = 0;
		pidx.resize(0);
	}
};

class GroupCatalog{
public:
	unsigned long group_number;	
	Group *group;
	GroupCatalog() {
		group_number = 0;
		group = nullptr;
	}
	~GroupCatalog() {
		delete[] group;
	}
};

typedef struct {
#ifdef SAVE_MEM
	float* pos;
	unsigned int pidx;
#else
	float pos[3];
	unsigned long pidx;
#endif
	int label_halo;
} grp_Body;

int compare_group(const void *h1, const void *h2);

template<typename Body_T>
class FOFHaloFinder {
private:
	int MAXBUFF_GROUP;
	int MAXBUFF_INDEX;
	fof_Group *grp;
	int* pidx_st;
	int* pidx_ed;
	int* pidx;

public:
	string data_path;
	int TASK_RANK;
	int PROC_RANK;
	int PROC_SIZE;
	
	grp_Body* part_intile;
	int NPART_intile;
	long long NPART_TOTAL;
	
	tNode* btree;
	int fnode, lnode;
	
	tPack* leaf;
	int fleaf, lleaf;
	int MAXinLEAF;
	
	float MASSPART;
	float BOXSIZE;
	double BOX_DOM_L[3];
	double BOX_DOM_H[3];
	
	int MINHALO_NPART;

	int nbnd;	
	int  direct_local_start;

	Group *Grps;
	
	double LINK_LENGTH;

	unsigned long num_group;

FOFHaloFinder(Body_T* prt, unsigned long npart, unsigned long long npart_tot, float boxsize, int nsidemesh, int proc_rank,
	 	int* MESH_L, int* MESH_H, int rank, int nbndmesh) {
	double t0 = dtime();
	NPART_TOTAL = npart_tot;
	BOXSIZE = boxsize;
	PROC_RANK = proc_rank;
	TASK_RANK = rank;
	nbnd = nbndmesh; 

	float dx = boxsize / nsidemesh; 
	for (int i=0;i<3;++i) BOX_DOM_L[i] = dx * (MESH_L[i] + nbnd);			
	for (int i=0;i<3;++i) BOX_DOM_H[i] = dx * (MESH_H[i] - nbnd);
	part_intile = (grp_Body*)mymalloc(sizeof(grp_Body) * npart, 200);
	LOG(3, "[%d] (%d) part_intile malloc (%lld, %lf GB) ", PROC_RANK, TASK_RANK, npart, (float)npart * sizeof(grp_Body) / 1.0e9);
	if (part_intile == NULL ) {
		LOG(0, "[%d] (%d) error part_intile malloc in FOFHaloFinder", PROC_RANK, TASK_RANK);
		exit(0);
	}
	int n = 0;
	for (unsigned long i=0;i<npart;++i) {
#ifndef INTXYZ
		if (prt[i].pos[0] >= dx * MESH_L[0] && prt[i].pos[0] < dx * MESH_H[0] 
			&& prt[i].pos[1] >= dx * MESH_L[1] && prt[i].pos[1] < dx * MESH_H[1] 
			&& prt[i].pos[2] >= dx * MESH_L[2] && prt[i].pos[2] < dx * MESH_H[2]) {

#ifdef SAVE_MEM
			part_intile[n].pos = prt[i].pos;
#else
			part_intile[n].pos[0] = prt[i].pos[0];
			part_intile[n].pos[1] = prt[i].pos[1];
			part_intile[n].pos[2] = prt[i].pos[2];
#endif
			part_intile[n].pidx = i;
			++n;
		}

#else
		float ppos[3];
		ppos[0] = (float) (prt[i].posi[0] * INT2POS );
		ppos[1] = (float) (prt[i].posi[1] * INT2POS );
		ppos[2] = (float) (prt[i].posi[2] * INT2POS );
		if (ppos[0] >= dx * MESH_L[0] && ppos[0] < dx * MESH_H[0] 
			&& ppos[1] >= dx * MESH_L[1] && ppos[1] < dx * MESH_H[1] 
			&& ppos[2] >= dx * MESH_L[2] && ppos[2] < dx * MESH_H[2]) {

			part_intile[n].pos[0] = ppos[0];
			part_intile[n].pos[1] = ppos[1];
			part_intile[n].pos[2] = ppos[2];
			part_intile[n].pidx = i;
			++n;
		}

#endif

	}
	NPART_intile = n;

	double t1 = dtime();
	LOG(3, "[%d] (%d) fetch local particles (%d) time %lf", PROC_RANK, TASK_RANK, NPART_intile, t1 - t0);
	LOG(3, "[%d] (%d) fetch domain (%lf %lf %lf)-(%lf %lf %lf)", PROC_RANK, TASK_RANK, dx * MESH_L[0], dx * MESH_L[1], dx * MESH_L[2], dx * MESH_H[0], dx * MESH_H[1], dx * MESH_H[2]);
		
	MAXinLEAF = 16;
	build_localtree();
	double t2 = dtime();
	LOG(3, "[%d] (%d) build localtree (n %d, l %d) time %lf", PROC_RANK, TASK_RANK, lnode - fnode, lleaf - fleaf, t2 - t1);
}

~FOFHaloFinder() {
	delete[] Grps;
}

void build_localtree() {
	MTKDtree<grp_Body, tNode, tPack> builder(part_intile, NPART_intile, MAXinLEAF);
	PartTree<grp_Body, tNode, tPack> tres = builder.build_kdtree();
	part_intile = tres.part;
	btree = tres.btree;
	fnode = tres.first_node;
	lnode = tres.last_node;
	leaf = tres.leaf;
	fleaf = tres.first_leaf;
	lleaf = tres.last_leaf;
}

void fof(double redshift, int min_npart_in_halo, double link_length, int rank_snap) {
	if ( 0 == PROC_RANK )
		LOG(3, "[%d] (%d) FOF Halo finder begin ( z = %lf )", PROC_RANK, TASK_RANK, redshift);

	fof_init(min_npart_in_halo, link_length);

	double tfof = dtime();
	fof_proc_new( rank_snap );
//	fof_proc( rank_snap );
	LOG(3, "[%d] (%d) fof halo finder time %lf", PROC_RANK, TASK_RANK, dtime() - tfof);

	fof_free();

	if ( 0 == PROC_RANK )
		LOG(3, "[%d] (%d) FOF Halo finder complete ", PROC_RANK, TASK_RANK);
}

inline void insert_id_group(int ipart, int target, int igroup) {
	if (target < fleaf || target >= fnode ) {
		LOG(0, "[%d] (%d) insert error !", PROC_RANK, TASK_RANK);
		exit(0);
	}

	int n, p, c, d;
	double dx[3], dr2;

	for (n=0; n<leaf[target].npart; n++) {
		p = leaf[target].ipart + n;		

		if (ipart == p)
			continue;

		if ( -1 != part_intile[p].label_halo )
			continue;

		dr2 = 0.0;
		for (d=0; d<NDIM; d++) {
			dx[d] = part_intile[p].pos[d] - part_intile[ipart].pos[d]; 
			dr2 += dx[d] * dx[d];
		}

		if ( dr2 < LINK_LENGTH * LINK_LENGTH  ) {
			// inser into pidx
			c = pidx_ed[igroup];
			pidx[c] = p;
			part_intile[p].label_halo = igroup;
			pidx_ed[igroup] ++;
		}
	}
}

void fopart(int direct, int iNode, int ipart, int igrp) {
	int new_direct = (direct + 1)%NDIM;
	int in, ou, in0, in1;
	int tp, d;
	double ds, dr, dr2;

	if ( iNode < fnode ) {
		insert_id_group(ipart, iNode, igrp);	
		return ;
	}

	in0 = btree[iNode].son[0];
	in1 = btree[iNode].son[1];

	ds = btree[iNode].split - part_intile[ipart].pos[direct]; 

	if ( ds >= 0.0 ) {
		in = in0;
		ou = in1;
	}
	else {
		in = in1;
		ou = in0;
	}


	if ( in >= fnode ) {
		fopart(new_direct, in, ipart, igrp);
	}
	else {
		insert_id_group(ipart, in, igrp);
	}

	if (  ou < fnode ) {
		insert_id_group(ipart, ou, igrp);	
	}
	else if ( fabs(ds) < LINK_LENGTH ) {
		fopart(new_direct, ou, ipart, igrp);
	}
}

int remove_bnd_particles(int grp_idx, int& npart_left) {

	double link_len = 1.01 * LINK_LENGTH;
	double wid_b = (double) FOF_WIDTH_BOUND;
	int close_bnd = 0;
	int m, n = grp_idx;
	double center_z[3] = {0, 0, 0};

	for (m=pidx_st[n]; m<=pidx_ed[n]; m++) {
		int r = pidx[m];

		double dx0 =  (double)part_intile[r].pos[0] - (BOX_DOM_L[0] - wid_b) ;
		double dx1 =  BOX_DOM_H[0] + wid_b - (double)part_intile[r].pos[0] ;

		double dy0 =  (double)part_intile[r].pos[1] - (BOX_DOM_L[1] - wid_b) ;
		double dy1 =  BOX_DOM_H[1] + wid_b - (double)part_intile[r].pos[1] ;

		double dz0 =  (double)part_intile[r].pos[2] - (BOX_DOM_L[2] - wid_b) ;
		double dz1 =  BOX_DOM_H[2] + wid_b - (double)part_intile[r].pos[2] ;

		if (dx0 < link_len || dx1 < link_len ) close_bnd = 1;
		if (dy0 < link_len || dy1 < link_len ) close_bnd = 1;
		if (dz0 < link_len || dz1 < link_len ) close_bnd = 1;

		if ( 1 == close_bnd) break;
	}

	if (close_bnd == 0)
		return 0;
	
	int n_zone = 0; 
	for (m=pidx_st[n]; m<=pidx_ed[n]; m++) {
		int r = pidx[m];

		double x_z =  (double)part_intile[r].pos[0];
		double y_z =  (double)part_intile[r].pos[1];
		double z_z =  (double)part_intile[r].pos[2];

		if (( BOX_DOM_L[0] - wid_b  <= x_z && x_z < BOX_DOM_L[0] + wid_b) ||  ( BOX_DOM_H[0] - wid_b  <= x_z && x_z < BOX_DOM_H[0] + wid_b)
     		|| ( BOX_DOM_L[1] - wid_b  <= y_z && y_z < BOX_DOM_L[1] + wid_b) ||  ( BOX_DOM_H[1] - wid_b  <= y_z && y_z < BOX_DOM_H[1] + wid_b)
		|| ( BOX_DOM_L[2] - wid_b  <= z_z && z_z < BOX_DOM_L[2] + wid_b) ||  ( BOX_DOM_H[2] - wid_b  <= z_z && z_z < BOX_DOM_H[2] + wid_b)){
			center_z[0] += x_z;	
			center_z[1] += y_z;	
			center_z[2] += z_z;	
			n_zone ++;	

		}
	}

	center_z[0] /= (double)n_zone;				
	center_z[1] /= (double)n_zone;				
	center_z[2] /= (double)n_zone;	
			
	if (BOX_DOM_L[0] <= center_z[0] && center_z[0] < BOX_DOM_H[0] 
	 && BOX_DOM_L[1] <= center_z[1] && center_z[1] < BOX_DOM_H[1] 
	 && BOX_DOM_L[2] <= center_z[2] && center_z[2] < BOX_DOM_H[2]) {
		return 0;
	} else {
		npart_left -= n_zone;
		return 1;
	}
}


int transit_bnd_particles(int grp_idx, int& npart_left) {

	double link_len = 1.01 * LINK_LENGTH;
	double wid_b = (double) FOF_WIDTH_BOUND;
	int close_bnd = 0;
	int m, n = grp_idx;
	double center_z[3] = {0, 0, 0};

	for (m=pidx_st[n]; m<=pidx_ed[n]; m++) {
		int r = pidx[m];

		double dx0 =  (double)part_intile[r].pos[0] - (BOX_DOM_L[0] - wid_b) ;
		double dx1 =  BOX_DOM_H[0] + wid_b - (double)part_intile[r].pos[0] ;

		double dy0 =  (double)part_intile[r].pos[1] - (BOX_DOM_L[1] - wid_b) ;
		double dy1 =  BOX_DOM_H[1] + wid_b - (double)part_intile[r].pos[1] ;

		double dz0 =  (double)part_intile[r].pos[2] - (BOX_DOM_L[2] - wid_b) ;
		double dz1 =  BOX_DOM_H[2] + wid_b - (double)part_intile[r].pos[2] ;

		if (dx0 < link_len || dx1 < link_len ) close_bnd = 1;
		if (dy0 < link_len || dy1 < link_len ) close_bnd = 1;
		if (dz0 < link_len || dz1 < link_len ) close_bnd = 1;

		if ( 1 == close_bnd) break;
	}

	if (close_bnd == 0)
		return 0;
	
	int n_zone = 0; 
	for (m=pidx_st[n]; m<=pidx_ed[n]; m++) {
		int r = pidx[m];

		double x_z =  (double)part_intile[r].pos[0];
		double y_z =  (double)part_intile[r].pos[1];
		double z_z =  (double)part_intile[r].pos[2];

		fprintf(fhabnd, "%lf %lf %lf\n", x_z, y_z, z_z);
		
	}

	return 0;
}

void fof_proc_new(int rank_snap) {
	int n, m;
	int igroup;

	double t0 = dtime();
	num_group = 0;
	int drt = direct_local_start;
	drt = 0;
	pidx_st[0] = 0;
	pidx_ed[0] = 0;

	int id;

	for (n=0; n<NPART_intile; n++) {
		if ( -1 ==  part_intile[n].label_halo   ) {

			igroup = num_group;

			part_intile[n].label_halo = igroup;

			pidx_ed[igroup] = pidx_st[igroup] + 1;
			pidx[pidx_st[igroup]] = n;

			fopart(drt, fnode, n, igroup);

			for (m=pidx_st[igroup]+1; m < pidx_ed[igroup]; m++) {
				id = pidx[m]; 	
				fopart(drt, fnode, id, igroup);
			}

			pidx_ed[igroup] -= 1;
			if (pidx_ed[igroup] > pidx_st[igroup]) {
				grp[igroup].npart = pidx_ed[igroup] - pidx_st[igroup] + 1;
				num_group++;
				pidx_st[igroup+1] = pidx_ed[igroup] + 1;
			}
			else {
				part_intile[n].label_halo = -1;
			}

		}

	}

	double t1 = dtime();
	LOG(3, "[%d] (%d) finding halos (NPART_intile %d) time %lf", PROC_RANK, TASK_RANK, NPART_intile, t1 - t0);

	Grps = new Group[num_group];
	int n_group_cnt = 0;
	double center[3];
	for (n=0; n<num_group; n++) {

		if (grp[n].npart >= MINHALO_NPART ) {

			center[0] = 0.0;
			center[1] = 0.0;
			center[2] = 0.0;

			for (m=pidx_st[n]; m<=pidx_ed[n]; m++) {
				int r = pidx[m];
				center[0] += (double)part_intile[r].pos[0];
				center[1] += (double)part_intile[r].pos[1];
				center[2] += (double)part_intile[r].pos[2];
			}

			center[0] /= (double)grp[n].npart;
			center[1] /= (double)grp[n].npart;
			center[2] /= (double)grp[n].npart;
			
			grp[n].center[0] = center[0];
			grp[n].center[1] = center[1];
			grp[n].center[2] = center[2];

			if (BOX_DOM_L[0] <= center[0] && center[0] < BOX_DOM_H[0] 
					&& BOX_DOM_L[1] <= center[1] && center[1] < BOX_DOM_H[1] 
					&& BOX_DOM_L[2] <= center[2] && center[2] < BOX_DOM_H[2]) 
			{

				double wid_b = (double) FOF_WIDTH_BOUND;
				int npart_left = grp[n].npart;

				if (transit_bnd_particles(n, npart_left)) {
					printf("[%d] remove bnd left %d BOUND %lf \n", PROC_RANK, npart_left, FOF_WIDTH_BOUND);
					if (npart_left < MINHALO_NPART)
						continue;

					double center_z[3] = {0, 0, 0};

					Grps[n_group_cnt].npart = npart_left;
					Grps[n_group_cnt].pidx.resize(npart_left);

					int ip = 0;
					for (m=pidx_st[n]; m<=pidx_ed[n]; m++) {
						int r = pidx[m];

						double x_z =  (double)part_intile[r].pos[0];
						double y_z =  (double)part_intile[r].pos[1];
						double z_z =  (double)part_intile[r].pos[2];

						if ( BOX_DOM_L[0] + wid_b  <= x_z && x_z < BOX_DOM_H[0] - wid_b 
					   	  && BOX_DOM_L[1] + wid_b  <= y_z && y_z < BOX_DOM_H[1] - wid_b
						  && BOX_DOM_L[2] + wid_b  <= z_z && z_z < BOX_DOM_H[2] - wid_b) {
							center_z[0] +=  x_z;
							center_z[1] +=  y_z;
							center_z[2] +=  z_z;
							Grps[n_group_cnt].pidx[ip++] = part_intile[r].pidx;
						}
					}
					center_z[0] /= (double)npart_left;
					center_z[1] /= (double)npart_left;
					center_z[2] /= (double)npart_left;

					Grps[n_group_cnt].center[0] = (float)center_z[0];
					Grps[n_group_cnt].center[1] = (float)center_z[1];
					Grps[n_group_cnt].center[2] = (float)center_z[2];

					if (Grps[n_group_cnt].pidx.size() != npart_left) {
						printf("[%d] grp pidx error.\n", PROC_RANK);
						exit(0);
					}	
					n_group_cnt ++;
					
				} 
				else 
				{

					Grps[n_group_cnt].npart = grp[n].npart;
					Grps[n_group_cnt].center[0] = (float)center[0];
					Grps[n_group_cnt].center[1] = (float)center[1];
					Grps[n_group_cnt].center[2] = (float)center[2];
					Grps[n_group_cnt].pidx.resize(grp[n].npart);
					int ip = 0;
					for (m=pidx_st[n]; m<=pidx_ed[n]; m++) {
						Grps[n_group_cnt].pidx[ip++] = part_intile[pidx[m]].pidx;
					}
					n_group_cnt ++;
				}

			}


		}
	}


	double t2 = dtime();

	num_group = n_group_cnt ;
	double t3 = dtime();
	LOG(3, "[%d] (%d) group(halo) sort time %lf", PROC_RANK, TASK_RANK, t3 - t2);
}


void fof_proc(int rank_snap) {
	int n, m;
	int igroup;

	double t0 = dtime();
	num_group = 0;
	int drt = direct_local_start;
	drt = 0;
	pidx_st[0] = 0;
	pidx_ed[0] = 0;

	int id;
//	long long nstart = 0; 

	for (n=0; n<NPART_intile; n++) {
		if ( -1 ==  part_intile[n].label_halo   ) {
//			printf("part n %d\n", n);

			igroup = num_group;

			part_intile[n].label_halo = igroup;

//			grp[igroup].id = nstart + (long long)n;

			pidx_ed[igroup] = pidx_st[igroup] + 1;
			pidx[pidx_st[igroup]] = n;

			fopart(drt, fnode, n, igroup);

			for (m=pidx_st[igroup]+1; m < pidx_ed[igroup]; m++) {
				id = pidx[m]; 	
//				printf("friend part n %d\n", id);
				fopart(drt, fnode, id, igroup);
			}

			pidx_ed[igroup] -= 1;
			if (pidx_ed[igroup] > pidx_st[igroup]) {
				grp[igroup].npart = pidx_ed[igroup] - pidx_st[igroup] + 1;
				num_group++;
				pidx_st[igroup+1] = pidx_ed[igroup] + 1;
			}
			else {
				part_intile[n].label_halo = -1;
			}

		}

	}

	double t1 = dtime();
	LOG(3, "[%d] (%d) finding halos (NPART_intile %d) time %lf", PROC_RANK, TASK_RANK, NPART_intile, t1 - t0);

	Grps = new Group[num_group];
	int n_group_cnt = 0;
	double center[3];
	for (n=0; n<num_group; n++) {

		if (grp[n].npart >= MINHALO_NPART ) {

			center[0] = 0.0;
			center[1] = 0.0;
			center[2] = 0.0;

			for (m=pidx_st[n]; m<=pidx_ed[n]; m++) {
				int r = pidx[m];
				center[0] += (double)part_intile[r].pos[0];
				center[1] += (double)part_intile[r].pos[1];
				center[2] += (double)part_intile[r].pos[2];
			}

			center[0] /= (double)grp[n].npart;
			center[1] /= (double)grp[n].npart;
			center[2] /= (double)grp[n].npart;
			
			grp[n].center[0] = center[0];
			grp[n].center[1] = center[1];
			grp[n].center[2] = center[2];

			if (BOX_DOM_L[0] <= center[0] && center[0] < BOX_DOM_H[0] 
			 && BOX_DOM_L[1] <= center[1] && center[1] < BOX_DOM_H[1] 
			 && BOX_DOM_L[2] <= center[2] && center[2] < BOX_DOM_H[2]) 
			{

				double wid_b = (double) FOF_WIDTH_BOUND;
				int npart_left = grp[n].npart;
#ifdef FOF_BND
				if (remove_bnd_particles(n, npart_left)) {
					printf("[%d] remove bnd left %d BOUND %lf \n", PROC_RANK, npart_left, FOF_WIDTH_BOUND);
					if (npart_left < MINHALO_NPART)
						continue;

					double center_z[3] = {0, 0, 0};

					Grps[n_group_cnt].npart = npart_left;
					Grps[n_group_cnt].pidx.resize(npart_left);

					int ip = 0;
					for (m=pidx_st[n]; m<=pidx_ed[n]; m++) {
						int r = pidx[m];

						double x_z =  (double)part_intile[r].pos[0];
						double y_z =  (double)part_intile[r].pos[1];
						double z_z =  (double)part_intile[r].pos[2];

						if ( BOX_DOM_L[0] + wid_b  <= x_z && x_z < BOX_DOM_H[0] - wid_b 
					   	  && BOX_DOM_L[1] + wid_b  <= y_z && y_z < BOX_DOM_H[1] - wid_b
						  && BOX_DOM_L[2] + wid_b  <= z_z && z_z < BOX_DOM_H[2] - wid_b) {
							center_z[0] +=  x_z;
							center_z[1] +=  y_z;
							center_z[2] +=  z_z;
							Grps[n_group_cnt].pidx[ip++] = part_intile[r].pidx;
						}
					}
					center_z[0] /= (double)npart_left;
					center_z[1] /= (double)npart_left;
					center_z[2] /= (double)npart_left;

					Grps[n_group_cnt].center[0] = (float)center_z[0];
					Grps[n_group_cnt].center[1] = (float)center_z[1];
					Grps[n_group_cnt].center[2] = (float)center_z[2];

					if (Grps[n_group_cnt].pidx.size() != npart_left) {
						printf("[%d] grp pidx error.\n", PROC_RANK);
						exit(0);
					}	
					n_group_cnt ++;
					
				} else
#endif 
				{

					Grps[n_group_cnt].npart = grp[n].npart;
					Grps[n_group_cnt].center[0] = (float)center[0];
					Grps[n_group_cnt].center[1] = (float)center[1];
					Grps[n_group_cnt].center[2] = (float)center[2];
					Grps[n_group_cnt].pidx.resize(grp[n].npart);
					int ip = 0;
					for (m=pidx_st[n]; m<=pidx_ed[n]; m++) {
						Grps[n_group_cnt].pidx[ip++] = part_intile[pidx[m]].pidx;
					}
					n_group_cnt ++;
				}

			}


		}
	}


	double t2 = dtime();

	num_group = n_group_cnt ;
	double t3 = dtime();
	LOG(3, "[%d] (%d) group(halo) sort time %lf", PROC_RANK, TASK_RANK, t3 - t2);
}

void fof_init( int min_npart_inhalo, double link_length) {
	int n;
	double mean_sep = BOXSIZE/pow((double)NPART_TOTAL, 0.333333333);
	LINK_LENGTH = link_length * mean_sep;
	if ( 0 == PROC_RANK && 0 == TASK_RANK)
		LOG(2, "[%d] (%d) FOF LINK_LENGTH %lf (BS %lf NP_TOT %lld link_length %lf)\n", PROC_RANK, TASK_RANK, LINK_LENGTH, BOXSIZE, NPART_TOTAL, 0.2);
	MAXBUFF_GROUP = NPART_intile/2;
	MAXBUFF_INDEX = NPART_intile;
	MINHALO_NPART = min_npart_inhalo;


	grp = (fof_Group*)mymalloc(sizeof(fof_Group)*MAXBUFF_GROUP, 201);
	if (grp == NULL ) {
		LOG(0, "[%d] (%d) error grp malloc in FOFHaloFinder", PROC_RANK, TASK_RANK);
		exit(0);
	}
	LOG(3, "[%d] (%d) grp malloc (%lld, %lf GB) ", PROC_RANK, TASK_RANK, MAXBUFF_GROUP, (float)MAXBUFF_GROUP * sizeof(fof_Group) / 1.0e9);
	pidx = (int*)mymalloc(sizeof(int)*MAXBUFF_INDEX, 202);
	LOG(3, "[%d] (%d) pidx malloc (%lld, %lf GB) \n", PROC_RANK, TASK_RANK, MAXBUFF_INDEX, (float)MAXBUFF_INDEX * sizeof(int) / 1.0e9);
	pidx_st = (int*)mymalloc(sizeof(int)*MAXBUFF_GROUP, 203);
	pidx_ed = (int*)mymalloc(sizeof(int)*MAXBUFF_GROUP, 204);
	LOG(3, "[%d] (%d) pidx_st/ed malloc (%lld, %lf GB) \n", PROC_RANK, TASK_RANK, MAXBUFF_GROUP, (float)MAXBUFF_GROUP * sizeof(int) / 1.0e9);

	for (n=0; n<MAXBUFF_GROUP; n++){
//		grp[n].id = -1LL;
		grp[n].npart = 0;
//		grp[n].pot = 0.0;

		grp[n].center[0] = 0.0;
		grp[n].center[1] = 0.0;
		grp[n].center[2] = 0.0;

	}

	for (n=0; n<NPART_intile; n++)
		part_intile[n].label_halo = -1;
}

void fof_free() {
	if (part_intile)
		myfree(part_intile, 200);
	if (grp)
		myfree(grp, 201);
	if (pidx)
		myfree(pidx, 202);
	if (pidx_st)
		myfree(pidx_st, 203);
	if (pidx_ed)
		myfree(pidx_ed, 204);
	leaf += fleaf;
	if (leaf)
		myfree(leaf, 250);
	btree += fnode;
	if (btree)
		myfree(btree, 251);
}

};

void mesh_decomp(int rank, int* npxyz, int nbnd, 
		int* MESH_L, int* MESH_H, int* mesh_l, int* mesh_h);

void halo_mass_func(Group* grps, unsigned long num_group, unsigned long* mcnt_tot, double* mrange, double mass_min, double scale, float MASSPART);

#endif
