#ifndef SUBFIND_H
#define SUBFIND_H

#define MyFloat float

//#define SAVE_MEM_EX

#ifdef HALO_FOF
#include "fof.h"
#endif
#include "MTkdtree.h"

using namespace std;

#define HIGHBIT (1<<30)

#ifdef ONEDIMS
#define NUMDIMS 1
#define KERNEL_COEFF_1 (4.0 / 3)
#define KERNEL_COEFF_2 (8.0)
#define KERNEL_COEFF_3 (24.0)
#define KERNEL_COEFF_4 (16.0)
#define KERNEL_COEFF_5 (8.0 / 3)
#define KERNEL_COEFF_6 (-8.0)
#define NORM_COEFF 2.0
#else
#ifndef TWODIMS
#define NUMDIMS 3                     /**< For 3D-normalized kernel */
#define KERNEL_COEFF_1 2.546479089470 /**< Coefficients for SPH spline kernel and its derivative */
#define KERNEL_COEFF_2 15.278874536822
#define KERNEL_COEFF_3 45.836623610466
#define KERNEL_COEFF_4 30.557749073644
#define KERNEL_COEFF_5 5.092958178941
#define KERNEL_COEFF_6 (-15.278874536822)
#define NORM_COEFF 4.188790204786 /**< Coefficient for kernel normalization. Note:  4.0/3 * PI = 4.188790204786 */
#else
#define NUMDIMS 2                                 /**< For 2D-normalized kernel */
#define KERNEL_COEFF_1 (5.0 / 7 * 2.546479089470) /**< Coefficients for SPH spline kernel and its derivative */
#define KERNEL_COEFF_2 (5.0 / 7 * 15.278874536822)
#define KERNEL_COEFF_3 (5.0 / 7 * 45.836623610466)
#define KERNEL_COEFF_4 (5.0 / 7 * 30.557749073644)
#define KERNEL_COEFF_5 (5.0 / 7 * 5.092958178941)
#define KERNEL_COEFF_6 (5.0 / 7 * (-15.278874536822))
#define NORM_COEFF M_PI /**< Coefficient for kernel normalization. */
#endif
#endif /* ONEDIMS */

#define SOFTFAC1 (32.0 / 3) /**< Coefficients for gravitational softening */
#define SOFTFAC2 32.0
#define SOFTFAC3 (-38.4)
#define SOFTFAC4 (-2.8)
#define SOFTFAC5 (16.0 / 3)
#define SOFTFAC6 6.4
#define SOFTFAC7 (-9.6)
#define SOFTFAC8 (64.0 / 3)
#define SOFTFAC9 (-48.0)
#define SOFTFAC10 38.4
#define SOFTFAC11 (-32.0 / 3)
#define SOFTFAC12 (-1.0 / 15)
#define SOFTFAC13 (-3.2)
#define SOFTFAC14 (1.0 / 15)
#define SOFTFAC15 (-16.0)
#define SOFTFAC16 9.6
#define SOFTFAC17 (-64.0 / 30)
#define SOFTFAC18 128.0
#define SOFTFAC19 (-115.2)
#define SOFTFAC20 (64.0 / 3)
#define SOFTFAC21 (-96.0)
#define SOFTFAC22 115.2
#define SOFTFAC23 (-128.0 / 3)
#define SOFTFAC24 (4.0 / 30)

#define SOFTFAC30 (32.0 / 3)
#define SOFTFAC31 (-576.0 / 5)
#define SOFTFAC32 (128.0)
#define SOFTFAC33 (-1152.0 / 5)
#define SOFTFAC34 (384.0)
#define SOFTFAC35 (2.0 * 384.0)

#define SOFTFAC40 (64.0 / 3)
#define SOFTFAC41 (2.0 / 15)
#define SOFTFAC42 (-96.0)
#define SOFTFAC43 (576.0 / 5)
#define SOFTFAC44 (-128.0 / 3)
#define SOFTFAC45 (-96.0)
#define SOFTFAC46 (-2.0 / 5)
#define SOFTFAC47 (1152.0 / 5)
#define SOFTFAC48 (-128.0)
#define SOFTFAC49 (8.0 / 5)
#define SOFTFAC50 (-256.0)
#define SOFTFAC51 (-8.0)

#define SQRT_PI 1.772453850906           /* sqrt(M_PI) */
#define FACT1 0.366025403785             /* FACT1 = 0.5 * (sqrt(3)-1) */
#define FACTSQRT3HALF 0.866025403785     /* sqrt(3)/2 */
#define FACTSQRT3 (2.0 * 0.866025403785) /* sqrt(3) */

// map gobal idx to in_grp idx
extern int *idxmap;

class subGroup {
public:
	int Len;
	MyFloat pos[3];
	MyFloat CM[3];
	MyFloat vel[3];
	MyFloat spin[3];
	MyFloat Mass;
	MyFloat SubVmax;
	MyFloat SubVmaxRad;
	MyFloat SubHalfMassRad;
	MyFloat SubMostBoundID;
	int SubRankInGr;
	int SubParentRank;

	int in_group;
	vector<int> pidx;
	subGroup () {
		Len = 0;
		pidx.resize(0);
	}
	~subGroup () {
		Len = 0;
		pidx.resize(0);
	}
};

class subGroupCatalog{
public:
	unsigned long subgrp_number;	
	vector<subGroup> subgrps;
	subGroupCatalog() {
		subgrp_number = 0;
		subgrps.resize(0);
	}
	~subGroupCatalog() {
		subgrp_number = 0;
		subgrps.resize(0);
	}
};

typedef struct subfind_part {
	unsigned long pidx;
#ifdef SAVE_MEM_EX
	float* pos;
#else
	float pos[3];
#endif
} SFPart;

typedef struct subfind_part_info {
	unsigned long pidx; // global
	int in_group;
 
	MyFloat DM_Hsml;
	union {
		MyFloat DM_Density;
		MyFloat DM_Potential;
	};
	MyFloat DM_BindingEnergy;
	int invidx;
} SFInfo;

typedef struct subfind_partinfo {
	int submark;
	int origintask;
} SFMark;

typedef struct subfind_density_part {
	unsigned long pidx;
	float pos[3];
	MyFloat DM_Hsml;
	MyFloat DM_Density;
	int nrst2_count;
	int nrst2_index[2];
	MyFloat nrst2_dist[2];	
} SFDPart;

typedef struct subfind_desity_sort {
	int locidx;
	int pidx;
	MyFloat DM_Density;
} SFDSort;

typedef struct subfind_linklist {
	int head;
	int len;
	int next;
	int tail;
//	int prevlen;
} SFLinkList;

typedef struct subfind_coll_cand {
	int head;
	int len;
	int rank;
	int parent;
	int nsub;
	int subnr;
	int bound_length;
} SFCandidate;

typedef struct subfind_submark_sort {
	int pidx;
	int submark;
} SFSMSort;

typedef struct gfactors {
	MyFloat fac0;
} GFac;

typedef struct sort_r2list {
	MyFloat r;
	MyFloat mass;
} R2ListSort;

template<typename Body_T>
class SIMInfo {
public:
	Body_T* part;
	unsigned long npart;
	long long NPART_TOTAL;
	int N_SIDEMESH;
	float MASSPART;
	float BOXSIZE;
	float Hubble0;
	float OmegaM0;
	float OmegaX0;
	float redshift;
	int fof_group_min_len;
	float fof_linklength;

	int PROC_RANK;
	int PROC_SIZE;
}; 

//#define DEBUG
#define DBGN 200000

#define MAX_DOUBLE_NUMBER 1e306
#define MAX_ITER_UNBIND 500

#define DesNumNgb 64
#define DesLinkNgb 20
#define MaxNumNgbDeviation 1

#define RECOMPUTE_ALL 0
#define UPDATE_ALL 1

#define ErrTolTheta 0.66
#define GravConst 43007.105732
#define Hubble 0.1

#define MAX_UNBOUND_FRAC_BEFORE_BULK_VELOCITY_UPDATE 0.02
#define MAX_UNBOUND_FRAC_BEFORE_POTENTIAL_UPDATE 0.20

int compare_density(const void *h1, const void *h2);
int compare_rank(const void *h1, const void *h2);
int compare_subnr(const void *h1, const void *h2);
int compare_nsubs(const void *h1, const void *h2);
int compare_PPS(const void *h1, const void *h2);
int compare_binding_energy(const void *h1, const void *h2);
int compare_bound_length(const void *h1, const void *h2);
int compare_len(const void *h1, const void *h2);
int compare_rotcurve(const void *h1, const void *h2);
int compare_subgroup(const void *h1, const void *h2);
void subhalo_mass_func(subGroup* subgrps, unsigned long num_subgrp, unsigned long* mcnt_tot, double* mrange, double mass_min, double scale, float MASSPART);

template <typename T>
bool sampling_print(T* pos) {
	return  (1 || pos[2] >= 6700 && pos[2] < 6740);
}

template <typename T>
void fwrite_pos(T* buf, int num, const char* fn) {
	FILE *fd;
	fd = fopen(fn, "w"); 
	for (int i=0; i<num; ++i)
		if (sampling_print(buf[i].pos))
			fprintf(fd, "%lf %lf %lf\n", buf[i].pos[0], buf[i].pos[1], buf[i].pos[2]); 
	fclose(fd);
}

template<typename Body_T>
class SubFind {
private:
	Body_T* part_inproc;
	unsigned long npart_inproc;
	int N_SIDEMESH;
	int PROC_RANK;
	int PROC_SIZE;
	long long NPART_TOTAL; 
	float MASSPART; 
	float BOXSIZE; 
	float Hubble0; 
	float OmegaM0; 
	float OmegaX0; 
	float redshift;
	float velFac;
	float ForceSoftening;

void subfind_find_nearesttwo(SFDPart* part, int npart, MyFloat LINKLENGTH) {
	MTKDtree<SFDPart, tNode, tPack> builder(part, npart, 16);
	PartTree<SFDPart, tNode, tPack> tree = builder.build_kdtree();
	KNN3DIndex<SFDPart, tNode, tPack> knnidx(&tree);
	MyFloat mean_sep = BOXSIZE/pow((MyFloat)NPART_TOTAL, 0.333333333);
	double estimate_h = DesLinkNgb * LINKLENGTH * mean_sep;
 
	int idx;
	MyFloat dist;
	for (unsigned long i=0; i<npart; ++i) {
		PartIndex* sidx;
		MyFloat maxh = min(part[i].DM_Hsml * 2., estimate_h);
		int Nngb = DesLinkNgb - 1; //exclude itself 
		sidx = knnidx.KNN_search(&part[i], 1, Nngb, maxh);
		while (sidx->matched_num < Nngb) {
			if (maxh >= estimate_h) {
				printf("[%d] part %lu matched_num %d in the radius of %lf (should not be less than %d)\n", PROC_RANK, i, sidx->matched_num, maxh, Nngb);
				exit(0);
			}
			maxh = min(maxh * 1.5, estimate_h*1.00001);
			sidx = knnidx.KNN_search(&part[i], 1, Nngb, maxh);
		}
		part[i].DM_Hsml = sidx->distance_list[Nngb-1];
		part[i].nrst2_count = 0;
		for (int j=0; j<Nngb; ++j) {
			idx = sidx->matched_list[j];
			if (part[idx].DM_Density > part[i].DM_Density) {
				dist = sidx->distance_list[j];
				if (part[i].nrst2_count < 2) {
					part[i].nrst2_index[part[i].nrst2_count] = idx;
					part[i].nrst2_dist[part[i].nrst2_count] = dist;
					part[i].nrst2_count++;
				} else {
					int k = part[i].nrst2_dist[0] > part[i].nrst2_dist[1] ? 0 : 1;
					if (dist < part[i].nrst2_dist[k]) {
						part[i].nrst2_dist[k] = dist;
						part[i].nrst2_index[k] = idx;
					}
				}
			}
		}
		delete[] sidx;
	}

	builder.MTKDtree_free();

}

int subfind_find_coll_candidates(SFDPart* part, int npart, SFCandidate* coll_candidates, int max_coll_candidates, SFLinkList* linklists) {
	SFDSort* psort = (SFDSort*)mymalloc(sizeof(SFDSort) * npart, 220);
	if (psort == NULL ) {
		LOG(0, "[%d] error psort malloc in subfind.cc", PROC_RANK);
		exit(0);
	}
	for (unsigned long i=0; i<npart; ++i) {
		psort[i].locidx = i; //idx in group 
		psort[i].pidx = part[i].pidx;
		psort[i].DM_Density = part[i].DM_Density;
	}
	
	// sort by the key of density
	qsort(psort, npart, sizeof(SFDSort), compare_density);
		
	for (unsigned long i=0; i<npart; ++i) {
		linklists[i].head = -1;
		linklists[i].next = -1;
		linklists[i].tail = -1;
		linklists[i].len = 0;
//		linklists[i].prevlen = 0;
	}

	int count_cand = 0;

	int old_tail, head, head_attach, tail, tail_attach, len, len_attach;
	for (unsigned long ii=0; ii<npart; ++ii) {
		int i = psort[ii].locidx;
		int cnt = part[i].nrst2_count;
		int idx1 = part[i].nrst2_index[0];
		int idx2 = part[i].nrst2_index[1];
		switch(cnt) {
			case 0: // a lonely maximum -> new group
				linklists[i].head = i;	
				linklists[i].tail = i;	
				linklists[i].len = 1;
				break;
			case 1: // the particle is attached to one group
				head = linklists[idx1].head;
				if (head == -1) {
					printf("[%d] head is not set, not possible!\n", PROC_RANK);
					exit(0);
				}
				old_tail = linklists[head].tail;
				linklists[head].tail = i;
				linklists[head].len += 1;
				linklists[i].head = head;
				linklists[old_tail].next = i;	
				break;
			case 2: // merge the two groups
				head = linklists[idx1].head;
				head_attach = linklists[idx2].head;
				if (head == -1 || head_attach == -1) {
					printf("[%d] any of heads is not set, not possible!\n", PROC_RANK);
					exit(0);
				}
				if (head != head_attach) {
					bool swap = false;
					len = linklists[head].len;
					tail = linklists[head].tail;
					len_attach = linklists[head_attach].len;
					tail_attach = linklists[head_attach].tail;
					if (len_attach > len || (len_attach == len && head_attach < head))
						swap = true;
					if (swap) {
						int tmp = head;
						head = head_attach;
						head_attach = tmp;

						tmp = tail;
						tail = tail_attach;
						tail_attach = tmp;

						tmp = len;
						len = len_attach;
						len_attach = tmp;
					}
					if (len_attach >= DesLinkNgb && len >= DesLinkNgb) {
						if (count_cand < max_coll_candidates) {
							coll_candidates[count_cand].len = len_attach;
							coll_candidates[count_cand].head = head_attach;
							count_cand++;
						}
					}
					// merge two group
					linklists[head].tail = tail_attach;
					linklists[head].len += len_attach;
					linklists[tail].next = head_attach;
					int ss = head_attach;
					do {
						linklists[ss].head = head;
						ss = linklists[ss].next;
					} while (ss >= 0); 	
				}
				// add this part
				old_tail = linklists[head].tail;
				linklists[head].tail = i;
				linklists[head].len += 1;
				linklists[i].head = head;
				linklists[old_tail].next = i;
				break;
		}
	}

	head = -1;
	int prev = -1, next;
	for (unsigned long i=0; i<npart; ++i) {
		if (linklists[i].head == i) {
			tail = linklists[i].tail;
			len = linklists[i].len;
			next = linklists[tail].next;
			if (next == -1) {
				if (prev < 0)
					head = i;
				if (prev >= 0)
					linklists[prev].next = i;
				prev = tail;
			}
		}
	}
	if (count_cand < max_coll_candidates) {
		coll_candidates[count_cand].len = npart;
		coll_candidates[count_cand].head = head;
		count_cand++;
	}

	int p = head;
	int rank = 0;
	while(p >= 0) {
		linklists[p].len = rank;
		rank++;
		p = linklists[p].next;
	}
	for (int k=0; k<count_cand; ++k) {
		coll_candidates[k].rank = linklists[coll_candidates[k].head].len;
	}

//	LOG(1, "this group (%lu) has %d subhalo candidates.", npart, count_cand);
//	for(int i=0; i<min(10, count_cand); ++i) {
//		LOG(2, "cand %d head %d len %d rank %d.", i, coll_candidates[i].head, coll_candidates[i].len, coll_candidates[i].rank);
//	}

	myfree(psort, 220);

	return count_cand;
}

void get_gfactors_potential(GFac &res, const MyFloat r, const MyFloat hmax, const MyFloat rinv) {
	if (r >= hmax) {
		res.fac0 = rinv;
	} else {
		MyFloat h_inv = 1 / hmax;
		MyFloat u = r * h_inv;
		if (u < 0.5) {
			MyFloat u2 = u * u;
			res.fac0 = -h_inv * (SOFTFAC4 + u2 * (SOFTFAC5 + u2 * (SOFTFAC6 * u + SOFTFAC7)));
		} else {
			MyFloat u2 = u * u;
			res.fac0 = -h_inv * (SOFTFAC13 + SOFTFAC14 / u + u2 * (SOFTFAC1 + u * (SOFTFAC15 + u *(SOFTFAC16 + SOFTFAC17 * u))));
		}
	}
}

void walk_tree(PartTree<SFPart, tNode, tPack>& tree, int inode, float* pos, MyFloat& potential) {
	int first_leaf = tree.first_leaf;
	int last_leaf = tree.last_leaf;
	int first_node = tree.first_node;
	SFPart* part = tree.part;
	tPack* leaf = tree.leaf;
	tNode* btree = tree.btree;

	if (inode == -1)
		return;
	if (inode < first_leaf) {
		printf("error\n");
		exit(0);
	}
	MyFloat dx[3];
	MyFloat hmax = ForceSoftening;
	MyFloat theta2 = ErrTolTheta * ErrTolTheta;

	// this node is leaf
	if (inode < last_leaf) {
		for (int n=0; n<leaf[inode].npart; n++) {
			int ip = leaf[inode].ipart + n;
			dx[0] = pos[0] - part[ip].pos[0];
			dx[1] = pos[1] - part[ip].pos[1];
			dx[2] = pos[2] - part[ip].pos[2];
			MyFloat r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
			MyFloat r = sqrt(r2);
			MyFloat rinv = (r > 0) ? 1.0 / r : 0;
			GFac gfac;
			get_gfactors_potential(gfac, r, hmax, rinv);
			potential -= gfac.fac0; 
			
		}
		return;
	}
	// this node is internal node
	if (inode >= first_node) {
		dx[0] = pos[0] - btree[inode].center[0];
		dx[1] = pos[1] - btree[inode].center[1];
		dx[2] = pos[2] - btree[inode].center[2];
		MyFloat r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
		MyFloat len = max(max(btree[inode].width[0], btree[inode].width[1]), btree[inode].width[2]);
		MyFloat len2 = len * len;
		if (r2 * theta2 >= len2) {
			MyFloat r = sqrt(r2);
			MyFloat rinv = (r > 0) ? 1.0 / r : 0;
			GFac gfac;
			get_gfactors_potential(gfac, r, hmax, rinv);
			potential -= gfac.fac0 * btree[inode].npart; 
		} else {
			// open the node
			walk_tree(tree, btree[inode].son[0], pos, potential);
			walk_tree(tree, btree[inode].son[1], pos, potential);
		}	
		return;
	}
	
}

void subfind_potential_compute(Body_T* prt, SFInfo* sfprt, int* d, int num, PartTree<SFPart, tNode, tPack>& tree) {
	int nth = min(omp_get_max_threads(), (num+49999)/50000);//omp_get_max_threads();
#pragma omp parallel for num_threads(nth) schedule(dynamic)
	for (int i=0; i<num; ++i) {
		MyFloat potential = 0.;
		MyFloat ppos[3];

#ifndef INTXYZ
		ppos[0] = prt[d[i]].pos[0];
		ppos[1] = prt[d[i]].pos[1];
		ppos[2] = prt[d[i]].pos[2];
#else

		ppos[0] = prt[d[i]].posi[0] * INT2POS;
		ppos[1] = prt[d[i]].posi[1] * INT2POS;
		ppos[2] = prt[d[i]].posi[2] * INT2POS;

#endif	

	
		walk_tree(tree, tree.first_node, ppos, potential);
		sfprt[idxmap[d[i]]].DM_Potential = potential * MASSPART;
	}
}

MyFloat hubble_function(MyFloat a) {
	MyFloat Omega0 = OmegaM0;
	MyFloat OmegaLambda = OmegaX0;
	MyFloat hubble_a = Omega0 / (a * a * a) + (1 - Omega0 - OmegaLambda) / (a * a) + OmegaLambda;
	hubble_a = Hubble * sqrt(hubble_a);
	return hubble_a;
}

void subfind_get_factors(MyFloat& fac_vel_to_phys, MyFloat& fac_hubbleflow, MyFloat& fac_comov_to_phys) {
	MyFloat Time = 1.0 / (1.0 + redshift); 
	fac_vel_to_phys = 1.0 / Time;
	fac_hubbleflow = hubble_function(Time);
	fac_comov_to_phys = Time;
}

int subfind_unbind (Body_T * prt, SFInfo * sfprt, int *d, int num, Group *host_group)
{
	int phaseflag = RECOMPUTE_ALL;
	int iter = 0;
	int count_unbind = 0;
	int num_removed = 0;
	long long removed;
	int *dremoved = (int *) mymalloc (sizeof (int) * num, 221);
	MyFloat *potold = (MyFloat *) mymalloc (sizeof (MyFloat) * num, 222);

	MyFloat fac_vel_to_phys, fac_comov_to_phys, fac_hubbleflow;
	subfind_get_factors (fac_vel_to_phys, fac_hubbleflow, fac_comov_to_phys);

	SFPart *part;
	PartTree < SFPart, tNode, tPack > tree;
	MTKDtree < SFPart, tNode, tPack > *builder;
	MyFloat left[3], right[3];

	count_unbind = (int) (10.0*( 1.0-log((double)num)/10.5) ) ;
	if (count_unbind < 0) count_unbind = 0;
//	printf(" qwang num = %d %d\n", num, count_unbind);
	
	do
	{
		if (phaseflag == RECOMPUTE_ALL)
		{
			part = (SFPart *) mymalloc (sizeof (SFPart) * num, 223);
			for (int i = 0; i < num; ++i)
			{
#ifdef SAVE_MEM_EX
				part[i].pos = prt[d[i]].pos;
#else

#ifndef INTXYZ
				part[i].pos[0] = prt[d[i]].pos[0];
				part[i].pos[1] = prt[d[i]].pos[1];
				part[i].pos[2] = prt[d[i]].pos[2];
#else

				part[i].pos[0] = prt[d[i]].posi[0] * INT2POS;
				part[i].pos[1] = prt[d[i]].posi[1] * INT2POS;
				part[i].pos[2] = prt[d[i]].posi[2] * INT2POS;

#endif

#endif
				part[i].pidx = d[i];
			}
			builder = new MTKDtree < SFPart, tNode, tPack > (part, num, 16);
			tree = builder->build_kdtree ();
			left[0] = BOXSIZE;
			left[1] = BOXSIZE;
			left[2] = BOXSIZE;
			right[0] = 0.0;
			right[1] = 0.0;
			right[2] = 0.0;
			for (int i = 0; i < num; ++i)
			{
				if (part[i].pos[0] < left[0])
					left[0] = part[i].pos[0];
				if (part[i].pos[1] < left[1])
					left[1] = part[i].pos[1];
				if (part[i].pos[2] < left[2])
					left[2] = part[i].pos[2];
				if (part[i].pos[0] > right[0])
					right[0] = part[i].pos[0];
				if (part[i].pos[1] > right[1])
					right[1] = part[i].pos[1];
				if (part[i].pos[2] > right[2])
					right[2] = part[i].pos[2];
			}
			builder->center_kdtree (0, tree.first_node, left, right);
		}
		else
		{
			part = (SFPart *) mymalloc (sizeof (SFPart) * num_removed, 223);
			for (int i = 0; i < num_removed; ++i)
			{
#ifdef SAVE_MEM_EX
				part[i].pos = prt[dremoved[i]].pos;
#else

#ifndef INTXYZ
				part[i].pos[0] = prt[dremoved[i]].pos[0];
				part[i].pos[1] = prt[dremoved[i]].pos[1];
				part[i].pos[2] = prt[dremoved[i]].pos[2];
#else

				part[i].pos[0] = prt[dremoved[i]].posi[0]*INT2POS ;
				part[i].pos[1] = prt[dremoved[i]].posi[1]*INT2POS ;
				part[i].pos[2] = prt[dremoved[i]].posi[2] *INT2POS ;
#endif

#endif
				part[i].pidx = dremoved[i];
			}
			builder =new MTKDtree < SFPart, tNode, tPack > (part, num_removed, 16);
			tree = builder->build_kdtree ();
			left[0] = BOXSIZE;
			left[1] = BOXSIZE;
			left[2] = BOXSIZE;
			right[0] = 0.0;
			right[1] = 0.0;
			right[2] = 0.0;
			for (int i = 0; i < num_removed; ++i)
			{
				if (part[i].pos[0] < left[0])
					left[0] = part[i].pos[0];
				if (part[i].pos[1] < left[1])
					left[1] = part[i].pos[1];
				if (part[i].pos[2] < left[2])
					left[2] = part[i].pos[2];
				if (part[i].pos[0] > right[0])
					right[0] = part[i].pos[0];
				if (part[i].pos[1] > right[1])
					right[1] = part[i].pos[1];
				if (part[i].pos[2] > right[2])
					right[2] = part[i].pos[2];
			}
			builder->center_kdtree (0, tree.first_node, left, right);
			for (int i = 0; i < num; ++i)
				potold[i] = sfprt[idxmap[d[i]]].DM_Potential;
		}
		if (part == NULL)
		{
			LOG (0, "error part malloc in subfind.cc");
			exit (0);
		}

		double ptime = dtime ();
		subfind_potential_compute (prt, sfprt, d, num, tree);
		//      pot_t += dtime() - ptime;
		builder->MTKDtree_free ();
		delete builder;

		if (phaseflag == RECOMPUTE_ALL)
		{
			for (int i = 0; i < num; ++i)
			{
				sfprt[idxmap[d[i]]].DM_Potential += MASSPART / (ForceSoftening / 2.8);
				sfprt[idxmap[d[i]]].DM_Potential *= GravConst / fac_comov_to_phys;
			}
		}
		else
		{
			for (int i = 0; i < num; ++i)
			{
				sfprt[idxmap[d[i]]].DM_Potential *= GravConst / fac_comov_to_phys;
				sfprt[idxmap[d[i]]].DM_Potential =
					potold[i] - sfprt[idxmap[d[i]]].DM_Potential;
			}

		}

		int minidx = -1;
		MyFloat minpot = MAX_DOUBLE_NUMBER;
		for (int i = 0; i < num; ++i)
		{
			if (sfprt[idxmap[d[i]]].DM_Potential < minpot || minidx == -1)
			{
				minpot = sfprt[idxmap[d[i]]].DM_Potential;
				minidx = d[i];
			}
		}


		num_removed = 0;
		long long unbound;
		do
		{
			MyFloat mass = 0, s[3] = { 0, 0, 0 }, v[3] ={0, 0, 0};
			for (int i = 0; i < num; ++i)
			{
				int pindex = d[i];
#ifndef INTXYZ
				for (int j = 0; j < 3; ++j)
					s[j] += MASSPART * (prt[pindex].pos[j] - prt[minidx].pos[j]);
#else
				for (int j = 0; j < 3; ++j)
					s[j] += MASSPART * (prt[pindex].posi[j] - prt[minidx].posi[j] ) *INT2POS;
#endif

				for (int j = 0; j < 3; ++j)
					v[j] += MASSPART * prt[pindex].vel[j];
			}
			mass = num * MASSPART;
			for (int j = 0; j < 3; ++j)
				s[j] /= mass;
			for (int j = 0; j < 3; ++j)
				v[j] /= mass;
			MyFloat cm[3];
#ifndef INTXYZ
			cm[0] = s[0] + prt[minidx].pos[0];
			cm[1] = s[1] + prt[minidx].pos[1];
			cm[2] = s[2] + prt[minidx].pos[2];
#else
			cm[0] = s[0] + prt[minidx].posi[0]*INT2POS;
			cm[1] = s[1] + prt[minidx].posi[1]*INT2POS;
			cm[2] = s[2] + prt[minidx].posi[2]*INT2POS;

#endif

			MyFloat *bnd_energy = (MyFloat *) mymalloc (sizeof (MyFloat) * num, 224);
			for (int i = 0; i < num; ++i)
			{
				int pindex = d[i];
				MyFloat dv[3];

				MyFloat ppos[3];

#ifndef INTXYZ
				ppos[0] = prt[pindex].pos[0];
				ppos[1] = prt[pindex].pos[1];
				ppos[2] = prt[pindex].pos[2];
#else

				ppos[0] = prt[pindex].posi[0]*INT2POS;
				ppos[1] = prt[pindex].posi[1]*INT2POS;
				ppos[2] = prt[pindex].posi[2] *INT2POS ;

#endif


				for (int j = 0; j < 3; ++j)
				{
					dv[j] =
						fac_vel_to_phys * (prt[pindex].vel[j] - v[j]);
					dv[j] +=
						fac_hubbleflow * fac_comov_to_phys *
						(ppos[j] - cm[j]);
				}
				sfprt[idxmap[pindex]].DM_BindingEnergy =
					sfprt[idxmap[pindex]].DM_Potential + 0.5 * (dv[0] * dv[0] +
							dv[1] * dv[1] +
							dv[2] * dv[2]);
				bnd_energy[i] = sfprt[idxmap[pindex]].DM_BindingEnergy;
			}
			qsort (bnd_energy, num, sizeof (bnd_energy[0]),
					compare_binding_energy);
			long long j =max ((long long) 5,(long long) (MAX_UNBOUND_FRAC_BEFORE_BULK_VELOCITY_UPDATE *num));
			MyFloat energy_limit = MAX_DOUBLE_NUMBER;
			if (j < num)
			{
				energy_limit = bnd_energy[j];
			}
			else
			{
				printf ("wrong j\n");
				exit (0);
			}
			int num_unbound = 0;
			for (int i = 0; i < num; ++i)
			{
				int p = d[i];
				if (sfprt[idxmap[p]].DM_BindingEnergy > 0 && sfprt[idxmap[p]].DM_BindingEnergy > energy_limit)
				{
					num_unbound++;
					dremoved[num_removed++] = d[i];
					d[i] = d[num - 1];
					num--;
					i--;
				}
			}

			myfree (bnd_energy, 224);

			unbound = num_unbound;
			removed = num_removed;

		}
		while ( count_unbind && unbound > 0 && num >= DesLinkNgb && removed < MAX_UNBOUND_FRAC_BEFORE_POTENTIAL_UPDATE * num);

		iter++;
		if (iter > MAX_ITER_UNBIND) {
			printf ("too many iterations\n");
			exit (0);
		}
		if (phaseflag == RECOMPUTE_ALL) {
			if (removed > 0)
				phaseflag = UPDATE_ALL;
		}else{
			if (removed == 0)
			{
				phaseflag = RECOMPUTE_ALL;
				removed = 1;
			}
		}
		myfree (part, 223);
	}
	while (count_unbind-- && removed > 0 && num >= DesLinkNgb);

	myfree (potold, 222);
	myfree (dremoved, 221);
	return num;
}
/*
int subfind_unbind_bak(Body_T* prt, SFInfo* sfprt, int* d, int num) {
	int phaseflag = RECOMPUTE_ALL;
	int iter = 0;

	int num_removed = 0;
	long long removed;
	int *dremoved = (int*)mymalloc(sizeof(int) * num, 221);
	MyFloat *potold = (MyFloat*)mymalloc(sizeof(MyFloat) * num, 222);

	MyFloat fac_vel_to_phys, fac_comov_to_phys, fac_hubbleflow;
	subfind_get_factors(fac_vel_to_phys, fac_hubbleflow, fac_comov_to_phys);

	SFPart* part;
	PartTree<SFPart, tNode, tPack> tree;
	MTKDtree<SFPart, tNode, tPack>* builder;
	MyFloat left[3], right[3];
	do {
		if (phaseflag == RECOMPUTE_ALL) {
			part = (SFPart*)mymalloc(sizeof(SFPart) * num, 223);
			for (int i=0; i<num; ++i) {
#ifdef SAVE_MEM_EX
				part[i].pos = prt[d[i]].pos;
#else
				part[i].pos[0] = prt[d[i]].pos[0];
				part[i].pos[1] = prt[d[i]].pos[1];
				part[i].pos[2] = prt[d[i]].pos[2];
#endif
				part[i].pidx = d[i];
			}
			builder = new MTKDtree<SFPart, tNode, tPack>(part, num, 16);
			tree = builder->build_kdtree();
			left[0] = BOXSIZE;
			left[1] = BOXSIZE;
			left[2] = BOXSIZE;
			right[0] = 0.0;
			right[1] = 0.0;
			right[2] = 0.0;
			for (int i=0; i<num; ++i) {
				if (part[i].pos[0] < left[0]) left[0] = part[i].pos[0];
				if (part[i].pos[1] < left[1]) left[1] = part[i].pos[1];
				if (part[i].pos[2] < left[2]) left[2] = part[i].pos[2];
				if (part[i].pos[0] > right[0]) right[0] = part[i].pos[0];
				if (part[i].pos[1] > right[1]) right[1] = part[i].pos[1];
				if (part[i].pos[2] > right[2]) right[2] = part[i].pos[2];
			}
			builder->center_kdtree(0, tree.first_node, left, right);
		} else {
			part = (SFPart*)mymalloc(sizeof(SFPart) * num_removed, 223);
			for (int i=0; i<num_removed; ++i) {
#ifdef SAVE_MEM_EX
				part[i].pos = prt[dremoved[i]].pos;
#else
				part[i].pos[0] = prt[dremoved[i]].pos[0];
				part[i].pos[1] = prt[dremoved[i]].pos[1];
				part[i].pos[2] = prt[dremoved[i]].pos[2];
#endif
				part[i].pidx = dremoved[i];
			}
			builder = new MTKDtree<SFPart, tNode, tPack>(part, num_removed, 16);
			tree = builder->build_kdtree();
			left[0] = BOXSIZE;
			left[1] = BOXSIZE;
			left[2] = BOXSIZE;
			right[0] = 0.0;
			right[1] = 0.0;
			right[2] = 0.0;
			for (int i=0; i<num_removed; ++i) {
				if (part[i].pos[0] < left[0]) left[0] = part[i].pos[0];
				if (part[i].pos[1] < left[1]) left[1] = part[i].pos[1];
				if (part[i].pos[2] < left[2]) left[2] = part[i].pos[2];
				if (part[i].pos[0] > right[0]) right[0] = part[i].pos[0];
				if (part[i].pos[1] > right[1]) right[1] = part[i].pos[1];
				if (part[i].pos[2] > right[2]) right[2] = part[i].pos[2];
			}
			builder->center_kdtree(0, tree.first_node, left, right);
			for (int i=0; i<num; ++i)
				potold[i] = sfprt[idxmap[d[i]]].DM_Potential;
		}
		if (part == NULL ) {
			LOG(0, "error part malloc in subfind.cc");
			exit(0);
		}

		double ptime = dtime();
		subfind_potential_compute(prt, sfprt, d, num, tree);
		builder->MTKDtree_free();
		delete builder;

		if (phaseflag == RECOMPUTE_ALL) {
			for (int i=0; i<num; ++i) {
				sfprt[idxmap[d[i]]].DM_Potential += MASSPART / (ForceSoftening / 2.8);
				sfprt[idxmap[d[i]]].DM_Potential *= GravConst / fac_comov_to_phys;
			}
		} else {
			for (int i=0; i<num; ++i) {
				sfprt[idxmap[d[i]]].DM_Potential *= GravConst / fac_comov_to_phys;
				sfprt[idxmap[d[i]]].DM_Potential = potold[i] - sfprt[idxmap[d[i]]].DM_Potential;
			}

		}

		int minidx = -1;
		MyFloat minpot = MAX_DOUBLE_NUMBER;
		for (int i=0; i<num; ++i) {
			if (sfprt[idxmap[d[i]]].DM_Potential < minpot || minidx == -1) {
				minpot = sfprt[idxmap[d[i]]].DM_Potential;
				minidx = d[i];
			}
		}
	

		num_removed = 0;
		long long unbound;
		do {
			MyFloat mass = 0, s[3] = {0, 0, 0}, v[3] = {0, 0, 0};
			for (int i=0; i<num; ++i) {
				int pindex = d[i];
				for (int j=0; j<3; ++j)
					s[j] += MASSPART * (prt[pindex].pos[j] - prt[minidx].pos[j]);  
				for (int j=0; j<3; ++j)
					v[j] += MASSPART * prt[pindex].vel[j];
			}
			mass = num * MASSPART;
			for (int j=0; j<3; ++j)
				s[j] /= mass;  
			for (int j=0; j<3; ++j)
				v[j] /= mass;
			MyFloat cm[3];
			cm[0] = s[0] + prt[minidx].pos[0];
			cm[1] = s[1] + prt[minidx].pos[1];
			cm[2] = s[2] + prt[minidx].pos[2];

			MyFloat* bnd_energy = (MyFloat*)mymalloc(sizeof(MyFloat) * num, 224);
			for (int i=0; i<num; ++i) {
				int pindex = d[i];
				MyFloat dv[3];
				for (int j=0; j<3; ++j) {
					dv[j] = fac_vel_to_phys * (prt[pindex].vel[j] - v[j]);
					dv[j] += fac_hubbleflow * fac_comov_to_phys * (prt[pindex].pos[j] - cm[j]);
				}
				sfprt[idxmap[pindex]].DM_BindingEnergy = sfprt[idxmap[pindex]].DM_Potential + 0.5 * (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
				bnd_energy[i] = sfprt[idxmap[pindex]].DM_BindingEnergy;
			}
			qsort(bnd_energy, num, sizeof(bnd_energy[0]), compare_binding_energy);
			long long j = max((long long)5, (long long)(MAX_UNBOUND_FRAC_BEFORE_BULK_VELOCITY_UPDATE * num));
			MyFloat energy_limit = MAX_DOUBLE_NUMBER;
			if (j < num) {
				energy_limit = bnd_energy[j];
			} else {
				printf("wrong j\n");
				exit(0);
			}
			int num_unbound = 0;
			for (int i=0; i<num; ++i) {
				int p = d[i];
				if (sfprt[idxmap[p]].DM_BindingEnergy > 0 && sfprt[idxmap[p]].DM_BindingEnergy > energy_limit) {
					num_unbound++;
					dremoved[num_removed++] = d[i];
					d[i] = d[num - 1];
					num--;
					i--;
				}
			}	
			
			myfree(bnd_energy, 224);

			unbound = num_unbound;
			removed = num_removed;

		}while(unbound>0 && num >= DesLinkNgb && removed < MAX_UNBOUND_FRAC_BEFORE_POTENTIAL_UPDATE * num);

		iter++;	
		if (iter > MAX_ITER_UNBIND) {
			printf("too many iterations\n");
			exit(0);
		}
		if (phaseflag == RECOMPUTE_ALL) {
			if (removed > 0)
				phaseflag = UPDATE_ALL;
		} else {
			if (removed == 0) {
				phaseflag = RECOMPUTE_ALL;
				removed = 1;
			}
		}

		myfree(part, 223);
	} while (removed > 0 && num >= DesLinkNgb);
	myfree(potold, 222);
	myfree(dremoved, 221);
	return num;
}
*/
void subfind_unbind_independent_ones(SFCandidate* coll_candidates, int count_cand, Body_T* prt, unsigned long totnpart, SFInfo* sfprt, SFDPart* part, SFSMSort* PPS, SFMark* sfmrk, unsigned long npart, Group *host_group) {
	int* unbind_list = (int*)mymalloc(sizeof(int) * npart, 225);
	if (unbind_list == NULL ) {
		LOG(0, "error unbind_list malloc in subfind.cc");
		exit(0);
	}
	qsort(coll_candidates, count_cand, sizeof(SFCandidate), compare_nsubs);
	for (int k=0,ii=0; k<count_cand; ++k) {
		if (coll_candidates[k].parent == 0) {
			int i = PPS[ii].pidx;
			while(sfmrk[i].submark < coll_candidates[k].nsub) {
				ii++;
				i = PPS[ii].pidx;
				if (i >= npart) {
					printf("i >= npart\n");
					exit(0);
				}
			}
			if (sfmrk[i].submark >=0 && sfmrk[i].submark < HIGHBIT) {
				int len = 0;
				int nsubs = sfmrk[i].submark;
				if (nsubs != coll_candidates[k].nsub) {
					printf("subs error\n");
					exit(0);
				}
				while (i < npart) {
					if (sfmrk[i].submark == nsubs) {
						sfmrk[i].submark = HIGHBIT;
						if ((sfmrk[i].origintask & HIGHBIT) == 0) {
							unbind_list[len] = i;
							len++;
						} 
						ii++;
						i = PPS[ii].pidx;
					}
					else
						break;
				}
				for (int i=0; i<len; ++i) {
					sfprt[idxmap[part[unbind_list[i]].pidx]].invidx = unbind_list[i];
					unbind_list[i] = part[unbind_list[i]].pidx;
					if (unbind_list[i] < 0 || unbind_list[i] >= totnpart) {
						printf("unbind_list element > npart of proc\n");
						exit(0);

					}
				}
				len = subfind_unbind(prt, sfprt, unbind_list, len, host_group);
				// go from particle indices to group indices
				for (int i=0; i<len; ++i) {
					unbind_list[i] = sfprt[idxmap[unbind_list[i]]].invidx;	
					if (unbind_list[i] < 0 || unbind_list[i] >= npart) {
						printf("unbind_list element > npart of group\n");
						exit(0);
					}
				}
				if (len >= DesLinkNgb) {
					coll_candidates[k].bound_length = len;
					for (int j=0; j<len;++j) {
						sfmrk[unbind_list[j]].submark = nsubs;
					}
				} else {
					coll_candidates[k].bound_length = 0;
				}
			}
		}
	}
	myfree(unbind_list, 225);
	
}	

void subfind_determine_sub_halo_properties(Body_T* prt, SFInfo* sfprt, int* d, int num, subGroupCatalog& sgrpcat) {
		MyFloat fac_vel_to_phys, fac_comov_to_phys, fac_hubbleflow;
		subfind_get_factors(fac_vel_to_phys, fac_hubbleflow, fac_comov_to_phys);

		int minidx = -1;
		MyFloat minpot = MAX_DOUBLE_NUMBER;
		for (int i=0; i<num; ++i) {
			if (sfprt[idxmap[d[i]]].DM_Potential < minpot || minidx == -1) {
				minpot = sfprt[idxmap[d[i]]].DM_Potential;
				minidx = d[i];
			}
		}
		MyFloat pos[3];

#ifndef INTXYZ
		pos[0] = prt[minidx].pos[0];
		pos[1] = prt[minidx].pos[1];
		pos[2] = prt[minidx].pos[2];
#else

		pos[0] = prt[minidx].posi[0]*INT2POS;
		pos[1] = prt[minidx].posi[1]*INT2POS;
		pos[2] = prt[minidx].posi[2] *INT2POS;

#endif
		minidx = -1;
		minpot = MAX_DOUBLE_NUMBER;
		for (int i=0; i<num; ++i) {
			if (sfprt[idxmap[d[i]]].DM_BindingEnergy < minpot || minidx == -1) {
				minpot = sfprt[idxmap[d[i]]].DM_BindingEnergy;
				minidx = d[i];
			}
		}
		int mostboundid = minidx;

		MyFloat mass = 0, s[3] = {0, 0, 0}, v[3] = {0, 0, 0};
		for (int i=0; i<num; ++i) {
			int pindex = d[i];
#ifndef INTXYZ
			for (int j=0; j<3; ++j)
				s[j] += MASSPART * (prt[pindex].pos[j] - pos[j]);  
#else
			for (int j=0; j<3; ++j)
				s[j] += MASSPART * (prt[pindex].posi[j]*INT2POS - pos[j]);  

#endif
			for (int j=0; j<3; ++j)
				v[j] += MASSPART * prt[pindex].vel[j];
		}
		mass = num * MASSPART;

		MyFloat vel[3];
		for (int j=0; j<3; ++j) {
			s[j] /= mass;  
			v[j] /= mass;
			vel[j] = fac_vel_to_phys * v[j];
		}
		MyFloat cm[3];
#ifndef INTXYZ	
		cm[0] = s[0] + prt[minidx].pos[0];
		cm[1] = s[1] + prt[minidx].pos[1];
		cm[2] = s[2] + prt[minidx].pos[2];
#else
		cm[0] = s[0] + prt[minidx].posi[0] * INT2POS;
		cm[1] = s[1] + prt[minidx].posi[1] * INT2POS;
		cm[2] = s[2] + prt[minidx].posi[2] * INT2POS;
#endif

  		MyFloat lx = 0, ly = 0, lz = 0;
		R2ListSort* rr_list = (R2ListSort*)mymalloc(sizeof(R2ListSort) * (num+1), 226);
		for (int i=0; i<num; ++i) {
			int pindex = d[i];
			MyFloat dx[3], dv[3];
#ifndef INTXYZ
			for (int j=0; j<3; ++j)
				dx[j] = (prt[pindex].pos[j] - cm[j]) * fac_comov_to_phys;
#else
			for (int j=0; j<3; ++j)
				dx[j] = (prt[pindex].posi[j] * INT2POS - cm[j]) * fac_comov_to_phys;
#endif
			for (int j=0; j<3; ++j)
        		{
          			dv[j] = fac_vel_to_phys * (prt[pindex].vel[j] - v[j]);
          			dv[j] += fac_hubbleflow * dx[j];
        		}

			lx += MASSPART * (dx[1] * dv[2] - dx[2] * dv[1]);
			ly += MASSPART * (dx[2] * dv[0] - dx[0] * dv[2]);
			lz += MASSPART * (dx[0] * dv[1] - dx[1] * dv[0]);

			MyFloat dxyz[3];
#ifndef INTXYZ
			for (int j=0; j<3; ++j)
				dxyz[j] = (prt[pindex].pos[j] - pos[j]) * fac_comov_to_phys;
#else
			for (int j=0; j<3; ++j)
				dxyz[j] = (prt[pindex].posi[j]*INT2POS - pos[j]) * fac_comov_to_phys;
#endif

			MyFloat r2 = dxyz[0] * dxyz[0] + dxyz[1] * dxyz[1] + dxyz[2] * dxyz[2];
			MyFloat r = sqrt(r2);
			rr_list[i].r = r;
			rr_list[i].mass = MASSPART; 
		}

		MyFloat spin[3];
		spin[0] = lx / mass; 
		spin[1] = ly / mass; 
		spin[2] = lz / mass; 

		qsort(rr_list, num, sizeof(R2ListSort), compare_rotcurve);
		for (int i=1; i<num; ++i)
			rr_list[i].mass = rr_list[i-1].mass + rr_list[i].mass;
		MyFloat halfmassrad = 0;
		MyFloat max = 0, maxrad = 0;
		rr_list[num].mass = 0;
		rr_list[num].r = 0;
		for (int i=num-1; i>=0; --i) {
			if (i > 5 && rr_list[i].mass > max * rr_list[i].r) {
				max = rr_list[i].mass / rr_list[i].r;
				maxrad = rr_list[i].r;
			}
			if (rr_list[i].mass < 0.5 * mass && rr_list[i+1].mass >= 0.5 * mass) {
				halfmassrad = 0.5 * (rr_list[i].r + rr_list[i+1].r);
			}
		}
		MyFloat vmax = sqrt(GravConst * max);
		MyFloat vmaxrad = maxrad;	
	
		subGroup subhalo;
		subhalo.Len = num;
		subhalo.pidx.resize(num);
		for (int i=0; i<num; ++i) {
			subhalo.pidx[i] = d[i];
		}
		sort(subhalo.pidx.begin(), subhalo.pidx.end(), [& sfprt](int i1, int i2){return sfprt[idxmap[i1]].DM_Potential < sfprt[idxmap[i2]].DM_Potential;});
		subhalo.in_group = sfprt[idxmap[d[0]]].in_group;
		subhalo.Mass = mass;
		subhalo.pos[0] = pos[0];
		subhalo.pos[1] = pos[1];
		subhalo.pos[2] = pos[2];
		subhalo.CM[0] = cm[0];
		subhalo.CM[1] = cm[1];
		subhalo.CM[2] = cm[2];
		subhalo.vel[0] = vel[0];
		subhalo.vel[1] = vel[1];
		subhalo.vel[2] = vel[2];
		subhalo.spin[0] = spin[0];
		subhalo.spin[1] = spin[1];
		subhalo.spin[2] = spin[2];
		subhalo.SubMostBoundID = mostboundid;
		subhalo.SubVmax = vmax;
		subhalo.SubVmaxRad = vmaxrad;
		subhalo.SubHalfMassRad = halfmassrad;
		sgrpcat.subgrps.push_back(subhalo);
		sgrpcat.subgrp_number++;
	
		myfree(rr_list, 226);
}

public:
SubFind(SIMInfo<Body_T>& sim) {
	part_inproc = sim.part;
	npart_inproc = sim.npart;
	N_SIDEMESH = sim.N_SIDEMESH;
	PROC_RANK = sim.PROC_RANK;
	PROC_SIZE = sim.PROC_SIZE;
	NPART_TOTAL = sim.NPART_TOTAL; 
	MASSPART = sim.MASSPART; 
	BOXSIZE = sim.BOXSIZE; 
	Hubble0 = sim.Hubble0; 
	OmegaM0 = sim.OmegaM0; 
	OmegaX0 = sim.OmegaX0; 
	redshift = sim.redshift;
	float Time = 1.0 / (1.0 + redshift);
	velFac = sqrt(Time) * Time;
	
	double SofteningMaxPhys = 0.0005, SofteningComoving = 0.001;
	double SofteningTable;
	if (SofteningComoving * Time > SofteningMaxPhys)
		SofteningTable = SofteningMaxPhys / Time;
	else
		SofteningTable = SofteningComoving;
	ForceSoftening = 2.8 * SofteningTable;
}
~SubFind() {}

void subfind_density(SFInfo* sfprt, int np_ingrp, MyFloat LINKLENGTH) {
	Body_T* prt = part_inproc;
	unsigned long npart = npart_inproc;
	double t0 = dtime();

	SFPart* part = (SFPart*)mymalloc(sizeof(SFPart) * npart, 227);
	if (part == NULL ) {
		LOG(0, "[%d] error part malloc in subfind.cc", PROC_RANK);
		exit(0);
	}

	for (unsigned long i=0; i<npart; ++i) {
#ifdef SAVE_MEM_EX
		part[i].pos = prt[i].pos;
#else


#ifndef INTXYZ
		part[i].pos[0] = prt[i].pos[0];
		part[i].pos[1] = prt[i].pos[1];
		part[i].pos[2] = prt[i].pos[2];
#else
		part[i].pos[0] = prt[i].posi[0] * INT2POS;
		part[i].pos[1] = prt[i].posi[1] * INT2POS;
		part[i].pos[2] = prt[i].posi[2] * INT2POS;
#endif

#endif
		part[i].pidx = i;
	}

	MTKDtree<SFPart, tNode, tPack> builder(part, npart, 16);
	PartTree<SFPart, tNode, tPack> tree = builder.build_kdtree();

	LOG(2, "[%d] finish building tree of in-group particles time %lf", PROC_RANK, dtime() - t0);
	double t1 = dtime();
	MyFloat mean_sep = BOXSIZE/pow((MyFloat)NPART_TOTAL, 0.333333333);
	MyFloat estimate_h = DesNumNgb * LINKLENGTH * mean_sep;

	// do it only for particles in a group
	int nth = omp_get_max_threads();
	LOG(1, "[%d] subfind density: openmp num of threads %d", PROC_RANK, nth);
#pragma omp parallel for num_threads(nth) schedule(dynamic)
	for (unsigned long i=0; i<npart; ++i) {
		// skip those particles not in  group
		if (idxmap[part[i].pidx] == -1) {
			continue;
		}
		PartIndex* sidx;
		MyFloat maxh = estimate_h/3.;
		KNN3DIndex<SFPart, tNode, tPack> knnidx(&tree);
		sidx = knnidx.KNN_search(&part[i], 1, DesNumNgb, maxh);
		
		while (abs(sidx->matched_num - DesNumNgb) > MaxNumNgbDeviation) {
			if (maxh > estimate_h) {
				printf("[%d] part %lu matched_num %d in the radius of %lf (should be around %d)\n", PROC_RANK, i, sidx->matched_num, maxh, DesNumNgb);
				exit(0);
			}
			maxh = maxh * 1.5;
			sidx = knnidx.KNN_search(&part[i], 1, DesNumNgb, maxh);
		}

		MyFloat h, r, u, hinv, hinv3, wk;

		MyFloat h0 = sidx->distance_list[sidx->matched_num-2];
		MyFloat h1 = sidx->distance_list[sidx->matched_num-1];
		h = min(h0 * 1.001, (h0+h1)/2.0);
		sfprt[idxmap[part[i].pidx]].DM_Hsml = h;

		hinv = 1.0 / h;
		hinv3 = hinv * hinv * hinv;
		// add DesNumNgb-1 particles' density
		for (int j=0; j<sidx->matched_num-1; ++j) {
			r = sidx->distance_list[j];
			u = r * hinv;
			if (u < 0.5)
				wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
			else
				wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
			sfprt[idxmap[part[i].pidx]].DM_Density += wk * MASSPART;
		}
		// add itself's density
		sfprt[idxmap[part[i].pidx]].DM_Density += hinv3 * KERNEL_COEFF_1 * MASSPART;
		delete []sidx;
	}
	builder.MTKDtree_free();

	LOG(2, "[%d] finish computing density of each in-group particle (%d of %lu) time %lf", PROC_RANK, np_ingrp, npart, dtime() - t1);
 
	myfree(part, 227);
}

 
void subfind_process_single_group(subGroupCatalog& sgrpcat, SFInfo* sfprt, vector<int>& pidx, MyFloat link_length, Group *host_group) {
	Body_T* prt = part_inproc;
	unsigned long totnpart = npart_inproc;
	int npart = pidx.size();
#ifdef DEBUG 
	FILE *fd;
	if (npart > DBGN) {
		fd = fopen("subfind.txt","w");
	}
#endif
	SFMark* sfmrk = (SFMark*)mymalloc(sizeof(SFMark) * npart, 228);
	if (sfmrk == NULL ) {
		LOG(0, "[%d] error sfmrk malloc in subfind.cc", PROC_RANK);
		exit(0);
	}

	SFDPart* part = (SFDPart*)mymalloc(sizeof(SFDPart) * npart, 229);
	if (part == NULL ) {
		LOG(0, "[%d] error part malloc in subfind.cc", PROC_RANK);
		exit(0);
	}
 
	for (int i=0; i<npart; ++i) {
#ifndef INTXYZ
		part[i].pos[0] = prt[pidx[i]].pos[0];
		part[i].pos[1] = prt[pidx[i]].pos[1];
		part[i].pos[2] = prt[pidx[i]].pos[2];
#else
		part[i].pos[0] = prt[pidx[i]].posi[0] * INT2POS;
		part[i].pos[1] = prt[pidx[i]].posi[1] * INT2POS;
		part[i].pos[2] = prt[pidx[i]].posi[2] * INT2POS;
#endif

		part[i].DM_Hsml = sfprt[idxmap[pidx[i]]].DM_Hsml;
		part[i].DM_Density = sfprt[idxmap[pidx[i]]].DM_Density;
		part[i].pidx = pidx[i];
	}

#ifdef DEBUG
	if (npart > DBGN) {
		fwrite_pos(part, npart, "part.txt");
	}
#endif

	subfind_find_nearesttwo(part, npart, link_length);

	int max_coll_candidates = max(npart / 50., 200.); 
	SFCandidate* coll_candidates = (SFCandidate*)mymalloc(sizeof(SFCandidate) * max_coll_candidates, 230);
	if (coll_candidates == NULL ) {
		LOG(0, "[%d] error coll_candidates malloc in subfind.cc", PROC_RANK);
		exit(0);
	}
	int count_cand;
	SFLinkList* linklists = (SFLinkList*)mymalloc(sizeof(SFLinkList) * npart, 231);
	if (linklists == NULL ) {
		LOG(0, "[%d] error linklists malloc in subfind.cc", PROC_RANK);
		exit(0);
	}

	count_cand = subfind_find_coll_candidates(part, npart, coll_candidates, max_coll_candidates, linklists);
	qsort(coll_candidates, count_cand, sizeof(SFCandidate), compare_len);

	for (int i=0; i<npart; ++i) {
		linklists[i].tail = -1;
	}
	for(int i=0; i<count_cand; ++i) {
		coll_candidates[i].parent = 0;
	}

	int max_length, count_leaves;
	long long nremaining = count_cand;

	do {
		// Let's see which coll_candidates can be unbound independent from each other.
		for (int k=0; k<count_cand; ++k) {
			coll_candidates[k].nsub = k;
			coll_candidates[k].subnr = k;
		}
		// serial rank in one group
		qsort(coll_candidates, count_cand, sizeof(SFCandidate), compare_rank);
		for (int k=0; k<count_cand; ++k) {
			if (coll_candidates[k].parent >=0) {
				coll_candidates[k].parent = 0;
				for (int j=k+1; j<count_cand; ++j) {
					//independent
					if (coll_candidates[j].rank > coll_candidates[k].rank + coll_candidates[k].len)
						break;
					//ignore
					if (coll_candidates[j].parent < 0)
						continue;
					//k enclose j
					if (coll_candidates[k].rank + coll_candidates[k].len >= coll_candidates[j].rank + coll_candidates[j].len)
						coll_candidates[k].parent++;
					else {
						printf("odd\n");
						exit(0);
					}
				}
			
			}
		}
		qsort(coll_candidates, count_cand, sizeof(SFCandidate), compare_subnr);
		count_leaves = 0;
		max_length = 0;
		for(int i=0; i<count_cand; ++i) {
			if (coll_candidates[i].parent == 0) {
				if (0 && coll_candidates[i].len > 0.2 * totnpart) {
					coll_candidates[i].parent++;
					printf("len seems large\n");
					exit(0);
				} else {
					if (coll_candidates[i].len > max_length)
						max_length = coll_candidates[i].len;
					count_leaves++;
				}
			}
		}
		if (count_leaves <= 0) {
		//	printf("two few, let's do the rest collectively\n");
		//	exit(0);
			break;
		}
		if (max_length > 0.5 * totnpart) {
		//	printf("two big coll_candidates, do the rest collectively\n");
		//	exit(0);
			break;
		}

		nremaining -= count_leaves;

		for (unsigned long i=0; i<npart; ++i) {
			sfmrk[i].origintask = 0;
			sfmrk[i].submark = HIGHBIT;
		}

		for(int i=0; i<npart; ++i) {
			// this part is already bound to a substructure
			if(linklists[i].tail >= 0)
				sfmrk[i].origintask |= HIGHBIT;
		}
		
		
		int nsubs = 0;
		for (int k=0; k<count_cand; ++k) {
			int len = coll_candidates[k].len;
			int parent = coll_candidates[k].parent;
			if (parent == 0) {
				int p = coll_candidates[k].head;
				for (int i=0; i<len; ++i) {
					if (sfmrk[p].submark != HIGHBIT) {
						printf("not HIGHBIT\n");
						exit(0);
					}
					if (p<0) {
						printf("Bummer\n");
						exit(0);
					}
					sfmrk[p].submark = nsubs;
					p = linklists[p].next;
				}
			}
			nsubs++;
		}

		SFSMSort* PPS = (SFSMSort*)mymalloc(sizeof(SFSMSort) * npart, 232);
		if (PPS == NULL ) {
			LOG(0, "[%d] error PPS malloc in subfind.cc", PROC_RANK);
			exit(0);
		}
		for (unsigned long i=0; i<npart; ++i) {
			PPS[i].pidx = i;
			PPS[i].submark = sfmrk[i].submark;
		}
		qsort(PPS, npart, sizeof(SFSMSort), compare_PPS);
		subfind_unbind_independent_ones(coll_candidates, count_cand, prt, totnpart, sfprt, part, PPS, sfmrk, npart, host_group);	

		myfree(PPS, 232);

		// now mark the bound particles
		for (int i=0; i<npart; ++i) {
			if (sfmrk[i].submark >=0 && sfmrk[i].submark < nsubs)
				linklists[i].tail = sfmrk[i].submark;
		}
		for (int i=0; i<count_cand; ++i) {
			if (coll_candidates[i].parent == 0)
				coll_candidates[i].parent = -1;
		}

		
	} while(count_leaves > 0);

	int* unbind_list = (int*)mymalloc(sizeof(int) * npart, 233);
	if (unbind_list == NULL ) {
		LOG(0, "error unbind_list malloc in subfind.cc");
		exit(0);
	}
	
	// Now we do the unbinding of the subhalo candidates that contain other subhalo candidates;
	for (int k=0; k<count_cand; ++k) {
		int len, parent, nsubs;
		len = coll_candidates[k].len;
		nsubs = coll_candidates[k].nsub;
		parent = coll_candidates[k].parent;
		if (parent >= 0) {
			int LocalLen = 0;
			int p = coll_candidates[k].head;
			for (int i=0; i<coll_candidates[k].len; ++i) {
				if (p<0) {
					printf("Bummer\n");
					exit(0);
				}
				if (linklists[p].tail < 0) {
					if (p >= npart) {
						printf("p index >= npart of group\n");
						exit(0);
					}
					unbind_list[LocalLen] = p;
					LocalLen++;
				}
				p = linklists[p].next;
			}
			if (LocalLen > npart) {
				printf("LocalLen > npart of group\n");
				exit(0);
			}
			for (int i=0; i<LocalLen; ++i) {
				sfprt[idxmap[part[unbind_list[i]].pidx]].invidx = unbind_list[i];
				unbind_list[i] = part[unbind_list[i]].pidx;
				if (unbind_list[i] < 0 || unbind_list[i] >= totnpart) {
					printf("unbind_list element > npart of proc\n");
					exit(0);

				}
			}
			
			LocalLen = subfind_unbind(prt, sfprt, unbind_list, LocalLen, host_group);

			// go from particle indices to group indices
			for (int i=0; i<LocalLen; ++i) {
				unbind_list[i] = sfprt[idxmap[unbind_list[i]]].invidx;	
				if (unbind_list[i] < 0 || unbind_list[i] >= npart) {
					printf("unbind_list element > npart of group\n");
					exit(0);
				}
			}
			int oldlen = len;
			len = LocalLen;
		
			if (len >= DesLinkNgb) {
				for (int i=0; i<LocalLen; i++)
					linklists[unbind_list[i]].tail = nsubs;
				coll_candidates[k].bound_length = len;
			} else {
				coll_candidates[k].bound_length = 0;

			}
		}
	}

	// get the total substructure count
	int countall = 0;	
	for (int k=0; k<count_cand; ++k) {
		if (coll_candidates[k].bound_length >= DesLinkNgb) {
			if (coll_candidates[k].len < DesLinkNgb) {
				printf("bound_length, len bummer\n");
				exit(0);
			}
			countall++;
		}
	}

	// now determine the parent subhalo for each candidate
	qsort(coll_candidates, count_cand, sizeof(SFCandidate), compare_bound_length);
	for (int k=0; k<count_cand; ++k) {
		coll_candidates[k].subnr = k;
		coll_candidates[k].parent = 0;
	}
	qsort(coll_candidates, count_cand, sizeof(SFCandidate), compare_rank);
	for (int k=0; k<count_cand; ++k) {
		for (int j=k+1; j<count_cand; ++j) {
			if (coll_candidates[j].rank > coll_candidates[k].rank + coll_candidates[k].len)
				break;
			if (coll_candidates[k].rank + coll_candidates[k].len >= coll_candidates[j].rank + coll_candidates[j].len) {
				if (coll_candidates[k].bound_length >= DesLinkNgb)
					coll_candidates[j].parent = coll_candidates[k].subnr;
			} else {
				printf("rank len wrong\n");
				exit(0);
			}
		}
	}
	qsort(coll_candidates, count_cand, sizeof(SFCandidate), compare_subnr);

	// Now let's save some properties of substructures
	int subnr = 0;
 	for (int k=0; k<count_cand; ++k) {
		int len, parent, nsubs;
		len = coll_candidates[k].bound_length;
		nsubs = coll_candidates[k].nsub;
		parent = coll_candidates[k].parent;
		if (len > 0) {
			int LocalLen = 0;
			int p = coll_candidates[k].head;
			for (int i=0; i<coll_candidates[k].len; ++i) {
				if (linklists[p].tail == nsubs) {
					unbind_list[LocalLen] = p;
					LocalLen++;
				}
				p = linklists[p].next;
			}
			for (int i=0; i<LocalLen; ++i) {
				unbind_list[i] = part[unbind_list[i]].pidx;
				if (unbind_list[i] < 0 || unbind_list[i] >= totnpart) {
					printf("index out of local npart\n");
					exit(0);
				}
			}
			subfind_determine_sub_halo_properties(prt, sfprt, unbind_list, LocalLen, sgrpcat);
			sgrpcat.subgrps[sgrpcat.subgrp_number-1].SubRankInGr = subnr;
			sgrpcat.subgrps[sgrpcat.subgrp_number-1].SubParentRank = parent;
			subnr++;
#ifdef DEBUG
			if (npart > DBGN) {
				subGroup* sh = &sgrpcat.subgrps[sgrpcat.subgrp_number-1];
				if (sampling_print(sh->pos))
					fprintf(fd, "%lf %lf %lf : %lf\n", sh->pos[0], sh->pos[1], sh->pos[2], sh->SubVmaxRad);
			}
#endif
			
		}
	}
	host_group->nsubhalo = subnr;
	LOG(3, "this group (%lu) has %d subhalo.", npart, subnr);
	myfree(unbind_list, 233);
	myfree(linklists, 231);
	myfree(coll_candidates, 230);
	myfree(part, 229);
	myfree(sfmrk, 228);
#ifdef DEBUG 
	if (npart > DBGN) {
		fclose(fd);
	}
#endif
}

};

#endif
