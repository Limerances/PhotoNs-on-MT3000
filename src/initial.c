#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <math.h>

#include "typedef.h"
#include "photoNs.h"
#include "setparam.h"
#include "snapshot.h"

void create_mpi_datatype() {
//	MPI_Datatype mpi_group_comp;
	MPI_Type_contiguous(sizeof(Body), MPI_BYTE, &Bodytype);
	MPI_Type_commit(&Bodytype);
	MPI_Type_contiguous(sizeof(BodyBnd), MPI_BYTE, &BodyBndtype);
	MPI_Type_commit(&BodyBndtype);
}

void free_mpi_datatype() {

	MPI_Type_free(&Bodytype);
	MPI_Type_free(&BodyBndtype);
}

void decomp2d_init_() {
	int local_xstart[3];
	int local_xend[3];
	int local_xsize[3];
	int nside[3] = { NSIDEMESH, NSIDEMESH, NSIDEMESH };
	int vproc[2] = {NSIDE0PROC, NSIDE1PROC};

	decomp_2d_init(nside, vproc, local_xstart, local_xend, local_xsize);
}

void decomp2d_free_() {
	decomp_2d_free();
}

void initialize_ic() {
	if ( 0 != NSIDE0PROC % NSIDE0SUDOM ) {
		printf(" error : nside 0 !\n");
		exit(0);
	}
	if ( 0 != NSIDE1PROC % NSIDE1SUDOM ) {
		printf(" error : nside 1 !\n");
		exit(0);
	}

	if ( 0 != NSIDEMESH % NSIDE0PROC ) {
		printf(" error : nside 0 p !\n");
		exit(0);
	}

	if ( 0 != NSIDEMESH % NSIDE1PROC ) {
		printf(" error : nside 1 p!\n");
		exit(0);
	}

	if ( 0 != NSIDEMESH % NSIDE0SUDOM ) {
		printf(" error : nside 0 m !\n");
		exit(0);
	}

	if ( 0 != NSIDEMESH % NSIDE1SUDOM ) {
		printf(" error : nside 1 m!\n");
		exit(0);
	}

	if ( NSIDE0PROC * NSIDE1PROC  !=  NSIDE0SUDOM * NSIDE1SUDOM * NDOMinSUDOM  ) {
		printf(" error : num NDOMinSUDOM !\n");
		exit(0);
	}
	if ( PROC_SIZE != NSIDE0PROC * NSIDE1PROC  ) {
		printf(" error : num proc !\n");
		exit(0);
	}

	setup_param();
}

void initial_restart(char file_path[], int snap_number) {

	if ( 0 != NSIDE0PROC % NSIDE0SUDOM ) {
		printf(" error : nside 0 !\n");
		exit(0);
	}
	if ( 0 != NSIDE1PROC % NSIDE1SUDOM ) {
		printf(" error : nside 1 !\n");
		exit(0);
	}

	if ( 0 != NSIDEMESH % NSIDE0PROC ) {
		printf(" error : nside 0 p !\n");
		exit(0);
	}

	if ( 0 != NSIDEMESH % NSIDE1PROC ) {
		printf(" error : nside 1 p!\n");
		exit(0);
	}

	if ( 0 != NSIDEMESH % NSIDE0SUDOM ) {
		printf(" error : nside 0 m !\n");
		exit(0);
	}

	if ( 0 != NSIDEMESH % NSIDE1SUDOM ) {
		printf(" error : nside 1 m!\n");
		exit(0);
	}

	if ( NSIDE0PROC * NSIDE1PROC  !=  NSIDE0SUDOM * NSIDE1SUDOM * NDOMinSUDOM  ) {
		printf(" error : num NDOMinSUDOM !\n");
		exit(0);
	}
	if ( PROC_SIZE != NSIDE0PROC * NSIDE1PROC  ) {
		printf(" error : num proc !\n");
		exit(0);
	}


	char file_path_name[128];
	sprintf(file_path_name, "%s/snapshot_%d", file_path, snap_number);
	input_snapshot(file_path_name, snap_number, PROC_RANK);

	setup_param();

	domain_restart() ;

	MPI_Barrier(MPI_COMM_WORLD);
}



void domain_free() {
	if (proc_rank_sudom != NULL )
		myfree(proc_rank_sudom, 1);
}

int domain_setup() {
	proc_rank_sudom = (int*) mymalloc(sizeof(int)*NSIDE0SUDOM * NSIDE1SUDOM, 1);

	int n, c, is_sudom;

	cpuType dx = BOXSIZE / NDOMinSUDOM;
	cpuType dy = BOXSIZE / NSIDE0SUDOM;
	cpuType dz = BOXSIZE / NSIDE1SUDOM;


	for (n=0; n<PROC_SIZE; n++) {
		int idx = n;

		int dside0 = NSIDE0PROC/NSIDE0SUDOM;
		int dside1 = NSIDE1PROC/NSIDE1SUDOM;

		int py = idx / NSIDE1PROC ; 
		int pz = idx % NSIDE1PROC ; 

		int sy = py / dside0;
		int sz = pz / dside1;

		int sudom_r  = sy * NSIDE1SUDOM + sz;

		if ( (0 == py % dside0) && (0 == pz % dside1) ) {
			is_sudom = 1;
			proc_rank_sudom[sudom_r] = idx;
		}
		else {
			is_sudom = 0;
		}
		if (PROC_RANK == idx) {
			SUDOM_RANK = sudom_r;
			isSUDOM = is_sudom;
		}
	}

	for (n=0, c=0; n<PROC_SIZE; n++) {
		int idx = n;

		int dside0 = NSIDE0PROC/NSIDE0SUDOM;
		int dside1 = NSIDE1PROC/NSIDE1SUDOM;

		int py = idx / NSIDE1PROC ; 
		int pz = idx % NSIDE1PROC ; 

		int sy = py / dside0;
		int sz = pz / dside1;

		int sudom_r  = sy * NSIDE1SUDOM + sz;

		if ( sudom_r == SUDOM_RANK ) {
			dom_grp[c] = n;
			c ++;
		}

		if (PROC_RANK == idx) {

			BOX0L = sy * dy;
			BOX0H = (sy + 1 )* dy;
			BOX1L = sz * dz;
			BOX1H = (sz + 1 )* dz;

			if (sy == NSIDE0SUDOM -1 )BOX0H = BOXSIZE;
			if (sz == NSIDE1SUDOM -1 )BOX1H = BOXSIZE;
		}

	}

	PROC_RANK_SUDOM = dom_grp[0];

	int first = NDOM_HEAD;
	int dmesh = NSIDEMESH / (NDOMinSUDOM - first);
	int pad = NSIDEMESH % (NDOMinSUDOM - first);
	int m0, m1, dm;
	m0 = m1 = 0;
	for (n=0; n<NDOMinSUDOM; n++){
		int idx = n-first;

		if ( n < first ) {
			//	COM_DOM = 0;
			dom_grp_mesh_start[n] = 0;
			dom_grp_mesh_end[n] = 0;
			dom_grp_mesh_size[n] = 0;
		}
		else {
			//	COM_DOM = 1;
			dm = dmesh;
			if (idx < pad)
				dm += 1;

			m0 = m1;
			m1 = m0 + dm;

			dom_grp_mesh_start[n] = m0;
			dom_grp_mesh_end[n] = m1;
			dom_grp_mesh_size[n] = dm;

		}
		if (PROC_RANK == dom_grp[n])
			dom_grp_rank = n;
	}
	if (dom_grp_rank < first)
		COM_DOM =0;
	else
		COM_DOM= 1;

	int mx0, mx1, mdx;	
	cpuType wid = BOXSIZE/ NSIDEMESH;

	mx0 = dom_grp_mesh_start[dom_grp_rank];
	mdx = dom_grp_mesh_size[dom_grp_rank];
	mx1 = dom_grp_mesh_end[dom_grp_rank];

	BOX_DOM_L[0] = mx0 * wid;
	BOX_DOM_L[1] = BOX0L;
	BOX_DOM_L[2] = BOX1L;

	BOX_DOM_H[0] = mx1 * wid;
	BOX_DOM_H[1] = BOX0H;
	BOX_DOM_H[2] = BOX1H;

//	printf(" > >[%d]  %d :  %lf %lf\n", PROC_RANK, dom_grp_rank, BOX_DOM_L[0], BOX_DOM_H[0]);
}


int domain_restart() {
	proc_rank_sudom = (int*) mymalloc(sizeof(int)*NSIDE0SUDOM * NSIDE1SUDOM, 1);

	int n, c, is_sudom;

	cpuType dx = BOXSIZE / NDOMinSUDOM;
	cpuType dy = BOXSIZE / NSIDE0SUDOM;
	cpuType dz = BOXSIZE / NSIDE1SUDOM;


	for (n=0; n<PROC_SIZE; n++) {
		int idx = n;

		int dside0 = NSIDE0PROC/NSIDE0SUDOM;
		int dside1 = NSIDE1PROC/NSIDE1SUDOM;

		int py = idx / NSIDE1PROC ; 
		int pz = idx % NSIDE1PROC ; 

		int sy = py / dside0;
		int sz = pz / dside1;

		int sudom_r  = sy * NSIDE1SUDOM + sz;

		if ( (0 == py % dside0) && (0 == pz % dside1) ) {
			is_sudom = 1;
			proc_rank_sudom[sudom_r] = idx;
		}
		else {
			is_sudom = 0;
		}
		if (PROC_RANK == idx) {
			SUDOM_RANK = sudom_r;
			isSUDOM = is_sudom;
		}
	}

	for (n=0, c=0; n<PROC_SIZE; n++) {
		int idx = n;

		int dside0 = NSIDE0PROC/NSIDE0SUDOM;
		int dside1 = NSIDE1PROC/NSIDE1SUDOM;

		int py = idx / NSIDE1PROC ; 
		int pz = idx % NSIDE1PROC ; 

		int sy = py / dside0;
		int sz = pz / dside1;

		int sudom_r  = sy * NSIDE1SUDOM + sz;

		if ( sudom_r == SUDOM_RANK ) {
			dom_grp[c] = n;
			c ++;
		}

		if (PROC_RANK == idx) {

			BOX0L = sy * dy;
			BOX0H = (sy + 1 )* dy;
			BOX1L = sz * dz;
			BOX1H = (sz + 1 )* dz;

			if (sy == NSIDE0SUDOM -1 )BOX0H = BOXSIZE;
			if (sz == NSIDE1SUDOM -1 )BOX1H = BOXSIZE;
		}

	}
	PROC_RANK_SUDOM = dom_grp[0];
	int mx0, mx1, mdx;	
	cpuType wid = BOXSIZE/ NSIDEMESH;

	mx0 = dom_grp_mesh_start[dom_grp_rank];
	mdx = dom_grp_mesh_size[dom_grp_rank];
	mx1 = dom_grp_mesh_end[dom_grp_rank];

	BOX_DOM_L[0] = mx0 * wid;
	BOX_DOM_L[1] = BOX0L;
	BOX_DOM_L[2] = BOX1L;

	BOX_DOM_H[0] = mx1 * wid;
	BOX_DOM_H[1] = BOX0H;
	BOX_DOM_H[2] = BOX1H;

	//	printf(" > >[%d]  %d :  %lf %lf\n", PROC_RANK, dom_grp_rank, BOX_DOM_L[0], BOX_DOM_H[0]);
}

