#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "photoNs.h"
#include "gravity.h"

#ifdef GPUP2P
#include "gpu_kernel.h"
#endif

void p2ppack_init() ;
void p2ppack_cln() ;

void particle_combine(cpuType boxL[3], cpuType boxH[3]) {
	int npart_tot;
	npart_tot = NPART + NPART_BND;



	if (npart_tot > MAXNPART) {
		printf(" error : need more buff MAXNPART: %d + %d > %d\n", NPART, NPART_BND, MAXNPART);
		exit(0);
	}


	boxL[0] = BOXSIZE;
	boxL[1] = BOXSIZE;
	boxL[2] = BOXSIZE;

	boxH[0] = 0.0;	
	boxH[1] = 0.0;	
	boxH[2] = 0.0;	

	int n, i ;
#if 0
	for ( n=0; n<NPART; n++ ) {
		part[n].tag=1;

		if (boxL[0] > part[n].pos[0]) boxL[0] = part[n].pos[0];
		if (boxL[1] > part[n].pos[1]) boxL[1] = part[n].pos[1];
		if (boxL[2] > part[n].pos[2]) boxL[2] = part[n].pos[2];

		if (boxH[0] < part[n].pos[0]) boxH[0] = part[n].pos[0];
		if (boxH[1] < part[n].pos[1]) boxH[1] = part[n].pos[1];
		if (boxH[2] < part[n].pos[2]) boxH[2] = part[n].pos[2];


	}

	for (i=0, n=NPART; i<NPART_BND; i++, n++) {
		part[n].tag = 0;
		
		part[n].pos[0] = part_b[i].pos[0];
		part[n].pos[1] = part_b[i].pos[1];
		part[n].pos[2] = part_b[i].pos[2];

		if (boxL[0] > part[n].pos[0]) boxL[0] = part[n].pos[0];
		if (boxL[1] > part[n].pos[1]) boxL[1] = part[n].pos[1];
		if (boxL[2] > part[n].pos[2]) boxL[2] = part[n].pos[2];

		if (boxH[0] < part[n].pos[0]) boxH[0] = part[n].pos[0];
		if (boxH[1] < part[n].pos[1]) boxH[1] = part[n].pos[1];
		if (boxH[2] < part[n].pos[2]) boxH[2] = part[n].pos[2];
	}

#endif

		float pos[3];
	for ( n=0; n<NPART; n++ ) {
		part[n].tag=1;
#ifndef INTXYZ
		pos[0] = part[n].pos[0];
		pos[1] = part[n].pos[1];
		pos[2] = part[n].pos[2];

#else
		pos[0] = (float)( part[n].posi[0] * INT2POS);
		pos[1] = (float)( part[n].posi[1] * INT2POS);
		pos[2] = (float)( part[n].posi[2] * INT2POS);
#endif
		if (boxL[0] > pos[0]) boxL[0] = pos[0];
		if (boxL[1] > pos[1]) boxL[1] = pos[1];
		if (boxL[2] > pos[2]) boxL[2] = pos[2];

		if (boxH[0] < pos[0]) boxH[0] = pos[0];
		if (boxH[1] < pos[1]) boxH[1] = pos[1];
		if (boxH[2] < pos[2]) boxH[2] = pos[2];


	}

	for (i=0, n=NPART; i<NPART_BND; i++, n++) {
		part[n].tag = 0;
	
#ifndef INTXYZ	
		part[n].pos[0] = part_b[i].pos[0];
		part[n].pos[1] = part_b[i].pos[1];
		part[n].pos[2] = part_b[i].pos[2];

		pos[0] = part[n].pos[0];
		pos[1] = part[n].pos[1];
		pos[2] = part[n].pos[2];


#else
		part[n].posi[0] = part_b[i].posi[0];
		part[n].posi[1] = part_b[i].posi[1];
		part[n].posi[2] = part_b[i].posi[2];

		pos[0] = (float)( part[n].posi[0] * INT2POS);
		pos[1] = (float)( part[n].posi[1] * INT2POS);
		pos[2] = (float)( part[n].posi[2] * INT2POS);
#endif
		if (boxL[0] > pos[0]) boxL[0] = pos[0];
		if (boxL[1] > pos[1]) boxL[1] = pos[1];
		if (boxL[2] > pos[2]) boxL[2] = pos[2];

		if (boxH[0] < pos[0]) boxH[0] = pos[0];
		if (boxH[1] < pos[1]) boxH[1] = pos[1];
		if (boxH[2] < pos[2]) boxH[2] = pos[2];
	}






	boxL[0] -= 0.000001 * BOXSIZE;
	boxL[1] -= 0.000001 * BOXSIZE;
	boxL[2] -= 0.000001 * BOXSIZE;

	boxH[0] += 0.000001 * BOXSIZE;
	boxH[1] += 0.000001 * BOXSIZE;
	boxH[2] += 0.000001 * BOXSIZE;


	NPART = npart_tot;

} 

void particle_reduce() {
	int n, c;

	Body tpart;
	int p0, p1;

	p0 = 0;
	p1 = NPART-1; 
	while (p0 < p1 ) {
		if ( 0==part[p0].tag ) {
			tpart = part[p0];
			while ( 0==part[p1].tag  ) {
				p1--; 
			}
			if ( p1 <= p0)
				break;

			part[p0] = part[p1];
			part[p1] = tpart;

		}
		p0++;
	}


	NPART = NPART - NPART_BND;

	for (n=0; n<NPART; n++ ) {
		if (part[n].tag != 1) {
			printf(" p erp\n");
			exit(0);
		}

	}
	for (; n<NPART+NPART_BND; n++ ) {
		if (part[n].tag != 0) {
			printf(" p erp 2\n");
			exit(0);
		}

	}


	//printf(" npr =  %d %d %d \n", part[NPART-1].act, part[NPART].act, NPART);
}

void p2ppack_init() {
	NPACK = 0;
	MAXNPACK = pow(2, packLevel) + 1;
	if (PROC_RANK == 0) printf("MAXNPACK %d\n", MAXNPACK);
	packlist = (P2PPack*)malloc(sizeof(P2PPack) * MAXNPACK);

}

void p2ppack_cln() {
#ifdef INTRA_LB
	int i;
	for (i=0; i<NPACK; ++i) {
		free(packlist[i].task_t);
		free(packlist[i].task_s);
	}
#endif	
	free(packlist);
	packlist = NULL;
	NPACK = 0;

}

#if 0
void gravity_act(int active_level) {

#ifdef ACT_SYNC
	acc_cln();
	gravity_cln();
	boundary();

	if (COM_DOM) {
		int direct = 0;
		cpuType bdl[3];
		cpuType bdr[3];
		particle_combine(bdl, bdr);
		fmm_prepare(direct, bdl, bdr);
	}
#endif


	if (COM_DOM) {
		int n;
		acc[0] = (cpuType*)mymalloc(sizeof(cpuType)*NPART, 14);
		acc[1] = (cpuType*)mymalloc(sizeof(cpuType)*NPART, 15);
		acc[2] = (cpuType*)mymalloc(sizeof(cpuType)*NPART, 16);

		for (n=0; n<NPART; n++) {
			acc[0][n] = 0.0;
			acc[1][n] = 0.0;
			acc[2][n] = 0.0;
		}

		fmm_task_act(active_level);

	}

}
#endif


void gravity_fmm_act(int level, int max_level ) {
	int n;
	double t0, t1, t2, t3, t4;
	t0 = dtime();
	
	int syn_act = 1;

#if 1
	if (level < max_level || level < 4) {
		syn_act = 1;
	}
	else 
		syn_act = 0;

#else


	if (level ==  max_level &&  leveltot[max_level] < 3000  ) {
		syn_act = 0;
	}
#endif


if ( syn_act ) {

//	if (1==PROC_RANK ) printf(" syn %d %d\n", level, max_level);
	acc_cln();
	gravity_cln();
	boundary();

	if (COM_DOM) {
		int direct = 0;
		cpuType bdl[3];
		cpuType bdr[3];

		particle_combine(bdl, bdr);

		fmm_prepare(direct, bdl, bdr);
	}
}


	if (COM_DOM) {
		p2ppack_init();
		fmm_task_act(level);
	}
	if (COM_DOM) {


if ( syn_act ) {
		acc[0] = (cpuType*)mymalloc(sizeof(cpuType)*NPART, 14);
		acc[1] = (cpuType*)mymalloc(sizeof(cpuType)*NPART, 15);
		acc[2] = (cpuType*)mymalloc(sizeof(cpuType)*NPART, 16);
}

		for (n=0; n<NPART; n++) {
			acc[0][n] = 0.0;
			acc[1][n] = 0.0;
			acc[2][n] = 0.0;
		}

		for (n = first_leaf; n < last_leaf; n++)
			if (leaf[n].active >= level)	
				l2p_act(n,level);	
	}
#if defined(GPUP2P)
	if (COM_DOM) {
		dtime_task += get_last_kntime(PROC_RANK);
		t1 = dtime();
		hostmalloc_out(NPART);
		devMemcpyDH_updateAcc(packlist, NPACK, -1);
		hostfree_out();
		p2ppack_cln();
		dtime_task += dtime() - t1;
	}
#endif
}

void gravity_fmm_pm_kick(cpuType dkick) {
	int n;
	double t0, t1, t2, t3, t4;
	t0 = dtime();

	if (37 == PROC_RANK) printf("gravity_fmm_pm_kcik %f\n", dkick);

	if (COM_DOM) {
		int direct = 0;
		cpuType bdl[3];
		cpuType bdr[3];

		particle_combine(bdl, bdr);

		fmm_prepare(direct, bdl, bdr);
	}

	
	if (COM_DOM ) {
		p2ppack_init();
		fmm_task();

	}

	int ALLOC_ACC_MEM = 1;
	partmesh_kick(dkick, ALLOC_ACC_MEM);

	if (COM_DOM) {
		for (n = first_leaf; n < last_leaf; n++) {
			int p;
			int ip = leaf[n].ipart;
			int np = leaf[n].npart;
			cpuType accp[3];
			int cntp = 0;
			accp[0] = accp[1] = accp[2] = 0.0;	
			for (p=ip; p<np+ip; p++) {			
				if (part[p].tag == 1) {
					accp[0] += acc[0][p];
					accp[1] += acc[1][p];
					accp[2] += acc[2][p];
					cntp++;	
				}
			}
			leaf[n].acc_pm[0] = accp[0]/(cntp);
			leaf[n].acc_pm[1] = accp[1]/(cntp);
			leaf[n].acc_pm[2] = accp[2]/(cntp);
		}

#if 0 
///////
		if ( 1== PROC_RANK ) {


#ifndef INTXYZ
			FILE *ff = fopen("pos_acc.txt", "w");
#else
			FILE *ff = fopen("posi_acc.txt", "w");
#endif
			int cnt = 0;
			for (n=0; n<NPART; n++) {
				float pp[3];

#ifndef INTXYZ
				pp[0] = part[n].pos[0];		
				pp[1] = part[n].pos[1];		
				pp[2] = part[n].pos[2];		
#else
				pp[0] = part[n].posi[0] * INT2POS;
				pp[1] = part[n].posi[1] * INT2POS;
				pp[2] = part[n].posi[2] * INT2POS;

#endif
				if ( -1000.0< pp[0] && pp[0] < 2000.0 ) 
				if ( -1000.0< pp[1] && pp[1] < 2000.0 ) 
				if ( -1000.0< pp[2] && pp[2] < 2000.0 ) 
				{
				fprintf(ff, "%f %f %f %e %e %e\n", pp[0], pp[1], pp[2], acc[0][n],acc[1][n],acc[2][n]);
				cnt++;
				}
			}
		}
//////////////////////			

#endif









		for (n=0; n<NPART; n++) {
			acc[0][n] = 0.0;
			acc[1][n] = 0.0;
			acc[2][n] = 0.0;
		}


		for (n = first_leaf; n < last_leaf; n++)
			if (leaf[n].active >= 0)	
				l2p(n);	








	}


#if defined(GPUP2P)
	if (COM_DOM) {
		dtime_task += get_last_kntime(PROC_RANK);
		t1 = dtime();
		hostmalloc_out(NPART);
		devMemcpyDH_updateAcc(packlist, NPACK, -1);
		hostfree_out();
		p2ppack_cln();
		dtime_task += dtime() - t1;



#if 0
///////
		if ( 1== PROC_RANK ) {


#ifndef INTXYZ
			FILE *ff = fopen("pos_acc.txt", "w");
#else
			FILE *ff = fopen("posi_acc.txt", "w");
#endif
			int cnt = 0;
			for (n=0; n<NPART; n++) {
				float pp[3];

#ifndef INTXYZ
				pp[0] = part[n].pos[0];		
				pp[1] = part[n].pos[1];		
				pp[2] = part[n].pos[2];		
#else
				pp[0] = part[n].posi[0] * INT2POS;
				pp[1] = part[n].posi[1] * INT2POS;
				pp[2] = part[n].posi[2] * INT2POS;

#endif
				if ( -1000.0< pp[0] && pp[0] < 2000.0 ) 
				if ( -1000.0< pp[1] && pp[1] < 2000.0 ) 
				if ( -1000.0< pp[2] && pp[2] < 2000.0 ) 
				{
				fprintf(ff, "%f %f %f %e %e %e\n", pp[0], pp[1], pp[2], acc[0][n],acc[1][n],acc[2][n]);
				cnt++;
				}
			}
		}
//////////////////////			

#endif

	}
//MPI_Barrier(MPI_COMM_WORLD);
//exit(0);

#endif
	TAG_FULL_GRAVITY = 1;
}


void acc_cln(){
	if (COM_DOM) {
		myfree(acc[0], 14);
		myfree(acc[1], 15);
		myfree(acc[2], 16);
	}
}

/*
void acc_pm_cln(){
	if (COM_DOM) {
		myfree(acc_pm[0], 110);
		myfree(acc_pm[1], 111);
		myfree(acc_pm[2], 112);
	}
}
*/

void gravity_cln(){

	if (COM_DOM) {
		fmm_deconstruct();
		particle_reduce();
	}
}

void gravity_cln_analysis(){

	if (COM_DOM) {
		particle_reduce();
	}
}

