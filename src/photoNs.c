#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "typedef.h"
#include "photoNs.h"
#include "setparam.h"
#include "step.h"
#include "snapshot.h"
#include "gravity.h"
#ifdef GPUP2P
#include "gpu_kernel.h"
#endif
#include "measure.h"
#include "schedule.h"

#ifdef HALO_FOF
void analysis_fof(double redshift, int min_npart_in_halo, double link_length, int timestamp) ;
#endif

void analysis(double redshift, int min_npart_in_halo, double link_length, int timestamp) ;

int grav_level_active(int step, int level_max) {
	int level = 0;
	int nstep = step ;
	nstep = step;
	int tra =  level_max -1 ;
	while (tra >= 0) {
		if (  0 == nstep % ( 1 << tra ) ) {
			level = level_max - tra;
			break;
		}
		tra --;		
	}
	if (nstep == 0 ) level = 0 ;

	return level;
}

void analysis_driver(int snap_idx, int nth, char path_snap[]){
	int npart_min = 20;

	if ( 0 == PROC_RANK)
		printf(" Read snapshot from %s/snapshot_%d \n\n", path_snap, snap_idx);

	initial_restart(path_snap, snap_idx) ;
	cpuType zz = 1.0/InitialTime - 1.0;

	MPI_Barrier(MPI_COMM_WORLD);

	if ( 0 == PROC_RANK ) {
		printf("  MASSPART = %lf [10^10Msun/h]\n", MASSPART);
		printf("  BOX SIZE = %lf [kpc/h]\n", BOXSIZE);
		printf("  OMEGA_M = %lf\n", OmegaM0);
		printf("  OMEGA_X = %lf\n", OmegaX0);
		printf("  HUBBLE0 = %lf\n", Hubble0);
		printf("  Rsoften = %lf [kpc/h]\n", SoftenLength);
		printf("* red-shift = %lf \n", zz);

		printf("  NSIDEMESH = %d\n", NSIDEMESH);
		printf("  NPART TOTAL = %lu NDOM_COM=%d\n", NPART_TOTAL, NDOM_COM);
		printf("  N Proc = %d (%d x %d)\n", PROC_SIZE, NSIDE0PROC, NSIDE1PROC);

	}


	domain_decomp();

	boundary();

	power_spec_daub( snap_idx);
		power_spec( snap_idx);

#ifdef HALO_FOF
	double linking_length = 0.2;
	double Redshift_Time;

	double tfof = dtime();
//		analysis_fof( 0.0, npart_min,  linking_length, snap_idx );
//	analysis( zz, npart_min,  linking_length, snap_idx );
//	analysis( zz, npart_min,  linking_length, nth );
	tfof = dtime() - tfof;

	printf("%5d analysis time %lf\n", PROC_RANK, tfof);

	MPI_Barrier(MPI_COMM_WORLD);
#else

	if (0==PROC_RANK)
		printf("\n WARNING! :  Turn on the OPTION 'HALO_FOF' !!! \n\n");
#endif

	if (COM_DOM ) {
		myfree(part, 9);
		myfree(part_b, 10);
	}


}




void driver(int Nstep, cpuType ai, cpuType af, int rank_snap, double walltime) ;
void driver_mem(int Nstep, cpuType ai, cpuType af, int loop_start, double walltime) ;
void make_title();

int comp(const void*a,const void*b)
{
	return ((Body*)a)->act-((Body*)b)->act;
}

int print_msg() {
	int nmsg = 57;
	return PROC_RANK % ((PROC_SIZE + nmsg - 1) / nmsg) == 0;
}


static double AI, AF;
static int step_total;
static int snap_idx;

int main(int argc, char* argv[]) 
{
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &PROC_SIZE);
	MPI_Comm_rank(MPI_COMM_WORLD, &PROC_RANK);

	make_title();

	create_mpi_datatype();

	MPI_Barrier(MPI_COMM_WORLD);

	NUM_OUTPUT_SCHEDUE = 1;

	MODE_IC = 2;
//		MODE_IC = 1;
		
	int Nstep = 3;

	sprintf(PathSnapshot ,"../data/");
	sprintf(InputSnapshot,"../data/");


	//	sprintf(InputSnapshot, "/public/hyper_test/PhotoNs367/");
	//	sprintf(InputSnapshot, "/public/hyper_test/PhotoNs367/");

	if (0 == PROC_RANK)
		printf(" PathSnapshot = %s\n", PathSnapshot);


	AI = 1.0/(127.0 + 1.0);
	AF = 1.0/(99.0 + 1.0);

	double ai = AI;
	double af = AF;

//	setup_plank_param(0.02242, 0.11933, 0.66771e-03, 0.0, 67.66, 0.8102, 0.9665,  0.6766 *  46700.0 ) ;
//	setup_plank_param(0.02242, 0.11933, 0.66771e-03, 0.0, 67.66, 0.8102, 0.9665,  2500000.0 ) ;
	setup_plank_param(0.02242, 0.11933, 0.66771e-03, 0.0, 67.66, 0.8119, 0.9665,  2500000.0*NSIDEMESH/16128.0 ) ;

	int loop_start = 0;
	double walltime_total=0.0;

#ifdef MEMLOG
	initialize_memlogfile();	
#endif

	if ( 2 ==  MODE_IC ) {

		double it0, it1;
		it0 = dtime();
		initialize_ic();
		domain_setup();
		loop_start = 0;
		step_total = 0;
		walltime_total = 0.0;

		ic_ZA(ai, NSIDEMESH);

		it1 = dtime();
		if ( 0 == PROC_RANK || 1 == PROC_RANK) {
			printf(" [%d] ic generation complete in %lf SEC !\n ", PROC_RANK, it1 - it0);
		}

	}

	if ( 1 ==  MODE_IC ) {
		//read  IC snapshot from gadget format

		initialize_ic();
		domain_setup();
		loop_start = 0;
		step_total = 0;
		double tf = dtime();
		char fpath[256];
		int nfile = 64;
//		sprintf(fpath, "/public/hyper_test/Gadget3/ICs/Scaling_Simu_n768_b125_pl13_IC_2lpt_as_prim_ind_1.0");
	//	sprintf(fpath, "/home/qwang/develop-3.7/ICs/ic_384");
		sprintf(fpath, "/public/cmp768/ICs/NS_Simu_125Mpc_768Part_Planck18_T1_ICs");

#ifdef GADGETFILE
		load_gadget_ics(fpath, nfile) ;
#else
		printf(" LACK OF GADGETFILE option !!!\n");
		exit(0);
#endif

		ai = InitialTime;
		MPI_Barrier(MPI_COMM_WORLD);

		if ( 0 == PROC_RANK)
			printf(" read Gadget IC snap from %s\n", fpath );

		cpuType zz = 1.0 /ai - 1.0;

	}

	if ( 0 ==  MODE_IC ){
//		loop_start = 2499;

		if ( 0 == PROC_RANK ) {
			FILE *fdr = fopen("./restart.txt","r");
			if (fdr == NULL){ printf(" canot open restart file\n"); exit(0);}
			fscanf(fdr, "%d %d %d %lf", &loop_start, &snap_idx, &step_total, &walltime_total);
			fclose(fdr);
			printf("\n\n restart from snapshot  %d\n\n", loop_start);

		}
		MPI_Bcast(&loop_start, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&snap_idx, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&step_total, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&walltime_total, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		//		walltime = walltime_total;

		double tf = dtime();
		initial_restart(InputSnapshot, snap_idx) ;
		ai = InitialTime;
		MPI_Barrier(MPI_COMM_WORLD);

		if ( 0 == PROC_RANK ) {
			printf("\n* PARAMETERS:\n");
			printf(" MASSPART = %e [Msun/h]\n", MASSPART * 1.0e10);
			printf(" BOX SIZE = %lf [kpc/h]\n", BOXSIZE);
			printf(" OMEGA_M = %lf\n", OmegaM0);
			printf(" OMEGA_X = %lf\n", OmegaX0);
			printf(" HUBBLE0 = %lf\n", Hubble0);
			printf(" Rcutoff = %lf\n", cutoffRadius);
			printf(" Rsoften = %lf\n", SoftenLength);

			printf(" NSIDEMESH = %d\n", NSIDEMESH);
			printf(" NPART TOTAL = %lu NDOM_COM=%d\n", NPART_TOTAL, NDOM_COM);
			printf(" N Proc = %d (%d x %d)\n", PROC_SIZE, NSIDE0PROC, NSIDE1PROC);
			printf(" init Z = %lf\n", 1.0/ai-1.0);
			printf(" ai = %lf\n", ai);
			printf(" af = %lf\n", af);

		}
		if ( 0 == PROC_RANK)
			printf(" restart from %s/snapshot_%d at z=%f time: %lf [SEC]\n", InputSnapshot, snap_idx, 1.0/ai - 1.0, dtime() - tf);


	}


	if ( 5 ==  MODE_IC ){

		int snap_idx = 2000;
		char path_snap[256];
		sprintf(path_snap, "/public/P7/i3716/");

		sprintf(PathSnapshot, "/public/P7/i3716/");

//		sprintf(path_snap, "/public/simulation_4T/");

		analysis_driver(snap_idx, snap_idx, path_snap);


		if ( 0 == PROC_RANK ) {
			printf("\n\n FOF halo find COMPLETE ! snap idx = %d \n\n", snap_idx );
		}
		MPI_Finalize();

		return 0;

	}



	if ( 6 ==  MODE_IC ){
		int n, Nsnap;
		int snap_idx[1000];

		sprintf(PathSnapshot, "/public/P7/f3710/");
		if ( 0 == PROC_RANK ) {
			char fname[256];
			sprintf(fname, "%s/output_schedue.txt", PathSnapshot);

			FILE *fd = fopen(fname, "r");
			if (fd == NULL) {
				printf(" error \n");
			}
	
			int idx, lp;
			float zd, ad, td;
			char c_idx[16], c_loop[16], c_z[16], c_a[16], c_t[16];
			fscanf(fd, "%s %s %s %s %s", c_idx, c_loop, c_z, c_a, c_t);
		
			n = 0;
			while ( ! feof(fd) ) {
				
				fscanf(fd, "%d %d %f %f %f", &idx, &lp, &zd, &ad, &td);
				printf("  %d %d %f %f %f\n", idx, lp, zd, ad, td);
				snap_idx[n] = lp;
				
				n++;
			}
			Nsnap = n-1;
			printf("  - %d - %d\n", PROC_RANK, Nsnap);
			fclose(fd);	
		}

		MPI_Bcast( &Nsnap, 1, MPI_INT, 0, MPI_COMM_WORLD );
		MPI_Bcast( snap_idx, Nsnap, MPI_INT, 0, MPI_COMM_WORLD );

		MPI_Barrier(MPI_COMM_WORLD);

		for (n=1; n<=Nsnap; n++) {	
	//		printf(" [%d] %d snap = %d\n", PROC_RANK,n, snap_idx[n]);	
			analysis_driver(snap_idx[n], n,  PathSnapshot);
		}
		if ( 0 == PROC_RANK ) {
			printf("\n\n FOF halo find COMPLETE ! snap idx = %d \n\n", snap_idx );
		}
		MPI_Finalize();

		return 0;

	}

	MPI_Barrier(MPI_COMM_WORLD);

#ifdef TASKLOG
	initialize_tasklogfile();	
#endif

	measure_len_task = 3 * NPART / 2 ; // 512
	measure_len_part = 3 * NPART / 2 ;

	if (0 == PROC_RANK) {
		printf("\n\n");
		printf(" MASSPART = %e  [Msun/h]\n", MASSPART *1.0e10);
		printf(" BOX SIZE = %lf [kpc/h]\n", BOXSIZE);
		printf(" OMEGA_M = %lf\n", OmegaM0);
		printf(" OMEGA_X = %lf\n", OmegaX0);
		printf(" HUBBLE0 = %lf\n", Hubble0);
		printf(" Rsoften = %lf [kpc/h]\n", SoftenLength);
		printf(" NPART TOTAL = %lu \n\n", NPART_TOTAL);
		printf(" NDOM_COM=%d\n", NDOM_COM);
		printf(" NSIDEMESH = %d\n", NSIDEMESH);
		printf(" Rcutoff = %lf\n", cutoffRadius);
		printf(" N Proc = %d (%d x %d)\n", PROC_SIZE, NSIDE0PROC, NSIDE1PROC);
		printf(" red-shift = %lf\n", 1.0/ai-1.0);
		printf(" ai = %lf\n", ai);
		printf(" af = %lf\n", af);
	}


	create_schedule(Nstep, AI, AF, 1.0/(1.0+30.0), AF, NUM_OUTPUT_SCHEDUE);

	driver_mem( Nstep, ai, af, loop_start, walltime_total);


	free_schedule();

	MPI_Barrier(MPI_COMM_WORLD);

	if (COM_DOM ) {
		myfree(part, 9);
		myfree(part_b, 10);
#ifdef GPUP2P
		devfree();
		dinfofree();
#endif
	}

#ifdef MEMLOG
	finalize_memlogfile();	
#endif

#ifdef TASKLOG
	finalize_tasklogfile();	
#endif

	MPI_Barrier(MPI_COMM_WORLD);

	if ( 0 == PROC_RANK ) {
		printf("\n\n  * * *  * * *  PHOTONS-3 COMPLETE !  * * *  * * *\n\n");
		printf("  *  MASSPART = %e [Msun/h]\n", MASSPART * 1.0e10  );
		printf("  *  BOX SIZE = %lf [kpc/h]\n", BOXSIZE);
		printf("  *  OMEGA_M  = %lf\n", OmegaM0);
		printf("  *  OMEGA_X  = %lf\n", OmegaX0);
		printf("  *  HUBBLE0  = %lf\n", Hubble0);
		printf("  *  Rsoften  = %lf [kpc/h]\n\n", SoftenLength);
	

	}

	MPI_Finalize();

	return 0;
}

void driver_mem(int Nstep, cpuType ai, cpuType af, int loop_start, double walltime) {

	TAG_FULL_GRAVITY = 0;
	double walltime_tot = walltime;
	int n, m;

	double tw0 ,tw1; 
	tw0= dtime();
	double tt0, tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8, tt3h;
	double dt_dom, dt_bnd, dt_gr, dt_adptv;
	tt0 = dtime();
	domain_decomp();
	MPI_Barrier(MPI_COMM_WORLD);
	tt1 = dtime();

	if (0==PROC_RANK)	
	printf(" domain_decomp \n");

	boundary();
	tt2 = dtime();
	dtime_task = 0.0;

	MPI_Barrier(MPI_COMM_WORLD);
	if (0==PROC_RANK)	
	printf(" boundary \n");
	walltime_tot += (dtime() - tw0);

#if 0
		if ( 1== PROC_RANK) {


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

				if ( 0.0< pp[0] && pp[0] < 2000.0 ) 
//				if ( 0.0< pp[1] && pp[1] < 2000.0 ) 
//				if ( 0.0< pp[2] && pp[2] < 2000.0 ) 
				{
					fprintf(ff, "%f %f %f\n", pp[0], pp[1], pp[2]);
					cnt++;
				}
			}
		}
#endif


	if (  0 == loop_start) {
		power_spectrum( loop_start );
	}



	tw0 = dtime();


	int loop, loop0, loop1, loop2;
	int level_max ;
	int nstep_active;
	loop = loop_start;

	loop0 = loop;	
	loop1 = loop + 1;
	loop2 = loop + 2;

	cpuType dloga =  (log(AF) - log(AI) ) / (Nstep  );
	cpuType loga_i, loga_m, loga_f, loga_2m, loga_2f;

	loga_i =  loop * dloga + log( AI);
	loga_m = (loop+0.5) * dloga + log( AI);
	loga_f = (loop+1.0) * dloga + log( AI);
	loga_2m= (loop+1.5) * dloga + log( AI);
	loga_2f= (loop+2.0) * dloga + log( AI);

	cpuType dk_im, dk_m2m, dk_mf, dk_f2m;

	dk_im = kick_loga(loga_i, loga_m) * GravConst;
	dk_mf = kick_loga(loga_m, loga_f) * GravConst;
	dk_f2m = kick_loga(loga_f, loga_2m) * GravConst;
	dk_m2m = kick_loga(loga_m, loga_2m) * GravConst;

	cpuType dKick1Half[MAX_STEP_LEVEL];
	cpuType dKick2Half[MAX_STEP_LEVEL];
	cpuType npart_mean = (cpuType)NPART_TOTAL/NDOM_COM;

	gravity_fmm_pm_kick(dk_im);

	level_max = active_label_ex(loga_i, loga_f);

	nstep_active = (1<<level_max);
	cpuType dloga_step = (loga_f - loga_i) / nstep_active;

	MPI_Barrier(MPI_COMM_WORLD);
	tt3 = dtime();

	dt_dom = tt1 - tt0;
	dt_bnd = tt2 - tt1;
	dt_gr = tt3 - tt2; 
	if (COM_DOM && print_msg())
		printf(" %5d GRAV dt : dom %lf, bnd %lf, gr %lf (pm %lf, prep %lf m2l %lf task %lf)\n", PROC_RANK, dt_dom, dt_bnd, dt_gr, dtime_pm, dtime_prep, dtime_m2l, dtime_task);
	else if (print_msg())
		printf(" %5d GRAV dt : dom %lf, bnd %lf, pm %lf\n", PROC_RANK, dt_dom, dt_bnd, dtime_pm);
	if( 0 == PROC_RANK ) printf("\n Reference Time of single step = %lf\n\n", tt3-tt0 );

	fflush(stdout);

	int step_cnt = step_total;
	for (loop = loop_start; loop < Nstep; loop++) {

		double ttl = dtime();

		loop1 = loop + 1;


		loga_i =  loop * dloga + log( AI);
		loga_m = (loop+0.5) * dloga + log( AI);
		loga_f = (loop+1.0) * dloga + log( AI);
		loga_2m= (loop+1.5) * dloga + log( AI);
		loga_2f= (loop+2.0) * dloga + log( AI);

#ifdef MESHOUT
		if ( analysis_required(loop1) ) {
			setup_mesh_output( analysis_output_index(loop1) );

			OUTPUT_DENSITY_MESH = 1;
		}
#endif

		dk_im  = kick_loga(loga_i, loga_m) * GravConst;
		dk_mf  = kick_loga(loga_m, loga_f) * GravConst;
		dk_f2m = kick_loga(loga_f, loga_2m) * GravConst;
		dk_m2m = kick_loga(loga_m, loga_2m) * GravConst;

		cpuType	dk  = kick_loga(loga_i, loga_f);
		cpuType	dd = drift_loga(loga_i, loga_f);

		cpuType aa = exp(loga_f);
		cpuType zz = 1.0 /aa - 1.0;

		if (PROC_RANK  == 0) 
			printf("\n\nLOOP         a=( %lf to %.10lf )        %5d\n\n",  exp(loga_i), exp(loga_f), loop);

		if (PROC_RANK  == 0) 
			printf(" dk = %f dd =%f \n zz = %lf\n\n", dk*GravConst, dd, zz);

		fflush(stdout);
		cpuType dkh0, dkh1;
		dkh0 = kick_loga(loga_i, loga_m) * GravConst;
		dkh1 = kick_loga(loga_m, loga_f) * GravConst;
		int current_level = 0;
		cpuType loga_i_step, loga_m_step, loga_f_step;

		dloga_step = (loga_f - loga_i) / nstep_active;

		for (n=0; n<nstep_active; n++) {
			loga_i_step = dloga_step * n      + loga_i;
			loga_f_step = dloga_step *(n+1.0) + loga_i;	
			cpuType d_step  = drift_loga(loga_i_step, loga_f_step);

			for (m=0; m<=level_max; m++ ) {
				cpuType la_i, la_m,  la_f;

				int wid = 1<< (level_max - m);
				int mi = (n / wid) ;
				la_i = loga_i + (mi    ) * wid * dloga_step;
				la_m = loga_i + (mi+0.5) * wid * dloga_step;
				la_f = loga_i + (mi+1.0) * wid * dloga_step;

				dKick1Half[m] = kick_loga(la_i, la_m) * GravConst ;
				dKick2Half[m] = kick_loga(la_m, la_f) * GravConst ;	

			}

			kick_half_active(dKick1Half, current_level);    

			drift_step(d_step);

			if (0 == PROC_RANK) {
				printf(" step = (%d / %d) level = (%d/%d) \n", n, nstep_active, current_level, level_max);
				fflush(stdout);
			}

			current_level = grav_level_active(n+1, level_max); 

			if (n < (nstep_active-1) ) {

				gravity_fmm_act(current_level, level_max);	

				kick_half_active(dKick2Half, current_level);

			}

		}

		acc_cln();

		gravity_cln();

		tt3 = dtime();

		domain_decomp();



		tt4 = dtime();

		boundary();



#if 0
		if (  analysis_required(loop1) || loop1 == 100 ) {
#ifdef HALO_FOF
	double linking_length = 0.2;
	double Redshift_Time;
	int npart_min = 20;
	double tfof = dtime();
//		analysis_fof( zz, npart_min,  linking_length, loop1 );
	analysis( zz, npart_min,  linking_length, snap_idx );
	tfof = dtime() - tfof;

	printf("%5d analysis time %lf\n", PROC_RANK, tfof);

	MPI_Barrier(MPI_COMM_WORLD);
#else

	if (0==PROC_RANK)
		printf("\n WARNING! :  Turn on the OPTION 'HALO_FOF' !!! \n\n");
#endif

}


#endif
		if (  analysis_required(loop1) ) {
			power_spectrum( loop1 );

			if (0 == PROC_RANK) printf("output powspec_%d zz = %lf\n", loop1, zz); 
		}

		tt5 = dtime();
		if (COM_DOM) {
			if (print_msg())
				printf(" [%d] NPART = %d (r=%f), NPART_BND = %d (r=%f)\n", PROC_RANK, NPART,( (float)(NPART+NPART_BND )/(float)MAXNPART), NPART_BND, ((float)NPART_BND)/(float)MAXNPART_BND);
			measure_len_part = NPART + NPART_BND;

			cpuType ratio = 1.16;

			if ((int)(measure_len_part*ratio) > MAXNPART) {
				MAXNPART = (int)( ratio * 1.1 * measure_len_part);
				printf(" [%d] realloc mem of part!\n", PROC_RANK);
				part = (Body*) myrealloc(part, sizeof(Body)* MAXNPART, 9);
			}
			if ((int)(NPART_BND*ratio) > MAXNPART_BND) {
				MAXNPART_BND = (int)( ratio * 1.1 * NPART_BND);
				printf(" [%d] realloc mem of part_bnd!\n", PROC_RANK);
				part_b= (BodyBnd*) myrealloc(part_b, sizeof(BodyBnd)* MAXNPART_BND, 10);

				MAXNPART_BND_SUD = MAXNPART_BND * NDOMinSUDOM;

				MAXNPART_DOM = MAXNPART_BND_SUD;
				MAXNPART_DOM_SUD = MAXNPART_BND_SUD;
			}
		}


		step_cnt += (1 << level_max );



		if (  analysis_required(loop1) ) {

			gravity_fmm_pm_kick( dk_mf );
		}
		else {
			gravity_fmm_pm_kick( dk_m2m );
		}

		kick_half_active(dKick2Half, 0);


		level_max = active_label_ex(loga_f, loga_2f);
		nstep_active = (1<<level_max);


		if (  analysis_required(loop1) ) {

			acc_cln();

			double tf = dtime();
			int snap_idx = loop1;

			walltime_tot += (dtime() - tw0);

			make_dir(PathSnapshot, "snapshot", snap_idx);	
			output_snapshot_new(zz, snap_idx, PROC_RANK) ;
			if (print_msg())
				printf("[%d] Disk IO write %d time: %lf [SEC]\n", PROC_RANK, loop1, dtime() - tf);



			MPI_Barrier(MPI_COMM_WORLD);
			if ( 0 == PROC_RANK) {
				FILE *fdw = fopen("./restart.txt","w");
				if (fdw == NULL){ printf(" canot open restart file\n"); exit(0);}
				fprintf(fdw,"%d %d %d %lf", loop1, snap_idx, step_cnt, walltime_tot  );
				fclose(fdw);
			}

			tw0 = dtime();

			// 83.2
			gravity_cln();

			boundary();
			// 74.6

			power_spectrum( loop1 );


			// 92.4 -88.9
			MPI_Barrier(MPI_COMM_WORLD);	
			gravity_fmm_pm_kick( dk_f2m );

			// 83.8
		}

		tt6 = dtime();

		dt_adptv = tt3 - tt2;
		dt_dom = tt4 - tt3;
		dt_bnd = tt5 - tt4;
		dt_gr  = tt6 - tt5;

		double bandwidth = NPART_BND * sizeof(BodyBnd) * 1.0e-6 / dt_bnd;
		if (COM_DOM && print_msg())
			printf(" %5d head_r %5d dt : dom %lf, bnd %lf, adptv %lf, gr %lf (pm %lf, prep %lf m2l %lf adptv_task %lf task %lf)\n", PROC_RANK, PROC_RANK_SUDOM, dt_dom, dt_bnd, dt_adptv, dt_gr, dtime_pm, dtime_prep, dtime_m2l, dtime_adptv_task, dtime_task);
		else if (print_msg())
			printf(" %5d head_r %5d dt : dom %lf, bnd %lf, pm %lf\n", PROC_RANK, PROC_RANK_SUDOM, dt_dom, dt_bnd, dtime_pm);

		MPI_Barrier(MPI_COMM_WORLD);	

		walltime_tot += (dtime() - tw0);
		tw0 = dtime();

		if ( 0== PROC_RANK) {
			printf(" [%d] LOOP TIME = %lf ( %lf ) SEC  \n", PROC_RANK, dtime() - ttl, walltime_tot);
			printf(" STEP %d   \n\n",step_cnt);
		}

		fflush(stdout);
	}


	double zf = 1.0/af - 1.0;	

	acc_cln();
	gravity_cln();

}



void make_title(){
	if ( 0 == PROC_RANK) {
		printf("\n\n\n\n   / PHOTONS-3 / \n\n");



		printf("\n   OPTIONS\n");
#ifdef MTHKDTree 
		printf("    + MTHKTree\n");
#endif

#ifdef ICPOTGRAD
		printf("    + ICPOTGRAD\n");
#endif

#ifdef INTRA_LB
		printf("    + INTRA_LB\n");
#endif

#ifdef GPUP2P
		printf("    + GPUP2P\n");
#endif

#ifdef FIXEDSTEP
		printf("    + FIXEDSTEP\n");
#endif

#ifdef SPI
		printf("    + SPI\n");
#endif

#ifdef SPI_T2
		printf("    + SPI_T2\n");
#endif

#ifdef ACT_OPT
		printf("    + ACT_OPT\n");
#endif

#ifdef GPUCUTOFF
		printf("    + GPUCUTOFF\n");
#endif

#ifdef JSPLIT
		printf("    + JSPLIT\n");
#endif 

#ifdef COMPRESS
		printf("    + COMPRESS\n");
#endif 

#ifdef GPU_DP
		printf("    + GPU_DP\n");
#endif

#ifdef CPU_DP
		printf("    + CPU_DP\n");
#endif

#ifdef SAVE_MEM
		printf("    + SAVE_MEM\n");
#endif

#ifdef MESHOUT
		printf("    + MESHOUT\n");
#endif 

#ifdef HALO_FOF
		printf("    + HALO_FOF\n");
#endif 

#ifdef HALO_FOF_SORT
		printf("    + HALO_FOF_SORT\n");
#endif 

#ifdef SUBHALO
		printf("    + SUBHALO\n");
#endif 

#ifdef SUBHALO_SORT
		printf("    + SUBHALO_SORT\n");
#endif 

#ifdef POW_SPEC_DAUB
		printf("    + POW_SPEC_DAUB\n");
#endif 

#ifdef IOBLK
		printf("    + IOBLK\n");
#endif 

#ifdef POW_SPEC
		printf("    + POW_SPEC\n");
#endif 

#ifdef INTXYZ
		printf("    + INTXYZ\n");
#endif 

#ifdef HDF5
		printf("    + HDF5\n");
#endif 
		printf("\n   BIT = %d\n" , (int) BITWIDTH);
	}
}


