#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "photoNs.h"
#include "step.h"

	
long long leveltot[MAX_STEP_LEVEL+1];

int active_label_ex(cpuType loga_i, cpuType loga_f){

	int gmax_level = 0;
	int n, lv, max_level;
	long long level[MAX_STEP_LEVEL+1];
//	long long leveltot[MAX_STEP_LEVEL+1];

	for (n=0; n<= MAX_STEP_LEVEL; n++) {
		level[n] = 0;
		leveltot[n] = 0;
	}

	max_level = 0;
	if (COM_DOM) {
		cpuType eta = 0.025;
		
		eta = 0.035;

		cpuType tp, ac, ax,ay,az;
		cpuType ai, af, a;
		ai = exp(loga_i);
		af = exp(loga_f);
		a = 0.5*(ai+af); 
		cpuType fac = GravConst / (a*a);

		cpuType dstep = t_flat_lcdm_a(af) - t_flat_lcdm_a(ai);

		max_level =0;

		//#pragma omp parallel for num_threads(NumThread)

		float tp_max = 0.0;
		int nl;
		for (nl = first_leaf; nl < last_leaf; nl++) { 

			int ip = leaf[nl].ipart; 
			int np = leaf[nl].npart;
			cpuType accp[3];

			accp[0] = leaf[nl].acc_pm[0];
			accp[1] = leaf[nl].acc_pm[1];
			accp[2] = leaf[nl].acc_pm[2];
			for (n=ip; n<ip+np; n++) {
				if (part[n].tag != 1) continue;
				ax = acc[0][n] + accp[0];
				ay = acc[1][n] + accp[1];
				az = acc[2][n] + accp[2];

				ac = fac*sqrt(ax*ax + ay*ay + az*az);

				tp = sqrt(2.0*eta*SoftenLength*a/ac);
				lv = 0;

	//			if (tp > tp_max) tp_max = tp;

				while (tp < dstep/(1<<lv) ) {
					lv ++;
				}
#ifdef FIXEDSTEP
				lv = 0;
#endif
				if (lv > MAX_STEP_LEVEL) {
					lv = MAX_STEP_LEVEL;
				}

	//			if (507 == n)printf(" %d part : %f %f %f, %f %f %f\n",n,part[n].posi[0]*INT2POS, part[n].posi[1]*INT2POS, part[n].posi[2]*INT2POS, part[n].pos[0], part[n].pos[1], part[n].pos[2]);
				//if (lv==1)printf(" %d part : %d %d %d, %f %f %f\n",n,part[n].posi[0], part[n].posi[1], part[n].posi[2], part[n].pos[0], part[n].pos[1], part[n].pos[2]);

				part[n].act = lv;
				if (part[n].act > max_level)
					max_level = part[n].act;
			}
		}
		//		printf(" label : [%d]  max_level = %d %f %f\n",PROC_RANK,  max_level, tp, dstep);
	//	printf("[%d] tp = %f\n", PROC_RANK, tp_max);
		for (n=0; n<NPART; n++) {
			if (part[n].tag != 1) continue;
			level[part[n].act] ++;
		}
	}

	MPI_Allreduce(&max_level, &gmax_level,1,MPI_INT, MPI_MAX, MPI_COMM_WORLD );
	MPI_Allreduce(level, leveltot, MAX_STEP_LEVEL+1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD );


	if (0==PROC_RANK) {
		int	n=0; 
		for(n=0; n<=gmax_level; n++) {
			if (n == gmax_level)			
				printf("   Finest level %d : %lld\n", n, leveltot[n]); 
			else
				printf("                %d : %lld\n", n, leveltot[n]);
			fflush(stdout);
		}

	}

	return gmax_level;
}


void kick_half_step(cpuType dkh_pm) {
	int n;

	if (COM_DOM) {
		for (n=0; n<NPART; n++) {
			if (part[n].tag == 1) {
				part[n].vel[0] += acc_pm[0][n]*dkh_pm;
				part[n].vel[1] += acc_pm[1][n]*dkh_pm;
				part[n].vel[2] += acc_pm[2][n]*dkh_pm;

			}
		}

	}
}

void kick_half_active(cpuType dkh[], int level) {
	int n;
	if (COM_DOM) {
		for (n=0; n<NPART; n++) {
			if (part[n].tag == 1) {

				if (part[n].act >= level) {	
					cpuType dkh_p = dkh[ part[n].act ] ;

					part[n].vel[0] += acc[0][n] * dkh_p;
					part[n].vel[1] += acc[1][n] * dkh_p;
					part[n].vel[2] += acc[2][n] * dkh_p;
				}

			}
		}
	}			
}






void kick_half_act(cpuType dkh[], int level) {
	int n;

	if (COM_DOM) {
		for (n=0; n<NPART; n++) {
			if (part[n].tag == 1) {
				if (0 == level && part[n].act >= 0 ) {	
					cpuType dkh_p = dkh[ 0 ] ;
					part[n].vel[0] += acc_pm[0][n]*dkh_p;
					part[n].vel[1] += acc_pm[1][n]*dkh_p;
					part[n].vel[2] += acc_pm[2][n]*dkh_p;
				}


#ifndef TEST_PM
				if (part[n].act >= level) {	
					cpuType dkh_p = dkh[ part[n].act ] ;

					part[n].vel[0] += acc[0][n] * dkh_p;
					part[n].vel[1] += acc[1][n] * dkh_p;
					part[n].vel[2] += acc[2][n] * dkh_p;
				}
#endif

			}
		}

	}
}

void drift_step(cpuType dd) {
	int n;

//	if( PROC_RANK == 2) printf( " d_step = %f \n" , dd  );

	int npart_disp[32];
	for (n=0; n<32; n++) npart_disp[n] = 0;

	
	if (COM_DOM) {
		for (n=0; n<NPART; n++) {
			if (part[n].tag == 1) {

#ifndef INTXYZ

#ifdef TEST_FLOAT
				float dispf = 2000000.0;
	
				part[n].pos[0] += dispf;	
				part[n].pos[1] += dispf;	
				part[n].pos[2] += dispf;	
#endif

				part[n].pos[0] +=  part[n].vel[0] * dd;
				part[n].pos[1] +=  part[n].vel[1] * dd;
				part[n].pos[2] +=  part[n].vel[2] * dd;
#ifdef TEST_FLOAT
				part[n].pos[0] -= dispf;	
				part[n].pos[1] -= dispf;	
				part[n].pos[2] -= dispf;	
#endif




#else
				int tmp[3];
				double dst[3];
				int m;
				for (m=0; m<3; m++) {
					dst[m] = part[n].vel[m] * dd * POS2INT;
					if (dst[m] > 0.0 ) 
						dst[m] += 0.5;
					else 
						dst[m] -= 0.5;
				
					part[n].posi[m] += (int) ( dst[m] );

					tmp[m] = (int) ( dst[m] );
					if (tmp[m] < 0) tmp[m] = -tmp[m];
				}
				int tmax = tmp[0] ;
				if (tmax < tmp[1]) tmax = tmp[1];
				if (tmax < tmp[2]) tmax = tmp[2];

					int ll = 0;
					if (tmax == 0 ) 
						npart_disp[0] ++;
					else {
						while ( tmax > 1<<ll ) ll++;
					 	npart_disp[ll] ++;
					}
//				if (100000<n && n<100010) printf("> %d %d,  %f %f %f, dd=%f\n", tmax, ll, dst[0], dst[1], dst[2], dd);


			
			//	if (n <5) printf(" dst(%d) %d, %d, %d \n",n, (int)(dst[0]), (int) dst[1], (int) dst[2] );

			//	part[n].posi[0] += (int) (part[n].vel[0]*dd*POS2INT);
			//	part[n].posi[1] += (int) (part[n].vel[1]*dd*POS2INT);
			//	part[n].posi[2] += (int) (part[n].vel[2]*dd*POS2INT);
#endif
			}
		}
//		if ( 1 ==PROC_RANK)
//		for (n=0; n<16; n++)  
//			printf("[%d] %d, %d\n", PROC_RANK, n, npart_disp[n]);
	}

}


