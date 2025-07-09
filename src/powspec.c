
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <limits.h>
#include "photoNs.h"

#define NUMPAD 3


#ifndef INTXYZ

void power_spec(int rank_snap){

	char fname[128];
	sprintf(fname, "%s/powspec_%d", PathSnapshot, rank_snap);
	double t0, t1, t2, t3, t4, t5;
	cpuType *smesh;
	cpuType *mbuff;
	MPI_Request req;
	MPI_Status  sta;
	MPI_Request request[NDOMinSUDOM];
	MPI_Status  status[NDOMinSUDOM];
	int n, ip;
	int meshsize_sud;
	int meshsize_proc;
	int meshsize ;
	int meshsize_pot ;

	int local_min[3], local_max[3], isize[3];
	int nrecv[NDOMinSUDOM];	

	if (  isSUDOM ) {
		int nx, ny, nz;
		nx = NSIDEMESH;
		ny = NSIDEMESH/NSIDE0SUDOM;
		nz = NSIDEMESH/NSIDE1SUDOM;
		meshsize_sud = nx * ny * nz;
		//		meshsize_proc = nx * ny * nz;
		smesh = (cpuType*)mymalloc(sizeof(cpuType)* meshsize_sud, 260) ;

		if (smesh == NULL) {
			printf(" error smesh \n");
			exit(0);
		}
		for (n=0; n<meshsize_sud; n++) 
			smesh[n] = 0.0;

	}

	if ( COM_DOM ) {	
		t0 = dtime();

		int  c ;

		local_min[0] = dom_grp_mesh_start[dom_grp_rank];
		local_max[0] = local_min[0] + dom_grp_mesh_size[dom_grp_rank];

		int idx = PROC_RANK;

		int dside0 = NSIDE0PROC/NSIDE0SUDOM;
		int dside1 = NSIDE1PROC/NSIDE1SUDOM;
		int y = idx / NSIDE1PROC ; 
		int z = idx % NSIDE1PROC ; 

		int sy = y / dside0;
		int sz = z / dside1;

		local_min[1] = sy * NSIDEMESH/NSIDE0SUDOM; 
		local_min[2] = sz * NSIDEMESH/NSIDE1SUDOM; 

		local_max[1] =(sy+1) * NSIDEMESH/NSIDE0SUDOM; 
		local_max[2] =(sz+1) * NSIDEMESH/NSIDE1SUDOM; 

		isize[0] = local_max[0] - local_min[0] ;			
		isize[1] = NSIDEMESH / NSIDE0SUDOM;
		isize[2] = NSIDEMESH / NSIDE1SUDOM;

		//		printf(" [%d]local_min (%d %d) <%d %d %d %d>\n", PROC_RANK, local_min[0], local_max[0],dom_grp[0], dom_grp[1], dom_grp[2], dom_grp[3]);

		meshsize = isize[0] * isize[1] * isize[2];
		meshsize_pot = (isize[0] + 2*NUMPAD)  * (isize[1]+2*NUMPAD) * (isize[2]+2*NUMPAD);
		pmesh = (cpuType*)mymalloc(sizeof(cpuType) * meshsize, 261);

		for (n=0; n<meshsize; n++) 
			pmesh[n] = 0.0;

		int i,j,k;
		int ii, jj, kk;
		cpuType norm = NSIDEMESH/BOXSIZE;
		cpuType delta = 1.0/norm;

		t1 = dtime();

		cpuType wi, wj , wk, win, wjn, wkn;
		c=0;	
		for (n=0; n<NPART; n++) {
			if (part[n].tag != 1) continue;

			i = (int) (part[n].pos[0] *norm) ;	
			j = (int) (part[n].pos[1] *norm) ;	
			k = (int) (part[n].pos[2] *norm) ;	

			wi = (part[n].pos[0] - (i+0.5)*delta)*norm;

			if (wi > 0) {
				ii = i + 1;
			}
			else {
				wi = -wi;
				ii = i - 1;
			}
			win =  1.0 - wi;

			wj = (part[n].pos[1] - (j+0.5)*delta)*norm;
			if (wj > 0) {
				jj = j + 1;
			}
			else {
				wj = -wj;
				jj = j - 1;
			}
			wjn =  1.0 - wj;

			wk = (part[n].pos[2] - (k+0.5)*delta)*norm;
			if (wk > 0) {
				kk = k + 1;
			}
			else {
				wk = -wk;
				kk = k - 1;
			}
			wkn =  1.0 - wk;

			i -= local_min[0] ;	
			j -= local_min[1] ;	
			k -= local_min[2] ;	

			ii -= local_min[0] ;	
			jj -= local_min[1] ;	
			kk -= local_min[2] ;	

			if (k >=0 && k< isize[2]) { 
				if (j>=0 && j<isize[1])	{
					if (i>=0 && i<isize[0]) {	
						idx = (i*isize[1] + j)*isize[2] + k;
						pmesh[idx] += MASSPART*win*wjn*wkn;	
					}

					if (ii>=0 && ii<isize[0]) {	
						idx = (ii*isize[1] + j)*isize[2] + k;
						pmesh[idx] += MASSPART*wi *wjn*wkn;
					}	
				}


				if (jj>=0 && jj<isize[1])	{
					if (i>=0 && i<isize[0]) {	
						idx = (i*isize[1] + jj)*isize[2] + k;
						pmesh[idx] += MASSPART*win*wj *wkn;	
					}
					if (ii>=0 && ii<isize[0]) {	
						idx = (ii*isize[1] + jj)*isize[2] + k;
						pmesh[idx] += MASSPART*wi *wj *wkn;	
					}
				}
			}

			if (kk>=0 && kk< isize[2]) {
				if (j>=0 && j<isize[1])	{
					if (i>=0 && i<isize[0]) {	
						idx = (i*isize[1] + j)*isize[2] + kk;
						pmesh[idx] += MASSPART*win*wjn*wk ;	
					}

					if (ii>=0 && ii<isize[0]) {	
						idx = (ii*isize[1] + j)*isize[2] + kk;
						pmesh[idx] += MASSPART*wi *wjn*wk ;	
					}
				}


				if (jj>=0 && jj<isize[1]){
					if (i>=0 && i<isize[0]) {	
						idx = (i*isize[1] + jj)*isize[2] + kk;
						pmesh[idx] += MASSPART*win*wj *wk ;	
					}
					if (ii>=0 && ii<isize[0]) {	
						idx = (ii*isize[1] + jj)*isize[2] + kk;
						pmesh[idx] += MASSPART* wi *wj *wk ;
					}
				}
			}
			c++;
		}

		t2 = dtime();

		c=0;	
		for (n=0; n<NPART_BND; n++) {
			i = (int) (part_b[n].pos[0] *norm) ;	
			j = (int) (part_b[n].pos[1] *norm) ;	
			k = (int) (part_b[n].pos[2] *norm) ;	

			if (part_b[n].pos[0] < 0.0 )
				i--;

			if (part_b[n].pos[1] < 0.0 )
				j--;

			if (part_b[n].pos[2] < 0.0 )
				k--;

			wi = (part_b[n].pos[0] - (i+0.5)*delta)*norm;

			if (wi > 0) {
				ii = i + 1;
			}
			else {
				wi = -wi;
				ii = i - 1;
			}
			win =  1.0 - wi;

			wj = (part_b[n].pos[1] - (j+0.5)*delta)*norm;
			if (wj > 0) {
				jj = j + 1;
			}
			else {
				wj = -wj;
				jj = j - 1;
			}
			wjn =  1.0 - wj;

			wk = (part_b[n].pos[2] - (k+0.5)*delta)*norm;
			if (wk > 0) {
				kk = k + 1;
			}
			else {
				wk = -wk;
				kk = k - 1;
			}
			wkn =  1.0 - wk;

			i -= local_min[0] ;	
			j -= local_min[1] ;	
			k -= local_min[2] ;	

			ii -= local_min[0] ;	
			jj -= local_min[1] ;	
			kk -= local_min[2] ;	

			//if	
			if (k >=0 && k< isize[2]) { 
				if (j>=0 && j<isize[1])	{
					if (i>=0 && i<isize[0]) {	
						idx = (i*isize[1] + j)*isize[2] + k;
						pmesh[idx] += MASSPART*win*wjn*wkn;	
					}

					if (ii>=0 && ii<isize[0]) {	
						idx = (ii*isize[1] + j)*isize[2] + k;
						pmesh[idx] += MASSPART*wi *wjn*wkn;
					}	
				}


				if (jj>=0 && jj<isize[1])	{
					if (i>=0 && i<isize[0]) {	
						idx = (i*isize[1] + jj)*isize[2] + k;
						pmesh[idx] += MASSPART*win*wj *wkn;	
					}
					if (ii>=0 && ii<isize[0]) {	
						idx = (ii*isize[1] + jj)*isize[2] + k;
						pmesh[idx] += MASSPART*wi *wj *wkn;	
					}
				}
			}

			if (kk>=0 && kk< isize[2]) {
				if (j>=0 && j<isize[1])	{
					if (i>=0 && i<isize[0]) {	
						idx = (i*isize[1] + j)*isize[2] + kk;
						pmesh[idx] += MASSPART*win*wjn*wk ;	
					}

					if (ii>=0 && ii<isize[0]) {	
						idx = (ii*isize[1] + j)*isize[2] + kk;
						pmesh[idx] += MASSPART*wi *wjn*wk ;	
					}
				}


				if (jj>=0 && jj<isize[1])	{
					if (i>=0 && i<isize[0]) {	
						idx = (i*isize[1] + jj)*isize[2] + kk;
						pmesh[idx] += MASSPART*win*wj *wk ;	
					}
					if (ii>=0 && ii<isize[0]) {	
						idx = (ii*isize[1] + jj)*isize[2] + kk;
						pmesh[idx] += MASSPART* wi *wj *wk ;
					}
				}
			}
			c++;
		}

		t3 = dtime();

		cpuType renormal = (NSIDEMESH/BOXSIZE);
		renormal = renormal*renormal*renormal;

		renormal = 1.0/MASSPART;

		c=0;
		for (i=0; i<isize[0]; i++)
			for (j=0; j<isize[1]; j++)
				for (k=0; k<isize[2]; k++) {

					pmesh[c] *= renormal;
					c++;
				}
		MPI_Isend(pmesh, meshsize, MPI_CPUTYPE, PROC_RANK_SUDOM, 100, MPI_COMM_WORLD, &req);

	}


	if ( isSUDOM ) {
		int local_size_y =  NSIDEMESH/NSIDE0SUDOM; 
		int local_size_z =  NSIDEMESH/NSIDE1SUDOM;
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			nrecv[n] = local_size_y  * local_size_z * dom_grp_mesh_size[n] ;	
		}
		ip = 0;
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			//		printf(" nrecv = %d\n", nrecv[n]);
			MPI_Recv( smesh + ip , nrecv[n], MPI_CPUTYPE, dom_grp[n], 100, MPI_COMM_WORLD, &(status[n]));
			ip += nrecv[n];
		}
	}

	if ( COM_DOM) {
		MPI_Wait(&req, &sta);		
		t4 = dtime();
		//	printf(" [%d] dtime = %lf ,mesh build %lf %lf , com=%lf tot=%lf \n", PROC_RANK, t1 -t0, t2-t1, t3-t2, t4-t3, t4-t0);
	}
	if ( COM_DOM ) {	

		myfree(pmesh, 261);
	}


	//////////////////

	if (  isSUDOM ) {
		mbuff = (cpuType*)mymalloc(sizeof(cpuType)*meshsize_sud, 262);
		int m;

		for (n=0; n<meshsize_sud; n++) 
			mbuff[n] = 0.0;

		int i,j,k;
		int ndisp[NDOMinSUDOM];
		ndisp[0] = 0;
		for (m=1; m<NDOMinSUDOM; m++) {
			ndisp[m] = ndisp[m-1] + nrecv[m-1];
		}
		int nzz = NSIDE1PROC/NSIDE1SUDOM;
		int y_n = NSIDEMESH/NSIDE0SUDOM;
		int z_n = NSIDEMESH/NSIDE1SUDOM;

		int dy = NSIDEMESH/NSIDE0PROC;
		int dz = NSIDEMESH/NSIDE1PROC;
		int nm = NSIDEMESH*dy*dz;

		cpuType *sm;
		for (m=0; m<NDOMinSUDOM; m++) {
			int x_n = dom_grp_mesh_size[m];
			// try		//	sm = smesh + m * x_n * y_n * z_n;
			sm = smesh + dom_grp_mesh_start[m] * y_n * z_n;
			int ii, jj, kk;
			int in, jn, kn;
			int ig, idx;
			for (i=0; i<x_n; i++) {
				for (j=0; j<y_n; j++) {
					for (k=0; k<z_n; k++) {
						jj = j/dy;
						kk = k/dz;
						ig = jj*nzz + kk;

						in = dom_grp_mesh_start[m] + i;
						jn = j%dy;
						kn = k%dz;
						idx = ( in * dy + jn ) * dz + kn;

						mbuff[ig*nm + idx] = sm[(i*y_n+j)*z_n+k];
					}
				}
			}
		}

		myfree(smesh, 260);
		for (m=0; m<NDOMinSUDOM; m++) {	
			MPI_Isend( mbuff + m*nm, nm, MPI_CPUTYPE, dom_grp[m], 200, MPI_COMM_WORLD, &(request[m]));
		}
	}


	int x_n = NSIDEMESH;
	int y_n = NSIDEMESH/NSIDE0PROC;
	int z_n = NSIDEMESH/NSIDE1PROC;
	meshsize_proc =  x_n  * y_n * z_n;
	mesh = (cpuType*)mymalloc(sizeof(cpuType) * meshsize_proc, 263);

	MPI_Recv(mesh, 	meshsize_proc, MPI_CPUTYPE, PROC_RANK_SUDOM, 200, MPI_COMM_WORLD, &sta); 
	if (  isSUDOM ) {
		int m;	
		for (m=0; m<NDOMinSUDOM; m++) {	
			MPI_Wait(&(request[m]), &(status[m]) );
		}
		myfree(mbuff, 262);
	}
	/// canbe free mbuff and smesh here once
	///////////////////////
	int nside[3] = {NSIDEMESH, NSIDEMESH, NSIDEMESH};
	cpuType param[2] = { splitRadius, BOXSIZE};

	FILE *fd;


	{
		int *psn = (int*)mymalloc(sizeof(int)*NSIDEMESH, 264);
		cpuType *ps = (cpuType*)mymalloc(sizeof(cpuType)*NSIDEMESH, 265);
		cpuType *ps_total;	
		for (n=0; n<NSIDEMESH; n++) {
			psn[n] = 0;
			ps[n] = 0.0;
		}
		powsk(mesh,nside,param, ps, psn);

		float *psn_total;	
		float *psn_rank = (float*)mymalloc(sizeof(float)*NSIDEMESH, 266);	
		for (n=0; n<NSIDEMESH; n++) {
			psn_rank[n] = (float)psn[n];

		}
		if (0 == PROC_RANK) {

			ps_total = (cpuType*)mymalloc(sizeof(cpuType)*NSIDEMESH, 267);
			psn_total = (float*)mymalloc(sizeof(float)*NSIDEMESH, 268);
		}

		MPI_Reduce(ps, ps_total, NSIDEMESH, MPI_CPUTYPE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(psn_rank, psn_total, NSIDEMESH, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

		if (0 == PROC_RANK ) {
			for (n=0; n<NSIDEMESH; n++) {
				if (psn_total[n] > 0)
					ps_total[n] /= (psn_total[n]);	
				else 
					ps_total[n] = 0.0;
			}

			fd = fopen(fname, "w");
			if (fd == NULL) {
				printf( " Error ! cannot open powspec file \n"  );
				exit(0);
			}
			float vol = BOXSIZE/1000.0;
			vol = vol * vol * vol;
			float Norm = vol/((float)NSIDEMESH*NSIDEMESH *NSIDEMESH);

			int output_num = (int)(NSIDEMESH*0.75);			
			output_num = NSIDEMESH/2;			
			for (n=1; n<output_num; n++) {
				fprintf(fd, "%e %e %lf\n",  (n+0.5)*2.0*M_PI/(BOXSIZE/1000.0), ps_total[n] * Norm *Norm /vol, psn_total[n]  );
			}

			fclose(fd);
		}



		if ( 0 == PROC_RANK){
			myfree(ps_total, 267);
			myfree(psn_total, 268);
		}
		myfree(psn, 264);
		myfree(psn_rank, 266);
		myfree(ps, 265);
	}
//	if (  isSUDOM ) {
//		myfree(mbuff, 262);
//	}

	//	MPI_Barrier(MPI_COMM_WORLD);
	myfree(mesh, 263);
}


#if 0
///////////////////////////// daubechies12 ////////////////////////////


#include "daub12hash.h"
void initial_daub12(void)
{
    int n, N;
    float x;
    FILE *fp;
    if (!(fp=fopen("../demo/daub12hash.dat","r")))
    {
        printf("can not open `daub12.dat'\n");
        exit(0);
    }
    N = 600001;
    DAUB12 = (float*)malloc((int)N*sizeof(float));
    for (n=0; n<N; n++)
        fscanf(fp, "%f %f", &x, &(DAUB12[n]) );
    
//    for (n=100000; n<100300; n++)
//        printf("%f  ", DAUB12[n]);
//    printf("\n");
    fclose(fp);
}

float daub12_scaling_function(float d)
{
//    if (d<0.0 || d>=6.0)
//    {
//        printf("error in daub12 %f \n", d);
//        exit(10);
//    }   
	if (0.0<d && d<6.0) 
		return DAUB12[(int)(d*100000.0)];
	else 
		return 0.0;
}


void free_daub12(void)
{
    free(DAUB12);
}

void power_spec_daub(int rank_snap)
{
	char fname[128];
	sprintf(fname, "%s/powspec_daub_%d", PathSnapshot, rank_snap);
	double t0, t1, t2, t3, t4, t5;
	cpuType *smesh;
	cpuType *mbuff;
	MPI_Request req;
	MPI_Status  sta;
	MPI_Request request[NDOMinSUDOM];
	MPI_Status  status[NDOMinSUDOM];
	int n, ip;
	int meshsize_sud;
	int meshsize_proc;
	int meshsize ;
	int meshsize_pot ;

	int local_min[3], local_max[3], isize[3];
	int nrecv[NDOMinSUDOM];	

	if (  isSUDOM ) {
		int nx, ny, nz;
		nx = NSIDEMESH;
		ny = NSIDEMESH/NSIDE0SUDOM;
		nz = NSIDEMESH/NSIDE1SUDOM;
		meshsize_sud = nx * ny * nz;
		//		meshsize_proc = nx * ny * nz;
		smesh = (cpuType*)mymalloc(sizeof(cpuType)* meshsize_sud, 260) ;

		if (smesh == NULL) {
			printf(" error smesh \n");
			exit(0);
		}

		for (n=0; n<meshsize_sud; n++) 
			smesh[n] = 0.0;
	}

	if ( COM_DOM ) {	
		t0 = dtime();

		int  c ;

		local_min[0] = dom_grp_mesh_start[dom_grp_rank];
		local_max[0] = local_min[0] + dom_grp_mesh_size[dom_grp_rank];

		int idx = PROC_RANK;

		int dside0 = NSIDE0PROC/NSIDE0SUDOM;
		int dside1 = NSIDE1PROC/NSIDE1SUDOM;
		int y = idx / NSIDE1PROC ; 
		int z = idx % NSIDE1PROC ; 

		int sy = y / dside0;
		int sz = z / dside1;

		local_min[1] = sy * NSIDEMESH/NSIDE0SUDOM; 
		local_min[2] = sz * NSIDEMESH/NSIDE1SUDOM; 

		local_max[1] =(sy+1) * NSIDEMESH/NSIDE0SUDOM; 
		local_max[2] =(sz+1) * NSIDEMESH/NSIDE1SUDOM; 

		isize[0] = local_max[0] - local_min[0] ;			
		isize[1] = NSIDEMESH / NSIDE0SUDOM;
		isize[2] = NSIDEMESH / NSIDE1SUDOM;

		//		printf(" [%d]local_min (%d %d) <%d %d %d %d>\n", PROC_RANK, local_min[0], local_max[0],dom_grp[0], dom_grp[1], dom_grp[2], dom_grp[3]);

		meshsize = isize[0] * isize[1] * isize[2];
		meshsize_pot = (isize[0] + 2*NUMPAD)  * (isize[1]+2*NUMPAD) * (isize[2]+2*NUMPAD);
		pmesh = (cpuType*)mymalloc(sizeof(cpuType) * meshsize, 261);

		for (n=0; n<meshsize; n++) 
			pmesh[n] = 0.0;

		int i,j,k;
		int ii, jj, kk;
		cpuType norm = NSIDEMESH/BOXSIZE;
		cpuType delta = 1.0/norm;



		t1 = dtime();

		cpuType wx[6], wy[6] , wz[6], wt;
		float gx, gy, gz;
		float px, py, pz;
		int Ic, Jc, Kc, Is, Js, Ks;
		int i_s, j_s, k_s;
	
		initial_daub12();

		int nside = NSIDEMESH;
		cpuType Dlt = 1.0/((cpuType) NSIDEMESH );
		c=0;	
		for (n=0; n<NPART; n++) {
			if (part[n].tag != 1) continue;
			px = part[n].pos[0] / BOXSIZE;
			py = part[n].pos[1] / BOXSIZE;
			pz = part[n].pos[2] / BOXSIZE;
			Ic = (int) (px * nside) ;	
			Jc = (int) (py * nside) ;	
			Kc = (int) (pz * nside) ;	
			gx = (float)(Ic+0.5)*Dlt;
			gy = (float)(Jc+0.5)*Dlt;
			gz = (float)(Kc+0.5)*Dlt;

			if (px>gx) Is = -2; else Is = -3;
			if (py>gy) Js = -2; else Js = -3;
			if (pz>gz) Ks = -2; else Ks = -3;

			for (i=0; i<6; i++)
				wx[i] = daub12_scaling_function( (float)( (gx-px)/Dlt + (i+Is+3) ));

			for (j=0; j<6; j++)
				wy[j] = daub12_scaling_function( (float)( (gy-py)/Dlt + (j+Js+3)) );

			for (k=0; k<6; k++)
				wz[k] = daub12_scaling_function( (float)( (gz-pz)/Dlt + (k+Ks+3)) );

			for (i=0; i<6; i++)
				for (j=0; j<6; j++)
					for (k=0; k<6; k++)  {
						i_s =(Ic+Is+i+nside)%nside;
						j_s =(Jc+Js+j+nside)%nside;
						k_s =(Kc+Ks+k+nside)%nside;          
						i_s -= local_min[0] ;	
						j_s -= local_min[1] ;	
						k_s -= local_min[2] ;	
						if (k_s >=0 && k_s< isize[2]) { 
							if (j_s>=0 && j_s<isize[1])	{
								if (i_s>=0 && i_s<isize[0]) {	

									//	grid[(i_s*nside + j_s) * nside + k_s] += (wx[i]*wy[j]*wz[k]);
									idx = (i_s*isize[1] + j_s)*isize[2] + k_s;
									pmesh[idx] += (wx[i]*wy[j]*wz[k]);
								}
							}
						}
					}
		}
		t2 = dtime();

		c=0;	
		for (n=0; n<NPART_BND; n++) {
			px = part_b[n].pos[0] / BOXSIZE;
			py = part_b[n].pos[1] / BOXSIZE;
			pz = part_b[n].pos[2] / BOXSIZE;
			Ic = (int) (px * NSIDEMESH) ;	
			Jc = (int) (py * NSIDEMESH) ;	
			Kc = (int) (pz * NSIDEMESH) ;	
		//	Ic = (int) (part[n].pos[0] *norm) ;	
		//	Jc = (int) (part[n].pos[1] *norm) ;	
		//	Kc = (int) (part[n].pos[2] *norm) ;	
			if (part_b[n].pos[0] < 0.0 )
				Ic--;

			if (part_b[n].pos[1] < 0.0 )
				Jc--;

			if (part_b[n].pos[2] < 0.0 )
				Kc--;

			gx = (float)(Ic+0.5)*Dlt;
			gy = (float)(Jc+0.5)*Dlt;
			gz = (float)(Kc+0.5)*Dlt;

			if (px>gx) Is = -2; else Is = -3;
			if (py>gy) Js = -2; else Js = -3;
			if (pz>gz) Ks = -2; else Ks = -3;

			for (i=0; i<6; i++)
				wx[i] = daub12_scaling_function( (float)( (gx-px)/Dlt + (i+Is+3) ));

			for (j=0; j<6; j++)
				wy[j] = daub12_scaling_function( (float)( (gy-py)/Dlt + (j+Js+3)) );

			for (k=0; k<6; k++)
				wz[k] = daub12_scaling_function( (float)( (gz-pz)/Dlt + (k+Ks+3)) );

			for (i=0; i<6; i++)
				for (j=0; j<6; j++)
					for (k=0; k<6; k++)  {
						i_s =(Ic+Is+i+nside)%nside;
						j_s =(Jc+Js+j+nside)%nside;
						k_s =(Kc+Ks+k+nside)%nside;          
						i_s -= local_min[0] ;	
						j_s -= local_min[1] ;	
						k_s -= local_min[2] ;	
						if (k_s >=0 && k_s< isize[2]) { 
							if (j_s>=0 && j_s<isize[1])	{
								if (i_s>=0 && i_s<isize[0]) {	

									//	grid[(i_s*nside + j_s) * nside + k_s] += (wx[i]*wy[j]*wz[k]);
									idx = (i_s*isize[1] + j_s)*isize[2] + k_s;
									pmesh[idx] += (wx[i]*wy[j]*wz[k]);
								}
							}
						}
					}

		}

//printf(" power 2 - daub12 \n");
	free_daub12();
//	rho = (double)npart/((double)nside*nside*nside);
		t3 = dtime();

		cpuType renormal = (NSIDEMESH/BOXSIZE);
		renormal = renormal*renormal*renormal;

		cpuType rho = 1.0/((float)NSIDEMESH );
		rho = rho * rho * rho;
		rho *= (float)NPART_TOTAL ;
		
		c=0;
		for (i=0; i<isize[0]; i++)
			for (j=0; j<isize[1]; j++)
				for (k=0; k<isize[2]; k++) {

					pmesh[c] -= rho;
					pmesh[c] /= rho;
					c++;
				}



		///////////////////////////
		MPI_Isend(pmesh, meshsize, MPI_CPUTYPE, PROC_RANK_SUDOM, 100, MPI_COMM_WORLD, &req);

	}


	if ( isSUDOM ) {
		int local_size_y =  NSIDEMESH/NSIDE0SUDOM; 
		int local_size_z =  NSIDEMESH/NSIDE1SUDOM;
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			nrecv[n] = local_size_y  * local_size_z * dom_grp_mesh_size[n] ;	
		}
		ip = 0;
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			//		printf(" nrecv = %d\n", nrecv[n]);
			MPI_Recv( smesh + ip , nrecv[n], MPI_CPUTYPE, dom_grp[n], 100, MPI_COMM_WORLD, &(status[n]));
			ip += nrecv[n];
		}
	}

	if ( COM_DOM) {
		MPI_Wait(&req, &sta);		
		t4 = dtime();
		//	printf(" [%d] dtime = %lf ,mesh build %lf %lf , com=%lf tot=%lf \n", PROC_RANK, t1 -t0, t2-t1, t3-t2, t4-t3, t4-t0);
	}

	if ( COM_DOM ) {	
		myfree(pmesh, 261);
	}

	//////////////////

	if (  isSUDOM ) {
		mbuff = (cpuType*)mymalloc(sizeof(cpuType)*meshsize_sud, 262);
		int m;

		if (mbuff == NULL) {
			printf(" error mbuff \n");
			exit(0);
		}

		for (n=0; n<meshsize_sud; n++) 
			mbuff[n] = 0.0;

		int i,j,k;
		int ndisp[NDOMinSUDOM];
		ndisp[0] = 0;
		for (m=1; m<NDOMinSUDOM; m++) {
			ndisp[m] = ndisp[m-1] + nrecv[m-1];
		}
		int nzz = NSIDE1PROC/NSIDE1SUDOM;
		int y_n = NSIDEMESH/NSIDE0SUDOM;
		int z_n = NSIDEMESH/NSIDE1SUDOM;

		int dy = NSIDEMESH/NSIDE0PROC;
		int dz = NSIDEMESH/NSIDE1PROC;
		int nm = NSIDEMESH*dy*dz;

		cpuType *sm;
		for (m=0; m<NDOMinSUDOM; m++) {
			int x_n = dom_grp_mesh_size[m];
			// try		//	sm = smesh + m * x_n * y_n * z_n;
			sm = smesh + dom_grp_mesh_start[m] * y_n * z_n;
			int ii, jj, kk;
			int in, jn, kn;
			int ig, idx;
			for (i=0; i<x_n; i++) {
				for (j=0; j<y_n; j++) {
					for (k=0; k<z_n; k++) {
						jj = j/dy;
						kk = k/dz;
						ig = jj*nzz + kk;

						in = dom_grp_mesh_start[m] + i;
						jn = j%dy;
						kn = k%dz;
						idx = ( in * dy + jn ) * dz + kn;

						mbuff[ig*nm + idx] = sm[(i*y_n+j)*z_n+k];
					}
				}
			}
		}

		myfree(smesh, 260);
		for (m=0; m<NDOMinSUDOM; m++) {	
			MPI_Isend( mbuff + m*nm, nm, MPI_CPUTYPE, dom_grp[m], 200, MPI_COMM_WORLD, &(request[m]));
		}
	}


	int x_n = NSIDEMESH;
	int y_n = NSIDEMESH/NSIDE0PROC;
	int z_n = NSIDEMESH/NSIDE1PROC;
	meshsize_proc =  x_n  * y_n * z_n;
	mesh = (cpuType*)mymalloc(sizeof(cpuType) * meshsize_proc, 263);

		if (mesh == NULL) {
			printf(" error mesh in power_spect \n");
			exit(0);
		}



	MPI_Recv(mesh, 	meshsize_proc, MPI_CPUTYPE, PROC_RANK_SUDOM, 200, MPI_COMM_WORLD, &sta); 
	if (  isSUDOM ) {
		int m;	
		for (m=0; m<NDOMinSUDOM; m++) {	
			MPI_Wait(&(request[m]), &(status[m]) );
		}
		myfree(mbuff, 262);
	}
	/// canbe free mbuff and smesh here once
	///////////////////////
	int nside[3] = {NSIDEMESH, NSIDEMESH, NSIDEMESH};
	cpuType param[2] = { splitRadius, BOXSIZE};

	FILE *fd;


	{
		int *psn = (int*)mymalloc(sizeof(int)*NSIDEMESH, 264);
		cpuType *ps = (cpuType*)mymalloc(sizeof(cpuType)*NSIDEMESH, 265);
		cpuType *ps_total;	
		for (n=0; n<NSIDEMESH; n++) {
			psn[n] = 0;
			ps[n] = 0.0;
		}
		powsk(mesh,nside,param, ps, psn);

		float *psn_total;	
		float *psn_rank = (float*)mymalloc(sizeof(float)*NSIDEMESH, 266);	
		for (n=0; n<NSIDEMESH; n++) {
			psn_rank[n] = (float)psn[n];

		}
		if (0 == PROC_RANK) {

			ps_total = (cpuType*)mymalloc(sizeof(cpuType)*NSIDEMESH, 267);
			psn_total = (float*)mymalloc(sizeof(float)*NSIDEMESH, 268);
		}

		MPI_Reduce(ps, ps_total, NSIDEMESH, MPI_CPUTYPE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(psn_rank, psn_total, NSIDEMESH, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

		if (0 == PROC_RANK ) {
			for (n=0; n<NSIDEMESH; n++) {
				if (psn_total[n] > 0)
					ps_total[n] /= (psn_total[n]);	
				else 
					ps_total[n] = 0.0;
			}

			fd = fopen(fname, "w");
			if (fd == NULL) {
				printf( " Error ! cannot open powspec file \n"  );
				exit(0);
			}
			float vol = BOXSIZE/1000.0;
			vol = vol * vol * vol;
			float Norm = vol/((float)NSIDEMESH*NSIDEMESH *NSIDEMESH);

			int output_num = (int)(NSIDEMESH*0.75);			
//			output_num = NSIDEMESH/2;			

			for (n=1; n<output_num; n++) {
				fprintf(fd, "%e %e %lf\n",  (n+0.5)*2.0*M_PI/(BOXSIZE/1000.0), ps_total[n] * Norm *Norm /vol, psn_total[n]  );
			}

			fclose(fd);
		}



		if ( 0 == PROC_RANK){
			myfree(ps_total, 267);
			myfree(psn_total, 268);
		}
		myfree(psn, 264);
		myfree(psn_rank, 266);
		myfree(ps, 265);
	}
///	if ( COM_DOM ) {	
//
//		myfree(pmesh, 261);
//	}

	if (  isSUDOM ) {
//		myfree(mbuff, 262);
//		myfree(smesh, 260);
	}

	//	MPI_Barrier(MPI_COMM_WORLD);
	myfree(mesh, 263);
}


void power_spectrum(int rank_snap) {

#ifdef POW_SPEC
	power_spec(rank_snap);
#endif

#ifdef POW_SPEC_DAUB
	power_spec_daub(rank_snap);
#endif

}
#endif

#else

void power_spec(int rank_snap){

	char fname[128];
	sprintf(fname, "%s/powspec_%04d", PathSnapshot, rank_snap);
	double t0, t1, t2, t3, t4, t5;
	cpuType *smesh;
	cpuType *mbuff;
	MPI_Request req;
	MPI_Status  sta;
	MPI_Request request[NDOMinSUDOM];
	MPI_Status  status[NDOMinSUDOM];
	int n, ip;
	int meshsize_sud;
	int meshsize_proc;
	int meshsize ;
	int meshsize_pot ;

	int local_min[3], local_max[3], isize[3];
	int nrecv[NDOMinSUDOM];	

	if (  isSUDOM ) {
		int nx, ny, nz;
		nx = NSIDEMESH;
		ny = NSIDEMESH/NSIDE0SUDOM;
		nz = NSIDEMESH/NSIDE1SUDOM;
		meshsize_sud = nx * ny * nz;
		//		meshsize_proc = nx * ny * nz;
		smesh = (cpuType*)mymalloc(sizeof(cpuType)* meshsize_sud, 260) ;

		if (smesh == NULL) {
			printf(" error smesh \n");
			exit(0);
		}
		for (n=0; n<meshsize_sud; n++) 
			smesh[n] = 0.0;

	}

	if ( COM_DOM ) {	
		t0 = dtime();

		int  c ;

		local_min[0] = dom_grp_mesh_start[dom_grp_rank];
		local_max[0] = local_min[0] + dom_grp_mesh_size[dom_grp_rank];

		int idx = PROC_RANK;

		int dside0 = NSIDE0PROC/NSIDE0SUDOM;
		int dside1 = NSIDE1PROC/NSIDE1SUDOM;
		int y = idx / NSIDE1PROC ; 
		int z = idx % NSIDE1PROC ; 

		int sy = y / dside0;
		int sz = z / dside1;

		local_min[1] = sy * NSIDEMESH/NSIDE0SUDOM; 
		local_min[2] = sz * NSIDEMESH/NSIDE1SUDOM; 

		local_max[1] =(sy+1) * NSIDEMESH/NSIDE0SUDOM; 
		local_max[2] =(sz+1) * NSIDEMESH/NSIDE1SUDOM; 

		isize[0] = local_max[0] - local_min[0] ;			
		isize[1] = NSIDEMESH / NSIDE0SUDOM;
		isize[2] = NSIDEMESH / NSIDE1SUDOM;

		//		printf(" [%d]local_min (%d %d) <%d %d %d %d>\n", PROC_RANK, local_min[0], local_max[0],dom_grp[0], dom_grp[1], dom_grp[2], dom_grp[3]);

		meshsize = isize[0] * isize[1] * isize[2];
		meshsize_pot = (isize[0] + 2*NUMPAD)  * (isize[1]+2*NUMPAD) * (isize[2]+2*NUMPAD);
		pmesh = (cpuType*)mymalloc(sizeof(cpuType) * meshsize, 261);

		for (n=0; n<meshsize; n++) 
			pmesh[n] = 0.0;

		int i,j,k;
		int ii, jj, kk;
		cpuType norm = NSIDEMESH/BOXSIZE;
		cpuType delta = 1.0/norm;

		t1 = dtime();

		cpuType wi, wj , wk, win, wjn, wkn;
		c=0;	
		for (n=0; n<NPART; n++) {
			if (part[n].tag != 1) continue;
	
			cpuType posf[3];
			posf[0] = (float) (part[n].posi[0] * INT2POS );
			posf[1] = (float) (part[n].posi[1] * INT2POS );
			posf[2] = (float) (part[n].posi[2] * INT2POS );

			i = (int) (posf[0] *norm) ;	
			j = (int) (posf[1] *norm) ;	
			k = (int) (posf[2] *norm) ;	

			wi = (posf[0] - (i+0.5)*delta)*norm;

			if (wi > 0) {
				ii = i + 1;
			}
			else {
				wi = -wi;
				ii = i - 1;
			}
			win =  1.0 - wi;

			wj = (posf[1] - (j+0.5)*delta)*norm;
			if (wj > 0) {
				jj = j + 1;
			}
			else {
				wj = -wj;
				jj = j - 1;
			}
			wjn =  1.0 - wj;

			wk = (posf[2] - (k+0.5)*delta)*norm;
			if (wk > 0) {
				kk = k + 1;
			}
			else {
				wk = -wk;
				kk = k - 1;
			}
			wkn =  1.0 - wk;

			i -= local_min[0] ;	
			j -= local_min[1] ;	
			k -= local_min[2] ;	

			ii -= local_min[0] ;	
			jj -= local_min[1] ;	
			kk -= local_min[2] ;	

			if (k >=0 && k< isize[2]) { 
				if (j>=0 && j<isize[1])	{
					if (i>=0 && i<isize[0]) {	
						idx = (i*isize[1] + j)*isize[2] + k;
						pmesh[idx] += MASSPART*win*wjn*wkn;	
					}

					if (ii>=0 && ii<isize[0]) {	
						idx = (ii*isize[1] + j)*isize[2] + k;
						pmesh[idx] += MASSPART*wi *wjn*wkn;
					}	
				}


				if (jj>=0 && jj<isize[1])	{
					if (i>=0 && i<isize[0]) {	
						idx = (i*isize[1] + jj)*isize[2] + k;
						pmesh[idx] += MASSPART*win*wj *wkn;	
					}
					if (ii>=0 && ii<isize[0]) {	
						idx = (ii*isize[1] + jj)*isize[2] + k;
						pmesh[idx] += MASSPART*wi *wj *wkn;	
					}
				}
			}

			if (kk>=0 && kk< isize[2]) {
				if (j>=0 && j<isize[1])	{
					if (i>=0 && i<isize[0]) {	
						idx = (i*isize[1] + j)*isize[2] + kk;
						pmesh[idx] += MASSPART*win*wjn*wk ;	
					}

					if (ii>=0 && ii<isize[0]) {	
						idx = (ii*isize[1] + j)*isize[2] + kk;
						pmesh[idx] += MASSPART*wi *wjn*wk ;	
					}
				}


				if (jj>=0 && jj<isize[1]){
					if (i>=0 && i<isize[0]) {	
						idx = (i*isize[1] + jj)*isize[2] + kk;
						pmesh[idx] += MASSPART*win*wj *wk ;	
					}
					if (ii>=0 && ii<isize[0]) {	
						idx = (ii*isize[1] + jj)*isize[2] + kk;
						pmesh[idx] += MASSPART* wi *wj *wk ;
					}
				}
			}
			c++;
		}

		t2 = dtime();

		c=0;	
		for (n=0; n<NPART_BND; n++) {

			cpuType posf[3];
			posf[0] = (float) (part_b[n].posi[0] * INT2POS );
			posf[1] = (float) (part_b[n].posi[1] * INT2POS );
			posf[2] = (float) (part_b[n].posi[2] * INT2POS );



			i = (int) (posf[0] *norm) ;	
			j = (int) (posf[1] *norm) ;	
			k = (int) (posf[2] *norm) ;	





			if (posf[0] < 0.0 )
				i--;

			if (posf[1] < 0.0 )
				j--;

			if (posf[2] < 0.0 )
				k--;

			wi = (posf[0] - (i+0.5)*delta)*norm;

			if (wi > 0) {
				ii = i + 1;
			}
			else {
				wi = -wi;
				ii = i - 1;
			}
			win =  1.0 - wi;

			wj = ( posf[1] - (j+0.5)*delta)*norm;
			if (wj > 0) {
				jj = j + 1;
			}
			else {
				wj = -wj;
				jj = j - 1;
			}
			wjn =  1.0 - wj;

			wk = ( posf[2] - (k+0.5)*delta)*norm;
			if (wk > 0) {
				kk = k + 1;
			}
			else {
				wk = -wk;
				kk = k - 1;
			}
			wkn =  1.0 - wk;

			i -= local_min[0] ;	
			j -= local_min[1] ;	
			k -= local_min[2] ;	

			ii -= local_min[0] ;	
			jj -= local_min[1] ;	
			kk -= local_min[2] ;	

			//if	
			if (k >=0 && k< isize[2]) { 
				if (j>=0 && j<isize[1])	{
					if (i>=0 && i<isize[0]) {	
						idx = (i*isize[1] + j)*isize[2] + k;
						pmesh[idx] += MASSPART*win*wjn*wkn;	
					}

					if (ii>=0 && ii<isize[0]) {	
						idx = (ii*isize[1] + j)*isize[2] + k;
						pmesh[idx] += MASSPART*wi *wjn*wkn;
					}	
				}


				if (jj>=0 && jj<isize[1])	{
					if (i>=0 && i<isize[0]) {	
						idx = (i*isize[1] + jj)*isize[2] + k;
						pmesh[idx] += MASSPART*win*wj *wkn;	
					}
					if (ii>=0 && ii<isize[0]) {	
						idx = (ii*isize[1] + jj)*isize[2] + k;
						pmesh[idx] += MASSPART*wi *wj *wkn;	
					}
				}
			}

			if (kk>=0 && kk< isize[2]) {
				if (j>=0 && j<isize[1])	{
					if (i>=0 && i<isize[0]) {	
						idx = (i*isize[1] + j)*isize[2] + kk;
						pmesh[idx] += MASSPART*win*wjn*wk ;	
					}

					if (ii>=0 && ii<isize[0]) {	
						idx = (ii*isize[1] + j)*isize[2] + kk;
						pmesh[idx] += MASSPART*wi *wjn*wk ;	
					}
				}


				if (jj>=0 && jj<isize[1])	{
					if (i>=0 && i<isize[0]) {	
						idx = (i*isize[1] + jj)*isize[2] + kk;
						pmesh[idx] += MASSPART*win*wj *wk ;	
					}
					if (ii>=0 && ii<isize[0]) {	
						idx = (ii*isize[1] + jj)*isize[2] + kk;
						pmesh[idx] += MASSPART* wi *wj *wk ;
					}
				}
			}
			c++;
		}

		t3 = dtime();

		cpuType renormal = (NSIDEMESH/BOXSIZE);
		renormal = renormal*renormal*renormal;

		renormal = 1.0/MASSPART;

		c=0;
		for (i=0; i<isize[0]; i++)
			for (j=0; j<isize[1]; j++)
				for (k=0; k<isize[2]; k++) {

					pmesh[c] *= renormal;
					c++;
				}
		MPI_Isend(pmesh, meshsize, MPI_CPUTYPE, PROC_RANK_SUDOM, 100, MPI_COMM_WORLD, &req);

	}


	if ( isSUDOM ) {
		int local_size_y =  NSIDEMESH/NSIDE0SUDOM; 
		int local_size_z =  NSIDEMESH/NSIDE1SUDOM;
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			nrecv[n] = local_size_y  * local_size_z * dom_grp_mesh_size[n] ;	
		}
		ip = 0;
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			//		printf(" nrecv = %d\n", nrecv[n]);
			MPI_Recv( smesh + ip , nrecv[n], MPI_CPUTYPE, dom_grp[n], 100, MPI_COMM_WORLD, &(status[n]));
			ip += nrecv[n];
		}
	}

	if ( COM_DOM) {
		MPI_Wait(&req, &sta);		
		t4 = dtime();
		//	printf(" [%d] dtime = %lf ,mesh build %lf %lf , com=%lf tot=%lf \n", PROC_RANK, t1 -t0, t2-t1, t3-t2, t4-t3, t4-t0);
	}
	if ( COM_DOM ) {	

		myfree(pmesh, 261);
	}


	//////////////////

	if (  isSUDOM ) {
		mbuff = (cpuType*)mymalloc(sizeof(cpuType)*meshsize_sud, 262);
		int m;

		for (n=0; n<meshsize_sud; n++) 
			mbuff[n] = 0.0;

		int i,j,k;
		int ndisp[NDOMinSUDOM];
		ndisp[0] = 0;
		for (m=1; m<NDOMinSUDOM; m++) {
			ndisp[m] = ndisp[m-1] + nrecv[m-1];
		}
		int nzz = NSIDE1PROC/NSIDE1SUDOM;
		int y_n = NSIDEMESH/NSIDE0SUDOM;
		int z_n = NSIDEMESH/NSIDE1SUDOM;

		int dy = NSIDEMESH/NSIDE0PROC;
		int dz = NSIDEMESH/NSIDE1PROC;
		int nm = NSIDEMESH*dy*dz;

		cpuType *sm;
		for (m=0; m<NDOMinSUDOM; m++) {
			int x_n = dom_grp_mesh_size[m];
			// try		//	sm = smesh + m * x_n * y_n * z_n;
			sm = smesh + dom_grp_mesh_start[m] * y_n * z_n;
			int ii, jj, kk;
			int in, jn, kn;
			int ig, idx;
			for (i=0; i<x_n; i++) {
				for (j=0; j<y_n; j++) {
					for (k=0; k<z_n; k++) {
						jj = j/dy;
						kk = k/dz;
						ig = jj*nzz + kk;

						in = dom_grp_mesh_start[m] + i;
						jn = j%dy;
						kn = k%dz;
						idx = ( in * dy + jn ) * dz + kn;

						mbuff[ig*nm + idx] = sm[(i*y_n+j)*z_n+k];
					}
				}
			}
		}

		myfree(smesh, 260);
		for (m=0; m<NDOMinSUDOM; m++) {	
			MPI_Isend( mbuff + m*nm, nm, MPI_CPUTYPE, dom_grp[m], 200, MPI_COMM_WORLD, &(request[m]));
		}
	}


	int x_n = NSIDEMESH;
	int y_n = NSIDEMESH/NSIDE0PROC;
	int z_n = NSIDEMESH/NSIDE1PROC;
	meshsize_proc =  x_n  * y_n * z_n;
	mesh = (cpuType*)mymalloc(sizeof(cpuType) * meshsize_proc, 263);

	MPI_Recv(mesh, 	meshsize_proc, MPI_CPUTYPE, PROC_RANK_SUDOM, 200, MPI_COMM_WORLD, &sta); 
	if (  isSUDOM ) {
		int m;	
		for (m=0; m<NDOMinSUDOM; m++) {	
			MPI_Wait(&(request[m]), &(status[m]) );
		}
		myfree(mbuff, 262);
	}
	/// canbe free mbuff and smesh here once
	///////////////////////
	int nside[3] = {NSIDEMESH, NSIDEMESH, NSIDEMESH};
	cpuType param[2] = { splitRadius, BOXSIZE};

	FILE *fd;


	{
		int *psn = (int*)mymalloc(sizeof(int)*NSIDEMESH, 264);
		cpuType *ps = (cpuType*)mymalloc(sizeof(cpuType)*NSIDEMESH, 265);
		cpuType *ps_total;	
		for (n=0; n<NSIDEMESH; n++) {
			psn[n] = 0;
			ps[n] = 0.0;
		}
		powsk(mesh,nside,param, ps, psn);

		float *psn_total;	
		float *psn_rank = (float*)mymalloc(sizeof(float)*NSIDEMESH, 266);	
		for (n=0; n<NSIDEMESH; n++) {
			psn_rank[n] = (float)psn[n];

		}
		if (0 == PROC_RANK) {

			ps_total = (cpuType*)mymalloc(sizeof(cpuType)*NSIDEMESH, 267);
			psn_total = (float*)mymalloc(sizeof(float)*NSIDEMESH, 268);
		}

		MPI_Reduce(ps, ps_total, NSIDEMESH, MPI_CPUTYPE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(psn_rank, psn_total, NSIDEMESH, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

		if (0 == PROC_RANK ) {
			for (n=0; n<NSIDEMESH; n++) {
				if (psn_total[n] > 0)
					ps_total[n] /= (psn_total[n]);	
				else 
					ps_total[n] = 0.0;
			}

			fd = fopen(fname, "w");
			if (fd == NULL) {
				printf( " Error ! cannot open powspec file \n"  );
				exit(0);
			}
			float vol = BOXSIZE/1000.0;
			vol = vol * vol * vol;
			float Norm = vol/((float)NSIDEMESH*NSIDEMESH *NSIDEMESH);

			int output_num = (int)(NSIDEMESH*0.75);			
			output_num = NSIDEMESH/2;			
//			float noise = sqrt(Norm);
			for (n=1; n<output_num; n++) {
			//	fprintf(fd, "%e %e %e %f\n",  (n+0.5)*2.0*M_PI/(BOXSIZE/1000.0) , ps_total[n] * Norm *Norm /vol,   ps_total[n] * Norm *Norm /vol - noise , psn_total[n]  );
				fprintf(fd, "%e %e %f\n",  (n+0.5)*2.0*M_PI/(BOXSIZE/1000.0), ps_total[n] * Norm *Norm /vol, psn_total[n]  );
			}

			fclose(fd);
		}



		if ( 0 == PROC_RANK){
			myfree(ps_total, 267);
			myfree(psn_total, 268);
		}
		myfree(psn, 264);
		myfree(psn_rank, 266);
		myfree(ps, 265);
	}
	myfree(mesh, 263);
}



///////////////////////////// daubechies12 ////////////////////////////


#include "daub12hash.h"
void initial_daub12(void)
{
    int n, N;
    float x;
    FILE *fp;
    if (!(fp=fopen("../demo/daub12hash.dat","r")))
    {
        printf("can not open `daub12.dat'\n");
        exit(0);
    }
    N = 600001;
    DAUB12 = (float*)malloc((int)N*sizeof(float));
    for (n=0; n<N; n++)
        fscanf(fp, "%f %f", &x, &(DAUB12[n]) );
    
//    for (n=100000; n<100300; n++)
//        printf("%f  ", DAUB12[n]);
//    printf("\n");
    fclose(fp);
}

float daub12_scaling_function(float d)
{
//    if (d<0.0 || d>=6.0)
//    {
//        printf("error in daub12 %f \n", d);
//        exit(10);
//    }   
	if (0.0<d && d<6.0) 
		return DAUB12[(int)(d*100000.0)];
	else 
		return 0.0;
}


void free_daub12(void)
{
    free(DAUB12);
}

void power_spec_daub(int rank_snap)
{
	char fname[128];
	sprintf(fname, "%s/powspec_daub_%d", PathSnapshot, rank_snap);
	double t0, t1, t2, t3, t4, t5;
	cpuType *smesh;
	cpuType *mbuff;
	MPI_Request req;
	MPI_Status  sta;
	MPI_Request request[NDOMinSUDOM];
	MPI_Status  status[NDOMinSUDOM];
	int n, ip;
	int meshsize_sud;
	int meshsize_proc;
	int meshsize ;
	int meshsize_pot ;

	int local_min[3], local_max[3], isize[3];
	int nrecv[NDOMinSUDOM];	

	if (  isSUDOM ) {
		int nx, ny, nz;
		nx = NSIDEMESH;
		ny = NSIDEMESH/NSIDE0SUDOM;
		nz = NSIDEMESH/NSIDE1SUDOM;
		meshsize_sud = nx * ny * nz;
		//		meshsize_proc = nx * ny * nz;
		smesh = (cpuType*)mymalloc(sizeof(cpuType)* meshsize_sud, 260) ;

		if (smesh == NULL) {
			printf(" error smesh \n");
			exit(0);
		}

		for (n=0; n<meshsize_sud; n++) 
			smesh[n] = 0.0;
	}

	if ( COM_DOM ) {	
		t0 = dtime();

		int  c ;

		local_min[0] = dom_grp_mesh_start[dom_grp_rank];
		local_max[0] = local_min[0] + dom_grp_mesh_size[dom_grp_rank];

		int idx = PROC_RANK;

		int dside0 = NSIDE0PROC/NSIDE0SUDOM;
		int dside1 = NSIDE1PROC/NSIDE1SUDOM;
		int y = idx / NSIDE1PROC ; 
		int z = idx % NSIDE1PROC ; 

		int sy = y / dside0;
		int sz = z / dside1;

		local_min[1] = sy * NSIDEMESH/NSIDE0SUDOM; 
		local_min[2] = sz * NSIDEMESH/NSIDE1SUDOM; 

		local_max[1] =(sy+1) * NSIDEMESH/NSIDE0SUDOM; 
		local_max[2] =(sz+1) * NSIDEMESH/NSIDE1SUDOM; 

		isize[0] = local_max[0] - local_min[0] ;			
		isize[1] = NSIDEMESH / NSIDE0SUDOM;
		isize[2] = NSIDEMESH / NSIDE1SUDOM;

		//		printf(" [%d]local_min (%d %d) <%d %d %d %d>\n", PROC_RANK, local_min[0], local_max[0],dom_grp[0], dom_grp[1], dom_grp[2], dom_grp[3]);

		meshsize = isize[0] * isize[1] * isize[2];
		meshsize_pot = (isize[0] + 2*NUMPAD)  * (isize[1]+2*NUMPAD) * (isize[2]+2*NUMPAD);
		pmesh = (cpuType*)mymalloc(sizeof(cpuType) * meshsize, 261);

		for (n=0; n<meshsize; n++) 
			pmesh[n] = 0.0;

		int i,j,k;
		int ii, jj, kk;
		cpuType norm = NSIDEMESH/BOXSIZE;
		cpuType delta = 1.0/norm;



		t1 = dtime();

		cpuType wx[6], wy[6] , wz[6], wt;
		float gx, gy, gz;
		float px, py, pz;
		int Ic, Jc, Kc, Is, Js, Ks;
		int i_s, j_s, k_s;
	
		initial_daub12();

		int nside = NSIDEMESH;
		cpuType Dlt = 1.0/((cpuType) NSIDEMESH );
		c=0;	
		for (n=0; n<NPART; n++) {
			if (part[n].tag != 1) continue;

			cpuType posf[3];
			posf[0] = (float) (part[n].posi[0] * INT2POS );
			posf[1] = (float) (part[n].posi[1] * INT2POS );
			posf[2] = (float) (part[n].posi[2] * INT2POS );

			px = posf[0] / BOXSIZE;
			py = posf[1] / BOXSIZE;
			pz = posf[2] / BOXSIZE;
			Ic = (int) (px * nside) ;	
			Jc = (int) (py * nside) ;	
			Kc = (int) (pz * nside) ;	
			gx = (float)(Ic+0.5)*Dlt;
			gy = (float)(Jc+0.5)*Dlt;
			gz = (float)(Kc+0.5)*Dlt;

			if (px>gx) Is = -2; else Is = -3;
			if (py>gy) Js = -2; else Js = -3;
			if (pz>gz) Ks = -2; else Ks = -3;

			for (i=0; i<6; i++)
				wx[i] = daub12_scaling_function( (float)( (gx-px)/Dlt + (i+Is+3) ));

			for (j=0; j<6; j++)
				wy[j] = daub12_scaling_function( (float)( (gy-py)/Dlt + (j+Js+3)) );

			for (k=0; k<6; k++)
				wz[k] = daub12_scaling_function( (float)( (gz-pz)/Dlt + (k+Ks+3)) );

			for (i=0; i<6; i++)
				for (j=0; j<6; j++)
					for (k=0; k<6; k++)  {
						i_s =(Ic+Is+i+nside)%nside;
						j_s =(Jc+Js+j+nside)%nside;
						k_s =(Kc+Ks+k+nside)%nside;          
						i_s -= local_min[0] ;	
						j_s -= local_min[1] ;	
						k_s -= local_min[2] ;	
						if (k_s >=0 && k_s< isize[2]) { 
							if (j_s>=0 && j_s<isize[1])	{
								if (i_s>=0 && i_s<isize[0]) {	

									//	grid[(i_s*nside + j_s) * nside + k_s] += (wx[i]*wy[j]*wz[k]);
									idx = (i_s*isize[1] + j_s)*isize[2] + k_s;
									pmesh[idx] += (wx[i]*wy[j]*wz[k]);
								}
							}
						}
					}
		}
		t2 = dtime();

		c=0;	
		for (n=0; n<NPART_BND; n++) {


			cpuType posf[3];
			posf[0] = (float) (part_b[n].posi[0] * INT2POS );
			posf[1] = (float) (part_b[n].posi[1] * INT2POS );
			posf[2] = (float) (part_b[n].posi[2] * INT2POS );


			px = posf[0] / BOXSIZE;
			py = posf[1] / BOXSIZE;
			pz = posf[2] / BOXSIZE;
			Ic = (int) (px * NSIDEMESH) ;	
			Jc = (int) (py * NSIDEMESH) ;	
			Kc = (int) (pz * NSIDEMESH) ;	
		//	Ic = (int) (part[n].pos[0] *norm) ;	
		//	Jc = (int) (part[n].pos[1] *norm) ;	
		//	Kc = (int) (part[n].pos[2] *norm) ;	
			if (posf[0] < 0.0 )
				Ic--;

			if (posf[1] < 0.0 )
				Jc--;

			if (posf[2] < 0.0 )
				Kc--;

			gx = (float)(Ic+0.5)*Dlt;
			gy = (float)(Jc+0.5)*Dlt;
			gz = (float)(Kc+0.5)*Dlt;

			if (px>gx) Is = -2; else Is = -3;
			if (py>gy) Js = -2; else Js = -3;
			if (pz>gz) Ks = -2; else Ks = -3;

			for (i=0; i<6; i++)
				wx[i] = daub12_scaling_function( (float)( (gx-px)/Dlt + (i+Is+3) ));

			for (j=0; j<6; j++)
				wy[j] = daub12_scaling_function( (float)( (gy-py)/Dlt + (j+Js+3)) );

			for (k=0; k<6; k++)
				wz[k] = daub12_scaling_function( (float)( (gz-pz)/Dlt + (k+Ks+3)) );

			for (i=0; i<6; i++)
				for (j=0; j<6; j++)
					for (k=0; k<6; k++)  {
						i_s =(Ic+Is+i+nside)%nside;
						j_s =(Jc+Js+j+nside)%nside;
						k_s =(Kc+Ks+k+nside)%nside;          
						i_s -= local_min[0] ;	
						j_s -= local_min[1] ;	
						k_s -= local_min[2] ;	
						if (k_s >=0 && k_s< isize[2]) { 
							if (j_s>=0 && j_s<isize[1])	{
								if (i_s>=0 && i_s<isize[0]) {	

									//	grid[(i_s*nside + j_s) * nside + k_s] += (wx[i]*wy[j]*wz[k]);
									idx = (i_s*isize[1] + j_s)*isize[2] + k_s;
									pmesh[idx] += (wx[i]*wy[j]*wz[k]);
								}
							}
						}
					}

		}

//printf(" power 2 - daub12 \n");
	free_daub12();
//	rho = (double)npart/((double)nside*nside*nside);
		t3 = dtime();

		cpuType renormal = (NSIDEMESH/BOXSIZE);
		renormal = renormal*renormal*renormal;

		cpuType rho = 1.0/((float)NSIDEMESH );
		rho = rho * rho * rho;
		rho *= (float)NPART_TOTAL ;
		
		c=0;
		for (i=0; i<isize[0]; i++)
			for (j=0; j<isize[1]; j++)
				for (k=0; k<isize[2]; k++) {

					pmesh[c] -= rho;
					pmesh[c] /= rho;
					c++;
				}



		///////////////////////////
		MPI_Isend(pmesh, meshsize, MPI_CPUTYPE, PROC_RANK_SUDOM, 100, MPI_COMM_WORLD, &req);

	}


	if ( isSUDOM ) {
		int local_size_y =  NSIDEMESH/NSIDE0SUDOM; 
		int local_size_z =  NSIDEMESH/NSIDE1SUDOM;
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			nrecv[n] = local_size_y  * local_size_z * dom_grp_mesh_size[n] ;	
		}
		ip = 0;
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			//		printf(" nrecv = %d\n", nrecv[n]);
			MPI_Recv( smesh + ip , nrecv[n], MPI_CPUTYPE, dom_grp[n], 100, MPI_COMM_WORLD, &(status[n]));
			ip += nrecv[n];
		}
	}

	if ( COM_DOM) {
		MPI_Wait(&req, &sta);		
		t4 = dtime();
		//	printf(" [%d] dtime = %lf ,mesh build %lf %lf , com=%lf tot=%lf \n", PROC_RANK, t1 -t0, t2-t1, t3-t2, t4-t3, t4-t0);
	}

	if ( COM_DOM ) {	
		myfree(pmesh, 261);
	}

	//////////////////

	if (  isSUDOM ) {
		mbuff = (cpuType*)mymalloc(sizeof(cpuType)*meshsize_sud, 262);
		int m;

		if (mbuff == NULL) {
			printf(" error mbuff \n");
			exit(0);
		}

		for (n=0; n<meshsize_sud; n++) 
			mbuff[n] = 0.0;

		int i,j,k;
		int ndisp[NDOMinSUDOM];
		ndisp[0] = 0;
		for (m=1; m<NDOMinSUDOM; m++) {
			ndisp[m] = ndisp[m-1] + nrecv[m-1];
		}
		int nzz = NSIDE1PROC/NSIDE1SUDOM;
		int y_n = NSIDEMESH/NSIDE0SUDOM;
		int z_n = NSIDEMESH/NSIDE1SUDOM;

		int dy = NSIDEMESH/NSIDE0PROC;
		int dz = NSIDEMESH/NSIDE1PROC;
		int nm = NSIDEMESH*dy*dz;

		cpuType *sm;
		for (m=0; m<NDOMinSUDOM; m++) {
			int x_n = dom_grp_mesh_size[m];
			// try		//	sm = smesh + m * x_n * y_n * z_n;
			sm = smesh + dom_grp_mesh_start[m] * y_n * z_n;
			int ii, jj, kk;
			int in, jn, kn;
			int ig, idx;
			for (i=0; i<x_n; i++) {
				for (j=0; j<y_n; j++) {
					for (k=0; k<z_n; k++) {
						jj = j/dy;
						kk = k/dz;
						ig = jj*nzz + kk;

						in = dom_grp_mesh_start[m] + i;
						jn = j%dy;
						kn = k%dz;
						idx = ( in * dy + jn ) * dz + kn;

						mbuff[ig*nm + idx] = sm[(i*y_n+j)*z_n+k];
					}
				}
			}
		}

		myfree(smesh, 260);
		for (m=0; m<NDOMinSUDOM; m++) {	
			MPI_Isend( mbuff + m*nm, nm, MPI_CPUTYPE, dom_grp[m], 200, MPI_COMM_WORLD, &(request[m]));
		}
	}


	int x_n = NSIDEMESH;
	int y_n = NSIDEMESH/NSIDE0PROC;
	int z_n = NSIDEMESH/NSIDE1PROC;
	meshsize_proc =  x_n  * y_n * z_n;
	mesh = (cpuType*)mymalloc(sizeof(cpuType) * meshsize_proc, 263);

		if (mesh == NULL) {
			printf(" error mesh in power_spect \n");
			exit(0);
		}



	MPI_Recv(mesh, 	meshsize_proc, MPI_CPUTYPE, PROC_RANK_SUDOM, 200, MPI_COMM_WORLD, &sta); 
	if (  isSUDOM ) {
		int m;	
		for (m=0; m<NDOMinSUDOM; m++) {	
			MPI_Wait(&(request[m]), &(status[m]) );
		}
		myfree(mbuff, 262);
	}
	/// canbe free mbuff and smesh here once
	///////////////////////
	int nside[3] = {NSIDEMESH, NSIDEMESH, NSIDEMESH};
	cpuType param[2] = { splitRadius, BOXSIZE};

	FILE *fd;


	{
		int *psn = (int*)mymalloc(sizeof(int)*NSIDEMESH, 264);
		cpuType *ps = (cpuType*)mymalloc(sizeof(cpuType)*NSIDEMESH, 265);
		cpuType *ps_total;	
		for (n=0; n<NSIDEMESH; n++) {
			psn[n] = 0;
			ps[n] = 0.0;
		}
		powsk(mesh,nside,param, ps, psn);

		float *psn_total;	
		float *psn_rank = (float*)mymalloc(sizeof(float)*NSIDEMESH, 266);	
		for (n=0; n<NSIDEMESH; n++) {
			psn_rank[n] = (float)psn[n];

		}
		if (0 == PROC_RANK) {

			ps_total = (cpuType*)mymalloc(sizeof(cpuType)*NSIDEMESH, 267);
			psn_total = (float*)mymalloc(sizeof(float)*NSIDEMESH, 268);
		}

		MPI_Reduce(ps, ps_total, NSIDEMESH, MPI_CPUTYPE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(psn_rank, psn_total, NSIDEMESH, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

		if (0 == PROC_RANK ) {
			for (n=0; n<NSIDEMESH; n++) {
				if (psn_total[n] > 0)
					ps_total[n] /= (psn_total[n]);	
				else 
					ps_total[n] = 0.0;
			}

			fd = fopen(fname, "w");
			if (fd == NULL) {
				printf( " Error ! cannot open powspec file \n"  );
				exit(0);
			}
			float vol = BOXSIZE/1000.0;
			vol = vol * vol * vol;
			float Norm = vol/((float)NSIDEMESH*NSIDEMESH *NSIDEMESH);

			int output_num = (int)(NSIDEMESH*0.75);			
//			output_num = NSIDEMESH/2;			

			for (n=1; n<output_num; n++) {
				fprintf(fd, "%e %e %lf\n",  (n+0.5)*2.0*M_PI/(BOXSIZE/1000.0), ps_total[n] * Norm *Norm /vol, psn_total[n]  );
			}

			fclose(fd);
		}



		if ( 0 == PROC_RANK){
			myfree(ps_total, 267);
			myfree(psn_total, 268);
		}
		myfree(psn, 264);
		myfree(psn_rank, 266);
		myfree(ps, 265);
	}
///	if ( COM_DOM ) {	
//
//		myfree(pmesh, 261);
//	}

	if (  isSUDOM ) {
//		myfree(mbuff, 262);
//		myfree(smesh, 260);
	}

	//	MPI_Barrier(MPI_COMM_WORLD);
	myfree(mesh, 263);
}


void power_spectrum(int rank_snap) {

#ifdef POW_SPEC
	power_spec(rank_snap);
#endif

#ifdef POW_SPEC_DAUB
	power_spec_daub(rank_snap);
#endif

}





#endif



