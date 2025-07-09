#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <limits.h>
#include "photoNs.h"

int OUTPUT_DENSITY_MESH;

#ifndef INTXYZ

void partmesh_kick(cpuType dkick , int mode_alloc ){
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

	
	if (100 == PROC_RANK) printf("partmesh_kcik %f\n", dkick);
	dtime_pm = dtime();
	if (  isSUDOM ) {
		int nx, ny, nz;
		nx = NSIDEMESH;
		ny = NSIDEMESH/NSIDE0SUDOM;
		nz = NSIDEMESH/NSIDE1SUDOM;
		meshsize_sud = nx * ny * nz;
		//		meshsize_proc = nx * ny * nz;
		smesh = (cpuType*)mymalloc(sizeof(cpuType)* meshsize_sud, 101) ;

		if (smesh == NULL) {
			printf(" ERROR : cannot open smesh  \n");
			fflush(stdout);
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
		meshsize_pot = (isize[0] + 2*NPAD)  * (isize[1]+2*NPAD) * (isize[2]+2*NPAD);
		pmesh = (cpuType*)mymalloc(sizeof(cpuType) * meshsize, 102);

		if (pmesh == NULL) {
			printf(" ERROR : cannot open pmesh  \n");
			fflush(stdout);
			exit(0);
		}

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


	myfree(pmesh, 102);
	}

	//////////////////

	if (  isSUDOM ) {
		mbuff = (cpuType*)mymalloc(sizeof(cpuType)*meshsize_sud, 103);
		int m;

		if (mbuff == NULL) {
			printf(" ERROR : cannot open mbuff  \n");
			fflush(stdout);
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


#ifdef MESHOUT
		if ( 1 ==  OUTPUT_DENSITY_MESH  ) {

			int ntg = 2;
			int nx = NSIDEMESH;
			int ny = NSIDEMESH/NSIDE0SUDOM;
			int nz = NSIDEMESH/NSIDE1SUDOM;

			float mean_dens = MASSPART / pow(BOXSIZE /((float)NSIDEMESH), 3.0) ;

			char fname_mesh[128];
			sprintf(fname_mesh, "%s.%d", PathMeshOut, SUDOM_RANK );
			FILE *fom = fopen(fname_mesh,"w");

			int y_rank = SUDOM_RANK / NSIDE1SUDOM;
			int z_rank = SUDOM_RANK % NSIDE1SUDOM;
			//	fprintf(fom, "%d %d %d ", z_rank * nz/ntg, y_rank * ny/ntg, 0);
			//	fprintf(fom, "%d %d %d ", (z_rank + 1) * nz/ntg, (y_rank + 1 ) * ny/ntg, nx/ntg );

			int sx, sy, sz;
			int ex, ey, ez;

			sx = 0;
			ex = nx / ntg;

			sy = y_rank * ny/ntg;
			ey =(y_rank+1) * ny/ntg;

			sz = z_rank * nz/ntg;
			ez =(z_rank+1) * nz/ntg;

			fwrite(&sx, sizeof(int), 1, fom);
			fwrite(&sy, sizeof(int), 1, fom);
			fwrite(&sz, sizeof(int), 1, fom);
			fwrite(&ex, sizeof(int), 1, fom);
			fwrite(&ey, sizeof(int), 1, fom);
			fwrite(&ez, sizeof(int), 1, fom);

			cpuType renormal = (NSIDEMESH/BOXSIZE);
			renormal = renormal*renormal*renormal;
			renormal *= MASSPART*((float) ntg *ntg *ntg);

			renormal = 1.0/renormal;


			for (i=0; i<nx; i+= ntg) {
				for (j=0; j<ny; j+= ntg) {
					for (k=0; k<nz; k+= ntg) {
						int di, dj, dk;
						float val = 0.0;
						for (dk=0; dk<ntg; dk++) {
							for (dj=0; dj<ntg; dj++) {	
								for (di=0; di<ntg; di++) {	

									m =0;
									while (i+di < ndisp[m]) m++;

									int im = i+di - ndisp[m]; 

									int idxm = ( im*ny + (j+dj) )*nz + k+dk;
									sm = smesh + dom_grp_mesh_start[m] * y_n * z_n;
									if (idxm < meshsize_sud)
										val +=(float) sm[idxm];
								}
							}
						}
						val *= renormal;


						//	fprintf(fom,"%e ", val);
						fwrite(&val, sizeof(float), 1, fom);
					}
				}
			}
			fclose(fom);
			OUTPUT_DENSITY_MESH = 0;
		}
#endif


		myfree(smesh, 101);


		for (m=0; m<NDOMinSUDOM; m++) {	
			MPI_Isend( mbuff + m*nm, nm, MPI_CPUTYPE, dom_grp[m], 200, MPI_COMM_WORLD, &(request[m]));
		}
	}


	int x_n = NSIDEMESH;
	int y_n = NSIDEMESH/NSIDE0PROC;
	int z_n = NSIDEMESH/NSIDE1PROC;
	meshsize_proc =  x_n  * y_n * z_n;
	mesh = (cpuType*)mymalloc(sizeof(cpuType) * meshsize_proc, 104);

	if (mesh == NULL) {
		printf(" ERROR : cannot open mesh  \n");
		fflush(stdout);
		exit(0);
	}


	MPI_Recv(mesh, 	meshsize_proc, MPI_CPUTYPE, PROC_RANK_SUDOM, 200, MPI_COMM_WORLD, &sta); 
	if (  isSUDOM ) {
		int m;	
		for (m=0; m<NDOMinSUDOM; m++) {	
			MPI_Wait(&(request[m]), &(status[m]) );
		}
		myfree(mbuff, 103);
	}
	/// canbe free mbuff and smesh here once
	///////////////////////
	int nside[3] = {NSIDEMESH, NSIDEMESH, NSIDEMESH};
	cpuType param[2] = { splitRadius, BOXSIZE};

	MPI_Barrier(MPI_COMM_WORLD);

	double tcnv1 = dtime();
	{	
		convolution(mesh,nside,param);
	}
	double tcnv2 = dtime();
	//	if (print_msg())

	if ( 0 == PROC_RANK)
		printf("[%d] Convolution time %lf.\n", PROC_RANK, tcnv2 - tcnv1);


	{
		MPI_Isend( mesh, meshsize_proc, MPI_CPUTYPE, PROC_RANK_SUDOM, 300, MPI_COMM_WORLD, &req )  ;
	}
	if ( isSUDOM ) {
		int m;
		mbuff = (cpuType*)mymalloc(sizeof(cpuType)*meshsize_sud, 103);
		for (m=0; m<NDOMinSUDOM; m++) {	
			MPI_Recv( mbuff + m*meshsize_proc, meshsize_proc, MPI_CPUTYPE, dom_grp[m], 300, MPI_COMM_WORLD, &(status[m]) )  ;
		}

	}
	{
		MPI_Wait(&req, &sta );
	}

	MPI_Barrier(MPI_COMM_WORLD);
	myfree(mesh, 104);
	//////////////////////////



	if (  isSUDOM ) {
		int m;

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

		smesh = (cpuType*)mymalloc(sizeof(cpuType)* meshsize_sud, 101) ;


		cpuType *sm;
		for (m=0; m<NDOMinSUDOM; m++) {
			int x_n = dom_grp_mesh_size[m];

			sm = smesh + dom_grp_mesh_start[m] * y_n * z_n;
			//	sm = smesh + m * x_n * y_n * z_n;
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

						sm[(i*y_n+j)*z_n+k] = mbuff[ig*nm + idx];


					}
				}
			}


		}

		myfree(mbuff, 103);


		int local_size_y =  NSIDEMESH/NSIDE0SUDOM; 
		int local_size_z =  NSIDEMESH/NSIDE1SUDOM;
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			nrecv[n] = local_size_y  * local_size_z * dom_grp_mesh_size[n] ;	
		}
		ip = 0;
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			MPI_Isend( smesh + ip , nrecv[n], MPI_CPUTYPE, dom_grp[n], 100, MPI_COMM_WORLD, &request[n]);
			ip += nrecv[n];
		}
	}
	if (COM_DOM) {
	
	pmesh = (cpuType*)mymalloc(sizeof(cpuType) * meshsize, 102);
		MPI_Recv( pmesh, meshsize, MPI_CPUTYPE, PROC_RANK_SUDOM, 100, MPI_COMM_WORLD, &sta);
	}

	if ( isSUDOM ){
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			MPI_Wait( &request[n], &status[n]);
		}

//		myfree(smesh, 101);
	}



	if (COM_DOM) {
		gmesh = (cpuType*)mymalloc(sizeof(cpuType) * meshsize_pot, 105);		

		if (gmesh == NULL) {
			printf(" ERROR : cannot open gmesh  \n");
			fflush(stdout);
			exit(0);
		}


		int i,j,k;
		int ii, jj, kk;
		int gsize[3];
		gsize[0] = isize[0] + 2*NPAD;
		gsize[1] = isize[1] + 2*NPAD;
		gsize[2] = isize[2] + 2*NPAD;
		int c =0 ;
		for (i=0; i<isize[0]; i++){
			ii = i + NPAD;
			for (j=0; j<isize[1]; j++){
				jj = j + NPAD;
				for (k=0; k<isize[2]; k++){
					kk = k+NPAD;
					gmesh[ (ii*gsize[1]+jj)*gsize[2] + kk] = pmesh[c];
					c++;
				}
			}
		}
	}

	if ( COM_DOM ) {	

		myfree(pmesh, 102);
	}

	if (  isSUDOM ) {
//		myfree(mbuff, 103);
	}



	if (isSUDOM ) {
		int m, i,j,k;
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

		int nmrecv[8];

		MeshBnd *cmesh_ex = (MeshBnd*)mymalloc(sizeof(MeshBnd) *( NSIDEMESH + 2*NPAD) * 2 * NPAD * ( z_n + y_n + 2*NPAD ), 106);	

		if (cmesh_ex == NULL) {
			printf(" ERROR : cannot open cmesh_ex  \n");
			fflush(stdout);
			exit(0);
		}


		///////////////////		
		cpuType *sm;
		nmrecv[0] = NPAD*NPAD*(NSIDEMESH+2*NPAD);
		nmrecv[2] = NPAD*NPAD*(NSIDEMESH+2*NPAD);
		nmrecv[5] = NPAD*NPAD*(NSIDEMESH+2*NPAD);
		nmrecv[7] = NPAD*NPAD*(NSIDEMESH+2*NPAD);
		//		nmrecv[1] = (NSIDEMESH+2*NPAD) * y_n * NPAD;
		//		nmrecv[6] = (NSIDEMESH+2*NPAD) * y_n * NPAD;
		//		nmrecv[3] = (NSIDEMESH+2*NPAD) * z_n * NPAD;
		//		nmrecv[4] = (NSIDEMESH+2*NPAD) * z_n * NPAD;

		nmrecv[1] = (NSIDEMESH+2*NPAD) * z_n * NPAD;
		nmrecv[6] = (NSIDEMESH+2*NPAD) * z_n * NPAD;
		nmrecv[3] = (NSIDEMESH+2*NPAD) * y_n * NPAD;
		nmrecv[4] = (NSIDEMESH+2*NPAD) * y_n * NPAD;

		int istart[8];
		istart[0] = 0;
		for (m=1; m<8; m++)
			istart[m] = istart[m-1] + nmrecv[m-1];

		int nmrecv_total = istart[7] + nmrecv[7];

		nmrecv[0] =0; 
		nmrecv[2] =0;
		nmrecv[5] =0;
		nmrecv[7] =0;
		nmrecv[1] =0;
		nmrecv[6] =0;
		nmrecv[3] =0;
		nmrecv[4] =0;


		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {
			int x_n = dom_grp_mesh_size[m];

			sm = smesh + dom_grp_mesh_start[m] * y_n * z_n;
			int ii, jj, kk;
			int in, jn, kn;
			int ig, idx;

			for (i=0; i<x_n; i++) {
				for (j=0; j<y_n; j++) {
					for (k=0; k<z_n; k++) {
						in = dom_grp_mesh_start[m] + i;
						jn = j;
						kn = k;
						if (j<NPAD && k<NPAD )  {
							idx = 0;
							ip = istart[idx] + nmrecv[idx];
							jn = j + y_n;
							kn = k + z_n;
							cmesh_ex[ ip ].x = in ; 
							cmesh_ex[ ip ].y = jn ; 
							cmesh_ex[ ip ].z = kn ; 
							cmesh_ex[ ip ].data = sm[(i*y_n+j)*z_n+k]; 
							nmrecv[idx]++;
						}

						if (k<NPAD  )  {
							idx = 3;
							ip = istart[idx] + nmrecv[idx];
							jn = j ;
							kn = k + z_n;
							cmesh_ex[ ip ].x = in ; 
							cmesh_ex[ ip ].y = jn ; 
							cmesh_ex[ ip ].z = kn ; 
							cmesh_ex[ ip ].data = sm[(i*y_n+j)*z_n+k]; 
							nmrecv[idx]++;
						}

						if (k<NPAD && j>= y_n -NPAD )  {
							idx = 5;
							ip = istart[idx] + nmrecv[idx];
							jn = j -  y_n;
							kn = k + z_n;
							cmesh_ex[ ip ].x = in ; 
							cmesh_ex[ ip ].y = jn ; 
							cmesh_ex[ ip ].z = kn ; 
							cmesh_ex[ ip ].data = sm[(i*y_n+j)*z_n+k]; 
							nmrecv[idx]++;
						}

						if (j<NPAD )  {
							idx = 1;
							ip = istart[idx] + nmrecv[idx];
							jn = j + y_n;
							kn = k ;
							cmesh_ex[ ip ].x = in ; 
							cmesh_ex[ ip ].y = jn ; 
							cmesh_ex[ ip ].z = kn ; 
							cmesh_ex[ ip ].data = sm[(i*y_n+j)*z_n+k]; 
							nmrecv[idx]++;
						}

						if (j>= y_n -NPAD )  {
							idx = 6;
							ip = istart[idx] + nmrecv[idx];
							jn = j - y_n;
							kn = k ;
							cmesh_ex[ ip ].x = in ; 
							cmesh_ex[ ip ].y = jn ; 
							cmesh_ex[ ip ].z = kn ; 
							cmesh_ex[ ip ].data = sm[(i*y_n+j)*z_n+k]; 
							nmrecv[idx]++;
						}

						if (j<NPAD && k>=z_n-NPAD )  {
							idx = 2;
							ip = istart[idx] + nmrecv[idx];
							jn = j + y_n;
							kn = k - z_n;
							cmesh_ex[ ip ].x = in ; 
							cmesh_ex[ ip ].y = jn ; 
							cmesh_ex[ ip ].z = kn ; 
							cmesh_ex[ ip ].data = sm[(i*y_n+j)*z_n+k]; 
							nmrecv[idx]++;
						}

						if (k>=z_n-NPAD )  {
							idx = 4;
							ip = istart[idx] + nmrecv[idx];
							jn = j ;
							kn = k - z_n;
							cmesh_ex[ ip ].x = in ; 
							cmesh_ex[ ip ].y = jn ; 
							cmesh_ex[ ip ].z = kn ; 
							cmesh_ex[ ip ].data = sm[(i*y_n+j)*z_n+k]; 
							nmrecv[idx]++;
						}

						if (k>=z_n-NPAD && j>= y_n -NPAD )  {
							idx = 7;
							ip = istart[idx] + nmrecv[idx];
							jn = j - y_n;
							kn = k - z_n;
							cmesh_ex[ ip ].x = in ; 
							cmesh_ex[ ip ].y = jn ; 
							cmesh_ex[ ip ].z = kn ; 
							cmesh_ex[ ip ].data = sm[(i*y_n+j)*z_n+k]; 
							nmrecv[idx]++;
						}

					}
				}
			}


		}



		for (n=0; n<8; n++) {
			int len = nmrecv[n];
			int im, in, i;
			for (i=0; i<len; i++) {
				im = istart[n] + i;
				if ( cmesh_ex[im].x < NPAD ) {
					in = cmesh_ex[im].x + NSIDEMESH;
					ip = istart[n] + nmrecv[n];
					cmesh_ex[ ip ].x = in;
					cmesh_ex[ ip ].y = cmesh_ex[im].y;
					cmesh_ex[ ip ].z = cmesh_ex[im].z;
					cmesh_ex[ ip ].data = cmesh_ex[im].data;
					nmrecv[n]++;
				}	
				if ( cmesh_ex[im].x >= NSIDEMESH-NPAD ) {
					in = cmesh_ex[im].x - NSIDEMESH;
					ip = istart[n] + nmrecv[n];
					cmesh_ex[ ip ].x = in;
					cmesh_ex[ ip ].y = cmesh_ex[im].y;
					cmesh_ex[ ip ].z = cmesh_ex[im].z;
					cmesh_ex[ ip ].data = cmesh_ex[im].data;

					nmrecv[n]++;
				}	

			}
		}


		//		for (m=0; m<8; m++)
		//		printf(" 2 check : [%d] (%d) nmrecv= %d\n", PROC_RANK, m, nmrecv[m]);

		MPI_Status  status2[8], status3[8];
		MPI_Request request3[8];

		int Nmesh_bnd_sud = (y_n+2*NPAD)*(z_n+2*NPAD)*(NSIDEMESH+2*NPAD) + y_n * z_n * 2*NPAD *(NDOMinSUDOM-NDOM_HEAD ) - y_n*z_n*NSIDEMESH;
		//	int Nmesh_bnd_sud2 = (y_n+2*NPAD)*(z_n+2*NPAD)*(NSIDEMESH+2*NPAD) - y_n*z_n*(2*NPAD+NSIDEMESH);
		//		int Nmesh_proc = Nmesh_bnd_sud / (NDOMinSUDOM - NRPOC_GAP);

		cmesh = (MeshBnd*)mymalloc(sizeof(MeshBnd) * Nmesh_bnd_sud, 107);	

		if (cmesh == NULL) {
			printf(" ERROR : cannot open cmesh  \n");
			fflush(stdout);
			exit(0);
		}


		for (m=0; m<8; m++) {
			if (nmrecv[m] > INT_MAX / sizeof(MeshBnd)) {
				printf("[%d] Error %s:%d! count of MPI_Isend overflow in INT32.\n", PROC_RANK, __FILE__, __LINE__);
				exit(0);
			}
			MPI_Isend(cmesh_ex+istart[m], sizeof(MeshBnd)*nmrecv[m], MPI_BYTE, PROC_RANK_SEND_NGB[m], 2, MPI_COMM_WORLD, &(request3[m]) );

		}

		for (m=0; m<8; m++) {
			MPI_Recv(cmesh + istart[m] , sizeof(MeshBnd)*nmrecv[m], MPI_BYTE, PROC_RANK_RECV_NGB[m], 2, MPI_COMM_WORLD, &(status2[m]));
			//			printf("%d:  nrecv[%d] = %d\n", SUDOM_RANK, m, nrecv[m]);
		}	
		for (m=0; m<8; m++) {
			MPI_Wait(&request3[m], &status3[m]);		
		}


		myfree(cmesh_ex, 106);
		///////////////////////return ;
		ip = 0;
		for (m=0; m<8; m++)
			ip += nmrecv[m];

		ip = nmrecv_total;
		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {

			int x_n = dom_grp_mesh_size[m];

			sm = smesh + dom_grp_mesh_start[m] * y_n * z_n;
			int ii, jj, kk;
			int in, jn, kn;
			int ig, idx;
			for (i=0; i<NPAD; i++) {
				for (j=0; j<y_n; j++) {
					for (k=0; k<z_n; k++) {
						in = dom_grp_mesh_start[m] + i;
						jn = j;
						kn = k;
						if (in < NPAD)
							in += NSIDEMESH;

						cmesh[ ip ].x = in ; 
						cmesh[ ip ].y = jn ; 
						cmesh[ ip ].z = kn ; 
						cmesh[ ip ].data = sm[(i*y_n+j)*z_n+k]; 

						ip ++;


					}
				}
			}
			for (i=dom_grp_mesh_size[m]-NPAD; i<dom_grp_mesh_size[m]; i++) {
				for (j=0; j<y_n; j++) {
					for (k=0; k<z_n; k++) {
						in = dom_grp_mesh_start[m] + i;
						jn = j;
						kn = k;
						if (in >= NSIDEMESH-NPAD)
							in -= NSIDEMESH;

						cmesh[ ip ].x = in ; 
						cmesh[ ip ].y = jn ; 
						cmesh[ ip ].z = kn ; 
						cmesh[ ip ].data = sm[(i*y_n+j)*z_n+k]; 
						ip ++;
					}
				}
			}



		}
		myfree(smesh, 101);
		int nmsend[NDOMinSUDOM] ;
		int imsend[NDOMinSUDOM+1] ;
		int len_bmesh_sud = 0;

		for (n=0; n<NDOMinSUDOM; n++ ) {
			nmsend[n] = 0;
			imsend[n] = 0;
		}
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++ ){	
			nmsend[n] = (dom_grp_mesh_size[n] + 2*NPAD ) * (y_n + 2*NPAD) *(z_n + 2*NPAD) - dom_grp_mesh_size[n]*y_n *z_n;	
			imsend[n+1] = imsend[n]+ nmsend[n];
		}
		len_bmesh_sud = imsend[NDOMinSUDOM];

		bmesh = (MeshBnd*)mymalloc( sizeof(MeshBnd) * len_bmesh_sud, 108);	

		if (bmesh == NULL) {
			printf(" ERROR : cannot open bmesh  \n");
			fflush(stdout);
			exit(0);
		}

		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {
			ip = imsend[m];
			for (n=0; n<Nmesh_bnd_sud; n++){
				if (cmesh[n].x >= dom_grp_mesh_start[m] - NPAD && cmesh[n].x <dom_grp_mesh_start[m] ) {
					bmesh[ip] = cmesh[n];
					bmesh[ip].x -= dom_grp_mesh_start[m];
					ip ++;
				}
				if (cmesh[n].x >= dom_grp_mesh_end[m] && cmesh[n].x < dom_grp_mesh_end[m]+NPAD ) {
					bmesh[ip] = cmesh[n];
					bmesh[ip].x -= dom_grp_mesh_start[m];
					ip ++;
				}

				if (cmesh[n].x >= dom_grp_mesh_start[m] && cmesh[n].x < dom_grp_mesh_end[m] ) {
					if (cmesh[n].y < 0 || cmesh[n].y >= y_n	|| cmesh[n].z < 0 || cmesh[n].z >= z_n ) 
					{

						bmesh[ip] = cmesh[n];
						bmesh[ip].x -= dom_grp_mesh_start[m];
						ip ++;
					}

				}
			}	

			if (ip - imsend[m] != nmsend[m])	
				printf("[%d] ip = %d, nmsend = %d\n", PROC_RANK, ip - imsend[m], nmsend[m]);	
			// MPI_Isend;
			//
		}

		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {

			MPI_Isend( bmesh + imsend[m], sizeof(MeshBnd)*nmsend[m], MPI_BYTE, dom_grp[m], 400, MPI_COMM_WORLD, &(request[m]) );
		}

	}
	////////////////////////////

	//fflush(stdout);
	//MPI_Barrier(MPI_COMM_WORLD);

	if (COM_DOM) {
		int y_n = NSIDEMESH/NSIDE0SUDOM;
		int z_n = NSIDEMESH/NSIDE1SUDOM;
		int n = dom_grp_rank;
		int nmrecv = (dom_grp_mesh_size[n] + 2*NPAD ) * (y_n + 2*NPAD) *(z_n + 2*NPAD) - dom_grp_mesh_size[n]*y_n *z_n;




		//		printf(" %d > %d %d \n",PROC_RANK, n, nmrecv  );
		bmeshc = (MeshBnd*)mymalloc(sizeof(MeshBnd) * nmrecv, 109);

		if (bmeshc == NULL) {
			printf(" ERROR : cannot open bmeshc  \n");
			fflush(stdout);
			exit(0);
		}





		MPI_Recv( bmeshc, sizeof(MeshBnd) * nmrecv, MPI_BYTE, PROC_RANK_SUDOM, 400, MPI_COMM_WORLD, &sta);
		int i,j,k;
		int ny = y_n + 2*NPAD;
		int nz = z_n + 2*NPAD;
		int idx;
		for (n=0; n<nmrecv; n++) {
			i = bmeshc[n].x + NPAD;			
			j = bmeshc[n].y + NPAD;			
			k = bmeshc[n].z + NPAD;
			idx = (i * ny + j) *nz + k;

			gmesh[idx] = bmeshc[n].data;				
		}
	}



	if (  isSUDOM ) {
		int m;
		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {
			MPI_Wait( &(request[m]), &(status[m]) );
		}
	}

	if (  isSUDOM ) {
		myfree(bmesh, 108);
		myfree(cmesh, 107);
	}

//MPI_Barrier(MPI_COMM_WORLD);
	if (COM_DOM) {
		myfree(bmeshc, 109);
	}

	if (COM_DOM) {
		cpuType *mesh = gmesh;
		int n;
		int i,j,k,idx;
		int ii, jj, kk;
		cpuType norm = NSIDEMESH/BOXSIZE;
		cpuType delta = 1.0/norm;

		cpuType invx = 0.5*NSIDEMESH/BOXSIZE;
		isize[0] += 2*NPAD   ;
		isize[1] += 2*NPAD  ;
		isize[2] += 2*NPAD  ;
		meshsize = isize[0] * isize[1] * isize[2];
		cpuType wi, wj, wk;
		cpuType win, wjn, wkn;	
		int nx, ny, nz;

		nx = isize[1]*isize[2];
		ny = isize[2];
		nz = 1;

		local_min[0] -= NPAD;
		local_min[1] -= NPAD;
		local_min[2] -= NPAD;

		cpuType dp[8];
		int idx1, idx2;

		if (mode_alloc) {

			acc[0] = (cpuType*)mymalloc(sizeof(cpuType)*NPART, 14);
			acc[1] = (cpuType*)mymalloc(sizeof(cpuType)*NPART, 15);
			acc[2] = (cpuType*)mymalloc(sizeof(cpuType)*NPART, 16);

			
			for (n=0; n<NPART; n++) {
				acc[0][n] = 0.0;
				acc[1][n] = 0.0;
				acc[2][n] = 0.0;
			}
		}
		cpuType acc_t[3];
		for (n=0; n<NPART; n++) {
			if (part[n].tag != 1) continue; // qwang
			i = (int) (part[n].pos[0] *norm) ;	
			j = (int) (part[n].pos[1] *norm) ;	
			k = (int) (part[n].pos[2] *norm) ;	

			acc_t[0] = acc_t[1] =acc_t[2]  = 0.0;

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

			idx = (i*isize[1] + j)*isize[2] + k;


			if (idx + nx > meshsize)
				printf(" partmesh 1 error : %d %d %d, %lf %lf %lf\n", i, j, k, part[n].pos[0], part[n].pos[1], part[n].pos[2]);

			if (idx - nx < 0)
				printf(" partmesh 2 error : %d %d %d, %lf %lf %lf\n", i, j, k, part[n].pos[0], part[n].pos[1], part[n].pos[2]);


			fflush(stdout);
			cpuType f1 = 4.0/3.0;
			cpuType f2 = 1.0/6.0;
			idx1 = ((i-1)*isize[1] + j)*isize[2] + k;
			idx2 = ((i+1)*isize[1] + j)*isize[2] + k;
			dp[0] =  f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-2)*isize[1] + j)*isize[2] + k;
			idx2 = ((i+2)*isize[1] + j)*isize[2] + k;
			dp[0] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-1)*isize[1] + j)*isize[2] + k;
			idx2 = ((ii+1)*isize[1] + j)*isize[2] + k;
			dp[1] = f1 *invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-2)*isize[1] + j)*isize[2] + k;
			idx2 = ((ii+2)*isize[1] + j)*isize[2] + k;
			dp[1] -= f2 * invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-1)*isize[1] + jj)*isize[2] + k;
			idx2 = ((i+1)*isize[1] + jj)*isize[2] + k;
			dp[2] = f1 * invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-2)*isize[1] + jj)*isize[2] + k;
			idx2 = ((i+2)*isize[1] + jj)*isize[2] + k;
			dp[2] -= f2 * invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-1)*isize[1] + jj)*isize[2] + k;
			idx2 = ((ii+1)*isize[1] + jj)*isize[2] + k;
			dp[3] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-2)*isize[1] + jj)*isize[2] + k;
			idx2 = ((ii+2)*isize[1] + jj)*isize[2] + k;
			dp[3] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-1)*isize[1] + j)*isize[2] + kk;
			idx2 = ((i+1)*isize[1] + j)*isize[2] + kk;
			dp[4] = f1* invx*(mesh[idx2] - mesh[idx1]);


			idx1 = ((i-2)*isize[1] + j)*isize[2] + kk;
			idx2 = ((i+2)*isize[1] + j)*isize[2] + kk;
			dp[4] -= f2 *invx*(mesh[idx2] - mesh[idx1]);
			idx1 = ((ii-1)*isize[1] + j)*isize[2] + kk;
			idx2 = ((ii+1)*isize[1] + j)*isize[2] + kk;
			dp[5] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-2)*isize[1] + j)*isize[2] + kk;
			idx2 = ((ii+2)*isize[1] + j)*isize[2] + kk;
			dp[5] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-1)*isize[1] + jj)*isize[2] + kk;
			idx2 = ((i+1)*isize[1] + jj)*isize[2] + kk;
			dp[6] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-2)*isize[1] + jj)*isize[2] + kk;
			idx2 = ((i+2)*isize[1] + jj)*isize[2] + kk;
			dp[6]-= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-1)*isize[1] + jj)*isize[2] + kk;
			idx2 = ((ii+1)*isize[1] + jj)*isize[2] + kk;
			dp[7] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-2)*isize[1] + jj)*isize[2] + kk;
			idx2 = ((ii+2)*isize[1] + jj)*isize[2] + kk;
			dp[7]-= f2*invx*(mesh[idx2] - mesh[idx1]);

			acc_t[0]= win*wjn*wkn*dp[0]
				+ wi *wjn*wkn*dp[1]
				+ win*wj *wkn*dp[2]
				+ wi *wj *wkn*dp[3]
				+ win*wjn*wk *dp[4]
				+ wi *wjn*wk *dp[5]
				+ win*wj *wk *dp[6]
				+ wi *wj *wk *dp[7];

			idx1 = (i*isize[1] + j-1)*isize[2] + k;
			idx2 = (i*isize[1] + j+1)*isize[2] + k;
			dp[0] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + j-2)*isize[2] + k;
			idx2 = (i*isize[1] + j+2)*isize[2] + k;
			dp[0] -= f2*invx*(mesh[idx2] - mesh[idx1]);


			idx1 = (ii*isize[1] + j-1)*isize[2] + k;
			idx2 = (ii*isize[1] + j+1)*isize[2] + k;
			dp[1] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j-2)*isize[2] + k;
			idx2 = (ii*isize[1] + j+2)*isize[2] + k;
			dp[1] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + jj-1)*isize[2] + k;
			idx2 = (i*isize[1] + jj+1)*isize[2] + k;
			dp[2] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + jj-2)*isize[2] + k;
			idx2 = (i*isize[1] + jj+2)*isize[2] + k;
			dp[2] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj-1)*isize[2] + k;
			idx2 = (ii*isize[1] + jj+1)*isize[2] + k;
			dp[3] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj-2)*isize[2] + k;
			idx2 = (ii*isize[1] + jj+2)*isize[2] + k;
			dp[3] -= f2*invx*(mesh[idx2] - mesh[idx1]);


			idx1 = (i*isize[1] + j-1)*isize[2] + kk;
			idx2 = (i*isize[1] + j+1)*isize[2] + kk;
			dp[4] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + j-2)*isize[2] + kk;
			idx2 = (i*isize[1] + j+2)*isize[2] + kk;
			dp[4] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j-1)*isize[2] + kk;
			idx2 = (ii*isize[1] + j+1)*isize[2] + kk;
			dp[5] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j-2)*isize[2] + kk;
			idx2 = (ii*isize[1] + j+2)*isize[2] + kk;
			dp[5] -= f2*invx*(mesh[idx2] - mesh[idx1]);


			idx1 = (i*isize[1] + jj-1)*isize[2] + kk;
			idx2 = (i*isize[1] + jj+1)*isize[2] + kk;
			dp[6] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + jj-2)*isize[2] + kk;
			idx2 = (i*isize[1] + jj+2)*isize[2] + kk;
			dp[6] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj-1)*isize[2] + kk;
			idx2 = (ii*isize[1] + jj+1)*isize[2] + kk;
			dp[7] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj-2)*isize[2] + kk;
			idx2 = (ii*isize[1] + jj+2)*isize[2] + kk;
			dp[7] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			acc_t[1]= win*wjn*wkn*dp[0]
				+ wi *wjn*wkn*dp[1]
				+ win*wj *wkn*dp[2]
				+ wi *wj *wkn*dp[3]
				+ win*wjn*wk *dp[4]
				+ wi *wjn*wk *dp[5]
				+ win*wj *wk *dp[6]
				+ wi *wj *wk *dp[7];



			idx1 = (i*isize[1] + j)*isize[2] + k-1;
			idx2 = (i*isize[1] + j)*isize[2] + k+1;
			dp[0] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + j)*isize[2] + k-2;
			idx2 = (i*isize[1] + j)*isize[2] + k+2;
			dp[0] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j)*isize[2] + k-1;
			idx2 = (ii*isize[1] + j)*isize[2] + k+1;
			dp[1] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j)*isize[2] + k-2;
			idx2 = (ii*isize[1] + j)*isize[2] + k+2;
			dp[1] -= f2*invx*(mesh[idx2] - mesh[idx1]);


			idx1 = (i*isize[1] + jj)*isize[2] + k-1;
			idx2 = (i*isize[1] + jj)*isize[2] + k+1;
			dp[2] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + jj)*isize[2] + k-2;
			idx2 = (i*isize[1] + jj)*isize[2] + k+2;
			dp[2] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj)*isize[2] + k-1;
			idx2 = (ii*isize[1] + jj)*isize[2] + k+1;
			dp[3] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj)*isize[2] + k-2;
			idx2 = (ii*isize[1] + jj)*isize[2] + k+2;
			dp[3] -= f2*invx*(mesh[idx2] - mesh[idx1]);


			idx1 = (i*isize[1] + j)*isize[2] + kk-1;
			idx2 = (i*isize[1] + j)*isize[2] + kk+1;
			dp[4] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + j)*isize[2] + kk-2;
			idx2 = (i*isize[1] + j)*isize[2] + kk+2;
			dp[4] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j)*isize[2] + kk-1;
			idx2 = (ii*isize[1] + j)*isize[2] + kk+1;
			dp[5] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j)*isize[2] + kk-2;
			idx2 = (ii*isize[1] + j)*isize[2] + kk+2;
			dp[5] -= f2*invx*(mesh[idx2] - mesh[idx1]);


			idx1 = (i*isize[1] + jj)*isize[2] + kk-1;
			idx2 = (i*isize[1] + jj)*isize[2] + kk+1;
			dp[6] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + jj)*isize[2] + kk-2;
			idx2 = (i*isize[1] + jj)*isize[2] + kk+2;
			dp[6] -= f2*invx*(mesh[idx2] - mesh[idx1]);


			idx1 = (ii*isize[1] + jj)*isize[2] + kk-1;
			idx2 = (ii*isize[1] + jj)*isize[2] + kk+1;
			dp[7] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj)*isize[2] + kk-2;
			idx2 = (ii*isize[1] + jj)*isize[2] + kk+2;
			dp[7] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			acc_t[2]= win*wjn*wkn*dp[0]
				+ wi *wjn*wkn*dp[1]
				+ win*wj *wkn*dp[2]
				+ wi *wj *wkn*dp[3]
				+ win*wjn*wk *dp[4]
				+ wi *wjn*wk *dp[5]
				+ win*wj *wk *dp[6]
				+ wi *wj *wk *dp[7];


			part[n].vel[0] += acc_t[0] * dkick ;	
			part[n].vel[1] += acc_t[1] * dkick ;	
			part[n].vel[2] += acc_t[2] * dkick ;	

			if (mode_alloc) {
				acc[0][n] = acc_t[0];
				acc[1][n] = acc_t[1];
				acc[2][n] = acc_t[2];
			}
		}




	}
	if ( COM_DOM ) {	
//		myfree(bmeshc, 109);

		myfree(gmesh, 105);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if ( 0 == PROC_RANK ) {
		printf( "\n partmesh complete !\n" );
		fflush(stdout);
	}
	dtime_pm = dtime() - dtime_pm;
}





#else

void partmesh_kick(cpuType dkick , int mode_alloc ){
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

	
	if (100 == PROC_RANK) printf("partmesh_kcik %f\n", dkick);
	dtime_pm = dtime();
	if (  isSUDOM ) {
		int nx, ny, nz;
		nx = NSIDEMESH;
		ny = NSIDEMESH/NSIDE0SUDOM;
		nz = NSIDEMESH/NSIDE1SUDOM;
		meshsize_sud = nx * ny * nz;
		//		meshsize_proc = nx * ny * nz;
		smesh = (cpuType*)mymalloc(sizeof(cpuType)* meshsize_sud, 101) ;

		if (smesh == NULL) {
			printf(" ERROR : cannot open smesh  \n");
			fflush(stdout);
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
		meshsize_pot = (isize[0] + 2*NPAD)  * (isize[1]+2*NPAD) * (isize[2]+2*NPAD);
		pmesh = (cpuType*)mymalloc(sizeof(cpuType) * meshsize, 102);

		if (pmesh == NULL) {
			printf(" ERROR : cannot open pmesh  \n");
			fflush(stdout);
			exit(0);
		}

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

			if ( posf[0] < 0.0 )
				i--;

			if ( posf[1] < 0.0 )
				j--;

			if ( posf[2] < 0.0 )
				k--;

			wi = ( posf[0] - (i+0.5)*delta)*norm;

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


	myfree(pmesh, 102);
	}

	//////////////////

	if (  isSUDOM ) {
		mbuff = (cpuType*)mymalloc(sizeof(cpuType)*meshsize_sud, 103);
		int m;

		if (mbuff == NULL) {
			printf(" ERROR : cannot open mbuff  \n");
			fflush(stdout);
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


#ifdef MESHOUT
		if ( 1 ==  OUTPUT_DENSITY_MESH  ) {

			int ntg = 2;
			int nx = NSIDEMESH;
			int ny = NSIDEMESH/NSIDE0SUDOM;
			int nz = NSIDEMESH/NSIDE1SUDOM;

			float mean_dens = MASSPART / pow(BOXSIZE /((float)NSIDEMESH), 3.0) ;

			char fname_mesh[128];
			sprintf(fname_mesh, "%s.%d", PathMeshOut, SUDOM_RANK );
			FILE *fom = fopen(fname_mesh,"w");

			int y_rank = SUDOM_RANK / NSIDE1SUDOM;
			int z_rank = SUDOM_RANK % NSIDE1SUDOM;
			//	fprintf(fom, "%d %d %d ", z_rank * nz/ntg, y_rank * ny/ntg, 0);
			//	fprintf(fom, "%d %d %d ", (z_rank + 1) * nz/ntg, (y_rank + 1 ) * ny/ntg, nx/ntg );

			int sx, sy, sz;
			int ex, ey, ez;

			sx = 0;
			ex = nx / ntg;

			sy = y_rank * ny/ntg;
			ey =(y_rank+1) * ny/ntg;

			sz = z_rank * nz/ntg;
			ez =(z_rank+1) * nz/ntg;

			fwrite(&sx, sizeof(int), 1, fom);
			fwrite(&sy, sizeof(int), 1, fom);
			fwrite(&sz, sizeof(int), 1, fom);
			fwrite(&ex, sizeof(int), 1, fom);
			fwrite(&ey, sizeof(int), 1, fom);
			fwrite(&ez, sizeof(int), 1, fom);

			cpuType renormal = (NSIDEMESH/BOXSIZE);
			renormal = renormal*renormal*renormal;
			renormal *= MASSPART*((float) ntg *ntg *ntg);

			renormal = 1.0/renormal;


			for (i=0; i<nx; i+= ntg) {
				for (j=0; j<ny; j+= ntg) {
					for (k=0; k<nz; k+= ntg) {
						int di, dj, dk;
						float val = 0.0;
						for (dk=0; dk<ntg; dk++) {
							for (dj=0; dj<ntg; dj++) {	
								for (di=0; di<ntg; di++) {	

									m =0;
									while (i+di < ndisp[m]) m++;

									int im = i+di - ndisp[m]; 

									int idxm = ( im*ny + (j+dj) )*nz + k+dk;
									sm = smesh + dom_grp_mesh_start[m] * y_n * z_n;
									if (idxm < meshsize_sud)
										val +=(float) sm[idxm];
								}
							}
						}
						val *= renormal;


						//	fprintf(fom,"%e ", val);
						fwrite(&val, sizeof(float), 1, fom);
					}
				}
			}
			fclose(fom);
			OUTPUT_DENSITY_MESH = 0;
		}
#endif


		myfree(smesh, 101);


		for (m=0; m<NDOMinSUDOM; m++) {	
			MPI_Isend( mbuff + m*nm, nm, MPI_CPUTYPE, dom_grp[m], 200, MPI_COMM_WORLD, &(request[m]));
		}
	}


	int x_n = NSIDEMESH;
	int y_n = NSIDEMESH/NSIDE0PROC;
	int z_n = NSIDEMESH/NSIDE1PROC;
	meshsize_proc =  x_n  * y_n * z_n;
	mesh = (cpuType*)mymalloc(sizeof(cpuType) * meshsize_proc, 104);

	if (mesh == NULL) {
		printf(" ERROR : cannot open mesh  \n");
		fflush(stdout);
		exit(0);
	}


	MPI_Recv(mesh, 	meshsize_proc, MPI_CPUTYPE, PROC_RANK_SUDOM, 200, MPI_COMM_WORLD, &sta); 
	if (  isSUDOM ) {
		int m;	
		for (m=0; m<NDOMinSUDOM; m++) {	
			MPI_Wait(&(request[m]), &(status[m]) );
		}
		myfree(mbuff, 103);
	}
	/// canbe free mbuff and smesh here once
	///////////////////////
	int nside[3] = {NSIDEMESH, NSIDEMESH, NSIDEMESH};
	cpuType param[2] = { splitRadius, BOXSIZE};

	MPI_Barrier(MPI_COMM_WORLD);

	double tcnv1 = dtime();
	{	
		convolution(mesh,nside,param);
	}
	double tcnv2 = dtime();
	//	if (print_msg())

	if ( 0 == PROC_RANK)
		printf("[%d] Convolution time %lf.\n", PROC_RANK, tcnv2 - tcnv1);


	{
		MPI_Isend( mesh, meshsize_proc, MPI_CPUTYPE, PROC_RANK_SUDOM, 300, MPI_COMM_WORLD, &req )  ;
	}
	if ( isSUDOM ) {
		int m;
		mbuff = (cpuType*)mymalloc(sizeof(cpuType)*meshsize_sud, 103);
		for (m=0; m<NDOMinSUDOM; m++) {	
			MPI_Recv( mbuff + m*meshsize_proc, meshsize_proc, MPI_CPUTYPE, dom_grp[m], 300, MPI_COMM_WORLD, &(status[m]) )  ;
		}

	}
	{
		MPI_Wait(&req, &sta );
	}

	MPI_Barrier(MPI_COMM_WORLD);
	myfree(mesh, 104);
	//////////////////////////



	if (  isSUDOM ) {
		int m;

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

		smesh = (cpuType*)mymalloc(sizeof(cpuType)* meshsize_sud, 101) ;


		cpuType *sm;
		for (m=0; m<NDOMinSUDOM; m++) {
			int x_n = dom_grp_mesh_size[m];

			sm = smesh + dom_grp_mesh_start[m] * y_n * z_n;
			//	sm = smesh + m * x_n * y_n * z_n;
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

						sm[(i*y_n+j)*z_n+k] = mbuff[ig*nm + idx];


					}
				}
			}


		}

		myfree(mbuff, 103);


		int local_size_y =  NSIDEMESH/NSIDE0SUDOM; 
		int local_size_z =  NSIDEMESH/NSIDE1SUDOM;
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			nrecv[n] = local_size_y  * local_size_z * dom_grp_mesh_size[n] ;	
		}
		ip = 0;
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			MPI_Isend( smesh + ip , nrecv[n], MPI_CPUTYPE, dom_grp[n], 100, MPI_COMM_WORLD, &request[n]);
			ip += nrecv[n];
		}
	}
	if (COM_DOM) {
	
	pmesh = (cpuType*)mymalloc(sizeof(cpuType) * meshsize, 102);
		MPI_Recv( pmesh, meshsize, MPI_CPUTYPE, PROC_RANK_SUDOM, 100, MPI_COMM_WORLD, &sta);
	}

	if ( isSUDOM ){
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++) {
			MPI_Wait( &request[n], &status[n]);
		}

//		myfree(smesh, 101);
	}



	if (COM_DOM) {
		gmesh = (cpuType*)mymalloc(sizeof(cpuType) * meshsize_pot, 105);		

		if (gmesh == NULL) {
			printf(" ERROR : cannot open gmesh  \n");
			fflush(stdout);
			exit(0);
		}


		int i,j,k;
		int ii, jj, kk;
		int gsize[3];
		gsize[0] = isize[0] + 2*NPAD;
		gsize[1] = isize[1] + 2*NPAD;
		gsize[2] = isize[2] + 2*NPAD;
		int c =0 ;
		for (i=0; i<isize[0]; i++){
			ii = i + NPAD;
			for (j=0; j<isize[1]; j++){
				jj = j + NPAD;
				for (k=0; k<isize[2]; k++){
					kk = k+NPAD;
					gmesh[ (ii*gsize[1]+jj)*gsize[2] + kk] = pmesh[c];
					c++;
				}
			}
		}
	}

	if ( COM_DOM ) {	

		myfree(pmesh, 102);
	}

	if (  isSUDOM ) {
//		myfree(mbuff, 103);
	}



	if (isSUDOM ) {
		int m, i,j,k;
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

		int nmrecv[8];

		MeshBnd *cmesh_ex = (MeshBnd*)mymalloc(sizeof(MeshBnd) *( NSIDEMESH + 2*NPAD) * 2 * NPAD * ( z_n + y_n + 2*NPAD ), 106);	

		if (cmesh_ex == NULL) {
			printf(" ERROR : cannot open cmesh_ex  \n");
			fflush(stdout);
			exit(0);
		}


		///////////////////		
		cpuType *sm;
		nmrecv[0] = NPAD*NPAD*(NSIDEMESH+2*NPAD);
		nmrecv[2] = NPAD*NPAD*(NSIDEMESH+2*NPAD);
		nmrecv[5] = NPAD*NPAD*(NSIDEMESH+2*NPAD);
		nmrecv[7] = NPAD*NPAD*(NSIDEMESH+2*NPAD);
		//		nmrecv[1] = (NSIDEMESH+2*NPAD) * y_n * NPAD;
		//		nmrecv[6] = (NSIDEMESH+2*NPAD) * y_n * NPAD;
		//		nmrecv[3] = (NSIDEMESH+2*NPAD) * z_n * NPAD;
		//		nmrecv[4] = (NSIDEMESH+2*NPAD) * z_n * NPAD;

		nmrecv[1] = (NSIDEMESH+2*NPAD) * z_n * NPAD;
		nmrecv[6] = (NSIDEMESH+2*NPAD) * z_n * NPAD;
		nmrecv[3] = (NSIDEMESH+2*NPAD) * y_n * NPAD;
		nmrecv[4] = (NSIDEMESH+2*NPAD) * y_n * NPAD;

		int istart[8];
		istart[0] = 0;
		for (m=1; m<8; m++)
			istart[m] = istart[m-1] + nmrecv[m-1];

		int nmrecv_total = istart[7] + nmrecv[7];

		nmrecv[0] =0; 
		nmrecv[2] =0;
		nmrecv[5] =0;
		nmrecv[7] =0;
		nmrecv[1] =0;
		nmrecv[6] =0;
		nmrecv[3] =0;
		nmrecv[4] =0;


		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {
			int x_n = dom_grp_mesh_size[m];

			sm = smesh + dom_grp_mesh_start[m] * y_n * z_n;
			int ii, jj, kk;
			int in, jn, kn;
			int ig, idx;

			for (i=0; i<x_n; i++) {
				for (j=0; j<y_n; j++) {
					for (k=0; k<z_n; k++) {
						in = dom_grp_mesh_start[m] + i;
						jn = j;
						kn = k;
						if (j<NPAD && k<NPAD )  {
							idx = 0;
							ip = istart[idx] + nmrecv[idx];
							jn = j + y_n;
							kn = k + z_n;
							cmesh_ex[ ip ].x = in ; 
							cmesh_ex[ ip ].y = jn ; 
							cmesh_ex[ ip ].z = kn ; 
							cmesh_ex[ ip ].data = sm[(i*y_n+j)*z_n+k]; 
							nmrecv[idx]++;
						}

						if (k<NPAD  )  {
							idx = 3;
							ip = istart[idx] + nmrecv[idx];
							jn = j ;
							kn = k + z_n;
							cmesh_ex[ ip ].x = in ; 
							cmesh_ex[ ip ].y = jn ; 
							cmesh_ex[ ip ].z = kn ; 
							cmesh_ex[ ip ].data = sm[(i*y_n+j)*z_n+k]; 
							nmrecv[idx]++;
						}

						if (k<NPAD && j>= y_n -NPAD )  {
							idx = 5;
							ip = istart[idx] + nmrecv[idx];
							jn = j -  y_n;
							kn = k + z_n;
							cmesh_ex[ ip ].x = in ; 
							cmesh_ex[ ip ].y = jn ; 
							cmesh_ex[ ip ].z = kn ; 
							cmesh_ex[ ip ].data = sm[(i*y_n+j)*z_n+k]; 
							nmrecv[idx]++;
						}

						if (j<NPAD )  {
							idx = 1;
							ip = istart[idx] + nmrecv[idx];
							jn = j + y_n;
							kn = k ;
							cmesh_ex[ ip ].x = in ; 
							cmesh_ex[ ip ].y = jn ; 
							cmesh_ex[ ip ].z = kn ; 
							cmesh_ex[ ip ].data = sm[(i*y_n+j)*z_n+k]; 
							nmrecv[idx]++;
						}

						if (j>= y_n -NPAD )  {
							idx = 6;
							ip = istart[idx] + nmrecv[idx];
							jn = j - y_n;
							kn = k ;
							cmesh_ex[ ip ].x = in ; 
							cmesh_ex[ ip ].y = jn ; 
							cmesh_ex[ ip ].z = kn ; 
							cmesh_ex[ ip ].data = sm[(i*y_n+j)*z_n+k]; 
							nmrecv[idx]++;
						}

						if (j<NPAD && k>=z_n-NPAD )  {
							idx = 2;
							ip = istart[idx] + nmrecv[idx];
							jn = j + y_n;
							kn = k - z_n;
							cmesh_ex[ ip ].x = in ; 
							cmesh_ex[ ip ].y = jn ; 
							cmesh_ex[ ip ].z = kn ; 
							cmesh_ex[ ip ].data = sm[(i*y_n+j)*z_n+k]; 
							nmrecv[idx]++;
						}

						if (k>=z_n-NPAD )  {
							idx = 4;
							ip = istart[idx] + nmrecv[idx];
							jn = j ;
							kn = k - z_n;
							cmesh_ex[ ip ].x = in ; 
							cmesh_ex[ ip ].y = jn ; 
							cmesh_ex[ ip ].z = kn ; 
							cmesh_ex[ ip ].data = sm[(i*y_n+j)*z_n+k]; 
							nmrecv[idx]++;
						}

						if (k>=z_n-NPAD && j>= y_n -NPAD )  {
							idx = 7;
							ip = istart[idx] + nmrecv[idx];
							jn = j - y_n;
							kn = k - z_n;
							cmesh_ex[ ip ].x = in ; 
							cmesh_ex[ ip ].y = jn ; 
							cmesh_ex[ ip ].z = kn ; 
							cmesh_ex[ ip ].data = sm[(i*y_n+j)*z_n+k]; 
							nmrecv[idx]++;
						}

					}
				}
			}


		}



		for (n=0; n<8; n++) {
			int len = nmrecv[n];
			int im, in, i;
			for (i=0; i<len; i++) {
				im = istart[n] + i;
				if ( cmesh_ex[im].x < NPAD ) {
					in = cmesh_ex[im].x + NSIDEMESH;
					ip = istart[n] + nmrecv[n];
					cmesh_ex[ ip ].x = in;
					cmesh_ex[ ip ].y = cmesh_ex[im].y;
					cmesh_ex[ ip ].z = cmesh_ex[im].z;
					cmesh_ex[ ip ].data = cmesh_ex[im].data;
					nmrecv[n]++;
				}	
				if ( cmesh_ex[im].x >= NSIDEMESH-NPAD ) {
					in = cmesh_ex[im].x - NSIDEMESH;
					ip = istart[n] + nmrecv[n];
					cmesh_ex[ ip ].x = in;
					cmesh_ex[ ip ].y = cmesh_ex[im].y;
					cmesh_ex[ ip ].z = cmesh_ex[im].z;
					cmesh_ex[ ip ].data = cmesh_ex[im].data;

					nmrecv[n]++;
				}	

			}
		}


		//		for (m=0; m<8; m++)
		//		printf(" 2 check : [%d] (%d) nmrecv= %d\n", PROC_RANK, m, nmrecv[m]);

		MPI_Status  status2[8], status3[8];
		MPI_Request request3[8];

		int Nmesh_bnd_sud = (y_n+2*NPAD)*(z_n+2*NPAD)*(NSIDEMESH+2*NPAD) + y_n * z_n * 2*NPAD *(NDOMinSUDOM-NDOM_HEAD ) - y_n*z_n*NSIDEMESH;
		//	int Nmesh_bnd_sud2 = (y_n+2*NPAD)*(z_n+2*NPAD)*(NSIDEMESH+2*NPAD) - y_n*z_n*(2*NPAD+NSIDEMESH);
		//		int Nmesh_proc = Nmesh_bnd_sud / (NDOMinSUDOM - NRPOC_GAP);

		cmesh = (MeshBnd*)mymalloc(sizeof(MeshBnd) * Nmesh_bnd_sud, 107);	

		if (cmesh == NULL) {
			printf(" ERROR : cannot open cmesh  \n");
			fflush(stdout);
			exit(0);
		}


		for (m=0; m<8; m++) {
			if (nmrecv[m] > INT_MAX / sizeof(MeshBnd)) {
				printf("[%d] Error %s:%d! count of MPI_Isend overflow in INT32.\n", PROC_RANK, __FILE__, __LINE__);
				exit(0);
			}
			MPI_Isend(cmesh_ex+istart[m], sizeof(MeshBnd)*nmrecv[m], MPI_BYTE, PROC_RANK_SEND_NGB[m], 2, MPI_COMM_WORLD, &(request3[m]) );

		}

		for (m=0; m<8; m++) {
			MPI_Recv(cmesh + istart[m] , sizeof(MeshBnd)*nmrecv[m], MPI_BYTE, PROC_RANK_RECV_NGB[m], 2, MPI_COMM_WORLD, &(status2[m]));
			//			printf("%d:  nrecv[%d] = %d\n", SUDOM_RANK, m, nrecv[m]);
		}	
		for (m=0; m<8; m++) {
			MPI_Wait(&request3[m], &status3[m]);		
		}


		myfree(cmesh_ex, 106);
		///////////////////////return ;
		ip = 0;
		for (m=0; m<8; m++)
			ip += nmrecv[m];

		ip = nmrecv_total;
		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {

			int x_n = dom_grp_mesh_size[m];

			sm = smesh + dom_grp_mesh_start[m] * y_n * z_n;
			int ii, jj, kk;
			int in, jn, kn;
			int ig, idx;
			for (i=0; i<NPAD; i++) {
				for (j=0; j<y_n; j++) {
					for (k=0; k<z_n; k++) {
						in = dom_grp_mesh_start[m] + i;
						jn = j;
						kn = k;
						if (in < NPAD)
							in += NSIDEMESH;

						cmesh[ ip ].x = in ; 
						cmesh[ ip ].y = jn ; 
						cmesh[ ip ].z = kn ; 
						cmesh[ ip ].data = sm[(i*y_n+j)*z_n+k]; 

						ip ++;


					}
				}
			}
			for (i=dom_grp_mesh_size[m]-NPAD; i<dom_grp_mesh_size[m]; i++) {
				for (j=0; j<y_n; j++) {
					for (k=0; k<z_n; k++) {
						in = dom_grp_mesh_start[m] + i;
						jn = j;
						kn = k;
						if (in >= NSIDEMESH-NPAD)
							in -= NSIDEMESH;

						cmesh[ ip ].x = in ; 
						cmesh[ ip ].y = jn ; 
						cmesh[ ip ].z = kn ; 
						cmesh[ ip ].data = sm[(i*y_n+j)*z_n+k]; 
						ip ++;
					}
				}
			}



		}
		myfree(smesh, 101);
		int nmsend[NDOMinSUDOM] ;
		int imsend[NDOMinSUDOM+1] ;
		int len_bmesh_sud = 0;

		for (n=0; n<NDOMinSUDOM; n++ ) {
			nmsend[n] = 0;
			imsend[n] = 0;
		}
		for (n=NDOM_HEAD; n<NDOMinSUDOM; n++ ){	
			nmsend[n] = (dom_grp_mesh_size[n] + 2*NPAD ) * (y_n + 2*NPAD) *(z_n + 2*NPAD) - dom_grp_mesh_size[n]*y_n *z_n;	
			imsend[n+1] = imsend[n]+ nmsend[n];
		}
		len_bmesh_sud = imsend[NDOMinSUDOM];

		bmesh = (MeshBnd*)mymalloc( sizeof(MeshBnd) * len_bmesh_sud, 108);	

		if (bmesh == NULL) {
			printf(" ERROR : cannot open bmesh  \n");
			fflush(stdout);
			exit(0);
		}

		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {
			ip = imsend[m];
			for (n=0; n<Nmesh_bnd_sud; n++){
				if (cmesh[n].x >= dom_grp_mesh_start[m] - NPAD && cmesh[n].x <dom_grp_mesh_start[m] ) {
					bmesh[ip] = cmesh[n];
					bmesh[ip].x -= dom_grp_mesh_start[m];
					ip ++;
				}
				if (cmesh[n].x >= dom_grp_mesh_end[m] && cmesh[n].x < dom_grp_mesh_end[m]+NPAD ) {
					bmesh[ip] = cmesh[n];
					bmesh[ip].x -= dom_grp_mesh_start[m];
					ip ++;
				}

				if (cmesh[n].x >= dom_grp_mesh_start[m] && cmesh[n].x < dom_grp_mesh_end[m] ) {
					if (cmesh[n].y < 0 || cmesh[n].y >= y_n	|| cmesh[n].z < 0 || cmesh[n].z >= z_n ) 
					{

						bmesh[ip] = cmesh[n];
						bmesh[ip].x -= dom_grp_mesh_start[m];
						ip ++;
					}

				}
			}	

			if (ip - imsend[m] != nmsend[m])	
				printf("[%d] ip = %d, nmsend = %d\n", PROC_RANK, ip - imsend[m], nmsend[m]);	
			// MPI_Isend;
			//
		}

		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {

			MPI_Isend( bmesh + imsend[m], sizeof(MeshBnd)*nmsend[m], MPI_BYTE, dom_grp[m], 400, MPI_COMM_WORLD, &(request[m]) );
		}

	}
	////////////////////////////

	//fflush(stdout);
	//MPI_Barrier(MPI_COMM_WORLD);

	if (COM_DOM) {
		int y_n = NSIDEMESH/NSIDE0SUDOM;
		int z_n = NSIDEMESH/NSIDE1SUDOM;
		int n = dom_grp_rank;
		int nmrecv = (dom_grp_mesh_size[n] + 2*NPAD ) * (y_n + 2*NPAD) *(z_n + 2*NPAD) - dom_grp_mesh_size[n]*y_n *z_n;




		//		printf(" %d > %d %d \n",PROC_RANK, n, nmrecv  );
		bmeshc = (MeshBnd*)mymalloc(sizeof(MeshBnd) * nmrecv, 109);

		if (bmeshc == NULL) {
			printf(" ERROR : cannot open bmeshc  \n");
			fflush(stdout);
			exit(0);
		}





		MPI_Recv( bmeshc, sizeof(MeshBnd) * nmrecv, MPI_BYTE, PROC_RANK_SUDOM, 400, MPI_COMM_WORLD, &sta);
		int i,j,k;
		int ny = y_n + 2*NPAD;
		int nz = z_n + 2*NPAD;
		int idx;
		for (n=0; n<nmrecv; n++) {
			i = bmeshc[n].x + NPAD;			
			j = bmeshc[n].y + NPAD;			
			k = bmeshc[n].z + NPAD;
			idx = (i * ny + j) *nz + k;

			gmesh[idx] = bmeshc[n].data;				
		}
	}



	if (  isSUDOM ) {
		int m;
		for (m=NDOM_HEAD; m<NDOMinSUDOM; m++) {
			MPI_Wait( &(request[m]), &(status[m]) );
		}
	}

	if (  isSUDOM ) {
		myfree(bmesh, 108);
		myfree(cmesh, 107);
	}

//MPI_Barrier(MPI_COMM_WORLD);
	if (COM_DOM) {
		myfree(bmeshc, 109);
	}

	if (COM_DOM) {
		cpuType *mesh = gmesh;
		int n;
		int i,j,k,idx;
		int ii, jj, kk;
		cpuType norm = NSIDEMESH/BOXSIZE;
		cpuType delta = 1.0/norm;

		cpuType invx = 0.5*NSIDEMESH/BOXSIZE;
		isize[0] += 2*NPAD   ;
		isize[1] += 2*NPAD  ;
		isize[2] += 2*NPAD  ;
		meshsize = isize[0] * isize[1] * isize[2];
		cpuType wi, wj, wk;
		cpuType win, wjn, wkn;	
		int nx, ny, nz;

		nx = isize[1]*isize[2];
		ny = isize[2];
		nz = 1;

		local_min[0] -= NPAD;
		local_min[1] -= NPAD;
		local_min[2] -= NPAD;

		cpuType dp[8];
		int idx1, idx2;

		if (mode_alloc) {

			acc[0] = (cpuType*)mymalloc(sizeof(cpuType)*NPART, 14);
			acc[1] = (cpuType*)mymalloc(sizeof(cpuType)*NPART, 15);
			acc[2] = (cpuType*)mymalloc(sizeof(cpuType)*NPART, 16);

			
			for (n=0; n<NPART; n++) {
				acc[0][n] = 0.0;
				acc[1][n] = 0.0;
				acc[2][n] = 0.0;
			}
		}
		cpuType acc_t[3];
		for (n=0; n<NPART; n++) {
			if (part[n].tag != 1) continue; // qwang

			cpuType posf[3];
			posf[0] = (float) (part[n].posi[0] * INT2POS );
			posf[1] = (float) (part[n].posi[1] * INT2POS );
			posf[2] = (float) (part[n].posi[2] * INT2POS );



			i = (int) (posf[0] *norm) ;	
			j = (int) (posf[1] *norm) ;	
			k = (int) (posf[2] *norm) ;	

			acc_t[0] = acc_t[1] =acc_t[2]  = 0.0;

			wi = ( posf[0] - (i+0.5)*delta)*norm;
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

			idx = (i*isize[1] + j)*isize[2] + k;


			if (idx + nx > meshsize)
				printf(" partmesh 1 error : %d %d %d, %lf %lf %lf\n", i, j, k, posf[0], posf[1], posf[2]);

			if (idx - nx < 0)
				printf(" partmesh 2 error : %d %d %d, %lf %lf %lf\n", i, j, k, posf[0], posf[1], posf[2]);


			fflush(stdout);
			cpuType f1 = 4.0/3.0;
			cpuType f2 = 1.0/6.0;
			idx1 = ((i-1)*isize[1] + j)*isize[2] + k;
			idx2 = ((i+1)*isize[1] + j)*isize[2] + k;
			dp[0] =  f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-2)*isize[1] + j)*isize[2] + k;
			idx2 = ((i+2)*isize[1] + j)*isize[2] + k;
			dp[0] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-1)*isize[1] + j)*isize[2] + k;
			idx2 = ((ii+1)*isize[1] + j)*isize[2] + k;
			dp[1] = f1 *invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-2)*isize[1] + j)*isize[2] + k;
			idx2 = ((ii+2)*isize[1] + j)*isize[2] + k;
			dp[1] -= f2 * invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-1)*isize[1] + jj)*isize[2] + k;
			idx2 = ((i+1)*isize[1] + jj)*isize[2] + k;
			dp[2] = f1 * invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-2)*isize[1] + jj)*isize[2] + k;
			idx2 = ((i+2)*isize[1] + jj)*isize[2] + k;
			dp[2] -= f2 * invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-1)*isize[1] + jj)*isize[2] + k;
			idx2 = ((ii+1)*isize[1] + jj)*isize[2] + k;
			dp[3] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-2)*isize[1] + jj)*isize[2] + k;
			idx2 = ((ii+2)*isize[1] + jj)*isize[2] + k;
			dp[3] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-1)*isize[1] + j)*isize[2] + kk;
			idx2 = ((i+1)*isize[1] + j)*isize[2] + kk;
			dp[4] = f1* invx*(mesh[idx2] - mesh[idx1]);


			idx1 = ((i-2)*isize[1] + j)*isize[2] + kk;
			idx2 = ((i+2)*isize[1] + j)*isize[2] + kk;
			dp[4] -= f2 *invx*(mesh[idx2] - mesh[idx1]);
			idx1 = ((ii-1)*isize[1] + j)*isize[2] + kk;
			idx2 = ((ii+1)*isize[1] + j)*isize[2] + kk;
			dp[5] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-2)*isize[1] + j)*isize[2] + kk;
			idx2 = ((ii+2)*isize[1] + j)*isize[2] + kk;
			dp[5] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-1)*isize[1] + jj)*isize[2] + kk;
			idx2 = ((i+1)*isize[1] + jj)*isize[2] + kk;
			dp[6] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((i-2)*isize[1] + jj)*isize[2] + kk;
			idx2 = ((i+2)*isize[1] + jj)*isize[2] + kk;
			dp[6]-= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-1)*isize[1] + jj)*isize[2] + kk;
			idx2 = ((ii+1)*isize[1] + jj)*isize[2] + kk;
			dp[7] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = ((ii-2)*isize[1] + jj)*isize[2] + kk;
			idx2 = ((ii+2)*isize[1] + jj)*isize[2] + kk;
			dp[7]-= f2*invx*(mesh[idx2] - mesh[idx1]);

			acc_t[0]= win*wjn*wkn*dp[0]
				+ wi *wjn*wkn*dp[1]
				+ win*wj *wkn*dp[2]
				+ wi *wj *wkn*dp[3]
				+ win*wjn*wk *dp[4]
				+ wi *wjn*wk *dp[5]
				+ win*wj *wk *dp[6]
				+ wi *wj *wk *dp[7];

			idx1 = (i*isize[1] + j-1)*isize[2] + k;
			idx2 = (i*isize[1] + j+1)*isize[2] + k;
			dp[0] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + j-2)*isize[2] + k;
			idx2 = (i*isize[1] + j+2)*isize[2] + k;
			dp[0] -= f2*invx*(mesh[idx2] - mesh[idx1]);


			idx1 = (ii*isize[1] + j-1)*isize[2] + k;
			idx2 = (ii*isize[1] + j+1)*isize[2] + k;
			dp[1] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j-2)*isize[2] + k;
			idx2 = (ii*isize[1] + j+2)*isize[2] + k;
			dp[1] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + jj-1)*isize[2] + k;
			idx2 = (i*isize[1] + jj+1)*isize[2] + k;
			dp[2] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + jj-2)*isize[2] + k;
			idx2 = (i*isize[1] + jj+2)*isize[2] + k;
			dp[2] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj-1)*isize[2] + k;
			idx2 = (ii*isize[1] + jj+1)*isize[2] + k;
			dp[3] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj-2)*isize[2] + k;
			idx2 = (ii*isize[1] + jj+2)*isize[2] + k;
			dp[3] -= f2*invx*(mesh[idx2] - mesh[idx1]);


			idx1 = (i*isize[1] + j-1)*isize[2] + kk;
			idx2 = (i*isize[1] + j+1)*isize[2] + kk;
			dp[4] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + j-2)*isize[2] + kk;
			idx2 = (i*isize[1] + j+2)*isize[2] + kk;
			dp[4] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j-1)*isize[2] + kk;
			idx2 = (ii*isize[1] + j+1)*isize[2] + kk;
			dp[5] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j-2)*isize[2] + kk;
			idx2 = (ii*isize[1] + j+2)*isize[2] + kk;
			dp[5] -= f2*invx*(mesh[idx2] - mesh[idx1]);


			idx1 = (i*isize[1] + jj-1)*isize[2] + kk;
			idx2 = (i*isize[1] + jj+1)*isize[2] + kk;
			dp[6] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + jj-2)*isize[2] + kk;
			idx2 = (i*isize[1] + jj+2)*isize[2] + kk;
			dp[6] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj-1)*isize[2] + kk;
			idx2 = (ii*isize[1] + jj+1)*isize[2] + kk;
			dp[7] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj-2)*isize[2] + kk;
			idx2 = (ii*isize[1] + jj+2)*isize[2] + kk;
			dp[7] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			acc_t[1]= win*wjn*wkn*dp[0]
				+ wi *wjn*wkn*dp[1]
				+ win*wj *wkn*dp[2]
				+ wi *wj *wkn*dp[3]
				+ win*wjn*wk *dp[4]
				+ wi *wjn*wk *dp[5]
				+ win*wj *wk *dp[6]
				+ wi *wj *wk *dp[7];



			idx1 = (i*isize[1] + j)*isize[2] + k-1;
			idx2 = (i*isize[1] + j)*isize[2] + k+1;
			dp[0] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + j)*isize[2] + k-2;
			idx2 = (i*isize[1] + j)*isize[2] + k+2;
			dp[0] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j)*isize[2] + k-1;
			idx2 = (ii*isize[1] + j)*isize[2] + k+1;
			dp[1] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j)*isize[2] + k-2;
			idx2 = (ii*isize[1] + j)*isize[2] + k+2;
			dp[1] -= f2*invx*(mesh[idx2] - mesh[idx1]);


			idx1 = (i*isize[1] + jj)*isize[2] + k-1;
			idx2 = (i*isize[1] + jj)*isize[2] + k+1;
			dp[2] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + jj)*isize[2] + k-2;
			idx2 = (i*isize[1] + jj)*isize[2] + k+2;
			dp[2] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj)*isize[2] + k-1;
			idx2 = (ii*isize[1] + jj)*isize[2] + k+1;
			dp[3] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj)*isize[2] + k-2;
			idx2 = (ii*isize[1] + jj)*isize[2] + k+2;
			dp[3] -= f2*invx*(mesh[idx2] - mesh[idx1]);


			idx1 = (i*isize[1] + j)*isize[2] + kk-1;
			idx2 = (i*isize[1] + j)*isize[2] + kk+1;
			dp[4] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + j)*isize[2] + kk-2;
			idx2 = (i*isize[1] + j)*isize[2] + kk+2;
			dp[4] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j)*isize[2] + kk-1;
			idx2 = (ii*isize[1] + j)*isize[2] + kk+1;
			dp[5] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + j)*isize[2] + kk-2;
			idx2 = (ii*isize[1] + j)*isize[2] + kk+2;
			dp[5] -= f2*invx*(mesh[idx2] - mesh[idx1]);


			idx1 = (i*isize[1] + jj)*isize[2] + kk-1;
			idx2 = (i*isize[1] + jj)*isize[2] + kk+1;
			dp[6] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (i*isize[1] + jj)*isize[2] + kk-2;
			idx2 = (i*isize[1] + jj)*isize[2] + kk+2;
			dp[6] -= f2*invx*(mesh[idx2] - mesh[idx1]);


			idx1 = (ii*isize[1] + jj)*isize[2] + kk-1;
			idx2 = (ii*isize[1] + jj)*isize[2] + kk+1;
			dp[7] = f1*invx*(mesh[idx2] - mesh[idx1]);

			idx1 = (ii*isize[1] + jj)*isize[2] + kk-2;
			idx2 = (ii*isize[1] + jj)*isize[2] + kk+2;
			dp[7] -= f2*invx*(mesh[idx2] - mesh[idx1]);

			acc_t[2]= win*wjn*wkn*dp[0]
				+ wi *wjn*wkn*dp[1]
				+ win*wj *wkn*dp[2]
				+ wi *wj *wkn*dp[3]
				+ win*wjn*wk *dp[4]
				+ wi *wjn*wk *dp[5]
				+ win*wj *wk *dp[6]
				+ wi *wj *wk *dp[7];


			part[n].vel[0] += acc_t[0] * dkick ;	
			part[n].vel[1] += acc_t[1] * dkick ;	
			part[n].vel[2] += acc_t[2] * dkick ;	

			if (mode_alloc) {
				acc[0][n] = acc_t[0];
				acc[1][n] = acc_t[1];
				acc[2][n] = acc_t[2];
			}
		}




	}
	if ( COM_DOM ) {	
//		myfree(bmeshc, 109);

		myfree(gmesh, 105);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if ( 0 == PROC_RANK ) {
		printf( "\n partmesh complete !\n" );
		fflush(stdout);
	}
	dtime_pm = dtime() - dtime_pm;
}

#endif

