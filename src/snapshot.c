#include "photoNs.h"
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>

#ifndef INTXYZ

void input_snapshot(char fname[], int timestamp, int rank ) {

	char fout[128];
	int n, proc;
	sprintf(fout, "%s/cfg_%d.%d", fname, timestamp, rank);

	FILE *fd = fopen(fout, "r");

	if ( NULL == fd) {
		printf("cannot find file %s\n", fout);
		exit(0);
	}
	cpuType redshift;

	fscanf(fd, "%d", &proc);
	fscanf(fd, "%d", &PROC_SIZE);
	fscanf(fd, "%lld", &NPART_TOTAL);
	fscanf(fd, "%e", &MASSPART);
	fscanf(fd, "%e", &BOXSIZE);
	fscanf(fd, "%e", &Hubble0);
	fscanf(fd, "%e", &OmegaM0);
	fscanf(fd, "%e", &OmegaX0);

	fscanf(fd, "%e", &redshift);
	fscanf(fd, "%d", &NPART);


	//	MAXNPART = fmax((double)NPART, (double)(NPART_TOTAL / PROC_SIZE)) * 1.5;
	//	MAXNPART_BND = MAXNPART * 0.3;
	fscanf(fd, "%d", &MAXNPART);
	fscanf(fd, "%d", &MAXNPART_BND);


	fscanf(fd, "%d", &COM_DOM);

	InitialTime = 1.0/(redshift + 1.0);

	if (COM_DOM == 0) 
		isSUDOM = 1;

	fscanf(fd, "%d", &dom_grp_rank);
	for (n=0; n<NDOMinSUDOM; n++)
		fscanf(fd, "%d\n", &dom_grp_mesh_start[n]);

	for (n=0; n<NDOMinSUDOM; n++)
		fscanf(fd, "%d\n", &dom_grp_mesh_size[n]);

	for (n=0; n<NDOMinSUDOM; n++)
		fscanf(fd, "%d\n", &dom_grp_mesh_end[n]);

	int nsidemesh;
	fscanf(fd, "%d", &nsidemesh);

#ifdef INTXYZ
	int bit_width;
	fscanf(fd, "%d", &bit_width);

	if ( bit_width != BITWIDTH ) {
		printf(" Warning! bit width \n");
		
	}
#endif 
	if (proc != PROC_RANK) {
		printf(" error ! read wrong rank \n");
		exit(0);
	}

	fclose(fd);

	MAXNPART_BND_SUD = MAXNPART_BND * NDOMinSUDOM;

	MAXNPART_DOM = MAXNPART_BND_SUD;
	MAXNPART_DOM_SUD = MAXNPART_BND_SUD;

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


	//printf("[%d] (%e %e) (%e %e) (%e %e).\n", PROC_RANK, BOX_DOM_L[0], BOX_DOM_H[0], BOX_DOM_L[1], BOX_DOM_H[1], BOX_DOM_L[2], BOX_DOM_H[2]);

	if ( COM_DOM) {
		sprintf(fout, "%s/dat_%d.%d", fname, timestamp, rank);

		part = (Body*)mymalloc(sizeof(Body) * MAXNPART, 9);

		if (part == NULL ){
			printf(" error part malloc %d\n", PROC_RANK);
			exit(0);
		}
		part_b = (BodyBnd*)mymalloc(sizeof(BodyBnd)*MAXNPART_BND, 10);

		fd = fopen(fout, "r");
		if ( NULL == fd) {
			printf("cannot open file %s\n", fout);
			exit(0);
		}
#ifdef IOBLK
		int BLK_N = 1024 * 256, m;
		float* vec = (float*)mymalloc(BLK_N * sizeof(float) * 3, 4);
		for (n=0; n<NPART; n+=BLK_N) {
			int np = (NPART-n < BLK_N) ? NPART-n : BLK_N;
			fread(vec, sizeof(float) * 3 * np, 1, fd);
			for (m=0; m<np; m++) {
				part[n+m].pos[0] = (float)vec[3*m];
				part[n+m].pos[1] = (float)vec[3*m+1];
				part[n+m].pos[2] = (float)vec[3*m+2];
			}
		}
		for (n=0; n<NPART; n+=BLK_N) {
			int np = (NPART-n < BLK_N) ? NPART-n : BLK_N;
			fread(vec, sizeof(float) * 3 * np, 1, fd);
			for (m=0; m<np; m++) {
				part[n+m].vel[0] = (cpuType)vec[3*m];
				part[n+m].vel[1] = (cpuType)vec[3*m+1];
				part[n+m].vel[2] = (cpuType)vec[3*m+2];
			}
		}
		myfree(vec, 4);
		//#ifdef HALO_FOF
		unsigned long* lvec = (unsigned long*)mymalloc(BLK_N * sizeof(unsigned long), 4);
		for (n=0; n<NPART; n+=BLK_N) {
			int np = (NPART-n < BLK_N) ? NPART-n : BLK_N;
			fread(lvec, sizeof(unsigned long) * np, 1, fd);
			for (m=0; m<np; m++) {
				part[n+m].id = lvec[m];
			}
		}
		myfree(lvec, 4);
		//#endif


#else

		for (n=0; n<NPART; n++) {

			float x[3];
			fread(x, sizeof(float), 3, fd);
			part[n].pos[0] = (float)x[0];
			part[n].pos[1] = (float)x[1];
			part[n].pos[2] = (float)x[2];

			if ( (part[n].pos[0]) < BOX_DOM_L[0] || (part[n].pos[0] ) > BOX_DOM_H[0])
				printf(" [%d]  er %lf (%lf %lf) \n", PROC_RANK, part[n].pos[0], BOX_DOM_L[0], BOX_DOM_H[0] );

		}
		for (n=0; n<NPART; n++) {
			float x[3];
			fread(x, sizeof(float), 3, fd);
			part[n].vel[0] = (cpuType)x[0];
			part[n].vel[1] = (cpuType)x[1];
			part[n].vel[2] = (cpuType)x[2];
		}

		//#ifdef HALO_FOF
		for (n=0; n<NPART; n++) {
			unsigned long ID;
			fread(&ID, sizeof(unsigned long), 1, fd);
			part[n].id = (unsigned long)ID;
		}
		//#endif

#endif

		fclose(fd);
	}

	if (print_msg())
		printf(" [%d] %d %d (%lf, %lf)\n", PROC_RANK, nsidemesh, NPART, BOX_DOM_L[1], BOX_DOM_H[1]);
}


void output_snapshot_new(cpuType redshift, int timestamp, int rank ) 
{
	char fout[128];
	char ofname[128];
	sprintf(ofname, "%s/snapshot_%d", PathSnapshot,timestamp);
	int n;
	FILE *fd;
	//	if (0 == PROC_RANK) {
	//	}

//	printf(" (%d) zz = %f\n", rank, redshift);

	int np = 0;
	if ( COM_DOM) {
		sprintf(fout, "%s/dat_%d.%d", ofname,timestamp, rank);

		fd = fopen(fout, "w");
		if ( NULL == fd) {
			printf("cannot open file %s\n", fout);
			exit(0);
		}

		int n;
		np = 0;
		for (n=0; n<NPART; n++) {
			if (part[n].tag == 1) {

				float x[3];
				x[0] = part[n].pos[0];
				x[1] = part[n].pos[1];
				x[2] = part[n].pos[2];
				fwrite(x, sizeof(float), 3, fd);
				np++;
			}
		}

		for (n=0; n<NPART; n++) {
			if (part[n].tag == 1) {
				float x[3];
				x[0] = part[n].vel[0];
				x[1] = part[n].vel[1];
				x[2] = part[n].vel[2];
				fwrite(x, sizeof(float), 3, fd);
			}
		}

		for (n=0; n<NPART; n++) {
			if (part[n].tag == 1) {
				unsigned long ID = part[n].id;
				fwrite(&ID, sizeof(unsigned long), 1, fd);
			}
		}

		fclose(fd);

	}

	sprintf(fout, "%s/cfg_%d.%d", ofname, timestamp, rank );

	fd = fopen(fout, "w");
	if ( NULL == fd) {
		printf("cannot open file %s\n", fout);
		exit(0);
	}

	fprintf(fd, "%d\n", PROC_RANK);
	fprintf(fd, "%d\n", PROC_SIZE);
	fprintf(fd, "%lld\n", NPART_TOTAL);
	fprintf(fd, "%e\n", MASSPART);
	fprintf(fd, "%e\n", BOXSIZE);
	fprintf(fd, "%f\n", Hubble0);
	fprintf(fd, "%f\n", OmegaM0);
	fprintf(fd, "%f\n", OmegaX0);

	fprintf(fd, "%e\n", redshift);
	if (COM_DOM)
		fprintf(fd, "%d\n", np);
	else
		fprintf(fd, "0\n");
	fprintf(fd, "%d\n", MAXNPART);
	fprintf(fd, "%d\n", MAXNPART_BND);
	fprintf(fd, "%d\n", COM_DOM);

	fprintf(fd, "%d\n", dom_grp_rank);

	for (n=0; n<NDOMinSUDOM; n++)
		fprintf(fd, "%d\n", dom_grp_mesh_start[n]);

	for (n=0; n<NDOMinSUDOM; n++)
		fprintf(fd, "%d\n", dom_grp_mesh_size[n]);

	for (n=0; n<NDOMinSUDOM; n++)
		fprintf(fd, "%d\n", dom_grp_mesh_end[n]);

	fprintf(fd, "%d\n", NSIDEMESH);
	fclose(fd);





	if (0 == PROC_RANK)
		printf("\n     Snapshot <<< %s/cfg+dat_%d (z=%2.3lf)\n\n", ofname,  timestamp, redshift);
}
#else

void input_snapshot(char fname[], int timestamp, int rank ) {

	char fout[128];
	int n, proc;
	sprintf(fout, "%s/cfg_%d.%d", fname, timestamp, rank);

	FILE *fd = fopen(fout, "r");

	if ( NULL == fd) {
		printf("cannot find file %s\n", fout);
		exit(0);
	}
	cpuType redshift;

	fscanf(fd, "%d", &proc);
	fscanf(fd, "%d", &PROC_SIZE);
	fscanf(fd, "%lld", &NPART_TOTAL);
	fscanf(fd, "%e", &MASSPART);
	fscanf(fd, "%e", &BOXSIZE);
	fscanf(fd, "%e", &Hubble0);
	fscanf(fd, "%e", &OmegaM0);
	fscanf(fd, "%e", &OmegaX0);

	fscanf(fd, "%e", &redshift);
	fscanf(fd, "%d", &NPART);


	//	MAXNPART = fmax((double)NPART, (double)(NPART_TOTAL / PROC_SIZE)) * 1.5;
	//	MAXNPART_BND = MAXNPART * 0.3;
	fscanf(fd, "%d", &MAXNPART);
	fscanf(fd, "%d", &MAXNPART_BND);


	fscanf(fd, "%d", &COM_DOM);

	InitialTime = 1.0/(redshift + 1.0);

	if (COM_DOM == 0) 
		isSUDOM = 1;

	fscanf(fd, "%d", &dom_grp_rank);
	for (n=0; n<NDOMinSUDOM; n++)
		fscanf(fd, "%d\n", &dom_grp_mesh_start[n]);

	for (n=0; n<NDOMinSUDOM; n++)
		fscanf(fd, "%d\n", &dom_grp_mesh_size[n]);

	for (n=0; n<NDOMinSUDOM; n++)
		fscanf(fd, "%d\n", &dom_grp_mesh_end[n]);

	int nsidemesh;
	fscanf(fd, "%d", &nsidemesh);


	int bit_width;
	fscanf(fd, "%d", &bit_width);


	if (proc != PROC_RANK) {
		printf(" error ! read wrong rank \n");
		exit(0);
	}


	fclose(fd);

	MAXNPART_BND_SUD = MAXNPART_BND * NDOMinSUDOM;

	MAXNPART_DOM = MAXNPART_BND_SUD;
	MAXNPART_DOM_SUD = MAXNPART_BND_SUD;

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


//	printf("[%d] (%e %e) (%e %e) (%e %e).\n", PROC_RANK, BOX_DOM_L[0], BOX_DOM_H[0], BOX_DOM_L[1], BOX_DOM_H[1], BOX_DOM_L[2], BOX_DOM_H[2]);
//return;

	if ( COM_DOM) {
		sprintf(fout, "%s/dat_%d.%d", fname, timestamp, rank);

		part = (Body*)mymalloc(sizeof(Body) * MAXNPART, 9);

		if (part == NULL ){
			printf(" error part malloc %d\n", PROC_RANK);
			exit(0);
		}
		part_b = (BodyBnd*)mymalloc(sizeof(BodyBnd)*MAXNPART_BND, 10);

		fd = fopen(fout, "r");
		if ( NULL == fd) {
			printf("cannot open file %s\n", fout);
			exit(0);
		}
#ifdef IOBLK
		int BLK_N = 1024 * 256, m;

		int* veci = (int*)mymalloc(BLK_N * sizeof(int) * 3, 4);
		for (n=0; n<NPART; n+=BLK_N) {
			int np = (NPART-n < BLK_N) ? NPART-n : BLK_N;
			fread(veci, sizeof(int) * 3 * np, 1, fd);
			for (m=0; m<np; m++) {
				part[n+m].posi[0] = veci[3*m];
				part[n+m].posi[1] = veci[3*m+1];
				part[n+m].posi[2] = veci[3*m+2];
			}
		}
		myfree(veci, 4);

		float* vec = (float*)mymalloc(BLK_N * sizeof(float) * 3, 4);
		for (n=0; n<NPART; n+=BLK_N) {
			int np = (NPART-n < BLK_N) ? NPART-n : BLK_N;
			fread(vec, sizeof(float) * 3 * np, 1, fd);
			for (m=0; m<np; m++) {
				part[n+m].vel[0] = (cpuType)vec[3*m];
				part[n+m].vel[1] = (cpuType)vec[3*m+1];
				part[n+m].vel[2] = (cpuType)vec[3*m+2];
			}
		}
		myfree(vec, 4);

		//#ifdef HALO_FOF
		unsigned long* lvec = (unsigned long*)mymalloc(BLK_N * sizeof(unsigned long), 4);
		for (n=0; n<NPART; n+=BLK_N) {
			int np = (NPART-n < BLK_N) ? NPART-n : BLK_N;
			fread(lvec, sizeof(unsigned long) * np, 1, fd);
			for (m=0; m<np; m++) {
				part[n+m].id = lvec[m];
			}
		}
		myfree(lvec, 4);
		//#endif


#else

		for (n=0; n<NPART; n++) {
			int x[3];
			fread(x, sizeof(int), 3, fd);
			part[n].posi[0] = x[0];
			part[n].posi[1] = x[1];
			part[n].posi[2] = x[2];
		//	if ( (float)(part[n].posi[0]*INT2POS) < BOX_DOM_L[0] ||(float) (part[n].posi[0] * INT2POS ) > BOX_DOM_H[0])
		//		printf(" [%d]  eri %f, %d (%f %f) \n", PROC_RANK,(float)( part[n].posi[0] *INT2POS), part[n].posi[0] , BOX_DOM_L[0], BOX_DOM_H[0] );


		}
		for (n=0; n<NPART; n++) {
			float x[3];
			fread(x, sizeof(float), 3, fd);
			part[n].vel[0] = (cpuType)x[0];
			part[n].vel[1] = (cpuType)x[1];
			part[n].vel[2] = (cpuType)x[2];
		}

		//#ifdef HALO_FOF
		for (n=0; n<NPART; n++) {
			unsigned long ID;
			fread(&ID, sizeof(unsigned long), 1, fd);
			part[n].id = (unsigned long)ID;
		}
		//#endif

#endif

		fclose(fd);
	}

	if (print_msg())
		printf(" [%d] %d %d (%lf, %lf)\n", PROC_RANK, nsidemesh, NPART, BOX_DOM_L[1], BOX_DOM_H[1]);
}

/*
void input_snapshot_compress(char fname[], int timestamp, int rank ) {

	char fout[128];
	int n, proc;
	sprintf(fout, "%s/cfg_%d.%d", fname, timestamp, rank);

	FILE *fd = fopen(fout, "r");

	if ( NULL == fd) {
		printf("cannot find file %s\n", fout);
		exit(0);
	}
	cpuType redshift;

	fscanf(fd, "%d", &proc);
	fscanf(fd, "%d", &PROC_SIZE);
	fscanf(fd, "%lld", &NPART_TOTAL);
	fscanf(fd, "%e", &MASSPART);
	fscanf(fd, "%e", &BOXSIZE);
	fscanf(fd, "%e", &Hubble0);
	fscanf(fd, "%e", &OmegaM0);
	fscanf(fd, "%e", &OmegaX0);

	fscanf(fd, "%e", &redshift);
	fscanf(fd, "%d", &NPART);
	fscanf(fd, "%d", &MAXNPART);
	fscanf(fd, "%d", &MAXNPART_BND);

	fscanf(fd, "%d", &COM_DOM);

	InitialTime = 1.0/(redshift + 1.0);

	if (COM_DOM == 0) 
		isSUDOM = 1;

	fscanf(fd, "%d", &dom_grp_rank);
	for (n=0; n<NDOMinSUDOM; n++)
		fscanf(fd, "%d\n", &dom_grp_mesh_start[n]);

	for (n=0; n<NDOMinSUDOM; n++)
		fscanf(fd, "%d\n", &dom_grp_mesh_size[n]);

	for (n=0; n<NDOMinSUDOM; n++)
		fscanf(fd, "%d\n", &dom_grp_mesh_end[n]);

	int nsidemesh;
	fscanf(fd, "%d", &nsidemesh);
	if (proc != PROC_RANK) {
		printf(" error ! read wrong rank \n");
		exit(0);
	}

	fclose(fd);

	MAXNPART_BND_SUD = MAXNPART_BND * NDOMinSUDOM;

	MAXNPART_DOM = MAXNPART_BND_SUD;
	MAXNPART_DOM_SUD = MAXNPART_BND_SUD;

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


	//printf("[%d] (%e %e) (%e %e) (%e %e).\n", PROC_RANK, BOX_DOM_L[0], BOX_DOM_H[0], BOX_DOM_L[1], BOX_DOM_H[1], BOX_DOM_L[2], BOX_DOM_H[2]);

	if ( COM_DOM) {
		sprintf(fout, "%s/pos_c_%d.%d", fname, timestamp, rank);

		part = (Body*)mymalloc(sizeof(Body) * MAXNPART, 9);

		if (part == NULL ){
			printf(" error part malloc %d\n", PROC_RANK);
			exit(0);
		}
		part_b = (BodyBnd*)mymalloc(sizeof(BodyBnd)*MAXNPART_BND, 10);

		int npart = cmprs_fread(part, NULL, fout);
		if (npart != NPART) {
			if (PROC_RANK == 1)
				printf("[%d] check npart %d != NPART %d\n", PROC_RANK, npart, NPART);
			exit(0);
		}
	}

	if (print_msg())
		printf(" [%d] %d %d (%lf, %lf)\n", PROC_RANK, nsidemesh, NPART, BOX_DOM_L[1], BOX_DOM_H[1]);
}

void output_snapshot_compress(cpuType redshift, int timestamp, int rank ){

	char fout[128];
	char ofname[128];
	sprintf(ofname, "%s/part_snap_%d", PathSnapshot,timestamp);
	int n;
	FILE *fd;
	sprintf(fout, "%s/cfg_%d.%d", ofname, timestamp, rank );

	fd = fopen(fout, "w");
	if ( NULL == fd) {
		printf("cannot open file %s\n", fout);
		exit(0);
	}



	fprintf(fd, "%d\n", PROC_RANK);
	fprintf(fd, "%d\n", PROC_SIZE);
	fprintf(fd, "%lld\n", NPART_TOTAL);
	fprintf(fd, "%e\n", MASSPART);
	fprintf(fd, "%e\n", BOXSIZE);
	fprintf(fd, "%e\n", Hubble0);
	fprintf(fd, "%e\n", OmegaM0);
	fprintf(fd, "%e\n", OmegaX0);

	fprintf(fd, "%e\n", redshift);
	if (COM_DOM)
		fprintf(fd, "%d\n", NPART);
	else
		fprintf(fd, "0\n");
	fprintf(fd, "%d\n", MAXNPART);
	fprintf(fd, "%d\n", MAXNPART_BND);
	fprintf(fd, "%d\n", COM_DOM);

	fprintf(fd, "%d\n", dom_grp_rank);

	for (n=0; n<NDOMinSUDOM; n++)
		fprintf(fd, "%d\n", dom_grp_mesh_start[n]);

	for (n=0; n<NDOMinSUDOM; n++)
		fprintf(fd, "%d\n", dom_grp_mesh_size[n]);

	for (n=0; n<NDOMinSUDOM; n++)
		fprintf(fd, "%d\n", dom_grp_mesh_end[n]);

	fprintf(fd, "%d\n", NSIDEMESH);
	fclose(fd);

	if ( COM_DOM) {
		sprintf(fout, "%s/pos_c_%d.%d", ofname,timestamp, rank);
		compress_fwrite(part, NPART, fout);
	}
	if (0 == PROC_RANK)
		printf("\n     Compress <<< %s/cfg+pos_c_%d (z=%2.3lf)\n\n", ofname,  timestamp, redshift);
}
*/
void output_snapshot_new(cpuType redshift, int timestamp, int rank ) 
{
	char fout[128];
	char ofname[128];
	sprintf(ofname, "%s/snapshot_%d", PathSnapshot,timestamp);
	int n;
	FILE *fd;
	//	if (0 == PROC_RANK) {
	//	}

//	printf(" (%d) zz = %f\n", rank, redshift);

	int np = 0;
	if ( COM_DOM) {
		sprintf(fout, "%s/dat_%d.%d", ofname,timestamp, rank);

		fd = fopen(fout, "w");
		if ( NULL == fd) {
			printf("cannot open file %s\n", fout);
			exit(0);
		}

		int n;
		np = 0;
		for (n=0; n<NPART; n++) {
			if (part[n].tag == 1) {

				int x[3];
				x[0] = part[n].posi[0];
				x[1] = part[n].posi[1];
				x[2] = part[n].posi[2];
				fwrite(x, sizeof(int), 3, fd);

				np++;
			}
		}

		for (n=0; n<NPART; n++) {
			if (part[n].tag == 1) {
				float x[3];
				x[0] = part[n].vel[0];
				x[1] = part[n].vel[1];
				x[2] = part[n].vel[2];
				fwrite(x, sizeof(float), 3, fd);
			}
		}

		for (n=0; n<NPART; n++) {
			if (part[n].tag == 1) {
				unsigned long ID = part[n].id;
				fwrite(&ID, sizeof(unsigned long), 1, fd);
			}
		}

		fclose(fd);

	}

	sprintf(fout, "%s/cfg_%d.%d", ofname, timestamp, rank );

	fd = fopen(fout, "w");
	if ( NULL == fd) {
		printf("cannot open file %s\n", fout);
		exit(0);
	}



	fprintf(fd, "%d\n", PROC_RANK);
	fprintf(fd, "%d\n", PROC_SIZE);
	fprintf(fd, "%lld\n", NPART_TOTAL);
	fprintf(fd, "%e\n", MASSPART);
	fprintf(fd, "%e\n", BOXSIZE);
	fprintf(fd, "%f\n", Hubble0);
	fprintf(fd, "%f\n", OmegaM0);
	fprintf(fd, "%f\n", OmegaX0);

	fprintf(fd, "%f\n", redshift);
	if (COM_DOM)
		fprintf(fd, "%d\n", np);
	else
		fprintf(fd, "0\n");
	fprintf(fd, "%d\n", MAXNPART);
	fprintf(fd, "%d\n", MAXNPART_BND);
	fprintf(fd, "%d\n", COM_DOM);

	fprintf(fd, "%d\n", dom_grp_rank);

	for (n=0; n<NDOMinSUDOM; n++)
		fprintf(fd, "%d\n", dom_grp_mesh_start[n]);

	for (n=0; n<NDOMinSUDOM; n++)
		fprintf(fd, "%d\n", dom_grp_mesh_size[n]);

	for (n=0; n<NDOMinSUDOM; n++)
		fprintf(fd, "%d\n", dom_grp_mesh_end[n]);

	fprintf(fd, "%d\n", NSIDEMESH);

	fprintf(fd, "%d\n", BITWIDTH );

	fclose(fd);

	if (0 == PROC_RANK)
		printf("\n     Snapshot <<< %s/cfg+dat_%d (z=%2.3lf)\n\n", ofname,  timestamp, redshift);
}

#endif

