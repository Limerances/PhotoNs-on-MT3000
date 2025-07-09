#include "photoNs.h"
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include "snapshot.h"

#ifdef GADGETFILE


typedef unsigned int uint4byte;
typedef int int4byte;
////
//typedef struct
//{
//	uint4byte npart[6];      /*!< npart[1] gives the number of particles in the present file, other particle types are ignored */
//	cpuType mass[6];          /*!< mass[1] gives the particle mass */
//	double time;             /*!< time (=cosmological scale factor) of snapshot */
//	cpuType redshift;         /*!< redshift of snapshot */
//	int4byte flag_sfr;       /*!< flags whether star formation is used (not available in L-Gadget2) */
//	int4byte flag_feedback;  /*!< flags whether feedback from star formation is included */
//	uint4byte npartTotal[6]; /*!< npart[1] gives the total number of particles in the run. If this number exceeds 2^32, the npartTotal[2] stores
//				   the result of a division of the particle number by 2^32, while npartTotal[1] holds the remainder. */
//	int4byte flag_cooling;   /*!< flags whether radiative cooling is included */
//	int4byte num_files;      /*!< determines the number of files that are used for a snapshot */
//	cpuType BoxSize;          /*!< Simulation box size (in code units) */
//	cpuType Omega0;           /*!< matter density */
//	cpuType OmegaLambda;      /*!< vacuum energy density */
//	cpuType HubbleParam;      /*!< little 'h' */
//	int4byte flag_stellarage;     /*!< flags whether the age of newly formed stars is recorded and saved */
//	int4byte flag_metals;         /*!< flags whether metal enrichment is included */
//	uint4byte npartTotalHighWord[6];   /*!< High word of the total number of particles of each type */
//	char fill[64];
//} GadgetHeader;
//
//
typedef struct
{
	uint4byte npart[6];      /*!< npart[1] gives the number of particles in the present file, other particle types are ignored */
	double mass[6];          /*!< mass[1] gives the particle mass */
	double time;             /*!< time (=cosmological scale factor) of snapshot */
	double redshift;         /*!< redshift of snapshot */
	int4byte flag_sfr;       /*!< flags whether star formation is used (not available in L-Gadget2) */
	int4byte flag_feedback;  /*!< flags whether feedback from star formation is included */
	uint4byte npartTotal[6]; /*!< npart[1] gives the total number of particles in the run. If this number exceeds 2^32, the npartTotal[2] stores
				   the result of a division of the particle number by 2^32, while npartTotal[1] holds the remainder. */
	int4byte flag_cooling;   /*!< flags whether radiative cooling is included */
	int4byte num_files;      /*!< determines the number of files that are used for a snapshot */
	double BoxSize;          /*!< Simulation box size (in code units) */
	double Omega0;           /*!< matter density */
	double OmegaLambda;      /*!< vacuum energy density */
	double HubbleParam;      /*!< little 'h' */
	int4byte flag_stellarage;     /*!< flags whether the age of newly formed stars is recorded and saved */
	int4byte flag_metals;         /*!< flags whether metal enrichment is included */
	uint4byte npartTotalHighWord[6];   /*!< High word of the total number of particles of each type */
	char fill[64];
} GadgetHeader;


static cpuType sp_unit = 1.0;

static cpuType partmass[6];
static cpuType Omega0;
static cpuType OmegaLambda;
static cpuType HubbleParam;


void npart_infile(const char filename[], int ith, int np[])
{
	int dummy, byte;
	FILE* fsnap;
	GadgetHeader head;

	char fname[200];
	sprintf(fname, "%s.%d",filename, ith);


	if ( !(fsnap=fopen(fname,"r")) )
	{
		printf("[read_GadgetHeader] cannot open file %s \n", filename);
		exit(0);
	}

	byte = fread(&dummy, sizeof(dummy),        1, fsnap);
	byte = fread(&head,  sizeof(GadgetHeader), 1, fsnap);
	byte = fread(&dummy, sizeof(dummy),        1, fsnap);

	BOXSIZE = head.BoxSize;

	np[0] = head.npart[0];
	np[1] = head.npart[1];
	np[2] = head.npart[2];
	np[3] = head.npart[3];
	np[4] = head.npart[4];
	np[5] = head.npart[5];

	fclose(fsnap);

}

void read_GadgetHeader(const char filename[], int *npart_file)
{
	int dummy, byte;
	FILE* fsnap;
	GadgetHeader head;

	if ( !(fsnap=fopen(filename,"r")) )
	{
		printf("[read_GadgetHeader] cannot open file %s \n", filename);
		exit(0);
	}

	byte = fread(&dummy, sizeof(dummy),        1, fsnap);
	byte = fread(&head,  sizeof(GadgetHeader), 1, fsnap);
	byte = fread(&dummy, sizeof(dummy),        1, fsnap);

	BOXSIZE = head.BoxSize;

	NPART_TOTAL = head.npartTotal[0];
	NPART_TOTAL+= head.npartTotal[1];
	NPART_TOTAL+= head.npartTotal[2];
	NPART_TOTAL+= head.npartTotal[3];
	NPART_TOTAL+= head.npartTotal[4];
	NPART_TOTAL+= head.npartTotal[5];

	MASSPART  = head.mass[1];

	partmass[0] = head.mass[0];
	partmass[1] = head.mass[1];
	partmass[2] = head.mass[2];
	partmass[3] = head.mass[3];
	partmass[4] = head.mass[4];
	partmass[5] = head.mass[5];
//printf(" 1  box  = %f\n", BOXSIZE);
	//#ifdef MEGAPARSEC
	sp_unit = 1000.0;
//	if ( fabs(BOXSIZE - head.BoxSize * sp_unit ) > 1.0e-5 ) {
////		printf("Error: setup box size \n\n\n");
//		exit(0);
//	}

	BOXSIZE =  (cpuType) (head.BoxSize);
	OmegaM0 =  (cpuType) (head.Omega0);
	OmegaX0 =  (cpuType) (head.OmegaLambda);
	Hubble0 =  (cpuType) (head.HubbleParam);
	cpuType	InitialRedshift = (cpuType) (head.redshift);
	InitialTime = 1.0/(1.0+InitialRedshift);

	*npart_file = head.npart[1];

	//#endif
	BOXSIZE *= sp_unit;
	if ( 0 &&  0==PROC_RANK) {
		printf("\n PARAMETERS from IC snapshots : \n");
		printf("            + Omega Matter =  %f\n", OmegaM0 );
		printf("            + Omega Lambda =  %f\n", OmegaX0 );
		printf("            + Hubble Const =  %f\n", Hubble0 );
		printf("            + Box Size = %f \n", BOXSIZE);   
		printf("            + Mass Part =  %f\n", MASSPART );
		printf("            + Num Particle =  %lld\n", NPART_TOTAL);   
		printf("            + IC Red-Shift =  %f\n", InitialRedshift);   
		printf("\n");   
		fflush(stdout);
	}
	dummy = byte;




	fclose(fsnap);

}

void read_Particle_Gadget2_domain(char filename[],int *part_start)
{
	int n, m, c, p;
	int n_count;
	int dummy, byte, ID;
	unsigned long disp, num_part;

	int DATA_SIZE = sizeof(float);
	int ID_SIZE = sizeof(int);

	float data[3];


	FILE* fsnap;
	GadgetHeader head;
	if ( !(fsnap=fopen(filename,"r")) ) {
		printf("[read_Particle] cannot open mfile %s \n", filename);
		exit(0);
	}


	byte = 0;
	byte += fread(&dummy, sizeof(dummy),        1, fsnap);
	byte += fread(&head,  sizeof(GadgetHeader), 1, fsnap);
	byte += fread(&dummy, sizeof(dummy),        1, fsnap);

//printf(" 3  box  = %f\n", BOXSIZE);

	///	printf(" reading particle box=%lf z=%lf\n", head.BoxSize * sp_unit, head.redshift);
	sp_unit = 1000.0;

	if ( COM_DOM  ) {
		int* pidx = (int*)malloc(sizeof(int)*MAXNPART);
		int ip = *part_start;

				//printf(" [%d] %f %f %f : %f %f %f\n", PROC_RANK, BOX_DOM_L[0], BOX_DOM_L[1], BOX_DOM_L[2], BOX_DOM_H[0], BOX_DOM_H[1], BOX_DOM_H[2]);
		cpuType px[3];

		byte += fread(&dummy, sizeof(int), 1, fsnap);   
		for (c=0, p=0, m=0; m<6; m++) {
			for (n=0; n<head.npart[m]; n++) {



#ifndef INTXYZ
				byte += fread(data, DATA_SIZE, 3, fsnap);

				px[0] = (cpuType)data[0] * sp_unit;
				px[1] = (cpuType)data[1] * sp_unit;
				px[2] = (cpuType)data[2] * sp_unit;

				if (px[0] >= BOXSIZE ) px[0] -= BOXSIZE;
				if (px[1] >= BOXSIZE ) px[1] -= BOXSIZE;
				if (px[2] >= BOXSIZE ) px[2] -= BOXSIZE;

				if (px[0] < 0.0 ) px[0] += BOXSIZE;
				if (px[1] < 0.0 ) px[1] += BOXSIZE;
				if (px[2] < 0.0 ) px[2] += BOXSIZE;

				if ( BOX_DOM_L[0] <= px[0] && px[0] < BOX_DOM_H[0]
					&& BOX_DOM_L[1] <= px[1] && px[1] < BOX_DOM_H[1]
					&& BOX_DOM_L[2] <= px[2] && px[2] < BOX_DOM_H[2] )
				{


					part[ip+p].pos[0] = px[0];
					part[ip+p].pos[1] = px[1];
					part[ip+p].pos[2] = px[2];

					part[ip+p].tag = 1;
					part[ip+p].act = 0;
					pidx[c] = 1;
					p++;

				}
				else 
					pidx[c] = 0;

#else
				byte += fread(data, DATA_SIZE, 3, fsnap);

				px[0] = (cpuType)data[0] * sp_unit;
				px[1] = (cpuType)data[1] * sp_unit;
				px[2] = (cpuType)data[2] * sp_unit;

				if (px[0] >= BOXSIZE ) px[0] -= BOXSIZE;
				if (px[1] >= BOXSIZE ) px[1] -= BOXSIZE;
				if (px[2] >= BOXSIZE ) px[2] -= BOXSIZE;

				if (px[0] < 0.0 ) px[0] += BOXSIZE;
				if (px[1] < 0.0 ) px[1] += BOXSIZE;
				if (px[2] < 0.0 ) px[2] += BOXSIZE;


				if ( BOX_DOM_L[0] <= px[0] && px[0] < BOX_DOM_H[0]
				&& BOX_DOM_L[1] <= px[1] && px[1] < BOX_DOM_H[1]
				&& BOX_DOM_L[2] <= px[2] && px[2] < BOX_DOM_H[2] )
				{


	//				part[ip+p].pos[0] = px[0];
	//				part[ip+p].pos[1] = px[1];
	//				part[ip+p].pos[2] = px[2];
					part[ip+p].posi[0] =(int) (px[0]*POS2INT);
					part[ip+p].posi[1] =(int) (px[1]*POS2INT);
					part[ip+p].posi[2] =(int) (px[2]*POS2INT);

					part[ip+p].tag = 1;
					part[ip+p].act = 0;
					pidx[c] = 1;
					p++;

				}
				else 
					pidx[c] = 0;



#endif


				c++;
			}
		}
		byte += fread(&dummy, sizeof(int), 1, fsnap);

		cpuType gdt2unit = pow( 1.0/(1.0+head.redshift), 1.5); //for gadget2 format
		byte += fread(&dummy, sizeof(int), 1, fsnap);
		for (c=0, p=0, m=0; m<6; m++) {
			for (n=0; n<head.npart[m]; n++) {
				byte += fread(data, DATA_SIZE, 3, fsnap);

				if ( 1 ==  pidx[c]) {
					part[ip+p].vel[0] = (cpuType)data[0] * gdt2unit ;
					part[ip+p].vel[1] = (cpuType)data[1] * gdt2unit ;
					part[ip+p].vel[2] = (cpuType)data[2] * gdt2unit ;
					p++;
				}
				c++;
			}
		}
		byte += fread(&dummy, sizeof(int), 1, fsnap);

		byte += fread(&dummy, sizeof(int), 1, fsnap);
		for (c=0, p=0, m=0; m<6; m++) {
			for (n=0; n<head.npart[m]; n++) {
				byte += fread(&ID, ID_SIZE, 1, fsnap);

				if ( 1 ==  pidx[c]) {
					part[ip+p].id = (long)ID;
					p++;
				}
				c++;
			}
		}
		byte += fread(&dummy, sizeof(int), 1, fsnap);

		free(pidx);
		*part_start += p;
	}

	fclose(fsnap);
}


void read_Particle_Gadget2_mfile(char filename[], int n_start, int n_end, int ip)
{
	int n, m, c, p;
	int n_count;
	int dummy, byte, ID;
	unsigned long disp, num_part;

	int DATA_SIZE = sizeof(float);
	int ID_SIZE = sizeof(int);

	float data[3];

	n_count = n_start - n_end;

	FILE* fsnap;
	GadgetHeader head;
	if ( !(fsnap=fopen(filename,"r")) ) {
		printf("[read_Particle] cannot open mfile %s \n", filename);
		exit(0);
	}

	//	printf(" < <  reading `%s'\n", filename);

	byte = 0;
	byte += fread(&dummy, sizeof(dummy),        1, fsnap);
	byte += fread(&head,  sizeof(GadgetHeader), 1, fsnap);
	byte += fread(&dummy, sizeof(dummy),        1, fsnap);


	Omega0 = head.Omega0;
	OmegaLambda = head.OmegaLambda;
	HubbleParam = head.HubbleParam;

	//    printf(" npart = %d %d %d \n", head.npart[0], head.npart[1], head.npart[2]);
	//    printf(" mass = %lf %lf %lf \n", head.mass[0], head.mass[1], head.mass[2]);
	//    printf(" time = %lf , redshift = %lf\n", head.time, head.redshift);
	//    printf(" OmegaM = %lf, %lf, Hubble= %lf \n", head.Omega0, head.OmegaLambda, head.HubbleParam);


	byte += fread(&dummy, sizeof(int), 1, fsnap);   
	for (c=0, p=0, m=0; m<6; m++) {
		for (n=0; n<head.npart[m]; n++) {
			byte += fread(data, DATA_SIZE, 3, fsnap);

			if (c>=n_start && c<n_end) {

#ifndef INTXYZ
				part[ip+p].pos[0] = (cpuType)data[0] * sp_unit;
				part[ip+p].pos[1] = (cpuType)data[1] * sp_unit;
				part[ip+p].pos[2] = (cpuType)data[2] * sp_unit;
#else

				part[ip+p].posi[0] = (int)data[0] * sp_unit * POS2INT;
				part[ip+p].posi[1] = (int)data[1] * sp_unit * POS2INT;
				part[ip+p].posi[2] = (int)data[2] * sp_unit * POS2INT;


#endif
				p++;
			}
			c++;
		}
	}
	byte += fread(&dummy, sizeof(int), 1, fsnap);

	cpuType gdt2unit = pow( 1.0/(1.0+head.redshift), 1.5); //for gadget2 format
	byte += fread(&dummy, sizeof(int), 1, fsnap);
	for (c=0, p=0, m=0; m<6; m++) {
		for (n=0; n<head.npart[m]; n++) {
			byte += fread(data, DATA_SIZE, 3, fsnap);

			if (c>=n_start && c<n_end) {
				part[ip+p].vel[0] = (cpuType)data[0] * gdt2unit;
				part[ip+p].vel[1] = (cpuType)data[1] * gdt2unit;
				part[ip+p].vel[2] = (cpuType)data[2] * gdt2unit;
				p++;
			}
			c++;
		}
	}
	byte += fread(&dummy, sizeof(int), 1, fsnap);

#ifdef PARTIDON



	byte += fread(&dummy, sizeof(int), 1, fsnap);
	for (c=0, p=0, m=0; m<6; m++) {
		for (n=0; n<head.npart[m]; n++) {
			byte += fread(&ID, ID_SIZE, 1, fsnap);

			if (c>=n_start && c<n_end) {
				part[ip+p].id = (long)ID;
				p++;
			}
			c++;
		}
	}
	byte += fread(&dummy, sizeof(int), 1, fsnap);
#endif

#ifdef LABEL_BODY
	for (p=0; p<NPART; p++) {
		part[p].id = (long)( p + n_start);
		//	printf(" %d %d %ld\n", PROC_RANK, p, part[p].id);
	}

#endif



	fclose(fsnap);
}

void read_Particle_Gadget2_s(char filename[], int n_start, int n_count)
{
	int n, m, c, p;
	int n_end;
	int dummy, byte, ID;
	unsigned long disp, num_part;

	int DATA_SIZE = sizeof(float);
	int ID_SIZE = sizeof(int);

	float data[3];

	n_end = n_start + n_count;

	FILE* fsnap;
	GadgetHeader head;
	if ( !(fsnap=fopen(filename,"r")) ) {
		printf("[read_Particle_Gadget2] cannot open ic %s \n", filename);
		exit(0);
	}

	printf(" > reading `%s'\n", filename);

	byte = 0;
	byte += fread(&dummy, sizeof(dummy),        1, fsnap);
	byte += fread(&head,  sizeof(GadgetHeader), 1, fsnap);
	byte += fread(&dummy, sizeof(dummy),        1, fsnap);

	OmegaM0 = head.Omega0;
	OmegaX0 = head.OmegaLambda;
	Hubble0 = head.HubbleParam;

	printf(" [%d] box_x %lf %lf\n", PROC_RANK, BOX_DOM_L[0], BOX_DOM_H[0]);
	printf(" [%d] box_y %lf %lf\n", PROC_RANK, BOX_DOM_L[1], BOX_DOM_H[1]);
	printf(" [%d] box_z %lf %lf\n", PROC_RANK, BOX_DOM_L[2], BOX_DOM_H[2]);

	MASSPART = head.mass[1];
	sp_unit = 1000.0;
	int *label = mymalloc(sizeof(int)*n_count, 17);	
	byte += fread(&dummy, sizeof(int), 1, fsnap);   
	for (c=0, p=0, m=0; m<6; m++) {
		for (n=0; n<head.npart[m]; n++) {
			byte += fread(data, DATA_SIZE, 3, fsnap);

			//	if (c>=n_start && c<n_end) {

			cpuType px = (cpuType)data[0] * sp_unit;
			cpuType py = (cpuType)data[1] * sp_unit;
			cpuType pz = (cpuType)data[2] * sp_unit;

			px += 1953.125;	
			py += 1953.125;	
			pz += 1953.125;	

			if (px < 0.0) px+=1000000.0;	
			if (py < 0.0) py+=1000000.0;	
			if (pz < 0.0) pz+=1000000.0;	

			if (px >= 1000000.0) px-=1000000.0;	
			if (py >= 1000000.0) py-=1000000.0;	
			if (pz >= 1000000.0) pz-=1000000.0;	

			if( BOX_DOM_L[0] <= px && px < BOX_DOM_H[0] 
					&&  BOX_DOM_L[1] <= py && py < BOX_DOM_H[1]
					&&  BOX_DOM_L[2] <= pz && pz < BOX_DOM_H[2] ) 
			{

#ifndef INTXYZ 
				part[p].pos[0] = px;
				part[p].pos[1] = py;
				part[p].pos[2] = pz;
#else
				part[p].posi[0] = (int)( px * POS2INT);
				part[p].posi[1] = (int)( py * POS2INT);
				part[p].posi[2] = (int)( pz * POS2INT);
#endif
				part[p].tag = 1;
				part[p].act = 0;
				label[c] =1;
				p++;
			}
			else
				label[c] = 0;

			//	}
			c++;
		}
	}
	NPART = p;
	//printf(" readin %d\n", NPART );
	byte += fread(&dummy, sizeof(int), 1, fsnap);

	cpuType gdt2unit = pow( 1.0/(1.0+head.redshift), 1.5); //for gadget2 format

	//	gdt2unit = 1.0;

	byte += fread(&dummy, sizeof(int), 1, fsnap);
	for (c=0, p=0, m=0; m<6; m++) {
		for (n=0; n<head.npart[m]; n++) {
			byte += fread(data, DATA_SIZE, 3, fsnap);

			if (label[c] == 1) {
				//		if (c>=n_start && c<n_end) {
				part[p].vel[0] = (cpuType)data[0] * gdt2unit;
				part[p].vel[1] = (cpuType)data[1] * gdt2unit;
				part[p].vel[2] = (cpuType)data[2] * gdt2unit;
				p++;
				//		}
			}
			c++;
		}
	}
	if (p != NPART)  {

		printf("error  \n");
		exit(0);
	}
	//	for (p=0; p< NPART; p++) {
	//		if (part[p].pos[0] < 2000.0 && part[p].pos[1] < 2000.0)
	//			printf(" %lf %lf %lf %lf %lf %lf\n", part[p].pos[0], part[p].pos[1], part[p].pos[2], part[p].vel[0], part[p].vel[1], part[p].vel[2]);

	//	}

	byte += fread(&dummy, sizeof(int), 1, fsnap);
	myfree(label, 17);
	fclose(fsnap);
}


void read_Particle_Gadget2(char filename[], int n_start, int n_count)
{
	int n, m, c, p;
	int n_end;
	int dummy, byte, ID;
	unsigned long disp, num_part;

	int DATA_SIZE = sizeof(float);
	int ID_SIZE = sizeof(int);

	float data[3];

	n_end = n_start + n_count;

	FILE* fsnap;
	GadgetHeader head;
	if ( !(fsnap=fopen(filename,"r")) ) {
		printf("[read_Particle_Gadget2] cannot open ic %s \n", filename);
		exit(0);
	}

	printf(" > reading `%s'\n", filename);

	byte = 0;
	byte += fread(&dummy, sizeof(dummy),        1, fsnap);
	byte += fread(&head,  sizeof(GadgetHeader), 1, fsnap);
	byte += fread(&dummy, sizeof(dummy),        1, fsnap);

	Omega0 = head.Omega0;
	OmegaLambda = head.OmegaLambda;
	HubbleParam = head.HubbleParam;
	sp_unit = 1000.0;
	byte += fread(&dummy, sizeof(int), 1, fsnap);   
	for (c=0, p=0, m=0; m<6; m++) {
		for (n=0; n<head.npart[m]; n++) {
			byte += fread(data, DATA_SIZE, 3, fsnap);

			if (c>=n_start && c<n_end) {

#ifndef INTXYZ 
				part[p].pos[0] = (cpuType)data[0] * sp_unit;
				part[p].pos[1] = (cpuType)data[1] * sp_unit;
				part[p].pos[2] = (cpuType)data[2] * sp_unit;
				//part[p].mass   =  head.mass[m];
#else

				part[p].posi[0] = (int) POS2INT * data[0] * sp_unit;
				part[p].posi[1] = (int) POS2INT * data[1] * sp_unit;
				part[p].posi[2] = (int) POS2INT * data[2] * sp_unit;


#endif

				p++;
			}
			c++;
		}
	}
	byte += fread(&dummy, sizeof(int), 1, fsnap);

	cpuType gdt2unit = pow( 1.0/(1.0+head.redshift), 1.5); //for gadget2 format
	byte += fread(&dummy, sizeof(int), 1, fsnap);
	for (c=0, p=0, m=0; m<6; m++) {
		for (n=0; n<head.npart[m]; n++) {
			byte += fread(data, DATA_SIZE, 3, fsnap);

			if (c>=n_start && c<n_end) {
				part[p].vel[0] = (cpuType)data[0] * gdt2unit;
				part[p].vel[1] = (cpuType)data[1] * gdt2unit;
				part[p].vel[2] = (cpuType)data[2] * gdt2unit;
				p++;
			}
			c++;
		}
	}
	byte += fread(&dummy, sizeof(int), 1, fsnap);



#ifdef PARTIDON
	byte += fread(&dummy, sizeof(int), 1, fsnap);
	for (c=0, p=0, m=0; m<6; m++) {
		for (n=0; n<head.npart[m]; n++) {
			byte += fread(&ID, ID_SIZE, 1, fsnap);

			if (c>=n_start && c<n_end) {
				part[p].id = (long)ID;
				p++;
			}
			c++;
		}
	}
	byte += fread(&dummy, sizeof(int), 1, fsnap);
#endif

#ifdef LABEL_BODY
	for (p=0; p<NPART; p++) {
		part[p].id = (long)( p + n_start);
		//	printf(" %d %d %ld\n", PROC_RANK, p, part[p].id);
	}

#endif

	fclose(fsnap);
}


void read_Particle(char filename[], int n_start, int n_count)
{
	int n;
	int dummy, byte;
	unsigned long disp, num_part;

	int DATA_SIZE = sizeof(float);
	float data[3];

	FILE* fsnap;
	GadgetHeader head;
	if ( !(fsnap=fopen(filename,"r")) ) {
		printf("[read_Particle] cannot open file %s \n", filename);
		exit(0);
	}

	printf(" > reading file `%s' \n", filename);

	byte = 0;
	byte += fread(&dummy, sizeof(dummy),        1, fsnap);
	byte += fread(&head,  sizeof(GadgetHeader), 1, fsnap);
	byte += fread(&dummy, sizeof(dummy),        1, fsnap);

	num_part = NPART_TOTAL;

	byte += fread(&dummy, sizeof(int), 1, fsnap);
	for (n=0; n<n_start; n++)
		byte += fread(data, DATA_SIZE, 3, fsnap);

	for (n=0; n<n_count; n++) {
		byte += fread(data, DATA_SIZE, 3, fsnap);

#ifndef INTXYZ
		part[n].pos[0] = (cpuType)data[0] * sp_unit;
		part[n].pos[1] = (cpuType)data[1] * sp_unit;
		part[n].pos[2] = (cpuType)data[2] * sp_unit;
#else
		part[n].posi[0] = (int)data[0] * sp_unit * POS2INT;
		part[n].posi[1] = (int)data[1] * sp_unit * POS2INT;
		part[n].posi[2] = (int)data[2] * sp_unit * POS2INT;
#endif
	}

	//    printf(" [%d] pos byte = %d\n", PROC_RANK, byte);
	for (n=0; n<num_part-n_count-n_start; n++)
		byte += fread(data, DATA_SIZE , 3, fsnap);

	byte += fread(&dummy, sizeof(int), 1, fsnap);

	byte += fread(&dummy, sizeof(int), 1, fsnap);
	for (n=0; n<n_start; n++)
		byte += fread(data, DATA_SIZE, 3, fsnap);

	for (n=0; n<n_count; n++) {
		byte += fread(&data[0], DATA_SIZE, 3, fsnap);
		part[n].vel[0] = (cpuType)data[0];
		part[n].vel[1] = (cpuType)data[1];
		part[n].vel[2] = (cpuType)data[2];
	}
	fclose(fsnap);

}

void write_Particle_Gadget2(char filename[], int n_start, int n_count, cpuType Redshift_Time)
{
	int n, m, c, p;
	int n_end;
	int dummy, byte, ID;
	unsigned long disp, num_part;

	int DATA_SIZE = sizeof(float);
	int ID_SIZE = sizeof(int);

	float data[3];

	n_end = n_start + n_count;

	FILE* fsnap;
	GadgetHeader head;
	if ( !(fsnap=fopen(filename,"w")) ) {
		printf("[read_Particle] cannot open file %s \n", filename);
		exit(0);
	}

	printf(" > writing `%s' \n", filename);

	head.num_files = 12;
	head.BoxSize = BOXSIZE;

	head.Omega0 = OmegaM0 ;
	head.OmegaLambda = OmegaX0;
	head.HubbleParam = Hubble0 ;

	head.npart[0] = 0; 
	head.npart[1] = n_count;
	head.npart[2] = 0; 
	head.npart[3] = 0; 
	head.npart[4] = 0; 
	head.npart[5] = 0; 

	head.mass[0]  = 0.0;
	head.mass[1]  = MASSPART;
	head.mass[2]  = 0.0;
	head.mass[3]  = 0.0;
	head.mass[4]  = 0.0;
	head.mass[5]  = 0.0;




	uint4byte npartTotalLowWord = (NPART_TOTAL & 0xFFFFFFFF);



	head.npartTotal[0] = 0; 
	head.npartTotal[1] = npartTotalLowWord; 
	head.npartTotal[2] = 0; 
	head.npartTotal[3] = 0; 
	head.npartTotal[4] = 0; 
	head.npartTotal[5] = 0; 

	head.npartTotalHighWord[1] = NPART_TOTAL >> 32;   /*!< High word of the total number of particles of each type */

	printf(" npart =  %u %u\n", npartTotalLowWord, head.npartTotalHighWord[1]);

	head.time = 1.0/(Redshift_Time+1);
	head.redshift = Redshift_Time;

	dummy  = 256;
	byte = 0;
	byte += fwrite(&dummy, sizeof(dummy),        1, fsnap);
	byte += fwrite(&head,  sizeof(GadgetHeader), 1, fsnap);
	byte += fwrite(&dummy, sizeof(dummy),        1, fsnap);

sp_unit = 0.001;

	byte += fwrite(&dummy, sizeof(int), 1, fsnap);   
	for (c=n_start; c<n_end; c++) {

#ifndef INTXYZ
		data[0] = (float)part[c].pos[0] * sp_unit;
		data[1] = (float)part[c].pos[1] * sp_unit;
		data[2] = (float)part[c].pos[2] * sp_unit;
#else

		data[0] = (float)part[c].posi[0] * INT2POS *  sp_unit;
		data[1] = (float)part[c].posi[1] * INT2POS * sp_unit;
		data[2] = (float)part[c].posi[2] * INT2POS * sp_unit;

#endif

		byte += fwrite(data, DATA_SIZE, 3, fsnap);
	}
	byte += fwrite(&dummy, sizeof(int), 1, fsnap);

	//printf(" size f= %d\n", sizeof(GadgetHeader));
	cpuType gdt2unit = pow( 1.0/(1.0+head.redshift), 1.5); 
	byte += fwrite(&dummy, sizeof(int), 1, fsnap);   
	for (c=n_start; c<n_end; c++) {

		data[0] = (float)part[c].vel[0] / gdt2unit;
		data[1] = (float)part[c].vel[1] / gdt2unit;
		data[2] = (float)part[c].vel[2] / gdt2unit;

		byte += fwrite(data, DATA_SIZE, 3, fsnap);
	}
	byte += fwrite(&dummy, sizeof(int), 1, fsnap);   

	fclose(fsnap);

	//MPI_Barrier(MPI_COMM_WORLD);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////

/*
void simple_snapshot(char fbase[]) {
	int n;
	FILE *fd;
	char fname[128];
	sprintf(fname, "%s.%d-%d", fbase, PROC_SIZE, PROC_RANK);

	if ( !(fd = fopen(fname,"w")) ) {
		printf("\n\nERROR! Cannot open/write snapshot file\n"); 
		printf("  - '%s' - (simple_snapshot)\n\n", fname);
		exit(0);
	}
	for (n=0; n<NPART; n++) {

		if (part[n].pos[2] < 0.1*BOXSIZE) {
			fprintf(fd, "%lf %lf %lf ", part[n].pos[0], part[n].pos[1], part[n].pos[2]);
			//		fprintf(fd, "%lf %lf %lf ", part[n].vel[0], part[n].vel[1], part[n].vel[2]);
			fprintf(fd, "\n");
		}
	}
	fclose(fd);
}
*/

void load_gadget_ics(char file_path[], int nfiles) {

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

	char file_path_name[128];

	int n;

	//printf(" BOX %d , %f %f %f : %f %f %f\n", PROC_RANK, BOX_DOM_L[0], BOX_DOM_L[1], BOX_DOM_L[2], BOX_DOM_H[0], BOX_DOM_H[1], BOX_DOM_H[2]);	

	NPART_TOTAL = (unsigned long long)(NSIDEMESH) * NSIDEMESH * NSIDEMESH;

	int npart_proc = (int) ( NPART_TOTAL /( PROC_SIZE - NDOM_HEAD * NSIDE0SUDOM * NSIDE1SUDOM ) );
	if (isSUDOM) {
		NPART = (int) ( NPART_TOTAL /( PROC_SIZE - NDOM_HEAD * NSIDE0SUDOM * NSIDE1SUDOM ) );
	}
	if (COM_DOM)
		NPART = (int) ( npart_proc );  

	MAXNPART = (int) (NPART * max_npart_ratio ) ; 





	if (COM_DOM )
		part = (Body*)mymalloc(sizeof(Body) * MAXNPART, 9);

	int start = 0;

	for (n=0; n<nfiles; n++) {
		sprintf(file_path_name, "%s.%d", file_path, n);
		int npart_infile;
		read_GadgetHeader(file_path_name,&npart_infile);	
		read_Particle_Gadget2_domain(file_path_name, &start);
	//	printf(" npart_infile[%d] = %d, %d\n", n, npart_infile, start);
//		printf("[%d] npart(%d) in %d\n",PROC_RANK, n, start);
	}
	NPART = start;

	

	int npart_gather;
	MPI_Reduce(&NPART, &npart_gather, 1, MPI_INT, MPI_SUM, 0,  MPI_COMM_WORLD );
	if ( 0 == PROC_RANK)
	printf(" npart_gather = %d, %lld\n", npart_gather, NPART_TOTAL);

	MAXNPART_BND = (int)  MAXNPART * max_bnd_ratio;
	MAXNPART_BND_SUD = MAXNPART_BND * NDOMinSUDOM;

	MAXNPART_DOM = MAXNPART_BND_SUD;
	MAXNPART_DOM_SUD = MAXNPART_BND_SUD;

	printf(" [%d] NPART = %d, MAXNPART = %d\n", PROC_RANK, NPART, MAXNPART);


	if ( COM_DOM  ) {
		part_b = (BodyBnd*)mymalloc(sizeof(BodyBnd)*MAXNPART_BND, 10);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (0 == PROC_RANK)
	printf(" compele read Gadget snap\n\n");
	//exit(0);
}


#endif

