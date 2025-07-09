#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "photoNs.h"
#include "typedef.h"
#include "setparam.h"



void setup_param(){

	open_angle = 0.4;
//	open_angle = 0.5;

	SoftenLength = 0.02*BOXSIZE/NSIDEMESH;
//	printf(" *** %f %f, %d\n", SoftenLength, BOXSIZE, NSIDEMESH);
//	SoftenLength = 3.0;
	
//	SoftenLength = 0.5;

	//SoftenScale = 2.5 * SoftenLength;
	SoftenScale = 1.5 * SoftenLength;


//	SoftenScale *= 2.5;
//	SoftenScale *= 3.0;

	splitRadius = 1.25*BOXSIZE/NSIDEMESH;
	cutoffRadius = 4.5*splitRadius;

	bitwidth_boxsize = (double) BOXSIZE + cutoffRadius * 1.1 ;
//	bitwidth_boxsize = (double) BOXSIZE ;
	POS2INT = BITWIDTH/bitwidth_boxsize;
        INT2POS = bitwidth_boxsize/BITWIDTH;

	double tmp = 1.0 * (1+PROC_RANK);
	int tmpi = (int)( tmp * POS2INT );
//	printf("P2I = %e, I2P = %e mul = %e i-%d [%d]\n", POS2INT, INT2POS, (tmpi * INT2POS)/tmp, tmpi , PROC_RANK);

	LEN_TASK = 80000000 * 1 * pow(NSIDEMESH / 512., 3) /  (PROC_SIZE - NSIDE0SUDOM * NSIDE1SUDOM * NDOM_HEAD);

	// For small side length, these values need to be
	// appropriately increased.  
	max_npart_ratio = 1.5;


	max_bnd_ratio = 0.35;
//	if (NSIDEMESH < 1000) max_bnd_ratio = 1.0;
//	max_bnd_ratio = 0.4; // modify on 419 test P369

	// for load-balance module
	packLevel = 3;

	MAXLEAF = 256;
//	MAXLEAF = 128;
//	MAXLEAF = 100;

	PIisq = 1.0/sqrt(M_PI);
	GravConst = 43007.105732;

	NDOM_COM = NSIDE0PROC*NSIDE1PROC - NDOM_HEAD *   NSIDE0SUDOM * NSIDE1SUDOM  ;

}


void set_cospar_omM_omX_h0_box(cpuType omegam, cpuType omegax, cpuType hubble, cpuType box) {

	double h2 = hubble;

	OmegaM0 = omegam;
	OmegaX0 = omegax;
	OmegaB0 = 0.0;
	Hubble0 = hubble;
	Sigma8 = 0.0;
	PrimordialIndex = 1.0;
	BOXSIZE = box;	

//	SEED_IC = -857291;
	SEED_IC = -1857291;
//	SEED_IC = -31857291;
}



void setup_plank_param(cpuType ombh2, cpuType omch2, cpuType omnuh2, cpuType omk, cpuType h0, cpuType s8, cpuType ns, cpuType box) {

	double h2 = h0*h0/10000.0;

	OmegaM0 = (ombh2+omnuh2+omch2)/h2;
	OmegaX0 = 1.0-OmegaM0-omk;
	OmegaB0 = ombh2/h2;
	Hubble0 = h0/100.0;
	Sigma8 = s8;
	PrimordialIndex = ns;
	BOXSIZE = box;	

	SEED_IC = -1857291; // simulation_4T
//	SEED_IC = -31857291;

//	SEED_IC = -3157291; /// o3
//	SEED_IC = -857291; /// o2
//	SEED_IC = -31857291; /// O1
}



void setup_cosmo_param(cpuType Om, cpuType Ob, cpuType Ox, cpuType h0, cpuType s8, cpuType ns, cpuType box) {

	OmegaM0 = Om;
	OmegaX0 = Ox;
	OmegaB0 = Ob;
	Hubble0 = h0;
	Sigma8 = s8;
	PrimordialIndex = ns;
	BOXSIZE = box;	

	SEED_IC = -3857291;
//	SEED_IC = -13857291;
}



