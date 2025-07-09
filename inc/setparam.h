#ifndef PARAM_H
#define PARAM_H



/// node 4
#define NDOMinSUDOM  9
#define NSIDE0PROC 6
#define NSIDE1PROC 6
#define NSIDE0SUDOM 2
#define NSIDE1SUDOM 2
#define NSIDEMESH 1620 
//#define NSIDEMESH 810  


/*
/// node 16
#define NDOMinSUDOM  9
#define NSIDE0PROC 12
#define NSIDE1PROC 12
#define NSIDE0SUDOM 4
#define NSIDE1SUDOM 4
#define NSIDEMESH 2556 
*/

/*
///// node 64
#define NDOMinSUDOM  9
#define NSIDE0PROC 24
#define NSIDE1PROC 24
#define NSIDE0SUDOM 8
#define NSIDE1SUDOM 8
#define NSIDEMESH 4104 
*/

#define NDOM_HEAD 1

#define MAX_STEP_LEVEL 0
//#define MAX_STEP_LEVEL 7


void setup_param();
void set_cospar_omM_omX_h0_box(cpuType omegam, cpuType omegax, cpuType hubble, cpuType box) ;
void setup_cosmo_param(cpuType Om, cpuType Ob, cpuType Ox, cpuType h0, cpuType s8, cpuType ns, cpuType box) ;
void setup_plank_param(cpuType omch2, cpuType ombh2, cpuType omnuh2, cpuType omk, cpuType h0, cpuType s8, cpuType ns, cpuType box) ;

void read_power_table(char FileWithInputSpectrum[]);

#endif
