#ifndef PHOTONS_H
#define PHOTONS_H
#include "typedef.h"
#include "setparam.h"
#include "utility.h"
#include "memory.h"
#include "tasklog.h"

///cosmos
extern long int SEED_IC;
extern cpuType OmegaM0;
extern cpuType OmegaB0;
extern cpuType OmegaX0;
extern cpuType Hubble0;
extern cpuType Sigma8;
extern cpuType PrimordialIndex;
extern cpuType InitialTime;
extern char OutputPath[128];
extern char OutputName[128];
//int PROC_TOTAL;
extern int PROC_SIZE;
extern int PROC_RANK;
extern int PROC_RANK_SUDOM;
extern int NEW_RANK;

extern int SUDOM_RANK;

extern unsigned long long NPART_TOTAL;
extern int  NPART;
extern int  NPART_EXP;
extern int  NPART_BND;
extern int  NPART_BND_SUD;
extern int  NPART_BND_COM;
extern int  NPART_DOM_SUD;
extern int  NPART_DOM_COM;

extern int  TAG_FULL_GRAVITY;

extern int MAXNPART;
extern int MAXNPART_BND; 
extern int MAXNPART_BND_SUD; 
extern int MAXNPART_DOM;
extern int MAXNPART_DOM_SUD;
//Sudomain *sudom;

extern int dom_grp[NDOMinSUDOM];
extern int dom_grp_mesh_start[NDOMinSUDOM];
extern int dom_grp_mesh_end[NDOMinSUDOM];
extern int dom_grp_mesh_size[NDOMinSUDOM];
extern int dom_grp_rank;

extern int PROC_RANK_SEND_NGB[8];
extern int PROC_RANK_RECV_NGB[8];

extern cpuType dom_bnd_l[NDOMinSUDOM];
extern cpuType dom_bnd_h[NDOMinSUDOM];
//int sudom_list[NSUDOM];

extern cpuType GravConst;

extern cpuType BOXSIZE;
extern cpuType MASSPART;

//extern double BOX0L, BOX0H;
//extern double BOX1L, BOX1H;
extern cpuType BOX0L, BOX0H;
extern cpuType BOX1L, BOX1H;
 
extern int COM_DOM;
extern int isSUDOM;
extern int *proc_rank_sudom;

//extern double BOX_DOM_L[3];
//extern double BOX_DOM_H[3];
extern cpuType BOX_DOM_L[3];
extern cpuType BOX_DOM_H[3];
 
//////////////////////

#ifdef OMPI_MPI_H
extern MPI_Datatype Bodytype;
extern MPI_Datatype BodyBndtype;

#ifdef CPU_DP
#define MPI_CPUTYPE MPI_DOUBLE
#else
#define MPI_CPUTYPE MPI_FLOAT
#endif

#endif

extern Body *part;
extern Body *part_dom; //recv
extern Body *part_dom_sud;
extern BodyBnd *part_b;
extern BodyBnd *part_bnd; //send
extern BodyBnd *part_com; //recv

extern cpuType *acc[3];
extern cpuType *acc_pm[3];

extern cpuType *pmesh;
extern cpuType *gmesh;
extern cpuType *mesh;
extern MeshBnd *bmesh; // sudom 2 dom
extern MeshBnd *bmeshc; // sudom 2 dom
extern MeshBnd *cmesh; // sudom 2 sudom

//Body* part_dom; // why again?
extern int NPART_IN;
//Node *btree; // why again?
extern TopNode *toptree;
extern LETList *let;

extern int first_node;
extern int first_leaf;
extern int last_node;
extern int last_leaf;
extern int MAXLEAF;
extern int MAXNLEAF;
extern int MAXNNODE;
extern int TOPNODE_LEVEL;
extern int TOPNODE_SIZE;
extern int LOCNODE_SIZE;
extern int TOPNODE_START;
 
extern int NumThread;

extern int NUM_OUTPUT_SCHEDUE;
 
extern cpuType cutoffRadius;
extern cpuType splitRadius;
 
extern cpuType SoftenScale; 
extern cpuType SoftenLength; 
extern cpuType PIisq;
 
extern int NBUCK;
extern int *np_buck;
extern int *ip_buck;
 
extern int *task_t;
extern int *task_s;
 
extern int idxP2P;
//int idxtask;
extern int idxM2L;

extern long long leveltot[MAX_STEP_LEVEL+1];
extern P2PPack* packlist;
extern int MAXNPACK;
extern int NPACK;
extern int packLevel;

extern double dtime_task;
extern double dtime_prep;
extern double dtime_m2l;
extern double dtime_fmm;
extern double dtime_p2p;
extern double dtime_pm;
extern double dtime_adptv_task;

extern cpuType open_angle;
extern int m2l_count;
extern int walk_m2l_count;

extern cpuType max_npart_ratio;
extern cpuType max_bnd_ratio;
extern int LEN_TASK;

extern int MODE_IC;
extern int OUTPUT_DENSITY_MESH;

extern char PathSnapshot[256];
extern char InputSnapshot[256];
extern char PathMeshOut[256];

cpuType t_flat_lcdm_a(cpuType a);
cpuType a_flat_lcdm_t(cpuType time);
cpuType drift_loga(cpuType loga_i, cpuType loga_f);
cpuType kick_loga(cpuType loga_i, cpuType loga_f);
float ran3(long *idum);
float ran2(long *idum);

void ic_ZA(cpuType, int);
void ic_ZA_nside(cpuType, int);

void gravity();
void gravity_act(int active_level);
void acc_cln();
void gravity_cln();
void acc_pm_cln();
void gravity_cln();

void boundary() ;
#define convolution convolution_
#define conv_pmonly  conv_pmonly_
#define decomp_2d_init get_local_size_
#define decomp_2d_free decomp_2d_free_
#define ic_grad_fft icgrad_
#define ic_fft icifft_

#define powsk powsk_
//#define ps_ic_k ps_ic_k_

#define NPAD 3
void convolution(cpuType *mesh, int nside[], cpuType param[]);

void rand_ic() ;
void rand_ic_2() ;


void write_particle_text(char filename[], int n_start, int n_count);

void checkpart(int label);
void checkbound(Body* pt, int np, int label);
void initial_restart(char file_path[], int snap_number);
void fmm_prepare(int direct, cpuType bdl[3], cpuType bdr[3]);
void fmm_task();
void fmm_task_act(int level);
void fmm();
int print_msg();
///////////
double growth(double a) ;

void setup_mesh_output( int r_snap);


void power_spectrum(int rank_snap) ;

void partmesh_kick(cpuType dkick , int mode_alloc );

//#define BITWIDTH 33554432
//#define BITWIDTH 134217728 
//#define BITWIDTH 268435456 // 28

#define BITWIDTH 1073741824 //30
//#define BITWIDTH 2147483647 //31

//#define BITWIDTH 16777216 // 24
//#define BITWIDTH 33554432 // 25
//#define BITWIDTH 67108864 // 26

//#define BITWIDTH 2147483640

extern double bitwidth_boxsize;
extern double POS2INT;
extern double INT2POS;

extern int domain_center_i[3];

#endif
