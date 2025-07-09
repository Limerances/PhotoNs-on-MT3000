#include <mpi.h>
#include "photoNs.h"
#include <stdio.h>
///cosmos
cpuType OmegaM0;
cpuType OmegaB0;
cpuType OmegaX0;
cpuType Hubble0;
cpuType Sigma8;
cpuType PrimordialIndex;
cpuType InitialTime;
char OutputPath[128];
char OutputName[128];
//int PROC_TOTAL;
int PROC_SIZE;
int PROC_RANK;
int PROC_RANK_SUDOM;
int NEW_RANK;

int SUDOM_RANK;

long int SEED_IC;

unsigned long long NPART_TOTAL;
int  NPART;
int  NPART_EXP;
int  NPART_BND;
int  NPART_BND_SUD;
int  NPART_BND_COM;
int  NPART_DOM_SUD;
int  NPART_DOM_COM;


int  TAG_FULL_GRAVITY;

int MAXNPART;
int MAXNPART_BND; 
int MAXNPART_BND_SUD; 
int MAXNPART_DOM;
int MAXNPART_DOM_SUD;
//Sudomain *sudom;

int dom_grp[NDOMinSUDOM];
int dom_grp_mesh_start[NDOMinSUDOM];
int dom_grp_mesh_end[NDOMinSUDOM];
int dom_grp_mesh_size[NDOMinSUDOM];
int dom_grp_rank;

int PROC_RANK_SEND_NGB[8];
int PROC_RANK_RECV_NGB[8];

cpuType dom_bnd_l[NDOMinSUDOM];
cpuType dom_bnd_h[NDOMinSUDOM];
//int sudom_list[NSUDOM];

cpuType GravConst;

cpuType BOXSIZE;
cpuType MASSPART;

//double BOX0L, BOX0H;
//double BOX1L, BOX1H;
cpuType BOX0L, BOX0H;
cpuType BOX1L, BOX1H;

int COM_DOM;
int isSUDOM;
int *proc_rank_sudom;

//double BOX_DOM_L[3];
//double BOX_DOM_H[3];

cpuType BOX_DOM_L[3];
cpuType BOX_DOM_H[3];

//////////////////////

MPI_Datatype Bodytype;
MPI_Datatype BodyBndtype;
Body *part;
Body *part_dom; //recv
Body *part_dom_sud;
BodyBnd *part_b;
BodyBnd *part_bnd; //send
BodyBnd *part_com; //recv

cpuType *acc[3];
cpuType *acc_pm[3];

cpuType *pmesh;
cpuType *gmesh;
cpuType *mesh;
MeshBnd *bmesh; // sudom 2 dom
MeshBnd *bmeshc; // sudom 2 dom
MeshBnd *cmesh; // sudom 2 sudom

//Body* part_dom; // why again?
int NPART_IN;
//Node *btree; // why again?
TopNode *toptree;
LETList *let;

int first_node;
int first_leaf;
int last_node;
int last_leaf;
int MAXLEAF;
int MAXNLEAF;
int MAXNNODE;
int TOPNODE_LEVEL;
int TOPNODE_SIZE;
int LOCNODE_SIZE;
int TOPNODE_START;

int NumThread;

cpuType cutoffRadius;
cpuType splitRadius;

cpuType SoftenScale; 
cpuType SoftenLength; 
cpuType PIisq;

int NBUCK;
int *np_buck;
int *ip_buck;

int *task_t;
int *task_s;

int idxP2P;
int idxM2L;

P2PPack* packlist;
int MAXNPACK;
int NPACK;
int packLevel;

cpuType max_npart_ratio;
cpuType max_bnd_ratio;
double dtime_task;
double dtime_prep;
double dtime_m2l;
double dtime_fmm;
double dtime_p2p;
double dtime_pm;
double dtime_adptv_task;

cpuType open_angle;
int m2l_count;
int walk_m2l_count;

int NDOM_COM;
Node *btree;
Pack *leaf;

int LEN_TASK;

char PathSnapshot[256];
char InputSnapshot[256];
char PathMeshOut[256];

int NUM_OUTPUT_SCHEDUE;


double bitwidth_boxsize;
double POS2INT;
double INT2POS;

int domain_center_i[3];

void setup_mesh_output( int r_snap) {
	make_dir(PathSnapshot, "mesh", r_snap);
	sprintf(PathMeshOut, "%s/mesh_%d/mesh_%03d", PathSnapshot, r_snap, r_snap);
}
