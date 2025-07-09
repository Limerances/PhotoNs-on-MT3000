#ifndef TYPEDEF_H
#define TYPEDEF_H

#ifdef CPU_DP
typedef double cpuType;
#else
typedef float cpuType;
#endif



#define QUADRUPOLE
#define  OCTUPOLE
#define WRITE_SNAP_ON

extern int NDOM_COM;

typedef unsigned long long INT64;
typedef struct {
	cpuType pos[3];
	int posi[3];
} BodyBnd;

typedef struct {
	int x,y,z;
	cpuType data;
} MeshBnd;

typedef struct {
	int tag;
	int act;
#ifndef INTXYZ
	cpuType pos[3];
#else
	int posi[3];
#endif



	cpuType vel[3];
//	cpuType acc_pm[3];
//#ifdef HALO_FOF
	unsigned long id;
//#endif

} Body;

#define NSON 2

#ifdef TREECODE
typedef struct {
	int  updated;
	int  npart;
	int  son[NSON];
	cpuType split;
	cpuType width[3];
	cpuType center[3];
} Node;
#endif

typedef struct {
	int  son[NSON];
	cpuType split;
	int npart;
	int ipart;
	int direct;
	cpuType center[3];
	cpuType width[3];
} TopNode;


typedef struct {
	int  npart;
	int  nnode;
	int  ipart;
	int  inode;
} LETList;



#define NMULTI 20
typedef struct {
	int active;
	int npart;
	int ipart;
	cpuType width[3];
	cpuType center[3];
	cpuType M[NMULTI];
	cpuType L[NMULTI];
	cpuType acc_pm[3];
} Pack;


typedef struct {
#ifdef INTRA_LB
	int level;
	int inner_p2p;
#endif
	int  updated;
	int  npart;
	int  son[NSON];
	cpuType split;
	cpuType width[3];
	cpuType center[3];
	cpuType M[NMULTI];
	cpuType L[NMULTI];
} Node;

typedef struct {
	int rankid;
	int ipart;
	int npart;
	int ileaf;
	int nleaf;
	int ntask; 

 	int* task_t;
	int* task_s;	
} P2PPack;

extern Node *btree;
extern Pack *leaf;
#endif
