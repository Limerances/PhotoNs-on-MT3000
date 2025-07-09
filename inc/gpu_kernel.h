#ifndef GPUKERNEL_H
#define GPUKERNEL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define WARP 32

#ifdef GPU_DP
typedef double gpuT;
#else
typedef float gpuT;
#endif

typedef unsigned char uint_8b;

//using in P2P host reorder
gpuT* h_part_pos;
gpuT* h_part_buf;
uint_8b* h_part_act;

int* h_leaf_ip;
int* h_leaf_np;


void *task_compute_p2p_gpu(int nt, int* ltask_t, int* ltask_s, int);
void hostmalloc_task(int LEN_TASK);
void hostrealloc_task(int LEN_TASK);
void hostfree_task();
void hostmalloc_in(int NPART, int NLEAF);
void hostfree_in();
void hostmalloc_out(int NPART);
void hostfree_out();
void dinfofree();
void devfree();

void updateAcc_fromHostBuff(P2PPack pack, int active_level);
void devMemcpyDH_updateAcc(P2PPack* packs, int npack, int active_level); 
void sync_with_dev();

double get_last_kntime(int proc_rank);
#endif
