/*
 *     photoNs-2
 *
 *     2017 - 8 - 21
 *	qwang@nao.cas.cn
 *
 *     modified by Q. Wang on 19 Apr 2019
 */	 
#ifdef INTRA_LB

#include <omp.h>
#include <math.h>
#include <pthread.h>
#include "photoNs.h"
#include <stdio.h>
#include <stdlib.h>
#include "operator.h"
static int ntask[2];

static int idxtask;

int *ts;
int *tt;

void add_one_p2ppack(int nt, int* ltask_t, int* ltask_s) {
	packlist[NPACK].ntask = nt;
	packlist[NPACK].task_s = (int*)malloc(sizeof(int) * nt);
	packlist[NPACK].task_t = (int*)malloc(sizeof(int) * nt);
	int i;
	for (i=0; i<nt; i++) {
		packlist[NPACK].task_s[i] = ltask_s[i];
		packlist[NPACK].task_t[i] = ltask_t[i];
	}
#ifdef MTHKDTree
	printf("LB module for multi-thread kdtree is not available.\n");
	exit(0);
#else
	int minleaf = last_leaf, maxleaf = first_leaf;
	for (i=0; i<nt; ++i) {
		minleaf = fmin(minleaf, ltask_t[i]);
		minleaf = fmin(minleaf, ltask_s[i]);
		maxleaf = fmax(maxleaf, ltask_t[i]);
		maxleaf = fmax(maxleaf, ltask_s[i]);
	}
	packlist[NPACK].ileaf = minleaf;
	packlist[NPACK].nleaf = maxleaf - minleaf + 1;
	packlist[NPACK].ipart = leaf[minleaf].ipart;
	packlist[NPACK].npart = leaf[maxleaf].ipart - leaf[minleaf].ipart + leaf[maxleaf].npart;
#endif
	packlist[NPACK].rankid = -1;
	NPACK++;
}

void walk_task_p2p_LB(int im, int jm)
{
	if ( -1 == im  || -1 == jm)
		return;

	if ( im < first_leaf || jm < first_leaf ) {
		printf("error\n");
		exit(0);
	}

	if ( im == jm ) {
		if ( im <  last_leaf) {
			if (leaf[im].active >= 0) {	
				*(ts + idxtask) = jm;
				*(tt + idxtask) = im;
				idxtask ++;
			}

			if ( idxtask == LEN_TASK ) {
				printf("[%d] len of p2ptasks is too large %d >= %d, LEN_TASK*1.2, realloc...\n", PROC_RANK, idxtask, LEN_TASK);
				LEN_TASK *= 1.2;
				hostrealloc_task(LEN_TASK);
				ts = task_s;
				tt = task_t;
			}
		}
		if ( im >= first_node ) {
			if (btree[im].son[0] >= first_node)
				btree[btree[im].son[0]].level = btree[im].level + 1;
			if (btree[im].son[1] >= first_node)
				btree[btree[im].son[1]].level = btree[im].level + 1;

			int org = idxtask;

			walk_task_p2p_LB(btree[im].son[0], btree[jm].son[0]);
			walk_task_p2p_LB(btree[im].son[0], btree[jm].son[1]);
			walk_task_p2p_LB(btree[im].son[1], btree[jm].son[0]);
			walk_task_p2p_LB(btree[im].son[1], btree[jm].son[1]);

			int np2p = idxtask - org;

			btree[im].inner_p2p = np2p;
			if (btree[im].level == packLevel && np2p > 0) {
				add_one_p2ppack(np2p, tt + org, ts + org);
				idxtask = org;
				idxP2P += np2p;
			}
		}
		return;
	}
	cpuType dx, dy, dz, r2, dr, wi, wj;
	cpuType dist[3];

	if ( im < last_leaf && jm < last_leaf ) {
		if (leaf[im].active >= 0) {	
			*(ts + idxtask) = jm;
			*(tt + idxtask) = im;
			idxtask ++;
		}

		if ( idxtask == LEN_TASK ) {
			printf("[%d] len of p2ptasks is too large %d >= %d, LEN_TASK*1.2, realloc...\n", PROC_RANK, idxtask, LEN_TASK);
			LEN_TASK *= 1.2;
			hostrealloc_task(LEN_TASK);
			ts = task_s;
			tt = task_t;
		}
		return;
	}

	if ( im < first_node && jm >= first_node ) 
	{
		dx = leaf[im].center[0] - btree[jm].center[0];
		dy = leaf[im].center[1] - btree[jm].center[1];
		dz = leaf[im].center[2] - btree[jm].center[2];
		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(leaf[im].width, btree[jm].width, dist);

		if ( 1 == flag ) {
		}
		else if (0 == flag) {
			walk_task_p2p_LB(im, btree[jm].son[0]);
			walk_task_p2p_LB(im, btree[jm].son[1]);
		}
		else if ( -1 == flag){
			return;
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}
		return ;
	}	

	if ( im >= first_node && jm < first_node ) {
		dx = btree[im].center[0] - leaf[jm].center[0];
		dy = btree[im].center[1] - leaf[jm].center[1];
		dz = btree[im].center[2] - leaf[jm].center[2];
		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(leaf[jm].width, btree[im].width, dist);

		if ( 1 == flag ) {
		}
		else if (0 == flag) {
			walk_task_p2p_LB(btree[im].son[0], jm);
			walk_task_p2p_LB(btree[im].son[1], jm);
		}
		else if ( -1 == flag){
			return;
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}
		return ;
	}

	if ( im >= first_node && jm >= first_node ) {
		dx = btree[im].center[0] - btree[jm].center[0];
		dy = btree[im].center[1] - btree[jm].center[1];
		dz = btree[im].center[2] - btree[jm].center[2];

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(btree[im].width, btree[jm].width, dist);

		if ( 1 == flag ) {

		}
		else if (0 == flag) {
			if ( btree[im].width[0]+btree[im].width[1] + btree[im].width[2] 
					> btree[jm].width[0]+btree[jm].width[1] + btree[jm].width[2] ) 
			{
				walk_task_p2p_LB(btree[im].son[0], jm);
				walk_task_p2p_LB(btree[im].son[1], jm);
			}
			else {
				walk_task_p2p_LB(im, btree[jm].son[0]);
				walk_task_p2p_LB(im, btree[jm].son[1]);
			}
		}
		else if ( -1 == flag){
			return;
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}
		return ;
	}
}


void task_prepare_p2p_LB() {
	idxtask = 0;
	ts = task_s;
	tt = task_t;

	btree[first_node].level = 0;
	walk_task_p2p_LB(first_node, first_node);

	// add the inter-nodes p2p list
	if (idxtask > 0) {
		add_one_p2ppack(idxtask, tt, ts);
		packlist[NPACK-1].rankid = 0;
		idxP2P += idxtask;
	}

//	int i, j;
//	for (i=0; i<MAXNPACK; ++i)
//		printf("PACK(%d) ipart %d-%d ileaf %d-%d p2p %d task_t %lld task_s %lld on g%d\n", i, packlist[i].ipart, packlist[i].ipart+packlist[i].npart, packlist[i].ileaf, packlist[i].ileaf+packlist[i].nleaf, packlist[i].ntask, packlist[i].task_t, packlist[i].task_s, packlist[i].rankid);
/*
	int maxnlvl = (int)(log2((double)(MAXNLEAF))) + 1;
	int* nlevel = (int*)malloc(sizeof(int) * maxnlvl);
	for (i=0; i<maxnlvl; ++i)
		nlevel[i] = 0;
	for (i=first_node; i<=last_node; i++) {
		if (btree[i].level == packLevel) {
			printf("inner_p2p of btree[%d] on level %d is %d\n", i, btree[i].level, btree[i].inner_p2p);
		}
		
		if (btree[i].level < maxnlvl)
			nlevel[btree[i].level] += btree[i].inner_p2p;
	}
	for (i=0; i<maxnlvl; ++i)
		printf("nlevel(%d) %d\n", i, nlevel[i]);
	free(nlevel);
*/
}

void walk_task_p2p_act_LB(int im, int jm, int act)
{
	if ( -1 == im  || -1 == jm)
		return;

	if ( im < first_leaf || jm < first_leaf ) {
		printf("error\n");
		exit(0);
	}

	if ( im == jm ) {
		if ( im <  last_leaf) {
			if (act == leaf[im].active) {
				*(ts + idxtask) = jm;
				*(tt + idxtask) = im;
				idxtask ++;
			}

			if ( idxtask == LEN_TASK ) {
				printf("[%d] len of p2ptasks is too large %d >= %d, LEN_TASK*1.2, realloc...\n", PROC_RANK, idxtask, LEN_TASK);
				LEN_TASK *= 1.2;
				hostrealloc_task(LEN_TASK);
				ts = task_s;
				tt = task_t;
			}
		}
		if ( im >= first_node ) {
			if (btree[im].son[0] >= first_node)
				btree[btree[im].son[0]].level = btree[im].level + 1;
			if (btree[im].son[1] >= first_node)
				btree[btree[im].son[1]].level = btree[im].level + 1;

			int org = idxtask;

			walk_task_p2p_act_LB(btree[im].son[0], btree[jm].son[0], act);
			walk_task_p2p_act_LB(btree[im].son[0], btree[jm].son[1], act);
			walk_task_p2p_act_LB(btree[im].son[1], btree[jm].son[0], act);
			walk_task_p2p_act_LB(btree[im].son[1], btree[jm].son[1], act);

			int np2p = idxtask - org;

			btree[im].inner_p2p = np2p;
			if (btree[im].level == packLevel && np2p > 0) {
				add_one_p2ppack(np2p, tt + org, ts + org);
				idxtask = org;
				idxP2P += np2p;
			}
		}
		return;
	}
	cpuType dx, dy, dz, r2, dr, wi, wj;
	cpuType dist[3];

	if ( im < last_leaf && jm < last_leaf ) {
		if (act == leaf[im].active) {
			*(ts + idxtask) = jm;
			*(tt + idxtask) = im;
			idxtask ++;
		}

		if ( idxtask == LEN_TASK ) {
			printf("[%d] len of p2ptasks is too large %d >= %d, LEN_TASK*1.2, realloc...\n", PROC_RANK, idxtask, LEN_TASK);
			LEN_TASK *= 1.2;
			hostrealloc_task(LEN_TASK);
			ts = task_s;
			tt = task_t;
		}
		return;
	}

	if ( im < first_node && jm >= first_node ) 
	{
		dx = leaf[im].center[0] - btree[jm].center[0];
		dy = leaf[im].center[1] - btree[jm].center[1];
		dz = leaf[im].center[2] - btree[jm].center[2];
		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(leaf[im].width, btree[jm].width, dist);

		if ( 1 == flag ) {
		}
		else if (0 == flag) {
			walk_task_p2p_act_LB(im, btree[jm].son[0], act);
			walk_task_p2p_act_LB(im, btree[jm].son[1], act);
		}
		else if ( -1 == flag){
			return;
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}
		return ;
	}	

	if ( im >= first_node && jm < first_node ) {
		dx = btree[im].center[0] - leaf[jm].center[0];
		dy = btree[im].center[1] - leaf[jm].center[1];
		dz = btree[im].center[2] - leaf[jm].center[2];
		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(leaf[jm].width, btree[im].width, dist);

		if ( 1 == flag ) {
		}
		else if (0 == flag) {
			walk_task_p2p_act_LB(btree[im].son[0], jm, act);
			walk_task_p2p_act_LB(btree[im].son[1], jm, act);
		}
		else if ( -1 == flag){
			return;
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}
		return ;
	}

	if ( im >= first_node && jm >= first_node ) {
		dx = btree[im].center[0] - btree[jm].center[0];
		dy = btree[im].center[1] - btree[jm].center[1];
		dz = btree[im].center[2] - btree[jm].center[2];

		dist[0] = dx;
		dist[1] = dy;
		dist[2] = dz;

		int flag = acceptance(btree[im].width, btree[jm].width, dist);

		if ( 1 == flag ) {

		}
		else if (0 == flag) {
			if ( btree[im].width[0]+btree[im].width[1] + btree[im].width[2] 
					> btree[jm].width[0]+btree[jm].width[1] + btree[jm].width[2] ) 
			{
				walk_task_p2p_act_LB(btree[im].son[0], jm, act);
				walk_task_p2p_act_LB(btree[im].son[1], jm, act);
			}
			else {
				walk_task_p2p_act_LB(im, btree[jm].son[0], act);
				walk_task_p2p_act_LB(im, btree[jm].son[1], act);
			}
		}
		else if ( -1 == flag){
			return;
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}
		return ;
	}
}


void task_prepare_p2p_act_LB() {
	idxtask = 0;
	ts = task_s;
	tt = task_t;

	btree[first_node].level = 0;
	walk_task_p2p_LB(first_node, first_node);

	// add the inter-nodes p2p list
	if (idxtask > 0) {
		add_one_p2ppack(idxtask, tt, ts);
		packlist[NPACK-1].rankid = 0;
		idxP2P += idxtask;
	}

	int i, j;
//	for (i=0; i<MAXNPACK; ++i)
//		printf("PACK(%d) ipart %d-%d ileaf %d-%d p2p %d task_t %lld task_s %lld on g%d\n", i, packlist[i].ipart, packlist[i].ipart+packlist[i].npart, packlist[i].ileaf, packlist[i].ileaf+packlist[i].nleaf, packlist[i].ntask, packlist[i].task_t, packlist[i].task_s, packlist[i].rankid);
/*
	int maxnlvl = (int)(log2((double)(MAXNLEAF))) + 1;
	int* nlevel = (int*)malloc(sizeof(int) * maxnlvl);
	for (i=0; i<maxnlvl; ++i)
		nlevel[i] = 0;
	for (i=first_node; i<=last_node; i++) {
		if (btree[i].level == packLevel) {
			printf("inner_p2p of btree[%d] on level %d is %d\n", i, btree[i].level, btree[i].inner_p2p);
		}
		
		if (btree[i].level < maxnlvl)
			nlevel[btree[i].level] += btree[i].inner_p2p;
	}
	for (i=0; i<maxnlvl; ++i)
		printf("nlevel(%d) %d\n", i, nlevel[i]);
	free(nlevel);
*/
}

#endif
