/*
 *     photoNs-2
 *
 *     2017 - 8 - 21
 *	qwang@nao.cas.cn
 *
 *     modified by Q. Wang on 19 Apr 2019
 */	 


#include <omp.h>
#include <math.h>
#include <pthread.h>
#include "photoNs.h"
#include <stdio.h>
#include <stdlib.h>
#include "operator.h"
static int ntask;

static int stask;
static int ttask;
static int idxtask;

static int Dir;
#ifndef INTXYZ

int cmpos(const void* p1, const void *p2) {
	if ( ((Body*)p1)->pos[Dir] > ((Body*)p2)->pos[Dir] )
		return 1;
	else
		return 0;
}

void sort_ratio(int D, Body *p, int length, int npart[2], cpuType *split, cpuType ratio)
{
	npart[0] = npart[1] = 0;
	if (length < 2)
	{
		npart[1] = length;
		return;
	}
	Body temp;
	if (length == 2)
	{
		npart[0] = npart[1] = 1;
		*split = 0.5 * (p[0].pos[D] + p[1].pos[D]);
		if (p[0].pos[D] > p[1].pos[D])
		{
			temp = p[0];
			p[0] = p[1];
			p[1] = temp;
		}
		return;
	}

	Dir = D;
	int but;
	qsort(p, length, sizeof(Body), cmpos);
	but = (int)(ratio * length);	

	npart[0] = but;
	npart[1] = length - npart[0];

	*split = 0.5*(p[but-1].pos[D] + p[but].pos[D]);
}



void bksort_inplace(int D, Body *p, int length, int npart[2], double *split)
{
	npart[0] = npart[1] = 0;
	if (length <2 ) {
		npart[1] = length;
		return;
	}
	Body temp;
	if (length == 2 )
	{
		npart[0] = npart[1] = 1;
		*split = 0.5*(p[0].pos[D] + p[1].pos[D]);
		if (p[0].pos[D] > p[1].pos[D]) {
			temp = p[0];
			p[0] = p[1];
			p[1] = temp;
		}
		return;
	}
	int n;

	double mean = 0.0;

	int mul =1  ;
	if (length > 100000)
		mul = length / 100000;

	for (n=0; n<length; n+=mul) {
		mean += p[n].pos[D];
	}
	mean /= (double)length;
	mean *= (double)mul;

	int but = length - 1;

	n=0;
	while ( n < but ) {
		if ( p[n].pos[D] > mean ) {
			temp = p[n];

			while ( p[but].pos[D] > mean && but > n )
				but--;

			p[n] = p[but];
			p[but] = temp;
		}
		n++;
	}

	npart[0] = but;
	npart[1] = length - npart[0];

	*split = mean;
}

#else

int cmpos(const void* p1, const void *p2) {
	if ( ((Body*)p1)->posi[Dir] > ((Body*)p2)->posi[Dir] )
		return 1;
	else
		return 0;
}

void sort_ratio(int D, Body *p, int length, int npart[2], cpuType *split, cpuType ratio)
{
	npart[0] = npart[1] = 0;
	if (length < 2)
	{
		npart[1] = length;
		return;
	}
	Body temp;
	if (length == 2)
	{
		npart[0] = npart[1] = 1;
	//	*split = 0.5*( (double) (INT2POS * p[0].posi[D]) + (double)(p[1].posi[D] * INT2POS  ) );
		*split = 0.5*( INT2POS * (p[0].posi[D] + p[1].posi[D] ) );
		if (p[0].posi[D] > p[1].posi[D])
		{
			temp = p[0];
			p[0] = p[1];
			p[1] = temp;
		}
		return;
	}

	Dir = D;
	int but;
	qsort(p, length, sizeof(Body), cmpos);
	but = (int)(ratio * length);	

	npart[0] = but;
	npart[1] = length - npart[0];

	*split = 0.5*( (float)(INT2POS * ( p[but-1].posi[D] + p[but].posi[D])) ) ;
//	*split = 0.5*( (float)(INT2POS * p[but-1].posi[D] + INT2POS * p[but].posi[D])) ;
}



void bksort_inplace(int D, Body *p, int length, int npart[2], double *split)
{
	npart[0] = npart[1] = 0;
	if (length <2 ) {
		npart[1] = length;
		return;
	}
	Body temp;
	if (length == 2 )
	{
		npart[0] = npart[1] = 1;
		*split = 0.5*( INT2POS * ( p[0].posi[D] + p[1].posi[D]   ) );
	//	*split = 0.5*( (double) (INT2POS * p[0].posi[D]) + (double)(p[1].posi[D] * INT2POS  ) );
		if (p[0].posi[D] > p[1].posi[D]) {
			temp = p[0];
			p[0] = p[1];
			p[1] = temp;
		}
		return;
	}
	int n;

	double mean = 0.0;

	int mul =1  ;
	if (length > 100000)
		mul = length / 100000;

	for (n=0; n<length; n+=mul) {
		mean += (double)( INT2POS *  p[n].posi[D]);
	}
	mean /= (double)length;
	mean *= (double)mul;

	int but = length - 1;

	n=0;
	while ( n < but ) {
		if ( (double) (INT2POS * p[n].posi[D]) > mean ) {
			temp = p[n];

			while ( (double) (INT2POS * p[but].posi[D] )> mean && but > n )
				but--;

			p[n] = p[but];
			p[but] = temp;
		}
		n++;
	}

	npart[0] = but;
	npart[1] = length - npart[0];

	*split = mean;
}




#endif

void build_kdtree(int direct, int iPart, int length, int iNode)
{
	if (length == 0)
		return;


	btree[iNode].npart = length;
	Body *p = part + iPart;

	int npart[2];
	double split;

	bksort_inplace(direct, p, length, npart, &split);

	btree[iNode].split = split;

	int n, ip;
	ip = iPart;
	int new_direct = (direct + 1)%3;

	for (n=0; n<NSON; n++) {

		if (npart[n]  <= MAXLEAF) {
			leaf[last_leaf].npart = npart[n];
			leaf[last_leaf].ipart = ip;

			btree[iNode].son[n] = last_leaf;

			last_leaf++;
		}
		else {
			last_node++;
			btree[iNode].son[n] = last_node;
			build_kdtree(new_direct, ip, npart[n],  last_node);
		}

		ip += npart[n];
	}

}


void center_kdtree(int direct, int iNode, cpuType left[3], cpuType right[3])
{
	int n, idx;
	btree[iNode].width[0] = right[0]-left[0];
	btree[iNode].width[1] = right[1]-left[1];
	btree[iNode].width[2] = right[2]-left[2];
	btree[iNode].center[0] = 0.5*(right[0]+left[0]);
	btree[iNode].center[1] = 0.5*(right[1]+left[1]);
	btree[iNode].center[2] = 0.5*(right[2]+left[2]);

	int newd = (direct + 1)%3;
	cpuType tmp;
	for (n=0; n<NSON; n++) {
		idx = btree[iNode].son[n];

		if (idx < last_leaf ) {

			leaf[idx].width[0] = btree[iNode].width[0];
			leaf[idx].width[1] = btree[iNode].width[1];
			leaf[idx].width[2] = btree[iNode].width[2];

			leaf[idx].center[0] = btree[iNode].center[0];
			leaf[idx].center[1] = btree[iNode].center[1];
			leaf[idx].center[2] = btree[iNode].center[2];

			if (0==n){
				leaf[idx].width[direct] = btree[iNode].split-left[direct];
				leaf[idx].center[direct] = 0.5*(left[direct]+btree[iNode].split);
			}

			if (1==n) {
				leaf[idx].width[direct] = right[direct] - btree[iNode].split;
				leaf[idx].center[direct] = 0.5*(right[direct]+btree[iNode].split);
			}


		} else {

			if (0==n) {
				tmp = right[direct];
				right[direct] = btree[iNode].split;
				center_kdtree(newd, idx, left, right);
				right[direct] =tmp;
			}
			if (1==n) {
				tmp = left[direct];
				left[direct] = btree[iNode].split;
				center_kdtree(newd, idx, left, right);
				left[direct] =tmp;

			}
		}
	}

}


void build_localtree(int direct, cpuType bdleft[3], cpuType bdright[3]) {

	dtime_p2p = 0.0;
	int d;

	MAXNNODE = (int)(2.0*((double)NPART)/((double)MAXLEAF));
	MAXNLEAF = (int)(2.0*((double)NPART)/((double)MAXLEAF));

	if (MAXNNODE > NPART)
		MAXNNODE = NPART +1 ;
	if (MAXNLEAF > NPART)
		MAXNLEAF = NPART +1;

	first_leaf = last_leaf = NPART;
	first_node = last_node = NPART+MAXNLEAF;

	// 0...NPART, fist_leaf...last_leaf-1, first_node...last_node

	leaf  = (Pack*)mymalloc(sizeof(Pack)*MAXNLEAF, 40);
	btree = (Node*)mymalloc(sizeof(Node)*MAXNNODE, 41);
#ifdef MOREINFO
	printf(" [%d] MEM Buff : leaf = %d MB, node = %d MB\n", PROC_RANK, sizeof(Pack)*MAXNLEAF/1024/1024, sizeof(Node)*MAXNNODE/1024/1024);
#endif
	int n, m;
	for (n=0; n<MAXNNODE; n++) {
		btree[n].npart = 0;
		btree[n].son[0] = btree[n].son[1] = -1;
		btree[n].split = 0.0;
		btree[n].width[0] = 0.0;
		btree[n].width[1] = 0.0;
		btree[n].width[2] = 0.0;

		btree[n].center[0] = 0.0;
		btree[n].center[1] = 0.0;
		btree[n].center[2] = 0.0;

		for (m=0; m<NMULTI; m++) {
			btree[n].M[m] = 0.0;
			btree[n].L[m] = 0.0;
		}
	}


	for (n=0; n<MAXNLEAF; n++) {
		leaf[n].npart = 0;
		leaf[n].ipart = 0;

		leaf[n].width[0] =0.0;
		leaf[n].width[1] =0.0;
		leaf[n].width[2] =0.0;

		leaf[n].center[0] = 0.0;
		leaf[n].center[1] = 0.0;
		leaf[n].center[2] = 0.0;

		for (m=0; m<NMULTI; m++) {
			leaf[n].M[m] = 0.0;
			leaf[n].L[m] = 0.0;
		}

	}

	btree -= first_node;
	leaf  -= first_leaf;

#ifdef MTHKDTree
	build_kdtree_threads();
#else
	build_kdtree(direct,0, NPART,  first_node);
#endif

#ifdef MOREINFO
	printf("[%d] nnode = %d (%d) nleaf=%d (%d)\n", PROC_RANK, last_node-first_node+1, MAXNNODE,last_leaf-first_leaf, MAXNLEAF);
#endif
	center_kdtree(direct, first_node, bdleft, bdright);

}


// 0==open, 1==accept, -1==abort
int acceptance(cpuType wi[3], cpuType wj[3], cpuType dist[3]);
inline int acceptance(cpuType wi[3], cpuType wj[3], cpuType dist[3]) 
{
	cpuType min[3];
	cpuType w[3];
	cpuType dd2, dm2, comp = 1.0;

	w[0] = ( wi[0] + wj[0] ) * 0.5;	
	w[1] = ( wi[1] + wj[1] ) * 0.5;	
	w[2] = ( wi[2] + wj[2] ) * 0.5;	

	min[0] = dist[0];
	min[1] = dist[1];
	min[2] = dist[2];


	dd2 = dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2];

	if (min[0] < 0.0)
		min[0] = - min[0];
	if (min[1] < 0.0)
		min[1] = - min[1];
	if (min[2] < 0.0)
		min[2] = - min[2];

	min[0] -= w[0];
	min[1] -= w[1];
	min[2] -= w[2];

	if (min[0] <= 0.0)
		min[0] = 0.0;
	if (min[1] <= 0.0)
		min[1] = 0.0;
	if (min[2] <= 0.0)
		min[2] = 0.0;

	// neighbour
	if (min[0] + min[1] + min[2] < 1.0e-15 )
		return 0;


	dm2 =  min[0]*min[0] + min[1]*min[1] + min[2]*min[2] ;

#ifdef LONGSHORT
	cpuType c2 = cutoffRadius * cutoffRadius;
	if ( dm2 > c2) {
		return -1;
	}
	else {
#if FIDUICAL 
		return 0;

	/// compare for fid
#endif
	}

#endif
//////////////////////////////////////////////
	cpuType wmax = w[0];
	if (w[1]> wmax)wmax = w[1];
	if (w[2]> wmax)wmax = w[2];
	wmax = fabs(wmax);
	wmax *= 2;
	if (  wmax*wmax < open_angle * open_angle * dd2 )  {
		return 1;
	}
	else
		return 0;

}
////////////////////////////////////////////////////////////////

///// need stort im, jm and for m2l or p2p list.

void *task_compute_p2p(void *arg) {
	if (print_msg())
		printf("Not implemented yet.\n");
}
void *task_compute_p2p_act(void *arg) {
	if (print_msg())
		printf("Not implemented yet.\n");
}

int argt[3];   //

int *ts;
int *tt;

pthread_t tid;

void turn2compute_p2p(){
	ntask = idxtask;

#ifdef GPUP2P
	task_compute_p2p_gpu(ntask, task_t, task_s, 0);
#else
	task_compute_p2p(argt);
#endif

	ts = task_s;
	tt = task_t;

	idxP2P +=idxtask;
	idxtask = 0;

}

void walk_task_p2p_level(int im, int jm, int level)
{
	if ( -1 == im  || -1 == jm)
		return;

	if ( im < first_leaf || jm < first_leaf ) {
		printf("Error in walk_task_p2p\n");
		exit(0);
	}

	if ( im == jm ) {
		if ( im <  last_leaf) {
			if (leaf[im].active >= level) {	
				*(ts + idxtask) = jm;
				*(tt + idxtask) = im;
				idxtask ++;
			}

			if ( idxtask == LEN_TASK ) {
				turn2compute_p2p();
			}
		}
		if ( im >= first_node ) {
			walk_task_p2p_level(btree[im].son[0], btree[jm].son[0], level);
			walk_task_p2p_level(btree[im].son[0], btree[jm].son[1], level);
			walk_task_p2p_level(btree[im].son[1], btree[jm].son[0], level);
			walk_task_p2p_level(btree[im].son[1], btree[jm].son[1], level);
		}
		return;
	}
	cpuType dx, dy, dz, r2, dr, wi, wj;
	cpuType dist[3];

	if ( im < last_leaf && jm < last_leaf ) {
		if (leaf[im].active >= level) {	
			*(ts + idxtask) = jm;
			*(tt + idxtask) = im;
			idxtask ++;
		}

		if ( idxtask == LEN_TASK ) {

			turn2compute_p2p();

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
			walk_task_p2p_level(im, btree[jm].son[0], level);
			walk_task_p2p_level(im, btree[jm].son[1], level);
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
			walk_task_p2p_level(btree[im].son[0], jm, level);
			walk_task_p2p_level(btree[im].son[1], jm, level);
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
				walk_task_p2p_level(btree[im].son[0], jm, level);
				walk_task_p2p_level(btree[im].son[1], jm, level);
			}
			else {
				walk_task_p2p_level(im, btree[jm].son[0], level);
				walk_task_p2p_level(im, btree[jm].son[1], level);
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


void walk_task_p2p(int im, int jm)
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
				turn2compute_p2p();
			}
		}
		if ( im >= first_node ) {
			walk_task_p2p(btree[im].son[0], btree[jm].son[0]);
			walk_task_p2p(btree[im].son[0], btree[jm].son[1]);
			walk_task_p2p(btree[im].son[1], btree[jm].son[0]);
			walk_task_p2p(btree[im].son[1], btree[jm].son[1]);
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

			turn2compute_p2p();

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
			walk_task_p2p(im, btree[jm].son[0]);
			walk_task_p2p(im, btree[jm].son[1]);
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
			walk_task_p2p(btree[im].son[0], jm);
			walk_task_p2p(btree[im].son[1], jm);
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
				walk_task_p2p(btree[im].son[0], jm);
				walk_task_p2p(btree[im].son[1], jm);
			}
			else {
				walk_task_p2p(im, btree[jm].son[0]);
				walk_task_p2p(im, btree[jm].son[1]);
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


void *task_compute_m2l(int nt);

void turn2compute_m2l(){
	ntask = idxtask;

	task_compute_m2l(ntask);

	ts = task_s;
	tt = task_t;

	idxM2L +=idxtask;
	idxtask = 0;

}



void walk_task_m2l(int im, int jm)
{
	if ( -1 == im  || -1 == jm)
		return;

	if ( (im < first_leaf || jm < first_leaf)
	|| ( im > last_node || jm > last_node )
	|| ( im >= last_leaf && im < first_node)
	|| ( jm >= last_leaf && jm < first_node)) {
		printf("[%d] error node im %d jm %d last_node %d\n", PROC_RANK, im, jm, last_node);
		exit(0);
	}

	if ( im == jm ) {

		if ( im >= first_node ) 
		{
			walk_task_m2l(btree[im].son[0], btree[jm].son[0]);
			walk_task_m2l(btree[im].son[0], btree[jm].son[1]);
			walk_task_m2l(btree[im].son[1], btree[jm].son[0]);
			walk_task_m2l(btree[im].son[1], btree[jm].son[1]);
		}

		return ;
	}
	cpuType dx, dy, dz, r2, dr, wi, wj;
	cpuType dist[3];

	if ( im < last_leaf && jm < last_leaf ) {
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

		if ( 1 == flag) {
			*(ts + idxtask) = jm;
			*(tt + idxtask) = im;
			idxtask ++;
			if ( idxtask == LEN_TASK ) {
				turn2compute_m2l();
			}
		}
		else if (0 == flag) {
			walk_task_m2l(im, btree[jm].son[0]);
			walk_task_m2l(im, btree[jm].son[1]);
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

		if ( 1 == flag) {
			*(ts + idxtask) = jm;
			*(tt + idxtask) = im;
			idxtask ++;
			if ( idxtask == LEN_TASK ) {
				turn2compute_m2l();
			}
		}
		else if (0 == flag) {
			walk_task_m2l(btree[im].son[0], jm);
			walk_task_m2l(btree[im].son[1], jm);
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

		if ( 1 == flag) {
			*(ts + idxtask) = jm;
			*(tt + idxtask) = im;
			idxtask ++;
			if ( idxtask == LEN_TASK ) {
				turn2compute_m2l();
			}

		}
		else if (0 == flag) {
			if ( btree[im].width[0]+btree[im].width[1] + btree[im].width[2] 
					> btree[jm].width[0]+btree[jm].width[1] + btree[jm].width[2] ) 
			{
				walk_task_m2l(btree[im].son[0], jm);
				walk_task_m2l(btree[im].son[1], jm);
			}
			else {
				walk_task_m2l(im, btree[jm].son[0]);
				walk_task_m2l(im, btree[jm].son[1]);
			}
		}
		else if ( -1 == flag){
			return;
		}
		else {
			printf(" error acceptance \n");
			exit(0);
		}
		return;
	}
}


static double dTlocal, dTremote, dTmirror;
static INT64 NTASKP2P;
static INT64 NTASKM2L;

void fmm_prepare(int direct, cpuType bdl[3], cpuType bdr[3]) 
{
	int rank, d;
	int n, mi, mj, mk;
	double t0, t1, t2, t3, t4, t5, t6, t7, t8;
	double dTlocal, dTremote, dTmirror;
	dtime_prep = dtime();

	t0 = dtime();
	m2l_count = 0;

	build_localtree(direct, bdl, bdr);

	for (n=first_leaf; n<last_leaf; n++)
		p2m(leaf[n].ipart, leaf[n].npart, leaf[n].center, leaf[n].M);


	walk_m2m(first_node);

	int idx[8];

	cpuType width_this_domain = bdr[0] - bdl[0];


	if ( width_this_domain < bdr[1]- bdl[1]  )
		width_this_domain = bdr[1] - bdl[1];

	if ( width_this_domain  < bdr[2] - bdl[2])
		width_this_domain = bdr[2] - bdl[2];


	dtime_prep = dtime() - dtime_prep;

	t1 = dtime();
	idxM2L = 0; 
	idxP2P = 0;
}

void task_prepare_p2p() {
	idxtask = 0;
	ts = task_s;
	tt = task_t;
	walk_task_p2p(first_node, first_node);
}


void task_prepare_p2p_level(int level) {
	idxtask = 0;
	ts = task_s;
	tt = task_t;
	walk_task_p2p_level(first_node, first_node, level);
}


void turn2compute_p2p_act(int active){
	ntask = idxtask;

#ifdef GPUP2P
//	task_compute_p2p_act_gpu(ntask, task_t, task_s, 0, active);
	task_compute_p2p_gpu(ntask, task_t, task_s, 0, active);
#else
	task_compute_p2p_act(argt);
#endif

	ts = task_s;
	tt = task_t;

	idxP2P +=idxtask;
	idxtask = 0;

}

void walk_task_p2p_act(int im, int jm, int act)
{
	if ( -1 == im  || -1 == jm)
		return;

	if ( im < first_leaf || jm < first_leaf ) {
		printf("error\n");
		exit(0);
	}

	if ( im == jm ) {
		if ( im <  last_leaf ) {
			if (leaf[im].active  == act ) {
				*(ts + idxtask) = jm;
				*(tt + idxtask) = im;
				idxtask ++;
			}
			if ( idxtask == LEN_TASK ) {
				turn2compute_p2p_act(act);
			}
		}
		if ( im >= first_node ) {
			walk_task_p2p_act(btree[im].son[0], btree[jm].son[0], act);
			walk_task_p2p_act(btree[im].son[0], btree[jm].son[1], act);
			walk_task_p2p_act(btree[im].son[1], btree[jm].son[0], act);
			walk_task_p2p_act(btree[im].son[1], btree[jm].son[1], act);
		}
		return;
	}
	cpuType dx, dy, dz, r2, dr, wi, wj;
	cpuType dist[3];

	if ( im < last_leaf && jm < last_leaf ) {
		if ( leaf[im].active == act) {
			*(ts + idxtask) = jm;
			*(tt + idxtask) = im;
			idxtask ++;
		}
		if ( idxtask == LEN_TASK ) {

			turn2compute_p2p_act(act);

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
			walk_task_p2p_act(im, btree[jm].son[0],act);
			walk_task_p2p_act(im, btree[jm].son[1],act);
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
			walk_task_p2p_act(btree[im].son[0], jm, act);
			walk_task_p2p_act(btree[im].son[1], jm, act);
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
				walk_task_p2p_act(btree[im].son[0], jm, act);
				walk_task_p2p_act(btree[im].son[1], jm, act);
			}
			else {
				walk_task_p2p_act(im, btree[jm].son[0], act);
				walk_task_p2p_act(im, btree[jm].son[1], act);
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

void task_prepare_p2p_act(int act) {
	idxtask = 0;
	ts = task_s;
	tt = task_t;
	walk_task_p2p_act(first_node, first_node, act);
}



void task_prepare_m2l() {
	idxtask = 0;

	ts = task_s;
	tt = task_t;

	walk_task_m2l(first_node, first_node);
}

void *task_compute_m2l(int nt) {
	int n;
	double t0 = dtime();

	for (n=0; n<nt; n++) {
		int im = task_t[n];
		int jm = task_s[n];
		cpuType dx, dy, dz;
		//   omp_set_lock(&writelock);
		if ( im < first_node ) {
			dx = leaf[im].center[0] - btree[jm].center[0];
			dy = leaf[im].center[1] - btree[jm].center[1];
			dz = leaf[im].center[2] - btree[jm].center[2];
			m2l(dx, dy ,dz, btree[jm].M, leaf[im].L);
		}
		else if ( jm < first_node ) {
			dx = btree[im].center[0] - leaf[jm].center[0];
			dy = btree[im].center[1] - leaf[jm].center[1];
			dz = btree[im].center[2] - leaf[jm].center[2];
			m2l(dx, dy ,dz, leaf[jm].M, btree[im].L);
		}
		else {
			dx = btree[im].center[0] - btree[jm].center[0];
			dy = btree[im].center[1] - btree[jm].center[1];
			dz = btree[im].center[2] - btree[jm].center[2];
			m2l(dx, dy ,dz, btree[jm].M, btree[im].L);
		}	
	}
	dtime_m2l = dtime() - t0;
}

#ifdef GPUP2P

void fmm_task()
{
	int n;
	double t0 = dtime();

	for (n=first_leaf; n<last_leaf; n++) {
		int i;
		int idx = leaf[n].ipart;
		leaf[n].active = -1;	
		for (i=0; i<leaf[n].npart; i++) {		
			if (part[i + idx ].tag == 1) {
				leaf[n].active = 0;
				break;
			}
		}
	}

#ifdef MOREINFO
	printf(" [%d] LEN_TASK : %d,  MEM Buff per proc = 2 * %d MB \n", PROC_RANK, LEN_TASK, sizeof(int)*LEN_TASK/1024/1024);
#endif    
	double t1 = dtime();

	hostmalloc_task(LEN_TASK);

#ifdef INTRA_LB
	task_prepare_p2p_LB();

	double tt1 = dtime();

	p2ppack_compute_gpu(packlist, NPACK);
#else
	task_prepare_p2p();

	idxP2P += idxtask;

	double tt1 = dtime();

	ntask = idxtask;
//	printf(" ntask_p2p = %d\n", ntask);
	if (ntask > LEN_TASK) {
		printf( " error LEN_TASK ! proc_rank = %d ntask =%d LEN_TASK = %d\n", PROC_RANK, ntask, LEN_TASK);
		exit(0);
	}

	//	printf( " k - - - %d %d\n", task_s[0][0]-first_leaf, task_t[0][0]-first_leaf);

	//int n;
	//	if ( 1 == PROC_RANK)
	//	for (n=0; n<ntask; n++) {
	//		if (task_s[0][n]-first_leaf ==100 ){
	//			printf(" task_t[0][%d] = %d\n",n, task_t[0][n] - first_leaf);
	//		}
	//	}

	task_compute_p2p_gpu(ntask, task_t, task_s, 1);
#endif

	double t2 = dtime();
	idxtask = 0;

	walk_task_m2l(first_node, first_node);
	double t3 = dtime();

	idxM2L += idxtask;

	ntask = idxtask;

	task_compute_m2l(ntask);

	hostfree_task();

	walk_l2l(first_node);

	dtime_task += dtime() - t0;

#ifdef MOREINFO
	printf("[%d] fmm time = %lf, memory_init=%lf, task_prep=%lf, walk/compute_p2p=%lf, walk_m2l:%lf, compute_m2l/l2p=%lf\n", PROC_RANK, dtime() - t0, t1 - t0, tt1 - t1, t2 - t1, t3 - t2, dtime() - t3);
#endif

}


void fmm_task_act(int level)
{
	int n;
	double t0 = dtime();

	for (n=first_leaf; n<last_leaf; n++) {
		int i;
		int idx = leaf[n].ipart;
		leaf[n].active = -1;	
		for (i=0; i<leaf[n].npart; i++) {		
			if (part[i + idx ].tag == 1) {
				if (part[i+idx].act >= level) {
				leaf[n].active = level;
					break;
				}
			}
		}
	}

#ifdef MOREINFO
	printf(" [%d] LEN_TASK : %d,  MEM Buff per proc = 2 * %d MB \n", PROC_RANK, LEN_TASK, sizeof(int)*LEN_TASK/1024/1024);
#endif    
	double t1 = dtime();

	hostmalloc_task(LEN_TASK);

#ifdef INTRA_LB
	task_prepare_p2p_LB();

	double tt1 = dtime();

	p2ppack_compute_gpu(packlist, NPACK);
#else
	task_prepare_p2p_act(level);

	idxP2P += idxtask;

	double tt1 = dtime();

	ntask = idxtask;
//	printf(" ntask_p2p = %d\n", ntask);
	if (ntask > LEN_TASK) {
		printf( " error LEN_TASK ! proc_rank = %d ntask =%d LEN_TASK = %d\n", PROC_RANK, ntask, LEN_TASK);
		exit(0);
	}
	task_compute_p2p_gpu(ntask, task_t, task_s, 1);
#endif

	double t2 = dtime();
	idxtask = 0;

	walk_task_m2l(first_node, first_node);
	double t3 = dtime();

	idxM2L += idxtask;

	ntask = idxtask;

	task_compute_m2l(ntask);

	hostfree_task();

	walk_l2l(first_node);

	dtime_task += dtime() - t0;

#ifdef MOREINFO
	printf("[%d] fmm time = %lf, memory_init=%lf, task_prep=%lf, walk/compute_p2p=%lf, walk_m2l:%lf, compute_m2l/l2p=%lf\n", PROC_RANK, dtime() - t0, t1 - t0, tt1 - t1, t2 - t1, t3 - t2, dtime() - t3);
#endif

}



#else
////  CPU Version

void fmm_task() {
	if (print_msg())
		printf("Not implemented yet.\n");
}
void fmm_task_act(int active) {
	if (print_msg())
		printf("Not implemented yet.\n");
}

#endif

void fmm_deconstruct() 
{
	btree += first_node;
	myfree(btree, 41);

	leaf  += first_leaf;
	myfree(leaf, 40);
}

