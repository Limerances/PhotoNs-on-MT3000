#ifndef MTKDTREE_H
#define MTKDTREE_H

#include <stdio.h>
#include <stdlib.h>
#include <queue>
#include <atomic>
#include <thread>
#include <omp.h>
#include <math.h>
#include "utility.h"

using namespace std;

#define NDIM 3
#define NSON 2

typedef struct {
	int npart;
	int ipart;
	MyFloat width[3];
	MyFloat center[3];
} tPack;

typedef struct {
	int  npart;
	MyFloat width[3];
	MyFloat center[3];
	int  son[NSON];
	MyFloat split;
} tNode;

template<typename Body_T, typename Node_T, typename Pack_T>
class PartTree {
public:
	Body_T* part;
	int NPART;
	Node_T* btree;
	int first_node, last_node;
	Pack_T* leaf;
	int first_leaf, last_leaf;
	PartTree():part(nullptr), NPART(0), 
		btree(nullptr), first_node(0), last_node(0), 
		leaf(nullptr), first_leaf(0), last_leaf(0) {};
	PartTree(Body_T *prt, int NP, Node_T* bt, int fnode, int lnode, 
		Pack_T *lf, int fleaf, int lleaf):part(prt), NPART(NP), 
		btree(bt), first_node(fnode), last_node(lnode), 
		leaf(lf), first_leaf(fleaf), last_leaf(lleaf) {};
	~PartTree() {}
};

class PartIndex {
public:
	int matched_num;
	vector<int> matched_list;
	vector<MyFloat> distance_list;
	PartIndex() {
		matched_num = 0;
		matched_list.resize(0);
		distance_list.resize(0);
	}
	~PartIndex() {
		matched_num = 0;
		matched_list.resize(0);
		distance_list.resize(0);
	}
};


struct cmp
{
	bool operator()(pair<int, MyFloat> left, pair<int, MyFloat> right)
	{
		if (left.second < right.second)
			return true;
		return false;
	}
};

template<typename Body_T, typename Node_T, typename Pack_T>
class MTKDtree {
private:
	atomic<int> num_threads;
	atomic<int> atm_lleaf, atm_lnode;
	int MAXinLEAF;
	int MAXNNode;
	int MAXNLeaf;
	
public:
	Body_T* tr_part;
	int tr_NPART;
	Node_T* btree;
	int fnode, lnode;
	Pack_T* leaf;
	int fleaf, lleaf;

	MTKDtree():tr_part(nullptr), tr_NPART(0), 
		btree(nullptr), fnode(0), lnode(0), 
		leaf(nullptr), fleaf(0), lleaf(0) {};
	MTKDtree(Body_T* prt, int NP, int maxleaf) {
		tr_part = prt;
		tr_NPART = NP;
		MAXinLEAF = maxleaf;

		MAXNNode = max(16, (int)ceil(2.0*((double)tr_NPART)/((double)MAXinLEAF)));
		MAXNLeaf = max(32, (int)ceil(2.0*((double)tr_NPART)/((double)MAXinLEAF)));
	
		if (MAXNNode > tr_NPART)
			MAXNNode = tr_NPART +1 ;
		if (MAXNLeaf > tr_NPART)
			MAXNLeaf = tr_NPART +1;
	
		fleaf = lleaf = tr_NPART;
		fnode = lnode = tr_NPART+MAXNLeaf;
	
		// 0...NPART-1, first_leaf...last_leaf-1, first_node...last_node
	
		leaf  = (Pack_T*)mymalloc(sizeof(Pack_T)*MAXNLeaf, 250);
		LOG(3, "leaf malloc (%lld, %lf GB) ", MAXNLeaf, (MyFloat)MAXNLeaf * sizeof(Pack_T) / 1.0e9);
		if (leaf == NULL ){
			LOG(0, "error leaf malloc");
			exit(0);
		}
		btree = (Node_T*)mymalloc(sizeof(Node_T)*MAXNNode, 251);
		LOG(3, "btree malloc (%lld, %lf GB) ", MAXNNode, (MyFloat)MAXNNode * sizeof(Node_T) / 1.0e9);
		if (btree == NULL ){
			LOG(0, "error btree malloc");
			exit(0);
		}
		int n, m;
		for (n=0; n<MAXNNode; n++) {
			btree[n].npart = 0;
			btree[n].son[0] = btree[n].son[1] = -1;
			btree[n].split = 0.0;
		}
	
	
		for (n=0; n<MAXNLeaf; n++) {
			leaf[n].npart = 0;
			leaf[n].ipart = 0;
		}
	
		btree -= fnode;
		leaf  -= fleaf;
	};

	~MTKDtree() {}

	void MTKDtree_free() {
		leaf += fleaf;
		if (leaf)
			myfree(leaf, 250);
		btree += fnode;
		if (btree)
			myfree(btree, 251);

	}
	
	void kdSplit(int direct, int iPart, int length, int iNode) {
		if (length == 0)
			return;
	
	
		btree[iNode].npart = length;
		Body_T *p = tr_part + iPart;
	
		int npart[2];
		double split;
	        // done by current single-thread
		bksort_inplace(direct, p, length, npart, &split);
	
		btree[iNode].split = split;
	
		int n, ip;
		ip = iPart;
		int new_direct = (direct + 1)%NDIM;
	
		thread* td_list[NSON-1];
		int td_cnt = 0;
		int leaf_id, node_id;
		for (n=0; n<NSON; n++) {
	
			if (npart[n]  <= MAXinLEAF) {
				leaf_id = (atm_lleaf++);
				leaf[leaf_id].npart = npart[n];
				leaf[leaf_id].ipart = ip;
	
				btree[iNode].son[n] = leaf_id;
			}
			else {
				node_id = (++atm_lnode);
				btree[iNode].son[n] = node_id;
				if (n != NSON - 1 && (--num_threads) >= 1) {
					td_list[td_cnt++] = new thread(&MTKDtree::kdSplit, this, new_direct, ip, npart[n], node_id);
				} else {
					kdSplit(new_direct, ip, npart[n],  node_id);
				}
			}
	
			ip += npart[n];
		}
		for (n=0; n<td_cnt; ++n)
			td_list[n]->join();
	}

	PartTree<Body_T, Node_T, Pack_T> build_kdtree(int dir=0) {
		num_threads = min((tr_NPART+8*MAXinLEAF-1)/(8*MAXinLEAF), omp_get_max_threads());
	
		atm_lleaf = fleaf;
		atm_lnode = fnode;
	
		kdSplit(dir, 0, tr_NPART, fnode);
	
		lleaf = atm_lleaf;
		lnode = atm_lnode;
		if (lnode - fnode + 1 > MAXNNode) {
			printf("MAXNNode (p %d n %d)is not large enough (< %d)! \n", tr_NPART, MAXNNode, lnode - fnode + 1);
			exit(0);
		}
		if (lleaf - fleaf > MAXNLeaf) {
			printf("MAXNLeaf (p %d l %d)is not large enough (< %d)! \n", tr_NPART, MAXNLeaf, lleaf - fleaf);
			exit(0);
		}
		return PartTree<Body_T, Node_T, Pack_T>(tr_part, tr_NPART, 
					btree, fnode, lnode,
					leaf, fleaf, lleaf);
	}

	void center_kdtree(int direct, int iNode, MyFloat left[3], MyFloat right[3]) {
		int n, idx;
		btree[iNode].width[0] = right[0]-left[0];
		btree[iNode].width[1] = right[1]-left[1];
		btree[iNode].width[2] = right[2]-left[2];
		btree[iNode].center[0] = 0.5*(right[0]+left[0]);
		btree[iNode].center[1] = 0.5*(right[1]+left[1]);
		btree[iNode].center[2] = 0.5*(right[2]+left[2]);
	
		int newd = (direct + 1)%3;
		MyFloat tmp;
		for (n=0; n<NSON; n++) {
			idx = btree[iNode].son[n];
	
			if (idx < lleaf ) {
	
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
};

template<typename Body_T>
void bksort_inplace(int D, Body_T *p, int length, int npart[2], double *split)
{
	npart[0] = npart[1] = 0;
	if (length <2 ) {
		npart[1] = length;
		return;
	}
	Body_T temp;
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
	mean *= mul;

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


template<typename Body_T, typename Node_T, typename Pack_T>
class KNN3DIndex {
private:
	priority_queue<pair<int, MyFloat>, vector<pair<int, MyFloat>>, cmp> maxque;
	double max_length;
	int min_NN;
public:
	Body_T* knn_part;
	int knn_NPART;
	Node_T* btree;
	int fnode, lnode;
	Pack_T* leaf;
	int fleaf, lleaf;
	
	KNN3DIndex() {};
	KNN3DIndex(PartTree<Body_T, Node_T, Pack_T>* prttree) {
		knn_part = prttree->part;
		knn_NPART = prttree->NPART;
		btree = prttree->btree;
		fnode = prttree->first_node;
		lnode = prttree->last_node;
		leaf = prttree->leaf;
		fleaf = prttree->first_leaf;
		lleaf = prttree->last_leaf;
	}
	KNN3DIndex(Body_T* prt, int np) {
		knn_part = prt;
		knn_NPART = np;
	}
	
	double distance(Body_T& ipart, int p) { 
		double dx, dr2 = 0;
		for (int d=0; d<NDIM; d++) {
			dx = knn_part[p].pos[d] - ipart.pos[d]; 
			dr2 += dx * dx;
		}
		return sqrt(dr2);
	}

	// min_N: at least minN Nearest Neighbours
	// max_r: max search radius
	PartIndex* KNN_search(Body_T *pos, int np, int min_N, double max_r) {

		PartIndex *pindex = new PartIndex[np];
		min_NN = min_N;
		for (int i=0; i<np; ++i) {
			max_length = max_r;
			que_clear();
			SEARCH(0, fnode, pos[i]);
			pindex[i].matched_num = maxque.size();
			pindex[i].matched_list.resize(pindex[i].matched_num);
			pindex[i].distance_list.resize(pindex[i].matched_num);
			for (int j=pindex[i].matched_num-1; j>=0; --j) {
				pindex[i].matched_list[j] = maxque.top().first;
				pindex[i].distance_list[j] = maxque.top().second;
				maxque.pop();
			}
		}
		return pindex;
	}

	PartIndex* KNN_search_traverse(Body_T *pos, int np, int min_N, double max_r) {
		PartIndex *pindex = new PartIndex[np];
		double dr;
		for (int i=0; i<np; ++i) {
			max_length = max_r;
			que_clear();
			for (int n=0; n<knn_NPART; n++) {
				if (pos[i].pidx == knn_part[n].pidx)
					continue;
		
				dr = distance(pos[i], n);
		
				if (dr <= max_length) {
					// inser into sHalos
					maxque.push(make_pair(n, dr));
					if (maxque.size() > min_NN) {
						maxque.pop();
						if (max_length > maxque.top().second) {
							max_length = maxque.top().second;
						}
					}
				}
			}
			pindex[i].matched_num = maxque.size();
			pindex[i].matched_list.resize(pindex[i].matched_num);
			pindex[i].distance_list.resize(pindex[i].matched_num);
			for (int j=pindex[i].matched_num-1; j>=0; --j) {
				pindex[i].matched_list[j] = maxque.top().first;
				pindex[i].distance_list[j] = maxque.top().second;
				maxque.pop();
			}
		}
		return pindex;
	}

	void insert_id_KNNque(Body_T& ipart, int ileaf) {
		if (ileaf < fleaf || ileaf >= fnode ) {
			LOG(0, "insert KNNque error !");
			exit(0);
		}
	
		int n, p;
		double dr;
	
		for (n=0; n<leaf[ileaf].npart; n++) {
			p = leaf[ileaf].ipart + n;		
	
			if (ipart.pidx == knn_part[p].pidx)
				continue;
	
			dr = distance(ipart, p);
	
			if (dr <= max_length) {
				// inser into sHalos
				maxque.push(make_pair(p, dr));
				if (maxque.size() > min_NN) {
					maxque.pop();
					if (max_length > maxque.top().second) {
						max_length = maxque.top().second;
					}
				}
			}
		}
	
	}

	void SEARCH(int direct, int iNode, Body_T& ipart) {
		int new_direct = (direct + 1)%NDIM;
		int in, out, in0, in1;
		double ds;
	
		// iNode is a leaf
		if ( iNode < fnode ) {
			insert_id_KNNque(ipart, iNode);	
			return ;
		}
	
		in0 = btree[iNode].son[0];
		in1 = btree[iNode].son[1];
	
		ds = btree[iNode].split - ipart.pos[direct]; 
	
		if ( ds >= 0.0 ) {
			in = in0;
			out = in1;
		}
		else {
			in = in1;
			out = in0;
		}
	
		if ( in >= fnode ) {
			// recursive search
			SEARCH(new_direct, in, ipart);
		}
		else {
			// in is a leaf
			insert_id_KNNque(ipart, in);
		}
	
		if (  out < fnode ) {
			// out is a leaf
			insert_id_KNNque(ipart, out);	
		}
		else if ( fabs(ds) < max_length ) {
			// recursive search
			SEARCH(new_direct, out, ipart);
		}
	
	}

	void que_clear() {
		if (maxque.size() != 0) {
			printf("not empty\n");
		}
//		priority_queue<pair<int, MyFloat>, vector<pair<int, MyFloat>>, cmp> empty;
//		swap(empty, maxque);
		while(!maxque.empty()) {
			maxque.pop();
		}
		if (maxque.size() != 0) {
			printf("clear fail\n");
		}
	}
	
};


#endif
