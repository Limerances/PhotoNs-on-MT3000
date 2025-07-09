#include "photoNs.h"
#include <atomic>
#include <thread>
#include <omp.h>

using namespace std;

atomic<int> num_threads;
atomic<int> atm_last_leaf, atm_last_node;

extern "C" void bksort_inplace(int D, Body *p, int length, int npart[2], double *split);

void kdSplit(int direct, int iPart, int length, int iNode) {
	if (length == 0)
		return;


	btree[iNode].npart = length;
	Body *p = part + iPart;

	int npart[2];
	double split;
        // done by current single-thread
	bksort_inplace(direct, p, length, npart, &split);

	btree[iNode].split = split;

	int n, ip;
	ip = iPart;
	int new_direct = (direct + 1)%3;

	thread* td_list[NSON-1];
	int td_cnt = 0;
	int leaf_id, node_id;
	for (n=0; n<NSON; n++) {

		if (npart[n]  <= MAXLEAF) {
			leaf_id = (atm_last_leaf++);
			leaf[leaf_id].npart = npart[n];
			leaf[leaf_id].ipart = ip;

			btree[iNode].son[n] = leaf_id;
		}
		else {
			node_id = (++atm_last_node);
			btree[iNode].son[n] = node_id;
			if (n != NSON - 1 && (--num_threads) >= 1) {
				td_list[td_cnt++] = new thread(kdSplit, new_direct, ip, npart[n], node_id);
			} else {
				kdSplit(new_direct, ip, npart[n],  node_id);
			}
		}

		ip += npart[n];
	}
	for (n=0; n<td_cnt; ++n)
		td_list[n]->join();
}

extern "C" void build_kdtree_threads() {
	num_threads = omp_get_max_threads();

	atm_last_leaf = first_leaf;
	atm_last_node = first_node;

	kdSplit(0, 0, NPART, first_node);

	last_leaf = atm_last_leaf;
	last_node = atm_last_node;

	if (last_node - first_node + 1 > MAXNNODE) {
		printf("[%d] MAXNNODE (p %d n %d)is not large enough (< %d)! \n", PROC_RANK, NPART, MAXNNODE, last_node - first_node + 1);
		exit(0);
	}
	if (last_leaf - first_leaf > MAXNLEAF) {
		printf("[%d] MAXNLEAF (p %d l %d)is not large enough (< %d)! \n", PROC_RANK, NPART, MAXNLEAF, last_leaf - first_leaf);
		exit(0);
	}
}
