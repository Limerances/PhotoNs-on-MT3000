#include "gpu_kernel.h"
#ifdef __CUDACC__

#define RETCODE cudaSuccess
#define devGetLastError cudaGetLastError 
#define devGetErrorString cudaGetErrorString
#define devHostGetDevicePointer cudaHostGetDevicePointer

#define devMalloc cudaMalloc
#define devFree cudaFree
#define devEvent_t cudaEvent_t
#define devEventCreate cudaEventCreate
#define devEventRecord cudaEventRecord
#define devEventElapsedTime cudaEventElapsedTime
#define devEventSynchronize cudaEventSynchronize
#define devEventDestroy cudaEventDestroy
#define devMemcpy cudaMemcpy
#define devMemcpyAsync cudaMemcpyAsync
#define devMemcpy2D cudaMemcpy2D
#define devMemcpy2DAsync cudaMemcpy2DAsync
#define devMemcpyHostToDevice cudaMemcpyHostToDevice
#define devMemcpyDeviceToHost cudaMemcpyDeviceToHost
#define devHostAlloc cudaHostAlloc
#define __WCB_FLAG__ cudaHostAllocWriteCombined
#define __ZCPY_FLAG__ cudaHostAllocMapped 
#define __ZCPY_FLAG_BAK__ cudaHostAllocWriteCombined|cudaHostAllocMapped 
#define __DEFAULT_FLAG__ cudaHostAllocDefault
#define devHostGetDevicePointer cudaHostGetDevicePointer
#define devHostFree cudaFreeHost
#define devMemset cudaMemset
#define devSetDevice cudaSetDevice
#define devDeviceSynchronize cudaDeviceSynchronize
#define devGetDeviceCount cudaGetDeviceCount

#else

#define RETCODE hipSuccess 
#define devGetLastError hipGetLastError 
#define devGetErrorString hipGetErrorString
#define devHostGetDevicePointer hipHostGetDevicePointer

#define devMalloc hipMalloc
#define devFree hipFree
#define devEvent_t hipEvent_t
#define devEventCreate hipEventCreate
#define devEventRecord hipEventRecord
#define devEventElapsedTime hipEventElapsedTime
#define devEventSynchronize hipEventSynchronize
#define devEventDestroy hipEventDestroy
#define devMemcpy hipMemcpy
#define devMemcpyAsync hipMemcpyAsync
#define devMemcpy2D hipMemcpy2D
#define devMemcpy2DAsync hipMemcpy2DAsync
#define devMemcpyHostToDevice hipMemcpyHostToDevice
#define devMemcpyDeviceToHost hipMemcpyDeviceToHost
#define devHostAlloc hipHostMalloc
#define __WCB_FLAG__ hipHostMallocWriteCombined
#define __ZCPY_FLAG__ hipHostMallocMapped 
#define __ZCPY_FLAG_BAK__ hipHostMallocWriteCombined|hipHostMallocMapped 
#define __DEFAULT_FLAG__ hipHostMallocDefault
#define devHostGetDevicePointer hipHostGetDevicePointer
#define devHostFree hipHostFree
#define devMemset hipMemset
#define devSetDevice hipSetDevice
#define devDeviceSynchronize hipDeviceSynchronize
#define devGetDeviceCount hipGetDeviceCount

#endif

#define checkDevErrors( a ) do { \
    if (RETCODE != (a)) { \
    fprintf(stderr, "Device runtime error in line %d of file %s \
    : %s \n", __LINE__, __FILE__, devGetErrorString(devGetLastError()) ); \
    exit(EXIT_FAILURE); \
    } \
} while(0);

#define checkKernelErrors(expr) do { \
    expr; \
    cudaError_t __err = devGetLastError(); \
    if (RETCODE != __err) { \
        fprintf(stderr, "Device Kernel error in line %d of file %s \
	: %s \n", __LINE__, __FILE__, devGetErrorString(__err)); \
    exit(EXIT_FAILURE); \
    } \
} while(0)

#ifndef HIP_RTM
#define NUMNODE 1
#endif

#define MERGE_PACK

int p2p_inited = 0;
int *max_npart, *max_nleaf, *max_ntask;
int *cur_npart, *cur_nleaf, *cur_ntask;
int *off_npart, *off_nleaf, *off_ntask;
float fac = 1.1, ifac = 0.9;
// kernel timing
devEvent_t* start_events, *stop_events;

void myhostmalloc(void** palloc, size_t sz, int idx, int flag=__DEFAULT_FLAG__) {
#ifdef GPUASYNC
	checkDevErrors(devHostAlloc((void**)palloc, sz, flag));
#ifdef MEMLOG
	malloc_update_log(sz, idx);
#endif
#else
	*palloc = mymalloc(sz, idx);
#endif
}

void myhostrealloc(void** palloc, size_t sz, int idx, int flag=__DEFAULT_FLAG__) {
#ifdef GPUASYNC
	printf("Not support for realloc host memory!\n");
	exit(0);
#else
	*palloc = myrealloc(*palloc, sz, idx);
#endif
}

void myhostfree(void* ptr, int idx) {
#ifdef GPUASYNC
	checkDevErrors(devHostFree(ptr));
#ifdef MEMLOG
	free_update_log(idx);
#endif
#else
	myfree(ptr, idx);
#endif
}

void sync_with_dev() {
	checkDevErrors(devDeviceSynchronize());
}

void hostmalloc_part_in(int NPART) {
#ifdef FASTTRANS
	checkDevErrors(devHostAlloc((void**)&h_part_pos, sizeof(gpuT) * NPART * 3, __DEFAULT_FLAG__ || __WCB_FLAG__));
#ifdef MEMLOG
	malloc_update_log(sizeof(gpuT) * NPART * 3, 52);
#endif
#else
	myhostmalloc((void**)&h_part_pos, sizeof(gpuT) * NPART * 3, 52);
#endif

#ifdef ACT_OPT
	myhostmalloc((void**)&h_part_act, sizeof(uint_8b) * NPART, 55);
#endif
}

void hostmalloc_part_out(int NPART) {
	myhostmalloc((void**)&h_part_buf, sizeof(gpuT) * NPART * 3, 54);
}

void hostfree_part_in() {
#ifdef FASTTRANS
	checkDevErrors(devHostFree(h_part_pos));
#ifdef MEMLOG
	free_update_log(52);
#endif
#else
	myhostfree(h_part_pos, 52);
#endif
#ifdef ACT_OPT
	myhostfree(h_part_act, 55);
#endif
}

void hostfree_part_out() {
	myhostfree(h_part_buf, 54);
}

void hostmalloc_leaf(int NLEAF) {
	myhostmalloc((void**)&h_leaf_ip ,sizeof(int) * NLEAF, 56);
	myhostmalloc((void**)&h_leaf_np ,sizeof(int) * NLEAF, 57);
}
 
void hostfree_leaf() {
	myhostfree(h_leaf_ip, 56);
	myhostfree(h_leaf_np, 57);
}

void hostmalloc_task(int LEN_TASK) {
	myhostmalloc((void**)&(task_s), sizeof(int) * LEN_TASK, 58);
	myhostmalloc((void**)&(task_t), sizeof(int) * LEN_TASK, 59);
}

void hostrealloc_task(int LEN_TASK) {
	myhostrealloc((void**)&(task_s), sizeof(int) * LEN_TASK, 58);
	myhostrealloc((void**)&(task_t), sizeof(int) * LEN_TASK, 59);
}

void hostfree_task() {
	myhostfree(task_s, 58);
	myhostfree(task_t, 59);
}

void hostmalloc_dinfo(int dev_cnt) {
	device_ids = (int*)mymalloc(sizeof(int) * dev_cnt, 51);
	max_npart = (int*)malloc(sizeof(int) * device_cnt);
	max_nleaf = (int*)malloc(sizeof(int) * device_cnt);
	max_ntask = (int*)malloc(sizeof(int) * device_cnt);
	cur_npart = (int*)malloc(sizeof(int) * device_cnt);
	cur_nleaf = (int*)malloc(sizeof(int) * device_cnt);
	cur_ntask = (int*)malloc(sizeof(int) * device_cnt);
	off_npart = (int*)malloc(sizeof(int) * MAXNPACK);
	off_nleaf = (int*)malloc(sizeof(int) * MAXNPACK);
	off_ntask = (int*)malloc(sizeof(int) * MAXNPACK);

	start_events = (devEvent_t*)malloc(sizeof(devEvent_t) * device_cnt);
	stop_events = (devEvent_t*)malloc(sizeof(devEvent_t) * device_cnt);

	d_keys.resize(device_cnt);	
	d_values.resize(device_cnt);	
	d_one.resize(device_cnt);

	d_part_pos.resize(device_cnt);	
	d_part_buf.resize(device_cnt);
	d_part_act.resize(device_cnt);
	
	d_leaf_ip.resize(device_cnt);	
	d_leaf_np.resize(device_cnt);
	d_inodes.resize(device_cnt);	
	d_in.resize(device_cnt);	
	d_nn.resize(device_cnt);	
}

void dinfofree() {
	myfree(device_ids, 51);
	free(max_npart);
	free(max_nleaf);
	free(max_ntask);
	free(cur_npart);
	free(cur_nleaf);
	free(cur_ntask);
	free(off_npart);
	free(off_nleaf);
	free(off_ntask);

	free(start_events);
	free(stop_events);

	d_keys.resize(0);	
	d_values.resize(0);	
	d_one.resize(0);

	d_part_pos.resize(0);	
	d_part_buf.resize(0);
	d_part_act.resize(0);
	
	d_leaf_ip.resize(0);	
	d_leaf_np.resize(0);
	d_inodes.resize(0);	
	d_in.resize(0);	
	d_nn.resize(0);	
}

void hostmalloc_in(int NPART, int NLEAF) {
	hostmalloc_part_in(NPART);
	hostmalloc_leaf(NLEAF);
}

void hostmalloc_out(int NPART) {
	hostmalloc_part_out(NPART);
}

void hostfree_in() {
	hostfree_part_in();
	hostfree_leaf();
}

void hostfree_out() {
	hostfree_part_out();
}

void devmalloc_task_core(int i, int len) {
	checkDevErrors(devSetDevice(device_ids[i]));
	d_keys[i] = thrust::device_vector<int>(len);
	d_values[i] = thrust::device_vector<int>(len);
	d_one[i] = thrust::device_vector<int>(len);
	thrust::fill(d_one[i].begin(), d_one[i].end(), 1);
}

void devfree_task_core(int i) {
	checkDevErrors(devSetDevice(device_ids[i]));
	d_keys[i].clear();
	d_keys[i].shrink_to_fit();
	d_values[i].clear();
	d_values[i].shrink_to_fit();
	d_one[i].clear();
	d_one[i].shrink_to_fit();
}

void devmalloc_task(int* lens) {
#pragma omp parallel for num_threads(device_cnt) 
    	for (int i=0; i<device_cnt; ++i) {
		if (lens[i] > max_ntask[i]) {
			max_ntask[i] = lens[i] * fac;
			devfree_task_core(i);
			devmalloc_task_core(i, max_ntask[i]);
		}
	}
}

void devfree_task() {
#pragma omp parallel for num_threads(device_cnt) 
    	for (int i=0; i<device_cnt; ++i) {
		devfree_task_core(i);
	}
}

void devmalloc_part_core(int i, int np) {
	checkDevErrors(devSetDevice(device_ids[i]));
        checkDevErrors(devMalloc((void**)&d_part_pos[i], sizeof(gpuT) * np * 3));
       	checkDevErrors(devMalloc((void**)&d_part_buf[i], sizeof(gpuT) * np * 3));
#ifdef ACT_OPT
        checkDevErrors(devMalloc((void**)&d_part_act[i], sizeof(uint_8b) * np));
#endif
}

void devfree_part_core(int i) {
	checkDevErrors(devSetDevice(device_ids[i]));
	checkDevErrors(devFree(d_part_pos[i]));
	checkDevErrors(devFree(d_part_buf[i]));
#ifdef ACT_OPT
	checkDevErrors(devFree(d_part_act[i]));
#endif

}

void devmalloc_part(int* nprts) {
#pragma omp parallel for num_threads(device_cnt) 
    	for (int i=0; i<device_cnt; ++i) {
		if (nprts[i] > max_npart[i]) {
			max_npart[i] = nprts[i] * fac;
			devfree_part_core(i);
			devmalloc_part_core(i, max_npart[i]);
		}
	}
}

void devfree_part() {
#pragma omp parallel for num_threads(device_cnt) 
    	for (int i=0; i<device_cnt; ++i) {
		devfree_part_core(i);
	}
}

void devmalloc_leaf_core(int i, int nleaf) {
	checkDevErrors(devSetDevice(device_ids[i]));
        checkDevErrors(devMalloc((void**)&d_leaf_ip[i], sizeof(int) * nleaf));
        checkDevErrors(devMalloc((void**)&d_leaf_np[i], sizeof(int) * nleaf));
	d_inodes[i] = thrust::device_vector<int>(nleaf);
	d_in[i] = thrust::device_vector<int>(nleaf);
	d_nn[i] = thrust::device_vector<int>(nleaf);
}

void devfree_leaf_core(int i) {
	checkDevErrors(devSetDevice(device_ids[i]));
	checkDevErrors(devFree(d_leaf_ip[i]));
	checkDevErrors(devFree(d_leaf_np[i]));
	d_inodes[i].clear();
	d_inodes[i].shrink_to_fit();
	d_in[i].clear();
	d_in[i].shrink_to_fit();
	d_nn[i].clear();
	d_nn[i].shrink_to_fit();

}

void devmalloc_leaf(int* nleaves) {
#pragma omp parallel for num_threads(device_cnt) 
    	for (int i=0; i<device_cnt; ++i) {
		if (nleaves[i] > max_nleaf[i]) {
			max_nleaf[i] = nleaves[i] * fac;
			devfree_leaf_core(i);
			devmalloc_leaf_core(i, max_nleaf[i]);
		}
	}
}

void devfree_leaf() {
#pragma omp parallel for num_threads(device_cnt) 
    	for (int i=0; i<device_cnt; ++i) {
		devfree_leaf_core(i);
	}
}

void devmalloc(int* lens, int* nprts, int* nleaves) {
	devmalloc_task(lens);
	devmalloc_part(nprts);
	devmalloc_leaf(nleaves);
}

void devfree() {
	devfree_task();
	devfree_part();
	devfree_leaf();
}

int get_new_rank(int proc_rank) {
	int new_rank;
#ifdef HIP_RTM
	int py = proc_rank / NSIDE1PROC ; 
	int pz = proc_rank % NSIDE1PROC ;
 
	int dside0 = NSIDE0PROC/NSIDE0SUDOM;
	int dside1 = NSIDE1PROC/NSIDE1SUDOM;

	int sy = py / dside0;
	int sz = pz / dside1;

	int sudom_r  = sy * NSIDE1SUDOM + sz;

	int iy = py -  sy*dside0 ;
	int iz = pz -  sz*dside1 ;

	int sudom_size = dside0 * dside1;

	new_rank = sudom_r * sudom_size + iy * dside1 + iz;
	return new_rank;
#else
	return proc_rank;
#endif
}

void initGPU(int proc_size, int proc_rank) {
	checkDevErrors(devGetDeviceCount(&device_cnt));
	if (device_cnt <= 0) {
		printf("No available device\n");
		exit(0);
	}
	NEW_RANK = get_new_rank(proc_rank);
#ifdef INTRA_LB
	LOG(3, "[%d] There is %d available devices\n", proc_rank, device_cnt);
	hostmalloc_dinfo(device_cnt);
	for (int i=0; i<device_cnt; ++i) {
		device_ids[i] = (NEW_RANK % device_cnt + i) % device_cnt;
	}
#else
	hostmalloc_dinfo(1);
	if (proc_size / device_cnt == 0)
		device_ids[0] = NEW_RANK;
	else {
#ifdef HIP_RTM
		device_ids[0] = 0;
#else
		device_ids[0] = (NEW_RANK % (proc_size / NUMNODE)) / (proc_size/ NUMNODE / device_cnt);
#endif
	}
	device_cnt = 1;
#endif
}

struct cmp
{
	bool operator()(std::pair<int, int> a, std::pair<int, int> b)
	{
		if (a.second > b.second)
			return true;
		return false;
	}
};

int compare_p2ppack_ntask(const void *p1, const void *p2) {
	if ( ((P2PPack*)p1)->ntask < ((P2PPack*)p2)->ntask )
		return 1;
	else
		return 0;
}

int compare_p2ppack_ipart(const void *p1, const void *p2) {
	if ( ((P2PPack*)p1)->rankid > ((P2PPack*)p2)->rankid )
		return 1;
	else if ( ((P2PPack*)p1)->rankid == ((P2PPack*)p2)->rankid ) {
		if ( ((P2PPack*)p1)->ipart > ((P2PPack*)p2)->ipart )
			return 1;
		else if ( ((P2PPack*)p1)->ipart == ((P2PPack*)p2)->ipart) {
			if ( ((P2PPack*)p1)->ntask < ((P2PPack*)p2)->ntask )
				return 1;
		}
	}
	return 0;
}
			
void assign_pack(P2PPack* pack, P2PPack* p) {
	if (pack == p)
		return;
	pack->ntask = p->ntask;
	pack->task_t = p->task_t;
	pack->task_s = p->task_s;
	pack->ileaf = p->ileaf;
	pack->nleaf = p->nleaf;
	pack->ipart = p->ipart;
	pack->npart = p->npart;
	pack->rankid = p->rankid;

	p->ntask = 0;
	p->task_t = NULL;
	p->task_s = NULL;
}

void schedule_on_devices(P2PPack* packs, int& npack, int dev_cnt) {
	std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, cmp> minque;
	std::pair<int, int> tmp;

	// make sure the wholepack being located on device 0
	if (packs[npack-1].rankid == 0) {
		minque.push(std::make_pair(0, packs[npack-1].ntask));
	}
	else {
		minque.push(std::make_pair(0, 0));
	}
	for (int i=1; i<dev_cnt; ++i)
		minque.push(std::make_pair(i, 0));

	// assign the packs under load banlance with greedy algo
	qsort(packs, npack, sizeof(P2PPack), compare_p2ppack_ntask);
	for (int i=0; i<npack; ++i) {
		if (packs[i].rankid != -1)
			continue;
		tmp = minque.top();
		minque.pop();
		packs[i].rankid = tmp.first;
		tmp.second += packs[i].ntask;
		minque.push(tmp);	
	}

	// merge overlapping packs on the same device to reduce the transfer overhead
#ifdef MERGE_PACK
	qsort(packs, npack, sizeof(P2PPack), compare_p2ppack_ipart);
	int tmp_ntask, id_task, maxleaf;
	int npck = 0, tpos = 0, spos = 1;
	while (tpos < npack) {
		if (spos >= npack 
			|| packs[spos].rankid != packs[tpos].rankid 
			|| packs[spos].ipart > packs[tpos].ipart + packs[tpos].npart) {
			// can not merge, update npck, tpos and spos
			if (spos - tpos > 1) {
				tmp_ntask = 0;
				for (int i=tpos; i<spos; ++i) {
					tmp_ntask += packs[i].ntask;
				}
				packs[tpos].task_t = (int*)realloc(packs[tpos].task_t, sizeof(int) * tmp_ntask);
				packs[tpos].task_s = (int*)realloc(packs[tpos].task_s, sizeof(int) * tmp_ntask);
				id_task = packs[tpos].ntask;
				for (int i=tpos+1; i<spos; ++i) {
					for (int j=0; j<packs[i].ntask; ++j) {
						packs[tpos].task_t[id_task] = packs[i].task_t[j];
						packs[tpos].task_s[id_task] = packs[i].task_s[j];
						id_task++;
					}
					free(packs[i].task_t);
					free(packs[i].task_s);
					packs[i].ntask = 0;
					packs[i].task_t = NULL;
					packs[i].task_s = NULL;
				}
				packs[tpos].ntask = tmp_ntask;
			}

			assign_pack(packs + npck, packs + tpos);
			npck ++;
			tpos = spos;
			spos = tpos + 1; 
		} else {
			// merge, update npart and nleaf
			maxleaf = max(packs[tpos].ileaf + packs[tpos].nleaf, 
						packs[spos].ileaf + packs[spos].nleaf);
			packs[tpos].nleaf = maxleaf - packs[tpos].ileaf;
			if (maxleaf < last_leaf) {
				packs[tpos].npart = leaf[maxleaf].ipart - packs[tpos].ipart;
			} else {	
				packs[tpos].npart = NPART - packs[tpos].ipart;
			}

			spos++;
		}
	}
	npack = npck;
	NPACK = npack;
	for (int i=0; i<NPACK; ++i)
		LOG(3, "[%d] after merging PACK(%d) ipart %d-%d ileaf %d-%d p2p %d tast_t %lld task_s %lld on g%d\n", PROC_RANK, i, packs[i].ipart, packs[i].ipart+packs[i].npart, packs[i].ileaf, packs[i].ileaf+packs[i].nleaf, packs[i].ntask, packs[i].task_t, packs[i].task_s, packs[i].rankid);
#endif	
}	

void devBuffPrep(P2PPack* packs, int& npack) {

	if (!p2p_inited) {
		initGPU(PROC_SIZE, PROC_RANK);
		for (int i=0; i<device_cnt; ++i) {
			max_npart[i] = 0;
			max_nleaf[i] = 0;
			max_ntask[i] = 0;
		}
		p2p_inited = 1;
	}

	for (int i=0; i<device_cnt; ++i) {
		cur_npart[i] = 0;
		cur_nleaf[i] = 0;
		cur_ntask[i] = 0;
	}
#ifdef INTRA_LB	
	schedule_on_devices(packs, npack, device_cnt);
#endif
	for (int i=0; i<npack; ++i) {
		off_npart[i] = cur_npart[packs[i].rankid];
		off_nleaf[i] = cur_nleaf[packs[i].rankid];
		off_ntask[i] = cur_ntask[packs[i].rankid];
		cur_npart[packs[i].rankid] += packs[i].npart;
		cur_nleaf[packs[i].rankid] += packs[i].nleaf;
		cur_ntask[packs[i].rankid] += packs[i].ntask;
	}

	devmalloc(cur_ntask, cur_npart, cur_nleaf);
}

void mydevMemcpyDH(void* tptr, void* sptr, size_t sz) {
//#ifdef GPUASYNC
//	checkDevErrors(devMemcpyAsync(tptr, sptr, sz, devMemcpyDeviceToHost));
//#else
	checkDevErrors(devMemcpy(tptr, sptr, sz, devMemcpyDeviceToHost));
//#endif
}

void mydevMemcpyHD(void* tptr, void* sptr, size_t sz) {
#ifdef GPUASYNC
	checkDevErrors(devMemcpyAsync(tptr, sptr, sz, devMemcpyHostToDevice));
#else
	checkDevErrors(devMemcpy(tptr, sptr, sz, devMemcpyHostToDevice));
#endif
}

void mydevMemcpy2DHD(void* tptr, size_t dpitch, void* sptr, size_t spitch, size_t width, size_t height) {
#ifdef GPUASYNC
	checkDevErrors(devMemcpy2DAsync(tptr, dpitch, sptr, spitch, width, height, devMemcpyHostToDevice));
#else
	checkDevErrors(devMemcpy2D(tptr, dpitch, sptr, spitch, width, height, devMemcpyHostToDevice));
#endif	
}

void updateAcc_fromHostBuff(P2PPack pack, int active_level=-1) {
	for (int i=0; i<pack.npart; ++i) {
		int j = pack.ipart + i;
		if (active_level == -1 || (active_level && part[j].act == active_level)) {
			acc[0][j] += h_part_buf[i];
	    		acc[1][j] += h_part_buf[i+pack.npart];
	    		acc[2][j] += h_part_buf[i+2*pack.npart];
		}
	}
}

void devMemcpyDH_updateAcc(P2PPack* packs, int npack, int active_level=-1) {
	for (int n=0; n<npack; ++n) {
		int i = packs[n].rankid;
		if (packs[n].ntask == 0)
			continue;
		checkDevErrors(devSetDevice(device_ids[i]));
		mydevMemcpyDH(h_part_buf, d_part_buf[i] + off_npart[n]*3, sizeof(gpuT) * packs[n].npart * 3);
		updateAcc_fromHostBuff(packs[n], active_level);
	}
#ifdef MOREINFO
	for (int j=0; j<NPART; j++) {
		if (abs(acc[0][j]) > 1.0e-6) {
			printf("acc(%d) %e %e %e\n", j, acc[0][j], acc[1][j], acc[2][j]);
			break;
		}
	}
/*
	char fname[64];
	sprintf(fname, "./lb_%d", packs[0].ntask);
	FILE* ftest = fopen(fname, "w");
	for (int j=0; j<NPART; j+=100) {
		fprintf(ftest, "acc(%d) %e %e %e\n", j, acc[0][j], acc[1][j], acc[2][j]);
	}
	fclose(ftest);
*/
#endif
}

void devMemcpyHD_part(P2PPack* packs, int pid=0) {
	P2PPack pack = packs[pid];
	int i = pack.rankid;
#ifdef INTRA_LB
	mydevMemcpy2DHD(d_part_pos[i] + off_npart[pid]*3, pack.npart * sizeof(gpuT), h_part_pos + pack.ipart, NPART * sizeof(gpuT), pack.npart * sizeof(gpuT), 3);
/*	mydevMemcpyHD(d_part_pos[i] + off_npart[pid]*3, h_part_pos + pack.ipart, sizeof(gpuT) * pack.npart);
	mydevMemcpyHD(d_part_pos[i] + off_npart[pid]*3 + pack.npart, h_part_pos + NPART + pack.ipart, sizeof(gpuT) * pack.npart);
	mydevMemcpyHD(d_part_pos[i] + off_npart[pid]*3 + 2*pack.npart, h_part_pos + 2*NPART + pack.ipart, sizeof(gpuT) * pack.npart);
*/
#else
	mydevMemcpyHD(d_part_pos[i] + off_npart[pid]*3, h_part_pos, sizeof(gpuT) * NPART * 3);
#endif
	checkDevErrors(devMemset(d_part_buf[i] + off_npart[pid]*3, 0, sizeof(gpuT) * pack.npart * 3));
#ifdef ACT_OPT
	mydevMemcpyHD(d_part_act[i] + off_npart[pid], h_part_act + pack.ipart, sizeof(uint_8b) * pack.npart);
#endif
	mydevMemcpyHD(d_leaf_ip[i] + off_nleaf[pid], h_leaf_ip + pack.ileaf - first_leaf, sizeof(int) * pack.nleaf);
	mydevMemcpyHD(d_leaf_np[i] + off_nleaf[pid], h_leaf_np + pack.ileaf - first_leaf, sizeof(int) * pack.nleaf);

}

void devMemcpyHD_task(P2PPack* packs, int pid=0) {
	P2PPack pack = packs[pid];
	int i = pack.rankid;
	for (int j=0; j<pack.ntask; ++j) {
		pack.task_t[j] -= pack.ileaf;
		pack.task_s[j] -= pack.ileaf;
	}
	mydevMemcpyHD(thrust::raw_pointer_cast(d_keys[i].data() + off_ntask[pid]), pack.task_t, sizeof(int) * pack.ntask);
	mydevMemcpyHD(thrust::raw_pointer_cast(d_values[i].data() + off_ntask[pid]), pack.task_s, sizeof(int) * pack.ntask);
}

double get_last_kntime(int proc_rank) {
	double elapsed_time = 0.;

	for (int i=0; i<device_cnt; ++i) {
		checkDevErrors(devSetDevice(device_ids[i]));
	        checkDevErrors(devEventSynchronize(stop_events[i]));
	        float time = 0;
	        checkDevErrors(devEventElapsedTime(&time, start_events[i], stop_events[i]));
	        LOG(3, "[%d] kernel time on dev %d: %lf ms\n", proc_rank, i, time);
		elapsed_time = max(elapsed_time, time);
		checkDevErrors ( devEventDestroy(start_events[i]) );
		checkDevErrors ( devEventDestroy(stop_events[i]) );
	}

	return elapsed_time / 1.0e3;
}
