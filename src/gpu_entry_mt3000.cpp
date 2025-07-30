#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <queue>
#include <sys/time.h>
#include "utility.h"

#ifdef __CUDACC__
#include <cuda_runtime.h>
#else
#include <hip/hip_runtime.h>
#include <cuda_wrappers/algorithm>
#include <cuda_wrappers/new>
#endif

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <thrust/generate.h>
#include <thrust/detail/type_traits.h>

extern "C" {
    #include "photoNs.h"
    #include "gpu_kernel.h" 
}
    #include "gpu_thrust.h" 

#include "gpu_kernel.impl.h"
#include "gpu_func.impl_mt3000.h"

struct P2P_kernel_Params
{
    gpuT SoftenScale;
    gpuT MASSPART;
    gpuT rrs,;
    gpuT coeff;
    int maxleaf;
    int numInodes;
    int *dev_inodes;
    int *dev_in;
    int *dev_nn;
    int *dev_jnodes;
    int ipart;
    int npart;
    gpuT *dev_part_pos;
    gpuT *dev_part_acc;
    int *dev_leaf_ip;
    int *dev_leaf_np;
    uint_8b *dev_part_act = NULL;
    uint_8b active_level = 0;
};


extern "C" void *task_compute_p2p_gpu(int ntask,
			int* ltask_t, int* ltask_s, int last_time) {
	double t0, t1;
	static int tot_tasks = 0, call_cnt = 0;
	bool first_time = (call_cnt == 0);
	tot_tasks += ntask;
	call_cnt ++;

	P2PPack* pack = packlist;
	pack->ntask = ntask;
	pack->task_t = ltask_t;
	pack->task_s = ltask_s;
	pack->ileaf = first_leaf;
	pack->nleaf = last_leaf - first_leaf;
	pack->ipart = 0;
	pack->npart = NPART;
	pack->rankid = 0;

#ifdef MOREINFO
	long long flopcnt = 0;
	int inode, jnode;
	for(int i = 0; i < pack->ntask; i++ ) 
	{
	    inode = pack->task_t[i];
	    jnode = pack->task_s[i];
	
	    flopcnt += leaf[inode].npart * leaf[jnode].npart;
//	    flopcnt += min(leaf[inode].npart, 128) * min(leaf[jnode].npart, 128);
	}
	LOG(0, "[%d] flop count: %lf (G)\n", PROC_RANK, flopcnt / 1.0e9);
#endif

	double ttest0, ttest1, ttest2, ttest3, ttest4;
	ttest0 = dtime();
	if (first_time) NPACK++;
	devBuffPrep(pack, NPACK);
	ttest1 = dtime();
	
	int nleaf = last_leaf - first_leaf;
	if (pack->nleaf != nleaf || pack->npart != NPART) {
		printf("Pack info error!\n");
		exit(0);
	}

	int rankid = pack->rankid;
	checkDevErrors(devSetDevice(device_ids[rankid]));

    checkDevErrors(hthread_dev_open(current_cluster_id));
    checkDevErrors(hthread_dat_load(current_cluster_id,"../obj/test.dat"));

	gpuT rrs = 1.0 / splitRadius;
	gpuT coeff = 2.0 * PIisq;

	if (first_time) {
		hostmalloc_in(NPART, nleaf);
		t1 = dtime();
		for (int i=0; i<NPART; ++i) {
//			h_part_pos[i] = part[i].pos[0];
//			h_part_pos[i+NPART] = part[i].pos[1];
//			h_part_pos[i+2*NPART] = part[i].pos[2];
#ifndef INTXYZ
		h_part_pos[i] = part[i].pos[0];
		h_part_pos[i+NPART] = part[i].pos[1];
		h_part_pos[i+2*NPART] = part[i].pos[2];
#else
		h_part_pos[i] = (float) (INT2POS * (part[i].posi[0] - domain_center_i[0]) );
		h_part_pos[i+NPART] = (float)(INT2POS * (part[i].posi[1] - domain_center_i[1]) );
		h_part_pos[i+2*NPART] = (float)(INT2POS * (part[i].posi[2] - domain_center_i[2]) );
#endif







		}
		LOG(3, "[%d] part (NPART %d) reorder time: %lf s\n", PROC_RANK, NPART, dtime()-t1);
		t1 = dtime();
		for (int i=0; i<nleaf; ++i) {
			h_leaf_ip[i] = leaf[i+first_leaf].ipart;
			h_leaf_np[i] = leaf[i+first_leaf].npart;
		}
		LOG(3, "[%d] leaf reorder time: %lf s\n", PROC_RANK, dtime()-t1);
		t1 = dtime();
		devMemcpyHD_part(pack);
		LOG(3, "[%d] copy input time: %lf s\n", PROC_RANK, dtime()-t1);
		hostfree_in();

        	checkDevErrors( devEventCreate(&start_events[rankid]) );
        	checkDevErrors( devEventCreate(&stop_events[rankid]) );
        	checkDevErrors( devEventRecord(start_events[rankid], 0));
        	checkDevErrors( devEventRecord(stop_events[rankid], 0));
	}
	
	if (pack->ntask > 0) {
		t0 = dtime();

		devMemcpyHD_task(pack);
		LOG(3, "[%d] copy task list time: %lf s\n", PROC_RANK, dtime()-t0);

		ttest2 = dtime();

		thrust::sort_by_key(d_keys[rankid].begin(), d_keys[rankid].begin() + pack->ntask, d_values[rankid].begin());
		// reorder <inode, jnode> pairs
		thrust::pair<thrust::host_vector<int, DevAllocator<int>>::iterator, thrust::host_vector<int, DevAllocator<int>>::iterator> new_end 
		= thrust::reduce_by_key(d_keys[rankid].begin(), d_keys[rankid].begin() + pack->ntask, d_one[rankid].begin(), d_inodes[rankid].begin(), d_nn[rankid].begin());
		thrust::exclusive_scan(d_nn[rankid].begin(), new_end.second, d_in[rankid].begin());

		int numInodes = new_end.second - (d_nn[rankid].begin());
		LOG(3, "[%d]:g%d numInodes: %d NLEAF: %d\n", PROC_RANK, device_ids[rankid], numInodes, nleaf);


		// p2p kernel
		int block = 24;
		int grid = numInodes; //80 * 4;
		// 3 gpuTs for pos (acc), (inode.pos, inode.acc, jnode.pos)
		int sharedMemSize = blockSize * sizeof(gpuT) * 3;

        uint64_t args[1];

        P2P_kernel_Params * params_ptr;
        checkDevErrors(devMalloc((void**)&params_ptr, sizeof(P2P_kernel_Params)));
        params_ptr->SoftenScale = SoftenScale;
        params_ptr->MASSPART = MASSPART;
        params_ptr->rrs = rrs;
        params_ptr->coeff = coeff;
        params_ptr->maxleaf = MAXLEAF;
        params_ptr->numInodes = numInodes;
        params_ptr->dev_inodes = thrust::raw_pointer_cast(d_inodes[rankid].data());
        params_ptr->dev_in = thrust::raw_pointer_cast(d_in[rankid].data());
        params_ptr->dev_nn = thrust::raw_pointer_cast(d_nn[rankid].data());
        params_ptr->dev_jnodes = thrust::raw_pointer_cast(d_values[rankid].data());
        params_ptr->ipart = pack->ipart;
        params_ptr->npart = pack->npart;
        params_ptr->dev_part_pos = d_part_pos[rankid];
        params_ptr->dev_part_acc = d_part_buf[rankid]; // acc
        params_ptr->dev_leaf_ip = d_leaf_ip[rankid];
        params_ptr->dev_leaf_np = d_leaf_np[rankid];
        params_ptr->dev_part_act = (uint_8b*)0;
        params_ptr->active_level = (uint_8b)0;

        args[0] = (uint64_t)params_ptr;

#ifdef JSPLIT
    int threads_group_id = hthread_group_create(current_cluster_id, 24, "P2P_kernel2",0,1,args);
#else
    int threads_group_id = hthread_group_create(current_cluster_id, 24, "P2P_kernel3",0,1,args);
#endif 

    checkDevErrors(hthread_group_wait(threads_group_id ));

	ttest3 = dtime();

	if (last_time) {
        	checkDevErrors(devEventRecord(stop_events[rankid], 0));
#ifdef TASKLOG
		task_record(0, tot_tasks);	
#endif
                LOG(2, " [%d] ntask = %d (call cnt: %d)\n", PROC_RANK, tot_tasks, call_cnt);
		tot_tasks = 0;
		call_cnt = 0;
	}
	ttest4 = dtime();
#ifdef MOREINFO
	printf("%lf %lf %lf %lf\n", ttest1-ttest0, ttest2-ttest1, ttest3-ttest2, ttest4-ttest3);
#endif
	return NULL;
}

extern "C" void *p2ppack_compute_gpu(P2PPack* packs, int npack) {
	int tot_tasks = 0;
	for (int i=0; i<npack; ++i)
		tot_tasks += packs[i].ntask;
#ifdef TASKLOG
	task_record(0, tot_tasks);	
#endif
        LOG(2, " [%d] ntask = %d (NPACK: %d)\n", PROC_RANK, tot_tasks, npack);
	if (npack == 0)
		return NULL;

	double t0, t1;

	double ttest0, ttest1, ttest2, ttest3, ttest4, ttest5;
	ttest0 = dtime();
	devBuffPrep(packs, npack);
	ttest1 = dtime();

	gpuT rrs = 1.0 / splitRadius;
	gpuT coeff = 2.0 * PIisq;

	int nleaf = last_leaf - first_leaf;
	hostmalloc_in(NPART, nleaf);
	t0 = dtime();
	for (int i=0; i<NPART; ++i) {
//		h_part_pos[i] = part[i].pos[0];
//		h_part_pos[i+NPART] = part[i].pos[1];
//		h_part_pos[i+2*NPART] = part[i].pos[2];
#ifndef INTXYZ
		h_part_pos[i] = part[i].pos[0];
		h_part_pos[i+NPART] = part[i].pos[1];
		h_part_pos[i+2*NPART] = part[i].pos[2];

#ifdef TEST_FLOAT
		float dispf = 2000000.0;
		
		h_part_pos[i] += dispf;
		h_part_pos[i+NPART] += dispf;
		h_part_pos[i+2*NPART] += dispf;
#endif

#else
	//	h_part_pos[i] = (float) (INT2POS * part[i].posi[0]);
	//	h_part_pos[i+NPART] = (float)(INT2POS * part[i].posi[1]);
	//	h_part_pos[i+2*NPART] = (float)(INT2POS *  part[i].posi[2]);

		h_part_pos[i] = (float) (INT2POS * (part[i].posi[0] - domain_center_i[0]) );
		h_part_pos[i+NPART] = (float)(INT2POS * (part[i].posi[1] - domain_center_i[1]) );
		h_part_pos[i+2*NPART] = (float)(INT2POS * (part[i].posi[2] - domain_center_i[2]) );

#endif




	}
	LOG(3, "[%d] part (NPART %d) reorder time: %lf s\n", PROC_RANK, NPART, dtime()-t0);
	t0 = dtime();
	for (int i=0; i<nleaf; ++i) {
		h_leaf_ip[i] = leaf[i+first_leaf].ipart;
		h_leaf_np[i] = leaf[i+first_leaf].npart;
	}
	LOG(3, "[%d] leaf reorder time: %lf s\n", PROC_RANK, dtime()-t0);

	int* numInodes_arr = (int*)malloc(sizeof(int) * npack);

#pragma omp parallel for num_threads(npack)
	for (int n=0; n<npack; ++n) {
		if (packs[n].ntask == 0)
			continue;
		int i = packs[n].rankid;
 
		checkDevErrors(devSetDevice(device_ids[i]));

		t1 = dtime();
		devMemcpyHD_task(packs, n);
		LOG(3, "[%d] copy task list %d time: %lf s\n", PROC_RANK, n, dtime()-t1);
	
		t1 = dtime();
		devMemcpyHD_part(packs, n);
		LOG(3, "[%d] copy input time: %lf s\n", PROC_RANK, dtime()-t1);

		thrust::sort_by_key(d_keys[i].begin() + off_ntask[n], d_keys[i].begin() + off_ntask[n] + packs[n].ntask, d_values[i].begin() + off_ntask[n]);
		// reorder <inode, jnode> pairs
		thrust::pair<thrust::host_vector<int, DevAllocator<int>>::iterator, thrust::host_vector<int, DevAllocator<int>>::iterator> new_end 
		= thrust::reduce_by_key(d_keys[i].begin() + off_ntask[n], d_keys[i].begin() + off_ntask[n] + packs[n].ntask, d_one[i].begin() + off_ntask[n], d_inodes[i].begin() + off_nleaf[n], d_nn[i].begin() + off_nleaf[n]);
		thrust::exclusive_scan(d_nn[i].begin() + off_nleaf[n], new_end.second, d_in[i].begin() + off_nleaf[n]);
		numInodes_arr[n] = new_end.second - (d_nn[i].begin() + off_nleaf[n]);
	}
	hostfree_in();
	ttest2 = dtime();

	for (int i=0; i<device_cnt; ++i) {
		checkDevErrors(devSetDevice(device_ids[i]));
        	checkDevErrors( devEventCreate(&start_events[i]) );
        	checkDevErrors( devEventCreate(&stop_events[i]) );
        	checkDevErrors( devEventRecord(start_events[i], 0));
        	checkDevErrors( devEventRecord(stop_events[i], 0));
	}

#pragma omp parallel for num_threads(npack)
	for (int n=0; n<npack; ++n) {
		if (packs[n].ntask == 0)
			continue;
		int i = packs[n].rankid;

		checkDevErrors(devSetDevice(device_ids[i]));

        checkDevErrors(hthread_dev_open(current_cluster_id));
        checkDevErrors(hthread_dat_load(current_cluster_id,"../obj/test.dat"));

		int numInodes = numInodes_arr[n]; 

		LOG(3, "[%d]:g%d numInodes: %d NLEAF: %d\n", PROC_RANK, device_ids[i], numInodes, packs[n].nleaf);
		// p2p kernel

		int blockSize = 256;

		int numBlocks = numInodes; //80 * 4;
		// 3 gpuTs for pos (acc), (inode.pos, inode.acc, jnode.pos)
		int sharedMemSize = blockSize * sizeof(gpuT) * 3;


        uint64_t args[1];

        P2P_kernel_Params * params_ptr;
        checkDevErrors(devMalloc((void**)&params_ptr, sizeof(P2P_kernel_Params)));
        params_ptr->SoftenScale = SoftenScale;
        params_ptr->MASSPART = MASSPART;
        params_ptr->rrs = rrs;
        params_ptr->coeff = coeff;
        params_ptr->maxleaf = MAXLEAF;
        params_ptr->numInodes = numInodes;
        params_ptr->dev_inodes = thrust::raw_pointer_cast(d_inodes[i].data() + off_nleaf[n]);
        params_ptr->dev_in = thrust::raw_pointer_cast(d_in[i].data() + off_nleaf[n]);
        params_ptr->dev_nn = thrust::raw_pointer_cast(d_nn[i].data() + off_nleaf[n]);
        params_ptr->dev_jnodes = thrust::raw_pointer_cast(d_values[i].data() + off_ntask[n]);
        params_ptr->ipart = packs[n].ipart;
        params_ptr->npart = packs[n].npart;
        params_ptr->dev_part_pos = d_part_pos[i] + off_npart[n]*3;
        params_ptr->dev_part_acc = d_part_buf[i] + off_npart[n]*3; // acc
        params_ptr->dev_leaf_ip = d_leaf_ip[i] + off_nleaf[n];
        params_ptr->dev_leaf_np = d_leaf_np[i] + off_nleaf[n];
        params_ptr->dev_part_act = (uint_8b*)0;
        params_ptr->active_level = (uint_8b)0;

        args[0] = (uint64_t)params_ptr;

    #ifdef JSPLIT
        int threads_group_id = hthread_group_create(current_cluster_id, 24, "P2P_kernel2",0,1,args);
    #else
        int threads_group_id = hthread_group_create(current_cluster_id, 24, "P2P_kernel3",0,1,args);
    #endif 

        checkDevErrors(hthread_group_wait(threads_group_id ));
	}
	ttest3 = dtime();

	for (int i=0; i<device_cnt; ++i) {
		checkDevErrors(devSetDevice(device_ids[i]));
        	checkDevErrors(devEventRecord(stop_events[i], 0));
	}

	free(numInodes_arr);
#ifdef MOREINFO
	printf("%lf %lf %lf\n", ttest1-ttest0, ttest2-ttest1, ttest3-ttest2);
#endif

	return NULL;
}

#include "gpu_active.impl.h"
