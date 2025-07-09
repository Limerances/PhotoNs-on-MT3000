extern "C" void *task_compute_p2p_act_gpu(int ntask, int* ltask_t, int* ltask_s, int last_time, int active_level) 
{
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

	if (first_time) NPACK++;
	devBuffPrep(pack, NPACK);

	int nleaf = last_leaf - first_leaf;
	if (pack->nleaf != nleaf || pack->npart != NPART) {
		printf("Pack info error!\n");
		exit(0);
	}

	int rankid = pack->rankid;
	checkDevErrors(devSetDevice(device_ids[rankid]));

	gpuT rrs = 1.0 / splitRadius;
	gpuT coeff = 2.0 * PIisq;

	if (first_time) {
		t0 = dtime();
		hostmalloc_in(NPART, nleaf);
		for (int i=0; i<NPART; ++i) {
#ifndef INTXYZ
		h_part_pos[i] = part[i].pos[0];
		h_part_pos[i+NPART] = part[i].pos[1];
		h_part_pos[i+2*NPART] = part[i].pos[2];
#else
		h_part_pos[i] = (float) (INT2POS * part[i].posi[0]);
		h_part_pos[i+NPART] = (float)(INT2POS * part[i].posi[1]);
		h_part_pos[i+2*NPART] = (float)(INT2POS *  part[i].posi[2]);
#endif





//
//			h_part_pos[i] = part[i].pos[0];
//			h_part_pos[i+NPART] = part[i].pos[1];
//			h_part_pos[i+2*NPART] = part[i].pos[2];
#ifdef ACT_OPT
			h_part_act[i] = (uint_8b)part[i].act;
#endif
		}
		LOG(3, "active [%d] part (NPART %d) reorder time: %lf s\n", PROC_RANK, NPART, dtime()-t0);
		t0 = dtime();
		for (int i=0; i<nleaf; ++i) {
			h_leaf_ip[i] = leaf[i+NPART].ipart;
			h_leaf_np[i] = leaf[i+NPART].npart;
		}
		LOG(3, "active [%d] leaf reorder time: %lf s\n", PROC_RANK, dtime()-t0);
		t0 = dtime();
		devMemcpyHD_part(pack);
		LOG(3, "active [%d] copy input time: %lf s\n", PROC_RANK, dtime()-t0);
		hostfree_in();

        	checkDevErrors( devEventCreate(&start_events[rankid]) );
        	checkDevErrors( devEventCreate(&stop_events[rankid]) );
        	checkDevErrors( devEventRecord(start_events[rankid], 0));
        	checkDevErrors( devEventRecord(stop_events[rankid], 0));
	}

	if (pack->ntask > 0) {
		t0 = dtime();

		devMemcpyHD_task(pack);
		LOG(3, "active [%d] copy task list time: %lf s\n", PROC_RANK, dtime()-t0);

		thrust::sort_by_key(d_keys[rankid].begin(), d_keys[rankid].begin() + pack->ntask, d_values[rankid].begin());
	// reorder <inode, jnode> pairs
		thrust::pair<thrust::device_vector<int>::iterator, thrust::device_vector<int>::iterator> new_end 
		= thrust::reduce_by_key(d_keys[rankid].begin(), d_keys[rankid].begin() + pack->ntask, d_one[rankid].begin(), d_inodes[rankid].begin(), d_nn[rankid].begin());
		thrust::exclusive_scan(d_nn[rankid].begin(), new_end.second, d_in[rankid].begin());

		int numInodes = new_end.second - d_nn[rankid].begin();
		LOG(3, "active [%d]:g%d numInodes: %d NLEAF: %d\n", PROC_RANK, device_ids[rankid], numInodes, nleaf);

		// p2p kernel
		int blockSize = MAXLEAF;
		if (blockSize % WARP != 0) {
			LOG(0, "Error! blockSize for GPU kernel should be mutiple of WARP.\n");
			exit(0);
		}
		int numBlocks = numInodes;
		// 3 gpuTs for pos (acc), (inode.pos, inode.acc, jnode.pos)
		int sharedMemSize = blockSize * 3 * sizeof(gpuT);
#ifdef __CUDACC__
		checkKernelErrors((P2P_kernel<<< numBlocks, blockSize, sharedMemSize>>>(SoftenScale, MASSPART, rrs, coeff,
		numInodes, 
		thrust::raw_pointer_cast(d_inodes[rankid].data()),
		thrust::raw_pointer_cast(d_in[rankid].data()),
		thrust::raw_pointer_cast(d_nn[rankid].data()),
		thrust::raw_pointer_cast(d_values[rankid].data()),
		pack->ipart,
		pack->npart,
		d_part_pos[rankid],
		d_part_buf[rankid], // acc
		d_leaf_ip[rankid],
		d_leaf_np[rankid],
		d_part_act[rankid], // active
		(uint_8b)active_level)));
#else
		hipLaunchKernelGGL(P2P_kernel, numBlocks, blockSize, sharedMemSize, 0, SoftenScale, MASSPART, rrs, coeff,
		numInodes, 
		thrust::raw_pointer_cast(d_inodes[rankid].data()),
		thrust::raw_pointer_cast(d_in[rankid].data()),
		thrust::raw_pointer_cast(d_nn[rankid].data()),
		thrust::raw_pointer_cast(d_values[rankid].data()),
		pack->ipart,
		pack->npart,
		d_part_pos[rankid],
		d_part_buf[rankid], // acc
		d_leaf_ip[rankid],
		d_leaf_np[rankid],
		d_part_act[rankid], // active
		(uint_8b)active_level);
#endif
	}

	if (last_time) {
		checkDevErrors(devEventRecord(stop_events[rankid], 0));
		get_last_kntime(PROC_RANK);
#ifdef TASKLOG
		task_record(active_level, tot_tasks);	
#endif
		LOG(2, " active (%d) [%d] ntask = %d (call cnt: %d)\n", active_level, PROC_RANK, tot_tasks, call_cnt);
		tot_tasks = 0;
		call_cnt = 0;

		t0 = dtime();
		hostmalloc_out(NPART);
		devMemcpyDH_updateAcc(pack, NPACK, active_level);
		hostfree_out();
		LOG(3, "active [%d] output (NPART %d) copy and reorder time: %lf s\n", PROC_RANK, NPART, dtime()-t0);
	}

	return NULL;
}

extern "C" void *p2ppack_compute_act_gpu(P2PPack* packs, int npack, int active_level) {
	int tot_tasks = 0;
	for (int i=0; i<npack; ++i)
		tot_tasks += packs[i].ntask;
#ifdef TASKLOG
	task_record(active_level, tot_tasks);	
#endif
	LOG(2, " active (%d) [%d] ntask = %d (NPACK: %d)\n", active_level, PROC_RANK, tot_tasks, npack);
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
///		h_part_pos[i+NPART] = part[i].pos[1];
//		h_part_pos[i+2*NPART] = part[i].pos[2];
#ifndef INTXYZ
		h_part_pos[i] = part[i].pos[0];
		h_part_pos[i+NPART] = part[i].pos[1];
		h_part_pos[i+2*NPART] = part[i].pos[2];
#else
		h_part_pos[i] = (float) (INT2POS * part[i].posi[0]);
		h_part_pos[i+NPART] = (float)(INT2POS * part[i].posi[1]);
		h_part_pos[i+2*NPART] = (float)(INT2POS *  part[i].posi[2]);
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
		thrust::pair<thrust::device_vector<int>::iterator, thrust::device_vector<int>::iterator> new_end 
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

		int numInodes = numInodes_arr[n]; 

		LOG(3, "[%d]:g%d numInodes: %d NLEAF: %d\n", PROC_RANK, device_ids[i], numInodes, packs[n].nleaf);
		// p2p kernel
#ifdef JSPLIT
		int blockSize = 256;
#else
		int blockSize = MAXLEAF;
#endif
		if (blockSize % WARP != 0) {
			LOG(0, "Error! blockSize for GPU kernel should be mutiple of WARP.\n");
			exit(0);
		}
		int numBlocks = numInodes; //80 * 4;
		// 3 gpuTs for pos (acc), (inode.pos, inode.acc, jnode.pos)
		int sharedMemSize = blockSize * sizeof(gpuT) * 3;
#ifdef __CUDACC__
#ifdef JSPLIT
		if (blockSize % MAXLEAF != 0 || blockSize > 1024) {
			LOG(0, "Error! blockSize (<=1024) for GPU p2p_kernel2 should be mutiple of MAXLEAF.\n");
			exit(0);
		}
		checkKernelErrors((P2P_kernel2<<< numBlocks, blockSize, sharedMemSize>>>(SoftenScale, MASSPART, rrs, coeff, MAXLEAF,
#else
		checkKernelErrors((P2P_kernel3<<< numBlocks, blockSize, sharedMemSize>>>(SoftenScale, MASSPART, rrs, coeff,
#endif 
			numInodes, 
			thrust::raw_pointer_cast(d_inodes[i].data() + off_nleaf[n]),
			thrust::raw_pointer_cast(d_in[i].data() + off_nleaf[n]),
			thrust::raw_pointer_cast(d_nn[i].data() + off_nleaf[n]),
			thrust::raw_pointer_cast(d_values[i].data() + off_ntask[n]),
			packs[n].ipart,
			packs[n].npart,
			d_part_pos[i] + off_npart[n]*3,
			d_part_buf[i] + off_npart[n]*3, // acc
			d_leaf_ip[i] + off_nleaf[n],
			d_leaf_np[i] + off_nleaf[n])));
#else
		hipLaunchKernelGGL(P2P_kernel, numBlocks, blockSize, sharedMemSize, 0, SoftenScale, MASSPART, rrs, coeff,
			numInodes, 
			thrust::raw_pointer_cast(d_inodes[i].data() + off_nleaf[n]),
			thrust::raw_pointer_cast(d_in[i].data() + off_nleaf[n]),
			thrust::raw_pointer_cast(d_nn[i].data() + off_nleaf[n]),
			thrust::raw_pointer_cast(d_values[i].data() + off_ntask[n]),
			packs[n].ipart,
			packs[n].npart,
			d_part_pos[i] + off_npart[n]*3,
			d_part_buf[i] + off_npart[n]*3, // acc
			d_leaf_ip[i] + off_nleaf[n],
			d_leaf_np[i] + off_nleaf[n],
			(uint_8b*)0,
			(uint_8b)0);
#endif
	}
	ttest3 = dtime();

	for (int i=0; i<device_cnt; ++i) {
		checkDevErrors(devSetDevice(device_ids[i]));
        	checkDevErrors(devEventRecord(stop_events[i], 0));
	}
	get_last_kntime(PROC_RANK);
	hostmalloc_out(NPART);
	devMemcpyDH_updateAcc(packs, NPACK, active_level);
	hostfree_out();

	free(numInodes_arr);
#ifdef MOREINFO
	printf("%lf %lf %lf\n", ttest1-ttest0, ttest2-ttest1, ttest3-ttest2);
#endif

	return NULL;
}
