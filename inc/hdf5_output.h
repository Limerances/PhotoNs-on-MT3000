#ifndef HALO_OUTPUT_H
#define HALO_OUTPUT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>
#include <math.h>
#include "stdint.h"
#include <unistd.h>
#include <algorithm>

#include "hdf5_utils.h"
#include "utility.h"
#include "memory.h"

#ifdef HALO_FOF
#include "fof.h"
#endif
#ifdef SUBHALO
#include "subfind.h"
#endif
using namespace std;

typedef struct halo_header {
	unsigned long ngrps;
	unsigned long nsubs;
	unsigned long nsubs_total;
	unsigned long ngrps_total;
	unsigned long subs_offset;
	unsigned long grps_offset;
	unsigned long npart_total;
	double MASSPART;
	double redshift;
	unsigned int fof_group_min_len;
	double fof_linklength;
	double BoxSize;
} Hheader;

typedef struct hdf5_dataset_info {
	char dname[128];
	hsize_t dims;
	hsize_t start;
	hsize_t count;
	hsize_t nvars_per_block;
	hid_t file_datatype;
	hid_t mem_datatype;
	void* buffer;
} H5Dinfo;

typedef struct halo_distribution_info {
	unsigned long intra_offset_halo;
	unsigned long intra_N_halo;
	unsigned long inter_offset_halo;
	unsigned long inter_N_halo;
	unsigned long intra_offset_subhalo;
	unsigned long intra_N_subhalo;
	unsigned long inter_offset_subhalo;
	unsigned long inter_N_subhalo;
	unsigned long intra_offset_part;
	unsigned long intra_N_part;
} HaloDistInfo;

int compare_ul(const void* a, const void* b);
void write_header_fields(hid_t handle, Hheader header);
void get_offsets(GroupCatalog* grpcat, subGroupCatalog* subgrpcat, MPI_Comm intra_file_comm, MPI_Comm inter_file_comm, HaloDistInfo& halodinfo);

template<typename Body_T>
void get_halo_dataset_info(int id, H5Dinfo& dinfo, HaloDistInfo halodinfo, Body_T* prt, GroupCatalog* grpcat) {

	Group* grps = &grpcat->group[0];
	unsigned long num_grp = grpcat->group_number;

	int* ibuf;
	float* fbuf;
	unsigned long* lbuf;
	void** buf = &dinfo.buffer;
	
	dinfo.dims = halodinfo.intra_N_halo;
	dinfo.start = halodinfo.intra_offset_halo;
	dinfo.count = num_grp;

	switch (id) {
		case 0:
			sprintf(dinfo.dname, "GroupLen");
			dinfo.file_datatype = H5T_NATIVE_UINT; 
			dinfo.mem_datatype = H5T_NATIVE_UINT; 
			dinfo.nvars_per_block = 1;
			*buf = (int*)mymalloc(sizeof(int) * num_grp * dinfo.nvars_per_block, 258);
			ibuf = (int*)(*buf);
			for (int i=0; i<num_grp; ++i)
				ibuf[i] = grps[i].npart; 
			break;	
		case 1:
			sprintf(dinfo.dname, "GroupFirstSub");
			dinfo.file_datatype = H5T_NATIVE_UINT64; 
			dinfo.mem_datatype = H5T_NATIVE_UINT64; 
			dinfo.nvars_per_block = 1;
			*buf = (int*)mymalloc(sizeof(unsigned long) * num_grp * dinfo.nvars_per_block, 258);
			lbuf = (unsigned long*)(*buf);
			lbuf[0] = halodinfo.inter_offset_subhalo + halodinfo.intra_offset_subhalo;
			for (int i=1; i<num_grp; ++i)
				lbuf[i] = lbuf[i-1] + grps[i-1].nsubhalo; 
			break;	
		case 2:
			sprintf(dinfo.dname, "GroupNsubs");
			dinfo.file_datatype = H5T_NATIVE_UINT; 
			dinfo.mem_datatype = H5T_NATIVE_UINT; 
			dinfo.nvars_per_block = 1;
			*buf = (int*)mymalloc(sizeof(int) * num_grp * dinfo.nvars_per_block, 258);
			ibuf = (int*)(*buf);
			for (int i=0; i<num_grp; ++i)
				ibuf[i] = grps[i].nsubhalo; 
			break;	
		case 3:
			sprintf(dinfo.dname, "GroupPos");
			dinfo.file_datatype = H5T_NATIVE_FLOAT; 
			dinfo.mem_datatype = H5T_NATIVE_FLOAT; 
			dinfo.nvars_per_block = 3;
			*buf = (float*)mymalloc(sizeof(float) * num_grp * dinfo.nvars_per_block, 258);
			fbuf = (float*)(*buf);
			for (int i=0; i<num_grp; ++i) {
				fbuf[3 * i] = grps[i].center[0]; 
				fbuf[3 * i + 1] = grps[i].center[1]; 
				fbuf[3 * i + 2] = grps[i].center[2]; 
			}
			break;
#ifdef HALOPROPERTY
		case 4:
			sprintf(dinfo.dname, "GroupMean200");
			dinfo.file_datatype = H5T_NATIVE_FLOAT; 
			dinfo.mem_datatype = H5T_NATIVE_FLOAT; 
			dinfo.nvars_per_block = 1;
			*buf = (float*)mymalloc(sizeof(float) * num_grp * dinfo.nvars_per_block, 258);
			fbuf = (float*)(*buf);
			for (int i=0; i<num_grp; ++i) {
				fbuf[i] = grps[i].Mmean200; 
			}		
			break;
		case 5:
			sprintf(dinfo.dname, "GroupCrit200");
			dinfo.file_datatype = H5T_NATIVE_FLOAT; 
			dinfo.mem_datatype = H5T_NATIVE_FLOAT; 
			dinfo.nvars_per_block = 1;
			*buf = (float*)mymalloc(sizeof(float) * num_grp * dinfo.nvars_per_block, 258);
			fbuf = (float*)(*buf);
			for (int i=0; i<num_grp; ++i) {
				fbuf[i] = grps[i].Mcrit200; 
			}		
			break;


#endif
	}
	return; 
}


template<typename Body_T>
void get_subhalo_dataset_info(int id, H5Dinfo& dinfo, HaloDistInfo halodinfo, int min_mostbound_npart, Body_T* prt, subGroupCatalog* subgrpcat) {

	subGroup* subgrp = &subgrpcat->subgrps[0];
	unsigned long num_subgrp = subgrpcat->subgrp_number;
	unsigned long tempN;

	int* ibuf;
	float* fbuf;
	unsigned long* lbuf;
	void** buf = &dinfo.buffer;
	
	dinfo.dims = halodinfo.intra_N_subhalo;
	dinfo.start = halodinfo.intra_offset_subhalo;
	dinfo.count = num_subgrp;

	switch (id) {
		case 0:
			sprintf(dinfo.dname, "SubhaloLen");
			dinfo.file_datatype = H5T_NATIVE_UINT; 
			dinfo.mem_datatype = H5T_NATIVE_UINT; 
			dinfo.nvars_per_block = 1;
			*buf = (int*)mymalloc(sizeof(int) * num_subgrp * dinfo.nvars_per_block, 258);
			ibuf = (int*)(*buf);
			for (int i=0; i<num_subgrp; ++i)
				ibuf[i] = subgrp[i].Len; 
			break;	
		case 1:
			sprintf(dinfo.dname, "SubhaloPos");
			dinfo.file_datatype = H5T_NATIVE_FLOAT; 
			dinfo.mem_datatype = H5T_NATIVE_FLOAT; 
			dinfo.nvars_per_block = 3;
			*buf = (float*)mymalloc(sizeof(float) * num_subgrp * dinfo.nvars_per_block, 258);
			fbuf = (float*)(*buf);
			for (int i=0; i<num_subgrp; ++i) {
				fbuf[3 * i] = subgrp[i].pos[0]; 
				fbuf[3 * i + 1] = subgrp[i].pos[1]; 
				fbuf[3 * i + 2] = subgrp[i].pos[2]; 
			}
			break;	
		case 2:
			sprintf(dinfo.dname, "SubhaloSpin");
			dinfo.file_datatype = H5T_NATIVE_FLOAT; 
			dinfo.mem_datatype = H5T_NATIVE_FLOAT; 
			dinfo.nvars_per_block = 3;
			*buf = (float*)mymalloc(sizeof(float) * num_subgrp * dinfo.nvars_per_block, 258);
			fbuf = (float*)(*buf);
			for (int i=0; i<num_subgrp; ++i) {
				fbuf[3 * i] = subgrp[i].spin[0]; 
				fbuf[3 * i + 1] = subgrp[i].spin[1]; 
				fbuf[3 * i + 2] = subgrp[i].spin[2]; 
			}
			break;	
		case 3:
			sprintf(dinfo.dname, "SubhaloVel");
			dinfo.file_datatype = H5T_NATIVE_FLOAT; 
			dinfo.mem_datatype = H5T_NATIVE_FLOAT; 
			dinfo.nvars_per_block = 3;
			*buf = (float*)mymalloc(sizeof(float) * num_subgrp * dinfo.nvars_per_block, 258);
			fbuf = (float*)(*buf);
			for (int i=0; i<num_subgrp; ++i) {
				fbuf[3 * i] = subgrp[i].vel[0]; 
				fbuf[3 * i + 1] = subgrp[i].vel[1]; 
				fbuf[3 * i + 2] = subgrp[i].vel[2]; 
			}
			break;	
		case 4:
			sprintf(dinfo.dname, "SubhaloVmax");
			dinfo.file_datatype = H5T_NATIVE_FLOAT; 
			dinfo.mem_datatype = H5T_NATIVE_FLOAT; 
			dinfo.nvars_per_block = 1;
			*buf = (float*)mymalloc(sizeof(float) * num_subgrp * dinfo.nvars_per_block, 258);
			fbuf = (float*)(*buf);
			for (int i=0; i<num_subgrp; ++i) {
				fbuf[i] = subgrp[i].SubVmax;
			}
			break;	
		case 5:
			sprintf(dinfo.dname, "SubhaloVmaxRad");
			dinfo.file_datatype = H5T_NATIVE_FLOAT; 
			dinfo.mem_datatype = H5T_NATIVE_FLOAT; 
			dinfo.nvars_per_block = 1;
			*buf = (float*)mymalloc(sizeof(float) * num_subgrp * dinfo.nvars_per_block, 258);
			fbuf = (float*)(*buf);
			for (int i=0; i<num_subgrp; ++i) {
				fbuf[i] = subgrp[i].SubVmaxRad; 
			}
			break;	
		case 6:
			sprintf(dinfo.dname, "SubhaloGroupNr");
			dinfo.file_datatype = H5T_NATIVE_UINT64;
			dinfo.mem_datatype = H5T_NATIVE_UINT64;
			dinfo.nvars_per_block = 1;
			*buf = (unsigned long*)mymalloc(sizeof(unsigned long) * num_subgrp * dinfo.nvars_per_block, 258);
			lbuf = (unsigned long*)(*buf);
			tempN = halodinfo.inter_offset_halo + halodinfo.intra_offset_halo;
			for (int i=0; i<num_subgrp; ++i) {
				lbuf[i] = tempN + (unsigned long)subgrp[i].in_group; 
			}
			break;	
		case 7:
			sprintf(dinfo.dname, "ParticleIDsMostbound");
			dinfo.file_datatype = H5T_NATIVE_UINT64;
			dinfo.mem_datatype = H5T_NATIVE_UINT64;
			dinfo.nvars_per_block = min_mostbound_npart;
			*buf = (unsigned long*)mymalloc(sizeof(unsigned long) * num_subgrp * dinfo.nvars_per_block, 258);
			lbuf = (unsigned long*)(*buf);
			for (int i=0; i<num_subgrp; ++i) {
				for (int j=0; j<dinfo.nvars_per_block;++j)
					lbuf[dinfo.nvars_per_block * i + j] = prt[subgrp[i].pidx[j]].id; 
			}
			break;
		case 8:
			sprintf(dinfo.dname, "SubhaloParticleIDs");
			dinfo.file_datatype = H5T_NATIVE_UINT64;
			dinfo.mem_datatype = H5T_NATIVE_UINT64;
			dinfo.nvars_per_block = 1;
			dinfo.dims = halodinfo.intra_N_part;
			dinfo.start = halodinfo.intra_offset_part;
			tempN = 0;
			for (int i=0; i<num_subgrp; ++i) {
				tempN += subgrp[i].Len; 
			}
			dinfo.count = tempN;
			
			*buf = (unsigned long*)mymalloc(sizeof(unsigned long) * tempN, 258);
			lbuf = (unsigned long*)(*buf);
			tempN = 0;
			for (int i=0; i<num_subgrp; ++i) {
				for (int j=0; j<subgrp[i].Len;++j)
					lbuf[tempN++] = prt[subgrp[i].pidx[j]].id; 
			}
			break;

	}
	return; 
}


template<typename Body_T>
void hdf5_output(SIMInfo<Body_T> um, GroupCatalog* grpcat, subGroupCatalog* subgrpcat, char* fn, int file_id, int file_rank, int min_mostbound_npart) {
	unsigned long n, m;

	Body_T* prt = um.part;
	int PROC_RANK = um.PROC_RANK;
	int PROC_SIZE = um.PROC_SIZE;
	float MASSPART = um.MASSPART;
	float BOXSIZE = um.BOXSIZE;
	float redshift = um.redshift;
	float npart_tot = um.NPART_TOTAL;

	MPI_Comm intra_file_comm, inter_file_comm;
	MPI_Comm_split(MPI_COMM_WORLD, file_id, PROC_RANK, &intra_file_comm);
	MPI_Comm_split(MPI_COMM_WORLD, file_rank, PROC_RANK, &inter_file_comm);
	int intra_size, intra_rank;
	MPI_Comm_size(intra_file_comm, &intra_size); 
	MPI_Comm_rank(intra_file_comm, &intra_rank);

	// get the info of halos and subhalos on each process
	HaloDistInfo halodinfo;
	get_offsets(grpcat, subgrpcat, intra_file_comm, inter_file_comm, halodinfo);

	char fsubhalo[128];
	sprintf(fsubhalo, "%s.%d.hdf5", fn, file_id);

	hid_t hdf5_dataset, hdf5_dataspace_in_file, hdf5_dataspace_memory;
	hsize_t dims[2], count[2], start[2];
	int rank;

	hid_t plist_id = my_H5Pcreate(H5P_FILE_ACCESS);
	my_H5Pset_fapl_mpio(plist_id, intra_file_comm, MPI_INFO_NULL);

	hid_t hdf5_file = my_H5Fcreate(fsubhalo, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

	hid_t plist_id_xfer = my_H5Pcreate(H5P_DATASET_XFER);
	my_H5Pset_dxpl_mpio(plist_id_xfer, H5FD_MPIO_COLLECTIVE);

	hid_t hdf5_header = my_H5Gcreate(hdf5_file, "/Header", 0);
	Hheader header;
	header.ngrps = halodinfo.intra_N_halo; 
	header.nsubs = halodinfo.intra_N_subhalo;
	header.ngrps_total = halodinfo.inter_N_halo;
	header.nsubs_total = halodinfo.inter_N_subhalo;
	header.grps_offset = halodinfo.inter_offset_halo;
	header.subs_offset = halodinfo.inter_offset_subhalo;
	header.npart_total = npart_tot;
	header.MASSPART = MASSPART;
	header.redshift = redshift;
	header.fof_group_min_len = um.fof_group_min_len;
	header.fof_linklength = um.fof_linklength;
	header.BoxSize = BOXSIZE;
	write_header_fields(hdf5_header, header);

        hid_t hdf5_grp = my_H5Gcreate(hdf5_file, "/FOFHalo", 0);

	H5Dinfo dinfo;
#ifdef HALOPROPERTY
	int N_Blocks_H = 6;
#else
	int N_Blocks_H = 4;
#endif
	for (int bid=0; bid < N_Blocks_H; ++bid) {
		get_halo_dataset_info(bid, dinfo, halodinfo, prt, grpcat);
	      	if (PROC_RANK == PROC_SIZE - 1) printf("halo block %d (%s)...\n", bid, dinfo.dname);
		dims[0] = dinfo.dims;
		dims[1] = dinfo.nvars_per_block;
		if(dims[1] == 1)
		  rank = 1;
		else
		  rank = 2;
		hdf5_dataspace_in_file = my_H5Screate_simple(rank, dims, NULL);
		hdf5_dataset = my_H5Dcreate(hdf5_grp, dinfo.dname, dinfo.file_datatype, hdf5_dataspace_in_file, H5P_DEFAULT);
	
	        start[0] = dinfo.start;
	        start[1] = 0;
	
	        count[0] = dinfo.count;
	        count[1] = dinfo.nvars_per_block;
	
	        my_H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET, start, NULL, count, NULL);
	
	        hdf5_dataspace_memory = my_H5Screate_simple(rank, count, NULL);

	        my_H5Dwrite(hdf5_dataset, dinfo.mem_datatype, hdf5_dataspace_memory, hdf5_dataspace_in_file,
	                    plist_id_xfer, dinfo.buffer, dinfo.dname);
		myfree(dinfo.buffer, 258);
	
		my_H5Sclose(hdf5_dataspace_memory, H5S_SIMPLE);
		my_H5Dclose(hdf5_dataset, dinfo.dname);
		my_H5Sclose(hdf5_dataspace_in_file, H5S_SIMPLE);
	
	}

        hid_t hdf5_subgrp = my_H5Gcreate(hdf5_file, "/Subhalo", 0);

	int N_Blocks_SH = 9;
	for (int bid=0; bid < N_Blocks_SH; ++bid) {
		get_subhalo_dataset_info(bid, dinfo, halodinfo, min_mostbound_npart, prt, subgrpcat);
              	if (PROC_RANK == PROC_SIZE - 1) printf("Subhalo block %d (%s)...\n", bid, dinfo.dname);
		dims[0] = dinfo.dims;
		dims[1] = dinfo.nvars_per_block;
		if(dims[1] == 1)
		  rank = 1;
		else
		  rank = 2;
		hdf5_dataspace_in_file = my_H5Screate_simple(rank, dims, NULL);
		hdf5_dataset = my_H5Dcreate(hdf5_subgrp, dinfo.dname, dinfo.file_datatype, hdf5_dataspace_in_file, H5P_DEFAULT);

                start[0] = dinfo.start;
                start[1] = 0;

                count[0] = dinfo.count;
                count[1] = dinfo.nvars_per_block;

                my_H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET, start, NULL, count, NULL);

                hdf5_dataspace_memory = my_H5Screate_simple(rank, count, NULL);

                my_H5Dwrite(hdf5_dataset, dinfo.mem_datatype, hdf5_dataspace_memory, hdf5_dataspace_in_file,
                            plist_id_xfer, dinfo.buffer, dinfo.dname);
		myfree(dinfo.buffer, 258);

		my_H5Sclose(hdf5_dataspace_memory, H5S_SIMPLE);
		my_H5Dclose(hdf5_dataset, dinfo.dname);
		my_H5Sclose(hdf5_dataspace_in_file, H5S_SIMPLE);

	}

        my_H5Gclose(hdf5_grp, "/FOFHalo");
        my_H5Gclose(hdf5_subgrp, "/Subhalo");
        my_H5Gclose(hdf5_header, "/Header");
        my_H5Pclose(plist_id_xfer);
        my_H5Pclose(plist_id);
        my_H5Fclose(hdf5_file, fsubhalo);

	MPI_Comm_free(&intra_file_comm);
	MPI_Comm_free(&inter_file_comm);
}

#endif
