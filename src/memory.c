#include "photoNs.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "memory.h"

#define FMM_MLOG_BASE 0
#define PM_MLOG_BASE 100
#define ANLYS_MLOG_BASE 200
#define UNDEF_MLOG_BASE 300

#define MAXN_MLOG 1000
#define MAXN_THREADS 128
#define MAXN_LOGPART 10

typedef struct mem_used_log {
	size_t cur_memused;
	size_t max_memused;
} MemLog;

size_t memused[MAXN_THREADS][MAXN_MLOG];
MemLog partial_mem[MAXN_THREADS][MAXN_LOGPART];
MemLog global_mem[MAXN_THREADS];

static FILE* flog[MAXN_THREADS];
int loop = -1;

void reset_memlog() {
	int i, j;
	int nth = omp_get_max_threads();
	for (i=0; i<nth; ++i) {
		global_mem[i].cur_memused = 0;
		global_mem[i].max_memused = 0;
		for (j=0; j<MAXN_LOGPART; ++j) {
			partial_mem[i][j].cur_memused = 0;
			partial_mem[i][j].max_memused = 0;
		}
	}
}

void initialize_memlogfile() {
	int i;
	int nth = omp_get_max_threads();
	int r_snap = 0;
	make_dir(PathSnapshot, "memlog", r_snap);
	char fname[256];
	for (i=0; i<nth; ++i) {
		sprintf(fname, "%s/memlog_%03d/%d_%d", PathSnapshot, r_snap, PROC_RANK, i);	
		flog[i] = fopen(fname, "a");
		if (0 == flog[i]) {
			printf("[%d] ith %d error memory log file %s!\n", PROC_RANK, i, fname);
			exit(0);
		}
		fprintf(flog[i], "\ntime loopid idx size logid partmem_cur partmem_max gmem_cur gmem_max\n");
	}

	reset_memlog();
}

void finalize_memlogfile() {
	int i;
	int nth = omp_get_max_threads();
	for (i=0; i<nth; ++i) {
		if (0 != flog[i]) {
			fclose(flog[i]);
		}
	}

}

void get_mem_info(int idx, char* mem_msg, int* logid) {
	if (idx < PM_MLOG_BASE) {
		sprintf(mem_msg, "fmm");
		*logid = 0;
	} else if (idx < ANLYS_MLOG_BASE) {
		sprintf(mem_msg, "pm");
		*logid = 1;
	} else if (idx < UNDEF_MLOG_BASE) {
		sprintf(mem_msg, "analysis");
		*logid = 2;
	} else {
		sprintf(mem_msg, "undefined");
		*logid = 3;
	}
}

void malloc_update_log(size_t size, int idx) {
	char mem_msg[64];
	int logid;
	get_mem_info(idx, mem_msg, &logid);
	int ith = omp_get_thread_num();
	
	if (memused[ith][idx] != 0) {
		printf("ith %d double malloc %d\n", ith, idx);
		exit(0);
	}
	memused[ith][idx] = size;	
	partial_mem[ith][logid].cur_memused += size;
	if (partial_mem[ith][logid].max_memused < partial_mem[ith][logid].cur_memused) {
		partial_mem[ith][logid].max_memused = partial_mem[ith][logid].cur_memused;
	}

	global_mem[ith].cur_memused += size;
	if (global_mem[ith].max_memused < global_mem[ith].cur_memused) {
		global_mem[ith].max_memused = global_mem[ith].cur_memused;
	}
	
	fprintf(flog[ith], "%lf %d %d %lld %d %lu %lu %lu %lu\n", dtime(), loop, idx, size, logid, partial_mem[ith][logid].cur_memused, partial_mem[ith][logid].max_memused, global_mem[ith].cur_memused, global_mem[ith].max_memused);
}

void free_update_log(int idx) {
	char mem_msg[64];
	int logid;
	get_mem_info(idx, mem_msg, &logid);

	int ith = omp_get_thread_num();

	size_t old_sz = memused[ith][idx];
	if (old_sz == 0) {
		printf("ith %d double free %d\n", ith, idx);
		exit(0);
	}
	memused[ith][idx] = 0;	
	partial_mem[ith][logid].cur_memused -= old_sz;
	global_mem[ith].cur_memused -= old_sz;

	fprintf(flog[ith], "%lf %d %d %lld %d %lu %lu %lu %lu\n", dtime(), loop, idx, -old_sz, logid, partial_mem[ith][logid].cur_memused, partial_mem[ith][logid].max_memused, global_mem[ith].cur_memused, global_mem[ith].max_memused);
}

void realloc_update_log(size_t size, int idx) {
	char mem_msg[64];
	int logid;
	get_mem_info(idx, mem_msg, &logid);
	int ith = omp_get_thread_num();
	
	size_t old_sz = memused[ith][idx];
	memused[ith][idx] = size;	
	partial_mem[ith][logid].cur_memused += size - old_sz;
	if (partial_mem[ith][logid].max_memused < partial_mem[ith][logid].cur_memused) {
		partial_mem[ith][logid].max_memused = partial_mem[ith][logid].cur_memused;
	}

	global_mem[ith].cur_memused += size - old_sz;
	if (global_mem[ith].max_memused < global_mem[ith].cur_memused) {
		global_mem[ith].max_memused = global_mem[ith].cur_memused;
	}

	fprintf(flog[ith], "%lf %d %d %lld %d %lu %lu %lu %lu\n", dtime(), loop, idx, size - old_sz, logid, partial_mem[ith][logid].cur_memused, partial_mem[ith][logid].max_memused, global_mem[ith].cur_memused, global_mem[ith].max_memused);
}

void* mymalloc(size_t size, int idx) {
	void* alloc = malloc(size);

	if (alloc == 0) {
		printf("[%d] malloc failed idx=%d size=%lu\n", PROC_RANK, idx, size);
		exit(0);
	}
	if (idx < 0 || idx >= MAXN_MLOG) {
		printf("wrong idx in mymalloc.\n");
		exit(0);
	}

	///////
	memset(alloc, 0, size);
	///////

#ifdef MEMLOG
	malloc_update_log(size, idx);
#endif
	return alloc;
}

void myfree(void* ptr, int idx) {
	if (idx >= MAXN_MLOG) {
		printf("wrong idx in myfree.\n");
		exit(0);
	}
	if (ptr)
		free(ptr);

#ifdef MEMLOG
	free_update_log(idx);
#endif
}

void* myrealloc(void* ptr, size_t size, int idx) {
	void* alloc = realloc(ptr, size);

	if (alloc == 0) {
		printf("[%d] realloc failed idx=%d size=%lu\n", PROC_RANK, idx, size);
		exit(0);
	}
	if (idx < 0 || idx >= MAXN_MLOG) {
		printf("wrong idx in mymalloc.\n");
		exit(0);
	}
#ifdef MEMLOG
	realloc_update_log(size, idx);
#endif
	return alloc;
}

void mock_mymalloc_(size_t* psize, int* idx) {
	malloc_update_log(*psize, *idx);
}

void mock_myfree_(void* ptr, int* idx) {
	free_update_log(*idx);
}

void mem_update_log(int loopid) {
	loop = loopid;
}
