#include "photoNs.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "tasklog.h"

static FILE* flog;

unsigned long ntasks_p2p[MAX_STEP_LEVEL+1];

void reset_record() {
	int i;
	for (i=0; i<=MAX_STEP_LEVEL; ++i)
		ntasks_p2p[i] = 0;
}

void initialize_tasklogfile() {
	int i;
	int r_snap = 0;
	make_dir(PathSnapshot, "tasklog", r_snap);
	char fname[256];
	sprintf(fname, "%s/tasklog_%03d/%d_%d", PathSnapshot, r_snap, PROC_RANK_SUDOM, PROC_RANK);	
	flog = fopen(fname, "a");
	if (0 == flog) {
		printf("[%d] error task log file %s!\n", PROC_RANK, fname);
		exit(0);
	}
	fprintf(flog, "\nloopid sudom ntasks timetasks ntasks_adptv timetasks_adptv level0 level1 level2 level3 level4 level5 ...\n");
	reset_record();
}

void finalize_tasklogfile() {
	if (0 != flog) {
		fclose(flog);
	}
}

void task_record(int active_level, int ntasks) {
	ntasks_p2p[active_level] += ntasks;	
}

void task_update_log(int loop) {
	int i;
	unsigned long ntasks = 0;
	unsigned long ntasks_adptv = 0;
	ntasks = ntasks_p2p[0];
	for (i=1; i<=MAX_STEP_LEVEL; ++i)
		ntasks_adptv +=ntasks_p2p[i];
	
	fprintf(flog, "%d %d", loop, PROC_RANK_SUDOM);
	fprintf(flog, " %lu %lf", ntasks, dtime_task);
	fprintf(flog, " %lu %lf", ntasks_adptv, dtime_adptv_task);
	for (i=0; i<=MAX_STEP_LEVEL; ++i)
		fprintf(flog, " %d", ntasks_p2p[i]);
	fprintf(flog, "\n");
	reset_record();
}
