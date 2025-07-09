#ifndef SCHEDULE_H
#define SCHEDULE_H

#include "typedef.h"

void create_output_list(int step_number, cpuType ai, cpuType af, cpuType aa_i, cpuType aa_f, int analysis_number);
void create_schedule(int step_number, cpuType ai, cpuType af, cpuType aa_i, cpuType aa_f, int analysis_number);

int analysis_output_index(int step);
int analysis_required(int step);

void free_schedule();


#endif
