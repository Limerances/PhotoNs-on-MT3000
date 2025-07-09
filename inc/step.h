#ifndef STEP_H
#define STEP_H

int active_label(cpuType loga_i, cpuType loga_f);
int active_label_ex(cpuType loga_i, cpuType loga_f);
void drift_step(cpuType dd) ;
void kick_half_step(cpuType dkh_pm) ;
void kick_half_active(cpuType dkh[], int level) ;

void kdk_level(int adaptive_level,int max, cpuType ddrift, cpuType dkickh);

#endif

