#ifndef GRAVITY_H
#define GRAVITY_H

void partmesh_kick(cpuType dkick , int mode_alloc );
void gravity_pm_fmm() ;
void gravity_fmm_act(int level, int max_level) ;

void gravity_fmm_pm_kick(cpuType dkick) ;

void acc_cln();

void gravity_cln();

#endif
