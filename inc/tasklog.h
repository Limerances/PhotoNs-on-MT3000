#ifndef TASKLOG_H
#define TASKLOG_H

#ifdef __cplusplus
extern "C" {
#endif
void task_record(int active_level, int ntasks);
void task_update_log(int loopid);
#ifdef __cplusplus
}
#endif

#endif
