#ifndef MEMORY_H
#define MEMORY_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif
void* mymalloc(size_t size, int idx);
void myfree(void* ptr, int idx);
void* myrealloc(void* ptr, size_t size, int idx);
void malloc_update_log(size_t size, int idx);
void free_update_log(int idx);
void mock_mymalloc_(size_t* psize, int* idx);
void mock_myfree_(void* ptr, int* idx);
void mem_update_log(int loopid);
#ifdef __cplusplus
}
#endif

#endif
