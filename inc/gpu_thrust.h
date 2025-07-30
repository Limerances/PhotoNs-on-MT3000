#ifndef P2PTHRUST_H
#define P2PTHRUST_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int device_cnt;
int* device_ids;

template <typename T>
struct DevAllocator {
    using value_type = T;

    DevAllocator() {}

    T* allocate(std::size_t n) {
        T* ptr = nullptr;
        checkDevErrors(devMalloc((void**)&ptr, sizeof(T) * n));
        return ptr;
    }

    void deallocate(T* p, std::size_t) {
        devFree(p);
    }

    template<typename U>
    struct rebind {
        typedef DevAllocator<U> other;
    };
};

// Copy data onto the GPU
// tasks
std::vector<thrust::host_vector<int, DevAllocator<int>> > d_keys;
std::vector<thrust::host_vector<int, DevAllocator<int>> > d_values;
std::vector<thrust::host_vector<int, DevAllocator<int>> > d_one;
// leaf
std::vector<thrust::host_vector<int, DevAllocator<int>> > d_inodes;
std::vector<thrust::host_vector<int, DevAllocator<int>> > d_nn;
std::vector<thrust::host_vector<int, DevAllocator<int>> > d_in;
std::vector<int*> d_leaf_ip;
std::vector<int*> d_leaf_np;
// particle
std::vector<gpuT*> d_part_pos;
std::vector<gpuT*> d_part_buf;
std::vector<uint_8b*> d_part_act;

#endif
