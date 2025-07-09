#ifndef P2PTHRUST_H
#define P2PTHRUST_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int device_cnt;
int* device_ids;

// Copy data onto the GPU
// tasks
std::vector<thrust::device_vector<int> > d_keys;
std::vector<thrust::device_vector<int> > d_values;
std::vector<thrust::device_vector<int> > d_one;
// leaf
std::vector<thrust::device_vector<int> > d_inodes;
std::vector<thrust::device_vector<int> > d_nn;
std::vector<thrust::device_vector<int> > d_in;
std::vector<int*> d_leaf_ip;
std::vector<int*> d_leaf_np;
// particle
std::vector<gpuT*> d_part_pos;
std::vector<gpuT*> d_part_buf;
std::vector<uint_8b*> d_part_act;

#endif
