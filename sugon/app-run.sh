#!/bin/bash

export LD_LIBRARY_PATH=../lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/public/home/wangqiao/local/fftw-3.3.8/lib:/public/home/wangqiao/local/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=8

#ulimit -s
#export UCX_TLS=rc_x,self,shm

APP="../run/PhotoNs-GPU"

local_rank=$(( $OMPI_COMM_WORLD_LOCAL_RANK % 4 ))

case ${local_rank} in
[0])
  export UCX_NET_DEVICES="mlx5_0:1"
  export UCX_IB_PCI_BW="mlx5_0:50Gbs"
  export HIP_VISIBLE_DEVICES=0
  numactl --cpunodebind=0 --membind=0 \
  $APP
  ;;
[1])
  export UCX_NET_DEVICES="mlx5_1:1"
  export UCX_IB_PCI_BW="mlx5_1:50Gbs"
  export HIP_VISIBLE_DEVICES=1
  numactl --cpunodebind=1 --membind=1 \
  $APP
  ;;
[2])
  export UCX_NET_DEVICES="mlx5_2:1"
  export UCX_IB_PCI_BW="mlx5_2:50Gbs"
  export HIP_VISIBLE_DEVICES=2
  numactl --cpunodebind=2 --membind=2 \
  $APP
  ;;
[3])
  export UCX_NET_DEVICES="mlx5_3:1"
  export UCX_IB_PCI_BW="mlx5_3:50Gbs"
  export HIP_VISIBLE_DEVICES=3
  numactl --cpunodebind=3 --membind=3 \
  $APP
  ;;
esac
