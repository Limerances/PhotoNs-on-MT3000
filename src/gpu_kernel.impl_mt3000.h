#define NTABINT 512
#ifdef GPU_DP
typedef double3 gpuT3;
#define HALF 0.5
#define ONE_3 (1.0/3)
#define FIVE_3 (5.0/3)
#define TWO_3 (2.0/3)
#define ONE_1 (1.0)
#define THREE_2 (3.0/2)
#define WIDTH (3.0/NTABINT)
#define WIDTH_INV (NTABINT/3.0)
#else
typedef float3 gpuT3;
#define HALF 0.5f
#define ONE_3 (1.0f/3)
#define FIVE_3 (5.0f/3)
#define TWO_3 (2.0f/3)
#define ONE_1 (1.0f)
#define THREE_2 (3.0f/2)
#define WIDTH (3.0f/NTABINT)
#define WIDTH_INV (NTABINT/3.0f)
#endif

//3 mul , 1 fma, 1 sub, 1 shr, 1 add
//mov.u64         %SPL, __local_depot7;
//ld.param.f32    %f1, [_Z7Q_rsqrtf_param_0];
//add.u64         %rd2, %SPL, 0;
//st.local.f32    [%rd2], %f1;
//ld.local.u64    %rd3, [%rd2];
//shr.s64         %rd4, %rd3, 1;
//mov.u64         %rd5, 1597463007;
//sub.s64         %rd6, %rd5, %rd4;
//mov.b64         {%r1, %r2}, %rd6;
//mov.b32         %f2, %r1;
//cvt.u32.u64     %r3, %rd3;
//mov.b32         %f3, %r3;
//mul.f32         %f4, %f3, 0fBF000000;
//mul.f32         %f5, %f4, %f2;
//fma.rn.f32      %f6, %f2, %f5, 0f3FC00000;
//mul.f32         %f7, %f2, %f6;
//st.param.f32    [func_retval0+0], %f7;
__device__ gpuT Q_rsqrt(gpuT number) {
	long i = *(long*)&number;
	i = 0x5F3759DF - (i >> 1);
	gpuT y = *(gpuT*)&i;
	y *= (THREE_2 - (number * HALF * y * y));
	return y;
}

#ifdef SPI
__device__ gpuT3 P2P_core(gpuT SoftenScale, gpuT MASSPART, gpuT rrs, gpuT coeff,
		gpuT3 ii, gpuT3 jj, gpuT3 acc, gpuT *tab0, gpuT *tab1)
#else
__device__ gpuT3 P2P_core(gpuT SoftenScale, gpuT MASSPART, gpuT rrs, gpuT coeff,
		gpuT3 ii, gpuT3 jj, gpuT3 acc)
#endif
{
	gpuT dr, ir3, half = 0.5;
	gpuT3 dx;
	dx.x = jj.x - ii.x;
	dx.y = jj.y - ii.y;
	dx.z = jj.z - ii.z;
	gpuT rdr, dr2;
	dr2 = dx.x * dx.x + dx.y * dx.y + dx.z * dx.z;

	if (dr2 > SoftenScale * SoftenScale) {
#ifdef HIP_RTM
#ifdef GPU_DP
		rdr = __llvm_amdgcn_rsq_f64(dr2);
#else
		rdr = __llvm_amdgcn_rsq_f32(dr2);//__frsqrt_rn(dr2);
#endif
#else
		rdr = rsqrt(dr2);
#endif
		dr = rdr * dr2;
	}
	else {
		dr = SoftenScale;
		rdr = 1.0 / SoftenScale;
	}
	gpuT drs = half * dr * rrs;

#ifdef GPUCUTOFF
//	if (drs > 3.0) {
	if (drs > 2.8125) {
		return acc;
	}
#endif

#ifdef SPI
	ir3 = MASSPART * (rdr * rdr * rdr);
#else
	ir3 = MASSPART / (dr * dr * dr);
#endif

#ifdef SPI
	gpuT t0, t1;
	gpuT widinv = NTABINT / 3.0;
	gpuT width = 3.0 / NTABINT;

#ifndef SPI_T2
	gpuT one_3 = 0.333333;
	gpuT five_3 = 1.66666667;
	gpuT two_3 = 0.6666667;
	gpuT one = 1.0;
	gpuT three_2 = 1.5;
#endif
	int idx = (int)(drs * widinv);
	gpuT xi = idx * width;
	gpuT eps = drs - xi;
	if (eps > half * width) {
		idx ++;
		xi = idx * width;
		eps = drs - xi;
	}
	if (idx >= NTABINT || idx <0) {
		t0 = 0.0;
		t1 = 0.0;
#ifdef GPUCUTOFF
		return acc;
#endif
	} else {
		t0 = tab0[idx];
		t1 = tab1[idx];
	}


#ifdef  SPI_T2
	gpuT a1 = xi* xi * eps;
	gpuT a2 = xi * eps * ( eps -  a1 );

	ir3 *= (t0 + (a1 + a2)* t1);
#else
	gpuT tdx = eps;
	gpuT x2  = xi * xi;
	gpuT a1  = x2 * tdx;
	tdx *= eps;
	gpuT a2  = xi * (one - x2) * tdx;
	tdx *= eps;
	gpuT a3  = (one_3 - five_3 * x2 + two_3 * x2 * x2) * tdx;
	tdx *= eps;
	gpuT a4  = xi * (three_2 * x2- one - one_3 * x2 * x2 ) * tdx; 

	ir3 *= (t0 + (a1 + a2 + a3 + a4)* t1);

#endif // T4

#else
	ir3 *= (erfc(drs) + coeff * drs * exp(-drs * drs));
#endif
	acc.x += dx.x * ir3;
	acc.y += dx.y * ir3;
	acc.z += dx.z * ir3;

	return acc;
}


#ifdef SPI
__device__ gpuT3 P2P_core_lite(gpuT SoftenScale, gpuT r_SoftenScale, gpuT rrs, gpuT coeff,
		gpuT3 ii, gpuT3 jj, gpuT3 acc, gpuT *tab0, gpuT *tab1)
#else
__device__ gpuT3 P2P_core_lite(gpuT SoftenScale, gpuT r_SoftenScale, gpuT rrs, gpuT coeff,
		gpuT3 ii, gpuT3 jj, gpuT3 acc)
#endif
{
	gpuT3 dx;
	dx.x = jj.x - ii.x;
	dx.y = jj.y - ii.y;
	dx.z = jj.z - ii.z;
	gpuT dr2 = dx.x * dx.x + dx.y * dx.y + dx.z * dx.z;

	gpuT dr = SoftenScale;
	gpuT rdr = r_SoftenScale;

	if (dr2 > SoftenScale * SoftenScale) {
#ifdef HIP_RTM
#ifdef GPU_DP
		rdr = __llvm_amdgcn_rsq_f64(dr2);
#else
		rdr = __llvm_amdgcn_rsq_f32(dr2);//__frsqrt_rn(dr2);
#endif
#else
		rdr = rsqrt(dr2);
#endif
		dr = rdr * dr2;
	}

	gpuT drs = dr * rrs;

#ifdef GPUCUTOFF
	if (drs > 2.8125) {
		return acc;
	}
#endif

#ifdef SPI
	gpuT ir3 = rdr * rdr * rdr;
#else
	gpuT ir3 = 1.0 / (dr * dr * dr);
#endif

#ifdef SPI
	int idx = (int)(drs * WIDTH_INV);
	gpuT xi = idx * WIDTH;
	gpuT eps = drs - xi;
	if (eps > HALF * WIDTH) {
		idx ++;
		xi = idx * WIDTH;
		eps = drs - xi;
	}
	gpuT t0 = 0., t1 = 0.;
	if (idx < NTABINT) {
		t0 = tab0[idx];
		t1 = tab1[idx];
	} 
#ifdef GPUCUTOFF
	else {
		return acc;
	}
#endif


#ifdef  SPI_T2
	gpuT a1 = xi* xi * eps;
	gpuT a2 = xi * eps * ( eps -  a1 );

	ir3 *= (t0 + (a1 + a2)* t1);
#else
	gpuT tdx = eps;
	gpuT x2  = xi * xi;
	gpuT a1  = x2 * tdx;
	tdx *= eps;
	gpuT a2 = xi * (ONE_1 - x2) * tdx;
	tdx *= eps;
	gpuT x4  = x2 * x2;
	gpuT a3 = (ONE_3 - FIVE_3 * x2 + TWO_3 * x4) * tdx;
	tdx *= eps;
	gpuT a4 = xi * (THREE_2 * x2- ONE_1 - ONE_3 * x4 ) * tdx; 

	ir3 *= (t0 + (a1 + a2 + a3 + a4)* t1);

#endif // T4

#else
	ir3 *= (erfc(drs) + coeff * drs * exp(-drs * drs));
#endif
//	acc[0] += dx.x * ir3;
//	acc[1] += dx.y * ir3;
//	acc[2] += dx.z * ir3;
	acc.x += dx.x * ir3;
	acc.y += dx.y * ir3;
	acc.z += dx.z * ir3;

	return acc;
}


#ifdef SPI
__device__ gpuT3 P2P_core_1(gpuT SoftenScale, gpuT MASSPART, gpuT rrs, gpuT coeff,
		gpuT3 ii, gpuT3 jj, gpuT3 acc, gpuT *tab0, gpuT *tab1)
#else
__device__ gpuT3 P2P_core_1(gpuT SoftenScale, gpuT MASSPART, gpuT rrs, gpuT coeff,
		gpuT3 ii, gpuT3 jj, gpuT3 acc)
#endif
{
	gpuT dr, ir3, half = 0.5;
	gpuT3 dx;
	dx.x = jj.x - ii.x;
	dx.y = jj.y - ii.y;
	dx.z = jj.z - ii.z;

#ifdef SPI
	gpuT rdr, dr2 = dx.x * dx.x + dx.y * dx.y + dx.z * dx.z + SoftenScale * SoftenScale;
	rdr = rsqrt(dr2);
	dr = rdr * dr2;
#else
	dr = sqrt(dx.x * dx.x + dx.y * dx.y + dx.z * dx.z + SoftenScale * SoftenScale);
#endif

	gpuT drs = half * dr * rrs;

#ifdef GPUCUTOFF
	if (drs > 2.8125) {
		return acc;
	}
#endif

#ifdef SPI
	ir3 = MASSPART * (rdr * rdr * rdr);
#else
	ir3 = MASSPART / (dr * dr * dr);
#endif

#ifdef SPI
	gpuT t0, t1;
	gpuT widinv = NTABINT / 3.0;
	gpuT width = 3.0 / NTABINT;

#ifndef SPI_T2
	gpuT one_3 = 0.333333;
	gpuT five_3 = 1.66666667;
	gpuT two_3 = 0.6666667;
	gpuT one = 1.0;
	gpuT three_2 = 1.5;
#endif

	int idx = (int)(drs * widinv);
	gpuT xi = idx * width;
	gpuT eps = drs - xi;
	if (eps > half * width) {
		idx ++;
		xi = idx * width;
		eps = drs - xi;
	}
	if (idx >= NTABINT) {
		t0 = 0.0;
		t1 = 0.0;
#ifdef GPUCUTOFF
		return acc;
#endif
	} else {
		t0 = tab0[idx];
		t1 = tab1[idx];
	}


#ifdef SPI_T2
	gpuT a1 = xi* xi * eps;
	gpuT a2 = xi * eps * ( eps -  a1 );

	ir3 *= (t0 + (a1 + a2)* t1);
#else

	gpuT tdx = eps;
	gpuT x2  = xi * xi;
	gpuT a1  = x2 * tdx;
	tdx *= eps;
	gpuT a2  = xi * (one - x2) * tdx;
	tdx *= eps;
	gpuT a3  = (one_3 - five_3 * x2 + two_3 * x2 * x2) * tdx;
	tdx *= eps;
	gpuT a4  = xi * (three_2 * x2- one - one_3 * x2 * x2 ) * tdx; 

	ir3 *= (t0 + (a1 + a2 + a3 + a4)* t1);

#endif // T4

#else
	ir3 *= (erfc(drs) + coeff * drs * exp(-drs * drs));
#endif
	if (drs > 2.8125) {
		ir3 = 0.0;
	}

	acc.x += dx.x * ir3;
	acc.y += dx.y * ir3;
	acc.z += dx.z * ir3;

	return acc;
}

__global__ void P2P_kernel(gpuT SoftenScale, gpuT MASSPART, gpuT rrs, gpuT coeff,
		int numInodes, int *dev_inodes, int *dev_in, int *dev_nn, int *dev_jnodes, 
		int ipart, int npart, gpuT *dev_part_pos, gpuT *dev_part_acc, int *dev_leaf_ip, int *dev_leaf_np, uint_8b *dev_part_act = NULL, uint_8b active_level = 0)
{
	extern __shared__ gpuT sh_jpos[];
#ifdef SPI
	__shared__ gpuT tab0[NTABINT]; 
	__shared__ gpuT tab1[NTABINT];
#endif
	int tid = threadIdx.x;
#ifdef SPI
	gpuT x, neg_two = -2.0;
	gpuT width = 3.0 / NTABINT;
	for (int i = tid; i < NTABINT; i += blockDim.x) {
		x = (gpuT)i * width;	
		tab0[i] = erfc(x) + x * exp(-x * x) * coeff; 
		tab1[i] = neg_two * coeff * exp(-x * x);
	}
	__syncthreads();
#endif
	int inode, in, nn, jnode;
	int iip, inp, jip, jnp;
	gpuT3 ipos, iacc, jpos;
	bool cmp;
	for (int i = blockIdx.x; i < numInodes; i += gridDim.x) {
		inode = dev_inodes[i];
		iip = dev_leaf_ip[inode] - ipart;
		inp = dev_leaf_np[inode];
#ifdef ACT_OPT
		cmp = (tid < inp && (!dev_part_act || dev_part_act[tid+iip] == active_level));
#else
		cmp = (tid < inp);
#endif
		if(cmp) {
			ipos.x = dev_part_pos[tid+iip];
			ipos.y = dev_part_pos[tid+iip+npart];
			ipos.z = dev_part_pos[tid+iip+2*npart];
			iacc.x = 0.0;
			iacc.y = 0.0;
			iacc.z = 0.0;
		}
		in = dev_in[i];
		nn = dev_nn[i];
		for (int j = in; j < in + nn; ++j) {
			jnode = dev_jnodes[j];
			jip = dev_leaf_ip[jnode] - ipart;
			jnp = dev_leaf_np[jnode];
			__syncthreads();
			if (tid < jnp) {
				sh_jpos[tid] = dev_part_pos[tid+jip];
				sh_jpos[tid+jnp] = dev_part_pos[tid+jip+npart];
				sh_jpos[tid+2*jnp] = dev_part_pos[tid+jip+2*npart];
			}
			__syncthreads();
			if(cmp) {
				for (int k = 0; k < jnp; ++k) {
					jpos.x = sh_jpos[k];
					jpos.y = sh_jpos[k+jnp];
					jpos.z = sh_jpos[k+2*jnp];
#ifdef SPI
					iacc = P2P_core(SoftenScale, MASSPART, rrs, coeff, ipos, jpos, iacc, tab0, tab1);
#else
					iacc = P2P_core(SoftenScale, MASSPART, rrs, coeff, ipos, jpos, iacc);
#endif
				}
			}
		}
		if(cmp) {
			dev_part_acc[tid+iip] += iacc.x;
			dev_part_acc[tid+iip+npart] += iacc.y;
			dev_part_acc[tid+iip+2*npart] += iacc.z;
		}
	}
}

__global__ void P2P_kernel3(gpuT SoftenScale, gpuT MASSPART, gpuT rrs, gpuT coeff,
		int numInodes, int *dev_inodes, int *dev_in, int *dev_nn, int *dev_jnodes, 
		int ipart, int npart, gpuT *dev_part_pos, gpuT *dev_part_acc, int *dev_leaf_ip, int *dev_leaf_np, uint_8b *dev_part_act = NULL, uint_8b active_level = 0)
{
	extern __shared__ gpuT sh_mem[];
#ifdef SPI
	__shared__ gpuT tab0[NTABINT]; 
	__shared__ gpuT tab1[NTABINT];
#endif
	unsigned int tid = threadIdx.x;
	gpuT (* volatile sh_jpos)[3] = (gpuT (*)[3])sh_mem;
//	gpuT (* volatile sh_ipos)[3] = (gpuT (*)[3])(sh_mem + bsize * 3);
//	gpuT (* volatile sh_iacc)[3] = (gpuT (*)[3])(sh_mem + bsize * 6);
#ifdef SPI
	gpuT x, neg_two = -2.0;
	for (int i = tid; i < NTABINT; i += blockDim.x) {
		x = (gpuT)i * WIDTH;	
		tab0[i] = erfc(x) + x * exp(-x * x) * coeff; 
		tab1[i] = neg_two * coeff * exp(-x * x);
	}
	__syncthreads();
#endif
	int inode, in, nn, jnode;
	int iip, inp, jip, jnp;
	bool cmp;
	gpuT3 ipos, jpos, iacc;
	gpuT hrrs = HALF * rrs, r_SoftenScale = 1.0 / SoftenScale;

	unsigned int bid = blockIdx.x;
	inode = dev_inodes[bid];
	iip = dev_leaf_ip[inode] - ipart;
	inp = dev_leaf_np[inode];
#ifdef ACT_OPT
	cmp = (tid < inp && (!dev_part_act || dev_part_act[tid+iip] == active_level));
#else
	cmp = (tid < inp);
#endif
	if(cmp) {
		ipos.x = dev_part_pos[tid+iip];
		ipos.y = dev_part_pos[tid+iip+npart];
		ipos.z = dev_part_pos[tid+iip+2*npart];
		iacc.x = 0.0;
		iacc.y = 0.0;
		iacc.z = 0.0;
	//	sh_ipos[tid][0] = dev_part_pos[tid+iip];
	//	sh_ipos[tid][1] = dev_part_pos[tid+iip+npart];
	//	sh_ipos[tid][2] = dev_part_pos[tid+iip+2*npart];
//		sh_iacc[tid][0] = 0.0;
//		sh_iacc[tid][1] = 0.0;
//		sh_iacc[tid][2] = 0.0;
	}
	in = dev_in[bid];
	nn = dev_nn[bid];
	for (int j = in; j < in + nn; ++j) {
		jnode = dev_jnodes[j];
		jip = dev_leaf_ip[jnode] - ipart;
		jnp = dev_leaf_np[jnode];
		__syncthreads();
		if (tid < jnp) {
			sh_jpos[tid][0] = dev_part_pos[tid+jip];
			sh_jpos[tid][1] = dev_part_pos[tid+jip+npart];
			sh_jpos[tid][2] = dev_part_pos[tid+jip+2*npart];
		}
		__syncthreads();
		if(cmp) {
			for (int k = 0; k < jnp; ++k) {
				jpos.x = sh_jpos[k][0];
				jpos.y = sh_jpos[k][1];
				jpos.z = sh_jpos[k][2];
#ifdef SPI
				iacc = P2P_core_lite(SoftenScale, r_SoftenScale, hrrs, coeff, ipos, jpos, iacc, tab0, tab1);
#else
				iacc = P2P_core_lite(SoftenScale, r_SoftenScale, hrrs, coeff, ipos, jpos, iacc);
#endif
			}
		}
	}
	if(cmp) {
		dev_part_acc[tid+iip] += iacc.x * MASSPART;
		dev_part_acc[tid+iip+npart] += iacc.y * MASSPART;
		dev_part_acc[tid+iip+2*npart] += iacc.z * MASSPART;
	//	dev_part_acc[tid+iip] += sh_iacc[tid][0] * MASSPART;
	//	dev_part_acc[tid+iip+npart] += sh_iacc[tid][1] * MASSPART;
	//	dev_part_acc[tid+iip+2*npart] += sh_iacc[tid][2] * MASSPART;
	}
}

inline __device__ gpuT warpReduce(int bs, gpuT t) {
#if CUDART_VERSION >= 9000
	if (bs >= 32) t += __shfl_xor_sync(0xffffffff,t,16);
	if (bs >= 16) t += __shfl_xor_sync(0xffffffff,t,8);
	if (bs >= 8) t += __shfl_xor_sync(0xffffffff,t,4);
	if (bs >= 4) t += __shfl_xor_sync(0xffffffff,t,2);
	if (bs >= 2) t += __shfl_xor_sync(0xffffffff,t,1);
#else
	if (bs >= 32) t += __shfl_xor(t,16);
	if (bs >= 16) t += __shfl_xor(t,8);
	if (bs >= 8) t += __shfl_xor(t,4);
	if (bs >= 4) t += __shfl_xor(t,2);
	if (bs >= 2) t += __shfl_xor(t,1);
#endif
	return t;
}

__device__ void warpReduceAndStore(int bs, gpuT * data, int tid, gpuT *result) {
	gpuT t = warpReduce(bs, data[tid]);
	if (tid==0) *result += t;
}

__device__ void warpReduceAndStore(int bs, int tid, gpuT t, gpuT *result) {
	t = warpReduce(bs, t);
	if (tid==0) *result = t;
}


__global__ void P2P_kernel2(gpuT SoftenScale, gpuT MASSPART, gpuT rrs, gpuT coeff, int maxleaf,
		int numInodes, int *dev_inodes, int *dev_in, int *dev_nn, int *dev_jnodes, 
		int ipart, int npart, gpuT *dev_part_pos, gpuT *dev_part_acc, int *dev_leaf_ip, int *dev_leaf_np, uint_8b *dev_part_act = NULL, uint_8b active_level = 0)
{
	extern __shared__ gpuT sh_ipos[];
	__shared__ gpuT sh_iacc[32 * 3]; // 1024 / 32 = 32
#ifdef SPI
	__shared__ gpuT tab0[NTABINT]; 
	__shared__ gpuT tab1[NTABINT]; 
#endif
	int tid = threadIdx.x;
#ifdef SPI
	gpuT x, neg_two = -2.0;
	gpuT width = 3.0 / NTABINT;
	for (int i = tid; i < NTABINT; i += blockDim.x) {
		x = (gpuT)i * width;	
		tab0[i] = erfc(x) + x * exp(-x * x) * coeff; 
		tab1[i] = neg_two * coeff * exp(-x * x);
	}
	__syncthreads();
#endif
	int inode, in, nn, jnode;
	int iip, inp, jip, jnp;
	gpuT3 ipos, iacc, jpos;
	bool cmp;
	int nWarp = blockDim.x / WARP;
	int iWarp = tid / WARP;
	int iInW = tid % WARP;
	int nLeaf = blockDim.x / maxleaf;
	int iLeaf = tid / maxleaf;
	int iInL = tid % maxleaf;
	for (int i = blockIdx.x; i < numInodes; i += gridDim.x) {
		inode = dev_inodes[i];
		iip = dev_leaf_ip[inode] - ipart;
		inp = dev_leaf_np[inode];
#ifdef ACT_OPT
		cmp = (!dev_part_act || dev_part_act[tid+iip] == active_level);
#else
		cmp = true;
#endif
		__syncthreads();
		if(tid < inp && cmp) {
			sh_ipos[tid] = dev_part_pos[tid+iip];
			sh_ipos[tid+inp] = dev_part_pos[tid+iip+npart];
			sh_ipos[tid+2*inp] = dev_part_pos[tid+iip+2*npart];
		}
		__syncthreads();
		in = dev_in[i];
		nn = dev_nn[i];
		for (int j = 0; j < inp; ++j) {
#ifdef ACT_OPT
			if (dev_part_act && dev_part_act[iip + j] != active_level) continue;
#endif
			ipos.x = sh_ipos[j];
			ipos.y = sh_ipos[j+inp];
			ipos.z = sh_ipos[j+2*inp];
			iacc.x = 0.0;
			iacc.y = 0.0;
			iacc.z = 0.0;
			for (int k = in; k < in + nn; k += nLeaf) {
				if (k + iLeaf >= in + nn) break;
				jnode = dev_jnodes[k+iLeaf];
				jip = dev_leaf_ip[jnode] - ipart;
				jnp = dev_leaf_np[jnode];
				if (iInL < jnp) {
					jpos.x = dev_part_pos[iInL+jip];
					jpos.y = dev_part_pos[iInL+jip+npart];
					jpos.z = dev_part_pos[iInL+jip+2*npart];
#ifdef SPI
					iacc = P2P_core(SoftenScale, MASSPART, rrs, coeff, ipos, jpos, iacc, tab0, tab1);
#else
					iacc = P2P_core(SoftenScale, MASSPART, rrs, coeff, ipos, jpos, iacc);
#endif
				}
			}
			__syncthreads();
			warpReduceAndStore(WARP, iInW, iacc.x, sh_iacc + iWarp);
			warpReduceAndStore(WARP, iInW, iacc.y, sh_iacc + iWarp + nWarp);
			warpReduceAndStore(WARP, iInW, iacc.z, sh_iacc + iWarp + nWarp*2);
			__syncthreads();
			warpReduceAndStore(nWarp, sh_iacc, 	     tid, dev_part_acc + iip + j);
			warpReduceAndStore(nWarp, sh_iacc + nWarp,   tid, dev_part_acc + iip + j + npart);
			warpReduceAndStore(nWarp, sh_iacc + nWarp*2, tid, dev_part_acc + iip + j + 2*npart);
		}
	}
}
