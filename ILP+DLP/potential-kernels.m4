/*
*  potential-kernels.mc
*  Part of MRAG/2d-treecode-potential
*
*  Created and authored by Diego Rossinelli on 2015-11-25.
*  Copyright 2015. All rights reserved.
*
*  Users are NOT authorized
*  to employ the present software for their own publications
*  before getting a written permission from the author of this file.
*/

define(NACC, 4)
include(unroll.m4) dnl

#include <cassert>
#define EPS (10 * __DBL_EPSILON__)
#define MAX(a,b) (((a)>(b))?(a):(b))

__device__ realtype potential_p2p(
	   const realtype * __restrict__ const xsrc,
	   const realtype * __restrict__ const ysrc,
	   const realtype * __restrict__ const vsrc,
	   const int nsources,
	   const realtype xt,
	   const realtype yt)
  {
    const int tid = threadIdx.x;
    assert(tid < 32 && blockDim.x == 32);
    
    realtype LUNROLL(`i', 0, eval(NACC - 1), `ifelse(i,0,,`,') TMP(s,i) = 0') ;

    const int nnice = eval(32 * NACC) * (nsources / eval(32 * NACC));

    for(int i = 0; i < nnice; i += eval(NACC * 32))
    {dnl
      LUNROLL(j, 0, eval(NACC - 1), `
      const realtype TMP(xr, j) = xt - xsrc[tid + i + eval(32 * j)];
      const realtype TMP(yr, j) = yt - ysrc[tid + i + eval(32 * j)];')
      dnl
      LUNROLL(j, 0, eval(NACC - 1), `
      TMP(s, j) += log(TMP(xr, j) * TMP(xr, j) + TMP(yr, j) * TMP(yr, j) + EPS) * vsrc[tid + i + eval(32 * j)];')
    }dnl

    REDUCE(`+=', LUNROLL(i, 0, eval(NACC - 1),`ifelse(i,0,,`,')TMP(s,i)'))
    
    
    for(int i = nnice + tid; i < nsources; i += 32)
    {
      const realtype xr = xt - xsrc[i];
      const realtype yr = yt - ysrc[i];

      TMP(s, 0) += log(xr * xr + yr * yr + EPS) * vsrc[i];
    }

    //#SEQ(`
    //#TMP(s, 0) += __shfl_xor(TMP(s, 0), L );', L, 16, 8, 4, 2, 1)
    
    return TMP(s, 0) / 2;
  }

define(`RARY', `scratch[$1]')
define(`IARY', `scratch[32 + $1]')

define(PASSES, `esyscmd(python z_to_k.py ORDER)')

__device__ realtype potential_e2p(const realtype mass,
  const realtype rz,
  const realtype iz,
  const realtype * __restrict__ const rxp,
  const realtype * __restrict__ const ixp,
  realtype * const scratch) // size of 32 * 2 * sizeof(realtype)
  {
    //volatile realtype * const scratch = _scratch;
    const int tid = threadIdx.x;
    assert(tid < 32 && blockDim.x == 32);
  
    const realtype r2 = rz * rz + iz * iz;

    scratch[0] = rz / r2;
    scratch[32] = -iz / r2;

     RARY(1) = RARY(0) * RARY(0) - IARY(0) * IARY(0);
     IARY(1) = 2 * IARY(0) * RARY(0);

    SEQ(`
    pushdef(`DELTAN', eval(n / 2))dnl
    pushdef(`NOLD', eval(n / 2))dnl

    if (tid < DELTAN)
    {
        RARY(tid + NOLD) = RARY(tid) * RARY(NOLD - 1) - IARY(tid) * IARY(NOLD - 1);
	IARY(tid + NOLD) = RARY(tid) * IARY(NOLD - 1) + IARY(tid) * RARY(NOLD - 1);
    }', n, PASSES)       

    realtype rsum = 0;
    if (tid < ORDER)
       rsum = rxp[tid] * RARY(tid) - ixp[tid] * IARY(tid);

    SEQ(`rsum += __shfl_xor(rsum, L );
    ', L, 16, 8, 4, 2, 1)
    
    return  mass * log(r2) / 2 + rsum;
  }

  __device__ realtype potential_e2p_individual(const realtype mass,
  const realtype rz, 
  const realtype iz,
  const realtype * __restrict__ const rxp,
  const realtype * __restrict__ const ixp)
  {
    const int tid = threadIdx.x;
    //assert(tid == 0);
  
    const realtype r2 = rz * rz + iz * iz;

    const realtype rinvz_1 = rz / r2;
    const realtype iinvz_1 = -iz / r2;

    LUNROLL(j, 2, ORDER, `
    const realtype TMP(rinvz, j) = TMP(rinvz, eval(j - 1)) * rinvz_1 - TMP(iinvz, eval(j - 1)) * iinvz_1;
    const realtype TMP(iinvz, j) = TMP(rinvz, eval(j - 1)) * iinvz_1 + TMP(iinvz, eval(j - 1)) * rinvz_1;')

    LUNROLL(j, 1, ORDER, `
    realtype TMP(rsum, eval(j - 1)) = rxp[eval(j - 1)] * TMP(rinvz, j) - ixp[eval(j - 1)] * TMP(iinvz, j);')

    REDUCE(`+=', LUNROLL(i, 0, eval(ORDER - 1),`ifelse(i,0,,`,')TMP(rsum,i)'))
        
    return  mass * log(r2) / 2 + TMP(rsum, 0);
  }
