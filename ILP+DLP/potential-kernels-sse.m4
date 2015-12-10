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

include(unroll.m4)
dnl
#include <emmintrin.h>
#include <math.h>

#if defined(__INTEL_COMPILER)
inline __m128d operator+(__m128d a, __m128d b){ return _mm_add_pd(a, b); }
inline __m128d operator/(__m128d a, __m128d b){ return _mm_div_pd(a, b); }
inline __m128d operator*(__m128d a, __m128d b){ return _mm_mul_pd(a, b); }
inline __m128d operator-(__m128d a, __m128d b){ return _mm_sub_pd(a, b); }
inline __m128d operator += (__m128d& a, __m128d b){ return a = _mm_add_pd(a, b); }
inline __m128d operator -= (__m128d& a, __m128d b){ return a = _mm_sub_pd(a, b); }
#endif


define(NACC, 32)
#define EPS (10 * __DBL_EPSILON__)
#define MAX(a,b) (((a)>(b))?(a):(b))

#ifdef __cplusplus
extern "C"
#endif
realtype potential_p2p(const realtype * __restrict__ const _xsrc,
  const realtype * __restrict__ const _ysrc,
  const realtype * __restrict__ const _vsrc,
  const int nsources,
  const realtype xt,
  const realtype yt)
  {
    realtype dummy LUNROLL(`i', 0, NACC, `, TMP(s,i) = 0') ;

    const int nnice = NACC * (nsources / NACC);

    for(int i = 0; i < nnice; i += NACC)
    {
      const realtype * __restrict__ const xsrc = _xsrc + i;
      const realtype * __restrict__ const ysrc = _ysrc + i;
      const realtype * __restrict__ const vsrc = _vsrc + i;

      LUNROLL(j, 0, eval(NACC - 1), `
      const realtype TMP(xr, j) = xt - xsrc[j];')
      LUNROLL(j, 0, eval(NACC - 1), `
      const realtype TMP(yr, j) = yt - ysrc[j];')
      LUNROLL(j, 0, eval(NACC - 1), `
      TMP(s, j) += log(TMP(xr, j) * TMP(xr, j) + TMP(yr, j) * TMP(yr, j) + EPS) * vsrc[j];')
    }

    realtype sum = 0;

    for(int i = nnice; i < nsources; ++i)
    {
      const realtype xr = xt - _xsrc[i];
      const realtype yr = yt - _ysrc[i];

      sum += log(xr * xr + yr * yr + EPS) * _vsrc[i];
    }

    LUNROLL(i, 0, eval(NACC - 1), `sum += TMP(s, i);')

    return sum / 2;
  }

#ifdef __cplusplus
extern "C"
#endif
realtype potential_e2p(const realtype mass,
  const realtype rz,
  const realtype iz,
  const realtype * __restrict__ const rxp,
  const realtype * __restrict__ const ixp)
  {
    const realtype r2 = rz * rz + iz * iz;

    const realtype rinvz_1 = rz / r2;
    const realtype iinvz_1 = -iz / r2;

    const realtype TMP(rinvz, 2) = TMP(rinvz, 1) * rinvz_1 - TMP(iinvz, 1) * iinvz_1;
    const realtype TMP(iinvz, 2) = TMP(rinvz, 1) * iinvz_1 + TMP(iinvz, 1) * rinvz_1;

    const __m128d rz2 = _mm_set1_pd(TMP(rinvz, 2));
    const __m128d iz2 = _mm_set1_pd(TMP(iinvz, 2));

    const __m128d rbase_0 = _mm_set_pd(TMP(rinvz, 2), TMP(rinvz, 1));
    const __m128d ibase_0 = _mm_set_pd(TMP(iinvz, 2), TMP(iinvz, 1));

    LUNROLL(j, 0, eval(ORDER/2 - 1), `
    	       ifelse(j, eval(ORDER/2 - 1), , 
	       const __m128d TMP(rbase, eval(j + 1)) = TMP(rbase, j) * rz2 - TMP(ibase, j) * iz2;)
    	       
	       ifelse(j, eval(ORDER/2 - 1), , 
	       const __m128d TMP(ibase, eval(j + 1)) = TMP(rbase, j) * iz2 + TMP(ibase, j) * rz2;)

    	       const __m128d TMP(rxp2, j) = _mm_loadu_pd(rxp + eval(2 * j));
    	       const __m128d TMP(ixp2, j) = _mm_loadu_pd(ixp + eval(2 * j));

    	       const __m128d TMP(rsum, j) = (TMP(rxp2, j) * TMP(rbase, j) - TMP(ixp2, j) * TMP(ibase, j)) 
	       	     	     	       ifelse(j, 0, ,+ TMP(rsum, eval(j - 1)));
    ')

    double partial;
    _mm_store_sd(&partial,  _mm_hadd_pd(TMP(rsum, eval(ORDER/2 - 1)), TMP(rsum, eval(ORDER/2 - 1))));

    return  mass * log(r2) / 2 + partial;
  }

