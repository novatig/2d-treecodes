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
divert(-1)

define(`forloop',
`pushdef(`$1', `$2')_forloop(`$1', `$2', `$3', `$4')popdef(`$1')')

define(`_forloop',
`$4`'ifelse($1, `$3', ,
`define(`$1', incr($1))_forloop(`$1', `$2', `$3', `$4')')')

define(`forrloop',
`pushdef(`$1', `$2')_forrloop(`$1', `$2', `$3', `$4')popdef(`$1')')

define(`_forrloop',
`$4`'ifelse($1, `$3', ,
`define(`$1', decr($1))_forrloop(`$1', `$2', `$3', `$4')')')

define(BINOMIAL, `syscmd(python binomial.py $1 $2)')

#USAGE LUNROLL
#$1 iteration variable
#$2 iteration start
#$3 iteration end
#$4 body

define(LUNROLL, `forloop($1, $2, $3,`$4')')
define(RLUNROLL, `forrloop($1, $2, $3, `$4')')
define(NACC, 32)
define(`TMP', $1_$2)
divert(0)
#include <math.h>
#include <immintrin.h>
#include "potential-kernels.h"

#define EPS (10 * __DBL_EPSILON__)
#define MAX(a,b) (((a)>(b))?(a):(b))

realtype treecode_p2p(const realtype * __restrict__ const _xsrc,
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

  realtype treecode_e2p(const realtype mass,
    const realtype rz,
    const realtype iz,
    const realtype * __restrict__ const rxp,
    const realtype * __restrict__ const ixp)
    {
      const realtype r2 = rz * rz + iz * iz;

      const realtype rinvz = rz / r2;
      const realtype iinvz = -iz / r2;

      const realtype rprod_0 = rinvz;
      const realtype iprod_0 = iinvz;

      realtype rs = mass * log(r2) / 2;

      rs += rprod_0 * rxp[0] - iprod_0 * ixp[0];

      LUNROLL(n, 1, eval(ORDER - 1), `
      const realtype TMP(rprod, n) = rinvz * TMP(rprod, eval(n - 1)) - iinvz * TMP(iprod, eval(n - 1));
      const realtype TMP(iprod, n) = iinvz * TMP(rprod, eval(n - 1)) + rinvz * TMP(iprod, eval(n - 1));

      rs += TMP(rprod, n) * rxp[n] - TMP(iprod, n) * ixp[n];
      ')

      return rs;
    }
