/*
 *  force-kernels.m4, force-kernels.ispc
 *  Part of 2d-treecodes
 *
 *  Created and authored by Diego Rossinelli on 2015-11-25.
 *  Copyright 2015. All rights reserved.
 *
 *  Users are NOT authorized
 *  to employ the present software for their own publications
 *  before getting a written permission from the author of this file.
 */

include(unroll.m4)
#define EPS (2 * __DBL_EPSILON__)
export void force_p2p_8x8(
			 uniform const realtype xsources[],
			 uniform const realtype ysources[],
			 uniform const realtype vsources[],
			 uniform const int nsources,
			 uniform const realtype x0,
			 uniform const realtype y0,
			 uniform const realtype h,
			 uniform realtype xresult[],
			 uniform realtype yresult[])
{
	uniform const double eps = EPS;

	LUNROLL(ix, 0, 1, `
	const realtype TMP(xt, ix) = x0 + (eval(ix * 4) + programIndex) * h;')
	
	LUNROLL(iy, 0, 7,`
	uniform const realtype TMP(yt, iy) = y0 + iy * h;') 

	for(uniform int s = 0; s < nsources; ++s)
	{
		uniform const realtype xsrc = xsources[s];
		uniform const realtype ysrc = ysources[s];
		uniform const realtype vsrc = vsources[s];
		
		LUNROLL(ix, 0, 1, `
		const realtype TMP(xr, ix) = TMP(xt, ix) - xsrc;
		const realtype TMP(xr2, ix) = TMP(xr, ix) * TMP(xr, ix);')
    		
		LUNROLL(iy, 0, 7, `
		uniform const realtype TMP(yr, iy) = TMP(yt, iy) - ysrc;
		uniform const realtype TMP(yr2, iy) = TMP(yr, iy) * TMP(yr, iy) + eps;')

		LUNROLL(iy, 0, 7, `
		LUNROLL(ix, 0, 1, `
    		const realtype TMP(factor, eval(4 * ix + 8 * iy)) = vsrc / (TMP(xr2, ix) + TMP(yr2, iy));
		xresult[programIndex + eval(4 * ix + 8 * iy)] += TMP(xr, ix) * TMP(factor, eval(4 * ix + 8 * iy));
    		yresult[programIndex + eval(4 * ix + 8 * iy)] += TMP(yr, iy) * TMP(factor, eval(4 * ix + 8 * iy));
		')')
	}
}

export void force_e2p_8x8(
			uniform const realtype mass,
			uniform const realtype x0,
			uniform const realtype y0,
			uniform const realtype h,
			uniform const realtype rxp[],
			uniform const realtype ixp[],
			uniform realtype xresult[],
			uniform realtype yresult[])
{
	LUNROLL(ix, 0, 1, `
	const realtype TMP(rz, ix) = x0 + (eval(ix * 4) + programIndex) * h;')
	
	for(uniform int iy = 0; iy < 8; ++iy)
	{
		uniform const realtype iz = y0 + iy * h;

		LUNROLL(ix, 0, 1, `
		const realtype TMP(r2, ix) = TMP(rz, ix) * TMP(rz, ix) + iz * iz;')
		LUNROLL(ix, 0, 1, `
		const realtype TMP(rinvz_1, ix) = TMP(rz, ix) / TMP(r2, ix);')
		LUNROLL(ix, 0, 1, `
		const realtype TMP(iinvz_1, ix) = -iz / TMP(r2, ix);')

		LUNROLL(ix, 0, 1, `
		realtype TMP(rsum, ix) = mass * TMP(rinvz_1, ix), TMP(isum, ix) = mass * TMP(iinvz_1, ix);')

		LUNROLL(ix, 0, 1, `
		realtype TMP(rprod, ix) = TMP(rinvz_1, ix), TMP(iprod, ix) = TMP(iinvz_1, ix);')

		LUNROLL(j, 0, eval(ORDER - 1),`
		{
			LUNROLL(ix, 0, 1, `
			const realtype TMP(rtmp, ix) = TMP(rprod, ix) * TMP(rinvz_1, ix) - TMP(iprod, ix) * TMP(iinvz_1, ix);
	    		const realtype TMP(itmp, ix) = TMP(rprod, ix) * TMP(iinvz_1, ix) + TMP(iprod, ix) * TMP(rinvz_1, ix);')

			LUNROLL(ix, 0, 1, `
			TMP(rprod, ix) = TMP(rtmp, ix);
	    		TMP(iprod, ix) = TMP(itmp, ix);')

			uniform const realtype prefactor = eval(j + 1);
			uniform const uniform realtype rrxp = rxp[j];
			uniform const uniform realtype iixp = ixp[j];
		
			LUNROLL(ix, 0, 1, `
			TMP(rsum, ix) -= prefactor * (rrxp * TMP(rprod, ix) - iixp * TMP(iprod, ix));
	    		TMP(isum, ix) -= prefactor * (rrxp * TMP(iprod, ix) + iixp * TMP(rprod, ix));')
		}')

		LUNROLL(ix, 0, 1, `
	    	xresult[8 * iy + eval(ix * 4) + programIndex] += TMP(rsum, ix);
		yresult[8 * iy + eval(ix * 4) + programIndex] -= TMP(isum, ix);')
	}
}

define(mysign, `ifelse(eval((-1)**($1)), -1,-,+)')

export void force_e2l(
	uniform const realtype x0s[],
	uniform const realtype y0s[],
	uniform const realtype masses[],
	uniform const realtype * uniform rxps[],
	uniform const realtype * uniform ixps[],
	uniform const int nexpansions,
	uniform realtype rlocal[],
	uniform realtype ilocal[])
{
	uniform const double eps = EPS;

	foreach(i = 0 ... nexpansions)
	{
		const realtype mass = masses[i];

		const realtype x0 = x0s[i];
		const realtype y0 = y0s[i];

		const realtype r2z0 = x0 * x0 + y0 * y0 + eps;
    		const realtype rinvz_1 = x0 / r2z0;
    		const realtype iinvz_1 = -y0 / r2z0; //5 FLOPs

		uniform const realtype * const rxp = rxps[i];
		uniform const realtype * const ixp = ixps[i];

		dnl
    		LUNROLL(j, 1, eval(ORDER),`
    		ifelse(j, 1, , `
      	  		  const realtype TMP(rinvz, j) = TMP(rinvz, eval(j - 1)) * rinvz_1 - TMP(iinvz, eval(j - 1)) * iinvz_1;
      	  		  const realtype TMP(iinvz, j) = TMP(rinvz, eval(j - 1)) * iinvz_1 + TMP(iinvz, eval(j - 1)) * rinvz_1;')

	  		  const realtype TMP(rcoeff, j) = rxp[eval(j - 1)] * TMP(rinvz, j) - ixp[eval(j - 1)] * TMP(iinvz, j);
      	  		  const realtype TMP(icoeff, j) = rxp[eval(j - 1)] * TMP(iinvz, j) + ixp[eval(j - 1)] * TMP(rinvz, j); //12 FLOPs
      		')//12 FLOPs * ORDER

		LUNROLL(l, 1, eval(ORDER),`
      		{
			realtype TMP(rtmp, l) = ifelse(l,1,` - mass', `mass * esyscmd(echo -1/eval(l) | bc --mathlib )');
			realtype TMP(itmp, l) = 0;// 2 FLOPs

			pushdef(`BINFAC', `BINOMIAL(eval(l + k - 1), eval(k - 1)).')dnl
			 LUNROLL(k, 1, eval(ORDER),`
			 TMP(rtmp, l) mysign(k)= ifelse(BINFAC,1.f,,`BINFAC *') TMP(rcoeff, k);
			 TMP(itmp, l) mysign(k)= ifelse(BINFAC,1.f,,`BINFAC *') TMP(icoeff, k);
			 ') //4 FLOPs * ORDER
			popdef(`BINFAC')dnl

       			realtype rpartial = TMP(rtmp, l) * TMP(rinvz, l) - TMP(itmp, l) * TMP(iinvz, l);
       			realtype ipartial = TMP(rtmp, l) * TMP(iinvz, l) + TMP(itmp, l) * TMP(rinvz, l); //6 FLOPs

			rlocal[l] += reduce_add(rpartial);
			ilocal[l] += reduce_add(ipartial); //2 FLOPs
		}')
	}

	//TOTAL FLOP count: nexpansions * eval(5 + ORDER * 12 + ORDER * (2 + ORDER * 4 + 6 + 2))
}

export void force_l2p_8x8(
       uniform const realtype x0,
       uniform const realtype y0,
       uniform const realtype h,
       uniform const realtype rlocal[],
       uniform const realtype ilocal[],
       uniform realtype xresult[],
       uniform realtype yresult[])
{
	foreach(d = 0 ... 64)
	{
		const int ix = d & 7;
		const int iy = d >> 3;

		const realtype rz_1 = x0 + ix * h;
		const realtype iz_1 = y0 + iy * h;

		const realtype TMP(rresult, 1) = rlocal[1];
          	const realtype TMP(iresult, 1) = ilocal[1];

		LUNROLL(l, 2, eval(ORDER),`
          	const realtype TMP(rz, l) = TMP(rz, eval(l - 1)) * rz_1 - TMP(iz, eval(l - 1)) * iz_1;
          	const realtype TMP(iz, l) = TMP(rz, eval(l - 1)) * iz_1 + TMP(iz, eval(l - 1)) * rz_1;

          	const realtype TMP(rresult, l) = TMP(rresult, eval(l - 1)) +
          	l * (rlocal[l] * TMP(rz, eval(l - 1)) - ilocal[l] * TMP(iz, eval(l - 1)));

          	const realtype TMP(iresult, l) = TMP(iresult, eval(l - 1)) +
          	l * (rlocal[l] * TMP(iz, eval(l - 1)) + ilocal[l] * TMP(rz, eval(l - 1)));
          	')

		xresult[d] += TMP(rresult, ORDER);
		yresult[d] -= TMP(iresult, ORDER);
	}
}
