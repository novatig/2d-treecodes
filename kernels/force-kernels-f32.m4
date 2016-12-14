include(unroll.m4)

export void force_p2p_8x8f(
			 uniform const float xsources[],
			 uniform const float ysources[],
			 uniform const float vsources[],
			 uniform const int nsources,
			 uniform const float x0,
			 uniform const float y0,
			 uniform const float h,
			 uniform float xresult[],
			 uniform float yresult[])
{
	uniform const double eps = 10 * __DBL_EPSILON__;
	
	const float xt = x0 + programIndex * h;
	
	LUNROLL(iy, 0, 7,`
	uniform const float TMP(yt, iy) = y0 + iy * h;') 

	LUNROLL(iy, 0, 7,`
	float TMP(xsum, iy) = 0, TMP(ysum, iy) = 0;')

	for(uniform int s = 0; s < nsources; ++s)
	{
		uniform const float xsrc = xsources[s];
		uniform const float ysrc = ysources[s];
		uniform const float vsrc = vsources[s];
		
		const float xr = xt - xsrc;
		const float xr2 = xr * xr;
    		
		LUNROLL(iy, 0, 7, `
		uniform const float TMP(yr, iy) = TMP(yt, iy) - ysrc;
		uniform const float TMP(yr2, iy) = TMP(yr, iy) * TMP(yr, iy);')

		LUNROLL(iy, 0, 7, `
    const float TMP(denom, iy) = xr2 + TMP(yr2, iy);
    const float TMP(factor, iy) = TMP(denom, iy) ? vsrc/ TMP(denom, iy) : 0;
		TMP(xsum, iy) += xr * TMP(factor, iy);
    		TMP(ysum, iy) += TMP(yr, iy) * TMP(factor, iy);')
	}
	
	LUNROLL(iy, 0, 7, `
	xresult[programIndex + eval(8 * iy)] += TMP(xsum, iy);
    	yresult[programIndex + eval(8 * iy)] += TMP(ysum, iy);')
}
