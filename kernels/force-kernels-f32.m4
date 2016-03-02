include(unroll.m4)

extern "C" void force_p2p_8x8f(
			  const float xsources[],
			  const float ysources[],
			  const float vsources[],
			  const int nsources,
			  const float x0,
			  const float y0,
			  const float h,
			  float xresult[],
			  float yresult[])
{
	const double eps = 10 * __DBL_EPSILON__;

	for(int d = 0; d < 64 ; ++d)
	{
		const int ix = d & 7;
		const int iy = d >> 3;

		const float xt = x0 + ix * h;
		const float yt = y0 + iy * h;

		float xsum = 0, ysum = 0;
		for(int s = 0; s < nsources; ++s)
		{
			const float xr = xt - xsources[s];
	    		const float yr = yt - ysources[s];
	    		const float factor = vsources[s] / (xr * xr + yr * yr + eps);

			xsum += xr * factor;
	    		ysum += yr * factor;
		}

		xresult[d] += xsum;
		yresult[d] += ysum;
	}
}
