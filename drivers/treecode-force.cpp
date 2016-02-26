/*
 *  treecode-force.cpp
 *  Part of 2d-treecodes
 *
 *  Created and authored by Diego Rossinelli on 2015-09-25.
 *  Copyright 2015. All rights reserved.
 *
 *  Users are NOT authorized
 *  to employ the present software for their own publications
 *  before getting a written permission from the author of this file.
 */

#include <cassert>
#include <cmath>
#include <cstring>
#include <tuple>
#include <algorithm>
#include <omp.h>
#include "upward.h"
#include "force-kernels.h"
#include "treecode-force.h"

#define MIXPREC 
#define INSTRUMENTATION

#ifndef INSTRUMENTATION
#define kernelcall(w, f, ...) f(__VA_ARGS__)
#else
static __inline__ unsigned long long rdtsc(void)
{
    unsigned hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return ((unsigned long long)lo) | (((unsigned long long)hi) << 32);
}

#define kernelcall(w, f, ...)			\
    do						\
    {						\
	const int64_t start = rdtsc();		\
	f(__VA_ARGS__);				\
	const int64_t stop = rdtsc();		\
	td.f##_cycles += stop - start;		\
	td.f##_calls += w;			\
    }						\
    while(0)
#endif

namespace EvaluateForce
{
#ifdef MIXPREC
    float * xdataf32 = NULL, *ydataf32 = NULL, *vdataf32 = NULL;
#endif
    
    struct TimeDistrib
    {
#ifdef MIXPREC
	int64_t force_p2p_8x8f_cycles, force_p2p_8x8f_calls;
#else
	int64_t force_p2p_8x8_cycles, force_p2p_8x8_calls;
#endif
	int64_t
	    force_e2p_8x8_cycles, force_e2p_8x8_calls,
	    force_e2l_cycles, force_e2l_calls,
	    force_l2p_8x8_cycles, force_l2p_8x8_calls;

	void init()
	    {
#ifdef INSTRUMENTATION
		memset(this, 0, sizeof(*this));
#endif
	    }

	void report()
	    {
#ifdef INSTRUMENTATION
		const double tot = std::max((int64_t)1,
#ifdef MIXPREC
					    force_p2p_8x8f_cycles +
#else
					    force_p2p_8x8_cycles +
#endif
					    force_e2p_8x8_cycles +
					    force_e2l_cycles +
					    force_l2p_8x8_cycles );
		
		const double e2pflops =
		    force_e2p_8x8_calls

		    *
		    (2 * 8 + 8 * (2 + 8 * (2 + 1 + 2 + 2) + 8 * ORDER * (6 + 10) + 8 * 2));

		const double e2pFPC = e2pflops / force_e2p_8x8_cycles;

		const double e2lflops = force_e2l_calls *
		    (5 + ORDER * 12 + ORDER * (2 + ORDER * 4 + 6 + 2));

		const double e2lFPC = e2lflops / force_e2l_cycles;
		    
		const double p2pCPP =
#ifdef MIXPREC
		    force_p2p_8x8f_cycles / 64. / force_p2p_8x8f_calls;
#else
		    force_p2p_8x8_cycles / 64. / force_p2p_8x8_calls;
#endif
		   

		printf("total: %.1e Mcycles p2p:%.1f%% (%.1f C/P) e2p:%.1f%% (%.1f F/C) e2l:%.1f%% (%.1f F/C) l2p:%.1f%%\n",  
		       tot * 1e-6,
#ifdef MIXPREC
		       force_p2p_8x8f_cycles * 100. / tot, p2pCPP,
#else
		       force_p2p_8x8_cycles * 100. / tot, p2pCPP,
#endif
		       force_e2p_8x8_cycles * 100. / tot, e2pFPC, 
		       force_e2l_cycles * 100. / tot, e2lFPC,
		       force_l2p_8x8_cycles * 100. / tot);
#endif
	    }

    } td;
    
#pragma omp threadprivate(td)

    template<int size>
    struct E2LWork 
    {
	int count;
	realtype *const rdst,  *const idst;

	realtype x0s[size], y0s[size], masses[size];
	const realtype * rxps[size], *ixps[size];

	E2LWork(realtype * const rlocal, realtype * const ilocal):
	    count(0), rdst(rlocal), idst(ilocal) { }

	void _flush()
	    {
		kernelcall(count, force_e2l, x0s, y0s, masses, rxps, ixps, count, rdst, idst);
		
		count = 0;
	    }

	void push(const realtype x0, const realtype y0, const realtype mass,
		  const realtype * const rxp, const realtype * const ixp)
	    {
		x0s[count] = x0;
		y0s[count] = y0;
		masses[count] = mass;
		rxps[count] = rxp;
		ixps[count] = ixp;

		if (++count >= size)
		    _flush();
	    }

	void finalize()
	    {
		if (count)
		    _flush();
	    }
    };

#define TILE 8
#define BRICKSIZE (TILE * TILE)
    
    void evaluate(realtype * const xresultbase, realtype * const yresultbase,
		  const realtype x0, const realtype y0, const realtype h,
		  const realtype theta)
    {
	int maxentry = 0;

	int stack[LMAX * 3];

	realtype result[2 * BRICKSIZE];

#ifdef MIXPREC
	float resultf[2 * BRICKSIZE];
#endif

	realtype rlocal[ORDER + 1], ilocal[ORDER + 1];

	E2LWork<32> e2lwork(rlocal, ilocal);

	const realtype rbrick = 1.4142135623730951 * h * (TILE - 1) * 0.5;

	for(int by = 0; by < BLOCKSIZE; by += TILE)
	    for(int bx = 0; bx < BLOCKSIZE; bx += TILE)
	    {
		const realtype x0brick = x0 + h * (bx + 0.5 * (TILE - 1));
		const realtype y0brick = y0 + h * (by + 0.5 * (TILE - 1));

		for(int i = 0; i <= ORDER; ++i)
		{
		    rlocal[i] = 0;
		    ilocal[i] = 0;
		}

		for(int i = 0; i < 2 * BRICKSIZE; ++i)
		    result[i] = 0;

#ifdef MIXPREC
		for(int i = 0; i < 2 * BRICKSIZE; ++i)
		    resultf[i] = 0;
#endif

		int stackentry = 0;
		stack[0] = 0;

		while(stackentry > -1)
		{
		    const int nodeid = stack[stackentry--];
		    const Tree::Node * const node = Tree::nodes + nodeid;

		    const realtype xcom = node->xcom;
		    const realtype ycom = node->ycom;

		    const realtype distance = sqrt(pow(x0brick - xcom, 2) + pow(y0brick - ycom, 2));

		    const bool localexpansion_converges = (distance / node->r - 1) > (1 / theta) && rbrick <= node->r;

		    if (localexpansion_converges)
			e2lwork.push(xcom - x0brick, ycom - y0brick, node->mass,
				     Tree::expansions + ORDER * (2 * nodeid + 0),
				     Tree::expansions + ORDER * (2 * nodeid + 1));
		    else
		    {
			const double xt = std::max(x0 + bx * h, std::min(x0 + (bx + TILE - 1) * h, xcom));
			const double yt = std::max(y0 + by * h, std::min(y0 + (by + TILE - 1) * h, ycom));

			const realtype r2 = pow(xt - xcom, 2) + pow(yt - ycom, 2);

			if (node->r * node->r < theta * theta * r2)
			    kernelcall(1, force_e2p_8x8, node->mass, x0 + (bx + 0) * h - xcom, y0 + (by + 0) * h - ycom, h,
				       Tree::expansions + ORDER * (2 * nodeid + 0),
				       Tree::expansions + ORDER * (2 * nodeid + 1),
				       result, result + BRICKSIZE);		    
			else
			{
			    if (!node->state.innernode)
			    {
				const int s = node->s;
				const int e = node->e;
#ifdef MIXPREC
				kernelcall(e - s, force_p2p_8x8f, &xdataf32[s], &ydataf32[s], &vdataf32[s],
					   e - s, x0 + (bx + 0) * h, y0 + (by + 0) * h, h,
					   resultf, resultf + BRICKSIZE);
#else
				kernelcall(e - s, force_p2p_8x8, &Tree::xdata[s], &Tree::ydata[s], &Tree::vdata[s],
					   e - s, x0 + (bx + 0) * h, y0 + (by + 0) * h, h,
					   result, result + BRICKSIZE);
#endif
			    }
			    else
			    {
				for(int c = 0; c < 4; ++c)
				    stack[++stackentry] = node->state.childbase + c;

				maxentry = std::max(maxentry, stackentry);
			    }
			}
		    }
		}
		
		e2lwork.finalize();
		
		kernelcall(1, force_l2p_8x8, h * (0 - 0.5 * (TILE - 1)),
			   h * (0 - 0.5 * (TILE - 1)),
			   h, rlocal, ilocal, result, result + BRICKSIZE);
		
#ifdef MIXPREC
		for(int i = 0; i < 2 * BRICKSIZE; ++i)
		    result[i] += resultf[i];
#endif
		
		for(int iy = 0; iy < TILE; ++iy)
		    for(int ix = 0; ix < TILE; ++ix)
			xresultbase[bx + ix + BLOCKSIZE * (by + iy)] = result[ix + TILE * iy];

		for(int iy = 0; iy < TILE; ++iy)
		    for(int ix = 0; ix < TILE; ++ix)
			yresultbase[bx + ix + BLOCKSIZE * (by + iy)] = result[BRICKSIZE + ix + TILE * iy];
	    }
    }

    extern "C"
    __attribute__ ((visibility ("default")))
    void treecode_force_mrag(const realtype theta,
			     const realtype * const xsrc,
			     const realtype * const ysrc,
			     const realtype * const vsrc,
			     const int nsrc,
			     const realtype * const x0s,
			     const realtype * const y0s,
			     const realtype * const hs,
			     const int nblocks,
			     realtype * const xdst,
			     realtype * const ydst)
    {
	const double t0 =  omp_get_wtime();
	
	Tree::build(xsrc, ysrc, vsrc, nsrc,
#ifdef MIXPREC
		    256
#else
		    96
#endif
	    );

#ifdef MIXPREC
	posix_memalign((void **)&xdataf32, 32, sizeof(*xdataf32) * nsrc);
	posix_memalign((void **)&ydataf32, 32, sizeof(*ydataf32) * nsrc);
	posix_memalign((void **)&vdataf32, 32, sizeof(*vdataf32) * nsrc);

#pragma omp for
	for(int i = 0; i < nsrc; ++i)
	{
	    xdataf32[i] = (float)Tree::xdata[i];
	    ydataf32[i] = (float)Tree::ydata[i];
	    vdataf32[i] = (float)Tree::vdata[i];
	}
#endif
	
	const double t1 = omp_get_wtime();
#pragma omp parallel
	{
	    td.init();

#pragma omp for schedule(dynamic,1)	
	    for(int i = 0; i < nblocks; ++i)
		evaluate(xdst + i * BLOCKSIZE * BLOCKSIZE, ydst + i * BLOCKSIZE * BLOCKSIZE, x0s[i], y0s[i], hs[i], theta);

	    td.report();
	}

	const double t2 = omp_get_wtime();
	
	Tree::dispose();
#ifdef MIXPREC
	free(xdataf32);
	free(ydataf32);
	free(vdataf32);
#endif

#ifdef INSTRUMENTATION
      	printf("UP: %.1es DW: %.1es (%.1f%%)\n", t1- t0, t2-t1, (t2 - t1) * 100. / (t2 - t0));
#endif
    }
}
