/*
 *  main.cpp
 *  Part of MRAG/2d-treecode-potential
 *
 *  Created and authored by Diego Rossinelli on 2015-09-25.
 *  Copyright 2015. All rights reserved.
 *
 *  Users are NOT authorized
 *  to employ the present software for their own publications
 *  before getting a written permission from the author of this file.
 */

#include <unistd.h>

#include <omp.h>

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <cstring>

#include <algorithm>
#include <limits>

#include "treecode.h"

void test_potential(realtype theta, double tol, FILE * f = NULL, bool verify = true)
{
    const realtype eps = std::numeric_limits<realtype>::epsilon() * 10;
    
    int NSRC = 1e5 * (0.1 + drand48());
    
    if (f)
	fread(&NSRC, sizeof(int), 1, f);
        
    realtype * xsrc = new realtype[NSRC];
    realtype * ysrc = new realtype[NSRC];
    realtype * sources = new realtype[NSRC];

    if (f)
    {
	fread(xsrc, sizeof(realtype), NSRC, f);
	fread(ysrc, sizeof(realtype), NSRC, f);
	fread(sources, sizeof(realtype), NSRC, f);
    }
    else
	for(int i = 0; i < NSRC; ++i)
	{
	    xsrc[i] = (drand48() - 0.25) * 10;
	    ysrc[i] = (drand48() - 0.75) * 20;
	    sources[i] = drand48() * 55;
	}
	
    int NDST = 3000 * (0.1 + drand48());
    
    if (f)
	fread(&NDST, sizeof(int), 1, f);
    
    realtype * xdst = new realtype[NDST];
    realtype * ydst = new realtype[NDST];
    realtype * ref = new realtype[NDST];

    if (f)
    {
	fread(xdst, sizeof(realtype), NDST, f);
	fread(ydst, sizeof(realtype), NDST, f);
	fread(ref, sizeof(realtype), NDST, f);

	for(int i = 0; i < NDST; ++i)
	    ref[i] *= 2.0 * M_PI;
    }
    else
    {
	for(int i = 0; i < NDST; ++i)
	{
	    xdst[i] = (drand48() - 0.25) * 10;
	    ydst[i] = (drand48() - 0.75) * 20;
	}
    }

#if 0 //while waiting clarifications from Sid, i overwrite his reference data with mine
    if (!f)
#endif
	if (verify)
#pragma omp parallel for
	    for(int i = 0; i < NDST; ++i)
	    {
		const realtype xd = xdst[i];
		const realtype yd = ydst[i];
	
		realtype s = 0;
	
		for(int j = 0; j < NSRC; ++j)
		{
		    const realtype xr = xd - xsrc[j];
		    const realtype yr = yd - ysrc[j];
		    const realtype r = sqrt(xr * xr + yr * yr + eps);
		    s += log(r) * sources[j]; 
		}
	
		ref[i] = s;
	    }

    realtype * targets = new realtype[NDST];

    printf("Testing POTENTIAL with %d sources and %d targets (theta %.3e)...\n", NDST, NSRC, theta);
    const double tstart = omp_get_wtime();
    treecode_potential(theta, xsrc, ysrc, sources, NSRC, xdst, ydst, NDST, targets);
    const double tend = omp_get_wtime();

    printf("\x1b[92msolved in %.2f ms\x1b[0m\n", (tend - tstart) * 1e3);
    
    if (verify)
    {
	double linf = 0, l1 = 0, linf_rel = 0, l1_rel = 0;
	for(int i = 0; i < NDST; ++i)
	{
	    assert(!std::isnan(ref[i]));
	    assert(!std::isnan(targets[i]));
	    
	    const double err = ref[i] - targets[i];
	    const double maxval = std::max(fabs(targets[i]), fabs(ref[i]));
	    const double relerr = err/std::max(1e-6, maxval); 
	    
	    if (fabs(relerr) >= tol && fabs(err) >= tol)
		printf("%d: %e ref: %e -> %e %e\n", i, targets[i], ref[i], err, relerr);
	    
	    assert(fabs(relerr) < tol || fabs(err) < tol);

	    l1 += fabs(err);
	    l1_rel += fabs(relerr);

	    linf = std::max(linf, fabs(err));
	    linf_rel = std::max(linf_rel, fabs(relerr));
	}
	
	printf("l-infinity errors: %.03e (absolute) %.03e (relative)\n", linf, linf_rel);
	printf("       l-1 errors: %.03e (absolute) %.03e (relative)\n", l1, l1_rel);
    }
	
    delete [] xdst;
    delete [] ydst;
    delete [] targets;
    delete [] ref;
	
    delete [] xsrc;
    delete [] ysrc;
    delete [] sources;

    printf("TEST PASSED.\n");
}

void test_force(realtype theta, double tol, FILE * f = NULL, bool verify = true)
{
    const realtype eps = std::numeric_limits<realtype>::epsilon() * 10;
    
    int NSRC = 1e5 * (0.1 + drand48());
    
    if (f)
	fread(&NSRC, sizeof(int), 1, f);
        
    realtype * xsrc = new realtype[NSRC];
    realtype * ysrc = new realtype[NSRC];
    realtype * sources = new realtype[NSRC];

    if (f)
    {
	fread(xsrc, sizeof(realtype), NSRC, f);
	fread(ysrc, sizeof(realtype), NSRC, f);
	fread(sources, sizeof(realtype), NSRC, f);
    }
    else
	for(int i = 0; i < NSRC; ++i)
	{
	    xsrc[i] = (drand48() - 0.25) * 10;
	    ysrc[i] = (drand48() - 0.75) * 20;
	    sources[i] = drand48() * 55;
	}
	
    int NDST = 3000 * (0.1 + drand48());
    
    if (f)
	fread(&NDST, sizeof(int), 1, f);
    
    realtype * xdst = new realtype[NDST];
    realtype * ydst = new realtype[NDST];
    realtype * xref = new realtype[NDST];
    realtype * yref = new realtype[NDST];

    if (f)
    {
	fread(xdst, sizeof(realtype), NDST, f);
	fread(ydst, sizeof(realtype), NDST, f);
	//fread(ref, sizeof(realtype), NDST, f);

	//for(int i = 0; i < NDST; ++i)
	//  ref[i] *= 2.0 * M_PI;
    }
    else
    {
	for(int i = 0; i < NDST; ++i)
	{
	    xdst[i] = (drand48() - 0.25) * 10;
	    ydst[i] = (drand48() - 0.75) * 20;
	}
    }

#if 0 //while waiting clarifications from Sid, i overwrite his reference data with mine
    if (!f)
#endif
	if (verify)
#pragma omp parallel for
	    for(int i = 0; i < NDST; ++i)
	    {
		const realtype xd = xdst[i];
		const realtype yd = ydst[i];
	
		realtype xs = 0, ys = 0;
	
		for(int j = 0; j < NSRC; ++j)
		{
		    const realtype xr = xd - xsrc[j];
		    const realtype yr = yd - ysrc[j];
		    const realtype factor = sources[j] / (xr * xr + yr * yr + eps);
		    xs += xr * factor;
		    ys += yr * factor;
		}
	
		xref[i] = xs;
		yref[i] = ys;
	    }

    realtype * xtargets = new realtype[NDST];
    realtype * ytargets = new realtype[NDST];

    printf("Testing FORCE with %d sources and %d targets (theta %.3e)...(tol is %e)\n", NDST, NSRC, theta, tol);
    const double tstart = omp_get_wtime();
    treecode_force(theta, xsrc, ysrc, sources, NSRC, xdst, ydst, NDST, xtargets, ytargets);
    const double tend = omp_get_wtime();

    printf("\x1b[94msolved in %.2f ms\x1b[0m\n", (tend - tstart) * 1e3);
    
    if (verify)
    {
	auto check = [&](const double * ref, const double * res)
	    {
		double linf = 0, l1 = 0, linf_rel = 0, l1_rel = 0;
		for(int i = 0; i < NDST; ++i)
		{
		    assert(!std::isnan(ref[i]));
		    assert(!std::isnan(res[i]));
		    
		    const double err = ref[i] - res[i];
		    const double maxval = std::max(fabs(res[i]), fabs(ref[i]));
		    const double relerr = err/std::max(1e-6, maxval); 
		    
		    if (fabs(relerr) >= tol && fabs(err) >= tol)
			printf("%d: %e ref: %e -> %e %e\n", i, res[i], ref[i], err, relerr);
		    
		    assert(fabs(relerr) < tol || fabs(err) < tol);
		    
		    l1 += fabs(err);
		    l1_rel += fabs(relerr);
		    
		    linf = std::max(linf, fabs(err));
		    linf_rel = std::max(linf_rel, fabs(relerr));
		}
		
		printf("l-infinity errors: %.03e (absolute) %.03e (relative)\n", linf, linf_rel);
		printf("       l-1 errors: %.03e (absolute) %.03e (relative)\n", l1, l1_rel);
	    };
	
	check(xref, xtargets);
	check(yref, ytargets);
    }
	
    delete [] xdst;
    delete [] ydst;
    
    delete [] xtargets;
    delete [] ytargets;
    
    delete [] xref;
    delete [] yref;
	
    delete [] xsrc;
    delete [] ysrc;
    delete [] sources;

    printf("TEST PASSED.\n");
}

int main(int argc, char ** argv)
{
    srand48(1451);

    double theta = 1, tol = 1e-8;
    bool verify = true;
    
    if (argc > 1)
	theta = atof(argv[1]);

    if (argc > 2)
	tol = atof(argv[2]);

    if (argc > 3)
	verify = strcmp(argv[3], "profile") != 0;
    
    auto file2test = [&] (const char * filename)
	{
	    if (access(filename, R_OK) == -1)
	    {
		printf("WARNING: reference file <%s> not found, i skip this test.\n", filename);
		return;
	    }
	    
	    FILE * fin = fopen(filename, "r");
	    assert(fin && sizeof(realtype) == sizeof(double));
	    test_potential(theta, tol, fin, verify);
	    
	    fseek(fin, 0, SEEK_SET);
	    
	    test_force(theta, tol * 100, fin, verify);
	    fclose(fin);
	};

    file2test("diegoBinaryN400");
    file2test("diegoBinaryN2000");
    //file2test("diegoBinaryN12000");

    for(int itest = 0; itest < 10; ++itest)
    {
	printf("test nr %d\n", itest);
	test_potential(theta, tol, NULL, verify);
	test_force(theta, tol * 100, NULL, verify);
    }
}
