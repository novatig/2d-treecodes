#pragma once

typedef REAL realtype; 

#include <cstdio>
#include <cuda_runtime.h>

#define CUDA_CHECK(ans) do { cudaAssert((ans), __FILE__, __LINE__); } while(0)
inline void cudaAssert(cudaError_t code, const char *file, int line)
{
    if (code != cudaSuccess)
    {
	fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);

	abort();
    }
}

#if !defined(__CUDA_ARCH__)
//#warning __CUDA_ARCH__ not defined! assuming 350
#define ACCESS(x) __ldg(&x)
#elif __CUDA_ARCH__ >= 350
#define ACCESS(x) __ldg(&x)
#else
#define ACCESS(x) (x)
#endif
