#
# drivers/Makefile
# Part of 2d-treecodes
#
# Created and authored by Diego Rossinelli on 2015-09-25.
# Copyright 2015. All rights reserved.
#
# Users are NOT authorized
# to employ the present software for their own publications
# before getting a written permission from the author of this file.
#

real ?= double
order ?= 12
mrag-blocksize ?= 32
config ?=release
#profiling = 1
compute ?= 35

TLPFLAGS = -I../kernels/ -std=c++11 -DREAL=$(real) -DBLOCKSIZE=$(mrag-blocksize) -DORDER=$(order) -DLMAX=15 -O4

NVCCFLAGS := $(TLPFLAGS) -arch=compute_$(compute) -code=sm_$(compute) \
	-Xcompiler -fopenmp \
	-Xcompiler '-fPIC' \
	-Xcompiler '-fvisibility=hidden' \
	-Xcompiler -O4  \
	-Xcompiler -mfpmath=sse

TLPFLAGS += -march=native -fopenmp -fvisibility=hidden -use_fast_math

ifeq "$(config)" "release"
	TLPFLAGS +=  -DNDEBUG
	NVCCFLAGS += -DNDEBUG
endif

ifeq "$(profiling)" "1"
	TLPFLAGS += -pg
	NVCCFLAGS += -lineinfo
endif

potential: potentialkernels treecode-potential.o order$(order)-upward.o sort-sources.o driver.o

force: forcekernels treecode-force.o order$(order)-upward.o sort-sources.o

potentialkernels:
	make -C ../kernels potential

forcekernels:
	make -C ../kernels force

treecode-potential.o: treecode-potential.cu treecode-potential.h Makefile
	nvcc $(NVCCFLAGS) -DHEADER_POTENTIAL_KERNEL="<order$(order)-potential-kernels.h>" -c $<

driver.o: driver.cpp
	g++ $(TLPFLAGS) -fPIC -c $<

treecode-force.o: treecode-force.cu treecode-force.h Makefile
	nvcc $(NVCCFLAGS) -DHEADER_FORCE_KERNEL="<order$(order)-force-kernels.h>" -c $<

order$(order)-upward.o: upward.cu upward.h Makefile
	nvcc $(NVCCFLAGS) -DHEADER_UPWARD_KERNEL="<order$(order)-upward-kernels.h>" -c $< -o $@

sort-sources.o: sort-sources.cu sort-sources.h Makefile
	nvcc $(NVCCFLAGS) -c $< -o $@

clean:
	rm -f *.o 

.PHONY = clean potential force

