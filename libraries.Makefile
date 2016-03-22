#
# libraries.Makefile
# Part of MRAG/2d-treecode-potential
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
compute ?= 35

OBJS = drivers/order$(order)-upward.o drivers/sort-sources.o
CUDARTPATH=$(CRAY_CUDATOOLKIT_POST_LINK_OPTS)
NVCCFLAGS = -code=sm_$(compute) -arch=compute_$(compute) -Xcompiler '-fPIC'

ifeq "$(MAKECMDGOALS)" "libtreecode-force.so"
	OBJS += drivers/treecode-force.o
	TARGET=force
else
	OBJS += drivers/treecode-potential.o drivers/driver.o
	TARGET=potential
endif

libtreecode-potential.so: drivers/treecode-potential.h alldrivers
	m4 -D realtype=$(real) drivers/treecode-potential.h | sed '/typedef/d' | sed '/attribute/d' > treecode-potential.h
	g++ -shared -o $@ $(OBJS)  -L/usr/local/cuda/lib64 $(CUDARTPATH) -lcudart

libtreecode-force.so: drivers/treecode-force.h alldrivers
	m4 -D realtype=$(real) drivers/treecode-force.h | sed '/typedef/d' | sed '/attribute/d' > treecode-force.h
	g++ -shared -o $@ $(OBJS) -L/usr/local/cuda/lib64 $(CUDARTPATH) -lcudart

alldrivers:
	make -C drivers $(TARGET)

clean:
	rm -f test *.so treecode-potential.h treecode-force.h
	make -C drivers clean
	make -C kernels clean

.PHONY = clean drivers kernels
