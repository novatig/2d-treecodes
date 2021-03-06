#
# Makefile
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
cpu ?= haswell
isa ?= avx2
masking ?= 64

OBJS = order$(order)-upward-kernels.o upward-common.o

ifneq "$(MAKECMDGOALS)" "force"
	OBJS += order$(order)-potential-kernels.o
else
	OBJS += order$(order)-force-kernels.o order$(order)-force-kernels-f32.o
endif

potential force: $(OBJS)

order$(order)-force-kernels-f32.o: order$(order)-force-kernels-f32.ispc Makefile
	ispc  $< -Drealtype=$(real) -O3 --arch=x86-64 --math-lib=fast \
	--opt=disable-assertions --opt=fast-math \
	--cpu=$(cpu) --target=$(isa)-i32x8 \
	--wno-perf --pic -o $@

%.o: %.ispc Makefile
	ispc  $< -Drealtype=$(real) -O3 --arch=x86-64 --math-lib=fast \
	--opt=disable-assertions --opt=fast-math \
	--cpu=$(cpu) --target=$(isa)-i$(masking)x4 \
	--wno-perf --pic -o $@

order$(order)-%.ispc: %.m4 Makefile
	m4 -D realtype=$(real) -D ORDER=$(order) $< > $@

clean:
	rm -f *.o *kernels*.c *kernels*.ispc

.PRECIOUS: order$(order)-%.ispc order$(order)-%.o

.PHONY = clean potential force
