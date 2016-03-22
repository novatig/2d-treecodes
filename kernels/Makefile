#
# kernels/Makefile
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
compute ?= 35

SRCS =order$(order)-upward-kernels.h

ifneq "$(MAKECMDGOALS)" "force"
	SRCS += order$(order)-potential-kernels.h
else
	SRCS += order$(order)-force-kernels.h
endif

potential force: $(SRCS)

order$(order)-%.h: %.m4 Makefile
	m4 -D realtype=$(real) -D ORDER=$(order) $< > $@

clean:
	rm -f *.o *kernels*.c *kernels*.h

.PRECIOUS: order$(order)-%.h

.PHONY = clean potential force

