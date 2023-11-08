## MAKEFILE FOR RascalC. This compiles the grid_covariance.cpp file into the ./cov exececutable.

CC = gcc
CFLAGS = -g -O3 -Wall -MMD
CXXFLAGS = -O3 -Wall -MMD -DOPENMP
#-DOPENMP  # use this to run multi-threaded with OPENMP
#-DPERIODIC # use this to enable periodic behavior
#-DLEGENDRE # use this to compute 2PCF covariances in Legendre bins (original mode, corresponding to direct accumulation into multipoles from pair counts)
#-DLEGENDRE_MIX # also compute 2PCF covariances in Legendre bins, but in other, "mixed" mode, corresponding to projection of s,µ bins into multipoles
# without either of the two Legendre flags above, the covariance is computed in s,µ bins
#-DJACKKNIFE # use this to compute (r,mu)-space 2PCF covariances and jackknife covariances. Incompatible with -DLEGENDRE but works with -DLEGENDRE_MIX
#-DTHREE_PCF # use this to compute 3PCF autocovariances
#-DPRINTPERCENTS # use this to print percentage of progress in each loop. This can be a lot of output

# Known OS-specific choices
ifeq ($(shell uname -s),Darwin)
# Here we use LLVM compiler to load the Mac OpenMP. Tested after installation commands:
# brew install llvm
# brew install libomp
# This may need to be modified with a different installation
ifndef HOMEBREW_PREFIX
HOMEBREW_PREFIX = /usr/local
endif
CXX = ${HOMEBREW_PREFIX}/opt/llvm/bin/clang++ -std=c++0x -fopenmp -ffast-math $(shell pkg-config --cflags gsl)
LD	= ${HOMEBREW_PREFIX}/opt/llvm/bin/clang++
LFLAGS	= $(shell pkg-config --libs gsl) -fopenmp -lomp
else
# default (Linux) case
CXX = g++ -fopenmp -lgomp -std=c++0x -ffast-math $(shell pkg-config --cflags gsl)
LD	= g++
LFLAGS	= -L/usr/local/lib -L/usr/lib/x86_64-linux-gnu $(shell pkg-config --libs gsl) -lgomp
endif

# The code below compiles all the valid variants for the 2PCF covariance

BASE_VARIANTS	= default default_jackknife legendre_orig legendre_mix legendre_mix_jackknife

DEFINES_FOR_default=
DEFINES_FOR_default_jackknife=-DJACKKNIFE
DEFINES_FOR_legendre_orig=-DLEGENDRE
DEFINES_FOR_legendre_mix=-DLEGENDRE_MIX
DEFINES_FOR_legendre_mix_jackknife=-DLEGENDRE_MIX -DJACKKNIFE

PERIODIC_VARIANTS = $(foreach variant,$(BASE_VARIANTS),$(variant)_periodic)

define add_periodic_defines
DEFINES_FOR_$(1)_periodic	= $$(DEFINES_FOR_$(1)) -DPERIODIC
endef

$(foreach variant,$(BASE_VARIANTS),$(eval $(call add_periodic_defines,$(variant))))

VARIANTS	= $(BASE_VARIANTS) $(PERIODIC_VARIANTS)

OBJDIR	= obj
EXECDIR	= bin

VAR_SRC_BASE	= grid_covariance
VAR_SRC	= $(VAR_SRC_BASE).cpp
VAR_OBJ_BASE	= $(OBJDIR)/$(VAR_SRC_BASE)

EXEC_BASE	= $(EXECDIR)/cov
COMMON_OBJS	= ./cubature/hcubature.o ./ransampl/ransampl.o

VAR_OBJS	= $(foreach variant,$(VARIANTS),$(VAR_OBJ_BASE).$(variant).o)
ALL_EXECS	= $(foreach variant,$(VARIANTS),$(EXEC_BASE).$(variant))

ALL_OBJS	= $(VAR_OBJS) $(COMMON_OBJS)
ALL_DEPS	= ${ALL_OBJS:.o=.d}

.PHONY: main clean

main: $(ALL_EXECS)

define make_targets
OBJ_$(1)	= $$(VAR_OBJ_BASE).$(1).o
EXEC_$(1)	= $$(EXEC_BASE).$(1)

$$(OBJ_$(1)): $$(VAR_SRC)
	mkdir -p $$(OBJDIR)
	$$(CXX) $$(CXXFLAGS) $$(DEFINES_FOR_$(1)) -c $$(VAR_SRC) -o $$@

$$(EXEC_$(1)): $$(OBJ_$(1)) $$(COMMON_OBJS)
	mkdir -p $$(EXECDIR)
	$$(LD) $$(OBJ_$(1)) $$(COMMON_OBJS) $$(LFLAGS) -o $$@

endef

$(foreach variant,$(VARIANTS),$(eval $(call make_targets,$(variant))))

clean:
	rm -f $(ALL_EXECS) $(ALL_OBJS) ${ALL_DEPS}
	rm -rf $(OBJDIR) $(EXECDIR)

$(ALL_OBJS): Makefile
-include ${ALL_DEPS}
