## MAKEFILE FOR RascalC. This compiles the grid_covariance.cpp file into the ./cov exececutable.

CC = gcc
CFLAGS = -O3 -Wall -MMD
CXXFLAGS	= -O3 -Wall -MMD -std=c++0x -fopenmp -ffast-math $(shell pkg-config --cflags gsl)
CXXFLAGS	+= -DOPENMP -DLEGENDRE_MIX -DJACKKNIFE -DPRINTPERCENTS
#-DOPENMP  # use this to run multi-threaded with OPENMP
#-DPERIODIC # use this to enable periodic behavior
#-DLEGENDRE # use this to compute 2PCF covariances in Legendre bins (original mode, corresponding to direct accumulation into multipoles from pair counts)
#-DLEGENDRE_MIX # also compute 2PCF covariances in Legendre bins, but in other, "mixed" mode, corresponding to projection of s,µ bins into multipoles
# without either of the two Legendre flags above, the covariance is computed in s,µ bins
#-DJACKKNIFE # use this to compute (r,mu)-space 2PCF covariances and jackknife covariances. Incompatible with -DLEGENDRE but works with -DLEGENDRE_MIX
#-DTHREE_PCF # use this to compute 3PCF autocovariances
#-DPRINTPERCENTS # use this to print percentage of progress in each loop. This can be a lot of output

LFLAGS	= $(shell pkg-config --libs gsl) # common part

# Known OS-specific choices
ifeq ($(shell uname -s),Darwin)
# Here we use LLVM compiler to load the Mac OpenMP. Tested after installation commands:
# brew install llvm
# brew install libomp
# This may need to be modified with a different installation
ifndef HOMEBREW_PREFIX
HOMEBREW_PREFIX = /usr/local
endif
CXX = ${HOMEBREW_PREFIX}/opt/llvm/bin/clang++
LFLAGS	+= -fopenmp -lomp
else
# default (Linux) case
CXX = g++
LFLAGS	+= -lgomp
endif

LD = $(CXX)

AUNTIE	= cov
AOBJS	= grid_covariance.o ./cubature/hcubature.o ./ransampl/ransampl.o
ADEPS   = ${AOBJS:.o=.d}

.PHONY: main clean

main: $(AUNTIE)

$(AUNTIE):	$(AOBJS) Makefile
	$(LD) $(AOBJS) $(LFLAGS) -o $(AUNTIE)

clean:
	rm -f $(AUNTIE) $(AOBJS) ${ADEPS}

$(AOBJS): Makefile
-include ${ADEPS}
