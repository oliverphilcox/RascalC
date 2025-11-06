## MAKEFILE FOR RascalC. This compiles the grid_covariance.cpp file into the ./cov exececutable.

CC = gcc
CFLAGS = -O3 -Wall -MMD
CXX = g++
CXXFLAGS	= -O3 -Wall -MMD -std=c++11 -ffast-math $(shell pkg-config --cflags gsl)
CXXFLAGS	+= -DOPENMP -DTHREE_PCF -DPRINTPERCENTS
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
# Here we load the Mac OpenMP. Tested after installation commands:
# brew install libomp
# This may need to be modified with a different installation
ifndef HOMEBREW_PREFIX
HOMEBREW_PREFIX = /usr/local
endif
CXXFLAGS += -I$(HOMEBREW_PREFIX)/opt/libomp/include
LFLAGS	+= -L$(HOMEBREW_PREFIX)/opt/libomp/lib -lomp
else
# default (Linux) case
CXXFLAGS += -fopenmp
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
