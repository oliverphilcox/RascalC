## MAKEFILE FOR RascalC. This compiles the grid_covariance.cpp file into the ./cov exececutable.

CC = gcc
CFLAGS = -g -O3 -Wall -MMD
CXXFLAGS = -O3 -Wall -MMD -DOPENMP -DLEGENDRE
#-DOPENMP  # use this to run multi-threaded with OPENMP
#-DPERIODIC # use this to enable periodic behavior
#-DLEGENDRE # use this to compute 2PCF covariances in Legendre bins
#-DJACKKNIFE # use this to compute (r,mu)-space 2PCF covariances and jackknife covariances
#-DTHREE_PCF # use this to compute 3PCF autocovariances

# Known OS-specific choices
ifeq ($(shell uname -s),Darwin)
# Here we use LLVM compiler to load the Mac OpenMP. Tested after installation commands:
# brew install llvm
# brew install libomp
# This may need to be modified with a different installation
ifndef $(HOMEBREW_PREFIX)
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
