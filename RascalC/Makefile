## MAKEFILE FOR RascalC. This compiles the grid_covariance.cpp file into the ./cov exececutable.

CC = gcc
CFLAGS = -O3 -Wall -MMD
CXX = g++
CXXFLAGS	= -O3 -Wall -MMD -std=c++11 -ffast-math $(shell pkg-config --cflags gsl)
CXXFLAGS	+= -DOPENMP
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

# The code below compiles all the valid variants for the 2PCF covariance

BASE_VARIANTS	= s_mu s_mu_jackknife legendre_accumulated legendre_projected legendre_projected_jackknife

DEFINES_FOR_s_mu=
DEFINES_FOR_s_mu_jackknife=-DJACKKNIFE
DEFINES_FOR_legendre_accumulated=-DLEGENDRE
DEFINES_FOR_legendre_projected=-DLEGENDRE_MIX
DEFINES_FOR_legendre_projected_jackknife=-DLEGENDRE_MIX -DJACKKNIFE

define add_periodic_variant
DEFINES_FOR_$(1)_periodic	= $$(DEFINES_FOR_$(1)) -DPERIODIC
VARIANTS_INT	+= $(1) $(1)_periodic
endef

$(foreach variant,$(BASE_VARIANTS),$(eval $(call add_periodic_variant,$(variant))))

define add_verbose_variant
DEFINES_FOR_$(1)_periodic	= $$(DEFINES_FOR_$(1)) -DPRINTPERCENTS
VARIANTS	+= $(1) $(1)_verbose
endef

$(foreach variant,$(VARIANTS_INT),$(eval $(call add_verbose_variant,$(variant))))

OBJDIR	= obj
EXECDIR	= bin

VAR_SRC_BASE	= grid_covariance
VAR_SRC	= ../$(VAR_SRC_BASE).cpp
VAR_OBJ_BASE	= $(OBJDIR)/$(VAR_SRC_BASE)

EXEC_BASE	= $(EXECDIR)/cov
COMMON_OBJS	= ../cubature/hcubature.o ../ransampl/ransampl.o

ALL_OBJS	= $(COMMON_OBJS)

define make_targets
OBJ_$(1)	= $$(VAR_OBJ_BASE).$(1).o
ALL_OBJS	+= $$(OBJ_$(1))
EXEC_$(1)	= $$(EXEC_BASE).$(1)
ALL_EXECS	+= $$(EXEC_$(1))

$$(OBJ_$(1)): $$(VAR_SRC)
	mkdir -p $$(OBJDIR)
	$$(CXX) $$(CXXFLAGS) $$(DEFINES_FOR_$(1)) -c $$(VAR_SRC) -o $$@

$$(EXEC_$(1)): $$(OBJ_$(1)) $$(COMMON_OBJS)
	mkdir -p $$(EXECDIR)
	$$(LD) $$(OBJ_$(1)) $$(COMMON_OBJS) $$(LFLAGS) -o $$@

endef

$(foreach variant,$(VARIANTS),$(eval $(call make_targets,$(variant))))

.PHONY: all clean

all: $(ALL_EXECS)
.DEFAULT_GOAL	:= all

clean:
	rm -f $(ALL_EXECS) $(ALL_OBJS) ${ALL_DEPS}
	rm -rf $(OBJDIR) $(EXECDIR)

$(ALL_OBJS): Makefile
ALL_DEPS	= ${ALL_OBJS:.o=.d}
-include ${ALL_DEPS}
