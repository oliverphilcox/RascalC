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

# The code below compiles all the valid variants for the 2PCF covariance, and one for 3PCF

BASE_VARIANTS	= s_mu s_mu_jackknife legendre_accumulated legendre_projected legendre_projected_jackknife 3pcf_legendre_accumulated

DEFINES_FOR_s_mu=
DEFINES_FOR_s_mu_jackknife=-DJACKKNIFE
DEFINES_FOR_legendre_accumulated=-DLEGENDRE
DEFINES_FOR_legendre_projected=-DLEGENDRE_MIX
DEFINES_FOR_legendre_projected_jackknife=-DLEGENDRE_MIX -DJACKKNIFE
DEFINES_FOR_3pcf_legendre_accumulated=-DTHREE_PCF

define add_periodic_variant
DEFINES_FOR_$(1)_periodic	= $$(DEFINES_FOR_$(1)) -DPERIODIC
VARIANTS_INT	+= $(1) $(1)_periodic
endef

$(foreach variant,$(BASE_VARIANTS),$(eval $(call add_periodic_variant,$(variant))))

define add_verbose_variant
DEFINES_FOR_$(1)_verbose	= $$(DEFINES_FOR_$(1)) -DPRINTPERCENTS
VARIANTS	+= $(1) $(1)_verbose
endef

$(foreach variant,$(VARIANTS_INT),$(eval $(call add_verbose_variant,$(variant))))

OBJDIR	= obj
EXECDIR	= bin

$(OBJDIR) $(EXECDIR):
	mkdir -p $@

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

$$(OBJ_$(1)): $$(VAR_SRC) | $$(OBJDIR)
	$$(CXX) $$(CXXFLAGS) $$(DEFINES_FOR_$(1)) -c $$(VAR_SRC) -o $$@

$$(EXEC_$(1)): $$(OBJ_$(1)) $$(COMMON_OBJS) | $$(EXECDIR)
	$$(LD) $$(OBJ_$(1)) $$(COMMON_OBJS) $$(LFLAGS) -o $$@
endef
# directories that need to exist are listed after "|" as order-only prerequisites, as suggested in https://stackoverflow.com/a/38525905/23322509
# if directories are listed as normal prerequisites, they are modified when objects/executables are made/remade and this alone triggers another remake of objects/executables unnecessarily

$(foreach variant,$(VARIANTS),$(eval $(call make_targets,$(variant))))

# end of variants for the covariance

# The code below compiles all the valid variants for the triple_counts

VARIANTS2 = s_mu s_mu_verbose s_mu_periodic s_mu_periodic_verbose

DEFINES2 = -DTRIPLE # triple_counts needs this to have mumin=-1 and not 0 to be consistent with the 3PCF covariance

VAR_SRC2_BASE	= triple_counts
VAR_SRC2	= ../$(VAR_SRC2_BASE)/$(VAR_SRC2_BASE).cpp
VAR_OBJ2_BASE	= $(OBJDIR)/$(VAR_SRC2_BASE)

EXEC2_BASE	= $(EXECDIR)/triple

define make_targets2
OBJ2_$(1)	= $$(VAR_OBJ2_BASE).$(1).o
ALL_OBJS	+= $$(OBJ2_$(1))
EXEC2_$(1)	= $$(EXEC2_BASE).$(1)
ALL_EXECS	+= $$(EXEC2_$(1))

$$(OBJ2_$(1)): $$(VAR_SRC2) | $$(OBJDIR)
	$$(CXX) $$(CXXFLAGS) $$(DEFINES2) $$(DEFINES_FOR_$(1)) -c $$(VAR_SRC2) -o $$@

$$(EXEC2_$(1)): $$(OBJ2_$(1)) $$(COMMON_OBJS) | $$(EXECDIR)
	$$(LD) $$(OBJ2_$(1)) $$(COMMON_OBJS) $$(LFLAGS) -o $$@
endef

$(foreach variant,$(VARIANTS2),$(eval $(call make_targets2,$(variant))))

# end of variants for the triple_counts

.PHONY: all clean

all: $(ALL_EXECS)
.DEFAULT_GOAL	:= all

clean:
	rm -f $(ALL_EXECS) $(ALL_OBJS) ${ALL_DEPS}
	rm -rf $(OBJDIR) $(EXECDIR)

$(ALL_OBJS): Makefile
ALL_DEPS	= ${ALL_OBJS:.o=.d}
-include ${ALL_DEPS}
