// utils.cpp - Useful functions for use in grid_covariance.cpp and subclasses

#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
//#include <random>
#include <complex>
#include <algorithm>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include "../threevector.hh"
#include <gsl/gsl_rng.h>
#include "../ransampl/ransampl.h"
#include "../STimer.cc"
#include "../cubature/cubature.h"
#include <limits>


#undef ALLOUT
//#undef PERIODIC
#define PAGE 4096     // To force some memory alignment.




// We need a vector floor3 function
Float3 floor3(float3 p) {
    return Float3(floor(p.x), floor(p.y), floor(p.z));
}
