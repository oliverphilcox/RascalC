// grid_covariance.cpp -- Alexander Wiegand, started Oct 14, 2016 based on grid_multipoles.cpp of Daniel Eisenstein.

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
#include "threevector.hh"
#include <gsl/gsl_rng.h>
#include "./ransampl/ransampl.h"
#include "STimer.cc"
#include "./cubature/cubature.h"
#include <limits>
#include <sys/stat.h>


//#define OPENMP

// For multi-threading:
#ifdef OPENMP
//#define _GNU_SOURCE
#include <omp.h>
#include <sched.h>
#endif

// In order to not print the output matrices
#define NOPRINT


//#undef PERIODIC
#define PAGE 4096     // To force some memory alignment.

typedef unsigned long long int uint64;

// Could swap between single and double precision here.
// Only double precision has been tested.
typedef double Float;
typedef double3 Float3;


// Define module files
#include "modules/parameters.h"
#include "modules/cell_utilities.h"
#include "modules/grid.h"
#include "modules/correlation_function.h"
#include "modules/random_draws2.h"
#include "modules/integrals.h"
#include "modules/driver.h"
#include "modules/jackknife_weights.h"

// Very ugly way to get the correlation function into the integrator, but hey, it works
CorrelationFunction * RandomDraws2::corr;


STimer TotalTime;

// ====================  Computing the integrals ==================


#include "modules/compute_integral3.h"



// ================================ main() =============================



int main(int argc, char *argv[]) {

	Parameters par=Parameters(argc,argv);
    
    // Read in jackknife weights and RR pair counts
    JK_weights weights(&par);
    
    Float3 shift;
    Particle *orig_p;
    if (!par.make_random){
        orig_p = read_particles(par.rescale, &par.np, par.fname, par.rstart, par.nmax, &weights);
        assert(par.np>0);
        par.perbox = compute_bounding_box(orig_p, par.np, par.rect_boxsize, par.rmax, shift, par.nside);
    } else {
    // If you want to just make random particles instead:
	assert(par.np>0);
	//assert(par.boxsize==par.rescale);    // Nonsense if not!
	orig_p = make_particles(par.rect_boxsize, par.np);
    // set as periodic if we make the random particles
    par.perbox = true;
    }
    
    if (par.qinvert) invert_weights(orig_p, par.np);
    if (par.qbalance) balance_weights(orig_p, par.np);

    // Now ready to compute!
    // Sort the particles into the grid.
    Grid grid(orig_p, par.np, par.rect_boxsize, par.nside, shift);

    Float grid_density = (double)par.np/grid.nf;//pow(par.nside,3);
    printf("Average number of particles per grid cell = %6.2f\n", grid_density);
    Float max_density = 16.0;
    if (grid_density>max_density){
        fprintf(stderr,"Average particle density exceeds maximum advised particle density (%.0f particles per cell) - exiting.\n",max_density);
        exit(1);
    }
    printf("Average number of particles per max_radius ball = %6.2f\n",
    		par.np*4.0*M_PI/3.0*pow(par.rmax,3.0)/(par.rect_boxsize.x*par.rect_boxsize.y*par.rect_boxsize.z));
    if (grid_density<2){
        printf("#\n# WARNING: grid appears inefficiently fine; exiting.\n#\n");
        exit(1);
    }

    printf("# Done gridding the particles\n");
    printf("# %d particles in use, %d with positive weight\n", grid.np, grid.np_pos);
    printf("# Weights: Positive particles sum to %f\n", grid.sumw_pos);
    printf("#          Negative particles sum to %f\n", grid.sumw_neg);
    free(orig_p); // Particles are now only stored in grid

    fflush(NULL);

    // Everything above here takes negligible time.  This line is nearly all of the work.
    compute_integral3(&grid,&par,&weights);
    
    rusage ru;
    getrusage(RUSAGE_SELF, &ru);

    fprintf(stderr,"# user CPU time used: %ld \n# system CPU time used: %ld \n# maximum resident set size: %ld \n# integral shared memory size: %ld \n# integral unshared data size: %ld \n# integral unshared stack size: %ld \n# page reclaims (soft page faults): %ld \n# page faults (hard page faults): %ld \n# swaps: %ld \n# block input operations: %ld \n# block output operations: %ld \n# IPC messages sent: %ld \n# IPC messages received: %ld \n# signals received: %ld \n# voluntary context switches: %ld \n# involuntary context switches: %ld \n",ru.ru_utime.tv_sec,ru.ru_stime.tv_sec,ru.ru_maxrss,ru.ru_ixrss,ru.ru_idrss,ru.ru_isrss,ru.ru_minflt,ru.ru_majflt,ru.ru_nswap,ru.ru_inblock,ru.ru_oublock,ru.ru_msgsnd,ru.ru_msgrcv,ru.ru_nsignals,ru.ru_nvcsw,ru.ru_nivcsw);

    return 0;
}


