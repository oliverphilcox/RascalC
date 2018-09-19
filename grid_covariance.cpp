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

//#define OPENMP

// For multi-threading:
#ifdef OPENMP
//#define _GNU_SOURCE
#include <omp.h>
#include <sched.h>
#endif


//#undef PERIODIC
#define PAGE 4096     // To force some memory alignment.

typedef unsigned long long int uint64;

// Could swap between single and double precision here.
// Only double precision has been tested.
typedef double Float;
typedef double3 Float3;


#include "modules/parameters.h"


// We need a vector floor3 function
Float3 floor3(float3 p) {
    return Float3(floor(p.x), floor(p.y), floor(p.z));
}

// =================== Particles ====================
// This is the info about each particle that we load in and store in the Grid.

class Particle {
  public:
    Float3 pos;
    Float w;  // The weight for each particle
};


// ====================  The Cell and Grid classes ==================

/* The Grid class holds a new copy of the particles.
These are sorted into cells, and all positions are referenced to the
cell center.  That way, we can handle periodic wrapping transparently,
simply by the cell indexing.

For simplicity, we opt to flatten the index of the cells into a 1-d number.
For example, this makes multi-threading over cells simpler.
*/

class Cell {
  public:
    int start;	// The starting index of the particle list
    int np;
};


// load modules
#include "modules/grid.h"
#include "modules/correlation_function.h"
#include "modules/random_draws.h"

// Very ugly way to get the correlation function into the integrator, but hey, it works

CorrelationFunction * RandomDraws::corr;

// ========================== Accumulate Integral ================

#include "modules/integrals.h"

/*--------------------------------------------------------------------*/

void kofnip(int* cellsel, int n, int k, gsl_rng* rng) {

	/* -------------------------------------------------------- */
	/* selects k points of an ensemble of n points              */
	/* -------------------------------------------------------- */

	int tmp, r;

	assert(k<=n);

	for (int i=0;i<k;i++) {
		r = i + (int) floor(gsl_rng_uniform(rng) * (n - i));
		tmp = cellsel[i];
		cellsel[i] = cellsel[r];
		cellsel[r]=tmp;
	}

	return;

};


/*--------------------------------------------------------------------*/

void kofn(int* cellsel, int n, int k,Integrals * integs, gsl_rng* rng) {

	/* -------------------------------------------------------- */
	/* selects k points of an ensemble of n points              */
	/* -------------------------------------------------------- */

	int* mapping;
	int max_rand, ind, r;
	int i;

	assert(k<=n);

	mapping = (int*)malloc(sizeof(int)*n);

	for (i = 0; i < n; i++) {
		mapping[i] = i;
	}

	max_rand = n - 1;


	for (i=0;i<k;i++) {
		r = (int) floor(gsl_rng_uniform(rng) * (max_rand + 1));
		ind = mapping[r];

		cellsel[i] = ind;

		mapping[r] = mapping[max_rand];
		max_rand = max_rand - 1;
	}

	free(mapping);

	return;

};

integer3 relid(int cin, int maxsep){

	integer3 res;
	res.z = cin%(2*maxsep+1)-maxsep;
	cin = cin/(2*maxsep+1);
	res.y = cin%(2*maxsep+1)-maxsep;
	res.x = cin/(2*maxsep+1)-maxsep;


	return res;
}

STimer TotalTime;

// ====================  Computing the integrals ==================


// Available sampling strategies:
#define ALL 0   // Use all cells in the grid up to a certain radius
#define FLAT 1	// Select cells equiprobably up to a certain radius
#define SEL 2	// Select cells weighted by the distance from the starting cell up to the radius defined by the input parameter xicutoff

#include "modules/compute_integral.h"

// ====================  The Driver ===========================

#include "modules/driver.h"

// ================================ main() =============================



int main(int argc, char *argv[]) {

	Parameters par=Parameters(argc,argv);

	Float3 shift;
    Particle *orig_p;
    if (par.make_random) {
	// If you want to just make random particles instead:
	assert(par.np>0);
	assert(par.boxsize==par.rescale);    // Nonsense if not!
	orig_p = make_particles(par.boxsize, par.np);
    } else {
	orig_p = read_particles(par.rescale, &par.np, par.fname, par.rstart, par.nmax);
	assert(par.np>0);
    }
    if (par.qinvert) invert_weights(orig_p, par.np);
    if (par.qbalance) balance_weights(orig_p, par.np);

    par.perbox=check_bounding_box(orig_p, par.np, par.boxsize, par.rmax,shift);

    // Now ready to compute!
    // Sort the particles into the grid.
    Grid grid(orig_p, par.np, par.boxsize, par.nside, shift);

    Float grid_density = (double)par.np/grid.nf;//pow(par.nside,3);
    printf("Average number of particles per grid cell = %6.2f\n", grid_density);
    printf("Average number of particles per max_radius ball = %6.2f\n",
    		par.np*4.0*M_PI/3.0*pow(par.rmax/par.boxsize,3.0));
    if (grid_density<1) printf("#\n# WARNING: grid appears inefficiently fine.\n#\n");


//   Count occupation of pos and neg cells
//    printf("BEG\n");
//    for(int ii=0;ii<grid.ncells;ii++){
//    	Cell fou=grid.c[ii];
//    	int sne=0;
//    	int spo=0;
//    	for (int l = fou.start; l<fou.start+fou.np; l++)
//    		grid.p[l].w>0?spo++:sne++;
//    	printf("%d %d\n",spo,sne);
//    }
//
//    return 0;

    printf("# Done gridding the particles\n");
    printf("# %d particles in use, %d with positive weight\n", grid.np, grid.np_pos);
    printf("# Weights: Positive particles sum to %f\n", grid.sumw_pos);
    printf("#          Negative particles sum to %f\n", grid.sumw_neg);
    free(orig_p); // Particles are now only stored in grid

    fflush(NULL);

    // Everything above here takes negligible time.  This line is nearly all of the work.
    compute_integral(&grid,&par);

    rusage ru;
    getrusage(RUSAGE_SELF, &ru);

    fprintf(stderr,"# user CPU time used: %ld \n# system CPU time used: %ld \n# maximum resident set size: %ld \n# integral shared memory size: %ld \n# integral unshared data size: %ld \n# integral unshared stack size: %ld \n# page reclaims (soft page faults): %ld \n# page faults (hard page faults): %ld \n# swaps: %ld \n# block input operations: %ld \n# block output operations: %ld \n# IPC messages sent: %ld \n# IPC messages received: %ld \n# signals received: %ld \n# voluntary context switches: %ld \n# involuntary context switches: %ld \n",ru.ru_utime.tv_sec,ru.ru_stime.tv_sec,ru.ru_maxrss,ru.ru_ixrss,ru.ru_idrss,ru.ru_isrss,ru.ru_minflt,ru.ru_majflt,ru.ru_nswap,ru.ru_inblock,ru.ru_oublock,ru.ru_msgsnd,ru.ru_msgrcv,ru.ru_nsignals,ru.ru_nvcsw,ru.ru_nivcsw);

    return 0;
}


