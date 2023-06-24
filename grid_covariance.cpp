// grid_covariance.cpp -- Oliver Philcox / Alexander Wiegand, started Oct 14, 2016 based on grid_multipoles.cpp of Daniel Eisenstein.

#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <complex>
#include <algorithm>
#include <random>
#include <gsl/gsl_sf.h>
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

// For multi-threading:
#ifdef OPENMP
#include <omp.h>
#include <sched.h>
#endif

// In order to not print the output matrices

#define PAGE 4096     // To force some memory alignment.

typedef unsigned long long int uint64;

// Could swap between single and double precision here.
typedef double Float;
typedef double3 Float3;


// Define module files
    #include "modules/parameters.h"
    #include "modules/cell_utilities.h"
    #include "modules/grid.h"
    #include "modules/correlation_function.h"
#ifdef THREE_PCF
    #include "modules/legendre_utilities.h"
    #include "modules/integrals_3pcf.h"
#elif defined LEGENDRE
    #include "modules/legendre_utilities.h"
    #include "modules/integrals_legendre.h"
#elif defined POWER
    #include "power_spectra/power_mod/kernel_interp.h"
    #include "modules/integrals_legendre_power.h"
#else
    #include "modules/jackknife_weights.h"
    #include "modules/integrals.h"
#endif
    #include "modules/random_draws.h"
    #include "modules/driver.h"

// Get the correlation function into the integrator
CorrelationFunction * RandomDraws::corr;

STimer TotalTime;

// ====================  Computing the integrals ==================

#ifdef THREE_PCF
    #include "modules/compute_integral_3pcf.h"
#else
    #include "modules/compute_integral.h"
#endif
#include "modules/rescale_correlation.h"

// ================================ main() =============================


int main(int argc, char *argv[]) {

	Parameters par=Parameters(argc,argv);

    int max_no_functions=1; // required number of xi / random_draws / jackknife_weight functions
    int no_fields=1; // number of different fields used
    if(par.multi_tracers==true){
        max_no_functions=3;
        no_fields=2;
    }
#if (defined LEGENDRE)||(defined THREE_PCF)||(defined POWER)
    // Define all possible survey correction functions
    SurveyCorrection all_survey[max_no_functions]; // create empty functions
#ifdef THREE_PCF
    SurveyCorrection tmp_sc(&par);
    all_survey[0].copy(&tmp_sc); // save into global memory
#else
    SurveyCorrection tmp_sc(&par,1,1);
    all_survey[0].copy(&tmp_sc); // save into global memory

    if(par.multi_tracers==true){
        SurveyCorrection tmp_sc12(&par,1,2), tmp_sc2(&par,2,2);
        all_survey[1].copy(&tmp_sc2);
        all_survey[2].copy(&tmp_sc12);
    }
#endif
#else
    // Read in jackknife weights and RR pair counts (if JACKKNIFE is not defined just includes RR counts)
    JK_weights all_weights[max_no_functions]; // create empty functions

    JK_weights tmp(&par,1,1);
    all_weights[0].copy(&tmp); // save into global memory

    if(par.multi_tracers==true){
        JK_weights tmp12(&par,1,2), tmp2(&par,2,2);
        all_weights[1].copy(&tmp2);
        all_weights[2].copy(&tmp12);
    }
#endif

    // Now read in particles
    Particle* all_particles[no_fields];
    int all_np[no_fields];

    for (int index = 0; index < no_fields; index++) {
        Float3 shift;
        if (!par.make_random) {
            char *filename;
            if (index == 0) filename = par.fname;
            else filename = par.fname2;
#ifdef JACKKNIFE
            all_particles[index] = read_particles(par.rescale, &all_np[index], filename, par.rstart, par.nmax, &all_weights[index]);
#else
            all_particles[index] = read_particles(par.rescale, &all_np[index], filename, par.rstart, par.nmax);
#endif
            assert(all_np[index] > 0);
        }
        else {
            // If you want to just make random particles instead:
            assert(par.np > 0);
            all_particles[index] = make_particles(par.rect_boxsize, all_np[index], index);
            all_np[index] = par.np;
        }
        if (par.qinvert) invert_weights(all_particles[index], all_np[index]);
        if (par.qbalance) balance_weights(all_particles[index], all_np[index]);
    }

    // Now put particles to grid(s)
    Grid all_grid[no_fields]; // create empty grids
    Float max_density = 16., min_density = 2.;
    int nside_attempts = 3; // number of attempts to meet the constraints by changing nside
    bool nside_global_failure = true; // assume failure until success

    for (int no_attempt = 0; no_attempt <= nside_attempts; no_attempt++) {
        // Compute bounding box using all particles. Do it inside the attempt loop because nside could change.
        Float3 shift; // default value is zero
        if (!par.make_random) {
            par.perbox = compute_bounding_box(all_particles, all_np, no_fields, par.rect_boxsize, par.cellsize, par.rmax, shift, par.nside);
#ifdef PERIODIC
            par.rect_boxsize = {par.boxsize, par.boxsize, par.boxsize}; // restore the given boxsize if periodic
            par.cellsize = par.boxsize / (Float)par.nside; // set cell size manually
            // keep the shift from compute_bounding_box, allowing for coordinate ranges other than [0, par.boxsize) but still of length par.boxsize - this is quite generic and precise at the same time.
#endif
        }
        else {
            // If randoms particles were made we keep the boxsize
            par.cellsize = par.boxsize / (Float)par.nside;
            // set as periodic if we make the random particles
            par.perbox = true;
        }

        // Create grid(s) and see if the particle density is acceptable
        bool nside_local_success = true; // assume this attempt succeeded be default, can be unset
        for (int index = 0; index < no_fields; index++) {
            // Now ready to compute!
            // Sort particles into grid(s)
            Float nofznorm = par.nofznorm;
            if (index == 1) nofznorm = par.nofznorm2;
            Grid tmp_grid(all_particles[index], all_np[index], par.rect_boxsize, par.cellsize, par.nside, shift, nofznorm);

            Float grid_density = (Float)tmp_grid.np/tmp_grid.nf;
            printf("\n RANDOM CATALOG %d DIAGNOSTICS:\n", index+1);
            printf("Average number of particles per grid cell = %6.2f\n", grid_density);
            if (grid_density > max_density) {
                nside_local_success = false; // unset this attempt's success flag
                Float aimed_density = cbrt(max_density * max_density * min_density); // aim for density between the limits but closer to max
                Float nside_approx = cbrt(grid_density/aimed_density) * par.nside; // approximate value of nside to reach this density
                par.nside = 2 * (int)round((nside_approx + 1)/2) - 1; // round to closest odd integer
                fprintf(stderr,"# WARNING: Average particle density exceeds maximum advised particle density (%.0f particles per cell). Setting nside=%d.\n", max_density, par.nside);
                break; // terminate the inner, tracer loop
            }
            if (grid_density < min_density) {
                nside_local_success = false; // unset this attempt's success flag
                Float aimed_density = cbrt(max_density * min_density * min_density); // aim for density between the limits but closer to min
                Float nside_approx = cbrt(grid_density/aimed_density) * par.nside; // approximate value of nside to reach this density
                par.nside = 2 * (int)round((nside_approx + 1)/2) - 1; // round to closest odd integer
                fprintf(stderr, "# WARNING: grid appears inefficiently fine (average density less than %.0f particles per cell). Setting nside=%d.\n", min_density, par.nside);
                break; // terminate the inner, tracer loop
            }
            printf("Average number of particles per max_radius ball = %6.2f\n",
                    tmp_grid.np*4.0*M_PI/3.0*pow(par.rmax,3.0)/(par.rect_boxsize.x*par.rect_boxsize.y*par.rect_boxsize.z));

            printf("# Done gridding the particles\n");
            printf("# %d particles in use, %d with positive weight\n", tmp_grid.np, tmp_grid.np_pos);
            printf("# Weights: Positive particles sum to %f\n", tmp_grid.sumw_pos);
            printf("#          Negative particles sum to %f\n", tmp_grid.sumw_neg);

            // Now save grid to global memory:
            all_grid[index].copy(&tmp_grid);

            free(all_particles[index]); // Particles are now only stored in grid

            fflush(NULL);
        }
        if (nside_local_success) {
            nside_global_failure = false; // unset global failure
            break; // terminate attempt loop
        }
    }
    if (nside_global_failure) { // report and terminate
        fprintf(stderr, "# ERROR: could not meet mean grid density constraints after %d additional attempts.\n", nside_attempts);
        exit(1);
    }

    // Print the resulting grid size to be sure the stderr messages are not missed
    printf("Final grid = %d\n", par.nside);
    // Print box size and max radius in grid units here, because they are adjusted while reading particles (non-periodic case)
    printf("Box Size = {%6.5e,%6.5e,%6.5e}\n", par.rect_boxsize.x, par.rect_boxsize.y, par.rect_boxsize.z);
    Float gridsize = par.rmax / par.cellsize;
    printf("Max Radius in Grid Units = %6.5e\n", gridsize);
    if (gridsize < 1) printf("#\n# WARNING: grid appears inefficiently coarse\n#\n");

#if (!defined LEGENDRE && !defined THREE_PCF && !defined POWER)
    // Now rescale weights based on number of particles (whether or not using jackknives)
    all_weights[0].rescale(all_grid[0].norm,all_grid[0].norm);
    if(par.multi_tracers==true){
        all_weights[1].rescale(all_grid[1].norm,all_grid[1].norm);
        all_weights[2].rescale(all_grid[0].norm,all_grid[1].norm);
    }
#else
    // Now rescale correction functions based on number of particles - similar to RR counts above
    all_survey[0].rescale(all_grid[0].norm,all_grid[0].norm);
    if(par.multi_tracers==true){
        all_survey[1].rescale(all_grid[1].norm,all_grid[1].norm);
        all_survey[2].rescale(all_grid[0].norm,all_grid[1].norm);
    }
#endif

    // Now define all possible correlation functions and random draws:
    CorrelationFunction all_cf[max_no_functions];
    RandomDraws all_rd[max_no_functions];

    CorrelationFunction tmp_cf(par.corname, par.nbin_cf, par.radial_bins_low_cf, par.radial_bins_high_cf, par.mbin_cf, par.mumax-par.mumin);
    all_cf[0].copy_function(&tmp_cf);
    RandomDraws tmp_rd(&tmp_cf,&par,NULL,0);
    all_rd[0].copy(&tmp_rd);

    if(par.multi_tracers==true){
         CorrelationFunction tmp_cf12(par.corname12, par.nbin_cf, par.radial_bins_low_cf, par.radial_bins_high_cf, par.mbin_cf, par.mumax-par.mumin), tmp_cf2(par.corname2, par.nbin_cf, par.radial_bins_low_cf, par.radial_bins_high_cf, par.mbin_cf, par.mumax-par.mumin);
        all_cf[1].copy_function(&tmp_cf2);
        all_cf[2].copy_function(&tmp_cf12);
        RandomDraws rd2(&tmp_cf2,&par,NULL,0), rd12(&tmp_cf12,&par,NULL,0);
        all_rd[1].copy(&rd2);
        all_rd[2].copy(&rd12);
    }

    // Rescale correlation functions
    rescale_correlation rescale(&par);
    rescale.refine_wrapper(&par, all_grid, all_cf, all_rd, max_no_functions);

#ifdef THREE_PCF
    // Compute threePCF integrals
    compute_integral(&all_grid[0],&par,&all_cf[0],&all_rd[0],&all_survey[0],1); // final digit is iteration number
    compute_integral(&all_grid[0],&par,&all_cf[0],&all_rd[0],&all_survey[0],2);
#elif (defined LEGENDRE || defined POWER)
    // Compute integrals
    compute_integral(all_grid,&par,all_cf,all_rd,all_survey,1,1,1,1,1); // final digit is iteration number

    if(par.multi_tracers==true){
        compute_integral(all_grid,&par,all_cf,all_rd,all_survey,1,2,1,1,2);
        compute_integral(all_grid,&par,all_cf,all_rd,all_survey,1,2,2,1,3);
        compute_integral(all_grid,&par,all_cf,all_rd,all_survey,1,2,1,2,4);
        compute_integral(all_grid,&par,all_cf,all_rd,all_survey,1,1,2,2,5);
        compute_integral(all_grid,&par,all_cf,all_rd,all_survey,2,1,2,2,6);
        compute_integral(all_grid,&par,all_cf,all_rd,all_survey,2,2,2,2,7);
    }
#else
    // Compute integrals
    compute_integral(all_grid,&par,all_weights,all_cf,all_rd,1,1,1,1,1); // final digit is iteration number

    if(par.multi_tracers==true){
        compute_integral(all_grid,&par,all_weights,all_cf,all_rd,1,2,1,1,2);
        compute_integral(all_grid,&par,all_weights,all_cf,all_rd,1,2,2,1,3);
        compute_integral(all_grid,&par,all_weights,all_cf,all_rd,1,2,1,2,4);
        compute_integral(all_grid,&par,all_weights,all_cf,all_rd,1,1,2,2,5);
        compute_integral(all_grid,&par,all_weights,all_cf,all_rd,2,1,2,2,6);
        compute_integral(all_grid,&par,all_weights,all_cf,all_rd,2,2,2,2,7);
    }
#endif

    rusage ru;
    getrusage(RUSAGE_SELF, &ru);

    fprintf(stderr,"# user CPU time used: %ld \n# system CPU time used: %ld \n# maximum resident set size: %ld \n# integral shared memory size: %ld \n# integral unshared data size: %ld \n# integral unshared stack size: %ld \n# page reclaims (soft page faults): %ld \n# page faults (hard page faults): %ld \n# swaps: %ld \n# block input operations: %ld \n# block output operations: %ld \n# IPC messages sent: %ld \n# IPC messages received: %ld \n# signals received: %ld \n# voluntary context switches: %ld \n# involuntary context switches: %ld \n",ru.ru_utime.tv_sec,ru.ru_stime.tv_sec,ru.ru_maxrss,ru.ru_ixrss,ru.ru_idrss,ru.ru_isrss,ru.ru_minflt,ru.ru_majflt,ru.ru_nswap,ru.ru_inblock,ru.ru_oublock,ru.ru_msgsnd,ru.ru_msgrcv,ru.ru_nsignals,ru.ru_nvcsw,ru.ru_nivcsw);

    return 0;
}
