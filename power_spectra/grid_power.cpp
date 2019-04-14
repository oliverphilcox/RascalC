// grid_power.cpp -- Oliver Philcox - based on Daniel Eisenstein's grid_multipoles.cpp code

#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "../threevector.hh"
#include "../STimer.cc"
#include <sys/stat.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

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
#include "power_mod/power_parameters.h"
#include "../modules/cell_utilities.h"
#include "../modules/grid.h"
#include "../modules/driver.h"
#include "power_mod/kernel_interp.h"
#include "power_mod/survey_correction_legendre.h"
#include "power_mod/pair_counter.h"
#include "power_mod/power_counts.h"

STimer TotalTime;

// =========================== main() ==================================


int main(int argc, char *argv[]) {
    
	Parameters par=Parameters(argc,argv);
        
    // Now read in particles to grid:
    Grid all_grid[2]; // create empty grids
    
    for(int index=0;index<2;index++){
        Float3 shift;
        Particle *orig_p;
        if (!par.make_random){
            char *filename;
            if(index==0) filename=par.fname;
            else filename=par.fname2;
            orig_p = read_particles(par.rescale, &par.np, filename, par.rstart, par.nmax);
            assert(par.np>0);
            par.perbox = compute_bounding_box(orig_p, par.np, par.rect_boxsize, par.rmax, shift, par.nside);
        } else {
        // If you want to just make random particles instead:
        assert(par.np>0);
        orig_p = make_particles(par.rect_boxsize, par.np);
        // set as periodic if we make the random particles
        par.perbox = true;
        }
        
        if (par.qinvert) invert_weights(orig_p, par.np);
        if (par.qbalance) balance_weights(orig_p, par.np);

        // Now ready to compute!
        // Sort the particles into the grid.
        Grid tmp_grid(orig_p, par.np, par.rect_boxsize, par.nside, shift, 1.);

        Float grid_density = (double)par.np/tmp_grid.nf;
        printf("\n RANDOM CATALOG %d DIAGNOSTICS:\n",index+1);
        printf("Grid cell-size = %.2fMpc/h\n", tmp_grid.cellsize);
        printf("Average number of particles per grid cell = %6.2f\n", grid_density);
        Float max_density = 200.0;
        if (grid_density>max_density){
            fprintf(stderr,"Average particle density exceeds maximum advised particle density (%.0f particles per cell) - exiting.\n",max_density);
            exit(1);
        }
        if (grid_density<0.01){
            printf("#\n# WARNING: grid appears inefficiently fine; exiting.\n#\n");
            exit(1);
        }

        printf("# Done gridding the particles\n");
        printf("# %d particles in use, %d with positive weight\n", tmp_grid.np, tmp_grid.np_pos);
        printf("# Weights: Positive particles sum to %f\n", tmp_grid.sumw_pos);
        printf("#          Negative particles sum to %f\n", tmp_grid.sumw_neg);

        // Now save grid to global memory:
        all_grid[index].copy(&tmp_grid);
        
        free(orig_p); // Particles are now only stored in grid
        
        fflush(NULL);
    }
    
    // Read in survey correction fucntion
    SurveyCorrection sc(&par);

    
    // Count number of second field cells enclosed by the maximum truncation radius
    Float cellsize = all_grid[1].cellsize;
    Float filled_vol = 4./3.*M_PI*pow(par.R0+2.*cellsize,3);
    int n_close = ceil(filled_vol/pow(cellsize,3)); // number of close cells
    
    // Define cell separations (dimensionless) within truncation radius
    Float3 cell_sep_close_tmp[n_close];
    int len = ceil((par.R0+cellsize/2.)/cellsize);
    int len_cell_sep_close=0.; // counter
    integer3 this_pos;
    
    for(int i=-1*len;i<=len;i++){
        for(int j=-1*len;j<=len;j++){
            for(int k=-1*len;k<=len;k++){
                if((sqrt(pow(i,2)+pow(j,2)+pow(k,2))*cellsize)<(par.R0+cellsize)){
                    this_pos = {i,j,k};
                    cell_sep_close_tmp[len_cell_sep_close] = this_pos;
                    len_cell_sep_close++;
                }
            }
        }
    }
    assert(len_cell_sep_close<=n_close);  
    
    integer3 cell_sep_close[len_cell_sep_close]; // proper array to house cell separations of correct length
    for(int i=0;i<len_cell_sep_close;i++) cell_sep_close[i] = cell_sep_close_tmp[i];
    
    // Compute kernel interpolation functions
    KernelInterp all_interp[par.nbin*par.mbin];
    for(int i=0;i<par.nbin;i++){
        printf("Creating kernel interpolator for k-bin %d of %d\n",i+1,par.nbin);
        for(int j=0;j<par.mbin;j++){
            KernelInterp tmp_interp(j*2,par.radial_bins_low[i],par.radial_bins_high[i],par.R0);
            all_interp[i*par.mbin+j].copy_function(tmp_interp);
            if((i==0)&&(j==0)) printf("%.2e\n",all_interp[i*par.mbin+j].kernel_vals[2000]);
        }
    }
        
    // RUN Pair Counter
    
    pair_counter(&all_grid[0],&all_grid[1],&par,&sc,all_interp,cell_sep_close,len_cell_sep_close);
    
    rusage ru;
    getrusage(RUSAGE_SELF, &ru);

    fprintf(stderr,"# user CPU time used: %ld \n# system CPU time used: %ld \n# maximum resident set size: %ld \n# integral shared memory size: %ld \n# integral unshared data size: %ld \n# integral unshared stack size: %ld \n# page reclaims (soft page faults): %ld \n# page faults (hard page faults): %ld \n# swaps: %ld \n# block input operations: %ld \n# block output operations: %ld \n# IPC messages sent: %ld \n# IPC messages received: %ld \n# signals received: %ld \n# voluntary context switches: %ld \n# involuntary context switches: %ld \n",ru.ru_utime.tv_sec,ru.ru_stime.tv_sec,ru.ru_maxrss,ru.ru_ixrss,ru.ru_idrss,ru.ru_isrss,ru.ru_minflt,ru.ru_majflt,ru.ru_nswap,ru.ru_inblock,ru.ru_oublock,ru.ru_msgsnd,ru.ru_msgrcv,ru.ru_nsignals,ru.ru_nvcsw,ru.ru_nivcsw);

    return 0;
}

