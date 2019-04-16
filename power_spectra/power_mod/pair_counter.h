// pair_counter.h - Module to run the pair counting for grid_power.cpp - Oliver Philcox 2019

#ifndef PAIR_COUNTER_H
#define PAIR_COUNTER_H

#include "power_counts.h"

class pair_counter{
    
private:
    int nbin,mbin;
    int one_grid; // boolean to check if grids are the same
    uint64 cnt2; // number of pairs counted
    
public:
    void check_threads(Parameters *par,int print){
        // Set up OPENMP and define which threads to use
#ifdef OPENMP
        cpu_set_t mask[par->nthread+1];
        int tnum=0;
        sched_getaffinity(0, sizeof(cpu_set_t), &mask[par->nthread]);
        if(print==1) fprintf(stderr, " CPUs used are: ");
        for(int ii=0;ii<64;ii++){
            if(CPU_ISSET(ii, &mask[par->nthread])){
                if(print==1) fprintf(stderr,"%d ", ii);
                CPU_ZERO(&mask[tnum]);
                CPU_SET(ii,&mask[tnum]);
                tnum++;
            }
        }
        fprintf(stderr,"\n");
#endif
    }

public:
    pair_counter(Grid *grid1, Grid *grid2, Parameters *par, SurveyCorrection *sc, KernelInterp *kernel_interp, integer3 *cell_sep, int len_cell_sep){
        
        // Define parameters
        nbin = par->nbin; // number of radial bins
        mbin = par->mbin; // number of Legendre bins
        Float percent_counter=0.,used_cells=0,cell_attempts=0;
        uint64 used_particles=0;
        
        one_grid=0;
        if(!strcmp(par->fname,par->fname2)) one_grid=1;
        
        STimer initial, TotalTime; // time initialization
        initial.Start();
        
        //-----------INITIALIZE OPENMP + CLASSES----------
        PowerCounts global_counts(par,sc,kernel_interp);

        check_threads(par,1); // define threads
        
        fprintf(stderr, "Init time: %g s\n",initial.Elapsed());
        printf("# 1st grid filled cells: %d\n",grid1->nf);
        printf("# All 1st grid points in use: %d\n",grid1->np);
        printf("# Max points in one cell in grid 1%d\n",grid1->maxnp);
        fflush(NULL);
        
        TotalTime.Start(); // Start timer
            
        
#ifdef OPENMP
        #pragma omp parallel firstprivate(par,grid1,grid2,kernel_interp) shared(global_counts,TotalTime) reduction(+:percent_counter,used_cells,used_particles)
        { // start parallel loop
        // Decide which thread we are in
        int thread = omp_get_thread_num();
        assert(omp_get_num_threads()<=par->nthread);
        if (thread==0) printf("# Starting power-spectrum pair-count computation on %d threads.\n", omp_get_num_threads());        
#else
        int thread = 0;
        printf("# Starting power-spectrum pair-count computation single threaded.\n");
        { // start loop
#endif
            // DEFINE LOCAL THREAD VARIABLES
            Particle *prim_list; // list of particles in first cell
            int* prim_ids; // list of particle IDs in first cell
            int prim_id_1D,sec_id_1D; // 1D cell ids
            integer3 prim_id,delta,sec_id; // particle positions in grid
            int mnp = grid1->maxnp; // max number of particles in grid1 cell
            Cell prim_cell,sec_cell; // cell objects
            
            PowerCounts loc_counts(par,sc,kernel_interp);
            
            // Assign memory for intermediate steps
            int ec=0;
            ec+=posix_memalign((void **) &prim_list, PAGE, sizeof(Particle)*mnp);
            ec+=posix_memalign((void **) &prim_ids, PAGE, sizeof(int)*mnp);
            assert(ec==0);
            
#ifdef OPENMP
#pragma omp for schedule(dynamic)
#endif
            // Loop over all filled n1 cells
            for(int n1=0;n1<grid1->nf;n1++){
                
                // Print time left
                if((float(n1)/float(grid1->nf)*100)>=percent_counter){
                    printf("Counting cell %d of %d on thread %d: %.0f percent complete\n",n1+1,grid1->nf,thread,percent_counter);
                        percent_counter+=5.;
                }
                
                // Pick first cell
                prim_id_1D = grid1-> filled[n1]; // 1d ID for cell i 
                prim_id = grid1->cell_id_from_1d(prim_id_1D); // define first cell
                prim_cell = grid1->c[prim_id_1D];
                if(prim_cell.np==0) continue; // skip if empty
                
                cell_attempts+=len_cell_sep; // total number of cells used
                
                // Iterate over all nearby second cells
                for(int n2=0;n2<len_cell_sep;n2++){ // second cell index
                    delta = cell_sep[n2]; // difference in cell positions
                    sec_id = prim_id + delta;
                    sec_id_1D = grid2->test_cell(sec_id);
                    if(sec_id_1D<0) continue; // if cell not in grid
                    sec_cell = grid2->c[sec_id_1D];
                    if(sec_cell.np==0) continue; // if empty cell
                    
                    used_cells++; // number of cells used
                    
                    // Now iterate over particles
                    for(int i=prim_cell.start;i<(prim_cell.start+prim_cell.np);i++){
                        for(int j=sec_cell.start;j<(sec_cell.start+sec_cell.np);j++){
                            used_particles++;
                            if((one_grid==1)&&(i==j)) continue; // skip if identical particles in same grids
                            loc_counts.count_pairs(grid1->p[i],grid2->p[j]);
                        }
                    }
                }
            }
#ifdef OPENMP
#pragma omp critical // only one processor at once
#endif
        {
            // Sum up power-sums
            global_counts.sum_counts(&loc_counts);
            loc_counts.reset();
        }
            
    } // end OPENMP loop
    
    // ----- REPORT AND SAVE OUTPUT ------------
    TotalTime.Stop();
    
    int runtime = TotalTime.Elapsed();
    fflush(NULL);
    printf("\nPAIR COUNTS COMPLETE\n\n"); 
    printf("\nTotal process time for %.2e particle pairs: %d s, i.e. %2.2d:%2.2d:%2.2d hms\n", double(global_counts.used_pairs),runtime, runtime/3600,runtime/60%60,runtime%60);
    printf("We tried %.2e pairs of cells and accepted %.2e pairs of cells.\n", double(cell_attempts),double(used_cells));
    printf("Cell acceptance ratio is %.3f.\n",(double)used_cells/cell_attempts);
    printf("We tried %.2e pairs of particles and accepted %.2e pairs of particles.\n", double(used_particles),double(global_counts.used_pairs));
    printf("Particle acceptance ratio is %.3f.\n",(double)global_counts.used_pairs/used_particles);
    printf("Average of %.2f pairs accepted per primary particle.\n\n",(Float)global_counts.used_pairs/grid1->np);
    
    printf("\nTrial speed: %.2e cell pairs per core per second\n",double(used_cells)/(runtime*double(par->nthread)));       printf("Acceptance speed: %.2e particle pairs per core per second\n\n",double(global_counts.used_pairs)/(runtime*double(par->nthread)));   
    
    char this_out[5];
    sprintf(this_out,"full");
    global_counts.save_counts(this_out);
    printf("Printed counts to file as %s%s_power_counts_n%d_m%d_full.txt\n", par->out_file,par->out_string,nbin, mbin);
    }
};

#endif
