// triple_counts.cpp -- Oliver Philcox started March 2019 based on grid_multipoles.cpp of Daniel Eisenstein.

#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
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
    #include "../modules/parameters.h"
    #include "../modules/cell_utilities.h"
    #include "../modules/grid.h"
    #include "../modules/correlation_function.h"
    #include "../modules/random_draws.h"
    #include "../modules/driver.h"

// Get the correlation function into the integrator
CorrelationFunction * RandomDraws::corr;

STimer TotalTime;

// ================================ main() =============================

// Class to hold triple counts
class TripleCounts{
private:
    int nbin,mbin;
    Float rmin,rmax,mumin,mumax,dmu; //Ranges in r and mu
    Float *r_high, *r_low; // Max and min of each radial bin
    Float *RRR; // Arrays to accumulate integrals
    char* out_file;
    bool box,rad=0; // Flags to decide whether we have a periodic box + if we have a radial correlation function only
    
public:
    TripleCounts(Parameters *par){
        // Initializer
        init(par);
    }
    
    void init(Parameters *par){
        nbin = par->nbin;
        mbin = par->mbin; // number of mu bins
        out_file = par->out_file; // output directory

        int ec=0;
        // Initialize the binning
        ec+=posix_memalign((void **) &RRR, PAGE, sizeof(double)*nbin*nbin*mbin);
        assert(ec==0);
        reset();
        
        box=par->perbox;

        rmax=par->rmax;
        rmin=par->rmin;

        mumax=par->mumax;
        mumin=par->mumin;
        
        r_high = par->radial_bins_high;
        r_low = par->radial_bins_low;
        
        dmu=(mumax-mumin)/mbin;

        rad=mbin==1&&dmu==1.;
    }
    
    ~TripleCounts(){
        free(RRR);
    }
    void reset(){
        for(int j=0;j<nbin*nbin*mbin;j++) RRR[j] = 0;
    }
    
    inline int get_radial_bin(Float r){
        // Computes radial bin
        
        int which_bin = -1; // default if outside bins
        for(int i=0;i<nbin;i++){
            if((r>r_low[i])&&(r<r_high[i])){
                which_bin=i;
                break;
            }
            if((i==nbin-1)&&(r>r_high[i])){
                which_bin=nbin; // if above top bin
            }
        }
        return which_bin;
        }
    inline int which_ang_bin(Float mu){
        // Computes angular bin
        return floor((mu-mumin)/dmu);
    }
    
     void triangle_bins(Float3 p1, Float3 p2, Float3 p3, Float* norm, Float* ang,Float norm_in, int preload){
        // Compute side lengths and internal angles of a triangle of positions
        // Indices are chosen such that norm[x] is opposite to ang[x]
        // NB we assume norm(p2,p3) [if preload==0] or norm(1,2) [if preload==2] is already computed here for speed.
        
        Float3 pos12 = p1-p2, pos23 = p2-p3, pos13 = p1-p3;
        
        // Compute lengths
        if(preload==0) norm[0] = norm_in;
        else norm[0]=pos23.norm();
        norm[1] = pos13.norm();
        if(preload==2) norm[2] = norm_in;
        else norm[2] = pos12.norm();
        
        // Compute internal angles
        ang[0] = pos12.dot(pos13)/norm[2]/norm[1];
        ang[1] = -pos12.dot(pos23)/norm[2]/norm[0];
        ang[2] = pos13.dot(pos23)/norm[1]/norm[0];
    }
    
    int all_bins(Float norm[3], Float bin_ab[6], const int bin_in,int preload){
        // Compute all triangle bin combinations and save as an array. We use the j-k bin as an input here
        
        int error=0; // to check if there are any correct bins
        Float bin_ij=0, bin_ik=0, bin_jk=0;
        
        // Load bins and check if correct bins are sampled
        if(preload==0) bin_jk=bin_in;
        else bin_jk = get_radial_bin(norm[0]);
        bin_ik = get_radial_bin(norm[1]);     
        if(preload==2) bin_ij=bin_in;
        else bin_ij = get_radial_bin(norm[2]);
        
        if ((bin_ij<0)||(bin_ij>=nbin)){
            bin_ij = -1;
            error++;
        }
        if((bin_jk<0)||(bin_jk>=nbin)){
            bin_jk = -1;
            error++;
        }
        if ((bin_ik<0)||(bin_ik>=nbin)){
            bin_ik = -1;
            error++;
        }
        
        // Now save all 3 combinations of these bins
        // Set 1
        bin_ab[0] = bin_ij;
        bin_ab[1] = bin_ik;
        // Set 2
        bin_ab[2] = bin_ij;
        bin_ab[3] = bin_jk;
        // Set 3
        bin_ab[4] = bin_ik;
        bin_ab[5] = bin_jk;
        
        return error;
    }
    
    
    inline void update_counts(const Particle* prim_list,const int* prim_ids, int pln, const Particle particle_j, const Particle particle_k, const int pid_j, const int pid_k, const double prob){
        
        Float norm_ijk[3], ang_ijk[3], bins_ijk[6];
        Particle pi;
        Float3 pjk;
        Float norm_jk;
        int bin_jk,out_bin,bin_1,bin_2,tmp_ang_bin;            
        
        if(pid_k==pid_j) return;
        
        // Define j-k norm and bin
        pjk = particle_j.pos-particle_k.pos;
        norm_jk = pjk.norm();
        bin_jk = get_radial_bin(norm_jk);
        
        // No self-counts here
        for(int i=0;i<pln;i++){ // Iterate over particle in pi_list

            // First check for any self-counts
            if((prim_ids[i]==pid_j)||(prim_ids[i]==pid_k)) continue;
            
            pi = prim_list[i]; // first particle
            
            // Define angles and separations
            triangle_bins(pi.pos,particle_j.pos,particle_k.pos,norm_ijk,ang_ijk,norm_jk,0);
            
            // Define radial bins
            int x = all_bins(norm_ijk,bins_ijk,bin_jk,0);
            if(x==3) continue; // skip if no correct bins     
            
            for(int bin_index=0;bin_index<3;bin_index++){
                bin_1 = bins_ijk[bin_index*2];
                if(bin_1==-1) continue;
                bin_2 = bins_ijk[bin_index*2+1];
                if(bin_2==-1) continue; // skip if bad bin
                
                // Now define full bin 
                tmp_ang_bin = which_ang_bin(ang_ijk[bin_index]);
                assert(tmp_ang_bin>=0);
                assert(tmp_ang_bin<mbin);
                out_bin = (bin_1*nbin+bin_2)*mbin+tmp_ang_bin;
                
                // Extra x2 is from symmetry in angle kernel
                RRR[out_bin]+=pi.w*particle_j.w*particle_k.w/prob*2.;
                
            }
        }
    }
    
    void normalize(Float sum_w, Float n_triple){
        for(int i=0;i<nbin*nbin*mbin;i++) RRR[i]/=(pow(sum_w,3)*n_triple);
    }
    
    void save_triples(char* suffix){
        char RRRname[1000];
        snprintf(RRRname, sizeof RRRname, "%s/RRR_counts_n%d_m%d_%s.txt", out_file,nbin, mbin,suffix);
        FILE * RRRFile = fopen(RRRname,"w"); 
        for(int j=0;j<nbin*nbin*mbin;j++) fprintf(RRRFile,"%le\n",RRR[j]);
        
        fflush(NULL);
        fclose(RRRFile);
    }
    
    inline void sum_counts(TripleCounts *cc){
        for(int i=0;i<nbin*nbin*mbin;i++) RRR[i]+=cc->RRR[i];
    }
        
};

class compute_triples{
    
    private:
        uint64 cnt2=0,cnt3=0;
        int nbin, mbin; 
     public:
        int particle_list(int id_1D, Particle* &part_list, int* &id_list, Grid *grid){
            // function updates a list of particles for a 1-dimensional ID. Output is number of particles in list.
            
            Cell cell = grid->c[id_1D]; // cell object 
            int no_particles = 0;
            // copy in list of particles into list
            for (int i = cell.start; i<cell.start+cell.np; i++, no_particles++){
                part_list[no_particles]=grid->p[i];
                id_list[no_particles]=i;
            }
            
        return no_particles;
        }
        
    public:
        int draw_particle(integer3 id_3D, Particle &particle, int &pid, Float3 shift, Grid *grid, int &n_particles, gsl_rng* locrng){
            // Draw a random particle from a cell given the cell ID.
            // This updates the particle and particle ID and returns 1 if error.
            // This is used for k,l cells (with no indication of particle random class)
            int id_1D = grid-> test_cell(id_3D); 
            if(id_1D<0) return 1; // error if cell not in grid
            Cell cell = grid->c[id_1D];
            if(cell.np==0) return 1; // error if empty cell
            pid = floor(gsl_rng_uniform(locrng)*cell.np) + cell.start; // draw random ID
            particle = grid->p[pid]; // define particle
            n_particles = cell.np; // no. of particles in cell 
    #ifdef PERIODIC
            particle.pos+=shift;
    #endif
            return 0;
        }
        
        
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
        compute_triples(Grid *grid, Parameters *par, CorrelationFunction *cf, RandomDraws *rd){
            
            // Define parameters
            nbin = par->nbin; // number of radial bins
            mbin = par->mbin; // number of mu bins
            
            STimer initial, TotalTime; // Time initialization
            initial.Start(); 
            
            //-----------INITIALIZE OPENMP + CLASSES----------
            std::random_device urandom("/dev/urandom");
            std::uniform_int_distribution<unsigned int> dist(1, std::numeric_limits<unsigned int>::max());
            unsigned long int steps = dist(urandom);        
            gsl_rng_env_setup(); // initialize gsl rng

            uint64 tot_pairs=0, tot_triples=0; // global number of particle pairs/triples/quads used (including those rejected for being in the wrong bins)
            uint64 cell_attempt2=0,cell_attempt3=0; // number of j,k,l cells attempted
            uint64 used_cell2=0,used_cell3=0; // number of used j,k,l cells
            
            TripleCounts global_counts(par);
            
            check_threads(par,1); // Define which threads we use

            initial.Stop();
            fprintf(stderr, "Init time: %g s\n",initial.Elapsed());
            printf("# 1st grid filled cells: %d\n",grid->nf);
            printf("# All 1st grid points in use: %d\n",grid->np);
            printf("# Max points in one cell in grid 1%d\n",grid->maxnp);
            fflush(NULL);
            
            TotalTime.Start(); // Start timer
            
#ifdef OPENMP      
    #pragma omp parallel firstprivate(steps,par,grid) shared(global_counts,TotalTime,gsl_rng_default,rd) reduction(+:cell_attempt2,cell_attempt3,used_cell2,used_cell3,tot_pairs,tot_triples)
            { // start parallel loop
            // Decide which thread we are in
            int thread = omp_get_thread_num();
            assert(omp_get_num_threads()<=par->nthread);
            if (thread==0) printf("# Starting integral computation on %d threads.\n", omp_get_num_threads());        
#else
            int thread = 0;
            printf("# Starting integral computation single threaded.\n");
            { // start loop
#endif
            
    //-----------DEFINE LOCAL THREAD VARIABLES
            Particle *prim_list; // list of particles in first cell
            int pln,sln,tln; // number of particles in each cell
            int pid_j, pid_k; // particle IDs particles drawn from j,k,l cell
            Particle particle_j, particle_k; // randomly drawn particle
            //int* bin; // a-b bins for particles
            int* prim_ids; // list of particle IDs in primary cell
            double p2,p3; // probabilities
            int mnp = grid->maxnp; // max number of particles in a grid cell
            Float percent_counter;
            int x, prim_id_1D;
            integer3 delta2, delta3, prim_id, sec_id, thi_id;
            Float3 cell_sep2,cell_sep3;
            
            TripleCounts local_counts(par); // holds triple counts for each thread
            
            gsl_rng* locrng = gsl_rng_alloc(gsl_rng_default); // one rng per thread
            gsl_rng_set(locrng, steps*(thread+1));
            
            // Assign memory for intermediate steps
            int ec=0;
            ec+=posix_memalign((void **) &prim_list, PAGE, sizeof(Particle)*mnp);
            ec+=posix_memalign((void **) &prim_ids, PAGE, sizeof(int)*mnp);
            assert(ec==0);
            
            uint64 loc_used_pairs,loc_used_triples; // local counts of used pairs/triples/quads
            
    //-----------START FIRST LOOP-----------
    #ifdef OPENMP
    #pragma omp for schedule(dynamic)
    #endif
           for (int n_loops = 0; n_loops<par->max_loops; n_loops++){
                percent_counter=0.;
                loc_used_pairs=0; loc_used_triples=0;  
                
                // LOOP OVER ALL FILLED I CELLS
                for (int n1=0; n1<grid->nf;n1++){
                    
                    // Print time left
                    if((float(n1)/float(grid->nf)*100)>=percent_counter){
                        printf("Run %d of %d on thread %d: Using cell %d of %d - %.0f percent complete\n",1+n_loops/par->nthread, int(ceil(float(par->max_loops)/(float)par->nthread)),thread, n1+1,grid->nf,percent_counter);
                        percent_counter+=5.;
                    }
                    
                    // Pick first particle
                    prim_id_1D = grid-> filled[n1]; // 1d ID for cell i 
                    prim_id = grid->cell_id_from_1d(prim_id_1D); // define first cell
                    pln = particle_list(prim_id_1D, prim_list, prim_ids, grid); // update list of particles and number of particles
                    
                    if(pln==0) continue; // skip if empty
                    
                    loc_used_pairs+=pln*par->N2;
                    loc_used_triples+=pln*par->N2*par->N3;
                    
                    // LOOP OVER N2 J CELLS
                    for (int n2=0; n2<par->N2; n2++){
                        cell_attempt2+=1; // new cell attempted
                        
                        // Draw second cell from i weighted by 1/r^2
                        delta2 = rd->random_cubedraw(locrng, &p2); 
                        // p2 is the ratio of sampling to true pair distribution here
                        sec_id = prim_id + delta2;
                        cell_sep2 = grid->cell_sep(delta2);
                        x = draw_particle(sec_id, particle_j, pid_j,cell_sep2, grid, sln, locrng);
                        if(x==1) continue; // skip if error
                        
                        used_cell2+=1; // new cell accepted
                        
                        // For all particles
                        p2*=1./(grid->np*(double)sln); // probability is divided by total number of i particles and number of particles in cell
                        
                        // LOOP OVER N3 K CELLS
                        for (int n3=0; n3<par->N3; n3++){
                            cell_attempt3+=1; // new third cell attempted
                            
                            // Draw third cell from i weighted by 1/r^2
                            delta3 = rd->random_cubedraw(locrng, &p3); 
                            thi_id = prim_id + delta3;
                            cell_sep3 = grid->cell_sep(delta3); 
                            x = draw_particle(thi_id,particle_k,pid_k,cell_sep3,grid,tln,locrng); // draw from third grid
                            if(x==1) continue; 
                            
                            used_cell3+=1; // new third cell used
                            
                            p3*=p2/(double)tln; // update probability
                            
                            local_counts.update_counts(prim_list,prim_ids,pln,particle_j,particle_k,pid_j,pid_k,p3);                           
                            
                        }
                    }
                }
                
                // Update used pair/triple counts
                tot_pairs+=loc_used_pairs;
                tot_triples+=loc_used_triples;
                
                #ifdef OPENMP
    #pragma omp critical // only one processor can access at once
    #endif
            {
                if ((n_loops+1)%par->nthread==0){ // Print every nthread loops
                    TotalTime.Stop(); // interrupt timing to access .Elapsed()
                    int current_runtime = TotalTime.Elapsed();
                    int remaining_time = current_runtime/((n_loops+1)/par->nthread)*(par->max_loops/par->nthread-(n_loops+1)/par->nthread);  // estimated remaining time
                    fprintf(stderr,"\nFinished integral loop %d of %d after %d s. Estimated time left:  %2.2d:%2.2d:%2.2d hms, i.e. %d s.\n",n_loops+1,par->max_loops, current_runtime,remaining_time/3600,remaining_time/60%60, remaining_time%60,remaining_time);
                    
                    TotalTime.Start(); // Restart the timer
                }
                
                // Sum up counts
                global_counts.sum_counts(&local_counts);
                
                // Reset local counts
                local_counts.reset();
            }
           } // end cycle loop
           
           // Free allocated memory
           free(prim_list);
           
            } // end OPENMP loop
                
        // ----- REPORT + SAVE OUTPUT -------
        
        // Normalize counts
        global_counts.normalize(grid->sum_weights,tot_triples);
        
        int runtime = TotalTime.Elapsed();
        printf("\n\nTRIPLE COUNTS COMPLETE\n"); 
        fprintf(stderr, "\nTotal process time for %.2e sets of cells and %.2e triples of particles: %d s, i.e. %2.2d:%2.2d:%2.2d hms\n", double(used_cell3),double(tot_triples),runtime, runtime/3600,runtime/60%60,runtime%60);
        printf("We tried %.2e pairs and %.2e triples.\n",double(cell_attempt2),double(cell_attempt3));
        printf("Of these, we accepted %.2e pairs and %.2e triples.\n",double(used_cell2),double(used_cell3));
        printf("We sampled %.2e pairs and %.2e triples.\n",double(tot_pairs),double(tot_triples));
        
        printf("\nTrial speed: %.2e triples per core per second\n",double(tot_triples)/(runtime*double(par->nthread)));
        
        char out_string[5];
        sprintf(out_string,"full");
        global_counts.save_triples(out_string);
        printf("Printed triple counts to file\n");
        fflush(NULL);
        return;
        }
};


int main(int argc, char *argv[]) {
    
	Parameters par=Parameters(argc,argv);
        
    int max_no_functions=1; // required number of xi / random_draws / jackknife_weight functions
    int no_fields=1; // number of different fields used
    if(par.multi_tracers==true){
        max_no_functions=3;
        no_fields=2;
    }
    
    // Now read in particles to grid:
    Grid all_grid[no_fields]; // create empty grids
    
    for(int index=0;index<no_fields;index++){
        Float3 shift;
        Particle *orig_p;
        if (!par.make_random){
            char *filename;
            if(index==0) filename=par.fname;
            else filename=par.fname2;
            orig_p = read_particles(par.rescale, &par.np, filename, par.rstart, par.nmax);
            assert(par.np>0);
            par.perbox = compute_bounding_box(orig_p, par.np, par.rect_boxsize, par.cellsize, par.rmax, shift, par.nside);
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
        Float nofznorm=par.nofznorm;
        if(index==1) nofznorm=par.nofznorm2;
        Grid tmp_grid(orig_p, par.np, par.rect_boxsize, par.cellsize, par.nside, shift, nofznorm);

        Float grid_density = (double)par.np/tmp_grid.nf;
        printf("\n RANDOM CATALOG %d DIAGNOSTICS:\n",index+1);
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
        printf("# %d particles in use, %d with positive weight\n", tmp_grid.np, tmp_grid.np_pos);
        printf("# Weights: Positive particles sum to %f\n", tmp_grid.sumw_pos);
        printf("#          Negative particles sum to %f\n", tmp_grid.sumw_neg);

        // Now save grid to global memory:
        all_grid[index].copy(&tmp_grid);
        
        free(orig_p); // Particles are now only stored in grid
        
        fflush(NULL);
    }
    
    // Now define all possible correlation functions and random draws:
    CorrelationFunction all_cf[max_no_functions];
    RandomDraws all_rd[max_no_functions];
    
    CorrelationFunction tmp_cf(par.corname,par.mbin,par.mumax-par.mumin);
    all_cf[0].copy_function(&tmp_cf);
    RandomDraws tmp_rd(&tmp_cf,&par,NULL,0);
    all_rd[0].copy(&tmp_rd);
    
    // Run main modules
    compute_triples(&all_grid[0],&par,&all_cf[0],&all_rd[0]);
   
    // END
    rusage ru;
    getrusage(RUSAGE_SELF, &ru);

    fprintf(stderr,"# user CPU time used: %ld \n# system CPU time used: %ld \n# maximum resident set size: %ld \n# integral shared memory size: %ld \n# integral unshared data size: %ld \n# integral unshared stack size: %ld \n# page reclaims (soft page faults): %ld \n# page faults (hard page faults): %ld \n# swaps: %ld \n# block input operations: %ld \n# block output operations: %ld \n# IPC messages sent: %ld \n# IPC messages received: %ld \n# signals received: %ld \n# voluntary context switches: %ld \n# involuntary context switches: %ld \n",ru.ru_utime.tv_sec,ru.ru_stime.tv_sec,ru.ru_maxrss,ru.ru_ixrss,ru.ru_idrss,ru.ru_isrss,ru.ru_minflt,ru.ru_majflt,ru.ru_nswap,ru.ru_inblock,ru.ru_oublock,ru.ru_msgsnd,ru.ru_msgrcv,ru.ru_nsignals,ru.ru_nvcsw,ru.ru_nivcsw);

    return 0;
}


