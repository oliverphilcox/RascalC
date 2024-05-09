// rescale_correlation script to rescale the input correlation function to correct for errors due to the interpolation function requiring the bin center value and the inputs giving the mean-bin positions only.

#ifndef RESCALE_CORRELATION_H
#define RESCALE_CORRELATION_H

#include <algorithm>
#include <inttypes.h>

class correlation_integral{
public:
    Float *cf_estimate; // estimated correlation function in each bin
    Float *rr_estimate; // RR bin count estimate
    uint64 *binct; // array to accumulate bin counts
    CorrelationFunction *old_cf; // input correlation function
    Integrals *integral; // integrals class
private:
    int nbin,mbin;
    bool rad;
    Float dmu,mumin,mumax,*r_high,*r_low; // parameter file values

public:
    correlation_integral(Parameters *par, CorrelationFunction *_cf){

        // Read in correlation function
        old_cf = new CorrelationFunction(_cf);

        // Empty integrals class to use class functions

        integral = new Integrals(); // integrals object

        // Important parameters
        nbin = par->nbin_cf;
        mbin = par->mbin_cf;
        mumin = par->mumin;
        mumax = par->mumax;
        r_high = par->radial_bins_high_cf;
        r_low = par->radial_bins_low_cf;
        dmu = (mumax-mumin)/mbin; // assume same mu ranges for correlation function and output covariance matrix
        rad=mbin==1&&dmu==1.;

        // Allocate memory;
        int ec=0;
        ec+=posix_memalign((void **) &cf_estimate, PAGE, sizeof(double)*nbin*mbin);
        ec+=posix_memalign((void **) &rr_estimate, PAGE, sizeof(double)*nbin*mbin);
        ec+=posix_memalign((void **) &binct, PAGE, sizeof(uint64)*nbin*mbin);
        assert(ec==0);

        // Reset functions
        reset();

    }

    ~correlation_integral(){
        free(cf_estimate);
        free(rr_estimate);
        free(binct);
    }
private:
    inline int getbin(Float r, Float mu){
        // Linearizes 2D indices - needs CORRELATION FUNCTION bins here, not covariance binning (i.e. should extend to zero)

        // First define which r bin we are in;
        Float* r_higher = std::upper_bound(r_high, r_high + nbin, r); // binary search for r_high element higher than r
        int which_bin = r_higher - r_high; // bin index is pointer difference; will be nbin if value not found, i.e. if we are above top bin
        if (which_bin < nbin) // safety check unless we are above top bin already
            if (r < r_low[which_bin]) // r < r_high[which_bin] is guaranteed above so only need to check that r >= r_low[which_bin]
                which_bin = -1; // if not then no bin fits the bill
        return which_bin*mbin + floor((mu-mumin)/dmu);
    }

public:
    void compile_integral(const Particle* pi_list, const int* prim_ids, int pln, const Particle pj, const int pj_id, const double prob){
        // Accumulates the correlation function integral xi

        Float tmp_weight, tmp_xi, rij_mag, rij_mu,xi_contrib,rr_contrib;
        Particle pi;
        int tmp_bin;

        for(int i=0;i<pln;i++){ // Iterate over particle in pi_list
            if(prim_ids[i]==pj_id){
                continue; // don't self-count
            }
            pi = pi_list[i]; // first particle
            integral->cleanup_l(pi.pos,pj.pos,rij_mag,rij_mu); // define |r_ij| and ang(r_ij)
            tmp_bin = getbin(rij_mag, rij_mu); // bin for each particle

            if ((tmp_bin<0)||(tmp_bin>=nbin*mbin)){
                continue; // if not in correct bin
            }

            tmp_weight = pi.w*pj.w; // product of weights
            tmp_xi = old_cf->xi(rij_mag, rij_mu); // correlation function for i-j

            // Now compute the integral:
            rr_contrib = tmp_weight / prob; // RR_a contribution
            xi_contrib = rr_contrib*tmp_xi; // c2 contribution

            // Add to local integral counts:
            cf_estimate[tmp_bin]+=xi_contrib;
            rr_estimate[tmp_bin]+=rr_contrib;
            binct[tmp_bin]++;
        }
    }

    void normalize(Float norm1, Float norm2, Float n_pairs){
        // Normalize the accumulated integrals and reweight by N_gal/N_rand
        double corrf2 = norm1*norm2; // correction factor
        bool fail_flag = false;
        uint64 min_binct = binct[0];
        int min_binct_bin = 0;

        for(int i = 0; i<nbin*mbin;i++){
            if (binct[i] < min_binct) {
                min_binct = binct[i];
                min_binct_bin = i;
            }
            if (binct[i] == 0) {
                fprintf(stderr, "No pairs sampled in correlation function bin %d (radial %d, angular %d)!\n", i, i / mbin, i % mbin);
                fail_flag = true;
            }
            else if (rr_estimate[i] == 0) {
                fprintf(stderr, "RR estimate in correlation function bin %d (radial %d, angular %d) is zero!\n", i, i / mbin, i % mbin);
                fail_flag = true;
            }
            rr_estimate[i]/=(n_pairs*corrf2);
            cf_estimate[i]/=(n_pairs*corrf2*rr_estimate[i]);
        }
        if (fail_flag) {
            fprintf(stderr, "Failed to obtain a correlation function estimate for rescaling.\nTry increasing N2 and/or number of integration loops, or using coarser correlation function binning, or providing denser randoms.\n");
            exit(1);
        }
        printf("Smallest number of pairs sampled is %llu in correlation function bin %d (radial %d, angular %d). If the number is not large, consider increasing N2 and/or number of integration loops, or using coarser correlation function binning, or providing denser randoms.\n", (unsigned long long) min_binct, min_binct_bin, min_binct_bin / mbin, min_binct_bin % mbin);
    }

    void sum(correlation_integral* corr){
        // Sum individual thread contributions
        for(int i=0;i<nbin*mbin;i++){
            cf_estimate[i]+=corr->cf_estimate[i];
            rr_estimate[i]+=corr->rr_estimate[i];
            binct[i] += corr->binct[i];
        }
    }

    void reset(){
        for(int i=0;i<nbin*mbin;i++){
            cf_estimate[i]=0;
            rr_estimate[i]=0;
            binct[i]=0;
        }
    }

};

class rescale_correlation{
    // class to rescale the correlation function to ensure that the binned correlation functions are reproduced

private:
    int nbin,mbin;
    Integrals *integral;
    compute_integral *compute;
    Float *r_centers, *mu_centers, *new_xi_array;
public:
    CorrelationFunction *new_cf;

public:
    rescale_correlation(){};

    ~rescale_correlation(){
        // Destructor function
        free(r_centers);
        free(mu_centers);
        free(new_xi_array);
    }

    rescale_correlation(Parameters *par){
        // Constructor function
        nbin = par->nbin_cf;
        mbin = par->mbin_cf;

        integral = new Integrals();
        compute = new compute_integral();

        // Create bin centers;
        r_centers = (Float *)malloc(sizeof(Float)*nbin);
        mu_centers = (Float *)malloc(sizeof(Float)*mbin);
        new_xi_array = (Float *)malloc(sizeof(Float)*nbin*mbin);
        for(int j=0;j<mbin;j++) mu_centers[j]=(0.5+j)/mbin;
        for(int i=0;i<nbin;i++) r_centers[i]= 0.5*(par->radial_bins_low_cf[i]+par->radial_bins_high_cf[i]);
        }


public:
    void refine_wrapper(Parameters *par, Grid all_grid[], CorrelationFunction all_cf[], RandomDraws all_rd[], int number_xi){
        // Wrapper to refine each individual correlation function. This updates the correlation function and random draw classes.

        CorrelationFunction true_cf;
        int grid1_index[3] = {0,1,0}, grid2_index[3] = {0,1,1};

        for(int index=0;index<number_xi;index++){ // iterate over correlation functions
            if (par->cf_loops>0){
                fprintf(stderr, "\nRefining correlation function %d of %d.\n",index+1,number_xi);
            }
                true_cf.copy_function(&all_cf[index]); // store initial correlation function
                for(int n_refine=0;n_refine<par->cf_loops;n_refine++){ // refine cf_loops times per correlation function
                    // Rescale correlation function
                    CorrelationFunction output = rescale_xi(par, &all_grid[grid1_index[index]], &all_grid[grid2_index[index]], &all_cf[index], &true_cf, &all_rd[index],n_refine);
                    // Update correlation function
                    all_cf[index].copy_function(&output);
                }
                // Only update random draws on the final iteration for speed
                RandomDraws tmp_rd(&all_cf[index],par,NULL,0);
                all_rd[index].copy(&tmp_rd);
        }
    }

    CorrelationFunction rescale_xi(Parameters *par, Grid *grid1, Grid *grid2, CorrelationFunction *old_cf, CorrelationFunction *true_cf, RandomDraws *rd, int index){
        // Rescale the xi function by computing the binned xi from some estimate and comparing it to the known value. Return a correlation function object with an updated estimate of the xi at the bin-centers.

        // gsl and random class setup
        unsigned long seed_step = (unsigned long)(std::numeric_limits<uint32_t>::max()) / (unsigned long)(par->max_loops); // restrict the seeds to 32 bits to avoid possible collisions, as most GSL random generators only accept 32-bit seeds and truncate the higher bits. unsigned long is at least 32 bits but can be more.
        unsigned long seed_shift = par->seed % seed_step; // for each loop, it will be seed = seed_step * n_loops + seed_shift, where n_loops goes from 0 to max_loops-1 => no overflow and guaranteed unique for each thread
        if (par->random_seed) {
            std::random_device urandom("/dev/urandom");
            std::uniform_int_distribution<unsigned long> dist(0, seed_step-1); // distribution of integers with these limits, inclusive
            seed_shift = dist(urandom); // for each thread, it will be seed = seed_step * n_loops + seed_shift, where n_loops goes from 0 to max_loops-1 => no overflow and guaranteed unique for each thread
        }
        gsl_rng_env_setup(); // initialize gsl rng
        correlation_integral full_xi_function(par,old_cf); // full correlation function class
        uint64 used_pairs=0;

#ifdef OPENMP
#pragma omp parallel firstprivate(steps,par,grid1, grid2, old_cf) shared(gsl_rng_default,rd) reduction(+:used_pairs)
        { // start parallel loop
        // Decide thread
        int thread = omp_get_thread_num();
        assert(omp_get_num_threads()<=par->nthread);
        if (thread==0) fprintf(stderr, "# Computing correlation function iteration %d of %d on %d threads.\n", index+1,par->cf_loops, omp_get_num_threads());
    #else
        { // start empty loop
        int thread = 0;
        fprintf(stderr, "# Computing correlation function iteration %d of %d single threaded.\n",index+1,par->cf_loops);
    #endif

        // Define local thread variables
        Particle *prim_list,particle_j;
        int *prim_ids,pln,sln,pid_j,prim_id_1D,x;
        integer3 delta2,prim_id,sec_id;
        double p2;
        Float3 cell_sep2;
        gsl_rng* locrng = gsl_rng_alloc(gsl_rng_default); // one rng per thread, seed will be set later inside the loop

        correlation_integral thread_xi_function(par, old_cf);

        // Allocate particle position memory
        int ec=0;
        ec+=posix_memalign((void **) &prim_list, PAGE, sizeof(Particle)*grid1->maxnp);
        ec+=posix_memalign((void **) &prim_ids, PAGE, sizeof(int)*grid1->maxnp);
        assert(ec==0);

    // start first loop
#ifdef OPENMP
#pragma omp for schedule(dynamic)
#endif
        for(int n_loops = 0; n_loops<par->max_loops; n_loops++){
            // Set/reset the RNG seed based on loop number instead of thread number to reproduce results with different number of threads but other parameters kept the same. Individual subsamples may differ because they are accumulated/written in order of loop completion which may depend on external factors at runtime, but the final integrals should be the same.
            gsl_rng_set(locrng, seed_step * (unsigned long)n_loops + seed_shift); // the second number, seed, will not overflow and will be unique for each loop in one run. Here, seed_shift is a random number between 0 and seed_step-1, inclusively.
            // Seed clashes between different runs seem less likely than for the old formula, seed = steps * (nthread+1), with steps being a random number between 1 and UINT_MAX (or ULONG_MAX / nthread) inclusively - there, steps of one run could be a multiple of steps of the other resulting in the same seeds for some threads, which is somewhat more likely than getting the same random number.
            for(int n1=0;n1<grid1->nf;n1++){
                // Pick first particle
                prim_id_1D = grid1-> filled[n1]; // 1d ID for cell i
                prim_id = grid1->cell_id_from_1d(prim_id_1D); // define first cell
                pln = compute->particle_list(prim_id_1D, prim_list, prim_ids, grid1); // update list of particles and number of particles

                if(pln==0) continue; // skip if empty

                used_pairs+=pln*par->N2;

                // Loop over N2 j cells;
                for (int n2=0;n2<par->N2;n2++){
                    delta2 = rd->random_xidraw(locrng,&p2);
                    sec_id = prim_id+delta2;
                    cell_sep2 = grid2->cell_sep(delta2);
#ifdef THREE_PCF
                    x = compute->draw_particle(sec_id,particle_j,pid_j,cell_sep2,grid2,sln,locrng);
#else
                    x = compute->draw_particle_without_class(sec_id,particle_j,pid_j,cell_sep2,grid2,sln,locrng);
#endif
                    if(x==1) continue; // skip if error
                    p2*=1./(grid1->np*(double)sln);
                    // Compute contributions to correlation function
                    thread_xi_function.compile_integral(prim_list, prim_ids, pln, particle_j, pid_j, p2);
                }
            }
#ifdef OPENMP
#pragma omp critical
#endif
{
        full_xi_function.sum(&thread_xi_function);
        thread_xi_function.reset();
}
        } // end of loops

        // Free up memory
        free(prim_list);
        free(prim_ids);

        } // end OPENMP loop

        // Normalize integral
        full_xi_function.normalize(grid1->norm,grid2->norm,(Float)used_pairs);

        // Compute ratios
        Float true_xi,old_xi;
        for(int i=0;i<nbin;i++){
          for(int j=0;j<mbin;j++){
                 old_xi = old_cf->xi(r_centers[i],mu_centers[j]);
                true_xi = true_cf->xi(r_centers[i],mu_centers[j]);
                // Rescale value
                new_xi_array[i*mbin+j]=true_xi/full_xi_function.cf_estimate[i*mbin+j]*old_xi;
            }
        }

        // Now write to cf function
        new_cf = new CorrelationFunction(new_xi_array, r_centers, mu_centers, nbin, mbin);
        return new_cf;
    }

};

#endif
