 // This class computes the contributions to the C_ab integral iterating over cells and particles. It is a heavily modified version of code by Alex Wiegand.

    #ifndef COMPUTE_INTEGRAL_H
    #define COMPUTE_INTEGRAL_H

// Choose relevant class to hold integrals
#ifdef LEGENDRE
    #include "integrals_legendre.h"
#elif defined POWER
    #include "integrals_legendre_power.h"
#else
    #include "integrals.h"
#endif
    class compute_integral{

    private:
        uint64 cnt2=0,cnt3=0,cnt4=0;
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

    private:
        int draw_particle(integer3 id_3D, Particle &particle, int &pid, Float3 shift, Grid *grid, int &n_particles, gsl_rng* locrng, int &n_particles1, int &n_particles2){
            // Draw a random particle from a cell given the cell ID.
            // This updates the particle and particle ID and returns 1 if error.

            int id_1D = grid-> test_cell(id_3D);
            if(id_1D<0) return 1; // error if cell not in grid
            Cell cell = grid->c[id_1D];
            if(cell.np==0) return 1; // error if empty cell
            pid = floor(gsl_rng_uniform(locrng)*cell.np) + cell.start; // draw random ID
            particle = grid->p[pid]; // define particle
            n_particles = cell.np; // no. of particles in cell
            n_particles1 = cell.np1; // no. particles in cell partition 1
            n_particles2 = cell.np2;
    #ifdef PERIODIC
            particle.pos+=shift;
    #endif
            return 0;
        }

    public:
        int draw_particle_without_class(integer3 id_3D, Particle &particle, int &pid, integer3 shift, Grid *grid, int &n_particles, gsl_rng* locrng){
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
    private:
        CorrelationFunction* which_cf(CorrelationFunction all_cf[], int Ia, int Ib){
            // Returns the relevant correlation function for two input field indices
            if((Ia==1)&(Ib==1)) return &all_cf[0];
            else if ((Ia==2)&&(Ib==2)) return &all_cf[1];
            else return &all_cf[2];
        }
        RandomDraws* which_rd(RandomDraws all_rd[], int Ia, int Ib){
            // Returns the relevant correlation function for two input field indices
            if((Ia==1)&(Ib==1)) return &all_rd[0];
            else if ((Ia==2)&&(Ib==2)) return &all_rd[1];
            else return &all_rd[2];
        }
        Grid* which_grid(Grid all_grid[], int Ia){
            // Returns the relevant correlation function for two input field indices
            if(Ia==1) return &all_grid[0];
            else return &all_grid[1];
        }

#if (defined LEGENDRE || defined POWER)
        SurveyCorrection* which_survey(SurveyCorrection all_survey[], int Ia,int Ib){
            // Returns relevant survey correction function for two input field indices
            if((Ia==1)&(Ib==1)) return &all_survey[0];
            else if ((Ia==2)&&(Ib==2)) return &all_survey[1];
            else return &all_survey[2];
        }
#else
        JK_weights* which_JK(JK_weights all_JK[], int Ia, int Ib){
            // Returns the relevant correlation function for two input field indices
            if((Ia==1)&(Ib==1)) return &all_JK[0];
            else if ((Ia==2)&&(Ib==2)) return &all_JK[1];
            else return &all_JK[2];
        }
#endif

    public:
        compute_integral(){};

#if (defined LEGENDRE || defined POWER)
        compute_integral(Grid all_grid[], Parameters *par, CorrelationFunction all_cf[], RandomDraws all_rd[], SurveyCorrection all_survey[], int I1, int I2, int I3, int I4, int iter_no){
#else
        compute_integral(Grid all_grid[], Parameters *par, JK_weights all_JK[], CorrelationFunction all_cf[], RandomDraws all_rd[], int I1, int I2, int I3, int I4, int iter_no){
#endif
            // MAIN FUNCTION TO COMPUTE INTEGRALS

            int tot_iter=1; // total number of iterations
            if(par->multi_tracers==true) tot_iter=7;

            // Define relevant grids
            Grid *grid1 = which_grid(all_grid,I1);
            Grid *grid2 = which_grid(all_grid,I2);
            Grid *grid3 = which_grid(all_grid,I3);
            Grid *grid4 = which_grid(all_grid,I4);

            // Define relevant correlation functions
            CorrelationFunction *cf12 = which_cf(all_cf,I1,I2);
            CorrelationFunction *cf13 = which_cf(all_cf,I1,I3);
            CorrelationFunction *cf24 = which_cf(all_cf,I2,I4);

            // Define relevant random draw classes:
            RandomDraws *rd13 = which_rd(all_rd,I1,I3);
            RandomDraws *rd24 = which_rd(all_rd,I2,I4);

#ifdef POWER
            // Define relevant survey correction factor
            SurveyCorrection *survey_corr_12 = which_survey(all_survey,I1,I2);
            SurveyCorrection *survey_corr_23 = which_survey(all_survey,I2,I3);
            SurveyCorrection *survey_corr_34 = which_survey(all_survey,I3,I4);
#elif defined LEGENDRE
            // Define relevant survey correction factor
            SurveyCorrection *survey_corr_12 = which_survey(all_survey,I1,I2);
            SurveyCorrection *survey_corr_23 = which_survey(all_survey,I2,I3);
            SurveyCorrection *survey_corr_34 = which_survey(all_survey,I3,I4);
#else
            // Define relevant jackknife JK_weights
            JK_weights *JK12 = which_JK(all_JK,I1,I2);
            JK_weights *JK23 = which_JK(all_JK,I2,I3);
            JK_weights *JK34 = which_JK(all_JK,I3,I4);
#endif

            nbin = par->nbin; // number of radial bins
            mbin = par->mbin; // number of mu bins

            STimer* LoopTimes = new STimer[par->max_loops];
            STimer initial, TotalTime; // Time initialization
            initial.Start();

#ifdef JACKKNIFE
    // --------Compute product weights----------------------
            int nbins=nbin*mbin;
            Float *product_weights12_12, *product_weights12_23, *product_weights12_34; // arrays to get products of jackknife weights to avoid recomputation
            product_weights12_12 = JK12->product_weights;

            // ------ w_12*w_34 weight:: -------

            // Read in product weights if already computed
            if(((I1==I3)&(I2==I4))||((I1==I4)&(I2==I3))) product_weights12_34 = JK12->product_weights;
            // Otherwise compute weights from scratch
            else{
                // Allocate memory first:
                int ec=0;
                ec+=posix_memalign((void **) &product_weights12_34, PAGE, sizeof(Float)*nbins*nbins);
                assert(ec==0);
                // Now compute products
                int partial_bin=0,partial_bin2;
                Float this_weight;
                for (int x=0;x<JK12->n_JK_filled;x++){
                    partial_bin2=0.;
                    for(int bin_a=0; bin_a<nbins;bin_a++){
                        this_weight = JK12->weights[partial_bin+bin_a];
                        for(int bin_b=0;bin_b<nbins;bin_b++){
                            product_weights12_34[partial_bin2+bin_b]+=this_weight*JK34->weights[partial_bin+bin_b];
                        }
                        partial_bin2+=nbins;
                    }
                    partial_bin+=nbins;
                }
            }
            // ---- w_12*w_23 weight:: ------------
            // Read in product weights if already computed
            if(I1==I3) product_weights12_23 = JK12->product_weights;
            else{
                // Allocate memory first
                int ec=0;
                ec+=posix_memalign((void **) &product_weights12_23, PAGE, sizeof(Float)*nbins*nbins);
                assert(ec==0);
                // Now compute products
                int partial_bin=0,partial_bin2;
                Float this_weight;
                for (int x=0;x<JK12->n_JK_filled;x++){
                    partial_bin2 = 0.;
                    for(int bin_a=0;bin_a<nbins;bin_a++){
                        this_weight = JK12->weights[partial_bin+bin_a];
                        for(int bin_b=0;bin_b<nbins;bin_b++){
                            product_weights12_23[partial_bin2+bin_b]+=this_weight*JK23->weights[partial_bin+bin_b];
                        }
                        partial_bin2+=nbins;
                    }
                    partial_bin+=nbins;
                }
            }

            printf("Computed relevant product weights\n");
#endif

    //-----------INITIALIZE OPENMP + CLASSES----------
            std::random_device urandom("/dev/urandom");
            unsigned long seed_step = std::numeric_limits<unsigned long>::max() / (unsigned long)(par->nthread);
            std::uniform_int_distribution<unsigned long> dist(0, seed_step-1); // distribution of integers with these limits, inclusive
            unsigned long seed_shift = dist(urandom); // for each thread, it will be seed = seed_step * thread + seed_shift, where thread goes from 0 to nthread-1 => no overflow and guaranteed unique for each thread

            gsl_rng_env_setup(); // initialize gsl rng
            int completed_loops = 0; // counter of completed loops, since may not finish in index order
            uint64 used_pairs_per_sample = 0, used_triples_per_sample = 0, used_quads_per_sample = 0; // pair, triple and quad counts for the normalization of the current output sample
#if (defined LEGENDRE || defined POWER)
            Integrals sumint(par, cf12, cf13, cf24, I1, I2, I3, I4,survey_corr_12,survey_corr_23,survey_corr_34); // total integral
            Integrals outint(par, cf12, cf13, cf24, I1, I2, I3, I4, survey_corr_12, survey_corr_23, survey_corr_34); // current output integral
#elif defined JACKKNIFE
            Integrals sumint(par, cf12, cf13, cf24, JK12, JK23, JK34, I1, I2, I3, I4,product_weights12_12,product_weights12_23,product_weights12_34); // total integral
            Integrals outint(par, cf12, cf13, cf24, JK12, JK23, JK34, I1, I2, I3, I4, product_weights12_12,product_weights12_23, product_weights12_34); // current output integral
#else
            Integrals sumint(par, cf12, cf13, cf24, JK12, JK23, JK34, I1, I2, I3, I4); // total integral
            Integrals outint(par, cf12, cf13, cf24, JK12, JK23, JK34, I1, I2, I3, I4); // current output integral
#endif
            uint64 tot_pairs=0, tot_triples=0, tot_quads=0; // global number of particle pairs/triples/quads used (including those rejected for being in the wrong bins)
            uint64 cell_attempt2=0,cell_attempt3=0,cell_attempt4=0; // number of j,k,l cells attempted
            uint64 used_cell2=0,used_cell3=0,used_cell4=0; // number of used j,k,l cells

            initial.Stop();
            fprintf(stderr, "Init time: %g s\n",initial.Elapsed());
            printf("# 1st grid filled cells: %d\n",grid1->nf);
            printf("# All 1st grid points in use: %d\n",grid1->np);
            printf("# Max points in one cell in grid 1%d\n",grid1->maxnp);
            fflush(NULL);

            TotalTime.Start(); // Start timer

#ifdef OPENMP

#if (defined LEGENDRE || defined POWER)
    #pragma omp parallel firstprivate(seed_step,seed_shift,par,grid1,grid2,grid3,grid4,cf12,cf13,cf24) shared(sumint,outint,completed_loops,used_pairs_per_sample,used_triples_per_sample,used_quads_per_sample,TotalTime,LoopTimes,gsl_rng_default,rd13,rd24) reduction(+:cell_attempt2,cell_attempt3,cell_attempt4,used_cell2,used_cell3,used_cell4,tot_pairs,tot_triples,tot_quads)
#elif defined JACKKNIFE
    #pragma omp parallel firstprivate(seed_step,seed_shift,par,grid1,grid2,grid3,grid4,cf12,cf13,cf24) shared(sumint,outint,completed_loops,used_triples_per_sample,used_quads_per_sample,TotalTime,LoopTimes,gsl_rng_default,rd13,rd24,JK12,JK23,JK34,product_weights12_12,product_weights12_23,product_weights12_34) reduction(+:cell_attempt2,cell_attempt3,cell_attempt4,used_cell2,used_cell3,used_cell4,tot_pairs,tot_triples,tot_quads)
#else
    #pragma omp parallel firstprivate(seed_step,seed_shift,par,grid1,grid2,grid3,grid4,cf12,cf13,cf24) shared(sumint,outint,completed_loops,used_triples_per_sample,used_quads_per_sample,TotalTime,LoopTimes,gsl_rng_default,rd13,rd24,JK12,JK23,JK34) reduction(+:cell_attempt2,cell_attempt3,cell_attempt4,used_cell2,used_cell3,used_cell4,tot_pairs,tot_triples,tot_quads)
#endif
            { // start parallel loop
            // Decide which thread we are in
            int thread = omp_get_thread_num();
            assert(omp_get_num_threads()<=par->nthread);
            if (thread==0) printf("# Starting integral computation %d of %d on %d threads.\n", iter_no, tot_iter, omp_get_num_threads());
#else
            int thread = 0;
            printf("# Starting integral computation %d of %d single threaded.\n",iter_no,tot_iter);
            { // start loop
#endif

    //-----------DEFINE LOCAL THREAD VARIABLES
            Particle *prim_list; // list of particles in first cell
            int pln,sln,tln,fln,sln1,sln2; // number of particles in each cell
            int pid_j, pid_k, pid_l; // particle IDs particles drawn from j,k,l cell
            Particle particle_j, particle_k, particle_l; // randomly drawn particle
            //int* bin; // a-b bins for particles
            int* prim_ids; // list of particle IDs in primary cell
            double p2,p3,p4; // probabilities
#if (!defined LEGENDRE && !defined POWER)
            Float p22,p21; // probabilities
#endif
            int *bin_ij; // i-j separation bin
            int mnp = grid1->maxnp; // max number of particles in a grid1 cell
            Float *xi_ik, *w_ijk, *w_ij; // arrays to store xi and weight values
#ifdef PRINTPERCENTS
            Float percent_counter;
#endif
            int x, prim_id_1D;
            integer3 delta2, delta3, delta4, prim_id, sec_id, thi_id;
            Float3 cell_sep2,cell_sep3;

#ifdef LEGENDRE
            Float *poly_ij, *factor_ij;
            Integrals locint(par, cf12, cf13, cf24, I1, I2, I3, I4, survey_corr_12,survey_corr_23,survey_corr_34); // Accumulates the integral contribution of each thread
            // Assign memory for useful quantities
            int ecL=0;
            ecL+=posix_memalign((void **) &factor_ij, PAGE, sizeof(Float)*mnp);
            ecL+=posix_memalign((void **) &poly_ij, PAGE, sizeof(Float)*mnp*mbin);
            assert(ecL==0);
#elif defined POWER
            Float *poly_ij;
            Integrals locint(par, cf12, cf13, cf24, I1, I2, I3, I4, survey_corr_12, survey_corr_23, survey_corr_34);
            // Assign memory
            int ecL=0;
            ecL+=posix_memalign((void **) &poly_ij, PAGE, sizeof(Float)*mnp*mbin);
            assert(ecL==0);
#elif defined JACKKNIFE
            Integrals locint(par, cf12, cf13, cf24, JK12, JK23, JK34, I1, I2, I3, I4, product_weights12_12, product_weights12_23, product_weights12_34); // Accumulates the integral contribution of each thread
#else
            Integrals locint(par, cf12, cf13, cf24, JK12, JK23, JK34, I1, I2, I3, I4); // Accumulates the integral contribution of each thread
#endif
            gsl_rng* locrng = gsl_rng_alloc(gsl_rng_default); // one rng per thread
            gsl_rng_set(locrng, seed_step * (unsigned long)thread + seed_shift); // the second number, seed, will not overflow and will be unique for each thread in one run. Here, seed_shift is a random number between 0 and seed_step-1, inclusively. Seed clashes between different runs seem less likely than for the old formula, seed = steps * (nthread+1), with steps being a random number between 1 and UINT_MAX (or ULONG_MAX / nthread) inclusively - there, steps of one run could be a multiple of steps of the other resulting in the same seeds for some threads, which is somewhat more likely than getting the same random number.

            // Assign memory for intermediate steps
            int ec=0;
            ec+=posix_memalign((void **) &prim_list, PAGE, sizeof(Particle)*mnp);
            ec+=posix_memalign((void **) &prim_ids, PAGE, sizeof(int)*mnp);
            ec+=posix_memalign((void **) &bin_ij, PAGE, sizeof(int)*mnp);
            ec+=posix_memalign((void **) &w_ij, PAGE, sizeof(Float)*mnp);
            ec+=posix_memalign((void **) &xi_ik, PAGE, sizeof(Float)*mnp);
            ec+=posix_memalign((void **) &w_ijk, PAGE, sizeof(Float)*mnp);
            assert(ec==0);

            uint64 loc_used_pairs,loc_used_triples, loc_used_quads; // local counts of used pairs/triples/quads
    //-----------START FIRST LOOP-----------
    #ifdef OPENMP
    #pragma omp for schedule(dynamic)
    #endif
            for (int n_loops = 0; n_loops<par->max_loops; n_loops++){
#ifdef PRINTPERCENTS
                percent_counter=0.;
#endif
                loc_used_pairs=0; loc_used_triples=0; loc_used_quads=0;
                LoopTimes[n_loops].Start();

                // LOOP OVER ALL FILLED I CELLS
                for (int n1=0; n1<grid1->nf;n1++){

#ifdef PRINTPERCENTS
                    // Print time left
                    if((float(n1)/float(grid1->nf)*100)>=percent_counter){
                        printf("Integral %d of %d, iteration %d of %d on thread %d: Using cell %d of %d - %.0f percent complete\n", iter_no, tot_iter, 1+n_loops, par->max_loops, thread, n1+1, grid1->nf, percent_counter);
                        percent_counter+=5.;
                    }
#endif

                    // Pick first particle
                    prim_id_1D = grid1-> filled[n1]; // 1d ID for cell i
                    prim_id = grid1->cell_id_from_1d(prim_id_1D); // define first cell
                    pln = particle_list(prim_id_1D, prim_list, prim_ids, grid1); // update list of particles and number of particles

                    if(pln==0) continue; // skip if empty

                    loc_used_pairs+=pln*par->N2;
                    loc_used_triples+=pln*par->N2*par->N3;
                    loc_used_quads+=pln*par->N2*par->N3*par->N4;

                    // LOOP OVER N2 J CELLS
                    for (int n2=0; n2<par->N2; n2++){
                        cell_attempt2+=1; // new cell attempted

                        // Draw second cell from i weighted by 1/r^2
                        delta2 = rd13->random_cubedraw(locrng, &p2); // can use any rd class here since drawing as 1/r^2
                        // p2 is the ratio of sampling to true pair distribution here
                        sec_id = prim_id + delta2;
                        cell_sep2 = grid2->cell_sep(delta2);
                        x = draw_particle(sec_id, particle_j, pid_j, cell_sep2, grid2, sln, locrng, sln1, sln2);
                        if (x==1) continue; // skip failed draws

                        used_cell2+=1; // new cell accepted

#if (!defined LEGENDRE && !defined POWER)
                        // Probabilities for two random-particle partitions
                        p21=p2/(grid1->np1*(double)sln1); // divide probability by total number of particles in 1st partition and number in cell
                        p22=p2/(grid1->np2*(double)sln2); // for partition 2
#endif
                        // For all particles
                        p2*=1./(grid1->np*(double)sln); // probability is divided by total number of i particles and number of particles in cell

                        // Compute C2 integral
#ifdef LEGENDRE
                        locint.second(prim_list, prim_ids, pln, particle_j, pid_j, bin_ij, w_ij, p2, factor_ij, poly_ij);
#elif defined POWER
                        locint.second(prim_list, prim_ids, pln, particle_j, pid_j, bin_ij, w_ij, p2, poly_ij);
#else
                        locint.second(prim_list, prim_ids, pln, particle_j, pid_j, bin_ij, w_ij, p2, p21, p22);
#endif

                        // LOOP OVER N3 K CELLS
                        for (int n3=0; n3<par->N3; n3++){
                            cell_attempt3+=1; // new third cell attempted

                            // Draw third cell from i weighted by xi(r)
                            delta3 = rd13->random_xidraw(locrng, &p3); // use 1-3 random draw class here for xi_13
                            thi_id = prim_id + delta3;
                            cell_sep3 = grid3->cell_sep(delta3);
                            x = draw_particle_without_class(thi_id,particle_k,pid_k,cell_sep3,grid3,tln,locrng); // draw from third grid
                            if (x==1) continue; // skip failed draws
                            if ((pid_j==pid_k) && (I2==I3)) continue; // skip jk self-counts

                            used_cell3+=1; // new third cell used

                            p3*=p2/(double)tln; // update probability

                            // Compute third integral
#ifdef LEGENDRE
                            locint.third(prim_list, prim_ids, pln, particle_j, particle_k, pid_j, pid_k, bin_ij, w_ij, xi_ik, w_ijk, p3,factor_ij,poly_ij);
#elif defined POWER
                            locint.third(prim_list, prim_ids, pln, particle_j, particle_k, pid_j, pid_k, bin_ij, w_ij, xi_ik, w_ijk, p3, poly_ij);
#else
                            locint.third(prim_list, prim_ids, pln, particle_j, particle_k, pid_j, pid_k, bin_ij, w_ij, xi_ik, w_ijk, p3);
#endif

                            // LOOP OVER N4 L CELLS
                            for (int n4=0; n4<par->N4; n4++){
                                cell_attempt4+=1; // new fourth cell attempted

                                // Draw fourth cell from j cell weighted by xi_24(r)
                                delta4 = rd24->random_xidraw(locrng,&p4);
                                x = draw_particle_without_class(sec_id+delta4,particle_l,pid_l,cell_sep2+grid4->cell_sep(delta4),grid4,fln,locrng); // draw from 4th grid
                                if (x==1) continue; // skip failed draws
                                if (((pid_l==pid_j) && (I4==I2)) || ((pid_l==pid_k) && (I4==I3))) continue; // skip jl and kl self-counts

                                used_cell4+=1; // new fourth cell used

                                p4*=p3/(double)fln;


                                // Now compute the four-point integral
#ifdef LEGENDRE
                                locint.fourth(prim_list, prim_ids, pln, particle_j, particle_k, particle_l, pid_j, pid_k, pid_l, bin_ij, w_ijk, xi_ik, p4, factor_ij, poly_ij);
#elif defined POWER
                                locint.fourth(prim_list, prim_ids, pln, particle_j, particle_k, particle_l, pid_j, pid_k, pid_l, bin_ij, w_ijk, xi_ik, p4, poly_ij);
#else
                                locint.fourth(prim_list, prim_ids, pln, particle_j, particle_k, particle_l, pid_j, pid_k, pid_l, bin_ij, w_ijk, xi_ik, p4);
#endif

                            }
                        }
                    }
                }

                // Update used pair/triple/quad counts
                tot_pairs+=loc_used_pairs;
                tot_triples+=loc_used_triples;
                tot_quads+=loc_used_quads;

    #ifdef OPENMP
    #pragma omp critical // only one processor can access at once
    #endif
            {
                printf("Integral %d of %d, iteration %d of %d on thread %d completed\n", iter_no, tot_iter, 1+n_loops, par->max_loops, thread);
                int subsample_index = completed_loops / par->loops_per_sample; // index of output subsample for this loop
                completed_loops++; // increment completed loops counter, since they may be done not according to n_loops order
                if (completed_loops % par->nthread == 0) { // Print every nthread completed loops
                    TotalTime.Stop(); // interrupt timing to access .Elapsed()
                    int current_runtime = TotalTime.Elapsed();
                    int remaining_time = current_runtime*(par->max_loops - completed_loops)/completed_loops;  // estimated remaining time
                    fprintf(stderr, "\nFinished %d integral loops of %d after %d s. Estimated time left:  %2.2d:%2.2d:%2.2d hms, i.e. %d s.\n", completed_loops, par->max_loops, current_runtime, remaining_time/3600, remaining_time/60%60, remaining_time%60, remaining_time);

                    TotalTime.Start(); // Restart the timer
                    Float frob_C2, frob_C3, frob_C4;
#ifndef JACKKNIFE
                    sumint.frobenius_difference_sum(&locint, subsample_index * par->loops_per_sample, frob_C2, frob_C3, frob_C4); // since sumint is only incremented every loops_per_sample iterations, subsample_index * loops_per_sample is exactly how many loops are stored in the sumint at the moment. Thus if loops_per_sample>1 the Frobenius percent difference may be an overestimate.
                    fprintf(stderr, "Frobenius percent difference after %d loops is %.3f (C2), %.3f (C3), %.3f (C4)\n", completed_loops, frob_C2, frob_C3, frob_C4);
#else
                    Float frob_C2j, frob_C3j, frob_C4j;
                    sumint.frobenius_difference_sum(&locint, subsample_index * par->loops_per_sample, frob_C2, frob_C3, frob_C4, frob_C2j, frob_C3j, frob_C4j); // since sumint is only incremented every loops_per_sample iterations, subsample_index * loops_per_sample is exactly how many loops are stored in the sumint at the moment. Thus if loops_per_sample>1 the Frobenius percent difference may be an overestimate.
                    fprintf(stderr, "Frobenius percent difference after %d loops is %.3f (C2), %.3f (C3), %.3f (C4)\n", completed_loops, frob_C2, frob_C3, frob_C4);
                    fprintf(stderr, "Frobenius jackknife percent difference after %d loops is %.3f (C2j), %.3f (C3j), %.3f (C4j)\n", completed_loops, frob_C2j, frob_C3j, frob_C4j);
#endif
                }

                // Sum up integrals
                outint.sum_ints(&locint); // subsample total

                // Sum up pairs, triples, quads for the group
                used_pairs_per_sample += loc_used_pairs;
                used_triples_per_sample += loc_used_triples;
                used_quads_per_sample += loc_used_quads;

                // Save output if the group is done
                if (completed_loops % par->loops_per_sample == 0) {
                    sumint.sum_ints(&outint); // add to grand total
                    std::string output_string = string_format("%d", subsample_index);
#ifndef POWER
                    outint.normalize(grid1->norm, grid2->norm, grid3->norm, grid4->norm, (Float)used_pairs_per_sample, (Float)used_triples_per_sample, (Float)used_quads_per_sample);
#else
                    outint.normalize(grid1->norm, grid2->norm, grid3->norm, grid4->norm, (Float)used_pairs_per_sample, (Float)used_triples_per_sample, (Float)used_quads_per_sample, par->power_norm);
#endif
                    outint.save_integrals(output_string.c_str(), 0);
#ifdef JACKKNIFE
                    outint.save_jackknife_integrals(output_string.c_str());
#endif
                    // Reset the current output sample variables
                    outint.reset();
                    used_pairs_per_sample = used_triples_per_sample = used_quads_per_sample = 0;
                }
                locint.sum_total_counts(cnt2, cnt3, cnt4);
                locint.reset();

                }

            LoopTimes[n_loops].Stop();
            } // end cycle loop

            // Free up allocated memory at end of process
            free(prim_list);
            free(xi_ik);
            free(bin_ij);
            free(w_ij);
            free(w_ijk);
    } // end OPENMP loop

    //-----------REPORT + SAVE OUTPUT---------------
        TotalTime.Stop();

        // Normalize the accumulated results, using the RR counts
#ifndef POWER
        sumint.normalize(grid1->norm,grid2->norm,grid3->norm,grid4->norm,(Float)tot_pairs, (Float)tot_triples,(Float)tot_quads);
#else
        sumint.normalize(grid1->norm,grid2->norm,grid3->norm,grid4->norm,(Float)tot_pairs, (Float)tot_triples,(Float)tot_quads, par->power_norm);
#endif

        int runtime = TotalTime.Elapsed();
        printf("\n\nINTEGRAL %d OF %d COMPLETE\n",iter_no,tot_iter);
        for (int n_loop = 0; n_loop < par->max_loops; n_loop++) {
            int loop_runtime = LoopTimes[n_loop].Elapsed();
            fprintf(stderr, "Loop %d time: %d s, i.e. %2.2d:%2.2d:%2.2d hms\n", n_loop, loop_runtime, loop_runtime/3600, loop_runtime/60%60, loop_runtime%60);
        }
        fprintf(stderr, "\nTotal process time for %.2e sets of cells and %.2e quads of particles: %d s, i.e. %2.2d:%2.2d:%2.2d hms\n", double(used_cell4),double(tot_quads),runtime, runtime/3600,runtime/60%60,runtime%60);
        printf("We tried %.2e pairs, %.2e triples and %.2e quads of cells.\n",double(cell_attempt2),double(cell_attempt3),double(cell_attempt4));
        printf("Of these, we accepted %.2e pairs, %.2e triples and %.2e quads of cells.\n",double(used_cell2),double(used_cell3),double(used_cell4));
        printf("We sampled %.2e pairs, %.2e triples and %.2e quads of particles.\n",double(tot_pairs),double(tot_triples),double(tot_quads));
        double real_cnt2 = cnt2, real_cnt3 = cnt3, real_cnt4 = cnt4;
#if (defined LEGENDRE || LEGENDRE_MIX || POWER)
        int n_l = par->max_l / 2 + 1; // general expression for the number of multipoles. This is equal to mbin for LEGENDRE and POWER but not LEGENDRE_MIX
        double cnt_norm = pow(n_l, 2); // in these modes, one pair contributes to all multipole "bins" so the real counts are n_l^2 times smaller
        real_cnt2 /= cnt_norm; real_cnt3 /= cnt_norm; real_cnt4 /= cnt_norm;
#endif
        printf("Of these, we have integral contributions from %.2e pairs, %.2e triples and %.2e quads of particles.\n", real_cnt2, real_cnt3, real_cnt4);
        printf("Cell acceptance ratios are %.3f for pairs, %.3f for triples and %.3f for quads.\n", (double)used_cell2/cell_attempt2, (double)used_cell3/cell_attempt3, (double)used_cell4/cell_attempt4);
        printf("Acceptance ratios are %.3f for pairs, %.3f for triples and %.3f for quads.\n", real_cnt2/tot_pairs, real_cnt3/tot_triples, real_cnt4/tot_quads);
        printf("Average of %.2f pairs accepted per primary particle.\n\n",(Float)cnt2/grid1->np);

        printf("\nTrial speed: %.2e quads per core per second\n",double(tot_quads)/(runtime*double(par->nthread)));
        printf("Acceptance speed: %.2e quads per core per second\n",double(cnt4)/(runtime*double(par->nthread)));

        const char out_string[5] = "full";
        sumint.save_integrals(out_string, 1); // save integrals to file
        sumint.save_counts(tot_pairs,tot_triples,tot_quads); // save total pair/triple/quads attempted to file
#ifdef POWER
        printf("Printed integrals to file in the %sPowerCovMatrices/ directory\n",par->out_file);
#elif defined THREE_PCF
        printf("Printed integrals to file in the %s3PCFCovMatricesAll/ directory\n",par->out_file);
#elif defined JACKKNIFE
        sumint.save_jackknife_integrals(out_string);
        printf("Printed jackknife integrals to file in the %sCovMatricesJack/ directory\n",par->out_file);
#endif
        fflush(NULL);
        return;
        }

    };

#endif
