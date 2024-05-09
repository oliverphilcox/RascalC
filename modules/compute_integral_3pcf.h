 // This class computes the contributions to the threePCF autocovariance C_ab integral iterating over cells and particles. 

    #ifndef COMPUTE_INTEGRAL_3PCF_H
    #define COMPUTE_INTEGRAL_3PCF_H
    
// Load class to hold integrals
#include "integrals_3pcf.h"

class compute_integral{
        
    private:
        uint64 cnt3=0,cnt4=0,cnt5=0,cnt6=0;
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
        
        inline int draw_particle(const integer3 id_3D, Particle& particle, int& pid, Float3 shift, Grid *grid, int& n_particles, gsl_rng* locrng){
            // Draw a random particle from a cell given the cell ID.
            // This updates the particle and particle ID and returns 1 if error.
            
            int id_1D = grid->test_cell(id_3D); 
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
        compute_integral(){};
        
        compute_integral(Grid *grid, Parameters *par, CorrelationFunction *cf, RandomDraws *rd, SurveyCorrection *survey_corr, int iter_no){

            // MAIN FUNCTION TO COMPUTE INTEGRALS
            
            int tot_iter=2; // total number of integral sets to compute
            
            nbin = par->nbin; // number of radial bins
            mbin = par->mbin; // number of Legendre bins
            
            STimer initial, TotalTime; // Time initialization
            initial.Start(); 
            
            int convergence_counter=0, printtime=0;// counter to stop loop early if convergence is reached.
              
    //-----------INITIALIZE OPENMP + CLASSES----------
            unsigned long seed_step = (unsigned long)(std::numeric_limits<uint32_t>::max()) / (unsigned long)(par->max_loops); // restrict the seeds to 32 bits to avoid possible collisions, as most GSL random generators only accept 32-bit seeds and truncate the higher bits. unsigned long is at least 32 bits but can be more.
            unsigned long seed_shift = par->seed % seed_step; // for each loop, it will be seed = seed_step * n_loops + seed_shift, where n_loops goes from 0 to max_loops-1 => no overflow and guaranteed unique for each thread
            if (par->random_seed) {
                std::random_device urandom("/dev/urandom");
                std::uniform_int_distribution<unsigned long> dist(0, seed_step-1); // distribution of integers with these limits, inclusive
                seed_shift = dist(urandom); // for each thread, it will be seed = seed_step * n_loops + seed_shift, where n_loops goes from 0 to max_loops-1 => no overflow and guaranteed unique for each thread
            }

            gsl_rng_env_setup(); // initialize gsl rng
            Integrals sumint(par, cf, survey_corr); // total integral

            uint64 tot_triples=0, tot_quads=0, tot_quints=0, tot_hexes=0; // global number of particle sets used (including those rejected for being in the wrong bins)
            uint64 cell_attempt3=0,cell_attempt4=0,cell_attempt5=0,cell_attempt6=0; // number of k,l,m,n cells attempted
            uint64 used_cell3=0,used_cell4=0,used_cell5=0,used_cell6=0; // number of used j,k,l cells
        
            initial.Stop();
            fprintf(stderr, "Init time: %g s\n",initial.Elapsed());
            printf("# 1st grid filled cells: %d\n",grid->nf);
            printf("# All 1st grid points in use: %d\n",grid->np);
            printf("# Max points in one cell in grid 1%d\n",grid->maxnp);
            fflush(NULL);
            
            TotalTime.Start(); // Start timer
            
#ifdef OPENMP       
    #pragma omp parallel firstprivate(seed_step,seed_shift,par,printtime,grid,cf) shared(sumint,TotalTime,gsl_rng_default,rd) reduction(+:convergence_counter,cell_attempt3,cell_attempt4,cell_attempt5,cell_attempt6,used_cell3,used_cell4,used_cell5,used_cell6,tot_triples,tot_quads,tot_quints,tot_hexes)
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
            int pln,sln,tln,fln,filn,siln; // number of particles in each cell
            int pid_j, pid_k, pid_l,pid_m,pid_n; // particle IDs particles drawn from cells
            Particle particle_j, particle_k, particle_l,particle_m,particle_n; // randomly drawn particle
            int *prim_ids; // list of particle IDs in primary cell
            double p2,p3,p4,p5,p6; // probabilities
            int mnp = grid->maxnp; // max number of particles in a grid cell
            Float *correction_ijk, *legendre_ijk, *xi_pass, *xi_pass2, xi_pass3=0, *w_ijk, *w_ijkl, *w_ijklm;
            Float norm_jk=0, norm_kl=0, norm_lm=0; // arrays to store xi and weight values
            int *bins_ijk, bin_jk=0, bin_kl=0, bin_lm=0; // i-j separation bin
            Float percent_counter;
            int x, prim_id_1D;
            integer3 delta2, delta3, delta4, delta5, delta6, prim_id, sec_id, thi_id, fou_id, fif_id, six_id;
            Float3 cell_sep2,cell_sep3,cell_sep4,cell_sep5,cell_sep6;
            
            Integrals locint(par, cf, survey_corr); // Accumulates the integral contribution of each thread
            
            gsl_rng* locrng = gsl_rng_alloc(gsl_rng_default); // one rng per thread, seed will be set later inside the loop
            
            // Assign memory for intermediate steps
            int ec=0;
            ec+=posix_memalign((void **) &prim_list, PAGE, sizeof(Particle)*mnp);
            ec+=posix_memalign((void **) &prim_ids, PAGE, sizeof(int)*mnp);
            ec+=posix_memalign((void **) &w_ijk, PAGE, sizeof(Float)*mnp);
            ec+=posix_memalign((void **) &w_ijkl, PAGE, sizeof(Float)*mnp);
            ec+=posix_memalign((void **) &w_ijklm, PAGE, sizeof(Float)*mnp);
            ec+=posix_memalign((void **) &bins_ijk, PAGE, sizeof(int)*mnp*6);
            ec+=posix_memalign((void **) &correction_ijk, PAGE, sizeof(Float)*mnp*3);
            ec+=posix_memalign((void **) &legendre_ijk, PAGE, sizeof(Float)*mbin*3*mnp);
            ec+=posix_memalign((void **) &xi_pass, PAGE, sizeof(Float)*mnp);
            ec+=posix_memalign((void **) &xi_pass2, PAGE, sizeof(Float)*mnp);
            assert(ec==0);
            
            uint64 loc_used_triples, loc_used_quads, loc_used_quints, loc_used_hexes; // local counts of used triples->hexes
            
    //-----------START FIRST LOOP-----------
    #ifdef OPENMP
    #pragma omp for schedule(dynamic)
    #endif
            for (int n_loops = 0; n_loops<par->max_loops; n_loops++){
                percent_counter=0.;
                loc_used_triples=0; loc_used_quads=0; loc_used_quints=0; loc_used_hexes=0;  
                
                // End loops early if convergence has been acheived
                if (convergence_counter==10){ 
                    if (printtime==0) printf("1 percent convergence acheived in C6 10 times, exiting.\n");
                    printtime++;
                    continue;
                    }
                
                // Set/reset the RNG seed based on loop number instead of thread number to reproduce results with different number of threads but other parameters kept the same. Individual subsamples may differ because they are accumulated/written in order of loop completion which may depend on external factors at runtime, but the final integrals should be the same.
                gsl_rng_set(locrng, seed_step * (unsigned long)n_loops + seed_shift); // the second number, seed, will not overflow and will be unique for each loop in one run. Here, seed_shift is a random number between 0 and seed_step-1, inclusively.
                // Seed clashes between different runs seem less likely than for the old formula, seed = steps * (nthread+1), with steps being a random number between 1 and UINT_MAX (or ULONG_MAX / nthread) inclusively - there, steps of one run could be a multiple of steps of the other resulting in the same seeds for some threads, which is somewhat more likely than getting the same random number.
                    
                // LOOP OVER ALL FILLED I CELLS
                for (int n1=0; n1<grid->nf;n1++){
                    
                    // Print time left
                    if((float(n1)/float(grid->nf)*100)>=percent_counter){
                        printf("Integral %d of %d, run %d of %d on thread %d: Using cell %d of %d - %.0f percent complete\n",iter_no+1,tot_iter,1+n_loops/par->nthread, int(ceil(float(par->max_loops)/(float)par->nthread)),thread, n1+1,grid->nf,percent_counter);
                        percent_counter+=5.;
                    }
                    
                    // Pick first particle
                    prim_id_1D = grid->filled[n1]; // 1d ID for cell i 
                    prim_id = grid->cell_id_from_1d(prim_id_1D); // define first cell
                    pln = particle_list(prim_id_1D, prim_list, prim_ids, grid); // update list of particles and number of particles
                    
                    if(pln==0) continue; // skip if empty
                    
                    loc_used_triples+=pln*par->N2*par->N3;
                    loc_used_quads+=pln*par->N2*par->N3*par->N4;
                    loc_used_quints+=pln*par->N2*par->N3*par->N4*par->N5;
                    loc_used_hexes+=pln*par->N2*par->N3*par->N4*par->N5*par->N6;
                    
                    // LOOP OVER N2 J CELLS
                    for (int n2=0; n2<par->N2; n2++){
                        
                        // Draw second cell from i                         
                        if(iter_no==0) delta2 = rd->random_xidraw(locrng, &p2); // weight by xi(r)
                        else delta2 = rd->random_cubedraw(locrng,&p2); // weight by 1/r^2
                        
                        sec_id = prim_id + delta2;
                        cell_sep2 = grid->cell_sep(delta2);
                        
                        x = draw_particle(sec_id, particle_j, pid_j,cell_sep2, grid, sln, locrng);
                        if(x==1) continue; // skip if error
                        
                        // For all particles
                        p2*=1./(grid->np*(double)sln); // probability is divided by total number of i particles and number of particles in cell
                        
                        // LOOP OVER N3 K CELLS
                        for (int n3=0; n3<par->N3; n3++){
                            cell_attempt3+=1; // new third cell attempted
                            
                            // Draw third cell from second cell
                            if(iter_no==0) delta3 = rd->random_cubedraw(locrng,&p3); // weight by 1/r^2
                            else delta3 = rd->random_xidraw(locrng,&p3); // xi(r) weighting
                            
                            thi_id = sec_id + delta3;
                            cell_sep3 = cell_sep2 + grid->cell_sep(delta3); 
                            x = draw_particle(thi_id,particle_k,pid_k,cell_sep3,grid,tln,locrng); 
                            if(x==1) continue; 
                            if(pid_k==pid_j) continue; // no self-counting
                            
                            used_cell3+=1; // new third cell used
                            
                            p3*=p2/(double)tln; // update probability
                            
                            // Compute third integral
                            locint.third(prim_list, prim_ids, pln, particle_j, particle_k, pid_j, pid_k, p3, w_ijk, bins_ijk, correction_ijk, legendre_ijk,xi_pass,norm_jk, bin_jk,iter_no);
                            
                            // LOOP OVER N4 L CELLS
                            for (int n4=0; n4<par->N4; n4++){
                                cell_attempt4+=1; // new fourth cell attempted
                                
                                // Draw fourth cell from k cell weighted by xi(r)
                                delta4 = rd->random_xidraw(locrng,&p4);
                                
                                if(iter_no==0){
                                    fou_id = thi_id + delta4;
                                    cell_sep4 = cell_sep3 + grid->cell_sep(delta4);
                                }
                                else{
                                    fou_id = prim_id + delta4;
                                    cell_sep4 = grid->cell_sep(delta4);
                                }
                                
                                x = draw_particle(fou_id,particle_l,pid_l,cell_sep4,grid,fln,locrng); // draw from grid
                                if(x==1) continue;
                                if((pid_l==pid_k)||(pid_l==pid_j)) continue;
                                
                                used_cell4+=1; // new fourth cell used
                                
                                p4*=p3/(double)fln;
                                
                                // Now compute the four-point integral
                                locint.fourth(prim_list, prim_ids, pln, particle_j, particle_k, particle_l, pid_j, pid_k, pid_l, p4, w_ijk, bins_ijk, correction_ijk, legendre_ijk, xi_pass, norm_jk, bin_jk, w_ijkl, xi_pass2,  norm_kl, bin_kl, iter_no);
                                
                                 // LOOP OVER N5 M CELLS
                                for (int n5=0; n5<par->N5; n5++){
                                    cell_attempt5+=1; // new fourth cell attempted
                                    
                                    // Draw fifth cell from l or i cell
                                    delta5 = rd->random_xidraw(locrng,&p5); 
                                    if(iter_no==0){
                                        fif_id = fou_id + delta5;
                                        cell_sep5 = cell_sep4 + grid->cell_sep(delta5);
                                    }
                                    else{
                                        fif_id = sec_id + delta5;
                                        cell_sep5 = cell_sep2 + grid->cell_sep(delta5);
                                    }
                                        
                                    x = draw_particle(fif_id,particle_m,pid_m,cell_sep5,grid,filn,locrng); // draw from grid
                                    if(x==1) continue;
                                    if((pid_m==pid_l)||(pid_m==pid_k)||(pid_m==pid_j)) continue;
                                    
                                    used_cell5+=1; // new fourth cell used
                                    
                                    p5*=p4/(double)filn;
                                    
                                    // Now compute the four-point integral
                                    locint.fifth(prim_list, prim_ids, pln, particle_j, particle_k, particle_l, particle_m, pid_j, pid_k, pid_l, pid_m, p5, w_ijkl, bins_ijk, correction_ijk, legendre_ijk, xi_pass, xi_pass2, norm_kl, bin_kl, w_ijklm, xi_pass3, norm_lm, bin_lm, iter_no); 
                                    
                                    // LOOP OVER N6 N CELLS
                                    for (int n6=0; n6<par->N6; n6++){
                                        cell_attempt6+=1; // new fourth cell attempted
                                        
                                        // Draw sixth cell from m cell weighted by xi(r)
                                        delta6 = rd->random_xidraw(locrng,&p6); 
                                        if (iter_no==0){
                                            six_id = fif_id + delta6;
                                            cell_sep6 = cell_sep5 + grid->cell_sep(delta6);
                                        }
                                        else{
                                            six_id = thi_id + delta6;
                                            cell_sep6 = cell_sep3 + grid->cell_sep(delta6);
                                        }
                                        
                                        x = draw_particle(six_id,particle_n,pid_n,cell_sep6,grid,siln,locrng); // draw from grid
                                        if(x==1) continue;
                                        if((pid_n==pid_m)||(pid_n==pid_l)||(pid_n==pid_k)||(pid_n==pid_j)) continue;
                                        
                                        used_cell6+=1; // new sixth cell used
                                        
                                        p6*=p5/(double)siln;
                                        
                                        // Now compute the four-point integral
                                        locint.sixth(prim_list, prim_ids, pln, particle_j, particle_k, particle_l, particle_m, particle_n, pid_j, pid_k, pid_l, pid_m, pid_n, p6, w_ijklm, bins_ijk, correction_ijk, legendre_ijk, xi_pass, xi_pass2, xi_pass3, norm_lm, bin_lm, iter_no); 
                                    }
                                }
                            }
                        }
                    }
                }
                
                // Update used pair/triple/quad counts
                tot_triples+=loc_used_triples;
                tot_quads+=loc_used_quads; 
                tot_quints+=loc_used_quints;
                tot_hexes+=loc_used_hexes;
               
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
                    Float frob_C3, frob_C4, frob_C5, frob_C6;
                    sumint.frobenius_difference_sum(&locint,n_loops, frob_C3, frob_C4, frob_C5, frob_C6);
                    if(frob_C6<0.01) convergence_counter++;
                    if (n_loops!=0){
                        fprintf(stderr,"Frobenius percent difference after loop %d is %.3f (C3), %.3f (C4), %.3f (C5), %.3f (C6) \n",n_loops,frob_C3, frob_C4, frob_C5, frob_C6);
                    }
                }
                    
                // Sum up integrals
                sumint.sum_ints(&locint); 
                
                // Save output after each loop
                std::string output_string = string_format("%d", n_loops);
                
                locint.normalize(grid->norm,(Float)loc_used_triples, (Float)loc_used_quads, (Float)loc_used_quints, (Float)loc_used_hexes);
                locint.save_integrals(output_string.c_str(), 0, iter_no);
                locint.sum_total_counts(cnt3, cnt4, cnt5, cnt6); 
                locint.reset();
                
                }
            
                
            } // end cycle loop
            
            // Free up allocated memory at end of process
            free(prim_list);
            free(prim_ids);
            free(w_ijk);
            free(w_ijkl);
            free(w_ijklm);
            free(bins_ijk);
            free(correction_ijk);
            free(legendre_ijk);
            free(xi_pass);
            free(xi_pass2);
            
    } // end OPENMP loop

    //-----------REPORT + SAVE OUTPUT---------------
        TotalTime.Stop();
        
        // Normalize the accumulated results, using the RR counts
        sumint.normalize(grid->norm,(Float)tot_triples,(Float)tot_quads,(Float)tot_quints,(Float)tot_hexes);
        
        // Normalize counts by expected number in each cell
        cnt3/=(9.*mbin*mbin);
        cnt4/=(9.*mbin*mbin);
        cnt5/=(9.*mbin*mbin);
        cnt6/=(9.*mbin*mbin);
        
        int runtime = TotalTime.Elapsed();
        printf("\n\nINTEGRAL %d OF %d COMPLETE\n",iter_no+1,tot_iter); 
        fprintf(stderr, "\nTotal process time for %.2e sets of cells and %.2e hexes of particles: %d s, i.e. %2.2d:%2.2d:%2.2d hms\n", double(used_cell6),double(tot_hexes),runtime, runtime/3600,runtime/60%60,runtime%60);
        printf("We tried %.2e triples, %.2e quads, %.2e quints and %.2e hexes of cells.\n",double(cell_attempt3),double(cell_attempt4),double(cell_attempt5),double(cell_attempt6));
        printf("Of these, we accepted %.2e triples, %.2e quads, %.2e quints and %.2e hexes of cells.\n",double(used_cell3),double(used_cell4),double(used_cell5),double(used_cell6));
        printf("We sampled %.2e triples, %.2e quads, %.2e quints and %.2e hexes of particles.\n",double(tot_triples),double(tot_quads),double(tot_quints),double(tot_hexes));
        printf("Of these, we have integral contributions from %.2e triples, %.2e quads, %.2e quints and %.2e hexes of particles.\n",double(cnt3),double(cnt4),double(cnt5),double(cnt6));
        printf("Cell acceptance ratios are %.3f for triples, %.3f for quads, %.3f for quints and %.3f for hexes.\n",(double)used_cell3/cell_attempt3,(double)used_cell4/cell_attempt4,(double)used_cell5/cell_attempt5,(double)used_cell6/cell_attempt6);
        printf("Acceptance ratios are %.3f for triples, %.3f for quads, %.3f for quints and %.3f for hexes.\n",(double)cnt3/tot_triples,(double)cnt4/tot_quads,(double)cnt5/tot_quints,(double)cnt6/tot_hexes);
        printf("Average of %.2f triples accepted per primary particle.\n\n",(double)cnt3/grid->np);
        
        printf("\nTrial speed: %.2e hexes per core per second\n",double(tot_hexes)/(runtime*double(par->nthread)));
        printf("Acceptance speed: %.2e hexes per core per second\n",double(cnt6)/(runtime*double(par->nthread)));
        
        const char out_string[5] = "full";
        sumint.save_integrals(out_string,1,iter_no); // save integrals to file
        sumint.save_counts(tot_triples,tot_quads,tot_quints,tot_hexes,iter_no); // save total pair/triple/quads attempted to file
        printf("Printed integrals to file in the %s3PCFCovMatricesAll/ directory\n",par->out_file);
        fflush(NULL);
        return;
        }
                
    };
            
#endif 
