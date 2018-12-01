 // This class computes the contributions to the C_ab integral iterating over cells and particles. It is a heavily modified version of code by Alex Wiegand.

    #ifndef COMPUTE_INTEGRAL3_H
    #define COMPUTE_INTEGRAL3_H

    #include "integrals2.h"

    class compute_integral3{
        
    private:
        Grid *grid;
        //CorrelationFunction *cf;
        JK_weights *JK;
        //RandomDraws2 *rd;
        uint64 cnt2=0,cnt3=0,cnt4=0;
        int nbin, mbin; 
        
    private:
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
        
        int draw_particle(integer3 id_3D, Particle &particle, int &pid, integer3 shift, Grid *grid, int &n_particles, gsl_rng* locrng, int &n_particles1, int &n_particles2){
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
        
        void check_threads(Parameters *par){
            // Set up OPENMP and define which threads ot use
    #ifdef OPENMP
            cpu_set_t mask[par->nthread+1];
            int tnum=0;
            sched_getaffinity(0, sizeof(cpu_set_t), &mask[par->nthread]);
            fprintf(stderr, " CPUs used are: ");
            for(int ii=0;ii<64;ii++){
                if(CPU_ISSET(ii, &mask[par->nthread])){
                    fprintf(stderr,"%d ", ii);
                    CPU_ZERO(&mask[tnum]);
                    CPU_SET(ii,&mask[tnum]);
                    tnum++;
                }
            }
            fprintf(stderr,"\n");
    #endif
        }
        
    public:    
        compute_integral3(Grid *grid, Parameters *par, JK_weights *JK){
            // MAIN FUNCTION TO COMPUTE INTEGRALS
            
            nbin = par->nbin; // number of radial bins
            mbin = par->mbin; // number of mu bins
            
            STimer initial, TotalTime; // Time initialization
            initial.Start(); 
            
            int convergence_counter=0, printtime=0;// counter to stop loop early if convergence is reached.
            
    //-----------INITIALIZE OPENMP + CLASSES----------
            std::random_device urandom("/dev/urandom");
            std::uniform_int_distribution<unsigned int> dist(1, std::numeric_limits<unsigned int>::max());
            unsigned long int steps = dist(urandom);        
            gsl_rng_env_setup(); // initialize gsl rng
            CorrelationFunction *cf=new CorrelationFunction(par->corname,par->mbin,par->mumax-par->mumin, par->r_cutoff);
            RandomDraws2 *rd=new RandomDraws2(cf,par, NULL, 0);
            Integrals2 sumint(par,cf,JK); // total integral

            
            uint64 tot_pairs=0, tot_triples=0, tot_quads=0; // global number of particle pairs/triples/quads used (including those rejected for being in the wrong bins)
            uint64 cell_attempt2=0,cell_attempt3=0,cell_attempt4=0; // number of j,k,l cells attempted
            uint64 used_cell2=0,used_cell3=0,used_cell4=0; // number of used j,k,l cells
            
            check_threads(par); // Define which threads we use
        
            initial.Stop();
            fprintf(stderr, "Init time: %g s\n",initial.Elapsed());
            printf("# Filled cells: %d\n",grid->nf);
            printf("# All points in use: %d\n",grid->np);
            printf("# Max points in one cell %d\n",grid->maxnp);
            fflush(NULL);
            
            TotalTime.Start(); // Start timer
            
    #ifdef OPENMP       
    #pragma omp parallel firstprivate(steps,grid,par,printtime) shared(sumint,TotalTime,JK,cf,rd,gsl_rng_default) reduction(+:convergence_counter,cell_attempt2,cell_attempt3,cell_attempt4,used_cell2,used_cell3,used_cell4,tot_pairs,tot_triples,tot_quads)
            { // start parallel loop
            // Decide which thread we are in
            int thread = omp_get_thread_num();
            assert(omp_get_num_threads()<=par->nthread);
            if (thread==0) printf("# Starting integral computation on %d threads.\n", omp_get_num_threads());        
    #else
            int thread = 0;
            printf("# Starting integral computation single threaded.\n");
    #endif
            //TODO: Make sure this works single-threaded
            
    //-----------DEFINE VARIABLES
            Particle *prim_list; // list of particles in first cell
            int pln,sln,tln,fln,sln1,sln2; // number of particles in each cell
            int pid_j, pid_k, pid_l; // particle IDs particles drawn from j,k,l cell
            Particle particle_j, particle_k, particle_l; // randomly drawn particle
            //int* bin; // a-b bins for particles
            int* prim_ids; // list of particle IDs in primary cell
            double p2,p3,p4,p22,p21; // probabilities
            int mnp = grid->maxnp; // max number of particles
            Float *xi_ik, *w_ijk, *w_ij; // arrays to store xi and weight values
            int *bin_ij; // i-j separation bin
            Float percent_counter;
            int x, prim_id_1D;
            integer3 delta2, delta3, delta4, prim_id, sec_id, thi_id;
            Float3 cell_sep2,cell_sep3;
            
            Integrals2 locint(par,cf,JK); // Accumulates the integral contribution of each thread
            
            gsl_rng* locrng = gsl_rng_alloc(gsl_rng_default); // one rng per thread
            gsl_rng_set(locrng, steps*(thread+1));
            
            // Assign memory
            int ec=0;
            ec+=posix_memalign((void **) &prim_list, PAGE, sizeof(Particle)*mnp);
            ec+=posix_memalign((void **) &prim_ids, PAGE, sizeof(int)*mnp);
            ec+=posix_memalign((void **) &bin_ij, PAGE, sizeof(int)*mnp);
            ec+=posix_memalign((void **) &w_ij, PAGE, sizeof(Float)*mnp);
            ec+=posix_memalign((void **) &xi_ik, PAGE, sizeof(Float)*mnp);
            ec+=posix_memalign((void **) &w_ijk, PAGE, sizeof(Float)*mnp);
            assert(ec==0);
            
            uint64 loc_used_pairs, loc_used_triples, loc_used_quads; // local counts of used pairs/triples/quads
            
    //-----------START FIRST LOOP-----------
    #ifdef OPENMP
    #pragma omp for schedule(dynamic)
    #endif
            for (int n_loops = 0; n_loops<par->max_loops; n_loops++){
                percent_counter=0.;
                loc_used_pairs=0; loc_used_triples=0; loc_used_quads=0;  
                
                
                RandomDraws2 *rd=new RandomDraws2(cf,par, NULL, 0);
            
                
                // End loops early if convergence has been acheived
                if (convergence_counter==5){ 
                    if (printtime==0) printf("\n0.1 percent convergence acheived in C4 50 times, exiting.\n");
                    printtime++;
                    continue;
                    }
                    
                    
                // LOOP OVER ALL FILLED I CELLS
                for (int n1=0; n1<grid->nf;n1++){
                    
                    if((float(n1)/float(grid->nf)*100)>=percent_counter){
                        printf("Using cell %d of %d on core %d on run %d of %d: %.0f percent complete\n",n1+1,grid->nf,thread,1+n_loops/par->nthread, int(ceil(float(par->max_loops)/(float)par->nthread)),percent_counter);
                        percent_counter+=5.;
                    }
                    
                    prim_id_1D = grid-> filled[n1]; // 1d ID for cell i 
                    prim_id = grid->cell_id_from_1d(prim_id_1D); // define first cell
                    pln = particle_list(prim_id_1D, prim_list, prim_ids, grid); // update list of particles and number of particles
                    
                    if(pln==0) continue; // skip if empty
                    
                    loc_used_pairs+=pln*par->N2;
                    loc_used_triples+=pln*par->N2*par->N3;
                    loc_used_quads+=pln*par->N2*par->N3*par->N4;
                    
                    // LOOP OVER N2 J CELLS
                    for (int n2=0; n2<par->N2; n2++){
                        cell_attempt2+=1; // new cell attempted
                        
                        // Draw second cell from i weighted by 1/r^2
                        delta2 = rd->random_cubedraw(locrng, &p2); 
                        // p2 is the ratio of sampling to true pair distribution here
                        sec_id = prim_id + delta2;
                        cell_sep2 = grid->cell_sep(delta2);
                        x = draw_particle(sec_id, particle_j, pid_j,cell_sep2, grid, sln, locrng, sln1, sln2);
                        if(x==1) continue; // skip if error
                        
                        used_cell2+=1; // new cell accepted
                        
                        // Probabilities for two random-particle partitions
                        p21=p2/(grid->np1*(double)sln1); // divide probability by total number of particles in 1st partition and number in cell
                        p22=p2/(grid->np2*(double)sln2); // for partition 2
                        
                        // For all particles
                        p2*=1./(grid->np*(double)sln); // probability is divided by total number of particles and number of particles in cell
                        
                        // Compute C2 integral
                        locint.second(prim_list, prim_ids, pln, particle_j, pid_j, bin_ij, w_ij, p2, p21, p22);
                        
                        // LOOP OVER N3 K CELLS
                        for (int n3=0; n3<par->N3; n3++){
                            cell_attempt3+=1; // new third cell attempted
                            
                            // Draw third cell from i weighted by xi(r)
                            delta3 = rd->random_xidraw(locrng, &p3);
                            thi_id = prim_id + delta3;
                            cell_sep3 = grid->cell_sep(delta3);
                            x = draw_particle_without_class(thi_id,particle_k,pid_k,cell_sep3,grid,tln,locrng);
                            if(x==1) continue; 
                            
                            used_cell3+=1; // new third cell used
                            
                            p3*=p2/(double)tln; // update probability
                            
                            // Compute third integral
                            locint.third(prim_list, prim_ids, pln, particle_j, particle_k, pid_j, pid_k, bin_ij, w_ij, xi_ik, w_ijk, p3);
                            
                            // LOOP OVER N4 L CELLS
                            for (int n4=0; n4<par->N4; n4++){
                                cell_attempt4+=1; // new fourth cell attempted
                                
                                // Draw fourth cell from j cell weighted by xi(r)
                                delta4 = rd->random_xidraw(locrng,&p4);
                                x = draw_particle_without_class(sec_id+delta4,particle_l,pid_l,cell_sep2+grid->cell_sep(delta4),grid,fln,locrng);
                                if(x==1) continue;
                                
                                used_cell4+=1; // new fourth cell used
                                
                                p4*=p3/(double)fln;
                                
                                // Now compute the four-point integral
                                locint.fourth(prim_list, prim_ids, pln, particle_j, particle_k, particle_l, pid_j, pid_k, pid_l, bin_ij, w_ijk, xi_ik, p4);
                                
                            }
                        }
                    }
                }
                
    #pragma omp critical // only one processor can access at once
            {
                // Update used pair/triple/quad counts
                tot_pairs+=loc_used_pairs;
                tot_triples+=loc_used_triples;
                tot_quads+=loc_used_quads; 
                
                if (n_loops%par->nthread==0){ // Print every nthread loops
                    TotalTime.Stop(); // interrupt timing to access .Elapsed()
                    int current_runtime = TotalTime.Elapsed();
                    int remaining_time = current_runtime/(n_loops/par->nthread+1)*(par->max_loops/par->nthread-n_loops/par->nthread);  // estimated remaining time
                    //TODO: Fix this calculation
                    fprintf(stderr,"\nFinished integral loop %d of %d after %d s. Estimated time left:  %2.2d:%2.2d:%2.2d hms, i.e. %d s.\n",n_loops,par->max_loops, current_runtime,remaining_time/3600,remaining_time/60%60, remaining_time%60,remaining_time);
                    
                    TotalTime.Start(); // Restart the timer
                    Float frob_C2, frob_C3, frob_C4, frob_C2j, frob_C3j, frob_C4j;
                    sumint.frobenius_difference_sum(&locint,n_loops, frob_C2, frob_C3, frob_C4, frob_C2j, frob_C3j, frob_C4j);
                    if((frob_C4<0.01)&&(frob_C4j<0.01)) convergence_counter++;
                    
                    // Report N_eff for the integrals:
                    Float N_eff, N_eff_jack;
                    sumint.compute_Neff((Float)tot_quads, N_eff, N_eff_jack);
                    
                    if (n_loops!=0){
                        fprintf(stderr,"Frobenius percent difference after loop %d is %.3f (C2), %.3f (C3), %.3f (C4)\n",n_loops,frob_C2, frob_C3, frob_C4);
                        fprintf(stderr,"Frobenius jackknife percent difference after loop %d is %.3f (C2j), %.3f (C3j), %.3f (C4j)\n",n_loops,frob_C2j, frob_C3j, frob_C4j);
                        fprintf(stderr, "N_eff estimate: %.3e (C4) %.3e (C4j) \n", N_eff, N_eff_jack);
                    }
                }
                else{
                    sumint.sum_ints(&locint); 
                }
                // Save output after each loop
                char output_string[50];
                sprintf(output_string,"%d", n_loops);
                
                locint.normalize(grid->np,par->nofznorm, (Float)loc_used_pairs, (Float)loc_used_triples, (Float)loc_used_quads, 1); // do normalize by RR here
                locint.save_integrals(output_string,0);
                locint.save_jackknife_integrals(output_string);
                
                locint.sum_total_counts(cnt2, cnt3, cnt4); 
                locint.reset();
                
                }
            
                
            } // end cycle loop
            
            // Free up allocated memory at end of process
            free(prim_list);
            //gsl_rng_free(locrng);
            free(xi_ik);
            free(bin_ij);
            free(w_ij);
            free(w_ijk);
            
            
            
    } // end OPENMP loop

    //-----------REPORT + SAVE OUTPUT---------------
        TotalTime.Stop();
        
        // Normalize the accumulated results, using the RR counts
        sumint.normalize(grid->np, par->nofznorm, (Float)tot_pairs, (Float)tot_triples,(Float)tot_quads,1);//(Float)n_triples, (Float)n_quads, 1);
        
        int runtime = TotalTime.Elapsed();
        fprintf(stderr, "\nTotal process time for %.2e sets of cells and %.2e quads of particles: %d s, i.e. %2.2d:%2.2d:%2.2d hms\n", double(used_cell4),double(tot_quads),runtime, runtime/3600,runtime/60%60,runtime%60);
        printf("We tried %.2e pairs, %.2e triples and %.2e quads of cells.\n",double(cell_attempt2),double(cell_attempt3),double(cell_attempt4));
        printf("Of these, we accepted %.2e pairs, %.2e triples and %.2e quads of cells.\n",double(used_cell2),double(used_cell3),double(used_cell4));
        printf("We sampled %.2e pairs, %.2e triples and %.2e quads of particles.\n",double(tot_pairs),double(tot_triples),double(tot_quads));
        printf("Of these, we have integral contributions from %.2e pairs, %.2e triples and %.2e quads of particles.\n",double(cnt2),double(cnt3),double(cnt4));
        printf("Cell acceptance ratios are %.3f for pairs, %.3f for triples and %.3f for quads.\n",(double)used_cell2/cell_attempt2,(double)used_cell3/cell_attempt3,(double)used_cell4/cell_attempt4);
        printf("Acceptance ratios are %.3f for pairs, %.3f for triples and %.3f for quads.\n",(double)cnt2/tot_pairs,(double)cnt3/tot_triples,(double)cnt4/tot_quads);
        printf("Average of %.2f pairs accepted per primary particle.\n\n",(Float)cnt2/grid->np);
        
        printf("\nTrial speed: %.2e quads per core per second\n",double(tot_quads)/(runtime*double(par->nthread)));
        printf("Acceptance speed: %.2e quads per core per second\n",double(cnt4)/(runtime*double(par->nthread)));
        
        char out_string[5];
        sprintf(out_string,"full");
        sumint.save_integrals(out_string,1); // save integrals to file
        sumint.save_counts(tot_pairs,tot_triples,tot_quads); // save total pair/triple/quads attempted to file
        printf("Printed integrals to file in the %sCovMatricesAll/ directory\n",par->out_file);
        sumint.save_jackknife_integrals(out_string);
        printf("Printed jackknife integrals to file in the %sCovMatricesJack/ directory\n",par->out_file);
        fflush(NULL);
        return;
        }
        
        
    };
            
#endif 
