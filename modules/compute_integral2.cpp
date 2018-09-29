#ifndef COMPUTE_INTEGRAL2_H
#define COMPUTE_INTEGRAL2_H

class compute_integral2{

public:
    CorrelationFunction *cf;
    RandomDraws *rd;
    
    // initialize particles and list of cells
    Particle *prim_list; // NEED
    Particle *sec_list;
    Particle *thi_list;
    Particle *fou_list;
    Cell primary, sec, thi, fou; // NEED
    int pln, sln, tln, fln;
    double p1,p2,p3,p4;    
    Grid *grid;
    gsl_rng* locrng; 


private:
    integer3 prim_id, sec_id, thi_id, fou_id;
    integer3 delta2, delta3, delta4; // DON'T NEED?
    int* bin;
    Float* xi;    
    int mnp;
    
    
public:
    
int particle_list(integer3 id_3D, Cell &cell, Particle* &part_list, int &no_particles, float3 shift){
    // function returns a cell and a list of particles for a 3-dimensional ID
    // shift is a 3-vector shift applied for periodic boxes, to ensure parrticle positions are relative to cell center
    int id_1D = grid->test_cell(id_3D); // compute the 1d cell index from the 3d index
    
    if (id_1D<0) return 0; // cell is not within bounding box so will be skipped
    
    cell = grid->c[id_1D]; // cell object 
    no_particles = 0;
    // copy in list of particles into list
    for (int i = cell.start; i<cell.start+cell.np; i++, no_particles++){
        part_list[no_particles]=grid->p[i];
#ifdef PERIODIC
        part_list[no_particles].pos+= shift; // In periodic case particles positions are relative to cell center
#endif
        }        
    return id_1D;
    }
    
void init(Grid *grid, Parameters *par){
    // Initialize compute integrals class and assign relevant memory
    
    mnp=grid->maxnp;

    // assign memory to bin, xi and lists of particles
    int ec=0;
    ec+=posix_memalign((void **) &prim_list, PAGE, sizeof(Particle)*mnp); 
    ec+=posix_memalign((void **) &sec_list, PAGE, sizeof(Particle)*mnp);
    ec+=posix_memalign((void **) &thi_list, PAGE, sizeof(Particle)*mnp);
    ec+=posix_memalign((void **) &fou_list, PAGE, sizeof(Particle)*mnp);
    ec+=posix_memalign((void **) &bin, PAGE, sizeof(int)*mnp*mnp);
    ec+=posix_memalign((void **) &xi, PAGE, sizeof(Float)*mnp*mnp);

    assert(ec==0);

    // SET THIS UP FOR OMP
    gsl_rng_env_setup();
    locrng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(locrng,0); // just one thread at the moment
    


    
}
    

compute_integral2(Grid *_grid, Parameters *par){
    grid=_grid;
			
    init(grid, par); // initizalize memory etc.

    // Set the factor between two seeds of the rng of the threads
    // For deterministic rng seeds, set "steps" manually
    // Could omit ("/dev/urandom") for portability, but there seem to be implementations where random_device is deterministic
    
    int sumct=0; //Track how many iterations went into intermediate result
    double begtime; // Define global variable for the start time
   

    
// ---- Initialize correlation function and selection grid ----
    STimer initial;
    initial.Start();

    cf=new CorrelationFunction(par->corname,par->mbin,par->mumax-par->mumin);
	rd=new RandomDraws(cf,par, NULL, 0);
    Integrals sumint(par,cf);
    
     for(int ll=0;ll<mnp*mnp;ll++){
    	bin[ll]=9999999;
    }
    
    initial.Stop();
    fprintf(stderr,"Init Time: %g s\n", initial.Elapsed());
    
// ---- ---- ----

    printf("# Filled cells: %d\n",grid->nf);
    printf("# All points in use: %d\n",grid->np);
    printf("# Max points in one cell %d\n",grid->maxnp);
    fflush(NULL);

    TotalTime.Start();
    
    int max_loops = par->max_loops; // copying locally for openmp
    
    // ---- Defines intermediate output ----
    int sumupdint=(int)ceil((double)max_loops/par->nthread/1000); //Update sum in each thread every 0.1% of the runtime; should be much smaller than sumrepint
    int sumrepint=(int)ceil((double)max_loops/20);	//Report sum every x% of the full runtime; set to > n1 to avoid intermediate results
    int interval=10; // Interval when to report estimates of the remaining runtime. Changed during evaluation
// ---- ---- ---- ---- ---- ---- ---- ----

    int maxsep = ceil(par->rmax/grid->cellsize);	// The maximal distance for the pairs is determined by the upper limit of the covariance matrix to be calculated

    uint64 cnt2=0,cnt3=0,cnt4=0;
    uint64 acc2=0,acc3=0,acc4=0;


#ifdef OPENMP
	cpu_set_t mask[par->nthread+1];
	int tnum=0;
	sched_getaffinity(0, sizeof(cpu_set_t),&mask[par->nthread]);
	fprintf(stderr,"# CPUs used are: ");
	for(int ii=0;ii<64;ii++){
		if(CPU_ISSET(ii, &mask[par->nthread])){
			fprintf(stderr,"%d ",ii);
			CPU_ZERO(&mask[tnum]);
			CPU_SET(ii,&mask[tnum]);
			tnum++;
		}
	}
	fprintf(stderr,"\n");

// Begin parallel section before omp for loop to define thread local variables
#pragma omp parallel firstprivate(sumupdint,sumrepint,maxsep,steps) reduction(+:cnt2,cnt3,cnt4) shared(TotalTime,begtime,sumint,sumct)
	{
	// Decide which thread we are in
	int thread = omp_get_thread_num();
	int loopct = 0;
    assert(omp_get_num_threads()<=par->nthread);
    if (thread==0) printf("# Running on %d threads.\n", omp_get_num_threads());
#else
	int thread = 0;
    printf("# Running single threaded.\n");
#endif
    
    Integrals locint(par,cf); // accumulates the integral contribution of each thread
    
    
    uint64 locacc2=0,locacc3=0,locacc4=0;
    
#ifdef OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int n_loops = 0; n_loops<max_loops; n_loops++){
        printf("Starting Integral Loop %d of %d\n",n_loops,max_loops);

        for (int n1=0; n1<grid->nf;n1++){
            // first loop over all the filled cells in the grid for the i cell
            
            // DO WE NEED TO RANORD CELLS I?
            
            int prim_id_1D = grid->filled[n1]; // 1d ID for cell i
            prim_id = grid->cell_id_from_1d(prim_id_1D); // define the 3D primary ID
            // define the cell (primary), the list of particles in the cell (prim_list) + the number of particles (pln)
            int x1 = particle_list(prim_id, primary, prim_list, pln, {0.,0.,0.}); 
            assert(x1>=0); // nonsense else!
            
            if(primary.np==0) continue; // skip if no particles in cell
            
            p1=1.;
            
            for (int n2=0; n2<par->N2; n2++){ // Determine second position
                // DRAW SECOND CELL (sec) FROM I CELL WEIGHTED BY 1/r^2 + up to max r only???
                delta2 = rd->random_cubedraw(locrng, &p2); // draw the second cell from a cube around the first (may draw the same cell again)
                sec_id = prim_id+delta2;
                // CHECK THIS FUNCTION
                
                
                // define list of particles
                int x2=particle_list(sec_id, sec, sec_list, sln, grid->cell_sep(delta2)); 
                
                if(x2<0) continue; // reject if outside the bounding box
                
                p2*=p1; // accumulate selection probability
                
                if(sec.np==0) continue;
                
                
                // COMPUTE C2 INTEGRAL FOR ALL PARTICLES IN prim + sec + thi
                
                for (int i = 0; i<pln; i++){    // loop over i cell particles
                    for (int j = 0; j<sln; j++){    // loop over j cell particles
                        if(bin[sln*i+j]<0) continue; // skip if current pair outside binning range
                        if (sec.start==primary.start&&i==j) continue; // exclude self count
                        locint.second(bin[sln*i+j],prim_list[i],sec_list[j],p2);
                        cnt2++;
                    }
                }
                
                // NB: Worth seeing which of these params are actually used by the integrand - can remove else.
                
                for (int n3=0; n3<par->N3; n3++){ // Determine third position
                    // DRAW THIRD CELL (thi) FROM J CELL WEIGHTED BY XI
                    delta3 = rd->random_xidraw(locrng,&p3); // draw third cell from a region around cell j weighted by xi
                    thi_id=sec_id+delta3;
                    int x3=particle_list(thi_id, thi, thi_list, tln, grid->cell_sep(delta2)+grid->cell_sep(delta3)); // define list of probabilities etc.
                
                    if(x3<0) continue; // reject if outside bounding box
                    
                    p3*=p2; // update probability
                    
                    if(thi.np==0) continue; // don't need empty cells here
                    // DOES THIS CHANGE THE PROBABILITIES??
                    
                    // COMPUTE C3 INTEGRAL FOR PARTICLES IN prim+sec+thi
                
                    for (int i = 0; i<pln; i++){    // loop over i cell particles
                        for (int j = 0; j<sln; j++){    // loop over j cell particles
                            
                            if(sec.start==primary.start&&i==j) continue; // exclude self count
                            if(bin[sln*i+j]<0) continue; // skip if current pair outside binning range
                            
                            for (int k = 0; k<tln; k++) {   // loop over k cell particles
                                if (thi.start==sec.start&&j==k) continue;   // exclude self-count
                                if (thi.start==primary.start&&i==k) continue;
                                locint.third(bin[sln*i + j],xi[tln*j + k],prim_list[i],sec_list[j],thi_list[k],p3); 
                                cnt3++;
                            }
                        }
                    }
                                    
                    
                    
                    for (int n4=0; n4<par->N4; n4++){ // Determine fourth position
                        // DRAW FOURTH CELL (fou) FROM K CELL WEIGHTED BY 1/r^2
                        
                        delta4 = rd->random_cubedraw(locrng,&p4);
                        int x4=particle_list(thi_id+delta4, fou, fou_list, fln, grid->cell_sep(delta4)); 
                        
                        if(x4<0) continue;
                        
                        
                        p4*=p3;
                        
                        if(fou.np==0) continue;
                        
                        
                        // CHECK PROBABILITIES HERE
                
                        // COMPUTE C4 INTEGRAL FOR PARTICLES IN prim+sec+thi+fou
                
                        for (int i = 0; i<pln; i++){    // loop over i cell particles
                            for (int j = 0; j<sln; j++){    // loop over j cell particles
                                if(sec.start==primary.start&&i==j) continue; // exclude self count
                                if(bin[sln*i+j]<0) continue; // skip if current pair outside binning range
                            
                                for (int k = 0; k<tln; k++) {   // loop over k cell particles
                                    if (thi.start==sec.start&&j==k) continue;   // exclude self-count
                                    if (thi.start==primary.start&&i==k) continue;
                                    
                                    for (int l = 0; l<fln; l++) {
                                        if (fou.start==primary.start&&i==l) continue;   // Exclude self-count
                                        if (fou.start==sec.start&&j==l) continue;
                                        if (fou.start==thi.start&&k==l) continue;

                                        locint.fourth(bin[sln*i + j],xi[tln*j + k],prim_list[i],sec_list[j],thi_list[k],fou_list[l],p4);
                                        cnt4++;
                                    }
                                }
                            }
                        } // end particles
                    } // end l cell
                } // end k cell
            } // end j cell
        } // end i cell
        //  Update global variables from the accumulated thread variables
	if (loopct%sumupdint==0&&loopct!=0) {
#pragma omp critical(accsumint)		// Critical section to prevent several threads trying to copy their results to the global variable simultaneously
    {
    	sumint.sum_ints(&locint);
        sumct+=sumupdint;
    }
	}

	//Print estimated remaining time and current speed
    if (n_loops%interval==0) {
#pragma omp critical
    	{

    	if(interval<0.05*max_loops)
    		interval*=5;
    	else
    		interval=(int)floor(0.1*max_loops);

    	TotalTime.Stop();  // Interrupt timing to access .Elapsed()

    	if (n_loops==0) begtime=TotalTime.Elapsed();

    	uint timeleft=(uint)(TotalTime.Elapsed()-begtime)*(1./n_loops*grid->nf-1.);
    	double runtime=TotalTime.Elapsed()-begtime;
    	uint64 tmplocacc2=locacc2, tmplocacc3=locacc3, tmplocacc4=locacc4;

    	locint.sum_total_counts(locacc2, locacc3, locacc4);

    	fprintf(stderr,"# Finished cycle %d of %d after %g s. "
    			"Estimated time left: %2.2d:%2.2d:%2.2d hms "
    			"i.e. %d s. Thread %d: %ld %e %e quads/sec.\n",
    			n_loops, max_loops,runtime,timeleft/3600,timeleft/60%60,
				timeleft%60,timeleft,thread,(long int) cnt4,locacc4/runtime,cnt4/runtime);

    	locacc2=tmplocacc2, locacc3=tmplocacc3, locacc4=tmplocacc4;

    	TotalTime.Start();
    }
    }

	if (n_loops%sumrepint==0&&n_loops!=0){
#pragma omp critical(accsumint)   // Reading sumint needs the same critical section as writing to sumint
		{
		printf("# Used cells: %d/%d\n",sumct,grid->nf);
		printf("UCF: %g\n",(double)sumct/grid->nf); // Used cell fraction, needed to normalize intermediate result
#ifndef NOPRINT
        sumint.report_integrals();
#endif
        
   fprintf(stderr,"Parallely evaluated cells: %d\nTotal cells: %d\n",sumct,grid->nf);

    TotalTime.Stop();

    fprintf(stderr,"\nTime for loop: %g s\n",TotalTime.Elapsed());

    // Normalization of the accumulated results. n3 and n4 should probably better go into the probabilities p3 and p4 instead of showing up here
	sumint.normalize(/*cnt2(long)grid->np*((long)grid->np-1)*/1.,grid->np,par->nofznorm,par->N3,par->N4);
	sumint.sum_total_counts(acc2, acc3, acc4);
    
    printf("\n\nCHECK NORMALIZATIONS FOR N3, N4\n\n");

//	Print some statistics
    printf("# We tried  %lld pairs, %lld triples and %lld quads within %f.\n", cnt2, cnt3, cnt4, par->rmax);
    printf("# We accepted %lld pairs, %lld triples and %lld quads.\n", acc2, acc3, acc4);
    printf("# Acceptance ratios are %f for pairs, %f for triples and %f for quads.\n",
    		(double)acc2/cnt2, (double)acc3/cnt3, (double)acc4/cnt4);
    printf("# Average of %f pairs per primary particle.\n",
    		(Float)cnt2/grid->np);
    float x1 = par->rmax/grid->rect_boxsize.x;
    float x2 = par->rmax/grid->rect_boxsize.y;
    float x3 = par->rmax/grid->rect_boxsize.z;
    float expected = grid->np * (4*M_PI/3.0)*grid->np*x1*x2*x3;
    printf("# In a periodic cuboid box we would expect %1.0f pairs, off by a factor %f.\n", expected, cnt2/expected);

//  Print the result
#ifndef NOPRINT
    sumint.report_integrals();
#endif
    // Save the integrals to file
    sumint.save_integrals();
    
    printf("\n# Time for initialization: %g s \n"
    		"# Time for main loop: %g s\n",initial.Elapsed(),TotalTime.Elapsed());
    fflush(NULL);
}

        
		}
	}

//	Reset local variable whenever contents had been written to sumint
	if (loopct%sumupdint==0&&loopct!=0) {
    	locint.sum_total_counts(locacc2, locacc3, locacc4);
        locint.reset();
	}

    loopct++;

     
    // UPDATE INTEGRALS HERE + CHECK CONVERGENCE
    } // end cycles
    
        free(prim_list);
	free(sec_list);
	free(thi_list);
	free(fou_list);
	free(bin);
	free(xi);
	gsl_rng_free(locrng);

   
}

};

#endif
