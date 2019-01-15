// compute_integral class for grid_covariance.cpp (originally from Alex Wiegand)

#ifndef COMPUTE_INTEGRAL_H
#define COMPUTE_INTEGRAL_H

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


void compute_integral(Grid *grid, Parameters *par) {

	/* -------------------------------------------------------- */
	/* Calculates the integral based on the random particles 	*/
	/* in the grid              								*/
	/* -------------------------------------------------------- */

   

 
  
    

// ---- Defines the sampling strategy ----

#define SECSAMP ALL		// Strategy for the selection of the second cell ...
#define THISAMP SEL		// the third cell ...
#define FOUSAMP SEL		// and the fourth cell

	// In the following define the number of cells selected at each stage
	uint n1=grid->nf; // Currently use all grid cells as first grid cell once. Could use less than all filled cells here
#if SECSAMP!=ALL
	uint n2=100;      // Sets the number of cells used as secondary cells. Not effective if strategy is ALL
#endif
	uint n3=1;		  // Number of cells selected as third cell for each pair of primary and secondary cell
	uint n4=1;		  // Number of cells selected as fourth cell for each triple
					  // Setting those to 1 currently also avoids selecting a quad twice

// ---- ---- ---- ---- ---- ---- ---- ----


// ---- Defines intermediate output ----
    int sumupdint=(int)ceil((double)n1/par->nthread/1000); //Update sum in each thread every 0.1% of the runtime; should be much smaller than sumrepint
    int sumrepint=(int)ceil((double)n1/20);	//Report sum every x% of the full runtime; set to > n1 to avoid intermediate results
    int interval=10; // Interval when to report estimates of the remaining runtime. Changed during evaluation
// ---- ---- ---- ---- ---- ---- ---- ----

	int maxsep = ceil(par->rmax/grid->cellsize);	// The maximal distance for the pairs is determined by the upper limit of the covariance matrix to be calculated
#if THISAMP != SEL
	int maxsep3=ceil(par->xicutoff/grid->cellsize);	// Not used for SEL. The maximal distance for the third cell is determined from the xicutoff
#endif
#if FOUSAMP !=SEL
	int maxsep4=maxsep;								// Not used for SEL.
#endif

#if SECAMP == FLAT
	int aneigh=(int)pow(2*maxsep+1,3);		// Number of all neighbour cells up to max radius
#endif
#if THISAMP == FLAT
    int aneigh3=(int)pow(2*maxsep3+1,3);
#endif
#if FOUSAMP == FLAT
	int aneigh4=(int)pow(2*maxsep4+1,3);
#endif
    
    uint64 cnt2=0,cnt3=0,cnt4=0;
    uint64 acc2=0,acc3=0,acc4=0;
//    uint64 mean=(uint64)floor((double)grid->np/grid->nf);

    int sumct=0; //Track how many iterations went into intermediate result
    double begtime; // Define global variable for the start time

    // Set the factor between two seeds of the rng of the threads
    // For deterministic rng seeds, set "steps" manually
    // Could omit ("/dev/urandom") for portability, but there seem to be implementations where random_device is deterministic
    std::random_device urandom("/dev/urandom");
    std::uniform_int_distribution<unsigned int> dist(1, std::numeric_limits<unsigned int>::max());
    unsigned long int steps = dist(urandom);

// ---- Initialize correlation function and selection grid ----
    STimer initial;
    initial.Start();
    
    gsl_rng_env_setup();
    CorrelationFunction *cf=new CorrelationFunction(par->corname,par->mbin,par->mumax-par->mumin);
	RandomDraws *rd=new RandomDraws(cf,par, NULL, 0);
    Integrals sumint(par,cf);

    // Shuffles the filled cells. Only needed if we do not use all cells in the outermost loop
    std::vector<int> cellsranord(grid->filled, grid->filled + grid->nf);
    std::random_shuffle ( cellsranord.begin(), cellsranord.end() );

    initial.Stop();
    fprintf(stderr,"Init time: %g s\n",initial.Elapsed());
// ---- ---- ----

    printf("# Filled cells: %d\n",grid->nf);
    printf("# All points in use: %d\n",grid->np);
    printf("# Max points in one cell %d\n",grid->maxnp);
    fflush(NULL);

    TotalTime.Start();

    
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
#pragma omp parallel firstprivate(sumupdint,sumrepint,maxsep,n1,steps) reduction(+:cnt2,cnt3,cnt4) shared(TotalTime,begtime,sumint,sumct)
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

    Integrals locint(par,cf); // Accumulates the integral contribution of each thread
	gsl_rng* locrng = gsl_rng_alloc( gsl_rng_default ); // One rng per thread
	gsl_rng_set(locrng,steps*(thread+1)); // seeds of the different rngs "steps" apart

    uint64 locacc2=0,locacc3=0,locacc4=0;

    int mnp=grid->maxnp;
    int pln,sln,tln,fln; // this is number of particles in the primary/secondary/etc. list 
    Particle *prim_list;
    Particle *sec_list;
    Particle *thi_list;
    Particle *fou_list;
    Cell zero({0,0});

    Float* xi;
    int* bin;

    bool changeij;
    bool changek;

    int ec=0;
    ec+=posix_memalign((void **) &prim_list, PAGE, sizeof(Particle)*mnp); // Not sure if memalign has a benefit over malloc here
    ec+=posix_memalign((void **) &sec_list, PAGE, sizeof(Particle)*mnp);
    ec+=posix_memalign((void **) &thi_list, PAGE, sizeof(Particle)*mnp);
    ec+=posix_memalign((void **) &fou_list, PAGE, sizeof(Particle)*mnp);

    ec+=posix_memalign((void **) &bin, PAGE, sizeof(int)*mnp*mnp);
    ec+=posix_memalign((void **) &xi, PAGE, sizeof(Float)*mnp*mnp);

    assert(ec==0);

    for(int ll=0;ll<mnp*mnp;ll++){
    	bin[ll]=9999999;
    }

#if SECSAMP == FLAT
		int* cellsel2 = (int*)malloc(sizeof(int)*aneigh);
		for(int ll=0;ll<aneigh;ll++){
			cellsel2[ll]=ll;
		}
#endif
#if THISAMP == FLAT
		int* cellsel3 = (int*)malloc(sizeof(int)*aneigh3);
		for(int ll=0;ll<aneigh3;ll++){
			cellsel3[ll]=ll;
		}
#endif
#if FOUSAMP == FLAT
		int* cellsel4 = (int*)malloc(sizeof(int)*aneigh4);
		for(int ll=0;ll<aneigh4;ll++){
			cellsel4[ll]=ll;
		}
#endif

//    for(int i=0;i<=1000;i++){
//    	for(int j=0;j<=1000;j++){
//    	fprintf(stdout,"%e ",locint.cf->xi(i,j*0.001));
//    	}
//    	fprintf(stdout,"\n");
//    }
//    fprintf(stdout,"\n");
//
//    exit(0);


// ---- Calculation starts here ----

// Loop over primary cells.
#ifdef OPENMP
#pragma omp	for schedule(dynamic)
#endif
   for (uint n=0; n<n1; n++) {
#ifdef OPENMP
	// Manually pin threads to CPUs (necessary on Odyssey)
//	if(CPU_COUNT(&mask[thread])!=0){
//		sched_setaffinity(0, sizeof(cpu_set_t),&mask[thread]);
//		CPU_ZERO(&mask[thread]);
//	}
	//fprintf(stderr,"# Thread %d on cpu %d calculates cell %d and mask is ",omp_get_thread_num(),sched_getcpu(),n);
	//sched_getaffinity(0, sizeof(cpu_set_t),&mask[thread]);
	//for(int ii=0;ii<64;ii++){
	//	if(CPU_ISSET(ii, &mask[thread]))
	//		fprintf(stderr,"%d ",ii);
	//}
	//fprintf(stderr,".\n");
#endif

        // Determine first cell; Cells randomly ordered so we can use only a limited number just by stopping the loop
		int m=cellsranord[n];
		integer3 prim_id = grid->cell_id_from_1d(m);
		Cell primary = grid->c[m];
		double p1=n1/grid->nf;

		if(primary.np==0) continue; // If cell empty skip (should not happen any more)

		// Copy global list to local list to enhance memory access pattern
		pln=0;
		for (int i = primary.start; i<primary.start+primary.np; i++,pln++) {
			prim_list[pln]=grid->p[i];
		}

//		Determine second position
		double p2; 
		integer3 delta;

// ---- Depending on sampling mode select second cell ----
#if SECSAMP == FLAT
        int fin2 = fmin(aneigh,n2);
		kofnip(cellsel2,aneigh,fin2,locrng); 	// Draw second cell equiprobably from all neighbouring cells
		p2=1./aneigh*fin2;						// fin2 points selected out of aneigh
        for (int ci=0;ci<fin2;ci++){
            printf("");
            delta=relid(cellsel2[ci],maxsep);
#elif SECSAMP == ALL
		for (delta.x = -maxsep; delta.x <= maxsep; delta.x++)	// Iterate through all other cells up to maxsep
		for (delta.y = -maxsep; delta.y <= maxsep; delta.y++)
		for (delta.z = -maxsep; delta.z <= maxsep; delta.z++){
			p2=1.; // All points selected, so probability 1
#elif SECSAMP == SEL
			for (uint ci=0;ci<n2;ci++){
	    		delta=rd->random_cubedraw(&p2); // Draw second cell from cube around first; May draw same cell twice
#endif

	    	integer3 sec_id = prim_id+delta;
	       	int id=grid->test_cell(sec_id);
			if(id<0) continue;					// If the cell is outside the bounding box
			Cell sec = grid->c[id];

			if(sec.np==0) continue; // If second cell empty there will be neither pairs, triples nor quads
									// Can't do that for the inner loops (thi and fou) as we need all pairs!

			// Copy global list to local list to enhance memory access pattern
			sln=0;
			for (int i = sec.start; i<sec.start+sec.np; i++,sln++) {
				sec_list[sln]=grid->p[i];
#ifdef PERIODIC
				sec_list[sln].pos+= grid->cell_sep(delta); // In periodic case particles positions are relative to cell center
#endif
			}

			p2*=p1;		   // Accumulates selection probability
			changeij=true; // Tracks if we changed the pair ij by selecting a new primary and/or secondary cell


//			Determine third position
			integer3 delta3;
	    	double p3;
// ---- Depending on sampling mode select third cell ----
#if THISAMP == SEL
				for (uint ck=0;ck<n3;ck++){
					delta3=rd->random_xidraw(locrng,&p3);
					integer3 thi_id=sec_id+delta3;			// select with respect to second cell
#elif THISAMP == FLAT
	    		int fin3 = fmin(aneigh3,n3);
	    		kofnip(cellsel3,aneigh3,fin3,locrng);
				for (int ck=0;ck<fin3;ck++){
		    		delta3=relid(cellsel3[ck],maxsep3);
		       		integer3 thi_id=sec_id+delta3;
					p3=1./aneigh3*fin3;
#elif THISAMP == ALL
				int ck=0;
				for (delta3.x = -maxsep3; delta3.x <= maxsep3; delta3.x++)
				for (delta3.y = -maxsep3; delta3.y <= maxsep3; delta3.y++)
				for (delta3.z = -maxsep3; delta3.z <= maxsep3; delta3.z++,ck++){
					integer3 thi_id = sec_id+delta3;
					p3=1.;
//					integer3 thi_id = prim_id+delta3;
#endif
					p3*=p2;
					id=grid->test_cell(thi_id);
		       		Cell thi;
		       		if(id<0)
		       			thi = zero;
		       		else
		       			thi = grid->c[id];

//					Can't skip empty cells at this stage because we need all pairs, independently of third draw
//					if(thi.np==0) continue;

					// Copy global list to local list to enhance memory access pattern
		       		tln=0;
					for (int i = thi.start; i<thi.start+thi.np; i++,tln++) {
						thi_list[tln]=grid->p[i];
#ifdef PERIODIC
//						thi_list[tln].pos += grid->cell_sep(delta3); //From cell 1
						thi_list[tln].pos += grid->cell_sep(delta) + grid->cell_sep(delta3); //From cell 2
#endif
					}
					changek=true; // Tracks if third cell has changed

//			Determine fourth position
					double p4;
					integer3 delta4;
// ---- Depending on sampling mode select fourth cell ----
#if FOUSAMP == SEL
					for (uint cl=0;cl<n4;cl++){
						// Determine fourth position
						delta4=rd->random_xidraw(locrng,&p4);
						integer3 fou_id = prim_id+delta4;
	//					delta4=rd->random_cubedraw(&p4);
	//					integer3 fou_id = thi_id+delta4;
#elif FOUSAMP == ALL
					int cl=0;

					for (delta4.x = -maxsep4; delta4.x <= maxsep4; delta4.x++)
					for (delta4.y = -maxsep4; delta4.y <= maxsep4; delta4.y++)
					for (delta4.z = -maxsep4; delta4.z <= maxsep4; delta4.z++,cl++){
						integer3 fou_id = thi_id+delta4;	// Samples around third cell!
						p4=1.;

#elif FOUSAMP == FLAT
					int fin4 = fmin(aneigh4,n4);
					kofnip(cellsel4,aneigh4,fin4,locrng);
					for (int cl=0;cl<fin4;cl++){
						delta4=relid(cellsel4[cl],maxsep4);
						integer3 fou_id = thi_id+delta4;	// Samples around third cell!
						p4=1./aneigh4*fin4;
#endif
					p4*=p3;
					id=grid->test_cell(fou_id);
					Cell fou;
					if(id<0)
						fou = zero;
					else
						fou = grid->c[id];

//					Not possible with current loopstructure (would impact 2pt part)
//					if(fou.np==0) continue;

					// Copy global list to local list to enhance memory access pattern
					fln=0;
					for (int i = fou.start; i<fou.start+fou.np; i++,fln++) {
						fou_list[fln]=grid->p[i];
#ifdef PERIODIC
#if FOUSAMP == SEL
						fou_list[fln].pos += grid->cell_sep(delta4); // delta4 from cell 1
#else
						fou_list[fln].pos += grid->cell_sep(delta) + grid->cell_sep(delta3) + grid->cell_sep(delta4); // delta4 from cell 3
#endif
#endif
					}

// 		---- Now we loop over all particles in the chosen cells ----
					for (int i = 0; i<pln; i++) {
//						if(sec.np==0) cnt2+=mean;
						for (int j = 0; j<sln; j++) {
		    				if (sec.start==primary.start&&i==j) continue;   // Exclude self-count

		    				if(changeij){
		    					locint.second(bin[sln*i + j],prim_list[i],sec_list[j],p2); // bin is io here!
		    					cnt2++;
		    				}
		    				if(bin[sln*i + j]<0) continue; // If current pair outside of binning range: skip
//		    				continue;

	    				for (int k = 0; k<tln; k++) {
							if (thi.start==sec.start&&j==k) continue;   // Exclude self-count
							if (thi.start==primary.start&&i==k) continue;

							if(changek){
								locint.third(bin[sln*i + j],xi[tln*j + k],prim_list[i],sec_list[j],thi_list[k],p3); // xi is io here!
								cnt3++;
							}
//		    				continue;

						for (int l = 0; l<fln; l++) {
							if (fou.start==primary.start&&i==l) continue;   // Exclude self-count
							if (fou.start==sec.start&&j==l) continue;
							if (fou.start==thi.start&&k==l) continue;

							locint.fourth(bin[sln*i + j],xi[tln*j + k],prim_list[i],sec_list[j],thi_list[k],fou_list[l],p4);
							cnt4++;

						} // Done with fourth cell
	    				} // Done with third cell
						} // Done with second cell
					} // Done with first cell
//			cnt4++;
			changeij=false;
			changek=false;

			if(tln==0) break; // Now that the pairs are taken into account, if the third cell is empty do no longer continue to select fourth cells to the empty third cell

			} // Done with fourth choice
		} // Done with third choice
		} // Done with second choice


//  Update global variables from the accumulated thread variables
	if (loopct%sumupdint==0&&loopct!=0) {
#pragma omp critical(accsumint)		// Critical section to prevent several threads trying to copy their results to the global variable simultaneously
    {
    	sumint.sum_ints(&locint);
        sumct+=sumupdint;
    }
	}

	//Print estimated remaining time and current speed
    if (n%interval==0) {
#pragma omp critical
    	{

    	if(interval<0.05*n1)
    		interval*=5;
    	else
    		interval=(int)floor(0.1*n1);

    	TotalTime.Stop();  // Interrupt timing to access .Elapsed()

    	if (n==0) begtime=TotalTime.Elapsed();

    	uint timeleft=(uint)(TotalTime.Elapsed()-begtime)*(1./n*grid->nf-1.);
    	double runtime=TotalTime.Elapsed()-begtime;
    	uint64 tmplocacc2=locacc2, tmplocacc3=locacc3, tmplocacc4=locacc4;

    	locint.sum_total_counts(locacc2, locacc3, locacc4);

    	fprintf(stderr,"# Finished cell %d of %d after %g s. "
    			"Estimated time left: %2.2d:%2.2d:%2.2d hms "
    			"i.e. %d s. Thread %d: %ld %e %e quads/sec.\n",
    			n, grid->nf,runtime,timeleft/3600,timeleft/60%60,
				timeleft%60,timeleft,thread,(long int) cnt4,locacc4/runtime,cnt4/runtime);

    	locacc2=tmplocacc2, locacc3=tmplocacc3, locacc4=tmplocacc4;

    	TotalTime.Start();
    }
    }

	if (n%sumrepint==0&&n!=0){
#pragma omp critical(accsumint)   // Reading sumint needs the same critical section as writing to sumint
		{
		printf("# Used cells: %d/%d\n",sumct,n1);
		printf("UCF: %g\n",(double)sumct/n1); // Used cell fraction, needed to normalize intermediate result
#ifndef NOPRINT
        sumint.report_integrals();
#endif
		}
	}

//	Reset local variable whenever contents had been written to sumint
	if (loopct%sumupdint==0&&loopct!=0) {
    	locint.sum_total_counts(locacc2, locacc3, locacc4);
        locint.reset();
	}

    loopct++;

	} // Done with first choice; End of omp for loop

	free(prim_list);
	free(sec_list);
	free(thi_list);
	free(fou_list);
	free(bin);
	free(xi);
	gsl_rng_free(locrng);

#if SECSAMP == FLAT
        free(cellsel2);
#endif
#if THISAMP == FLAT
        free(cellsel3);
#endif
#if FOUSAMP == FLAT
        free(cellsel4);
#endif

// After the end of the loop add all remaining contents of the local threads
#ifdef OPENMP
#pragma omp critical(accsumint)
    {
		sumint.sum_ints(&locint);
// TODO: sumct sometimes wrong. Maybe because of loopct++ after  sumct+=sumupdint; ? Should not be a problem anyway because sumct not used outside the loop
	    sumct+=loopct%sumupdint;
    }
	} // End of omp parallel pragma
#endif

    fprintf(stderr,"Parallely evaluated cells: %d\nTotal cells: %d\n",sumct,n1);

    TotalTime.Stop();

    fprintf(stderr,"\nTime for loop: %g s\n",TotalTime.Elapsed());

    // Normalization of the accumulated results. n3 and n4 should probably better go into the probabilities p3 and p4 instead of showing up here
	sumint.normalize(/*cnt2(long)grid->np*((long)grid->np-1)*/1.,grid->np,par->nofznorm,n3,n4);
	sumint.sum_total_counts(acc2, acc3, acc4);

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
    return;
}

#endif
