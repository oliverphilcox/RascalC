// integrals class for grid_covariance.cpp (originally from Alex Wiegand)

#ifndef INTEGRALS_H
#define INTEGRALS_H

class Integrals {

	public:
//	RandomDraws *rd;
	CorrelationFunction *cf;

	private:
	int nbin,mbin;
	Float rmin,rmax,mumin,mumax,dr,dmu; // Ranges in r and mu
	Float *Ra, *cx, *c2, *c3, *c4; // Arrays to accumulate the integrals
	Float *Raerr, *cxerr, *c2err, *c3err, *c4err; // Arrays to accumulate the errors; however, currently used for other purposes!
	Float xiij; //Only temporary for xi^2

	uint64 *binct,*binct3,*binct4; // Arrays to accumulate the counts of the bins
	bool box,rad; // Flags to decide whether we have a periodic box and if we have a radial correlation function only

	double empty[8];   // Leftover from grid_multipoles. Probably no longer needed

	  public:
		Integrals(Parameters *par, CorrelationFunction *_cf){
			cf=new CorrelationFunction(_cf);
			init(par);
		}

//		void init_cf_rd(double*x,long n,Parameters *par) {
//			cf=new CorrelationFunction(par->corname,par->mbin,par->mumax-par->mumin);
//			if(n==0)
//				rd=new RandomDraws(cf,par,x,n);
//			else
//				rd=new RandomDraws(cf,par,x,n);
//		}

	    void init(Parameters *par) {

	    	nbin=par->nbin;
	    	mbin=par->mbin;

			int ec=0;
			// Initialize the binning
			ec+=posix_memalign((void **) &Ra, PAGE, sizeof(double)*nbin*mbin);
			ec+=posix_memalign((void **) &cx, PAGE, sizeof(double)*nbin*mbin);
			ec+=posix_memalign((void **) &c2, PAGE, sizeof(double)*nbin*mbin);
			ec+=posix_memalign((void **) &c3, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);
			ec+=posix_memalign((void **) &c4, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);

			ec+=posix_memalign((void **) &Raerr, PAGE, sizeof(double)*nbin*mbin);
			ec+=posix_memalign((void **) &cxerr, PAGE, sizeof(double)*nbin*mbin);
			ec+=posix_memalign((void **) &c2err, PAGE, sizeof(double)*nbin*mbin);
			ec+=posix_memalign((void **) &c3err, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);
			ec+=posix_memalign((void **) &c4err, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);

			ec+=posix_memalign((void **) &binct, PAGE, sizeof(uint64)*nbin*mbin);
			ec+=posix_memalign((void **) &binct3, PAGE, sizeof(uint64)*nbin*mbin*nbin*mbin);
			ec+=posix_memalign((void **) &binct4, PAGE, sizeof(uint64)*nbin*mbin*nbin*mbin);

			assert(ec==0);

			reset();

			box=par->perbox;

			rmax=par->rmax;
			rmin=par->rmin;

			mumax=par->mumax;
			mumin=par->mumin;

			dr=(rmax-rmin)/nbin;
			dmu=(mumax-mumin)/mbin;

	    	rad=mbin==1&&dmu==1.;

			empty[0] = 0.0;   // To avoid a warning


		}

		~Integrals() {
			free(Ra);
			free(cx);
			free(c2);
			free(c3);
			free(c4);
			free(c3err);
			free(c4err);
			free(binct);
			free(binct3);
			free(binct4);

	    }

	    void reset(){
	    	for (int j=0; j<nbin*mbin; j++) {
				Ra[j] = 0;
				cx[j] = 0;
				c2[j] = 0;
				binct[j] = 0;
				Raerr[j] = 0;
				cxerr[j] = 0;
				c2err[j] = 0;
			}
			for (int j=0; j<nbin*mbin*nbin*mbin; j++) {
				c3[j]=0;
				c4[j]=0;
				binct3[j] = 0;
				binct4[j] = 0;
				c3err[j]=0;
				c4err[j]=0;
			}
	    }

	    inline int getbin(Float r, Float mu){
	    	// Linearizes 2D indices
	       	return floor((r-rmin)/dr) *mbin + floor((mu-mumin)/dmu);
	    }

	    inline void second(int& bin, const Particle& pi, const Particle& pj, Float pro) {
	    	// Accumulates the two point integral C2
	    	// Note that the value of bin is returned
			Float rij_mag,rij_mu;
			Float rav,c2v,cxv, weight;

//			Bin will be reused later
	        cleanup_l(pi.pos, pj.pos, rij_mag, rij_mu);
	        bin=getbin(rij_mag,rij_mu);

	        if ((bin>=0) && (bin<nbin*mbin)){
	        	weight=pi.w*pj.w;
		        xiij=cf->xi(rij_mag,rij_mu);
	        	rav=weight / pro;
	        	c2v= weight*weight *(1+xiij) / pro;
	        	cxv = weight*xiij / pro;

	        	Ra[bin]+= rav;
	        	c2[bin]+= c2v;
	        	cx[bin]+= cxv;
	        	binct[bin]++;
	        	Raerr[bin]+= pow(rav,2);
//	        	c2err[bin]+= pow(c2v,2);
	        	c2err[bin]+= weight*weight / pro * pow(xiij,2); //Misuse to accumulate xi^2 term
	        	cxerr[bin]+= pow(cxv,2);
	        }else{
	        	bin=-1;
	        }

	    }

	    inline void third(const int& bin, Float& xijk, const Particle& pi, const Particle& pj, const Particle& pk, Float p3) {
	    	// Accumulates the three-point integral C3
	    	// Returns the correlation function between point j and k
	    	int binik;
	    	Float c3l;
			Float rik_mag,rik_mu, rjk_mag,rjk_mu, xiik;

			cleanup_l(pi.pos, pk.pos, rik_mag, rik_mu);
			binik=getbin(rik_mag, rik_mu);
			xiik=cf->xi(rik_mag, rik_mu); //Only temporary for xi^2; can be commented out if c3err is used for error as intended

//			Correlation function will be reused in fourth
			cleanup_l(pj.pos, pk.pos, rjk_mag, rjk_mu);
	        xijk=cf->xi(rjk_mag, rjk_mu);

//			Checking if ik bin is in chosen range and assuming that "bin" is the correct ij bin
			if ((binik<nbin*mbin) && (binik>=0)){
				c3l=pi.w*pi.w*pj.w*pk.w/p3 * xijk;
				c3[bin*nbin*mbin+binik]+=c3l;
//				c3err[bin*nbin*mbin+binik]+=pow(c3l,2);
				c3err[bin*nbin*mbin+binik]+=pi.w*pi.w*pj.w*pk.w/p3 * 0.5*(xijk*xiik+xijk*xiij); //Misuse to accumulate xi^2 term
				binct3[bin*nbin*mbin+binik]++;
			}


		}

	    inline void fourth(const int& bin, Float& xijk, const Particle& pi, const Particle& pj, const Particle& pk, const Particle& pl, Float p4) {
	    	// Accumulates the four-point integral C4
	    	// Input: The ij bin, the correlation function between jk, the four particle positions
	    	int binkl;
	    	Float c4l,xiil;
	    	Float ril_mag,ril_mu, rkl_mag,rkl_mu;

	    	cleanup_l(pk.pos, pl.pos, rkl_mag, rkl_mu);
	    	binkl=getbin(rkl_mag,rkl_mu);

//			Checking if kl bin is in chosen range and assuming that "bin" is the correct ij bin
	    	if ((binkl<nbin*mbin) && (binkl>=0)){

		    	cleanup_l(pi.pos, pl.pos, ril_mag, ril_mu);
		    	xiil=cf->xi(ril_mag, ril_mu);

	    		c4l=pi.w*pj.w*pk.w*pl.w/p4  * xijk  * xiil;
	    		c4[bin*nbin*mbin+binkl]+=c4l;
//	    		c4err[bin*nbin*mbin+binkl]+=pow(c4l,2);
	    		c4err[bin*nbin*mbin+binkl]+=pi.w*pj.w*pk.w*pl.w/p4; //Misuse of c4err to test normalization against 2pt randoms
	    		binct4[bin*nbin*mbin+binkl]++;
	    	}

	    }

	    void sum_ints(Integrals* ints) {
	    	// Add the values accumulated in ints to the corresponding internal sums
			for(int i=0;i<nbin*mbin;i++){
				Ra[i]+=ints->Ra[i];
				c2[i]+=ints->c2[i];
				cx[i]+=ints->cx[i];
				binct[i]+=ints->binct[i];
				Raerr[i]+=ints->Raerr[i];
				c2err[i]+=ints->c2err[i];
				cxerr[i]+=ints->cxerr[i];
			}

			for(int i=0;i<nbin*mbin;i++){
				for(int j=0;j<nbin*mbin;j++){
					c3[i*nbin*mbin+j]+=ints->c3[i*nbin*mbin+j];
					c4[i*nbin*mbin+j]+=ints->c4[i*nbin*mbin+j];
					c3err[i*nbin*mbin+j]+=ints->c3err[i*nbin*mbin+j];
					c4err[i*nbin*mbin+j]+=ints->c4err[i*nbin*mbin+j];
					binct3[i*nbin*mbin+j]+=ints->binct3[i*nbin*mbin+j];
					binct4[i*nbin*mbin+j]+=ints->binct4[i*nbin*mbin+j];
				}
			}
		}

	    void sum_total_counts(uint64& acc2, uint64& acc3, uint64& acc4){
	    	// Add local counts to total bin counts in acc2-4
//	    	acc2=0,acc3=0,acc4=0;
	    	for (int i=0; i<nbin*mbin; i++) {
	    		acc2+=binct[i];
				for (int j=0; j<nbin*mbin; j++) {
					acc3+=binct3[i*nbin*mbin+j];
					acc4+=binct4[i*nbin*mbin+j];
				}
			}
	    }

	    void normalize(long n,int np,int ngal, int n3, int n4) {
	    	// Normalize the acccumulated integrals; partially already done by normalization of the proposal function of the cubes selected
	    	// np is the number of random particles used, ngal is the number of galaxies in the survey
	    	double corrf=pow((double)np/ngal,2); // Correction factor for the different densities of randoms and points

	    	for(int i=0;i<nbin*mbin;i++){

//	    		printf("Predefined Integrals %2d %le +- %le   %le +- %le   %le +- %le\n",
//						i, Ra[i], Raerr[i], cx[i], cxerr[i], c2[i], c2err[i]);

	    		Ra[i]/=((Float)n*corrf);
	    		c2[i]/=(Ra[i]*Ra[i]*(Float)n*corrf);
				cx[i]/=(Ra[i]*(Float)n*corrf);

// TODO: Figure out how to normalize the quadratic counts of err to arrive at an estimate for the integration error
//				Raerr[i]=sqrt((Raerr[i]/((Float)n)-pow(Ra[i]*corrf,2))/(Float)(n-1))/corrf;
///*wrong->*/		c2err[i]=sqrt((c2err[i]/(pow((Ra[i]*Ra[i]*corrf*corrf),2)*(Float)n)-pow(c2[i]/corrf,2))/(Float)(n-1))/corrf;
///*wrong->*/		cxerr[i]=sqrt((cxerr[i]/(pow((Ra[i]*pow((double)np/ngal,2)),2)*(Float)n)-pow(cx[i],2))/(Float)(n-1));

//				printf("Postdefined Integrals %2d %le +- %le   %le +- %le   %le +- %le\n",
//						i, Ra[i], Raerr[i], cx[i], cxerr[i], c2[i], c2err[i]);


	    	}

	    	for(int i=0;i<nbin*mbin;i++){
	    		for(int j=0;j<nbin*mbin;j++){
	    			c3[i*nbin*mbin+j]/=(n3*Ra[i]*Ra[j]*(Float)n*pow((double)np/ngal,3));
	    			c4[i*nbin*mbin+j]/=(n3*n4*Ra[i]*Ra[j]*(Float)n*pow((double)np/ngal,4));
	    		}
			}
	    }
        
        void save_integrals() {
        /* Print integral outputs to file. 
         * In txt files {c2,c3,c4,RR}_n{nbin}_m{mbin}.txt there are lists of the outputs of c2,c3,c4 and RR_a that are already normalized and multiplied by combinatoric factors. The n and m strings specify the number of n and m bins present.
         */
            // Create output files
            char c2name[1000];
            snprintf(c2name, sizeof c2name, "CovMatrices/c2_n%d_m%d.txt", nbin, mbin);
            char c3name[1000];
            snprintf(c3name, sizeof c3name, "CovMatrices/c3_n%d_m%d.txt", nbin, mbin);
            char c4name[1000];
            snprintf(c4name, sizeof c4name, "CovMatrices/c4_n%d_m%d.txt", nbin, mbin);
            char RRname[1000];
            snprintf(RRname, sizeof RRname, "CovMatrices/RR_n%d_m%d.txt", nbin, mbin);
            FILE * C2File = fopen(c2name,"w"); // for c2 part of integral
            FILE * C3File = fopen(c3name,"w"); // for c3 part of integral
            FILE * C4File = fopen(c4name,"w"); // for c4 part of integral
            FILE * RRFile = fopen(RRname,"w"); // for RR part of integral
            
            for (int j=0;j<nbin*mbin;j++){
                fprintf(C2File,"%le\n",c2[j]);
                fprintf(RRFile,"%le\n",Ra[j]);
            }
            for(int i=0;i<nbin*mbin;i++){
                for(int j=0;j<nbin*mbin;j++){
                    fprintf(C3File,"%le\n",c3[i*nbin*mbin+j]);
                    fprintf(C4File,"%le\n",c4[i*nbin*mbin+j]);
                }
            }
            printf("Printed output to file in the CovMatrices/ directory");            
        }
	    void report_integrals() {
	///*
		for (int j=0; j<nbin*mbin; j++) {
		    printf("C2 %2d %le +- %le   %le +- %le   %le +- %le\n",
				j, Ra[j], Raerr[j], cx[j], cxerr[j], 2*c2[j], 2*c2err[j]);
		}
		printf("\n");


		for(int i=0;i<nbin*mbin;i++){
			printf("C3 %2d ",i);
			for(int j=0;j<nbin*mbin;j++){
				printf("%e ", 4*c3[i*nbin*mbin+j]);
			}
			printf("\n");
		}
		printf("\n");

		for(int i=0;i<nbin*mbin;i++){
			printf("C4 %2d ",i);
			for(int j=0;j<nbin*mbin;j++){
				printf("%e ",2*c4[i*nbin*mbin+j]);
			}
			printf("\n");
		}
		printf("\n");

		for(int i=0;i<nbin*mbin;i++){
			printf("E3 %2d ",i);
			for(int j=0;j<nbin*mbin;j++){
				printf("%e ",4*c3err[i*nbin*mbin+j]);
			}
			printf("\n");
		}
		printf("\n");

		for(int i=0;i<nbin*mbin;i++){
			printf("E4 %2d ",i);
			for(int j=0;j<nbin*mbin;j++){
				printf("%e ",2*c4err[i*nbin*mbin+j]);
			}
			printf("\n");
		}
		printf("\n");

		for (int i=0; i<nbin; i++) {
			printf("N2 %2d ",i);
			for (int j=0; j<mbin; j++) {
				printf("%lld ", binct[i*mbin+j]);
			}
			printf("\n");
		}
		printf("\n");

		for (int i=0; i<nbin*mbin; i++) {
			printf("N3 %2d ",i);
			for (int j=0; j<nbin*mbin; j++) {
				printf("%lld ", binct3[i*nbin*mbin+j]);
			}
			printf("\n");
		}
		printf("\n");

		for (int i=0; i<nbin*mbin; i++) {
			printf("N4 %2d ",i);
			for (int j=0; j<nbin*mbin; j++) {
				printf("%lld ", binct4[i*nbin*mbin+j]);
			}
			printf("\n");
		}
		printf("\n");

		fflush(NULL);

	//*/
	/*
		printf("\n{");
		printf("%.0f", xi0[0]);
		for (int j=1; j<nbin; j++) {
		    printf(",%.0f", xi0[j]);
		}
		printf("}\n");
	*/

	    }

	  private:
	    void cleanup_l(Float3 p1,Float3 p2,Float& norm,Float& mu){
	    	Float3 pos=p1-p2;
	    	norm = pos.norm();
	    	// If the input correlation function had only radial information and no mu bins
	    	// fill the dummy mu bins by 0.5
	    	if(rad){
	    		mu=0.5;
	    	}
	    	else{
#ifndef PERIODIC
				Float3 los=p1+p2; // No 1/2 as normalized anyway below
	//	    	if(!box){
					//Absolute value as mu bins defined in terms of |mu|
					mu = fabs(pos.dot(los)/norm/los.norm());
	//	    	}
#else
	    		// In the periodic case use z-direction for mu
				mu = fabs(pos.z/norm);
#endif
	    	}
	    }
	};
#endif
