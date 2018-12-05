// random draws class for grid_covariance.cpp (originally from Alex Wiegand, rewritten by Oliver Philcox)
#include "correlation_function.h"
#include "parameters.h"
#include <gsl/gsl_sf_dawson.h>

#ifndef RANDOM_DRAWS2_H
#define RANDOM_DRAWS2_H





class RandomDraws2{

	// Class that handles the selection of neighbouring boxes

public:
	static CorrelationFunction *corr;
	int nside;     // Number of cells in each direction of large draw
	int nsidecube; // Number of cells in each direction of maxsep cube
	double boxside;
	double *x; // Probability grid for 1/r^2 kernel
	double *xcube; // Probability grid for xi(r) kernel

	private:
		// Sampling of long distance
		ransampl_ws* ws;

		// Sampling of short distance
		ransampl_ws* cube;


	public:

		RandomDraws2(CorrelationFunction *fun,Parameters *par,const double *xin, long np){ 

		long n=np;
		corr=fun;

        // Define grid of nside up to the maximum xi cut-off scale (as a cubic grid)
        Float box_max = fmax(fmax(par->rect_boxsize.x,par->rect_boxsize.y),par->rect_boxsize.z);
        boxside=box_max/par->nside;
        nside=2*ceil(par->xicutoff/boxside)+1; 
        
		// If precalculated grid has been saved, load it
		if (par->loadname!=NULL&&(xin==NULL||np==0)){
			int len = strlen(par->loadname);
			if(len>4&&!strcmp(&par->loadname[len-4],".bin")){
				// Read binary format that was written by writeDataBin
                printf("\n# Loading pre-computed binary-format probability grid\n");
				readDataBin(&x,&n,par->loadname);
			}
			else{
				// Read ascii format possibly created externally
                printf("\n# Loading pre-computed ascii-format probability grid\n");
				readData(&x,&n,par->loadname);
			}
			if(n==0){
                printf("\n# Computing the probability grid\n");
				integData(&x,&n,f_corr,nside,boxside,0);
                printf("# Probability grid computation complete\n");
            }
		}
		else
			// If no correlation function is given to copy, do integration
			if (xin==NULL||np==0){
				printf("\n# Computing the probability grid");
                //TODO: Change this line
                printf("\n\nStarting new integration:\n");
                integData2(&x,&n,xi_integrand,nside,boxside); // could change to 1/r^2 here by using f instead of f_corr and 1 instead of 0
                printf("# Probability grid computation complete\n");
                exit(1);
            }
			else
				copyData(&x,&n,xin);

		// Save grid of precalculated probabilities
		if (par->savename!=NULL){
			int len = strlen(par->savename);
			if(len>4&&!strcmp(&par->savename[len-4],".bin"))
				writeDataBin(&x,&n,par->savename);
			else
				fprintf(stderr,"Save file %s does not end in \".bin\". No output written.\n",par->savename);
		}

		// Set up actual sampler
		ws = ransampl_alloc( n );
		ransampl_set( ws, x );

		// Normalize grid probabilities to one
		double sum=0.;
		for(int i=0;i<n;i++){
			sum+=x[i];
		}
		for(int i=0;i<n;i++){
			x[i]/=sum;
		}


//		Initialize second sampler

		int maxsep = ceil(par->rmax/boxside);
		nsidecube = 2 * maxsep + 1;
		long nn=0;

        printf("\nUsing xi(r) sampling for j-k cells\n");
        printf("Using 1/r^2 sampling for i-j and k-l cells\n");
        //TODO: Change this
        integData(&xcube,&nn,f,nsidecube,boxside,1.);// (previously 1./2.);

		// Set up actual sampler
		cube = ransampl_alloc( nn );
		ransampl_set( cube, xcube );

		// Normalize grid probabilities to one
		sum=0.;
		for(int i=0;i<nn;i++){
			sum+=xcube[i];
		}
		for(int i=0;i<nn;i++){
			xcube[i]/=sum;
		}


		}


~RandomDraws2() {
		ransampl_free( ws );
		ransampl_free( cube );
		free(x);
		free(xcube);
	}

		integer3 random_xidraw(gsl_rng* rng, double* p){
			// Draws the index of a box at some distance which is weighted by the correlation function
			int n=ransampl_draw( ws, gsl_rng_uniform(rng), gsl_rng_uniform(rng) );
			*p=x[n];
			return cubifyindex(nside,n);
		}

		integer3 random_cubedraw(gsl_rng* rng, double* p){
			// Can be used to draw only a subset of the boxes within maxsep
			int n=ransampl_draw( cube, gsl_rng_uniform(rng), gsl_rng_uniform(rng) );
			*p=xcube[n];
			return cubifyindex(nsidecube,n);
		}

		// Undo 1d back to 3-d indexing
		integer3 cubifyindex(int nside,int n){
            
			assert(n>=0&&n<pow(nside,3));

			integer3 cid;
			cid.z = n%nside-((nside-1)/2);
			n = n/nside;
			cid.y = n%nside-((nside-1)/2);
			cid.x = n/nside-((nside-1)/2);
		
			return cid;
		}

	private:
		void readData(double **x,long *np,const char *filename){

			char line[10000];
			int  n=0;
			FILE *fp;
			fp = fopen(filename, "r");
			if (fp==NULL) {
				fprintf(stderr,"File %s not found\n", filename); abort();
			}

			//Count lines
			while (fgets(line,10000,fp)!=NULL) {
			    if (line[0]=='#') continue;
				if (line[0]=='\n') continue;
				n++;
			}
			rewind(fp);
			*np = n;

			if(n!=pow(nside,3)){
				fprintf(stderr,"File %s does not contain the correct probability grid.\n", filename); abort();
			}


			//Read content of lines and columns
			*x = (double *)malloc(sizeof(double)*n);
			printf("# Found %d lines in %s\n", n, filename);


			int ct=0;
			while (fgets(line,10000,fp)!=NULL) {
				if (line[0]=='#') continue;
				if (line[0]=='\n') continue;
				sscanf(line, "%lf",  &((*x)[ct++]) );
			}

			assert(ct==n);

			fclose(fp);

		}
		
		void integData2(double **x, long *np, integrand xi_fun, int nside, double boxside){
            // Implements the 1d integration of the xi_integrand over all possible distances between 
            // points in two boxes in a grid (one of them being the central box)
            // of nside*nside*nside with the boxes having the sidelength boxside
            
            // Number of points in grid
            (*np)=(int)pow(nside,3);
            
            // Array to house probabilities
			*x = (double *)malloc(sizeof(double)*(*np));

			printf("\nNumber of Boxes in Probability Grid: %ld\n",(*np));
			fflush(NULL);

#pragma omp parallel
        {            
            int len=(nside-1)/2; // This works because nside has been required to be odd
            Float R = boxside/2, n;

            // Define integration limits (using some large upper limit)
            double xmin[1]={0}, xmax[1]={2*boxside*len}, val, err, param[2]={R,0};
            
#ifdef OPENMP
#pragma omp for schedule(dynamic,32)
#endif
            for(int i = 0; i<=len;i++){
                for(int k = 0; k<=len;k++){
                    for(int l = 0;l<=len;l++){
                        //NB: This for-loop recomputes a lot of values but its speed makes this irrelevant
                        n = sqrt(pow(i,2)+pow(k,2)+pow(l,2))*boxside; // distance from origin
                        param[1]=n;
                        
                        hcubature(1, xi_fun, &param[0], 1, xmin, xmax, 0, 0, 1e-5, ERROR_INDIVIDUAL, &val, &err);
 
                        // Insert the calculated value to its grid position
                        for(int is=-1;is<=1;is+=2){
                            for(int ks=-1;ks<=1;ks+=2){
                                for(int ls=-1;ls<=1;ls+=2){
                                    int ic=is*i+len;
                                    int kc=ks*k+len;
                                    int lc=ls*l+len;
                                    (*x)[nside*nside*ic+nside*kc+lc]=val;
                                }
                            }
                        }
                    }
                }
            }
        }
        }
		

		void integData(double **x, long *np, integrand fun, int nside, double boxside, double power){
			// Implements the integration of a function fun over all possible distance between
			// points in two boxes in a grid (one of them being the central box)
			// of nside*nside*nside and with the boxes having sidelength boxside

			(*np)=(int)pow(nside,3);

			*x = (double *)malloc(sizeof(double)*(*np));

			printf("\nNumber of Boxes in Probability Grid: %ld\n",(*np));
			fflush(NULL);

//			TODO: assure it compiles without OPENMP defed
			// parallel environment before for loop to have local integrator variables
			#pragma omp	parallel
			{
			int len=(nside-1)/2; // This works because nside has been required to be odd

			// Parameters for hcubature basically integration limits for the three x dimensions
			// and parameters for the integrand in n
			double xmin[3] = {0,0,0}, xmax[3] = {1,1,1}, n[5] = {0,0,0,boxside/2., power}, val, err;
#ifdef OPENMP
			#pragma omp for schedule(dynamic,32)
#endif
			// Iterates over the grid
				for(int i=0;i<=len;i++){
					for(int k=0;k<=len;k++){
						for(int l=0;l<=len;l++){
							n[0]=i;
							n[1]=k;
							n[2]=l;
							hcubature(1, fun, &n[0], 3, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &val, &err);

							// Insert the calculated value to its grid position
							for(int is=-1;is<=1;is+=2){
								for(int ks=-1;ks<=1;ks+=2){
									for(int ls=-1;ls<=1;ls+=2){
										int ic=is*i+len;
										int kc=ks*k+len;
										int lc=ls*l+len;
                                        
                                        //printf("Just using 1/r^2 probabilities of centers of cells for both samplers\n");
                                        (*x)[nside*nside*ic+nside*kc+lc]=val;//1./(0.1+pow(i,2.)+pow(k,2.)+pow(l,2.));//val;
		//								if(i==len&&k==len)
		//									printf("%d %f, ",side*side*ic+side*kc+lc,y[side*side*ic+side*kc+lc]/val-1);
									}
								}
							}
		//					if(i==len&&k==len)
		//						printf("\n");



		//	    			hcubature_v(1, f_v, &n[0], 3, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &val, &err);
		//		    			printf("Computed integral = %0.10g +/- %g\n", val, err);

						}
					}
				}
			}
		}
		void copyData(double **x, long *n,const double *xin){
			// Copy probability grid
			(*x) = (double *)malloc(sizeof(double)*(*n));
			for(int i=0;i<(*n);i++){
				(*x)[i]=xin[i];
			}
		}

		void readDataBin(double **x, long *n, const char *filename){

			FILE *fp;
			int _nside;
			double _boxside;
			int stat=0;
			fp = fopen(filename, "rb");
			if (fp==NULL) {
				fprintf(stderr,"# File %s not found\n", filename); return;
			}
			stat+=fread(&_nside, sizeof(int), 1, fp);
			stat+=fread(&_boxside, sizeof(double), 1, fp);
			if(!(_nside==nside&&_boxside==boxside)){
				fprintf(stderr,"# Size has changed. Recalculating probability grid.\n");
				fflush(NULL);
				*n=0;
				fclose(fp);
				return;
			}
			stat+=fread(n, sizeof(long), 1, fp);
			*x = (double *)malloc(sizeof(double)*(*n));
			if(*x==NULL){
				fprintf(stderr,"Allocation error.\n");
				fflush(NULL);
				abort();
			}

			if(stat!=3||*n!=(long)fread(*x, sizeof(double), *n, fp)){
	        	fprintf(stderr,"# Error reading file %s. Recalculating probability grid.\n", filename);
	        	fflush(NULL);
	        	*n=0;
	        	free(*x);
	        	fclose(fp);
	        	return;
	        }

	        assert(*n==pow(nside,3));

	        fclose(fp);

		}



		void writeDataBin(double **x, long *n, const char *filename){

			FILE *fp;
			fp = fopen(filename, "wb");
			if (fp==NULL) {
				fprintf(stderr,"# File %s not found\n", filename);
				return;
			}
			fwrite(&nside, 1, sizeof(int), fp);
			fwrite(&boxside, 1, sizeof(double), fp);
			fwrite(n, 1, sizeof(long), fp);
	        if(*n*sizeof(double)!=fwrite(*x, 1, *n*sizeof(double), fp)){
	        	fprintf(stderr,"# Error writing to file %s.\n", filename);
	        	return;
	        }
	        printf("# Probability grid written to file %s.\n", filename);

	        fclose(fp);

		}

	public:
        
    static int xi_integrand(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
        // Integrand for the 1d probability grid integral weighted by xi(r) (in correct cubature format).
        // NB: This must be multiplied by a prefactor independent of the integration variable x.
        // This uses the transform x = u/1-u to transform to a finite integration range
        double* param = (double *) fdata;
        
        // Read in parameters
        const double R = param[0];
        const double n = param[1];
        
        // Compute integrand
        Float factor_1 = pow((x[0]+n)/(2*R),2);
        Float factor_2 = pow((x[0]-n)/(2*R),2);
        
        fval[0] = x[0]*corr->xi(x[0]) / (R*n) * (exp(-factor_1)+exp(-factor_2));
        return 0;
    }
    
    Float xi_integrand_eval(double *x, const double R, const double n){
        Float factor_1 = pow((x[0]+n)/(2*R),2);
        Float factor_2 = pow((x[0]-n)/(2*R),2);
        Float fval = x[0]*corr->xi(x[0]) / (R*n) * (exp(-factor_1)+exp(-factor_2));
        return fval;
    }
        
    Float r2prob(const Float n, const Float a){
        // Return the (unnormalized) 1d probability for picking a cell at distance n from the center of the grid
        // This uses a 1/r^2 kernel giving a semi-analytic form involving a Dawson function.
        // a is the cubic-boxsize. (Here used as a spherical gaussian window)
        return gsl_sf_dawson(n/a)/(a*n);
    }
    
	static int f_corr(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
	    // integrand using the correlation function to weight the cells as done in Rascal
		double* n = (double *) fdata;
	    double sum = 0;
	    int i,l,k;

	    for(i=-1;i<=1;i+=2){
	    	for(k=-1;k<=1;k+=2){
	    		for(l=-1;l<=1;l+=2){
	    		    sum+=corr->xi(2*n[3]*sqrt((pow(x[0] + i*n[0],2) + pow(x[1] + k*n[1],2) + pow(x[2] + l*n[2],2))));
	    		}
	    	}
	    }
	    fval[0] = fabs(pow(2*n[3],6)*(1 - x[0])*(1 - x[1])*(1 - x[2])*sum);
		return 0;
	}

	static int f(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
		// integrand implementing 1/r^2n where n can be an arbitrary double
	    double* n = (double *) fdata;
	    double sum = 0;
	    int i,l,k;

	    for(i=-1;i<=1;i+=2){
	    	for(k=-1;k<=1;k+=2){
	    		for(l=-1;l<=1;l+=2){
	    		    sum+=1./(1./pow(2*n[3],2.*n[4]) + pow(pow(x[0] + i*n[0],2) + pow(x[1] + k*n[1],2) + pow(x[2] + l*n[2],2),n[4]));
	    		}
	    	}
	    }
	    fval[0] = pow(2*n[3],6.-2.*n[4])*(1 - x[0])*(1 - x[1])*(1 - x[2])*sum;
	    return 0;
	}

	static int f_cube(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
	    // cubic integrand
		double* n = (double *) fdata;
	    double sum = 0;
	    int i,l,k;

	    for(i=-1;i<=1;i+=2){
	    	for(k=-1;k<=1;k+=2){
	    		for(l=-1;l<=1;l+=2){
	    		    sum+=1./(1./pow(2*n[3],6) + pow(pow(x[0] + i*n[0],2) + pow(x[1] + k*n[1],2) + pow(x[2] + l*n[2],2),3));
	    		}
	    	}
	    }
	    fval[0] = pow(2*n[3],0)*(1 - x[0])*(1 - x[1])*(1 - x[2])*sum;
	    return 0;
	}


	// Test of different parallelization of the integration used for setting up the probability grid
	// However, no significant performance gain observed
	// Also, the setup of the grid already takes only a very small fraction of the full calculation
	// so further optimization not urgent

	#pragma omp declare simd uniform(n,x,ndim) linear(j:1)
	double eval(const double* n,const double * x, unsigned ndim,unsigned j){
	    int i,l,k;
		double sum = 0;
			for(i=-1;i<=1;i+=2){
				for(k=-1;k<=1;k+=2){
					for(l=-1;l<=1;l+=2){
						sum+=1./(1./pow(2*n[3],4) + pow(pow(x[j*ndim+0] + i*n[0],2) + pow(x[j*ndim+1] + k*n[1],2) + pow(x[j*ndim+2] + l*n[2],2),2));
					}
				}
			}
			return pow(2*n[3],2)*(1 - x[j*ndim+0])*(1 - x[j*ndim+1])*(1 - x[j*ndim+2])*sum;
	}


	int f_v(unsigned ndim, unsigned long npts, const double *x, void *fdata, unsigned fdim, double *fval) {

	    double* n = (double *) fdata; // we can pass σ via fdata argument
	    unsigned j;

	#pragma omp simd
	    for (j = 0; j < npts; ++j) {
		    fval[j] = eval(n,x,ndim,j);
	    }
	    return 0; // success
	}


};


#endif
