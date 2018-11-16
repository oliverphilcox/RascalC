// correlation function class for grid_covariance.cpp file (originally from Alex Wiegand)

#ifndef CORRELATION_FUNCTION_H
#define CORRELATION_FUNCTION_H

class CorrelationFunction{
	/* Reads and stores a 2d correlation function. Values in between the grid positions of the input
	 * are interpolated using gsl_interp2d
	 */
	private:
		int xsize, ysize;
		double *x,*y,*z;
		double rmin,rmax,mumin,mumax;
		bool mudim;
		gsl_interp_accel *xa, *ya, *x1a;
        
        
		gsl_interp2d* interp_2d;
		gsl_spline* corfu1d;
        
        double r_cutoff;

	public:
		double xi(double r, double mu){
			// 2D correlation function in radius and angular bins
			// xi values beyond the maximal radius in the correlation function file read in are extrapolated
			// as a simple r^-4 power law for each mu bin independently (may lead to different signs for
			// neighbouring mu bins if the correlation function is noisy and fluctuates around 0
			if(mudim){
				if(mu<mumin||mu>mumax||r<rmin||r>rmax){
					double tmu = fmin(fmax(mu,mumin),mumax);
					double tr = fmin(fmax(r,rmin),rmax);
					if(r>rmax){
						return gsl_interp2d_eval_extrap(interp_2d,y,x,z,tmu,tr, xa, ya)/pow(tr,2)*pow(r/rmax,-4);
					}
					else if(r<r_cutoff){
                        return 0.;
                    }
					else
						return gsl_interp2d_eval_extrap(interp_2d,y,x,z,tmu,tr, xa, ya)/pow(tr,2);
				}
				else{
					return gsl_interp2d_eval_extrap(interp_2d,y,x,z,mu,r, xa, ya)/pow(r,2);
				}
			}
			else
				return xi(r);
		}
		double xi2(double r12, double r2) {
			// Simple test correlation function to compare to sobol
			double r0 = 8.0;
		    return r0*r0/(r12*r12+1.0);
		    	// Simple case, r_0 = 8 Mpc, softening at 1 Mpc
		}

		double xi(double r){
			// Radial correlation function
			// xi values beyond the maximal radius in the correlation function file read in are extrapolated
			// as a simple r^-4 power law
			if(r>rmax){
				return gsl_spline_eval(corfu1d,rmax, x1a)/pow(rmax,2)*pow(r/rmax,-4);
			}
			else{
				if(r<rmin){
					return 0.;
				}
				else
					return gsl_spline_eval(corfu1d,r, x1a)/pow(r,2);
			}
		}


	private:

		void readLine(char *line,double **z, int *ct){
			char * pch;
			pch = strtok (line," \t\n");
			sscanf(pch, "%lf",  &((*z)[(*ct)++]) );
//				printf("\n%f ",(*z)[ct-1]);
			while (pch != NULL)
			{
			  pch = strtok (NULL, " \t\n");
			  if(pch != NULL){
				  sscanf(pch, "%lf",  &((*z)[(*ct)++]) );
//				  printf("%f ", (*z)[ct-1]);
			  }
			}
		}

		void readData(const char *filename,double **x,double **y, double **z,int *np,int *mp){

			char line[10000];
			char * pch;
			int n=0,m=0;
			double x0;
			FILE *fp;
			fp = fopen(filename, "r");
			if (fp==NULL) {
				fprintf(stderr,"File %s not found\n", filename); abort();
			}

			//Count lines and columns
			while (fgets(line,10000,fp)!=NULL) {
			    if (line[0]=='#') continue;
				if (line[0]=='\n') continue;
				if(n==0){
					pch = strtok (line," \t\n");
					sscanf(pch, "%lf",  &x0 );
				}
				if(n==2){
					m=0;
					pch = strtok (line," \t\n");
					while (pch != NULL)
					{
					  m++;
					  pch = strtok (NULL, " \t\n");
					}
				}

				n++;
			}
			rewind(fp);
			n-=2; //Remove the first two lines that should contain r and mu bin centers

			// If there is no r=0 entry, add one
			if(x0!=0.){
				n++;
			}

			*np = n;
			*mp = m;

			*x = (double *)malloc(sizeof(double)*n);
			*y = (double *)malloc(sizeof(double)*m);
			*z = (double *)malloc(sizeof(double)*n*m);

			printf("# Found %d radial and %d mu bins in %s\n", n,m, filename);

			int nx=0, ny=0, nz=0, lnnr=0;

			// Add a first entry at r=0 assuming that we shall multiply by r^2 anyway
			if(x0!=0.){
				(*x)[0]=0.;
				nx++;
				for(int i=0;i<m;i++){
					(*z)[i]=0.;
					nz++;
				}
			}

			//Read content of lines and columns
			while (fgets(line,10000,fp)!=NULL) {
				if (line[0]=='#') continue;
				if (line[0]=='\n') continue;
				if(lnnr==0)	readLine(line, x, &nx);
				if(lnnr==1)	readLine(line, y, &ny);
				if(lnnr>1)	readLine(line, z, &nz);
				lnnr++;
			}

			if(n!=nx){
				fprintf(stderr,"Found %d lines but %d r-values in the first line. Aborting.\n", n, nx); abort();
			}
			if(m!=ny){
				fprintf(stderr,"Found %d columns but %d mu-values in the second line. Aborting.\n", m, ny); abort();
			}

//			for(int i=0;i<n*m;i++){
//				printf("%f ",(*z)[i]);
//			}

		}

	private:
		void copy(int& _xsize, int& _ysize, double ** _x, double ** _y, double ** _z, double& _rmin,
				double& _rmax, double& _mumin, double& _mumax, bool& _mudim){
			// Copy correlation function object

			_xsize=xsize;
			_ysize=ysize;
			_rmin=rmin;
			_rmax=rmax;
			_mumin=mumin;
			_mumax=mumax;
			_mudim=mudim;

			*_x = (double *)malloc(sizeof(double)*xsize);
			*_y = (double *)malloc(sizeof(double)*ysize);
			*_z = (double *)malloc(sizeof(double)*xsize*ysize);

			int ct=0;
			for(int i=0;i<xsize;i++){
				(*_x)[i]=x[i];
				for(int j=0;j<ysize;j++){
					(*_z)[ct]=z[ct];
					ct++;
				}
			}

			for(int j=0;j<ysize;j++){
				(*_y)[j]=y[j];
			}
		}

		void interpolate(){
			double *y1;

			y1 = (double *)malloc(sizeof(double)*xsize);

			// Sum to obtain radial correlation function
			int col=0;
			int ct=0;
			for(int i=0;i<xsize;i++){
				y1[i]=0;
				col=0;
				for(int j=0;j<ysize;j++){
	//				if((j==0&&y[0]==0.)||(j==ysize-1&&y[ysize-1]==1.))
	//					ct++;
	//				else{
						y1[i]+=z[ct++];
						col++;
	//				}

				}
				y1[i]/=col;
	//			fprintf(stderr,"%f ",y1[i]);
			}

			if(mudim){
				interp_2d=gsl_interp2d_alloc(gsl_interp2d_bicubic, ysize, xsize);
				gsl_interp2d_init(interp_2d, y, x, z, ysize, xsize);
				xa = gsl_interp_accel_alloc();
				ya = gsl_interp_accel_alloc();
			}

			corfu1d=gsl_spline_alloc(gsl_interp_cspline, xsize);
			gsl_spline_init(corfu1d, x, y1, xsize);
			x1a = gsl_interp_accel_alloc();
		}
	public:

		CorrelationFunction(CorrelationFunction* corr){
			// Copy constructor

			corr->copy(xsize,ysize,&x,&y,&z,rmin,rmax,mumin,mumax,mudim);
			interpolate();
		}

	CorrelationFunction(const char *filename, int mbin, double dmu, double _r_cutoff){
		// Construct from input file
        
        r_cutoff=_r_cutoff; // read in from input

		readData(filename,&x,&y,&z,&xsize,&ysize);

		if(ysize<=1&&mbin>1){
			fprintf(stderr,"Requested %d mu-bins but correlation function file only has %d column.",mbin,ysize);
			abort();
		}

		mudim=!(mbin==1&&dmu==1.);

//		mudim=true;

		rmin=x[0];
		rmax=x[xsize-1];
		// Mu range hardcoded to be between 0 and 1
		mumin=0.0;
		mumax=1.;

		// Multiply correlation function by r^2 for smoother interpolation
		// If one wants to change this also the readData function needs to be changed
		int ct=0;
		for(int i=0;i<xsize;i++){
			for(int j=0;j<ysize;j++){
				z[ct++]*=pow(x[i],2); //Multiply by r^2
			}
		}

		interpolate();

	}

	~CorrelationFunction() {
		// Destructor

		if(mudim){
			gsl_interp2d_free(interp_2d);
			gsl_interp_accel_free(xa);
			gsl_interp_accel_free(ya);
		}
		gsl_spline_free(corfu1d);
		gsl_interp_accel_free(x1a);
	}

};
#endif
