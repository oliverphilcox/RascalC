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
        bool mudim = 0;
        gsl_interp_accel *xa, *ya, *x1a;
        gsl_interp2d* interp_2d;
        gsl_spline* corfu1d;
        bool interp_setup = 0;

        double smooth_transition(double x, double x1, double x2, double f1, double f2) {
            // Smooth transition for f(x) given f(x1) = f1 and f(x2) = f2, assuming x2 > x1
            if (x >= x2) return f2;
            if (x <= x1) return f1;
            double y = (x2-x)/(x2-x1); // rescaled variable that ranges from 0 for x=x2 to 1 for x=x1
            double y2 = y*y;
            double y3 = y2*y;
            double weight1 = y3 * (10. - 15. * y + 6. * y2); // this polynomial is 0 at y=0 and 1 at y=1, and has zero 1st and 2nd derivatives at both points. Thus this transition is continuous, smooth and has continuous 2nd derivative, just like the (bi)cubic spline used below. Moreover, the polynomial values are contained in the [0, 1] interval for 0<=y<=1
            double weight2 = 1. - weight1;
            return weight1 * f1 + weight2 * f2;
        }

    public:
        double xi(double r, double mu){
            // 2D correlation function in radius and angular bins
            // xi values beyond the maximal radius in the correlation function file read in are extrapolated
            // as a simple r^-4 power law for each mu bin independently (may lead to different signs for
            // neighbouring mu bins if the correlation function is noisy and fluctuates around 0)
            if(mudim){
                assert(mu>=0);
                if(r>rmax){
                    double tmu = fmin(fmax(mu,mumin),mumax);
                    // Choose r-max to ensure function integrates to zero
                    return gsl_interp2d_eval_extrap(interp_2d,y,x,z,tmu,rmax, xa, ya)/pow(rmax,2)*pow(r/rmax,-4.);
                }
                else if(r<rmin){
                    double tmu = fmin(fmax(mu,mumin),mumax);
                    return smooth_transition(r, rmin/2., rmin, gsl_interp2d_eval_extrap(interp_2d, y, x, z, tmu, rmin/2., xa, ya)/pow(rmin/2., 2), gsl_interp2d_eval_extrap(interp_2d, y, x, z, tmu, r, xa, ya)/pow(r, 2)); // this ensures xi goes as constant below rmin/2 (rmin is the middle of smallest bin) and not as 1/r which might cause problems, and the transition is smooth
                }
                if(mu<mumin){
                    double tr = fmin(fmax(r,rmin),rmax);
                    return gsl_interp2d_eval_extrap(interp_2d,y,x,z,mumin,tr, xa, ya)/pow(tr,2);
                }
                else if(mu>mumax){
                    double tr = fmin(fmax(r,rmin),rmax);
                    return gsl_interp2d_eval_extrap(interp_2d,y,x,z,mumax,tr, xa, ya)/pow(tr,2);
                    }
                else{
                    return gsl_interp2d_eval_extrap(interp_2d,y,x,z,mu,r, xa, ya)/pow(r,2);
                }
            }
            else
                return xi(r);
        }
        double xi(double r){
            // Radial correlation function
            // xi values beyond the maximal radius in the correlation function file read in are extrapolated as a simple 1/r^4 power law
            if(r>rmax){
                return gsl_spline_eval(corfu1d,rmax, x1a)/pow(rmax,2)*pow(r/rmax,-4.);
            }
            else{
                if(r<rmin){
                    return smooth_transition(r, rmin/2., rmin, gsl_spline_eval(corfu1d, rmin/2., x1a)/pow(rmin/2., 2), gsl_spline_eval(corfu1d, r, x1a)/pow(r, 2)); // this ensures xi goes as some constant below rmin/2 (rmin is the middle of smallest bin) and not as 1/r which might cause problems, and the transition is smooth
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
            while (pch != NULL)
            {
            pch = strtok (NULL, " \t\n");
            if(pch != NULL){
                sscanf(pch, "%lf",  &((*z)[(*ct)++]) );
            }
            }
        }

        void readData(const char *filename,double **x,double **y, double **z,int *np,int *mp){

            char line[100000];
            char * pch;
            int n=0,m=0;
            double x0;
            FILE *fp;
            fp = fopen(filename, "r");
            if (fp==NULL) {
                fprintf(stderr,"File %s not found\n", filename); abort();
            }

            //Count lines and columns
            while (fgets(line,100000,fp)!=NULL) {
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

            printf("# Found %d radial and %d mu bins in %s\n", n-1,m, filename);

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
            while (fgets(line,100000,fp)!=NULL) {
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

        }
    public:
        void copy_function(CorrelationFunction *cf){
        // Copy a preexisting correlation function into this object
            xsize=cf->xsize;
            ysize=cf->ysize;
            rmin=cf->rmin;
            rmax=cf->rmax;
            mumin=cf->mumin;
            mumax=cf->mumax;
            mudim=cf->mudim;

            // Allocate memory
            x = (double *)malloc(sizeof(double)*xsize);
            y = (double *)malloc(sizeof(double)*ysize);
            z = (double *)malloc(sizeof(double)*xsize*ysize);

            // Read in x,y,z grids
            int ct=0;
            for(int i=0;i<xsize;i++){
                x[i]=cf->x[i];
                for(int j=0;j<ysize;j++){
                    z[ct]=cf->z[ct];
                    ct++;
                }
            }
            for(int j=0;j<ysize;j++) y[j]=cf->y[j];
            // activate the interpolator function here
            interpolate();
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
                        y1[i]+=z[ct++];
                        col++;
                }
                y1[i]/=col;
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
            interp_setup = 1;
        }
    public:
        CorrelationFunction(){
            // empty constructor
        }
        CorrelationFunction(CorrelationFunction* corr){
            // Copy constructor

            corr->copy(xsize,ysize,&x,&y,&z,rmin,rmax,mumin,mumax,mudim);
            interpolate();
        }
    CorrelationFunction(const char *filename, int nbin, Float *r_low, Float *r_high, int mbin, Float dmu){
        // Construct from input file

        readData(filename,&x,&y,&z,&xsize,&ysize);

        if (xsize != nbin+1) {
          fprintf(stderr,"%d r-bins found in correlation function file but %d specified in parameters.\n", xsize-1, nbin);
          abort();
        }

        for (int i = 0; i < nbin; i++) {
            if ((x[i+1] < r_low[i]) || (x[i+1] > r_high[i])) {
                fprintf(stderr,"%d'th r-bin found in correlation function file is %le but expected between %le and %le from binning file.\n", i, x[i+1], r_low[i], r_high[i]);
                abort();
            }
        }

        if(ysize!=mbin){
          fprintf(stderr,"%d mu-bins found in correlation function file but %d specified in parameters.\n",ysize,mbin);
          abort();
        }

        if(ysize<=1&&mbin>1){
            fprintf(stderr,"Requested %d mu-bins but correlation function file only has %d columns.\n",mbin,ysize);
            abort();
        }

        mudim=!(mbin==1&&dmu==1.);
        rmin=x[1];
        rmax=x[xsize-1];
        // Mu range hardcoded to be between 0 and 1
        mumin=0.;
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
    CorrelationFunction(Float* xi_array, Float * r_array, Float* mu_array, int nbin, int mbin){
        // Construct from input array

        xsize=nbin+1;
        ysize=mbin;
        mumin=mu_array[0];
        mumax=mu_array[mbin-1];
        mudim=mu_array[1]-mu_array[0];

        // Allocate memory
        x = (double *)malloc(sizeof(double)*xsize);
        y = (double *)malloc(sizeof(double)*ysize);
        z = (double *)malloc(sizeof(double)*xsize*ysize);

        // Read in x,y,z grids
        x[0]=0;
        for(int j=0;j<ysize;j++) z[j] = 0.; // set initial value to zero (since we multiply by r^2 anyway)
        int ct=0;
        for(int i=1;i<xsize;i++){
            x[i]=r_array[i-1];
            for(int j=0;j<ysize;j++){
                z[ct+ysize]=xi_array[ct]*pow(x[i],2.); // also multiply by r^2
                ct++;
            }
        }
        assert(ct==nbin*mbin);

        // Define r ranges (excluding zero)
        rmin=x[1];
        rmax=x[xsize-1];

        for(int j=0;j<ysize;j++) y[j]=mu_array[j];

        // activate the interpolator function here
        interpolate();
    }

    ~CorrelationFunction() {
        // Destructor

        if (mudim) {
            gsl_interp2d_free(interp_2d);
            gsl_interp_accel_free(xa);
            gsl_interp_accel_free(ya);
        }
        if (interp_setup) {
            gsl_spline_free(corfu1d);
            gsl_interp_accel_free(x1a);
        }
    }

};
#endif
