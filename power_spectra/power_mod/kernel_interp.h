// kernel_interp.h - Interpolator for kernel functions for power spectrum computation

#ifndef KERNEL_INTERP_H
#define KERNEL_INTERP_H

class KernelInterp{
    // Compute a interpolation kernel function
private:
    Float R0, k_min, k_max;

public:
    int npoint=100000;
    Float min_val,max_val;
    double *sep, *kernel_vals_0, *kernel_vals_2, *kernel_vals_4, *kernel_vals_6;
    gsl_interp_accel *sep_a;
    gsl_spline *interp_kernel_0, *interp_kernel_2, *interp_kernel_4,*interp_kernel_6;

public:
    inline double kernel(int ell, double this_sep){
        // Kernel function accelerated by interpolation
        if(this_sep<=min_val) this_sep=min_val;
        if (ell==0) return gsl_spline_eval(interp_kernel_0,this_sep,sep_a);
        else if (ell==2) return gsl_spline_eval(interp_kernel_2,this_sep,sep_a);
        else if (ell==4) return gsl_spline_eval(interp_kernel_4,this_sep,sep_a);
        else if (ell==6) return gsl_spline_eval(interp_kernel_6,this_sep,sep_a);
        else{
            printf("Multipoles greater than ell = 6 not yet implemented");
            exit(1);
        }
    }

public:
    void copy_function(KernelInterp *kern){
        // Copy pre-existing kernel function into object
        npoint = kern->npoint;
        min_val = kern->min_val;

        // Allocate memory
        sep = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_0 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_2 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_4 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_6 = (double *)malloc(sizeof(double)*npoint);

        // Read in grids
        for(int i=0;i<npoint;i++){
            sep[i]=kern->sep[i];
            kernel_vals_0[i]=kern->kernel_vals_0[i];
            kernel_vals_2[i]=kern->kernel_vals_2[i];
            kernel_vals_4[i]=kern->kernel_vals_4[i];
            kernel_vals_6[i]=kern->kernel_vals_6[i];
        }
        // activate interpolator
        interpolate();
    }

private:
    void copy(int& _npoint, Float& _min_val, double **_sep, double **_kernel_vals_0, double **_kernel_vals_2, double **_kernel_vals_4, double **_kernel_vals_6){
        _npoint = npoint;
        _min_val = min_val;
        *_sep = (double *)malloc(sizeof(double)*npoint);
        *_kernel_vals_0 = (double *)malloc(sizeof(double)*npoint);
        *_kernel_vals_2 = (double *)malloc(sizeof(double)*npoint);
        *_kernel_vals_4 = (double *)malloc(sizeof(double)*npoint);
        *_kernel_vals_6 = (double *)malloc(sizeof(double)*npoint);

        for(int i=0;i<npoint;i++){
            (*_sep)[i]=sep[i];
            (*_kernel_vals_0)[i]=kernel_vals_0[i];
            (*_kernel_vals_2)[i]=kernel_vals_2[i];
            (*_kernel_vals_4)[i]=kernel_vals_4[i];
            (*_kernel_vals_6)[i]=kernel_vals_6[i];
        }
    }


    void interpolate(){
        interp_kernel_0 = gsl_spline_alloc(gsl_interp_cspline,npoint);
        interp_kernel_2 = gsl_spline_alloc(gsl_interp_cspline,npoint);
        interp_kernel_4 = gsl_spline_alloc(gsl_interp_cspline,npoint);
        interp_kernel_6 = gsl_spline_alloc(gsl_interp_cspline,npoint);
        gsl_spline_init(interp_kernel_0,sep,kernel_vals_0,npoint);
        gsl_spline_init(interp_kernel_2,sep,kernel_vals_2,npoint);
        gsl_spline_init(interp_kernel_4,sep,kernel_vals_4,npoint);
        gsl_spline_init(interp_kernel_6,sep,kernel_vals_6,npoint);
        sep_a = gsl_interp_accel_alloc();
    }

public:
    KernelInterp(){
        // Empty constructor
    }

    KernelInterp(KernelInterp *kern){
        // Copy constructor

        kern->copy(npoint,min_val,&sep,&kernel_vals_0, &kernel_vals_2, &kernel_vals_4, &kernel_vals_6);
        interpolate();
    }

    KernelInterp(Parameters *par){

        R0 = par->R0;
        k_min = par->kmin;
        k_max = par->kmax;

        // Construct interpolation object
        Float tmp_kw,Si_int,tmp_bessel, tmp_sin, tmp_cos, tmp_kernel;

        // Allocate memory
        sep = (double *)malloc(sizeof(double)*npoint); // k*|r_i-r_j|
        kernel_vals_0 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_2 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_4 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_6 = (double *)malloc(sizeof(double)*npoint);

        // This minimum value should be sufficient for multipoles up to ell = 10
        // Going to too-small ell gives nasty divergences at high ell
        min_val = 0.01*(R0*k_min)/double(npoint); // minimum interpolation value
        max_val = 2.01*(R0*k_max); // maximum interpolation value

        for(int i=0;i<npoint;i++){
            tmp_kw = min_val + double(i)/double(npoint-1)*(max_val-min_val);
            sep[i] = tmp_kw;
            tmp_sin = sin(tmp_kw);
            tmp_cos = cos(tmp_kw);

            Si_int = gsl_sf_Si(tmp_kw);

            // Insert kw = 0 limit for low k to avoid blow-up
            if(tmp_kw<0.01){
                for(int ell=0;ell<=6;ell+=2){
                    if(ell==0) kernel_vals_0[i] = 0;
                    else if(ell==2) kernel_vals_2[i] = 0;
                    else if(ell==4) kernel_vals_4[i] = 0;
                    else if(ell==6) kernel_vals_6[i] = 0;
                    else{
                        printf("\nOnly ell = 0,2,4,6 implemented\n");
                        exit(1);
                    }
                }
            }
            else{
                // Compute multipole contributions
                for(int ell=0;ell<=6;ell+=2){
                    tmp_kernel = 3.*Float(pow(-1,ell/2)*(2*ell+1.));
                    if(ell==0){
                        tmp_bessel = tmp_sin - tmp_kw*tmp_cos;
                        kernel_vals_0[i] = tmp_bessel*tmp_kernel;
                    }
                    else if(ell==2){
                        tmp_bessel = tmp_kw*tmp_cos - 4*tmp_sin+3*Si_int;
                        kernel_vals_2[i] = tmp_bessel*tmp_kernel;
                    }
                    else if(ell==4){
                        tmp_bessel = 0.5*((105./tmp_kw-2.*tmp_kw)*tmp_cos+(22.-105./pow(tmp_kw,2))*tmp_sin+15*Si_int);
                        kernel_vals_4[i] = tmp_bessel*tmp_kernel;
                    }
                    else if(ell==6){
                        tmp_bessel = 1./(8.*pow(tmp_kw,4))*(tmp_kw*(20790. - 1575.*pow(tmp_kw,2) + 8.*pow(tmp_kw,4.))*tmp_cos + (-20790. + 8505.*pow(tmp_kw,2) - 176.*pow(tmp_kw,4.))*tmp_sin + 105.*pow(tmp_kw,4.)*Si_int);
                        kernel_vals_6[i] = tmp_bessel*tmp_kernel;
                    }

                    else{
                        printf("\nOnly ell = 0,2,4,6 implemented\n");
                        exit(1);
                    }
                }
            }
        }

        // activate interpolator function
        interpolate();
    }

    ~KernelInterp() {
        //Destructor
        gsl_spline_free(interp_kernel_0);
        gsl_spline_free(interp_kernel_2);
        gsl_spline_free(interp_kernel_4);
        gsl_spline_free(interp_kernel_6);
        gsl_interp_accel_free(sep_a);
        free(sep);
        free(kernel_vals_0);
        free(kernel_vals_2);
        free(kernel_vals_4);
        free(kernel_vals_6);
    }
};

#endif
