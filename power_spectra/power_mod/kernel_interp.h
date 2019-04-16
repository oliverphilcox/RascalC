// kernel_interp.h - Interpolator for kernel functions for power spectrum computation

#ifndef KERNEL_INTERP_H
#define KERNEL_INTERP_H

class KernelInterp{
    // Compute a interpolation kernel function
    
public:
    int npoint=100000;
    Float min_val;
    double *sep, *kernel_vals_0, *kernel_vals_2, *kernel_vals_4;
    gsl_interp_accel *sep_a;
    gsl_spline *interp_kernel_0, *interp_kernel_2, *interp_kernel_4;
    
public:
    inline double kernel(int ell, double this_sep){
        // Kernel function accelerated by interpolation
        if(this_sep<=min_val) return 0.;
        else{
            if (ell==0) return gsl_spline_eval(interp_kernel_0,this_sep,sep_a);
            else if (ell==2) return gsl_spline_eval(interp_kernel_2,this_sep,sep_a);
            else if (ell==4) return gsl_spline_eval(interp_kernel_4,this_sep,sep_a);
            else{
                printf("Multipoles greater than hexadecapole not yet implemented");
                exit(1);
            }
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
        
        // Read in grids
        for(int i=0;i<npoint;i++){
            sep[i]=kern->sep[i];
            kernel_vals_0[i]=kern->kernel_vals_0[i];
            kernel_vals_2[i]=kern->kernel_vals_2[i];
            kernel_vals_4[i]=kern->kernel_vals_4[i];
        }
        // activate interpolator
        interpolate();
    }
    
private:
    void copy(int& _npoint, Float& _min_val, double **_sep, double **_kernel_vals_0, double **_kernel_vals_2, double **_kernel_vals_4){
        _npoint = npoint;
        _min_val = min_val;
        *_sep = (double *)malloc(sizeof(double)*npoint);
        *_kernel_vals_0 = (double *)malloc(sizeof(double)*npoint);
        *_kernel_vals_2 = (double *)malloc(sizeof(double)*npoint);
        *_kernel_vals_4 = (double *)malloc(sizeof(double)*npoint);
        
        for(int i=0;i<npoint;i++){
            (*_sep)[i]=sep[i];
            (*_kernel_vals_0)[i]=kernel_vals_0[i];
            (*_kernel_vals_2)[i]=kernel_vals_2[i];
            (*_kernel_vals_4)[i]=kernel_vals_4[i];
        }
    }
    
    
    void interpolate(){
        interp_kernel_0 = gsl_spline_alloc(gsl_interp_cspline,npoint);
        interp_kernel_2 = gsl_spline_alloc(gsl_interp_cspline,npoint);
        interp_kernel_4 = gsl_spline_alloc(gsl_interp_cspline,npoint);
        gsl_spline_init(interp_kernel_0,sep,kernel_vals_0,npoint);
        gsl_spline_init(interp_kernel_2,sep,kernel_vals_2,npoint);
        gsl_spline_init(interp_kernel_4,sep,kernel_vals_4,npoint);
        sep_a = gsl_interp_accel_alloc();
    }
        
public:
    KernelInterp(){
        // Empty constructor
    }
    
    KernelInterp(KernelInterp *kern){
        // Copy constructor
        
        kern->copy(npoint,min_val,&sep,&kernel_vals_0, &kernel_vals_2, &kernel_vals_4);
        interpolate();
    }
    
    KernelInterp(Float R0, Float k_min){
        // Construct interpolation object
        Float tmp_kw,Si_int,tmp_bessel, tmp_sin, tmp_cos, tmp_kernel;
        
        // Allocate memory
        sep = (double *)malloc(sizeof(double)*npoint); // k*|r_i-r_j|
        kernel_vals_0 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_2 = (double *)malloc(sizeof(double)*npoint);
        kernel_vals_4 = (double *)malloc(sizeof(double)*npoint);
        
        min_val = (R0*k_min)/double(npoint); // minimum interpolation value
        
        for(int i=0;i<npoint;i++){
            tmp_kw = double(i+1)/double(npoint)*(R0*k_min); // must be greater than zero here
            sep[i] = tmp_kw;
            tmp_sin = sin(tmp_kw);
            tmp_cos = cos(tmp_kw);
            
            Si_int = gsl_sf_Si(tmp_kw);
            
            // Compute multipole contributions
            for(int ell=0;ell<=4;ell+=2){
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
                else{
                    printf("\nOnly ell = 0,2,4 implemented\n");
                    exit(1);
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
        gsl_interp_accel_free(sep_a);
        free(sep);
        free(kernel_vals_0);
        free(kernel_vals_2);
        free(kernel_vals_4);
    }
};

#endif
        
