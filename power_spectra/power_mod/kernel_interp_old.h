// kernel_interp.h - Interpolator for kernel functions for power spectrum computation

#ifndef KERNEL_INTERP_H
#define KERNEL_INTERP_H

class KernelInterp{
    // Compute a interpolation kernel function
    
public:
    int npoint=10000;
    double *sep, *kernel_vals;
    gsl_interp_accel *sep_a;
    gsl_spline *interp_kernel;
    
private:
    inline Float pair_weight(Float sep, Float R0){
        // Compute weight function W(r;R_0)
        if(sep<R0/2) return 1.;
        else{
            Float x = sep/R0;
            if(sep<3*R0/4) return 1.-8.*pow(2*x-1,3)+8.*pow(2*x-1,4);
            else return -64.*pow(x-1,3)-32*pow(x-1,4);
        }
    }
    
public:
    inline double kernel(double this_sep){
        // Kernel function accelerated by interpolation
        return gsl_spline_eval(interp_kernel,this_sep,sep_a);
    }
    
public:
    void copy_function(KernelInterp *kern){
        // Copy pre-existing kernel function into object
        npoint = kern->npoint;
        
        // Allocate memory
        sep = (double *)malloc(sizeof(double)*npoint);
        kernel_vals = (double *)malloc(sizeof(double)*npoint);
        
        // Read in grids
        for(int i=0;i<npoint;i++){
            sep[i]=kern->sep[i];
            kernel_vals[i]=kern->kernel_vals[i];
        }
        // activate interpolator
        interpolate();
    }
    
private:
    void copy(int& _npoint, double **_sep, double **_kernel_vals){
        _npoint = npoint;
        *_sep = (double *)malloc(sizeof(double)*npoint);
        *_kernel_vals = (double *)malloc(sizeof(double)*npoint);
        
        for(int i=0;i<npoint;i++){
            (*_sep)[i]=sep[i];
            (*_kernel_vals)[i]=kernel_vals[i];
        }
    }
    
    
    void interpolate(){
        interp_kernel = gsl_spline_alloc(gsl_interp_cspline,npoint);
        gsl_spline_init(interp_kernel,sep,kernel_vals,npoint);
        sep_a = gsl_interp_accel_alloc();
    }
        
public:
    KernelInterp(){
        // Empty constructor
    }
    
    KernelInterp(KernelInterp *kern){
        // Copy constructor
        
        kern->copy(npoint,&sep,&kernel_vals);
        interpolate();
    }
    
    KernelInterp(int ell, Float kw, Float R0){
        // Construct interpolation object
        
        Float tmp_kernel,tmp_sep,kr_lo,kr_hi,k_sin_lo,k_cos_lo,k_sin_hi,k_cos_hi,k_Si_lo=0,k_Si_hi=0,tmp_bessel;
        
        // Allocate memory
        sep = (double *)malloc(sizeof(double)*npoint);
        kernel_vals = (double *)malloc(sizeof(double)*npoint);
        
        tmp_kernel = 3.*Float(pow(-1,ell/2)*(2.*ell+1))/(pow(k2,3)-pow(k1,3));
        
        // Compute x,z grids
        for(int i=0;i<npoint;i++){
            tmp_sep = double(i+1)/double(npoint+1)*R0; // must be greater than zero here
            sep[i] = tmp_sep;
            
            kr_lo = tmp_sep*k1;
            kr_hi = tmp_sep*k2;
            
            k_sin_lo = sin(kr_lo);
            k_sin_hi = sin(kr_hi);
            k_cos_lo = cos(kr_lo);
            k_cos_hi = cos(kr_hi);
            
            // Define sine integrals
            if(ell>1){
                k_Si_lo = gsl_sf_Si(kr_lo);
                k_Si_hi = gsl_sf_Si(kr_hi);
            }
            
            // Now compute multipole contributions
            if(ell==0) tmp_bessel = (k_sin_hi-kr_hi*k_cos_hi)-(k_sin_lo-kr_lo*k_cos_lo);
            else if(ell==2) tmp_bessel = (kr_hi*k_cos_hi - 4.*k_sin_hi+3.*k_Si_hi) - (kr_lo*k_cos_lo - 4.*k_sin_lo + 3.*k_Si_lo);
            else if(ell==4) tmp_bessel = 0.5*(((105/kr_hi-2.*kr_hi)*k_cos_hi + (22.-105./pow(kr_hi,2))*k_sin_hi + 15.*k_Si_hi)-((105/kr_lo-2.*kr_lo)*k_cos_lo+(22.-105./pow(kr_lo,2))*k_sin_lo + 15.*k_Si_lo));
            else{
                printf("\nOnly ell = 0,2,4 implemented\n");
                exit(1);
            }
            
            kernel_vals[i] = pair_weight(tmp_sep,R0)/pow(tmp_sep,3)*tmp_bessel*tmp_kernel;
        }
        
        // activate interpolator function
        interpolate();
    }
    
    ~KernelInterp() {
        //Destructor
        gsl_spline_free(interp_kernel);
        gsl_interp_accel_free(sep_a);
        free(sep);
        free(kernel_vals);
    }
};

#endif
        
