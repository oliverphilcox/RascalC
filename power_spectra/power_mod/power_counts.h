// power_counts.h for the grid_power.cpp code - Oliver Philcox 2019

#ifndef POWER_COUNTS_H
#define POWER_COUNTS_H

class PowerCounts{
    
private:
    int nbin,mbin;
    Float rmin,rmax,R0; //Ranges in r and truncation radius
    Float *r_high, *r_low; // Max and min of each radial bin
    Float *power_counts; // Power counts
    bool box; // to decide if we have a periodic box
    char *out_file, *out_string;
    int max_legendre; // maximum order of Legendre polynomials needed
    SurveyCorrection *sc; // survey correction function
    KernelInterp *kernel_interp;
public:
    uint64 used_pairs; // total number of pairs used
    
    void cleanup_l(Float3 p1,Float3 p2,Float& norm,Float& mu){
        Float3 pos=p1-p2;
        norm = pos.norm();
#ifndef PERIODIC
        Float3 los=p1+p2; // No 1/2 as normalized anyway below
        mu = fabs(pos.dot(los)/norm/los.norm());
#else
        // In the periodic case use z-direction for mu
        mu = fabs(pos.z/norm);
#endif
        }
    
public:
    void sum_counts(PowerCounts *pc){
        // Add counts accumulated in different threads
        for(int i=0;i<nbin*mbin;i++) power_counts[i]+=pc->power_counts[i];
        used_pairs+=pc->used_pairs;
    }
    
public:
    PowerCounts(Parameters *par, SurveyCorrection *_sc, KernelInterp *_kernel_interp){
        nbin = par->nbin; // number of radial bins
        mbin = par->mbin; // number of Legendre bins
        out_file = par->out_file; // output directory
        out_string = par->out_string; // type specifier string
        R0 = par-> R0; // truncation radius
        
        kernel_interp = new KernelInterp(_kernel_interp);
        
        sc = _sc;
        max_legendre = fmax((sc->l_bins-1)*2,par->max_l);
        
        int ec=0;
        ec+=posix_memalign((void **) &power_counts, PAGE, sizeof(double)*nbin*mbin);
        assert(ec==0);
        
        reset();
        
        box=par->perbox;

        rmax=par->rmax;
        rmin=par->rmin;
        
        r_high = par->radial_bins_high;
        r_low = par->radial_bins_low;
    }

    void reset(){
        for(int j=0;j<nbin*mbin;j++){
            power_counts[j]=0.;
        }
        used_pairs=0.;
    }
    
    ~PowerCounts(){
        free(power_counts);
    }
    
    inline void count_pairs(Particle pi, Particle pj){
        // Count pairs of particles with some weighting function        
        Float r_ij, mu_ij, w_ij, tmp_phi_inv=0, new_kr, old_kr,new_kernel,diff_kr,new_kr3, old_kr3;
        Float legendre[max_legendre/2+1]; // Even-ell Legendre coefficients
        Float old_kernel[mbin]; // saved kernel functions
        
        cleanup_l(pi.pos,pj.pos,r_ij,mu_ij); // define distance and angle
        
        if(r_ij>R0) return; // outside correct size
        
        used_pairs++; // update number of pairs used
        
        // First define all required legendre moments
        legendre[0]=1.;
        if(max_legendre>1){
            Float mu2 = mu_ij*mu_ij;
            legendre[1]=0.5*(3.*mu2-1.);
            if(max_legendre>3){
                Float mu4 = mu2*mu2;
                legendre[2]=1./8.*(35.*mu4-30.*mu2+3.);
                if(max_legendre>5){
                    Float mu6 = mu4*mu2;
                    legendre[3] = 1./16.*(231.*mu6-315.*mu4+105.*mu2-5.);
                    if(max_legendre>7){
                        Float mu8 = mu6*mu2;
                        legendre[4] = 1./128.*(6435.*mu8-12012.*mu6+6930.*mu4-1260.*mu2+35.);
                        if(max_legendre>9){
                            legendre[5] = 1./256.*(46189.*mu8*mu2-109395.*mu8+90090.*mu6-30030.*mu4+3465.*mu2-63.);
                        }
                    }
                }
            }
        }
        
        for(int l_i=0;l_i<sc->l_bins;l_i++) tmp_phi_inv+=legendre[l_i]*sc->inv_correction_function(l_i*2,r_ij);
        
        w_ij = pi.w*pj.w*pair_weight(r_ij);///tmp_phi_inv
        
        // NB: This assumes bins are contiguous
        for(int i=0;i<nbin;i++){
            // Now compute multipole contributions
            if(i==0){
                old_kr = r_ij*r_low[0];
                old_kr3 = pow(old_kr,3);
            }
            new_kr = r_ij*r_high[i];
            new_kr3 = pow(new_kr,3);
            diff_kr = new_kr3 - old_kr3;
            for(int j=0;j<mbin;j++){
                if(i==0) old_kernel[j] = kernel_interp->kernel(j*2,old_kr);
                new_kernel = kernel_interp->kernel(j*2,new_kr);
                power_counts[i*mbin+j]+=w_ij*legendre[j]*(new_kernel-old_kernel[j])/diff_kr;
                old_kernel[j] = new_kernel;
            }
            // Update values for next time
            old_kr = new_kr;
            old_kr3 = new_kr3;
        }
    }
    
    inline void count_pairs_old(Particle pi, Particle pj){
        // Count pairs of particles with some weighting function        
        Float r_ij, mu_ij, kr_ij_lo, kr_ij_hi,w_ij;
        Float legendre[max_legendre/2+1]; // Even-ell Legendre coefficients
        Float tmp_phi_inv=0.,tmp_bessel;
//#define BINNED
#ifdef BINNED
        Float k_sin_lo, k_cos_lo, k_sin_hi, k_cos_hi,k_Si_lo=0,k_Si_hi=0;
        Float diff_k3, tmp_kernel;
#else
        Float kr_pow[2*(mbin+1)];
        Float kr_mean,k_sin,k_cos;
#endif
        
        cleanup_l(pi.pos,pj.pos,r_ij,mu_ij); // define distance and angle
        
        if(r_ij>R0) return; // outside correct size
        
        used_pairs++; // update number of pairs used
        
        // First define all required legendre moments
        legendre[0]=1.;
        if(max_legendre>1){
            Float mu2 = mu_ij*mu_ij;
            legendre[1]=0.5*(3.*mu2-1.);
            if(max_legendre>3){
                Float mu4 = mu2*mu2;
                legendre[2]=1./8.*(35.*mu4-30.*mu2+3.);
                if(max_legendre>5){
                    Float mu6 = mu4*mu2;
                    legendre[3] = 1./16.*(231.*mu6-315.*mu4+105.*mu2-5.);
                    if(max_legendre>7){
                        Float mu8 = mu6*mu2;
                        legendre[4] = 1./128.*(6435.*mu8-12012.*mu6+6930.*mu4-1260.*mu2+35.);
                        if(max_legendre>9){
                            legendre[5] = 1./256.*(46189.*mu8*mu2-109395.*mu8+90090.*mu6-30030.*mu4+3465.*mu2-63.);
                        }
                    }
                }
            }
        }
        
        for(int l_i=0;l_i<sc->l_bins;l_i++) tmp_phi_inv+=legendre[l_i]*sc->inv_correction_function(l_i*2,r_ij);
        
        //TODO: expand pair weight to just use pre-computed r2,r3,r4??
#ifdef BINNED
        w_ij = pi.w*pj.w*pair_weight(r_ij)/tmp_phi_inv/pow(r_ij,3)*3;
#else
        w_ij = pi.w*pj.w*pair_weight(r_ij)/tmp_phi_inv;
#endif
        for(int i=0;i<nbin;i++){
            
            kr_ij_lo = r_ij*r_low[i];
            kr_ij_hi = r_ij*r_high[i];
#ifdef BINNED
            diff_k3 = pow(r_high[i],3)-pow(r_low[i],3);
            tmp_kernel = w_ij/diff_k3;
            k_sin_lo = sin(kr_ij_lo);
            k_sin_hi = sin(kr_ij_hi);
            k_cos_lo = cos(kr_ij_lo);
            k_cos_hi = cos(kr_ij_hi);
            // Define sine integrals
            if(max_legendre>1){
                k_Si_lo = gsl_sf_Si(kr_ij_lo);
                k_Si_hi = gsl_sf_Si(kr_ij_hi);
            }
            
            // Now compute multipole contributions
            for(int j=0;j<mbin;j++){
                if(j==0) tmp_bessel = (k_sin_hi-kr_ij_hi*k_cos_hi)-(k_sin_lo-kr_ij_lo*k_cos_lo);
                else if(j==1) tmp_bessel = (kr_ij_hi*k_cos_hi - 4.*k_sin_hi+3.*k_Si_hi) - (kr_ij_lo*k_cos_lo - 4.*k_sin_lo + 3.*k_Si_lo);
                else if(j==2) tmp_bessel = 0.5*(((105/kr_ij_hi-2.*kr_ij_hi)*k_cos_hi + (22.-105./pow(kr_ij_hi,2))*k_sin_hi + 15.*k_Si_hi)-((105/kr_ij_lo-2.*kr_ij_lo)*k_cos_lo+(22.-105./pow(kr_ij_lo,2))*k_sin_lo + 15.*k_Si_lo));
                else{
                    printf("\nOnly ell = 0,2,4 implemented\n");
                    exit(1);
                }
                power_counts[i*mbin+j]+= tmp_kernel*Float(pow(-1,j)*(4.*j+1))*legendre[j]*tmp_bessel;
            }
#else
            kr_mean = 0.5*(kr_ij_hi+kr_ij_lo);
            k_sin = sin(kr_mean);
            k_cos = cos(kr_mean);
            
            // Define powers of kr
            kr_pow[0]=1;
            for(int j=1;j<2*(mbin+1);j++) kr_pow[j] = kr_pow[j-1]*kr_mean;
            k_sin = sin(kr_mean);
            k_cos = cos(kr_mean);
            
            // Now compute multipole contributions
            for(int j=0;j<mbin;j++){
                if(j==0) tmp_bessel = k_sin/kr_pow[1];
                else if(j==1) tmp_bessel = (3.-kr_pow[2])*k_sin/kr_pow[3]-3.*k_cos/kr_pow[2];
                else if(j==2) tmp_bessel = 5.*(2*kr_pow[2]-21)*k_cos/kr_pow[4]+(kr_pow[4]-45*kr_pow[2]+105)/kr_pow[5]*k_sin;
                else{
                    printf("\nOnly ell = 0,2,4 implemented\n");
                    exit(1);
                }
                power_counts[i*mbin+j]+=w_ij*Float(pow(-1,j)*(4*j+1))*legendre[j]*tmp_bessel;
            }
#endif
        }
    }
    
    inline Float pair_weight(Float sep){
        // Compute weight function W(r;R_0)
        if(sep<R0/2) return 1.;
        else{
            Float x = sep/R0;
            if(sep<3*R0/4) return 1.-8.*pow(2*x-1,3)+8.*pow(2*x-1,4);
            else return -64.*pow(x-1,3)-32*pow(x-1,4);
        }
    }
    
    void save_counts(char* suffix) {
        // Print power-count output to file. 
        // Create output files
        
        char pow_name[1000];
        snprintf(pow_name, sizeof pow_name, "%s/%s_power_counts_n%d_m%d_%s.txt", out_file,out_string,nbin, mbin,suffix);
        FILE * PowFile = fopen(pow_name,"w"); 
        
        for (int i=0;i<nbin;i++){
            for(int j=0;j<mbin;j++){
                fprintf(PowFile,"%le\t",power_counts[i*mbin+j]);
            }
            fprintf(PowFile,"\n");
        }
                
        fflush(NULL);
        
        // Close open files
        fclose(PowFile);
    }
};

#endif
