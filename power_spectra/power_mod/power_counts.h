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
    char* out_file;
    SurveyCorrection *sc; // survey correction function
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
    PowerCounts(Parameters *par, SurveyCorrection *_sc){
        nbin = par->nbin; // number of radial bins
        mbin = par->mbin; // number of Legendre bins
        out_file = par->out_file; // output directory
        R0 = par-> R0; // truncation radius
        
        sc = _sc;
        
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

    ~PowerCounts(){
        free(power_counts);
    }
    
    void reset(){
        for(int j=0;j<nbin*mbin;j++){
            power_counts[j]=0.;
        }
        used_pairs=0.;
    }
    
    inline void count_pairs(Particle pi, Particle pj){
        // Count pairs of particles with some weighting function        
        Float r_ij, mu_ij, kr_ij_lo, kr_ij_hi,w_ij,k_sin,k_cos;
        
        cleanup_l(pi.pos,pj.pos,r_ij,mu_ij); // define distance and angle
        
        if(r_ij>R0) return; // outside correct size
        
        used_pairs++; // update number of pairs used
        
        for(int i=0;i<nbin;i++){
            w_ij = pi.w*pj.w*pair_weight(r_ij)/sc->inv_correction_function(r_ij,mu_ij);
            
            kr_ij_lo = r_ij*r_high[i];
            kr_ij_hi = r_ij*r_low[i];
            
            k_sin = sin(kr_ij_lo);
            k_cos = cos(kr_ij_hi);
            
            for(int j=0;j<mbin;j++){
                power_counts[i*mbin+j]+=w_ij*(k_sin+k_cos);
            }
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
        snprintf(pow_name, sizeof pow_name, "%s/power_counts_n%d_m%d_%s.txt", out_file,nbin, mbin,suffix);
        FILE * PowFile = fopen(pow_name,"w"); 
        
        for (int j=0;j<nbin*mbin;j++) fprintf(PowFile,"%le\n",power_counts[j]);
        
        fflush(NULL);
        
        // Close open files
        fclose(PowFile);
    }
};

#endif
