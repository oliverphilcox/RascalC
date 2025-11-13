// legendre_utilities.h - This contains utility functions for computing matrices binned in Legendre moments

#ifndef LEGENDRE_UTILITIES_H
#define LEGENDRE_UTILITIES_H

#ifdef THREE_PCF
const int max_no_multipoles = 11; // for array allocation on stack - no variable length arrays in C++ standards. needs to be consistent with (not less than) the number of multipoles written below
inline void legendre_polynomials_preload(Float mu, int max_l, Float *poly_out){
    // Declare Legendre polynomials using pre-computed powers of mu from correction function
    
    // Compute all Legendre polynomials in a list
    poly_out[0] = 1.;
    
    if(max_l>0){
        poly_out[1] = mu;
        if(max_l>1){
            Float mu2 = mu*mu;
            poly_out[2] = 0.5*(3.*mu2-1.);
            if(max_l>2){
                Float mu3 = mu2*mu;
                poly_out[3] = 0.5*(5.*mu3-3.*mu);
                if(max_l>3){
                    Float mu4 = mu3*mu;
                    poly_out[4] = 1./8.*(35.*mu4-30.*mu2+3.);
                    if(max_l>4){
                        Float mu5 = mu4*mu;
                        poly_out[5] = 1./8.*(63.*mu5-70.*mu3+15.*mu);
                        if(max_l>5){
                            Float mu6 = mu5*mu;
                            poly_out[6] = 1./16.*(231.*mu6-315.*mu4+105.*mu2-5.);
                            if(max_l>6){
                                Float mu7 = mu6*mu;
                                poly_out[7] = 1./16.*(429.*mu7-693.*mu5+315.*mu3-35.*mu);
                                if(max_l>7){
                                    Float mu8 = mu7*mu;
                                    poly_out[8] = 1./128.*(6435.*mu8-12012.*mu6+6930.*mu4-1260.*mu2+35.);
                                    if(max_l>8){
                                        Float mu9 = mu8*mu;
                                        poly_out[9] = 1./128.*(12155.*mu9-25740.*mu7+18018.*mu5-4620.*mu3+315.*mu);
                                        if(max_l>9){
                                            Float mu10 = mu9*mu;
                                            poly_out[10] = 1./256.*(46189.*mu10-109395.*mu8+90090.*mu6-30030.*mu4+3465.*mu2-63.);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
#else
const int max_no_multipoles = 6; // for array allocation on stack - no variable length arrays in C++ standards. needs to be consistent with (not less than) the number of multipoles written below
inline void legendre_polynomials(Float mu,int max_l, Float *poly_out){
    // Declare Legendre polynomials here
    // now compute relevant Legendre polynomials
    poly_out[0]=1;
    if(max_l>1){
        Float mu2 = mu*mu; // mu^2
        poly_out[1] = 0.5*(3.*mu2-1.);
        if(max_l>3){
            Float mu4 = pow(mu2,2.); // mu^4
            poly_out[2] = 1./8.*(35.*mu4-30.*mu2+3.);
            if(max_l>5){
                Float mu6 = mu4*mu2;
                poly_out[3] = 1./16.*(231.*mu6-315.*mu4+105.*mu2-5.);
                if(max_l>7){
                    Float mu8 = mu6*mu2; // mu^8
                    poly_out[4] = 1./128.*(6435.*mu8-12012.*mu6+6930.*mu4-1260.*mu2+35.);
                    if(max_l>9){
                        Float mu10 = mu8*mu2; // mu^10
                        poly_out[5] = 1./256.*(46189.*mu10-109395.*mu8+90090.*mu6-30030.*mu4+3465.*mu2-63.);
                    }
                }
            }
        }
    }
}
#endif

#ifdef POWER
    // Use the inverse Phi definition in this case
    #include "../power_spectra/power_mod/survey_correction_legendre.h"
#else
    // Create new survey correction definition

class SurveyCorrection{
    // this class stores the correction functions for each bin, giving the difference between the true and estimated RR counts. It is created by reading in coefficients to recompute smooth Phi(mu) functions for each radial bin.
    
public:
#ifdef THREE_PCF
    Float* inv_phi_multipoles = nullptr; // houses multipole moments of the inverse correction function
    int n_param = 7;
#else
    Float* phi_coeffs = nullptr; // houses polynomial coefficients for the correction function
    int n_param = 7; // number of polynomial coefficients in fitting model (only 5 are independent)
    Float mu_crit = 0.75; // critical mu to change between polynomial models
#endif
private:
    int nbin; // number of radial bins (equal to par->nbin)
    
public:
    void copy(SurveyCorrection *sc){
        // Copy survey correction object
        n_param=sc->n_param;
        nbin=sc->nbin;
#ifdef THREE_PCF
        inv_phi_multipoles = (Float *)malloc(sizeof(Float)*nbin*nbin*n_param);
        for(int i=0; i<nbin*nbin*n_param; i++) inv_phi_multipoles[i] = sc->inv_phi_multipoles[i];
#else
        mu_crit=sc->mu_crit;
        phi_coeffs=(Float *)malloc(sizeof(Float)*nbin*n_param);
        for(int i=0;i<nbin*n_param;i++) phi_coeffs[i]=sc->phi_coeffs[i];
#endif
    }
    
#ifdef THREE_PCF
    void rescale(Float norm1, Float norm2, Float norm3){
        // Rescale the survey correction function by a factor (N_rand_1/N_gal_1)*(N_rand2/N_gal2)*(N_rand3/N_gal3) - inverse of RRR counts rescaling (if it was implemented)
        Float rescale_factor = norm1*norm2*norm3;
        printf("Rescaling survey correction function by a factor (N_rand_1/N_gal_1)*(N_rand2/N_gal2)*(N_rand3/N_gal3) = %.1e\n", rescale_factor);
        for(int i = 0; i < nbin * nbin * n_param; i++) inv_phi_multipoles[i] /= rescale_factor; // need to divide because the coefficients are actually for the inverse of the correction function
    }
#else
    void rescale(Float norm1, Float norm2){
        // Rescale the survey correction function by a factor (N_rand_1/N_gal_1)*(N_rand2/N_gal2) - inverse of RR counts rescaling
        Float rescale_factor = norm1*norm2;
        printf("Rescaling survey correction function by a factor (N_rand_1/N_gal_1)*(N_rand2/N_gal2) = %.1e\n", rescale_factor);
        for(int i = 0; i < nbin * n_param; i++) phi_coeffs[i] *= rescale_factor;
    }
#endif
    
    // Empty operator
    SurveyCorrection(){};
    
    // Assignment operator creation
    SurveyCorrection& operator=(const SurveyCorrection& survey_corr);
    
#ifdef THREE_PCF
    SurveyCorrection(Parameters *par){
        // This initializes the function and reads in the relevant polynomial coefficients for each radial bin. 
        // NB: coefficients are indexed as INDEX = RADIAL_BIN_1*N_COEFF^2. + RADIAL_BIN_2*N_COEFF + COEFF_ID where N_COEFF is the total number of coefficients for each model; here 7.
#else
    SurveyCorrection(Parameters *par, int index1, int index2){
        // This initializes the function and reads in the relevant polynomial coefficients for each radial bin
        // NB: coefficients are indexed as INDEX = RADIAL_BIN*N_COEFF + COEFF_ID where N_COEFF is the total number of coefficients for each model; here 7. 
#endif
        
        nbin = par->nbin; // number of radial bins 
        
        // READ IN FILE
        char line[1000000], *phi_file;
        int line_no = 0;
        FILE *fp;
        
#if (defined THREE_PCF || defined POWER)
        phi_file = par->phi_file;
#else
        if((index1==1)&&(index2==1)) phi_file=par->phi_file;
        else if((index1==2)&&(index2==2)) phi_file=par->phi_file2;
        else phi_file=par->phi_file12;
#endif  
        fp = fopen(phi_file,"r");
        if (fp==NULL){
            fprintf(stderr,"Survey correction function coefficient file %s not found\n",phi_file);
            abort();
        }
                
        fprintf(stderr,"\nReading survey correction function coefficient file '%s'\n",phi_file);
        
        // Count lines to construct the correct size
        while (fgets(line,1000000,fp)!=NULL){
            if (line[0]=='#') continue; // comment line
            if (line[0]=='\n') continue;
            line_no++;
        }
        rewind(fp); // restart file
        
        int ec=0;
        // Now allocate memory to the weights array
#ifdef THREE_PCF
        assert(line_no==nbin*nbin); // need correct number of functions
        ec+=posix_memalign((void **) &inv_phi_multipoles, PAGE, sizeof(Float)*n_param*nbin*nbin);
#else
        assert(line_no==nbin); // need correct number of functions
        ec+=posix_memalign((void **) &phi_coeffs, PAGE, sizeof(Float)*n_param*nbin);
#endif
        assert(ec==0);
        int line_count=0; // line counter
        int index=0; // indexes array
            
        // Read in values to file
        while (fgets(line,1000000,fp)!=NULL) {
            // Select required lines in file
        
            if (line[0]=='#') continue;
            if (line[0]=='\n') continue;
                
            // Split into variables
            char * split_string;
            split_string = strtok(line, " \t");
                
            // Iterate over line
            while (split_string!=NULL){
#ifdef THREE_PCF
                inv_phi_multipoles[index] = atof(split_string);
#else
                phi_coeffs[index] = atof(split_string);
#endif
                split_string = strtok(NULL, " \t");
                index++;
                }
            line_count++;
        }
            
#ifdef THREE_PCF
        assert(line_count==nbin*nbin);
        assert(index==nbin*nbin*n_param);
#else
        assert(line_count==nbin);
        assert(index==nbin*n_param); // to ensure we get all particles
#endif       
        printf("Read in survey correction function coefficients successfully.\n\n"); 
            
        }
        
public:
    ~SurveyCorrection(){
        // The destructor
#ifdef THREE_PCF
        free(inv_phi_multipoles);
#else
        free(phi_coeffs);
#endif
        return;
    }
    
#ifdef THREE_PCF
    Float inv_correction_function_3pcf(int radial_bin1, int radial_bin2, int ell){
        // Create function to output polynomial RRR correction model given two radial bin numbers and a Legendre polynomial index
        // This gives the ell-th multipole of 1/Phi for convenience
        int base_bin = (radial_bin1*nbin+radial_bin2)*n_param;
        //if(ell==0) return 1;
        //else return 0;
        
        return inv_phi_multipoles[base_bin+ell];
        
    }
#else
    Float correction_function(int radial_bin, Float mu){
        // Create function to output polynomial RR correction model given the radial bin number
        int base_bin = radial_bin*n_param;
        if(mu<mu_crit) return phi_coeffs[base_bin]+phi_coeffs[base_bin+1]*mu+phi_coeffs[base_bin+2]*mu*mu;
        else{
            Float mu2 = mu*mu;
            return phi_coeffs[base_bin+3]+phi_coeffs[base_bin+4]*mu+phi_coeffs[base_bin+5]*mu2+phi_coeffs[base_bin+6]*mu2*mu;
        }
    }
#endif
};

#endif
#endif
