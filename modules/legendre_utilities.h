// legendre_utilities.h - This contains utility functions for computing matrices binned in Legendre moments

#ifndef LEGENDRE_UTILITIES_H
#define LEGENDRE_UTILITIES_H

inline void legendre_polynomials(Float mu,int max_l, Float *poly_out){
    // Declare Legendre polynomials here
    Float mu2,mu4,mu6,mu8,mu10; // declare all powers of mu here
    
    // now compute relevant Legendre polynomials
    poly_out[0]=1;
    if(max_l>1){
        mu2 = mu*mu;
        poly_out[1] = 0.5*(3*mu2-1);
    }
    if(max_l>3){
        mu4 = mu2*mu2;
        poly_out[2] = 1/8*(35*mu4-30*mu2+3);
    }
    if(max_l>5){
        mu6 = mu4*mu2;
        poly_out[3] = 1/16*(231*mu6-315*mu4+105*mu2-5);
    }
    if(max_l>7){
        mu8 = mu6*mu2;
        poly_out[4] = 1/128*(6435*mu8-12012*mu6+6930*mu4-1260*mu2+35);
    }
    if(max_l>9){
        mu10 = mu8*mu2;
        poly_out[5] = 1/256*(46189*mu10-109395*mu8+90090*mu6-30030*mu4+3465*mu2-63);
    }
}

class SurveyCorrection{
    // this class stores the correction functions for each bin, giving the difference between the true and estimated RR counts. It is created by reading in coefficients to recompute smooth Phi(mu) functions for each radial bin.
    
public:        
    Float* phi_coeffs; // houses polynomial coefficients for the correction function
    int n_param = 7; // number of polynomial coefficients in fitting model (only 5 are independent)
private:
    int nbin; // number of radial bins (equal to par->nbin)
    Float mu_crit = 0.75; // critical mu to change between polynomial models
    
public:
    void copy(SurveyCorrection *sc){
        // Copy survey correction object
        n_param=sc->n_param;
        nbin=sc->nbin;
        mu_crit=sc->mu_crit;
        phi_coeffs=(Float *)malloc(sizeof(Float)*nbin*n_param);
        for(int i=0;i<nbin*n_param;i++) phi_coeffs[i]=sc->phi_coeffs[i];
    }
    
    // Empty operator
    SurveyCorrection(){};
    
    // Assignment operator creation
    SurveyCorrection& operator=(const SurveyCorrection& survey_corr);
    
    SurveyCorrection(Parameters *par, int index1, int index2){
        // This initializes the function and reads in the relevant polynomial coefficients for each radial bin
        // NB: coefficients are indexed as INDEX = RADIAL_BIN*N_COEFF + COEFF_ID where N_COEFF is the total number of coefficients for each model; here 7. 
        
        nbin = par->nbin; // number of radial bins 
        
        // READ IN FILE
        char line[1000000], *phi_file;
        int line_no = 0;
        FILE *fp;
        
        if((index1==1)&&(index2==1)) phi_file=par->phi_file;
        else if((index1==2)&&(index2==2)) phi_file=par->phi_file2;
        else phi_file=par->phi_file12;
        
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
        
        assert(line_no==nbin); // must have same number of coefficients as radial bins
        
        // Now allocate memory to the weights array
        int ec=0;
        ec+=posix_memalign((void **) &phi_coeffs, PAGE, sizeof(Float)*n_param*nbin);
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
            split_string = strtok(line, "\t");
                
            // Iterate over line
            while (split_string!=NULL){
                phi_coeffs[index]=atof(split_string);
                split_string = strtok(NULL,"\t");
                index++;
                }
            line_count++;
        }
            
        assert(line_count==nbin);
        assert(index==nbin*n_param); // to ensure we get all particles
            
        printf("Read in survey correction function coefficients successfully.\n\n"); 
            
        }
        
public:
    ~SurveyCorrection(){
        // The destructor
        free(phi_coeffs);
        return;
    }
    
    Float correction_function(int radial_bin, Float mu){
        // Create function to output polynomial model given the radial bin number
        int base_bin = radial_bin*n_param;
        if(mu<mu_crit) return phi_coeffs[base_bin]+phi_coeffs[base_bin+1]*mu+phi_coeffs[base_bin+2]*mu*mu;
        else{
            Float mu2 = mu*mu;
            return phi_coeffs[base_bin+3]+phi_coeffs[base_bin+4]*mu+phi_coeffs[base_bin+5]*mu2+phi_coeffs[base_bin+6]*mu2*mu;
        }
    }
};

#endif
