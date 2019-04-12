// survey_correction_legendre.h - This contains utility functions for computing Legendre moments etc

#ifndef SURVEY_CORRECTION_LEGENDRE_H
#define SURVEY_CORRECTION_LEGENDRE_H

/*inline void legendre_polynomials_preload(Float* even_mu, int max_l, Float *poly_out){
    // Declare Legendre polynomials using pre-computed (even) powers of mu from correction function
    
    // Compute all Legendre polynomials in a list
    poly_out[0] = 1;
    
    if(max_l>1) poly_out[1] = 0.5*(3.*even_mu[0]-1.);
    if(max_l>3) poly_out[2] = 1./8.*(35.*even_mu[1]-30.*even_mu[0]+3.);
    if(max_l>5) poly_out[3] = 1./16.*(231.*even_mu[2]-315.*even_mu[1]+105.*even_mu[0]-5.);
    if(max_l>7) poly_out[4] = 1./128.*(6435.*even_mu[3]-12012.*even_mu[2]+6930.*even_mu[1]-1260.*even_mu[0]+35.);
    if(max_l>9) poly_out[5] = 1./256.*(46189.*even_mu[4]-109395.*even_mu[3]+90090.*even_mu[2]-30030.*even_mu[1]+3465.*even_mu[0]-63.);

}

inline void legendre_polynomials(Float mu,int max_l, Float *poly_out){
    // Declare Legendre polynomials here
    Float all_mu[5]; // declare all powers of mu here
    _
    // now compute relevant Legendre polynomials
    poly_out[0]=1;
    if(max_l>1){
        all_mu[0] = mu*mu; // mu^4
        poly_out[1] = 0.5*(3.*all_mu[0]-1.);
    }
    if(max_l>3){
        all_mu[1] = pow(all_mu[0],2.); // mu^4
        poly_out[2] = 1./8.*(35.*all_mu[1]-30.*all_mu[0]+3.);
    }
    if(max_l>5){
        all_mu[2] = all_mu[0]*all_mu[1]; // mu^6
        poly_out[3] = 1./16.*(231.*all_mu[2]-315.*all_mu[1]+105.*all_mu[0]-5.);
    }
    if(max_l>7){
        all_mu[3] = all_mu[0]*all_mu[2]; // mu^8
        poly_out[4] = 1./128.*(6435.*all_mu[3]-12012.*all_mu[2]+6930.*all_mu[1]-1260.*all_mu[0]+35.);
    }
    if(max_l>9){
        all_mu[4] = all_mu[0]*all_mu[3]; // mu^10
        poly_out[5] = 1./256.*(46189.*all_mu[4]-109395.*all_mu[3]+90090.*all_mu[2]-30030.*all_mu[1]+3465.*all_mu[0]-63.);
    }
}*/

class SurveyCorrection{
    // this class stores the correction functions for each bin, giving the difference between the true and estimated RR counts. It is created by reading in coefficients to  compute smooth Phi(r,mu).
    // For convenience we read in coefficients describing the monopoles of 1/Phi(r,mu)
    
public:        
    Float* phi_coeffs; // houses polynomial coefficients for the correction function
    int n_param = 3; // number of coefficients per multipole
    int l_bins = 3; // number of multipoles
    
public:
    void copy(SurveyCorrection *sc){
        // Copy survey correction object
        n_param=sc->n_param;
        l_bins=sc->l_bins;
        phi_coeffs=(Float *)malloc(sizeof(Float)*n_param*l_bins);
        for(int i=0;i<l_bins*n_param;i++) phi_coeffs[i]=sc->phi_coeffs[i];
    }
    
    // Empty operator
    SurveyCorrection(){};
    
    // Assignment operator creation
    SurveyCorrection& operator=(const SurveyCorrection& survey_corr);
        
    SurveyCorrection(Parameters *par){
        // This initializes the function and reads in the relevant polynomial coefficients for each radial bin. 
        // NB: coefficients are indexed as INDEX = MULTIPOLE_INDEX*N_COEFF. + COEFF_ID where N_COEFF is the total number of coefficients for each multipole; here 3. 
        
        // READ IN FILE
        char line[1000000], *phi_file;
        int line_no = 0;
        FILE *fp;
        
        phi_file = par->phi_file;
        fp = fopen(phi_file,"r");
        if (fp==NULL){
            fprintf(stderr,"Survey correction function coefficient file %s not found\n",phi_file);
            abort();
        }
                
        printf("\nReading survey correction function coefficient file '%s'\n",phi_file);
        
        // Count lines to construct the correct size
        while (fgets(line,1000000,fp)!=NULL){
            if (line[0]=='#') continue; // comment line
            if (line[0]=='\n') continue;
            line_no++;
        }
        rewind(fp); // restart file
        
        int ec=0;
        // Now allocate memory to the weights array
        assert(line_no==l_bins); // need correct number of functions
        ec+=posix_memalign((void **) &phi_coeffs, PAGE, sizeof(Float)*n_param*l_bins);
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
            
        assert(line_count==l_bins);
        assert(index==n_param*l_bins);
        printf("Read in survey correction function coefficients successfully.\n\n"); 
            
        }
        
public:
    ~SurveyCorrection(){
        // The destructor
        free(phi_coeffs);
        return;
    }
    
    Float inv_correction_function(int ell, Float r){
        int l_index = ell/2,base_bin = l_index*n_param;
        return phi_coeffs[base_bin]+phi_coeffs[base_bin+1]*r+phi_coeffs[base_bin+2]*pow(r,2);
    }
};

#endif
