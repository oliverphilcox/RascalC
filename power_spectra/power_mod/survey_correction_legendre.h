// survey_correction_legendre.h - This contains utility functions for computing Legendre moments etc

#ifndef SURVEY_CORRECTION_LEGENDRE_H
#define SURVEY_CORRECTION_LEGENDRE_H

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
        
    SurveyCorrection(Parameters *par, int I1, int I2){
        // This initializes the function and reads in the relevant polynomial coefficients for each radial bin. 
        // NB: coefficients are indexed as INDEX = MULTIPOLE_INDEX*N_COEFF. + COEFF_ID where N_COEFF is the total number of coefficients for each multipole; here 3. 
        
        // READ IN FILE
        char line[1000000], *phi_file;
        int line_no = 0;
        FILE *fp;
        
        if((I1==1)&&(I2==1))    phi_file = par->inv_phi_file;
        else if ((I1==2)&&(I2==2)) phi_file = par->inv_phi_file2;
        else phi_file = par->inv_phi_file12;
        
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
        //if(ell==0) return 1;
        //else return 0;
        return phi_coeffs[base_bin]+phi_coeffs[base_bin+1]*r+phi_coeffs[base_bin+2]*pow(r,2);
    }
};

#endif
