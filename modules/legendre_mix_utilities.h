// legendre_mix_utilities.h - This contains utility functions for computing matrices binned in Legendre moments projected from s/r and mu binned correlation function

#ifndef LEGENDRE_MIX_UTILITIES_H
#define LEGENDRE_MIX_UTILITIES_H

class MuBinLegendreFactors{
    // this class stores the (2*ell+1) times integral of Legendre multipoles for each mu bin.
    // Those are the weights with which radially and angularly binned 2PCF values contribute to Legendre moments
    
public:
    Float* data_array; // houses polynomial coefficients for the correction function
private:
    int mbin; // number of angular bins (equal to par->mbin)
    int n_l; // number of ells
    
public:
    void copy(MuBinLegendreFactors *mblf){
        // Copy survey correction object
        mbin = mblf->mbin;
        n_l = mblf->n_l;
        data_array = (Float *)malloc(sizeof(Float)*n_l*mbin);
        for(int i = 0; i < n_l*mbin; i++) data_array[i] = mblf->data_array[i];
    }
    
    // Empty operator
    MuBinLegendreFactors(){};
    
    MuBinLegendreFactors(char* filename, int _mbin, int max_l){
        // This initializes the function and reads in the relevant polynomial coefficients for each radial bin
        // NB: data are indexed as INDEX = MU_BIN*N_L + ELL_INDEX
        
        mbin = _mbin; // number of angular bins
        n_l = max_l/2+1; // number of multipoles, even only
        
        // READ IN FILE
        char line[1000000];
        int line_no = 0;
        FILE *fp;
        
        fp = fopen(filename, "r");
        if (fp == NULL) {
            fprintf(stderr, "Mu bin Legendre factors file %s not found\n", filename);
            abort();
        }
                
        fprintf(stderr, "\nMu bin Legendre factors file '%s'\n", filename);
        
        // Count lines to construct the correct size
        while (fgets(line, 1000000, fp) != NULL) {
            if (line[0] == '#') continue; // comment line
            if (line[0] == '\n') continue;
            line_no++;
        }
        rewind(fp); // restart file
        
        int ec = 0;
        // Now allocate memory to the array
        assert(line_no == mbin); // need correct number of rows
        ec += posix_memalign((void **) &data_array, PAGE, sizeof(Float)*n_l*mbin);
        assert(ec == 0);
        int line_count = 0; // line counter
        int index = 0; // indexes array
            
        // Read in values from file
        while (fgets(line, 1000000, fp) != NULL) {
            // Select required lines in file
        
            if (line[0] == '#') continue;
            if (line[0] == '\n') continue;
                
            // Split into variables
            char * split_string;
            split_string = strtok(line, " \t");
                
            // Iterate over line
            while (split_string != NULL){
                data_array[index] = atof(split_string);
                split_string = strtok(NULL, " \t");
                index++;
                }
            line_count++;
        }

        assert(line_count==mbin);
        assert(index==n_l*mbin); // to ensure we get all the data

        printf("Read in mu bin Legendre factors successfully.\n\n"); 
            
        }

    ~MuBinLegendreFactors() {
        // The destructor
        free(data_array);
        return;
    }
    
    inline Float* get_factors(int mu_bin) {
        // Simple function to return the indexable pointer to the n_l factors for a given angular bin
        return data_array + mu_bin * n_l;
    }
};

#endif
