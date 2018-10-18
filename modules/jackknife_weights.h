// jackknife_weihts.h - this contains c++ functions to read-in jackknife weights and pair counts from files.

#ifndef JACKKNIFE_WEIGHTS_H
#define JACKKNIFE_WEIGHTS_H

// This class stores the RR count weights for a given jackknife
class JK_weights{
public:
    Float* weights; // houses the weights for each bin for this jackknife
    int* filled_JKs; // houses indices for the filled jackknife arrays
    Float* RR_pair_counts; // houses the weighted pair counts summed over jackknife regions.
    int n_JK_filled; // number of non-empty jackknife regions
    int nbins; // total number of bins
    Float* product_weights; // houses a matrix of SUM_A{w_aA*w_bA} terms for later use with indexing bin_a*nbins+bin_b
    
public:
    ~JK_weights() {
        // The destructor
        free(weights);
        free(RR_pair_counts);
        free(product_weights);
        free(filled_JKs);
        return;
    }
    
    JK_weights(Parameters *par){
        
        // This reads in weights for each jackknife region for each bin from file.
        // File should have bins in space-separated columns and bins in rows, with the jackknife number in the first column.
        // NB: The indexing is defined as INDEX = JACKKNIFE_ID * NBINS + BIN_ID
        
        nbins = par->nbin*par->mbin; // define number of bins in total
        char line[50000];
        n_JK_filled = 0; // line number
        FILE *fp;
        FILE *fp2;
        fp = fopen(par->jk_weight_file,"r");
        fp2 = fopen(par->RR_bin_file,"r");
        if (fp==NULL){
            fprintf(stderr,"Jackknife file %s not found\n",par->jk_weight_file);
            abort();
        }
        if (fp2==NULL){
            fprintf(stderr,"RR bin count file %s not found\n",par->RR_bin_file);
            abort();
        }
            
        fprintf(stderr,"\nReading jackknife file '%s' and RR bin count file '%s'\n",par->jk_weight_file,par->RR_bin_file);
        // Count lines to construct the correct size
        while (fgets(line,50000,fp)!=NULL){
            if (line[0]=='#') continue; // comment line
            if (line[0]=='\n') continue;
            n_JK_filled++;
        }
        printf("\n# Found %d non-empty jackknives in the file\n",n_JK_filled);
        rewind(fp); // restart file
        //TODO: Check this is consistent with expected number of jackknives
        
        // Now allocate memory to the weights array
        int ec=0;
        ec+=posix_memalign((void **) &weights, PAGE, sizeof(Float)*nbins*n_JK_filled);
        ec+=posix_memalign((void **) &filled_JKs, PAGE, sizeof(int)*n_JK_filled);
        ec+=posix_memalign((void **) &product_weights, PAGE, sizeof(Float)*nbins*nbins);
        ec+=posix_memalign((void **) &RR_pair_counts, PAGE, sizeof(Float)*nbins);
        assert(ec==0);
        
        int index=0; // index for jackknives and bins
        int line_count=0; // line counter
        int counter; // counts which element in line
        
        // Read in values to file
        while (fgets(line,50000,fp)!=NULL) {
            // Select required lines in file
            if (line[0]=='#') continue;
            if (line[0]=='\n') continue;
            
            // Split into variables
            char * split_string;
            split_string = strtok(line, "\t");
            counter=0;
            
            // Iterate over line
            while (split_string!=NULL){
                if(counter==0){
                    filled_JKs[line_count]=atoi(split_string);
                    }
                else{
                    weights[index]=atof(split_string);
                    index++;
                }
                split_string = strtok(NULL,"\t");
                counter++;
            }
            line_count++;
        }
        
        assert(line_count==n_JK_filled);
        assert(index==nbins*n_JK_filled); // to ensure we get all particles
        
        printf("Read in jackknife weights successfully.\n"); 
        
        index=0;
        while (fgets(line,1000,fp2)!=NULL){
            // Select required lines in file
            if (line[0]=='#') continue;
            if (line[0]=='\n') continue;
            RR_pair_counts[index]=atof(line);
            index++;
        }
        assert(index==nbins);
        printf("Read in RR pair counts successfully.\n");
        
        // Compute SUM_A(w_aA*w_bA) for all jackknives
        for (int bin_a=0;bin_a<nbins;bin_a++){
            for (int bin_b=0;bin_b<nbins;bin_b++){
                Float tmp_product=0.;
                for (int x=0;x<n_JK_filled;x++){ // sum over jackknives
                    tmp_product+=weights[x*nbins+bin_a]*weights[x*nbins+bin_b];
                product_weights[bin_a*nbins+bin_b]=tmp_product;
                }
            }
        }
        
    }
};
  
#endif
