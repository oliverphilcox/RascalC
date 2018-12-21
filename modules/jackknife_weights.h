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
    
    void copy(JK_weights *JK){
        // Copy JK_weights object
        n_JK_filled=JK->n_JK_filled;
        nbins=JK->nbins;
        weights = (Float *)malloc(sizeof(Float)*nbins*n_JK_filled);
        filled_JKs = (int *)malloc(sizeof(int)*n_JK_filled);
        RR_pair_counts = (Float *)malloc(sizeof(Float)*nbins);
        product_weights = (Float *)malloc(sizeof(Float)*nbins*nbins);

        int ct=0;
        for(int i=0;i<nbins;i++){
            RR_pair_counts[i]=JK->RR_pair_counts[i];
            for(int j=0;j<nbins;j++){
                product_weights[ct]=JK->product_weights[ct];
                ct++;
            }
        }
        ct=0;
        for(int i=0;i<n_JK_filled;i++){
            filled_JKs[i]=JK->filled_JKs[i];
            for(int j=0;j<nbins;j++){
                weights[ct]=JK->weights[ct];
                ct++;
            }
        }
    }
    
    void rescale(Float n_gal1, Float n_gal2, int n_rand1, int n_rand2){
        // Rescale the RR pair counts by a factor (N_gal1/N_rand1)*(N_gal2/N_rand2)
        Float rescale_factor = (n_gal1*n_gal2)/(float(n_rand1)*float(n_rand2));
        printf("Rescaling RR pair counts by a factor (N_gal_1/N_rand_1)*(N_gal2/N_rand2) = %.1e\n",rescale_factor);
        for(int i=0;i<nbins;i++){
            RR_pair_counts[i]*=rescale_factor;
        }
    }
    
public:
    ~JK_weights() {
        // The destructor
        free(weights);
        free(RR_pair_counts);
        free(product_weights);
        free(filled_JKs);
        return;
    }
    
    JK_weights(){};
    
    //JK_weights(JK_weights* JK){
    //    // Constructor to copy an existing file:
    //    JK->copy(&weights, &filled_JKs, &RR_pair_counts, n_JK_filled, nbins, &product_weights);
    //}
    
    // Assignment operator creation
    JK_weights& operator=(const JK_weights& jk);
    
    JK_weights(Parameters *par, int index1, int index2){
        
        // This reads in weights for each jackknife region for each bin from file.
        // File should have bins in space-separated columns and bins in rows, with the jackknife number in the first column.
        // NB: The indexing is defined as INDEX = JACKKNIFE_ID * NBINS + BIN_ID
        // index1 and index2 define which random set of particles to use here
        
        nbins = par->nbin*par->mbin; // define number of bins in total
        char line[1000000];
        n_JK_filled = 0; // line number
        FILE *fp;
        FILE *fp2;
        char *jk_file;
        char *RR_file;
        
        if((index1==1)&&(index2==1)){
            jk_file=par->jk_weight_file;
            RR_file = par->RR_bin_file;
        }
        else if((index1==2)&&(index2==2)){
            jk_file = par->jk_weight_file2;
            RR_file = par->RR_bin_file2;
        }
        else{
            jk_file = par->jk_weight_file12;
            RR_file = par->RR_bin_file12;
        }
        
        fp = fopen(jk_file,"r");
        fp2 = fopen(RR_file,"r");
        if (fp==NULL){
            fprintf(stderr,"Jackknife file %s not found\n",jk_file);
            abort();
        }
        if (fp2==NULL){
            fprintf(stderr,"RR bin count file %s not found\n",RR_file);
            abort();
        }
            
        fprintf(stderr,"\nReading jackknife file '%s' and RR bin count file '%s'\n",jk_file,RR_file);
        // Count lines to construct the correct size
        while (fgets(line,1000000,fp)!=NULL){
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
        while (fgets(line,1000000,fp)!=NULL) {
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
        while (fgets(line,5000,fp2)!=NULL){
            // Select required lines in file
            if (line[0]=='#') continue;
            if (line[0]=='\n') continue;
            RR_pair_counts[index]=atof(line);
            index++;
        }
        assert(index==nbins);
        printf("Read in RR pair counts successfully.\n");
        
        // Compute SUM_A(w_aA*w_bA) for all jackknives
        int partial_bin=0,partial_bin2;
        Float this_weight;
        for (int x=0;x<n_JK_filled;x++){
            printf("Computing product weights for jackknife %d of %d\n",x+1,n_JK_filled);
            partial_bin2=0;
            for(int bin_a=0;bin_a<nbins;bin_a++){
                this_weight=weights[partial_bin+bin_a];
                for(int bin_b=0;bin_b<nbins;bin_b++){
                    product_weights[partial_bin2+bin_b]+=this_weight*weights[partial_bin+bin_b];
                }
                partial_bin2+=nbins;
            }
            partial_bin+=nbins;
        }
        
        /*Float tmp_product;
        int partial_bin;
        for (int bin_a=0;bin_a<nbins;bin_a++){
            for (int bin_b=0;bin_b<nbins;bin_b++){
                tmp_product=0.;
                for (int x=0;x<n_JK_filled;x++){ // sum over jackknives
                    tmp_product+=weights[x*nbins+bin_a]*weights[x*nbins+bin_b];
                product_weights[bin_a*nbins+bin_b]=tmp_product;
                }
            }
        }
        */
        printf("Computed product weights successfully.\n");        
    }
};
  
#endif
