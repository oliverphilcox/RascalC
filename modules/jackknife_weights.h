// jackknife_weights.h - this contains c++ functions to read-in jackknife weights and pair counts from files.

#ifndef JACKKNIFE_WEIGHTS_H
#define JACKKNIFE_WEIGHTS_H

// This class stores the RR count weights for a given jackknife
class JK_weights{
public:
    Float* RR_pair_counts; // houses the weighted pair counts summed over jackknife regions.
    int nbins; // total number of bins
#ifdef JACKKNIFE
    Float* weights; // houses the weights for each bin for this jackknife
    int* filled_JKs; // houses indices for the filled jackknife arrays
    int n_JK_filled; // number of non-empty jackknife regions
    Float* product_weights; // houses a matrix of SUM_A{w_aA*w_bA} terms for later use with indexing bin_a*nbins+bin_b
#endif
    
public: 
    
    void copy(JK_weights *JK){
        // Copy JK_weights object
        nbins=JK->nbins;
        RR_pair_counts = (Float *)malloc(sizeof(Float)*nbins);
#ifdef JACKKNIFE
        n_JK_filled=JK->n_JK_filled;
        weights = (Float *)malloc(sizeof(Float)*nbins*n_JK_filled);
        filled_JKs = (int *)malloc(sizeof(int)*n_JK_filled);
        product_weights = (Float *)malloc(sizeof(Float)*nbins*nbins);
        int ct=0;
#endif
        for(int i=0;i<nbins;i++){
            RR_pair_counts[i]=JK->RR_pair_counts[i];
#ifdef JACKKNIFE
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
#endif
        }
    }
    
    void rescale(Float norm1, Float norm2){
        // Rescale the RR pair counts by a factor (N_gal1/N_rand1)*(N_gal2/N_rand2)
        Float rescale_factor = norm1*norm2;
        printf("Rescaling RR pair counts by a factor (N_gal_1/N_rand_1)*(N_gal2/N_rand2) = %.1e\n",1./rescale_factor);
        for(int i=0;i<nbins;i++){
            RR_pair_counts[i]/=rescale_factor;
        }
    }
    
public:
    ~JK_weights() {
        // The destructor
        free(RR_pair_counts);
#ifdef JACKKNIFE
        free(weights);
        free(product_weights);
        free(filled_JKs);
#endif
        return;
    }
    
    JK_weights(){};
    
    // Assignment operator creation
    JK_weights& operator=(const JK_weights& jk);
    
    JK_weights(Parameters *par, int index1, int index2){
        
        // This reads in weights for each jackknife region for each bin from file.
        // File should have bins in space-separated columns and bins in rows, with the jackknife number in the first column.
        // NB: The indexing is defined as INDEX = JACKKNIFE_ID * NBINS + BIN_ID
        // index1 and index2 define which random set of particles to use here
        // If Jackknife directive is not set, we only read in RR pair counts here
        
        nbins = par->nbin*par->mbin; // define number of bins in total
        char line[1000000];
        FILE *fp2;
        char *RR_file;
        
        if (!(par->make_random && par->perbox)) { // if particles are made random in periodic box RR counts can (and will) be done analytically
        if((index1==1)&&(index2==1)) RR_file = par->RR_bin_file;
        else if((index1==2)&&(index2==2)) RR_file = par->RR_bin_file2;
        else RR_file = par->RR_bin_file12;
        fp2 = fopen(RR_file,"r");
        
        if (fp2==NULL){
            fprintf(stderr,"RR bin count file %s not found\n",RR_file);
            abort();
        }
        fprintf(stderr,"\nReading RR bin count file '%s'\n",RR_file);
        }
        
        int ec2=0;
        ec2+=posix_memalign((void **) &RR_pair_counts, PAGE, sizeof(Float)*nbins);
        assert(ec2==0);
        
        if (par->make_random && par->perbox) {
            for (int i = 0; i < par->nbin; i++) {
                Float counts = pow(par->np, 2) * 4*M_PI/3 * (pow(par->radial_bins_high[i], 3) - pow(par->radial_bins_low[i], 3)) / pow(par->boxsize, 3) / par->mbin;
                for (int j = 0; j < par->mbin; j++) RR_pair_counts[i * par->mbin + j] = counts;
            }
        }
        else {
        int index=0;
        while (fgets(line,5000,fp2)!=NULL){
            // Select required lines in file
            if (line[0]=='#') continue;
            if (line[0]=='\n') continue;
            RR_pair_counts[index]=atof(line);
            index++;
        }
        assert(index==nbins);
        printf("Read in RR pair counts successfully.\n");
        }
        
#ifdef JACKKNIFE
        char *jk_file;
        FILE *fp;
        n_JK_filled = 0; // line number
        
        if((index1==1)&&(index2==1)) jk_file=par->jk_weight_file;
        else if((index1==2)&&(index2==2)) jk_file = par->jk_weight_file2;
        else jk_file = par->jk_weight_file12;        
        fp = fopen(jk_file,"r");
        
        if (fp==NULL){
            fprintf(stderr,"Jackknife file %s not found\n",jk_file);
            abort();
        }
        fprintf(stderr,"\nReading jackknife file '%s'\n",jk_file);
        
        // Count lines to construct the correct size
        while (fgets(line,1000000,fp)!=NULL){
            if (line[0]=='#') continue; // comment line
            if (line[0]=='\n') continue;
            n_JK_filled++;
        }
        printf("\n# Found %d non-empty jackknives in the file\n",n_JK_filled);
        rewind(fp); // restart file
        // Now allocate memory to the weights array
        int ec=0;
        ec+=posix_memalign((void **) &weights, PAGE, sizeof(Float)*nbins*n_JK_filled);
        ec+=posix_memalign((void **) &filled_JKs, PAGE, sizeof(int)*n_JK_filled);
        ec+=posix_memalign((void **) &product_weights, PAGE, sizeof(Float)*nbins*nbins);
        assert(ec==0);
            
        index=0; // index for jackknives and bins
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
        
        // Compute SUM_A(w_aA*w_bA) for all jackknives
        int partial_bin=0,partial_bin2;
        Float this_weight;
        for (int x=0;x<n_JK_filled;x++){
            //printf("Computing product weights for jackknife %d of %d\n",x+1,n_JK_filled);
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
        
        printf("Computed product weights successfully.\n");        
#endif
    }
};
  
#endif
