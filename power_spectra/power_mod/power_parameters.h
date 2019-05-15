
// parameter function file for grid_power.cpp 

#ifndef PARAMETERS_H
#define PARAMETERS_H

class Parameters{

public:
	// Important variables to set!  Here are the defaults:
	
    //---------- ESSENTIAL PARAMETERS -----------------
    
    // Name of the first particle field
    char *fname = NULL;
    const char default_fname[500] = "/mnt/store1/oliverphilcox/QPM_proc/qpm_galaxy_0001.xyzw";//PowerSpectra/qpm_galaxy_1.xyzwj";//randoms_10x.xyzwj"; 
    
    // Name of second particle field 
    char *fname2 = NULL;
	const char default_fname2[500] = "/mnt/store1/oliverphilcox/QPM_proc/qpm_galaxy_0001.xyzw";//PowerSpectra/qpm_galaxy_1.xyzwj";// "/mnt/store1/oliverphilcox/3PCF/qpm_galaxy_1.xyzwj";
    
    // Optional File Prefix for output
    char *out_string = NULL;
    const char default_out_string[500] = "DD"; 
    
    // Name of the radial binning .csv file in k-space
    char *radial_bin_file = NULL;
    const char default_radial_bin_file[500] = "/mnt/store1/oliverphilcox/PowerSpectra/k_binning2.csv";
    
    // Output directory 
    char *out_file = NULL;
    const char default_out_file[500] = "/mnt/store1/oliverphilcox/PowerQPM/";
    
    // The number of threads to run on
	int nthread = 20;

    // The grid size, which should be tuned to match boxsize. 
	// This uses the maximum width of the cuboidal box.
	int nside = 151;
    
    // Whether or not we are using a periodic box
	bool perbox = false;

    int max_l = 4; // max Legendre moment (must be even)
    
    Float R0 = 100; // kernel truncation radius (in Mpc/h)
    
    char *inv_phi_file = NULL; // Survey correction function coefficient file 
    const char default_inv_phi_file[500] = "/mnt/store1/oliverphilcox/PowerSpectra/InvPhiCoeff_DR12.txt";
    
    //-------- OTHER PARAMETERS ----------------------------------------------
    
	// The periodicity of the position-space cube.
	Float boxsize = 200; // this is only used if the input particles are made randomly
    
	// The particles will be read from the unit cube, but then scaled by boxsize.
	Float rescale = 1.;   // If left zero or negative, set rescale=boxsize

	// The maximum number of points to read
	uint64 nmax = 1000000000000;
	
    // The location and name of a integrated grid of probabilities to be saved
	char *savename = NULL;
    // The location and name of a integrated grid of probabilities to be loaded
	char *loadname = NULL; //	
	
	// Whether to balance the weights or multiply them by -1
	int qinvert = 0, qbalance = 0;
    
	// If set, we'll just throw random periodic points instead of reading the file
	int make_random = 0;

	// Will be number of particles in a random distribution, but gets overwritten if reading from a file.
	int np = -1; // NB: This is only used for grid creation so we don't need a separate variable for the second set of randoms

	// The index from which on to invert the sign of the weights
	int rstart = 0;

	//---------------- INTERNAL PARAMETERS -----------------------------------
    // (no more user defined parameters below this line)
    
    // For consistency with other modules    
    char *inv_phi_file2 = NULL; // Survey correction function coefficient file 
    char *inv_phi_file12 = NULL; // Survey correction function coefficient file 
    
	// The periodicity of the position-space cuboid in 3D. 
    Float3 rect_boxsize = {boxsize,boxsize,boxsize}; // this is overwritten on particle read-in
    
    // Radial binning parameters (will be set from file)
    int nbin=0,mbin;
    Float rmin, rmax;
    Float * radial_bins_low;
    Float * radial_bins_high;
    
    // Variable to decide if we are using multiple tracers:
    bool multi_tracers;
    
    // Constructor
	Parameters(int argc, char *argv[]){
        
	    if (argc==1) usage();
	    int i=1;
	    while (i<argc) {
            if (!strcmp(argv[i],"-boxsize")){
                 // set cubic boxsize by default
                Float tmp_box=atof(argv[++i]);
                rect_boxsize = {tmp_box,tmp_box,tmp_box};
                }
        else if (!strcmp(argv[i],"-fname")) fname = argv[++i];
        else if (!strcmp(argv[i],"-fname2")) fname2 = argv[++i];
        else if (!strcmp(argv[i],"-rescale")) rescale = atof(argv[++i]);
	    else if (!strcmp(argv[i],"-R0")) R0 = atof(argv[++i]);
        else if (!strcmp(argv[i],"-nside")) nside = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-in")) fname = argv[++i];
        else if (!strcmp(argv[i],"-in2")) fname2 = argv[++i];
		else if (!strcmp(argv[i],"-rs")) rstart = atoi(argv[++i]);
		else if (!strcmp(argv[i],"-nmax")) nmax = atoll(argv[++i]);
		else if (!strcmp(argv[i],"-nthread")) nthread = atoi(argv[++i]);
		else if (!strcmp(argv[i],"-save")) savename = argv[++i];
		else if (!strcmp(argv[i],"-load")) loadname = argv[++i];
		else if (!strcmp(argv[i],"-balance")) qbalance = 1;
		else if (!strcmp(argv[i],"-invert")) qinvert = 1;
        else if (!strcmp(argv[i],"-output")) out_file = argv[++i];
        else if (!strcmp(argv[i],"-out_string")) out_string = argv[++i];
        else if (!strcmp(argv[i],"-binfile")) radial_bin_file=argv[++i];
        else if (!strcmp(argv[i],"-max_l")) max_l=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-inv_phi_file")) inv_phi_file=argv[++i];
        else if (!strcmp(argv[i],"-perbox")) perbox = 1;
        else if (!strcmp(argv[i],"-np")) {
			double tmp;
			if (sscanf(argv[++i],"%lf", &tmp)!=1) {
			    fprintf(stderr, "Failed to read number in %s %s\n",
			    	argv[i-1], argv[i]);
			    usage();
			}
			np = tmp;
			make_random=1;
		    }
		else if (!strcmp(argv[i],"-def")) { fname = NULL; }
		else {
		    fprintf(stderr, "Don't recognize %s\n", argv[i]);
		    usage();
		}
		i++;
	    }
	    
	        
#ifdef PERIODIC
        if (perbox!=true){
            printf("\nC++ code compiled with periodic flag, but periodic box parameter is not set! Exiting.\n\n");
            exit(1);
        }
#else
        if (perbox==true){
            printf("\nC++ code not compiled with periodic flag, but periodic box parameter is set! Exiting.\n\n");
            exit(1);
        }
#endif
        
        if(R0<40){
            printf("\nTruncation radius (%.0f Mpc/h) is too small for accurate power computation. Exiting.\n\n",R0);
            exit(1);
        }
        if(R0>400){
            printf("\nTruncation radius (%.0f Mpc/h) is too large for efficient power computation. Exiting.\n\n",R0);
            exit(1);
        }

	    // compute smallest and largest boxsizes
	    Float box_min = fmin(fmin(rect_boxsize.x,rect_boxsize.y),rect_boxsize.z);
	    Float box_max = fmax(fmax(rect_boxsize.x,rect_boxsize.y),rect_boxsize.z);
	    
	    assert(i==argc);  // For example, we might have omitted the last argument, causing disaster.

        assert(max_l%2==0); // check maximum ell is even
        assert(max_l<=10); // ell>10 not yet implemented!
        if (inv_phi_file==NULL) {inv_phi_file = (char *) default_inv_phi_file;} // no phi file specified
        assert(max_l%2==0); // check maximum ell is even
        assert(max_l<=10); // ell>10 not yet implemented!
        mbin = max_l/2+1; // number of angular bins is set to number of Legendre bins
        if (rescale<=0.0) rescale = box_max;   // This would allow a unit cube to fill the periodic volume
	    if (out_file==NULL) out_file = (char *) default_out_file; // no output savefile
	    if (out_string==NULL) out_string = (char *) default_out_string; // no output string
	    if (radial_bin_file==NULL) {radial_bin_file = (char *) default_radial_bin_file;} // No radial binning 
	    
	    if (fname==NULL) fname = (char *) default_fname;   // No name was given
	    if (fname2==NULL) fname2 = (char *) default_fname2;   // No name was given
	    
	    create_directory();
        
	    // Read in the radial binning
	    read_radial_binning(radial_bin_file);
        printf("Read in %d radial k-space bins in range (%.0f, %.0f) successfully.\n",nbin,rmin,rmax);
        
	    assert(box_min>0.0);
	    assert(rmax>0.0);
	    assert(nside>0);

#ifdef OPENMP
		omp_set_num_threads(nthread);
#else
		nthread=1;
#endif

		// Output for posterity
		printf("Box Size = {%6.5e,%6.5e,%6.5e}\n", rect_boxsize.x,rect_boxsize.y,rect_boxsize.z);
		printf("Grid = %d\n", nside);
		printf("Radial Bins = %d\n", nbin);
		printf("Radial k-space binning = {%6.5f, %6.5f} over %d bins (user-defined bin widths) \n",rmin,rmax,nbin);
		printf("Output directory: '%s'\n",out_file);

	}
private:
	void usage() {
	    fprintf(stderr, "\nUsage for grid_covariance:\n\n");
        fprintf(stderr, "   -def: This allows one to accept the defaults without giving other entries.\n");
	    fprintf(stderr, "   -in <file>: The input random particle file for particle-set 1 (space-separated x,y,z,w).\n");
	    fprintf(stderr, "   -in2 <file>: The input random particle file for particle-set 2 (space-separated x,y,z,w).\n");
        fprintf(stderr, "   -binfile <filename>: File containing the desired radial bins\n");
        fprintf(stderr, "   -output: (Pre-existing) directory to save output covariance matrices into\n"); 
        fprintf(stderr, "   -out_string: (Optional) String to add to file name to specify field type (e.g. RR)\n");
        fprintf(stderr, "   -nside <nside>: The grid size for accelerating the pair count.  Default 250.\n");
	    fprintf(stderr, "          There are {nside} cells along the longest dimension of the periodic box.\n");
	    fprintf(stderr, "   -nthread <nthread>: The number of CPU threads ot use for parallelization.\n");
        fprintf(stderr, "   -perbox <perbox>: Boolean, whether the box is periodic is not\n");
        fprintf(stderr, "\n");
	    fprintf(stderr, "   -max_l <max_l>: Maximum legendre multipole (must be even)\n");
        fprintf(stderr, "   -R0 <R0>: Truncation radius for pair-separation window (in Mpc/h)\n");
        fprintf(stderr, "   -inv_phi_file <filename>: Survey inverse correction function multipole coefficient file\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "   -boxsize <boxsize> : If creating particles randomly, this is the periodic size of the cubic computational domain.\n");
        fprintf(stderr, "           Default 400. If reading from file, this is reset dynamically creating a cuboidal box.\n");
	    fprintf(stderr, "   -rescale <rescale>: How much to dilate the input positions by.  Default 1.\n");
        fprintf(stderr, "            Zero or negative value causes =boxsize, rescaling unit cube to full periodicity\n");
        fprintf(stderr, "   -nmax <nmax>: The maximum number of particles to read in from the random particle files. Default 1000000000000\n");
	    fprintf(stderr, "   -save <filename>: Triggers option to store probability grid. <filename> has to end on \".bin\"\n");
	    fprintf(stderr, "      For advanced use, there is an option store the grid of probabilities used for sampling.\n");
	    fprintf(stderr, "      The file can then be reloaded on subsequent runs\n");
	    fprintf(stderr, "   -load <filename>: Triggers option to load the probability grid\n");
	    fprintf(stderr, "   -invert: Multiply all the weights by -1.\n");
	    fprintf(stderr, "   -balance: Rescale the negative weights so that the total weight is zero.\n");
        fprintf(stderr, "   -np <np>: Ignore any file and use np random perioidic points instead.\n");
        fprintf(stderr, "   -rs <rstart>:  If inverting particle weights, this sets the index from which to start weight inversion. Default 0\n");
	    fprintf(stderr, "\n");
	    fprintf(stderr, "\n");

	    exit(1);
	}
	
	void create_directory(){
        // Initialize output directory:
	    if (mkdir(out_file,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)==0){
            printf("\nCreating output directory\n");
        }
    }

    void read_radial_binning(char* binfile_name){
        // Read the radial binning file and determine the number of bins
        char line[100000];
    
        FILE *fp;
        fp = fopen(binfile_name,"r");
        if (fp==NULL){
            fprintf(stderr,"Radial binning file %s not found\n",binfile_name);
            abort();
        }
        fprintf(stderr,"\nReading radial binning file '%s'\n",binfile_name);
        
        // Count lines to construct the correct size
        while (fgets(line,10000,fp)!=NULL){
            if (line[0]=='#') continue; // comment line
            if (line[0]=='\n') continue;
                nbin++;
            }
            printf("\n# Found %d radial bins in the file\n",nbin);
            rewind(fp); // restart file
            
            // Now allocate memory to the weights array
            int ec=0;
            ec+=posix_memalign((void **) &radial_bins_low, PAGE, sizeof(Float)*nbin);
            ec+=posix_memalign((void **) &radial_bins_high, PAGE, sizeof(Float)*nbin);
            assert(ec==0);
            
            int line_count=0; // line counter
            int counter=0; // counts which element in line
            
            // Read in values to file
            while (fgets(line,100000,fp)!=NULL) {
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
                        radial_bins_low[line_count]=atof(split_string);
                        }
                    if(counter==1){
                        radial_bins_high[line_count]=atof(split_string);
                        }
                    if(counter>1){
                        fprintf(stderr,"Incorrect file format");
                        abort();
                    }
                    split_string = strtok(NULL,"\t");
                    counter++;
                }
                line_count++;
            }
            
            rmin = radial_bins_low[0];
            rmax = radial_bins_high[line_count-1];
            assert(line_count==nbin);
    }
};
#endif
