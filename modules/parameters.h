// parameter function file for grid_covariance.cpp (originally from Alex Wiegand)

#ifndef PARAMETERS_H
#define PARAMETERS_H

class Parameters{

public:
	// Important variables to set!  Here are the defaults:
	// The periodicity of the position-space cube.
	Float boxsize = 400;
    
    // The periodicity of the position-space cuboid in 3D. 
    Float3 rect_boxsize = {400,400,400};
    
	// The particles will be read from the unit cube, but then scaled by boxsize.
	Float rescale = 1.;   // If left zero or negative, set rescale=boxsize

	// The maximum mu of the largest bin.
	Float mumax = 1.0;

	// The minimum mu of the smallest bin.
	Float mumin = 0.0;

	// The radius beyond which the correlation function is set to zero
	Float xicutoff = 400.0;
    
    // Cut-off radius below which correlation function is set to zero:
    Float r_cutoff = 0.;    

	Float nofznorm=1198006; //1198006 - for BOSS DR12

	// The grid size, which should be tuned to match boxsize and rmax. 
	// This uses the maximum width of the cuboidal box.
	// Don't forget to adjust this if changing boxsize!
	int nside = 301;

	// If set, we'll just throw random periodic points instead of reading the file
	int make_random = 0;

	// Will be number of particles in a random distribution,
	// but gets overwritten if reading from a file.
	int np = -1;

	// Whether to balance the weights or multiply them by -1
	int qbalance = 0, qinvert = 0;

	// The maximum number of points to read
	uint64 nmax = 1000000000000;

    // The number of mu bins
	int mbin = 10;
    
    // The number of threads to run on
	int nthread=10;

	// The location and name of a integrated grid of probabilities to be saved
	char *savename = NULL;
	// The location and name of a integrated grid of probabilities to be loaded
	char *loadname = NULL; //	"ProbListGrid100box10-corfu.dat"

	// The index from which on to invert the sign of the weights
	int rstart=0;

	// Whether or not we are using a periodic box
	bool perbox=false;

	// The name of the input file
	char *fname = NULL;
	const char default_fname[500] = "../random_particles/test_particles.txt";//test_particles_mid.txt";//test_particles_small.txt";

	// The name of the correlation function file
	char *corname = NULL;
	const char default_corname[500] = "xi_functions/QPM_Mash.xi";
    
    // Name of the radial binning .csv file
    char *radial_bin_file = NULL;
    const char default_radial_bin_file[500] = "python/binfile_linear.csv";
    
    // Name of the jackknife weight file
    char *jk_weight_file = NULL;
    const char default_jk_weight_file[500] = "weight_files/jackknife_weights_n36_m10_j169.dat";//169.dat";
    
    // Name of the RR bin file
    char *RR_bin_file = NULL;
    const char default_RR_bin_file[500] = "weight_files/binned_pair_counts_n36_m10_j169.dat";//;169.dat";
    
    // Maximum number of iterations to compute the C_ab integrals over
    int max_loops=10;//10; 
    int N2 = 5;//20; // number of j cells per i cell
    int N3 = 10;//25; // number of k cells per j cell
    int N4 = 0;//50; // number of l cells per k cell
    
    // Radial binning parameters (will be set from file)
    int nbin=0;
    Float rmin, rmax;
    Float * radial_bins_low;
    Float * radial_bins_high;
    
    char *out_file = NULL;
    const char default_out_file[500] = "/mnt/store1/oliverphilcox/DR12/";
    
	// Constructor
	Parameters(int argc, char *argv[]){

	    if (argc==1) usage();
	    int i=1;
	    while (i<argc) {
		     if (!strcmp(argv[i],"-boxsize")||!strcmp(argv[i],"-box")){
                 // set cubic boxsize by default
                Float tmp_box=atof(argv[++i]);
                rect_boxsize = {tmp_box,tmp_box,tmp_box};
                }
        else if (!strcmp(argv[i],"-maxloops")||!strcmp(argv[i],"-loops")) max_loops = atof(argv[++i]);
        else if (!strcmp(argv[i],"-rescale")||!strcmp(argv[i],"-scale")) rescale = atof(argv[++i]);
		else if (!strcmp(argv[i],"-mumax")) mumax = atof(argv[++i]);
		else if (!strcmp(argv[i],"-mumin")) mumin = atof(argv[++i]);
		else if (!strcmp(argv[i],"-xicut")) xicutoff = atof(argv[++i]);
		else if (!strcmp(argv[i],"-norm")) nofznorm = atof(argv[++i]);
		else if (!strcmp(argv[i],"-nside")||!strcmp(argv[i],"-ngrid")||!strcmp(argv[i],"-grid")) nside = atoi(argv[++i]);
		else if (!strcmp(argv[i],"-in")) fname = argv[++i];
		else if (!strcmp(argv[i],"-cor")) corname = argv[++i];
		else if (!strcmp(argv[i],"-rs")) rstart = atoi(argv[++i]);
		else if (!strcmp(argv[i],"-nmax")) nmax = atoll(argv[++i]);
		else if (!strcmp(argv[i],"-nthread")||!strcmp(argv[i],"-nthreads")) nthread = atoi(argv[++i]);
		else if (!strcmp(argv[i],"-mbin")) mbin = atoi(argv[++i]);
		else if (!strcmp(argv[i],"-save")||!strcmp(argv[i],"-store")) savename = argv[++i];
		else if (!strcmp(argv[i],"-load")) loadname = argv[++i];
		else if (!strcmp(argv[i],"-balance")) qbalance = 1;
		else if (!strcmp(argv[i],"-invert")) qinvert = 1;
        else if (!strcmp(argv[i],"-output")) out_file = argv[++i];
        else if (!strcmp(argv[i],"-binfile")||!strcmp(argv[i],"-radial_binfile")) radial_bin_file=argv[++i];
        else if (!strcmp(argv[i],"-jackknife")||!strcmp(argv[i],"-jackknife_weights")||!strcmp(argv[i],"-weights")) jk_weight_file=argv[++i];
        else if (!strcmp(argv[i],"-RRbin")||!strcmp(argv[i],"-RRbin_file")) RR_bin_file=argv[++i];
		else if (!strcmp(argv[i],"-ran")||!strcmp(argv[i],"-np")) {
			double tmp;
			if (sscanf(argv[++i],"%lf", &tmp)!=1) {
			    fprintf(stderr, "Failed to read number in %s %s\n",
			    	argv[i-1], argv[i]);
			    usage();
			}
			np = tmp;
			make_random=1;
		    }
		else if (!strcmp(argv[i],"-def")||!strcmp(argv[i],"-default")) { fname = NULL; }
		else {
		    fprintf(stderr, "Don't recognize %s\n", argv[i]);
		    usage();
		}
		i++;
	    }
	    // compute smallest and largest boxsizes
	    Float box_min = fmin(fmin(rect_boxsize.x,rect_boxsize.y),rect_boxsize.z);
	    Float box_max = fmax(fmax(rect_boxsize.x,rect_boxsize.y),rect_boxsize.z);
	    
	    assert(i==argc);  // For example, we might have omitted the last argument, causing disaster.

	    assert(nside%2!=0); // The probability integrator needs an odd grid size

	    assert(mumin>=0); // We take the absolte value of mu
	    assert(mumax<=1); // mu > 1 makes no sense
        
        if (rescale<=0.0) rescale = box_max;   // This would allow a unit cube to fill the periodic volume
	    if (fname==NULL) fname = (char *) default_fname;   // No name was given
	    if (jk_weight_file==NULL) jk_weight_file = (char *) default_jk_weight_file; // No jackknife name was given
	    if (RR_bin_file==NULL) RR_bin_file = (char *) default_RR_bin_file; // no binning file was given
	    if (corname==NULL) { corname = (char *) default_corname; }// No name was given
	    if (radial_bin_file==NULL) {radial_bin_file = (char *) default_radial_bin_file;} // No radial binning 
	    if (out_file==NULL) out_file = (char *) default_out_file; // no output savefile
	    
	    create_directory();
        
	    // Read in the radial binning
	    //TODO: Initialize nbin, rmin, rmax somewhere - not an input
        read_radial_binning(radial_bin_file);
        printf("Read in %d radial bins in range (%.0f, %.0f) successfully.\n",nbin,rmin,rmax);
        
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
		printf("Maximum Radius = %6.5e\n", rmax);
		Float gridsize = rmax/(box_max/nside);
		printf("Max Radius in Grid Units = %6.5e\n", gridsize);
		if (gridsize<1) printf("#\n# WARNING: grid appears inefficiently coarse\n#\n");
        printf("Cutoff radius = %.2f\n", r_cutoff);
		printf("Radial Bins = %d\n", nbin);
		printf("Radial Binning = {%6.5f, %6.5f} over %d bins (user-defined bin widths) \n",rmin,rmax,nbin);
		printf("Mu Bins = %d\n", mbin);
		printf("Mu Binning = {%6.5f, %6.5f, %6.5f}\n",mumin,mumax,(mumax-mumin)/mbin);
		printf("Density Normalization = %6.5e\n",nofznorm);
        printf("Maximum number of integration loops = %d\n",max_loops);
        printf("Output directory: '%s'\n",out_file);

	}
private:
	void usage() {
	    fprintf(stderr, "\nUsage for grid_covariance:\n");
	    fprintf(stderr, "   -box <boxsize> : If creating particles randomly, this is the periodic size of the cubic computational domain.  Default 400. If reading from file, this is reset dynamically creating a cuboidal box.\n");
	    fprintf(stderr, "   -scale <rescale>: How much to dilate the input positions by.  Default 0.\n");
	    fprintf(stderr, "             Zero or negative value causes =boxsize, rescaling unit cube to full periodicity\n");
	    fprintf(stderr, "   -xicut <xicutoff>: The radius beyond which xi is set to zero.  Default 1000.\n");
	    fprintf(stderr, "   -nside <nside>: The grid size for accelerating the pair count.  Default 50.\n");
	    fprintf(stderr, "             Recommend having several grid cells per rmax. There are {nside} cells along the longest dimension of the periodic box.\n");
	    fprintf(stderr, "   -in <file>: The input file (space-separated x,y,z,w).  Default sample.dat.\n");
	    fprintf(stderr, "   -ran <np>: Ignore any file and use np random perioidic points instead.\n");
	    fprintf(stderr, "   -norm <nofznorm>: Number of galaxies in the survey. Used to normalize n(z).\n");
	    fprintf(stderr, "   -def: This allows one to accept the defaults without giving other entries.\n");
	    fprintf(stderr, "\n");
        fprintf(stderr, "   -cor <file>: File location of inut correlation function file.\n");
	    fprintf(stderr, "   -mbin:  The number of mu bins (spaced linearly).\n");
        fprintf(stderr, "   -mumin / -mumax: Minimum / maximum mu binning to use.\n");
	    fprintf(stderr, "\n");
	    fprintf(stderr, "For advanced use, there is an option store the grid of probabilities used for sampling.\n");
	    fprintf(stderr, "    -save <filename>: Triggers option to store probability grid. <filename> has to end on \".bin\"\n");
	    fprintf(stderr, "The file can then be reloaded on subsequent runs\n");
	    fprintf(stderr, "    -load <filename>: Triggers option to load the probability grid\n");
        fprintf(stderr, "    -binfile <filename>: File containing the desired radial bins\n");
        fprintf(stderr, "    -weights <filename>: File containing the jackknife weights (normally computed from Corrfunc\n)");
        fprintf(stderr, "    -RRbin_file <filename>: File containing the jackknife RR bin counts (computed from Corrfunc\n)");
	    fprintf(stderr, "    -balance: Rescale the negative weights so that the total weight is zero.\n");
	    fprintf(stderr, "    -invert: Multiply all the weights by -1.\n");
        fprintf(stderr, "    -loops: Maximum number of integral loops\n");
        fprintf(stderr, "    -output: (Pre-existing) directory to save output covariance matrices into\n");   


	    exit(1);
	}
	
	void create_directory(){
        // Initialize output directory:
	    // First create whole directory if it doesn't exist:
	    if (mkdir(out_file,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)==0){
            printf("\nCreating output directory\n");
        }
	    char cname[1000],cjname[1000];
        snprintf(cname, sizeof cname, "%sCovMatricesAll/",out_file);
        snprintf(cjname, sizeof cjname, "%sCovMatricesJack/",out_file);
	    if (mkdir(cname,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0){
            }
        if (mkdir(cjname,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)!=0){
        }
    }
	    

    void read_radial_binning(char* binfile_name){//, int nbin, Float min_r, Float max_r){
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
