// parameter function file for grid_covariance.cpp (originally from Alex Wiegand)

#ifndef PARAMETERS_H
#define PARAMETERS_H

class Parameters{

public:
	// Important variables to set!  Here are the defaults:
	
    //---------- ESSENTIAL PARAMETERS -----------------
    
    // The name of the input random particle files (first set)
	char *fname = NULL;
	const char default_fname[500] = "/mnt/store1/oliverphilcox/Legendre2PCF/qpm_randoms_10x.xyzwj"; 
    
    // Name of the radial binning .csv file
    char *radial_bin_file = NULL;
    const char default_radial_bin_file[500] = "/mnt/store1/oliverphilcox/Legendre2PCF/radial_binning_cov.csv";
    
    // The name of the correlation function file for the first set of particles
	char *corname = NULL;
	const char default_corname[500] = "/mnt/store1/oliverphilcox/Legendre2PCF/mean_xi.xi";
    
    // Name of the correlation function radial binning .csv file
    char *radial_bin_file_cf = NULL;
    const char default_radial_bin_file_cf[500] = "/mnt/store1/oliverphilcox/Legendre2PCF/radial_binning_corr.csv";
    
    // Number of galaxies in first dataset
    Float nofznorm = 642051;
    
    // Output directory 
    char *out_file = NULL;
    const char default_out_file[500] = "/mnt/store1/oliverphilcox/Legendre2PCF/MeanXi/";
    
    // The number of mu bins in the correlation function
    int mbin_cf = 10;
    
    // The number of threads to run on
	int nthread = 20;

    // The grid size, which should be tuned to match boxsize and rmax. 
	// This uses the maximum width of the cuboidal box.
	int nside = 251;
    
    // Whether or not we are using a periodic box
	bool perbox = false;

    //---------- (r,mu) PARAMETERS ------------------------------------------
    
	// The number of mu bins
	int mbin = 30;
    
     // Name of the RR bin file
    char *RR_bin_file = NULL; // RR_{aA}^{11} file
    const char default_RR_bin_file[500] = "/mnt/store1/oliverphilcox/Legendre2PCF/RR_counts_n35_m30_11.txt";
    
    //---------- JACKKNIFE PARAMETERS ---------------------------------------
    
    // Name of the jackknife weight file
    char *jk_weight_file = NULL; // w_{aA}^{11} weights
    const char default_jk_weight_file[500] = "";
    
    //-------- LEGENDRE PARAMETERS -------------------------------------------
    
    int max_l = 10; // max Legendre moment (must be even)
    
    char *phi_file = NULL; // Survey correction function coefficient file 
    const char default_phi_file[500] = "/mnt/store1/oliverphilcox/Legendre2PCF/BinCorrectionFactor_n35_m1000_11.txt";
    
    //---------- PRECISION PARAMETERS ---------------------------------------
	
    // Maximum number of iterations to compute the C_ab integrals over
    int max_loops=20;
    
    // Number of random cells to draw at each stage
    int N2 = 40; // number of j cells per i cell
    int N3 = 40; // number of k cells per j cell
    int N4 = 40; // number of l cells per k cell
    
    //------------------ GENERAL MULTI-FIELD PARAMETERS ----------------------
    
    // Second set of random particles
    char *fname2 = NULL;
    const char default_fname2[500] = ""; 
    
    // Correlation functions
    char *corname2 = NULL; // xi_22 file
    const char default_corname2[500] = "";
    
    char *corname12 = NULL; // xi_12 file
    const char default_corname12[500] = "";
    
    // Number of galaxies in second dataset
    Float nofznorm2=0; // 
    
    //---------- (r,mu) MULTI-FIELD PARAMETERS ------------------------------
    
    // Summed pair count files 
    char *RR_bin_file12 = NULL; // RR_{aA}^{12} file
    const char default_RR_bin_file12[500] = "";
    
    char *RR_bin_file2 = NULL; // RR_{aA}^{22} file
    const char default_RR_bin_file2[500] = "";
    
    //-------- JACKKNIFE MULTI-FIELD PARAMETERS ------------------------------
    
    // Jackknife weight files
    char *jk_weight_file12 = NULL; // w_{aA}^{12} weights
    const char default_jk_weight_file12[500] = "";
    
    char *jk_weight_file2 = NULL; // w_{aA}^{22} weights
    const char default_jk_weight_file2[500] = "";
    
    //-------- LEGENDRE MULTI-FIELD PARAMETERS -------------------------------
    
    char *phi_file2 = NULL; // (Normalized) survey correction function survey_22
    const char default_phi_file2[500] = "";
    
    char *phi_file12 = NULL; // (Normalized) survey correction function survey_12
    const char default_phi_file12[500] = "";
    
    //-------- OTHER PARAMETERS ----------------------------------------------
    
	// The minimum mu of the smallest bin.
	Float mumin = 0.0;

	// The maximum mu of the largest bin.
	Float mumax = 1.0;
    
    // Number of loops over which to refine the correlation function
    int cf_loops = 20;
    
    // The periodicity of the position-space cube.
	Float boxsize = 200; // this is only used if the input particles are made randomly
    
	// The particles will be read from the unit cube, but then scaled by boxsize.
	Float rescale = 1.;   // If left zero or negative, set rescale=boxsize

	// The radius beyond which the correlation function is set to zero
	Float xicutoff = 400.0;
    
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
    
	// The periodicity of the position-space cuboid in 3D. 
    Float3 rect_boxsize = {boxsize,boxsize,boxsize}; // this is overwritten on particle read-in
    
    
    // Radial binning parameters (will be set from file)
    int nbin=0,nbin_cf=0;
    Float rmin, rmax, rmin_cf,rmax_cf;
    Float * radial_bins_low, * radial_bins_low_cf;
    Float * radial_bins_high, * radial_bins_high_cf;
    
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
        else if (!strcmp(argv[i],"-maxloops")) max_loops = atof(argv[++i]);
        else if (!strcmp(argv[i],"-rescale")) rescale = atof(argv[++i]);
		else if (!strcmp(argv[i],"-mumax")) mumax = atof(argv[++i]);
		else if (!strcmp(argv[i],"-mumin")) mumin = atof(argv[++i]);
        else if (!strcmp(argv[i],"-cf_loops")) cf_loops = atoi(argv[++i]);
		else if (!strcmp(argv[i],"-xicut")) xicutoff = atof(argv[++i]);
		else if (!strcmp(argv[i],"-norm")) nofznorm = atof(argv[++i]);
		else if (!strcmp(argv[i],"-norm2")) nofznorm2 = atof(argv[++i]);
        else if (!strcmp(argv[i],"-nside")) nside = atoi(argv[++i]);
		else if (!strcmp(argv[i],"-in")) fname = argv[++i];
        else if (!strcmp(argv[i],"-in2")) fname2 = argv[++i];
		else if (!strcmp(argv[i],"-cor")) corname = argv[++i];
		else if (!strcmp(argv[i],"-cor12")) corname12 = argv[++i];
		else if (!strcmp(argv[i],"-cor2")) corname2 = argv[++i];
        else if (!strcmp(argv[i],"-rs")) rstart = atoi(argv[++i]);
		else if (!strcmp(argv[i],"-nmax")) nmax = atoll(argv[++i]);
		else if (!strcmp(argv[i],"-nthread")) nthread = atoi(argv[++i]);
		else if (!strcmp(argv[i],"-mbin_cf")) mbin_cf = atoi(argv[++i]);
		else if (!strcmp(argv[i],"-save")) savename = argv[++i];
		else if (!strcmp(argv[i],"-load")) loadname = argv[++i];
		else if (!strcmp(argv[i],"-balance")) qbalance = 1;
		else if (!strcmp(argv[i],"-invert")) qinvert = 1;
        else if (!strcmp(argv[i],"-output")) out_file = argv[++i];
        else if (!strcmp(argv[i],"-binfile")) radial_bin_file=argv[++i];
        else if (!strcmp(argv[i],"-binfile_cf")) radial_bin_file_cf=argv[++i];
        else if (!strcmp(argv[i],"-N2")) N2=atof(argv[++i]);
        else if (!strcmp(argv[i],"-N3")) N3=atof(argv[++i]);
        else if (!strcmp(argv[i],"-N4")) N4=atof(argv[++i]);
#ifdef LEGENDRE
        else if (!strcmp(argv[i],"-max_l")) max_l=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-phi_file")) phi_file=argv[++i];
        else if (!strcmp(argv[i],"-phi_file12")) phi_file12=argv[++i];
        else if (!strcmp(argv[i],"-phi_file2")) phi_file2=argv[++i];
#elif defined JACKKNIFE
        else if (!strcmp(argv[i],"-jackknife")) jk_weight_file=argv[++i];
        else if (!strcmp(argv[i],"-jackknife12")) jk_weight_file12=argv[++i];
        else if (!strcmp(argv[i],"-jackknife2")) jk_weight_file2=argv[++i];
#else
		else if (!strcmp(argv[i],"-mbin")) mbin = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-RRbin")) RR_bin_file=argv[++i];
		else if (!strcmp(argv[i],"-RRbin12")) RR_bin_file12=argv[++i];
		else if (!strcmp(argv[i],"-RRbin2")) RR_bin_file2=argv[++i];
#endif
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

        if((!multi_tracers)&&(make_random==1)){
            printf("\nRunning for multiple tracers but creating particles at random; this is not yet supported. Exiting.\n\n");
            exit(1);
        }

	    // compute smallest and largest boxsizes
	    Float box_min = fmin(fmin(rect_boxsize.x,rect_boxsize.y),rect_boxsize.z);
	    Float box_max = fmax(fmax(rect_boxsize.x,rect_boxsize.y),rect_boxsize.z);
	    
	    assert(i==argc);  // For example, we might have omitted the last argument, causing disaster.

	    assert(nside%2!=0); // The probability integrator needs an odd grid size
	
	    assert(nofznorm>0); // need some galaxies!

	    assert(mumin>=0); // We take the absolte value of mu
	    assert(mumax<=1); // mu > 1 makes no sense
        
#ifdef LEGENDRE
        assert(max_l%2==0); // check maximum ell is even
        assert(max_l<=10); // ell>10 not yet implemented!
        if (phi_file==NULL) {phi_file = (char *) default_phi_file;} // no phi file specified
        if (phi_file2==NULL) {phi_file2 = (char *) default_phi_file2;}
        if (phi_file12==NULL) {phi_file12 = (char *) default_phi_file12;} 
        mbin = max_l/2+1; // number of angular bins is set to number of Legendre bins
#elif defined JACKKNIFE
	    if (jk_weight_file==NULL) jk_weight_file = (char *) default_jk_weight_file; // No jackknife name was given
	    if (jk_weight_file12==NULL) jk_weight_file12 = (char *) default_jk_weight_file12; // No jackknife name was given
	    if (jk_weight_file2==NULL) jk_weight_file2 = (char *) default_jk_weight_file2; // No jackknife name was given
#else
	    if (RR_bin_file==NULL) RR_bin_file = (char *) default_RR_bin_file; // no binning file was given
	    if (RR_bin_file12==NULL) RR_bin_file12 = (char *) default_RR_bin_file12; // no binning file was given
	    if (RR_bin_file2==NULL) RR_bin_file2 = (char *) default_RR_bin_file2; // no binning file was given
#endif
        if (rescale<=0.0) rescale = box_max;   // This would allow a unit cube to fill the periodic volume
	    if (corname==NULL) { corname = (char *) default_corname; }// No name was given
	    if (out_file==NULL) out_file = (char *) default_out_file; // no output savefile
	    if (radial_bin_file==NULL) {radial_bin_file = (char *) default_radial_bin_file;} // No radial binning 
	    if (radial_bin_file_cf==NULL) {radial_bin_file_cf = (char *) default_radial_bin_file_cf;} // No radial binning 
	    
	    if (fname==NULL) fname = (char *) default_fname;   // No name was given
	    if (fname2==NULL) fname2 = (char *) default_fname2;   // No name was given
	    if (corname2==NULL) { corname2 = (char *) default_corname2; }// No name was given
	    if (corname12==NULL) { corname12 = (char *) default_corname12; }// No name was given
	    
	    // Decide if we are using multiple tracers:
	    if (strlen(fname2)!=0){
#ifdef LEGENDRE
            if ((strlen(phi_file12)==0)||(strlen(phi_file2)==0)){
                printf("Two random particle sets input but not enough survey correction function files! Exiting.");
                exit(1);
            }
            else if ((strlen(corname2)==0)||(strlen(corname12)==0)){
                printf("Two random particle sets input but not enough correlation function files! Exiting.");
                exit(1);
            }
            else if (nofznorm2==0){
                printf("Two random particle sets input but only one galaxy number provided; exiting.");
                exit(1);
            }
            else{
                // Set multi tracers parameter since we have all required data
                printf("Using two sets of tracer particles");
                multi_tracers=true;
            }
        }
        else{
            printf("\nUsing a single set of tracer particles\n");
            multi_tracers=false;
            // set variables for later use
            phi_file12=phi_file;
            phi_file2=phi_file;
            nofznorm2=nofznorm;
            corname12=corname;
            corname2=corname;
        }            
#else
            if ((strlen(RR_bin_file12)==0)||(strlen(RR_bin_file2)==0)){
                printf("Two random particle sets input but not enough RR pair count files! Exiting.");
                exit(1);
            }
#ifdef JACKKNIFE
            else if ((strlen(jk_weight_file12)==0)||(strlen(jk_weight_file2)==0)){
                printf("Two random particle sets input but not enough jackknife weight files! Exiting.");
                exit(1);
            }
#endif
            else if ((strlen(corname2)==0)||(strlen(corname12)==0)){
                printf("Two random particle sets input but not enough correlation function files! Exiting.");
                exit(1);
            }
            else if (nofznorm2==0){
                printf("Two random particle sets input but only one galaxy number provided; exiting.");
                exit(1);
            }
            else{
                // Set multi_tracers parameter since we have all required data
                printf("Using two sets of tracer particles");
                multi_tracers=true;
            }
        }
        else{
            printf("\nUsing a single set of tracer particles\n");
            multi_tracers=false;
            // set variables for later use
            RR_bin_file12=RR_bin_file;
            RR_bin_file2=RR_bin_file;
            nofznorm2=nofznorm;
            corname12=corname;
            corname2=corname;
#ifdef JACKKNIFE
            jk_weight_file12=jk_weight_file;
            jk_weight_file2=jk_weight_file;
#endif
        }
#endif
	    
	    create_directory();
        
	    // Read in the radial binning
	    read_radial_binning(radial_bin_file);
        printf("Read in %d radial bins in range (%.0f, %.0f) successfully.\n",nbin,rmin,rmax);
        
        read_radial_binning_cf(radial_bin_file_cf);
        printf("Read in %d radial bins in range (%.0f, %.0f) successfully.\n",nbin_cf,rmin_cf,rmax_cf);
        
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
	    fprintf(stderr, "\nUsage for grid_covariance:\n\n");
        fprintf(stderr, "   -def: This allows one to accept the defaults without giving other entries.\n");
	    fprintf(stderr, "   -in <file>: The input random particle file for particle-set 1 (space-separated x,y,z,w).\n");
        fprintf(stderr, "   -binfile <filename>: File containing the desired radial bins\n");
        fprintf(stderr, "   -cor <file>: File location of input xi_1 correlation function file.\n");
	    fprintf(stderr, "   -binfile_cf <filename>: File containing the desired radial bins for the correlation function.\n");
        fprintf(stderr, "   -norm <nofznorm>: Number of galaxies in the first tracer set.\n");
        fprintf(stderr, "   -jackknife <filename>: File containing the {1,1} jackknife weights (normally computed from Corrfunc)\n");
        fprintf(stderr, "   -RRbin <filename>: File containing the {1,1} jackknife RR bin counts (computed from Corrfunc)\n");
        fprintf(stderr, "   -output: (Pre-existing) directory to save output covariance matrices into\n");   
        fprintf(stderr, "   -mbin <mbin>:  The number of mu bins (spaced linearly).\n");
	    fprintf(stderr, "   -mbin_cf <mbin_cf>:  The number of mu bins in the correlation function (spaced linearly).\n");
	    fprintf(stderr, "   -nside <nside>: The grid size for accelerating the pair count.  Default 250.\n");
	    fprintf(stderr, "          Recommend having several grid cells per rmax.\n");
        fprintf(stderr, "          There are {nside} cells along the longest dimension of the periodic box.\n");
	    fprintf(stderr, "   -nthread <nthread>: The number of CPU threads ot use for parallelization.\n");
        fprintf(stderr, "   -perbox <perbox>: Boolean, whether the box is periodic is not\n");
        fprintf(stderr, "\n");
        
	    fprintf(stderr, "   -in2 <file>: (Optional) The input random particle file for particle-set 2 (space-separated x,y,z,w).\n");
	    fprintf(stderr, "   -cor12 <file>: (Optional) File location of input xi_{12} cross-correlation function file.\n");
	    fprintf(stderr, "   -cor2 <file>: (Optional) File location of input xi_2 correlation function file.\n");
	    fprintf(stderr, "   -norm2 <nofznorm2>: (Optional) Number of galaxies in the survey for the second tracer set.\n");
        
#ifdef LEGENDRE
        fprintf(stderr, "   -max_l <max_l>: Maximum legendre multipole (must be even)");
        fprintf(stderr, "   -phi_file <filename>: Survey correction function coefficient file\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "   -phi_file12 <filename>: (Optional) Survey correction function coefficient file for fields 1 x 2\n");
        fprintf(stderr, "   -phi_file2 <filename>: (Optional) Survey correction function coefficent file for field 2\n");
        fprintf(stderr, "\n");
#else
        fprintf(stderr, "   -RRbin12 <filename>: (Optional) File containing the {1,2} jackknife RR bin counts (computed from Corrfunc)\n");
	    fprintf(stderr, "   -RRbin2 <filename>: (Optional) File containing the {2,2} jackknife RR bin counts (computed from Corrfunc)\n");
        fprintf(stderr, "\n");
#ifdef JACKKNIFE
        fprintf(stderr, "   -jackknife12 <filename>: (Optional) File containing the {1,2} jackknife weights (normally computed from Corrfunc)\n");
        fprintf(stderr, "   -jackknife2 <filename>: (Optional) File containing the {2,2} jackknife weights (normally computed from Corrfunc)\n");
#endif
#endif
        fprintf(stderr, "   -maxloops <max_loops>: Maximum number of integral loops\n");
        fprintf(stderr, "   -N2 <N2>: Number of secondary particles to choose per primary particle\n");
        fprintf(stderr, "   -N3 <N3>: Number of tertiary particles to choose per secondary particle\n");
        fprintf(stderr, "   -N4 <N4>: Number of quaternary particles to choose per tertiary particle\n");
        fprintf(stderr, "\n");
        
        fprintf(stderr, "   -mumin <mumin> : Minimum mu binning to use.\n");
        fprintf(stderr, "   -mumax <mumax> : Maximum mu binning to use.\n");
        fprintf(stderr, "   -cf_loops <cf_loops>: Number of iterations over which to refine the correlation functions.\n");
        fprintf(stderr, "   -boxsize <boxsize> : If creating particles randomly, this is the periodic size of the cubic computational domain.\n");
        fprintf(stderr, "           Default 400. If reading from file, this is reset dynamically creating a cuboidal box.\n");
	    fprintf(stderr, "   -rescale <rescale>: How much to dilate the input positions by.  Default 1.\n");
        fprintf(stderr, "            Zero or negative value causes =boxsize, rescaling unit cube to full periodicity\n");
	    fprintf(stderr, "   -xicut <xicutoff>: The radius beyond which xi is set to zero.  Default 400.\n");
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
	    // First create whole directory if it doesn't exist:
	    if (mkdir(out_file,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)==0){
            printf("\nCreating output directory\n");
        }
	    char cname[1000];
        snprintf(cname, sizeof cname, "%sCovMatricesAll/",out_file);
	    if (mkdir(cname,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0){
            }
        // Check if this was successful:
        struct stat info;
        
        if( stat( cname, &info ) != 0 ){
            printf( "\nCreation of directory %s failed\n", cname);
            exit(1);
        }
#ifdef JACKKNIFE
        char cjname[1000];
        snprintf(cjname, sizeof cjname, "%sCovMatricesJack/",out_file);
        if (mkdir(cjname,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)!=0){
        }

        if(stat(cjname,&info)!=0){
            printf("\nCreation of directory %s fiailed\n",cjname);
            exit(1);
        }
#endif
    }
    
    void read_radial_binning_cf(char* binfile_name){
        // Read the radial binning file for a correlation function
        char line[100000];
    
        FILE *fp;
        fp = fopen(binfile_name,"r");
        if (fp==NULL){
            fprintf(stderr,"Radial correlation function binning file %s not found\n",binfile_name);
            abort();
        }
        fprintf(stderr,"\nReading radial correlation function binning file '%s'\n",binfile_name);
        
        // Count lines to construct the correct size
        while (fgets(line,10000,fp)!=NULL){
            if (line[0]=='#') continue; // comment line
            if (line[0]=='\n') continue;
                nbin_cf++;
            }
            printf("\n# Found %d radial bins in the correlation function binning file\n",nbin_cf);
            rewind(fp); // restart file
            
            // Now allocate memory to the weights array
            int ec=0;
            ec+=posix_memalign((void **) &radial_bins_low_cf, PAGE, sizeof(Float)*nbin_cf);
            ec+=posix_memalign((void **) &radial_bins_high_cf, PAGE, sizeof(Float)*nbin_cf);
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
                        radial_bins_low_cf[line_count]=atof(split_string);
                        }
                    if(counter==1){
                        radial_bins_high_cf[line_count]=atof(split_string);
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
            
            rmin_cf = radial_bins_low_cf[0];
            rmax_cf = radial_bins_high_cf[line_count-1];
            assert(line_count==nbin_cf);
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
