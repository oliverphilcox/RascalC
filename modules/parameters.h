
// parameter function file for grid_covariance.cpp (originally from Alex Wiegand)

#ifndef PARAMETERS_H
#define PARAMETERS_H

#ifdef LEGENDRE_MIX
#include "legendre_mix_utilities.h"
#endif

class Parameters{

public:
	// Important variables to set!  Here are the defaults:

    //---------- ESSENTIAL PARAMETERS -----------------

    // The name of the input random particle files (first set)
	char *fname = NULL;
	const char default_fname[500] = "/mnt/store1/oliverphilcox/Yuting/LRG_randoms_v2_10x.txt";

    // Name of the radial binning .csv file
    char *radial_bin_file = NULL;
    const char default_radial_bin_file[500] = "/home/oliverphilcox/eBOSS_MockChallenge/radial_binning_cov.csv";

    // The name of the correlation function file for the first set of particles
	char *corname = NULL;
	const char default_corname[500] = "/home/oliverphilcox/eBOSS_MockChallenge/v3_low/xi_n100_m10_periodic_11.dat";

    // Name of the correlation function radial binning .csv file
    char *radial_bin_file_cf = NULL;
    const char default_radial_bin_file_cf[500] = "/home/oliverphilcox/eBOSS_MockChallenge/v3_low/radial_binning_corr_low.csv";

    // Number of galaxies in first dataset
    Float nofznorm =156800;

    // Output directory
    char *out_file = NULL;
    const char default_out_file[500] = "/home/oliverphilcox/eBOSS_MockChallenge/v7/";

    // The number of mu bins in the correlation function
    int mbin_cf = 10;

    // The number of threads to run on
	int nthread = 30;

    // The grid size, which should be tuned to match boxsize and rmax.
	// This uses the maximum width of the cuboidal box.
	int nside = 71;

    // Whether or not we are using a periodic box
	bool perbox = false;

    //---------- (r,mu) PARAMETERS ------------------------------------------

	// The number of mu bins
	int mbin = 20;

     // Name of the RR bin file
    char *RR_bin_file = NULL; // RR_{aA}^{11} file
    const char default_RR_bin_file[500] = "";

    //---------- JACKKNIFE PARAMETERS ---------------------------------------

    // Name of the jackknife weight file
    char *jk_weight_file = NULL; // w_{aA}^{11} weights
    const char default_jk_weight_file[500] = "";

    //-------- LEGENDRE PARAMETERS -------------------------------------------

#if (defined LEGENDRE || defined LEGENDRE_MIX)
    int max_l = 2; // max Legendre moment (must be even unless computing 3PCF)
#endif

#ifdef LEGENDRE_MIX
    char *mu_bin_legendre_file = NULL; // Mu bin Legendre factors file
    const char default_mu_bin_legendre_file[500] = "mu_bin_legendre_factors_m20_l2.txt";

    MuBinLegendreFactors mu_bin_legendre_factors;
#endif

#ifdef LEGENDRE
    char *phi_file = NULL; // Survey correction function coefficient file
    const char default_phi_file[500] = "/home/oliverphilcox/eBOSS_MockChallenge/BinCorrectionFactor_n25_periodic_11.txt";
#endif

    //-------- POWER PARAMETERS (not yet publicly released) ------------------

    Float R0 = 100; // truncation radius in Mpc/h
    Float power_norm = 4.; // normalization of the power spectrum

    char *inv_phi_file = NULL; // Inverse survey correction function multipole coefficient file
    const char default_inv_phi_file[500] = "";

    //---------- PRECISION PARAMETERS ---------------------------------------

    // Maximum number of iterations to compute the C_ab integrals over
    int max_loops = 40;
    // Number of loops to output into each subsample/file
    int loops_per_sample = 1;
    // Number of output subsamples/files
    int no_subsamples = 120;

    // Number of random cells to draw at each stage
    int N2 = 20; // number of j cells per i cell
    int N3 = 40; // number of k cells per j cell
    int N4 = 80; // number of l cells per k cell

    //------------------ EXTRA 3PCF AUTOCOVARIANCE PARAMETERS ----------------

    int N5 = 10; // number of m cells per l cell
    int N6 = 10; // number of n cells per m cell

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
    Float nofznorm2=3398430; //

    //---------- (r,mu) MULTI-FIELD PARAMETERS ------------------------------

#if (!defined LEGENDRE && !defined POWER && !defined THREE_PCF)
    // Summed pair count files
    char *RR_bin_file12 = NULL; // RR_{aA}^{12} file
    const char default_RR_bin_file12[500] = "";

    char *RR_bin_file2 = NULL; // RR_{aA}^{22} file
    const char default_RR_bin_file2[500] = "";
#endif

    //-------- JACKKNIFE MULTI-FIELD PARAMETERS ------------------------------

#ifdef JACKKNIFE
    // Jackknife weight files
    char *jk_weight_file12 = NULL; // w_{aA}^{12} weights
    const char default_jk_weight_file12[500] = "";

    char *jk_weight_file2 = NULL; // w_{aA}^{22} weights
    const char default_jk_weight_file2[500] = "";
#endif

    //-------- LEGENDRE MULTI-FIELD PARAMETERS -------------------------------

#ifdef LEGENDRE
    const char default_phi_file12[500] = "";
    char *phi_file12 = NULL; // (Normalized) survey correction function survey_12

    char *phi_file2 = NULL; // (Normalized) survey correction function survey_22
    const char default_phi_file2[500] = "";
#endif

    // ------- POWER MULTI-FIELD PARAMETERS ----------------------------------

#ifdef POWER

    char *inv_phi_file2 = NULL; // (Normalized) inverse survey correction function multipoles survey_22
    const char default_inv_phi_file2[500] = "";

    char *inv_phi_file12 = NULL; // (Normalized) inverse survey correction function multipoles survey_12
    const char default_inv_phi_file12[500] = "";

    Float power_norm12 = 0; // power spectrum normalization for field 1 x 2

    Float power_norm2 = 0; // power spectrum normalization for field 2 x 2

#endif

    //-------- OTHER PARAMETERS ----------------------------------------------

	// The minimum mu of the smallest bin.
#ifdef THREE_PCF
    Float mumin = -1.0;
#else
	Float mumin = 0.0;
#endif

	// The maximum mu of the largest bin.
	Float mumax = 1.0;

    // Number of loops over which to refine the correlation function
    int cf_loops = 10;

    // The periodicity of the position-space cube.
	Float boxsize = 200.; // this is only used if the input particles are made randomly

	// The particles will be read from the unit cube, but then scaled by boxsize.
	Float rescale = 1.;   // If left zero or negative, set rescale=boxsize

	// The radius beyond which the correlation function is set to zero
	Float xicutoff = 250.;

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
	int np = 3000000; // NB: This is only used for grid creation so we don't need a separate variable for the second set of randoms

	// The index from which on to invert the sign of the weights
	int rstart = 0;

	//---------------- INTERNAL PARAMETERS -----------------------------------
    // (no more user defined parameters below this line)

	// The periodicity of the position-space cuboid in 3D.
    Float3 rect_boxsize = {boxsize,boxsize,boxsize}; // this is overwritten on particle read-in

    Float cellsize;

    // Radial binning parameters (will be set from file)
    int nbin=0,nbin_cf=0;
    Float rmin, rmax, rmin_cf,rmax_cf;
    Float * radial_bins_low, * radial_bins_low_cf;
    Float * radial_bins_high, * radial_bins_high_cf;

    // Variable to decide if we are using multiple tracers:
    bool multi_tracers;

    // Options for a deterministic/reproducible run
    bool random_seed = true;
    unsigned long seed = 0;

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
        else if (!strcmp(argv[i],"-maxloops")) max_loops = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-loopspersample")) loops_per_sample = atoi(argv[++i]);
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
#if (defined LEGENDRE || defined LEGENDRE_MIX)
        else if (!strcmp(argv[i],"-max_l")) max_l=atoi(argv[++i]);
#endif
#ifdef JACKKNIFE
        else if (!strcmp(argv[i],"-jackknife")) jk_weight_file=argv[++i];
        else if (!strcmp(argv[i],"-jackknife12")) jk_weight_file12=argv[++i];
        else if (!strcmp(argv[i],"-jackknife2")) jk_weight_file2=argv[++i];
#endif
#ifdef LEGENDRE_MIX
        else if (!strcmp(argv[i],"-mu_bin_legendre_file")) mu_bin_legendre_file=argv[++i];
#elif defined LEGENDRE
        else if (!strcmp(argv[i],"-phi_file")) phi_file=argv[++i];
        else if (!strcmp(argv[i],"-phi_file12")) phi_file12=argv[++i];
        else if (!strcmp(argv[i],"-phi_file2")) phi_file2=argv[++i];
#elif defined POWER
        else if (!strcmp(argv[i],"-max_l")) max_l=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-inv_phi_file")) inv_phi_file=argv[++i];
        else if (!strcmp(argv[i],"-inv_phi_file12")) inv_phi_file12=argv[++i];
        else if (!strcmp(argv[i],"-inv_phi_file2")) inv_phi_file2=argv[++i];
        else if (!strcmp(argv[i],"-R0")) R0 = atof(argv[++i]);
        else if (!strcmp(argv[i],"-power_norm")) power_norm = atof(argv[++i]);
        else if (!strcmp(argv[i],"-power_norm12")) power_norm12 = atof(argv[++i]);
        else if (!strcmp(argv[i],"-power_norm")) power_norm2 = atof(argv[++i]);
#elif defined THREE_PCF
        else if (!strcmp(argv[i],"-max_l")) max_l=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-phi_file")) phi_file=argv[++i];
        else if (!strcmp(argv[i],"-N5")) N5=atof(argv[++i]);
        else if (!strcmp(argv[i],"-N6")) N6=atof(argv[++i]);
#endif
#if (!defined LEGENDRE && !defined POWER && !defined THREE_PCF)
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
		else if (!strcmp(argv[i],"-seed")) { random_seed = false; sscanf(argv[++i], "%lu", &seed); }
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
	    // compute smallest and largest boxsizes
	    Float box_min = fmin(fmin(rect_boxsize.x,rect_boxsize.y),rect_boxsize.z);
	    Float box_max = fmax(fmax(rect_boxsize.x,rect_boxsize.y),rect_boxsize.z);

	    assert(i==argc);  // For example, we might have omitted the last argument, causing disaster.

	    assert(nside%2!=0); // The probability integrator needs an odd grid size

	    assert(nofznorm>0); // need some galaxies!
        assert(max_loops % loops_per_sample == 0); // group size need to divide the number of loops
        no_subsamples = max_loops / loops_per_sample;
#ifndef THREE_PCF
	    //assert(mumin>=0); // We take the absolte value of mu
#endif
	    assert(mumax<=1); // mu > 1 makes no sense

#if (defined LEGENDRE || defined LEGENDRE_MIX)
        assert(max_l%2==0); // check maximum ell is even
#endif
#ifdef LEGENDRE_MIX
        if (mu_bin_legendre_file == NULL) mu_bin_legendre_file = (char *) default_mu_bin_legendre_file; // no mu bin Legendre file specified
        new (&mu_bin_legendre_factors) MuBinLegendreFactors(mu_bin_legendre_file, mbin, max_l); // construct in place
#elif defined LEGENDRE
        assert(max_l<=10); // ell>10 not yet implemented!
        mbin = max_l/2+1; // number of angular bins is set to number of Legendre bins
        if (phi_file==NULL) {phi_file = (char *) default_phi_file;} // no phi file specified
        if (phi_file2==NULL) {phi_file2 = (char *) default_phi_file2;}
        if (phi_file12==NULL) {phi_file12 = (char *) default_phi_file12;}
#elif defined POWER
        assert(max_l%2==0); // check maximum ell is even
        assert(max_l<=10); // ell>10 not yet implemented!
        if (inv_phi_file==NULL) {inv_phi_file = (char *) default_inv_phi_file;} // no phi file specified
        if (inv_phi_file2==NULL) {inv_phi_file2 = (char *) default_inv_phi_file2;}
        if (inv_phi_file12==NULL) {inv_phi_file12 = (char *) default_inv_phi_file12;}
        mbin = max_l/2+1; // number of angular bins is set to number of Legendre bins
        if(R0<40){
            printf("\nTruncation radius (%.0f Mpc/h) is too small for accurate power computation. Exiting.\n\n",R0);
            exit(1);
        }
        if(R0>400){
            printf("\nTruncation radius (%.0f Mpc/h) is too large for efficient power computation. Exiting.\n\n",R0);
            exit(1);
        }
#elif defined THREE_PCF
        assert(max_l<=10); // ell>10 not yet implemented!
        if (phi_file==NULL) {phi_file = (char *) default_phi_file;} // no phi file specified
        mbin = max_l+1; // number of angular bins is set to number of Legendre bins (including odd bins)
#elif defined JACKKNIFE
	    if (jk_weight_file==NULL) jk_weight_file = (char *) default_jk_weight_file; // No jackknife name was given
	    if (jk_weight_file12==NULL) jk_weight_file12 = (char *) default_jk_weight_file12; // No jackknife name was given
	    if (jk_weight_file2==NULL) jk_weight_file2 = (char *) default_jk_weight_file2; // No jackknife name was given

	    if (RR_bin_file==NULL) RR_bin_file = (char *) default_RR_bin_file; // no binning file was given
	    if (RR_bin_file12==NULL) RR_bin_file12 = (char *) default_RR_bin_file12; // no binning file was given
	    if (RR_bin_file2==NULL) RR_bin_file2 = (char *) default_RR_bin_file2; // no binning file was given
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
#if (defined LEGENDRE || defined POWER)
#ifdef LEGENDRE
            if ((strlen(phi_file12)==0)||(strlen(phi_file2)==0)){
                printf("Two random particle sets input but not enough survey correction function files! Exiting.");
                exit(1);
            }
#else
            if ((strlen(inv_phi_file12)==0)||(strlen(inv_phi_file2)==0)){
                printf("Two random particle sets input but not enough survey correction function files! Exiting.");
                exit(1);
            }
            else if ((power_norm12==0)||(power_norm2==0)){
                printf("Two random particle sets input but not enough power normalizations provided; exiting.");
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
                // Set multi tracers parameter since we have all required data
                printf("Using two sets of tracer particles");
                multi_tracers=true;
            }
        }
        else{
            printf("\nUsing a single set of tracer particles\n");
            multi_tracers=false;
            // set variables for later use
#ifdef LEGENDRE
            phi_file12=phi_file;
            phi_file2=phi_file;
#else
            inv_phi_file12=inv_phi_file;
            inv_phi_file2=inv_phi_file;
            power_norm2 = power_norm;
            power_norm12 = power_norm;
#endif
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

#ifdef THREE_PCF
        if(multi_tracers){
            printf("\nSupport for multi-tracer threePCF covariance matrices not yet available. Exiting.\n\n");
            exit(1);
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
		printf("Grid = %d\n", nside);
		printf("Maximum Radius = %6.5e\n", rmax);
        printf("Radial Bins = %d\n", nbin);
		printf("Radial Binning = {%6.5f, %6.5f} over %d bins (user-defined bin widths) \n",rmin,rmax,nbin);
#if (!defined LEGENDRE && !defined THREE_PCF && !defined POWER)
		printf("Mu Bins = %d\n", mbin);
		printf("Mu Binning = {%6.5f, %6.5f, %6.5f}\n",mumin,mumax,(mumax-mumin)/mbin);
#endif
		printf("Number of galaxies = %6.5e\n",nofznorm);
        printf("Maximum number of integration loops = %d\n",max_loops);
        printf("Number of output subsamples = %d\n", no_subsamples);
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
#ifdef JACKKNIFE
        fprintf(stderr, "   -jackknife <filename>: File containing the {1,1} jackknife weights (normally computed from Corrfunc)\n");
#endif
#if (!defined LEGENDRE && !defined POWER && !defined THREE_PCF)
        fprintf(stderr, "   -RRbin <filename>: File containing the {1,1} jackknife RR bin counts (computed from Corrfunc)\n");
        fprintf(stderr, "   -mbin <mbin>:  The number of mu bins (spaced linearly).\n");
#endif
	    fprintf(stderr, "   -mbin_cf <mbin_cf>:  The number of mu bins in the correlation function (spaced linearly).\n");
        fprintf(stderr, "   -output: (Pre-existing) directory to save output covariance matrices into\n");
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

#if (defined LEGENDRE || defined LEGENDRE_MIX)
        fprintf(stderr, "   -max_l <max_l>: Maximum legendre multipole (must be even)\n");
#endif
#ifdef LEGENDRE_MIX
        fprintf(stderr, "   -mu_bin_legendre_file <filename>: Mu bin Legendre factors file\n");
#elif defined LEGENDRE
        fprintf(stderr, "   -phi_file <filename>: Survey correction function coefficient file\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "   -phi_file12 <filename>: (Optional) Survey correction function coefficient file for fields 1 x 2\n");
        fprintf(stderr, "   -phi_file2 <filename>: (Optional) Survey correction function coefficent file for field 2\n");
        fprintf(stderr, "\n");
#elif defined POWER
        fprintf(stderr, "   -inv_phi_file <filename>: Inverse survey correction function multipole coefficient file\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "   -inv_phi_file12 <filename>: (Optional) Inverse survey correction function coefficient file for fields 1 x 2\n");
        fprintf(stderr, "   -inv_phi_file2 <filename>: (Optional) Inverse survey correction function coefficent file for field 2\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "   -R0 <R0>: Truncation radius for pair-wise separation window function in Mpc/h. Default: R0 = 100\n");
        fprintf(stderr, "   -power_norm: Power spectrum normalization = V*<(nw)^2> = Sum(nw^2)\n");
        fprintf(stderr, "   -power_norm2: Power spectrum normalization = V*<(nw)^2> = Sum(nw^2) for field 1 x 2\n");
        fprintf(stderr, "   -power_norm12: Power spectrum normalization = V*<(nw)^2> = Sum(nw^2) for field 2\n");
#endif
#if (!defined LEGENDRE && !defined POWER && !defined THREE_PCF)
        fprintf(stderr, "   -RRbin12 <filename>: (Optional) File containing the {1,2} jackknife RR bin counts (computed from Corrfunc)\n");
	    fprintf(stderr, "   -RRbin2 <filename>: (Optional) File containing the {2,2} jackknife RR bin counts (computed from Corrfunc)\n");
        fprintf(stderr, "\n");
#endif
#ifdef JACKKNIFE
        fprintf(stderr, "   -jackknife12 <filename>: (Optional) File containing the {1,2} jackknife weights (normally computed from Corrfunc)\n");
        fprintf(stderr, "   -jackknife2 <filename>: (Optional) File containing the {2,2} jackknife weights (normally computed from Corrfunc)\n");
#endif
        fprintf(stderr, "   -maxloops <max_loops>: Maximum number of integral loops\n");
        fprintf(stderr, "   -loopspersample <loops_per_sample>: Number of loops to collapse into each subsample. Default 1.\n");
        fprintf(stderr, "   -N2 <N2>: Number of secondary particles to choose per primary particle\n");
        fprintf(stderr, "   -N3 <N3>: Number of tertiary particles to choose per secondary particle\n");
        fprintf(stderr, "   -N4 <N4>: Number of quaternary particles to choose per tertiary particle\n");
        fprintf(stderr, "\n");
#ifdef THREE_PCF
        fprintf(stderr, "   -N5 <N5>: Number of fifth particles to choose per quaternary particle\n");
        fprintf(stderr, "   -N6 <N6>: Number of sixth particles to choose per fifth particle\n");
        fprintf(stderr, "\n");
#endif
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
        fprintf(stderr, "   -seed <seed>:  If given, allows to reproduce the results with the same settings, except the number of threads.\n");
        fprintf(stderr, "      Individual subsamples may differ because they are accumulated/written in order of loop completion which may depend on external factors at runtime, but the final integrals should be the same.\n");
        fprintf(stderr, "      Needs to be a non-negative integer fitting into 32 bits (max 4294967295). If not given (default), the seed will be created randomly.\n");
	    fprintf(stderr, "\n");
	    fprintf(stderr, "\n");

	    exit(1);
	}

	void create_directory(){
        // Initialize output directory:
        size_t out_file_len = strlen(out_file);
        if (out_file[out_file_len - 1] != '/') { // append the slash if absent
            out_file_len += 2; // add one for the slash, and another for the terminating zero which was not included
            char * tmp_out_file = (char *) malloc(sizeof(char) * out_file_len); // one more character for the slash and another for zero terminator
            snprintf(tmp_out_file, out_file_len, "%s/", out_file);
            out_file = tmp_out_file;
            // one might think the old out_file should be freed, but actually not, because out_file is either pointing to default_out_file (a constant member of this Parameters class) or an element of argv[], neither can and should be freed. There is still some additional memory usage, but should be small in realistic cases.
        }
	    // First create whole directory if it doesn't exist:
	    if (mkdir(out_file,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)==0){
            printf("\nCreating output directory\n");
        }
#ifdef THREE_PCF
        std::string cname = string_format("%s3PCFCovMatricesAll/", out_file);
#elif defined POWER
        std::string cname = string_format("%sPowerCovMatrices/", out_file);
#else
	    std::string cname = string_format("%sCovMatricesAll/", out_file);
#endif
        mkdir(cname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        // Check if this was successful:
        struct stat info;

        if( stat( cname.c_str(), &info ) != 0 ){
            printf( "\nCreation of directory %s failed\n", cname.c_str());
            exit(1);
        }
#ifdef JACKKNIFE
        std::string cjname = string_format("%sCovMatricesJack/",out_file);
        mkdir(cjname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        if(stat(cjname.c_str(), &info) != 0) {
            printf("\nCreation of directory %s failed\n", cjname.c_str());
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
                split_string = strtok(line, " \t");
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
                    split_string = strtok(NULL, " \t");
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
                split_string = strtok(line, " \t");
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
                    split_string = strtok(NULL, " \t");
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
