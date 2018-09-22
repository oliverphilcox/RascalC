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

	// The maximum radius of the largest bin.
	Float rmax = 200.0;

	// The minimum radius of the smallest bin.
	Float rmin = 50.0;

	// The maximum mu of the largest bin.
	Float mumax = 1.0;

	// The minimum mu of the smallest bin.
	Float mumin = 0.0;

	// The radius beyond which the correlation function is set to zero
	Float xicutoff = 400.0;

	Float nofznorm=626798;//6684485//681013//672940//674847 629310

	// The grid size, which should be tuned to match boxsize and rmax.
	// Don't forget to adjust this if changing boxsize!
	int nside = 51;

	// If set, we'll just throw random periodic points instead of reading the file
	int make_random = 0;

	// Will be number of particles in a random distribution,
	// but gets overwritten if reading from a file.
	int np = -1;

	// Whether to balance the weights or multiply them by -1
	int qbalance = 0, qinvert = 0;

	// The maximum number of points to read
	uint64 nmax = 1000000000000;

    // The number of radial bins
    // NB: This doesn't need to be equal to the number in the xi file
	int nbin = 10;

    // The number of mu bins
	int mbin = 6;

	// The number of threads to run on
	int nthread=4;

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
	const char default_fname[500] = "../random_particles/test_particles_small.txt";

	// The name of the correlation function file
	char *corname = NULL;
	const char default_corname[500] = "xi_functions/QPM_Mash.xi";
	//"../grid_multipoles_own/PatchySkyCorrSingle361.xi"//PatchySkyCorrMean.xi//QPMCorrMean.xi//QPMExtrapolated.xi//"QPM_D_ngc_rsd_fix3.xi"

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
        else if (!strcmp(argv[i],"-rescale")||!strcmp(argv[i],"-scale")) rescale = atof(argv[++i]);
		else if (!strcmp(argv[i],"-rmax")||!strcmp(argv[i],"-max")) rmax = atof(argv[++i]);
		else if (!strcmp(argv[i],"-rmin")) rmin = atof(argv[++i]);
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
		else if (!strcmp(argv[i],"-nbin")) nbin = atoi(argv[++i]);
		else if (!strcmp(argv[i],"-mbin")) mbin = atoi(argv[++i]);
		else if (!strcmp(argv[i],"-save")||!strcmp(argv[i],"-store")) savename = argv[++i];
		else if (!strcmp(argv[i],"-load")) loadname = argv[++i];
		else if (!strcmp(argv[i],"-balance")) qbalance = 1;
		else if (!strcmp(argv[i],"-invert")) qinvert = 1;
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

	    assert(box_min>0.0);
	    assert(rmax>0.0);
	    assert(nside>0);
	    if (rescale<=0.0) rescale = box_max;   // This would allow a unit cube to fill the periodic volume
	    if (fname==NULL) fname = (char *) default_fname;   // No name was given
	    if (corname==NULL) { corname = (char *) default_corname; }//fprintf(stderr,"No correlation file."); return 1;}// No name was given

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
		printf("Radial Binning = {%6.5f, %6.5f, %6.5f}\n",rmin,rmax,(rmax-rmin)/nbin);
		printf("Mu Bins = %d\n", mbin);
		printf("Mu Binning = {%6.5f, %6.5f, %6.5f}\n",mumin,mumax,(mumax-mumin)/mbin);
		printf("Density Normalization = %6.5e\n",nofznorm);

	}
private:
	void usage() {
	    fprintf(stderr, "\nUsage for grid_covariance:\n");
	    fprintf(stderr, "   -box <boxsize> : If creating particles randomly, this is the periodic size of the cubic computational domain.  Default 400. If reading from file, this is reset dynamically creating a cuboidal box.\n");
	    fprintf(stderr, "   -scale <rescale>: How much to dilate the input positions by.  Default 0.\n");
	    fprintf(stderr, "             Zero or negative value causes =boxsize, rescaling unit cube to full periodicity\n");
	    fprintf(stderr, "   -rmax <rmax>: The maximum radius of the largest pair bin.  Default 200.\n");
	    fprintf(stderr, "   -xicut <xicutoff>: The radius beyond which xi is set to zero.  Default 1000.\n");
	    fprintf(stderr, "   -nside <nside>: The grid size for accelerating the pair count.  Default 8.\n");
	    fprintf(stderr, "             Recommend having several grid cells per rmax.\n");
	    fprintf(stderr, "   -in <file>: The input file (space-separated x,y,z,w).  Default sample.dat.\n");
	    fprintf(stderr, "   -ran <np>: Ignore any file and use np random perioidic points instead.\n");
	    fprintf(stderr, "   -norm <nofznorm>: Number of galaxies in the survey. Used to normalize n(z).\n");
	    fprintf(stderr, "   -def: This allows one to accept the defaults without giving other entries.\n");
	    fprintf(stderr, "\n");
        fprintf(stderr, "   -cor <file>: File location of inut correlation function file.\n");
	    fprintf(stderr, "   -nbin:  The number of radial bins.\n");
	    fprintf(stderr, "   -mbin:  The number of mu bins.\n");
	    fprintf(stderr, "The radial bin spacing (currently linear) is hard-coded.\n");
	    fprintf(stderr, "\n");
	    fprintf(stderr, "For advanced use, there is an option store the grid of probabilities used for sampling.\n");
	    fprintf(stderr, "    -save <filename>: Triggers option to store probability grid. <filename> has to end on \".bin\"\n");
	    fprintf(stderr, "The file can then be reloaded on subsequent runs\n");
	    fprintf(stderr, "    -load <filename>: Triggers option to load the probability grid\n");
	    fprintf(stderr, "    -balance: Rescale the negative weights so that the total weight is zero.\n");
	    fprintf(stderr, "    -invert: Multiply all the weights by -1.\n");


	    exit(1);
	}

};

#endif
