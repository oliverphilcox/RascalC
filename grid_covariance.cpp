// grid_covariance.cpp -- Alexander Wiegand, started Oct 14, 2016 based on grid_multipoles.cpp of Daniel Eisenstein.

#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
//#include <random>
#include <complex>
#include <algorithm>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include "threevector.hh"
#include <gsl/gsl_rng.h>
#include "./ransampl/ransampl.h"
#include "STimer.cc"
#include "./cubature/cubature.h"
#include <limits>

#define OPENMP

// For multi-threading:
#ifdef OPENMP
//#define _GNU_SOURCE
#include <omp.h>
#include <sched.h>
#endif

#undef ALLOUT
//#undef PERIODIC
#define PAGE 4096     // To force some memory alignment.

typedef unsigned long long int uint64;

// Could swap between single and double precision here.
// Only double precision has been tested.
typedef double Float;
typedef double3 Float3;

class Parameters{

public:
	// Important variables to set!  Here are the defaults:
	// The periodicity of the position-space cube.
	Float boxsize = 400;

	// The particles will be read from the unit cube, but then scaled by boxsize.
	Float rescale = 0.0;   // If left zero or negative, set rescale=boxsize

	// The maximum radius of the largest bin.
	Float rmax = 200.0;

	// The minimum radius of the smallest bin.
	Float rmin = 0.0;

	// The maximum mu of the largest bin.
	Float mumax = 1.0;

	// The minimum mu of the smallest bin.
	Float mumin = 0.0;

	// The radius beyond which the correlation function is set to zero
	Float xicutoff = 1000.0;

	Float nofznorm=626798;//6684485//681013//672940//674847 629310

	// The grid size, which should be tuned to match boxsize and rmax.
	// Don't forget to adjust this if changing boxsize!
	int nside = 8;

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
	int nbin = 25;

    // The number of mu bins
	int mbin = 1;

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
	const char default_fname[500] = "sample.dat";

	// The name of the correlation function file
	char *corname = NULL;
	const char default_corname[500] = "../grid_multipoles_own/PatchySkyCorrSingle361.xi";
	//"../grid_multipoles_own/PatchySkyCorrSingle361.xi"//PatchySkyCorrMean.xi//QPMCorrMean.xi//QPMExtrapolated.xi//"QPM_D_ngc_rsd_fix3.xi"

	// Constructor
	Parameters(int argc, char *argv[]){

	    if (argc==1) usage();
	    int i=1;
	    while (i<argc) {
		     if (!strcmp(argv[i],"-boxsize")||!strcmp(argv[i],"-box")) boxsize = atof(argv[++i]);
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
	    assert(i==argc);  // For example, we might have omitted the last argument, causing disaster.

	    assert(nside%2!=0); // The probability integrator needs an odd grid size

	    assert(mumin>=0); // We take the absolte value of mu
	    assert(mumax<=1); // mu > 1 makes no sense

	    assert(boxsize>0.0);
	    assert(rmax>0.0);
	    assert(nside>0);
	    if (rescale<=0.0) rescale = boxsize;   // This would allow a unit cube to fill the periodic volume
	    if (fname==NULL) fname = (char *) default_fname;   // No name was given
	    if (corname==NULL) { corname = (char *) default_corname; }//fprintf(stderr,"No correlation file."); return 1;}// No name was given

#ifdef OPENMP
		omp_set_num_threads(nthread);
#else
		nthread=1;
#endif


		// Output for posterity
		printf("Box Size = %6.5e\n", boxsize);
		printf("Grid = %d\n", nside);
		printf("Maximum Radius = %6.5e\n", rmax);
		Float gridsize = rmax/(boxsize/nside);
		printf("Radius in Grid Units = %6.5e\n", gridsize);
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
	    fprintf(stderr, "   -box <boxsize> : The periodic size of the computational domain.  Default 400.\n");
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




// We need a vector floor3 function
Float3 floor3(float3 p) {
    return Float3(floor(p.x), floor(p.y), floor(p.z));
}

// =================== Particles ====================
// This is the info about each particle that we load in and store in the Grid.

class Particle {
  public:
    Float3 pos;
    Float w;  // The weight for each particle
};


// ====================  The Cell and Grid classes ==================

/* The Grid class holds a new copy of the particles.
These are sorted into cells, and all positions are referenced to the
cell center.  That way, we can handle periodic wrapping transparently,
simply by the cell indexing.

For simplicity, we opt to flatten the index of the cells into a 1-d number.
For example, this makes multi-threading over cells simpler.
*/

class Cell {
  public:
    int start;	// The starting index of the particle list
    int np;
};

class Grid {
  public:
    Float boxsize;   // Size of the periodic volume
    int nside, ncells;       // Grid size (per linear and per volume)
    Cell *c;		// The list of cells
    Float cellsize;   // Size of one cell
    Particle *p;	// Pointer to the list of particles
    int np;		// Number of particles
    int np_pos;		// Number of particles
    int *pid;		// The original ordering
    int *filled; //List of filled cells
    int nf;      //Number of filled cells
    int maxnp;   //Max number of particles in a single cell
    Float sumw_pos, sumw_neg; // Summing the weights

    int test_cell(integer3 cell){
    	// returns -1 if cell is outside the grid or wraps around for the periodic grid
#ifndef PERIODIC
    	if(nside<cell.x||cell.x<0||nside<cell.y||cell.y<0||nside<cell.z||cell.z<0)
    		return -1;
    	else
#endif
    		return wrap_cell(cell);
    }

    int wrap_cell(integer3 cell) {
        // Return the 1-d cell number, after wrapping
	// We apply a very large bias, so that we're
	// guaranteed to wrap any reasonable input.
	int cx = (cell.x+ncells)%nside;
	int cy = (cell.y+ncells)%nside;
	int cz = (cell.z+ncells)%nside;
	// return (cx*nside+cy)*nside+cz;
	int answer = (cx*nside+cy)*nside+cz;
	assert(answer<ncells&&answer>=0);
	/* printf("Cell: %d %d %d -> %d %d %d -> %d\n",
	    cell.x, cell.y, cell.z, cx, cy, cz, answer); */
	return answer;
    }

    integer3 cell_id_from_1d(int n) {
	// Undo 1d back to 3-d indexing
        assert(n>=0&&n<ncells);
	// printf("Cell ID: %d ", n);
	integer3 cid;
	cid.z = n%nside;
	n = n/nside;
	cid.y = n%nside;
	cid.x = n/nside;
	// printf("-> %d %d %d\n", cid.x, cid.y, cid.z);
	return cid;
    }

    int pos_to_cell(Float3 pos) {
        // Return the 1-d cell number for this position, properly wrapped
	// We assume the first cell is centered at cellsize/2.0
	// return wrap_cell( floor3(pos/cellsize+Float3(0.5,0.5,0.5)));
	return wrap_cell( floor3(pos/cellsize));
    }

    Float3 cell_centered_pos(Float3 pos) {
        // Subtract off the cell center from the given position.
	// This is safe for positions not in the primary box.
	return pos-cellsize*(floor3(pos/cellsize)+Float3(0.5,0.5,0.5));
    }

    Float3 cell_sep(integer3 sep) {
	// Return the position difference corresponding to a cell separation
        return cellsize*sep;
    }

    ~Grid() {
	// The destructor
        free(p);
	free(pid);
	free(c);
	return;
    }

    Grid(Particle *input, int _np, Float _boxsize, int _nside, Float3 shift) {
	// The constructor: the input set of particles is copied into a
	// new list, which is ordered by cell.
	// After this, Grid is self-sufficient; one could discard *input
    boxsize = _boxsize;
	nside = _nside;
	assert(nside<1025);   // Can't guarantee won't spill int32 if bigger
	ncells = nside*nside*nside;
	np = _np;
	np_pos = 0;
	assert(boxsize>0&&nside>0&&np>=0);
	cellsize = boxsize/nside;

	p = (Particle *)malloc(sizeof(Particle)*np);
	pid = (int *)malloc(sizeof(int)*np);
	c = (Cell *)malloc(sizeof(Cell)*ncells);

	// Now we want to copy the particles, but do so into grid order.

	// First, figure out the cell for each particle
	// Shift them to the primary volume first
	int *cell = (int *)malloc(sizeof(int)*np);
	for (int j=0; j<np; j++) cell[j] = pos_to_cell(input[j].pos - shift);

	// Histogram the number of particles in each cell
	int *incell = (int *)malloc(sizeof(int)*ncells);
	for (int j=0; j<ncells; j++) incell[j] = 0.0;
	for (int j=0; j<np; j++) incell[cell[j]]++;

	// Create list of filled cells
	nf=0;
	for (int j=0; j<ncells; j++) if(incell[j]>0) nf++;
	filled = (int *)malloc(sizeof(int)*nf);
	for (int j=0,k=0; j<ncells; j++) if(incell[j]>0) filled[k++]=j;

	// Count the number of positively weighted particles
	sumw_pos = sumw_neg = 0.0;
	for (int j=0; j<np; j++)
	    if (input[j].w>=0) {
	    	np_pos++;
		sumw_pos += input[j].w;
	    } else {
		sumw_neg += input[j].w;
	    }

	// Cumulate the histogram, so we know where to start each cell
	for (int j=0, tot=0; j<ncells; tot+=incell[j], j++) {
	    c[j].start = tot;
	    c[j].np = 0;  // We'll count these as we add the particles
	}

	// Copy the particles into the cell-ordered list
	for (int j=0; j<np; j++) {
	    Cell *thiscell = c+cell[j];
	    int index = thiscell->start+thiscell->np;
	    p[index] = input[j];
#ifdef PERIODIC
	    p[index].pos = cell_centered_pos(input[j].pos);
	    	// Switch to cell-centered positions
#endif
	    pid[index] = j;	 // Storing the original index

	    // Diagnostics:
	    /*
	    integer3 cid = cell_id_from_1d(cell[j]);
	    printf("P->C: %d %7.4f %7.4f %7.4f -> %d (%d %d %d) %7.4f %7.4f %7.4f %d\n",
		j, input[j].pos.x, input[j].pos.y, input[j].pos.z,
		cell[j], cid.x, cid.y, cid.z,
		p[index].pos.x, p[index].pos.y, p[index].pos.z, index);
	    */

	    thiscell->np += 1;
	}

	// Checking that all is well.
	int tot = 0;
	maxnp=0;
	for (int j=0; j<ncells; j++) {
	    assert(c[j].start == tot);
	    assert(c[j].np == incell[j]);
	    if(c[j].np>maxnp) maxnp=c[j].np;
	    tot += c[j].np;
	}
	free(incell);
	assert(tot == np);

	free(cell);
	return;
    }

};   // End Grid class


class CorrelationFunction{
	/* Reads and stores a 2d correlation function. Values in between the grid positions of the input
	 * are interpolated using gsl_interp2d
	 */
	private:
		int xsize, ysize;
		double *x,*y,*z;
		double rmin,rmax,mumin,mumax;
		bool mudim;
		gsl_interp_accel *xa, *ya, *x1a;
		gsl_interp2d* interp_2d;
		gsl_spline* corfu1d;

	public:
		double xi(double r, double mu){
			// 2D correlation function in radius and angular bins
			// xi values beyond the maximal radius in the correlation function file read in are extrapolated
			// as a simple r^-4 power law for each mu bin independently (may lead to different signs for
			// neighbouring mu bins if the correlation function is noisy and fluctuates around 0
			if(mudim){
				if(mu<mumin||mu>mumax||r<rmin||r>rmax){
					double tmu = fmin(fmax(mu,mumin),mumax);
					double tr = fmin(fmax(r,rmin),rmax);
					if(r>rmax){
						return gsl_interp2d_eval_extrap(interp_2d,y,x,z,tmu,tr, xa, ya)/pow(tr,2)*pow(r/rmax,-4);
					}
					else
						return gsl_interp2d_eval_extrap(interp_2d,y,x,z,tmu,tr, xa, ya)/pow(tr,2);
				}
				else{
					return gsl_interp2d_eval_extrap(interp_2d,y,x,z,mu,r, xa, ya)/pow(r,2);
				}
			}
			else
				return xi(r);
		}
		double xi2(double r12, double r2) {
			// Simple test correlation function to compare to sobol
			double r0 = 8.0;
		    return r0*r0/(r12*r12+1.0);
		    	// Simple case, r_0 = 8 Mpc, softening at 1 Mpc
		}

		double xi(double r){
			// Radial correlation function
			// xi values beyond the maximal radius in the correlation function file read in are extrapolated
			// as a simple r^-4 power law
			if(r>rmax){
				return gsl_spline_eval(corfu1d,rmax, x1a)/pow(rmax,2)*pow(r/rmax,-4);
			}
			else{
				if(r<rmin){
					return 0.;
				}
				else
					return gsl_spline_eval(corfu1d,r, x1a)/pow(r,2);
			}
		}


	private:

		void readLine(char *line,double **z, int *ct){
			char * pch;
			pch = strtok (line," \t\n");
			sscanf(pch, "%lf",  &((*z)[(*ct)++]) );
//				printf("\n%f ",(*z)[ct-1]);
			while (pch != NULL)
			{
			  pch = strtok (NULL, " \t\n");
			  if(pch != NULL){
				  sscanf(pch, "%lf",  &((*z)[(*ct)++]) );
//				  printf("%f ", (*z)[ct-1]);
			  }
			}
		}

		void readData(const char *filename,double **x,double **y, double **z,int *np,int *mp){

			char line[10000];
			char * pch;
			int n=0,m=0;
			double x0;
			FILE *fp;
			fp = fopen(filename, "r");
			if (fp==NULL) {
				fprintf(stderr,"File %s not found\n", filename); abort();
			}

			//Count lines and columns
			while (fgets(line,10000,fp)!=NULL) {
			    if (line[0]=='#') continue;
				if (line[0]=='\n') continue;
				if(n==0){
					pch = strtok (line," \t\n");
					sscanf(pch, "%lf",  &x0 );
				}
				if(n==2){
					m=0;
					pch = strtok (line," \t\n");
					while (pch != NULL)
					{
					  m++;
					  pch = strtok (NULL, " \t\n");
					}
				}

				n++;
			}
			rewind(fp);
			n-=2; //Remove the first two lines that should contain r and mu bin centers

			// If there is no r=0 entry, add one
			if(x0!=0.){
				n++;
			}

			*np = n;
			*mp = m;

			*x = (double *)malloc(sizeof(double)*n);
			*y = (double *)malloc(sizeof(double)*m);
			*z = (double *)malloc(sizeof(double)*n*m);

			printf("# Found %d radial and %d mu bins in %s\n", n,m, filename);

			int nx=0, ny=0, nz=0, lnnr=0;

			// Add a first entry at r=0 assuming that we shall multiply by r^2 anyway
			if(x0!=0.){
				(*x)[0]=0.;
				nx++;
				for(int i=0;i<m;i++){
					(*z)[i]=0.;
					nz++;
				}
			}

			//Read content of lines and columns
			while (fgets(line,10000,fp)!=NULL) {
				if (line[0]=='#') continue;
				if (line[0]=='\n') continue;
				if(lnnr==0)	readLine(line, x, &nx);
				if(lnnr==1)	readLine(line, y, &ny);
				if(lnnr>1)	readLine(line, z, &nz);
				lnnr++;
			}

			if(n!=nx){
				fprintf(stderr,"Found %d lines but %d r-values in the first line. Aborting.\n", n, nx); abort();
			}
			if(m!=ny){
				fprintf(stderr,"Found %d columns but %d mu-values in the second line. Aborting.\n", m, ny); abort();
			}

//			for(int i=0;i<n*m;i++){
//				printf("%f ",(*z)[i]);
//			}

		}

	private:
		void copy(int& _xsize, int& _ysize, double ** _x, double ** _y, double ** _z, double& _rmin,
				double& _rmax, double& _mumin, double& _mumax, bool& _mudim){
			// Copy correlation function object

			_xsize=xsize;
			_ysize=ysize;
			_rmin=rmin;
			_rmax=rmax;
			_mumin=mumin;
			_mumax=mumax;
			_mudim=mudim;

			*_x = (double *)malloc(sizeof(double)*xsize);
			*_y = (double *)malloc(sizeof(double)*ysize);
			*_z = (double *)malloc(sizeof(double)*xsize*ysize);

			int ct=0;
			for(int i=0;i<xsize;i++){
				(*_x)[i]=x[i];
				for(int j=0;j<ysize;j++){
					(*_z)[ct]=z[ct];
					ct++;
				}
			}

			for(int j=0;j<ysize;j++){
				(*_y)[j]=y[j];
			}
		}

		void interpolate(){
			double *y1;

			y1 = (double *)malloc(sizeof(double)*xsize);

			// Sum to obtain radial correlation function
			int col=0;
			int ct=0;
			for(int i=0;i<xsize;i++){
				y1[i]=0;
				col=0;
				for(int j=0;j<ysize;j++){
	//				if((j==0&&y[0]==0.)||(j==ysize-1&&y[ysize-1]==1.))
	//					ct++;
	//				else{
						y1[i]+=z[ct++];
						col++;
	//				}

				}
				y1[i]/=col;
	//			fprintf(stderr,"%f ",y1[i]);
			}

			if(mudim){
				interp_2d=gsl_interp2d_alloc(gsl_interp2d_bicubic, ysize, xsize);
				gsl_interp2d_init(interp_2d, y, x, z, ysize, xsize);
				xa = gsl_interp_accel_alloc();
				ya = gsl_interp_accel_alloc();
			}

			corfu1d=gsl_spline_alloc(gsl_interp_cspline, xsize);
			gsl_spline_init(corfu1d, x, y1, xsize);
			x1a = gsl_interp_accel_alloc();
		}
	public:

		CorrelationFunction(CorrelationFunction* corr){
			// Copy constructor

			corr->copy(xsize,ysize,&x,&y,&z,rmin,rmax,mumin,mumax,mudim);
			interpolate();
		}

	CorrelationFunction(const char *filename, int mbin, double dmu){
		// Construct from input file

		readData(filename,&x,&y,&z,&xsize,&ysize);

		if(ysize<=1&&mbin>1){
			fprintf(stderr,"Requested %d mu-bins but correlation function file only has %d column.",mbin,ysize);
			abort();
		}

		mudim=!(mbin==1&&dmu==1.);

//		mudim=true;

		rmin=x[0];
		rmax=x[xsize-1];
		// Mu range hardcoded to be between 0 and 1
		mumin=0.0;
		mumax=1.;

		// Multiply correlation function by r^2 for smoother interpolation
		// If one wants to change this also the readData function needs to be changed
		int ct=0;
		for(int i=0;i<xsize;i++){
			for(int j=0;j<ysize;j++){
				z[ct++]*=pow(x[i],2); //Multiply by r^2
			}
		}

		interpolate();

//		Old tests of the correlation function and its interpolation
//		printf("\n");
//		printf("\n");
//		for(int i=0;i<xsize;i++){
//			for(int j=0;j<ysize;j++){
//			printf("%e ", xi(x[i],y[j]-0.001*0));
//			}
//			printf("\n");
//		}

//		printf("\n");
//		printf("\n");
//		for(int i=0;i<2000;i++){
//			for(int j=0;j<1000;j++){
//			printf("%e ", xi(0.25+i/2.,j/1000.));
//			}
//			printf("\n");
//		}

//		for(int i=0;i<xsize-1;i++){
//			printf("%f=%f\n",x[i+1],xs[i]);
//		}
//
//		for(int i=0;i<ysize-2;i++){
//			printf("%f=%f\n",y[i+1],ys[i]);
//		}
//		for(int i=0;i<ysize*xsize;i++){
//			printf("%f ",z[i]);
//		}

//    gsl_rng_env_setup();
//    gsl_rng* rng = gsl_rng_alloc( gsl_rng_default );
//
//    double* binxi = new double[40];
//    int* bincount = new int[40];
//    for(int i=0;i<40;i++){
//    	binxi[i]=0.;
//    	bincount[i]=0;
//    }
//    double ha,hi;
//    for(int i=0;i<200000000;i++){
//    	ha=160*gsl_rng_uniform(rng);
//    	hi=gsl_rng_uniform(rng);
//    	binxi[(int)floor(ha/4)]+=xi(ha,hi);
//    	bincount[(int)floor(ha/4)]++;
//    }
//
//    for(int i=0;i<40;i++){
//		printf("%f %d %f\n",binxi[i],bincount[i],binxi[i]/bincount[i]);
//	}
//		exit(0);
	}

	~CorrelationFunction() {
		// Destructor

		if(mudim){
			gsl_interp2d_free(interp_2d);
			gsl_interp_accel_free(xa);
			gsl_interp_accel_free(ya);
		}
		gsl_spline_free(corfu1d);
		gsl_interp_accel_free(x1a);
	}

};



class RandomDraws{

	// Class that handles the selection of neighbouring boxes

public:
	static CorrelationFunction *corr;
	int nside;     // Number of cells in each direction of large draw
	int nsidecube; // Number of cells in each direction of maxsep cube
	double boxside;
	double *x;
	double *xcube;

	private:
		// Sampling of long distance
		ransampl_ws* ws;

		// Sampling of short distance
		ransampl_ws* cube;

//		std::default_random_engine genran;
//		std::uniform_real_distribution<double> celln;
//	    gsl_rng_env_setup();
//  	gsl_rng* rng = gsl_rng_alloc( gsl_rng_default );
//  	gsl_rng_uniform(rng)


	public:

		RandomDraws(CorrelationFunction *fun,Parameters *par,const double *xin, long np){ //,Float boxsize_,int nside_

		long n=np;
		corr=fun;

		//nside=par->nside;
		//boxside=par->boxsize/nside;
		boxside=par->boxsize/par->nside;
//TODO: Fix nside calculation
		nside=2*ceil(par->xicutoff/boxside)+1;


//		celln=std::uniform_real_distribution<double>(0.,1.);

		// If precalculated grid has been saved, load it
		if (par->loadname!=NULL&&(xin==NULL||np==0)){
			int len = strlen(par->loadname);
			if(len>4&&!strcmp(&par->loadname[len-4],".bin")){
				// Read binary format that was written by writeDataBin
				readDataBin(&x,&n,par->loadname);
			}
			else{
				// Read ascii format possibly created externally
				readData(&x,&n,par->loadname);
			}
			if(n==0)
				integData(&x,&n,f_corr,nside,boxside,0);
		}
		else
			// If no correlation function is given to copy, do integration
			if (xin==NULL||np==0)
				integData(&x,&n,f_corr,nside,boxside,0); // could change to 1/r^2 here by using f instead of f_corr and 1 instead of 0
			else
				copyData(&x,&n,xin);

		// Save grid of precalculated probabilities
		if (par->savename!=NULL){
			int len = strlen(par->savename);
			if(len>4&&!strcmp(&par->savename[len-4],".bin"))
				writeDataBin(&x,&n,par->savename);
			else
				fprintf(stderr,"Save file %s does not end on \".bin\". No output written.\n",par->savename);
		}

		// Set up actual sampler
		ws = ransampl_alloc( n );
		ransampl_set( ws, x );

		// Normalize grid probabilities to one
		double sum=0.;
		for(int i=0;i<n;i++){
			sum+=x[i];
		}
		for(int i=0;i<n;i++){
			x[i]/=sum;
		}


//		Initialize second sampler

		int maxsep = ceil(par->rmax/boxside);
		nsidecube = 2 * maxsep + 1;
		long nn=0;

		integData(&xcube,&nn,f,nsidecube,boxside,1./2.);

		// Set up actual sampler
		cube = ransampl_alloc( nn );
		ransampl_set( cube, xcube );

		// Normalize grid probabilities to one
		sum=0.;
		for(int i=0;i<nn;i++){
			sum+=xcube[i];
		}
		for(int i=0;i<nn;i++){
			xcube[i]/=sum;
		}


		}


~RandomDraws() {
		ransampl_free( ws );
		ransampl_free( cube );
		free(x);
		free(xcube);
	}

		integer3 random_xidraw(gsl_rng* rng, double* p){
			// Draws the index of a box at some distance which is currently weighted by the correlation function
			int n=ransampl_draw( ws, gsl_rng_uniform(rng), gsl_rng_uniform(rng) );
			*p=x[n];
			return cubifyindex(nside,n);
		}

		integer3 random_cubedraw(gsl_rng* rng, double* p){
			// Can be used to draw only a subset of the boxes within maxsep
			int n=ransampl_draw( cube, gsl_rng_uniform(rng), gsl_rng_uniform(rng) );
			*p=xcube[n];
			return cubifyindex(nsidecube,n);
		}

		// Undo 1d back to 3-d indexing
		integer3 cubifyindex(int nside,int n){

			assert(n>=0&&n<pow(nside,3));

			// printf("Cell ID: %d ", n);
			integer3 cid;
			cid.z = n%nside-((nside-1)/2);
			n = n/nside;
			cid.y = n%nside-((nside-1)/2);
			cid.x = n/nside-((nside-1)/2);
			// printf("-> %d %d %d\n", cid.x, cid.y, cid.z);

			return cid;
		}

	private:
		void readData(double **x,long *np,const char *filename){

			char line[10000];
			int  n=0;
			FILE *fp;
//			int stat;
//			double tmp[6];
			fp = fopen(filename, "r");
			if (fp==NULL) {
				fprintf(stderr,"File %s not found\n", filename); abort();
			}

			//Count lines
			while (fgets(line,10000,fp)!=NULL) {
			    if (line[0]=='#') continue;
				if (line[0]=='\n') continue;
				n++;
			}
			rewind(fp);
			*np = n;

			if(n!=pow(nside,3)){
				fprintf(stderr,"File %s does not contain the correct probability grid.\n", filename); abort();
			}


			//Read content of lines and columns
			*x = (double *)malloc(sizeof(double)*n);
			printf("# Found %d lines in %s\n", n, filename);


			int ct=0;
			while (fgets(line,10000,fp)!=NULL) {
				if (line[0]=='#') continue;
				if (line[0]=='\n') continue;
				sscanf(line, "%lf",  &((*x)[ct++]) );
			}

			assert(ct==n);

			fclose(fp);

		}

		void integData(double **x, long *np, integrand fun, int nside, double boxside, double power){
			// Implements the integration of a function fun over all possible distance between
			// points in two boxes in a grid (one of them being the central box)
			// of nside*nside*nside and with the boxes having sidelength boxside

			(*np)=(int)pow(nside,3);

			*x = (double *)malloc(sizeof(double)*(*np));

			printf("\nSize: %ld\n",(*np));
			fflush(NULL);

//			TODO: assure it compiles without OPENMP defed
			// parallel environment before for loop to have local integrator variables
			#pragma omp	parallel
			{
			int len=(nside-1)/2; // This works because nside has been required to be odd

			// Parameters for hcubature basically integration limits for the three x dimensions
			// and parameters for the integrand in n
			double xmin[3] = {0,0,0}, xmax[3] = {1,1,1}, n[5] = {0,0,0,boxside/2., power}, val, err;
#ifdef OPENMP
			#pragma omp for schedule(dynamic,32)
#endif
			// Iterates over the grid
				for(int i=0;i<=len;i++){
					for(int k=0;k<=len;k++){
						for(int l=0;l<=len;l++){
							n[0]=i;
							n[1]=k;
							n[2]=l;
							hcubature(1, fun, &n[0], 3, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &val, &err);

							// Insert the calculated value to its grid position
							for(int is=-1;is<=1;is+=2){
								for(int ks=-1;ks<=1;ks+=2){
									for(int ls=-1;ls<=1;ls+=2){
										int ic=is*i+len;
										int kc=ks*k+len;
										int lc=ls*l+len;

										(*x)[nside*nside*ic+nside*kc+lc]=val;
		//								if(i==len&&k==len)
		//									printf("%d %f, ",side*side*ic+side*kc+lc,y[side*side*ic+side*kc+lc]/val-1);
									}
								}
							}
		//					if(i==len&&k==len)
		//						printf("\n");



		//	    			hcubature_v(1, f_v, &n[0], 3, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &val, &err);
		//		    			printf("Computed integral = %0.10g +/- %g\n", val, err);

						}
					}
				}
			}
		}
		void copyData(double **x, long *n,const double *xin){
			// Copy probability grid
			(*x) = (double *)malloc(sizeof(double)*(*n));
			for(int i=0;i<(*n);i++){
				(*x)[i]=xin[i];
			}
		}

		void readDataBin(double **x, long *n, const char *filename){

			FILE *fp;
			int _nside;
			double _boxside;
			int stat=0;
			fp = fopen(filename, "rb");
			if (fp==NULL) {
				fprintf(stderr,"# File %s not found\n", filename); return;
			}
			stat+=fread(&_nside, sizeof(int), 1, fp);
			stat+=fread(&_boxside, sizeof(double), 1, fp);
			if(!(_nside==nside&&_boxside==boxside)){
				fprintf(stderr,"# Size has changed. Recalculating probability grid.\n");
				fflush(NULL);
				*n=0;
				fclose(fp);
				return;
			}
			stat+=fread(n, sizeof(long), 1, fp);
			*x = (double *)malloc(sizeof(double)*(*n));
			if(*x==NULL){
				fprintf(stderr,"Allocation error.\n");
				fflush(NULL);
				abort();
			}

			if(stat!=3||*n!=(long)fread(*x, sizeof(double), *n, fp)){
	        	fprintf(stderr,"# Error reading file %s. Recalculating probability grid.\n", filename);
	        	fflush(NULL);
	        	*n=0;
	        	free(*x);
	        	fclose(fp);
	        	return;
	        }

	        assert(*n==pow(nside,3));

	        fclose(fp);

		}



		void writeDataBin(double **x, long *n, const char *filename){

			FILE *fp;
			fp = fopen(filename, "wb");
			if (fp==NULL) {
				fprintf(stderr,"# File %s not found\n", filename);
				return;
			}
			fwrite(&nside, 1, sizeof(int), fp);
			fwrite(&boxside, 1, sizeof(double), fp);
			fwrite(n, 1, sizeof(long), fp);
	        if(*n*sizeof(double)!=fwrite(*x, 1, *n*sizeof(double), fp)){
	        	fprintf(stderr,"# Error writing to file %s.\n", filename);
	        	return;
	        }

	        fclose(fp);

		}

	public:
	static int f_corr(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
	    // integrand using the correlation function to weight the cells as done in Rascal
		double* n = (double *) fdata;
	    double sum = 0;
	    int i,l,k;

	    for(i=-1;i<=1;i+=2){
	    	for(k=-1;k<=1;k+=2){
	    		for(l=-1;l<=1;l+=2){
	    		    sum+=corr->xi(2*n[3]*sqrt((pow(x[0] + i*n[0],2) + pow(x[1] + k*n[1],2) + pow(x[2] + l*n[2],2))));
	    		}
	    	}
	    }
	    fval[0] = fabs(pow(2*n[3],6)*(1 - x[0])*(1 - x[1])*(1 - x[2])*sum);
		return 0;
	}

	static int f(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
		// integrand implementing 1/r^2n where n can be an arbitrary double
	    double* n = (double *) fdata;
	    double sum = 0;
	    int i,l,k;

	    for(i=-1;i<=1;i+=2){
	    	for(k=-1;k<=1;k+=2){
	    		for(l=-1;l<=1;l+=2){
	    		    sum+=1./(1./pow(2*n[3],2.*n[4]) + pow(pow(x[0] + i*n[0],2) + pow(x[1] + k*n[1],2) + pow(x[2] + l*n[2],2),n[4]));
	    		}
	    	}
	    }
	    fval[0] = pow(2*n[3],6.-2.*n[4])*(1 - x[0])*(1 - x[1])*(1 - x[2])*sum;
	    return 0;
	}

	static int f_cube(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
	    // cubic integrand
		double* n = (double *) fdata;
	    double sum = 0;
	    int i,l,k;

	    for(i=-1;i<=1;i+=2){
	    	for(k=-1;k<=1;k+=2){
	    		for(l=-1;l<=1;l+=2){
	    		    sum+=1./(1./pow(2*n[3],6) + pow(pow(x[0] + i*n[0],2) + pow(x[1] + k*n[1],2) + pow(x[2] + l*n[2],2),3));
	    		}
	    	}
	    }
	    fval[0] = pow(2*n[3],0)*(1 - x[0])*(1 - x[1])*(1 - x[2])*sum;
	    return 0;
	}


	// Test of different parallelization of the integration used for setting up the probability grid
	// However, no significant performance gain observed
	// Also, the setup of the grid already takes only a very small fraction of the full calculation
	// so further optimization not urgent

	#pragma omp declare simd uniform(n,x,ndim) linear(j:1)
	double eval(const double* n,const double * x, unsigned ndim,unsigned j){
	    int i,l,k;
		double sum = 0;
			for(i=-1;i<=1;i+=2){
				for(k=-1;k<=1;k+=2){
					for(l=-1;l<=1;l+=2){
						sum+=1./(1./pow(2*n[3],4) + pow(pow(x[j*ndim+0] + i*n[0],2) + pow(x[j*ndim+1] + k*n[1],2) + pow(x[j*ndim+2] + l*n[2],2),2));
					}
				}
			}
			return pow(2*n[3],2)*(1 - x[j*ndim+0])*(1 - x[j*ndim+1])*(1 - x[j*ndim+2])*sum;
	}


	int f_v(unsigned ndim, unsigned long npts, const double *x, void *fdata, unsigned fdim, double *fval) {

	    double* n = (double *) fdata; // we can pass σ via fdata argument
	    unsigned j;

	#pragma omp simd
	    for (j = 0; j < npts; ++j) {
		    fval[j] = eval(n,x,ndim,j);
	    }
	    return 0; // success
	}


};

// Very ugly way to get the correlation function into the integrator, but hey, it works

CorrelationFunction * RandomDraws::corr;

// ========================== Accumulate Integral ================

class Integrals {

	public:
//	RandomDraws *rd;
	CorrelationFunction *cf;

#ifdef ALLOUT
	Float rij_magt;
	Float rij_mut;
	Float rjk_magt;
	Float rjk_mut;
#endif

	private:
	int nbin,mbin;
	Float rmin,rmax,mumin,mumax,dr,dmu; // Ranges in r and mu
	Float *Ra, *cx, *c2, *c3, *c4; // Arrays to accumulate the integrals
	Float *Raerr, *cxerr, *c2err, *c3err, *c4err; // Arrays to accumulate the errors; however, currently used for other purposes!
	Float xiij; //Only temporary for xi^2

	uint64 *binct,*binct3,*binct4; // Arrays to accumulate the counts of the bins
	bool box,rad; // Flags to decide whether we have a periodic box and if we have a radial correlation function only

#ifdef ALLOUT
	FILE**rij_magf;
	FILE**rij_muf;
	FILE**rik_magf;
	FILE**rjk_magf;
	FILE**rik_muf;
	FILE**rjk_muf;
#endif

	double empty[8];   // Leftover from grid_multipoles. Probably no longer needed

	  public:
		Integrals(Parameters *par, CorrelationFunction *_cf){
			cf=new CorrelationFunction(_cf);
			init(par);
		}

//		void init_cf_rd(double*x,long n,Parameters *par) {
//			cf=new CorrelationFunction(par->corname,par->mbin,par->mumax-par->mumin);
//			if(n==0)
//				rd=new RandomDraws(cf,par,x,n);
//			else
//				rd=new RandomDraws(cf,par,x,n);
//		}

	    void init(Parameters *par) {

	    	nbin=par->nbin;
	    	mbin=par->mbin;

			int ec=0;
			// Initialize the binning
			ec+=posix_memalign((void **) &Ra, PAGE, sizeof(double)*nbin*mbin);
			ec+=posix_memalign((void **) &cx, PAGE, sizeof(double)*nbin*mbin);
			ec+=posix_memalign((void **) &c2, PAGE, sizeof(double)*nbin*mbin);
			ec+=posix_memalign((void **) &c3, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);
			ec+=posix_memalign((void **) &c4, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);

			ec+=posix_memalign((void **) &Raerr, PAGE, sizeof(double)*nbin*mbin);
			ec+=posix_memalign((void **) &cxerr, PAGE, sizeof(double)*nbin*mbin);
			ec+=posix_memalign((void **) &c2err, PAGE, sizeof(double)*nbin*mbin);
			ec+=posix_memalign((void **) &c3err, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);
			ec+=posix_memalign((void **) &c4err, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);

			ec+=posix_memalign((void **) &binct, PAGE, sizeof(uint64)*nbin*mbin);
			ec+=posix_memalign((void **) &binct3, PAGE, sizeof(uint64)*nbin*mbin*nbin*mbin);
			ec+=posix_memalign((void **) &binct4, PAGE, sizeof(uint64)*nbin*mbin*nbin*mbin);

			assert(ec==0);

			reset();

			box=par->perbox;

			rmax=par->rmax;
			rmin=par->rmin;

			mumax=par->mumax;
			mumin=par->mumin;

			dr=(rmax-rmin)/nbin;
			dmu=(mumax-mumin)/mbin;

	    	rad=mbin==1&&dmu==1.;

			empty[0] = 0.0;   // To avoid a warning

#ifdef ALLOUT
			rij_magf=(FILE**)malloc(sizeof(FILE*)*nbin*mbin);

	//		rij_muf=(FILE**)malloc(sizeof(FILE*)*nbin*mbin);
	//		rik_magf=(FILE**)malloc(sizeof(FILE*)*nbin*mbin);
	//		rjk_magf=(FILE**)malloc(sizeof(FILE*)*nbin*mbin);
	//		rik_muf=(FILE**)malloc(sizeof(FILE*)*nbin*mbin);
	//		rjk_muf=(FILE**)malloc(sizeof(FILE*)*nbin*mbin);

			char buf[2000];

			for (int j=0; j<nbin*mbin; j++) {

				snprintf(buf, sizeof buf, "/mnt/store1/awiegand/CovarianceCalculation/Data/own_xiright_xisampling/own_10Mpc_all2pt_rij_mag-%d-%d.txt", j/mbin, j%mbin);
				rij_magf[j]=fopen(buf,"w");
	//			snprintf(buf, sizeof buf, "/mnt/store1/awiegand/CovarianceCalculation/Data/own_all2pt_rij_mu-%d-%d.txt", j/mbin, j%mbin);
	//			rij_muf[j]=fopen(buf,"w");
	//			snprintf(buf, sizeof buf, "/mnt/store1/awiegand/CovarianceCalculation/Data/own_all2pt_rik_mag-%d-%d.txt", j/mbin, j%mbin);
	//			rik_magf[j]=fopen(buf,"w");
	//			snprintf(buf, sizeof buf, "/mnt/store1/awiegand/CovarianceCalculation/Data/own_all2pt_rjk_mag-%d-%d.txt", j/mbin, j%mbin);
	//			rjk_magf[j]=fopen(buf,"w");
	//			snprintf(buf, sizeof buf, "/mnt/store1/awiegand/CovarianceCalculation/Data/own_all2pt_rik_mu-%d-%d.txt", j/mbin, j%mbin);
	//			rik_muf[j]=fopen(buf,"w");
	//			snprintf(buf, sizeof buf, "/mnt/store1/awiegand/CovarianceCalculation/Data/own_all2pt_rjk_mu-%d-%d.txt", j/mbin, j%mbin);
	//			rjk_muf[j]=fopen(buf,"w");
			}
#endif

		}

		~Integrals() {
			free(Ra);
			free(cx);
			free(c2);
			free(c3);
			free(c4);
			free(c3err);
			free(c4err);
			free(binct);
			free(binct3);
			free(binct4);

#ifdef ALLOUT
			for (int j=0; j<nbin*mbin; j++) {
				fclose(rij_magf[j]);
//				fclose(rij_muf[j]);
//				fclose(rik_magf[j]);
//				fclose(rjk_magf[j]);
//				fclose(rik_muf[j]);
//				fclose(rjk_muf[j]);
			}
#endif
//			rd.~RandomDraws();
	    }

	    void reset(){
	    	for (int j=0; j<nbin*mbin; j++) {
				Ra[j] = 0;
				cx[j] = 0;
				c2[j] = 0;
				binct[j] = 0;
				Raerr[j] = 0;
				cxerr[j] = 0;
				c2err[j] = 0;
			}
			for (int j=0; j<nbin*mbin*nbin*mbin; j++) {
				c3[j]=0;
				c4[j]=0;
				binct3[j] = 0;
				binct4[j] = 0;
				c3err[j]=0;
				c4err[j]=0;
			}
	    }

	    inline int getbin(Float r, Float mu){
	    	// Linearizes 2D indices
	       	return floor((r-rmin)/dr) *mbin + floor((mu-mumin)/dmu);
	    }

	    inline void second(int& bin, const Particle& pi, const Particle& pj, Float pro) {
	    	// Accumulates the two point integral C2
	    	// Note that the value of bin is returned
			Float rij_mag,rij_mu;
			Float rav,c2v,cxv, weight;

//			Bin will be reused later
	        cleanup_l(pi.pos, pj.pos, rij_mag, rij_mu);
	        bin=getbin(rij_mag,rij_mu);

	        if ((bin>=0) && (bin<nbin*mbin)){
	        	weight=pi.w*pj.w;
		        xiij=cf->xi(rij_mag,rij_mu);
	        	rav=weight / pro;
	        	c2v= weight*weight *(1+xiij) / pro;
	        	cxv = weight*xiij / pro;

	        	Ra[bin]+= rav;
	        	c2[bin]+= c2v;
	        	cx[bin]+= cxv;
	        	binct[bin]++;
	        	Raerr[bin]+= pow(rav,2);
//	        	c2err[bin]+= pow(c2v,2);
	        	c2err[bin]+= weight*weight / pro * pow(xiij,2); //Misuse to accumulate xi^2 term
	        	cxerr[bin]+= pow(cxv,2);
	        }else{
	        	bin=-1;
	        }
#ifdef ALLOUT
	        rij_magt=rij_mag;
	        rij_mut=rij_mu;
#endif
	    }

	    inline void third(const int& bin, Float& xijk, const Particle& pi, const Particle& pj, const Particle& pk, Float p3) {
	    	// Accumulates the three-point integral C3
	    	// Returns the correlation function between point j and k
	    	int binik;
	    	Float c3l;
			Float rik_mag,rik_mu, rjk_mag,rjk_mu, xiik;

			cleanup_l(pi.pos, pk.pos, rik_mag, rik_mu);
			binik=getbin(rik_mag, rik_mu);
			xiik=cf->xi(rik_mag, rik_mu); //Only temporary for xi^2; can be commented out if c3err is used for error as intended

//			Correlation function will be reused in fourth
			cleanup_l(pj.pos, pk.pos, rjk_mag, rjk_mu);
	        xijk=cf->xi(rjk_mag, rjk_mu);

//			Checking if ik bin is in chosen range and assuming that "bin" is the correct ij bin
			if ((binik<nbin*mbin) && (binik>=0)){
				c3l=pi.w*pi.w*pj.w*pk.w/p3 * xijk;
				c3[bin*nbin*mbin+binik]+=c3l;
//				c3err[bin*nbin*mbin+binik]+=pow(c3l,2);
				c3err[bin*nbin*mbin+binik]+=pi.w*pi.w*pj.w*pk.w/p3 * 0.5*(xijk*xiik+xijk*xiij); //Misuse to accumulate xi^2 term
				binct3[bin*nbin*mbin+binik]++;
			}

#ifdef ALLOUT
			rjk_magt=rjk_mag;
	        rjk_mut=rjk_mu;
#endif
		}

	    inline void fourth(const int& bin, Float& xijk, const Particle& pi, const Particle& pj, const Particle& pk, const Particle& pl, Float p4) {
	    	// Accumulates the four-point integral C4
	    	// Input: The ij bin, the correlation function between jk, the four particle positions
	    	int binkl;
	    	Float c4l,xiil;
	    	Float ril_mag,ril_mu, rkl_mag,rkl_mu;

	    	cleanup_l(pk.pos, pl.pos, rkl_mag, rkl_mu);
	    	binkl=getbin(rkl_mag,rkl_mu);

//			Checking if kl bin is in chosen range and assuming that "bin" is the correct ij bin
	    	if ((binkl<nbin*mbin) && (binkl>=0)){

		    	cleanup_l(pi.pos, pl.pos, ril_mag, ril_mu);
		    	xiil=cf->xi(ril_mag, ril_mu);

	    		c4l=pi.w*pj.w*pk.w*pl.w/p4  * xijk  * xiil;
	    		c4[bin*nbin*mbin+binkl]+=c4l;
//	    		c4err[bin*nbin*mbin+binkl]+=pow(c4l,2);
	    		c4err[bin*nbin*mbin+binkl]+=pi.w*pj.w*pk.w*pl.w/p4; //Misuse of c4err to test normalization against 2pt randoms
	    		binct4[bin*nbin*mbin+binkl]++;
	    	}

#ifdef ALLOUT
#pragma omp critical
	    		fprintf(rij_magf[bin],"%d %e %e %e %e %e %e %e %e %e %e\n",binkl,rij_mut,rkl_mu,rjk_mut,ril_mu,rij_magt,rkl_mag,rjk_magt,ril_mag, pi.w*pj.w*pk.w*pl.w/p4 * xijk  * xiil, pi.w*pj.w*pk.w*pl.w/p4);
#endif
	    }

	    void sum_ints(Integrals* ints) {
	    	// Add the values accumulated in ints to the corresponding internal sums
			for(int i=0;i<nbin*mbin;i++){
				Ra[i]+=ints->Ra[i];
				c2[i]+=ints->c2[i];
				cx[i]+=ints->cx[i];
				binct[i]+=ints->binct[i];
				Raerr[i]+=ints->Raerr[i];
				c2err[i]+=ints->c2err[i];
				cxerr[i]+=ints->cxerr[i];
			}

			for(int i=0;i<nbin*mbin;i++){
				for(int j=0;j<nbin*mbin;j++){
					c3[i*nbin*mbin+j]+=ints->c3[i*nbin*mbin+j];
					c4[i*nbin*mbin+j]+=ints->c4[i*nbin*mbin+j];
					c3err[i*nbin*mbin+j]+=ints->c3err[i*nbin*mbin+j];
					c4err[i*nbin*mbin+j]+=ints->c4err[i*nbin*mbin+j];
					binct3[i*nbin*mbin+j]+=ints->binct3[i*nbin*mbin+j];
					binct4[i*nbin*mbin+j]+=ints->binct4[i*nbin*mbin+j];
				}
			}
		}

	    void sum_total_counts(uint64& acc2, uint64& acc3, uint64& acc4){
	    	// Add local counts to total bin counts in acc2-4
//	    	acc2=0,acc3=0,acc4=0;
	    	for (int i=0; i<nbin*mbin; i++) {
	    		acc2+=binct[i];
				for (int j=0; j<nbin*mbin; j++) {
					acc3+=binct3[i*nbin*mbin+j];
					acc4+=binct4[i*nbin*mbin+j];
				}
			}
	    }

	    void normalize(long n,int np,int ngal, int n3, int n4) {
	    	// Normalize the acccumulated integrals; partially already done by normalization of the proposal function of the cubes selected
	    	// np is the number of random particles used, ngal is the number of galaxies in the survey
	    	double corrf=pow((double)np/ngal,2); // Correction factor for the different densities of randoms and points

	    	for(int i=0;i<nbin*mbin;i++){

//	    		printf("Predefined Integrals %2d %le +- %le   %le +- %le   %le +- %le\n",
//						i, Ra[i], Raerr[i], cx[i], cxerr[i], c2[i], c2err[i]);

	    		Ra[i]/=((Float)n*corrf);
	    		c2[i]/=(Ra[i]*Ra[i]*(Float)n*corrf);
				cx[i]/=(Ra[i]*(Float)n*corrf);

// TODO: Figure out how to normalize the quadratic counts of err to arrive at an estimate for the integration error
//				Raerr[i]=sqrt((Raerr[i]/((Float)n)-pow(Ra[i]*corrf,2))/(Float)(n-1))/corrf;
///*wrong->*/		c2err[i]=sqrt((c2err[i]/(pow((Ra[i]*Ra[i]*corrf*corrf),2)*(Float)n)-pow(c2[i]/corrf,2))/(Float)(n-1))/corrf;
///*wrong->*/		cxerr[i]=sqrt((cxerr[i]/(pow((Ra[i]*pow((double)np/ngal,2)),2)*(Float)n)-pow(cx[i],2))/(Float)(n-1));

//				printf("Postdefined Integrals %2d %le +- %le   %le +- %le   %le +- %le\n",
//						i, Ra[i], Raerr[i], cx[i], cxerr[i], c2[i], c2err[i]);


	    	}

	    	for(int i=0;i<nbin*mbin;i++){
	    		for(int j=0;j<nbin*mbin;j++){
	    			c3[i*nbin*mbin+j]/=(n3*Ra[i]*Ra[j]*(Float)n*pow((double)np/ngal,3));
	    			c4[i*nbin*mbin+j]/=(n3*n4*Ra[i]*Ra[j]*(Float)n*pow((double)np/ngal,4));
	    		}
			}
	    }

	    void report_integrals() {
	///*
		for (int j=0; j<nbin*mbin; j++) {
		    printf("C2 %2d %le +- %le   %le +- %le   %le +- %le\n",
				j, Ra[j], Raerr[j], cx[j], cxerr[j], 2*c2[j], 2*c2err[j]);
		}
		printf("\n");


		for(int i=0;i<nbin*mbin;i++){
			printf("C3 %2d ",i);
			for(int j=0;j<nbin*mbin;j++){
				printf("%e ", 4*c3[i*nbin*mbin+j]);
			}
			printf("\n");
		}
		printf("\n");

		for(int i=0;i<nbin*mbin;i++){
			printf("C4 %2d ",i);
			for(int j=0;j<nbin*mbin;j++){
				printf("%e ",2*c4[i*nbin*mbin+j]);
			}
			printf("\n");
		}
		printf("\n");

		for(int i=0;i<nbin*mbin;i++){
			printf("E3 %2d ",i);
			for(int j=0;j<nbin*mbin;j++){
				printf("%e ",4*c3err[i*nbin*mbin+j]);
			}
			printf("\n");
		}
		printf("\n");

		for(int i=0;i<nbin*mbin;i++){
			printf("E4 %2d ",i);
			for(int j=0;j<nbin*mbin;j++){
				printf("%e ",2*c4err[i*nbin*mbin+j]);
			}
			printf("\n");
		}
		printf("\n");

		for (int i=0; i<nbin; i++) {
			printf("N2 %2d ",i);
			for (int j=0; j<mbin; j++) {
				printf("%lld ", binct[i*mbin+j]);
			}
			printf("\n");
		}
		printf("\n");

		for (int i=0; i<nbin*mbin; i++) {
			printf("N3 %2d ",i);
			for (int j=0; j<nbin*mbin; j++) {
				printf("%lld ", binct3[i*nbin*mbin+j]);
			}
			printf("\n");
		}
		printf("\n");

		for (int i=0; i<nbin*mbin; i++) {
			printf("N4 %2d ",i);
			for (int j=0; j<nbin*mbin; j++) {
				printf("%lld ", binct4[i*nbin*mbin+j]);
			}
			printf("\n");
		}
		printf("\n");

		fflush(NULL);

	//*/
	/*
		printf("\n{");
		printf("%.0f", xi0[0]);
		for (int j=1; j<nbin; j++) {
		    printf(",%.0f", xi0[j]);
		}
		printf("}\n");
	*/

	    }

	  private:
	    void cleanup_l(Float3 p1,Float3 p2,Float& norm,Float& mu){
	    	Float3 pos=p1-p2;
	    	norm = pos.norm();
	    	// If the input correlation function had only radial information and no mu bins
	    	// fill the dummy mu bins by 0.5
	    	if(rad){
	    		mu=0.5;
	    	}
	    	else{
#ifndef PERIODIC
				Float3 los=p1+p2; // No 1/2 as normalized anyway below
	//	    	if(!box){
					//Absolute value as mu bins defined in terms of |mu|
					mu = fabs(pos.dot(los)/norm/los.norm());
	//	    	}
#else
	    		// In the periodic case use z-direction for mu
				mu = fabs(pos.z/norm);
#endif
	    	}
	    }
	};

/*--------------------------------------------------------------------*/

void kofnip(int* cellsel, int n, int k, gsl_rng* rng) {

	/* -------------------------------------------------------- */
	/* selects k points of an ensemble of n points              */
	/* -------------------------------------------------------- */

	int tmp, r;

	assert(k<=n);

	for (int i=0;i<k;i++) {
		r = i + (int) floor(gsl_rng_uniform(rng) * (n - i));
		tmp = cellsel[i];
		cellsel[i] = cellsel[r];
		cellsel[r]=tmp;
	}

	return;

};


/*--------------------------------------------------------------------*/

void kofn(int* cellsel, int n, int k,Integrals * integs, gsl_rng* rng) {

	/* -------------------------------------------------------- */
	/* selects k points of an ensemble of n points              */
	/* -------------------------------------------------------- */

	int* mapping;
	int max_rand, ind, r;
	int i;

	assert(k<=n);

	mapping = (int*)malloc(sizeof(int)*n);

	for (i = 0; i < n; i++) {
		mapping[i] = i;
	}

	max_rand = n - 1;


	for (i=0;i<k;i++) {
		r = (int) floor(gsl_rng_uniform(rng) * (max_rand + 1));
		ind = mapping[r];

		cellsel[i] = ind;

		mapping[r] = mapping[max_rand];
		max_rand = max_rand - 1;
	}

	free(mapping);

	return;

};

integer3 relid(int cin, int maxsep){

	integer3 res;
	res.z = cin%(2*maxsep+1)-maxsep;
	cin = cin/(2*maxsep+1);
	res.y = cin%(2*maxsep+1)-maxsep;
	res.x = cin/(2*maxsep+1)-maxsep;


	return res;
}

STimer TotalTime;

// ====================  Computing the integrals ==================


// Available sampling strategies:
#define ALL 0   // Use all cells in the grid up to a certain radius
#define FLAT 1	// Select cells equiprobably up to a certain radius
#define SEL 2	// Select cells weighted by the distance from the starting cell up to the radius defined by the input parameter xicutoff

void compute_integral(Grid *grid, Parameters *par) {

	/* -------------------------------------------------------- */
	/* Calculates the integral based on the random particles 	*/
	/* in the grid              								*/
	/* -------------------------------------------------------- */


// ---- Defines the sampling strategy ----

#define SECSAMP ALL		// Strategy for the selection of the second cell ...
#define THISAMP SEL		// the third cell ...
#define FOUSAMP SEL		// and the fourth cell

	// In the following define the number of cells selected at each stage
	uint n1=grid->nf; // Currently use all grid cells as first grid cell once. Could use less than all filled cells here
	uint n2=100;      // Sets the number of cells used as secondary cells. Not effective if strategy is ALL
	uint n3=1;		  // Number of cells selected as third cell for each pair of primary and secondary cell
	uint n4=1;		  // Number of cells selected as fourth cell for each triple
					  // Setting those to 1 currently also avoids selecting a quad twice

// ---- ---- ---- ---- ---- ---- ---- ----


// ---- Defines intermediate output ----
    int sumupdint=(int)ceil((double)n1/par->nthread/1000); //Update sum in each thread every 0.1% of the runtime; should be much smaller than sumrepint
    int sumrepint=(int)ceil((double)n1/20);	//Report sum every x% of the full runtime; set to > n1 to avoid intermediate results
    int interval=10; // Interval when to report estimates of the remaining runtime. Changed during evaluation
// ---- ---- ---- ---- ---- ---- ---- ----

	int maxsep = ceil(par->rmax/grid->cellsize);	// The maximal distance for the pairs is determined by the upper limit of the covariance matrix to be calculated
	int maxsep3=ceil(par->xicutoff/grid->cellsize);	// Not used for SEL. The maximal distance for the third cell is determined from the xicutoff
	int maxsep4=maxsep;								// Not used for SEL.

	int aneigh=(int)pow(2*maxsep+1,3);		// Number of all neighbour cells up to max radius
	int aneigh3=(int)pow(2*maxsep3+1,3);
	int aneigh4=(int)pow(2*maxsep4+1,3);

    uint64 cnt2=0,cnt3=0,cnt4=0;
    uint64 acc2=0,acc3=0,acc4=0;
//    uint64 mean=(uint64)floor((double)grid->np/grid->nf);

    int sumct=0; //Track how many iterations went into intermediate result
    double begtime; // Define global variable for the start time

    // Set the factor between two seeds of the rng of the threads
    // For deterministic rng seeds, set "steps" manually
    // Could omit ("/dev/urandom") for portability, but there seem to be implementations where random_device is deterministic
    std::random_device urandom("/dev/urandom");
    std::uniform_int_distribution<unsigned int> dist(1, std::numeric_limits<unsigned int>::max());
    unsigned long int steps = dist(urandom);

// ---- Initialize correlation function and selection grid ----
    STimer initial;
    initial.Start();

    CorrelationFunction *cf=new CorrelationFunction(par->corname,par->mbin,par->mumax-par->mumin);
	RandomDraws *rd=new RandomDraws(cf,par, NULL, 0);
    Integrals sumint(par,cf);

    // Shuffles the filled cells. Only needed if we do not use all cells in the outermost loop
    std::vector<int> cellsranord(grid->filled, grid->filled + grid->nf);
    std::random_shuffle ( cellsranord.begin(), cellsranord.end() );

    initial.Stop();
    fprintf(stderr,"Init time: %g s\n",initial.Elapsed());
// ---- ---- ----

    printf("# Filled cells: %d\n",grid->nf);
    printf("# All points in use: %d\n",grid->np);
    printf("# Max points in one cell %d\n",grid->maxnp);
    fflush(NULL);

    TotalTime.Start();

    gsl_rng_env_setup();

#ifdef OPENMP
	cpu_set_t mask[par->nthread+1];
	int tnum=0;
	sched_getaffinity(0, sizeof(cpu_set_t),&mask[par->nthread]);
	fprintf(stderr,"# CPUs used are: ");
	for(int ii=0;ii<64;ii++){
		if(CPU_ISSET(ii, &mask[par->nthread])){
			fprintf(stderr,"%d ",ii);
			CPU_ZERO(&mask[tnum]);
			CPU_SET(ii,&mask[tnum]);
			tnum++;
		}
	}
	fprintf(stderr,"\n");

// Begin parallel section before omp for loop to define thread local variables
#pragma omp parallel firstprivate(sumupdint,sumrepint,maxsep,n1,steps) reduction(+:cnt2,cnt3,cnt4) shared(TotalTime,begtime,sumint,sumct)
	{
	// Decide which thread we are in
	int thread = omp_get_thread_num();
	int loopct = 0;
    assert(omp_get_num_threads()<=par->nthread);
    if (thread==0) printf("# Running on %d threads.\n", omp_get_num_threads());
#else
	int thread = 0;
    printf("# Running single threaded.\n");
#endif

    Integrals locint(par,cf); // Accumulates the integral contribution of each thread
	gsl_rng* locrng = gsl_rng_alloc( gsl_rng_default ); // One rng per thread
	gsl_rng_set(locrng,steps*(thread+1)); // seeds of the different rngs "steps" apart

    uint64 locacc2=0,locacc3=0,locacc4=0;

    int mnp=grid->maxnp;
    int pln,sln,tln,fln;
    Particle *prim_list;
    Particle *sec_list;
    Particle *thi_list;
    Particle *fou_list;
    Cell zero({0,0});

    Float* xi;
    int* bin;

    bool changeij;
    bool changek;

    int ec=0;
    ec+=posix_memalign((void **) &prim_list, PAGE, sizeof(Particle)*mnp); // Not sure if memalign has a benefit over malloc here
    ec+=posix_memalign((void **) &sec_list, PAGE, sizeof(Particle)*mnp);
    ec+=posix_memalign((void **) &thi_list, PAGE, sizeof(Particle)*mnp);
    ec+=posix_memalign((void **) &fou_list, PAGE, sizeof(Particle)*mnp);

    ec+=posix_memalign((void **) &bin, PAGE, sizeof(int)*mnp*mnp);
    ec+=posix_memalign((void **) &xi, PAGE, sizeof(Float)*mnp*mnp);

    assert(ec==0);

    for(int ll=0;ll<mnp*mnp;ll++){
    	bin[ll]=9999999;
    }

#if SECSAMP == FLAT
		int* cellsel2 = (int*)malloc(sizeof(int)*aneigh);
		for(int ll=0;ll<aneigh;ll++){
			cellsel2[ll]=ll;
		}
#endif
#if THISAMP == FLAT
		int* cellsel3 = (int*)malloc(sizeof(int)*aneigh3);
		for(int ll=0;ll<aneigh3;ll++){
			cellsel3[ll]=ll;
		}
#endif
#if FOUSAMP == FLAT
		int* cellsel4 = (int*)malloc(sizeof(int)*aneigh4);
		for(int ll=0;ll<aneigh4;ll++){
			cellsel4[ll]=ll;
		}
#endif

//    for(int i=0;i<=1000;i++){
//    	for(int j=0;j<=1000;j++){
//    	fprintf(stdout,"%e ",locint.cf->xi(i,j*0.001));
//    	}
//    	fprintf(stdout,"\n");
//    }
//    fprintf(stdout,"\n");
//
//    exit(0);


// ---- Calculation starts here ----

// Loop over primary cells.
#ifdef OPENMP
#pragma omp	for schedule(dynamic)
#endif
   for (int n=0; n<n1; n++) {
#ifdef OPENMP
	// Manually pin threads to CPUs (necessary on Odyssey)
//	if(CPU_COUNT(&mask[thread])!=0){
//		sched_setaffinity(0, sizeof(cpu_set_t),&mask[thread]);
//		CPU_ZERO(&mask[thread]);
//	}
	//fprintf(stderr,"# Thread %d on cpu %d calculates cell %d and mask is ",omp_get_thread_num(),sched_getcpu(),n);
	//sched_getaffinity(0, sizeof(cpu_set_t),&mask[thread]);
	//for(int ii=0;ii<64;ii++){
	//	if(CPU_ISSET(ii, &mask[thread]))
	//		fprintf(stderr,"%d ",ii);
	//}
	//fprintf(stderr,".\n");
#endif

        // Determine first cell; Cells randomly ordered so we can use only a limited number just by stopping the loop
		int m=cellsranord[n];
		integer3 prim_id = grid->cell_id_from_1d(m);
		Cell primary = grid->c[m];
		double p1=n1/grid->nf;

		if(primary.np==0) continue; // If cell empty skip (should not happen any more)

		// Copy global list to local list to enhance memory access pattern
		pln=0;
		for (int i = primary.start; i<primary.start+primary.np; i++,pln++) {
			prim_list[pln]=grid->p[i];
		}

//		Determine second position
		double p2;
		integer3 delta;

// ---- Depending on sampling mode select second cell ----
#if SECSAMP == FLAT
		int fin2 = fmin(aneigh,n2);
		kofnip(cellsel2,aneigh,fin2,locrng); 	// Draw second cell equiprobably from all neighbouring cells
		p2=1./aneigh*fin2;						// fin2 points selected out of aneigh
        for (int ci=0;ci<fin2;ci++){
    		delta=relid(cellsel2[ci],maxsep);
#elif SECSAMP == ALL
		for (delta.x = -maxsep; delta.x <= maxsep; delta.x++)	// Iterate through all other cells up to maxsep
		for (delta.y = -maxsep; delta.y <= maxsep; delta.y++)
		for (delta.z = -maxsep; delta.z <= maxsep; delta.z++){
			p2=1.; // All points selected, so probability 1
#elif SECSAMP == SEL
			for (uint ci=0;ci<n2;ci++){
	    		delta=rd->random_cubedraw(&p2); // Draw second cell from cube around first; May draw same cell twice
#endif

	    	integer3 sec_id = prim_id+delta;
	       	int id=grid->test_cell(sec_id);
			if(id<0) continue;					// If the cell is outside the bounding box
			Cell sec = grid->c[id];

			if(sec.np==0) continue; // If second cell empty there will be neither pairs, triples nor quads
									// Can't do that for the inner loops (thi and fou) as we need all pairs!

			// Copy global list to local list to enhance memory access pattern
			sln=0;
			for (int i = sec.start; i<sec.start+sec.np; i++,sln++) {
				sec_list[sln]=grid->p[i];
#ifdef PERIODIC
				sec_list[sln].pos+= grid->cell_sep(delta); // In periodic case particles positions are relative to cell center
#endif
			}

			p2*=p1;		   // Accumulates selection probability
			changeij=true; // Tracks if we changed the pair ij by selecting a new primary and/or secondary cell


//			Determine third position
			integer3 delta3;
	    	double p3;
// ---- Depending on sampling mode select third cell ----
#if THISAMP == SEL
				for (uint ck=0;ck<n3;ck++){
					delta3=rd->random_xidraw(locrng,&p3);
					integer3 thi_id=sec_id+delta3;			// select with respect to second cell
#elif THISAMP == FLAT
	    		int fin3 = fmin(aneigh3,n3);
	    		kofnip(cellsel3,aneigh3,fin3,locrng);
				for (int ck=0;ck<fin3;ck++){
		    		delta3=relid(cellsel3[ck],maxsep3);
		       		integer3 thi_id=sec_id+delta3;
					p3=1./aneigh3*fin3;
#elif THISAMP == ALL
				int ck=0;
				for (delta3.x = -maxsep3; delta3.x <= maxsep3; delta3.x++)
				for (delta3.y = -maxsep3; delta3.y <= maxsep3; delta3.y++)
				for (delta3.z = -maxsep3; delta3.z <= maxsep3; delta3.z++,ck++){
					integer3 thi_id = sec_id+delta3;
					p3=1.;
//					integer3 thi_id = prim_id+delta3;
#endif
					p3*=p2;
					id=grid->test_cell(thi_id);
		       		Cell thi;
		       		if(id<0)
		       			thi = zero;
		       		else
		       			thi = grid->c[id];

//					Can't skip empty cells at this stage because we need all pairs, independently of third draw
//					if(thi.np==0) continue;

					// Copy global list to local list to enhance memory access pattern
		       		tln=0;
					for (int i = thi.start; i<thi.start+thi.np; i++,tln++) {
						thi_list[tln]=grid->p[i];
#ifdef PERIODIC
//						thi_list[tln].pos += grid->cell_sep(delta3); //From cell 1
						thi_list[tln].pos += grid->cell_sep(delta) + grid->cell_sep(delta3); //From cell 2
#endif
					}
					changek=true; // Tracks if third cell has changed

//			Determine fourth position
					double p4;
					integer3 delta4;
// ---- Depending on sampling mode select fourth cell ----
#if FOUSAMP == SEL
					for (uint cl=0;cl<n4;cl++){
						// Determine fourth position
						delta4=rd->random_xidraw(locrng,&p4);
						integer3 fou_id = prim_id+delta4;
	//					delta4=rd->random_cubedraw(&p4);
	//					integer3 fou_id = thi_id+delta4;
#elif FOUSAMP == ALL
					int cl=0;

					for (delta4.x = -maxsep4; delta4.x <= maxsep4; delta4.x++)
					for (delta4.y = -maxsep4; delta4.y <= maxsep4; delta4.y++)
					for (delta4.z = -maxsep4; delta4.z <= maxsep4; delta4.z++,cl++){
						integer3 fou_id = thi_id+delta4;	// Samples around third cell!
						p4=1.;

#elif FOUSAMP == FLAT
					int fin4 = fmin(aneigh4,n4);
					kofnip(cellsel4,aneigh4,fin4,locrng);
					for (int cl=0;cl<fin4;cl++){
						delta4=relid(cellsel4[cl],maxsep4);
						integer3 fou_id = thi_id+delta4;	// Samples around third cell!
						p4=1./aneigh4*fin4;
#endif
					p4*=p3;
					id=grid->test_cell(fou_id);
					Cell fou;
					if(id<0)
						fou = zero;
					else
						fou = grid->c[id];

//					Not possible with current loopstructure (would impact 2pt part)
//					if(fou.np==0) continue;

					// Copy global list to local list to enhance memory access pattern
					fln=0;
					for (int i = fou.start; i<fou.start+fou.np; i++,fln++) {
						fou_list[fln]=grid->p[i];
#ifdef PERIODIC
#if FOUSAMP == SEL
						fou_list[fln].pos += grid->cell_sep(delta4); // delta4 from cell 1
#else
						fou_list[fln].pos += grid->cell_sep(delta) + grid->cell_sep(delta3) + grid->cell_sep(delta4); // delta4 from cell 3
#endif
#endif
					}

// 		---- Now we loop over all particles in the chosen cells ----
					for (int i = 0; i<pln; i++) {
//						if(sec.np==0) cnt2+=mean;
						for (int j = 0; j<sln; j++) {
		    				if (sec.start==primary.start&&i==j) continue;   // Exclude self-count

		    				if(changeij){
		    					locint.second(bin[sln*i + j],prim_list[i],sec_list[j],p2); // bin is io here!
		    					cnt2++;
		    				}
		    				if(bin[sln*i + j]<0) continue; // If current pair outside of binning range: skip
//		    				continue;

	    				for (int k = 0; k<tln; k++) {
							if (thi.start==sec.start&&j==k) continue;   // Exclude self-count
							if (thi.start==primary.start&&i==k) continue;

							if(changek){
								locint.third(bin[sln*i + j],xi[tln*j + k],prim_list[i],sec_list[j],thi_list[k],p3); // xi is io here!
								cnt3++;
							}
//		    				continue;

						for (int l = 0; l<fln; l++) {
							if (fou.start==primary.start&&i==l) continue;   // Exclude self-count
							if (fou.start==sec.start&&j==l) continue;
							if (fou.start==thi.start&&k==l) continue;

							locint.fourth(bin[sln*i + j],xi[tln*j + k],prim_list[i],sec_list[j],thi_list[k],fou_list[l],p4);
							cnt4++;

						} // Done with fourth cell
	    				} // Done with third cell
						} // Done with second cell
					} // Done with first cell
//			cnt4++;
			changeij=false;
			changek=false;

			if(tln==0) break; // Now that the pairs are taken into account, if the third cell is empty do no longer continue to select fourth cells to the empty third cell

			} // Done with fourth choice
		} // Done with third choice
		} // Done with second choice


//  Update global variables from the accumulated thread variables
	if (loopct%sumupdint==0&&loopct!=0) {
#pragma omp critical(accsumint)		// Critical section to prevent several threads trying to copy their results to the global variable simultaneously
    {
    	sumint.sum_ints(&locint);
        sumct+=sumupdint;
    }
	}

	//Print estimated remaining time and current speed
    if (n%interval==0) {
#pragma omp critical
    	{

    	if(interval<0.05*n1)
    		interval*=5;
    	else
    		interval=(int)floor(0.1*n1);

    	TotalTime.Stop();  // Interrupt timing to access .Elapsed()

    	if (n==0) begtime=TotalTime.Elapsed();

    	uint timeleft=(uint)(TotalTime.Elapsed()-begtime)*(1./n*grid->nf-1.);
    	double runtime=TotalTime.Elapsed()-begtime;
    	uint64 tmplocacc2=locacc2, tmplocacc3=locacc3, tmplocacc4=locacc4;

    	locint.sum_total_counts(locacc2, locacc3, locacc4);

    	fprintf(stderr,"# Finished cell %d of %d after %g s. "
    			"Estimated time left: %2.2d:%2.2d:%2.2d hms "
    			"i.e. %d s. Thread %d: %ld %e %e quads/sec.\n",
    			n, grid->nf,runtime,timeleft/3600,timeleft/60%60,
				timeleft%60,timeleft,thread,cnt4,locacc4/runtime,cnt4/runtime);

    	locacc2=tmplocacc2, locacc3=tmplocacc3, locacc4=tmplocacc4;

    	TotalTime.Start();
    }
    }

	if (n%sumrepint==0&&n!=0){
#pragma omp critical(accsumint)   // Reading sumint needs the same critical section as writing to sumint
		{
		printf("# Used cells: %d/%d\n",sumct,n1);
		printf("UCF: %g\n",(double)sumct/n1); // Used cell fraction, needed to normalize intermediate result
		sumint.report_integrals();
		}
	}

//	Reset local variable whenever contents had been written to sumint
	if (loopct%sumupdint==0&&loopct!=0) {
    	locint.sum_total_counts(locacc2, locacc3, locacc4);
        locint.reset();
	}

    loopct++;

	} // Done with first choice; End of omp for loop

	free(prim_list);
	free(sec_list);
	free(thi_list);
	free(fou_list);
	free(bin);
	free(xi);
	gsl_rng_free(locrng);

#if SECSAMP == FLAT
        free(cellsel2);
#endif
#if THISAMP == FLAT
        free(cellsel3);
#endif
#if FOUSAMP == FLAT
        free(cellsel4);
#endif

// After the end of the loop add all remaining contents of the local threads
#ifdef OPENMP
#pragma omp critical(accsumint)
    {
		sumint.sum_ints(&locint);
// TODO: sumct sometimes wrong. Maybe because of loopct++ after  sumct+=sumupdint; ? Should not be a problem anyway because sumct not used outside the loop
	    sumct+=loopct%sumupdint;
    }
	} // End of omp parallel pragma
#endif

    fprintf(stderr,"Parallely evaluated cells: %d\nTotal cells: %d\n",sumct,n1);

    TotalTime.Stop();

    fprintf(stderr,"\nTime for loop: %g s\n",TotalTime.Elapsed());

    // Normalization of the accumulated results. n3 and n4 should probably better go into the probabilities p3 and p4 instead of showing up here
	sumint.normalize(/*cnt2(long)grid->np*((long)grid->np-1)*/1.,grid->np,par->nofznorm,n3,n4);
	sumint.sum_total_counts(acc2, acc3, acc4);

//	Print some statistics
    printf("# We tried  %lld pairs, %lld triples and %lld quads within %f.\n", cnt2, cnt3, cnt4, par->rmax);
    printf("# We accepted %lld pairs, %lld triples and %lld quads.\n", acc2, acc3, acc4);
    printf("# Acceptance ratios are %f for pairs, %f for triples and %f for quads.\n",
    		(double)acc2/cnt2, (double)acc3/cnt3, (double)acc4/cnt4);
    printf("# Average of %f pairs per primary particle.\n",
    		(Float)cnt2/grid->np);
    float x = par->rmax/grid->boxsize;
    float expected = grid->np * (4*M_PI/3.0)*grid->np*x*x*x;
    printf("# In a periodic box we would expect %1.0f pairs, off by %f.\n", expected, cnt2/expected);

//  Print the result
    sumint.report_integrals();

    printf("\n# Time for initialization: %g s \n"
    		"# Time for main loop: %g s\n",initial.Elapsed(),TotalTime.Elapsed());
    fflush(NULL);
    return;
}



// ====================  The Driver ===========================

Particle *make_particles(Float boxsize, int np) {
    // Make np random particles
    srand48(1);      // For reproducibility
    Particle *p = (Particle *)malloc(sizeof(Particle)*np);
    for (int j=0; j<np; j++) {
        p[j].pos.x = drand48()*boxsize;
        p[j].pos.y = drand48()*boxsize;
        p[j].pos.z = drand48()*boxsize;
        p[j].w = 1.0;
    }
    printf("# Done making %d random particles, periodically distributed.\n", np);
    return p;
}

Particle *read_particles(Float rescale, int *np, const char *filename, const int rstart, uint64 nmax) {
    // This will read particles from a file, space-separated x,y,z,w
    // Particle positions will be rescaled by the variable 'rescale'.
    // For example, if rescale==boxsize, then inputing the unit cube will cover the periodic volume
    char line[1000];
    int j=0,n=0;
    FILE *fp;
    int stat;
    double tmp[6];
    fp = fopen(filename, "r");
    if (fp==NULL) {
        fprintf(stderr,"File %s not found\n", filename); abort();
    }
    // Count lines to construct the correct size
    while (fgets(line,1000,fp)!=NULL&&(uint)n<nmax) {
        if (line[0]=='#') continue;
        if (line[0]=='\n') continue;
        n++;
    }
    rewind(fp);
    *np = n;
    Particle *p = (Particle *)malloc(sizeof(Particle)*n);
    printf("# Found %d particles from %s\n", n, filename);
    printf("# Rescaling input positions by factor %f\n", rescale);
    while (fgets(line,1000,fp)!=NULL&&j<n) {
        if (line[0]=='#') continue;
        if (line[0]=='\n') continue;
        stat=sscanf(line, "%lf %lf %lf %lf %lf %lf", tmp, tmp+1, tmp+2, tmp+3, tmp+4, tmp+5);

        if (stat<3) {
        	fprintf(stderr,"Particle %d has bad format\n", j); // Not enough coordinates
        	abort();
        }

        p[j].pos.x = tmp[0]*rescale;
        p[j].pos.y = tmp[1]*rescale;
        p[j].pos.z = tmp[2]*rescale;

        // If there are 4 or 6 entries per line get the weights from line 4 or 6
        // Otherwise fill the weights with 1 or -1 depending on the value of rstart
        // For grid_covariance rstart is typically not used
        if(stat!=6&&stat!=4)
		   if(rstart>0&&j>=rstart)
			   p[j].w = -1.;
		   else
			   p[j].w = 1.;
		else{
		   if(rstart>0&&j>=rstart)
			   p[j].w = -tmp[stat-1];
		   else
			   p[j].w = tmp[stat-1];
		}
		j++;
    }
    fclose(fp);
    printf("# Done reading the particles\n");
    return p;
}

bool check_bounding_box(Particle *p, int np, Float boxsize, Float rmax, Float3& pmin) {
    // Check that the bounding box is reasonable and return minimal position
    Float3 pmax;
    bool box=false;
    pmin.x = pmin.y = pmin.z = 1e30;
    pmax.x = pmax.y = pmax.z = -1e30;
    for (int j=0; j<np; j++) {
        pmin.x = fmin(pmin.x, p[j].pos.x);
        pmin.y = fmin(pmin.y, p[j].pos.y);
        pmin.z = fmin(pmin.z, p[j].pos.z);
        pmax.x = fmax(pmax.x, p[j].pos.x);
        pmax.y = fmax(pmax.y, p[j].pos.y);
        pmax.z = fmax(pmax.z, p[j].pos.z);
    }
    printf("# Range of x positions are %6.2f to %6.2f\n", pmin.x, pmax.x);
    printf("# Range of y positions are %6.2f to %6.2f\n", pmin.y, pmax.y);
    printf("# Range of z positions are %6.2f to %6.2f\n", pmin.z, pmax.z);
    Float3 prange = pmax-pmin;
    Float         biggest = prange.x;
    biggest = fmax(biggest, prange.y);
    biggest = fmax(biggest, prange.z);
    printf("# Biggest range is %6.2f\n", biggest);
    if (biggest>boxsize*1.001)
	printf("#\n# WARNING: particles will overlap on period wrapping!\n#\n");
    if (biggest+rmax<boxsize*0.6)
	printf("#\n# WARNING: box periodicity seems too generous, will hurt grid efficiency!\n#\n");

    if (prange.x>0.99*biggest && prange.y>0.99*biggest && prange.z>0.99*biggest) {
        // Probably using a cube of inputs, intended for a periodic box
    	box=true;
#ifndef PERIODIC
    	fprintf(stderr,"#\n# WARNING: cubic input detected but you have not compiled with PERIODIC flag!\n#\n");
    	printf("#\n# WARNING: cubic input detected but you have not compiled with PERIODIC flag!\n#\n");
#endif
	if (biggest<0.99*boxsize)
	    printf("#\n# WARNING: cubic input detected, but smaller than periodicity!\n#\n");
    } else {
        // Probably a non-periodic input
    	box=false;
#ifdef PERIODIC
    	fprintf(stderr,"#\n# WARNING: non-cubic input detected but you have compiled with PERIODIC flag!\n#\n");
    	printf("#\n# WARNING: non-cubic input detected but you have compiled with PERIODIC flag!\n#\n");
#endif
	if (biggest+rmax > boxsize)
	    printf("#\n# WARNING: non-cubic input detected, but could overlap periodically!\n#\n");
    }
    return box;
}

void invert_weights(Particle *p, int np) {
    for (int j=0; j<np; j++) p[j].w *= -1.0;
    printf("# Multiplying all weights by -1\n");
}

void balance_weights(Particle *p, int np) {
    Float sumpos = 0.0, sumneg = 0.0;
    for (int j=0; j<np; j++)
	if (p[j].w>=0.0) sumpos += p[j].w;
	    else sumneg += p[j].w;
    if (sumneg==0.0 || sumpos==0.0) {
	fprintf(stderr,"Asked to rebalance weights, but there are not both positive and negative weights\n");
	abort();
    }
    Float rescale = sumpos/(-sumneg);
    printf("# Rescaling negative weights by %f\n", rescale);
    for (int j=0; j<np; j++)
	if (p[j].w<0.0) p[j].w *= rescale;
    return;
}

// ================================ main() =============================



int main(int argc, char *argv[]) {

	Parameters par=Parameters(argc,argv);

	Float3 shift;
    Particle *orig_p;
    if (par.make_random) {
	// If you want to just make random particles instead:
	assert(par.np>0);
	assert(par.boxsize==par.rescale);    // Nonsense if not!
	orig_p = make_particles(par.boxsize, par.np);
    } else {
	orig_p = read_particles(par.rescale, &par.np, par.fname, par.rstart, par.nmax);
	assert(par.np>0);
    }
    if (par.qinvert) invert_weights(orig_p, par.np);
    if (par.qbalance) balance_weights(orig_p, par.np);

    par.perbox=check_bounding_box(orig_p, par.np, par.boxsize, par.rmax,shift);

    // Now ready to compute!
    // Sort the particles into the grid.
    Grid grid(orig_p, par.np, par.boxsize, par.nside, shift);

    Float grid_density = (double)par.np/grid.nf;//pow(par.nside,3);
    printf("Average number of particles per grid cell = %6.2f\n", grid_density);
    printf("Average number of particles per max_radius ball = %6.2f\n",
    		par.np*4.0*M_PI/3.0*pow(par.rmax/par.boxsize,3.0));
    if (grid_density<1) printf("#\n# WARNING: grid appears inefficiently fine.\n#\n");


//   Count occupation of pos and neg cells
//    printf("BEG\n");
//    for(int ii=0;ii<grid.ncells;ii++){
//    	Cell fou=grid.c[ii];
//    	int sne=0;
//    	int spo=0;
//    	for (int l = fou.start; l<fou.start+fou.np; l++)
//    		grid.p[l].w>0?spo++:sne++;
//    	printf("%d %d\n",spo,sne);
//    }
//
//    return 0;

    printf("# Done gridding the particles\n");
    printf("# %d particles in use, %d with positive weight\n", grid.np, grid.np_pos);
    printf("# Weights: Positive particles sum to %f\n", grid.sumw_pos);
    printf("#          Negative particles sum to %f\n", grid.sumw_neg);
    free(orig_p); // Particles are now only stored in grid

    fflush(NULL);

    // Everything above here takes negligible time.  This line is nearly all of the work.
    compute_integral(&grid,&par);

    rusage ru;
    getrusage(RUSAGE_SELF, &ru);

    fprintf(stderr,"# user CPU time used: %ld \n# system CPU time used: %ld \n# maximum resident set size: %ld \n# integral shared memory size: %ld \n# integral unshared data size: %ld \n# integral unshared stack size: %ld \n# page reclaims (soft page faults): %ld \n# page faults (hard page faults): %ld \n# swaps: %ld \n# block input operations: %ld \n# block output operations: %ld \n# IPC messages sent: %ld \n# IPC messages received: %ld \n# signals received: %ld \n# voluntary context switches: %ld \n# involuntary context switches: %ld \n",ru.ru_utime.tv_sec,ru.ru_stime.tv_sec,ru.ru_maxrss,ru.ru_ixrss,ru.ru_idrss,ru.ru_isrss,ru.ru_minflt,ru.ru_majflt,ru.ru_nswap,ru.ru_inblock,ru.ru_oublock,ru.ru_msgsnd,ru.ru_msgrcv,ru.ru_nsignals,ru.ru_nvcsw,ru.ru_nivcsw);

    return 0;
}
