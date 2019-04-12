// driver.h - this contains various c++ functions to create particles in random positions / read them in from file. Based on code by Alex Wiegand.
#include "cell_utilities.h"
#ifndef LEGENDRE
#ifndef POWER
    #include "jackknife_weights.h"
#endif
#endif

#ifndef DRIVER_H
#define DRIVER_H

// ====================  The Driver ===========================

Particle *make_particles(Float3 rect_boxsize, int np) {
    // Make np random particles
    srand48(1);      // For reproducibility
    Particle *p = (Particle *)malloc(sizeof(Particle)*np);
    for (int j=0; j<np; j++) {
        p[j].pos.x = drand48()*rect_boxsize.x;
        p[j].pos.y = drand48()*rect_boxsize.y;
        p[j].pos.z = drand48()*rect_boxsize.z;
        p[j].w = 1.0;
        p[j].JK = 0.;
        p[j].rand_class = rand()%2; // assign random class to each particle
    }
    fprintf(stderr,"# WARNING - Jackknife regions are not yet implemented for particles made at random. All particles are assigned to the same region.\n");
    printf("# Done making %d random particles, periodically distributed.\n", np);
    return p;
}

#ifdef JACKKNIFE
Particle *read_particles(Float rescale, int *np, const char *filename, const int rstart, uint64 nmax, const JK_weights *JK) {
#else
Particle *read_particles(Float rescale, int *np, const char *filename, const int rstart, uint64 nmax) {
#endif
    // This will read particles from a file, space-separated x,y,z,w,JK for weight w, (jackknife region JK)
    // Particle positions will be rescaled by the variable 'rescale'.
    // For example, if rescale==boxsize, then inputting the unit cube will cover the periodic volume
    char line[1000];
    int j=0,n=0;
    FILE *fp;
    int stat;
    double tmp[5];

    fp = fopen(filename, "r");
    if (fp==NULL) {
        fprintf(stderr,"File %s not found\n", filename); abort();
    }

#ifdef JACKKNIFE
    // Store filled jackknives in local memory to avoid file corruption
    int tmp_n_JK = JK->n_JK_filled;
    int tmp_filled_JK[tmp_n_JK];
    for(int ii=0;ii<tmp_n_JK;ii++) tmp_filled_JK[ii]=JK->filled_JKs[ii];
#endif
    
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
        stat=sscanf(line, "%lf %lf %lf %lf %lf", tmp, tmp+1, tmp+2, tmp+3, tmp+4);

        if (stat<4) {
        	fprintf(stderr,"Particle %d has bad format\n", j); // Not enough coordinates
        	abort();
        }

        p[j].pos.x = tmp[0]*rescale;
        p[j].pos.y = tmp[1]*rescale;
        p[j].pos.z = tmp[2]*rescale;
        p[j].rand_class = rand()%2;
        
        // Get the weights from line 4 if present, else fill with +1/-1 depending on the value of rstart
        // For grid_covariance rstart is typically not used
        if(stat!=5)
		   if(rstart>0&&j>=rstart)
			   p[j].w = -1.;
		   else
			   p[j].w = 1.;
		else{
		   if(rstart>0&&j>=rstart)
			   p[j].w = -tmp[stat-2]; //read in weights
		   else
			   p[j].w = tmp[stat-2]; 
#ifdef JACKKNIFE
        int tmp_JK = tmp[stat-1]; // read in JK region
		
		// Collapse jacknife indices to only include filled JKs:
		p[j].JK=-1;
        
        for (int x=0;x<tmp_n_JK;x++){
            if (tmp_filled_JK[x]==tmp_JK) p[j].JK=x;
        }
        assert(p[j].JK!=-1); // ensure we find jackknife index		    
#endif
        }
		
		j++;
    }
    fclose(fp);
    printf("# Done reading the particles\n");
    
    return p;
}

bool compute_bounding_box(Particle *p, int np, Float3 &rect_boxsize, Float rmax, Float3& pmin, int nside) {
    // Compute the boxsize of the bounding cuboid box and determine whether we are periodic
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
    if (prange.x>0.99*biggest && prange.y>0.99*biggest && prange.z>0.99*biggest) {
        // Probably using a cube of inputs, intended for a periodic box
    	box=true;
#ifndef PERIODIC
    	fprintf(stderr,"#\n# WARNING: cubic input detected but you have not compiled with PERIODIC flag!\n#\n");
    	printf("#\n# WARNING: cubic input detected but you have not compiled with PERIODIC flag!\n#\n");
#endif
        // Set boxsize to be the biggest dimension which allows for periodic overlap
        rect_boxsize= {biggest,biggest,biggest};
        printf("# Setting periodic box-size to %6.2f\n", biggest);
    } else {
        // Probably a non-periodic input (e.g. a real dataset)
        box=false;
#ifdef PERIODIC
    	fprintf(stderr,"#\n# WARNING: non-cubic input detected but you have compiled with PERIODIC flag!\n#\n");
    	printf("#\n# WARNING: non-cubic input detected but you have compiled with PERIODIC flag!\n#\n");
#endif
        // set max_boxsize to just enclose the biggest dimension plus r_max 
        // NB: We natively wrap the grid (to allow for any position of the center of the grid)
        // Must add rmax to biggest to ensure there is no periodic overlap in this case.
        Float max_boxsize = 1.05*(biggest+rmax);
        Float cellsize = max_boxsize/nside; // compute the width of each cell
        // Now compute the size of the box in every dimension
        rect_boxsize = ceil3(prange/cellsize)*cellsize; // to ensure we fit an integer number of cells in each direction
        printf("# Setting non-periodic box-size to {%6.2f,%6.2f,%6.2f}\n", rect_boxsize.x,rect_boxsize.y,rect_boxsize.z);
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

#endif
