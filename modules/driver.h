// driver.h - this contains various c++ functions to create particles in random positions / read them in from file. Based on code by Alex Wiegand.

#ifndef DRIVER_H
#define DRIVER_H

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
    // This will read particles from a file, space-separated x,y,z,w,JK for weight w, jackknife region JK
    // Particle positions will be rescaled by the variable 'rescale'.
    // For example, if rescale==boxsize, then inputing the unit cube will cover the periodic volume
    char line[1000];
    int j=0,n=0;
    FILE *fp;
    int stat;
    double tmp[5];
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
        stat=sscanf(line, "%lf %lf %lf %lf %lf", tmp, tmp+1, tmp+2, tmp+3, tmp+4);// %lf %lf", tmp, tmp+1, tmp+2, tmp+3, tmp+4, tmp+5, tmp+6);

        if (stat<4) {
        	fprintf(stderr,"Particle %d has bad format\n", j); // Not enough coordinates
        	abort();
        }

        p[j].pos.x = tmp[0]*rescale;
        p[j].pos.y = tmp[1]*rescale;
        p[j].pos.z = tmp[2]*rescale;

        // NB: Only works for 5 entries per line now.
        // If there are 5 or 7 entries per line get the weights from line 4 or 6
        // Otherwise fill the weights with 1 or -1 depending on the value of rstart
        // For grid_covariance rstart is typically not used
        if(stat!=5)//7&&stat!=5)
		   if(rstart>0&&j>=rstart)
			   p[j].w = -1.;
		   else
			   p[j].w = 1.;
		else{
		   if(rstart>0&&j>=rstart)
			   p[j].w = -tmp[stat-2]; //read in weights
		   else
			   p[j].w = tmp[stat-2]; 
        p[j].JK = tmp[stat-1]; // read in JK region
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

#endif
