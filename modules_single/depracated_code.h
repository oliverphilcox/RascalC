
bool check_bounding_box(Particle *p, int np, Float boxsize, Float rmax, Float3& pmin) {
    // From driver.h. Now replaced by compute_bounding_box function.
    
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


	    void report_integrals() {
            // From the integrals class - now replaced with SAVE integral.
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
