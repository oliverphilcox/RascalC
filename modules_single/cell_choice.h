
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
