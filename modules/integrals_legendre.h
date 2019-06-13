// Rewritten integrals.h code for legendre multipoles
#include "parameters.h"
#include "correlation_function.h"
#include "cell_utilities.h"
#include "legendre_utilities.h"

#ifndef INTEGRALS_LEGENDRE_H
#define INTEGRALS_LEGENDRE_H

class Integrals{
public:
    CorrelationFunction *cf12, *cf13, *cf24;
    
private:
    int nbin, mbin, max_l;
    Float rmin,rmax,mumin,mumax; //Ranges in r and mu
    Float *r_high, *r_low; // Max and min of each radial bin
    Float *c2, *c3, *c4; // Arrays to accumulate integrals
    char* out_file;
    bool box; // Flags to decide whether we have a periodic box
    int I1, I2, I3, I4; // indices for which fields to use for each particle
    
    uint64 *binct, *binct3, *binct4; // Arrays to accumulate bin counts
    SurveyCorrection *sc12,*sc23,*sc34; // survey correction function
    
public:
    Integrals(){};

    Integrals(Parameters *par, CorrelationFunction *_cf12, CorrelationFunction *_cf13, CorrelationFunction *_cf24, int _I1, int _I2, int _I3, int _I4, SurveyCorrection *_sc12, SurveyCorrection *_sc23, SurveyCorrection *_sc34){
        sc12 = _sc12;
        sc23 = _sc23;
        sc34 = _sc34;
        max_l = par->max_l;
        cf12 = new CorrelationFunction(_cf12);
        cf13 = new CorrelationFunction(_cf13);
        cf24 = new CorrelationFunction(_cf24);
        I1=_I1;
        I2=_I2;
        I3=_I3;
        I4=_I4;
        init(par);
    }
    
    void init(Parameters *par){
        nbin = par->nbin; // number of radial bins
        mbin = par->mbin; // number of mu bins
        out_file = par->out_file; // output directory
        max_l = par->max_l; // maximum Legendre polynomial
        
        int ec=0;
        // Initialize the binning
        ec+=posix_memalign((void **) &c2, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);
        ec+=posix_memalign((void **) &c3, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);
        ec+=posix_memalign((void **) &c4, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);

        ec+=posix_memalign((void **) &binct, PAGE, sizeof(uint64)*nbin*mbin*nbin*mbin);
        ec+=posix_memalign((void **) &binct3, PAGE, sizeof(uint64)*nbin*mbin*nbin*mbin);
        ec+=posix_memalign((void **) &binct4, PAGE, sizeof(uint64)*nbin*mbin*nbin*mbin);
        
        assert(ec==0);
        reset();
        
        box=par->perbox;

        rmax=par->rmax;
        rmin=par->rmin;

        mumax=par->mumax;
        mumin=par->mumin;
        
        r_high = par->radial_bins_high;
        r_low = par->radial_bins_low;
        
    }

    ~Integrals() {
        free(c2);
        free(c3);
        free(c4);
        free(binct);
        free(binct3);
        free(binct4);
    }

    void reset(){
        for (int j=0; j<nbin*mbin*nbin*mbin; j++) {
            c2[j]=0;
            c3[j]=0;
            c4[j]=0;
            binct[j] = 0;
            binct3[j] = 0;
            binct4[j] = 0;
        }
    }

    inline int get_radial_bin(Float r){
        // Computes radial bin
        
        int which_bin = -1; // default if outside bins
        for(int i=0;i<nbin;i++){
            if((r>r_low[i])&&(r<r_high[i])){
                which_bin=i;
                break;
            }
            if((i==nbin-1)&&(r>r_high[i])){
                which_bin=nbin; // if above top bin
            }
        }
        return which_bin;
        }
    
    inline void second(const Particle* pi_list, const int* prim_ids, int pln, const Particle pj, const int pj_id, int* &bin, Float* &wij, const double prob, Float* &factor_ij, Float* &poly_ij){
        // Accumulates the two point integral C2. Also outputs an array of bin values for later reuse.
        // Prob. here is defined as g_ij / f_ij where g_ij is the sampling PDF and f_ij is the true data PDF for picking pairs (equal to n_i/N n_j/N for N particles)
        // Prob1/2 are for when we divide the random particles into two subsets 1 and 2.
        
        Float tmp_weight, tmp_xi, rij_mag, rij_mu, c2v;
        Particle pi;
        int tmp_bin;
        Float correction_factor;
        int max_bin=nbin;
        Float polynomials[mbin];

        for(int i=0;i<pln;i++){ // Iterate over particle in pi_list
                if(prim_ids[i]==pj_id){
                    wij[i]=-1;
                    continue; // don't self-count
                }
                pi = pi_list[i]; // first particle
                cleanup_l(pi.pos,pj.pos,rij_mag,rij_mu); // define |r_ij| and ang(r_ij)
                
                tmp_bin = get_radial_bin(rij_mag); // radial bin
                if ((tmp_bin<0)||(tmp_bin>=max_bin)){
                    wij[i]=-1;
                    continue;
                }
                 
                tmp_weight = pi.w*pj.w; // product of weights
                tmp_xi = cf12->xi(rij_mag, rij_mu); // correlation function for i-j
                
                // Save into arrays for later
                bin[i] = tmp_bin;
                wij[i] = tmp_weight;
                
                
                // Now compute the integral:
                c2v = tmp_weight*tmp_weight*(1.+tmp_xi) / prob*2.; // c2 contribution with symmetry factor
                correction_factor = sc12->correction_function(tmp_bin,rij_mu);
                factor_ij[i] = correction_factor;
                
                // load all legendre polynomials
                legendre_polynomials(rij_mu, max_l,polynomials);
                
                // Now add to relevant bins
                int out_bin,tmp_out_bin;
                for(int p_bin=0;p_bin<mbin;p_bin++){ // iterate over all legendre polynomials
                    // add to output list
                    poly_ij[i*mbin+p_bin]=polynomials[p_bin];
                    tmp_out_bin = (tmp_bin*mbin+p_bin)*nbin*mbin+tmp_bin*mbin;
                    for(int q_bin=0;q_bin<mbin;q_bin++){ // second polynomial
                        out_bin = tmp_out_bin+q_bin; // output bin

                        // Now add to integral with correct kernel
                        c2[out_bin]+=c2v*pow(correction_factor,2)*polynomials[p_bin]*polynomials[q_bin];
                        
                        binct[out_bin]++; // only count actual contributions to bin
                    }
                }
        }
    }

    inline void third(const Particle* pi_list, const int* prim_ids, const int pln, const Particle pj, const Particle pk, const int pj_id, const int pk_id, const int* bin_ij, const Float* wij, Float* &xi_ik, Float* wijk, const double prob, const Float* factor_ij, const Float* poly_ij){
        // Accumulates the three point integral C3. Also outputs an array of xi_ik and bin_ik values for later reuse.
        
        // First define variables:
        Particle pi;
        Float rik_mag, rik_mu, c3v, rjk_mag, rjk_mu, tmp_weight, xi_ik_tmp;
        int tmp_bin, tmp_full_bin;
        
        cleanup_l(pj.pos,pk.pos,rjk_mag,rjk_mu);         
        tmp_bin = get_radial_bin(rjk_mag); // radial bin
        int max_bin = nbin,out_bin;
        Float correction_factors,polynomials_jk[mbin];
        
        // load all legendre polynomials
        legendre_polynomials(rjk_mu, max_l, polynomials_jk);            
        
        for(int i=0;i<pln;i++){ // Iterate over particle in pi_list
            if((pk_id==pj_id)||(wij[i]==-1)||(prim_ids[i]==pk_id)){
                wijk[i]=-1;
                continue; // skip incorrect bins / ij,jk self counts
            }
            pi = pi_list[i];
            cleanup_l(pi.pos,pk.pos,rik_mag,rik_mu); // define angles/lengths
            xi_ik_tmp = cf13->xi(rik_mag, rik_mu); 
            
            tmp_weight = wij[i]*pk.w; // product of weights, w_iw_jw_k 
            
            // save arrays for later
            xi_ik[i]=xi_ik_tmp;
            wijk[i]=tmp_weight;
            
            if ((tmp_bin<0)||(tmp_bin>=max_bin)){
                // Don't add contributions of this to C3 - but STILL save xi_ik etc. for later
                continue; // if not in correct bin
            }
            
            // Now compute the integral;
            c3v = tmp_weight*pj.w/prob*xi_ik_tmp*4.; // include symmetry factor
            correction_factors = factor_ij[i]*sc23->correction_function(tmp_bin,rjk_mu);
                
            // Now add to relevant bins
            for(int p_bin=0;p_bin<mbin;p_bin++){ // iterate over all legendre polynomials
                tmp_full_bin = (bin_ij[i]*mbin+p_bin)*nbin*mbin+tmp_bin*mbin;
                for(int q_bin=0;q_bin<mbin;q_bin++){ // second polynomial
                    out_bin = tmp_full_bin+q_bin; // output bin

                    // Now add to integral with correct kernel
                    c3[out_bin]+=c3v*correction_factors*poly_ij[p_bin]*polynomials_jk[q_bin];
                        
                    binct3[out_bin]++; // only count actual contributions to bin
                }
            }
        }
    }

    inline void fourth(const Particle* pi_list, const int* prim_ids, const int pln, const Particle pj, const Particle pk, const Particle pl, const int pj_id, const int pk_id, const int pl_id, const int* bin_ij, const Float* wijk, const Float* xi_ik, const double prob, const Float* factor_ij, const Float* poly_ij){
        // Accumulates the four point integral C4.
        
        // First define variables
        Particle pi;
        Float rjl_mag, rjl_mu, rkl_mag, rkl_mu, c4v, xi_jl, tmp_weight;
        int tmp_bin, tmp_full_bin;
        cleanup_l(pl.pos,pk.pos,rkl_mag,rkl_mu); 
        
        int max_bin = nbin, out_bin;
        Float correction_factors, polynomials_kl[mbin];
        tmp_bin = get_radial_bin(rkl_mag); // radial kl bin
        
        if ((tmp_bin<0)||(tmp_bin>=max_bin)) return; // if not in correct bin                    
        cleanup_l(pl.pos,pj.pos,rjl_mag,rjl_mu); 
        xi_jl = cf24->xi(rjl_mag, rjl_mu); // j-l correlation
        
        // load all legendre polynomials
        legendre_polynomials(rkl_mu, max_l, polynomials_kl);
        
        for(int i=0;i<pln;i++){ // Iterate over particle in pi_list
            if(wijk[i]==-1) continue; // skip incorrect bins / ij self counts
            if(prim_ids[i]==pl_id) continue; // don't self-count
            
            pi = pi_list[i];
            tmp_weight = wijk[i]*pl.w; // product of weights, w_i*w_j*w_k*w_l
            
            // Now compute the integral;
            c4v = tmp_weight/prob*2.*xi_ik[i]*xi_jl; // with xi_ik*xi_jl = xi_il*xi_jk symmetry factor
            correction_factors = factor_ij[i]*sc34->correction_function(tmp_bin,rkl_mu);
            
            // Now add to relevant bins
            for(int p_bin=0;p_bin<mbin;p_bin++){ // iterate over all legendre polynomials
                tmp_full_bin = (bin_ij[i]*mbin+p_bin)*nbin*mbin+tmp_bin*mbin;
                for(int q_bin=0;q_bin<mbin;q_bin++){ // second polynomial
                    out_bin = tmp_full_bin+q_bin; // output bin

                    // Now add to integral with correct kernel
                    c4[out_bin]+=c4v*correction_factors*poly_ij[p_bin]*polynomials_kl[q_bin];
                    binct4[out_bin]++; // only count actual contributions to bin
                }
            }            
        }
    }
        
public:
    void cleanup_l(Float3 p1,Float3 p2,Float& norm,Float& mu){
        Float3 pos=p1-p2;
        norm = pos.norm();
#ifndef PERIODIC
        Float3 los=p1+p2; // No 1/2 as normalized anyway below
        mu = fabs(pos.dot(los)/norm/los.norm());
#else
        // In the periodic case use z-direction for mu
        mu = fabs(pos.z/norm);
#endif
    }
public:
    void sum_ints(Integrals* ints) {
        // Add the values accumulated in ints to the corresponding internal sums
        for(int i=0;i<nbin*mbin*nbin*mbin;i++){
            c2[i]+=ints->c2[i];
            c3[i]+=ints->c3[i];
            c4[i]+=ints->c4[i];
            binct[i]+=ints->binct[i];
            binct3[i]+=ints->binct3[i];
            binct4[i]+=ints->binct4[i];
        }
    }
    
    void frobenius_difference_sum(Integrals* ints, int n_loop, Float &frobC2, Float &frobC3, Float &frobC4){
        // Add the values accumulated in ints to the corresponding internal sums and compute the Frobenius norm difference between integrals
        Float n_loops = (Float)n_loop;
        Float self_c2=0, diff_c2=0;
        Float self_c3=0, diff_c3=0;
        Float self_c4=0, diff_c4=0;
        
        // Compute Frobenius norms and sum integrals
        for(int i=0;i<nbin*mbin;i++){
            for(int j=0;j<nbin*mbin;j++){
                self_c2+=pow(c2[i*nbin*mbin+j]/n_loops,2.);
                diff_c2+=pow(c2[i*nbin*mbin+j]/n_loops-(c2[i*nbin*mbin+j]+ints->c2[i*nbin*mbin+j])/(n_loops+1.),2.);            
                self_c3+=pow(c3[i*nbin*mbin+j]/n_loops,2.);
                diff_c3+=pow(c3[i*nbin*mbin+j]/n_loops-(c3[i*nbin*mbin+j]+ints->c3[i*nbin*mbin+j])/(n_loops+1.),2.);
                self_c4+=pow(c4[i*nbin*mbin+j]/n_loops,2.);
                diff_c4+=pow(c4[i*nbin*mbin+j]/n_loops-(c4[i*nbin*mbin+j]+ints->c4[i*nbin*mbin+j])/(n_loops+1.),2.);
                }
            }
        
        self_c2=sqrt(self_c2);
        diff_c2=sqrt(diff_c2);
        diff_c3=sqrt(diff_c3);
        diff_c4=sqrt(diff_c4);
        self_c3=sqrt(self_c3);
        self_c4=sqrt(self_c4);
        
        // Return percent difference
        frobC2=100.*(diff_c2/self_c2);
        frobC3=100.*(diff_c3/self_c3);
        frobC4=100.*(diff_c4/self_c4);
        }

    void sum_total_counts(uint64& acc2, uint64& acc3, uint64& acc4){
        // Add local counts to total bin counts in acc2-4
        for (int i=0; i<nbin*mbin*nbin*mbin; i++) {
            acc2+=binct[i];
            acc3+=binct3[i];
            acc4+=binct4[i];
        }
    }
    
    void normalize(Float norm1, Float norm2, Float norm3, Float norm4, Float n_pairs, Float n_triples, Float n_quads){
        // Normalize the accumulated integrals (partly done by the normalising probabilities used from the selected cubes)
        // n_pair etc. are the number of PARTICLE pairs etc. attempted (not including rejected cells, but including pairs which don't fall in correct bin ranges)
        // To avoid recomputation
        double corrf2 = norm1*norm2; // correction factor for densities of random points
        double corrf3 = corrf2*norm3;
        double corrf4 = corrf3*norm4;;
        
        for(int i = 0; i<nbin*mbin*nbin*mbin;i++){
            c2[i]/=(n_pairs*corrf2);
            c3[i]/=(n_triples*corrf3);
            c4[i]/=(n_quads*corrf4);
        }
        
        // Further normalize by pre-factor
        int legendre_p,legendre_q;
        Float v_a,v_b,normalization;
        
        for(int i=0; i<nbin*mbin;i++){
            legendre_p = (i%mbin)*2.; // first legendre index
            v_a = 4.*M_PI/3.*(pow(r_high[i/mbin],3)-pow(r_low[i/mbin],3));
            
            for(int j=0;j<nbin*mbin;j++){
                legendre_q = (j%mbin)*2.; // second legendre index
                v_b = 4.*M_PI/3.*(pow(r_high[j/mbin],3)-pow(r_low[j/mbin],3));
            
                normalization = float((2*legendre_p+1)*(2*legendre_q+1))/(v_a*v_b);                
                c2[i*nbin*mbin+j]*=normalization;
                c3[i*nbin*mbin+j]*=normalization;
                c4[i*nbin*mbin+j]*=normalization;
            }
        }
    }
    
    void save_counts(uint64 pair_counts,uint64 triple_counts,uint64 quad_counts){
        // Print the counts for each integral (used for combining the estimates outside of C++)
        // This is the number of counts used in each loop [always the same]
        char counts_file[1000];
        snprintf(counts_file, sizeof counts_file, "%sCovMatricesAll/total_counts_n%d_l%d_%d%d,%d%d.txt",out_file,nbin,max_l,I1,I2,I3,I4);
        FILE * CountsFile = fopen(counts_file,"w");
        fprintf(CountsFile,"%llu\n",pair_counts);
        fprintf(CountsFile,"%llu\n",triple_counts);
        fprintf(CountsFile,"%llu\n",quad_counts);
        
        fflush(NULL);
        fclose(CountsFile);
    }
        
    
    void save_integrals(char* suffix, bool save_all) {
    /* Print integral outputs to file. 
        * In txt files {c2,c3,c4}_leg_n{nbin}_m{mbin}.txt there are lists of the outputs of c2,c3,c4 that are already normalized and multiplied by combinatoric factors. The n and m strings specify the number of n and m bins present.
        */
        // Create output files
        
        char c2name[1000];
        snprintf(c2name, sizeof c2name, "%sCovMatricesAll/c2_n%d_l%d_%d%d_%s.txt", out_file,nbin, max_l,I1,I2,suffix);
        char c3name[1000];
        snprintf(c3name, sizeof c3name, "%sCovMatricesAll/c3_n%d_l%d_%d,%d%d_%s.txt", out_file, nbin, max_l,I2,I1,I3,suffix);
        char c4name[1000];
        snprintf(c4name, sizeof c4name, "%sCovMatricesAll/c4_n%d_l%d_%d%d,%d%d_%s.txt", out_file, nbin, max_l, I1,I2,I3,I4,suffix);
        FILE * C2File = fopen(c2name,"w"); // for c2 part of integral
        FILE * C3File = fopen(c3name,"w"); // for c3 part of integral
        FILE * C4File = fopen(c4name,"w"); // for c4 part of integral
        
        for(int i=0;i<nbin*mbin;i++){
            for(int j=0;j<nbin*mbin;j++){
                fprintf(C2File,"%le\t",c2[i*nbin*mbin+j]);
                fprintf(C3File,"%le\t",c3[i*nbin*mbin+j]);
                fprintf(C4File,"%le\t",c4[i*nbin*mbin+j]);
            }
            fprintf(C2File,"\n");
            fprintf(C3File,"\n"); // new line each end of row
            fprintf(C4File,"\n");
        }
        
        fflush(NULL);
        
        // Close open files
        fclose(C2File);
        fclose(C3File);
        fclose(C4File);
        
        if(save_all==1){
            char binname4[1000];
            snprintf(binname4,sizeof binname4, "%sCovMatricesAll/binct_c4_n%d_l%d_%d%d,%d%d_%s.txt",out_file, nbin,max_l,I1,I2,I3,I4,suffix);
            FILE * BinFile4 = fopen(binname4,"w");
            
            char binname3[1000];
            snprintf(binname3,sizeof binname3, "%sCovMatricesAll/binct_c3_n%d_l%d_%d,%d%d_%s.txt",out_file, nbin,max_l,I2,I1,I3,suffix);
            FILE * BinFile3 = fopen(binname3,"w");
            
            char binname2[1000];
            snprintf(binname2,sizeof binname2, "%sCovMatricesAll/binct_c2_n%d_l%d_%d%d_%s.txt",out_file, nbin,max_l,I1,I2,suffix);
            FILE * BinFile2 = fopen(binname2,"w");
        
            for(int i=0;i<nbin*mbin;i++){
                for(int j=0;j<nbin*mbin;j++){
                    fprintf(BinFile2,"%llu\n",binct[i*nbin*mbin+j]);
                    fprintf(BinFile3,"%llu\t",binct3[i*nbin*mbin+j]);
                    fprintf(BinFile4,"%llu\t",binct4[i*nbin*mbin+j]);
                    }
                fprintf(BinFile2,"\n");
                fprintf(BinFile3,"\n");
                fprintf(BinFile4,"\n");
            }
        
            fclose(BinFile2);
            fclose(BinFile3);
            fclose(BinFile4);
        }
    }

};

#endif
