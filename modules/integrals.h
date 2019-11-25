// Rewritten integrals.h code for grid_covariance.cpp (originally from Alex Wiegand) to parallelize and compute integrands to a given quad of particles
#include "parameters.h"
#include "correlation_function.h"
#include "cell_utilities.h"
#include "jackknife_weights.h"
#include "legendre_utilities.h"

#ifndef INTEGRALS_H
#define INTEGRALS_H

class Integrals{
private:
    CorrelationFunction *cf12, *cf13, *cf24;
    int nbin, mbin;
    Float rmin,rmax,mumin,mumax,dmu; //Ranges in r and mu
    Float *r_high, *r_low; // Max and min of each radial bin
    Float *Ra, *c2, *c3, *c4; // Arrays to accumulate integrals
    JK_weights *JK12, *JK23, *JK34; // RR counts and jackknife weightts
#ifdef JACKKNIFE
    int n_jack;
    Float *c2j, *c3j, *c4j; // Arrays to accumulate jackknife integrals
    Float *EEaA1, *EEaA2; // Array to accumulate the two-independent xi-weighted pair counts
    Float *RRaA1, *RRaA2; // Array to accumulate the two-independent pair count estimates
    Float *product_weights12_12, *product_weights12_23, *product_weights12_34; // arrays to get products of jackknife weights to avoid recomputation
#endif
    char* out_file;
    bool box,rad=0; // Flags to decide whether we have a periodic box + if we have a radial correlation function only
    int I1, I2, I3, I4; // indices for which fields to use for each particle

    uint64 *binct, *binct3, *binct4; // Arrays to accumulate bin counts

public:
    Integrals(){};
#ifdef JACKKNIFE
    Integrals(Parameters *par, CorrelationFunction *_cf12, CorrelationFunction *_cf13, CorrelationFunction *_cf24, JK_weights *_JK12, JK_weights *_JK23, JK_weights *_JK34, int _I1, int _I2, int _I3, int _I4, Float* _product_weights12_12, Float* _product_weights12_23, Float* _product_weights12_34){
        product_weights12_12=_product_weights12_12;
        product_weights12_23=_product_weights12_23;
        product_weights12_34=_product_weights12_34;
#else
    Integrals(Parameters *par, CorrelationFunction *_cf12, CorrelationFunction *_cf13, CorrelationFunction *_cf24, JK_weights *_JK12, JK_weights *_JK23, JK_weights *_JK34, int _I1, int _I2, int _I3, int _I4){
#endif
        cf12 = new CorrelationFunction(_cf12);
        cf13 = new CorrelationFunction(_cf13);
        cf24 = new CorrelationFunction(_cf24);
        JK12 = _JK12;
        JK23 = _JK23;
        JK34 = _JK34;
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

        int ec=0;
        // Initialize the binning
        ec+=posix_memalign((void **) &Ra, PAGE, sizeof(double)*nbin*mbin);
        ec+=posix_memalign((void **) &c2, PAGE, sizeof(double)*nbin*mbin);
        ec+=posix_memalign((void **) &c3, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);
        ec+=posix_memalign((void **) &c4, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);
        ec+=posix_memalign((void **) &binct, PAGE, sizeof(uint64)*nbin*mbin);
        ec+=posix_memalign((void **) &binct3, PAGE, sizeof(uint64)*nbin*mbin*nbin*mbin);
        ec+=posix_memalign((void **) &binct4, PAGE, sizeof(uint64)*nbin*mbin*nbin*mbin);
#ifdef JACKKNIFE
        n_jack = fmax(fmax(JK12->n_JK_filled,JK23->n_JK_filled),JK34->n_JK_filled); // number of non-empty jackknives
        ec+=posix_memalign((void **) &c2j, PAGE, sizeof(double)*nbin*mbin);
        ec+=posix_memalign((void **) &c3j, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);
        ec+=posix_memalign((void **) &c4j, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);

        ec+=posix_memalign((void **) &EEaA1, PAGE, sizeof(double)*nbin*mbin*n_jack);
        ec+=posix_memalign((void **) &EEaA2, PAGE, sizeof(double)*nbin*mbin*n_jack);
        ec+=posix_memalign((void **) &RRaA1, PAGE, sizeof(double)*nbin*mbin*n_jack);
        ec+=posix_memalign((void **) &RRaA2, PAGE, sizeof(double)*nbin*mbin*n_jack);
#endif

        assert(ec==0);
        reset();

        box=par->perbox;

        rmax=par->rmax;
        rmin=par->rmin;

        mumax=par->mumax;
        mumin=par->mumin;

        r_high = par->radial_bins_high;
        r_low = par->radial_bins_low;

        dmu=(mumax-mumin)/mbin;

        rad=mbin==1&&dmu==1.;
    }

    ~Integrals() {
        free(Ra);
        free(c2);
        free(c3);
        free(c4);
        free(binct);
        free(binct3);
        free(binct4);
#ifdef JACKKNIFE
        free(c2j);
        free(c3j);
        free(c4j);
        free(EEaA1);
        free(EEaA2);
        free(RRaA1);
        free(RRaA2);
#endif
    }

    void reset(){
        for (int j=0; j<nbin*mbin; j++) {
            Ra[j] = 0;
            c2[j] = 0;
#ifdef JACKKNIFE
            c2j[j] = 0;
#endif
            binct[j] = 0;
        }
        for (int j=0; j<nbin*mbin*nbin*mbin; j++) {
            c3[j]=0;
            c4[j]=0;
            binct3[j] = 0;
            binct4[j] = 0;
#ifdef JACKKNIFE
            c3j[j]=0;
            c4j[j]=0;
        }
        for (int j=0;j<nbin*mbin*n_jack;j++){
            EEaA1[j] = 0;
            EEaA2[j] = 0;
            RRaA1[j] = 0;
            RRaA2[j] = 0;
#endif
        }
    }

    inline int getbin(Float r, Float mu){
        // Linearizes 2D indices
        // First define which r bin we are in;
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
        return which_bin*mbin + floor((mu-mumin)/dmu);
    }

    inline void second(const Particle* pi_list, const int* prim_ids, int pln, const Particle pj, const int pj_id, int* &bin, Float* &wij, const double prob, const double prob1, const double prob2){
        // Accumulates the two point integral C2. Also outputs an array of bin values for later reuse.
        // Prob. here is defined as g_ij / f_ij where g_ij is the sampling PDF and f_ij is the true data PDF for picking pairs (equal to n_i/N n_j/N for N particles)
        // Prob1/2 are for when we divide the random particles into two subsets 1 and 2.
        Float tmp_weight, tmp_xi, rij_mag, rij_mu, c2v,rav;
        Particle pi;
        int tmp_bin;
#ifdef JACKKNIFE
        Float c2vj,JK_weight;
        int jk_bin_i, jk_bin_j;
#endif
        for(int i=0;i<pln;i++){ // Iterate over particle in pi_list
                if((prim_ids[i]==pj_id)&&(I1==I2)){
                    wij[i]=-1;
                    continue; // don't self-count
                }

                pi = pi_list[i]; // first particle
                cleanup_l(pi.pos,pj.pos,rij_mag,rij_mu); // define |r_ij| and ang(r_ij)
                tmp_bin = getbin(rij_mag,rij_mu); // define i-j bin

                if ((tmp_bin<0)||(tmp_bin>=nbin*mbin)){
                    wij[i]=-1;
                    continue;
                }

                tmp_weight = pi.w*pj.w; // product of weights
                tmp_xi = cf12->xi(rij_mag, rij_mu); // correlation function for i-j

                // Save into arrays for later
                bin[i]=tmp_bin;
                wij[i] = tmp_weight;

                // Now compute the integral:
                c2v = tmp_weight*tmp_weight*(1.+tmp_xi) / prob*2.; // c2 contribution with symmetry factor
                rav = tmp_weight / prob; // RR_a contribution
                // Add to local integral counts:
                Ra[tmp_bin]+=rav;
                c2[tmp_bin]+=c2v;
                binct[tmp_bin]++; // only count actual contributions to bin
#ifdef JACKKNIFE
                // Compute jackknife weight tensor:
                JK_weight=weight_tensor(int(pi.JK),int(pj.JK),int(pi.JK),int(pj.JK),tmp_bin,tmp_bin,JK12,JK12, product_weights12_12);
                c2vj = c2v*JK_weight;
                c2j[tmp_bin]+=c2vj;
                // Now add EEaA bin counts:
                // If both in random set-0
                if ((pi.rand_class==0)&&(pj.rand_class==0)){
                    jk_bin_i = int(pi.JK)*nbin*mbin+tmp_bin; // EEaA bin for i particle
                    jk_bin_j = int(pj.JK)*nbin*mbin+tmp_bin; // EEaA bin for j particle
                    EEaA1[jk_bin_i]+=tmp_weight/prob1*tmp_xi/2.; // add half contribution to each jackknife
                    EEaA1[jk_bin_j]+=tmp_weight/prob1*tmp_xi/2.;
                    RRaA1[jk_bin_i]+=tmp_weight/prob1/2;
                    RRaA1[jk_bin_j]+=tmp_weight/prob1/2;
                }
                // If both in random set-1
                if((pi.rand_class==1)&&(pj.rand_class==1)){
                    jk_bin_i = int(pi.JK)*nbin*mbin+tmp_bin; // EEaA bin for i particle
                    jk_bin_j = int(pj.JK)*nbin*mbin+tmp_bin; // EEaA bin for j particle
                    EEaA2[jk_bin_i]+=tmp_weight/prob2*tmp_xi/2; // add half contribution to each jackknife
                    EEaA2[jk_bin_j]+=tmp_weight/prob2*tmp_xi/2;
                    RRaA2[jk_bin_i]+=tmp_weight/prob2/2;
                    RRaA2[jk_bin_j]+=tmp_weight/prob2/2;
                }
#endif
        }
    }
    inline void third(const Particle* pi_list, const int* prim_ids, const int pln, const Particle pj, const Particle pk, const int pj_id, const int pk_id, const int* bin_ij, const Float* wij, Float* &xi_ik, Float* wijk, const double prob){
        // Accumulates the three point integral C3. Also outputs an array of xi_ik and bin_ik values for later reuse.
        // First define variables:
        Particle pi;
        Float rik_mag, rik_mu, c3v, rjk_mag, rjk_mu, tmp_weight, xi_ik_tmp;
        int tmp_bin, tmp_full_bin, max_bin = nbin*mbin;
#ifdef JACKKNIFE
        Float c3vj, JK_weight;
#endif
        cleanup_l(pj.pos,pk.pos,rjk_mag,rjk_mu);
        tmp_bin = getbin(rjk_mag,rjk_mu); // bin for each particle

        for(int i=0;i<pln;i++){ // Iterate over particle in pi_list
            if(((pk_id==pj_id)&&(I2==I3))||((prim_ids[i]==pk_id)&&(I1==I3))||(wij[i]==-1)){
              wijk[i]=-1;
              continue; // skip incorrect bins / ij,jk self counts
            }
            pi = pi_list[i];
            cleanup_l(pi.pos,pk.pos,rik_mag,rik_mu); // define angles/lengths
            xi_ik_tmp = cf13->xi(rik_mag, rik_mu);

            tmp_weight = wij[i]*pk.w; // product of weights, w_iw_jw_k

            if(rik_mag<1e-4){
              printf("Particle separation of %.2e Mpc/h found between random particle files %d and %d. This is unusually small and will cause errors.\n",rik_mag,I2,I3);
              printf("Are the random particle files independent? The code will now exit.");
              exit(1);
            }

            // save arrays for later
            xi_ik[i]=xi_ik_tmp;
            wijk[i]=tmp_weight;
            if ((tmp_bin<0)||(tmp_bin>=max_bin)){
                // Don't add contributions of this to C3 - but STILL save xi_ik etc. for later
                continue; // if not in correct bin
            }
            // Now compute the integral;
            c3v = tmp_weight*pj.w/prob*xi_ik_tmp*4.; // include symmetry factor
            tmp_full_bin = bin_ij[i]*mbin*nbin+tmp_bin;
            // Add to local counts
            c3[tmp_full_bin]+=c3v;
            binct3[tmp_full_bin]++;
#ifdef JACKKNIFE
            // Compute jackknife weight tensor:
            JK_weight=weight_tensor(int(pi.JK),int(pj.JK),int(pj.JK),int(pk.JK),bin_ij[i],tmp_bin, JK12, JK23, product_weights12_23);
            c3vj = c3v*JK_weight;
            c3j[tmp_full_bin]+=c3vj;
#endif
        }
    }
    inline void fourth(const Particle* pi_list, const int* prim_ids, const int pln, const Particle pj, const Particle pk, const Particle pl, const int pj_id, const int pk_id, const int pl_id, const int* bin_ij, const Float* wijk, const Float* xi_ik, const double prob){
        // Accumulates the four point integral C4.
        // First define variables
        Particle pi;
        Float rjl_mag, rjl_mu, rkl_mag, rkl_mu, c4v, xi_jl, tmp_weight;
        int tmp_bin, tmp_full_bin, max_bin = nbin*mbin;
        cleanup_l(pl.pos,pk.pos,rkl_mag,rkl_mu);
#ifdef JACKKNIFE
        Float c4vj,JK_weight;
#endif
        tmp_bin = getbin(rkl_mag,rkl_mu);

        if ((tmp_bin<0)||(tmp_bin>=max_bin)) return; // if not in correct bin
        cleanup_l(pl.pos,pj.pos,rjl_mag,rjl_mu);
        xi_jl = cf24->xi(rjl_mag, rjl_mu); // j-l correlation

        for(int i=0;i<pln;i++){ // Iterate over particle in pi_list
            if(wijk[i]==-1) continue; // skip incorrect bins / ij self counts
            if(((prim_ids[i]==pl_id)&&(I1==I4))||((pj_id==pl_id)&&(I2==I4))||((pk_id==pl_id)&&(I3==I4))) continue; // don't self-count

            pi = pi_list[i];
            tmp_weight = wijk[i]*pl.w; // product of weights, w_i*w_j*w_k*w_l

            // Now compute the integral;
            c4v = tmp_weight/prob*2.*xi_ik[i]*xi_jl; // with xi_ik*xi_jl = xi_il*xi_jk symmetry factor
            // Compute jackknife weight tensor:
            tmp_full_bin = bin_ij[i]*mbin*nbin+tmp_bin;
            // Add to local counts
            c4[tmp_full_bin]+=c4v;
            binct4[tmp_full_bin]++;
#ifdef JACKKNIFE
            JK_weight=weight_tensor(int(pi.JK),int(pj.JK),int(pk.JK),int(pl.JK),bin_ij[i],tmp_bin, JK12, JK34, product_weights12_34);
            c4vj = c4v*JK_weight;
            c4j[tmp_full_bin]+=c4vj;
#endif
        }
    }

#ifdef JACKKNIFE
private:
    inline Float weight_tensor(const int Ji, const int Jj, const int Jk, const int Jl, const int bin_a, const int bin_b, JK_weights *JK_XY, JK_weights *JK_ZW, Float *product_weights){
        // Compute the jackknife weight tensor for jackknife regions J_i, J_j, J_k, J_l and w_aA weight matrix w_matrix. NB: J_x are collapsed jackknife indices here - only using the non-empty regions.
        // JK_weights class holds list of filled JKs, number of filled JKs and the weight matrix.
        // bin_a, bin_b specify the binning.
        // product_weights is a (flattened) matrix of Sum_A(w_aA * w_bA) for the relevant jackknife weights

        int nbins = JK_XY->nbins;

        // Compute q_ij^A q_kl^A term
        int first_weight = 0;
        if (Ji==Jk) first_weight++;
        if (Jj==Jk) first_weight++;
        if (Ji==Jl) first_weight++;
        if (Jj==Jl) first_weight++;

        // Compute q_ijw_bA _ q_klw_aA part
        Float second_weight=JK_ZW->weights[Ji*nbins+bin_b]+JK_ZW->weights[Jj*nbins+bin_b]+JK_XY->weights[Jk*nbins+bin_a]+JK_XY->weights[Jl*nbins+bin_a];

        // Compute total weight (third part is precomputed)
        Float total_weight = (Float)first_weight/4. - 0.5*second_weight + product_weights[bin_a*nbins+bin_b];

        return total_weight;
    }
#endif

public:
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
            mu = fabs(pos.dot(los)/norm/los.norm());
#else
            // In the periodic case use z-direction for mu
            mu = fabs(pos.z/norm);
#endif
        }
    }
public:
    void sum_ints(Integrals* ints) {
        // Add the values accumulated in ints to the corresponding internal sums
        for(int i=0;i<nbin*mbin;i++){
            Ra[i]+=ints->Ra[i];
            c2[i]+=ints->c2[i];
            binct[i]+=ints->binct[i];
#ifdef JACKKNIFE
            c2j[i]+=ints->c2j[i];
#endif
        }
        for(int i=0;i<nbin*mbin;i++){
            for(int j=0;j<nbin*mbin;j++){
                c3[i*nbin*mbin+j]+=ints->c3[i*nbin*mbin+j];
                c4[i*nbin*mbin+j]+=ints->c4[i*nbin*mbin+j];
                binct3[i*nbin*mbin+j]+=ints->binct3[i*nbin*mbin+j];
                binct4[i*nbin*mbin+j]+=ints->binct4[i*nbin*mbin+j];
#ifdef JACKKNIFE
                c3j[i*nbin*mbin+j]+=ints->c3j[i*nbin*mbin+j];
                c4j[i*nbin*mbin+j]+=ints->c4j[i*nbin*mbin+j];
            }
        }
        for(int i=0;i<n_jack;i++){
            for(int j=0;j<nbin*mbin;j++){
                EEaA1[i*nbin*mbin+j]+=ints->EEaA1[i*nbin*mbin+j];
                EEaA2[i*nbin*mbin+j]+=ints->EEaA2[i*nbin*mbin+j];
                RRaA1[i*nbin*mbin+j]+=ints->RRaA1[i*nbin*mbin+j];
                RRaA2[i*nbin*mbin+j]+=ints->RRaA2[i*nbin*mbin+j];
#endif
            }
        }
    }
#ifdef JACKKNIFE
    void frobenius_difference_sum(Integrals* ints, int n_loop, Float &frobC2, Float &frobC3, Float &frobC4, Float &frobC2j, Float &frobC3j, Float &frobC4j){
        Float self_c2j=0, diff_c2j=0;
        Float self_c3j=0, diff_c3j=0;
        Float self_c4j=0, diff_c4j=0;
#else
    void frobenius_difference_sum(Integrals* ints, int n_loop, Float &frobC2, Float &frobC3, Float &frobC4){
#endif
        // Add the values accumulated in ints to the corresponding internal sums and compute the Frobenius norm difference between integrals
        Float n_loops = (Float)n_loop;
        Float self_c2=0, diff_c2=0;
        Float self_c3=0, diff_c3=0;
        Float self_c4=0, diff_c4=0;
        // Compute Frobenius norms and sum integrals
        for(int i=0;i<nbin*mbin;i++){
            self_c2+=pow(c2[i]/n_loops,2.);
            diff_c2+=pow(c2[i]/n_loops-(c2[i]+ints->c2[i])/(n_loops+1.),2.);
#ifdef JACKKNIFE
            self_c2j+=pow(c2j[i]/n_loops,2.);
            diff_c2j+=pow(c2j[i]/n_loops-(c2j[i]+ints->c2j[i])/(n_loops+1.),2.);
#endif
            for(int j=0;j<nbin*mbin;j++){
                self_c4+=pow(c4[i*nbin*mbin+j]/n_loops,2.);
                diff_c4+=pow(c4[i*nbin*mbin+j]/n_loops-(c4[i*nbin*mbin+j]+ints->c4[i*nbin*mbin+j])/(n_loops+1.),2.);
                self_c3+=pow(c3[i*nbin*mbin+j]/n_loops,2.);
                diff_c3+=pow(c3[i*nbin*mbin+j]/n_loops-(c3[i*nbin*mbin+j]+ints->c3[i*nbin*mbin+j])/(n_loops+1.),2.);
#ifdef JACKKNIFE
                self_c4j+=pow(c4j[i*nbin*mbin+j]/n_loops,2.);
                diff_c4j+=pow(c4j[i*nbin*mbin+j]/n_loops-(c4j[i*nbin*mbin+j]+ints->c4j[i*nbin*mbin+j])/(n_loops+1.),2.);
                self_c3j+=pow(c3j[i*nbin*mbin+j]/n_loops,2.);
                diff_c3j+=pow(c3j[i*nbin*mbin+j]/n_loops-(c3j[i*nbin*mbin+j]+ints->c3j[i*nbin*mbin+j])/(n_loops+1.),2.);
#endif
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
#ifdef JACKKNIFE
        self_c2j=sqrt(self_c2j);
        diff_c2j=sqrt(diff_c2j);
        diff_c3j=sqrt(diff_c3j);
        diff_c4j=sqrt(diff_c4j);
        self_c3j=sqrt(self_c3j);
        self_c4j=sqrt(self_c4j);
        // Return percent difference
        frobC2j=100.*(diff_c2j/self_c2j);
        frobC3j=100.*(diff_c3j/self_c3j);
        frobC4j=100.*(diff_c4j/self_c4j);
#endif
        }

    void sum_total_counts(uint64& acc2, uint64& acc3, uint64& acc4){
        // Add local counts to total bin counts in acc2-4
        for (int i=0; i<nbin*mbin; i++) {
            acc2+=binct[i];
            for (int j=0; j<nbin*mbin; j++) {
                acc3+=binct3[i*nbin*mbin+j];
                acc4+=binct4[i*nbin*mbin+j];
            }
        }
    }
    void normalize(Float norm1, Float norm2, Float norm3, Float norm4, Float n_pairs, Float n_triples, Float n_quads){
        // Normalize the accumulated integrals (partly done by the normalising probabilities used from the selected cubes)
        // n_pair etc. are the number of PARTICLE pairs etc. attempted (not including rejected cells, but including pairs which don't fall in correct bin ranges)
        // To avoid recomputation
        double corrf2 = norm1*norm2; // correction factor for densities of random points
        double corrf3 = corrf2*norm3;
        double corrf4 = corrf3*norm4;;
#ifdef JACKKNIFE
        for(int i=0;i<n_jack;i++){
            for(int j=0;j<nbin*mbin;j++){
                // We get 1/4 of total pair counts in each EE bin
                EEaA1[i*nbin*mbin+j]/=(n_pairs*corrf2*0.25);
                EEaA2[i*nbin*mbin+j]/=(n_pairs*corrf2*0.25);
                RRaA1[i*nbin*mbin+j]/=(n_pairs*corrf2*0.25);
                RRaA2[i*nbin*mbin+j]/=(n_pairs*corrf2*0.25);
            }
        }
#endif
        for(int i = 0; i<nbin*mbin;i++){
            Ra[i]/=(n_pairs*corrf2);
            c2[i]/=(n_pairs*corrf2);
#ifdef JACKKNIFE
            c2j[i]/=(n_pairs*corrf2);
#endif
            for(int j=0;j<nbin*mbin;j++){
                c3[i*nbin*mbin+j]/=(n_triples*corrf3);
                c4[i*nbin*mbin+j]/=(n_quads*corrf4);
#ifdef JACKKNIFE
                c3j[i*nbin*mbin+j]/=(n_triples*corrf3);
                c4j[i*nbin*mbin+j]/=(n_quads*corrf4);
#endif
            }
        }

        // Further normalize by RR counts from corrfunc
        for(int i=0; i<nbin*mbin;i++){
            Float Ra_i = JK12->RR_pair_counts[i];
            c2[i]/=pow(Ra_i,2.); // must normalize by galaxy number here
#ifdef JACKKNIFE
            c2j[i]/=pow(Ra_i,2.)*(1.-product_weights12_12[i*nbin*mbin+i]);
#endif
            for(int j=0;j<nbin*mbin;j++){
                Float Rab3=Ra_i*JK23->RR_pair_counts[j];
                Float Rab4=Ra_i*JK34->RR_pair_counts[j];
                c3[i*nbin*mbin+j]/=Rab3;
                c4[i*nbin*mbin+j]/=Rab4;
#ifdef JACKKNIFE
                Float Rab_jk3 = Rab3*(1.-product_weights12_23[i*nbin*mbin+j]);
                Float Rab_jk4 = Rab4*(1.-product_weights12_34[i*nbin*mbin+j]);
                c3j[i*nbin*mbin+j]/=Rab_jk3;
                c4j[i*nbin*mbin+j]/=Rab_jk4;
#endif
            }
        }
    }
    void save_counts(uint64 pair_counts,uint64 triple_counts,uint64 quad_counts){
        // Print the counts for each integral (used for combining the estimates outside of C++)
        // This is the number of counts used in each loop [always the same]
        char counts_file[1000];
        snprintf(counts_file, sizeof counts_file, "%sCovMatricesAll/total_counts_n%d_m%d_%d%d,%d%d.txt",out_file,nbin,mbin,I1,I2,I3,I4);
        FILE * CountsFile = fopen(counts_file,"w");
        fprintf(CountsFile,"%llu\n",pair_counts);
        fprintf(CountsFile,"%llu\n",triple_counts);
        fprintf(CountsFile,"%llu\n",quad_counts);

        fflush(NULL);
        fclose(CountsFile);
    }


    void save_integrals(char* suffix, bool save_all) {
    /* Print integral outputs to file.
        * In txt files {c2,c3,c4,RR}_n{nbin}_m{mbin}.txt there are lists of the outputs of c2,c3,c4 and RR_a that are already normalized and multiplied by combinatoric factors. The n and m strings specify the number of n and m bins present.
        */
        // Create output files

        char c2name[1000];
        snprintf(c2name, sizeof c2name, "%sCovMatricesAll/c2_n%d_m%d_%d%d_%s.txt", out_file,nbin, mbin,I1,I2,suffix);
        char c3name[1000];
        snprintf(c3name, sizeof c3name, "%sCovMatricesAll/c3_n%d_m%d_%d,%d%d_%s.txt", out_file, nbin, mbin,I2,I1,I3,suffix);
        char c4name[1000];
        snprintf(c4name, sizeof c4name, "%sCovMatricesAll/c4_n%d_m%d_%d%d,%d%d_%s.txt", out_file, nbin, mbin, I1,I2,I3,I4,suffix);
        char RRname[1000];
        snprintf(RRname, sizeof RRname, "%sCovMatricesAll/RR_n%d_m%d_%d%d_%s.txt", out_file, nbin, mbin, I1,I2,suffix);
        FILE * C2File = fopen(c2name,"w"); // for c2 part of integral
        FILE * C3File = fopen(c3name,"w"); // for c3 part of integral
        FILE * C4File = fopen(c4name,"w"); // for c4 part of integral
        FILE * RRFile = fopen(RRname,"w"); // for RR part of integral

        for (int j=0;j<nbin*mbin;j++){
            fprintf(C2File,"%le\n",c2[j]);
            fprintf(RRFile,"%le\n",Ra[j]);
        }

        for(int i=0;i<nbin*mbin;i++){
            for(int j=0;j<nbin*mbin;j++){
                fprintf(C3File,"%le\t",c3[i*nbin*mbin+j]);
                fprintf(C4File,"%le\t",c4[i*nbin*mbin+j]);
            }
            fprintf(C3File,"\n"); // new line each end of row
            fprintf(C4File,"\n");
        }

        fflush(NULL);

        // Close open files
        fclose(C2File);
        fclose(C3File);
        fclose(C4File);
        fclose(RRFile);
        if(save_all==1){
            char binname[1000];
            snprintf(binname,sizeof binname, "%sCovMatricesAll/binct_c4_n%d_m%d_%d%d,%d%d_%s.txt",out_file, nbin,mbin,I1,I2,I3,I4,suffix);
            FILE * BinFile = fopen(binname,"w");

            char bin3name[1000];
            snprintf(bin3name,sizeof bin3name, "%sCovMatricesAll/binct_c3_n%d_m%d_%d,%d%d_%s.txt",out_file, nbin,mbin,I2,I1,I3,suffix);
            FILE * Bin3File = fopen(bin3name,"w");

            char bin2name[1000];
            snprintf(bin2name,sizeof bin2name, "%sCovMatricesAll/binct_c2_n%d_m%d_%d%d_%s.txt",out_file, nbin,mbin,I1,I2,suffix);
            FILE * Bin2File = fopen(bin2name,"w");

            for (int j=0;j<nbin*mbin;j++){
                fprintf(Bin2File,"%llu\n",binct[j]);
            }

            for(int i=0;i<nbin*mbin;i++){
                for(int j=0;j<nbin*mbin;j++){
                    fprintf(BinFile,"%llu\t",binct4[i*nbin*mbin+j]);
                    fprintf(Bin3File,"%llu\t",binct3[i*nbin*mbin+j]);
                    }
                fprintf(BinFile,"\n");
                fprintf(Bin3File,"\n");
            }
            fclose(BinFile);
            fclose(Bin2File);
            fclose(Bin3File);
        }
    }

#ifdef JACKKNIFE
    void save_jackknife_integrals(char* suffix) {
    /* Print jackknife integral outputs to file.
        * In txt files {c2,c3,c4,RR}_n{nbin}_m{mbin}.txt there are lists of the outputs of c2,c3,c4 and RR_a that are already normalized and multiplied by combinatoric factors. The n and m strings specify the number of n and m bins present.
        */
        // Create output files

        char c2name[1000];
        snprintf(c2name, sizeof c2name, "%sCovMatricesJack/c2_n%d_m%d_%d%d_%s.txt", out_file,nbin, mbin,I1,I2,suffix);
        char c3name[1000];
        snprintf(c3name, sizeof c3name, "%sCovMatricesJack/c3_n%d_m%d_%d,%d%d_%s.txt", out_file, nbin, mbin,I2,I1,I3,suffix);
        char c4name[1000];
        snprintf(c4name, sizeof c4name, "%sCovMatricesJack/c4_n%d_m%d_%d%d,%d%d_%s.txt", out_file, nbin, mbin, I1,I2,I3,I4,suffix);
        FILE * C2File = fopen(c2name,"w"); // for c2 part of integral
        FILE * C3File = fopen(c3name,"w"); // for c3 part of integral
        FILE * C4File = fopen(c4name,"w"); // for c4 part of integral
        char RR1name[1000];
        snprintf(RR1name,sizeof RR1name, "%sCovMatricesJack/RR1_n%d_m%d_%d%d_%s.txt", out_file, nbin, mbin, I1,I2,suffix);
        char RR2name[1000];
        snprintf(RR2name,sizeof RR2name, "%sCovMatricesJack/RR2_n%d_m%d_%d%d_%s.txt", out_file, nbin, mbin, I1,I2,suffix);
        char EE1name[1000];
        snprintf(EE1name,sizeof EE1name, "%sCovMatricesJack/EE1_n%d_m%d_%d%d_%s.txt", out_file, nbin, mbin, I1,I2,suffix);
        char EE2name[1000];
        snprintf(EE2name,sizeof EE2name, "%sCovMatricesJack/EE2_n%d_m%d_%d%d_%s.txt", out_file, nbin, mbin, I1,I2,suffix);
        FILE * EE1File = fopen(EE1name, "w"); // for EE1 integral
        FILE * EE2File = fopen(EE2name,"w"); // for EE2 integral
        FILE * RR1File = fopen(RR1name, "w"); // for RR1 integral
        FILE * RR2File = fopen(RR2name,"w"); // for RR2 integral
        for (int j=0;j<nbin*mbin;j++){
            fprintf(C2File,"%le\n",c2j[j]);
        }
        for(int i=0;i<nbin*mbin;i++){
            for(int j=0;j<nbin*mbin;j++){
                fprintf(C3File,"%le\t",c3j[i*nbin*mbin+j]);
                fprintf(C4File,"%le\t",c4j[i*nbin*mbin+j]);
            }
            fprintf(C3File,"\n"); // new line each end of row
            fprintf(C4File,"\n");
        }
        for(int i=0;i<n_jack;i++){
            for(int j=0;j<nbin*mbin;j++){
                fprintf(EE1File,"%le\t",EEaA1[i*nbin*mbin+j]);
                fprintf(EE2File,"%le\t",EEaA2[i*nbin*mbin+j]);
                fprintf(RR1File,"%le\t",RRaA1[i*nbin*mbin+j]);
                fprintf(RR2File,"%le\t",RRaA2[i*nbin*mbin+j]);
            }
            fprintf(EE1File,"\n");
            fprintf(EE2File,"\n");
            fprintf(RR1File,"\n");
            fprintf(RR2File,"\n");
        }

        fflush(NULL);

        // Close open files
        fclose(C2File);
        fclose(C3File);
        fclose(C4File);
        fclose(EE1File);
        fclose(EE2File);
        fclose(RR1File);
        fclose(RR2File);
    }
#endif

};

#endif
