// Rewritten integrals.h code for grid_covariance.cpp (originally from Alex Wiegand) to parallelize and compute integrands to a given quad of particles
#include "parameters.h"
#include "correlation_function.h"
#include "cell_utilities.h"

#ifndef INTEGRALS_2_H
#define INTEGRALS_2_H

class Integrals2{
public:
    CorrelationFunction *cf;
    
private:
    int nbin, mbin;
    Float rmin,rmax,mumin,mumax,dr,dmu; //Ranges in r and mu
    Float *Ra, *c2, *c3, *c4; // Arrays to accumulate integrals
    Float *cxj, *c2j, *c3j, *c4j; // Arrays to accumulate jackknife integrals
    Float *errc4; // Integral to house the variance in C4;
    
    bool box,rad; // Flags to decide whether we have a periodic box and if we have a radial correlation function only
    
    uint64 *binct, *binct3, *binct4; // Arrays to accumulate bin counts
    
public:
    Integrals2(Parameters *par, CorrelationFunction *_cf){
        cf = new CorrelationFunction(_cf);
        init(par);
    }
    
    void init(Parameters *par){
        nbin = par->nbin; // number of radial bins
        mbin = par->mbin; // number of mu bins
        
        int ec=0;
        // Initialize the binning
        ec+=posix_memalign((void **) &Ra, PAGE, sizeof(double)*nbin*mbin);
        ec+=posix_memalign((void **) &c2, PAGE, sizeof(double)*nbin*mbin);
        ec+=posix_memalign((void **) &c3, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);
        ec+=posix_memalign((void **) &c4, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);
        
        ec+=posix_memalign((void **) &c2j, PAGE, sizeof(double)*nbin*mbin);
        ec+=posix_memalign((void **) &c3j, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);
        ec+=posix_memalign((void **) &c4j, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);
        ec+=posix_memalign((void **) &cxj, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);
        ec+=posix_memalign((void **) &errc4, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);
        
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

    }

    ~Integrals2() {
        free(Ra);
        free(c2);
        free(c3);
        free(c4);
        free(cxj);
        free(c2j);
        free(c3j);
        free(c4j);
        free(errc4);
        free(binct);
        free(binct3);
        free(binct4);

    }

    void reset(){
        for (int j=0; j<nbin*mbin; j++) {
            Ra[j] = 0;
            c2[j] = 0;
            c2j[j] = 0;
            binct[j] = 0;
        }
        for (int j=0; j<nbin*mbin*nbin*mbin; j++) {
            c3[j]=0;
            c4[j]=0;
            c3j[j]=0;
            c4j[j]=0;
            cxj[j] = 0;
            errc4[j] = 0;
            binct3[j] = 0;
            binct4[j] = 0;
        }
    }

    inline int getbin(Float r, Float mu){
        // Linearizes 2D indices
        return floor((r-rmin)/dr) *mbin + floor((mu-mumin)/dmu);
    }
    
    inline void second(const Particle* pi_list, const int* prim_ids, int pln, const Particle pj, const int pj_id, Float* &xi_ij, int* &bin, Float* &wij, const double prob){
        // Accumulates the two point integral C2. Also outputs an array of xi_ij and bin values for later reuse.
        
        for(int i=0;i<pln;i++){ // Iterate over particle in pi_list
                Float rij_mag,rij_mu, rav, c2v, c2vj;
                Particle pi = pi_list[i]; // first particle
                
                if(prim_ids[i]==pj_id){
                    wij[i]=-1;
                    continue; // don't self-count
                }
                
                cleanup_l(pi.pos,pj.pos,rij_mag,rij_mu); // define |r_ij| and ang(r_ij)
                int tmp_bin = getbin(rij_mag, rij_mu); // bin for each particle
                
                if ((tmp_bin<0)||(tmp_bin>=nbin*mbin)){
                    wij[i]=-1;
                    continue; // if not in correct bin
                }
                
                Float tmp_weight = pi.w*pj.w; // product of weights
                Float tmp_xi = cf->xi(rij_mag, rij_mu); // correlation function for i-j
                
                // Save into arrays for later
                bin[i]=tmp_bin;
                xi_ij[i] = tmp_xi;
                wij[i] = tmp_weight;
                
                // Now compute the integral:
                c2vj = tmp_weight*tmp_weight*(1.+tmp_xi)/prob;
                c2v = tmp_weight*tmp_weight*(1.+tmp_xi) / prob; // c2 contribution
                rav = tmp_weight / prob; // RR_a contribution
                
                // Add to local integral counts:
                Ra[tmp_bin]+=rav;
                c2[tmp_bin]+=c2v;
                c2j[tmp_bin]+=c2vj;
                binct[tmp_bin]++;
        }
    }
    
    inline void third(const Particle* pi_list, const int* prim_ids, const int pln, const Particle pj, const Particle pk, const int pj_id, const int pk_id, const int* bin_ij, const Float* wij, Float* &xi_jk, Float* &xi_ik, Float* wijk, const double prob){
        // Accumulates the three point integral C3. Also outputs an array of xi_ik and bin_ik values for later reuse.
        
        for(int i=0;i<pln;i++){ // Iterate ovr particle in pi_list
            Particle pi = pi_list[i];
            Float rik_mag,rik_mu,c3v,c3vj,rjk_mag, rjk_mu;
            if(wij[i]==-1){
                wijk[i]=-1;
                continue; // skip incorrect bins / ij self counts
            }
            
            if((prim_ids[i]==pk_id)||(pk_id==pj_id)){
                wijk[i]=-1; // re-read to skip this later
                continue; // don't self-count
            }
            
            cleanup_l(pi.pos,pk.pos,rik_mag,rik_mu); // define angles/lengths
            
            int tmp_bin = getbin(rik_mag,rik_mu); // bin for each particle
            
            if ((tmp_bin<0)||(tmp_bin>=nbin*mbin)){
                wijk[i] = -1;
                continue; // if not in correct bin
            }
            
            cleanup_l(pj.pos,pk.pos,rjk_mag,rjk_mu); 
            
            Float tmp_weight = wij[i]*pk.w; // product of weights, w_iw_jw_k
            Float xi_jk_tmp = cf->xi(rjk_mag, rjk_mu); // correlation function for j-k
            Float xi_ik_tmp = cf->xi(rik_mag, rik_mu); // not used here but used later
            
            // save arrays for later
            xi_jk[i]=xi_jk_tmp;
            xi_ik[i]=xi_ik_tmp;
            wijk[i]=tmp_weight;
            
            // Now compute the integral;
            c3v = tmp_weight*pi.w/prob*xi_jk_tmp;
            c3vj = tmp_weight*pi.w/prob*xi_jk_tmp;
            
            // Add to local counts
            int tmp_full_bin = bin_ij[i]*mbin*nbin+tmp_bin;
            c3[tmp_full_bin]+=c3v;
            c3j[tmp_full_bin]+=c3vj;
            binct3[tmp_full_bin]++;
        }
    }
        
    inline void fourth(const Particle* pi_list, const int* prim_ids, const int pln, const Particle pj, const Particle pk, const Particle pl, const int pj_id, const int pk_id, const int pl_id, const int* bin_ij, const Float* wijk, const Float* xi_ik, const Float* xi_jk, const Float* xi_ij, const double prob){
        // Accumulates the three point integral C3. Also outputs an array of xi_ik and bin_ik values for later reuse.
        
        for(int i=0;i<pln;i++){ // Iterate ovr particle in pi_list
            Particle pi = pi_list[i];
            Float ril_mag,ril_mu,rjl_mag, rjl_mu, rkl_mag, rkl_mu, c4v, c4vj, cxvj;
            if(wijk[i]==-1) continue; // skip incorrect bins / ij self counts
            
            if((prim_ids[i]==pl_id)||(pj_id==pl_id)||(pk_id==pl_id)) continue; // don't self-count
            
            cleanup_l(pl.pos, pk.pos, rkl_mag, rkl_mu); // define angles/lengths
            int tmp_bin = getbin(rkl_mag,rkl_mu); // kl bin for each particle
            
            if ((tmp_bin<0)||(tmp_bin>=nbin*mbin)) continue; // if not in correct bin
            
            cleanup_l(pl.pos,pi.pos,ril_mag,ril_mu); 
            cleanup_l(pl.pos,pj.pos,rjl_mag,rjl_mu); 
            
            Float tmp_weight = wijk[i]*pl.w; // product of weights, w_iw_jw_kw_l
            Float xi_il = cf->xi(ril_mag, ril_mu); // correlation function for i-l
            Float xi_jl = cf->xi(rjl_mag, rjl_mu); // j-l correlation
            Float xi_kl = cf->xi(rkl_mag, rkl_mu); // k-l correlation
            
            // Now compute the integral;
            c4v = tmp_weight/prob*(xi_il*xi_jk[i]+xi_jl*xi_ik[i]);
            c4vj = tmp_weight/prob*(xi_il*xi_jk[i]+xi_jl*xi_ik[i]);
            cxvj = tmp_weight/prob*(xi_ij[i]*xi_kl);
            
            // Add to local counts
            int tmp_full_bin = bin_ij[i]*mbin*nbin+tmp_bin;
            c4[tmp_full_bin]+=c4v;
            c4j[tmp_full_bin]+=c4vj;
            cxj[tmp_full_bin]+=cxvj;
            errc4[tmp_full_bin]+=pow(c4v,2.);
            binct4[tmp_full_bin]++;
        }
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
            mu = fabs(pos.dot(los)/norm/los.norm());
#else
            // In the periodic case use z-direction for mu
            mu = fabs(pos.z/norm);
#endif
        }
    }
public:
    void sum_ints(Integrals2* ints) {
        // Add the values accumulated in ints to the corresponding internal sums
        for(int i=0;i<nbin*mbin;i++){
            Ra[i]+=ints->Ra[i];
            c2[i]+=ints->c2[i];
            c2j[i]+=ints->c2j[i];
            binct[i]+=ints->binct[i];
        }
        for(int i=0;i<nbin*mbin;i++){
            for(int j=0;j<nbin*mbin;j++){
                c3[i*nbin*mbin+j]+=ints->c3[i*nbin*mbin+j];
                c3j[i*nbin*mbin+j]+=ints->c3j[i*nbin*mbin+j];
                c4[i*nbin*mbin+j]+=ints->c4[i*nbin*mbin+j];
                c4j[i*nbin*mbin+j]+=ints->c4j[i*nbin*mbin+j];
                cxj[i*nbin*mbin+j]+=ints->cxj[i*nbin*mbin+j];
                binct3[i*nbin*mbin+j]+=ints->binct3[i*nbin*mbin+j];
                binct4[i*nbin*mbin+j]+=ints->binct4[i*nbin*mbin+j];
            }
        }
    }
    
    void frobenius_difference_sum(Integrals2* ints, int n_loop, Float &frobC2, Float &frobC3, Float &frobC4){
        // Add the values accumulated in ints to the corresponding internal sums and compute the Frobenius norm difference between integrals
        Float n_loops = (Float)n_loop;
        Float self_c2=0, diff_c2=0;
        Float self_c3=0, diff_c3=0;
        Float self_c4=0, diff_c4=0;
        
        // Compute Frobenius norms
        for(int i=0;i<nbin*mbin;i++){
            self_c2+=pow(c2[i]/n_loops,2.);
            diff_c2+=pow(c2[i]/n_loops-(c2[i]+ints->c2[i])/(n_loops+1.),2.);
            
            for(int j=0;j<nbin*mbin;j++){
                self_c4+=pow(c4[i*nbin*mbin+j]/n_loops,2.);
                diff_c4+=pow(c4[i*nbin*mbin+j]/n_loops-(c4[i*nbin*mbin+j]+ints->c4[i*nbin*mbin+j])/(n_loops+1.),2.);
                self_c3+=pow(c3[i*nbin*mbin+j]/n_loops,2.);
                diff_c3+=pow(c3[i*nbin*mbin+j]/n_loops-(c3[i*nbin*mbin+j]+ints->c3[i*nbin*mbin+j])/(n_loops+1.),2.);
                c3[i*nbin*mbin+j]+=ints->c3[i*nbin*mbin+j];
                c4[i*nbin*mbin+j]+=ints->c4[i*nbin*mbin+j];
                binct3[i*nbin*mbin+j]+=ints->binct3[i*nbin*mbin+j];
                binct4[i*nbin*mbin+j]+=ints->binct4[i*nbin*mbin+j];
            }
            c2[i]+=ints->c2[i];
            binct[i]+=ints->binct[i];
        }
        
        // Now update Ra values (must be separate from above calculation)
        for(int i=0;i<nbin*mbin;i++){
            Ra[i]+=ints->Ra[i];
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
        for (int i=0; i<nbin*mbin; i++) {
            acc2+=binct[i];
            for (int j=0; j<nbin*mbin; j++) {
                acc3+=binct3[i*nbin*mbin+j];
                acc4+=binct4[i*nbin*mbin+j];
            }
        }
    }
    
    void normalize(int np, int ngal, Float n_pairs, Float n_triples, Float n_quads, bool use_RR){
        // Normalize the accumulated integrals (partly done by the normalising probabilities used from the selected cubes)
        // np is the number of random particles used, ngal is the number of galaxies in the survey
        // n_pair etc. are the number of PARTICLE pairs etc. attempted (not including rejected cells, but including pairs which don't fall in correct bin ranges)
        // If use_RR=True, we normalize by RR_a, RR_b also
        double corrf = (double)np/ngal; // Correction factor for the different densities of random points
        
        for(int i = 0; i<nbin*mbin;i++){
            Ra[i]/=(n_pairs*corrf*corrf);
            c2[i]/=(n_pairs*corrf*corrf);
            c2j[i]/=(n_pairs*corrf*corrf);
            for(int j=0;j<nbin*mbin;j++){
                c3[i*nbin*mbin+j]/=(n_triples*pow(corrf,3));
                c4[i*nbin*mbin+j]/=(n_quads*pow(corrf,4));
                c3j[i*nbin*mbin+j]/=(n_triples*pow(corrf,3));
                c4j[i*nbin*mbin+j]/=(n_quads*pow(corrf,4));
                cxj[i*nbin*mbin+j]/=(n_quads*pow(corrf,4));
            }
        }
        if(use_RR==1){
            // Also normalize by RR counts
            for(int i=0; i<nbin*mbin;i++){
                c2[i]/=(Ra[i]*Ra[i]);
                c2j[i]/=(Ra[i]*Ra[i]);
                for(int j=0;j<nbin*mbin;j++){
                    Float Ra_ij = Ra[i]*Ra[j];
                    c3[i*nbin*mbin+j]/=Ra_ij;
                    c4[i*nbin*mbin+j]/=Ra_ij;
                    c3j[i*nbin*mbin+j]/=Ra_ij;
                    c4j[i*nbin*mbin+j]/=Ra_ij;
                    cxj[i*nbin*mbin+j]/=Ra_ij;
                }
            }
        }
    }
        
    
    void save_integrals(char* suffix) {
    /* Print integral outputs to file. 
        * In txt files {c2,c3,c4,RR}_n{nbin}_m{mbin}.txt there are lists of the outputs of c2,c3,c4 and RR_a that are already normalized and multiplied by combinatoric factors. The n and m strings specify the number of n and m bins present.
        */
        // Create output files
        char c2name[1000];
        snprintf(c2name, sizeof c2name, "CovMatricesAll/c2_n%d_m%d_%s.txt", nbin, mbin,suffix);
        char c3name[1000];
        snprintf(c3name, sizeof c3name, "CovMatricesAll/c3_n%d_m%d_%s.txt", nbin, mbin, suffix);
        char c4name[1000];
        snprintf(c4name, sizeof c4name, "CovMatricesAll/c4_n%d_m%d_%s.txt", nbin, mbin, suffix);
        char RRname[1000];
        snprintf(RRname, sizeof RRname, "CovMatricesAll/RR_n%d_m%d_%s.txt", nbin, mbin, suffix);
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
                    
    }

    
    void save_jackknife_integrals(char* suffix) {
    /* Print jackknife integral outputs to file. 
        * In txt files {c2,c3,c4,RR}_n{nbin}_m{mbin}.txt there are lists of the outputs of c2,c3,c4 and RR_a that are already normalized and multiplied by combinatoric factors. The n and m strings specify the number of n and m bins present.
        */
        // Create output files
        char c2name[1000];
        snprintf(c2name, sizeof c2name, "CovMatricesJack/c2j_n%d_m%d_%s.txt", nbin, mbin,suffix);
        char c3name[1000];
        snprintf(c3name, sizeof c3name, "CovMatricesJack/c3j_n%d_m%d_%s.txt", nbin, mbin, suffix);
        char c4name[1000];
        snprintf(c4name, sizeof c4name, "CovMatricesJack/c4j_n%d_m%d_%s.txt", nbin, mbin, suffix);
        char cxname[1000];
        snprintf(cxname, sizeof cxname, "CovMatricesJack/cxj_n%d_m%d_%s.txt", nbin, mbin, suffix);
        FILE * C2File = fopen(c2name,"w"); // for c2 part of integral
        FILE * C3File = fopen(c3name,"w"); // for c3 part of integral
        FILE * C4File = fopen(c4name,"w"); // for c4 part of integral
        FILE * CXFile = fopen(cxname,"w"); // for RR part of integral
        
        for (int j=0;j<nbin*mbin;j++){
            fprintf(C2File,"%le\n",c2j[j]);
        }
        for(int i=0;i<nbin*mbin;i++){
            for(int j=0;j<nbin*mbin;j++){
                fprintf(C3File,"%le\t",c3j[i*nbin*mbin+j]);
                fprintf(C4File,"%le\t",c4j[i*nbin*mbin+j]);
                fprintf(CXFile,"%le\t",cxj[i*nbin*mbin+j]);
            }
            fprintf(C3File,"\n"); // new line each end of row
            fprintf(C4File,"\n");
            fprintf(CXFile,"\n");
        }
        fflush(NULL);
    }
    
    void compute_Neff(Float n_quads){
    // Compute the effective N from the data given the number of sets of quads used 
        for(int i=0;i<nbin*mbin;i++){
            for(int j=0;j<nbin*mbin;j++){
                Float var_c4; // variance in this bin
                var_c4=errc4[i*nbin*mbin+j]n_quads-pow(c4[i*nbin*mbin+j],2.)
                Float neff;
                neff = pow(c4[i*nbin*mbin+j],2.)+c4[i*nbin*mbin+i]*c4[j*nbin*mbin+j];
        
    }

};

#endif
