// Rewritten integrals.h code for legendre multipoles
#include "parameters.h"
#include "correlation_function.h"
#include "cell_utilities.h"
#include "legendre_utilities.h"

#ifndef INTEGRALS_POWER_H
#define INTEGRALS_POWER_H

class Integrals{
public:
    CorrelationFunction *cf12, *cf13, *cf24;
<<<<<<< HEAD

=======
    
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
private:
    int nbin, mbin, max_legendre;
    Float R0; // truncation radius
    Float rmin,rmax,mumin,mumax; //Ranges in r and mu
    Float *r_high, *r_low; // Max and min of each radial bin
    Float *c2, *c3, *c4; // Arrays to accumulate integrals
    char* out_file;
    bool box; // Flags to decide whether we have a periodic box
    int I1, I2, I3, I4; // indices for which fields to use for each particle
<<<<<<< HEAD

    SurveyCorrection *sc12,*sc23,*sc34; // survey correction function
    KernelInterp *kernel_interp; // kernel interpolation function

public:
    uint64 used_pairs, used_triples, used_quads; // total number of particles used

=======
    
    SurveyCorrection *sc12,*sc23,*sc34; // survey correction function
    KernelInterp *kernel_interp; // kernel interpolation function
    
public:
    uint64 used_pairs, used_triples, used_quads; // total number of particles used
    
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
public:
    Integrals(){};

    Integrals(Parameters *par, CorrelationFunction *_cf12, CorrelationFunction *_cf13, CorrelationFunction *_cf24, int _I1, int _I2, int _I3, int _I4, SurveyCorrection *_sc12, SurveyCorrection *_sc23, SurveyCorrection *_sc34, KernelInterp *_kernel_interp){
        sc12 = _sc12;
        sc23 = _sc23;
        sc34 = _sc34;
<<<<<<< HEAD

=======
        
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
        cf12 = new CorrelationFunction(_cf12);
        cf13 = new CorrelationFunction(_cf13);
        cf24 = new CorrelationFunction(_cf24);
        kernel_interp = new KernelInterp(_kernel_interp);
        I1=_I1;
        I2=_I2;
        I3=_I3;
        I4=_I4;
        init(par);
        reset();
    }
<<<<<<< HEAD

=======
    
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
    void init(Parameters *par){
        nbin = par->nbin; // number of radial bins
        mbin = par->mbin; // number of mu bins
        R0 = par->R0; // truncation radius
        out_file = par->out_file; // output directory
        max_legendre = fmax((sc12->l_bins-1)*2,par->max_l); // maximum Legendre polynomial used

        int ec=0;
        // Initialize the binning
        ec+=posix_memalign((void **) &c2, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);
        ec+=posix_memalign((void **) &c3, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);
        ec+=posix_memalign((void **) &c4, PAGE, sizeof(double)*nbin*mbin*nbin*mbin);

        assert(ec==0);
        reset();
<<<<<<< HEAD

=======
        
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
        box=par->perbox;

        rmax=par->rmax;
        rmin=par->rmin;

        mumax=par->mumax;
        mumin=par->mumin;
<<<<<<< HEAD

        r_high = par->radial_bins_high;
        r_low = par->radial_bins_low;

=======
        
        r_high = par->radial_bins_high;
        r_low = par->radial_bins_low;
        
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
    }

    ~Integrals() {
        free(c2);
        free(c3);
        free(c4);
    }

    void reset(){
        for (int j=0; j<nbin*mbin*nbin*mbin; j++) {
            c2[j]=0;
            c3[j]=0;
            c4[j]=0;
        }
        used_pairs = 0;
        used_triples = 0;
        used_quads = 0;
    }
<<<<<<< HEAD

=======
    
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c

    inline void second(const Particle* pi_list, const int* prim_ids, int pln, const Particle pj, const int pj_id, Float* &wij, const double prob, Float* &kernel_ij){
        // Accumulates the two point integral C2.
        // Prob. here is defined as g_ij / f_ij where g_ij is the sampling PDF and f_ij is the true data PDF for picking pairs (equal to n_i/N n_j/N for N particles)
<<<<<<< HEAD

=======
        
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
        Float tmp_weight, tmp_xi, rij_mag, rij_mu, c2v, tmp_phi_inv, this_kernel[nbin*mbin], old_kernel[mbin];
        Float new_kr, old_kr, new_kernel, legendre[max_legendre/2+1];
        Particle pi;
        int tmp_bin;
<<<<<<< HEAD

        for(int i=0;i<pln;i++){ // Iterate over particle in pi_list
            if((prim_ids[i]==pj_id)&&(I1==I2)){
=======
        
        for(int i=0;i<pln;i++){ // Iterate over particle in pi_list
            if(prim_ids[i]==pj_id){
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
                wij[i]=-1;
                continue; // don't self-count
            }
            pi = pi_list[i]; // first particle
            cleanup_l(pi.pos,pj.pos,rij_mag,rij_mu); // define |r_ij| and ang(r_ij)
<<<<<<< HEAD

=======
            
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
            if(rij_mag>R0){
                wij[i] = -1;
                continue; // outside correct size
            }
<<<<<<< HEAD

            used_pairs++; // count a pair

            tmp_xi = cf12->xi(rij_mag, rij_mu); // correlation function for i-j

            // load all relevant legendre polynomials
            legendre_polynomials(rij_mu, max_legendre, legendre);

            // Compute 1/Phi function
            tmp_phi_inv=0.;
            for(int l_i=0;l_i<sc12->l_bins;l_i++) tmp_phi_inv+=legendre[l_i]*sc12->inv_correction_function(l_i*2,rij_mu);

=======
            
            used_pairs++; // count a pair  
            
            tmp_xi = cf12->xi(rij_mag, rij_mu); // correlation function for i-j
            
            // load all relevant legendre polynomials
            legendre_polynomials(rij_mu, max_legendre, legendre);
            
            // Compute 1/Phi function
            tmp_phi_inv=0.;
            for(int l_i=0;l_i<sc12->l_bins;l_i++) tmp_phi_inv+=legendre[l_i]*sc12->inv_correction_function(l_i*2,rij_mu);
            
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
#define UNBINNED
#ifdef UNBINNED
            tmp_weight = pi.w*pj.w*pair_weight(rij_mag)/tmp_phi_inv; // w_i*w_j*W(r_ij; R_0)*Phi_ij
#else
            tmp_weight = pi.w*pj.w*pair_weight(rij_mag)/tmp_phi_inv/pow(rij_mag,3); // w_i*w_j*W(r_ij; R_0)*Phi_ij
#endif
            // Save into arrays for later
<<<<<<< HEAD
            wij[i] = tmp_weight;

            // Now add to relevant bins
            c2v = tmp_weight*tmp_weight*(1.+tmp_xi) / prob*2.; // c2 contribution with symmetry factor

=======
            wij[i] = tmp_weight;                
            
            // Now add to relevant bins
            c2v = tmp_weight*tmp_weight*(1.+tmp_xi) / prob*2.; // c2 contribution with symmetry factor
            
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
            // Now compute multipole row elements
#ifdef UNBINNED
            for(int ii=0;ii<nbin;ii++){
                old_kr = rij_mag*0.5*(r_low[ii]+r_high[ii]);
                new_kr = sin(old_kr)/old_kr;
                for(int jj=0;jj<mbin;jj++){
                    this_kernel[ii*mbin+jj] = new_kr*legendre[jj];
                }
            }
#else
            for(int ii=0;ii<nbin;ii++){
                if(ii==0) old_kr = rij_mag*r_low[ii];
                new_kr = rij_mag*r_high[ii];
                for(int jj=0;jj<mbin;jj++){
                    if(ii==0) old_kernel[jj] = kernel_interp->kernel(jj*2.,old_kr);
                    new_kernel = kernel_interp->kernel(jj*2,new_kr);
                    this_kernel[ii*mbin+jj] = legendre[jj]*(new_kernel - old_kernel[jj]);
                    kernel_ij[i*nbin*mbin+ii*mbin+jj] = this_kernel[ii*mbin+jj];
                    old_kernel[jj] = new_kernel;
                }
                old_kr = new_kr;
            }
#endif
            // Fill up matrix
            for(int a=0;a<nbin;a++){
                for(int p=0;p<mbin;p++){
                    tmp_bin = (a*mbin+p)*nbin*mbin;
                    for(int b=0;b<nbin;b++){
                        for(int q=0;q<mbin;q++){
                            c2[tmp_bin+b*mbin+q]+=c2v*this_kernel[a*mbin+p]*this_kernel[b*mbin+q];
                        }
                    }
                }
            }
        }
    }

    inline void third(const Particle* pi_list, const int* prim_ids, const int pln, const Particle pj, const Particle pk, const int pj_id, const int pk_id, const Float* wij, Float* &xi_ik, Float* wijk, const double prob, const Float* kernel_ij){
        // Accumulates the three point integral C3. Also outputs an array of xi_ik and bin_ik values for later reuse.
<<<<<<< HEAD

=======
        
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
        // First define variables:
        Particle pi;
        Float rik_mag, rik_mu, rjk_mag, rjk_mu, c3v, tmp_phi_inv=0, this_kernel[nbin*mbin], old_kernel[mbin];
        Float tmp_weight, tmp_kernel, xi_ik_tmp, new_kr, old_kr, new_kernel, legendre_jk[max_legendre/2+1];
        int tmp_bin;
<<<<<<< HEAD

        // Define jk distance and angle
        cleanup_l(pj.pos,pk.pos,rjk_mag,rjk_mu);

        if(rjk_mag<R0){
            // Only do this for usable jk separations
            // but don't exit script yet, since xi_ik is still needed for later

            // load legendre polynomials
            legendre_polynomials(rjk_mu,max_legendre,legendre_jk);

            // Compute 1/Phi_jk function
            for(int l_i=0;l_i<sc23->l_bins;l_i++) tmp_phi_inv+=legendre_jk[l_i]*sc23->inv_correction_function(l_i*2,rjk_mu);

=======
        
        // Define jk distance and angle
        cleanup_l(pj.pos,pk.pos,rjk_mag,rjk_mu);
        
        if(rjk_mag<R0){
            // Only do this for usable jk separations
            // but don't exit script yet, since xi_ik is still needed for later
            
            // load legendre polynomials
            legendre_polynomials(rjk_mu,max_legendre,legendre_jk);
        
            // Compute 1/Phi_jk function
            for(int l_i=0;l_i<sc23->l_bins;l_i++) tmp_phi_inv+=legendre_jk[l_i]*sc23->inv_correction_function(l_i*2,rjk_mu);
            
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
            // Now compute multipole row elements
#ifdef UNBINNED
            for(int ii=0;ii<nbin;ii++){
                old_kr = rjk_mag*0.5*(r_low[ii]+r_high[ii]);
                new_kr = sin(old_kr)/old_kr;
                for(int jj=0;jj<mbin;jj++){
                    this_kernel[ii*mbin+jj] = new_kr*legendre_jk[jj];
                }
            }
        }
            tmp_kernel = pj.w*pair_weight(rjk_mag)*4./(tmp_phi_inv*prob);
#else
            for(int ii=0;ii<nbin;ii++){
                if(ii==0) old_kr = rjk_mag*r_low[ii];
                new_kr = rjk_mag*r_high[ii];
                for(int jj=0;jj<mbin;jj++){
                    if(ii==0) old_kernel[jj] = kernel_interp->kernel(jj*2.,old_kr);
                    new_kernel = kernel_interp->kernel(jj*2,new_kr);
                    this_kernel[ii*mbin+jj] = legendre_jk[jj]*(new_kernel - old_kernel[jj]);
                    old_kernel[jj] = new_kernel;
                }
                old_kr = new_kr;
            }
        }
        tmp_kernel = pj.w*pair_weight(rjk_mag)*4./(tmp_phi_inv*pow(rjk_mag,3)*prob);
<<<<<<< HEAD
#endif

        for(int i=0;i<pln;i++){ // Iterate over particle in pi_list
            if(((pk_id==pj_id)&&(I2==I3))||((prim_ids[i]==pk_id)&&(I1==I3))||(wij[i]==-1)){
=======
#endif        
        
        for(int i=0;i<pln;i++){ // Iterate over particle in pi_list
            if((wij[i]==-1)||(prim_ids[i]==pk_id)){
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
                wijk[i]=-1;
                continue; // skip incorrect bins / ij,jk self counts
            }
            pi = pi_list[i];
            cleanup_l(pi.pos,pk.pos,rik_mag,rik_mu); // define angles/lengths
<<<<<<< HEAD
            xi_ik_tmp = cf13->xi(rik_mag, rik_mu);

            tmp_weight = wij[i]*pk.w; // product of weights, w_iw_jw_k *W(r_ij,R0) * Phi_ij/r_ij^3

            // save arrays for later
            xi_ik[i]=xi_ik_tmp;
            wijk[i]=tmp_weight;

            if(rjk_mag>R0) continue;

            used_triples++; // count a triple

            // Now compute the integral;
            c3v = tmp_kernel*tmp_weight*xi_ik_tmp; // include symmetry factor

=======
            xi_ik_tmp = cf13->xi(rik_mag, rik_mu); 
            
            tmp_weight = wij[i]*pk.w; // product of weights, w_iw_jw_k *W(r_ij,R0) * Phi_ij/r_ij^3
            
            // save arrays for later
            xi_ik[i]=xi_ik_tmp;
            wijk[i]=tmp_weight;
        
            if(rjk_mag>R0) continue;       
            
            used_triples++; // count a triple
            
            // Now compute the integral;
            c3v = tmp_kernel*tmp_weight*xi_ik_tmp; // include symmetry factor
                
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
            // Fill up matrix
            for(int a=0;a<nbin;a++){
                for(int p=0;p<mbin;p++){
                    tmp_bin = (a*mbin+p)*nbin*mbin;
                    for(int b=0;b<nbin;b++){
                        for(int q=0;q<mbin;q++){
                            c3[tmp_bin+b*mbin+q]+=c3v*kernel_ij[i*mbin*nbin+a*mbin+p]*this_kernel[b*mbin+q];
                        }
                    }
                }
            }
        }
    }

    inline void fourth(const Particle* pi_list, const int* prim_ids, const int pln, const Particle pj, const Particle pk, const Particle pl, const int pj_id, const int pk_id, const int pl_id, const Float* wijk, const Float* xi_ik, const double prob, const Float* kernel_ij){
        // Accumulates the four point integral C4.
<<<<<<< HEAD

=======
        
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
        // First define variables
        Particle pi;
        Float rjl_mag, rjl_mu, rkl_mag, rkl_mu, c4v, xi_jl, tmp_phi_inv=0,this_kernel[nbin*mbin],tmp_weight;
        Float old_kernel[mbin], new_kr, old_kr, new_kernel, legendre_kl[max_legendre/2+1];
        int tmp_bin;
<<<<<<< HEAD

        cleanup_l(pl.pos,pk.pos,rkl_mag,rkl_mu);
        if(rkl_mag>R0) return; // if k-l separation too large

        // Load legendre polynomials
        legendre_polynomials(rkl_mu, max_legendre, legendre_kl);

        // Compute 1/Phi_kl function
        for(int l_i=0;l_i<sc34->l_bins;l_i++) tmp_phi_inv+=legendre_kl[l_i]*sc34->inv_correction_function(l_i*2,rkl_mu);

=======
        
        cleanup_l(pl.pos,pk.pos,rkl_mag,rkl_mu); 
        if(rkl_mag>R0) return; // if k-l separation too large
        
        // Load legendre polynomials
        legendre_polynomials(rkl_mu, max_legendre, legendre_kl);
        
        // Compute 1/Phi_kl function
        for(int l_i=0;l_i<sc34->l_bins;l_i++) tmp_phi_inv+=legendre_kl[l_i]*sc34->inv_correction_function(l_i*2,rkl_mu);
        
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
        // Now compute multipole matrix elements
#ifdef UNBINNED
            for(int ii=0;ii<nbin;ii++){
                old_kr = rkl_mag*0.5*(r_low[ii]+r_high[ii]);
                new_kr = sin(old_kr)/old_kr;
                for(int jj=0;jj<mbin;jj++){
                    this_kernel[ii*mbin+jj] = new_kr*legendre_kl[jj];
                }
            }
#else
        for(int ii=0;ii<nbin;ii++){
            if(ii==0) old_kr = rkl_mag*r_low[ii];
            new_kr = rkl_mag*r_high[ii];
            for(int jj=0;jj<mbin;jj++){
                if(ii==0) old_kernel[jj] = kernel_interp->kernel(jj*2.,old_kr);
                new_kernel = kernel_interp->kernel(jj*2,new_kr);
                this_kernel[ii*mbin+jj] = legendre_kl[jj]*(new_kernel - old_kernel[jj]);
                old_kernel[jj] = new_kernel;
            }
            old_kr = new_kr;
        }
<<<<<<< HEAD
#endif
        // Define 2PCF
        cleanup_l(pl.pos,pj.pos,rjl_mag,rjl_mu);
        xi_jl = cf24->xi(rjl_mag, rjl_mu); // j-l correlation

=======
#endif 
        // Define 2PCF
        cleanup_l(pl.pos,pj.pos,rjl_mag,rjl_mu); 
        xi_jl = cf24->xi(rjl_mag, rjl_mu); // j-l correlation
        
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
#ifdef UNBINNED
        tmp_weight = pl.w*pair_weight(rkl_mag)/(tmp_phi_inv*prob*pow(rkl_mag,3))*2*xi_jl;
#else
        tmp_weight = pl.w*pair_weight(rkl_mag)/(tmp_phi_inv*prob*pow(rkl_mag,3))*2*xi_jl;
<<<<<<< HEAD
#endif
        for(int i=0;i<pln;i++){ // Iterate over particle in pi_list
            if(wijk[i]==-1) continue; // skip incorrect bins / ij self counts
            if(((prim_ids[i]==pl_id)&&(I1==I4))||((pj_id==pl_id)&&(I2==I4))||((pk_id==pl_id)&&(I3==I4))) continue; // don't self-count
            pi = pi_list[i];

            used_quads++; // count a quad

            // Now compute the integral;
            c4v = wijk[i]*tmp_weight*xi_ik[i]; // with xi_ik*xi_jl = xi_il*xi_jk symmetry factor

=======
#endif   
        for(int i=0;i<pln;i++){ // Iterate over particle in pi_list
            if((wijk[i]==-1)||(pl_id==prim_ids[i])) continue;
            pi = pi_list[i];
            
            used_quads++; // count a quad
            
            // Now compute the integral;
            c4v = wijk[i]*tmp_weight*xi_ik[i]; // with xi_ik*xi_jl = xi_il*xi_jk symmetry factor
                
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
            // Fill up matrix
            for(int a=0;a<nbin;a++){
                for(int p=0;p<mbin;p++){
                    tmp_bin = (a*mbin+p)*nbin*mbin;
                    for(int b=0;b<nbin;b++){
                        for(int q=0;q<mbin;q++){
                            c4[tmp_bin+b*mbin+q]+=c4v*kernel_ij[i*mbin*nbin+a*mbin+p]*this_kernel[b*mbin+q];
                        }
                    }
                }
            }
        }
    }
<<<<<<< HEAD


=======
        
        
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
    inline Float pair_weight(Float sep){
        // Compute weight function W(r;R_0)
        if(sep<R0/2) return 1.;
        else if(sep>R0) return 0.;
        else{
            Float x = sep/R0;
            if(sep<3*R0/4) return 1.-8.*pow(2*x-1,3)+8.*pow(2*x-1,4);
            else return -64.*pow(x-1,3)-32*pow(x-1,4);
        }
    }

<<<<<<< HEAD

=======
        
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
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
        }
    }
<<<<<<< HEAD

=======
    
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
    void frobenius_difference_sum(Integrals* ints, int n_loop, Float &frobC2, Float &frobC3, Float &frobC4){
        // Add the values accumulated in ints to the corresponding internal sums and compute the Frobenius norm difference between integrals
        Float n_loops = (Float)n_loop;
        Float self_c2=0, diff_c2=0;
        Float self_c3=0, diff_c3=0;
        Float self_c4=0, diff_c4=0;
<<<<<<< HEAD

=======
        
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
        // Compute Frobenius norms and sum integrals
        for(int i=0;i<nbin*mbin;i++){
            for(int j=0;j<nbin*mbin;j++){
                self_c2+=pow(c2[i*nbin*mbin+j]/n_loops,2.);
<<<<<<< HEAD
                diff_c2+=pow(c2[i*nbin*mbin+j]/n_loops-(c2[i*nbin*mbin+j]+ints->c2[i*nbin*mbin+j])/(n_loops+1.),2.);
=======
                diff_c2+=pow(c2[i*nbin*mbin+j]/n_loops-(c2[i*nbin*mbin+j]+ints->c2[i*nbin*mbin+j])/(n_loops+1.),2.);            
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
                self_c3+=pow(c3[i*nbin*mbin+j]/n_loops,2.);
                diff_c3+=pow(c3[i*nbin*mbin+j]/n_loops-(c3[i*nbin*mbin+j]+ints->c3[i*nbin*mbin+j])/(n_loops+1.),2.);
                self_c4+=pow(c4[i*nbin*mbin+j]/n_loops,2.);
                diff_c4+=pow(c4[i*nbin*mbin+j]/n_loops-(c4[i*nbin*mbin+j]+ints->c4[i*nbin*mbin+j])/(n_loops+1.),2.);
                }
            }
<<<<<<< HEAD

=======
        
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
        self_c2=sqrt(self_c2);
        diff_c2=sqrt(diff_c2);
        diff_c3=sqrt(diff_c3);
        diff_c4=sqrt(diff_c4);
        self_c3=sqrt(self_c3);
        self_c4=sqrt(self_c4);
<<<<<<< HEAD

=======
        
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
        // Return percent difference
        frobC2=100.*(diff_c2/self_c2);
        frobC3=100.*(diff_c3/self_c3);
        frobC4=100.*(diff_c4/self_c4);
        }

    void sum_total_counts(Integrals *ints){
        // Add counts accumulated in different threads
        used_pairs+=ints->used_pairs;
        used_triples+=ints->used_triples;
        used_quads+=ints->used_quads;
    }
<<<<<<< HEAD

=======
    
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
    void normalize(Float norm1, Float norm2, Float norm3, Float norm4, Float n_pairs, Float n_triples, Float n_quads, Float norm_factor){
        // Normalize the accumulated integrals (partly done by the normalising probabilities used from the selected cubes)
        // n_pair etc. are the number of PARTICLE pairs etc. attempted (not including rejected cells, but including pairs which don't fall in correct bin ranges)
        // NB: norm_factor is V*<(nw)^2> or Sum nw^2 for the galaxies
<<<<<<< HEAD

=======
        
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
        // To avoid recomputation
        double corrf2 = norm1*norm2; // correction factor for densities of random points
        double corrf3 = corrf2*norm3;
        double corrf4 = corrf3*norm4;;
<<<<<<< HEAD

=======
        
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
        for(int i = 0; i<nbin*mbin*nbin*mbin;i++){
            c2[i]/=(n_pairs*corrf2);
            c3[i]/=(n_triples*corrf3);
            c4[i]/=(n_quads*corrf4);
        }
<<<<<<< HEAD

=======
        
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
        // Further normalize by pre-factor
        int legendre_p,legendre_q;
        Float normalization,k_norm;
        for(int i=0; i<nbin*mbin;i++){
            legendre_p = (i%mbin)*2; // first legendre index
            for(int j=0;j<nbin*mbin;j++){
                legendre_q = (j%mbin)*2.; // second legendre index
<<<<<<< HEAD

=======
                
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
                // Correction for bin sizes
#ifdef UNBINNED
                k_norm = 1;
#else
                k_norm = (pow(r_high[i],3)-pow(r_low[i],3))*(pow(r_high[j],3)-pow(r_low[j],3));
<<<<<<< HEAD
#endif
                normalization = float((2*legendre_p+1)*(2*legendre_q+1))/pow(norm_factor,2);
=======
#endif           
                normalization = float((2*legendre_p+1)*(2*legendre_q+1))/pow(norm_factor,2);                
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
                c2[i*nbin*mbin+j]*=normalization/k_norm;
                c3[i*nbin*mbin+j]*=normalization/k_norm;
                c4[i*nbin*mbin+j]*=normalization/k_norm;
            }
        }
    }
<<<<<<< HEAD

=======
    
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
    void save_counts(uint64 pair_counts,uint64 triple_counts,uint64 quad_counts){
        // Print the counts for each integral (used for combining the estimates outside of C++)
        // This is the number of counts used in each loop [always the same]
        char counts_file[1000];
        snprintf(counts_file, sizeof counts_file, "%sPowerCovMatrices/total_counts_n%d_l%d_%d%d,%d%d.txt",out_file,nbin,(mbin-1)*2,I1,I2,I3,I4);
        FILE * CountsFile = fopen(counts_file,"w");
        fprintf(CountsFile,"%llu\n",pair_counts);
        fprintf(CountsFile,"%llu\n",triple_counts);
        fprintf(CountsFile,"%llu\n",quad_counts);
<<<<<<< HEAD

        fflush(NULL);
        fclose(CountsFile);
    }


    void save_integrals(char* suffix, bool save_all) {
    /* Print integral outputs to file.
        * In txt files {c2,c3,c4}_leg_n{nbin}_m{mbin}.txt there are lists of the outputs of c2,c3,c4 that are already normalized and multiplied by combinatoric factors. The n and m strings specify the number of n and m bins present.
        */
        // Create output files

=======
        
        fflush(NULL);
        fclose(CountsFile);
    }
        
    
    void save_integrals(char* suffix, bool save_all) {
    /* Print integral outputs to file. 
        * In txt files {c2,c3,c4}_leg_n{nbin}_m{mbin}.txt there are lists of the outputs of c2,c3,c4 that are already normalized and multiplied by combinatoric factors. The n and m strings specify the number of n and m bins present.
        */
        // Create output files
        
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
        char c2name[1000];
        snprintf(c2name, sizeof c2name, "%sPowerCovMatrices/c2_leg_n%d_l%d_%d%d_%s.txt", out_file,nbin, (mbin-1)*2,I1,I2,suffix);
        char c3name[1000];
        snprintf(c3name, sizeof c3name, "%sPowerCovMatrices/c3_leg_n%d_l%d_%d,%d%d_%s.txt", out_file, nbin, (mbin-1)*2,I2,I1,I3,suffix);
        char c4name[1000];
        snprintf(c4name, sizeof c4name, "%sPowerCovMatrices/c4_leg_n%d_l%d_%d%d,%d%d_%s.txt", out_file, nbin, (mbin-1)*2, I1,I2,I3,I4,suffix);
        FILE * C2File = fopen(c2name,"w"); // for c2 part of integral
        FILE * C3File = fopen(c3name,"w"); // for c3 part of integral
        FILE * C4File = fopen(c4name,"w"); // for c4 part of integral
<<<<<<< HEAD

=======
        
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
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
<<<<<<< HEAD

        fflush(NULL);

=======
        
        fflush(NULL);
        
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
        // Close open files
        fclose(C2File);
        fclose(C3File);
        fclose(C4File);
<<<<<<< HEAD

=======
        
>>>>>>> 1822107b9340b8e7d0ce82303f62cf7b35d03f8c
    }

};

#endif
