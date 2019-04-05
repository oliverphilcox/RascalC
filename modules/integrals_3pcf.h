// Rewritten integrals.h code for legendre multipoles
#include "parameters.h"
#include "correlation_function.h"
#include "cell_utilities.h"
#include "legendre_utilities.h"

#ifndef INTEGRALS_3PCF_H
#define INTEGRALS_3PCF_H

class Integrals{
public:
    CorrelationFunction *cf;
    
private:
    int nbin, mbin, max_l, array_len;
    Float rmin,rmax,mumin,mumax; //Ranges in r and mu
    Float *r_high, *r_low; // Max and min of each radial bin
    Float *c3, *c4, *c5, *c6; // Arrays to accumulate integrals
    char* out_file;
    bool box; // Flags to decide whether we have a periodic box
    
    uint64 *binct3, *binct4, *binct5, *binct6; // Arrays to accumulate bin counts
    SurveyCorrection *sc; // survey correction function
    
public:
    Integrals(){};

    Integrals(Parameters *par, CorrelationFunction *_cf, SurveyCorrection *_sc){
        sc = _sc;
        max_l = par->max_l;
        cf = new CorrelationFunction(_cf);
        init(par);
    }
    
    void init(Parameters *par){
        nbin = par->nbin; // number of radial bins
        mbin = par->mbin; // number of mu bins
        out_file = par->out_file; // output directory
        max_l = par->max_l; // maximum Legendre polynomial
        
        int ec=0;
        array_len = pow(nbin,2)*mbin;
        // Initialize the binning
        ec+=posix_memalign((void **) &c3, PAGE, sizeof(double)*array_len*array_len);
        ec+=posix_memalign((void **) &c4, PAGE, sizeof(double)*array_len*array_len);
        ec+=posix_memalign((void **) &c5, PAGE, sizeof(double)*array_len*array_len);
        ec+=posix_memalign((void **) &c6, PAGE, sizeof(double)*array_len*array_len);

        ec+=posix_memalign((void **) &binct3, PAGE, sizeof(uint64)*array_len*array_len);
        ec+=posix_memalign((void **) &binct4, PAGE, sizeof(uint64)*array_len*array_len);
        ec+=posix_memalign((void **) &binct5, PAGE, sizeof(uint64)*array_len*array_len);
        ec+=posix_memalign((void **) &binct6, PAGE, sizeof(uint64)*array_len*array_len);
        
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
        free(c3);
        free(c4);
        free(c5);
        free(c6);
        free(binct3);
        free(binct4);
        free(binct5);
        free(binct6);
    }

    void reset(){
        for (int j=0; j<array_len*array_len; j++) {
            c3[j]=0;
            c4[j]=0;
            c5[j]=0;
            c6[j]=0;
            binct3[j] = 0;
            binct4[j] = 0;
            binct5[j] = 0;
            binct6[j] = 0;
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
    
    Float compute_los(Float3 p1, Float3 p2, Float norm){
        // Compute line of sight angle for two vectors p1, p2
        Float3 pos=p1-p2;
#ifndef PERIODIC
        Float3 los=p1+p2; // No 1/2 as normalized anyway below
        return fabs(pos.dot(los)/norm/los.norm());
#else
        // In the periodic case use z-direction for mu
        return fabs(pos.z/norm);
#endif
        }
        
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
        
    void triangle_bins(Float3 p1, Float3 p2, Float3 p3, Float* norm, Float* ang,Float norm_in, int preload){
        // Compute side lengths and internal angles of a triangle of positions
        // Indices are chosen such that norm[x] is opposite to ang[x]
        // NB we assume norm(p2,p3) [if preload==0] or norm(1,2) [if preload==2] is already computed here for speed.
        
        Float3 pos12 = p1-p2, pos23 = p2-p3, pos13 = p1-p3;
        
        // Compute lengths
        if(preload==0) norm[0] = norm_in;
        else norm[0]=pos23.norm();
        norm[1] = pos13.norm();
        if(preload==2) norm[2] = norm_in;
        else norm[2] = pos12.norm();
        
        // Compute internal angles
        ang[0] = fabs(pos12.dot(pos13)/norm[2]/norm[1]);
        ang[1] = fabs(pos12.dot(pos23)/norm[2]/norm[0]);
        ang[2] = fabs(pos13.dot(pos23)/norm[1]/norm[0]);
    }
    
    int all_bins(Float norm[3], Float bin_ab[6], const int bin_in,int preload){
        // Compute all triangle bin combinations and save as an array. We use the j-k bin as an input here
        
        int error=0; // to check if there are any correct bins
        Float bin_ij=0, bin_ik=0, bin_jk=0;
        
        // Load bins and check if correct bins are sampled
        if(preload==0) bin_jk=bin_in;
        else bin_jk = get_radial_bin(norm[0]);
        bin_ik = get_radial_bin(norm[1]);     
        if(preload==2) bin_ij=bin_in;
        else bin_ij = get_radial_bin(norm[2]);
        
        if ((bin_ij<0)||(bin_ij>=nbin)){
            bin_ij = -1;
            error++;
        }
        if((bin_jk<0)||(bin_jk>=nbin)){
            bin_jk = -1;
            error++;
        }
        if ((bin_ik<0)||(bin_ik>=nbin)){
            bin_ik = -1;
            error++;
        }
        
        // Now save all 3 combinations of these bins
        // Set 1
        bin_ab[0] = bin_ij;
        bin_ab[1] = bin_ik;
        // Set 2
        bin_ab[2] = bin_ij;
        bin_ab[3] = bin_jk;
        // Set 3
        bin_ab[4] = bin_ik;
        bin_ab[5] = bin_jk;
        
        return error;
    }
    
    inline void third(const Particle* pi_list, const int* prim_ids, int pln, const Particle pj, const Particle pk, const int pj_id, const int pk_id, const double prob, Float* &wijk, int* &all_bins_ijk, Float* &all_correction_factor, Float* &all_legendre, Float* &xi_pass, Float &norm_jk, int &bin_jk, int index){
        // Accumulates the three point integral C3. Also outputs several arrays for later reuse.
        // Prob. here is defined as g_ij / f_ij where g_ij is the sampling PDF and f_ij is the true data PDF for picking sets of particles
        
        Float tmp_weight, tmp_xi=0, this_poly, los_tmp,c3v;
        Float norm_ijk[3], ang_ijk[3], bins_ijk[6];
        Particle pi;
        Float3 pjk;
        int tmp_radial_bin, tmp_radial_bin2, out_bin, bin_1, bin_2,bin_3,bin_4;
        Float polynomials_tmp[mbin],correction_factor1,correction_factor2;
        
        // Check to avoid self-counts
        if((pj_id==pk_id)){
            wijk[0]=-2; // special value to be checked later
            return;
        }
        
        // Define j-k norm and bin
        pjk = pj.pos-pk.pos;
        norm_jk = pjk.norm();
        bin_jk = get_radial_bin(norm_jk);
            
        // Compute xi_jk if needed
        if(index==1){
            los_tmp = compute_los(pj.pos,pk.pos,norm_jk);
            tmp_xi = cf->xi(norm_jk,los_tmp);
            xi_pass[0] = tmp_xi; // save for next integrator
        }
        
        // No self-counts here
        for(int i=0;i<pln;i++){ // Iterate over particle in pi_list
            
            // First check for any self-counts
            if((prim_ids[i]==pj_id)||(prim_ids[i]==pk_id)){
                wijk[i]=-1;
                continue;
            }
            
            pi = pi_list[i]; // first particle
            
            // Define angles and separations
            triangle_bins(pi.pos,pj.pos,pk.pos,norm_ijk,ang_ijk,norm_jk,0);
            
            // Define radial bins
            int x = all_bins(norm_ijk,bins_ijk,bin_jk,0);
            if(x==3){
                wijk[i]=-1;
                continue; // skip if no correct bins
            }
            
            // Define weights and desired correlation function
            tmp_weight = pi.w*pj.w*pk.w; 
            // Save into arrays for later
            wijk[i] = tmp_weight; 
            
            if(index==0){
                los_tmp = compute_los(pi.pos, pj.pos, norm_ijk[2]);
                tmp_xi = cf->xi(norm_ijk[2], los_tmp); //xi_ij
                xi_pass[i] = tmp_xi; // save for next integrator
            }
            
            // Now compute the integral contribution
            c3v = 3.*tmp_weight*tmp_weight*(1.+3.*tmp_xi) / prob;
            
            // Pre-load relevant Legendre polynomials and correction factors;
            for(int bin_index=0;bin_index<3;bin_index++){
                bin_1 = bins_ijk[bin_index*2];
                if(bin_1==-1){
                    all_correction_factor[i*3+bin_index]=-1;
                    continue;
                }
                bin_2 = bins_ijk[bin_index*2+1];
                if(bin_2==-1){
                    all_correction_factor[i*3+bin_index]=-1;
                    continue; // skip if bad bin
                }
                correction_factor1 = sc->correction_function_3pcf(bin_1,bin_2,ang_ijk[bin_index]);
                all_correction_factor[i*3+bin_index] = correction_factor1;
                legendre_polynomials(ang_ijk[bin_index],max_l,polynomials_tmp);
                for(int p_bin=0;p_bin<mbin;p_bin++){
                    all_legendre[(i*3+bin_index)*mbin+p_bin]=polynomials_tmp[p_bin];
                }
            }
            
            for(int bin_index=0;bin_index<3;bin_index++){
                // Load correction factor
                correction_factor1 = all_correction_factor[i*3+bin_index];
                if(correction_factor1==-1) continue; // skip if bad bin combination
                
                bin_1 = bins_ijk[bin_index*2];
                bin_2 = bins_ijk[bin_index*2+1];
                
                // Save bins for later
                all_bins_ijk[i*6+bin_index*2]=bin_1;
                all_bins_ijk[i*6+bin_index*2+1]=bin_2;
                
                tmp_radial_bin = (bin_1*nbin+bin_2)*mbin;
                
                for(int bin_index2=0;bin_index2<3;bin_index2++){
                    // Load correction factor
                    correction_factor2 = all_correction_factor[i*3+bin_index2];
                    if(correction_factor2==-1) continue; // skip if bad bin combination
                    
                    bin_3 = bins_ijk[bin_index2*2];
                    bin_4 = bins_ijk[bin_index2*2+1];
                
                    // Add to relevant Legendre bins
                    tmp_radial_bin2 = (bin_3*nbin+bin_4)*mbin;
                    for(int p_bin=0;p_bin<mbin;p_bin++){
                        this_poly = all_legendre[(i*3+bin_index)+p_bin];
                        for(int q_bin=0;q_bin<mbin;q_bin++){
                            
                            out_bin = (tmp_radial_bin+p_bin)*array_len+tmp_radial_bin2+q_bin;
                        
                            // Add to integral (extra 4 is from 2x symmetry in each angle kernel)
                            c3[out_bin]+=c3v*this_poly*all_legendre[(i*3+bin_index2)+q_bin]*correction_factor1*correction_factor2;
                            binct3[out_bin]++;
                        }
                    }
                }
            }
        }
    }
    
    inline void fourth(const Particle* pi_list, const int* prim_ids, const int pln, const Particle pj, const Particle pk, const Particle pl, const int pj_id, const int pk_id, const int pl_id, const double prob, const Float* wijk, const int* bins_ijk, const Float* all_correction_factor_ijk, const Float* all_legendre_ijk, const Float* xi_pass, const Float norm_jk, const int bin_jk, Float* &wijkl, Float* &xi_pass2, Float &norm_kl, int &bin_kl, int index){
        // Accumulates the four point integral C4. Also outputs an array of xi_ik and bin_ik values for later reuse.
        
        Float norm_jkl[3],ang_jkl[3], bins_jkl[6], tmp_weight, tmp_xi1,tmp_xi2,c4v,this_poly,los_tmp;
        Particle pi;
        int bin_1,bin_2,bin_3,bin_4,out_bin, tmp_radial_bin, tmp_radial_bin2;
        Float polynomials_tmp1[mbin], correction_factor1,correction_factor2;
        Float all_correction_factor_jkl[3],all_legendre_jkl[3*mbin];
        
        // Check to avoid self-counts
        if((pl_id==pk_id)||(pl_id==pj_id)||(wijk[0]==-2)){
            wijkl[0]=-2;
            return;
        }
        
        // Define triangle sides independent of i
        triangle_bins(pj.pos,pk.pos,pl.pos,norm_jkl,ang_jkl,norm_jk,2);
        
        // Define first xi function
        if(index==0){
            los_tmp = compute_los(pk.pos,pl.pos,norm_jkl[0]);
            tmp_xi1 = cf->xi(norm_jkl[0],los_tmp); //xi_kl
            xi_pass2[0] = tmp_xi1;
        }
        else tmp_xi1 = xi_pass[0]; //xi_jk
        
        // Save lengths and bins for later
        norm_kl = norm_jkl[0];
        
        // Compute radial bins for this triangle
        int x = all_bins(norm_jkl,bins_jkl,bin_jk,2);
        if(x==3){
            wijkl[0]=-2;
            return; // if no correct radial bins
        }
        bin_kl = bins_jkl[5];
        
        // Preload correction factors and Legendre polynomials
        for(int bin_index=0;bin_index<3;bin_index++){
            bin_1 = bins_jkl[bin_index*2];
            if(bin_1==-1){ // if incorrect bin
                all_correction_factor_jkl[bin_index]=-1;
                continue;
            }
            bin_2 = bins_jkl[bin_index*2+1];
            if(bin_2==-1){
                all_correction_factor_jkl[bin_index]=-1;
                continue;
            }
            all_correction_factor_jkl[bin_index]=sc->correction_function_3pcf(bin_1,bin_2,ang_jkl[bin_index]);
            legendre_polynomials(ang_jkl[bin_index],max_l,polynomials_tmp1);
            for(int p_bin=0;p_bin<mbin;p_bin++) all_legendre_jkl[bin_index*mbin+p_bin]=polynomials_tmp1[p_bin];
        }
        
        for(int i=0;i<pln;i++){ // Iterate over particle in pi_list
            
            // First check for any self-counts;
            if((prim_ids[i]==pl_id)||(wijk[i]==-1)){
                wijkl[i]=-1;
                continue;
            }
            pi = pi_list[i]; // first particle
            
            // Define weights and correlation function
            tmp_weight = wijk[i]*pl.w; // 
            wijkl[i] = tmp_weight; // save for later
            
            if(index==0) tmp_xi2 = xi_pass[i]; // xi_ij
            else{
                // Compute xi_il
                Float norm_tmp;
                cleanup_l(pi.pos,pl.pos,norm_tmp,los_tmp);
                tmp_xi2 = cf->xi(norm_tmp,los_tmp); // xi_il
                xi_pass2[i] = tmp_xi2; // save for later
            }
            
            // Now compute integral contribution
            if(index==0) c4v = 36.*tmp_weight*pj.w*pk.w/prob * tmp_xi1 * tmp_xi2;
            else c4v = 18.*tmp_weight*pj.w*pk.w/prob * tmp_xi2*(tmp_xi1+1.);
            
            for(int bin_index=0;bin_index<3;bin_index++){
                correction_factor1 = all_correction_factor_ijk[i*3+bin_index];
                if(correction_factor1==-1) continue;
                bin_1 = bins_ijk[bin_index*2];
                bin_2 = bins_ijk[bin_index*2+1];
                
                // Load correction factor
                tmp_radial_bin = (bin_1*nbin+bin_2)*mbin;
                
                for(int bin_index2=0;bin_index2<3;bin_index2++){
                    correction_factor2 = all_correction_factor_jkl[bin_index2];
                    if(correction_factor2==-1) continue;
                    bin_3 = bins_jkl[bin_index2*2];
                    bin_4 = bins_jkl[bin_index2*2+1];
                    
                    // Add to relevant Legendre bins
                    tmp_radial_bin2 = (bin_3*nbin+bin_4)*mbin;
                    for(int p_bin=0;p_bin<mbin;p_bin++){
                        this_poly = all_legendre_ijk[(i*3+bin_index)*mbin+p_bin];
                        for(int q_bin=0;q_bin<mbin;q_bin++){
                            out_bin = (tmp_radial_bin+p_bin)*array_len+tmp_radial_bin2+q_bin;
                        
                            // Add to integral (extra 4 is from 2x symmetry in each angle kernel)
                            c4[out_bin] +=c4v*this_poly*all_legendre_jkl[bin_index2*mbin+q_bin]*correction_factor1*correction_factor2*4.;
                            binct4[out_bin]++;
                        }
                    }
                }
            }
        }
    }

    inline void fifth(const Particle* pi_list, const int* prim_ids, const int pln, const Particle pj, const Particle pk, const Particle pl, Particle pm, const int pj_id, const int pk_id, const int pl_id, const int pm_id, const double prob, const Float* wijkl, const int* bins_ijk, const Float* all_correction_factor_ijk, const Float* all_legendre_ijk, const Float* xi_pass, const Float* xi_pass2, const Float norm_kl, const int bin_kl, Float* &wijklm, Float &xi_pass3, Float &norm_lm, int &bin_lm, int index){
        // Accumulates the five point integral C5.
        
        Float norm_klm[3],ang_klm[3],bins_klm[6],los_tmp, tmp_xi1,tmp_xi2, polynomials_tmp[mbin];
        Float all_correction_factor_klm[3],all_legendre_klm[3*mbin],tmp_weight,c5v, this_poly, correction_factor1, correction_factor2;
        int bin_1,bin_2,bin_3,bin_4,tmp_radial_bin, tmp_radial_bin2, out_bin;
        Particle pi;
        
        // Check to avoid self-counts
        if((pm_id==pl_id)||(pm_id==pk_id)||(pm_id==pj_id)||(wijkl[0]==-2)){
            wijklm[0]=-2;
            return;
        }
        
        // Define triangle sides independent of i
        triangle_bins(pk.pos,pl.pos,pm.pos,norm_klm,ang_klm,norm_kl,2);
        
        norm_lm = norm_klm[0];
        
        // Define first xi function
        if(index==0){
            los_tmp = compute_los(pl.pos,pm.pos,norm_klm[0]);
            tmp_xi1 = cf->xi(norm_klm[0],los_tmp); // xi_lm
        }
        else{
            Float norm_tmp;
            cleanup_l(pj.pos,pm.pos,norm_tmp,los_tmp);
            tmp_xi1 = cf->xi(norm_tmp,los_tmp); // xi_jm
            xi_pass3 = tmp_xi1; // save for next integrator
        }
        
        // Compute radial bins for this triangle
        int x = all_bins(norm_klm,bins_klm,bin_kl,2);
        bin_lm = bins_klm[0];
        if(x==3){
            wijklm[0]=-2;
            return; // if no correct radial bins
        } 
        
        // Preload correction factors and Legendre polynomials
        for(int bin_index=0;bin_index<3;bin_index++){
            bin_1 = bins_klm[bin_index*2];
            if(bin_1==-1){
                all_correction_factor_klm[bin_index]=-1;
                continue;
            }
            bin_2 = bins_klm[bin_index*2+1];
            if(bin_2 ==-1){
                all_correction_factor_klm[bin_index]=-1;
                continue;
            }
            all_correction_factor_klm[bin_index]=sc->correction_function_3pcf(bin_1,bin_2,ang_klm[bin_index]);
            legendre_polynomials(ang_klm[bin_index],max_l,polynomials_tmp);
            for(int p_bin=0;p_bin<mbin;p_bin++) all_legendre_klm[bin_index*mbin+p_bin]=polynomials_tmp[p_bin];
        }
        
        for(int i=0;i<pln;i++){ // Iterate over particle in pi_list
            
            // First check for any self-counts;
            if((prim_ids[i]==pm_id)||(wijkl[i]==-1)){
                wijklm[i]=-1;
                continue;
            }
            pi = pi_list[i]; // first particle
            
            // Define weights and correlation function
            tmp_weight = wijkl[i]*pm.w; // 
            wijklm[i] = tmp_weight; // save for later
            
            //TODO: maybe rearrange xi-pass to get less smaller arrays?
            if(index==0) tmp_xi2 = xi_pass[i]; // xi_ij
            else tmp_xi2 = xi_pass2[i]; // xi_il
            
            // Now compute integral contribution
            if(index==0) c5v = 9.*tmp_weight*pk.w/prob * tmp_xi1 * tmp_xi2;
            else c5v = 18.*tmp_weight*pk.w/prob * tmp_xi1 * tmp_xi2;
            
            for(int bin_index=0;bin_index<3;bin_index++){
                // Load correction factor
                correction_factor1 = all_correction_factor_ijk[i*3+bin_index];
                if(correction_factor1==-1) continue;
                
                bin_1 = bins_ijk[bin_index*2];
                bin_2 = bins_ijk[bin_index*2+1];
                tmp_radial_bin = (bin_1*nbin+bin_2)*mbin;
                
                for(int bin_index2=0;bin_index2<3;bin_index2++){
                    correction_factor2 = all_correction_factor_klm[bin_index2];
                    if(correction_factor2==-1) continue;
           
                    bin_3 = bins_klm[bin_index2*2];
                    bin_4 = bins_klm[bin_index2*2+1];
                    
                    // Add to relevant Legendre bins
                    tmp_radial_bin2 = (bin_3*nbin+bin_4)*mbin;
                    for(int p_bin=0;p_bin<mbin;p_bin++){
                        this_poly = all_legendre_ijk[(i*3+bin_index)*mbin+p_bin];
                        for(int q_bin=0;q_bin<mbin;q_bin++){
                            out_bin = (tmp_radial_bin+p_bin)*array_len+tmp_radial_bin2+q_bin;
                        
                            // Add to integral (extra 4 is from 2x symmetry in each angle kernel)
                            c5[out_bin] += c5v*this_poly*all_legendre_klm[bin_index2*mbin+q_bin]*correction_factor1*correction_factor2*4.;
                            binct5[out_bin]++;
                        }
                    }
                }
            }
        }
    }
    
    inline void sixth(const Particle* pi_list, const int* prim_ids, const int pln, const Particle pj, const Particle pk, const Particle pl, const Particle pm, const Particle pn, const int pj_id, const int pk_id, const int pl_id, const int pm_id, const int pn_id, const double prob, const Float* wijklm, const int* bins_ijk,  const Float* all_correction_factor_ijk, const Float* all_legendre_ijk, const Float* xi_pass, const Float* xi_pass2, const Float xi_pass3, const Float norm_lm, const int bin_lm, int index){
        // Accumulates the six point integral C6.
         
        Float norm_lmn[3],ang_lmn[3],bins_lmn[6],norm_tmp,los_tmp, tmp_xi1,tmp_xi2,tmp_xi3, polynomials_tmp[mbin],correction_factor1,correction_factor2,this_poly;
        Float all_correction_factor_lmn[3],all_legendre_lmn[3*mbin],c6v;
        int bin_1,bin_2,bin_3,bin_4,tmp_radial_bin,tmp_radial_bin2,out_bin;
        Particle pi;
        
        // Check to avoid self-counts
        if((pn_id==pm_id)||(pn_id==pl_id)||(pn_id==pk_id)||(pn_id==pj_id)||(wijklm[0]==-2)) return;
        
        // Define triangle sides independent of i
        triangle_bins(pl.pos,pm.pos,pn.pos,norm_lmn,ang_lmn,norm_lm,2);
        
        // Compute radial bins for this triangle
        int x = all_bins(norm_lmn,bins_lmn,bin_lm,2);
        if(x==3) return; // if no correct radial bins
        
        // Define first two xi functions
        if(index==0){
            cleanup_l(pm.pos,pn.pos,norm_tmp,los_tmp);
            tmp_xi1 = cf->xi(norm_tmp,los_tmp); // xi_mn
            tmp_xi2 = xi_pass2[0]; // xi_kl
        }
        else{
            tmp_xi1 = xi_pass3; // xi_jm
            cleanup_l(pk.pos,pn.pos,norm_tmp,los_tmp);
            tmp_xi2 = cf->xi(norm_tmp,los_tmp); // xi_kn
        }
        
        // Preload correction factors and Legendre polynomials
        for(int bin_index=0;bin_index<3;bin_index++){
            bin_1 = bins_lmn[bin_index*2];
            if(bin_1==-1){
                all_correction_factor_lmn[bin_index]=-1;
                continue;
            }
            bin_2 = bins_lmn[bin_index*2+1];
            if(bin_2==-1){
                all_correction_factor_lmn[bin_index]=-1;
                continue;
            }
            all_correction_factor_lmn[bin_index]=sc->correction_function_3pcf(bin_1,bin_2,ang_lmn[bin_index]);
            legendre_polynomials(ang_lmn[bin_index],max_l,polynomials_tmp);
            for(int p_bin=0;p_bin<mbin;p_bin++) all_legendre_lmn[bin_index*mbin+p_bin]=polynomials_tmp[p_bin];
        }
        
        for(int i=0;i<pln;i++){ // Iterate over particle in pi_list
            
            // First check for any self-counts;
            if((prim_ids[i]==pm_id)||(wijklm[i]==-1)){
                continue;
            }
            pi = pi_list[i]; // first particle
            
            // Define correlation functions
            if(index==0) tmp_xi3 = xi_pass[i]; // xi_ij
            else tmp_xi3 = xi_pass2[i]; // xi_il
            
            // Now compute integral contribution
            if(index==0) c6v = 9.*wijklm[i]*pn.w/prob * tmp_xi1 * tmp_xi2 * tmp_xi3;
            else c6v = 6.*wijklm[i]*pn.w/prob * tmp_xi1 * tmp_xi2 * tmp_xi3;
            
            for(int bin_index=0;bin_index<3;bin_index++){
                // Load correction factor
                correction_factor1 = all_correction_factor_ijk[i*3+bin_index];
                if(correction_factor1==-1) continue;
                
                bin_1 = bins_ijk[bin_index*2];
                bin_2 = bins_ijk[bin_index*2+1];
                tmp_radial_bin = (bin_1*nbin+bin_2)*mbin;
                
                for(int bin_index2=0;bin_index2<3;bin_index2++){
                    correction_factor2 = all_correction_factor_lmn[bin_index2];
                    if(correction_factor2==-1) continue;
           
                    bin_3 = bins_lmn[bin_index2*2];
                    bin_4 = bins_lmn[bin_index2*2+1];
                
                    // Add to relevant Legendre bins
                    tmp_radial_bin2 = (bin_3*nbin+bin_4)*mbin;
                    for(int p_bin=0;p_bin<mbin;p_bin++){
                        this_poly = all_legendre_ijk[(i*3+bin_index)*mbin+p_bin];
                        for(int q_bin=0;q_bin<mbin;q_bin++){
                            out_bin = (tmp_radial_bin+p_bin)*array_len+tmp_radial_bin2+q_bin;
                        
                            // Add to integral (extra 4 is from 2x symmetry in each angle kernel)
                            c6[out_bin] +=c6v*this_poly*all_legendre_lmn[bin_index2*mbin+q_bin]*correction_factor1*correction_factor2*4.;
                            binct6[out_bin]++;
                        }
                    }
                }
            }
        }
    }

public:
    void sum_ints(Integrals* ints) {
        // Add the values accumulated in ints to the corresponding internal sums
        for(int i=0;i<pow(array_len,2);i++){
            c3[i]+=ints->c3[i];
            c4[i]+=ints->c4[i];
            c5[i]+=ints->c5[i];
            c6[i]+=ints->c6[i];
            binct3[i]+=ints->binct3[i];
            binct4[i]+=ints->binct4[i];
            binct5[i]+=ints->binct5[i];
            binct6[i]+=ints->binct6[i];
        }
    }
    
    void frobenius_difference_sum(Integrals* ints, int n_loop, Float &frobC3, Float &frobC4, Float &frobC5, Float &frobC6){
        // Add the values accumulated in ints to the corresponding internal sums and compute the Frobenius norm difference between integrals
        Float n_loops = (Float)n_loop;
        Float self_c3=0, diff_c3=0;
        Float self_c4=0, diff_c4=0;
        Float self_c5=0, diff_c5=0;
        Float self_c6=0, diff_c6=0;
        
        // Compute Frobenius norms and sum integrals
        for(int i=0;i<array_len;i++){
            for(int j=0;j<array_len;j++){
                self_c3+=pow(c3[i*array_len+j]/n_loops,2.);
                diff_c3+=pow(c3[i*array_len+j]/n_loops-(c3[i*array_len+j]+ints->c3[i*array_len+j])/(n_loops+1.),2.);
                self_c4+=pow(c4[i*array_len+j]/n_loops,2.);
                diff_c4+=pow(c4[i*array_len+j]/n_loops-(c4[i*array_len+j]+ints->c4[i*array_len+j])/(n_loops+1.),2.);
                self_c5+=pow(c5[i*array_len+j]/n_loops,2.);
                diff_c5+=pow(c5[i*array_len+j]/n_loops-(c5[i*array_len+j]+ints->c5[i*array_len+j])/(n_loops+1.),2.);
                self_c6+=pow(c6[i*array_len+j]/n_loops,2.);
                diff_c6+=pow(c6[i*array_len+j]/n_loops-(c6[i*array_len+j]+ints->c6[i*array_len+j])/(n_loops+1.),2.);
            }
            }
        
        diff_c3=sqrt(diff_c3);
        diff_c4=sqrt(diff_c4);
        self_c3=sqrt(self_c3);
        self_c4=sqrt(self_c4);
        self_c5=sqrt(self_c5);
        diff_c5=sqrt(diff_c5);
        self_c6=sqrt(self_c6);
        diff_c6=sqrt(diff_c6);
        
        // Return percent difference
        frobC3=100.*(diff_c3/self_c3);
        frobC4=100.*(diff_c4/self_c4);
        frobC5=100.*(diff_c5/self_c5);
        frobC6=100.*(diff_c6/self_c6);
    }

    void sum_total_counts(uint64& acc3, uint64& acc4, uint64& acc5, uint64& acc6){
        // Add local counts to total bin counts in acc2-4
        // Divide by 9*mbin*mbin since we add every set of particles to this many bins
        for (int i=0; i<pow(array_len,2); i++) {
            acc3+=binct3[i]/(9.*mbin*mbin);
            acc4+=binct4[i]/(9.*mbin*mbin);
            acc5+=binct5[i]/(9.*mbin*mbin);
            acc6+=binct6[i]/(9.*mbin*mbin);
        }
    }
    
    void normalize(Float norm, Float n_triples, Float n_quads, Float n_quints, Float n_hexes){
        // Normalize the accumulated integrals (partly done by the normalising probabilities used from the selected cubes)
        // n_pair etc. are the number of PARTICLE pairs etc. attempted (not including rejected cells, but including pairs which don't fall in correct bin ranges)
        // To avoid recomputation
        double corrf2 = norm*norm; // correction factor for densities of random points
        double corrf3 = corrf2*norm;
        double corrf4 = corrf3*norm;
        double corrf5 = corrf4*norm;
        double corrf6 = corrf5*norm;
        
        for(int i = 0; i<pow(array_len,2);i++){
            c3[i]/=(n_triples*corrf3);
            c4[i]/=(n_quads*corrf4);
            c5[i]/=(n_quints*corrf5);
            c6[i]/=(n_hexes*corrf6);
        }
        
        // Further normalize by pre-factor
        int legendre_p,legendre_q,ind1,ind2;
        Float r_a,r_b,r_c,r_d,delta_a,delta_b,delta_c,delta_d,normalization;
        
        for(int i=0; i<array_len;i++){
            legendre_p = (i%mbin)*2; // first legendre index
            ind1 = i/mbin;
            r_a = 0.5*(r_low[ind1/nbin]+r_high[ind1/nbin]); // mid-points of bin
            r_b = 0.5*(r_low[ind1%nbin]+r_high[ind1%nbin]);
            delta_a = r_high[ind1/nbin]-r_low[ind1/nbin]; // width of bin
            delta_b = r_high[ind1%nbin]-r_low[ind1%nbin];
            
            for(int j=0;j<array_len;j++){
                legendre_q = (j%mbin)*2.; // second legendre index
                ind2 = j/mbin;
                r_c = 0.5*(r_low[ind2/nbin]+r_high[ind2/nbin]); // mid-point of bin
                r_d = 0.5*(r_low[ind2%nbin]+r_high[ind2%nbin]);
                delta_c = r_high[ind2/nbin]-r_low[ind2/nbin];  // width of bin
                delta_d = r_high[ind2%nbin]-r_low[ind2%nbin];
                
                normalization = float((2*legendre_p+1)*(2*legendre_q+1))/(pow((r_a*r_b*r_c*r_d),2)*delta_a*delta_b*delta_c*delta_d);                
                c3[i*array_len+j]*=normalization;
                c4[i*array_len+j]*=normalization;
                c5[i*array_len+j]*=normalization;
                c6[i*array_len+j]*=normalization;
            }
        }
    }
    
    void save_counts(uint64 triple_counts,uint64 quad_counts, uint64 quint_counts, uint64 hex_counts, int index){
        // Print the counts for each integral (used for combining the estimates outside of C++)
        // This is the number of counts used in each loop [always the same]
        char counts_file[1000];
        snprintf(counts_file, sizeof counts_file, "%s3PCFCovMatricesAll/total_counts_n%d_m%d_%d.txt",out_file,nbin,mbin,index);
        FILE * CountsFile = fopen(counts_file,"w");
        fprintf(CountsFile,"%llu\n",triple_counts);
        fprintf(CountsFile,"%llu\n",quad_counts);
        fprintf(CountsFile,"%llu\n",quint_counts);
        fprintf(CountsFile,"%llu\n",hex_counts);
        
        fflush(NULL);
        fclose(CountsFile);
    }
        
    
    void save_integrals(char* suffix, bool save_all,int index) {
    /* Print integral outputs to file. 
        * In txt files {c2,c3,c4}_leg_n{nbin}_m{mbin}.txt there are lists of the outputs of c2,c3,c4 that are already normalized and multiplied by combinatoric factors. The n and m strings specify the number of n and m bins present.
        */
        // Create output files
        
        char c3name[1000];
        snprintf(c3name, sizeof c3name, "%s3PCFCovMatricesAll/c3_n%d_l%d_%d_%s.txt", out_file, nbin, max_l,index,suffix);
        char c4name[1000];
        snprintf(c4name, sizeof c4name, "%s3PCFCovMatricesAll/c4_n%d_l%d_%d_%s.txt", out_file, nbin, max_l, index, suffix);
        char c5name[1000];
        snprintf(c5name, sizeof c5name, "%s3PCFCovMatricesAll/c5_n%d_l%d_%d_%s.txt", out_file,nbin, max_l,index, suffix);
        char c6name[1000];
        snprintf(c6name, sizeof c6name, "%s3PCFCovMatricesAll/c6_n%d_l%d_%d_%s.txt", out_file,nbin, max_l,index,suffix);
        
        FILE * C3File = fopen(c3name,"w"); // for c3 part of integral
        FILE * C4File = fopen(c4name,"w"); // for c4 part of integral
        FILE * C5File = fopen(c5name,"w"); // for c5 part of integral
        FILE * C6File = fopen(c6name,"w"); // for c6 part of integral
        
        for(int i=0;i<array_len;i++){
            for(int j=0;j<array_len;j++){
                fprintf(C3File,"%le\t",c3[i*array_len+j]);
                fprintf(C4File,"%le\t",c4[i*array_len+j]);
                fprintf(C5File,"%le\t",c5[i*array_len+j]);
                fprintf(C6File,"%le\t",c6[i*array_len+j]);
            }
            fprintf(C3File,"\n"); // new line each end of row
            fprintf(C4File,"\n");
            fprintf(C5File,"\n");
            fprintf(C6File,"\n");
        }
        
        fflush(NULL);
        
        // Close open files
        fclose(C3File);
        fclose(C4File);
        fclose(C5File);
        fclose(C6File);
        
        if(save_all==1){
            char binname6[1000];
            snprintf(binname6,sizeof binname6, "%s3PCFCovMatricesAll/binct_c6_n%d_l%d_%d_%s.txt",out_file, nbin,max_l,index, suffix);
            FILE * BinFile6 = fopen(binname6,"w");
        
            char binname5[1000];
            snprintf(binname5,sizeof binname5, "%s3PCFCovMatricesAll/binct_c5_n%d_l%d_%d_%s.txt",out_file, nbin,max_l,index, suffix);
            FILE * BinFile5 = fopen(binname5,"w");
            
            char binname4[1000];
            snprintf(binname4,sizeof binname4, "%s3PCFCovMatricesAll/binct_c4_n%d_l%d_%d_%s.txt",out_file, nbin,max_l,index, suffix);
            FILE * BinFile4 = fopen(binname4,"w");
            
            char binname3[1000];
            snprintf(binname3,sizeof binname3, "%s3PCFCovMatricesAll/binct_c3_n%d_l%d_%d_%s.txt",out_file, nbin,max_l,index,suffix);
            FILE * BinFile3 = fopen(binname3,"w");
            
            for(int i=0;i<array_len;i++){
                for(int j=0;j<array_len;j++){
                    fprintf(BinFile3,"%llu\t",binct3[i*array_len+j]);
                    fprintf(BinFile4,"%llu\t",binct4[i*array_len+j]);
                    fprintf(BinFile5,"%llu\t",binct5[i*array_len+j]);
                    fprintf(BinFile6,"%llu\t",binct6[i*array_len+j]);
                    }
                fprintf(BinFile3,"\n");
                fprintf(BinFile4,"\n");
                fprintf(BinFile5,"\n");
                fprintf(BinFile6,"\n");
            }
        
            fclose(BinFile3);
            fclose(BinFile4);
            fclose(BinFile5);
            fclose(BinFile6);
        }
    }

};

#endif
