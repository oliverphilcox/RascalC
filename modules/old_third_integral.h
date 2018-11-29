    
    inline void third(const Particle* pi_list, const int* prim_ids, const int pln, const Particle pj, const Particle pk, const int pj_id, const int pk_id, const int* bin_ij, const Float* wij, Float* &xi_jk, Float* &xi_ik, Float* wijk, const double prob){
        // Accumulates the three point integral C3. Also outputs an array of xi_ik and bin_ik values for later reuse.
        
        // First define variables:
        Particle pi;
        Float rik_mag, rik_mu, c3v, c3vj, rjk_mag, rjk_mu, tmp_weight, xi_jk_tmp, xi_ik_tmp, JK_weight;
        int tmp_bin, tmp_full_bin;
        
        cleanup_l(pj.pos,pk.pos,rjk_mag,rjk_mu); 
        xi_jk_tmp = cf->xi(rjk_mag, rjk_mu); // correlation function for j-k
        
        if (pk_id==pj_id) return;
        
        for(int i=0;i<pln;i++){ // Iterate ovr particle in pi_list
            pi = pi_list[i];
            if(wij[i]==-1){
                wijk[i]=-1;
                continue; // skip incorrect bins / ij self counts
            }
            
            if(prim_ids[i]==pk_id){
                wijk[i]=-1; // re-read to skip this later
                continue; // don't self-count
            }
            
            cleanup_l(pi.pos,pk.pos,rik_mag,rik_mu); // define angles/lengths
            
            tmp_bin = getbin(rik_mag,rik_mu); // bin for each particle
            
            tmp_weight = wij[i]*pk.w; // product of weights, w_iw_jw_k
            xi_ik_tmp = cf->xi(rik_mag, rik_mu); // not used here but used later
            
            // save arrays for later
            xi_jk[i]=xi_jk_tmp;
            xi_ik[i]=xi_ik_tmp;
            wijk[i]=tmp_weight;
            
            if ((tmp_bin<0)||(tmp_bin>=nbin*mbin)){
                // Don't add contributions of this to C3 - but STILL save xi_ik etc. for later
                continue; // if not in correct bin
            }
            // Compute jackknife weight tensor:
            JK_weight=weight_tensor(int(pi.JK),int(pj.JK),int(pk.JK),int(pi.JK),bin_ij[i],tmp_bin);
            
            // Now compute the integral (using symmetry factor of 4);
            c3v = 4.*tmp_weight*pi.w/prob*xi_jk_tmp;
            c3vj = c3v*JK_weight;
            
            // Add to local counts
            tmp_full_bin = bin_ij[i]*mbin*nbin+tmp_bin;
            c3[tmp_full_bin]+=c3v;
            c3j[tmp_full_bin]+=c3vj;
            binct3[tmp_full_bin]++;
            errc3[tmp_full_bin]+=pow(c3v,2.);
            errc3j[tmp_full_bin]+=pow(c3vj,2.);
        }
    }
     
