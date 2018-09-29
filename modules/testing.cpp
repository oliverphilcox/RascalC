class TestClass{

public:
    
    // initialize particles and list of cells
    Particle *prim_list; // NEED
    Particle *sec_list;
    Particle *thi_list;
    Particle *fou_list;
    Cell primary, sec, thi, fou; // NEED
    int pln, sln, tln, fln;
    double p1,p2,p3,p4;    
    Grid *grid;

private:
    integer3 prim_id;
    integer3 delta2, delta3, delta4; // DON'T NEED?
    
    
public:
    
int particle_list(integer3 id_3D, Cell &cell, Particle* &part_list, int &no_particles, float3 shift){
    // function returns a cell and a list of particles for a 3-dimensional ID
    // shift is a 3-vector shift applied for periodic boxes, to ensure parrticle positions are relative to cell center
    int id_1D = grid->test_cell(id_3D); // compute the 1d cell index from the 3d index
    
    if (id_1D<0) return 0; // cell is not within bounding box so will be skipped
    
    cell = grid->c[id_1D]; // cell object 
    no_particles = 0;
    // copy in list of particles into list
    for (int i = cell.start; i<cell.start+cell.np; i++, no_particles++){
        part_list[no_particles]=grid->p[i];
#ifdef PERIODIC
        part_list[no_particles].pos+= shift; // In periodic case particles positions are relative to cell center
#endif
        }        
    return id_1D;
    }
    

TestClass(Grid *_grid, Parameters *par){
    grid=_grid;
			
    CorrelationFunction *cf=new CorrelationFunction(par->corname,par->mbin,par->mumax-par->mumin);
	RandomDraws *rd=new RandomDraws(cf,par, NULL, 0);
    Integrals sumint(par,cf);
    
    int* bin;
    Float* xi;
    int mnp=grid->maxnp;
    // assign memory to bin, xi and lists of particles
    int ec=0;
    ec+=posix_memalign((void **) &prim_list, PAGE, sizeof(Particle)*mnp); // Not sure if memalign has a benefit over malloc here
    ec+=posix_memalign((void **) &sec_list, PAGE, sizeof(Particle)*mnp);
    ec+=posix_memalign((void **) &thi_list, PAGE, sizeof(Particle)*mnp);
    ec+=posix_memalign((void **) &fou_list, PAGE, sizeof(Particle)*mnp);

    ec+=posix_memalign((void **) &bin, PAGE, sizeof(int)*mnp*mnp);
    ec+=posix_memalign((void **) &xi, PAGE, sizeof(Float)*mnp*mnp);

    assert(ec==0);

    
void c2_contributions(){
    for (int i = 0; i<pln; i++){    // loop over i cell particles
        for (int j = 0; j<sln; j++){    // loop over j cell particles
            if(bin[sln*i+j]<0) continue; // skip if current pair outside binning range
            if (sec.start==primary.start&&i==j) continue; // exclude self count
            locint.second(bin[sln*i+j],prim_list[i],sec_list[j],p2);
            cnt2++;
        }
    }
}
    
    for(int ll=0;ll<mnp*mnp;ll++){
    	bin[ll]=9999999;
    }
    
    // SET THIS UP FOR OMP
    gsl_rng_env_setup();
    gsl_rng* locrng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(locrng,0); // just one thread at the moment
    Integrals locint(par,cf); // accumulates the integral contribution of each thread
    uint64 cnt2=0,cnt3=0,cnt4=0;

    for (int n_loops = 0; n_loops<par->max_loops; n_loops++){
        printf("Starting Integral Loop %d of %d\n",n_loops,par->max_loops);

        for (int n1=0; n1<grid->nf;n1++){
            // first loop over all the filled cells in the grid for the i cell
            
            // DO WE NEED TO RANORD CELLS I?
            
            int prim_id_1D = grid->filled[n1]; // 1d ID for cell i
            prim_id = grid->cell_id_from_1d(prim_id_1D); // define the 3D primary ID
            // define the cell (primary), the list of particles in the cell (prim_list) + the number of particles (pln)
            int x1 = particle_list(prim_id, primary, prim_list, pln, {0.,0.,0.}); 
            assert(x1>=0); // nonsense else!
            
            if(primary.np==0) continue; // skip if no particles in cell
            
            p1=1.;
            
            for (int n2=0; n2<par->N2; n2++){ // Determine second position
                // DRAW SECOND CELL (sec) FROM I CELL WEIGHTED BY 1/r^2 + up to max r only???
                delta2 = rd->random_cubedraw(locrng, &p2); // draw the second cell from a cube around the first (may draw the same cell again)
                // CHECK THIS FUNCTION
                
                
                // define list of particles
                int x2=particle_list(prim_id+delta2, sec, sec_list, sln, grid->cell_sep(delta2)); 
                
                if(x2<0) continue; // reject if outside the bounding box
                
                p2*=p1; // accumulate selection probability
                
                if(sec.np==0) continue;
                
                
                // COMPUTE C2 INTEGRAL FOR ALL PARTICLES IN prim + sec + thi
                
                for (int i = 0; i<pln; i++){    // loop over i cell particles
                    for (int j = 0; j<sln; j++){    // loop over j cell particles
                        if(bin[sln*i+j]<0) continue; // skip if current pair outside binning range
                        if (sec.start==primary.start&&i==j) continue; // exclude self count
                        locint.second(bin[sln*i+j],prim_list[i],sec_list[j],p2);
                        cnt2++;
                    }
                }
                
                // NB: Worth seeing which of these params are actually used by the integrand - can remove else.
                
                for (int n3=0; n3<par->N3; n3++){ // Determine third position
                    // DRAW THIRD CELL (thi) FROM J CELL WEIGHTED BY XI
                    
                    delta3 = rd->random_xidraw(locrng,&p3); // draw third cell from a region around cell j weighted by xi
                    int x3=particle_list(prim_id+delta2+delta3, thi, thi_list, tln, grid->cell_sep(delta2)+grid->cell_sep(delta3)); // define list of probabilities etc.
                
                    if(x3<0) continue; // reject if outside bounding box
                    
                    p3*=p2; // update probability
                    
                    if(thi.np==0) continue; // don't need empty cells here
                    // DOES THIS CHANGE THE PROBABILITIES??
                    
                    // COMPUTE C3 INTEGRAL FOR PARTICLES IN prim+sec+thi
                
                    for (int i = 0; i<pln; i++){    // loop over i cell particles
                        for (int j = 0; j<sln; j++){    // loop over j cell particles
                            
                            if(sec.start==primary.start&&i==j) continue; // exclude self count
                            if(bin[sln*i+j]<0) continue; // skip if current pair outside binning range
                            
                            for (int k = 0; k<tln; k++) {   // loop over k cell particles
                                if (thi.start==sec.start&&j==k) continue;   // exclude self-count
                                if (thi.start==primary.start&&i==k) continue;
                                locint.third(bin[sln*i + j],xi[tln*j + k],prim_list[i],sec_list[j],thi_list[k],p3); 
                                cnt3++;
                            }
                        }
                    }
                                    
                    
                    
                    for (int n4=0; n4<par->N4; n4++){ // Determine fourth position
                        // DRAW FOURTH CELL (fou) FROM I CELL WEIGHTED BY XI
                        
                        delta4 = rd->random_xidraw(locrng,&p4);
                        int x4=particle_list(prim_id+delta4, fou, fou_list, fln, grid->cell_sep(delta4)); 
                        
                        if(x4<0) continue;
                        
                        
                        p4*=p3;
                        
                        if(fou.np==0) continue;
                        
                        
                        // CHECK PROBABILITIES HERE
                
                        // COMPUTE C4 INTEGRAL FOR PARTICLES IN prim+sec+thi+fou
                
                        for (int i = 0; i<pln; i++){    // loop over i cell particles
                            for (int j = 0; j<sln; j++){    // loop over j cell particles
                                if(sec.start==primary.start&&i==j) continue; // exclude self count
                                if(bin[sln*i+j]<0) continue; // skip if current pair outside binning range
                            
                                for (int k = 0; k<tln; k++) {   // loop over k cell particles
                                    if (thi.start==sec.start&&j==k) continue;   // exclude self-count
                                    if (thi.start==primary.start&&i==k) continue;
                                    
                                    for (int l = 0; l<fln; l++) {
                                        if (fou.start==primary.start&&i==l) continue;   // Exclude self-count
                                        if (fou.start==sec.start&&j==l) continue;
                                        if (fou.start==thi.start&&k==l) continue;

                                        locint.fourth(bin[sln*i + j],xi[tln*j + k],prim_list[i],sec_list[j],thi_list[k],fou_list[l],p4);
                                        cnt4++;
                                    }
                                }
                            }
                        } // end particles
                    } // end l cell
                } // end k cell
            } // end j cell
        } // end i cell
    // UPDATE INTEGRALS HERE + CHECK CONVERGENCE
    } // end cycles
}
};

