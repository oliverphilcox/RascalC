class TestClass{

	public:
    
    // initialize particles and list of cells
    Particle *prim_list;
    Particle *sec_list;
    Particle *thi_list;
    Particle *fou_list;
    Cell primary, sec, thi, fou;
    int pln, sln, tln, fln;
    integer3 prim_id, sec_id, thi_id, fou_id;
    integer3 delta2, delta3, delta4;
    double p1,p2,p3,p4;
    
    Grid *grid;

    
    

void particle_list(Cell &cell, Particle* &part_list, int id_1D, int &no_particles, float3 shift){
    // function returns a cell and a list of particles for a 1-dimensional ID
    // shift is a 3-vector shift applied for periodic boxes, to ensure parrticle positions are relative to cell center
    cell = grid->c[id_1D]; // cell object 
    no_particles = 0;
    // copy in list of particles into list
    for (int i = cell.start; i<cell.start+cell.np; i++, no_particles++){
        part_list[no_particles]=grid->p[i];
#ifdef PERIODIC
        part_list[no_particles].pos+= shift; // In periodic case particles positions are relative to cell center
#endif
        }        
    }

TestClass(Grid *_grid, Parameters *par){
    grid=_grid;
			
    CorrelationFunction *cf=new CorrelationFunction(par->corname,par->mbin,par->mumax-par->mumin);
	RandomDraws *rd=new RandomDraws(cf,par, NULL, 0);
    // SET THIS UP FOR OMP
    gsl_rng_env_setup();
    gsl_rng* locrng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(locrng,0); // just one thread at the moment

    for (int n_loops = 0; n_loops<par->max_loops; n_loops++){

        for (int n1=0; n1<grid->nf;n1++){
            // first loop over all the filled cells in the grid for the i cell
            
            // DO WE NEED TO RANORD CELLS I?
            
            int prim_id_1D = grid->filled[n1]; // 1d ID for cell i
            prim_id = grid->cell_id_from_1d(prim_id_1D); // define the 3D primary ID
            // define the cell (primary), the list of particles in the cell (prim_list) + the number of particles (pln)
            particle_list(primary, prim_list, prim_id_1D, pln, {0.,0.,0.}); 
            
            if(primary.np==0) continue; // no particles in cell
            
            p1=1.;
            
            for (int n2=0; n2<par->N2; n2++){ // Determine second position
                // DRAW SECOND CELL (sec) FROM I CELL WEIGHTED BY 1/r^2 + up to max r only???
                delta2 = rd->random_cubedraw(locrng, &p2);
                sec_id = prim_id + delta2; // draw the second cell from a cube around the first (may draw the same cell again)
                // CHECK THIS FUNCTION
                int sec_id_1D=grid->test_cell(sec_id);
                if(sec_id_1D<0) continue; // reject if outside the bounding box
                
                // define list of particles
                particle_list(sec, sec_list, sec_id_1D, sln, grid->cell_sep(delta2)); 
                
                p2*=p1; // accumulate selection probability
                
                if(sec.np==0) continue;
                
                // COMPUTE C2 INTEGRAL FOR ALL PARTICLES IN prim + sec + thi
                // exclude self-counts for particles
                // UPDATE COUNTS cnt2++
                
                // NB: Worth seeing which of these params are actually used by the integrand - can remove else.
                
                for (int n3=0; n3<par->N3; n3++){ // Determine third position
                    // DRAW THIRD CELL (thi) FROM J CELL WEIGHTED BY XI
                    
                    delta3 = rd->random_xidraw(locrng,&p3); // draw third cell from a region around cell j weighted by xi
                    thi_id = sec_id+delta3; // 3D cell ID
                    int thi_id_1D=grid->test_cell(thi_id);
                    if(thi_id_1D<0) continue; // reject if outside bounding box
                    
                    particle_list(thi, thi_list, thi_id_1D, tln, grid->cell_sep(delta2)+grid->cell_sep(delta3)); // define list of probabilities etc.
                
                    p3*=p2; // update probability
                    
                    if(thi.np==0) continue; // don't need empty cells here
                    // DOES THIS CHANGE THE PROBABILITIES??
                    
                    // COMPUTE C3 INTEGRAL FOR PARTICLES IN prim+sec+thi
                    // exclude self-counts for particles
                    // cnt3++
                    
                    for (int n4=0; n4<par->N4; n4++){ // Determine fourth position
                        // DRAW FOURTH CELL (fou) FROM I CELL WEIGHTED BY XI
                        
                        delta3 = rd->random_xidraw(locrng,&p4);
                        fou_id = prim_id+delta4; // define l cell ID from i cell
                        int fou_id_1D = grid->test_cell(fou_id);
                        if(fou_id_1D<0) continue;
                        
                        particle_list(fou, fou_list, fou_id_1D, fln, grid->cell_sep(delta4)); 
                
                        p4*=p3;
                        
                        if(fou.np==0) continue;
                        
                        
                        // CHECK PROBABILITIES HERE
                        // COMPUTE C4 INTEGRAL FOR PARTICLES IN prim+sec+thi+fou
                        // exclude self-counts for particles
                        // cnt4++
                    }
                }
            }
        }
        }
    // UPDATE INTEGRALS HERE + CHECK CONVERGENCE
    }


};

