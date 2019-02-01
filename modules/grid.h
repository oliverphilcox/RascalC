//grid class for grid_covariance.cpp file (modified from Alex Wiegand)
#include "cell_utilities.h"

#ifndef GRID_H
#define GRID_H

class Grid {
  public:
    Float3 rect_boxsize; // 3D dimensions of the periodic volume
    int nside, ncells;       // Grid size (per linear and per volume)
    Cell *c;		// The list of cells
    Float cellsize;   // Size of one cell
    Float max_boxsize; // largest dimension of the cuboid box
    Particle *p;	// Pointer to the list of particles
    int np,np1,np2;		// Number of particles (total and number in each partition
    integer3 nside_cuboid; // number of cells along each dimension of cuboidal box
    int np_pos;		// Number of particles
    int *pid;		// The original ordering
    int *filled; //List of filled cells
    int nf;      //Number of filled cells
    int maxnp;   //Max number of particles in a single cell
    Float norm; // sum_weights randoms / sum_weights galaxies for normalization
    Float sum_weights; // total summed weights
    Float sumw_pos, sumw_neg; // Summing the weights

    int test_cell(integer3 cell){
    	// returns -1 if cell is outside the grid or wraps around for the periodic grid
#ifndef PERIODIC
    	if(nside_cuboid.x<cell.x||cell.x<0||nside_cuboid.y<cell.y||cell.y<0||nside_cuboid.z<cell.z||cell.z<0)
            return -1;
    	else
#endif
    		return wrap_cell(cell);
    }

    int wrap_cell(integer3 cell) {
        // Return the 1-d cell number, after wrapping
        // We apply a very large bias, so that we're
        // guaranteed to wrap any reasonable input.
        int cx = (cell.x+ncells)%nside_cuboid.x;
        int cy = (cell.y+ncells)%nside_cuboid.y;
        int cz = (cell.z+ncells)%nside_cuboid.z;
        // return (cx*nside+cy)*nside+cz;
        int answer = (cx*nside_cuboid.y+cy)*nside_cuboid.z+cz;
        assert(answer<ncells&&answer>=0);

        return answer;
    }

    integer3 cell_id_from_1d(int n) {
	// Undo 1d back to 3-d indexing
        assert(n>=0&&n<ncells);
        
        integer3 cid;
        cid.z = n%nside_cuboid.z;
        n = n/nside_cuboid.z;
        cid.y = n%nside_cuboid.y;
        cid.x = n/nside_cuboid.y;
        return cid;
    }

    int pos_to_cell(Float3 pos) {
        // Return the 1-d cell number for this position, properly wrapped
        // We assume the first cell is centered at cellsize/2.0
        // return wrap_cell( floor3(pos/cellsize+Float3(0.5,0.5,0.5)));
        return wrap_cell( floor3(pos/cellsize));
    }

    Float3 cell_centered_pos(Float3 pos) {
        // Subtract off the cell center from the given position.
        // This is safe for positions not in the primary box.
        return pos-cellsize*(floor3(pos/cellsize)+Float3(0.5,0.5,0.5));
    }

    Float3 cell_sep(integer3 sep) {
        // Return the position difference corresponding to a cell separation
        return cellsize*sep;
    }
    
    void copy(Grid *g){
        // Copy grid object
        rect_boxsize=g->rect_boxsize; 
        nside=g->nside;
        ncells=g->ncells;
        cellsize=g->cellsize;
        max_boxsize=g->max_boxsize;
        np=g->np;
        np1=g->np1;
        np2=g->np2;
        nside_cuboid = g->nside_cuboid;
        np_pos = g->np_pos;
        norm = g->norm;
        nf=g->nf;
        maxnp=g->maxnp;
        sumw_pos=g->sumw_pos;
        sumw_neg=g->sumw_neg;
        sum_weights=g->sum_weights;
        
        // Allocate memory:
        p = (Particle *)malloc(sizeof(Particle)*np);
        pid = (int *)malloc(sizeof(int)*np);
        c  = (Cell *)malloc(sizeof(Cell)*ncells);
        filled = (int *)malloc(sizeof(int)*nf);
	
        // Copy in lists elementwise
        for(int j=0;j<ncells;j++) c[j]=g->c[j];
        for(int j=0;j<np;j++) p[j]=g->p[j];
        for(int j=0;j<np;j++) pid[j]=g->pid[j];
        for(int j=0;j<nf;j++) filled[j]=g->filled[j];
    }

    ~Grid() {
	// The destructor
        free(p);
        free(pid);
        free(c);
        free(filled);
        return;
    }
    
    Grid(){
       //empty constructor
    }

    Grid(Particle *input, int _np, Float3 _rect_boxsize, int _nside, Float3 shift, Float nofznorm) {
        // The constructor: the input set of particles is copied into a
        // new list, which is ordered by cell.
        // After this, Grid is self-sufficient; one could discard *input
        rect_boxsize = _rect_boxsize;
        nside = _nside;
        assert(nside<1025);   // Can't guarantee won't spill int32 if bigger
        np = _np;
        np_pos = 0;
        max_boxsize=fmax(rect_boxsize.x,fmax(rect_boxsize.y,rect_boxsize.z));
        assert(max_boxsize>0&&nside>0&&np>=0);
        cellsize = max_boxsize/nside;
        nside_cuboid = integer3(ceil3(rect_boxsize/cellsize));
        ncells = nside_cuboid.x*nside_cuboid.y*nside_cuboid.z;
            
        p = (Particle *)malloc(sizeof(Particle)*np);
        pid = (int *)malloc(sizeof(int)*np);
        c = (Cell *)malloc(sizeof(Cell)*ncells);

        // Now we want to copy the particles, but do so into grid order.
        // First, figure out the cell for each particle
        // Shift them to the primary volume first
        int *cell = (int *)malloc(sizeof(int)*np);
        for (int j=0; j<np; j++) cell[j] = pos_to_cell(input[j].pos - shift);
        
        // Histogram the number of particles in each cell
        int *incell = (int *)malloc(sizeof(int)*ncells);
        for (int j=0; j<ncells; j++) incell[j] = 0.0;
        for (int j=0; j<np; j++) incell[cell[j]]++;

        // Create list of filled cells
        nf=0;
        for (int j=0; j<ncells; j++) if(incell[j]>0) nf++;
        filled = (int *)malloc(sizeof(int)*nf);
        for (int j=0,k=0; j<ncells; j++) if(incell[j]>0) filled[k++]=j;

        printf("\nThere are %d filled cells compared with %d total cells.\n",nf,ncells);
        
        // Count the number of positively weighted particles + total weights
        sumw_pos = sumw_neg = sum_weights = 0.0;
        for (int j=0; j<np; j++){
            sum_weights+=input[j].w;
            if (input[j].w>=0) {
                np_pos++;
            sumw_pos += input[j].w;
            } else {
            sumw_neg += input[j].w;
            }
	}

        // Cumulate the histogram, so we know where to start each cell
        for (int j=0, tot=0; j<ncells; tot+=incell[j], j++) {
            c[j].start = tot;
            c[j].np = 0;  // We'll count these as we add the particles
            c[j].np1 = 0;
            c[j].np2 = 0;
        }

        // Initialize number of particles in each partition
        np1=0;
        np2=0;
        
        // Copy the particles into the cell-ordered list
        for (int j=0; j<np; j++) {
            Cell *thiscell = c+cell[j];
            int index = thiscell->start+thiscell->np;
            p[index] = input[j];
#ifdef PERIODIC
            p[index].pos = cell_centered_pos(input[j].pos);
                // Switch to cell-centered positions
#endif
            pid[index] = j;	 // Storing the original index

        
            thiscell->np += 1;
            if(p[index].rand_class==0){
                thiscell->np1+=1;
                np1++;
            }
            if(p[index].rand_class==1){
                thiscell->np2+=1;
                np2++;
            }
        }

        // Checking that all is well.
        int tot = 0;
        maxnp=0;
        for (int j=0; j<ncells; j++) {
            assert(c[j].start == tot);
            assert(c[j].np == incell[j]);
            if(c[j].np>maxnp) maxnp=c[j].np;
            tot += c[j].np;
        }
        free(incell);
        assert(tot == np);

        // compute normalization
        norm = sum_weights/nofznorm;

        free(cell);
        
        return;
        }

};   // End Grid class

#endif
