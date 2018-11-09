// cell_utilities.h - This contains many cell-related functions including the cell and particle classes. Based on code by Alex Wiegand.

#ifndef CELL_UTILITIES_H
#define CELL_UTILITIES_H

// We need a vector floor3 function
Float3 floor3(float3 p) {
    return Float3(floor(p.x), floor(p.y), floor(p.z));
}

// we need a vector ceil3 function
Float3 ceil3(float3 p) {
    return Float3(ceil(p.x), ceil(p.y), ceil(p.z));
}

// =================== Particles ====================
// This is the info about each particle that we load in and store in the Grid.

class Particle {
  public:
    Float3 pos;
    Float w;  // The weight for each particle
    Float JK; // The Jackknife region ID for each particle (stored as float for ease)
    int rand_class; // Integer 0 or 1 which defines which set of random particles this is in. 
    // This is set at read-in at random and used for the EE computation to avoid diagonal non-cancellation.
};


// ====================  The Cell and Grid classes ==================

/* The Grid class holds a new copy of the particles.
These are sorted into cells, and all positions are referenced to the
cell center.  That way, we can handle periodic wrapping transparently,
simply by the cell indexing.

For simplicity, we opt to flatten the index of the cells into a 1-d number.
For example, this makes multi-threading over cells simpler.
*/

class Cell {
  public:
    int start;	// The starting index of the particle list
    int np;
    int np1; // number of particles in cell in random-partition 1
    int np2;
};


#endif

