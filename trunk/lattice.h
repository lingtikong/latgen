#ifndef LATTICE_H
#define LATTICE_H

#include "memory.h"

class lattice {
public:
  
  lattice();
  ~lattice();

  Memory *memory;

  char *name;          // name of the lattice
 
  int initialized;     // flag to indicate if initialization successful (1) or not (0)
  int flag_z_perp_xy;  // flag to indicate that z is perpendicular to xy plane

  int nucell, ntype;   // number of atoms and atomic types in a unit cell
  int    *attyp;       // array to store atomic types
  double alat;         // lattice constant of lattice
  double latvec[3][3]; // lattice vectors, in unit of alat
  double **atpos;      // fractional coordinate for atoms in unit cell

  int nlayer;          // # of layers per unit cell; will count automatically in setup
  int *layer;          // layer ID for each atom
  int *numlayer;       // # of atoms per layer; will count automatically in setup
  double *h;           // height above each layer; will be set automatically in setup

  void display();      // method to display lattice info

  void OrientLattice();// to orient the lattice, following the rule of LAMMPS.

  int count_words(const char *);

private:
  void setup();        // to setup "numlayer", "h"
};
#endif
