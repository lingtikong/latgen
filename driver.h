#ifndef DRIVER_H
#define DRIVER_H

#include "lists.h"
#include "random.h"
#include "input.h"
#include "memory.h"
#include "elements.h"
#include <map>

class Driver{
public:
  Driver(int);
  ~Driver();
  
  void generate();  // method to generate the atomic configuration
  void modify();    // to modify the resultant model

  void write(int);  // method to write the atomic configuration (xyz) and mapping info

  void FormLayers();
  void Interstitial();

#ifdef Poly
  void PolyCrystal();
#endif

  Memory *memory;
  UserInput *uin;

private:
  int ShowMenu(int);
  void MainMenu();
  void ShowVersion();

  lattice *latt;

  char *name;
  double alat;
  int nx, ny, nz, nucell;         // size in three dimension and # of atoms per unit cell
  int natom, ntype;               // otal number of atoms and atom types
  int *attyp, *numtype, *typeID;  // array to store the atomic types for all
  int *xmap, *ymap, *zmap, *umap; // arrays to store the mapping info
  double **atpos, latvec[3][3];   // arrays to store the atomic positions and lattice info
  int flag_orient;                // > 0, ask for orient; <= 0, not ask.
  void typescan();                // to scan the total number of atomic types in system
  int lookup(int);                // to find the ID of an atomic type

  RanPark *random;                // class object to create random numbers

  // private modification methods
  void solidsol(void);            // method to create subsutitutional solid solution
  void ResetTypeID(void);         // method to reset the atomic type ID
  void MapElement(void);          // method to map atomic type to real elements

  ChemElements * element;
  std::map<int,int> type2num;

};

#endif
