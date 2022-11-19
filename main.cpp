#include "driver.h"
#include "string.h"
#include "version.h"

void help()
{
  printf("\nlatgen  v 1.%d\nCode to generate atomic configuration for MD simulations.\n\n", VERSION);
  printf("The code is driven by menu, please follow the instructions step by step.\n");
  printf("By default, a LAMMPS data file, a compatible map file for fix-phonon, and\n");
  printf("an xyz format file will be written.\n");
  printf("\nUsage:\n    latgen [-poscar] [-h]\n");
  printf("\nOptions:\n      -poscar : to output VASP POSCAR instead of LAMMPS data file.\n");
  printf("      -q      : to output atomic_style of charge with all q = 0.\n");
  printf("      -h      : to print this hlep info.\n");
  printf("\nAuthor: Lingti Kong, konglt@sjtu.edu.cn\n(C) 2021\n\n");

  exit(0);
}

int main (int narg, char **arg)
{
  // analyse command line options
  int iarg = 1;
  int fmt = 1; // lammps-atomic-data

  while (narg > iarg){
    if (strcmp(arg[iarg],"-poscar") == 0 || strcmp(arg[iarg],"-p") == 0){
      fmt = 2;

    } else if (strcmp(arg[iarg],"-h") == 0){
      help();

    } else if (strcmp(arg[iarg],"-q") == 0){
      fmt = 3;

    }
    ++iarg;
  }

  Driver *driver = new Driver;

  driver->modify();
  driver->write(fmt);

  delete driver;

return 0;
}
