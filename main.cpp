#include "driver.h"
#include "string.h"
#include "version.h"

void help()
{
  printf("\nlatgen  v 1.%d\nCode to generate atomic configuration for MD simulations.\n\n", VERSION);
  printf("The code is driven by menu, please follow the instructions step by step.\n");
  printf("By default, a LAMMPS data file, a compatible map file for fix-phonon, and\n");
  printf("an xyz format file will be written.\n");
  printf("\nUsage:\n    latgen [options]\n");
  printf("\nAvailable options:\n      -poscar : to output VASP POSCAR instead of LAMMPS data file.\n");
  printf("      -q      : to output atomic_style of charge with all q = 0.\n");
  printf("      -s      : to write all user input to script.inp, facilitating scripting.\n");
  printf("      -h      : to print this hlep info.\n");
  printf("\nAuthor: Lingti Kong, konglt@sjtu.edu.cn\n(C) 2026\n\n");

  exit(0);
}

int main (int narg, char **arg)
{
  // analyse command line options
  int iarg = 1;
  int fmt = 1; // lammps-atomic-data
  int save = 0; // to save user input or not.

  while (narg > iarg){
    if (strcmp(arg[iarg],"-poscar") == 0 || strcmp(arg[iarg],"-p") == 0){
      fmt = 2;

    } else if (strcmp(arg[iarg],"-q") == 0){
      fmt = 3;

    } else if (strcmp(arg[iarg],"-save") == 0 || strcmp(arg[iarg],"-s") == 0){
      save = 1;

    } else if (strcmp(arg[iarg],"-h") == 0){
      help();

    }
    ++iarg;
  }

  Driver *driver = new Driver(save);

  driver->modify();
  driver->write(fmt);

  delete driver;

return 0;
}
