#include "driver.h"
#include "string.h"

int main (int narg, char **arg)
{
  // analyse command line options
  int iarg = 1;
  int fmt = 1; // lammps-atomic-data

  while (narg > iarg){
    if (strcmp(arg[iarg],"-f") == 0){
       ++iarg;
       if (iarg >= narg) exit(0);
       if (strcmp(arg[iarg],"poscar") == 0) fmt = 2;
    }

    ++iarg;
  }

  Driver *driver = new Driver;

  driver->modify();
  driver->write(fmt);

  delete driver;

return 0;
}
