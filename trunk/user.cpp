#include "user.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#define MAXLINE 256

using namespace std;

/* -----------------------------------------------------------------------------
 * To select the orientation of the lattice
 * -------------------------------------------------------------------------- */
USER::USER() : lattice()
{
  initialized = 0;
  memory->create(name,10,"user:name");
  strcpy(name, "USER_LATT");
  char str[MAXLINE];

  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");
  printf("Input the file name if you want to read the unit cell information from\n");
  printf("a POSCAR-like file, or simply ENTER to read from stdin: ");
  fgets(str, MAXLINE, stdin);

  int flag = 1;
  if (count_words(str) > 0) flag = read_file(strtok(str," \t\n\r\f"));
  if (flag) flag = read_stdin();
  if (flag == 0) initialized = 1;

  for (int i = 0; i < 14; ++i) printf("====="); printf("\n");

return;
}

/* -----------------------------------------------------------------------------
 * Deconstructor does nothing
 * -------------------------------------------------------------------------- */
USER::~USER()
{

}

/* -----------------------------------------------------------------------------
 * To read the unit cell info from file
 * -------------------------------------------------------------------------- */
int USER::read_file(const char *fname)
{
  FILE *fp = fopen(fname,"r");
  // check file existence
  if (fp == NULL){
    printf("Error: File %s not found! Read from stdin instead.\n", fname);
    return 1;
  }

  // to read file
  // if read from file, the file would be similar to vasp POSCAR (direct)
  //  alat
  //  xx xy xz
  //  yx yy yz
  //  zx zy zz
  //  ntype1 ntype2 ntype3 ...
  //  sx1 sy1 sz1
  //  ...

  char str[MAXLINE];
  // scaling factor (lattice constant)
  fgets(str,MAXLINE,fp); if (feof(fp)){fclose(fp); return 2;}
  alat = numeric(strtok(str, " \t\n\r\f"));
  if (alat <= 0.) return 2;

  // basis vectors
  for (int i = 0; i < 3; ++i){
    fgets(str,MAXLINE,fp); if (feof(fp)){fclose(fp); return i+3;}
    latvec[i][0] = numeric(strtok(str,  " \t\n\r\f"));
    latvec[i][1] = numeric(strtok(NULL, " \t\n\r\f"));
    latvec[i][2] = numeric(strtok(NULL, " \t\n\r\f"));
  }
  fgets(str,MAXLINE,fp); if (feof(fp)){fclose(fp); return 6;}
  ntype = count_words(str);

  int *ntm = new int[ntype];
  char *ptr = strtok(str," \t\n\r\f");
  nucell = ntm[0] = inumeric(ptr);

  for (int i = 1; i < ntype; ++i){
    ptr = strtok(NULL," \t\n\r\f");
    ntm[i] = inumeric(ptr);
    nucell += ntm[i];
  }

  memory->create(atpos, nucell, 3, "USER_atpos");
  memory->create(attyp, nucell, "USER:attyp");

  int iatom =0;
  for (int ip = 0; ip < ntype;   ++ip)
  for (int i  = 0; i  < ntm[ip]; ++i){
    fgets(str,MAXLINE,fp); if (feof(fp)){fclose(fp); return 7;}
    atpos[iatom][0] = numeric(strtok(str,  " \t\n\r\f"));
    atpos[iatom][1] = numeric(strtok(NULL, " \t\n\r\f"));
    atpos[iatom][2] = numeric(strtok(NULL, " \t\n\r\f"));
    attyp[iatom++] = ip+1;
  }
  fclose(fp);
  delete []ntm;

  return 0;
}

/* -----------------------------------------------------------------------------
 * To read the unit cell info from input
 * -------------------------------------------------------------------------- */
int USER::read_stdin()
{
  char str[MAXLINE];
  // ask for lattice constant
  alat = 1.;
  printf("\nPlease input the lattice constant of the USER lattice [%g]: ", alat);
  if (count_words(fgets(str,MAXLINE,stdin))>0) alat = numeric(strtok(str, " \t\n\r\f"));
  if (alat <= 0.) alat = 1.;

  // ask for lattice vectors
  for (int i = 0; i < 3; ++i){
    while ( 1 ){
      printf("Please input the lattice vector A%d: ", i+1);
      if (count_words(fgets(str,MAXLINE,stdin)) < 3) continue;
      latvec[i][0] = numeric(strtok(str,  " \t\n\r\f"));
      latvec[i][1] = numeric(strtok(NULL, " \t\n\r\f"));
      latvec[i][2] = numeric(strtok(NULL, " \t\n\r\f"));

      break;
    }
  }

  // ask for # of atoms and # of atom types
  nucell = ntype = 1;
  printf("Please input the number of atoms per unit cell [1]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) nucell = inumeric(strtok(str, " \t\n\r\f"));
  if (nucell < 1) return 1;

  if (nucell != 1){
    printf("Please input the number of atom  types in cell [1]: ");
    if (count_words(fgets(str,MAXLINE,stdin))>0) ntype = inumeric(strtok(str, " \t\n\r\f"));
    if (ntype < 1) return 2;
    if (ntype > nucell) ntype = nucell;
  }
    
  memory->create(atpos, nucell, 3, "USER_atpos");
  memory->create(attyp, nucell, "USER:attyp");
  // ask for atom coordinates and types
  for (int i = 0; i < nucell; ++i){
    do printf("Please input [type xs ys zs] for atom %d: ", i+1);
    while (count_words(fgets(str,MAXLINE,stdin)) < 4);
    attyp[i]    = inumeric(strtok(str,  " \t\n\r\f"));
    atpos[i][0] = numeric(strtok(NULL, " \t\n\r\f"));
    atpos[i][1] = numeric(strtok(NULL, " \t\n\r\f"));
    atpos[i][2] = numeric(strtok(NULL, " \t\n\r\f"));
  }

return 0;
}
/* ------------------------------------------------------------------- */
