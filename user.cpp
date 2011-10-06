#include "user.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#define MAX_LINE_LENGTH 256

using namespace std;

/* ----------------------------------------------------------------------
   To select the orientation of the lattice
------------------------------------------------------------------------- */
USER::USER() : lattice()
{
  initialized = 0;
  name = new char[10];
  strcpy(name, "USER_LATT");
  char str[MAX_LINE_LENGTH];
  printf("\n"); for (int i=0; i<70; i++) printf("="); printf("\n");
  printf("Would you like to read the lattice info from a file?(y/n)[n]: ");
  gets(str);
  if (strcmp(str,"y")==0 || strcmp(str,"Y")==0){
    FILE *fp;
    do printf("\nPlease input the name of the file that contains the lattice info: ");
    while (count_words(gets(str)) < 1);
    fp = fopen(str,"r");
    if (fp == NULL){
      printf("Error: File %s not found!\n", str);
      return;
    }
    // if read from file, the file would be similar to vasp POSCAR
    //  alat
    //  xx xy xz
    //  yx yy yz
    //  zx zy zz
    //  ntype1 ntype2 ntype3 ...
    //  sx1 sy1 sz1 layer-ID
    //  ...
    fgets(str,MAX_LINE_LENGTH,fp); if (feof(fp)){fclose(fp); return;}
    sscanf(str,"%lg",&alat);
    for (int i=0; i<3; i++){
      fgets(str,MAX_LINE_LENGTH,fp); if (feof(fp)){fclose(fp); return;}
      sscanf(str,"%lg %lg %lg",&latvec[i][0], &latvec[i][1], &latvec[i][2]);
    }
    fgets(str,MAX_LINE_LENGTH,fp); if (feof(fp)){fclose(fp); return;}
    ntype = count_words(str);

    int *ntm = new int[ntype];
    char *ptr = strtok(str," \t\n\r\f");
    nucell = ntm[0] = atoi(ptr);

    for (int i=1; i<ntype; i++){
      ptr = strtok(NULL," \t\n\r\f");
      ntm[i] = atoi(ptr);
      nucell += ntm[i];
    }

    atpos = memory->create(atpos,nucell, 3, "USER_atpos");
    attyp = new int[nucell];
    layer = new int[nucell];

    int iatom =0;
    for (int ip=0; ip<ntype; ip++){
      for (int i=0; i<ntm[ip]; i++){
        fgets(str,MAX_LINE_LENGTH,fp); if (feof(fp)){fclose(fp); return;}
        sscanf(str,"%lg %lg %lg", &atpos[iatom][0],&atpos[iatom][1],&atpos[iatom][2], &layer[iatom]);
        attyp[iatom++] = ip+1;
      }
    }
    fclose(fp);
    delete []ntm;
  } else {
    // ask for lattice constant
    alat = 1.;
    printf("\nPlease input the lattice constant of the USER lattice [1.0]: ");
    if (count_words(gets(str))>0) sscanf(str,"%lg",&alat);

    // ask for lattice vectors
    for (int i=0; i<3; i++){
      while (1){
        printf("Please input the lattice vector A%d: ", i+1);
        if (count_words(gets(str))<3) continue;
        int nr = sscanf(str,"%lg %lg %lg", &latvec[i][0], &latvec[i][1], &latvec[i][2]);
        if (nr == 3) break;
      }
    }
    // ask for # of atoms and # of atom types
    nucell = ntype = 1;
    printf("Please input the number of atoms per unit cell [1]: ");
    if (count_words(gets(str)) > 0) sscanf(str,"%d",&nucell);
    if (nucell < 1) return;
  
    if (nucell != 1){
      printf("Please input the number of atom types in cell [1]: ");
      if (count_words(gets(str))>0) sscanf(str,"%d",&ntype);
      if (ntype < 1) return;
      if (ntype > nucell) ntype = nucell;
    }
      
    atpos = memory->create(atpos,nucell, 3, "USER_atpos");
    attyp = new int[nucell];
    layer = new int[nucell];
    // ask for atom coordinates and types
    for (int i=0; i<nucell; i++){
      int nr;
      while (1){
        printf("Please input [type xs ys zs layerID] for atom %d: ", i+1);
        if (count_words(gets(str)) < 4) continue;
        nr = sscanf(str,"%d %lg %lg %lg %d", &attyp[i],&atpos[i][0], &atpos[i][1], &atpos[i][2], &layer[i]);
        if (nr == 5) break;
      }
    }
  }
  printf(""); for (int i=0; i<70; i++) printf("="); printf("\n");
  initialized = 1;

return;
}

/* ----------------------------------------------------------------------
   Deconstructor does nothing
------------------------------------------------------------------------- */
USER::~USER()
{

}

