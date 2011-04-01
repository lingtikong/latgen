#include "driver.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <time.h>

#define MAX_LINE_LENGTH 256
/* ----------------------------------------------------------------------
   constructor to initialize
------------------------------------------------------------------------- */
Driver::Driver()
{
  char str[MAX_LINE_LENGTH];
  int lattyp = 1;
  // print out the menu
  printf("\n===========================================================\n");
  printf("Please select the lattice type of your system:\n");
  printf(" 1. FCC/NaCl/Diamond;      |  4. A3B;\n");
  printf(" 2. BCC;                   |  5. A2B;\n");
  printf(" 3. HCP;                   |  6. AB;\n");
  printf("---------------------------+-------------------------------\n");
  printf(" 7. User defined;          |  0. Exit.\n");
  printf("---------------------------+-------------------------------\n");
  printf("Your  choice[1]: ");
  if (strlen(gets(str)) > 0) sscanf(str,"%d", &lattyp);
  printf("You selected: %d", lattyp);
  printf("\n===========================================================\n");
  // select the job
  switch (lattyp){
  case 1: latt = new FCC; break;
  case 2: latt = new BCC; break;
  case 3: latt = new HCP; break;
  case 4: latt = new A3B; break;
  case 5: latt = new A2B; break;
  case 6: latt = new AB; break;
  case 7: latt = new USER; break;
  default: exit(1);}
  // re-orient the lattice
  printf("Would you like to re-orient the unit cell? (y/n)[n]: ");
  if (strlen(gets(str)) > 0){
    if (strcmp(str,"y")==0 || strcmp(str,"Y")==0) latt->OrientLattice();
  }

  latt->display();
}

/* ----------------------------------------------------------------------
   deconstructor to free memory
------------------------------------------------------------------------- */
Driver::~Driver()
{
  if (atpos!= NULL) latt->memory->destroy_2d_double_array(atpos);
  if (attyp!= NULL) delete []attyp;
  if (xmap != NULL) delete []xmap;
  if (ymap != NULL) delete []ymap;
  if (zmap != NULL) delete []zmap;
  if (umap != NULL) delete []umap;
  if (latt != NULL) delete latt;

  if (typeID!= NULL) delete []typeID;
  if (numtype!= NULL) delete []numtype;
  if (random != NULL) delete random;
}

/* ----------------------------------------------------------------------
   method to generate atomic configurations
------------------------------------------------------------------------- */
void Driver::generate()
{
  char str[MAX_LINE_LENGTH];
  int leading_dir = 1;
  printf("\n===========================================================\n");
  while (1){
    printf("Please input the dimensions in x, y, and z directions: ");
    if (latt->count_words(gets(str)) < 3) continue;
    sscanf(str,"%d %d %d", &nx, &ny, &nz);
    natom = nx*ny*nz*latt->nucell;
    if (natom > 0) break;
  }

  printf("Your system would be %d x %d x %d with %d atoms.\n",nx,ny,nz,natom);
  printf("Please indicate which direction should goes fast(1:x; other: z)[1]: ");
  if (latt->count_words(gets(str)) > 0) sscanf(str,"%d", &leading_dir);
  printf("===========================================================\n");

  atpos = latt->memory->create_2d_double_array(natom, 3, "driver->generate:atpos");
  xmap = new int[natom];
  ymap = new int[natom];
  zmap = new int[natom];
  umap = new int[natom];
  attyp = new int[natom];

  int iatom = 0;
  if ( leading_dir == 1){
    for (int k=0; k<nz; k++){
      for (int j=0; j<ny; j++){
        for (int i=0; i<nx; i++){
          for (int u=0; u<latt->nucell; u++){
            xmap[iatom] = i;
            ymap[iatom] = j;
            zmap[iatom] = k;
            umap[iatom] = u;
            attyp[iatom] = latt->attyp[u];
            atpos[iatom][0] = latt->atpos[u][0] + double(i);
            atpos[iatom][1] = latt->atpos[u][1] + double(j);
            atpos[iatom][2] = latt->atpos[u][2] + double(k);
            iatom++;
          }
        }
      }
    }
  }else{
    for (int i=0; i<nx; i++){
      for (int j=0; j<ny; j++){
        for (int k=0; k<nz; k++){
          for (int u=0; u<latt->nucell; u++){
            xmap[iatom] = i;
            ymap[iatom] = j;
            zmap[iatom] = k;
            umap[iatom] = u;
            attyp[iatom] = latt->attyp[u];
            atpos[iatom][0] = latt->atpos[u][0] + double(i);
            atpos[iatom][1] = latt->atpos[u][1] + double(j);
            atpos[iatom][2] = latt->atpos[u][2] + double(k);
            iatom++;
          }
        }
      }
    }
  }
  // convert fractional coordinate to cartesian
  double tmp[3];
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = latt->latvec[i][j]*latt->alat;
  }
  for (int i=0; i<natom; i++){
    for (int idim=0; idim<3; idim++) tmp[idim] = atpos[i][idim];
    atpos[i][0] = tmp[0]*latvec[0][0] + tmp[1]*latvec[1][0] + tmp[2]*latvec[2][0];
    atpos[i][1] = tmp[0]*latvec[0][1] + tmp[1]*latvec[1][1] + tmp[2]*latvec[2][1];
    atpos[i][2] = tmp[0]*latvec[0][2] + tmp[1]*latvec[1][2] + tmp[2]*latvec[2][2];
  }
  for (int i=0; i<3; i++){
    latvec[0][i] *= double(nx);
    latvec[1][i] *= double(ny);
    latvec[2][i] *= double(nz);
  }
  // find the total # of types and # of atoms for each type
  typescan();

return;
}

/* ----------------------------------------------------------------------
   method to find the total # of types and # of atoms for each type
------------------------------------------------------------------------- */
void Driver::typescan()
{
  // allocate memory
  int typmax = 10;
  if (typeID != NULL) delete []typeID;
  if (numtype!= NULL) delete []numtype;
  typeID  = new int[typmax];
  numtype = new int[typmax];
  for (int i=0; i<typmax; i++) numtype[i] = 0;

  ntype = 0;
  // now to identify the total number of types
  for (int i=0; i<natom; i++){
    int id = lookup(attyp[i]);
    if (id < 0){
      if (ntype == typmax){
        int *int1 = new int[ntype];
        int *int2 = new int[ntype];
        for (int k=0; k<ntype; k++){int1[k] = typeID[k]; int2[k] = numtype[k];}
        delete []typeID;
        delete []numtype;

        typmax += 5;
        typeID  = new int[typmax];
        numtype = new int[typmax];
        for (int k=0; k<ntype; k++){typeID[k]=int1[k]; numtype[k]=int2[k];}
        for (int k=ntype; k<typmax; k++) numtype[k] = 0;
        delete []int1; delete []int2;
      }
      typeID[ntype] = attyp[i];
      id            = ntype;
      ntype++;
    }
    numtype[id]++;
  }
return;
}

/* ----------------------------------------------------------------------
   method to find the ID of the current atomic type
------------------------------------------------------------------------- */
int Driver::lookup(int ip)
{
  for (int i=0; i<ntype; i++){if (ip == typeID[i]) return i;}

  return -1;
}

/* ----------------------------------------------------------------------
   method to write out atomic configuraton and mapping info
------------------------------------------------------------------------- */
void Driver::write()
{
  FILE *fp;
  char str[MAX_LINE_LENGTH], *posfile, *mapfile;
  printf("\n===========================================================\n");
  printf("Please input the filename of the output xyz file [atomcfg.xyz]: ");
  if (latt->count_words(gets(str)) > 0){
    int n = strlen(str) + 1;
    posfile = new char[n];
    strcpy(posfile, str);
  } else {
    posfile = new char[12];
    strcpy(posfile, "atomcfg.xyz");
  }

  printf("Please input the filename of the output map file [map.in]: ");
  if (latt->count_words(gets(str)) > 0){
    int n = strlen(str) + 1;
    mapfile = new char[n];
    strcpy(mapfile, str);
  } else {
    mapfile = new char[7];
    strcpy(mapfile, "map.in");
  } 

  printf("\nThe atomic configuration will be written to file: %s\n", posfile);
  printf("The FFT map information  will be written to file: %s", mapfile);
  printf("\n===========================================================\n");
  // write the xyz position file 
  fp = fopen(posfile, "w");
  fprintf(fp, "%d\n", natom);
  fprintf(fp, "%s cell with dimension %d x %d x %d and a = %g\n",latt->name, nx, ny, nz, latt->alat);
  int nr = 3;
  if (natom < nr) nr = natom;
  for (int i=0; i<nr; i++){
    fprintf(fp,"%d %16.16e %16.16e %16.16e crystal_vector %d %16.16e %16.16e %16.16e\n", attyp[i], atpos[i][0],
    atpos[i][1], atpos[i][2], i+1, latvec[i][0], latvec[i][1], latvec[i][2]);
  }
  for (int i=nr; i<natom; i++)
    fprintf(fp,"%d %16.16e %16.16e %16.16e\n", attyp[i], atpos[i][0], atpos[i][1], atpos[i][2]);
  fclose(fp);
  // write the map file, useful to fix_phonon only.
  fp = fopen(mapfile, "w");
  fprintf(fp,"%d %d %d %d\n", nx, ny, nz, latt->nucell);
  fprintf(fp,"Map file for %dx%dx%d %s cell.\n",nx,ny,nz,latt->name);
  for (int i=0; i<natom; i++)
    fprintf(fp,"%d %d %d %d %d\n", xmap[i], ymap[i], zmap[i], umap[i], i+1);
  fclose(fp);

  delete []posfile;
  delete []mapfile;
return;
}

/* ----------------------------------------------------------------------
   method to modify the resultant model
------------------------------------------------------------------------- */
void Driver::modify()
{
  char str[MAX_LINE_LENGTH];
  int ncycle = 1;
  while (ncycle){
    int job=0;
    // to display the menu for modification
    printf("\n===========================================================\n");
    printf("Please select the modification you want to do:\n");
    printf("  1. Create substitutional solid solution;\n");

    if (ncycle == 1) printf("  0. Nothing.\n");
    else printf("  0. Done.\n");
    printf("Your choice[0]: ");

    if (latt->count_words(gets(str)) >0) sscanf(str,"%d",&job);

    printf("You selected: %d", job);
    printf("\n===========================================================\n");
   
    switch (job){ 
    case 1: solidsol(); break;
    default: return;
    }

    ncycle++;
  }
return;
}

/* ----------------------------------------------------------------------
   private method to create substitutional solid solution
------------------------------------------------------------------------- */
void Driver::solidsol()
{
  char str[MAX_LINE_LENGTH];
  printf("\n===========================================================\n");
  int lrange = 0, nrange;
  printf("Limit the solid solution within a region? (y/n)[n]: ");
  if (strlen(gets(str)) > 0) if (strcmp(str,"y") == 0 || strcmp(str,"Y")== 0) lrange = 1;
  int dir = 2, outside = 0;
  double lo, hi;
  if (lrange){
    printf("Please input the normal direction of the region [z]: ");
    if (strlen(gets(str)) > 0){
      if (strcmp(str,"x") == 0 || strcmp(str,"X")== 0) dir = 0;
      else if (strcmp(str,"y") == 0 || strcmp(str,"Y")== 0) dir = 1;
      else dir = 2;
    }
    lo = 0.; hi = latvec[dir][dir];
    printf("The lattice vector along the selected direction: %g %g %g\n", latvec[dir][0], latvec[dir][1], latvec[dir][2]);
    printf("Please input the lower and upper limit of the region [0 %lg]: ", latvec[dir][dir]);
    if (strlen(gets(str)) > 0){
      char *ptr;
      ptr = strtok(str," \t\n\r\f");  if (ptr != NULL) lo = atof(ptr);
      ptr = strtok(NULL," \t\n\r\f"); if (ptr != NULL) hi = atof(ptr);
    }
    printf("The desired region will be along the %d-th direction within [%g %g].\n", dir+1, lo, hi);
    printf("Will you limit the solid solution (0) inside or (1) outside the region?[0]: ");
    if (strlen(gets(str)) > 0) sscanf(str,"%d", &outside);
    
    nrange = 0;
    for (int i=0; i<natom; i++){
      if (outside){
        if (atpos[i][dir] < lo || atpos[i][dir] > hi) nrange++;
      } else {
        if (atpos[i][dir] >= lo && atpos[i][dir] <= hi) nrange++;
      }
    }
  }
  printf("Total number of atoms in the system: %d\n", natom);
  if (lrange) printf("Total number of atoms in the region: %d\n", nrange);
  printf("Total number of atomimc types      : %d\n", ntype);
  printf("Atomic type number for each type   :");
  for (int i=0; i<ntype; i++) printf(" %d", typeID[i]);
  printf("\nNumber of atoms for each  type     :");
  for (int i=0; i<ntype; i++) printf(" %d", numtype[i]); printf("\n");

  int ipsrc, idsrc, numsub, ipdes, iddes;
  do {
    printf("\nPlease input the atomic type to be substituted: ");
    scanf("%d", &ipsrc); while (getchar() != '\n');
    idsrc = lookup(ipsrc);
  } while (idsrc <0);
  printf("Total # of atoms with type %d is %d.\n", ipsrc, numtype[idsrc]);
  double frac;
  do {
    printf("Please input the fraction or total # of atoms to be replaced: ");
    scanf("%lg", &frac); while (getchar() != '\n');
  } while (frac<0. || int(frac) >= numtype[idsrc]);
  if (frac < 1.) numsub = int(frac*double(numtype[idsrc]));
  else numsub = int(frac);
  printf("There will be %d atoms with type %d to be replaced.\n", numsub, ipsrc);
  if (numsub < 1) return;

  do {
    printf("Please input the atomic type to be assigned: ");
    scanf("%d", &ipdes); while (getchar() != '\n');
    iddes = lookup(ipdes);
    if (iddes >= 0){
      printf("***Note: assigned type already exist, continue? (0= no, 1=yes):");
      int flag;
      scanf("%d", &flag); while (getchar() != '\n');
      if (flag == 1) iddes = -2;
    }
  } while (iddes >= 0);

  // create random number generator if not created; current time as seed;
  if (random == NULL){
    time_t ctime;
    time(&ctime);
    random = new RanPark((int) ctime);
  }

  // now to generate the solid solution
  int isub =0;
  while (isub < numsub){
    int i = int(random->uniform()*double(natom));
    if (attyp[i] == ipsrc){
      if (lrange){
        if (outside){
          if (atpos[i][dir] >= lo && atpos[i][dir] <= hi) continue;
        } else {
          if (atpos[i][dir] < lo || atpos[i][dir] > hi) continue;
        }
      }

      attyp[i] = ipdes;
      isub++;
    }
  }

  // reset type info
  typescan();
return;
}
