#include "driver.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <time.h>
#include <math.h>

#define MAXLINE 512
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
/* ----------------------------------------------------------------------
   constructor to initialize
------------------------------------------------------------------------- */
Driver::Driver()
{
  nx = ny = nz = natom = ntype = 0;

  alat = 0.;
  latt = NULL;
  atpos = NULL;
  attyp = NULL;
  random = NULL;
  typeID = numtype = NULL;
  xmap = ymap = zmap = umap = NULL;

  memory = new Memory();

  for (int i=0; i<3; i++)
  for (int j=0; j<3; j++) latvec[i][j] = 0.;

  MainMenu();

return;
}

int Driver::ShowMenu(const int flag)
{
  int ltype = 1;
  char str[MAXLINE];

  // print out the menu
  printf("\n"); for (int i=0; i<70; i++) printf("="); printf("\n");
  if (flag>0) printf("Please select the lattice type for lattice: %c\n", flag+'A'-1);
  else printf("Please select the lattice type of your system:\n");
  printf(" 1. FCC/NaCl/Diamond;          |  4. A3B;\n");
  printf(" 2. BCC;                       |  5. A2B;\n");
  printf(" 3. HCP/Graphene;              |  6. AB;\n");
  for (int i=0; i<31; i++) printf("-"); printf("+"); for (int i=0; i<38; i++) printf("-"); printf("\n");
  if (flag != 0){
    printf(" 7. User defined;              |  0. Exit.\n");
  } else {
    printf(" 7. User defined;              |  8. Multi-layer.\n");
    for (int i=0; i<70; i++) printf("-");
#ifdef Poly
    printf("\n 9. Polycrystal;               |  0. Exit.\n");
#else
    printf("\n 0. Exit.\n");
#endif
  }
  for (int i=0; i<70; i++) printf("-");
  printf("\nYour choice [1]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) ltype = atoi(strtok(str," \t\n\r\f"));
  printf("You selected: %d", ltype);
  printf("\n"); for (int i=0; i<70; i++) printf("="); printf("\n");

  switch (ltype){
  case 1: latt = new FCC(); break;
  case 2: latt = new BCC(); break;
  case 3: latt = new HCP(); break;
  case 4: latt = new A3B(); break;
  case 5: latt = new A2B(); break;
  case 6: latt = new AB(); break;
  case 7: latt = new USER(); break;
  case 8: if (flag == 0) FormLayers(); break;
#ifdef Poly
  case 9: if (flag == 0) PolyCrystal(); break;
#endif
  default: exit(1);}

  int rflag = 8 - ltype;
  if (flag == 0 && rflag > 0){
    // re-orient the lattice
    printf("Would you like to re-orient the unit cell to comply with LAMMPS? (y/n)[n]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      char *ptr = strtok(str," \t\n\r\f");
      if (strcmp(ptr,"y")==0 || strcmp(ptr,"Y")==0) latt->OrientLattice();
    }
  }

return rflag;
}

void Driver::MainMenu()
{
  if ( ShowMenu(0) > 0 ){
    name = memory->create(name, strlen(latt->name)+1, "driver->MainMenu:name");
    strcpy(name, latt->name);
    alat = latt->alat;
   
    latt->display();

    generate();
  }
return;
}

/* ----------------------------------------------------------------------
   deconstructor to free memory
------------------------------------------------------------------------- */
Driver::~Driver()
{
  memory->destroy(atpos);
  memory->destroy(attyp);
  memory->destroy(xmap);
  memory->destroy(ymap);
  memory->destroy(zmap);
  memory->destroy(umap);
  memory->destroy(name);

  if (latt != NULL) delete latt;

  memory->destroy(typeID);
  memory->destroy(numtype);
  if (random != NULL) delete random;
}

/* ----------------------------------------------------------------------
   method to generate atomic configurations
------------------------------------------------------------------------- */
void Driver::generate()
{
  char str[MAXLINE];
  int leading_dir = 1;
  printf("\n"); for (int i=0; i<70; i++) printf("="); printf("\n");
  while (1){
    printf("Please input the extensions in x, y, and z directions: ");
    if (count_words(fgets(str,MAXLINE,stdin)) < 3) continue;
    nx = atoi(strtok(str,  " \t\n\r\f"));
    ny = atoi(strtok(NULL, " \t\n\r\f"));
    nz = atoi(strtok(NULL, " \t\n\r\f"));
    natom = nx*ny*nz*latt->nucell;
    if (natom > 0) break;
  }

  printf("Your system would be %d x %d x %d with %d atoms.\n",nx,ny,nz,natom);
  printf("Please indicate which direction should goes fast(1:x; other: z)[1]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) leading_dir = atoi(strtok(str, " \t\n\r\f"));
  for (int i=0; i<70; i++) printf("="); printf("\n");

  atpos = memory->create(atpos, natom, 3, "driver->generate:atpos");
  attyp = memory->create(attyp, natom, "driver->generate:attyp");
  xmap = memory->create(xmap, natom, "driver->generate:xmap");
  ymap = memory->create(ymap, natom, "driver->generate:ymap");
  zmap = memory->create(zmap, natom, "driver->generate:zmap");
  umap = memory->create(umap, natom, "driver->generate:umap");

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
  if (typeID != NULL) memory->destroy(typeID);
  if (numtype!= NULL) memory->destroy(numtype);
  typeID  = memory->create(typeID,  typmax, "driver->typescan:typeID");
  numtype = memory->create(numtype, typmax, "driver->typescan:numtype");
  for (int i=0; i<typmax; i++) numtype[i] = 0;

  ntype = 0;
  // now to identify the total number of types
  for (int i=0; i<natom; i++){
    int id = lookup(attyp[i]);
    if (id < 0){
      if (ntype == typmax){
        typmax += 5;
        typeID  = memory->grow(typeID, typmax, "driver->typescan:typeID");
        numtype = memory->grow(numtype,typmax, "driver->typescan:numtype");
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
   method to reset the atomic type ID
------------------------------------------------------------------------- */
void Driver::ResetTypeID()
{
  char str[MAXLINE];
  printf("\n"); for (int i=0; i<70; i++) printf("=");
  printf("\nThere are %d atomic types in system, and their IDs are:\nIndex : ", ntype);
  for (int i=0; i<ntype; i++) printf("%4d", i+1); printf("\nTypeID: ");
  for (int i=0; i<ntype; i++) printf("%4d", typeID[i]);
  printf("\nPlease input the new type IDs in sequence, enter to skip: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0){
    int newID[ntype], num=0;
    char *ptr;
    ptr = strtok(str, " \t\n\r\f");
    while (ptr != NULL){
      newID[num] = atoi(ptr);
      if (++num >= ntype) break;
      ptr = strtok(NULL, " \t\n\r\f");
    }
    if (num == ntype){
      printf("\nThe new type IDs are:\nIndex : ");
      for (int i=0; i<ntype; i++) printf("%4d", i+1); printf("\nTypeID: ");
      for (int i=0; i<ntype; i++) printf("%4d", newID[i]);
      printf("\nIs this what you want? (y/n)[y]: ");
      int nw = count_words(fgets(str,MAXLINE,stdin));
      char *flag = strtok(str, " \t\n\r\f");
      if ( nw < 1 || ((strcmp(flag,"n") != 0) && (strcmp(flag,"N") != 0))){
        for (int i=0; i<natom; i++){
          int id = lookup(attyp[i]);
          attyp[i] = newID[id];
        }

        typescan();
      }
    }
  }
  for (int i=0; i<70; i++) printf("="); printf("\n");

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
  if (natom < 1) return;
  FILE *fp;
  char str[MAXLINE], *posfile, *mapfile, *lmpfile;
  int flag_lmp_data = 1;
  if (latvec[0][1]*latvec[0][1]+latvec[0][2]*latvec[0][2]+latvec[1][2]*latvec[1][2] > 1.e-6) flag_lmp_data = 0;

  printf("\n"); for (int i=0; i<70; i++) printf("="); printf("\n");
  printf("Please input the filename of the output xyz file [atomcfg.xyz]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0){
    int n = strlen(str) + 1;
    posfile = new char[n];
    strcpy(posfile, strtok(str," \t\n\r\f"));
  } else {
    posfile = new char[12];
    strcpy(posfile, "atomcfg.xyz");
  }
  if (flag_lmp_data){
    printf("Please input the filename of the lammps atomic file [data.pos]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      int n = strlen(str) + 1;
      lmpfile = new char[n];
      strcpy(lmpfile, strtok(str," \t\n\r\f"));
    } else {
      lmpfile = new char[9];
      strcpy(lmpfile, "data.pos");
    }
  }

  if (xmap){
    printf("Please input the filename of the output map file [map.in]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      int n = strlen(str) + 1;
      mapfile = new char[n];
      strcpy(mapfile, strtok(str," \t\n\r\f"));
    } else {
      mapfile = new char[7];
      strcpy(mapfile, "map.in");
    } 
  }

  printf("\nThe atomic configuration will be written to files %s", posfile);
  if (flag_lmp_data) printf(" and %s\n", lmpfile); else printf("\n");
  if (xmap) printf("The FFT map information  will be written to file: %s\n", mapfile);
  for (int i=0; i<70; i++) printf("="); printf("\n");

  // write the xyz position file 
  fp = fopen(posfile, "w");
  fprintf(fp, "%d\n", natom);
  fprintf(fp, "%s cell with dimension %dx%dx%d and a= %g\n", name, nx, ny, nz, alat);
  int nr = 3;
  if (natom < nr) nr = natom;
  for (int i=0; i<nr; i++){
    fprintf(fp,"%d %16.16e %16.16e %16.16e crystal_vector %d %16.16e %16.16e %16.16e\n", attyp[i], atpos[i][0],
    atpos[i][1], atpos[i][2], i+1, latvec[i][0], latvec[i][1], latvec[i][2]);
  }
  for (int i=nr; i<natom; i++)
    fprintf(fp,"%d %16.16e %16.16e %16.16e\n", attyp[i], atpos[i][0], atpos[i][1], atpos[i][2]);
  fclose(fp);
  delete []posfile;

  // write the lammps atomic style file
  if (flag_lmp_data){
    fp = fopen(lmpfile,"w");
    fprintf(fp, "# %s cell with dimension %d x %d x %d and a = %g\n", name, nx, ny, nz, alat);
    fprintf(fp, "%10d  atoms\n", natom);
    fprintf(fp, "%10d  atom types\n\n", ntype);
    fprintf(fp, " 0. %20.14f  xlo xhi\n", latvec[0][0]);
    fprintf(fp, " 0. %20.14f  ylo yhi\n", latvec[1][1]);
    fprintf(fp, " 0. %20.14f  zlo zhi\n", latvec[2][2]);
    if ( latvec[1][0]*latvec[1][0] + latvec[2][0]*latvec[2][0] + latvec[2][1]*latvec[2][1] > 1.e-8 )
      fprintf(fp, "%20.14f %20.14f %20.14f xy xz yz\n", latvec[1][0], latvec[2][0], latvec[2][1]);
    fprintf(fp, "\nAtoms\n\n");
  
    for (int i=0; i<natom; i++) fprintf(fp,"%d %d %20.14f %20.14f %20.14f\n", i+1, attyp[i], atpos[i][0], atpos[i][1], atpos[i][2]);
    fclose(fp);
    delete []lmpfile;
  }

  // write the map file, useful to fix_phonon only.
  if (xmap){
     fp = fopen(mapfile, "w");
     fprintf(fp,"%d %d %d %d\n", nx, ny, nz, latt->nucell);
     fprintf(fp,"Map file for %dx%dx%d %s cell.\n",nx,ny,nz, name);
     for (int i=0; i<natom; i++)
       fprintf(fp,"%d %d %d %d %d\n", xmap[i], ymap[i], zmap[i], umap[i], i+1);
     fclose(fp);

     delete []mapfile;
  }

return;
}

/* ----------------------------------------------------------------------
   method to modify the resultant model
------------------------------------------------------------------------- */
void Driver::modify()
{
  if (natom < 1) return;
  char str[MAXLINE];
  int ncycle = 1;
  while (ncycle){
    int job=0;
    // to display the menu for modification
    printf("\n"); for (int i=0; i<70; i++) printf("="); printf("\n");
    printf("Please select the modification you want to do:\n");
    printf("  1. Create substitutional solid solution;\n");
    printf("  2. Reset atomic types;\n");

    if (ncycle == 1) printf("  0. Nothing.\n");
    else printf("  0. Done.\n");
    printf("Your choice [0]: ");

    if (count_words(fgets(str,MAXLINE,stdin)) >0) job = atoi(strtok(str, " \t\n\r\f"));

    printf("You selected: %d", job);
    printf("\n"); for (int i=0; i<70; i++) printf("="); printf("\n");
   
    switch (job){ 
    case 1: solidsol(); break;
    case 2: ResetTypeID(); break;
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
  char str[MAXLINE];
  printf("\n"); for (int i=0; i<70; i++) printf("="); printf("\n");
  int lrange = 0, nrange;
  printf("Limit the solid solution within a region? (y/n)[n]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0){
    char *flag = strtok(str, " \r\n\t\f");
    if (strcmp(flag,"y") == 0 || strcmp(flag,"Y")== 0) lrange = 1;
  }
  int dir = 2, outside = 0;
  double lo, hi;
  if (lrange){
    printf("Please input the normal direction of the region [z]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      char *flag = strtok(str," \t\n\r\f");
      if (strcmp(flag,"x") == 0 || strcmp(flag,"X")== 0) dir = 0;
      else if (strcmp(flag,"y") == 0 || strcmp(flag,"Y")== 0) dir = 1;
      else dir = 2;
    }
    lo = 0.; hi = latvec[dir][dir];
    printf("The lattice vector along the selected direction: %g %g %g\n", latvec[dir][0], latvec[dir][1], latvec[dir][2]);
    printf("Please input the lower and upper limit of the region [0 %lg]: ", latvec[dir][dir]);
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      char *ptr;
      ptr = strtok(str," \t\n\r\f");  if (ptr != NULL) lo = atof(ptr);
      ptr = strtok(NULL," \t\n\r\f"); if (ptr != NULL) hi = atof(ptr);
    }
    printf("The desired region will be along the %d-th direction within [%g %g].\n", dir+1, lo, hi);
    printf("Will you limit the solid solution (0) inside or (1) outside the region?[0]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0) outside = atoi(strtok(str, " \t\n\r\f"));
    
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
  while (1){
    printf("\nPlease input the atomic type to be substituted: ");
    while (count_words(fgets(str,MAXLINE,stdin)) < 1);
    ipsrc = atoi(strtok(str, " \t\n\r\f"));
    idsrc = lookup(ipsrc);

    if (idsrc >= 0) break;
  }
  printf("Total # of atoms with type %d is %d.\n", ipsrc, numtype[idsrc]);
  double frac;
  while (1){
    printf("Please input the fraction or total # of atoms to be replaced: ");
    while (count_words(fgets(str,MAXLINE,stdin)) < 1);
    frac = atof(strtok(str, " \t\n\r\f"));

    if (frac >= 0. && int(frac) <= numtype[idsrc]) break;
  }
  if (frac < 1.) numsub = int(frac*double(numtype[idsrc]));
  else numsub = int(frac);
  printf("There will be %d atoms with type %d to be replaced.\n", numsub, ipsrc);
  if (numsub < 1) return;

  while (1){
    printf("Please input the atomic type to be assigned: ");
    while (count_words(fgets(str,MAXLINE,stdin)) < 1);
    ipdes = atoi(strtok(str, " \t\n\r\f"));
    iddes = lookup(ipdes);
    if (iddes >= 0){
      printf("***Note: assigned type already exist, continue? (0= no, 1=yes)[1]:");
      if (count_words(fgets(str,MAXLINE,stdin)) > 0){
        int flag = atoi(strtok(str, " \t\n\r\f"));
        if (flag == 1) iddes = -2;
      }
    }

    if (iddes < 0) break;
  }

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

/*------------------------------------------------------------------------------
 * Method to form multilayers
 *----------------------------------------------------------------------------*/
void Driver::FormLayers()
{
  int nlat = 1, idum;
  int *mynx, *myny;
  char str[MAXLINE];

  printf("\n\n>>>>>>======  To form multilayers with multiple lattices  ======<<<<<<\n");
  printf("NOTE: The 3rd axis of these lattices must be perpendicular to the other 2!\n");
  printf("\nPlease input the number of lattices in your multi-layer system: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) nlat = atoi(strtok(str," \t\n\r\f"));
  if (nlat < 1) return;

  lattice *latts[nlat];
  for (int i=0; i<nlat; i++){
    if (! ShowMenu(i+1) ){nlat = i; break;}

    latts[i] = latt;
    latt = NULL;
  }

  for (int i=0; i<nlat; i++){
    printf("\n>>>>>>   Lattice info for lattice: %c - %s    <<<<<", 'A'+i, latts[i]->name);
    if (latts[i]->perp_z == 0)
      printf("\nWARNING: A3 is not perpendicular to A1 and A2, this lattice cannot be used to form layers!\n");

    latts[i]->display();
  }

  idum = 0;
  printf("\nYou have defined %d lattices: ", nlat);
  for (int i=0; i<nlat; i++) printf(" %c = %s;", 'A'+i, latts[i]->name); printf("\n");
  
  mynx = memory->create(mynx, nlat, "mynx");
  myny = memory->create(myny, nlat, "myny");
  for (int ilat=0; ilat<nlat; ilat++){
    printf("Please input the lateral extensions (nx & ny) for lattice %c: ", 'A'+ilat);
    while (1){
      if ( count_words(fgets(str,MAXLINE,stdin)) == 2 ){
        mynx[ilat] = atoi(strtok(str, " \t\n\r\f"));
        myny[ilat] = atoi(strtok(NULL," \t\n\r\f"));
       if (mynx[ilat] > 0 && myny[ilat] > 0) break;
      }
    }
  }

  nx = ny = nz = 0;
  for (int j=0; j<3; j++){
    latvec[0][j] = latvec[1][j] = 0.;
  }

  printf("\nThe surface vectors for each lattice will be:\n");
  for (int i=0; i<nlat; i++){
    printf("  %c: [", i+'A');
    for (int j=0; j<3; j++){
      double xi = mynx[i]*latts[i]->latvec[0][j]*latts[i]->alat;
      printf("%g ", xi);
      latvec[0][j] += xi;
    }
    printf("] [");
    for (int j=0; j<3; j++){
      double yi = myny[i]*latts[i]->latvec[1][j]*latts[i]->alat;
      printf("%g ", yi);
      latvec[1][j] += yi;
    }
    printf("]\n");
  }
  for (int j=0; j<3; j++){
    latvec[0][j] /= double(nlat);
    latvec[1][j] /= double(nlat);
  }
  printf("Please input your desired surface vectors [%g %g %g, %g %g %g]: ",
    latvec[0][0], latvec[0][1], latvec[0][2], latvec[1][0], latvec[1][1], latvec[1][2]);
  if ( count_words(fgets(str,MAXLINE,stdin)) == 6 ){
    char *ptr = strtok(str," \n\r\t\f");
    for (int i=0; i<2; i++)
    for (int j=0; j<3; j++){ latvec[i][j] = atof(ptr); ptr = strtok(NULL, " \n\r\t\f");}
  }

  double lx, ly, lx0=0., ly0=0.;
  for (int j=0; j<3; j++){
    lx0 += latvec[0][j]*latvec[0][j];
    ly0 += latvec[1][j]*latvec[1][j];
  }
  lx0 = sqrt(lx0); ly0 = sqrt(ly0);

  for (int i=0; i<nlat; i++){
    lx = ly = 0.;
    for (int j=0; j<3; j++){
      double xi = mynx[i]*latts[i]->latvec[0][j]*latts[i]->alat;
      double yi = myny[i]*latts[i]->latvec[1][j]*latts[i]->alat;
      lx += xi*xi;
      ly += yi*yi;
    }
    lx = sqrt(lx); ly = sqrt(ly);
    printf("Lateral misfit for lattice %c is: [%lg %lg]\n", 'A'+i, (lx-lx0)/lx0, (ly-ly0)/ly0);
  }

  double H = 0.;
  int zprev[nlat], ntprev[nlat], zflag = 0, iatom = 0;
  for (int i=0; i<nlat; i++) zprev[i] = 0;

  ntprev[0] = 0;
  for (int i=1; i<nlat; i++) ntprev[i] = ntprev[i-1] + latts[i-1]->ntype;

  char realized[MAXLINE]; strcpy(realized, "");

  int first = 1, flag_no_interlayer = 0;
  double Hlast = 0., Hfirst = 0., Hextra;

  printf("\nPlease input the layer sequences, for example, if you have two lattices: A and B,\n");
  printf("and you want to have 4 layers A, 5 layers B and then 3 layers A, input A4 B5 A3.\n");
  printf("If extra distance between different lattices is needed, just insert a number between\n");
  printf("them, for example: A4 0.5 B5 0.4 A4 B5...; multiple numbers will add multiple distances.\n");
  printf("If you want to form the 2nd A layers from its first layer in the unit cell, use\n");
  printf("lower case 'a' instead of 'A'. Now, input your sequences: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 1) {
    char *ptr = strtok(str," \n\r\t\f");
    while (ptr != NULL){
      zflag = 0; Hextra = 0.;
      int ilat = ptr[0] - 'A';
      if (ilat > nlat){ ilat = ptr[0] - 'a'; zflag = 1;}
      if (ilat >=0 && ilat < nlat){
        latt = latts[ilat];

        strcat(realized," ");strcat(realized, ptr);
        if (zflag) zprev[ilat] = 0;
        ptr[0] = ' ';
        int nl_new  = atoi(ptr);
        int ntm_new = 0;
        for (int i=0; i<nl_new; i++) ntm_new += latt->numlayer[(i+zprev[ilat])%latt->nlayer];

        double Hbelow = latt->h[(latt->nlayer-1+zprev[ilat])%latt->nlayer];
        if (first){
          first = 0;
          Hlast = 0.;
          Hfirst = Hbelow;
          H = -Hfirst;
        }
        if (nl_new > 0){
          if (flag_no_interlayer) flag_no_interlayer = 0;
          else H += MAX(Hlast,Hbelow);
        }

        ntm_new *= (mynx[ilat]*myny[ilat]);
        natom += ntm_new;
        atpos = memory->grow(atpos, natom, 3, "atpos");
        attyp = memory->grow(attyp, natom, "attyp");

        for (int k=0; k<nl_new; k++){
          int il = (k+zprev[ilat])%latt->nlayer;
          for (int ia=0; ia<latt->nucell; ia++){
            if ( latt->layer[ia] == il ){
              for (int i=0; i<mynx[ilat]; i++)
              for (int j=0; j<myny[ilat]; j++){
                atpos[iatom][0] = (double(i)+latt->atpos[ia][0])/double(mynx[ilat]);
                atpos[iatom][1] = (double(j)+latt->atpos[ia][1])/double(myny[ilat]);
                atpos[iatom][2] = H;
                attyp[iatom++]  = latt->attyp[ia] + ntprev[ilat];
              }
            }
          }
          nz++;

          Hlast = latt->h[il];
          H += Hlast;
        }
        if (nl_new > 0) H -= Hlast;
        zprev[ilat] += nl_new%latt->nlayer;

      } else if (strcmp(ptr, "-z") == 0){
        flag_no_interlayer = 1;
        strcat(realized," ");strcat(realized, ptr);

      } else {
        Hextra = atof(ptr);
        strcat(realized," ");strcat(realized, ptr);
        H += Hextra;
      }

      ptr = strtok(NULL, " \n\r\t\f");
    }
  }

  printf("\nThe layer sequences realized is: %s\n", realized);
  printf("In total, %d layers and %d atoms are created.\n", nz, iatom);

  if (flag_no_interlayer == 0) H += MAX(Hlast,Hfirst);
  latt = NULL; alat = 1.;
  latvec[2][0] = latvec[2][1] = 0.; latvec[2][2] = H;

  strcpy(str,"Multilayer: ");
  for (int i=0; i<nlat; i++){
    char info[MAXLINE];
    sprintf(info, "%dx%d-%s ", mynx[i], myny[i], latts[i]->name);
    strcat(str,info);
  }
  name = memory->create(name, strlen(str)+1, "driver->FormLayers"); strcpy(name, str);

  double tmp[2];
  for (int i=0; i<natom; i++){
    for (int idim=0; idim<2; idim++) tmp[idim] = atpos[i][idim];
    atpos[i][0] = tmp[0]*latvec[0][0] + tmp[1]*latvec[1][0];
    atpos[i][1] = tmp[0]*latvec[0][1] + tmp[1]*latvec[1][1];
  }

  for (int i=0; i<nlat; i++) delete latts[i];

  memory->destroy(mynx);
  memory->destroy(myny);
  // find the total # of types and # of atoms for each type
  typescan();

return;
}

/*------------------------------------------------------------------------------
 * Method to count # of words in a string, without destroying the string
 *----------------------------------------------------------------------------*/
int Driver::count_words(const char *line)
{
  int n = strlen(line) + 1;
  char *copy = (char *) memory->smalloc(n*sizeof(char),"copy");
  strcpy(copy,line);

  char *ptr;
  if (ptr = strchr(copy,'#')) *ptr = '\0';

  if (strtok(copy," \t\n\r\f") == NULL) {
    memory->sfree(copy);
    return 0;
  }
  n = 1;
  while (strtok(NULL," \t\n\r\f")) n++;

  memory->sfree(copy);
  return n;
}