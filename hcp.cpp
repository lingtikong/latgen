#include "hcp.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#define MAXLINE 256

using namespace std;

/* -----------------------------------------------------------------------------
 * To select the orientation of the lattice
 * -------------------------------------------------------------------------- */
HCP::HCP() : lattice()
{
  char str[MAXLINE];
  alat = 1.; ca = sqrt(8./3.);
  // print out the menu
  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");
  printf("Please input the lattice constant of the HCP lattice [1]:");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) alat = numeric(strtok(str, " \t\n\r\f"));
  if (alat <= 0.) alat = 1.;

  printf("Please input the value of c/a ratio or c (negative) [1.633]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) ca = numeric(strtok(str, " \t\n\r\f"));
  if (ca == 0.) ca = sqrt(8./3.);
  if (ca <  0.) ca = -ca/alat;
  printf("The lattice constants of your HCP: a = %g, c/a = %g.\n\n", alat, ca);

  int orient = 1;
  printf("Please select the orientation of the HCP lattice:\n");
  printf("   1. (001)/(0001);         5. Graphene;\n");
  printf("   2. (100)/(2-1-10);       6. Graphite;\n");
  printf("   3. (110)/(11-20);        7. (101)/(2-1-13);\n");
  printf("   4.(-110)/(-1100);        8. (112)/(11-26);\n");
  printf("Your choice [%d]: ", orient);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) orient = inumeric(strtok(str, " \t\n\r\f"));
  printf("Your selection : %d", orient);
  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");
  
  // initialize according to orientation
  initialized = 0;
  switch (orient){
  case 1:
    HCP001();
    break;
  case 2:
    HCP100();
    break;
  case 3:
    HCP110();
    break;
  case 4:
    HCPm10();
    break;
  case 5:
    Graphene();
    break;
  case 6:
    Graphite();
    break;
  case 7:
    HCP10m1();
    break;
  case 8:
    HCP112();
    break;
  default:
    break;
  }

}

/* -----------------------------------------------------------------------------
 * Deconstructor does nothing
 * -------------------------------------------------------------------------- */
HCP::~HCP()
{

}

/* -----------------------------------------------------------------------------
 * Initialize for (001) orientation
 * -------------------------------------------------------------------------- */
void HCP::HCP001()
{
  char str[MAXLINE];
  int surftype = 5;
  // print out the menu
  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");
  printf("Please select the type of HCP(001) surface:\n");
  printf("   1. U along x,  60 degree;\n");
  printf("   2. V along y,  60 degree;\n");
  printf("   3. U along x, 120 degree;\n");
  printf("   4. V along y, 120 degree;\n");
  printf("   5. Rectangle, long along x;\n");
  printf("   6. Rectangle, long along y;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = inumeric(strtok(str, " \t\n\r\f"));
  printf("Your selection : %d", surftype);
  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");

  for (int i = 0; i < 3; ++i)
  for (int j = 0; j < 3; ++j) latvec[i][j] = 0.;

  memory->create(name,9,"HCP:name");
  strcpy(name, "HCP(001)");

  // initialize according to surface type
  switch (surftype){
  case 1:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] = 1.;
    latvec[1][0] = 0.5;
    latvec[1][1] = sqrt(0.75);
    latvec[2][2] = ca;

    memory->create(atpos, nucell, 3, "HCP001_atpos");
    memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.5;

    initialized = 1;
    break;
  case 2:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] = sqrt(0.75);
    latvec[0][1] = 0.5;
    latvec[1][1] = 1.;
    latvec[2][2] = ca;

    memory->create(atpos, nucell, 3, "HCP001_atpos");
    memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.5;

    initialized = 1;
    break;
  case 3:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] =  1.;
    latvec[1][0] = -0.5;
    latvec[1][1] =  sqrt(0.75);
    latvec[2][2] = ca;

    memory->create(atpos, nucell, 3, "HCP001_atpos");
    memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 2./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.5;

    initialized = 1;
    break;
  case 4:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] = sqrt(0.75);
    latvec[0][1] = -0.5;
    latvec[1][1] = 1.;
    latvec[2][2] = ca;

    memory->create(atpos, nucell, 3, "HCP001_atpos");
    memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 2./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.5;

    initialized = 1;
    break;
  case 5:
    nucell = 4;
    ntype  = 1;

    latvec[0][0] = sqrt(3.0);
    latvec[1][1] = 1.;
    latvec[2][2] = ca;

    memory->create(atpos, nucell, 3, "HCP001_atpos");
    memory->create(attyp, nucell, "HCP:attyp");

    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;

    atpos[2][0] = 1./6.;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.5;

    atpos[3][0] = 2./3.;
    atpos[3][1] = 0.0;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  case 6:
    nucell = 4;
    ntype  = 1;

    latvec[0][0] = 1.;
    latvec[1][1] = sqrt(3.0);
    latvec[2][2] = ca;

    memory->create(atpos, nucell, 3, "HCP001_atpos");
    memory->create(attyp, nucell, "HCP:attyp");

    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;

    atpos[2][0] = 0.5;
    atpos[2][1] = 1./6.;
    atpos[2][2] = 0.5;

    atpos[3][0] = 0.0;
    atpos[3][1] = 2./3.;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* -----------------------------------------------------------------------------
 * Initialize for (110) orientation
 * -------------------------------------------------------------------------- */
void HCP::HCP100()
{
  for (int i = 0; i < 3; ++i)
  for (int j = 0; j < 3; ++j) latvec[i][j] = 0.;

  memory->create(name,9,"HCP:name");
  strcpy(name, "HCP(100)");

  // initialize
  nucell = 4;
  ntype  = 1;
    
  latvec[0][0] = 1.;
  latvec[1][1] = ca;
  latvec[2][2] = sqrt(3.);

  memory->create(atpos, nucell, 3, "HCP100_atpos");
  memory->create(attyp, nucell, "HCP:attyp");
    
  for (int i = 0; i < nucell; ++i) attyp[i] = 1;
  atpos[0][0] = 0.50;
  atpos[0][1] = 0.75;
  atpos[0][2] = 0.00;

  atpos[1][0] = 0.50;
  atpos[1][1] = 0.25;
  atpos[1][2] = 1./3.;

  atpos[2][0] = 0.00;
  atpos[2][1] = 0.75;
  atpos[2][2] = 0.50;

  atpos[3][0] = 0.00;
  atpos[3][1] = 0.25;
  atpos[3][2] = 5./6.;

  initialized = 1;

return;
}

/* -----------------------------------------------------------------------------
 * Initialize for (110) orientation
 * -------------------------------------------------------------------------- */
void HCP::HCP110()
{
  for (int i = 0; i < 3; ++i)
  for (int j = 0; j < 3; ++j) latvec[i][j] = 0.;

  memory->create(name,9,"HCP:name");
  strcpy(name, "HCP(110)");

  // initialize
  nucell = 4;
  ntype  = 1;
    
  latvec[0][0] = ca;
  latvec[1][1] = sqrt(3.);
  latvec[2][2] = 1.;

  memory->create(atpos, nucell, 3, "HCP110_atpos");
  memory->create(attyp, nucell, "HCP:attyp");
    
  for (int i = 0; i < nucell; ++i) attyp[i] = 1;
  atpos[0][0] = 0.25;
  atpos[0][1] = 1./3.;
  atpos[0][2] = 0.00;

  atpos[1][0] = 0.75;
  atpos[1][1] = 2./3.;
  atpos[1][2] = 0.;

  atpos[2][0] = 0.75;
  atpos[2][1] = 1./6.;
  atpos[2][2] = 0.50;

  atpos[3][0] = 0.25;
  atpos[3][1] = 5./6.;
  atpos[3][2] = 0.50;

  initialized = 1;

return;
}

/* -----------------------------------------------------------------------------
 * Initialize for (-110) orientation
 * -------------------------------------------------------------------------- */
void HCP::HCPm10()
{
  for (int i = 0; i < 3; ++i)
  for (int j = 0; j < 3; ++j) latvec[i][j] = 0.;

  memory->create(name,10,"HCP:name");
  strcpy(name, "HCP(-110)");

  // initialize
  nucell = 4;
  ntype  = 1;
    
  latvec[0][0] = 1.;
  latvec[1][1] = ca;
  latvec[2][2] = sqrt(3.);

  memory->create(atpos, nucell, 3, "HCPm10_atpos");
  memory->create(attyp, nucell, "HCP:attyp");
    
  for (int i = 0; i < nucell; ++i) attyp[i] = 1;
  atpos[0][0] = 0.50;
  atpos[0][1] = 0.25;
  atpos[0][2] = 0.00;

  atpos[1][0] = 0.00;
  atpos[1][1] = 0.75;
  atpos[1][2] = 1./6.;

  atpos[2][0] = 0.00;
  atpos[2][1] = 0.25;
  atpos[2][2] = 0.50;

  atpos[3][0] = 0.50;
  atpos[3][1] = 0.75;
  atpos[3][2] = 2./3.;

  initialized = 1;

return;
}

/* -----------------------------------------------------------------------------
 * Graphene with different orientations
 * -------------------------------------------------------------------------- */
void HCP::Graphene()
{
  char str[MAXLINE];
  int surftype = 5;
  // print out the menu
  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");
  printf("   1. Primitive, U along x,  60 degree;\n");
  printf("   2. Primitive, V along y,  60 degree;\n");
  printf("   3. Primitive, U along x, 120 degree;\n");
  printf("   4. Primitive, V along y, 120 degree;\n");
  printf("   5. Rectangle, long along x;\n");
  printf("   6. Rectangle, long along y;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = inumeric(strtok(str, " \t\n\r\f"));
  printf("Your selection : %d", surftype);
  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");

  for (int i = 0; i < 3; ++i)
  for (int j = 0; j < 3; ++j) latvec[i][j] = 0.;

  memory->create(name,9,"HCP:name");
  strcpy(name, "Graphene");

  // initialize according to surface type
  switch (surftype){
  case 1:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] = 1.;
    latvec[1][0] = 0.5;
    latvec[1][1] = sqrt(0.75);
    latvec[2][2] = ca;

    memory->create(atpos, nucell, 3, "Graphene_atpos");
    memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.;

    initialized = 1;
    break;
  case 2:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] = sqrt(0.75);
    latvec[0][1] = 0.5;
    latvec[1][1] = 1.;
    latvec[2][2] = ca;

    memory->create(atpos, nucell, 3, "Graphene_atpos");
    memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.;

    initialized = 1;
    break;
  case 3:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] =  1.;
    latvec[1][0] = -0.5;
    latvec[1][1] =  sqrt(0.75);
    latvec[2][2] = ca;

    memory->create(atpos, nucell, 3, "Graphene_atpos");
    memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 2./3.;
    atpos[1][2] = 0.;

    initialized = 1;
    break;
  case 4:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] = sqrt(0.75);
    latvec[0][1] = -0.5;
    latvec[1][1] = 1.;
    latvec[2][2] = ca;

    memory->create(atpos, nucell, 3, "Graphene_atpos");
    memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 2./3.;
    atpos[1][2] = 0.;

    initialized = 1;
    break;
  case 5:
    nucell = 4;
    ntype  = 1;

    latvec[0][0] = sqrt(3.0);
    latvec[1][1] = 1.;
    latvec[2][2] = ca;

    memory->create(atpos, nucell, 3, "Graphene_atpos");
    memory->create(attyp, nucell, "HCP:attyp");

    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 0.;
    atpos[1][2] = 0.;

    atpos[2][0] = 0.5;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.;

    atpos[3][0] = 5./6.;
    atpos[3][1] = 0.5;
    atpos[3][2] = 0.;

    initialized = 1;
    break;
  case 6:
    nucell = 4;
    ntype  = 1;

    latvec[0][0] = 1.;
    latvec[1][1] = sqrt(3.0);
    latvec[2][2] = ca;

    memory->create(atpos, nucell, 3, "Graphene_atpos");
    memory->create(attyp, nucell, "HCP:attyp");

    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 0.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.;

    atpos[2][0] = 0.5;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.;

    atpos[3][0] = 0.5;
    atpos[3][1] = 5./6.;
    atpos[3][2] = 0.;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* -----------------------------------------------------------------------------
 * Graphite with different orientations
 * -------------------------------------------------------------------------- */
void HCP::Graphite()
{
  char str[MAXLINE];
  int surftype = 5;
  // print out the menu
  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");
  printf("   1. Primitive, U along x,  60 degree;\n");
  printf("   2. Primitive, V along y,  60 degree;\n");
  printf("   3. Primitive, U along x, 120 degree;\n");
  printf("   4. Primitive, V along y, 120 degree;\n");
  printf("   5. Rectangle, long along x;\n");
  printf("   6. Rectangle, long along y;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = inumeric(strtok(str, " \t\n\r\f"));
  printf("Your selection : %d", surftype);
  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");

  for (int i = 0; i < 3; ++i)
  for (int j = 0; j < 3; ++j) latvec[i][j] = 0.;

  memory->create(name,9,"HCP:name");
  strcpy(name, "Graphite");

  if (fabs(ca-2.725) > 0.1){
    printf("\nThe experimental c/a for graphite is ~2.725, while yours is %g, if you\n", ca);
    printf("want to redefine it, input now, enter to keep c/a = %g: ", ca);
    if (count_words(fgets(str,MAXLINE,stdin)) > 0) ca = numeric(strtok(str, " \t\n\r\f"));
    printf("The adopted c/a will be: %g\n\n", ca);
  }

  // initialize according to surface type
  switch (surftype){
  case 1:
    nucell = 4;
    ntype  = 1;
    
    latvec[0][0] = 1.;
    latvec[1][0] = 0.5;
    latvec[1][1] = sqrt(0.75);
    latvec[2][2] = ca;

    memory->create(atpos, nucell, 3, "Graphite_atpos");
    memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.;

    atpos[2][0] = 0.;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.5;

    atpos[3][0] = 2./3.;
    atpos[3][1] = 2./3.;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  case 2:
    nucell = 4;
    ntype  = 1;
    
    latvec[0][0] = sqrt(0.75);
    latvec[0][1] = 0.5;
    latvec[1][1] = 1.;
    latvec[2][2] = ca;

    memory->create(atpos, nucell, 3, "Graphite_atpos");
    memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.;

    atpos[2][0] = 0.;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.5;

    atpos[3][0] = 2./3.;
    atpos[3][1] = 2./3.;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  case 3:
    nucell = 4;
    ntype  = 1;
    
    latvec[0][0] =  1.;
    latvec[1][0] = -0.5;
    latvec[1][1] =  sqrt(0.75);
    latvec[2][2] = ca;

    memory->create(atpos, nucell, 3, "Graphite_atpos");
    memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 2./3.;
    atpos[1][2] = 0.;

    atpos[2][0] = 0.;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.5;

    atpos[3][0] = 2./3.;
    atpos[3][1] = 1./3.;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  case 4:
    nucell = 4;
    ntype  = 1;
    
    latvec[0][0] = sqrt(0.75);
    latvec[0][1] = -0.5;
    latvec[1][1] = 1.;
    latvec[2][2] = ca;

    memory->create(atpos, nucell, 3, "Graphite_atpos");
    memory->create(attyp, nucell, "HCP:attyp");
    
    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 2./3.;
    atpos[1][2] = 0.;

    atpos[2][0] = 0.;
    atpos[2][1] = 0.;
    atpos[2][2] = 0.5;

    atpos[3][0] = 2./3.;
    atpos[3][1] = 1./3.;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  case 5:
    nucell = 8;
    ntype  = 1;

    latvec[0][0] = sqrt(3.0);
    latvec[1][1] = 1.;
    latvec[2][2] = ca;

    memory->create(atpos, nucell, 3, "Graphite_atpos");
    memory->create(attyp, nucell, "HCP:attyp");

    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 0.;
    atpos[1][2] = 0.;

    atpos[2][0] = 0.5;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.;

    atpos[3][0] = 5./6.;
    atpos[3][1] = 0.5;
    atpos[3][2] = 0.;

    atpos[4][0] = 0.;
    atpos[4][1] = 0.;
    atpos[4][2] = 0.5;

    atpos[5][0] = 2./3.;
    atpos[5][1] = 0.;
    atpos[5][2] = 0.5;

    atpos[6][0] = 0.5;
    atpos[6][1] = 0.5;
    atpos[6][2] = 0.5;

    atpos[7][0] = 1./6.;
    atpos[7][1] = 0.5;
    atpos[7][2] = 0.5;

    initialized = 1;
    break;
  case 6:
    nucell = 8;
    ntype  = 1;

    latvec[0][0] = 1.;
    latvec[1][1] = sqrt(3.0);
    latvec[2][2] = ca;

    memory->create(atpos, nucell, 3, "Graphite_atpos");
    memory->create(attyp, nucell, "HCP:attyp");

    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 0.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 0.;

    atpos[2][0] = 0.5;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.;

    atpos[3][0] = 0.5;
    atpos[3][1] = 5./6.;
    atpos[3][2] = 0.;

    atpos[4][0] = 0.;
    atpos[4][1] = 0.;
    atpos[4][2] = 0.5;

    atpos[5][0] = 0.;
    atpos[5][1] = 2./3.;
    atpos[5][2] = 0.5;

    atpos[6][0] = 0.5;
    atpos[6][1] = 0.5;
    atpos[6][2] = 0.5;

    atpos[7][0] = 0.5;
    atpos[7][1] = 1./6.;
    atpos[7][2] = 0.5;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* -----------------------------------------------------------------------------
 * HCP (10 -1)
 * -------------------------------------------------------------------------- */
void HCP::HCP10m1()
{
  char str[MAXLINE];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");
  printf("   1. U = [010]/[-12-10] along x, V = [101]/[2-1-1-3];\n");
  printf("   2. U = [010]/[-12-10], V = [101]/[3-1-13] along y;\n");
  printf("   3. U = [212]/[10-1-2] along x, V = [010/[-12-10] along y;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = inumeric(strtok(str, " \t\n\r\f"));
  printf("Your selection : %d", surftype);
  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");

  for (int i = 0; i < 3; ++i)
  for (int j = 0; j < 3; ++j) latvec[i][j] = 0.;

  memory->create(name,10,"HCP:name");
  strcpy(name, "HCP(10-1)");
 
  // initialize according to surface type
  double lb = sqrt(1. + ca*ca);
  double c_ang = sqrt((ca*ca + 0.75)/(ca*ca + 1.));
  switch (surftype){
  case 1:
    latvec[0][0] = 1.;
    latvec[1][0] = -0.5;
    latvec[1][1] = sqrt(ca*ca + 0.75);

    latvec[2][2] = 41./sqrt(4./3. + 1./(ca*ca));

    nucell = 82;
    ntype  = 1;
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
  
    attyp[ 0] =  1;
    atpos[ 0][0] =     0.0000000000;
    atpos[ 0][1] =     0.0000000000;
    atpos[ 0][2] =     0.0000000000;
    
    attyp[ 1] =  1;
    atpos[ 1][0] =     0.3414600000;
    atpos[ 1][1] =     0.6829300000;
    atpos[ 1][2] =     0.0203300000;
    
    attyp[ 2] =  1;
    atpos[ 2][0] =     0.6097600000;
    atpos[ 2][1] =     0.2195100000;
    atpos[ 2][2] =     0.0243900000;
    
    attyp[ 3] =  1;
    atpos[ 3][0] =     0.9512200000;
    atpos[ 3][1] =     0.9024400000;
    atpos[ 3][2] =     0.0447200000;
    
    attyp[ 4] =  1;
    atpos[ 4][0] =     0.2195100000;
    atpos[ 4][1] =     0.4390200000;
    atpos[ 4][2] =     0.0487800000;
    
    attyp[ 5] =  1;
    atpos[ 5][0] =     0.5609800000;
    atpos[ 5][1] =     0.1219500000;
    atpos[ 5][2] =     0.0691100000;
    
    attyp[ 6] =  1;
    atpos[ 6][0] =     0.8292700000;
    atpos[ 6][1] =     0.6585400000;
    atpos[ 6][2] =     0.0731700000;
    
    attyp[ 7] =  1;
    atpos[ 7][0] =     0.1707300000;
    atpos[ 7][1] =     0.3414600000;
    atpos[ 7][2] =     0.0935000000;
    
    attyp[ 8] =  1;
    atpos[ 8][0] =     0.4390200000;
    atpos[ 8][1] =     0.8780500000;
    atpos[ 8][2] =     0.0975600000;
    
    attyp[ 9] =  1;
    atpos[ 9][0] =     0.7804900000;
    atpos[ 9][1] =     0.5609800000;
    atpos[ 9][2] =     0.1178900000;
    
    attyp[10] =  1;
    atpos[10][0] =     0.0487800000;
    atpos[10][1] =     0.0975600000;
    atpos[10][2] =     0.1219500000;
    
    attyp[11] =  1;
    atpos[11][0] =     0.3902400000;
    atpos[11][1] =     0.7804900000;
    atpos[11][2] =     0.1422800000;
    
    attyp[12] =  1;
    atpos[12][0] =     0.6585400000;
    atpos[12][1] =     0.3170700000;
    atpos[12][2] =     0.1463400000;
    
    attyp[13] =  1;
    atpos[13][0] =     0.0000000000;
    atpos[13][1] =     0.0000000000;
    atpos[13][2] =     0.1666700000;
    
    attyp[14] =  1;
    atpos[14][0] =     0.2682900000;
    atpos[14][1] =     0.5365900000;
    atpos[14][2] =     0.1707300000;
    
    attyp[15] =  1;
    atpos[15][0] =     0.6097600000;
    atpos[15][1] =     0.2195100000;
    atpos[15][2] =     0.1910600000;
    
    attyp[16] =  1;
    atpos[16][0] =     0.8780500000;
    atpos[16][1] =     0.7561000000;
    atpos[16][2] =     0.1951200000;
    
    attyp[17] =  1;
    atpos[17][0] =     0.2195100000;
    atpos[17][1] =     0.4390200000;
    atpos[17][2] =     0.2154500000;
    
    attyp[18] =  1;
    atpos[18][0] =     0.4878000000;
    atpos[18][1] =     0.9756100000;
    atpos[18][2] =     0.2195100000;
    
    attyp[19] =  1;
    atpos[19][0] =     0.8292700000;
    atpos[19][1] =     0.6585400000;
    atpos[19][2] =     0.2398400000;
    
    attyp[20] =  1;
    atpos[20][0] =     0.0975600000;
    atpos[20][1] =     0.1951200000;
    atpos[20][2] =     0.2439000000;
    
    attyp[21] =  1;
    atpos[21][0] =     0.4390200000;
    atpos[21][1] =     0.8780500000;
    atpos[21][2] =     0.2642300000;
    
    attyp[22] =  1;
    atpos[22][0] =     0.7073200000;
    atpos[22][1] =     0.4146300000;
    atpos[22][2] =     0.2682900000;
    
    attyp[23] =  1;
    atpos[23][0] =     0.0487800000;
    atpos[23][1] =     0.0975600000;
    atpos[23][2] =     0.2886200000;
    
    attyp[24] =  1;
    atpos[24][0] =     0.3170700000;
    atpos[24][1] =     0.6341500000;
    atpos[24][2] =     0.2926800000;
    
    attyp[25] =  1;
    atpos[25][0] =     0.6585400000;
    atpos[25][1] =     0.3170700000;
    atpos[25][2] =     0.3130100000;
    
    attyp[26] =  1;
    atpos[26][0] =     0.9268300000;
    atpos[26][1] =     0.8536600000;
    atpos[26][2] =     0.3170700000;
    
    attyp[27] =  1;
    atpos[27][0] =     0.2682900000;
    atpos[27][1] =     0.5365900000;
    atpos[27][2] =     0.3374000000;
    
    attyp[28] =  1;
    atpos[28][0] =     0.5365900000;
    atpos[28][1] =     0.0731700000;
    atpos[28][2] =     0.3414600000;
    
    attyp[29] =  1;
    atpos[29][0] =     0.8780500000;
    atpos[29][1] =     0.7561000000;
    atpos[29][2] =     0.3617900000;
    
    attyp[30] =  1;
    atpos[30][0] =     0.1463400000;
    atpos[30][1] =     0.2926800000;
    atpos[30][2] =     0.3658500000;
    
    attyp[31] =  1;
    atpos[31][0] =     0.4878000000;
    atpos[31][1] =     0.9756100000;
    atpos[31][2] =     0.3861800000;
    
    attyp[32] =  1;
    atpos[32][0] =     0.7561000000;
    atpos[32][1] =     0.5122000000;
    atpos[32][2] =     0.3902400000;
    
    attyp[33] =  1;
    atpos[33][0] =     0.0975600000;
    atpos[33][1] =     0.1951200000;
    atpos[33][2] =     0.4105700000;
    
    attyp[34] =  1;
    atpos[34][0] =     0.3658500000;
    atpos[34][1] =     0.7317100000;
    atpos[34][2] =     0.4146300000;
    
    attyp[35] =  1;
    atpos[35][0] =     0.7073200000;
    atpos[35][1] =     0.4146300000;
    atpos[35][2] =     0.4349600000;
    
    attyp[36] =  1;
    atpos[36][0] =     0.9756100000;
    atpos[36][1] =     0.9512200000;
    atpos[36][2] =     0.4390200000;
    
    attyp[37] =  1;
    atpos[37][0] =     0.3170700000;
    atpos[37][1] =     0.6341500000;
    atpos[37][2] =     0.4593500000;
    
    attyp[38] =  1;
    atpos[38][0] =     0.5853700000;
    atpos[38][1] =     0.1707300000;
    atpos[38][2] =     0.4634100000;
    
    attyp[39] =  1;
    atpos[39][0] =     0.9268300000;
    atpos[39][1] =     0.8536600000;
    atpos[39][2] =     0.4837400000;
    
    attyp[40] =  1;
    atpos[40][0] =     0.1951200000;
    atpos[40][1] =     0.3902400000;
    atpos[40][2] =     0.4878000000;
    
    attyp[41] =  1;
    atpos[41][0] =     0.5365900000;
    atpos[41][1] =     0.0731700000;
    atpos[41][2] =     0.5081300000;
    
    attyp[42] =  1;
    atpos[42][0] =     0.8048800000;
    atpos[42][1] =     0.6097600000;
    atpos[42][2] =     0.5122000000;
    
    attyp[43] =  1;
    atpos[43][0] =     0.1463400000;
    atpos[43][1] =     0.2926800000;
    atpos[43][2] =     0.5325200000;
    
    attyp[44] =  1;
    atpos[44][0] =     0.4146300000;
    atpos[44][1] =     0.8292700000;
    atpos[44][2] =     0.5365900000;
    
    attyp[45] =  1;
    atpos[45][0] =     0.7561000000;
    atpos[45][1] =     0.5122000000;
    atpos[45][2] =     0.5569100000;
    
    attyp[46] =  1;
    atpos[46][0] =     0.0243900000;
    atpos[46][1] =     0.0487800000;
    atpos[46][2] =     0.5609800000;
    
    attyp[47] =  1;
    atpos[47][0] =     0.3658500000;
    atpos[47][1] =     0.7317100000;
    atpos[47][2] =     0.5813000000;
    
    attyp[48] =  1;
    atpos[48][0] =     0.6341500000;
    atpos[48][1] =     0.2682900000;
    atpos[48][2] =     0.5853700000;
    
    attyp[49] =  1;
    atpos[49][0] =     0.9756100000;
    atpos[49][1] =     0.9512200000;
    atpos[49][2] =     0.6056900000;
    
    attyp[50] =  1;
    atpos[50][0] =     0.2439000000;
    atpos[50][1] =     0.4878000000;
    atpos[50][2] =     0.6097600000;
    
    attyp[51] =  1;
    atpos[51][0] =     0.5853700000;
    atpos[51][1] =     0.1707300000;
    atpos[51][2] =     0.6300800000;
    
    attyp[52] =  1;
    atpos[52][0] =     0.8536600000;
    atpos[52][1] =     0.7073200000;
    atpos[52][2] =     0.6341500000;
    
    attyp[53] =  1;
    atpos[53][0] =     0.1951200000;
    atpos[53][1] =     0.3902400000;
    atpos[53][2] =     0.6544700000;
    
    attyp[54] =  1;
    atpos[54][0] =     0.4634100000;
    atpos[54][1] =     0.9268300000;
    atpos[54][2] =     0.6585400000;
    
    attyp[55] =  1;
    atpos[55][0] =     0.8048800000;
    atpos[55][1] =     0.6097600000;
    atpos[55][2] =     0.6788600000;
    
    attyp[56] =  1;
    atpos[56][0] =     0.0731700000;
    atpos[56][1] =     0.1463400000;
    atpos[56][2] =     0.6829300000;
    
    attyp[57] =  1;
    atpos[57][0] =     0.4146300000;
    atpos[57][1] =     0.8292700000;
    atpos[57][2] =     0.7032500000;
    
    attyp[58] =  1;
    atpos[58][0] =     0.6829300000;
    atpos[58][1] =     0.3658500000;
    atpos[58][2] =     0.7073200000;
    
    attyp[59] =  1;
    atpos[59][0] =     0.0243900000;
    atpos[59][1] =     0.0487800000;
    atpos[59][2] =     0.7276400000;
    
    attyp[60] =  1;
    atpos[60][0] =     0.2926800000;
    atpos[60][1] =     0.5853700000;
    atpos[60][2] =     0.7317100000;
    
    attyp[61] =  1;
    atpos[61][0] =     0.6341500000;
    atpos[61][1] =     0.2682900000;
    atpos[61][2] =     0.7520300000;
    
    attyp[62] =  1;
    atpos[62][0] =     0.9024400000;
    atpos[62][1] =     0.8048800000;
    atpos[62][2] =     0.7561000000;
    
    attyp[63] =  1;
    atpos[63][0] =     0.2439000000;
    atpos[63][1] =     0.4878000000;
    atpos[63][2] =     0.7764200000;
    
    attyp[64] =  1;
    atpos[64][0] =     0.5122000000;
    atpos[64][1] =     0.0243900000;
    atpos[64][2] =     0.7804900000;
    
    attyp[65] =  1;
    atpos[65][0] =     0.8536600000;
    atpos[65][1] =     0.7073200000;
    atpos[65][2] =     0.8008100000;
    
    attyp[66] =  1;
    atpos[66][0] =     0.1219500000;
    atpos[66][1] =     0.2439000000;
    atpos[66][2] =     0.8048800000;
    
    attyp[67] =  1;
    atpos[67][0] =     0.4634100000;
    atpos[67][1] =     0.9268300000;
    atpos[67][2] =     0.8252000000;
    
    attyp[68] =  1;
    atpos[68][0] =     0.7317100000;
    atpos[68][1] =     0.4634100000;
    atpos[68][2] =     0.8292700000;
    
    attyp[69] =  1;
    atpos[69][0] =     0.0731700000;
    atpos[69][1] =     0.1463400000;
    atpos[69][2] =     0.8495900000;
    
    attyp[70] =  1;
    atpos[70][0] =     0.3414600000;
    atpos[70][1] =     0.6829300000;
    atpos[70][2] =     0.8536600000;
    
    attyp[71] =  1;
    atpos[71][0] =     0.6829300000;
    atpos[71][1] =     0.3658500000;
    atpos[71][2] =     0.8739800000;
    
    attyp[72] =  1;
    atpos[72][0] =     0.9512200000;
    atpos[72][1] =     0.9024400000;
    atpos[72][2] =     0.8780500000;
    
    attyp[73] =  1;
    atpos[73][0] =     0.2926800000;
    atpos[73][1] =     0.5853700000;
    atpos[73][2] =     0.8983700000;
    
    attyp[74] =  1;
    atpos[74][0] =     0.5609800000;
    atpos[74][1] =     0.1219500000;
    atpos[74][2] =     0.9024400000;
    
    attyp[75] =  1;
    atpos[75][0] =     0.9024400000;
    atpos[75][1] =     0.8048800000;
    atpos[75][2] =     0.9227600000;
    
    attyp[76] =  1;
    atpos[76][0] =     0.1707300000;
    atpos[76][1] =     0.3414600000;
    atpos[76][2] =     0.9268300000;
    
    attyp[77] =  1;
    atpos[77][0] =     0.5122000000;
    atpos[77][1] =     0.0243900000;
    atpos[77][2] =     0.9471500000;
    
    attyp[78] =  1;
    atpos[78][0] =     0.7804900000;
    atpos[78][1] =     0.5609800000;
    atpos[78][2] =     0.9512200000;
    
    attyp[79] =  1;
    atpos[79][0] =     0.1219500000;
    atpos[79][1] =     0.2439000000;
    atpos[79][2] =     0.9715400000;
    
    attyp[80] =  1;
    atpos[80][0] =     0.3902400000;
    atpos[80][1] =     0.7804900000;
    atpos[80][2] =     0.9756100000;
    
    attyp[81] =  1;
    atpos[81][0] =     0.7317100000;
    atpos[81][1] =     0.4634100000;
    atpos[81][2] =     0.9959300000;
    
    initialized = 1;
    break;

  case 2:
    latvec[0][0] =  c_ang;
    latvec[0][1] = -sqrt(1. - c_ang*c_ang);
    latvec[1][1] = sqrt(ca*ca + 1.);

    latvec[2][2] = 41./sqrt(4./3. + 1./(ca*ca));

    nucell = 82;
    ntype  = 1;
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
  
    attyp[ 0] =  1;
    atpos[ 0][0] =     0.0000000000;
    atpos[ 0][1] =     0.0000000000;
    atpos[ 0][2] =     0.0000000000;
    
    attyp[ 1] =  1;
    atpos[ 1][0] =     0.3414600000;
    atpos[ 1][1] =     0.6829300000;
    atpos[ 1][2] =     0.0203300000;
    
    attyp[ 2] =  1;
    atpos[ 2][0] =     0.6097600000;
    atpos[ 2][1] =     0.2195100000;
    atpos[ 2][2] =     0.0243900000;
    
    attyp[ 3] =  1;
    atpos[ 3][0] =     0.9512200000;
    atpos[ 3][1] =     0.9024400000;
    atpos[ 3][2] =     0.0447200000;
    
    attyp[ 4] =  1;
    atpos[ 4][0] =     0.2195100000;
    atpos[ 4][1] =     0.4390200000;
    atpos[ 4][2] =     0.0487800000;
    
    attyp[ 5] =  1;
    atpos[ 5][0] =     0.5609800000;
    atpos[ 5][1] =     0.1219500000;
    atpos[ 5][2] =     0.0691100000;
    
    attyp[ 6] =  1;
    atpos[ 6][0] =     0.8292700000;
    atpos[ 6][1] =     0.6585400000;
    atpos[ 6][2] =     0.0731700000;
    
    attyp[ 7] =  1;
    atpos[ 7][0] =     0.1707300000;
    atpos[ 7][1] =     0.3414600000;
    atpos[ 7][2] =     0.0935000000;
    
    attyp[ 8] =  1;
    atpos[ 8][0] =     0.4390200000;
    atpos[ 8][1] =     0.8780500000;
    atpos[ 8][2] =     0.0975600000;
    
    attyp[ 9] =  1;
    atpos[ 9][0] =     0.7804900000;
    atpos[ 9][1] =     0.5609800000;
    atpos[ 9][2] =     0.1178900000;
    
    attyp[10] =  1;
    atpos[10][0] =     0.0487800000;
    atpos[10][1] =     0.0975600000;
    atpos[10][2] =     0.1219500000;
    
    attyp[11] =  1;
    atpos[11][0] =     0.3902400000;
    atpos[11][1] =     0.7804900000;
    atpos[11][2] =     0.1422800000;
    
    attyp[12] =  1;
    atpos[12][0] =     0.6585400000;
    atpos[12][1] =     0.3170700000;
    atpos[12][2] =     0.1463400000;
    
    attyp[13] =  1;
    atpos[13][0] =     0.0000000000;
    atpos[13][1] =     0.0000000000;
    atpos[13][2] =     0.1666700000;
    
    attyp[14] =  1;
    atpos[14][0] =     0.2682900000;
    atpos[14][1] =     0.5365900000;
    atpos[14][2] =     0.1707300000;
    
    attyp[15] =  1;
    atpos[15][0] =     0.6097600000;
    atpos[15][1] =     0.2195100000;
    atpos[15][2] =     0.1910600000;
    
    attyp[16] =  1;
    atpos[16][0] =     0.8780500000;
    atpos[16][1] =     0.7561000000;
    atpos[16][2] =     0.1951200000;
    
    attyp[17] =  1;
    atpos[17][0] =     0.2195100000;
    atpos[17][1] =     0.4390200000;
    atpos[17][2] =     0.2154500000;
    
    attyp[18] =  1;
    atpos[18][0] =     0.4878000000;
    atpos[18][1] =     0.9756100000;
    atpos[18][2] =     0.2195100000;
    
    attyp[19] =  1;
    atpos[19][0] =     0.8292700000;
    atpos[19][1] =     0.6585400000;
    atpos[19][2] =     0.2398400000;
    
    attyp[20] =  1;
    atpos[20][0] =     0.0975600000;
    atpos[20][1] =     0.1951200000;
    atpos[20][2] =     0.2439000000;
    
    attyp[21] =  1;
    atpos[21][0] =     0.4390200000;
    atpos[21][1] =     0.8780500000;
    atpos[21][2] =     0.2642300000;
    
    attyp[22] =  1;
    atpos[22][0] =     0.7073200000;
    atpos[22][1] =     0.4146300000;
    atpos[22][2] =     0.2682900000;
    
    attyp[23] =  1;
    atpos[23][0] =     0.0487800000;
    atpos[23][1] =     0.0975600000;
    atpos[23][2] =     0.2886200000;
    
    attyp[24] =  1;
    atpos[24][0] =     0.3170700000;
    atpos[24][1] =     0.6341500000;
    atpos[24][2] =     0.2926800000;
    
    attyp[25] =  1;
    atpos[25][0] =     0.6585400000;
    atpos[25][1] =     0.3170700000;
    atpos[25][2] =     0.3130100000;
    
    attyp[26] =  1;
    atpos[26][0] =     0.9268300000;
    atpos[26][1] =     0.8536600000;
    atpos[26][2] =     0.3170700000;
    
    attyp[27] =  1;
    atpos[27][0] =     0.2682900000;
    atpos[27][1] =     0.5365900000;
    atpos[27][2] =     0.3374000000;
    
    attyp[28] =  1;
    atpos[28][0] =     0.5365900000;
    atpos[28][1] =     0.0731700000;
    atpos[28][2] =     0.3414600000;
    
    attyp[29] =  1;
    atpos[29][0] =     0.8780500000;
    atpos[29][1] =     0.7561000000;
    atpos[29][2] =     0.3617900000;
    
    attyp[30] =  1;
    atpos[30][0] =     0.1463400000;
    atpos[30][1] =     0.2926800000;
    atpos[30][2] =     0.3658500000;
    
    attyp[31] =  1;
    atpos[31][0] =     0.4878000000;
    atpos[31][1] =     0.9756100000;
    atpos[31][2] =     0.3861800000;
    
    attyp[32] =  1;
    atpos[32][0] =     0.7561000000;
    atpos[32][1] =     0.5122000000;
    atpos[32][2] =     0.3902400000;
    
    attyp[33] =  1;
    atpos[33][0] =     0.0975600000;
    atpos[33][1] =     0.1951200000;
    atpos[33][2] =     0.4105700000;
    
    attyp[34] =  1;
    atpos[34][0] =     0.3658500000;
    atpos[34][1] =     0.7317100000;
    atpos[34][2] =     0.4146300000;
    
    attyp[35] =  1;
    atpos[35][0] =     0.7073200000;
    atpos[35][1] =     0.4146300000;
    atpos[35][2] =     0.4349600000;
    
    attyp[36] =  1;
    atpos[36][0] =     0.9756100000;
    atpos[36][1] =     0.9512200000;
    atpos[36][2] =     0.4390200000;
    
    attyp[37] =  1;
    atpos[37][0] =     0.3170700000;
    atpos[37][1] =     0.6341500000;
    atpos[37][2] =     0.4593500000;
    
    attyp[38] =  1;
    atpos[38][0] =     0.5853700000;
    atpos[38][1] =     0.1707300000;
    atpos[38][2] =     0.4634100000;
    
    attyp[39] =  1;
    atpos[39][0] =     0.9268300000;
    atpos[39][1] =     0.8536600000;
    atpos[39][2] =     0.4837400000;
    
    attyp[40] =  1;
    atpos[40][0] =     0.1951200000;
    atpos[40][1] =     0.3902400000;
    atpos[40][2] =     0.4878000000;
    
    attyp[41] =  1;
    atpos[41][0] =     0.5365900000;
    atpos[41][1] =     0.0731700000;
    atpos[41][2] =     0.5081300000;
    
    attyp[42] =  1;
    atpos[42][0] =     0.8048800000;
    atpos[42][1] =     0.6097600000;
    atpos[42][2] =     0.5122000000;
    
    attyp[43] =  1;
    atpos[43][0] =     0.1463400000;
    atpos[43][1] =     0.2926800000;
    atpos[43][2] =     0.5325200000;
    
    attyp[44] =  1;
    atpos[44][0] =     0.4146300000;
    atpos[44][1] =     0.8292700000;
    atpos[44][2] =     0.5365900000;
    
    attyp[45] =  1;
    atpos[45][0] =     0.7561000000;
    atpos[45][1] =     0.5122000000;
    atpos[45][2] =     0.5569100000;
    
    attyp[46] =  1;
    atpos[46][0] =     0.0243900000;
    atpos[46][1] =     0.0487800000;
    atpos[46][2] =     0.5609800000;
    
    attyp[47] =  1;
    atpos[47][0] =     0.3658500000;
    atpos[47][1] =     0.7317100000;
    atpos[47][2] =     0.5813000000;
    
    attyp[48] =  1;
    atpos[48][0] =     0.6341500000;
    atpos[48][1] =     0.2682900000;
    atpos[48][2] =     0.5853700000;
    
    attyp[49] =  1;
    atpos[49][0] =     0.9756100000;
    atpos[49][1] =     0.9512200000;
    atpos[49][2] =     0.6056900000;
    
    attyp[50] =  1;
    atpos[50][0] =     0.2439000000;
    atpos[50][1] =     0.4878000000;
    atpos[50][2] =     0.6097600000;
    
    attyp[51] =  1;
    atpos[51][0] =     0.5853700000;
    atpos[51][1] =     0.1707300000;
    atpos[51][2] =     0.6300800000;
    
    attyp[52] =  1;
    atpos[52][0] =     0.8536600000;
    atpos[52][1] =     0.7073200000;
    atpos[52][2] =     0.6341500000;
    
    attyp[53] =  1;
    atpos[53][0] =     0.1951200000;
    atpos[53][1] =     0.3902400000;
    atpos[53][2] =     0.6544700000;
    
    attyp[54] =  1;
    atpos[54][0] =     0.4634100000;
    atpos[54][1] =     0.9268300000;
    atpos[54][2] =     0.6585400000;
    
    attyp[55] =  1;
    atpos[55][0] =     0.8048800000;
    atpos[55][1] =     0.6097600000;
    atpos[55][2] =     0.6788600000;
    
    attyp[56] =  1;
    atpos[56][0] =     0.0731700000;
    atpos[56][1] =     0.1463400000;
    atpos[56][2] =     0.6829300000;
    
    attyp[57] =  1;
    atpos[57][0] =     0.4146300000;
    atpos[57][1] =     0.8292700000;
    atpos[57][2] =     0.7032500000;
    
    attyp[58] =  1;
    atpos[58][0] =     0.6829300000;
    atpos[58][1] =     0.3658500000;
    atpos[58][2] =     0.7073200000;
    
    attyp[59] =  1;
    atpos[59][0] =     0.0243900000;
    atpos[59][1] =     0.0487800000;
    atpos[59][2] =     0.7276400000;
    
    attyp[60] =  1;
    atpos[60][0] =     0.2926800000;
    atpos[60][1] =     0.5853700000;
    atpos[60][2] =     0.7317100000;
    
    attyp[61] =  1;
    atpos[61][0] =     0.6341500000;
    atpos[61][1] =     0.2682900000;
    atpos[61][2] =     0.7520300000;
    
    attyp[62] =  1;
    atpos[62][0] =     0.9024400000;
    atpos[62][1] =     0.8048800000;
    atpos[62][2] =     0.7561000000;
    
    attyp[63] =  1;
    atpos[63][0] =     0.2439000000;
    atpos[63][1] =     0.4878000000;
    atpos[63][2] =     0.7764200000;
    
    attyp[64] =  1;
    atpos[64][0] =     0.5122000000;
    atpos[64][1] =     0.0243900000;
    atpos[64][2] =     0.7804900000;
    
    attyp[65] =  1;
    atpos[65][0] =     0.8536600000;
    atpos[65][1] =     0.7073200000;
    atpos[65][2] =     0.8008100000;
    
    attyp[66] =  1;
    atpos[66][0] =     0.1219500000;
    atpos[66][1] =     0.2439000000;
    atpos[66][2] =     0.8048800000;
    
    attyp[67] =  1;
    atpos[67][0] =     0.4634100000;
    atpos[67][1] =     0.9268300000;
    atpos[67][2] =     0.8252000000;
    
    attyp[68] =  1;
    atpos[68][0] =     0.7317100000;
    atpos[68][1] =     0.4634100000;
    atpos[68][2] =     0.8292700000;
    
    attyp[69] =  1;
    atpos[69][0] =     0.0731700000;
    atpos[69][1] =     0.1463400000;
    atpos[69][2] =     0.8495900000;
    
    attyp[70] =  1;
    atpos[70][0] =     0.3414600000;
    atpos[70][1] =     0.6829300000;
    atpos[70][2] =     0.8536600000;
    
    attyp[71] =  1;
    atpos[71][0] =     0.6829300000;
    atpos[71][1] =     0.3658500000;
    atpos[71][2] =     0.8739800000;
    
    attyp[72] =  1;
    atpos[72][0] =     0.9512200000;
    atpos[72][1] =     0.9024400000;
    atpos[72][2] =     0.8780500000;
    
    attyp[73] =  1;
    atpos[73][0] =     0.2926800000;
    atpos[73][1] =     0.5853700000;
    atpos[73][2] =     0.8983700000;
    
    attyp[74] =  1;
    atpos[74][0] =     0.5609800000;
    atpos[74][1] =     0.1219500000;
    atpos[74][2] =     0.9024400000;
    
    attyp[75] =  1;
    atpos[75][0] =     0.9024400000;
    atpos[75][1] =     0.8048800000;
    atpos[75][2] =     0.9227600000;
    
    attyp[76] =  1;
    atpos[76][0] =     0.1707300000;
    atpos[76][1] =     0.3414600000;
    atpos[76][2] =     0.9268300000;
    
    attyp[77] =  1;
    atpos[77][0] =     0.5122000000;
    atpos[77][1] =     0.0243900000;
    atpos[77][2] =     0.9471500000;
    
    attyp[78] =  1;
    atpos[78][0] =     0.7804900000;
    atpos[78][1] =     0.5609800000;
    atpos[78][2] =     0.9512200000;
    
    attyp[79] =  1;
    atpos[79][0] =     0.1219500000;
    atpos[79][1] =     0.2439000000;
    atpos[79][2] =     0.9715400000;
    
    attyp[80] =  1;
    atpos[80][0] =     0.3902400000;
    atpos[80][1] =     0.7804900000;
    atpos[80][2] =     0.9756100000;
    
    attyp[81] =  1;
    atpos[81][0] =     0.7317100000;
    atpos[81][1] =     0.4634100000;
    atpos[81][2] =     0.9959300000;
    
    initialized = 1;
    break;

  case 3:
    latvec[0][0] = 1.;
    latvec[1][1] = sqrt(3.+4.*ca*ca);
    latvec[2][2] = 41./sqrt(4./3. + 1./(ca*ca));
    
    nucell =  164;
    ntype  =  1;
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[  0] =  1;
    atpos[  0][0] =     0.0000000000;
    atpos[  0][1] =     0.0000000000;
    atpos[  0][2] =     0.0000000000;
    
    attyp[  1] =  1;
    atpos[  1][0] =     0.5000000000;
    atpos[  1][1] =     0.5000000000;
    atpos[  1][2] =     0.0000000000;
    
    attyp[  2] =  1;
    atpos[  2][0] =     0.0000000000;
    atpos[  2][1] =     0.3414600000;
    atpos[  2][2] =     0.0203300000;
    
    attyp[  3] =  1;
    atpos[  3][0] =     0.5000000000;
    atpos[  3][1] =     0.8414600000;
    atpos[  3][2] =     0.0203300000;
    
    attyp[  4] =  1;
    atpos[  4][0] =     0.5000000000;
    atpos[  4][1] =     0.1097600000;
    atpos[  4][2] =     0.0243900000;
    
    attyp[  5] =  1;
    atpos[  5][0] =     0.0000000000;
    atpos[  5][1] =     0.6097600000;
    atpos[  5][2] =     0.0243900000;
    
    attyp[  6] =  1;
    atpos[  6][0] =     0.5000000000;
    atpos[  6][1] =     0.4512200000;
    atpos[  6][2] =     0.0447200000;
    
    attyp[  7] =  1;
    atpos[  7][0] =     0.0000000000;
    atpos[  7][1] =     0.9512200000;
    atpos[  7][2] =     0.0447200000;
    
    attyp[  8] =  1;
    atpos[  8][0] =     0.0000000000;
    atpos[  8][1] =     0.2195100000;
    atpos[  8][2] =     0.0487800000;
    
    attyp[  9] =  1;
    atpos[  9][0] =     0.5000000000;
    atpos[  9][1] =     0.7195100000;
    atpos[  9][2] =     0.0487800000;
    
    attyp[ 10] =  1;
    atpos[ 10][0] =     0.0000000000;
    atpos[ 10][1] =     0.5609800000;
    atpos[ 10][2] =     0.0691100000;
    
    attyp[ 11] =  1;
    atpos[ 11][0] =     0.5000000000;
    atpos[ 11][1] =     0.0609800000;
    atpos[ 11][2] =     0.0691100000;
    
    attyp[ 12] =  1;
    atpos[ 12][0] =     0.5000000000;
    atpos[ 12][1] =     0.3292700000;
    atpos[ 12][2] =     0.0731700000;
    
    attyp[ 13] =  1;
    atpos[ 13][0] =     0.0000000000;
    atpos[ 13][1] =     0.8292700000;
    atpos[ 13][2] =     0.0731700000;
    
    attyp[ 14] =  1;
    atpos[ 14][0] =     0.0000000000;
    atpos[ 14][1] =     0.1707300000;
    atpos[ 14][2] =     0.0935000000;
    
    attyp[ 15] =  1;
    atpos[ 15][0] =     0.5000000000;
    atpos[ 15][1] =     0.6707300000;
    atpos[ 15][2] =     0.0935000000;
    
    attyp[ 16] =  1;
    atpos[ 16][0] =     0.0000000000;
    atpos[ 16][1] =     0.4390200000;
    atpos[ 16][2] =     0.0975600000;
    
    attyp[ 17] =  1;
    atpos[ 17][0] =     0.5000000000;
    atpos[ 17][1] =     0.9390200000;
    atpos[ 17][2] =     0.0975600000;
    
    attyp[ 18] =  1;
    atpos[ 18][0] =     0.5000000000;
    atpos[ 18][1] =     0.2804900000;
    atpos[ 18][2] =     0.1178900000;
    
    attyp[ 19] =  1;
    atpos[ 19][0] =     0.0000000000;
    atpos[ 19][1] =     0.7804900000;
    atpos[ 19][2] =     0.1178900000;
    
    attyp[ 20] =  1;
    atpos[ 20][0] =     0.5000000000;
    atpos[ 20][1] =     0.5487800000;
    atpos[ 20][2] =     0.1219500000;
    
    attyp[ 21] =  1;
    atpos[ 21][0] =     0.0000000000;
    atpos[ 21][1] =     0.0487800000;
    atpos[ 21][2] =     0.1219500000;
    
    attyp[ 22] =  1;
    atpos[ 22][0] =     0.0000000000;
    atpos[ 22][1] =     0.3902400000;
    atpos[ 22][2] =     0.1422800000;
    
    attyp[ 23] =  1;
    atpos[ 23][0] =     0.5000000000;
    atpos[ 23][1] =     0.8902400000;
    atpos[ 23][2] =     0.1422800000;
    
    attyp[ 24] =  1;
    atpos[ 24][0] =     0.5000000000;
    atpos[ 24][1] =     0.1585400000;
    atpos[ 24][2] =     0.1463400000;
    
    attyp[ 25] =  1;
    atpos[ 25][0] =     0.0000000000;
    atpos[ 25][1] =     0.6585400000;
    atpos[ 25][2] =     0.1463400000;
    
    attyp[ 26] =  1;
    atpos[ 26][0] =     0.5000000000;
    atpos[ 26][1] =     0.5000000000;
    atpos[ 26][2] =     0.1666700000;
    
    attyp[ 27] =  1;
    atpos[ 27][0] =     0.0000000000;
    atpos[ 27][1] =     0.0000000000;
    atpos[ 27][2] =     0.1666700000;
    
    attyp[ 28] =  1;
    atpos[ 28][0] =     0.0000000000;
    atpos[ 28][1] =     0.2682900000;
    atpos[ 28][2] =     0.1707300000;
    
    attyp[ 29] =  1;
    atpos[ 29][0] =     0.5000000000;
    atpos[ 29][1] =     0.7682900000;
    atpos[ 29][2] =     0.1707300000;
    
    attyp[ 30] =  1;
    atpos[ 30][0] =     0.5000000000;
    atpos[ 30][1] =     0.1097600000;
    atpos[ 30][2] =     0.1910600000;
    
    attyp[ 31] =  1;
    atpos[ 31][0] =     0.0000000000;
    atpos[ 31][1] =     0.6097600000;
    atpos[ 31][2] =     0.1910600000;
    
    attyp[ 32] =  1;
    atpos[ 32][0] =     0.5000000000;
    atpos[ 32][1] =     0.3780500000;
    atpos[ 32][2] =     0.1951200000;
    
    attyp[ 33] =  1;
    atpos[ 33][0] =     0.0000000000;
    atpos[ 33][1] =     0.8780500000;
    atpos[ 33][2] =     0.1951200000;
    
    attyp[ 34] =  1;
    atpos[ 34][0] =     0.0000000000;
    atpos[ 34][1] =     0.2195100000;
    atpos[ 34][2] =     0.2154500000;
    
    attyp[ 35] =  1;
    atpos[ 35][0] =     0.5000000000;
    atpos[ 35][1] =     0.7195100000;
    atpos[ 35][2] =     0.2154500000;
    
    attyp[ 36] =  1;
    atpos[ 36][0] =     0.0000000000;
    atpos[ 36][1] =     0.4878000000;
    atpos[ 36][2] =     0.2195100000;
    
    attyp[ 37] =  1;
    atpos[ 37][0] =     0.5000000000;
    atpos[ 37][1] =     0.9878000000;
    atpos[ 37][2] =     0.2195100000;
    
    attyp[ 38] =  1;
    atpos[ 38][0] =     0.5000000000;
    atpos[ 38][1] =     0.3292700000;
    atpos[ 38][2] =     0.2398400000;
    
    attyp[ 39] =  1;
    atpos[ 39][0] =     0.0000000000;
    atpos[ 39][1] =     0.8292700000;
    atpos[ 39][2] =     0.2398400000;
    
    attyp[ 40] =  1;
    atpos[ 40][0] =     0.5000000000;
    atpos[ 40][1] =     0.5975600000;
    atpos[ 40][2] =     0.2439000000;
    
    attyp[ 41] =  1;
    atpos[ 41][0] =     0.0000000000;
    atpos[ 41][1] =     0.0975600000;
    atpos[ 41][2] =     0.2439000000;
    
    attyp[ 42] =  1;
    atpos[ 42][0] =     0.0000000000;
    atpos[ 42][1] =     0.4390200000;
    atpos[ 42][2] =     0.2642300000;
    
    attyp[ 43] =  1;
    atpos[ 43][0] =     0.5000000000;
    atpos[ 43][1] =     0.9390200000;
    atpos[ 43][2] =     0.2642300000;
    
    attyp[ 44] =  1;
    atpos[ 44][0] =     0.5000000000;
    atpos[ 44][1] =     0.2073200000;
    atpos[ 44][2] =     0.2682900000;
    
    attyp[ 45] =  1;
    atpos[ 45][0] =     0.0000000000;
    atpos[ 45][1] =     0.7073200000;
    atpos[ 45][2] =     0.2682900000;
    
    attyp[ 46] =  1;
    atpos[ 46][0] =     0.5000000000;
    atpos[ 46][1] =     0.5487800000;
    atpos[ 46][2] =     0.2886200000;
    
    attyp[ 47] =  1;
    atpos[ 47][0] =     0.0000000000;
    atpos[ 47][1] =     0.0487800000;
    atpos[ 47][2] =     0.2886200000;
    
    attyp[ 48] =  1;
    atpos[ 48][0] =     0.0000000000;
    atpos[ 48][1] =     0.3170700000;
    atpos[ 48][2] =     0.2926800000;
    
    attyp[ 49] =  1;
    atpos[ 49][0] =     0.5000000000;
    atpos[ 49][1] =     0.8170700000;
    atpos[ 49][2] =     0.2926800000;
    
    attyp[ 50] =  1;
    atpos[ 50][0] =     0.5000000000;
    atpos[ 50][1] =     0.1585400000;
    atpos[ 50][2] =     0.3130100000;
    
    attyp[ 51] =  1;
    atpos[ 51][0] =     0.0000000000;
    atpos[ 51][1] =     0.6585400000;
    atpos[ 51][2] =     0.3130100000;
    
    attyp[ 52] =  1;
    atpos[ 52][0] =     0.5000000000;
    atpos[ 52][1] =     0.4268300000;
    atpos[ 52][2] =     0.3170700000;
    
    attyp[ 53] =  1;
    atpos[ 53][0] =     0.0000000000;
    atpos[ 53][1] =     0.9268300000;
    atpos[ 53][2] =     0.3170700000;
    
    attyp[ 54] =  1;
    atpos[ 54][0] =     0.0000000000;
    atpos[ 54][1] =     0.2682900000;
    atpos[ 54][2] =     0.3374000000;
    
    attyp[ 55] =  1;
    atpos[ 55][0] =     0.5000000000;
    atpos[ 55][1] =     0.7682900000;
    atpos[ 55][2] =     0.3374000000;
    
    attyp[ 56] =  1;
    atpos[ 56][0] =     0.0000000000;
    atpos[ 56][1] =     0.5365900000;
    atpos[ 56][2] =     0.3414600000;
    
    attyp[ 57] =  1;
    atpos[ 57][0] =     0.5000000000;
    atpos[ 57][1] =     0.0365900000;
    atpos[ 57][2] =     0.3414600000;
    
    attyp[ 58] =  1;
    atpos[ 58][0] =     0.5000000000;
    atpos[ 58][1] =     0.3780500000;
    atpos[ 58][2] =     0.3617900000;
    
    attyp[ 59] =  1;
    atpos[ 59][0] =     0.0000000000;
    atpos[ 59][1] =     0.8780500000;
    atpos[ 59][2] =     0.3617900000;
    
    attyp[ 60] =  1;
    atpos[ 60][0] =     0.0000000000;
    atpos[ 60][1] =     0.1463400000;
    atpos[ 60][2] =     0.3658500000;
    
    attyp[ 61] =  1;
    atpos[ 61][0] =     0.5000000000;
    atpos[ 61][1] =     0.6463400000;
    atpos[ 61][2] =     0.3658500000;
    
    attyp[ 62] =  1;
    atpos[ 62][0] =     0.0000000000;
    atpos[ 62][1] =     0.4878000000;
    atpos[ 62][2] =     0.3861800000;
    
    attyp[ 63] =  1;
    atpos[ 63][0] =     0.5000000000;
    atpos[ 63][1] =     0.9878000000;
    atpos[ 63][2] =     0.3861800000;
    
    attyp[ 64] =  1;
    atpos[ 64][0] =     0.5000000000;
    atpos[ 64][1] =     0.2561000000;
    atpos[ 64][2] =     0.3902400000;
    
    attyp[ 65] =  1;
    atpos[ 65][0] =     0.0000000000;
    atpos[ 65][1] =     0.7561000000;
    atpos[ 65][2] =     0.3902400000;
    
    attyp[ 66] =  1;
    atpos[ 66][0] =     0.5000000000;
    atpos[ 66][1] =     0.5975600000;
    atpos[ 66][2] =     0.4105700000;
    
    attyp[ 67] =  1;
    atpos[ 67][0] =     0.0000000000;
    atpos[ 67][1] =     0.0975600000;
    atpos[ 67][2] =     0.4105700000;
    
    attyp[ 68] =  1;
    atpos[ 68][0] =     0.0000000000;
    atpos[ 68][1] =     0.3658500000;
    atpos[ 68][2] =     0.4146300000;
    
    attyp[ 69] =  1;
    atpos[ 69][0] =     0.5000000000;
    atpos[ 69][1] =     0.8658500000;
    atpos[ 69][2] =     0.4146300000;
    
    attyp[ 70] =  1;
    atpos[ 70][0] =     0.5000000000;
    atpos[ 70][1] =     0.2073200000;
    atpos[ 70][2] =     0.4349600000;
    
    attyp[ 71] =  1;
    atpos[ 71][0] =     0.0000000000;
    atpos[ 71][1] =     0.7073200000;
    atpos[ 71][2] =     0.4349600000;
    
    attyp[ 72] =  1;
    atpos[ 72][0] =     0.5000000000;
    atpos[ 72][1] =     0.4756100000;
    atpos[ 72][2] =     0.4390200000;
    
    attyp[ 73] =  1;
    atpos[ 73][0] =     0.0000000000;
    atpos[ 73][1] =     0.9756100000;
    atpos[ 73][2] =     0.4390200000;
    
    attyp[ 74] =  1;
    atpos[ 74][0] =     0.0000000000;
    atpos[ 74][1] =     0.3170700000;
    atpos[ 74][2] =     0.4593500000;
    
    attyp[ 75] =  1;
    atpos[ 75][0] =     0.5000000000;
    atpos[ 75][1] =     0.8170700000;
    atpos[ 75][2] =     0.4593500000;
    
    attyp[ 76] =  1;
    atpos[ 76][0] =     0.0000000000;
    atpos[ 76][1] =     0.5853700000;
    atpos[ 76][2] =     0.4634100000;
    
    attyp[ 77] =  1;
    atpos[ 77][0] =     0.5000000000;
    atpos[ 77][1] =     0.0853700000;
    atpos[ 77][2] =     0.4634100000;
    
    attyp[ 78] =  1;
    atpos[ 78][0] =     0.5000000000;
    atpos[ 78][1] =     0.4268300000;
    atpos[ 78][2] =     0.4837400000;
    
    attyp[ 79] =  1;
    atpos[ 79][0] =     0.0000000000;
    atpos[ 79][1] =     0.9268300000;
    atpos[ 79][2] =     0.4837400000;
    
    attyp[ 80] =  1;
    atpos[ 80][0] =     0.0000000000;
    atpos[ 80][1] =     0.1951200000;
    atpos[ 80][2] =     0.4878000000;
    
    attyp[ 81] =  1;
    atpos[ 81][0] =     0.5000000000;
    atpos[ 81][1] =     0.6951200000;
    atpos[ 81][2] =     0.4878000000;
    
    attyp[ 82] =  1;
    atpos[ 82][0] =     0.0000000000;
    atpos[ 82][1] =     0.5365900000;
    atpos[ 82][2] =     0.5081300000;
    
    attyp[ 83] =  1;
    atpos[ 83][0] =     0.5000000000;
    atpos[ 83][1] =     0.0365900000;
    atpos[ 83][2] =     0.5081300000;
    
    attyp[ 84] =  1;
    atpos[ 84][0] =     0.5000000000;
    atpos[ 84][1] =     0.3048800000;
    atpos[ 84][2] =     0.5122000000;
    
    attyp[ 85] =  1;
    atpos[ 85][0] =     0.0000000000;
    atpos[ 85][1] =     0.8048800000;
    atpos[ 85][2] =     0.5122000000;
    
    attyp[ 86] =  1;
    atpos[ 86][0] =     0.0000000000;
    atpos[ 86][1] =     0.1463400000;
    atpos[ 86][2] =     0.5325200000;
    
    attyp[ 87] =  1;
    atpos[ 87][0] =     0.5000000000;
    atpos[ 87][1] =     0.6463400000;
    atpos[ 87][2] =     0.5325200000;
    
    attyp[ 88] =  1;
    atpos[ 88][0] =     0.0000000000;
    atpos[ 88][1] =     0.4146300000;
    atpos[ 88][2] =     0.5365900000;
    
    attyp[ 89] =  1;
    atpos[ 89][0] =     0.5000000000;
    atpos[ 89][1] =     0.9146300000;
    atpos[ 89][2] =     0.5365900000;
    
    attyp[ 90] =  1;
    atpos[ 90][0] =     0.5000000000;
    atpos[ 90][1] =     0.2561000000;
    atpos[ 90][2] =     0.5569100000;
    
    attyp[ 91] =  1;
    atpos[ 91][0] =     0.0000000000;
    atpos[ 91][1] =     0.7561000000;
    atpos[ 91][2] =     0.5569100000;
    
    attyp[ 92] =  1;
    atpos[ 92][0] =     0.5000000000;
    atpos[ 92][1] =     0.5243900000;
    atpos[ 92][2] =     0.5609800000;
    
    attyp[ 93] =  1;
    atpos[ 93][0] =     0.0000000000;
    atpos[ 93][1] =     0.0243900000;
    atpos[ 93][2] =     0.5609800000;
    
    attyp[ 94] =  1;
    atpos[ 94][0] =     0.0000000000;
    atpos[ 94][1] =     0.3658500000;
    atpos[ 94][2] =     0.5813000000;
    
    attyp[ 95] =  1;
    atpos[ 95][0] =     0.5000000000;
    atpos[ 95][1] =     0.8658500000;
    atpos[ 95][2] =     0.5813000000;
    
    attyp[ 96] =  1;
    atpos[ 96][0] =     0.5000000000;
    atpos[ 96][1] =     0.1341500000;
    atpos[ 96][2] =     0.5853700000;
    
    attyp[ 97] =  1;
    atpos[ 97][0] =     0.0000000000;
    atpos[ 97][1] =     0.6341500000;
    atpos[ 97][2] =     0.5853700000;
    
    attyp[ 98] =  1;
    atpos[ 98][0] =     0.5000000000;
    atpos[ 98][1] =     0.4756100000;
    atpos[ 98][2] =     0.6056900000;
    
    attyp[ 99] =  1;
    atpos[ 99][0] =     0.0000000000;
    atpos[ 99][1] =     0.9756100000;
    atpos[ 99][2] =     0.6056900000;
    
    attyp[100] =  1;
    atpos[100][0] =     0.0000000000;
    atpos[100][1] =     0.2439000000;
    atpos[100][2] =     0.6097600000;
    
    attyp[101] =  1;
    atpos[101][0] =     0.5000000000;
    atpos[101][1] =     0.7439000000;
    atpos[101][2] =     0.6097600000;
    
    attyp[102] =  1;
    atpos[102][0] =     0.0000000000;
    atpos[102][1] =     0.5853700000;
    atpos[102][2] =     0.6300800000;
    
    attyp[103] =  1;
    atpos[103][0] =     0.5000000000;
    atpos[103][1] =     0.0853700000;
    atpos[103][2] =     0.6300800000;
    
    attyp[104] =  1;
    atpos[104][0] =     0.5000000000;
    atpos[104][1] =     0.3536600000;
    atpos[104][2] =     0.6341500000;
    
    attyp[105] =  1;
    atpos[105][0] =     0.0000000000;
    atpos[105][1] =     0.8536600000;
    atpos[105][2] =     0.6341500000;
    
    attyp[106] =  1;
    atpos[106][0] =     0.0000000000;
    atpos[106][1] =     0.1951200000;
    atpos[106][2] =     0.6544700000;
    
    attyp[107] =  1;
    atpos[107][0] =     0.5000000000;
    atpos[107][1] =     0.6951200000;
    atpos[107][2] =     0.6544700000;
    
    attyp[108] =  1;
    atpos[108][0] =     0.0000000000;
    atpos[108][1] =     0.4634100000;
    atpos[108][2] =     0.6585400000;
    
    attyp[109] =  1;
    atpos[109][0] =     0.5000000000;
    atpos[109][1] =     0.9634100000;
    atpos[109][2] =     0.6585400000;
    
    attyp[110] =  1;
    atpos[110][0] =     0.5000000000;
    atpos[110][1] =     0.3048800000;
    atpos[110][2] =     0.6788600000;
    
    attyp[111] =  1;
    atpos[111][0] =     0.0000000000;
    atpos[111][1] =     0.8048800000;
    atpos[111][2] =     0.6788600000;
    
    attyp[112] =  1;
    atpos[112][0] =     0.5000000000;
    atpos[112][1] =     0.5731700000;
    atpos[112][2] =     0.6829300000;
    
    attyp[113] =  1;
    atpos[113][0] =     0.0000000000;
    atpos[113][1] =     0.0731700000;
    atpos[113][2] =     0.6829300000;
    
    attyp[114] =  1;
    atpos[114][0] =     0.0000000000;
    atpos[114][1] =     0.4146300000;
    atpos[114][2] =     0.7032500000;
    
    attyp[115] =  1;
    atpos[115][0] =     0.5000000000;
    atpos[115][1] =     0.9146300000;
    atpos[115][2] =     0.7032500000;
    
    attyp[116] =  1;
    atpos[116][0] =     0.5000000000;
    atpos[116][1] =     0.1829300000;
    atpos[116][2] =     0.7073200000;
    
    attyp[117] =  1;
    atpos[117][0] =     0.0000000000;
    atpos[117][1] =     0.6829300000;
    atpos[117][2] =     0.7073200000;
    
    attyp[118] =  1;
    atpos[118][0] =     0.5000000000;
    atpos[118][1] =     0.5243900000;
    atpos[118][2] =     0.7276400000;
    
    attyp[119] =  1;
    atpos[119][0] =     0.0000000000;
    atpos[119][1] =     0.0243900000;
    atpos[119][2] =     0.7276400000;
    
    attyp[120] =  1;
    atpos[120][0] =     0.0000000000;
    atpos[120][1] =     0.2926800000;
    atpos[120][2] =     0.7317100000;
    
    attyp[121] =  1;
    atpos[121][0] =     0.5000000000;
    atpos[121][1] =     0.7926800000;
    atpos[121][2] =     0.7317100000;
    
    attyp[122] =  1;
    atpos[122][0] =     0.5000000000;
    atpos[122][1] =     0.1341500000;
    atpos[122][2] =     0.7520300000;
    
    attyp[123] =  1;
    atpos[123][0] =     0.0000000000;
    atpos[123][1] =     0.6341500000;
    atpos[123][2] =     0.7520300000;
    
    attyp[124] =  1;
    atpos[124][0] =     0.5000000000;
    atpos[124][1] =     0.4024400000;
    atpos[124][2] =     0.7561000000;
    
    attyp[125] =  1;
    atpos[125][0] =     0.0000000000;
    atpos[125][1] =     0.9024400000;
    atpos[125][2] =     0.7561000000;
    
    attyp[126] =  1;
    atpos[126][0] =     0.0000000000;
    atpos[126][1] =     0.2439000000;
    atpos[126][2] =     0.7764200000;
    
    attyp[127] =  1;
    atpos[127][0] =     0.5000000000;
    atpos[127][1] =     0.7439000000;
    atpos[127][2] =     0.7764200000;
    
    attyp[128] =  1;
    atpos[128][0] =     0.0000000000;
    atpos[128][1] =     0.5122000000;
    atpos[128][2] =     0.7804900000;
    
    attyp[129] =  1;
    atpos[129][0] =     0.5000000000;
    atpos[129][1] =     0.0122000000;
    atpos[129][2] =     0.7804900000;
    
    attyp[130] =  1;
    atpos[130][0] =     0.5000000000;
    atpos[130][1] =     0.3536600000;
    atpos[130][2] =     0.8008100000;
    
    attyp[131] =  1;
    atpos[131][0] =     0.0000000000;
    atpos[131][1] =     0.8536600000;
    atpos[131][2] =     0.8008100000;
    
    attyp[132] =  1;
    atpos[132][0] =     0.0000000000;
    atpos[132][1] =     0.1219500000;
    atpos[132][2] =     0.8048800000;
    
    attyp[133] =  1;
    atpos[133][0] =     0.5000000000;
    atpos[133][1] =     0.6219500000;
    atpos[133][2] =     0.8048800000;
    
    attyp[134] =  1;
    atpos[134][0] =     0.0000000000;
    atpos[134][1] =     0.4634100000;
    atpos[134][2] =     0.8252000000;
    
    attyp[135] =  1;
    atpos[135][0] =     0.5000000000;
    atpos[135][1] =     0.9634100000;
    atpos[135][2] =     0.8252000000;
    
    attyp[136] =  1;
    atpos[136][0] =     0.5000000000;
    atpos[136][1] =     0.2317100000;
    atpos[136][2] =     0.8292700000;
    
    attyp[137] =  1;
    atpos[137][0] =     0.0000000000;
    atpos[137][1] =     0.7317100000;
    atpos[137][2] =     0.8292700000;
    
    attyp[138] =  1;
    atpos[138][0] =     0.5000000000;
    atpos[138][1] =     0.5731700000;
    atpos[138][2] =     0.8495900000;
    
    attyp[139] =  1;
    atpos[139][0] =     0.0000000000;
    atpos[139][1] =     0.0731700000;
    atpos[139][2] =     0.8495900000;
    
    attyp[140] =  1;
    atpos[140][0] =     0.0000000000;
    atpos[140][1] =     0.3414600000;
    atpos[140][2] =     0.8536600000;
    
    attyp[141] =  1;
    atpos[141][0] =     0.5000000000;
    atpos[141][1] =     0.8414600000;
    atpos[141][2] =     0.8536600000;
    
    attyp[142] =  1;
    atpos[142][0] =     0.5000000000;
    atpos[142][1] =     0.1829300000;
    atpos[142][2] =     0.8739800000;
    
    attyp[143] =  1;
    atpos[143][0] =     0.0000000000;
    atpos[143][1] =     0.6829300000;
    atpos[143][2] =     0.8739800000;
    
    attyp[144] =  1;
    atpos[144][0] =     0.5000000000;
    atpos[144][1] =     0.4512200000;
    atpos[144][2] =     0.8780500000;
    
    attyp[145] =  1;
    atpos[145][0] =     0.0000000000;
    atpos[145][1] =     0.9512200000;
    atpos[145][2] =     0.8780500000;
    
    attyp[146] =  1;
    atpos[146][0] =     0.0000000000;
    atpos[146][1] =     0.2926800000;
    atpos[146][2] =     0.8983700000;
    
    attyp[147] =  1;
    atpos[147][0] =     0.5000000000;
    atpos[147][1] =     0.7926800000;
    atpos[147][2] =     0.8983700000;
    
    attyp[148] =  1;
    atpos[148][0] =     0.0000000000;
    atpos[148][1] =     0.5609800000;
    atpos[148][2] =     0.9024400000;
    
    attyp[149] =  1;
    atpos[149][0] =     0.5000000000;
    atpos[149][1] =     0.0609800000;
    atpos[149][2] =     0.9024400000;
    
    attyp[150] =  1;
    atpos[150][0] =     0.5000000000;
    atpos[150][1] =     0.4024400000;
    atpos[150][2] =     0.9227600000;
    
    attyp[151] =  1;
    atpos[151][0] =     0.0000000000;
    atpos[151][1] =     0.9024400000;
    atpos[151][2] =     0.9227600000;
    
    attyp[152] =  1;
    atpos[152][0] =     0.0000000000;
    atpos[152][1] =     0.1707300000;
    atpos[152][2] =     0.9268300000;
    
    attyp[153] =  1;
    atpos[153][0] =     0.5000000000;
    atpos[153][1] =     0.6707300000;
    atpos[153][2] =     0.9268300000;
    
    attyp[154] =  1;
    atpos[154][0] =     0.0000000000;
    atpos[154][1] =     0.5122000000;
    atpos[154][2] =     0.9471500000;
    
    attyp[155] =  1;
    atpos[155][0] =     0.5000000000;
    atpos[155][1] =     0.0122000000;
    atpos[155][2] =     0.9471500000;
    
    attyp[156] =  1;
    atpos[156][0] =     0.5000000000;
    atpos[156][1] =     0.2804900000;
    atpos[156][2] =     0.9512200000;
    
    attyp[157] =  1;
    atpos[157][0] =     0.0000000000;
    atpos[157][1] =     0.7804900000;
    atpos[157][2] =     0.9512200000;
    
    attyp[158] =  1;
    atpos[158][0] =     0.0000000000;
    atpos[158][1] =     0.1219500000;
    atpos[158][2] =     0.9715400000;
    
    attyp[159] =  1;
    atpos[159][0] =     0.5000000000;
    atpos[159][1] =     0.6219500000;
    atpos[159][2] =     0.9715400000;
    
    attyp[160] =  1;
    atpos[160][0] =     0.0000000000;
    atpos[160][1] =     0.3902400000;
    atpos[160][2] =     0.9756100000;
    
    attyp[161] =  1;
    atpos[161][0] =     0.5000000000;
    atpos[161][1] =     0.8902400000;
    atpos[161][2] =     0.9756100000;
    
    attyp[162] =  1;
    atpos[162][0] =     0.5000000000;
    atpos[162][1] =     0.2317100000;
    atpos[162][2] =     0.9959300000;
    
    attyp[163] =  1;
    atpos[163][0] =     0.0000000000;
    atpos[163][1] =     0.7317100000;
    atpos[163][2] =     0.9959300000;
    
    initialized = 1;
  default:
    break;
  }
return;
}
/* -----------------------------------------------------------------------------
 * HCP(112)
 * -------------------------------------------------------------------------- */
void HCP::HCP112()
{
  char str[MAXLINE];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");
  printf("   1. [1 -1 0] along x;\n");
  printf("   2. [1 -1 0] along y;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = inumeric(strtok(str, " \t\n\r\f"));
  printf("Your selection : %d", surftype);
  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");

  for (int i = 0; i < 3; ++i)
  for (int j = 0; j < 3; ++j) latvec[i][j] = 0.;

  memory->create(name,9,"HCP:name");
  strcpy(name, "HCP(112)");
 
  // initialize according to surface type
  switch (surftype){
  case 1:
    nucell = 44;
    ntype  = 1;
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");

    latvec[0][0] = sqrt(3.);
    latvec[1][1] = sqrt(ca*ca + 1.);
    latvec[2][2] =  11./sqrt(1.+1./(ca*ca));
    
    attyp[ 0] =  1;
    atpos[ 0][0] =     0.0000000000;
    atpos[ 0][1] =     0.0000000000;
    atpos[ 0][2] =     0.0000000000;
    
    attyp[ 1] =  1;
    atpos[ 1][0] =     0.8333300000;
    atpos[ 1][1] =     0.5000000000;
    atpos[ 1][2] =     0.0000000000;
    
    attyp[ 2] =  1;
    atpos[ 2][0] =     0.5000000000;
    atpos[ 2][1] =     0.1363600000;
    atpos[ 2][2] =     0.0454500000;
    
    attyp[ 3] =  1;
    atpos[ 3][0] =     0.3333300000;
    atpos[ 3][1] =     0.6363600000;
    atpos[ 3][2] =     0.0454500000;
    
    attyp[ 4] =  1;
    atpos[ 4][0] =     0.0000000000;
    atpos[ 4][1] =     0.2727300000;
    atpos[ 4][2] =     0.0909100000;
    
    attyp[ 5] =  1;
    atpos[ 5][0] =     0.8333300000;
    atpos[ 5][1] =     0.7727300000;
    atpos[ 5][2] =     0.0909100000;
    
    attyp[ 6] =  1;
    atpos[ 6][0] =     0.5000000000;
    atpos[ 6][1] =     0.4090900000;
    atpos[ 6][2] =     0.1363600000;
    
    attyp[ 7] =  1;
    atpos[ 7][0] =     0.3333300000;
    atpos[ 7][1] =     0.9090900000;
    atpos[ 7][2] =     0.1363600000;
    
    attyp[ 8] =  1;
    atpos[ 8][0] =     0.8333300000;
    atpos[ 8][1] =     0.0454500000;
    atpos[ 8][2] =     0.1818200000;
    
    attyp[ 9] =  1;
    atpos[ 9][0] =     0.0000000000;
    atpos[ 9][1] =     0.5454500000;
    atpos[ 9][2] =     0.1818200000;
    
    attyp[10] =  1;
    atpos[10][0] =     0.5000000000;
    atpos[10][1] =     0.6818200000;
    atpos[10][2] =     0.2272700000;
    
    attyp[11] =  1;
    atpos[11][0] =     0.3333300000;
    atpos[11][1] =     0.1818200000;
    atpos[11][2] =     0.2272700000;
    
    attyp[12] =  1;
    atpos[12][0] =     0.0000000000;
    atpos[12][1] =     0.8181800000;
    atpos[12][2] =     0.2727300000;
    
    attyp[13] =  1;
    atpos[13][0] =     0.8333300000;
    atpos[13][1] =     0.3181800000;
    atpos[13][2] =     0.2727300000;
    
    attyp[14] =  1;
    atpos[14][0] =     0.5000000000;
    atpos[14][1] =     0.9545500000;
    atpos[14][2] =     0.3181800000;
    
    attyp[15] =  1;
    atpos[15][0] =     0.3333300000;
    atpos[15][1] =     0.4545500000;
    atpos[15][2] =     0.3181800000;
    
    attyp[16] =  1;
    atpos[16][0] =     0.0000000000;
    atpos[16][1] =     0.0909100000;
    atpos[16][2] =     0.3636400000;
    
    attyp[17] =  1;
    atpos[17][0] =     0.8333300000;
    atpos[17][1] =     0.5909100000;
    atpos[17][2] =     0.3636400000;
    
    attyp[18] =  1;
    atpos[18][0] =     0.5000000000;
    atpos[18][1] =     0.2272700000;
    atpos[18][2] =     0.4090900000;
    
    attyp[19] =  1;
    atpos[19][0] =     0.3333300000;
    atpos[19][1] =     0.7272700000;
    atpos[19][2] =     0.4090900000;
    
    attyp[20] =  1;
    atpos[20][0] =     0.0000000000;
    atpos[20][1] =     0.3636400000;
    atpos[20][2] =     0.4545500000;
    
    attyp[21] =  1;
    atpos[21][0] =     0.8333300000;
    atpos[21][1] =     0.8636400000;
    atpos[21][2] =     0.4545500000;
    
    attyp[22] =  1;
    atpos[22][0] =     0.3333300000;
    atpos[22][1] =     0.0000000000;
    atpos[22][2] =     0.5000000000;
    
    attyp[23] =  1;
    atpos[23][0] =     0.5000000000;
    atpos[23][1] =     0.5000000000;
    atpos[23][2] =     0.5000000000;
    
    attyp[24] =  1;
    atpos[24][0] =     0.0000000000;
    atpos[24][1] =     0.6363600000;
    atpos[24][2] =     0.5454500000;
    
    attyp[25] =  1;
    atpos[25][0] =     0.8333300000;
    atpos[25][1] =     0.1363600000;
    atpos[25][2] =     0.5454500000;
    
    attyp[26] =  1;
    atpos[26][0] =     0.5000000000;
    atpos[26][1] =     0.7727300000;
    atpos[26][2] =     0.5909100000;
    
    attyp[27] =  1;
    atpos[27][0] =     0.3333300000;
    atpos[27][1] =     0.2727300000;
    atpos[27][2] =     0.5909100000;
    
    attyp[28] =  1;
    atpos[28][0] =     0.0000000000;
    atpos[28][1] =     0.9090900000;
    atpos[28][2] =     0.6363600000;
    
    attyp[29] =  1;
    atpos[29][0] =     0.8333300000;
    atpos[29][1] =     0.4090900000;
    atpos[29][2] =     0.6363600000;
    
    attyp[30] =  1;
    atpos[30][0] =     0.5000000000;
    atpos[30][1] =     0.0454500000;
    atpos[30][2] =     0.6818200000;
    
    attyp[31] =  1;
    atpos[31][0] =     0.3333300000;
    atpos[31][1] =     0.5454500000;
    atpos[31][2] =     0.6818200000;
    
    attyp[32] =  1;
    atpos[32][0] =     0.0000000000;
    atpos[32][1] =     0.1818200000;
    atpos[32][2] =     0.7272700000;
    
    attyp[33] =  1;
    atpos[33][0] =     0.8333300000;
    atpos[33][1] =     0.6818200000;
    atpos[33][2] =     0.7272700000;
    
    attyp[34] =  1;
    atpos[34][0] =     0.5000000000;
    atpos[34][1] =     0.3181800000;
    atpos[34][2] =     0.7727300000;
    
    attyp[35] =  1;
    atpos[35][0] =     0.3333300000;
    atpos[35][1] =     0.8181800000;
    atpos[35][2] =     0.7727300000;
    
    attyp[36] =  1;
    atpos[36][0] =     0.0000000000;
    atpos[36][1] =     0.4545500000;
    atpos[36][2] =     0.8181800000;
    
    attyp[37] =  1;
    atpos[37][0] =     0.8333300000;
    atpos[37][1] =     0.9545500000;
    atpos[37][2] =     0.8181800000;
    
    attyp[38] =  1;
    atpos[38][0] =     0.3333300000;
    atpos[38][1] =     0.0909100000;
    atpos[38][2] =     0.8636400000;
    
    attyp[39] =  1;
    atpos[39][0] =     0.5000000000;
    atpos[39][1] =     0.5909100000;
    atpos[39][2] =     0.8636400000;
    
    attyp[40] =  1;
    atpos[40][0] =     0.0000000000;
    atpos[40][1] =     0.7272700000;
    atpos[40][2] =     0.9090900000;
    
    attyp[41] =  1;
    atpos[41][0] =     0.8333300000;
    atpos[41][1] =     0.2272700000;
    atpos[41][2] =     0.9090900000;
    
    attyp[42] =  1;
    atpos[42][0] =     0.5000000000;
    atpos[42][1] =     0.8636400000;
    atpos[42][2] =     0.9545500000;
    
    attyp[43] =  1;
    atpos[43][0] =     0.3333300000;
    atpos[43][1] =     0.3636400000;
    atpos[43][2] =     0.9545500000;
    
    initialized = 1;
    break;

  case 2:
    nucell = 44;
    ntype  = 1;
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");

    latvec[0][0] = sqrt(ca*ca + 1.);
    latvec[1][1] = sqrt(3.);
    latvec[2][2] =  11./sqrt(1.+1./(ca*ca));
    
    attyp[ 0] =  1;
    atpos[ 0][1] =     0.0000000000;
    atpos[ 0][0] =     0.0000000000;
    atpos[ 0][2] =     0.0000000000;
    
    attyp[ 1] =  1;
    atpos[ 1][1] =     0.8333300000;
    atpos[ 1][0] =     0.5000000000;
    atpos[ 1][2] =     0.0000000000;
    
    attyp[ 2] =  1;
    atpos[ 2][1] =     0.5000000000;
    atpos[ 2][0] =     0.1363600000;
    atpos[ 2][2] =     0.0454500000;
    
    attyp[ 3] =  1;
    atpos[ 3][1] =     0.3333300000;
    atpos[ 3][0] =     0.6363600000;
    atpos[ 3][2] =     0.0454500000;
    
    attyp[ 4] =  1;
    atpos[ 4][1] =     0.0000000000;
    atpos[ 4][0] =     0.2727300000;
    atpos[ 4][2] =     0.0909100000;
    
    attyp[ 5] =  1;
    atpos[ 5][1] =     0.8333300000;
    atpos[ 5][0] =     0.7727300000;
    atpos[ 5][2] =     0.0909100000;
    
    attyp[ 6] =  1;
    atpos[ 6][1] =     0.5000000000;
    atpos[ 6][0] =     0.4090900000;
    atpos[ 6][2] =     0.1363600000;
    
    attyp[ 7] =  1;
    atpos[ 7][1] =     0.3333300000;
    atpos[ 7][0] =     0.9090900000;
    atpos[ 7][2] =     0.1363600000;
    
    attyp[ 8] =  1;
    atpos[ 8][1] =     0.8333300000;
    atpos[ 8][0] =     0.0454500000;
    atpos[ 8][2] =     0.1818200000;
    
    attyp[ 9] =  1;
    atpos[ 9][1] =     0.0000000000;
    atpos[ 9][0] =     0.5454500000;
    atpos[ 9][2] =     0.1818200000;
    
    attyp[10] =  1;
    atpos[10][1] =     0.5000000000;
    atpos[10][0] =     0.6818200000;
    atpos[10][2] =     0.2272700000;
    
    attyp[11] =  1;
    atpos[11][1] =     0.3333300000;
    atpos[11][0] =     0.1818200000;
    atpos[11][2] =     0.2272700000;
    
    attyp[12] =  1;
    atpos[12][1] =     0.0000000000;
    atpos[12][0] =     0.8181800000;
    atpos[12][2] =     0.2727300000;
    
    attyp[13] =  1;
    atpos[13][1] =     0.8333300000;
    atpos[13][0] =     0.3181800000;
    atpos[13][2] =     0.2727300000;
    
    attyp[14] =  1;
    atpos[14][1] =     0.5000000000;
    atpos[14][0] =     0.9545500000;
    atpos[14][2] =     0.3181800000;
    
    attyp[15] =  1;
    atpos[15][1] =     0.3333300000;
    atpos[15][0] =     0.4545500000;
    atpos[15][2] =     0.3181800000;
    
    attyp[16] =  1;
    atpos[16][1] =     0.0000000000;
    atpos[16][0] =     0.0909100000;
    atpos[16][2] =     0.3636400000;
    
    attyp[17] =  1;
    atpos[17][1] =     0.8333300000;
    atpos[17][0] =     0.5909100000;
    atpos[17][2] =     0.3636400000;
    
    attyp[18] =  1;
    atpos[18][1] =     0.5000000000;
    atpos[18][0] =     0.2272700000;
    atpos[18][2] =     0.4090900000;
    
    attyp[19] =  1;
    atpos[19][1] =     0.3333300000;
    atpos[19][0] =     0.7272700000;
    atpos[19][2] =     0.4090900000;
    
    attyp[20] =  1;
    atpos[20][1] =     0.0000000000;
    atpos[20][0] =     0.3636400000;
    atpos[20][2] =     0.4545500000;
    
    attyp[21] =  1;
    atpos[21][1] =     0.8333300000;
    atpos[21][0] =     0.8636400000;
    atpos[21][2] =     0.4545500000;
    
    attyp[22] =  1;
    atpos[22][1] =     0.3333300000;
    atpos[22][0] =     0.0000000000;
    atpos[22][2] =     0.5000000000;
    
    attyp[23] =  1;
    atpos[23][1] =     0.5000000000;
    atpos[23][0] =     0.5000000000;
    atpos[23][2] =     0.5000000000;
    
    attyp[24] =  1;
    atpos[24][1] =     0.0000000000;
    atpos[24][0] =     0.6363600000;
    atpos[24][2] =     0.5454500000;
    
    attyp[25] =  1;
    atpos[25][1] =     0.8333300000;
    atpos[25][0] =     0.1363600000;
    atpos[25][2] =     0.5454500000;
    
    attyp[26] =  1;
    atpos[26][1] =     0.5000000000;
    atpos[26][0] =     0.7727300000;
    atpos[26][2] =     0.5909100000;
    
    attyp[27] =  1;
    atpos[27][1] =     0.3333300000;
    atpos[27][0] =     0.2727300000;
    atpos[27][2] =     0.5909100000;
    
    attyp[28] =  1;
    atpos[28][1] =     0.0000000000;
    atpos[28][0] =     0.9090900000;
    atpos[28][2] =     0.6363600000;
    
    attyp[29] =  1;
    atpos[29][1] =     0.8333300000;
    atpos[29][0] =     0.4090900000;
    atpos[29][2] =     0.6363600000;
    
    attyp[30] =  1;
    atpos[30][1] =     0.5000000000;
    atpos[30][0] =     0.0454500000;
    atpos[30][2] =     0.6818200000;
    
    attyp[31] =  1;
    atpos[31][1] =     0.3333300000;
    atpos[31][0] =     0.5454500000;
    atpos[31][2] =     0.6818200000;
    
    attyp[32] =  1;
    atpos[32][1] =     0.0000000000;
    atpos[32][0] =     0.1818200000;
    atpos[32][2] =     0.7272700000;
    
    attyp[33] =  1;
    atpos[33][1] =     0.8333300000;
    atpos[33][0] =     0.6818200000;
    atpos[33][2] =     0.7272700000;
    
    attyp[34] =  1;
    atpos[34][1] =     0.5000000000;
    atpos[34][0] =     0.3181800000;
    atpos[34][2] =     0.7727300000;
    
    attyp[35] =  1;
    atpos[35][1] =     0.3333300000;
    atpos[35][0] =     0.8181800000;
    atpos[35][2] =     0.7727300000;
    
    attyp[36] =  1;
    atpos[36][1] =     0.0000000000;
    atpos[36][0] =     0.4545500000;
    atpos[36][2] =     0.8181800000;
    
    attyp[37] =  1;
    atpos[37][1] =     0.8333300000;
    atpos[37][0] =     0.9545500000;
    atpos[37][2] =     0.8181800000;
    
    attyp[38] =  1;
    atpos[38][1] =     0.3333300000;
    atpos[38][0] =     0.0909100000;
    atpos[38][2] =     0.8636400000;
    
    attyp[39] =  1;
    atpos[39][1] =     0.5000000000;
    atpos[39][0] =     0.5909100000;
    atpos[39][2] =     0.8636400000;
    
    attyp[40] =  1;
    atpos[40][1] =     0.0000000000;
    atpos[40][0] =     0.7272700000;
    atpos[40][2] =     0.9090900000;
    
    attyp[41] =  1;
    atpos[41][1] =     0.8333300000;
    atpos[41][0] =     0.2272700000;
    atpos[41][2] =     0.9090900000;
    
    attyp[42] =  1;
    atpos[42][1] =     0.5000000000;
    atpos[42][0] =     0.8636400000;
    atpos[42][2] =     0.9545500000;
    
    attyp[43] =  1;
    atpos[43][1] =     0.3333300000;
    atpos[43][0] =     0.3636400000;
    atpos[43][2] =     0.9545500000;
    
    initialized = 1;
    break;

  default:
    break;
  }
return;
}
/* -------------------------------------------------------------------------- */
