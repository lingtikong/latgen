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
  printf("   1. (001);         5. Graphene;\n");
  printf("   2. (100);         6. Graphite;\n");
  printf("   3. (110);         7. (10 -1);\n");
  printf("   4. (-110);        8. (112);\n");
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
  printf("   1. [010] along x;          2. [101] along y;\n");
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
    break;

  case 2:
    latvec[0][0] =  c_ang;
    latvec[0][1] = -sqrt(1. - c_ang*c_ang);
    latvec[1][1] = sqrt(ca*ca + 1.);

    latvec[2][2] = 41./sqrt(4./3. + 1./(ca*ca));
    break;

  default:
    break;
  }

  if (surftype >= 1 && surftype <= 2){
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
