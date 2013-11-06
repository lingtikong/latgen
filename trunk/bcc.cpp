#include "bcc.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#define MAXLINE 256

using namespace std;

/* -----------------------------------------------------------------------------
 * To select the orientation of the lattice
 * -------------------------------------------------------------------------- */
BCC::BCC() : lattice()
{
  char str[MAXLINE];
  // print out the menu
  alat = 1.;
  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");
  printf("Please input the lattice constant of the BCC lattice [1.]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) alat = numeric(strtok(str, " \t\n\r\f"));

  int orient = 1;
  printf("Please selection the orientation of the BCC lattice:\n");
  printf("   1. (001);\n");
  printf("   2. (110);\n");
  printf("   3. (111);\n");
  printf("   4. (112);\n");
  printf("   5. primitive cell;\n");
  printf("Your choice [%d]: ", orient);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) orient = inumeric(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", orient);
  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");
  
  // initialize according to orientation
  initialized = 0;
  switch (orient){
  case 1:
    BCC001();
    break;
  case 2:
    BCC110();
    break;
  case 3:
    BCC111();
    break;
  case 4:
    BCC112();
    break;
  case 5:
    Primitive();
    break;
  default:
    break;
  }

}

/* -----------------------------------------------------------------------------
 * Deconstructor does nothing
 * -------------------------------------------------------------------------- */
BCC::~BCC()
{

}

/* -----------------------------------------------------------------------------
 * Initialize for (001) orientation
 * -------------------------------------------------------------------------- */
void BCC::BCC001()
{
  char str[MAXLINE];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");
  printf("Please selection the type of BCC(001) surface:\n");
  printf("   1. conventional orientation;\n");
  printf("   2. B2 structure;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = inumeric(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", surftype);
  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");

  for (int i = 0; i < 3; ++i)
  for (int j = 0; j < 3; ++j) latvec[i][j] = 0.;

  memory->create(name,9,"BCC:name");
  strcpy(name, "BCC(001)");

  // initialize according to surface type
  switch (surftype){
  case 1:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1.;
    latvec[2][2] = 1.;

    memory->create(atpos, nucell, 3, "BCC001_atpos");
    memory->create(attyp, nucell, "BCC:attyp");
    
    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.5;

    initialized = 1;
    break;
  case 2:
    nucell = 2;
    ntype  = 2;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1.;
    latvec[2][2] = 1.;

    memory->create(atpos, nucell, 3, "BCC001_atpos");
    memory->create(attyp, nucell, "BCC:attyp");
    
    attyp[0] = 1;
    attyp[1] = 2;

    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.5;

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
void BCC::BCC110()
{
  char str[MAXLINE];
  int surftype=1;
  // print out the menu
  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");
  printf("Please selection the type of BCC(110) surface:\n");
  printf("   1. orthogonal, long side along x\n");
  printf("   2. orthogonal, long side along y\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = inumeric(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", surftype);
  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");
  
  for (int i = 0; i < 3; ++i)
  for (int j = 0; j < 3; ++j) latvec[i][j] = 0.;

  memory->create(name,9,"BCC:name");
  strcpy(name, "BCC(110)");

  // initialize according to surface type
  switch (surftype){
  case 1:
    nucell = 4;
    ntype  = 1;
    
    latvec[0][0] = sqrt(2.);
    latvec[1][1] = 1.;
    latvec[2][2] = sqrt(2.);

    memory->create(atpos, nucell, 3, "BCC110_atpos");
    memory->create(attyp, nucell, "BCC:attyp");
    
    for (int i = 0; i < nucell; ++i) attyp[i] = 1;

    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;

    atpos[2][0] = 0.5;
    atpos[2][1] = 0.0;
    atpos[2][2] = 0.5;

    atpos[3][0] = 0.0;
    atpos[3][1] = 0.5;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  case 2:
    nucell = 4;
    ntype  = 1;
    
    latvec[0][0] = 1.;
    latvec[1][1] = sqrt(2.);
    latvec[2][2] = sqrt(2.);

    memory->create(atpos, nucell, 3, "BCC110_atpos");
    memory->create(attyp, nucell, "BCC:attyp");
    
    for (int i = 0; i < nucell; ++i) attyp[i] = 1;

    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;

    atpos[2][0] = 0.;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.5;

    atpos[3][0] = 0.5;
    atpos[3][1] = 0.0;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* -----------------------------------------------------------------------------
 * Initialize for (111) orientation
 * -------------------------------------------------------------------------- */
void BCC::BCC111()
{
  char str[MAXLINE];
  int surftype =1;
  // print out the menu
  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");
  printf("Please selection the type of BCC(111) surface:\n");
  printf("   1. U = [1-10], V = [10-1]; U // x;\n");
  printf("   2. U = [1-10], V = [10-1]; V // y;\n");
  printf("   3. U = [10-1], V = [1-21]; U // x;\n");
  printf("   4. U = [1-21], V = [10-1]; U // x;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = inumeric(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", surftype);
  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");
  
  for (int i = 0; i < 3; ++i)
  for (int j = 0; j < 3; ++j) latvec[i][j] = 0.;

  memory->create(name,9,"BCC:name");
  strcpy(name, "BCC(111)");

  // initialize according to surface type
  switch (surftype){
  case 1:
    nucell = 3;
    ntype  = 1;
    
    latvec[0][0] =  sqrt(2.);
    latvec[1][0] = -sqrt(0.5);
    latvec[1][1] =  sqrt(1.5);
    latvec[2][2] =  sqrt(0.75);
    
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[0] =  1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] =  1;
    atpos[1][0] = 1./3.;
    atpos[1][1] = 2./3.;
    atpos[1][2] = 1./3.;
    
    attyp[2] =  1;
    atpos[2][0] = 2./3.;
    atpos[2][1] = 1./3.;
    atpos[2][2] = 2./3.;
    
    initialized = 1;
    break;

   case 2:
    nucell = 3;
    ntype  = 1;
    
    latvec[0][0] =  sqrt(1.5);
    latvec[0][1] = -sqrt(0.5);
    latvec[1][1] =  sqrt(2.0);
    latvec[2][2] =  sqrt(0.75);
    
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[0] =  1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] =  1;
    atpos[1][0] = 1./3.;
    atpos[1][1] = 2./3.;
    atpos[1][2] = 1./3.;
    
    attyp[2] =  1;
    atpos[2][0] = 2./3.;
    atpos[2][1] = 1./3.;
    atpos[2][2] = 2./3.;
    
    initialized = 1;
    break;

  case 3:
    nucell =    6;
    ntype  =  1;
    
    latvec[0][0] = sqrt(2.);
    latvec[1][1] = sqrt(6.);
    latvec[2][2] = sqrt(0.75);
    
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[0] =  1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] =  1;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] =  1;
    atpos[2][0] = 0.5;
    atpos[2][1] = 1./6.;
    atpos[2][2] = 1./3.;
    
    attyp[3] =  1;
    atpos[3][0] = 0.;
    atpos[3][1] = 2./3.;
    atpos[3][2] = 1./3.;
    
    attyp[4] =  1;
    atpos[4][0] = 0.;
    atpos[4][1] = 1./3.;
    atpos[4][2] = 2./3.;
    
    attyp[5] =  1;
    atpos[5][0] = 0.5;
    atpos[5][1] = 5./6.;
    atpos[5][2] = 2./3.;
    
    initialized = 1;
    break;

  case 4:
    nucell =    6;
    ntype  =  1;
    
    latvec[0][0] = sqrt(6.);
    latvec[1][1] = sqrt(2.);
    latvec[2][2] = sqrt(0.75);
    
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[0] =  1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] =  1;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;
    
    attyp[2] =  1;
    atpos[2][0] = 1./6.;
    atpos[2][1] = 0.5;
    atpos[2][2] = 1./3.;
    
    attyp[3] =  1;
    atpos[3][0] = 2./3.;
    atpos[3][1] = 0.;
    atpos[3][2] = 1./3.;
    
    attyp[4] =  1;
    atpos[4][0] = 1./3.;
    atpos[4][1] = 0.;
    atpos[4][2] = 2./3.;
    
    attyp[5] =  1;
    atpos[5][0] = 5./6.;
    atpos[5][1] = 0.5;
    atpos[5][2] = 2./3.;
    
    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* -----------------------------------------------------------------------------
 * Initialize for (112) orientation
 * -------------------------------------------------------------------------- */
void BCC::BCC112()
{
  char str[MAXLINE];
  int surftype =1;
  // print out the menu
  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");
  printf("Please selection the type of BCC(111) surface:\n");
  printf("   1. U = [1-10], V = [.5,.5,-.5]; U // x\n");
  printf("   2. U = [.5,.5,-.5], V = [-110]; U // x\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = inumeric(strtok(str, " \t\n\r\f"));
  printf("You   selected : %d", surftype);
  printf("\n"); for (int i = 0; i < 14; ++i) printf("====="); printf("\n");
  
  for (int i = 0; i < 3; ++i)
  for (int j = 0; j < 3; ++j) latvec[i][j] = 0.;

  memory->create(name,9,"BCC:name");
  strcpy(name, "BCC(112)");

  // initialize according to surface type
  switch (surftype){
  case 1:
    nucell = 6;
    ntype  = 1;
    
    latvec[0][0] = sqrt(2.);
    latvec[1][1] = sqrt(0.75);
    latvec[2][2] = sqrt(6.);
    
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[0] =  1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] =  1;
    atpos[1][0] = 0.5;
    atpos[1][1] = 2./3.;
    atpos[1][2] = 1./6.;
    
    attyp[2] =  1;
    atpos[2][0] = 0.;
    atpos[2][1] = 1./3.;
    atpos[2][2] = 1./3.;
    
    attyp[3] =  1;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.0;
    atpos[3][2] = 0.5;
    
    attyp[4] =  1;
    atpos[4][0] = 0.;
    atpos[4][1] = 2./3.;
    atpos[4][2] = 2./3.;
    
    attyp[5] =  1;
    atpos[5][0] = 0.5;
    atpos[5][1] = 1./3.;
    atpos[5][2] = 5./6.;
    
    initialized = 1;
    break;

  case 2:
    nucell = 6;
    ntype  = 1;
    
    latvec[0][0] = sqrt(0.75);
    latvec[1][1] = sqrt(2.);
    latvec[2][2] = sqrt(6.);
    
    atpos = memory->create(atpos,nucell,3,"atpos");
    attyp = memory->create(attyp,nucell,"attyp");
    
    attyp[0] =  1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;
    
    attyp[1] =  1;
    atpos[1][0] = 2./3.;
    atpos[1][1] = 0.5;
    atpos[1][2] = 1./6.;
    
    attyp[2] =  1;
    atpos[2][0] = 1./3.;
    atpos[2][1] = 0.;
    atpos[2][2] = 1./3.;
    
    attyp[3] =  1;
    atpos[3][0] = 0.0;
    atpos[3][1] = 0.5;
    atpos[3][2] = 0.5;
    
    attyp[4] =  1;
    atpos[4][0] = 2./3.;
    atpos[4][1] = 0.;
    atpos[4][2] = 2./3.;
    
    attyp[5] =  1;
    atpos[5][0] = 1./3.;
    atpos[5][1] = 0.5;
    atpos[5][2] = 5./6.;
    
    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* -----------------------------------------------------------------------------
 * Initialize for Primitive cell
 * -------------------------------------------------------------------------- */
void BCC::Primitive()
{
  memory->create(name,9,"BCC:name");
  strcpy(name, "BCC-prim");

  nucell = 1;
  ntype  = 1;
  
  latvec[0][0] = -0.5;
  latvec[0][1] =  0.5;
  latvec[0][2] =  0.5;
  latvec[1][0] =  0.5;
  latvec[1][1] = -0.5;
  latvec[1][2] =  0.5;
  latvec[2][0] =  0.5;
  latvec[2][1] =  0.5;
  latvec[2][2] = -0.5;

  memory->create(atpos, nucell, 3, "BCC001_atpos");
  memory->create(attyp, nucell, "BCC:attyp");
  
  for (int i = 0; i < nucell; ++i) attyp[i] = 1;
  atpos[0][0] = 0.;
  atpos[0][1] = 0.;
  atpos[0][2] = 0.;

  initialized = 1;
return;
}

/* --------------------------------------------------------------------------- */
