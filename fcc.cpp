#include "fcc.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#define MAX_LINE_LENGTH 256

using namespace std;

/* ----------------------------------------------------------------------
   To select the orientation of the lattice
------------------------------------------------------------------------- */
FCC::FCC() : lattice()
{
  char str[MAX_LINE_LENGTH];
  alat = 1.;
  // print out the menu
  printf("\n"); for (int i=0; i<70; i++) printf("="); printf("\n");
  printf("Please input the lattice constant of the FCC lattice [1.]: ");
  if (count_words(gets(str)) > 0) alat = atof(strtok(str, " \t\n\r\f"));

  int orient = 3;
  printf("Please selection the orientation of the FCC lattice:\n");
  printf("   1. (001);                5. Diamond primitive cell;\n");
  printf("   2. (110);                6. Diamond conventional cell;\n");
  printf("   3. (111);                7. NaCl primitive cell;\n");
  printf("   4. Primitive cell;       8. NaCl convertion cell;\n");
  for (int i=0; i<70; i++) printf("-");
  printf("\nYour  choice [3]: ");
  if (count_words(gets(str)) > 0) orient = atoi(strtok(str, " \t\n\r\f"));
  printf("You selected: %d", orient);
  printf("\n"); for (int i=0; i<70; i++) printf("="); printf("\n");
  
  // initialize according to orientation
  initialized = 0;
  switch (orient){
  case 1:
    FCC001();
    break;
  case 2:
    FCC110();
    break;
  case 3:
    FCC111();
    break;
  case 4:
    Primitive();
    break;
  case 5:
    DiamondPrim();
    break;
  case 6:
    DiamondConv();
    break;
  case 7:
    NaClPrim();
    break;
  case 8:
    NaClConv();
    break;
  default:
    break;
  }

}

/* ----------------------------------------------------------------------
   Deconstructor does nothing
------------------------------------------------------------------------- */
FCC::~FCC()
{

}

/* ----------------------------------------------------------------------
   Initialize for (001) orientation
------------------------------------------------------------------------- */
void FCC::FCC001()
{
  char str[MAX_LINE_LENGTH];
  int surftype = 2;
  // print out the menu
  printf("\n"); for (int i=0; i<70; i++) printf("="); printf("\n");
  printf("Please selection the type of FCC(001) surface:\n");
  printf("   1. primitive, horizental orientation;\n");
  printf("   2. conventional orientation;\n");
  printf("Your  choice [2]: ");
  if (count_words(gets(str)) > 0) surftype = atoi(strtok(str, " \t\n\r\f"));
  printf("You selected: %d", surftype);
  printf("\n"); for (int i=0; i<70; i++) printf("="); printf("\n");
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = new char[9];
  strcpy(name, "FCC(001)");

  // initialize according to surface type
  switch (surftype){
  case 1:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] = 1./sqrt(2.);
    latvec[1][1] = 1./sqrt(2.);
    latvec[2][2] = 1.;

    atpos = memory->create(atpos,nucell, 3, "FCC001_atpos");
    attyp = new int[nucell];
    layer = new int[nucell];
  
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;
    layer[0]    = 0;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    layer[1]    = 1;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.5;

    initialized = 1;
    break;
  case 2:
    nucell = 4;
    ntype  = 1;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1.;
    latvec[2][2] = 1.;

    atpos = memory->create(atpos,nucell, 3, "FCC001_atpos");
    attyp = new int[nucell];
    layer = new int[nucell];
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;
    layer[0]    = 0;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    layer[1]    = 0;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;

    layer[2]    = 1;
    atpos[2][0] = 0.;
    atpos[2][1] = 0.5;
    atpos[2][2] = 0.5;

    layer[3]    = 1;
    atpos[3][0] = 0.5;
    atpos[3][1] = 0.;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* ----------------------------------------------------------------------
   Initialize for (110) orientation
------------------------------------------------------------------------- */
void FCC::FCC110()
{
  char str[MAX_LINE_LENGTH];
  int surftype = 1;
  // print out the menu
  printf("\n"); for (int i=0; i<70; i++) printf("="); printf("\n");
  printf("Please selection the type of FCC(110) surface:\n");
  printf("   1. orthogonal, long side along x\n");
  printf("   2. orthogonal, long side along y\n");
  printf("Your  choice [1]: ");
  if (count_words(gets(str)) > 0) surftype = atoi(strtok(str, " \t\n\r\f"));
  printf("You selected: %d", surftype);
  printf("\n"); for (int i=0; i<70; i++) printf("="); printf("\n");
  
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = new char[9];
  strcpy(name, "FCC(110)");

  // initialize according to surface type
  switch (surftype){
  case 2:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] = 0.5*sqrt(2.);
    latvec[1][1] = 1.;
    latvec[2][2] = 0.5*sqrt(2.);

    atpos = memory->create(atpos,nucell, 3, "FCC110_atpos");
    attyp = new int[nucell];
    layer = new int[nucell];
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;
    layer[0]    = 0;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    layer[1]    = 1;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.5;

    initialized = 1;
    break;
  case 1:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] = 1.;
    latvec[1][1] = sqrt(2.)*0.5;
    latvec[2][2] = sqrt(2.)*0.5;

    atpos = memory->create(atpos,nucell, 3, "FCC110_atpos");
    attyp = new int[nucell];
    layer = new int[nucell];
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;

    layer[0]    = 0;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    layer[1]    = 1;
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

/* ----------------------------------------------------------------------
   Initialize for (111) orientation
------------------------------------------------------------------------- */
void FCC::FCC111()
{
  char str[MAX_LINE_LENGTH];
  int surftype = 5;
  // print out the menu
  printf("\n"); for (int i=0; i<70; i++) printf("="); printf("\n");
  printf("Please selection the type of FCC(111) surface:\n");
  printf("   1. hexgonal U along x, 60 deg;\n");
  printf("   2. hexgonal V along y, 60 deg;\n");
  printf("   3. hexgonal U along x, 120 deg;\n");
  printf("   4. hexgonal V along y, 120 deg;\n");
  printf("   5. orthogonal long side along x;\n");
  printf("   6. orthogonal long side along y;\n");
  printf("Your  choice [5]: ");
  if (count_words(gets(str)) > 0) surftype = atoi(strtok(str, " \t\n\r\f"));
  printf("You selected: %d", surftype);
  printf("\n"); for (int i=0; i<70; i++) printf("="); printf("\n");
  
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = new char[9];
  strcpy(name, "FCC(111)");

  // initialize according to surface type
  switch (surftype){
  case 1:
    nucell = 3;
    ntype  = 1;
    
    latvec[0][0] = sqrt(0.5);
    latvec[1][0] = sqrt(0.125);
    latvec[1][1] = sqrt(0.375);
    latvec[2][2] = sqrt(3.);

    atpos = memory->create(atpos,nucell, 3, "FCC111_atpos");
    attyp = new int[nucell];
    layer = new int[nucell];
    
    for (int i=0; i<nucell; i++) attyp[i] = 1;

    layer[0]    = 0;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    layer[1]    = 1;
    atpos[1][0] = 1./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 1./3.;

    layer[2]    = 2;
    atpos[2][0] = 2./3.;
    atpos[2][1] = 2./3.;
    atpos[2][2] = 2./3.;

    initialized = 1;
    break;
  case 2:
    nucell = 3;
    ntype  = 1;

    latvec[0][0] = sqrt(0.375);
    latvec[0][1] = sqrt(0.125);
    latvec[1][1] = sqrt(0.5);
    latvec[2][2] = sqrt(3.);

    atpos = memory->create(atpos,nucell, 3, "FCC111_atpos");
    attyp = new int[nucell];
    layer = new int[nucell];

    for (int i=0; i<nucell; i++) attyp[i] = 1;

    layer[0]    = 0;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    layer[1]    = 1;
    atpos[1][0] = 1./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 1./3.;

    layer[2]    = 2;
    atpos[2][0] = 2./3.;
    atpos[2][1] = 2./3.;
    atpos[2][2] = 2./3.;

    initialized = 1;
    break;
  case 3:
    nucell = 3;
    ntype  = 1;

    latvec[0][0] =  sqrt(0.5);
    latvec[1][0] = -sqrt(0.125);
    latvec[1][1] =  sqrt(0.375);
    latvec[2][2] =  sqrt(3.);

    atpos = memory->create(atpos,nucell, 3, "FCC111_atpos");
    attyp = new int[nucell];
    layer = new int[nucell];

    for (int i=0; i<nucell; i++) attyp[i] = 1;

    layer[0]    = 0;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    layer[1]    = 1;
    atpos[1][0] = 2./3.;
    atpos[1][1] = 1./3.;
    atpos[1][2] = 1./3.;

    layer[2]    = 2;
    atpos[2][0] = 1./3.;
    atpos[2][1] = 2./3.;
    atpos[2][2] = 2./3.;

    initialized = 1;
    break;
  case 4:
    nucell = 3;
    ntype  = 1;

    latvec[0][0] =  sqrt(0.375);
    latvec[0][1] = -sqrt(0.125);
    latvec[1][1] =  sqrt(0.5);
    latvec[2][2] =  sqrt(3.);

    atpos = memory->create(atpos,nucell, 3, "FCC111_atpos");
    attyp = new int[nucell];
    layer = new int[nucell];

    for (int i=0; i<nucell; i++) attyp[i] = 1;

    layer[0]    = 0;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    layer[1]    = 1;
    atpos[1][0] = 1./3.;
    atpos[1][1] = 2./3.;
    atpos[1][2] = 1./3.;

    layer[2]    = 2;
    atpos[2][0] = 2./3.;
    atpos[2][1] = 1./3.;
    atpos[2][2] = 2./3.;

    initialized = 1;
    break;
  case 5:
    nucell = 6;
    ntype  = 1;

    latvec[0][0] =  sqrt(1.5);
    latvec[1][1] =  sqrt(0.5);
    latvec[2][2] =  sqrt(3.);

    atpos = memory->create(atpos,nucell, 3, "FCC111_atpos");
    attyp = new int[nucell];
    layer = new int[nucell];

    for (int i=0; i<nucell; i++) attyp[i] = 1;

    layer[0]    = 0;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    layer[1]    = 0;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;

    layer[2]    = 1;
    atpos[2][0] = 1./6.;
    atpos[2][1] = 0.5;
    atpos[2][2] = 1./3.;

    layer[3]    = 1;
    atpos[3][0] = 2./3.;
    atpos[3][1] = 1.;
    atpos[3][2] = 1./3.;

    layer[4]    = 2;
    atpos[4][0] = 1./3.;
    atpos[4][1] = 0.;
    atpos[4][2] = 2./3.;

    layer[5]    = 2;
    atpos[5][0] = 5./6.;
    atpos[5][1] = 0.5;
    atpos[5][2] = 2./3.;

    initialized = 1;
    break;
  case 6:
    nucell = 6;
    ntype  = 1;

    latvec[0][0] = sqrt(0.5);
    latvec[1][1] = sqrt(1.5);
    latvec[2][2] = sqrt(3.);

    atpos = memory->create(atpos,nucell, 3, "FCC111_atpos");
    attyp = new int[nucell];
    layer = new int[nucell];

    for (int i=0; i<nucell; i++) attyp[i] = 1;

    layer[0]    = 0;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    layer[1]    = 0;
    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.;

    layer[2]    = 1;
    atpos[2][0] = 0.5;
    atpos[2][1] = 1./6.;
    atpos[2][2] = 1./3.;

    layer[3]    = 1;
    atpos[3][0] = 1.;
    atpos[3][1] = 2./3.;
    atpos[3][2] = 1./3.;

    layer[4]    = 2;
    atpos[4][0] = 0.;
    atpos[4][1] = 1./3.;
    atpos[4][2] = 2./3.;

    layer[5]    = 2;
    atpos[5][0] = 0.5;
    atpos[5][1] = 5./6.;
    atpos[5][2] = 2./3.;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* ----------------------------------------------------------------------
   Initialize for (001) orientation
------------------------------------------------------------------------- */
void FCC::Primitive()
{
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = new char[9];
  strcpy(name, "FCC-prim");

  nucell = 1;
  ntype  = 1;
  
  latvec[0][1] =  0.5;
  latvec[0][2] =  0.5;
  latvec[1][0] =  0.5;
  latvec[1][2] =  0.5;
  latvec[2][0] =  0.5;
  latvec[2][1] =  0.5;

  atpos = memory->create(atpos,nucell, 3, "Primitive_atpos");
  attyp = new int[nucell]; layer = new int[nucell];
  
  for (int i=0; i<nucell; i++) attyp[i] = 1;
  layer[0]    = 0;
  atpos[0][0] = 0.;
  atpos[0][1] = 0.;
  atpos[0][2] = 0.;

  initialized = 1;
return;
}

/* ----------------------------------------------------------------------
   Initialize Diamond primitive cell
------------------------------------------------------------------------- */
void FCC::DiamondPrim()
{
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = new char[9];
  strcpy(name, "Dia-prim");

  nucell = 2;
  ntype  = 1;
  
  latvec[0][1] =  0.5;
  latvec[0][2] =  0.5;
  latvec[1][0] =  0.5;
  latvec[1][2] =  0.5;
  latvec[2][0] =  0.5;
  latvec[2][1] =  0.5;

  atpos = memory->create(atpos,nucell, 3, "Primitive_atpos");
  attyp = new int[nucell]; layer = new int[nucell];
  
  for (int i=0; i<nucell; i++) attyp[i] = 1;
  layer[0]    = 0;
  atpos[0][0] = 0.;
  atpos[0][1] = 0.;
  atpos[0][2] = 0.;

  layer[1]    = 1;
  atpos[1][0] = 0.25;
  atpos[1][1] = 0.25;
  atpos[1][2] = 0.25;

  initialized = 1;
return;
}

/* ----------------------------------------------------------------------
   Initialize Diamond Conventional cell
------------------------------------------------------------------------- */
void FCC::DiamondConv()
{
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = new char[9];
  strcpy(name, "Dia-conv");

  nucell = 8;
  ntype  = 1;
  
  latvec[0][0] = 1.;
  latvec[1][1] = 1.;
  latvec[2][2] = 1.;

  atpos = memory->create(atpos,nucell, 3, "Primitive_atpos");
  attyp = new int[nucell]; layer = new int[nucell];
  
  for (int i=0; i<nucell; i++) attyp[i] = 1;
  layer[0]    = 0;
  atpos[0][0] = 0.;
  atpos[0][1] = 0.;
  atpos[0][2] = 0.;

  layer[1]    = 0;
  atpos[1][0] = 0.5;
  atpos[1][1] = 0.5;
  atpos[1][2] = 0.;

  layer[2]    = 1;
  atpos[2][0] = 0.25;
  atpos[2][1] = 0.25;
  atpos[2][2] = 0.25;

  layer[3]    = 1;
  atpos[3][0] = 0.75;
  atpos[3][1] = 0.75;
  atpos[3][2] = 0.25;

  layer[4]    = 2;
  atpos[4][0] = 0.5;
  atpos[4][1] = 0.0;
  atpos[4][2] = 0.5;

  layer[5]    = 2;
  atpos[5][0] = 0.0;
  atpos[5][1] = 0.5;
  atpos[5][2] = 0.5;

  layer[6]    = 3;
  atpos[6][0] = 0.75;
  atpos[6][1] = 0.25;
  atpos[6][2] = 0.75;

  layer[7]    = 3;
  atpos[7][0] = 0.25;
  atpos[7][1] = 0.75;
  atpos[7][2] = 0.75;

  initialized = 1;
return;
}

/* ----------------------------------------------------------------------
   Initialize Diamond primitive cell
------------------------------------------------------------------------- */
void FCC::NaClPrim()
{
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = new char[10];
  strcpy(name, "NaCl-prim");

  nucell = 2;
  ntype  = 2;
  
  latvec[0][1] =  0.5;
  latvec[0][2] =  0.5;
  latvec[1][0] =  0.5;
  latvec[1][2] =  0.5;
  latvec[2][0] =  0.5;
  latvec[2][1] =  0.5;

  atpos = memory->create(atpos,nucell, 3, "NaCl_prim_atpos");
  attyp = new int[nucell]; layer = new int[nucell];
  
  attyp[0] = 1;
  attyp[1] = 2;

  layer[0]    = 0;
  atpos[0][0] = 0.;
  atpos[0][1] = 0.;
  atpos[0][2] = 0.;

  layer[1]    = 1;
  atpos[1][0] = 0.5;
  atpos[1][1] = 0.5;
  atpos[1][2] = 0.5;

  initialized = 1;
return;
}

/* ----------------------------------------------------------------------
   Initialize NaCl Conventional cell
------------------------------------------------------------------------- */
void FCC::NaClConv()
{
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = new char[10];
  strcpy(name, "NaCl-conv");

  nucell = 8;
  ntype  = 2;
  
  latvec[0][0] = 1.;
  latvec[1][1] = 1.;
  latvec[2][2] = 1.;

  atpos = memory->create(atpos,nucell, 3, "NaCl_conv_atpos");
  attyp = new int[nucell]; layer = new int[nucell];
  
  for (int i=0; i<nucell; i+=2) attyp[i] = 1;
  for (int i=1; i<nucell; i+=2) attyp[i] = 2;

  layer[0]    = 0;
  atpos[0][0] = 0.;
  atpos[0][1] = 0.;
  atpos[0][2] = 0.;

  layer[1]    = 0;
  atpos[1][0] = 0.5;
  atpos[1][1] = 0.;
  atpos[1][2] = 0.;

  layer[2]    = 0;
  atpos[2][0] = 0.5;
  atpos[2][1] = 0.5;
  atpos[2][2] = 0.;

  layer[3]    = 0;
  atpos[3][0] = 0.;
  atpos[3][1] = 0.5;
  atpos[3][2] = 0.;

  layer[4]    = 1;
  atpos[4][0] = 0.5;
  atpos[4][1] = 0.0;
  atpos[4][2] = 0.5;

  layer[5]    = 1;
  atpos[5][0] = 0.0;
  atpos[5][1] = 0.0;
  atpos[5][2] = 0.5;

  layer[6]    = 1;
  atpos[6][0] = 0.0;
  atpos[6][1] = 0.5;
  atpos[6][2] = 0.5;

  layer[7]    = 1;
  atpos[7][0] = 0.5;
  atpos[7][1] = 0.5;
  atpos[7][2] = 0.5;

  initialized = 1;
return;
}
