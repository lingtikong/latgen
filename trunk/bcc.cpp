#include "bcc.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#define MAX_LINE_LENGTH 256

using namespace std;

/* ----------------------------------------------------------------------
   To select the orientation of the lattice
------------------------------------------------------------------------- */
BCC::BCC() : lattice()
{
  char str[MAX_LINE_LENGTH];
  // print out the menu
  alat = 1.;
  printf("\n===========================================================\n");
  printf("Please input the lattice constant of the BCC lattice [1.]: ");
  if (count_words(gets(str)) > 0) sscanf(str,"%lg",&alat);

  int orient = 1;
  printf("Please selection the orientation of the BCC lattice:\n");
  printf("   1. (001);\n");
  printf("   2. (110);\n");
  printf("   3. (111);\n");
  printf("   4. primitive cell;\n");
  printf("Your  choice [1]: ");
  if (count_words(gets(str)) > 0) sscanf(str,"%d",&orient);
  printf("You selected: %d", orient);
  printf("\n===========================================================\n");
  
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
    Primitive();
    break;
  default:
    break;
  }

}

/* ----------------------------------------------------------------------
   Deconstructor does nothing
------------------------------------------------------------------------- */
BCC::~BCC()
{

}

/* ----------------------------------------------------------------------
   Initialize for (001) orientation
------------------------------------------------------------------------- */
void BCC::BCC001()
{
  char str[MAX_LINE_LENGTH];
  int surftype = 1;
  // print out the menu
  printf("\n===========================================================\n");
  printf("Please selection the type of BCC(001) surface:\n");
  printf("   1. conventional orientation;\n");
  printf("   2. B2 structure;\n");
  printf("Your  choice [1]: ");
  if (count_words(gets(str)) > 0) sscanf(str,"%d",&surftype);
  printf("You selected: %d", surftype);
  printf("\n===========================================================\n");
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = new char[9];
  strcpy(name, "BCC(001)");

  // initialize according to surface type
  switch (surftype){
  case 1:
    nucell = 2;
    ntype  = 1;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1.;
    latvec[2][2] = 1.;

    atpos = memory->create(atpos,nucell, 3, "BCC001_atpos");
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
    nucell = 2;
    ntype  = 2;
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1.;
    latvec[2][2] = 1.;

    atpos = memory->create(atpos,nucell, 3, "BCC001_atpos");
    attyp = new int[nucell];
    layer = new int[nucell];
    
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
    break;
  default:
    break;
  }
return;
}

/* ----------------------------------------------------------------------
   Initialize for (110) orientation
------------------------------------------------------------------------- */
void BCC::BCC110()
{
  char str[MAX_LINE_LENGTH];
  int surftype=1;
  // print out the menu
  printf("\n===========================================================\n");
  printf("Please selection the type of BCC(110) surface:\n");
  printf("   1. orthogonal, long side along x\n");
  printf("   2. orthogonal, long side along y\n");
  printf("Your  choice [1]: ");
  if (count_words(gets(str)) > 0) sscanf(str,"%d",&surftype);
  printf("You selected: %d", surftype);
  printf("\n===========================================================\n");
  
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = new char[9];
  strcpy(name, "BCC(110)");

  // initialize according to surface type
  switch (surftype){
  case 1:
    nucell = 4;
    ntype  = 1;
    
    latvec[0][0] = sqrt(2.);
    latvec[1][1] = 1.;
    latvec[2][2] = sqrt(2.);

    atpos = memory->create(atpos,nucell, 3, "BCC110_atpos");
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
    atpos[2][1] = 0.0;
    atpos[2][2] = 0.5;

    layer[3]    = 1;
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

    atpos = memory->create(atpos,nucell, 3, "BCC110_atpos");
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
    atpos[3][1] = 0.0;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* ----------------------------------------------------------------------
   Initialize for (111) orientation // this part is not corrected yet
------------------------------------------------------------------------- */
void BCC::BCC111()
{
  char str[MAX_LINE_LENGTH];
  int surftype =1;
  // print out the menu
  printf("\n===========================================================\n");
  printf("Please selection the type of BCC(111) surface:\n");
  printf("   1. orthogonal, long side along x\n");
  printf("   2. orthogonal, long side along y\n");
  printf("Your  choice [1]: ");
  if (count_words(gets(str)) > 0) sscanf(str,"%d",&surftype);
  printf("You selected: %d", surftype);
  printf("\n===========================================================\n");
  
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++) latvec[i][j] = 0.;
  }
  name = new char[9];
  strcpy(name, "BCC(111)");

  // initialize according to surface type
  switch (surftype){
  case 1:
    nucell = 4;
    ntype  = 1;
    
    latvec[0][0] = sqrt(2.);
    latvec[1][1] = 1.;
    latvec[2][2] = sqrt(2.);

    atpos = memory->create(atpos,nucell, 3, "BCC111_atpos");
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
    atpos[2][1] = 0.0;
    atpos[2][2] = 0.5;

    layer[3]    = 1;
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

    atpos = memory->create(atpos,nucell, 3, "BCC111_atpos");
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
    atpos[3][1] = 0.0;
    atpos[3][2] = 0.5;

    initialized = 1;
    break;
  default:
    break;
  }
return;
}

/* ----------------------------------------------------------------------
   Initialize for Primitive cell
------------------------------------------------------------------------- */
void BCC::Primitive()
{
  name = new char[9];
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

  atpos = memory->create(atpos,nucell, 3, "BCC001_atpos");
  attyp = new int[nucell];
  layer = new int[nucell];
  
  for (int i=0; i<nucell; i++) attyp[i] = 1;
  layer[0]    = 0;
  atpos[0][0] = 0.;
  atpos[0][1] = 0.;
  atpos[0][2] = 0.;

  initialized = 1;
return;
}
