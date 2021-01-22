#include "bcc.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "common.h"

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
  if (alat <= 0.) alat = 1.;

  int orient = 1;
  printf("Please select the orientation/type of the BCC lattice:\n");
  printf("   1. [001] along z;\n");
  printf("   2. [110] along z;\n");
  printf("   3. [111] along z;\n");
  printf("   4. [112] along z;\n");
  printf("   5. primitive cell;\n");
  printf("Your choice [%d]: ", orient);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) orient = inumeric(strtok(str, " \t\n\r\f"));
  printf("Your selection : %d", orient);
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
  printf("Please select the type of BCC(001) cell:\n");
  printf("   1. conventional, [100] along x;\n");
  printf("   2. conventional, B2 structure;\n");
  printf("   3. supercell, [110] along x;\n");
  printf("   4. B2 supercell, [110] along x;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = inumeric(strtok(str, " \t\n\r\f"));
  printf("Your selection : %d", surftype);
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
    noct   = nucell * 3;
    ntetra = nucell * 6;

    memory->create(atpos, nucell + noct + ntetra, 3, "BCC001:atpos");
    memory->create(attyp, nucell + noct + ntetra, "BCC001:attyp");
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1.;
    latvec[2][2] = 1.;

    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.5;

    for (int i = nucell; i < nucell + noct + ntetra; ++i) attyp[i] = 2;
    // octohedra site
    atpos[2][0] = 0.500;
    atpos[2][1] = 0.000;
    atpos[2][2] = 0.000;

    atpos[3][0] = 0.000;
    atpos[3][1] = 0.500;
    atpos[3][2] = 0.000;

    atpos[4][0] = 0.500;
    atpos[4][1] = 0.500;
    atpos[4][2] = 0.000;

    atpos[5][0] = 0.000;
    atpos[5][1] = 0.000;
    atpos[5][2] = 0.500;

    atpos[6][0] = 0.000;
    atpos[6][1] = 0.500;
    atpos[6][2] = 0.500;

    atpos[7][0] = 0.500;
    atpos[7][1] = 0.000;
    atpos[7][2] = 0.500;

    // tetrahedral sites
    atpos[8][0] = 0.250;
    atpos[8][1] = 0.500;
    atpos[8][2] = 0.000;

    atpos[9][0] = 0.750;
    atpos[9][1] = 0.500;
    atpos[9][2] = 0.000;

    atpos[10][0] = 0.500;
    atpos[10][1] = 0.250;
    atpos[10][2] = 0.000;

    atpos[11][0] = 0.500;
    atpos[11][1] = 0.750;
    atpos[11][2] = 0.000;

    atpos[12][0] = 0.500;
    atpos[12][1] = 0.000;
    atpos[12][2] = 0.250;

    atpos[13][0] = 0.000;
    atpos[13][1] = 0.500;
    atpos[13][2] = 0.250;

    atpos[14][0] = 0.000;
    atpos[14][1] = 0.250;
    atpos[14][2] = 0.500;

    atpos[15][0] = 0.000;
    atpos[15][1] = 0.750;
    atpos[15][2] = 0.500;

    atpos[16][0] = 0.250;
    atpos[16][1] = 0.000;
    atpos[16][2] = 0.500;

    atpos[17][0] = 0.750;
    atpos[17][1] = 0.000;
    atpos[17][2] = 0.500;

    atpos[18][0] = 0.000;
    atpos[18][1] = 0.500;
    atpos[18][2] = 0.750;

    atpos[19][0] = 0.500;
    atpos[19][1] = 0.000;
    atpos[19][2] = 0.750;

    initialized = 1;
    break;

  case 2:
    nucell = 2;
    ntype  = 2;
    noct   = nucell * 3;
    ntetra = nucell * 6;

    memory->create(atpos, nucell + noct + ntetra, 3, "BCC001:atpos");
    memory->create(attyp, nucell + noct + ntetra, "BCC001:attyp");
    
    latvec[0][0] = 1.;
    latvec[1][1] = 1.;
    latvec[2][2] = 1.;

    attyp[0] = 1;
    attyp[1] = 2;

    atpos[0][0] = 0.;
    atpos[0][1] = 0.;
    atpos[0][2] = 0.;

    atpos[1][0] = 0.5;
    atpos[1][1] = 0.5;
    atpos[1][2] = 0.5;

    for (int i = nucell; i < nucell + noct + ntetra; ++i) attyp[i] = 3;
    // octohedra site
    atpos[2][0] = 0.500;
    atpos[2][1] = 0.000;
    atpos[2][2] = 0.000;

    atpos[3][0] = 0.000;
    atpos[3][1] = 0.500;
    atpos[3][2] = 0.000;

    atpos[4][0] = 0.500;
    atpos[4][1] = 0.500;
    atpos[4][2] = 0.000;

    atpos[5][0] = 0.000;
    atpos[5][1] = 0.000;
    atpos[5][2] = 0.500;

    atpos[6][0] = 0.000;
    atpos[6][1] = 0.500;
    atpos[6][2] = 0.500;

    atpos[7][0] = 0.500;
    atpos[7][1] = 0.000;
    atpos[7][2] = 0.500;

    // tetrahedral sites
    atpos[8][0] = 0.250;
    atpos[8][1] = 0.500;
    atpos[8][2] = 0.000;

    atpos[9][0] = 0.750;
    atpos[9][1] = 0.500;
    atpos[9][2] = 0.000;

    atpos[10][0] = 0.500;
    atpos[10][1] = 0.250;
    atpos[10][2] = 0.000;

    atpos[11][0] = 0.500;
    atpos[11][1] = 0.750;
    atpos[11][2] = 0.000;

    atpos[12][0] = 0.500;
    atpos[12][1] = 0.000;
    atpos[12][2] = 0.250;

    atpos[13][0] = 0.000;
    atpos[13][1] = 0.500;
    atpos[13][2] = 0.250;

    atpos[14][0] = 0.000;
    atpos[14][1] = 0.250;
    atpos[14][2] = 0.500;

    atpos[15][0] = 0.000;
    atpos[15][1] = 0.750;
    atpos[15][2] = 0.500;

    atpos[16][0] = 0.250;
    atpos[16][1] = 0.000;
    atpos[16][2] = 0.500;

    atpos[17][0] = 0.750;
    atpos[17][1] = 0.000;
    atpos[17][2] = 0.500;

    atpos[18][0] = 0.000;
    atpos[18][1] = 0.500;
    atpos[18][2] = 0.750;

    atpos[19][0] = 0.500;
    atpos[19][1] = 0.000;
    atpos[19][2] = 0.750;

    initialized = 1;
    break;

  case 3:
    nucell = 4;
    ntype  = 1;
    noct   = nucell * 3;
    ntetra = nucell * 6;

    memory->create(atpos, nucell + noct + ntetra, 3, "BCC001:atpos");
    memory->create(attyp, nucell + noct + ntetra, "BCC001:attyp");
    
    latvec[0][0] = sqrt(2.);
    latvec[1][1] = sqrt(2.);
    latvec[2][2] = 1.;

    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    atpos[0][0] = 0.000;
    atpos[0][1] = 0.000;
    atpos[0][2] = 0.000;

    atpos[1][0] = 0.500;
    atpos[1][1] = 0.500;
    atpos[1][2] = 0.000;

    atpos[2][0] = 0.000;
    atpos[2][1] = 0.500;
    atpos[2][2] = 0.500;

    atpos[3][0] = 0.500;
    atpos[3][1] = 0.000;
    atpos[3][2] = 0.500;

    for (int i = nucell; i < nucell + noct + ntetra; ++i) attyp[i] = 2;
    // octohedra site
    atpos[4][0] = 0.250;
    atpos[4][1] = 0.750;
    atpos[4][2] = 0.000;

    atpos[5][0] = 0.750;
    atpos[5][1] = 0.250;
    atpos[5][2] = 0.000;

    atpos[6][0] = 0.250;
    atpos[6][1] = 0.250;
    atpos[6][2] = 0.000;

    atpos[7][0] = 0.750;
    atpos[7][1] = 0.750;
    atpos[7][2] = 0.000;

    atpos[8][0] = 0.000;
    atpos[8][1] = 0.500;
    atpos[8][2] = 0.000;

    atpos[9][0] = 0.500;
    atpos[9][1] = 0.000;
    atpos[9][2] = 0.000;

    atpos[10][0] = 0.000;
    atpos[10][1] = 0.000;
    atpos[10][2] = 0.500;

    atpos[11][0] = 0.500;
    atpos[11][1] = 0.500;
    atpos[11][2] = 0.500;

    atpos[12][0] = 0.250;
    atpos[12][1] = 0.250;
    atpos[12][2] = 0.500;

    atpos[13][0] = 0.750;
    atpos[13][1] = 0.750;
    atpos[13][2] = 0.500;

    atpos[14][0] = 0.250;
    atpos[14][1] = 0.750;
    atpos[14][2] = 0.500;

    atpos[15][0] = 0.750;
    atpos[15][1] = 0.250;
    atpos[15][2] = 0.500;

    // tetrahedral site
    atpos[16][0] = 0.375;
    atpos[16][1] = 0.125;
    atpos[16][2] = 0.000;

    atpos[17][0] = 0.875;
    atpos[17][1] = 0.625;
    atpos[17][2] = 0.000;

    atpos[18][0] = 0.125;
    atpos[18][1] = 0.375;
    atpos[18][2] = 0.000;

    atpos[19][0] = 0.625;
    atpos[19][1] = 0.875;
    atpos[19][2] = 0.000;

    atpos[20][0] = 0.375;
    atpos[20][1] = 0.875;
    atpos[20][2] = 0.000;

    atpos[21][0] = 0.875;
    atpos[21][1] = 0.375;
    atpos[21][2] = 0.000;

    atpos[22][0] = 0.125;
    atpos[22][1] = 0.625;
    atpos[22][2] = 0.000;

    atpos[23][0] = 0.625;
    atpos[23][1] = 0.125;
    atpos[23][2] = 0.000;

    atpos[24][0] = 0.250;
    atpos[24][1] = 0.750;
    atpos[24][2] = 0.250;

    atpos[25][0] = 0.750;
    atpos[25][1] = 0.250;
    atpos[25][2] = 0.250;

    atpos[26][0] = 0.250;
    atpos[26][1] = 0.250;
    atpos[26][2] = 0.250;

    atpos[27][0] = 0.750;
    atpos[27][1] = 0.750;
    atpos[27][2] = 0.250;

    atpos[28][0] = 0.125;
    atpos[28][1] = 0.125;
    atpos[28][2] = 0.500;

    atpos[29][0] = 0.625;
    atpos[29][1] = 0.625;
    atpos[29][2] = 0.500;

    atpos[30][0] = 0.375;
    atpos[30][1] = 0.375;
    atpos[30][2] = 0.500;

    atpos[31][0] = 0.875;
    atpos[31][1] = 0.875;
    atpos[31][2] = 0.500;

    atpos[32][0] = 0.125;
    atpos[32][1] = 0.875;
    atpos[32][2] = 0.500;

    atpos[33][0] = 0.625;
    atpos[33][1] = 0.375;
    atpos[33][2] = 0.500;

    atpos[34][0] = 0.375;
    atpos[34][1] = 0.625;
    atpos[34][2] = 0.500;

    atpos[35][0] = 0.875;
    atpos[35][1] = 0.125;
    atpos[35][2] = 0.500;

    atpos[36][0] = 0.250;
    atpos[36][1] = 0.250;
    atpos[36][2] = 0.750;

    atpos[37][0] = 0.750;
    atpos[37][1] = 0.750;
    atpos[37][2] = 0.750;

    atpos[38][0] = 0.250;
    atpos[38][1] = 0.750;
    atpos[38][2] = 0.750;

    atpos[39][0] = 0.750;
    atpos[39][1] = 0.250;
    atpos[39][2] = 0.750;

    initialized = 1;
    break;

  case 4:
    nucell = 4;
    ntype  = 2;
    noct   = nucell * 3;
    ntetra = nucell * 6;

    memory->create(atpos, nucell + noct + ntetra, 3, "BCC001:atpos");
    memory->create(attyp, nucell + noct + ntetra, "BCC001:attyp");
    
    latvec[0][0] = sqrt(2.);
    latvec[1][1] = sqrt(2.);
    latvec[2][2] = 1.;

    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    for (int i = 2; i < nucell; ++i) attyp[i] = 2;

    atpos[0][0] = 0.000;
    atpos[0][1] = 0.000;
    atpos[0][2] = 0.000;

    atpos[1][0] = 0.500;
    atpos[1][1] = 0.500;
    atpos[1][2] = 0.000;

    atpos[2][0] = 0.000;
    atpos[2][1] = 0.500;
    atpos[2][2] = 0.500;

    atpos[3][0] = 0.500;
    atpos[3][1] = 0.000;
    atpos[3][2] = 0.500;

    for (int i = nucell; i < nucell + noct + ntetra; ++i) attyp[i] = 3;
    // octohedra site
    atpos[4][0] = 0.250;
    atpos[4][1] = 0.750;
    atpos[4][2] = 0.000;

    atpos[5][0] = 0.750;
    atpos[5][1] = 0.250;
    atpos[5][2] = 0.000;

    atpos[6][0] = 0.250;
    atpos[6][1] = 0.250;
    atpos[6][2] = 0.000;

    atpos[7][0] = 0.750;
    atpos[7][1] = 0.750;
    atpos[7][2] = 0.000;

    atpos[8][0] = 0.000;
    atpos[8][1] = 0.500;
    atpos[8][2] = 0.000;

    atpos[9][0] = 0.500;
    atpos[9][1] = 0.000;
    atpos[9][2] = 0.000;

    atpos[10][0] = 0.000;
    atpos[10][1] = 0.000;
    atpos[10][2] = 0.500;

    atpos[11][0] = 0.500;
    atpos[11][1] = 0.500;
    atpos[11][2] = 0.500;

    atpos[12][0] = 0.250;
    atpos[12][1] = 0.250;
    atpos[12][2] = 0.500;

    atpos[13][0] = 0.750;
    atpos[13][1] = 0.750;
    atpos[13][2] = 0.500;

    atpos[14][0] = 0.250;
    atpos[14][1] = 0.750;
    atpos[14][2] = 0.500;

    atpos[15][0] = 0.750;
    atpos[15][1] = 0.250;
    atpos[15][2] = 0.500;

    // tetrahedral site
    atpos[16][0] = 0.375;
    atpos[16][1] = 0.125;
    atpos[16][2] = 0.000;

    atpos[17][0] = 0.875;
    atpos[17][1] = 0.625;
    atpos[17][2] = 0.000;

    atpos[18][0] = 0.125;
    atpos[18][1] = 0.375;
    atpos[18][2] = 0.000;

    atpos[19][0] = 0.625;
    atpos[19][1] = 0.875;
    atpos[19][2] = 0.000;

    atpos[20][0] = 0.375;
    atpos[20][1] = 0.875;
    atpos[20][2] = 0.000;

    atpos[21][0] = 0.875;
    atpos[21][1] = 0.375;
    atpos[21][2] = 0.000;

    atpos[22][0] = 0.125;
    atpos[22][1] = 0.625;
    atpos[22][2] = 0.000;

    atpos[23][0] = 0.625;
    atpos[23][1] = 0.125;
    atpos[23][2] = 0.000;

    atpos[24][0] = 0.250;
    atpos[24][1] = 0.750;
    atpos[24][2] = 0.250;

    atpos[25][0] = 0.750;
    atpos[25][1] = 0.250;
    atpos[25][2] = 0.250;

    atpos[26][0] = 0.250;
    atpos[26][1] = 0.250;
    atpos[26][2] = 0.250;

    atpos[27][0] = 0.750;
    atpos[27][1] = 0.750;
    atpos[27][2] = 0.250;

    atpos[28][0] = 0.125;
    atpos[28][1] = 0.125;
    atpos[28][2] = 0.500;

    atpos[29][0] = 0.625;
    atpos[29][1] = 0.625;
    atpos[29][2] = 0.500;

    atpos[30][0] = 0.375;
    atpos[30][1] = 0.375;
    atpos[30][2] = 0.500;

    atpos[31][0] = 0.875;
    atpos[31][1] = 0.875;
    atpos[31][2] = 0.500;

    atpos[32][0] = 0.125;
    atpos[32][1] = 0.875;
    atpos[32][2] = 0.500;

    atpos[33][0] = 0.625;
    atpos[33][1] = 0.375;
    atpos[33][2] = 0.500;

    atpos[34][0] = 0.375;
    atpos[34][1] = 0.625;
    atpos[34][2] = 0.500;

    atpos[35][0] = 0.875;
    atpos[35][1] = 0.125;
    atpos[35][2] = 0.500;

    atpos[36][0] = 0.250;
    atpos[36][1] = 0.250;
    atpos[36][2] = 0.750;

    atpos[37][0] = 0.750;
    atpos[37][1] = 0.750;
    atpos[37][2] = 0.750;

    atpos[38][0] = 0.250;
    atpos[38][1] = 0.750;
    atpos[38][2] = 0.750;

    atpos[39][0] = 0.750;
    atpos[39][1] = 0.250;
    atpos[39][2] = 0.750;

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
  printf("Please select the type of BCC(110) cell:\n");
  printf("   1. orthogonal, [110] along x\n");
  printf("   2. orthogonal, [110] along y\n");
  printf("   3. non-orthogonal, [1-11] along x\n");
  printf("   4. non-orthogonal, [1-1-1] along y\n");
  printf("   5. non-orthogonal, [-111] along x\n");
  printf("   6. non-orthogonal, [1-11] along y\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = inumeric(strtok(str, " \t\n\r\f"));
  printf("Your selection : %d", surftype);
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
    noct   = nucell * 3;
    ntetra = nucell * 6;

    memory->create(atpos, nucell + noct + ntetra, 3, "BCC110:atpos");
    memory->create(attyp, nucell + noct + ntetra, "BCC110:attyp");
    
    latvec[0][0] = sqrt(2.);
    latvec[1][1] = 1.;
    latvec[2][2] = sqrt(2.);

    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    for (int i = nucell; i < nucell + noct + ntetra; ++i) attyp[i] = 2;

    atpos[0][0] = 0.000;
    atpos[0][1] = 0.000;
    atpos[0][2] = 0.000;

    atpos[1][0] = 0.500;
    atpos[1][1] = 0.500;
    atpos[1][2] = 0.000;

    atpos[2][0] = 0.500;
    atpos[2][1] = 0.000;
    atpos[2][2] = 0.500;

    atpos[3][0] = 0.000;
    atpos[3][1] = 0.500;
    atpos[3][2] = 0.500;

    atpos[4][0] = 0.000;
    atpos[4][1] = 0.500;
    atpos[4][2] = 0.000;

    atpos[5][0] = 0.500;
    atpos[5][1] = 0.000;
    atpos[5][2] = 0.000;

    atpos[6][0] = 0.250;
    atpos[6][1] = 0.000;
    atpos[6][2] = 0.250;

    atpos[7][0] = 0.750;
    atpos[7][1] = 0.500;
    atpos[7][2] = 0.250;

    atpos[8][0] = 0.750;
    atpos[8][1] = 0.000;
    atpos[8][2] = 0.250;

    atpos[9][0] = 0.250;
    atpos[9][1] = 0.500;
    atpos[9][2] = 0.250;

    atpos[10][0] = 0.500;
    atpos[10][1] = 0.500;
    atpos[10][2] = 0.500;

    atpos[11][0] = 0.000;
    atpos[11][1] = 0.000;
    atpos[11][2] = 0.500;

    atpos[12][0] = 0.750;
    atpos[12][1] = 0.500;
    atpos[12][2] = 0.750;

    atpos[13][0] = 0.250;
    atpos[13][1] = 0.500;
    atpos[13][2] = 0.750;

    atpos[14][0] = 0.250;
    atpos[14][1] = 0.000;
    atpos[14][2] = 0.750;

    atpos[15][0] = 0.750;
    atpos[15][1] = 0.000;
    atpos[15][2] = 0.750;

    atpos[16][0] = 0.375;
    atpos[16][1] = 0.000;
    atpos[16][2] = 0.125;

    atpos[17][0] = 0.125;
    atpos[17][1] = 0.500;
    atpos[17][2] = 0.125;

    atpos[18][0] = 0.625;
    atpos[18][1] = 0.000;
    atpos[18][2] = 0.125;

    atpos[19][0] = 0.875;
    atpos[19][1] = 0.500;
    atpos[19][2] = 0.125;

    atpos[20][0] = 0.750;
    atpos[20][1] = 0.250;
    atpos[20][2] = 0.250;

    atpos[21][0] = 0.750;
    atpos[21][1] = 0.750;
    atpos[21][2] = 0.250;

    atpos[22][0] = 0.250;
    atpos[22][1] = 0.750;
    atpos[22][2] = 0.250;

    atpos[23][0] = 0.250;
    atpos[23][1] = 0.250;
    atpos[23][2] = 0.250;

    atpos[24][0] = 0.125;
    atpos[24][1] = 0.000;
    atpos[24][2] = 0.375;

    atpos[25][0] = 0.625;
    atpos[25][1] = 0.500;
    atpos[25][2] = 0.375;

    atpos[26][0] = 0.375;
    atpos[26][1] = 0.500;
    atpos[26][2] = 0.375;

    atpos[27][0] = 0.875;
    atpos[27][1] = 0.000;
    atpos[27][2] = 0.375;

    atpos[28][0] = 0.375;
    atpos[28][1] = 0.500;
    atpos[28][2] = 0.625;

    atpos[29][0] = 0.875;
    atpos[29][1] = 0.000;
    atpos[29][2] = 0.625;

    atpos[30][0] = 0.125;
    atpos[30][1] = 0.000;
    atpos[30][2] = 0.625;

    atpos[31][0] = 0.625;
    atpos[31][1] = 0.500;
    atpos[31][2] = 0.625;

    atpos[32][0] = 0.250;
    atpos[32][1] = 0.250;
    atpos[32][2] = 0.750;

    atpos[33][0] = 0.250;
    atpos[33][1] = 0.750;
    atpos[33][2] = 0.750;

    atpos[34][0] = 0.750;
    atpos[34][1] = 0.750;
    atpos[34][2] = 0.750;

    atpos[35][0] = 0.750;
    atpos[35][1] = 0.250;
    atpos[35][2] = 0.750;

    atpos[36][0] = 0.875;
    atpos[36][1] = 0.500;
    atpos[36][2] = 0.875;

    atpos[37][0] = 0.625;
    atpos[37][1] = 0.000;
    atpos[37][2] = 0.875;

    atpos[38][0] = 0.375;
    atpos[38][1] = 0.000;
    atpos[38][2] = 0.875;

    atpos[39][0] = 0.125;
    atpos[39][1] = 0.500;
    atpos[39][2] = 0.875;

    initialized = 1;
    break;

  case 2:
    nucell = 4;
    ntype  = 1;
    noct   = nucell * 3;
    ntetra = nucell * 6;

    memory->create(atpos, nucell + noct + ntetra, 3, "BCC110:atpos");
    memory->create(attyp, nucell + noct + ntetra, "BCC110:attyp");
    
    latvec[0][0] = 1.;
    latvec[1][1] = sqrt(2.);
    latvec[2][2] = sqrt(2.);

    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    for (int i = nucell; i < nucell + noct + ntetra; ++i) attyp[i] = 2;

    atpos[0][0] = 0.000;
    atpos[0][1] = 0.000;
    atpos[0][2] = 0.000;

    atpos[1][0] = 0.500;
    atpos[1][1] = 0.500;
    atpos[1][2] = 0.000;

    atpos[2][0] = 0.000;
    atpos[2][1] = 0.500;
    atpos[2][2] = 0.500;

    atpos[3][0] = 0.500;
    atpos[3][1] = 0.000;
    atpos[3][2] = 0.500;

    atpos[4][0] = 0.500;
    atpos[4][1] = 0.000;
    atpos[4][2] = 0.000;

    atpos[5][0] = 0.000;
    atpos[5][1] = 0.500;
    atpos[5][2] = 0.000;

    atpos[6][0] = 0.000;
    atpos[6][1] = 0.750;
    atpos[6][2] = 0.250;

    atpos[7][0] = 0.500;
    atpos[7][1] = 0.250;
    atpos[7][2] = 0.250;

    atpos[8][0] = 0.000;
    atpos[8][1] = 0.250;
    atpos[8][2] = 0.250;

    atpos[9][0] = 0.500;
    atpos[9][1] = 0.750;
    atpos[9][2] = 0.250;

    atpos[10][0] = 0.500;
    atpos[10][1] = 0.500;
    atpos[10][2] = 0.500;

    atpos[11][0] = 0.000;
    atpos[11][1] = 0.000;
    atpos[11][2] = 0.500;

    atpos[12][0] = 0.500;
    atpos[12][1] = 0.250;
    atpos[12][2] = 0.750;

    atpos[13][0] = 0.500;
    atpos[13][1] = 0.750;
    atpos[13][2] = 0.750;

    atpos[14][0] = 0.000;
    atpos[14][1] = 0.750;
    atpos[14][2] = 0.750;

    atpos[15][0] = 0.000;
    atpos[15][1] = 0.250;
    atpos[15][2] = 0.750;

    atpos[16][0] = 0.000;
    atpos[16][1] = 0.625;
    atpos[16][2] = 0.125;

    atpos[17][0] = 0.500;
    atpos[17][1] = 0.875;
    atpos[17][2] = 0.125;

    atpos[18][0] = 0.000;
    atpos[18][1] = 0.375;
    atpos[18][2] = 0.125;

    atpos[19][0] = 0.500;
    atpos[19][1] = 0.125;
    atpos[19][2] = 0.125;

    atpos[20][0] = 0.250;
    atpos[20][1] = 0.250;
    atpos[20][2] = 0.250;

    atpos[21][0] = 0.750;
    atpos[21][1] = 0.250;
    atpos[21][2] = 0.250;

    atpos[22][0] = 0.750;
    atpos[22][1] = 0.750;
    atpos[22][2] = 0.250;

    atpos[23][0] = 0.250;
    atpos[23][1] = 0.750;
    atpos[23][2] = 0.250;

    atpos[24][0] = 0.000;
    atpos[24][1] = 0.875;
    atpos[24][2] = 0.375;

    atpos[25][0] = 0.500;
    atpos[25][1] = 0.375;
    atpos[25][2] = 0.375;

    atpos[26][0] = 0.500;
    atpos[26][1] = 0.625;
    atpos[26][2] = 0.375;

    atpos[27][0] = 0.000;
    atpos[27][1] = 0.125;
    atpos[27][2] = 0.375;

    atpos[28][0] = 0.500;
    atpos[28][1] = 0.625;
    atpos[28][2] = 0.625;

    atpos[29][0] = 0.000;
    atpos[29][1] = 0.125;
    atpos[29][2] = 0.625;

    atpos[30][0] = 0.000;
    atpos[30][1] = 0.875;
    atpos[30][2] = 0.625;

    atpos[31][0] = 0.500;
    atpos[31][1] = 0.375;
    atpos[31][2] = 0.625;

    atpos[32][0] = 0.250;
    atpos[32][1] = 0.750;
    atpos[32][2] = 0.750;

    atpos[33][0] = 0.750;
    atpos[33][1] = 0.750;
    atpos[33][2] = 0.750;

    atpos[34][0] = 0.750;
    atpos[34][1] = 0.250;
    atpos[34][2] = 0.750;

    atpos[35][0] = 0.250;
    atpos[35][1] = 0.250;
    atpos[35][2] = 0.750;

    atpos[36][0] = 0.500;
    atpos[36][1] = 0.125;
    atpos[36][2] = 0.875;

    atpos[37][0] = 0.000;
    atpos[37][1] = 0.375;
    atpos[37][2] = 0.875;

    atpos[38][0] = 0.000;
    atpos[38][1] = 0.625;
    atpos[38][2] = 0.875;

    atpos[39][0] = 0.500;
    atpos[39][1] = 0.875;
    atpos[39][2] = 0.875;

    initialized = 1;
    break;

  case 3:
    nucell = 2;
    ntype  = 1;
    noct   = nucell * 3;
    ntetra = nucell * 6;

    memory->create(atpos, nucell + noct + ntetra, 3, "BCC110:atpos");
    memory->create(attyp, nucell + noct + ntetra, "BCC110:attyp");
    
    latvec[0][0] = sqrt(3.)/2.;
    latvec[1][0] = sqrt(3.)/6.;
    latvec[1][1] = sqrt(6.)/3.;
    latvec[2][2] = sqrt(2.);

    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    for (int i = nucell; i < nucell + noct + ntetra; ++i) attyp[i] = 2;

    atpos[0][0] = 0.000;
    atpos[0][1] = 0.000;
    atpos[0][2] = 0.000;

    atpos[1][0] = 0.500;
    atpos[1][1] = 0.500;
    atpos[1][2] = 0.500;

    atpos[2][0] = 0.500;
    atpos[2][1] = 0.500;
    atpos[2][2] = 0.000;

    atpos[3][0] = 0.750;
    atpos[3][1] = 0.750;
    atpos[3][2] = 0.250;

    atpos[4][0] = 0.250;
    atpos[4][1] = 0.250;
    atpos[4][2] = 0.250;

    atpos[5][0] = 0.000;
    atpos[5][1] = 0.000;
    atpos[5][2] = 0.500;

    atpos[6][0] = 0.250;
    atpos[6][1] = 0.250;
    atpos[6][2] = 0.750;

    atpos[7][0] = 0.750;
    atpos[7][1] = 0.750;
    atpos[7][2] = 0.750;

    atpos[8][0] = 0.625;
    atpos[8][1] = 0.625;
    atpos[8][2] = 0.125;

    atpos[9][0] = 0.375;
    atpos[9][1] = 0.375;
    atpos[9][2] = 0.125;

    atpos[10][0] = 0.500;
    atpos[10][1] = 0.000;
    atpos[10][2] = 0.250;

    atpos[11][0] = 0.000;
    atpos[11][1] = 0.500;
    atpos[11][2] = 0.250;

    atpos[12][0] = 0.125;
    atpos[12][1] = 0.125;
    atpos[12][2] = 0.375;

    atpos[13][0] = 0.875;
    atpos[13][1] = 0.875;
    atpos[13][2] = 0.375;

    atpos[14][0] = 0.125;
    atpos[14][1] = 0.125;
    atpos[14][2] = 0.625;

    atpos[15][0] = 0.875;
    atpos[15][1] = 0.875;
    atpos[15][2] = 0.625;

    atpos[16][0] = 0.000;
    atpos[16][1] = 0.500;
    atpos[16][2] = 0.750;

    atpos[17][0] = 0.500;
    atpos[17][1] = 0.000;
    atpos[17][2] = 0.750;

    atpos[18][0] = 0.375;
    atpos[18][1] = 0.375;
    atpos[18][2] = 0.875;

    atpos[19][0] = 0.625;
    atpos[19][1] = 0.625;
    atpos[19][2] = 0.875;

    initialized = 1;
    break;

  case 4:
    nucell = 2;
    ntype  = 1;
    noct   = nucell * 3;
    ntetra = nucell * 6;

    memory->create(atpos, nucell + noct + ntetra, 3, "BCC110:atpos");
    memory->create(attyp, nucell + noct + ntetra, "BCC110:attyp");
    
    latvec[0][0] = sqrt(6.)/3.;
    latvec[0][1] = sqrt(3.)/6.;
    latvec[1][1] = sqrt(3.)/2.;
    latvec[2][2] = sqrt(2.);

    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    for (int i = nucell; i < nucell + noct + ntetra; ++i) attyp[i] = 2;

    atpos[0][0] = 0.000;
    atpos[0][1] = 0.000;
    atpos[0][2] = 0.000;

    atpos[1][0] = 0.500;
    atpos[1][1] = 0.500;
    atpos[1][2] = 0.500;

    atpos[2][0] = 0.500;
    atpos[2][1] = 0.500;
    atpos[2][2] = 0.000;

    atpos[3][0] = 0.750;
    atpos[3][1] = 0.750;
    atpos[3][2] = 0.250;

    atpos[4][0] = 0.250;
    atpos[4][1] = 0.250;
    atpos[4][2] = 0.250;

    atpos[5][0] = 0.000;
    atpos[5][1] = 0.000;
    atpos[5][2] = 0.500;

    atpos[6][0] = 0.250;
    atpos[6][1] = 0.250;
    atpos[6][2] = 0.750;

    atpos[7][0] = 0.750;
    atpos[7][1] = 0.750;
    atpos[7][2] = 0.750;

    atpos[8][0] = 0.625;
    atpos[8][1] = 0.625;
    atpos[8][2] = 0.125;

    atpos[9][0] = 0.375;
    atpos[9][1] = 0.375;
    atpos[9][2] = 0.125;

    atpos[10][0] = 0.500;
    atpos[10][1] = 0.000;
    atpos[10][2] = 0.250;

    atpos[11][0] = 0.000;
    atpos[11][1] = 0.500;
    atpos[11][2] = 0.250;

    atpos[12][0] = 0.125;
    atpos[12][1] = 0.125;
    atpos[12][2] = 0.375;

    atpos[13][0] = 0.875;
    atpos[13][1] = 0.875;
    atpos[13][2] = 0.375;

    atpos[14][0] = 0.125;
    atpos[14][1] = 0.125;
    atpos[14][2] = 0.625;

    atpos[15][0] = 0.875;
    atpos[15][1] = 0.875;
    atpos[15][2] = 0.625;

    atpos[16][0] = 0.000;
    atpos[16][1] = 0.500;
    atpos[16][2] = 0.750;

    atpos[17][0] = 0.500;
    atpos[17][1] = 0.000;
    atpos[17][2] = 0.750;

    atpos[18][0] = 0.375;
    atpos[18][1] = 0.375;
    atpos[18][2] = 0.875;

    atpos[19][0] = 0.625;
    atpos[19][1] = 0.625;
    atpos[19][2] = 0.875;

    initialized = 1;
    break;

  case 5:
    nucell = 2;
    ntype  = 1;
    noct   = nucell * 3;
    ntetra = nucell * 6;

    memory->create(atpos, nucell + noct + ntetra, 3, "BCC110:atpos");
    memory->create(attyp, nucell + noct + ntetra, "BCC110:attyp");
    
    latvec[0][0] =  sqrt(3.)/2.;
    latvec[1][0] = -sqrt(3.)/6.;
    latvec[1][1] =  sqrt(6.)/3.;
    latvec[2][2] =  sqrt(2.);

    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    for (int i = nucell; i < nucell + noct + ntetra; ++i) attyp[i] = 2;
    atpos[0][0] = 0.000;
    atpos[0][1] = 0.000;
    atpos[0][2] = 0.000;

    atpos[1][0] = 0.500;
    atpos[1][1] = 0.500;
    atpos[1][2] = 0.500;

    atpos[2][0] = 0.500;
    atpos[2][1] = 0.500;
    atpos[2][2] = 0.000;

    atpos[3][0] = 0.250;
    atpos[3][1] = 0.750;
    atpos[3][2] = 0.250;

    atpos[4][0] = 0.750;
    atpos[4][1] = 0.250;
    atpos[4][2] = 0.250;

    atpos[5][0] = 1.000;
    atpos[5][1] = 0.000;
    atpos[5][2] = 0.500;

    atpos[6][0] = 0.750;
    atpos[6][1] = 0.250;
    atpos[6][2] = 0.750;

    atpos[7][0] = 0.250;
    atpos[7][1] = 0.750;
    atpos[7][2] = 0.750;

    atpos[8][0] = 0.625;
    atpos[8][1] = 0.375;
    atpos[8][2] = 0.125;

    atpos[9][0] = 0.375;
    atpos[9][1] = 0.625;
    atpos[9][2] = 0.125;

    atpos[10][0] = 0.000;
    atpos[10][1] = 0.500;
    atpos[10][2] = 0.250;

    atpos[11][0] = 0.500;
    atpos[11][1] = 0.000;
    atpos[11][2] = 0.250;

    atpos[12][0] = 0.875;
    atpos[12][1] = 0.125;
    atpos[12][2] = 0.375;

    atpos[13][0] = 0.125;
    atpos[13][1] = 0.875;
    atpos[13][2] = 0.375;

    atpos[14][0] = 0.875;
    atpos[14][1] = 0.125;
    atpos[14][2] = 0.625;

    atpos[15][0] = 0.125;
    atpos[15][1] = 0.875;
    atpos[15][2] = 0.625;

    atpos[16][0] = 0.000;
    atpos[16][1] = 0.500;
    atpos[16][2] = 0.750;

    atpos[17][0] = 0.500;
    atpos[17][1] = 0.000;
    atpos[17][2] = 0.750;

    atpos[18][0] = 0.375;
    atpos[18][1] = 0.625;
    atpos[18][2] = 0.875;

    atpos[19][0] = 0.625;
    atpos[19][1] = 0.375;
    atpos[19][2] = 0.875;

    initialized = 1;
    break;

  case 6:
    nucell = 2;
    ntype  = 1;
    noct   = nucell * 3;
    ntetra = nucell * 6;

    memory->create(atpos, nucell + noct + ntetra, 3, "BCC110:atpos");
    memory->create(attyp, nucell + noct + ntetra, "BCC110:attyp");
    
    latvec[0][0] =  sqrt(6.)/3.;
    latvec[0][1] = -sqrt(3.)/6.;
    latvec[1][1] =  sqrt(3.)/2.;
    latvec[2][2] =  sqrt(2.);

    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    for (int i = nucell; i < nucell + noct + ntetra; ++i) attyp[i] = 2;
    atpos[0][0] = 0.000;
    atpos[0][1] = 0.000;
    atpos[0][2] = 0.000;

    atpos[1][0] = 0.500;
    atpos[1][1] = 0.500;
    atpos[1][2] = 0.500;

    atpos[2][0] = 0.500;
    atpos[2][1] = 0.500;
    atpos[2][2] = 0.000;

    atpos[3][0] = 0.250;
    atpos[3][1] = 0.750;
    atpos[3][2] = 0.250;

    atpos[4][0] = 0.750;
    atpos[4][1] = 0.250;
    atpos[4][2] = 0.250;

    atpos[5][0] = 1.000;
    atpos[5][1] = 0.000;
    atpos[5][2] = 0.500;

    atpos[6][0] = 0.750;
    atpos[6][1] = 0.250;
    atpos[6][2] = 0.750;

    atpos[7][0] = 0.250;
    atpos[7][1] = 0.750;
    atpos[7][2] = 0.750;

    atpos[8][0] = 0.625;
    atpos[8][1] = 0.375;
    atpos[8][2] = 0.125;

    atpos[9][0] = 0.375;
    atpos[9][1] = 0.625;
    atpos[9][2] = 0.125;

    atpos[10][0] = 0.000;
    atpos[10][1] = 0.500;
    atpos[10][2] = 0.250;

    atpos[11][0] = 0.500;
    atpos[11][1] = 0.000;
    atpos[11][2] = 0.250;

    atpos[12][0] = 0.875;
    atpos[12][1] = 0.125;
    atpos[12][2] = 0.375;

    atpos[13][0] = 0.125;
    atpos[13][1] = 0.875;
    atpos[13][2] = 0.375;

    atpos[14][0] = 0.875;
    atpos[14][1] = 0.125;
    atpos[14][2] = 0.625;

    atpos[15][0] = 0.125;
    atpos[15][1] = 0.875;
    atpos[15][2] = 0.625;

    atpos[16][0] = 0.000;
    atpos[16][1] = 0.500;
    atpos[16][2] = 0.750;

    atpos[17][0] = 0.500;
    atpos[17][1] = 0.000;
    atpos[17][2] = 0.750;

    atpos[18][0] = 0.375;
    atpos[18][1] = 0.625;
    atpos[18][2] = 0.875;

    atpos[19][0] = 0.625;
    atpos[19][1] = 0.375;
    atpos[19][2] = 0.875;

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
  printf("Please select the type of BCC(111) cell:\n");
  printf("   1. U = [1-10], V = [10-1]; U // x;\n");
  printf("   2. U = [1-10], V = [10-1]; V // y;\n");
  printf("   3. U = [10-1], V = [1-21]; U // x;\n");
  printf("   4. U = [1-21], V = [10-1]; U // x;\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = inumeric(strtok(str, " \t\n\r\f"));
  printf("Your selection : %d", surftype);
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
    noct   = nucell * 3;
    ntetra = nucell * 6;

    memory->create(atpos, nucell + noct + ntetra, 3, "BCC110:atpos");
    memory->create(attyp, nucell + noct + ntetra, "BCC110:attyp");
    
    latvec[0][0] =  sqrt(2.);
    latvec[1][0] = -sqrt(0.5);
    latvec[1][1] =  sqrt(1.5);
    latvec[2][2] =  sqrt(0.75);
    
    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    for (int i = nucell; i < nucell + noct + ntetra; ++i) attyp[i] = 2;

    atpos[0][0] = 0.000;
    atpos[0][1] = 0.000;
    atpos[0][2] = 0.000;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 2./3.;
    atpos[1][2] = 1./3.;

    atpos[2][0] = 2./3.;
    atpos[2][1] = 1./3.;
    atpos[2][2] = 2./3.;

    atpos[3][0] = 0.000;
    atpos[3][1] = 0.500;
    atpos[3][2] = 0.000;

    atpos[4][0] = 0.500;
    atpos[4][1] = 0.500;
    atpos[4][2] = 0.000;

    atpos[5][0] = 0.500;
    atpos[5][1] = 0.000;
    atpos[5][2] = 0.000;

    atpos[6][0] = 5./6.;
    atpos[6][1] = 2./3.;
    atpos[6][2] = 1./3.;

    atpos[7][0] = 5./6.;
    atpos[7][1] = 1./6.;
    atpos[7][2] = 1./3.;

    atpos[8][0] = 1./3.;
    atpos[8][1] = 1./6.;
    atpos[8][2] = 1./3.;

    atpos[9][0] = 1./6.;
    atpos[9][1] = 1./3.;
    atpos[9][2] = 2./3.;

    atpos[10][0] = 2./3.;
    atpos[10][1] = 5./6.;
    atpos[10][2] = 2./3.;

    atpos[11][0] = 1./6.;
    atpos[11][1] = 5./6.;
    atpos[11][2] = 2./3.;

    atpos[12][0] = 2./3.;
    atpos[12][1] = 1./12.;
    atpos[12][2] = 1./6.;

    atpos[13][0] = 11./12.;
    atpos[13][1] = 7./12.;
    atpos[13][2] = 1./6.;

    atpos[14][0] = 5./12.;
    atpos[14][1] = 1./3.;
    atpos[14][2] = 1./6.;

    atpos[15][0] = 5./12.;
    atpos[15][1] = 1./12.;
    atpos[15][2] = 1./6.;

    atpos[16][0] = 2./3.;
    atpos[16][1] = 7./12.;
    atpos[16][2] = 1./6.;

    atpos[17][0] = 11./12.;
    atpos[17][1] = 1./3.;
    atpos[17][2] = 1./6.;

    atpos[18][0] = 0.750;
    atpos[18][1] = 0.750;
    atpos[18][2] = 0.500;

    atpos[19][0] = 0.000;
    atpos[19][1] = 0.250;
    atpos[19][2] = 0.500;

    atpos[20][0] = 0.000;
    atpos[20][1] = 0.750;
    atpos[20][2] = 0.500;

    atpos[21][0] = 0.250;
    atpos[21][1] = 0.250;
    atpos[21][2] = 0.500;

    atpos[22][0] = 0.250;
    atpos[22][1] = 0.000;
    atpos[22][2] = 0.500;

    atpos[23][0] = 0.750;
    atpos[23][1] = 0.000;
    atpos[23][2] = 0.500;

    atpos[24][0] = 1./3.;
    atpos[24][1] = 5./12.;
    atpos[24][2] = 5./6.;

    atpos[25][0] = 7./12.;
    atpos[25][1] = 2./3.;
    atpos[25][2] = 5./6.;

    atpos[26][0] = 1./12.;
    atpos[26][1] = 5./12.;
    atpos[26][2] = 5./6.;

    atpos[27][0] = 7./12.;
    atpos[27][1] = 11./12.;
    atpos[27][2] = 5./6.;

    atpos[28][0] = 1./12.;
    atpos[28][1] = 2./3.;
    atpos[28][2] = 5./6.;

    atpos[29][0] = 1./3.;
    atpos[29][1] = 11./12.;
    atpos[29][2] = 5./6.;

    initialized = 1;
    break;

   case 2:
    nucell = 3;
    ntype  = 1;
    noct   = nucell * 3;
    ntetra = nucell * 6;

    memory->create(atpos, nucell + noct + ntetra, 3, "BCC110:atpos");
    memory->create(attyp, nucell + noct + ntetra, "BCC110:attyp");
    
    latvec[0][0] =  sqrt(1.5);
    latvec[0][1] = -sqrt(0.5);
    latvec[1][1] =  sqrt(2.0);
    latvec[2][2] =  sqrt(0.75);
    
    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    for (int i = nucell; i < nucell + noct + ntetra; ++i) attyp[i] = 2;

    atpos[0][0] = 0.000;
    atpos[0][1] = 0.000;
    atpos[0][2] = 0.000;

    atpos[1][0] = 1./3.;
    atpos[1][1] = 2./3.;
    atpos[1][2] = 1./3.;

    atpos[2][0] = 2./3.;
    atpos[2][1] = 1./3.;
    atpos[2][2] = 2./3.;

    atpos[3][0] = 0.000;
    atpos[3][1] = 0.500;
    atpos[3][2] = 0.000;

    atpos[4][0] = 0.500;
    atpos[4][1] = 0.500;
    atpos[4][2] = 0.000;

    atpos[5][0] = 0.500;
    atpos[5][1] = 0.000;
    atpos[5][2] = 0.000;

    atpos[6][0] = 5./6.;
    atpos[6][1] = 2./3.;
    atpos[6][2] = 1./3.;

    atpos[7][0] = 5./6.;
    atpos[7][1] = 1./6.;
    atpos[7][2] = 1./3.;

    atpos[8][0] = 1./3.;
    atpos[8][1] = 1./6.;
    atpos[8][2] = 1./3.;

    atpos[9][0] = 1./6.;
    atpos[9][1] = 1./3.;
    atpos[9][2] = 2./3.;

    atpos[10][0] = 2./3.;
    atpos[10][1] = 5./6.;
    atpos[10][2] = 2./3.;

    atpos[11][0] = 1./6.;
    atpos[11][1] = 5./6.;
    atpos[11][2] = 2./3.;

    atpos[12][0] = 2./3.;
    atpos[12][1] = 1./12.;
    atpos[12][2] = 1./6.;

    atpos[13][0] = 11./12.;
    atpos[13][1] = 7./12.;
    atpos[13][2] = 1./6.;

    atpos[14][0] = 5./12.;
    atpos[14][1] = 1./3.;
    atpos[14][2] = 1./6.;

    atpos[15][0] = 5./12.;
    atpos[15][1] = 1./12.;
    atpos[15][2] = 1./6.;

    atpos[16][0] = 2./3.;
    atpos[16][1] = 7./12.;
    atpos[16][2] = 1./6.;

    atpos[17][0] = 11./12.;
    atpos[17][1] = 1./3.;
    atpos[17][2] = 1./6.;

    atpos[18][0] = 0.750;
    atpos[18][1] = 0.750;
    atpos[18][2] = 0.500;

    atpos[19][0] = 0.000;
    atpos[19][1] = 0.250;
    atpos[19][2] = 0.500;

    atpos[20][0] = 0.000;
    atpos[20][1] = 0.750;
    atpos[20][2] = 0.500;

    atpos[21][0] = 0.250;
    atpos[21][1] = 0.250;
    atpos[21][2] = 0.500;

    atpos[22][0] = 0.250;
    atpos[22][1] = 0.000;
    atpos[22][2] = 0.500;

    atpos[23][0] = 0.750;
    atpos[23][1] = 0.000;
    atpos[23][2] = 0.500;

    atpos[24][0] = 1./3.;
    atpos[24][1] = 5./12.;
    atpos[24][2] = 5./6.;

    atpos[25][0] = 7./12.;
    atpos[25][1] = 2./3.;
    atpos[25][2] = 5./6.;

    atpos[26][0] = 1./12.;
    atpos[26][1] = 5./12.;
    atpos[26][2] = 5./6.;

    atpos[27][0] = 7./12.;
    atpos[27][1] = 11./12.;
    atpos[27][2] = 5./6.;

    atpos[28][0] = 1./12.;
    atpos[28][1] = 2./3.;
    atpos[28][2] = 5./6.;

    atpos[29][0] = 1./3.;
    atpos[29][1] = 11./12.;
    atpos[29][2] = 5./6.;
    
    initialized = 1;
    break;

  case 3:
    nucell =  6;
    ntype  =  1;
    noct   = nucell * 3;
    ntetra = nucell * 6;

    memory->create(atpos, nucell + noct + ntetra, 3, "BCC110:atpos");
    memory->create(attyp, nucell + noct + ntetra, "BCC110:attyp");
    
    latvec[0][0] = sqrt(2.);
    latvec[1][1] = sqrt(6.);
    latvec[2][2] = sqrt(0.75);
    
    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    for (int i = nucell; i < nucell + noct + ntetra; ++i) attyp[i] = 2;

    atpos[0][0] = 0.000;
    atpos[0][1] = 0.000;
    atpos[0][2] = 0.000;

    atpos[1][0] = 0.500;
    atpos[1][1] = 0.500;
    atpos[1][2] = 0.000;

    atpos[2][0] = 0.000;
    atpos[2][1] = 2./3.;
    atpos[2][2] = 1./3.;

    atpos[3][0] = 0.500;
    atpos[3][1] = 1./6.;
    atpos[3][2] = 1./3.;

    atpos[4][0] = 0.000;
    atpos[4][1] = 1./3.;
    atpos[4][2] = 2./3.;

    atpos[5][0] = 0.500;
    atpos[5][1] = 5./6.;
    atpos[5][2] = 2./3.;

    atpos[6][0] = 0.250;
    atpos[6][1] = 0.250;
    atpos[6][2] = 0.000;

    atpos[7][0] = 0.750;
    atpos[7][1] = 0.750;
    atpos[7][2] = 0.000;

    atpos[8][0] = 0.000;
    atpos[8][1] = 0.500;
    atpos[8][2] = 0.000;

    atpos[9][0] = 0.500;
    atpos[9][1] = 0.000;
    atpos[9][2] = 0.000;

    atpos[10][0] = 0.250;
    atpos[10][1] = 0.750;
    atpos[10][2] = 0.000;

    atpos[11][0] = 0.750;
    atpos[11][1] = 0.250;
    atpos[11][2] = 0.000;

    atpos[12][0] = 0.250;
    atpos[12][1] = 5./12.;
    atpos[12][2] = 1./3.;

    atpos[13][0] = 0.750;
    atpos[13][1] = 11./12.;
    atpos[13][2] = 1./3.;

    atpos[14][0] = 0.000;
    atpos[14][1] = 1./6.;
    atpos[14][2] = 1./3.;

    atpos[15][0] = 0.500;
    atpos[15][1] = 2./3.;
    atpos[15][2] = 1./3.;

    atpos[16][0] = 0.250;
    atpos[16][1] = 11./12.;
    atpos[16][2] = 1./3.;

    atpos[17][0] = 0.750;
    atpos[17][1] = 5./12.;
    atpos[17][2] = 1./3.;

    atpos[18][0] = 0.250;
    atpos[18][1] = 1./12.;
    atpos[18][2] = 2./3.;

    atpos[19][0] = 0.750;
    atpos[19][1] = 7./12.;
    atpos[19][2] = 2./3.;

    atpos[20][0] = 0.500;
    atpos[20][1] = 1./3.;
    atpos[20][2] = 2./3.;

    atpos[21][0] = 0.750;
    atpos[21][1] = 1./12.;
    atpos[21][2] = 2./3.;

    atpos[22][0] = 0.000;
    atpos[22][1] = 5./6.;
    atpos[22][2] = 2./3.;

    atpos[23][0] = 0.250;
    atpos[23][1] = 7./12.;
    atpos[23][2] = 2./3.;

    atpos[24][0] = 0.375;
    atpos[24][1] = 17./24.;
    atpos[24][2] = 1./6.;

    atpos[25][0] = 0.875;
    atpos[25][1] = 5./24.;
    atpos[25][2] = 1./6.;

    atpos[26][0] = 0.250;
    atpos[26][1] = 1./3.;
    atpos[26][2] = 1./6.;

    atpos[27][0] = 0.750;
    atpos[27][1] = 5./6.;
    atpos[27][2] = 1./6.;

    atpos[28][0] = 0.375;
    atpos[28][1] = 23./24.;
    atpos[28][2] = 1./6.;

    atpos[29][0] = 0.875;
    atpos[29][1] = 11./24.;
    atpos[29][2] = 1./6.;

    atpos[30][0] = 0.250;
    atpos[30][1] = 5./6.;
    atpos[30][2] = 1./6.;

    atpos[31][0] = 0.750;
    atpos[31][1] = 1./3.;
    atpos[31][2] = 1./6.;

    atpos[32][0] = 0.125;
    atpos[32][1] = 11./24.;
    atpos[32][2] = 1./6.;

    atpos[33][0] = 0.625;
    atpos[33][1] = 23./24.;
    atpos[33][2] = 1./6.;

    atpos[34][0] = 0.125;
    atpos[34][1] = 5./24.;
    atpos[34][2] = 1./6.;

    atpos[35][0] = 0.625;
    atpos[35][1] = 17./24.;
    atpos[35][2] = 1./6.;

    atpos[36][0] = 0.250;
    atpos[36][1] = 0.500;
    atpos[36][2] = 0.500;

    atpos[37][0] = 0.750;
    atpos[37][1] = 0.000;
    atpos[37][2] = 0.500;

    atpos[38][0] = 0.125;
    atpos[38][1] = 0.125;
    atpos[38][2] = 0.500;

    atpos[39][0] = 0.625;
    atpos[39][1] = 0.625;
    atpos[39][2] = 0.500;

    atpos[40][0] = 0.375;
    atpos[40][1] = 0.375;
    atpos[40][2] = 0.500;

    atpos[41][0] = 0.875;
    atpos[41][1] = 0.875;
    atpos[41][2] = 0.500;

    atpos[42][0] = 0.250;
    atpos[42][1] = 0.000;
    atpos[42][2] = 0.500;

    atpos[43][0] = 0.750;
    atpos[43][1] = 0.500;
    atpos[43][2] = 0.500;

    atpos[44][0] = 0.125;
    atpos[44][1] = 0.875;
    atpos[44][2] = 0.500;

    atpos[45][0] = 0.625;
    atpos[45][1] = 0.375;
    atpos[45][2] = 0.500;

    atpos[46][0] = 0.375;
    atpos[46][1] = 0.625;
    atpos[46][2] = 0.500;

    atpos[47][0] = 0.875;
    atpos[47][1] = 0.125;
    atpos[47][2] = 0.500;

    atpos[48][0] = 0.375;
    atpos[48][1] = 1./24.;
    atpos[48][2] = 5./6.;

    atpos[49][0] = 0.875;
    atpos[49][1] = 13./24.;
    atpos[49][2] = 5./6.;

    atpos[50][0] = 0.125;
    atpos[50][1] = 13./24.;
    atpos[50][2] = 5./6.;

    atpos[51][0] = 0.625;
    atpos[51][1] = 1./24.;
    atpos[51][2] = 5./6.;

    atpos[52][0] = 0.250;
    atpos[52][1] = 1./6.;
    atpos[52][2] = 5./6.;

    atpos[53][0] = 0.750;
    atpos[53][1] = 2./3.;
    atpos[53][2] = 5./6.;

    atpos[54][0] = 0.250;
    atpos[54][1] = 2./3.;
    atpos[54][2] = 5./6.;

    atpos[55][0] = 0.750;
    atpos[55][1] = 1./6.;
    atpos[55][2] = 5./6.;

    atpos[56][0] = 0.375;
    atpos[56][1] = 7./24.;
    atpos[56][2] = 5./6.;

    atpos[57][0] = 0.875;
    atpos[57][1] = 19./24.;
    atpos[57][2] = 5./6.;

    atpos[58][0] = 0.125;
    atpos[58][1] = 19./24.;
    atpos[58][2] = 5./6.;

    atpos[59][0] = 0.625;
    atpos[59][1] = 7./24.;
    atpos[59][2] = 5./6.;
    
    initialized = 1;
    break;

  case 4:
    nucell =    6;
    ntype  =  1;
    noct   = nucell * 3;
    ntetra = nucell * 6;

    memory->create(atpos, nucell + noct + ntetra, 3, "BCC110:atpos");
    memory->create(attyp, nucell + noct + ntetra, "BCC110:attyp");
    
    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    for (int i = nucell; i < nucell + noct + ntetra; ++i) attyp[i] = 2;
    
    latvec[0][0] = sqrt(6.);
    latvec[1][1] = sqrt(2.);
    latvec[2][2] = sqrt(0.75);
    
    atpos[0][0] = 0.000;
    atpos[0][1] = 0.000;
    atpos[0][2] = 0.000;

    atpos[1][0] = 0.500;
    atpos[1][1] = 0.500;
    atpos[1][2] = 0.000;

    atpos[2][0] = 1./3.;
    atpos[2][1] = 0.000;
    atpos[2][2] = 1./3.;

    atpos[3][0] = 5./6.;
    atpos[3][1] = 0.500;
    atpos[3][2] = 1./3.;

    atpos[4][0] = 2./3.;
    atpos[4][1] = 0.000;
    atpos[4][2] = 2./3.;

    atpos[5][0] = 1./6.;
    atpos[5][1] = 0.500;
    atpos[5][2] = 2./3.;

    atpos[6][0] = 0.750;
    atpos[6][1] = 0.250;
    atpos[6][2] = 0.000;

    atpos[7][0] = 0.250;
    atpos[7][1] = 0.750;
    atpos[7][2] = 0.000;

    atpos[8][0] = 0.500;
    atpos[8][1] = 0.000;
    atpos[8][2] = 0.000;

    atpos[9][0] = 0.000;
    atpos[9][1] = 0.500;
    atpos[9][2] = 0.000;

    atpos[10][0] = 0.250;
    atpos[10][1] = 0.250;
    atpos[10][2] = 0.000;

    atpos[11][0] = 0.750;
    atpos[11][1] = 0.750;
    atpos[11][2] = 0.000;

    atpos[12][0] = 7./12.;
    atpos[12][1] = 0.250;
    atpos[12][2] = 1./3.;

    atpos[13][0] = 1./12.;
    atpos[13][1] = 0.750;
    atpos[13][2] = 1./3.;

    atpos[14][0] = 5./6.;
    atpos[14][1] = 0.000;
    atpos[14][2] = 1./3.;

    atpos[15][0] = 1./3.;
    atpos[15][1] = 0.500;
    atpos[15][2] = 1./3.;

    atpos[16][0] = 1./12.;
    atpos[16][1] = 0.250;
    atpos[16][2] = 1./3.;

    atpos[17][0] = 7./12.;
    atpos[17][1] = 0.750;
    atpos[17][2] = 1./3.;

    atpos[18][0] = 5./12.;
    atpos[18][1] = 0.250;
    atpos[18][2] = 2./3.;

    atpos[19][0] = 11./12.;
    atpos[19][1] = 0.750;
    atpos[19][2] = 2./3.;

    atpos[20][0] = 11./12.;
    atpos[20][1] = 0.250;
    atpos[20][2] = 2./3.;

    atpos[21][0] = 5./12.;
    atpos[21][1] = 0.750;
    atpos[21][2] = 2./3.;

    atpos[22][0] = 1./6.;
    atpos[22][1] = 0.000;
    atpos[22][2] = 2./3.;

    atpos[23][0] = 2./3.;
    atpos[23][1] = 0.500;
    atpos[23][2] = 2./3.;

    atpos[24][0] = 7./24.;
    atpos[24][1] = 0.375;
    atpos[24][2] = 1./6.;

    atpos[25][0] = 19./24.;
    atpos[25][1] = 0.875;
    atpos[25][2] = 1./6.;

    atpos[26][0] = 2./3.;
    atpos[26][1] = 0.250;
    atpos[26][2] = 1./6.;

    atpos[27][0] = 1./6.;
    atpos[27][1] = 0.750;
    atpos[27][2] = 1./6.;

    atpos[28][0] = 1./24.;
    atpos[28][1] = 0.375;
    atpos[28][2] = 1./6.;

    atpos[29][0] = 13./24.;
    atpos[29][1] = 0.875;
    atpos[29][2] = 1./6.;

    atpos[30][0] = 1./6.;
    atpos[30][1] = 0.250;
    atpos[30][2] = 1./6.;

    atpos[31][0] = 2./3.;
    atpos[31][1] = 0.750;
    atpos[31][2] = 1./6.;

    atpos[32][0] = 13./24.;
    atpos[32][1] = 0.125;
    atpos[32][2] = 1./6.;

    atpos[33][0] = 1./24.;
    atpos[33][1] = 0.625;
    atpos[33][2] = 1./6.;

    atpos[34][0] = 19./24.;
    atpos[34][1] = 0.125;
    atpos[34][2] = 1./6.;

    atpos[35][0] = 7./24.;
    atpos[35][1] = 0.625;
    atpos[35][2] = 1./6.;

    atpos[36][0] = 0.500;
    atpos[36][1] = 0.250;
    atpos[36][2] = 0.500;

    atpos[37][0] = 0.000;
    atpos[37][1] = 0.750;
    atpos[37][2] = 0.500;

    atpos[38][0] = 0.875;
    atpos[38][1] = 0.125;
    atpos[38][2] = 0.500;

    atpos[39][0] = 0.375;
    atpos[39][1] = 0.625;
    atpos[39][2] = 0.500;

    atpos[40][0] = 0.625;
    atpos[40][1] = 0.375;
    atpos[40][2] = 0.500;

    atpos[41][0] = 0.125;
    atpos[41][1] = 0.875;
    atpos[41][2] = 0.500;

    atpos[42][0] = 0.000;
    atpos[42][1] = 0.250;
    atpos[42][2] = 0.500;

    atpos[43][0] = 0.500;
    atpos[43][1] = 0.750;
    atpos[43][2] = 0.500;

    atpos[44][0] = 0.125;
    atpos[44][1] = 0.125;
    atpos[44][2] = 0.500;

    atpos[45][0] = 0.625;
    atpos[45][1] = 0.625;
    atpos[45][2] = 0.500;

    atpos[46][0] = 0.375;
    atpos[46][1] = 0.375;
    atpos[46][2] = 0.500;

    atpos[47][0] = 0.875;
    atpos[47][1] = 0.875;
    atpos[47][2] = 0.500;

    atpos[48][0] = 23./24.;
    atpos[48][1] = 0.375;
    atpos[48][2] = 5./6.;

    atpos[49][0] = 11./24.;
    atpos[49][1] = 0.875;
    atpos[49][2] = 5./6.;

    atpos[50][0] = 11./24.;
    atpos[50][1] = 0.125;
    atpos[50][2] = 5./6.;

    atpos[51][0] = 23./24.;
    atpos[51][1] = 0.625;
    atpos[51][2] = 5./6.;

    atpos[52][0] = 5./6.;
    atpos[52][1] = 0.250;
    atpos[52][2] = 5./6.;

    atpos[53][0] = 1./3.;
    atpos[53][1] = 0.750;
    atpos[53][2] = 5./6.;

    atpos[54][0] = 1./3.;
    atpos[54][1] = 0.250;
    atpos[54][2] = 5./6.;

    atpos[55][0] = 5./6.;
    atpos[55][1] = 0.750;
    atpos[55][2] = 5./6.;

    atpos[56][0] = 17./24.;
    atpos[56][1] = 0.375;
    atpos[56][2] = 5./6.;

    atpos[57][0] = 5./24.;
    atpos[57][1] = 0.875;
    atpos[57][2] = 5./6.;

    atpos[58][0] = 5./24.;
    atpos[58][1] = 0.125;
    atpos[58][2] = 5./6.;

    atpos[59][0] = 17./24.;
    atpos[59][1] = 0.625;
    atpos[59][2] = 5./6.;

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
  printf("Please select the type of BCC(111) cell:\n");
  printf("   1. U = [1-10], V = [.5,.5,-.5]; U // x\n");
  printf("   2. U = [.5,.5,-.5], V = [-110]; U // x\n");
  printf("Your choice [%d]: ", surftype);
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) surftype = inumeric(strtok(str, " \t\n\r\f"));
  printf("Your selection : %d", surftype);
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
    noct   = nucell * 3;
    ntetra = nucell * 6;

    memory->create(atpos, nucell + noct + ntetra, 3, "BCC110:atpos");
    memory->create(attyp, nucell + noct + ntetra, "BCC110:attyp");
    
    for (int i = 0; i < nucell; ++i) attyp[i] = 1;
    for (int i = nucell; i < nucell + noct + ntetra; ++i) attyp[i] = 2;
    
    latvec[0][0] = sqrt(2.);
    latvec[1][1] = sqrt(0.75);
    latvec[2][2] = sqrt(6.);
    
    atpos[0][0] = 0.000;
    atpos[0][1] = 0.000;
    atpos[0][2] = 0.000;

    atpos[1][0] = 0.500;
    atpos[1][1] = 2./3.;
    atpos[1][2] = 1./6.;

    atpos[2][0] = 0.000;
    atpos[2][1] = 1./3.;
    atpos[2][2] = 1./3.;

    atpos[3][0] = 0.500;
    atpos[3][1] = 0.000;
    atpos[3][2] = 0.500;

    atpos[4][0] = 0.000;
    atpos[4][1] = 2./3.;
    atpos[4][2] = 2./3.;

    atpos[5][0] = 0.500;
    atpos[5][1] = 1./3.;
    atpos[5][2] = 5./6.;

    atpos[6][0] = 0.500;
    atpos[6][1] = 0.000;
    atpos[6][2] = 0.000;

    atpos[7][0] = 0.750;
    atpos[7][1] = 1./3.;
    atpos[7][2] = 1./12.;

    atpos[8][0] = 0.250;
    atpos[8][1] = 1./3.;
    atpos[8][2] = 1./12.;

    atpos[9][0] = 0.000;
    atpos[9][1] = 2./3.;
    atpos[9][2] = 1./6.;

    atpos[10][0] = 0.750;
    atpos[10][1] = 0.000;
    atpos[10][2] = 0.250;

    atpos[11][0] = 0.250;
    atpos[11][1] = 0.000;
    atpos[11][2] = 0.250;

    atpos[12][0] = 0.500;
    atpos[12][1] = 1./3.;
    atpos[12][2] = 1./3.;

    atpos[13][0] = 0.750;
    atpos[13][1] = 2./3.;
    atpos[13][2] = 5./12.;

    atpos[14][0] = 0.250;
    atpos[14][1] = 2./3.;
    atpos[14][2] = 5./12.;

    atpos[15][0] = 0.000;
    atpos[15][1] = 0.000;
    atpos[15][2] = 0.500;

    atpos[16][0] = 0.250;
    atpos[16][1] = 1./3.;
    atpos[16][2] = 7./12.;

    atpos[17][0] = 0.750;
    atpos[17][1] = 1./3.;
    atpos[17][2] = 7./12.;

    atpos[18][0] = 0.500;
    atpos[18][1] = 2./3.;
    atpos[18][2] = 2./3.;

    atpos[19][0] = 0.250;
    atpos[19][1] = 0.000;
    atpos[19][2] = 0.750;

    atpos[20][0] = 0.750;
    atpos[20][1] = 0.000;
    atpos[20][2] = 0.750;

    atpos[21][0] = 0.000;
    atpos[21][1] = 1./3.;
    atpos[21][2] = 5./6.;

    atpos[22][0] = 0.750;
    atpos[22][1] = 2./3.;
    atpos[22][2] = 11./12.;

    atpos[23][0] = 0.250;
    atpos[23][1] = 2./3.;
    atpos[23][2] = 11./12.;

    atpos[24][0] = 0.250;
    atpos[24][1] = 0.500;
    atpos[24][2] = 0.000;

    atpos[25][0] = 0.750;
    atpos[25][1] = 0.500;
    atpos[25][2] = 0.000;

    atpos[26][0] = 0.375;
    atpos[26][1] = 1./6.;
    atpos[26][2] = 1./24.;

    atpos[27][0] = 0.625;
    atpos[27][1] = 1./6.;
    atpos[27][2] = 1./24.;

    atpos[28][0] = 0.125;
    atpos[28][1] = 0.500;
    atpos[28][2] = 0.125;

    atpos[29][0] = 0.875;
    atpos[29][1] = 0.500;
    atpos[29][2] = 0.125;

    atpos[30][0] = 0.250;
    atpos[30][1] = 1./6.;
    atpos[30][2] = 1./6.;

    atpos[31][0] = 0.750;
    atpos[31][1] = 1./6.;
    atpos[31][2] = 1./6.;

    atpos[32][0] = 0.125;
    atpos[32][1] = 5./6.;
    atpos[32][2] = 5./24.;

    atpos[33][0] = 0.875;
    atpos[33][1] = 5./6.;
    atpos[33][2] = 5./24.;

    atpos[34][0] = 0.625;
    atpos[34][1] = 1./6.;
    atpos[34][2] = 7./24.;

    atpos[35][0] = 0.375;
    atpos[35][1] = 1./6.;
    atpos[35][2] = 7./24.;

    atpos[36][0] = 0.750;
    atpos[36][1] = 5./6.;
    atpos[36][2] = 1./3.;

    atpos[37][0] = 0.250;
    atpos[37][1] = 5./6.;
    atpos[37][2] = 1./3.;

    atpos[38][0] = 0.375;
    atpos[38][1] = 0.500;
    atpos[38][2] = 0.375;

    atpos[39][0] = 0.625;
    atpos[39][1] = 0.500;
    atpos[39][2] = 0.375;

    atpos[40][0] = 0.875;
    atpos[40][1] = 5./6.;
    atpos[40][2] = 11./24.;

    atpos[41][0] = 0.125;
    atpos[41][1] = 5./6.;
    atpos[41][2] = 11./24.;

    atpos[42][0] = 0.750;
    atpos[42][1] = 0.500;
    atpos[42][2] = 0.500;

    atpos[43][0] = 0.250;
    atpos[43][1] = 0.500;
    atpos[43][2] = 0.500;

    atpos[44][0] = 0.875;
    atpos[44][1] = 1./6.;
    atpos[44][2] = 13./24.;

    atpos[45][0] = 0.125;
    atpos[45][1] = 1./6.;
    atpos[45][2] = 13./24.;

    atpos[46][0] = 0.375;
    atpos[46][1] = 0.500;
    atpos[46][2] = 0.625;

    atpos[47][0] = 0.625;
    atpos[47][1] = 0.500;
    atpos[47][2] = 0.625;

    atpos[48][0] = 0.750;
    atpos[48][1] = 1./6.;
    atpos[48][2] = 2./3.;

    atpos[49][0] = 0.250;
    atpos[49][1] = 1./6.;
    atpos[49][2] = 2./3.;

    atpos[50][0] = 0.625;
    atpos[50][1] = 5./6.;
    atpos[50][2] = 17./24.;

    atpos[51][0] = 0.375;
    atpos[51][1] = 5./6.;
    atpos[51][2] = 17./24.;

    atpos[52][0] = 0.875;
    atpos[52][1] = 1./6.;
    atpos[52][2] = 19./24.;

    atpos[53][0] = 0.125;
    atpos[53][1] = 1./6.;
    atpos[53][2] = 19./24.;

    atpos[54][0] = 0.250;
    atpos[54][1] = 5./6.;
    atpos[54][2] = 5./6.;

    atpos[55][0] = 0.750;
    atpos[55][1] = 5./6.;
    atpos[55][2] = 5./6.;

    atpos[56][0] = 0.125;
    atpos[56][1] = 0.500;
    atpos[56][2] = 0.875;

    atpos[57][0] = 0.875;
    atpos[57][1] = 0.500;
    atpos[57][2] = 0.875;

    atpos[58][0] = 0.625;
    atpos[58][1] = 5./6.;
    atpos[58][2] = 23./24.;

    atpos[59][0] = 0.375;
    atpos[59][1] = 5./6.;
    atpos[59][2] = 23./24.;
    
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
  noct   = nucell * 3;
  ntetra = nucell * 6;

  memory->create(atpos, nucell + noct + ntetra, 3, "BCC:atpos");
  memory->create(attyp, nucell + noct + ntetra, "BCC:attyp");
  
  latvec[0][0] = -0.5;
  latvec[0][1] =  0.5;
  latvec[0][2] =  0.5;
  latvec[1][0] =  0.5;
  latvec[1][1] = -0.5;
  latvec[1][2] =  0.5;
  latvec[2][0] =  0.5;
  latvec[2][1] =  0.5;
  latvec[2][2] = -0.5;

  for (int i = 0; i < nucell; ++i) attyp[i] = 1;
  for (int i = nucell; i < nucell + noct + ntetra; ++i) attyp[i] = 2;

  atpos[0][0] = 0.000;
  atpos[0][1] = 0.000;
  atpos[0][2] = 0.000;

  atpos[1][0] = 0.500;
  atpos[1][1] = 0.500;
  atpos[1][2] = 0.000;

  atpos[2][0] = 0.500;
  atpos[2][1] = 0.000;
  atpos[2][2] = 0.500;

  atpos[3][0] = 0.000;
  atpos[3][1] = 0.500;
  atpos[3][2] = 0.500;

  atpos[4][0] = 0.750;
  atpos[4][1] = 0.250;
  atpos[4][2] = 0.500;

  atpos[5][0] = 0.750;
  atpos[5][1] = 0.500;
  atpos[5][2] = 0.250;

  atpos[6][0] = 0.250;
  atpos[6][1] = 0.750;
  atpos[6][2] = 0.500;

  atpos[7][0] = 0.250;
  atpos[7][1] = 0.500;
  atpos[7][2] = 0.750;

  atpos[8][0] = 0.500;
  atpos[8][1] = 0.250;
  atpos[8][2] = 0.750;

  atpos[9][0] = 0.500;
  atpos[9][1] = 0.750;
  atpos[9][2] = 0.250;

  initialized = 1;
return;
}

/* --------------------------------------------------------------------------- */
