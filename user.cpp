#include "user.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "common.h"

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
  printf("a POSCAR file, or simply ENTER to read from stdin: ");
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
  // if read from file, the file would be in the format of vasp5 POSCAR (direct)
  //  Comment
  //  alat  (must be positive)
  //  xx xy xz
  //  yx yy yz
  //  zx zy zz
  //  Element-1 Element-2 ... (read but not used)
  //  ntype1 ntype2 ntype3 ...
  //  Direct (must be direct)
  //  sx1 sy1 sz1
  //  ...

  char str[MAXLINE];
  // scaling factor (lattice constant)
  fgets(str,MAXLINE,fp); if (feof(fp)){fclose(fp); return 2;} // Comment Line
  fgets(str,MAXLINE,fp); if (feof(fp)){fclose(fp); return 2;} // Scaling factor
  alat = numeric(strtok(str, " \t\n\r\f"));
  if (alat <= 0.) return 2;

  // basis vectors
  for (int i = 0; i < 3; ++i){
    fgets(str,MAXLINE,fp); if (feof(fp)){fclose(fp); return i+3;}
    latvec[i][0] = numeric(strtok(str,  " \t\n\r\f"));
    latvec[i][1] = numeric(strtok(NULL, " \t\n\r\f"));
    latvec[i][2] = numeric(strtok(NULL, " \t\n\r\f"));
  }
  // Elements names will be read but not used/stored.
  fgets(str,MAXLINE,fp); if (feof(fp)){fclose(fp); return 6;}
  // # of atoms for each type
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

  // Direct or Cartesian
  fgets(str,MAXLINE,fp); if (feof(fp)){fclose(fp); return 6;}
  int cartesian = 0;
  if (str[0] != 'd' && str[0] != 'D') cartesian = 1;

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

  if (cartesian) car2dir();

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

/*------------------------------------------------------------------------------
 * Private method to do matrix inversion
 *----------------------------------------------------------------------------*/
void USER::GaussJordan(const int n, const double *MatA, double *Mat)
{
  int i,icol,irow,j,k,l,ll,idr,idc;
  int indxc[n],indxr[n],ipiv[n];
  double big, dum, pivinv;

  for (int i=0; i<n*n; i++) Mat[i] = MatA[i];

  for (i=0; i<n; i++) ipiv[i] = 0;
  for (i=0; i<n; i++){
    big = 0.;
    for (j=0; j<n; j++){
      if (ipiv[j] != 1){
        for (k=0; k<n; k++){
          if (ipiv[k] == 0){
            idr = j*n+k;
            if (fabs(Mat[idr]) >= big){
              big  = fabs(Mat[idr]);
              irow = j;
              icol = k;
            }
          }else if (ipiv[k] >1){
            printf("\nError: Singular matrix in double GaussJordan!\n");
          }
        }
      }
    }
    ipiv[icol] += 1;
    if (irow != icol){
      for (l=0; l<n; l++){
        idr  = irow*n+l;
        idc  = icol*n+l;
        dum  = Mat[idr];
        Mat[idr] = Mat[idc];
        Mat[idc] = dum;
      }
    }
    indxr[i] = irow;
    indxc[i] = icol;
    idr = icol*n+icol;
    if (Mat[idr] == 0.) printf("\nError: Singular matrix in double GaussJordan!\n");
    
    pivinv = 1./ Mat[idr];
    Mat[idr] = 1.;
    idr = icol*n;
    for (l=0; l<n; l++) Mat[idr+l] *= pivinv;
    for (ll=0; ll<n; ll++){
      if (ll != icol){
        idc = ll*n+icol;
        dum = Mat[idc];
        Mat[idc] = 0.;
        idc -= icol;
        for (l=0; l<n; l++) Mat[idc+l] -= Mat[idr+l]*dum;
      }
    }
  }
  for (l=n-1; l>=0; l--){
    int rl = indxr[l];
    int cl = indxc[l];
    if (rl != cl){
      for (k=0; k<n; k++){
        idr = k*n+rl;
        idc = k*n+cl;
        dum = Mat[idr];
        Mat[idr] = Mat[idc];
        Mat[idc] = dum;
      }
    }
  }

return;
}

/*------------------------------------------------------------------------------
 * Method to convert cartesian coordinate into fractional
 *----------------------------------------------------------------------------*/
void USER::car2dir()
{
  double invaxis[3][3];
  GaussJordan(3,&latvec[0][0],&invaxis[0][0]);
  double x[nucell][3], **s = atpos;
  for (int i=0; i<nucell; i++)
  for (int idim=0; idim<3; idim++) x[i][idim] = atpos[i][idim];

  for (int i=0; i<nucell; i++){
    for (int idim=0; idim<3; idim++){
      s[i][idim] = 0.;
      for (int jdim=0; jdim<3; jdim++) s[i][idim] += x[i][jdim]*invaxis[jdim][idim];

      while (s[i][idim] >= 1.) s[i][idim] -= 1.;
      while (s[i][idim] <  0.) s[i][idim] += 1.;
    }
  }

return;
}
/* ------------------------------------------------------------------- */
