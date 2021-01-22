#include "driver.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <time.h>
#include <math.h>
#include <vector>
#include <set>
#include "common.h"

/* -----------------------------------------------------------------------------
 * Method to create crystals with interstitial solutes
 * -------------------------------------------------------------------------- */
int Driver::Interstitial()
{
  char str[MAXLINE], *ptr;
  printf("\n\n>>>>>>=========== Crystal with Interstitial Solutes ============<<<<<<\n");
  // ask for lattice info
  while ( ShowMenu(-1) < 1 );
  if (latt->noct + latt->ntetra < 1){
     printf("The lattice chosen does not define potential interstitial sites.\n");
     return 0;
  }
  latt->display();
  memory->create(name,strlen(latt->name),"interstitial:name");
  strcpy(name, latt->name);

  // ask for crystal size
  int leading_dir = 1;
  printf("Please input the extensions of your unit cell in x, y, and z directions: ");
  if (count_words(fgets(str,MAXLINE,stdin)) < 3) exit(2);
  nx = latt->inumeric(strtok(str,  " \t\n\r\f"));
  ny = latt->inumeric(strtok(NULL, " \t\n\r\f"));
  nz = latt->inumeric(strtok(NULL, " \t\n\r\f"));
  nucell = latt->nucell;

  int ncell = nx*ny*nz;
  natom = ncell*nucell;
  if (natom < 1) exit(3);

  // show interstitial sites info
  for (int i = 0; i < 14; ++i) printf("-----");
  printf("\nThere are %d octahedral interstitial sites within lattice %s:\n", latt->noct, name);
  for (int i = 0; i < latt->noct; ++i){
    int idx = i + latt->nucell;
    printf("  %d) [%lg %lg %lg]\n", i+1, latt->atpos[idx][0], latt->atpos[idx][1], latt->atpos[idx][2]);
  }
  printf("\nThere are %d tetrahedral interstitial sites within lattice %s:\n", latt->ntetra, name);
  for (int i = 0; i < latt->ntetra; ++i){
    int idx = i + latt->nucell + latt->noct;
    int idy = i + latt->noct + 1;
    printf("  %d) [%lg %lg %lg]\n", idy, latt->atpos[idx][0], latt->atpos[idx][1], latt->atpos[idx][2]);
  }
  for (int i = 0; i < 14; ++i) printf("-----");

  vector<int> sites; sites.clear();
  printf("\nPlease input the indices of the interstitial sites to add solute atoms: ");
  while (count_words(fgets(str,MAXLINE,stdin)) < 1) continue;
  ptr = strtok(str, " \t\n\r\f");
  while (ptr){
    sites.push_back(atoi(ptr));
    ptr = strtok(NULL, " \t\n\r\f");
  }
  if (sites.size() < 1){
    printf("No interstitial sites defined! Bye~\n\n");
    natom = 0;
    return 0;
  }

  vector<double> fraction; fraction.clear();
  printf("Please input the fraction (positive)  or number(negative) of solutes\nat each site in sequence: ");
  while (count_words(fgets(str,MAXLINE,stdin)) < 1) continue;
  ptr = strtok(str, " \t\n\r\f");
  while (ptr){
    fraction.push_back(atof(ptr));
    ptr = strtok(NULL, " \t\n\r\f");
  }
  if (fraction.size() != sites.size()){
    printf("\nIncomplete info provided for the fraction/number of solutes at each site. Bye~\n\n");
    natom = 0; return 0;
  }
  double tf = 0.;
  for (int i = 0; i < fraction.size(); ++i) tf += fabs(fraction[i]);
  if (tf <= ZERO){
    printf("\nIt seems that no interstitial solutes (%lg) are required. Bye~\n\n", tf);
    natom = 0; return 0;
  }
  
  int noct = latt->noct;
  int ntet = latt->ntetra;
  int *ninter;
  memory->create(ninter, noct+ntet, "ninter");
  for (int i = 0; i < noct+ntet; ++i) ninter[i] = 0;
  for (int i = 0; i < sites.size(); ++i){
    int ii = sites[i] - 1;
    if (ii < 0 || ii >= noct+ntet) continue;

    double fi = fraction[i];
    if (fi >= 0.) ninter[ii] = int(ncell * fi);
    else ninter[ii] = int(-fi);
  }
  int n_interstitial = 0;
  for (int i = 0; i < noct+ntet; ++i) n_interstitial += ninter[i];
  if (n_interstitial < 1){
    printf("\nIt seems that no interstitial solutes are required. Bye~\n\n");
    natom = 0; return 0;
  }
  
  natom += n_interstitial;
  printf("Your system would be of size %d x %d x %d with %d atoms, including %d interstitial solutes.\n",nx,ny,nz,natom, n_interstitial);
  printf("Please indicate which direction should goes fast (1:x; other: z)[1]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0) leading_dir = latt->inumeric(strtok(str, " \t\n\r\f"));

  memory->create(atpos, natom, 3, "interstitial:atpos");
  memory->create(attyp, natom, "interstitial:attyp");
  memory->create(xmap, natom, "interstitial:xmap");
  memory->create(ymap, natom, "interstitial:ymap");
  memory->create(zmap, natom, "interstitial:zmap");
  memory->create(umap, natom, "interstitial:umap");

  double **int_pos;
  int *int_type, *int_site, *int_occ;
  int ni_total = ncell * (noct + ntet);
  memory->create(int_pos,  ni_total, 3, "interstitial:int_pos");
  memory->create(int_type, ni_total, "interstitial:int_pos");
  memory->create(int_site, ni_total, "interstitial:int_pos");
  memory->create(int_occ,  ni_total, "interstitial:int_pos");

  int iatom = 0, iint = 0;
  if ( leading_dir == 1){
    for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i){
      for (int u = 0; u < nucell; ++u){
        xmap[iatom] = i;
        ymap[iatom] = j;
        zmap[iatom] = k;
        umap[iatom] = u;
        attyp[iatom] = latt->attyp[u];
        atpos[iatom][0] = latt->atpos[u][0] + double(i);
        atpos[iatom][1] = latt->atpos[u][1] + double(j);
        atpos[iatom][2] = latt->atpos[u][2] + double(k);
        ++iatom;
      }
      for (int u = nucell; u < nucell+noct+ntet; ++u){
        int_type[iint] = latt->attyp[u];
        int_site[iint] = u - nucell;
        int_pos[iint][0] = latt->atpos[u][0] + double(i);
        int_pos[iint][1] = latt->atpos[u][1] + double(j);
        int_pos[iint][2] = latt->atpos[u][2] + double(k);
        ++iint;
      }
    }

  } else {
    for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
    for (int k = 0; k < nz; ++k){
      for (int u = 0; u < latt->nucell; ++u){
        xmap[iatom] = i;
        ymap[iatom] = j;
        zmap[iatom] = k;
        umap[iatom] = u;
        attyp[iatom] = latt->attyp[u];
        atpos[iatom][0] = latt->atpos[u][0] + double(i);
        atpos[iatom][1] = latt->atpos[u][1] + double(j);
        atpos[iatom][2] = latt->atpos[u][2] + double(k);
        ++iatom;
      }
      for (int u = nucell; u < nucell+noct+ntet; ++u){
        int_type[iint] = latt->attyp[u];
        int_site[iint] = u - nucell;
        int_pos[iint][0] = latt->atpos[u][0] + double(i);
        int_pos[iint][1] = latt->atpos[u][1] + double(j);
        int_pos[iint][2] = latt->atpos[u][2] + double(k);
        ++iint;
      }
    }
  }

  // create random number generator if not created; current time as seed;
  if (random == NULL){
    time_t ctime;
    time(&ctime);
    random = new RanPark((int) ctime);
  }

  // now to define the interstitials
  for (int i = 0; i < ni_total; ++i) int_occ[i] = 0;

  for (int ii = 0; ii < noct+ntet; ++ii){
    if (ninter[ii] < 1) continue;
    sites.clear();
    for (int i = 0; i < ni_total;  ++i)
      if (int_site[i] == ii) sites.push_back(i);

    set<int> list; list.clear();
    while (list.size() < ninter[ii]){
      int idx = random->uniform()*ncell;
      if (idx >= ncell) continue;
      list.insert(idx);
    }
    for (std::set<int>::iterator it = list.begin(); it != list.end(); ++it){
      int id = sites[*it];
      int_occ[id] = 1;
    }
    list.clear();
  }

  for (int i = 0; i < ni_total; ++i){
    if (int_occ[i] == 0) continue;
    attyp[iatom] = int_type[i];
    atpos[iatom][0] = int_pos[i][0];
    atpos[iatom][1] = int_pos[i][1];
    atpos[iatom][2] = int_pos[i][2];
    ++iatom;
  }

  // free temporary memory
  memory->destroy(ninter);
  memory->destroy(int_pos);
  memory->destroy(int_type);
  memory->destroy(int_site);
  memory->destroy(int_occ);

  sites.clear();
  fraction.clear();

  // convert fractional coordinate to cartesian
  double tmp[3];
  for (int i = 0; i < 3; ++i)
  for (int j = 0; j < 3; ++j) latvec[i][j] = latt->latvec[i][j]*latt->alat;

  for (int i = 0; i < natom; ++i){
    for (int idim = 0; idim < 3; ++idim) tmp[idim] = atpos[i][idim];
    atpos[i][0] = tmp[0]*latvec[0][0] + tmp[1]*latvec[1][0] + tmp[2]*latvec[2][0];
    atpos[i][1] = tmp[0]*latvec[0][1] + tmp[1]*latvec[1][1] + tmp[2]*latvec[2][1];
    atpos[i][2] = tmp[0]*latvec[0][2] + tmp[1]*latvec[1][2] + tmp[2]*latvec[2][2];
  }
  for (int i = 0; i < 3; ++i){
    latvec[0][i] *= double(nx);
    latvec[1][i] *= double(ny);
    latvec[2][i] *= double(nz);
  }
  // find the total # of types and # of atoms for each type
  typescan();

  printf(">>>>>>============= End of Xtal with interstitial ==============<<<<<<\n\n");
  return -1;
}

/* ------------------------------------------------------------------- */
