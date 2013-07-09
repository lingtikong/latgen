#ifdef Poly
#include "driver.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <time.h>
#include <math.h>
#include "random.h"

#include "voro++.hh"
using namespace voro;

#define MAXLINE 512
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/*------------------------------------------------------------------------------
 * Method to create polycrystals
 *----------------------------------------------------------------------------*/
void Driver::PolyCrystal()
{
  char str[MAXLINE], *ptr;
  int ngrain = 0, type_by_grain = 0;
  int nmax = 1000;
  double lo[3], hi[3], box[3];
  double **grains, **orient;
  int pbc[3];
  bool bpbc[3];
  FILE *fp = NULL;

  printf("\n\n>>>>>>================ Polycrystal generation ==================<<<<<<\n");
  // ask for lattice info
  while ( ShowMenu(-1) < 1 );
  latt->display();
  memory->create(name,strlen(latt->name)+20,"polycrystal:name");
  sprintf(name,"Polycrystal of %s", latt->name);

  if (latt->perp_x != 1 || latt->perp_y != 1 || latt->perp_z != 1){
    printf("WARNING: your lattice is not othorgonal, might result in wrong config!\n");
  }

  // create random number generator if not created; current time as seed;
  if (random == NULL){
    time_t ctime;
    time(&ctime);
    random = new RanPark((int) ctime);
  }

  grains = orient = NULL;

  // Polycrystal info can be read from file
  printf("\nThe necessary info for polycrystal info can be provided either by stdin, or\n");
  printf("read from a file. In the previous case, the position and the orientation of\n");
  printf("each grain are assigned randomly. While in the latter case, such information\n");
  printf("can be defined explicitly. If this is preferred, input the file name now: ");
  if (count_words(fgets(str,MAXLINE,stdin)) > 0){
    ptr = strtok(str, " \n\t\r\f");
    fp = fopen(ptr, "r");
  }

  // read related info from file
  if (fp){
    int idir = 0;
    fgets(str, MAXLINE, fp);
    while ( !feof(fp) ) {  // first three effective lines: xlo, xhi, xpbc
      if (count_words(str) >= 3){
        lo[idir]  = atof(strtok(str, " \n\t\r\f"));
        hi[idir]  = atof(strtok(NULL," \n\t\r\f"));
        pbc[idir] = atoi(strtok(NULL," \n\t\r\f"));

        box[idir] = hi[idir] - lo[idir];

        if (++idir >= 3) break;
      }

      fgets(str, MAXLINE, fp);
    }

    fgets(str, MAXLINE, fp);
    while ( !feof(fp) ) {   // fourth effective line: ngrain, type-by-grain
      if (count_words(str) >= 2){
        ngrain = atoi(strtok(str, " \n\t\r\f"));
        type_by_grain = atoi(strtok(NULL," \n\t\r\f"));

        break;
      }
      fgets(str, MAXLINE, fp);
    }

    if (ngrain > 0){
      memory->create(grains, ngrain, 3, "grains");
      memory->create(orient, ngrain, 3, "orient");

      int igrain = 0;
      // remaining ngrain effective lines: x y z alpha beta gamma
      // alpha, beta, gamma are in unit of 2*PI; but can be 'ran' for random number within [0, 1]
      fgets(str, MAXLINE, fp);
      while ( !feof(fp) ) {
        if (count_words(str) >= 6){
          grains[igrain][0] = atof(strtok(str, " \n\t\r\f"));
          grains[igrain][1] = atof(strtok(NULL," \n\t\r\f"));
          grains[igrain][2] = atof(strtok(NULL," \n\t\r\f"));
          for (int ii=0; ii<3; ii++){
            ptr = strtok(NULL," \n\t\r\f");
            if (strcmp(ptr, "ran") == 0) orient[igrain][ii] = random->uniform();
            else orient[igrain][ii] = atof(ptr);
          }
  
          if (++igrain >= ngrain) break;
        }

        fgets(str, MAXLINE, fp);
      }
  
      if (igrain < ngrain){
        memory->destroy(grains); memory->destroy(orient);
        grains = orient = NULL;
      }
    }

    fclose(fp); fp = NULL;
  }

  // either asked to input from stdin or read file failed
  if (grains == NULL){
    printf("\nThe necessary information will be read from stdin.\n");
    // ask for box info
    for (int i=0; i<3; i++){
      printf("Please input the lower and upper bound of your box along %c: ", 'X'+i);
      while (1) if (count_words(fgets(str,MAXLINE,stdin)) == 2){
        lo[i] = atof(strtok(str, " \n\t\r\f"));
        hi[i] = atof(strtok(NULL," \n\t\r\f"));
  
        box[i] = hi[i] - lo[i];
        if (box[i] > 0.) break;
      }
    }
    pbc[0] = pbc[1] = pbc[2] = 1;
    printf("Please indicate where you want pbc in x, y, and z (0/1)[1 1 1]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) == 3){
      pbc[0] = atoi(strtok(str, " \t\n\r\f"));
      pbc[1] = atoi(strtok(NULL," \n\t\r\f"));
      pbc[2] = atoi(strtok(NULL," \n\t\r\f"));
    }

    // ask for number of grains in box
    printf("Please input the desired number of grains in your box: ");
    while (1) if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      ngrain = atoi(strtok(str, " \t\n\r\f"));
      if (ngrain > 1) break;
    }

    // create the grains and set their orientations randomly
    memory->create(grains, ngrain, 3, "grains");
    memory->create(orient, ngrain, 3, "orient");

    for (int ii=0; ii<ngrain; ii++){
      for (int idim=0; idim< 3; idim++){
        grains[ii][idim] = lo[idim] + random->uniform()*box[idim];
        orient[ii][idim] = random->uniform();
      }
    }

    // atoms in different grain can be assigned as different type
    printf("Would you like to assign different atomic types for different grains? (y/n)[n]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      ptr = strtok(str, " \n\t\r\f");
      if (strcmp(ptr,"y")==0 || strcmp(ptr,"Y")==0) type_by_grain = 1;
    }
  }

  double cvol = box[0]*box[1]*box[2];
  alat = 1.;
  latvec[0][0] = box[0];
  latvec[1][1] = box[1];
  latvec[2][2] = box[2];

  bpbc[0] = bpbc[1] = bpbc[2] = false;
  for (int idim=0; idim<3; idim++) if (pbc[idim] != 0) bpbc[idim] = true;

  int ngd = MAX(2,ngrain/2);
  // now create the box (container)
  container con(lo[0],hi[0],lo[1],hi[1],lo[2],hi[2],ngd,ngd,ngd,bpbc[0],bpbc[1],bpbc[2],8);
  
  // insert grain centers into the container
  for (int i=1; i<=ngrain; i++){
    int ii= i-1;
    con.put(i, grains[ii][0], grains[ii][1], grains[ii][2]);
  }
    
  // output generating info
  printf("\nYour box is bounded by        : [%g %g], [%g %g], [%g %g]\n",lo[0],hi[0],lo[1],hi[1],lo[2],hi[2]);
  printf("Periodicity in each direction : %d %d %d\n", pbc[0], pbc[1], pbc[2]);
  printf("Number of grains in the box   : %d\n", ngrain);
  printf("The crystalline grain info will be written to grain_part_info and\ngrain_cell_info, respectively.\n\n");
  printf("Now to generate the polycrystal ... "); fflush(stdout);

  fp = fopen("grain_part_info","w");
  fprintf(fp,"# box bounds : %g %g %g %g %g %g\n", lo[0],hi[0],lo[1],hi[1],lo[2],hi[2]);
  fprintf(fp,"# # of grains: %d\n", ngrain);
  fprintf(fp,"# 1  2  3  4  5   6        7    8\n");
  fprintf(fp,"# id xc yc zc vol surfarea rmax natom\n");

  // prepare array for inserting atoms
  memory->grow(atpos,nmax,3,"atpos");
  memory->grow(attyp,nmax,"attyp");

  // Create a loop class to iterate over all of the grains in the box
  c_loop_all cl(con);
  voronoicell c;

  // now to generate the polycrystal
  if (cl.start()) do if (con.compute_cell(c,cl)){
    // get the radius of the grain
    double rmax2 = 0.25*c.max_radius_squared();
    double rmax = sqrt(rmax2);
    double cvol = c.volume();
    double area = c.surface_area();

    int nxmax = int(1.5*rmax/latt->hx+0.5)+1;
    int nymax = int(1.5*rmax/latt->hy+0.5)+1;
    int nzmax = int(1.5*rmax/latt->hz+0.5)+1;

    // get the id and center position of the grain
    double xc, yc, zc, xp[3], rx, ry, rz;
    int id = cl.pid();
    cl.pos(xc,yc,zc);
    
    // rotated the lattice randomly
    double rotated[3][3], ang[3];
    for (int i=0; i<3; i++) ang[i] = orient[id-1][i];
    latt->RotateLattice(&ang[0], rotated);

    int n_in_grain = 0;

    // now to create local atoms
    for (int ix=-nxmax; ix<=nxmax; ix++)
    for (int iy=-nymax; iy<=nymax; iy++)
    for (int iz=-nzmax; iz<=nzmax; iz++)
    for (int iu=0; iu<latt->nucell; iu++){
      rx = double(ix) + latt->atpos[iu][0];
      ry = double(iy) + latt->atpos[iu][1];
      rz = double(iz) + latt->atpos[iu][2];

      xp[0] = (rx*rotated[0][0]+ry*rotated[1][0]+rz*rotated[2][0])*latt->alat;
      xp[1] = (rx*rotated[0][1]+ry*rotated[1][1]+rz*rotated[2][1])*latt->alat;
      xp[2] = (rx*rotated[0][2]+ry*rotated[1][2]+rz*rotated[2][2])*latt->alat;

      double r2 = xp[0]*xp[0] + xp[1]*xp[1] + xp[2]*xp[2];
      if (r2 > rmax2) continue;

      xp[0] += xc; xp[1] += yc; xp[2] += zc;

      int inside = 1;
      for (int idim=0; idim<3; idim++){
        if (pbc[idim]){
          if (xp[idim] < lo[idim]) xp[idim] += box[idim];
          if (xp[idim] >=hi[idim]) xp[idim] -= box[idim];
        }
        if (xp[idim] < lo[idim] || xp[idim] > hi[idim]) inside = 0;
      }
      if (inside != 1) continue;

      int ifnd;
      if (!con.find_voronoi_cell(xp[0],xp[1],xp[2],rx,ry,rz,ifnd)) continue;

      if (ifnd == id){
        if (natom >= nmax){
          nmax += 1000;
          atpos = memory->grow(atpos,nmax,3,"polycrystal:atpos");
          attyp = memory->grow(attyp,nmax,"polycrystal:attyp");
        }

        for (int idim=0; idim<3; idim++) atpos[natom][idim] = xp[idim];

        if (type_by_grain) attyp[natom] = latt->attyp[iu] + (id-1)*latt->ntype;
        else attyp[natom] = latt->attyp[iu];

        natom++; n_in_grain++;
      }
    }

    fprintf(fp,"%d %g %g %g %g %g %g %d\n", id, xc, yc, zc, cvol, area, rmax, n_in_grain);

  } while (cl.inc());

  fclose(fp);
  printf("Done! Total # of atoms: %d\n", natom);
  printf("WARNING: there might be overlaped atoms, check by yourself!!!\n\n");

  // Output files for diagnosic purposes
  con.draw_cells_gnuplot("grain_cell_info");

  // find the total # of types and # of atoms for each type
  typescan();

  // free memory
  memory->destroy(grains);
  memory->destroy(orient);

  printf(">>>>>>============= End of Polycrystal generation ==============<<<<<<\n");
return;
}

/* ------------------------------------------------------------------- */

#endif
