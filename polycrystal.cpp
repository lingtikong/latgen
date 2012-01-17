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

/*------------------------------------------------------------------------------
 * Method to create polycrystals
 *----------------------------------------------------------------------------*/
void Driver::PolyCrystal()
{
  char str[MAXLINE];
  int ngrain = 0;
  int nmax = 1000;
  double lo[3], hi[3], box[3];
  int pbc[3];
  bool bpbc[3];

  printf("\n\n>>>>>>================ Polycrystal generation ==================<<<<<<\n");
  // ask for lattice info
  while ( ShowMenu(-1) < 1 );
  latt->display();
  name = memory->create(name,strlen(latt->name)+20,"polycrystal:name");
  sprintf(name,"Polycrystal of %s", latt->name);

  if (latt->perp_x != 1 || latt->perp_y != 1 || latt->perp_z != 1){
    printf("WARNING: your lattice is not othorgonal, might result in wrong config!\n");
  }

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
  double cvol = box[0]*box[1]*box[2];
  pbc[0] = pbc[1] = pbc[2] = 1;
  printf("Please indicate where you want pbc in x, y, and z (0/1)[1 1 1]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) == 3){
    pbc[0] = atoi(strtok(str, " \t\n\r\f"));
    pbc[1] = atoi(strtok(NULL," \n\t\r\f"));
    pbc[2] = atoi(strtok(NULL," \n\t\r\f"));
  }
  alat = 1.;
  latvec[0][0] = box[0];
  latvec[1][1] = box[1];
  latvec[2][2] = box[2];

  // ask for number of grains in box
  printf("Please input the desired number of grains in your box: ");
  while (1) if (count_words(fgets(str,MAXLINE,stdin)) > 0){
    ngrain = atoi(strtok(str, " \t\n\r\f"));
    if (ngrain > 0) break;
  }
  
  bpbc[0] = bpbc[1] = bpbc[2] = false;
  for (int idim=0; idim<3; idim++) if (pbc[idim]) bpbc[idim] = true;
  // now create the box (container)
  container con(lo[0], hi[0], lo[1], hi[1], lo[2], hi[2],10,10,10,bpbc[0],bpbc[1],bpbc[2],8);
  
  // create random number generator if not created; current time as seed;
  if (random == NULL){
    time_t ctime;
    time(&ctime);
    random = new RanPark((int) ctime);
  }

  // insert grain centers into the container
  for (int i=1; i<=ngrain; i++){
    double x = lo[0] + random->uniform()*box[0];
    double y = lo[1] + random->uniform()*box[1];
    double z = lo[2] + random->uniform()*box[2];

    con.put(i,x,y,z);
  }
    
  // output generating info
  printf("\nYour box is bounded by        : [%g %g], [%g %g], [%g %g]\n",lo[0],hi[0],lo[1],hi[1],lo[2],hi[2]);
  printf("Periodicity in each direction : %d %d %d\n", pbc[0], pbc[1], pbc[2]);
  printf("Number of grains in the box   : %d\n", ngrain);
  printf("The crystalline grain info will be written to grain_part_info and\ngrain_cell_info, respectively.\n\n");
  printf("Now to generate the polycrystal ... "); fflush(stdout);

  FILE *fp = fopen("grain_part_info","w");
  fprintf(fp,"# box bounds : %g %g %g %g %g %g\n", lo[0],hi[0],lo[1],hi[1],lo[2],hi[2]);
  fprintf(fp,"# # of grains: %d\n", ngrain);
  fprintf(fp,"# 1  2  3  4  5   6        7    8\n");
  fprintf(fp,"# id xc yc zc vol surfarea rmax natom\n");

  // prepare array for inserting atoms
  atpos = memory->grow(atpos,nmax,3,"atpos");
  attyp = memory->grow(attyp,nmax,"attyp");

  // Create a loop class to iterate over all of the grains in the box
  c_loop_all cl(con);
  voronoicell c;

  // now to generate the polycrystal
  if (cl.start()) do if (con.compute_cell(c,cl)){
    // rotated the lattice randomly
    double **latvec = memory->create(latvec,3,3,"polycrystal:latvec");
    double ang[3];
    for (int i=0; i<3; i++) ang[i] = random->uniform();
    latt->RotateLattice(&ang[0], latvec);

    // get the radius of the grain
    double rmax = sqrt(0.25*c.max_radius_squared());
    double cvol = c.volume();
    double area = c.surface_area();

    int nx = int(rmax/latt->hx)+1;
    int ny = int(rmax/latt->hy)+1;
    int nz = int(rmax/latt->hz)+1;

    // get the id and center position of the grain
    double xc, yc, zc, xp[3], rx, ry, rz;
    int id = cl.pid();
    cl.pos(xc,yc,zc);
    
    int n_in_grain = 0;
    // now to create local atoms
    for (int i=-nx; i<nx; i++)
    for (int j=-ny; j<ny; j++)
    for (int k=-nz; k<nz; k++)
    for (int m=0; m<latt->nucell; m++){
      rx = double(i) + latt->atpos[m][0];
      ry = double(j) + latt->atpos[m][1];
      rz = double(k) + latt->atpos[m][2];

      xp[0] = xc + (rx*latvec[0][0]+ry*latvec[1][0]+rz*latvec[2][0])*latt->alat;
      xp[1] = yc + (rx*latvec[0][1]+ry*latvec[1][1]+rz*latvec[2][1])*latt->alat;
      xp[2] = zc + (rx*latvec[0][2]+ry*latvec[1][2]+rz*latvec[2][2])*latt->alat;

      for (int ii=0; ii<3; ii++){
        if (pbc[ii]){
          if (xp[ii] < lo[ii]) xp[ii] += box[ii];
          if (xp[ii] >=hi[ii]) xp[ii] -= box[ii];
        }
      }

      int ifnd;
      con.find_voronoi_cell(xp[0],xp[1],xp[2],rx,ry,rz,ifnd);

      if (ifnd == id){
        if (natom >= nmax){
          nmax += 1000;
          atpos = memory->grow(atpos,nmax,3,"polycrystal:atpos");
          attyp = memory->grow(attyp,nmax,"polycrystal:attyp");
        }

        for (int idim=0; idim<3; idim++) atpos[natom][idim] = xp[idim];
        attyp[natom++] = latt->attyp[m];
        n_in_grain++;
      }
    }

    fprintf(fp,"%d %g %g %g %g %g %g %d\n", id, xc, yc, zc, cvol, area, rmax, n_in_grain);

  } while (cl.inc());

  fclose(fp);
  printf("Done! Total # of atoms: %d\n", natom);
  printf("WARNING: there might be overlaped atoms, check by yourself!\n\n");

  // Output files for diagnosic purposes
  con.draw_cells_gnuplot("grain_cell_info");

  // find the total # of types and # of atoms for each type
  typescan();

return;
}

#endif
