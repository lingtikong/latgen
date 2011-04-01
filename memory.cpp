#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "memory.h"

/* ----------------------------------------------------------------------
   safe malloc 
------------------------------------------------------------------------- */

void *Memory::smalloc(int n, const char *name)
{
  if (n == 0) return NULL;
  void *ptr = malloc(n);
  if (ptr == NULL) {
    printf("Failed to allocate %d bytes for array %s",n,name);
    exit(1);
  }
  return ptr;
}

/* ----------------------------------------------------------------------
   safe free 
------------------------------------------------------------------------- */
void Memory::sfree(void *ptr)
{
  if (ptr == NULL) return;
  free(ptr);
}

/* ----------------------------------------------------------------------
   create a 2d double array 
------------------------------------------------------------------------- */

double **Memory::create_2d_double_array(int n1, int n2, const char *name)
{
  double *data = (double *) smalloc(n1*n2*sizeof(double),name);
  double **array = (double **) smalloc(n1*sizeof(double *),name);

  int n = 0;
  for (int i = 0; i < n1; i++) {
    array[i] = &data[n];
    n += n2;
  }

  return array;
}

/* ----------------------------------------------------------------------
   free a 2d double array 
------------------------------------------------------------------------- */
void Memory::destroy_2d_double_array(double **array)
{
  if (array == NULL) return;
  sfree(array[0]);
  sfree(array);
}

/* ----------------------------------------------------------------------
   create a 2d int array
   if either dim is 0, return NULL 
------------------------------------------------------------------------- */
int **Memory::create_2d_int_array(int n1, int n2, const char *name)
{
  if (n1 == 0 || n2 == 0) return NULL;

  int *data = (int *) smalloc(n1*n2*sizeof(int),name);
  int **array = (int **) smalloc(n1*sizeof(int *),name);

  int n = 0;
  for (int i = 0; i < n1; i++) {
    array[i] = &data[n];
    n += n2;
  }

  return array;
}

/* ----------------------------------------------------------------------
   free a 2d int array 
------------------------------------------------------------------------- */
void Memory::destroy_2d_int_array(int **array)
{
  if (array == NULL) return;
  sfree(array[0]);
  sfree(array);
}

/* ----------------------------------------------------------------------
   create a 3d double array 
------------------------------------------------------------------------- */
double ***Memory::create_3d_double_array(int n1, int n2, int n3,
					 const char *name)
{
  int i,j;

  double *data = (double *) smalloc(n1*n2*n3*sizeof(double),name);
  double **plane = (double **) smalloc(n1*n2*sizeof(double *),name);
  double ***array = (double ***) smalloc(n1*sizeof(double **),name);

  int n = 0;
  for (i = 0; i < n1; i++) {
    array[i] = &plane[i*n2];
    for (j = 0; j < n2; j++) {
      plane[i*n2+j] = &data[n];
      n += n3;
    }
  }

  return array;
}

/* ----------------------------------------------------------------------
   free a 3d double array 
------------------------------------------------------------------------- */
void Memory::destroy_3d_double_array(double ***array)
{
  if (array == NULL) return;
  sfree(array[0][0]);
  sfree(array[0]);
  sfree(array);
}

/* ----------------------------------------------------------------------
   create a 3d int array 
------------------------------------------------------------------------- */

int ***Memory::create_3d_int_array(int n1, int n2, int n3, const char *name)
{
  int i,j;

  int *data = (int *) smalloc(n1*n2*n3*sizeof(int),name);
  int **plane = (int **) smalloc(n1*n2*sizeof(int *),name);
  int ***array = (int ***) smalloc(n1*sizeof(int **),name);

  int n = 0;
  for (i = 0; i < n1; i++) {
    array[i] = &plane[i*n2];
    for (j = 0; j < n2; j++) {
      plane[i*n2+j] = &data[n];
      n += n3;
    }
  }

  return array;
}

/* ----------------------------------------------------------------------
   free a 3d int array 
------------------------------------------------------------------------- */
void Memory::destroy_3d_int_array(int ***array)
{
  if (array == NULL) return;
  sfree(array[0][0]);
  sfree(array[0]);
  sfree(array);
}

/* ----------------------------------------------------------------------
   create a 2d string array
   if either dim is 0, return NULL 
------------------------------------------------------------------------- */
char **Memory::create_2d_string_array(int n1, int n2, const char *name)
{
  if (n1 == 0 || n2 == 0) return NULL;

  char *data = (char *) smalloc(n1*n2*sizeof(char),name);
  char **array = (char **) smalloc(n1*sizeof(char *),name);

  int n = 0;
  for (int i = 0; i < n1; i++) {
    array[i] = &data[n];
    n += n2;
  }

  return array;
}

/* ----------------------------------------------------------------------
   free a 2d string array 
------------------------------------------------------------------------- */
void Memory::destroy_2d_string_array(char **array)
{
  if (array == NULL) return;
  sfree(array[0]);
  sfree(array);
}
