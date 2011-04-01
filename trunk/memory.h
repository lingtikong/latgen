#ifndef MEMORY_H
#define MEMORY_H

class Memory{
 public:

  void *smalloc(int n, const char *);
  void sfree(void *);

  double **create_2d_double_array(int, int, const char *);
  void destroy_2d_double_array(double **);

  int **create_2d_int_array(int, int, const char *);
  void destroy_2d_int_array(int **);

  double ***create_3d_double_array(int, int, int, const char *);
  void destroy_3d_double_array(double ***);

  double ***create_3d_double_array(int, int, int, int, const char *);
  void destroy_3d_double_array(double ***, int);

  int ***create_3d_int_array(int, int, int, const char *);
  void destroy_3d_int_array(int ***);

  char **create_2d_string_array(int, int, const char *);
  void destroy_2d_string_array(char **);
};

#endif
