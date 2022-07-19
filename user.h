#ifndef USER_H
#define USER_H

#include "lattice.h"

class USER : public lattice {
public:
  
  USER();
  ~USER();

private:
  int read_file(const char *);
  int read_stdin();
  void GaussJordan(const int n, const double *MatA, double *Mat);
  void car2dir();
};
#endif
