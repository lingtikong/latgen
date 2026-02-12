#ifndef USER_H
#define USER_H

#include "lattice.h"

class USER : public lattice {
public:
  
  USER(UserInput *);
  ~USER();

private:
  int from_file(const char *);
  int from_stdin();
  void GaussJordan(const int n, const double *MatA, double *Mat);
  void car2dir();
};
#endif
