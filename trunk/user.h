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
};
#endif
