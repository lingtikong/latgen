#ifndef INPUT_H
#define INPUT_H

#include <cstdio>
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "memory.h"

class UserInput {
public:
   UserInput(int, Memory *);
   ~UserInput();

   int read_stdin(char *);

   Memory *memory;
   int count_words(const char *);

private:
   FILE *fp;

};
#endif
