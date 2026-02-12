#include "input.h"
#include "common.h"

/* -----------------------------------------------------------------------------
 * Constructor. If flag = 1, output user inputs as script.inp
 * -------------------------------------------------------------------------- */ 
UserInput::UserInput(int flag, Memory *m)
{
   if (flag) fp = fopen("script.inp", "w");
   memory = m;

   return;
}

/* -----------------------------------------------------------------------------
 * Deconstructor. Output user inputs as required and clear workspace.
 * -------------------------------------------------------------------------- */ 
UserInput::~UserInput()
{
   if (fp){
      fclose(fp);
      fp = NULL;
   }
}

/* -----------------------------------------------------------------------------
 * Read stdin and keep a record of it.
 * -------------------------------------------------------------------------- */ 
int UserInput::read_stdin(char *str)
{
   fgets(str, MAXLINE, stdin);
   if (fp) fprintf(fp, "%s", str);

   return count_words(str);
}

/* -----------------------------------------------------------------------------
 * Method to count # of words in a string, without destroying the string
 * -------------------------------------------------------------------------- */
int UserInput::count_words(const char *line)
{
  int n = strlen(line) + 1;
  char *copy;
  memory->create(copy, n,"copy");
  strcpy(copy,line);

  char *ptr;
  if (ptr = strchr(copy,'#')) *ptr = '\0';

  if (strtok(copy," \t\n\r\f") == NULL) {
    memory->sfree(copy);
    return 0;
  }
  n = 1;
  while (strtok(NULL," \t\n\r\f")) n++;

  memory->sfree(copy);
  return n;
}

/* -------------------------------------------------------------------------- */
