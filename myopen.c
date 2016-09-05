/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: myopen.c,v 1.2 2001/03/08 03:04:49 becker Exp $ */
#include "myio.h"

FILE *myopen(char *name,char *modus)
{
  FILE *tmp;
  if( (tmp = (FILE *)fopen(name,modus)) == NULL)
    {	fprintf(stderr,"error opening \"%s\"\n",name);
    exit(-1);};
  return ((FILE *)tmp);
}	
FILE *myopen_wn(char *name,char *modus, char *program)
{
  FILE *tmp;
  if( (tmp = (FILE *)fopen(name,modus)) == NULL)
    {	fprintf(stderr,"%s: myopen: error opening \"%s\"\n",program,name);
    exit(-1);};
  return ((FILE *)tmp);
}	
