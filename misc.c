#include "shansyn.h"


/* safe callocs */
void mycalloc_cp(COMP_PRECISION **x,int n)
{
  *x = (COMP_PRECISION *)calloc((size_t)n,sizeof(COMP_PRECISION));
  if(! (*x)){
    fprintf(stderr,"memerror\n");
    exit(-1);
  }
}
void mycalloc_dp(DATA_PRECISION **x,int n)
{
  *x = (DATA_PRECISION *)calloc((size_t)n,sizeof(DATA_PRECISION));
  if(! (*x)){
    fprintf(stderr,"memerror\n");
    exit(-1);
  }
}
/* safe mallocs */
void mymalloc_dp(DATA_PRECISION **x,int n)
{
  *x = (DATA_PRECISION *)malloc((size_t)n*sizeof(DATA_PRECISION));
  if(! (*x)){
    fprintf(stderr,"memerror\n");
    exit(-1);
  }
}
void mymalloc_cp(COMP_PRECISION **x,int n)
{
  *x = (COMP_PRECISION *)malloc((size_t)n*sizeof(COMP_PRECISION));
  if(! (*x)){
    fprintf(stderr,"memerror\n");
    exit(-1);
  }
}
/* safe mallocs */
void myrealloc_dp(DATA_PRECISION **x,int n)
{
  *x = (DATA_PRECISION *)realloc(*x,(size_t)n*sizeof(DATA_PRECISION));
  if(! (*x)){
    fprintf(stderr,"memerror\n");
    exit(-1);
  }
}
void myrealloc_cp(COMP_PRECISION **x,int n)
{
  *x = (COMP_PRECISION *)realloc(*x,(size_t)n*sizeof(COMP_PRECISION));
  if(! (*x)){
    fprintf(stderr,"memerror\n");
    exit(-1);
  }
}

void zero_cp(COMP_PRECISION *x,int n)
{
  int i;
  for(i=0;i<n;)
    x[i++] = 0.0;
}
void zero_dp(DATA_PRECISION *x,int n)
{
  int i;
  for(i=0;i<n;)
    x[i++] = 0.0;
}

