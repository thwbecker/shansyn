/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: spear.h,v 1.3 2009/04/16 00:55:53 becker Exp $ */
/* 
   spear.c  related numerical recipes routines

*/

COMP_PRECISION corr_sub(COMP_PRECISION *, COMP_PRECISION *, int , COMP_PRECISION *, int );

void pearsn(COMP_PRECISION [],COMP_PRECISION [], unsigned long ,
	    COMP_PRECISION *,COMP_PRECISION *, COMP_PRECISION *,int);


void spear(COMP_PRECISION [], COMP_PRECISION [], unsigned long, COMP_PRECISION *, COMP_PRECISION *, COMP_PRECISION *, COMP_PRECISION *, COMP_PRECISION *);
void crank(unsigned long, COMP_PRECISION [], COMP_PRECISION *);
void sort2(unsigned long, COMP_PRECISION [], COMP_PRECISION []);
void nrerror(char []);
COMP_PRECISION *vector(long, long);
int *ivector(long, long);
unsigned char *cvector(long, long);
unsigned long *lvector(long, long);
double *dvector(long, long);
COMP_PRECISION **matrix(long, long, long, long);
double **dmatrix(long, long, long, long);
int **imatrix(long, long, long, long);
COMP_PRECISION **submatrix(COMP_PRECISION **, long, long, long, long, long, long);
COMP_PRECISION **convert_matrix(COMP_PRECISION *, long, long, long, long);
COMP_PRECISION ***f3tensor(long, long, long, long, long, long);
void free_vector(COMP_PRECISION *, long, long);
void free_ivector(int *, long, long);
void free_cvector(unsigned char *, long, long);
void free_lvector(unsigned long *, long, long);
void free_dvector(double *, long, long);
void free_matrix(COMP_PRECISION **, long, long, long, long);
void free_dmatrix(double **, long, long, long, long);
void free_imatrix(int **, long, long, long, long);
void free_submatrix(COMP_PRECISION **, long, long, long, long);
void free_convert_matrix(COMP_PRECISION **, long, long, long, long);
void free_f3tensor(COMP_PRECISION ***, long, long, long, long, long, long);
COMP_PRECISION erfcc(double);
COMP_PRECISION betai(double, double, double);
COMP_PRECISION gammln(double);
COMP_PRECISION betacf(double, double, double);
static COMP_PRECISION __spear__sqrarg;
#define SPEAR_SQR(a) ((__spear__sqrarg=(a)) == 0.0 ? 0.0 : __spear__sqrarg*__spear__sqrarg)


