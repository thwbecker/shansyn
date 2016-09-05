#include "shansyn.h"
/* 

   compute how many coefficients are needed

   is_gsh:

   0: regular scalar
   1: GSH scalar
   3: GSH 2 phi
   5: GSH 4 phi
   
   other inputs: lmax
   
   output: lms : number of coefficients 
   vector_field: GSH code
   izero: l start 
   

 */
void select_lms(int is_gsh,int lmax, int *lms, 
		int *vector_field, int *izero)
{
  switch(is_gsh){
  case 0:			/* regular scalar */
    *lms= (lmax+1)*(lmax+2)/2;
    *izero = 0;
    *vector_field=0;
    break;
  case 1:			/* GSH scalar */
    *lms= (lmax+1)*(lmax+1);
    *izero = 0;
    *vector_field=0;
    break;
  case 3:			/* GSH 2phi */
    *lms= (lmax-1)*(2*lmax+6);
    *izero = 2;
    *vector_field=2;
    break;
  case 5:			/* GSH 4phi */
    *lms= (lmax-3)*(2*lmax+10);
    *izero = 4;
    *vector_field=2;
    break;
  default:
    fprintf(stderr,"select_lms: ERROR: is_gsh code should be 0, 1, 3, or 5, but is %i\n",is_gsh);
    exit(-1);
    break;
  }
}
