#include "shansyn.h"
/* 

   write all coefficients of a GSH expansion to out, rescaling to BW
   convention

*/
void write_gsh_coeff_set(int ialpha, int lmax, 
			 COMP_PRECISION *a, COMP_PRECISION *b, /* for scalars */
			 COMP_PRECISION *ar, COMP_PRECISION *ai, /* for  re(t++) and im(t++) */
			 COMP_PRECISION *br, COMP_PRECISION *bi ,/* for  re(t++) and im(t++) */
			 COMP_PRECISION scale, FILE *out)
{
  COMP_PRECISION fac;
  int os,l,m,nout,izero;
  fprintf(stdout,"%i %i\n",lmax,ialpha);
  nout=0;
  /* prefactors to scale back to BW GSH convection */
  if(ialpha == 0){
    fac = SQRT_TWO;
    izero = 0;
  }else if(ialpha == 2){
    fac = 1.0/SQRT_THREE_FOURTHS;
    izero = 2;
  }else if(ialpha == 4){
    izero = 4;
    fac = 1/ SQRT_THREE;
  }else{
    fprintf(stderr,"write_gsh_coeff: error: %i undefined as ialpha\n",
	    ialpha);
    exit(-1);
  }
  /* start at the non-zero terms only  */
  for(l=izero;l <= lmax;l++){
    for(m=0;m<=l;m++){
      os = POSLM(l,m);
      switch(ialpha){
      case 0:
	if(m == 0){		/* no rescaling */
	  gsh_print(scale*a[os],&nout,out);
	}else{
	  gsh_print(scale*a[os]*fac,&nout,out);
	  gsh_print(scale*b[os]*fac,&nout,out);
	}
	break;
      case 2:
      case 4:			/* all coeffs rescaled */
	gsh_print(scale*ar[os]*fac,&nout,out);
	gsh_print(scale*ai[os]*fac,&nout,out);
	if(m != 0){
	  gsh_print(scale*br[os]*fac,&nout,out);
	  gsh_print(scale*bi[os]*fac,&nout,out);
	}
	break;
      default:
	fprintf(stderr,"%s: ialpha %i undefined\n","write_gsh_coeff",ialpha);
	break;
      }
    }
  }
  if(nout != 0)
    fprintf(out,"\n");
}

/* print a single coeff and increment counter, make 
   five values in a row */
void gsh_print(COMP_PRECISION a,int *nout,FILE *out)
{
  fprintf(out,"%15.7e ",a);
  *nout += 1;
  if(*nout > 4){
    fprintf(out,"\n");
    *nout = 0;
  }
}












/* 

   read generalized spherical harmoncis in the the ylmv4 format of the
   bwgsh packages, and rescale to physical normalization
   
   this is tested for scalar expansions with ialpha == 0

   the routine will also ensure that the coefficients are othonormalized

   this will increment iread by the number of items read
*/
void read_gsh_coeff(FILE *in,
		    int ialpha, int l, int m, 
		    COMP_PRECISION *a, COMP_PRECISION *b, /* for scalars */
		    COMP_PRECISION *ar, COMP_PRECISION *ai, /* for  re(t++) and im(t++) */
		    COMP_PRECISION *br, COMP_PRECISION *bi, /* for  re(t++) and im(t++) */
		    int *iread)
{
  COMP_PRECISION tmp,fac;
  int os;
  
  os = POSLM(l,m);
  switch(ialpha){
  case 0:
    /* scalar */
    if(fscanf(in,DATA_FSCAN_FORMAT,&tmp) != 1){
      /* read A */
      fprintf(stderr,"%s: gsh read error ia: %i l: %i m: %i A\n","read_gsh_coeff",ialpha,l,m);
      exit(-1);
    }
    *iread += 1;
    a[os] = tmp;		/* no m == 0 rescaling */
    if(m != 0){
      /* read B */
      if(fscanf(in,DATA_FSCAN_FORMAT,&tmp) != 1){
	fprintf(stderr,"%s: gsh read error ia: %i l: %i m: %i B\n","read_gsh_coeff",ialpha,l,m);
	exit(-1);
      }
      *iread += 1;
      b[os] = tmp;
      /* 
	 rescale m !=0 terms to physical, DT convention 
      */
      a[os] /= SQRT_TWO;	/* for m != 0 */
      b[os] /= SQRT_TWO;
    }else{
      b[os] = 0.0;
    }
    break;
  case 2:
  case 4:
    if(ialpha == 2){		/* 2phi, for all terms */
      fac = SQRT_THREE_FOURTHS;
    }else{			/* 4phi  */
      fac = SQRT_THREE;
    }
    /* re/im */

    /* read A */
    if(fscanf(in,DATA_FSCAN_FORMAT,&tmp) != 1){
      fprintf(stderr,"%s: gsh read error ia: %i l: %i m: %i AR\n","read_gsh_coeff",ialpha,l,m);
      exit(-1);
    }   
    *iread += 1;ar[os] = tmp * fac;
    if(fscanf(in,DATA_FSCAN_FORMAT,&tmp) != 1){
      fprintf(stderr,"%s: gsh read error ia: %i l: %i m: %i AI\n","read_gsh_coeff",ialpha,l,m);
      exit(-1);
    }   
    *iread += 1;ai[os] = tmp * fac;
    if(m != 0){
      /* read B */
      /* m != 0 terms */
      if(fscanf(in,DATA_FSCAN_FORMAT,&tmp) != 1){
	fprintf(stderr,"%s: gsh read error ia: %i l: %i m: %i BR\n","read_gsh_coeff",ialpha,l,m);
	exit(-1);
      }
      *iread += 1;br[os] = tmp * fac;
      if(fscanf(in,DATA_FSCAN_FORMAT,&tmp) != 1){
	fprintf(stderr,"%s: gsh read error ia: %i l: %i m: %i BI\n","read_gsh_coeff",ialpha,l,m);
	exit(-1);
      }
      *iread += 1;bi[os] = tmp * fac;
      /* 
	 rescale m !=0 terms to physical, DT convention 
      */
    }else{
      br[os] = bi[os] = 0.0;
    }
    break;
  default:
    fprintf(stderr,"%s: ialpha %i undefined\n","read_gsh_coeff",ialpha);
    break;
  }

}
