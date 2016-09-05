/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: write_coefficients.c,v 1.8 2009/04/16 00:55:29 becker Exp becker $ */

#include "shansyn.h"


void write_nonzero_coefficients(COMP_PRECISION **a,
				COMP_PRECISION **b,int lmax,
				COMP_PRECISION sfac,
				int print_lmax,
				FILE *out, int nexp)
{
  int l,m,i,os;
  if(print_lmax == 1)		/* short format */
    fprintf(out,"%i\n",lmax);
  else if(print_lmax == 2)	/* long format */
    fprintf(out,"%i %i %g %i %i %i\n",lmax,0,0.,1,nexp,0);
  for(l=0;l <= lmax;l++)
    for(m=0;m <= l;m++){
      os = POSLM(l,m);
      if(m == 0){
	/* don't write the zero B term for m == 0 */
	for(i=0;i < nexp;i++)
	  fprintf(out,ASCII_DATA_FORMAT,*(a[i]+os));
      }else{
	for(i=0;i < nexp;i++)
	  fprintf(out,COEFF_DATA_FORMAT,*(a[i]+os),*(b[i]+os));
      }
      fprintf(out,"\n");
    }

}
void write_coefficients(COMP_PRECISION **a,
			COMP_PRECISION **b,int lmax,
			COMP_PRECISION sfac,int print_lmax,
			FILE *out, int nexp)
{
  int l,m,i,os;
  /* print header */
  if(print_lmax == 1)
    fprintf(out,"%i\n",lmax);
  else if(print_lmax == 2)	/* long format */
    fprintf(out,"%i %i %g %i %i %i\n",lmax,0,0.,1,nexp,0);
  
  for(l=0;l<=lmax;l++)
    for(m=0;m<=l;m++){
      os = POSLM(l,m);
      for(i=0;i < nexp;i++){
	if((*(a[i]+os) != 0)||(*(b[i]+os)!=0))
	  fprintf(out,COEFF_DATA_FORMAT,*(a[i]+os)*sfac,*(b[i]+os)*sfac);
	else
	  fprintf(out,"              0               0");
      }
      fprintf(out,"\n");
    }
}
void write_vector_coefficients(COMP_PRECISION **amp,
			       COMP_PRECISION **amt,
			       COMP_PRECISION **bmp,
			       COMP_PRECISION **bmt,
			       int lmax,
			       COMP_PRECISION sfac,
			       int print_lmax,
			       FILE *out, int nexp)
{
  int l,m,i,os;
  char sout[200];
  /* make the output string */
  for(i=0;i<4;i++)
    if(i==0)
      sprintf(sout,"%s",ASCII_DATA_FORMAT);
    else
      sprintf(sout,"%s %s",sout,ASCII_DATA_FORMAT);
  
  /* print a header? */
  if(print_lmax == 1)
    fprintf(out,"%i\n",lmax);
  else if(print_lmax == 2)	/* long format */
    fprintf(out,"%i %i %g %i %i %i\n",lmax,0,0.,1,nexp*2,0);

  for(l=0;l <= lmax;l++)
    for(m=0;m <= l;m++){
      os = POSLM(l,m);
      for(i=0;i < nexp;i++){	/*                1                3                   2                4 */
	fprintf(out,sout,
		*(amp[i]+os)*sfac,*(bmp[i]+os)*sfac,
		*(amt[i]+os)*sfac,*(bmt[i]+os)*sfac);
      }
      fprintf(out,"\n");
    }
}


void write_model(struct mod *model,COMP_PRECISION fac,
		 int nexp,FILE *out)
{
  int i,j;


  switch(model[0].radial_type){
  case DISCRETE:{
    printf("%i\n",model[0].n);// number of layers
    for(i=0;i < model[0].n;i++){
      printf("%.8e\n",model[0].d[i]);// depth
      write_model_layer(model,i,fac,TRUE,nexp,out);
    }
    break;
  }
  case CHEBYSHEV:{
    printf("%i\n",-model[0].n);// use negative order as layer number
    printf("%.8e %.8e\n",model[0].dmin,model[0].dmax);// boundaries
    for(i=0;i<model[0].n;i++){// coefficients
      write_model_layer(model,i,fac,TRUE,nexp,out);
    }
    break;
  }
  case SPLINES:{
    printf("%i\n",model[0].n+1000);// use order + 1000 as layer number
    printf("%.8e %.8e\n",model[0].dmin,model[0].dmax);// boundaries
    for(i=0;i<model[0].n;i++){// coefficients
      write_model_layer(model,i,fac,TRUE,nexp,out);
    }
    break;
  }
  case HC_TYPE:{
    for(i=0;i < model[0].n;i++){
      fprintf(stdout,"%i %i %.8e %i %i %i\n",
	      model[0].lmax,i,model[0].d[i],model[0].n,1,0);
      write_model_layer(model,i,fac,TRUE,nexp,out);
    }
    break;
  }
  default:{
    fprintf(stderr,"write_model: type %i undefined\n",
	    model[0].radial_type);
    exit(-1);
  }}
}
void write_model_layer(struct mod *model, int layer,
		       COMP_PRECISION fac,int write_lmax,
		       int nexp, FILE *out)
{
  COMP_PRECISION **a,**b,**amp,**amt,**bmp,**bmt;
  int j;
  /*  */
  a=(COMP_PRECISION **)malloc(sizeof(COMP_PRECISION *)*nexp);
  b=(COMP_PRECISION **)malloc(sizeof(COMP_PRECISION *)*nexp);
  amt=(COMP_PRECISION **)malloc(sizeof(COMP_PRECISION *)*nexp);
  amp=(COMP_PRECISION **)malloc(sizeof(COMP_PRECISION *)*nexp);
  bmt=(COMP_PRECISION **)malloc(sizeof(COMP_PRECISION *)*nexp);
  bmp=(COMP_PRECISION **)malloc(sizeof(COMP_PRECISION *)*nexp);
  for(j=0;j < nexp;j++){
    if(!model[j].vector_field){
      a[j] = model[j].a[layer];
      b[j] = model[j].b[layer];
    }else{
      amp[j] = model[j].amp[layer];
      bmp[j] = model[j].bmp[layer];
      amt[j] = model[j].amt[layer];
      bmt[j] = model[j].bmt[layer];
    }
  }
  if(model[0].is_gsh){		/* GSH type */
    for(j=0;j < nexp;j++)
      write_gsh_coeff_set((model[j].is_gsh-1), model[j].lmax, 
			  a[j],b[j],amp[j],amt[j],
			  bmp[j],bmt[j],fac,out);
  }else{			/* regular */
    if(model[0].vector_field){
      // coefficients
      write_vector_coefficients(amp,amt,bmp,bmt,
				model[0].lmax,fac,
				write_lmax,out,nexp);
    }else{
      write_coefficients(a,b,model[0].lmax,fac,
			 write_lmax,out,nexp);
      
      }
  }

  free(a);free(b);
  free(amp);free(bmp);
  free(amt);free(bmt);
}
