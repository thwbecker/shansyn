/* 

compute radial correlation in a new way, there is also cradialcorr


*/
#include "shansyn.h"

int main(int argc, char **argv)
{
  COMP_PRECISION z1,z2,zlim[2],dz,r,zrange,zleft,zright,zcorr,rt;
  int lmax,zlfound,zrfound,l1,l2,dl,l,z1lim1,z1lim2,z1dz,cmode,i;

  int expect_gsh = 0;

  struct mod *model,out_model[2];
  static int lmin  = 0;		/* for correlations */
  static int extrapolate = 0;
  static int nsteps = 500;
  static int mode = 2;		/* 0: plot a few correlation curves 

                                   1: find correlation lengths for all
				   degrees 
				   
				   2: find correlation length for some degrees
				   
				*/
  rt = 0.7;			/* find the drop to which correlation? */	
  cmode = 1;			/* type of correlation */
  //rt = 1/exp(1.0);			
  
  switch(argc){
  case 2:
    ;
    break;
  default:
    fprintf(stderr,"%s model_file\n",argv[0]);
    break;
  }
  /* read in model */

  model = (struct mod *)calloc(1,sizeof(struct mod));

  read_she_model(argv[1],model,-1,1,expect_gsh);
  for(i=0;i < 2;i++){
    copy_model_par(model,(out_model+i));
    out_model[i].n = 1;
    allocate_model_coefficients((out_model+i));
  }

  lmax = model->lmax;
  
  switch(mode){
  case 0:			/* plot a few correlation curves at
				   different depths */
    zrange = 1500;
    dz = 2*zrange/(COMP_PRECISION)nsteps;
    for(z1=500;z1 <= 2500;z1 += 1000){
      interpolate_she_model(out_model,model,z1,lmax,extrapolate);
      
      zlim[0] = z1-zrange;zlim[1] = z1+zrange+1e-6;
      for(z2 = zlim[0];z2 <= zlim[1];z2 += dz){
	interpolate_she_model((out_model+1),model,z2,lmax,extrapolate);
	r = correlation(*out_model[0].a,*out_model[0].b,
			*out_model[1].a,*out_model[1].b,
			-lmax,0,lmin,cmode);
	fprintf(stdout,"%11g %11g %11g\n",z1,z2-z1,r);
      }
      printf("\n\n");
    }
    break;
  case 1:		
  case 2:
    if(mode == 1){
      z1lim1 = 50;z1lim2 = 2850+1e-6; z1dz = 100;/* loop through many depth layers */
      l1 = -lmax;l2=-lmax;dl=lmax;	/* output for all degrees */
    }else{
      l1 = 1; l2 = 31;dl=1;		/* output for some degrees */
      z1lim1 = 500;z1lim2 = 2500+1e-6 ; z1dz = 1000;; /* loop through few depth layers */
    }
    /* output for some degrees */
    /* find the correlation length to a certain drop  */

    dz = .1;
    for(z1=z1lim1;z1 <= z1lim2;z1 += z1dz){
      interpolate_she_model(out_model,model,z1,lmax,extrapolate);
      for(l=l1;l <= l2;l+=dl){	/* loop through degrees */

	/* 
	   search to the top 
	*/
	zleft = z1;r=1;
	while(finite(r) && (r > rt)){
	  zleft -= dz;
	  interpolate_she_model((out_model+1),model,zleft,lmax,extrapolate);
	  r = correlation(*out_model[0].a,*out_model[0].b,
			  *out_model[1].a,*out_model[1].b,
			  l,0,lmin,cmode);
	}
	if(!finite(r) && (zleft>=0))zlfound = 0;else zlfound = 1; /* could we find a correlation? */
	/* search to the bottom */
	zright = z1;r=1;
	while(finite(r) && (r > rt)){
	  zright += dz;
	  interpolate_she_model((out_model+1),model,zright,lmax,extrapolate);
	  r = correlation(*out_model[0].a,*out_model[0].b,
			  *out_model[1].a,*out_model[1].b,
			  l,0,lmin,cmode);
	}
	if(!finite(r) && (zright <= 2891))zrfound = 0;else zrfound = 1;
	/* pick either mean, left or right value, or set to nan if tr
	   is not reached within model */
	if(zlfound && zrfound)
	  zcorr = (-zleft+zright)/2.0;
	else if(zlfound && (!zrfound))
	  zcorr = -zleft+z1;
	else if(zrfound && (!zlfound))
	  zcorr = zright-z1;
	else
	  zcorr = my_make_nan();
	//printf("%11g\t%7.1f %i\t%7.1f %i\t%7.1f\n",z1,zleft-z1,zlfound,zright-z1,zrfound,zcorr);
	printf("%11g %3i %7.1f\n",z1,l,zcorr);
      }	/* end l loop */
      if(mode == 2)
	printf("\n\n");
    }

    break;
  }


  return 0;
}
