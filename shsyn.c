/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: shsyn.c,v 1.19 2004/03/10 18:05:19 becker Exp becker $ */
/* 
   routine to expand spherical harmonic coefficients into
   a field
   
   the coefficients are assumed to be in the physical normalization, that 
   is normalized as in Dahlen and Tromp

   (see abconvert.c for a conversion between geodetic and physical convention)

   uses plgndr to compute the associated Legendre polynomials by the stable recurrence
   formula (see e.g. Numerical Rec.) with included normalization factors

   function value at theta, phi is calculated by summing up l and m contributions
   or doing the m sum as an inverse FFT



*/
#include "shsyn.h"

#define FFT_CAST unsigned long

int main(int argc, char **argv )
{
  int lmax,i,j,l,m,nlon,nlon1,wrap_around,
    nlat,sum_mode,out_mode,rsize,lmax1,nn=0,os1,
    lmsize,rslmsize,nnt2,nnt4,m2,j2,rc,
    writeparameterfile=WPF,nr_of_coeff_sets;
  COMP_PRECISION a,b,dx,dy,dtheta,thetamin,thetamax,
    dincx,dincy,*y,min,max,avg,*r,xmin,xmax,
    ymin,ymax,one_over_l_fac,phi,
    *x,*sinarr,*cosarr,theta,value,p_array_size,
    sin_fac,fac,tmp,*tmparr,*ap,*bp,*phival,*dptheta,
    dummy,*p;
  BOOLEAN use_p_array,calcder,verbose,changed_bounds=FALSE,grd_out;
  GMT_PRECISION *gmtval;
  char grdfilename[2000]="";
  FILE *out;
  void *API;   
#ifdef BE_VERBOSE
  verbose=TRUE;
#else 
  verbose=FALSE;
#endif
  /* 
     handle command line arguments 
  */
  if((argc > 1) && (strcmp(argv[1],"-h")==0))
    argc=-9;
  /*
    default settings
  */
  sum_mode=DEF_SUM_MODE;
  out_mode=DEF_OUT_MODE;
  dincx=DEF_INC;
  dincy=-1;
  xmin=XMINS;
  xmax=XMAXS;
  ymin=YMINS;
  ymax=YMAXS;
  //
  switch(argc){
  case 1:{
    break;}
  case 2:{
    sscanf(argv[1],DATA_FSCAN_FORMAT,&dincx);
    break;}
  case 3:{
    sscanf(argv[1],DATA_FSCAN_FORMAT,&dincx);
    sscanf(argv[2],"%i",&sum_mode);
    break;}
  case 4:{
    sscanf(argv[1],DATA_FSCAN_FORMAT,&dincx);
    sscanf(argv[2],"%i",&sum_mode);
    check_out_mode(argv[3],&out_mode,grdfilename);
    break;}
  case 8:{
    sscanf(argv[1],DATA_FSCAN_FORMAT,&dincx);
    sscanf(argv[2],"%i",&sum_mode);
    check_out_mode(argv[3],&out_mode,grdfilename);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&xmin);
    sscanf(argv[5],DATA_FSCAN_FORMAT,&xmax);
    sscanf(argv[6],DATA_FSCAN_FORMAT,&ymin);
    sscanf(argv[7],DATA_FSCAN_FORMAT,&ymax);
    changed_bounds=TRUE;
    break;
  }
  case 9:{
    sscanf(argv[1],DATA_FSCAN_FORMAT,&dincx);
    sscanf(argv[2],"%i",&sum_mode);
    check_out_mode(argv[3],&out_mode,grdfilename);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&xmin);
    sscanf(argv[5],DATA_FSCAN_FORMAT,&xmax);
    sscanf(argv[6],DATA_FSCAN_FORMAT,&ymin);
    sscanf(argv[7],DATA_FSCAN_FORMAT,&ymax);
    changed_bounds=TRUE;
    sscanf(argv[8],DATA_FSCAN_FORMAT,&dincy);
    break;
  }
  default:{
    phelp(argv[0]);
    exit(-1);}
  }
  if(dincx <= 0.0){
    fprintf(stderr,"%s: spatial increment should be >= 0, dx: %g\n",
	    argv[0],dincx);
    exit(-1);
  }
  if(ymax < ymin){tmp=ymin;ymin=ymax;ymax=tmp;}
  if(xmax < xmin){tmp=xmin;xmin=xmax;xmax=tmp;}
  if((ymin< -90)||(ymax > 90)){
    fprintf(stderr,"%s: latitude should be between -90 and 90, ymin: %g ymax: %g\n",
	    argv[0],ymin,ymax);
    exit(-1);
  }
  if(changed_bounds && (sum_mode != SUM_OVER_L_AND_M)){
    if(verbose)
      fprintf(stderr,"%s: sum_mode reset to loops since limits changed from defaults\n",
	      argv[0]);
    sum_mode=SUM_OVER_L_AND_M;
  }
  if(dincy == -1)// same y spacing as x spacing
    dincy=dincx;
  // within limits?
  if((sum_mode < SUM_OVER_L_AND_M)  || 
     (sum_mode >  FFT_SUM) || 
     (out_mode < ASCII_STDOUT)|| 
     (out_mode > GRADIENT)){
    fprintf(stderr,"%s: sum_mode %i and out_mode %i out of limits\n",
	    argv[0],sum_mode,out_mode);
    phelp(argv[0]);
    exit(-1);
  }else{
    if(verbose)
      fprintf(stderr,"%s: sum_mode: %i and out_mode: %i\n",
	      argv[0],sum_mode,out_mode);
  }
  if((out_mode != TWO_GRDS) && 
     (out_mode != GRADIENT)){
    // no derivatives
    calcder=FALSE;
    nr_of_coeff_sets=1;
  }else{
    /* 
       have to calculate vector fields, 
       poloidal and tororidal part 
       or gradient of scalar fields
    */
    calcder=TRUE;
    if(out_mode != GRADIENT)
      nr_of_coeff_sets=2;
    else
      nr_of_coeff_sets=1;
    if(sum_mode != FFT_SUM){
      sum_mode=FFT_SUM;
      fprintf(stderr,"%s: switching to Fourier sum mode for vector fields\n",argv[0]);
    }
  }
  if((out_mode==ONE_GRD) || (out_mode==TWO_GRDS) || (out_mode==GRADIENT))
    grd_out = TRUE;
  else
    grd_out = FALSE;
  if(grd_out){
#ifdef USE_GMT4
    GMT_begin (argc, argv);
#else
    API = GMT_Create_Session (argv[0], 2U, 0, NULL);
#endif
  }

  /* 
     read in the maximum order of expansion 
  */
  rc=fscanf(IN_FILE,"%i",&lmax);
  if(!rc){
    fprintf(stderr,"%s: header read error\n",argv[0]);
    exit(-1);
  }
  if(lmax<0){ 
    phelp(argv[0]);
    exit(-1);
  }
  // limits and ranges
  lmax1=lmax+1;
  rslmsize=lmax1*(lmax1+1);
  lmsize=(int)((COMP_PRECISION)rslmsize/2.0);
  nlat=(int)((ymax-ymin)/dincy)+1;
  nlon=(int)((xmax-xmin)/dincx)+1;
  /* 
     for FFT, change the number of points in longitude to multiple of the FFT 2^n
     further, adjust the phi values
     
     since we are using vector harmonics, avoid poles for now and adjust latitudes
  */
  if(sum_mode == FFT_SUM){
    i=MAX((int)nlon/2,lmax1);
    nn=nextpwrtwo(i);
    nlon=nn*2;
    nlat=(nlat>=nlon/2)?(nlat):(nlon/2);
    if(verbose)
      fprintf(stderr,"%s: readjusting nlon, nlat and xrange for FFT\n",argv[0]);
    xmin=0.0;
    dx = 360.0/(COMP_PRECISION)(nlon);
    // the output will be from 0...360, really
    // even though we only calculate up to 360-dx
    xmax=xmin+dx*((COMP_PRECISION)(nlon-1));

    dy=(ymax-ymin)/(COMP_PRECISION)(nlat-1);
    ymax =  90.0-(((out_mode==TWO_GRDS)||(out_mode==GRADIENT))?(dy/2.0):(0.0));
    ymin = -90.0+(((out_mode==TWO_GRDS)||(out_mode==GRADIENT))?(dy/2.0):(0.0));

    wrap_around=1;
  }else{
    // no changes necessary
    dx = dincx;
    dy = dincy;
    wrap_around=0;
  }
  if(verbose){
    fprintf(stderr,"%s: lmax=%i nlon=%i nlat=%i\n",
	    argv[0],lmax,nlon+wrap_around,nlat);
    fprintf(stderr,"%s: E=%g dx=%g W=%g S=%g dy=%g N=%g\n",
	    argv[0],xmin,dx,xmax+wrap_around*(COMP_PRECISION)dx,
	    ymin,dy,ymax);
  }
  if(writeparameterfile){// write size of field to file
    out=fopen("expansion.par","w");
    fprintf(out,"%20.10lf %20.10lf %20.10lf %20.10lf %20.10lf %20.10lf %i %i\n",
	    xmin,dx,xmax+wrap_around*dx,
	    ymin,dy,ymax,nlon+wrap_around,nlat);
    fclose(out);
  }
  /* 
     lookup tables for Legendre function values 
     and coordinates 
  */ 
  if(verbose)
    fprintf(stderr,"%s: constructing lookup tables (Legendre functions) ...\n",argv[0]);
  /* 
     y coordinates, evenly spaced in z-units, 
     i.e. cos(theta) 
  */
  if((y=(COMP_PRECISION *)
      calloc(nlat,sizeof(COMP_PRECISION)))==NULL)
    MEMERROR;
  thetamin= LATITUDE2THETA(ymin); 
  thetamax= LATITUDE2THETA(ymax);
  dtheta=(thetamax-thetamin)/
    ((COMP_PRECISION)(nlat-1));
  for(theta=thetamin,i=0;i<nlat;i++,theta += dtheta)
    y[i] = cos(theta);
  
  // decide if we want to keep all Legendre functions in memory
  p_array_size=(double)(lmsize*nlat*sizeof(COMP_PRECISION))/ONE_MEGABYTE;
  if(p_array_size < P_MEM_LIMIT){
    if((p=(COMP_PRECISION *)
	calloc(lmsize*nlat,sizeof(COMP_PRECISION)))==NULL){
      fprintf(stderr,"%s: memerror for p array: %g MB, limit is set to %g\n",
	      argv[0],p_array_size,P_MEM_LIMIT);exit(-1);
    }
    else{
      if(verbose)
	fprintf(stderr,"%s: using %6g MB for P array, limit is set to %g\n",
		argv[0],p_array_size,P_MEM_LIMIT);
    }
    use_p_array=TRUE;
    plgndr(y,nlat,p,lmax);
  }else{
    if(verbose){
	fprintf(stderr,"%s: since P array with %6g MB is higher than limit of %g\n",
		argv[0],p_array_size,P_MEM_LIMIT);
	fprintf(stderr,"%s: resorting to (slow) method of individual Legendre functions\n",
		argv[0]);
    }
    use_p_array=FALSE;
    p = NULL;
  }
  if(calcder){// calculate derivatives
    if(!use_p_array){
      fprintf(stderr,"%s: attempting to calculate Legendre derivatives without P array\n",
	      argv[0]);
      fprintf(stderr,"%s: not implemented yet\n",argv[0]);
      exit(-1);
    }
    if(verbose)
      fprintf(stderr,"%s: constructing their derivatives....\n",argv[0]);
    if((dptheta=(COMP_PRECISION *)calloc(lmsize*nlat,sizeof(COMP_PRECISION)))==NULL)
      {fprintf(stderr,"memerror\n");exit(-1);}
    else
      if(verbose)
	fprintf(stderr,"%s: using %6g MB for DPTHETA array\n",
		argv[0],(float)lmsize*(float)sizeof(COMP_PRECISION)*(float)nlat/ONE_MEGABYTE);
    /* 
       create the dptheta array which will 
       hold d_theta X_lm 
    */
    pdtheta_lgndr(y,nlat,p,dptheta,lmax,&dummy,FALSE);
    /* 
       overwrite P with d_phi X_lm which is 
       m/sin(theta) X_lm WITH A AND B FLIPPED 
       FOR SYNTHESIS 
       
    */
    if(verbose)
      fprintf(stderr,"%s: changing P and DPTHETA\n",argv[0]);
     if(out_mode != GRADIENT){
       /* 
	  calculate derivatives plus 1/sqrt(l(l+1)) 
	  factor for spheroidal/poloidal expansion 
       */
       for(j=0;j<nlat;j++){
	 // sin(theta)
	 sin_fac = sqrt((1.0+EPS_SMALL)-y[j]*y[j]);
	 for(l=0;l<=lmax;l++){
	   one_over_l_fac=((l==0)?(1.0):
		       (1.0/sqrt((COMP_PRECISION)(l*(l+1)))));
	   P(l,0,j) = 0.0; 
	   DPTHETA(l,0,j) *= one_over_l_fac;
	   for(m=1;m<=l;m++){
	     // d Y_lm/d phi
	     P(l,m,j)       *= (COMP_PRECISION)m;
	     P(l,m,j)       /= sin_fac;
	     // scaling for d Y_lm/d phi 
	     // and d Y_lm/d theta
	     P(l,m,j)       *= one_over_l_fac;
	     DPTHETA(l,m,j) *= one_over_l_fac;
	   }
	 }
       }
     }else{
       /* 
	  need the derivatives for the calculation 
	  of the gradient of a scalar field only, 
	  hence do not include the 1/sqrt(l(l+1)) 
	  factor in the Ylm terms
       */
       for(j=0;j<nlat;j++){
	 sin_fac = sqrt((1.0+EPS_SMALL)-y[j]*y[j]);
	 for(l=0;l<=lmax;l++){
	   P(l,0,j) = 0.0; 
	   for(m=1;m<=l;m++){
	     P(l,m,j)       *= (COMP_PRECISION)m;
	     P(l,m,j)       /= sin_fac;
	   }
	 }
       }
     }
  }else{
    dptheta = NULL;
  }
  if(sum_mode == SUM_OVER_L_AND_M){// sum up all contributions with loops
    if( out_mode != TWO_GRDS ){
      /* for non FFT and scalar field
	 
	 create R array with coefficients
	 
	 TP(l, m, theta, size) = A_{lm} * P_{lm}(theta)
	 and
	 UP(l, m, theta, size) = B_{lm} * P_{lm}(theta)
	 
	 for faster summing
	 
      */
      rsize = nlat * rslmsize;
      if(verbose)
	fprintf(stderr,"%s: constructing R, S and trig. factors...\n",argv[0]);
      
      if((r=(COMP_PRECISION *)calloc(rsize,sizeof(COMP_PRECISION)))==NULL)
	{fprintf(stderr,"memerror, rsize=%iMB\n",
		 (int)(rsize*(COMP_PRECISION)sizeof(COMP_PRECISION)/
		 ONE_MEGABYTE));exit(-1);}
      else {
	if(verbose)fprintf(stderr,"%s: allocated %10g MB for R\n",
			   argv[0],(double)(rsize*sizeof(COMP_PRECISION)/ONE_MEGABYTE));}
      if(use_p_array){ // Legendre in loopup table
	for(l=0;l<=lmax;l++)
	  for(m=0;m<=l;m++){
	    if((fscanf(IN_FILE,TWO_DATA_FSCAN_FORMAT,&a,&b))!=2){
	      fprintf(stderr,"%s: AB coefficient read error, l=%i m=%i\n",
		      argv[0],l,m);exit(-1);}
	    else{
	      for(i=0;i<nlat;i++){
		TP(l, m, i)= a * P(l, m, i);
		UP(l, m, i)= b * P(l, m, i);
	      }
	    }
	  }
      }else{
	for(l=0;l<=lmax;l++)
	  for(m=0;m<=l;m++){
	    if((fscanf(IN_FILE,TWO_DATA_FSCAN_FORMAT,&a,&b))!=2){
	      fprintf(stderr,"%s: AB coefficient read error, l=%i m=%i\n",
		      argv[0],l,m);exit(-1);}
	    else{
	      for(i=0;i<nlat;i++){
		TP(l, m, i)= a * 
		  (tmp=slgndr(l,m,y[i]));
		UP(l, m, i)= b * tmp;
	      }
	    }
	  }
      }
    }else{
      r = NULL;
    }
    /* 
       x coordinates 
    */
    if((x=(COMP_PRECISION *)calloc(nlon,sizeof(COMP_PRECISION)))==NULL)
      {fprintf(stderr,"%s: memerror for x coordinates\n",argv[0]);exit(-1);};
    for(phi=LONGITUDE2PHI(xmin),j=0,tmp=DEG2RAD(dx); 
	j < nlon ; 
	j++,phi += tmp){
      x[j]=phi;
    }
    /* 
       cos and sin factors for R l m  and S l m 
    */
    if((cosarr=(COMP_PRECISION *)calloc(nlon*lmax1,sizeof(COMP_PRECISION)))==NULL)
      {fprintf(stderr,"memerror cosarr\n");exit(-1);};
    if((sinarr=(COMP_PRECISION *)calloc(nlon*lmax1,sizeof(COMP_PRECISION)))==NULL)
      {fprintf(stderr,"memerror sinarr\n");exit(-1);};
    for(j=os1=0; j < nlon ;j++,os1+=lmax1){
      *(cosarr+os1)=1.0;
      *(sinarr+os1)=0.0;}
    for(m=1;m<=lmax;m++){
      for(j=os1=0; j < nlon ;j++,os1+=lmax1){
	/* leave it like that, otherwise had 
	   to care abut sign */
	*(cosarr+ os1 + m)=
	  cos(((COMP_PRECISION)m)*x[j]);
	*(sinarr+ os1 + m)=
	  sin(((COMP_PRECISION)m)*x[j]); 
      }
    }
  }else{
    cosarr = sinarr = r = NULL;
  }
  
  switch(sum_mode){
    /* 
       
       sum formula  
       
    */
  case SUM_OVER_L_AND_M:{
    switch(out_mode){
    case ASCII_STDOUT:{/* ASCII lines output */
      /* loop over all locations */
      for(min=FLT_MAX,max=FLT_MIN,avg=0.0,theta=ymin,i=0;
	  i<nlat;
	  i++,theta+=dy){
	for(j=os1=0;j<nlon;j++,os1+=lmax1){
	  /* sum contributions from l */
	  value=0.0;
	  for(l=0;l<=lmax;l++){
	    m=0;
	    value  += TP(l, m, i)*  *(cosarr+os1+m);
	    /* sum contributions from m */
	    for(fac=0.0,m=1;m<=l;m++){
	      fac  += TP(l, m, i)*  *(cosarr+os1+m);
	      fac  += UP(l, m, i)*  *(sinarr+os1+m);
	    }
	    value += SQRT_TWO*fac;
	  }
	  fprintf(OUT_FILE,ASCII_DATA_FORMAT,value); 
	  fprintf(OUT_FILE,"\n");
	  if(value < min)
	    min=value;
	  if(value > max)
	    max=value;
	  avg += value;
	}
	if(verbose)
	  fprintf(stderr,"%s: lat = %10g\r",argv[0],theta);
      }
      if(verbose)
      	fprintf(stderr,"\n%s: min=%10g avg=%10g(avg*4pi=%10g) max=%10g\n",
		argv[0],min,avg/((COMP_PRECISION)(nlon*nlat)),
		avg/((COMP_PRECISION)(nlon*nlat))*2.0*TWOPI,max);
      break;
    }
    case ONE_GRD:
    case BINARY_STDOUT:{/* grdfile or binary block outputs */
      /* sum formula */
      if((gmtval=(GMT_PRECISION *)malloc(sizeof(GMT_PRECISION)*nlon*nlat))==NULL){
	fprintf(stderr,"%s: memerror for gmtval \n",argv[0]);exit(-1);}
      if(verbose)
	fprintf(stderr,"%s: creating the data array for block output\n",argv[0]);
      /* loop over all locations */
      for(i=0,theta=ymin;i<nlat;i++,theta+=dy){
	if(verbose)fprintf(stderr,"%s: %11g\r",argv[0],theta);
	for(os1=j=0;j<nlon;j++,os1+=lmax1){
	  /* sum contributions from l */
	  value=0.0;
	  for(l=0;l<=lmax;l++){
	    m=0;
	    value  += TP(l, m, i)*  *(cosarr+os1+m);
	    /* sum contributions from m */
	    for(fac=0.0,m=1;m<=l;m++){
	      fac  += TP(l, m, i)*  *(cosarr+os1+m);
	      fac  += UP(l, m, i)*  *(sinarr+os1+m);
	    }
	    value += SQRT_TWO*fac;
	  }
	  *(gmtval+(nlat-i-1)*nlon+j)=(GMT_PRECISION)value;
	}
      }
      grid_output(out_mode,grdfilename,gmtval,nlon,nlat,
		  xmin,xmax,ymin,ymax,dx,dy,
		  argc,argv,lmax,verbose,&API);
      break;
    }
    default:{
      fprintf(stderr,"inconsistency, out_mode %i undefined for sum rule\n",
	      out_mode);
      exit(-1);
      break;
    }}
    break;
  }// end of sum over l and m
  /* 
     begin FFT for m 
  */
  case FFT_SUM:{
    // add one for wrapping phi around (0...360)
    nlon1=nlon+wrap_around;
    nnt2=nn*2;
    nnt4=nnt2*2;
    //
    switch(out_mode){
    case ASCII_STDOUT:{/* ASCII line output */
      /* 
	 read the coefficients 
      */
      if((ap=(COMP_PRECISION *)malloc(nr_of_coeff_sets*lmsize*sizeof(COMP_PRECISION)))==NULL ||
	 (bp=(COMP_PRECISION *)malloc(nr_of_coeff_sets*lmsize*sizeof(COMP_PRECISION)))==NULL){
	fprintf(stderr,"%s: A B memerror, lmax=%i lmsize=%i\n",argv[0],lmax,lmsize); exit(-1);}
      read_ab(ap,bp,lmax, argv[0],IN_FILE,nr_of_coeff_sets,lmsize,verbose);
      // FFT array
      if((tmparr=(COMP_PRECISION *)malloc(nnt4*sizeof(COMP_PRECISION)))==NULL)
	{fprintf(stderr,"\nmemerror %i\n\n",(int)(nnt4*sizeof(COMP_PRECISION)));exit(-1);}
      /* function values along phi */
      if((phival=(COMP_PRECISION *)malloc(nlon1*sizeof(COMP_PRECISION)))==NULL)
	{fprintf(stderr,"\nmemerror phival one line %i\n\n",
		 (int)(nlon1*sizeof(COMP_PRECISION)));exit(-1);}
      for(min=FLT_MAX,max=FLT_MIN,avg=0.0,theta=ymin,i=0;
	  i<nlat;
	  i++,theta+=dy){ /* loop over all latitudes */
	if(verbose)
	  fprintf(stderr,"%s: Using FFT line by line, lat=%11g   \r",
		  argv[0],theta);
	for(j=0;j<nlon1;j++)
	  phival[j]= 0.0;
	if(use_p_array){
	  for(l=0;l<=lmax;l++){ /* sum up all contributions from different ls */
	    /* 
	       initialize array for FFT 
	       first, constant m=0 term
	    */
	    /* m=0 */
	    *(tmparr) =  *(tmparr+nnt4-1) = AP(l, 0) * P(l, 0, i); 
	    /* now higher frequencies, m != 0 */
	    for(m=1,m2=2;m<=l;m++,m2+=2){
	      tmp = FOURIER_FACTOR *  P(l, m, i);
	      *(tmparr+m2  ) = AP(l, m) * tmp; 
	      *(tmparr+m2+1) = BP(l, m) * tmp;
	      /* complex conjugates */
	      *(tmparr+nnt4-m2  ) =   *(tmparr+m2  );
	      *(tmparr+nnt4-m2+1) = - *(tmparr+m2+1);
	    }
	    /* set the higher frequencies to zero (have to do that each time
	       since the array is overwritten by the inverse FFT with function 
	       values*/
	    for(m2=2*l+2;m2<=nnt2;m2+=2){
	      *(tmparr+m2  )=0.0;
	      *(tmparr+m2+1)=0.0;
	      *(tmparr+nnt4-m2  )=0.0;
	      *(tmparr+nnt4-m2+1)= 0.0;
	    }
	    /* do inverse FFT  to get function values at different phis */
	    dfour1(tmparr-1,(FFT_CAST)(nnt2),-1); 
	    /* assign these to the phival array for all l contributions */
	    for(j=0;j<nlon;j++){
	      phival[j] += *(tmparr+j*2);
	    }
	  }// end l loop
	}else{/* have to calculate Plmi */
	  for(l=0;l<=lmax;l++){ // 
	    *(tmparr) =  *(tmparr+nnt4-1) = AP(l, 0) * slgndr(l,0,y[i]); 
	    for(m=1,m2=2;m<=l;m++,m2+=2){
	      tmp = FOURIER_FACTOR *  slgndr(l,m,y[i]);
	      *(tmparr+m2  ) = AP(l, m) * tmp; 
	      *(tmparr+m2+1) = BP(l, m) * tmp;
	      *(tmparr+nnt4-m2  ) =   *(tmparr+m2  );
	      *(tmparr+nnt4-m2+1) = - *(tmparr+m2+1);
	    }
	    for(m2=2*l+2;m2<=nnt2;m2+=2)
	      *(tmparr+m2  )=  *(tmparr+m2+1)= 
		*(tmparr+nnt4-m2  )=  *(tmparr+nnt4-m2+1)= 0.0;
	    dfour1(tmparr-1,(FFT_CAST)(nnt2),-1); 
	    for(j=0;j<nlon;j++) 
	      phival[j] += *(tmparr+j*2);
	  }
	}
	// add end node
	if(wrap_around)
	  phival[nlon] = phival[0];
	/* print row  */
	for(j=0;j<nlon1;j++){
	  fprintf(OUT_FILE,ASCII_DATA_FORMAT,phival[j]);
	  fprintf(OUT_FILE,"\n");
	  if(phival[j] < min)
	    min= phival[j];
	  if(phival[j] > max)
	    max= phival[j];
	  avg += phival[j];
	}
	
      }// end latitude (y) loop
      free(tmparr);free(phival);
      if(verbose){
	fprintf(stderr,"\n%s: min= %10g avg= %10g (avg*4pi=%10g) max= %10g\n",
		argv[0],min,avg/((COMP_PRECISION)(nlon1*nlat)),
		avg/((COMP_PRECISION)(nlon*nlat))*TWOPI*2.0,max);
      }
      break;
    }// end ASCII stdout
    case ONE_GRD:;
    case BINARY_STDOUT:{
      /* 
	 grdfile or binary 
	 block outputs with FFT 
      */
      /* read the coefficients */
      if((ap=(COMP_PRECISION *)malloc((size_t)nr_of_coeff_sets*lmsize*sizeof(COMP_PRECISION)))==NULL ||
	 (bp=(COMP_PRECISION *)malloc((size_t)nr_of_coeff_sets*lmsize*sizeof(COMP_PRECISION)))==NULL){
	fprintf(stderr,"%s: A B memerror, lmax=%i lmsize=%i\n",argv[0],lmax,lmsize); exit(-1);}
      read_ab(ap,bp,lmax, argv[0],IN_FILE,nr_of_coeff_sets,lmsize,verbose);
      // FFT array
      if((tmparr=(COMP_PRECISION *)malloc(nnt4*sizeof(COMP_PRECISION)))==NULL)
	{fprintf(stderr,"\nmemerror %i tmparr\n\n",
		 (int)(nnt4*sizeof(COMP_PRECISION)));exit(-1);}
      /* function values along phi and theta */
      if((gmtval=(GMT_PRECISION *)calloc((size_t)nr_of_coeff_sets*nlon1*nlat,sizeof(GMT_PRECISION)))==NULL)
	{fprintf(stderr,"\nmemerror gmtval %i\n\n",
		 (int)(nr_of_coeff_sets*nlon1*nlat*sizeof(GMT_PRECISION)));exit(-1);}
      for(theta=ymin,
	    i=0;i<nlat;i++,
	    theta+=dy){ /* loop over all latitudes, y[i] */
	if(verbose)
	  fprintf(stderr,"%s: using FFT block method, lat= %11g\r",
		  argv[0],theta);
	if(use_p_array){
	  for(l=0;l<=lmax;l++){ 
	    /* sum up all contributions from different ls */
	    /* initialize array for FFT 
	       first, constant m=0 term
	    */ 
	    /*  m=0 */
	    *(tmparr) =  *(tmparr+nnt4-1) = AP(l, 0) * P(l, 0, i); 
	    /* now higher frequencies, m != 0 */
	    for(m=1,m2=2;m<=l;m++,m2+=2){
	      tmp =  FOURIER_FACTOR * P(l, m, i);
	      *(tmparr+m2  ) = AP(l, m) * tmp;
	      *(tmparr+m2+1) = BP(l, m) * tmp;
	      /* complex conjugates */
	      *(tmparr+nnt4-m2  ) =   *(tmparr+m2  );
	      *(tmparr+nnt4-m2+1) = - *(tmparr+m2+1);
	    }
	    /* set the higher frequencies to zero (have to do that each time
	       since the array is overwritten by the inverse FFT with function 
	       values*/
	    for(m2=2*l+2;m2<=nnt2;m2+=2)
	      *(tmparr+m2  )=  *(tmparr+m2+1)= 
		*(tmparr+nnt4-m2  )=  *(tmparr+nnt4-m2+1)= 0.0;
	    /* do inverse FFT  to get function values at different phis */
	    dfour1(tmparr-1,(FFT_CAST)(nnt2),-1); 
	    /* assign these to the gmtval array for all l contributions */
	    for(os1=(nlat-i-1)*nlon1,j=0;j<nlon;j++) 
	      *(gmtval+os1+j) += (GMT_PRECISION) *(tmparr+j*2);
	    if(wrap_around)
	      *(gmtval+os1+nlon) += (GMT_PRECISION) *(tmparr);
	  }// end l loop
	}else{
	  for(l=0;l<=lmax;l++){ 
	    *(tmparr) =  *(tmparr+nnt4-1) = AP(l, 0) * slgndr(l,0,y[i]); 
	    for(m=1,m2=2;m<=l;m++,m2+=2){
	      tmp =  FOURIER_FACTOR * slgndr(l,m,y[i]);
	      *(tmparr+m2  ) = AP(l, m) * tmp;
	      *(tmparr+m2+1) = BP(l, m) * tmp;
	      *(tmparr+nnt4-m2  ) =   *(tmparr+m2  );
	      *(tmparr+nnt4-m2+1) = - *(tmparr+m2+1);
	    }
	    for(m2=2*l+2;m2<=nnt2;m2+=2)
	      *(tmparr+m2  )=  *(tmparr+m2+1)= 
		*(tmparr+nnt4-m2  )=  *(tmparr+nnt4-m2+1)= 0.0;
	    dfour1(tmparr-1,(FFT_CAST)(nnt2),-1); 
	    for(os1=(nlat-i-1)*nlon1,j2=j=0;j<nlon;j++,j2+=2) 
	      *(gmtval+os1+j) += (GMT_PRECISION) *(tmparr+j2);
	    if(wrap_around)
	      *(gmtval+os1+nlon) += (GMT_PRECISION) *(tmparr);
	  }
	}
      }// end theta loop
      if(verbose)
	fprintf(stderr,"\n");
      free(tmparr);free(ap);free(bp);
      grid_output(out_mode,grdfilename,gmtval,nlon1,nlat,xmin,xmax+wrap_around*dx,
		  ymin,ymax,dx,dy,argc,argv,lmax,verbose,&API);
      free(gmtval);
      break;
    }// end binary stdout and one_grd
    /* 
       grdfile with FFT for vectors 
    */
    case TWO_GRDS:{
      /* read the coefficients, first set is polodial, second toroidal */
      if((ap=(COMP_PRECISION *)malloc(nr_of_coeff_sets*lmsize*sizeof(COMP_PRECISION)))==NULL ||
	 (bp=(COMP_PRECISION *)malloc(nr_of_coeff_sets*lmsize*sizeof(COMP_PRECISION)))==NULL){
	fprintf(stderr,"%s: A B memerror, lmax=%i lmsize=%i\n",argv[0],lmax,lmsize); exit(-1);}
      read_ab(ap,bp,lmax, argv[0],IN_FILE,nr_of_coeff_sets,lmsize,verbose);
      if((tmparr=(COMP_PRECISION *)malloc(nnt4*sizeof(COMP_PRECISION)))==NULL)
	{fprintf(stderr,"\nmemerror %i\n\n",(int)(nnt4*sizeof(COMP_PRECISION)));exit(-1);}
      /* 
	 function values along phi and theta
	 the first set will hold v_x, the second v_y
       */
      if((gmtval=(GMT_PRECISION *)calloc(nr_of_coeff_sets*nlon1*nlat,sizeof(GMT_PRECISION)))==NULL)
	{fprintf(stderr,"\nmemerror %i\n\n",(int)(nr_of_coeff_sets*nlon1*nlat*sizeof(GMT_PRECISION)));exit(-1);}
      if(verbose)
	fprintf(stderr,"%s: using Dahlen and Tromp  convention\n",argv[0]);
      for(theta=ymin,i=0;i<nlat;i++,theta+=dy){ // begin y loop
	if(verbose)
	  fprintf(stderr,"%s: using FFT block method for velocities, lat=%11g\r",
		  argv[0],theta);
	for(l=0;l<=lmax;l++){ /* l loop */
	  /* 
	     (mind: the 1/sqrt(l(l+1)) factor is already included in the derivatives)

	     u_phi = -(BP DPHI + AT DPTHETA) cos m phi + ( -AP DPHI + BT DPTHETA) sin m phi
	  */
	  //m=0;
	  *(tmparr) =  *(tmparr+nnt4-1) = -AT(l,0) * DPTHETA(l,0,i); 
	  for(m=1,m2=2;m<=l;m++,m2+=2){
	    tmp = FOURIER_FACTOR * P(l, m, i);
	    *(tmparr+    m2  ) =    ( BP(l,m) * P(l,m,i) - AT(l,m) * DPTHETA(l,m,i))*FOURIER_FACTOR; /* cosine part */
	    *(tmparr+    m2+1) =    (-AP(l,m) * P(l,m,i) - BT(l,m) * DPTHETA(l,m,i))*FOURIER_FACTOR; /* sine part */
	    *(tmparr+nnt4-m2  ) =   *(tmparr+m2  );
	    *(tmparr+nnt4-m2+1) = - *(tmparr+m2+1);
	  }
	  for(m2=2*l+2;m2<=nnt2;m2+=2)
	    *(tmparr+m2  )=  *(tmparr+m2+1)= 
	      *(tmparr+nnt4-m2  )=  *(tmparr+nnt4-m2+1)= 0.0;
	  dfour1(tmparr-1,(FFT_CAST)(nnt2),-1); 
	  for(os1=(nlat-i-1)*nlon1,j=0,j2=0;j<nlon;j++,j2+=2) 
	    *(gmtval+os1+j) += (GMT_PRECISION) *(tmparr+j2);
	  if(wrap_around)
	    *(gmtval+os1+nlon) += (GMT_PRECISION) *(tmparr);
	  /* 
	     u_theta = (AP DPTHETA - BT P ) cos (m phi) - (BP DPTHETA + AT P) sin ( m phi )
	  */
	  //m=0;
	  *(tmparr) =  *(tmparr+nnt4-1) = AP(l,0) * DPTHETA(l,0,i); 
	  for(m=1,m2=2;m<=l;m++,m2+=2){
	    *(tmparr+    m2  ) =   (AP(l,m) * DPTHETA(l,m,i) + BT(l,m) * P(l,m,i)) *  FOURIER_FACTOR;/* cosine part */
	    *(tmparr+    m2+1) =   (BP(l,m) * DPTHETA(l,m,i) - AT(l,m) * P(l,m,i)) *  FOURIER_FACTOR;/* sine part */
	    *(tmparr+nnt4-m2  ) =   *(tmparr+m2  );
	    *(tmparr+nnt4-m2+1) = - *(tmparr+m2+1);
	  }
	  for(m2=2*l+2;m2<=nnt2;m2+=2)
	    *(tmparr+m2  )=  *(tmparr+m2+1)= 
	      *(tmparr+nnt4-m2  )=  *(tmparr+nnt4-m2+1)= 0.0;
	  dfour1(tmparr-1,(FFT_CAST)(nnt2),-1); 
	  for(os1=(nlon1*nlat)+(nlat-i-1)*nlon1,j=j2=0;j<nlon;j++,j2+=2) 
	    *(gmtval+os1+j) += (GMT_PRECISION) *(tmparr+j2);
	  if(wrap_around)
	    *(gmtval+os1+nlon) += (GMT_PRECISION) *(tmparr);
	}
      }
      if(verbose)fprintf(stderr,"\n");
      if((nr_of_coeff_sets != 2)&&(out_mode != GRADIENT)){
	fprintf(stderr,"%s: output mode inconsistency\n",argv[0]);
	exit(-1);
      }
      free(tmparr);free(ap);free(bp);
      my_gmt_write_grd(gmtval, verbose, argc,argv,
		       FIRST_VEL_OUT, nlon1,nlat, 
		       xmin, xmax+dx*wrap_around,
		       ymin, ymax, dx,  dy,&API);
      my_gmt_write_grd((gmtval+nlon1*nlat), verbose, argc,argv,
		       SECOND_VEL_OUT, nlon1,nlat, xmin, 
		       xmax+dx*wrap_around,
		       ymin, ymax, dx,  dy,&API);
      free(gmtval);
      break;
    }// end vector field out
    /* 
       calculate the gradient of the expanded field
       and write binary to two grid files 
    */
    case GRADIENT:{
      /* read the coefficients */
      if((ap=(COMP_PRECISION *)malloc(nr_of_coeff_sets*lmsize*sizeof(COMP_PRECISION)))==NULL ||
	 (bp=(COMP_PRECISION *)malloc(nr_of_coeff_sets*lmsize*sizeof(COMP_PRECISION)))==NULL){
	fprintf(stderr,"%s: A B memerror, lmax=%i lmsize=%i\n",argv[0],lmax,lmsize); exit(-1);}
      read_ab(ap,bp,lmax, argv[0],IN_FILE,nr_of_coeff_sets,lmsize,verbose);
      if((tmparr=(COMP_PRECISION *)malloc((int)nnt4*sizeof(COMP_PRECISION)))==NULL)
	{fprintf(stderr,"\nmemerror %i\n\n",(int)(nnt4*sizeof(COMP_PRECISION)));exit(-1);}
      /* 
	 function values along phi and theta, two arrays for 1/sin(t) d/dp and d/dt
      */
      if((gmtval=(GMT_PRECISION *)calloc(2*nlon1*nlat,sizeof(GMT_PRECISION)))==NULL)
	{fprintf(stderr,"\nmemerror %i\n\n",(int)(2*nlon1*nlat*sizeof(GMT_PRECISION)));exit(-1);}
      for(theta=ymin,i=0;i<nlat;i++,theta+=dy){ 
	if(verbose)
	  fprintf(stderr,"%s: using FFT block method for gradients, lat= %11g\r",
		  argv[0],theta);
	for(l=0;l<=lmax;l++){ /* l loop */
	  /* 
	     grad_phi f = 1/sin(t) df/d_phi = 
	     
	     -BP DPHI cos m phi +  -AP DPHI  sin m phi

	  */
	  //m=0;
	  *(tmparr) =  *(tmparr+nnt4-1) = 0.0; 
	  for(m=1,m2=2;m<=l;m++,m2+=2){
	    tmp = FOURIER_FACTOR * P(l, m, i);
	    *(tmparr+    m2  ) =    BP(l,m) * P(l,m,i) * FOURIER_FACTOR; /* cosine part */
	    *(tmparr+    m2+1) =   -AP(l,m) * P(l,m,i) * FOURIER_FACTOR; /* sine part */
	    *(tmparr+nnt4-m2  ) =   *(tmparr+m2  );
	    *(tmparr+nnt4-m2+1) = - *(tmparr+m2+1);
	  }
	  for(m2=2*l+2;m2<=nnt2;m2+=2)
	    *(tmparr+m2  )=  *(tmparr+m2+1)= 
	      *(tmparr+nnt4-m2  )=  *(tmparr+nnt4-m2+1)= 0.0;
	  dfour1(tmparr-1,(FFT_CAST)(nnt2),-1); 
	  for(os1=(nlat-i-1)*nlon1,j=j2=0;j<nlon;j++,j2+=2) 
	    *(gmtval+os1+j) += (GMT_PRECISION) *(tmparr+j2);
	  if(wrap_around)
	    *(gmtval+os1+nlon) += (GMT_PRECISION) *(tmparr);
	  /* 
	     grad_theta f  = df/d_theta = 
	     AP DPTHETA cos (m phi) - BP DPTHETA  sin ( m phi )
	  */
	  //m=0;
	  *(tmparr) =  *(tmparr+nnt4-1) = AP(l,0) * DPTHETA(l,0,i); 
	  for(m=1,m2=2;m<=l;m++,m2+=2){
	    *(tmparr+    m2  ) =   AP(l,m) * DPTHETA(l,m,i) *  FOURIER_FACTOR;/* cosine part */
	    *(tmparr+    m2+1) =   BP(l,m) * DPTHETA(l,m,i) *  FOURIER_FACTOR;/* sine part */
	    *(tmparr+nnt4-m2  ) =   *(tmparr+m2  );
	    *(tmparr+nnt4-m2+1) = - *(tmparr+m2+1);
	  }
	  for(m2=2*l+2;m2<=nnt2;m2+=2)
	    *(tmparr+m2  )=  *(tmparr+m2+1)= 
	      *(tmparr+nnt4-m2  )=  *(tmparr+nnt4-m2+1)= 0.0;
	  dfour1(tmparr-1,(FFT_CAST)(nnt2),-1); 
	  for(os1=nlon1*nlat+(nlat-i-1)*nlon1,j=j2=0;j<nlon;j++,j2+=2) 
	    *(gmtval+os1+j) += (GMT_PRECISION) *(tmparr+j2);
	  if(wrap_around)
	    *(gmtval+os1+nlon) += (GMT_PRECISION) *(tmparr);
	}
      }
      if(verbose)fprintf(stderr,"\n");
      if((nr_of_coeff_sets != 2)&&(out_mode != GRADIENT)){
	fprintf(stderr,"%s: output mode inconsistency in gradient\n",argv[0]);
	exit(-1);
      }
      free(tmparr);free(ap);free(bp);
      my_gmt_write_grd(gmtval, verbose, argc,argv,
		       FIRST_VEL_OUT, nlon1,nlat, 
		       xmin, xmax+dx*wrap_around,
		       ymin, ymax, dx,  dy,&API);
      my_gmt_write_grd((gmtval+nlon1*nlat), verbose,
		       argc,argv,
		       SECOND_VEL_OUT, nlon1,nlat, 
		       xmin, xmax+dx*wrap_around,
		       ymin, ymax, dx,  dy,&API);
      free(gmtval);
      break;
    }
    default:{
      fprintf(stderr,"inconsistency, out_mode %i undefined \n",out_mode);
      exit(-1);
      break;
    }
    }
    break;
  }
  default:{
    fprintf(stderr,"inconsistency, sum_mode %i undefined \n",sum_mode);
    exit(-1);
    break;
  }
  }
  
  if(verbose)
    fprintf(stderr,"%s: Done.\n",argv[0]);
  if(grd_out){
#ifdef USE_GMT4
    GMT_end (argc, argv);
#else
    GMT_Destroy_Session (API);
#endif
  }
  return 0;
}


int nextpwrtwo(int i)
{
  COMP_PRECISION y;
  y = ceil(log((COMP_PRECISION)i)*1.44269504088896341);
  return (int)pow(2,y);
}
void phelp(char *name)
{
  fprintf(stderr,"%s [dincx(%g)] [sum_mode(%i)] [out_mode(%i)] [ xmin(%g) xmax(%g) ymin(%g) ymax(%g) ] [dincy (dincx)]\n\n",
	  name,DEF_INC,DEF_SUM_MODE,DEF_OUT_MODE,
	  XMINS,XMAXS,YMINS,YMAXS);
  fprintf(stderr,"\tExpands spherical harmonics coefficients into spatial domain, the\n");
  fprintf(stderr,"\tinput is Thorsten's ASCII AB format as produced by \"shana\"\n");
  fprintf(stderr,"\t(fully normalized physical convention as in Dahlen and Tromp:\n");
  fprintf(stderr,"\tl_max\n");
  fprintf(stderr,"\tA_00 B_00\n\tA_10 B_10\n\tA_11 B_11\n\t...)\n");
  fprintf(stderr,"\t%s was compiled with %s precision.\n",name,COMPILE_PREC);
  fprintf(stderr,"\tArguments have to be given in the order indicated above.\n");
  fprintf(stderr,"\tdinc:     increment in longitudinal and latitudinal direction\n");
  fprintf(stderr,"\t          (nlon and nlat will be reset if the FFT mode is chosen)\n");     
  fprintf(stderr,"\tsum_mode: %i      for sum over l and m\n",
	  SUM_OVER_L_AND_M);
  fprintf(stderr,"\t          %i      for sum over l and inverse FFT for m\n",
	  FFT_SUM);
  fprintf(stderr,"\tout_mode: %i      for ASCII output in one column\n",
	  ASCII_STDOUT);
  fprintf(stderr,"\t                 bottom to top slow, left to right fast index\n");
  fprintf(stderr,"\t                 can be read, eg., by piping to xyz2grd -ZBLa\n");
  fprintf(stderr,"\t          string for GMT grd output in file \"string.grd\"\n");
  fprintf(stderr,"\t          %i      for GMT style binary output to stdout\n",
	  BINARY_STDOUT);
  fprintf(stderr,"\t                 (scanline, ie. top to bottom, left to right)\n");
  fprintf(stderr,"\t                 can be read, eg., by piping to xyz2grd -ZTLf\n");
  fprintf(stderr,"\t          %i      expects AB for pol. and tor. coefficients\n",
	  TWO_GRDS);
  fprintf(stderr,"\t                 for GMT style binary output to %s (v_phi)and\n",FIRST_VEL_OUT);
  fprintf(stderr,"\t                 %s (v_theta)\n",SECOND_VEL_OUT);
  fprintf(stderr,"\t                 in this case, the format of the expansion file is:\n");
  fprintf(stderr,"\t                 l_max\n");
  fprintf(stderr,"\t                 AP_00 BP_00 AT_00 BT_00\n");
  fprintf(stderr,"\t                 AP_10 BP_10 AT_10 BT_10\n");
  fprintf(stderr,"\t                 AP_11 BP_11 AT_11 BT_11 ...)\n");
  fprintf(stderr,"\t          %i      output of GRADIENT of AB expanded field to %s and\n",
	  GRADIENT,FIRST_VEL_OUT);
  fprintf(stderr,"\t                 %s\n",SECOND_VEL_OUT);
  fprintf(stderr,"\t          This scheme implies that you can not create grid files with\n");
  fprintf(stderr,"\t          names of %i.grd, %i.grd, %i.grd, or %i.grd.\n",
	  ASCII_STDOUT,BINARY_STDOUT,TWO_GRDS,GRADIENT);
  fprintf(stderr,"\txmin...:  range for expansion in geographical coordinates\n");
  fprintf(stderr,"\tdinxy:    increment in y direction, if not set equal to increment in x\n");
  fprintf(stderr,"\n");
}

/* read in the spherical harmonic coefficients en block */

void read_ab(COMP_PRECISION *ap, COMP_PRECISION *bp, 
	     int lmax, char *mainname, FILE *in, 
	     int nr_of_coeff_sets,int lmsize,int verbose)
{
  int l,m,i;
  COMP_PRECISION *total_power,tmp_pow;
  total_power=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*nr_of_coeff_sets);
  for(i=0;i<nr_of_coeff_sets;i++)
    total_power[i]=0.0;
  
  for(l=0;l<=lmax;l++){
    for(m=0;m<=l;m++){
      for(i=0;i<nr_of_coeff_sets;i++){
	if((fscanf(IN_FILE,TWO_DATA_FSCAN_FORMAT,
		   (ap+i*lmsize+POSLM(l, m)),
		   (bp+i*lmsize+POSLM(l, m))))!=2){
	  fprintf(stderr,"%s: A B read error, set=%i l=%i m=%i\n",mainname,i,l,m);
	  fprintf(stderr,"%s: expecting %i sets of coefficients\n",mainname,nr_of_coeff_sets);
	  exit(-1);}
      }
    }
  }
  if(verbose){
    fprintf(stderr,"%s: read %i sets of coefficients with lmax=%i \n",
	   mainname,nr_of_coeff_sets,lmax);
    for(l=0;l<=lmax && l<5;l++){
      fprintf(stderr,"%s: ",mainname);
      for(i=0;i<nr_of_coeff_sets;i++){
	tmp_pow=degree_power((ap+i*lmsize),(bp+i*lmsize),l,TRUE);
	total_power[i] += tmp_pow;
	fprintf(stderr,"set:%3i l:%4i pow^2:%10g ",i+1,l,tmp_pow);
      }
      fprintf(stderr,"\n");
    }
    for(l=5;l<=lmax;l++)
      for(i=0;i<nr_of_coeff_sets;i++)
	total_power[i] += 
	  degree_power((ap+i*lmsize),(bp+i*lmsize),l,TRUE);
    for(i=0;i<nr_of_coeff_sets;i++)
      fprintf(stderr,"%s: set %i: total power^2 per lmax: %20g\n",
	      mainname,i+1,
	      total_power[i]/((lmax!=0)?
			      ((COMP_PRECISION)lmax):(1.0)));
  }
  free(total_power);
}
void check_out_mode(char *argument, int *out_mode,char *grdfilename)
{
  if(strcmp(argument,"0") == 0)
    *out_mode=0;
  else if(strcmp(argument,"2") == 0)
    *out_mode=2;
  else if(strcmp(argument,"3") == 0)
    *out_mode=3;
  else if(strcmp(argument,"4") == 0)
    *out_mode=4;
  else{
    sprintf(grdfilename,"%s.grd",argument);
    *out_mode=1;
  }
}
