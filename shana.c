/* part of the shansyn spherical harmonics package, 

   see COPYRIGHT for license 

   $Id: shana.c,v 1.30 2009/04/16 00:55:29 becker Exp becker $ 

*/
#include "shana.h"
int finite(double );		/* why doesn't this get included with math.h? */
#define myfinite(x) (finite((double)(x)))

void spline(COMP_PRECISION *,COMP_PRECISION *,
	    int ,COMP_PRECISION  ,COMP_PRECISION ,
	    COMP_PRECISION *);
/* 
   
   shana

   calculates a real spherical harmonic expansion using theoretical
   physics-type fully normalized spherical harmonics as specified in
   Dahlen and Tromp (Theoretical Global Seismology, 1998, Appendix B)

   uses 

   - binary or ASCII block gridded data on the surface of the Earth
   (can be in GMT/netcdf format) and solves by straightforward
   integration
   
   - irregularly distributed data and performs least squares fit

   - contour lines along which Delta functions are integrated, or

   - individual points whose Delta function contributions are added up

   type shana -h for information on the parameter settings


   TODO:

   clean up coding
   
   implement inverse FFT in longitude


   
   (c) Thorsten Becker 1999 - 2009 twb@usc.edu


*/
// check for floating exceptions in velocity routines
//
//#define CHECK_FLOATING_EXCEPTION

int main(int argc, char *argv[] )
{
  int lmax,lmc,lmax1,lmsize,nlon,nlat,i,j,opmode,nrp=0,intmode,verbose=BE_VERBOSE,flip=0,pixelreg,
    calculate_derivatives,vectors,nlontimesnlat,os1,os2;
  DATA_PRECISION *func=NULL,avg,min,max,*cloc,minphi,maxphi,
    mintheta,maxtheta,dtheta=0.,dphi=0.,*tmp;
  COMP_PRECISION damping,theta,phi;
  double tmpd;
  char *filename;
#ifdef USE_GMT4
  int ret_code;
  struct GRD_HEADER *header;
  GMT_LONG dummy[4]={0,0,0,0};
  header=(struct GRD_HEADER *)calloc(1,sizeof(struct GRD_HEADER));
  double wesn[6];
#else  /* new GMT API */
  void *API;                        /* The API control structure */
  struct GMT_GRID *G = NULL;        /* Structure to hold output grid */
  struct GMT_GRID_HEADER *header = NULL;
#endif

  filename=(char *)malloc(sizeof(char)*STRING_LENGTH*2);

#ifdef USE_GMT4  
  /*  */
  GMT_begin (argc, argv);
  GMT_grd_init (header, argc, argv, FALSE);
  GMT_program = argv[0];
#else
  API = GMT_Create_Session (argv[0], 2U, 0, NULL);
#endif
  
  if(argc > 1 && strcmp(argv[1],"-h") == 0 )argc=-9;

  calculate_derivatives=vectors=FALSE;
  lmax=DEF_LMAX;
  opmode=DEF_OPMODE;
  intmode=DEF_INT_MODE;
  damping = DEF_DAMPING;
  /* 
     argument assignment , check some of the flags
  */
  switch(argc){
  case 1:{
    break;
  }
  case 2:{
    sscanf(argv[1],"%i",&lmax);
    break;}
  case 3:{
    sscanf(argv[1],"%i",&lmax);
    check_opmode(argv[2],&opmode,filename,
		 &vectors,&calculate_derivatives);
    break;
  }
  case 4:{
    sscanf(argv[1],"%i",&lmax);
    check_opmode(argv[2],&opmode,filename,
		 &vectors,&calculate_derivatives);
    sscanf(argv[3],"%i",&intmode);
    break;
  }
  case 5:{
    sscanf(argv[1],"%i",&lmax);
    check_opmode(argv[2],&opmode,filename,
		 &vectors,&calculate_derivatives);
    sscanf(argv[3],"%i",&intmode);
    sscanf(argv[4],"%lf",&tmpd);damping = (COMP_PRECISION) tmpd;
    break;
  }
  default:{
    phelp(argv[0]);
    exit(-1);}}

  lmsize= (int)((((COMP_PRECISION)lmax)+1.0)*(((COMP_PRECISION)lmax)+2)/2.0);
  lmax1 = lmax+1;
  lmc = lmsize*2 - lmax1;
  if(lmax >= 0){
    if(verbose)fprintf(stderr,"%s: expansion to l_max=%i, DOF %i-1\n",argv[0],lmax,lmc);
  }else {
    fprintf(stderr,"%s: expansion to l_max=%i doesn't make sense\n",argv[0],lmax);
    phelp(argv[0]);exit(-1);
  }
  if(opmode < 0){
    flip=1;
    opmode*= -1;
    if(verbose)fprintf(stderr,"%s: latitudes are flipped\n",argv[0]);
  }
  if((lmax < 30) && (intmode == GAUSSIAN))
    fprintf(stderr,"%s: WARNING: Gaussian integration for low order expansion\n",argv[0]);
  

  
  /* 

     DATA INPUT 
     
  */
  // flags and default settings for input
  pixelreg=FALSE;
  switch(opmode){
    /* 

       read in binary block data of a global scalar function

    */
  case BLOCK:{
    if(verbose)fprintf(stderr,"%s: expecting binary block on offset (grid) coordinates\n",argv[0]);
    fread(&nlat,sizeof(int),1,stdin);
    fread(&nlon,sizeof(int),1,stdin);
    if(verbose)fprintf(stderr,"%s: read dimensions expecting nlon:%i times nlat:%i values\n",
		       argv[0],nlon,nlat);
    myrealloc_dp(&func,(nlontimesnlat=nlon*nlat));zero_dp(func,nlontimesnlat);
    fread(func,sizeof(DATA_PRECISION)*nlontimesnlat,1,stdin);
    minmax(&min,&max,&avg,func,nlontimesnlat);
    //
    // this are the expected input data coordinates
    // dx/2 <= x <= 360-dx/2 and -90+dy/2 <= y <= 90-dy/2
    //
    dphi=      TWOPI/((COMP_PRECISION)nlon);
    minphi=    0.5*dphi;
    maxphi=    minphi+(COMP_PRECISION)(nlon-1)*dphi;
    //
    dtheta=    PI/((COMP_PRECISION)nlat);
    mintheta=  PI - 0.5*dtheta;
    maxtheta=  mintheta-(COMP_PRECISION)(nlat-1)*dtheta;
    if(verbose)
      field_message(minphi,dphi,maxphi,mintheta,dtheta,maxtheta,
		    nlon,nlat,min,max,avg,argv[0]);
    break;
  }
    /* 

       lon lat value data

       CONTOUR:

       read in ASCII contours for integration of Delta functions along
       contours

       POINTS_FOR_AB:
       POINTS_FOR_A_MATRIX:
     
       read in ASCII lon lat coordinates for creation of an A matrix for
       least squares solution of fitting problem OR for spherical
       harmonic analysis of spotted delta function type points that have
       a scalar value

     
    */
  case POINTS_FOR_AB:
  case POINTS_FOR_A_MATRIX:
  case POINTS_FOR_A_MATRIX_ASCII:
  case POINTS_FOR_AVEC_MATRIX:
  case POINTS_FOR_AVEC_MATRIX_ASCII:
  case CONTOUR:{ 
    if(verbose){
      if(intmode == LEAST_SQUARES)
	fprintf(stderr,"%s: reading lon lat z for least squares fit of spherical harmonics\n",
		argv[0]);
      else{
	if(opmode == CONTOUR)
	  fprintf(stderr,"%s: reading contour lon lat z tripels along which to integrate\n",argv[0]);
	else if(opmode == POINTS_FOR_A_MATRIX)
	  fprintf(stderr,"%s: expecting lon lat z for A matrix (binary, single precision) output (z unused)\n",
		  argv[0]);
	else if(opmode == POINTS_FOR_AVEC_MATRIX)
	  fprintf(stderr,"%s: expecting lon lat z for A matrix vector field (binary, single precision) output (z unused)\n",
		  argv[0]);
	else if(opmode == POINTS_FOR_A_MATRIX_ASCII)
	  fprintf(stderr,"%s: expecting lon lat z for A matrix (ascii) output (z unused)\n",
		  argv[0]);
	else if(opmode == POINTS_FOR_AVEC_MATRIX_ASCII)
	  fprintf(stderr,"%s: expecting lon lat z for A matrix vector field (ascii) output (z unused)\n",
		  argv[0]);
	else if(opmode == POINTS_FOR_AB)
	  fprintf(stderr,"%s: expecting lon lat z for sum over Delta points\n",argv[0]);
      }
    }
    nrp=1;
    mymalloc_dp(&cloc,2*nrp);
    myrealloc_dp(&func,nrp);
    os1 = 0;
    while(fscanf(stdin,THREE_GMTDATA_FSCAN_FORMAT,
		 (cloc+os1),(cloc+os1+1),(func+(nrp-1))) == 3){
      /* reformat coordinates */
      *(cloc+os1)  =LONGITUDE2PHI( *(cloc+os1));
      *(cloc+os1+1)=LATITUDE2THETA(*(cloc+os1+1));
      nrp++;
      os1 += 2;
      myrealloc_dp(&cloc,2*nrp);myrealloc_dp(&func,nrp);
    }
    nrp--;
    if(verbose)fprintf(stderr,"%s: read %i points\n",argv[0],nrp);
    if(!nrp){
      fprintf(stderr,"%s: Exiting, no data.\n",argv[0]);exit(-1);}
    minmax(&min,&max,&avg,func,nrp);
    if((opmode == CONTOUR)&&(intmode != LEAST_SQUARES))	/* contour input */
      if(nrp!=1)
	for(i=0;i<nrp;i++)
	  if(dist_rad(cloc,i, (i<nrp-1)?(i+1):(i-1)) > 0.5*PI/(COMP_PRECISION)lmax){
	    fprintf(stderr,"%s: pts. %i and %i have low spacing of %g degrees\n",
		    argv[0],i,i+1,dist_rad(cloc,i, (i<nrp-1)?(i+1):(i-1))/PIOVERONEEIGHTY);
	    fprintf(stderr,"%s: maximum should be around %g for lmax=%i\n",
		    argv[0],90.0/(COMP_PRECISION)lmax,lmax);
	  }
    
    if(verbose)fprintf(stderr,"%s: min=%10g avg=%10g max=%10g\n",
		       argv[0],min,avg/((COMP_PRECISION)(nrp)),max);
    break;
  } 
    /* 

       read in ASCII block data of a global scalar function
  
    */
  case ASCII_BLOCK:{
    if(verbose)
      fprintf(stderr,"%s: reading ascii block on grid registration coordinates\n",argv[0]);
    nlon=1;

    myrealloc_dp(&func,nlon);
    while((fscanf(stdin,"%f ",(func+nlon-1)))==1){
      nlon++;
      myrealloc_dp(&func,nlon);
    }
    nlon--;
    nlat=(int)sqrt((COMP_PRECISION)nlon/2.);
    nlon=2*nlat;

    // swap ordering
    mymalloc_dp(&tmp,nlon*nlat);
    for(os1=i=0;i<nlat;i++,os1+=nlon)
      for(os2=j=0;j<nlon;j++,os2+=nlat)
	*(tmp+os2+i)= *(func+os1+j);
    j=nlon*nlat;
    for(i=0;i<j;i++)
      func[i]= tmp[i];
    free(tmp);

    // dx/2 <= x <= 360-dx/2 and -90+dy/2 <= y <= 90-dy/2
    if(verbose)
      fprintf(stderr,"%s: nlon=%i, nlat=%i\n",argv[0],nlon,nlat);
    minmax(&min,&max,&avg,func,nlat*nlon);
    dphi=  TWOPI/((COMP_PRECISION)nlon);
    minphi=0.5*dphi;
    maxphi=minphi+(COMP_PRECISION)(nlon-1)*dphi;
    dtheta=  PI/((COMP_PRECISION)nlat);
    mintheta=PI - 0.5*dtheta;
    maxtheta=mintheta - (COMP_PRECISION)(nlat-1)*dtheta;
    if(verbose)
      field_message(minphi,dphi,maxphi,mintheta,dtheta,maxtheta,nlon,nlat,min,max,avg,argv[0]);
    break;
  }
    /* 
       
       read in ASCII block data of a global scalar function with header
       
    */
  case ASCII_BLOCK_HEADER:{
    if(verbose)
      fprintf(stderr,"%s: reading ascii block with header\n",argv[0]);
    
    if(fscanf(stdin,"%i %i %f %f %f %f\n",&nlon,&nlat,&minphi,&maxphi,&mintheta,&maxtheta)!=6){
      fprintf(stderr,"%s: need header line with nlon nlat minphi maxphi mintheta maxtheta\n\n",argv[0]);
      exit(-1);
    }
    minphi=  LONGITUDE2PHI(minphi);
    maxphi=  LONGITUDE2PHI(maxphi);
    mintheta=LATITUDE2THETA(mintheta);
    maxtheta=LATITUDE2THETA(maxtheta);
    dphi=  (maxphi-minphi)    /((COMP_PRECISION)(nlon-1));
    dtheta=fabs((maxtheta-mintheta)/((COMP_PRECISION)(nlat-1)));
    myrealloc_dp(&func,(nlontimesnlat=nlon*nlat));
    /* read in data */
    for(i=0;i<nlat;i++)
      for(os1=j=0;j<nlon;j++,os1+=nlat)
	if(fscanf(stdin,"%f ",(func+os1+i))!=1){
	  fprintf(stderr,"%s: read error at  %i %i\n",argv[0],i,j);
	  exit(-1);
	}
    if(flip)/* flip the latitudes */
      gmt2myconvention_rotate4(func,nlon,nlat,1);
    minmax(&min,&max,&avg,func,nlontimesnlat);
    if(verbose)
      field_message(minphi,dphi,maxphi,mintheta,dtheta,maxtheta,nlon,nlat,min,max,avg,argv[0]);
    break;
  }
    /* 
     
       read in GMT grd file of scalar data 

    */
  case GMT_BLOCK:{
#ifdef USE_GMT4
    if((ret_code = GMT_read_grd_info (filename,header))< 0){
      /* return code has changed for GMT4, i presume? */
      fprintf(stderr,"%s: error opening %s for grdinfo\n",argv[0],filename);
      exit(-1);
    }else{
      pixelreg=(header->node_offset?TRUE:FALSE);
      minphi=LONGITUDE2PHI(header->x_min+(pixelreg?header->x_inc/2.0:0.0));	
      maxphi=LONGITUDE2PHI(header->x_max-(pixelreg?header->x_inc/2.0:0.0));
      mintheta=LATITUDE2THETA(header->y_min+(pixelreg?header->y_inc/2.0:0.0));	
      maxtheta=LATITUDE2THETA(header->y_max-(pixelreg?header->y_inc/2.0:0.0));
      dphi=  DEG2RAD( header->x_inc);
      dtheta=DEG2RAD( header->y_inc);
      // for grid line registered nx_i = (x_i^max-x_i^min)/dx_i + 1
      // for pixel reg            nx_i = (x_i^max-x_i^min)/dx_i
      nlon=header->nx;
      nlat=header->ny;
    }
#else
    /* read header info */
    if((G = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE,
			   GMT_IS_SURFACE, GMT_CONTAINER_ONLY, NULL, filename, NULL))==NULL)
      return (-1);
    header = G->header;

    pixelreg=(header->registration == GMT_GRID_PIXEL_REG)?TRUE:FALSE;
    minphi=LONGITUDE2PHI(header->wesn[XLO]+(pixelreg?header->inc[GMT_X]/2.0:0.0));	
    maxphi=LONGITUDE2PHI(header->wesn[XHI]-(pixelreg?header->inc[GMT_X]/2.0:0.0));
    mintheta=LATITUDE2THETA(header->wesn[YLO]+(pixelreg?header->inc[GMT_Y]/2.0:0.0));	
    maxtheta=LATITUDE2THETA(header->wesn[YHI]-(pixelreg?header->inc[GMT_Y]/2.0:0.0));
    dphi=  DEG2RAD( header->inc[GMT_X]);
    dtheta=DEG2RAD( header->inc[GMT_Y]);
    // for grid line registered nx_i = (x_i^max-x_i^min)/dx_i + 1
    // for pixel reg            nx_i = (x_i^max-x_i^min)/dx_i
    nlon=header->n_columns;
    nlat=header->n_rows;
    
#endif
    if(!vectors){ 
      /* 
	 read in scalar data field from GMT grd file 
      */
      fprintf(stderr,"%s: reading from grd-file %s\n",argv[0],filename);
      myrealloc_dp(&func,nlon*nlat);
#ifdef USE_GMT4
      GMT_read_grd(filename,header, func, 0.0, 0.0, 0.0, 0.0, dummy, 0);
      gmt2myconvention_rotate4(func,nlon,nlat,1.0);
#else
      /* read data */
      if((G = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE,
			     GMT_IS_SURFACE, GMT_DATA_ONLY, NULL, filename, G))==NULL)
	return (-1);

      gmt2myconvention_rotate(func,nlon,nlat,1,G);
#endif
      minmax(&min,&max,&avg,func,(nlontimesnlat=nlat*nlon));
      if(verbose)
	field_message(minphi,dphi,maxphi,mintheta,dtheta,maxtheta,nlon,nlat,min,max,avg,argv[0]);
      if((minphi == maxphi) || (mintheta == maxtheta)){
	fprintf(stderr,"%s: use arrays with non-overlapping borders: %g <= x <= %g, %g <= y <= %g\n",
		argv[0],PHI2LONGITUDE(minphi),PHI2LONGITUDE(maxphi),
		THETA2LATITUDE(mintheta),THETA2LATITUDE(maxtheta));
	exit(-1);
      }
    }else{ 
      /* 
	 read in vectors from two grd files 
	 use dimensions specs from first one and check if second has identical
	 bounds
      */
      if(verbose)
	fprintf(stderr,"%s: reading phi component from grd-file %s\n",argv[0],filename);
      myrealloc_dp(&func,(nlontimesnlat=nlon*nlat)*2);
      /* read from grid */
#ifdef USE_GMT4
      GMT_read_grd (filename,header, func, 0.0, 0.0, 0.0, 0.0, dummy, 0);
      gmt2myconvention_rotate4(func,nlon,nlat,1.0);
#else  /* read datA */
      if((G = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE,
			     GMT_IS_SURFACE, GMT_DATA_ONLY, NULL, filename, G))==NULL)
	return (-1);
      gmt2myconvention_rotate(func,nlon,nlat,1.0,G);
#endif
      minmax(&min,&max,&avg,func,nlontimesnlat);
      if(verbose)
	field_message(minphi,dphi,maxphi,mintheta,dtheta,maxtheta,nlon,nlat,min,max,avg,argv[0]);
      /* second velocity file */
      if(verbose)
	fprintf(stderr,"%s:     and theta component from grd-file %s\n",argv[0],(filename+STRING_LENGTH));

#ifdef USE_GMT4
      if(GMT_read_grd_info ((filename+STRING_LENGTH),header)== -1){
	fprintf(stderr,"%s: error opening %s for grdinfo\n",argv[0],filename);
	exit(-1);
      }
      if(SIGNIFICANTLY_DIFFERENT(header->x_min, PHI2LONGITUDE(minphi)) ||
	 SIGNIFICANTLY_DIFFERENT(header->y_min, THETA2LATITUDE(mintheta)) ||
	 SIGNIFICANTLY_DIFFERENT(header->x_max, PHI2LONGITUDE(maxphi)) ||
	 SIGNIFICANTLY_DIFFERENT(header->y_max, THETA2LATITUDE(maxtheta)) ||
	 SIGNIFICANTLY_DIFFERENT(header->nx, nlon) || 
	 SIGNIFICANTLY_DIFFERENT(header->ny, nlat)){
	fprintf(stderr,"%s: %s and %s have to have identical size\n",
		argv[0],(filename),(filename+STRING_LENGTH));
	fprintf(stderr,"%s: delta xmin: %g\n",
		argv[0],header->x_min-PHI2LONGITUDE(minphi));
	fprintf(stderr,"%s: delta ymin: %g\n",
		argv[0],header->y_min-THETA2LATITUDE(mintheta));
	fprintf(stderr,"%s: delta xmax: %g\n",
		argv[0],header->x_max-PHI2LONGITUDE(maxphi));
	fprintf(stderr,"%s: delta ymax: %g\n",
		argv[0],header->y_max-THETA2LATITUDE(maxtheta));
	fprintf(stderr,"%s: delta nx:   %i\n",
		argv[0],header->nx -nlon );
	fprintf(stderr,"%s: delta ny:   %i\n",
		argv[0],header->ny -nlat );
	exit(-1);
      }
      GMT_read_grd ((filename+STRING_LENGTH),header,
		    (func+nlontimesnlat), 0.0, 0.0, 0.0, 0.0, dummy, 0);
      gmt2myconvention_rotate4((func+nlontimesnlat),nlon,nlat,1.0);
      
#else
      if((G = GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_CONTAINER_ONLY, NULL, (filename + STRING_LENGTH), NULL))==NULL)
       return (-1);
      header = G->header;
      if(SIGNIFICANTLY_DIFFERENT(header->wesn[XLO], PHI2LONGITUDE(minphi)) ||
	 SIGNIFICANTLY_DIFFERENT(header->wesn[YLO], THETA2LATITUDE(mintheta)) ||
	 SIGNIFICANTLY_DIFFERENT(header->wesn[XHI], PHI2LONGITUDE(maxphi)) ||
	 SIGNIFICANTLY_DIFFERENT(header->wesn[YHI], THETA2LATITUDE(maxtheta)) ||
	 SIGNIFICANTLY_DIFFERENT(header->n_columns, nlon) || 
	 SIGNIFICANTLY_DIFFERENT(header->n_rows, nlat)){
	fprintf(stderr,"%s: %s and %s have to have identical size\n",
		argv[0],(filename),(filename+STRING_LENGTH));
	fprintf(stderr,"%s: delta xmin: %g\n",
		argv[0],header->wesn[XLO]-PHI2LONGITUDE(minphi));
	fprintf(stderr,"%s: delta ymin: %g\n",
		argv[0],header->wesn[YLO]-THETA2LATITUDE(mintheta));
	fprintf(stderr,"%s: delta xmax: %g\n",
		argv[0],header->wesn[XHI]-PHI2LONGITUDE(maxphi));
	fprintf(stderr,"%s: delta ymax: %g\n",
		argv[0],header->wesn[YHI]-THETA2LATITUDE(maxtheta));
	fprintf(stderr,"%s: delta nx:   %i\n",
		argv[0],header->n_columns -nlon );
	fprintf(stderr,"%s: delta ny:   %i\n",
		argv[0],header->n_rows -nlat );
	exit(-1);
      }
      if((G = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE,
			     GMT_IS_SURFACE, GMT_DATA_ONLY, NULL, (filename+STRING_LENGTH), G))==NULL)
	return (-1);
      gmt2myconvention_rotate((func+nlontimesnlat),nlon,nlat,1.0,G);
#endif
      minmax(&min,&max,&avg,(func+nlontimesnlat),nlontimesnlat);
      if(verbose)
	field_message(minphi,dphi,maxphi,mintheta,
		      dtheta,maxtheta,nlon,nlat,min,max,avg,argv[0]);
    } /* end vectors branch */
    if(intmode == LEAST_SQUARES){
      /* 
	 need to initialize cloc array, my convention
      */
      mymalloc_dp(&cloc,nlontimesnlat*2);
      for(j=0,phi=minphi;j < nlon;j++,phi += dphi){ /*  */
	for(i=0,theta=mintheta;i<nlat;i++,theta -= dtheta){ /* theta is the fast index */
	  *(cloc+(j*nlat+i)*2)  =phi;
	  *(cloc+(j*nlat+i)*2+1)=theta;
	}
      }
      nrp = nlontimesnlat;
    }
    break;
  }
  default:{
    fprintf(stderr,"%s: opmode %i is not supported\n",argv[0],opmode);exit(-1);
  }}
  if((opmode!=CONTOUR) && (opmode!=POINTS_FOR_A_MATRIX) && (opmode!=POINTS_FOR_AVEC_MATRIX) 
     && (opmode!=POINTS_FOR_A_MATRIX_ASCII)&& (opmode!=POINTS_FOR_AVEC_MATRIX_ASCII)	\
     && (opmode!=POINTS_FOR_AB) 
     && (nlat < lmax+1)){
    fprintf(stderr,"%s: nlat (%i) is lower than l_{max}+1 (%i).\n",argv[0],nlat,lmax+1);
    fprintf(stderr,"%s: you might want to use lmax=%i instead\n",argv[0],nlat-1);
    if(intmode == GAUSSIAN){
      fprintf(stderr,"%s: for Gaussian integration mode, we need at least lmax+1 points\n",
	      argv[0]);
      exit(-1);
    }
  }
  if(verbose)fprintf(stderr,"%s: input done, computing coefficients...\n",argv[0]);
  calc_coeff(&func,lmax,nlon,nlat,argv[0],opmode,cloc,nrp,intmode,verbose,
	     minphi,maxphi,dphi,mintheta,maxtheta,dtheta,
	     vectors, calculate_derivatives,pixelreg,damping);
  if(verbose)
    fprintf(stderr,"%s: Done.\n",argv[0]);

#ifdef USE_GMT4
  // this would be nice, but it saometimes gave a double free error, need to check
  //GMT_end (argc, argv);
#else
  GMT_Destroy_Session (API);
#endif

	
  return 0;
}

/*************************************************************/
/* 

   routine to do the SPHERICAL HARMONIC EXPANSION

*/
/*************************************************************/
void calc_coeff(DATA_PRECISION **func,int lmax,
		int nlon,int nlat, char *program,
		int opmode,DATA_PRECISION *cloc, int nrp,
		int intmode,int verbose,DATA_PRECISION minphi, 
		DATA_PRECISION maxphi,
		DATA_PRECISION dphi,DATA_PRECISION mintheta,
		DATA_PRECISION maxtheta, 
		DATA_PRECISION dtheta, int vectors, 
		int calculate_derivatives, 
		int pixelreg, COMP_PRECISION damping)
{
  COMP_PRECISION *y=NULL,*x=NULL,
    *p=NULL,theta,phi,f[3][2],*sinfac=NULL,*dptheta=NULL,tmp,tmp2,
    *cosarr=NULL,*sinarr=NULL,atmp[2],btmp[2],ctmp[2],dist,
    *oneoverl=NULL,sinphi,cosphi,
    *tmp_func,*absc=NULL,sintheta,costheta,dummy,*sndder,phirange,
    thetarange,norm,chi2,sum,dnorm;
  static COMP_PRECISION normfac,sqrt2normfac;
  int i=0,j,k,l,m,n,lmsize,lmax1,memneeded,lookup_p=1,coeff,nlontimesnlat,klim,
    lmc,nr_gauss_pts,lmax_max,os1,os2,ic1,ic2,os3,nsize[3],nrp2,
    jc1,jc2,phisumlim,thetasumlim,ysumn,warned,damping_mode;
  float *aflt = NULL,*sol=NULL, *sigma = NULL,*aflt_copy=NULL;
  DATA_PRECISION *tmp_dfunc,*tmp_cloc;
  /* 
     y is cos(theta) 
     x is phi
     cosarr and sinarr hold sin(m phi) and cos(m phi)
     p is the legendre function arrar
     DPTHETA and DP the derivatives of the legendre functions (factors)
     
  */
  y=x=sinarr=cosarr=NULL;
#ifdef REDUCE_LMAX
  /*
    this is some sort of check if we have enough data points
    for trapeziodal integration: degree l would normally need
    2l data points to be exact;
  */
  lmax_max=(int)(nlat-1)/2;
#else
  lmax_max=lmax;
#endif
  // some size and length variables
  lmsize= (int)((((COMP_PRECISION)lmax)+1.0)*(((COMP_PRECISION)lmax)+2)/2.0);
  lmax1=lmax+1;
  // number of Gauss integration points, should be lmax+1
  nr_gauss_pts=lmax1;
  nlontimesnlat=nlon*nlat;

  if(intmode == LEAST_SQUARES){
    /* 

       
       LEAST SQUARES ESTIMATE


    */
    /* switch damping mode */
    if(fabs(damping) > 0)
      if(damping < 0){
	damping_mode = 2;		/* gradient */
	damping = -damping;
      }
      else
	damping_mode = 1;		/* norm */
    else
      damping_mode = 0;
    if(verbose){
      fprintf(stderr,"%s: least squares inversion, ",program);
      if(damping_mode)
	fprintf(stderr,"%s damping with %g\n",(damping_mode == 1)?("norm"):("roughness"),damping);
      else
	fprintf(stderr,"no damping\n");
    }
    /* 

       start inversion part 

    */
    switch(opmode){
    case POINTS_FOR_A_MATRIX:    
    case POINTS_FOR_A_MATRIX_ASCII:
    case POINTS_FOR_AVEC_MATRIX:    
    case POINTS_FOR_AVEC_MATRIX_ASCII:
    case POINTS_FOR_AB:
    case CONTOUR:
    case GMT_BLOCK:
      /* 
	 weights for every data point
      */
      nrp2 = vectors ? nrp * 2 : nrp;
      myrealloc_dp(&sigma,nrp2);
      /* assign sigma = 1 / weight */
      if(opmode != GMT_BLOCK){
	for(i=0;i < nrp2;i++)
	  sigma[i] = 1.0;
      }else{			/* GMT block */
	/* 
	   block data, weigh by latitude 
	*/
	k = n = 0;
	tmp_dfunc = tmp_cloc = NULL;
	for(i=0;i < nrp;i++){
	  sum = sin(cloc[i*2+1]); /* sin(theta) */
	  if(fabs(sum) > EPS_DATA_PREC){
	    /* use this point */
	    sigma[n] = 1/sum;
	    /* location */
	    myrealloc_dp(&tmp_cloc,2*(k+1));
	    tmp_cloc[k*2] = cloc[i*2];
	    tmp_cloc[k*2+1]=cloc[i*2+1];
	    /* first function value */
	    myrealloc_dp(&tmp_dfunc,(n+1));
	    tmp_dfunc[n] = *(*func+i);
	    n++;
	    if(vectors){	/*  */
	      myrealloc_dp(&tmp_dfunc,(n+1));
	      tmp_dfunc[n] = *(*func+nrp+i);
	      sigma[n] = 1/sum;
	      n++;
	    }
	    k++;
	  }
	}
	if(n != nrp2){fprintf(stderr,"logic error 001, %i %i\n",n,nrp2);exit(-1);}

	fprintf(stderr,"%s: used %i out of %i total points (%i total values), rejected %i polar latitude points in GMT block\n",
		program,k,nrp,n,nrp-k);
	nrp = k;		/* shrink */
	nrp2 = vectors ? nrp * 2 : nrp;
	myrealloc_dp(&sigma,nrp2);
	myrealloc_dp(&cloc,nrp*2);
	for(i=0;i < nrp;i++){	/* assign */
	  cloc[i*2]   = tmp_cloc[i*2];
	  cloc[i*2+1] = tmp_cloc[i*2+1];
	}
	if(vectors){
	  if(n != nrp2){
	    fprintf(stderr,"logic error 12\n");exit(-1);}
	}
	myrealloc_dp(func,nrp2);
	for(i=0;i < nrp2;i++)	/* assign */
	  *(*func+i) = tmp_dfunc[i];
	free(tmp_dfunc);free(tmp_cloc);
      }
      /*

	 get the A matrix, fortran style, and use sigma

      */
      make_a(lmax,lmsize,nrp,&aflt,cloc,verbose,
	     program,damping_mode,nsize,damping,TRUE,TRUE,sigma,vectors);
      /* 
	 save, for variance reduction (choose speed over mem usage) 
      */
      myrealloc_dp(&aflt_copy,nsize[2]);
      memcpy(aflt_copy,aflt,(size_t)(nsize[2])*sizeof(float));
      /* 
	 assign data to solution vector which will get overwritten 
      */
      myrealloc_dp(&sol,nsize[0]);
      for(i=0;i < nrp2;i++)	/* function values */
	sol[i] = *(*func+i)/sigma[i];
      for(;i < nsize[0];i++)	/* for damping */
	sol[i] = 0.0;
      /* solve and overwrite function values with solution
	 coefficients and a matrix with something else */
#ifdef USE_LAPACK
      solver_ab_lls(aflt,nsize[0],nsize[1],sol,program);
#else
      fprintf(stderr,"%s: recompile with LAPACK support, need solver_ab_lls for least squares \n",program);
      exit(-1);
#endif
      /* don't need that output */
      free(aflt);
      /* 
	 assess solution 
      */
      /* solution norm */
      norm = chi2 = 0.0;
      for(i=0;i < nsize[1];i++)
	norm += sol[i]*sol[i];
      norm = sqrt(norm);
      /* 
	 misfit 
      */
      chi2 = dnorm = 0.0;
      klim = vectors ? 2 : 1;
      n = 0;
      for(i=0;i < nrp;i++){
	for(k = 0; k < klim;k++){
	  sum = 0.0;
	  for(j=0;j < nsize[1];j++) /* aflt is fortran style and
				       weighted, need to undo */
	    sum += aflt_copy[j*nsize[0]+n] * sol[j];
	  
	  sum *= sigma[n];	/* undo 1/sigma scaling for A(i,j) */
	  /* misfit/sigma */
	  sum = (sum - (*(*func+k*nrp+i)))/sigma[n];
	  sum = sum * sum;
	  chi2 += sum;
	  dnorm += (*(*func+k*nrp+i)) * (*(*func+k*nrp+i));
	  n++;
	}
      }
      dnorm = sqrt(dnorm);
      fprintf(stderr,"%s: solution DOF: %i solution norm: %.5e VR(1-chi/|x|): %.1f%%\n",
	      program,nsize[1],norm,(1-sqrt(chi2)/dnorm)*100);
      free(aflt_copy);free(sigma);		/* no use anymore  */

      /* 
	 output of solution 
      */
      fprintf(stdout,"%i\n",lmax);
      if(vectors){
	fprintf(stdout,COEFF_DATA_FORMAT,0.0,0.0);fprintf(stdout,"\t");fprintf(stdout,COEFF_DATA_FORMAT,0.0,0.0); /* no degree zero term */
	fprintf(stdout,"\n");
	
	for(j=0,l=1;l <= lmax;l++){
	  for(m=0;m <= l;m++){
	    if(m == 0){
	      fprintf(stdout,COEFF_DATA_FORMAT,sol[j],0.0); /* A */
	      j++;
	      fprintf(stdout,"\t");	/* C */
	      fprintf(stdout,COEFF_DATA_FORMAT,sol[j],0.0);
	      j++;
	    }else{
	      fprintf(stdout,COEFF_DATA_FORMAT,sol[j],sol[j+1]); /* AB */
	      j += 2;
	      fprintf(stdout,"\t");
	      fprintf(stdout,COEFF_DATA_FORMAT,sol[j],sol[j+1]); /* CD */
	      j += 2;
	    }
	    fprintf(stdout,"\n");
	  }	    
	}
      }else{			/* scalar */
	for(j=l=0;l <= lmax;l++){
	  for(m=0;m <= l;m++){
	    /* scalar */
	    if(m == 0){
	      fprintf(stdout,COEFF_DATA_FORMAT,sol[j],0.0);
	      j++;
	    }else{
	      fprintf(stdout,COEFF_DATA_FORMAT,sol[j],sol[j+1]);
	      j += 2;
	    }
	    fprintf(stdout,"\n");
	  }
	}
      }
      free(sol);
      break;
    default:
      fprintf(stderr,"opmode %i not implemented for least squares yet\n",opmode);
      exit(-1);
    }
    /* 
       END LEAST SQUARES
       
     */
  }else{
    /* 
       
       ALL OTHER STRAIGHTFORWARD INTEGRATION MODES

    */
    //
    // decide on integration mode depending on which opmode is chosen
    // and determine summing bounds and alike
    //
    switch(opmode){
    case GMT_BLOCK:
    case ASCII_BLOCK_HEADER:
    case ASCII_BLOCK:
    case BLOCK:{
      /* use block data */ 
      phirange=fabs(maxphi-minphi);
      thetarange=fabs(maxtheta-mintheta);
      // check the boundaries in longitude to see if full geographical coverage
      ic1=(fabs(phirange-TWOPI)<1e-6);
      ic2=(fabs(phirange-TWOPI+dphi)<1e-6);
      if(!ic1 && !ic2){
	fprintf(stderr,"%s: expecting full longitudinal coverage, ics: %i %i, range: %g dx: %g\n",
		program,ic1,ic2,PHI2LONGITUDE(phirange),PHI2LONGITUDE(dphi));
	exit(-1);
      }else{
	if(ic1){
	  /* 
	     full 360 range with redundant repetition 
	     of phi=0 value. limit summation to 0...360-dx
	  */
	  phisumlim=nlon-1;
	  // check if values at Greenwich are same
	  for(os1=(nlon-1)*nlat,i=warned=0;i<nlat;i++)
	    if((*func)[i] != ((*func)[os1+i])){
	      if(!warned){
		fprintf(stderr,
			"%s: at least one mismatch at the Greenwich, averaging\n",
			program);
		warned=1;
	      }
	      // average in this case
	      (*func)[i]=((*func)[i] + ((*func)[os1+i]))/2.0;
	    }
	}else{ // Greenwich only in field once, use all points in x direction
	  phisumlim=nlon;
	}
      }
      // now check latitudes
      jc1=(fabs(thetarange-PI)<1e-2);
      jc2=(fabs(thetarange-PI+dtheta)<1e-2);
      if(!jc1 && !jc2){
	fprintf(stderr,"%s: expecting full latitudinal coverage, range: %g (%g), %g(%g) %g(%g)\n",
		program,RAD2DEG(thetarange),RAD2DEG(PI-dtheta),mintheta,THETA2LATITUDE(mintheta),
		maxtheta,THETA2LATITUDE(maxtheta));
	exit(-1);
      }else{
	thetasumlim=nlat;
	if(jc1){
	  /* 
	     full 180 range with both North and South pole values 
	     use only half values at poles if we use sum rule, for Gauss, it's OK
	  */
	  if(intmode == TRAPEZOIDAL)
	    for(os1=nlat-1,j=i=0;i<nlon;i++,j+=nlat){
	      ((*func)[j])     /= 2.0;
	      ((*func)[j+os1]) /= 2.0;
	    }
	  // we should count both polar entries only as a half box, thus 
	  ysumn=nlat-1;
	  // NOTE: this is different from the integration loop end which is 
	  // thetasumlim and always equal to nlat
	}else{ 
	  // check if whole lat range covered, ie. if offset like -90+dy/2 and 90-dy/2
	  if((fabs(mintheta-PI+dtheta/2.0)>1e-6)||
	     (fabs(maxtheta-dtheta/2.0)>1e-6)){
	    fprintf(stderr,"%s: error: use latitude range from %g to %g for %g spacing\n",
		    program,THETA2LATITUDE(PI-dtheta/2.0),
		    THETA2LATITUDE(dtheta/2.0),RAD2DEG(dtheta));
	    exit(-1);
	  }
	  // all boxes count in this case to divide PI
	  ysumn=nlat;
	}
      }
      /* 
	 x coordinates (0 to 2 pi) if pixel registration, shift from original values 
	 to mid of cell
      */    
      myrealloc_cp(&x,nlon);zero_cp(x,nlon);
      for(j=0,phi=minphi; j < nlon;j++,phi += dphi)
	x[j]= phi;
      /* 
	 y coordinates, form 0 (north) to pi (south), y is cos(theta)
      */
      myrealloc_cp(&y,nlat);zero_cp(y,nlat);
      // shouldn't use dtheta since this is positive by definition
      tmp=(COMP_PRECISION)(maxtheta-mintheta)/(COMP_PRECISION)(nlat-1);
      // y coords from now on are in cos(theta)
      for(i=0,theta=mintheta; i<nlat;i++,theta += tmp){
	y[i] = cos(theta);
	/* fprintf(stderr,"i: %i nlat: %i y[i]: %12g\n",i,nlat,y[i]); */
      }
#ifdef DEBUG
      fprintf(stderr,"%s: integration coord range: x: %g to %g y: %g(%g) to %g(%g)\n",
	      program,
	      PHI2LONGITUDE(x[0]), PHI2LONGITUDE(x[phisumlim-1]),
	      THETA2LATITUDE(acos(y[0])),acos(y[0]),
	      THETA2LATITUDE(acos(y[thetasumlim-1])),
	      acos(y[thetasumlim-1]));
#endif
      if((intmode == TRAPEZOIDAL) || (intmode == GAUSSIAN)){
	/* 
	   cos and sin factors for R l m  and S l m, 
	   cos(m phi) and sin(m phi)  
	   would not need that for FFT method
	*/
	myrealloc_cp(&cosarr,nlon*lmax1);zero_cp(cosarr,nlon*lmax1);
	myrealloc_cp(&sinarr,nlon*lmax1);zero_cp(sinarr,nlon*lmax1);

	// m=0 terms
	for(j=0; j < nlon ;j++){
	  cosarr[j]=1.0; sinarr[j]=0.0;
	}
	// 1 <= m <= lmax
	for(os1=nlon,m=1;m<=lmax;m++,os1+=nlon){
	  for(j=0; j < nlon ;j++){
	    tmp = ((COMP_PRECISION)m)*x[j];
	    cosarr[os1+j]=cos(tmp);
	    sinarr[os1+j]=sin(tmp); 
	  }
	}
      } 
      /* 
	 Legendre functions, assign P_lm(theta) to array P 
	 using y=cos(theta)[j] 
      */
      if(verbose)
	fprintf(stderr,
		"%s: constructing Assoc. legendre functions...\n",program);
      myrealloc_cp(&p,lmsize*nlat);zero_cp(p,lmsize*nlat);
      if(verbose)
	fprintf(stderr,"%s: using %6g MBs for P array\n",
		program,
		lmsize*sizeof(COMP_PRECISION)*nlat/ONE_MEGABYTE);
      if(intmode != GAUSSIAN){
	//
	// obtain Legendre functions at the latitudinal points of the input data
	//
	plgndr(y,nlat,p,lmax);
	// set these such that derivatives and sinfac etc. 
	// work for Gauss points and original input which is given on nlat points
	//
	nr_gauss_pts=nlat;
	absc=y;// absc values are the y coordinates
      }else{
	/* 
	   obtain Legendre functions at Gauss integration points
	   (should be lmax+1), the P_lm will then also be already  
	   multiplied by the weight factor for quadrature

	   we furthermore obtain the abscissae values from this 
	   routine, to replace the y[j] later (still need y for interpolation later)
	*/
	myrealloc_cp(&absc,nr_gauss_pts);
	plgndr_g(p,lmax,nr_gauss_pts,absc);
      }
      /* 
	 integretation area sin (theta) d theta d phi, sine factor, 
	 sin(arccos(y))=sqrt(1-y^2) 
	 also used for derivatives, therefore calculate here
	 for real abscissae values, not necessary y
      */
      myrealloc_cp(&sinfac,nr_gauss_pts);zero_cp(sinfac,nr_gauss_pts);
      for(i=0;i<nr_gauss_pts;i++)
	sinfac[i] = (((tmp=(1.0-absc[i]*absc[i]))>=0.0)?(sqrt(tmp)):(0.0));
      if(calculate_derivatives){
	/* 
	   to expand vector fields calculate derivative of P with respect to theta 
	   --> DPTHETA, 
	   and overwrite P with m/sin(theta) P which is 
	   d_phi P, sort of 
	*/
	if(verbose)
	  fprintf(stderr,"%s: constructing derivatives....\n",program);
	myrealloc_cp(&dptheta,lmsize*nr_gauss_pts);zero_cp(dptheta,lmsize*nr_gauss_pts);
	if(verbose)
	  fprintf(stderr,"%s: using %6g MBs for DPTHETA array\n",
		  program,
		  lmsize*sizeof(COMP_PRECISION)*nr_gauss_pts/ONE_MEGABYTE);
	/* 
	   create the dptheta array which will hold d_theta X_lm 
	   based on the X_lm which we calculated before (weights might
	   have been included for Gauss integration, that should be 
	   fine)
	*/
	pdtheta_lgndr(absc,nr_gauss_pts,p,dptheta,lmax,&dummy,FALSE);
	/* 
	   overwrite P with d_phi X_lm which is 
	   m/sin(theta) X_lm WITH A AND B FLIPPED FOR SYNTHESIS 
	*/
	if(verbose){
	  fprintf(stderr,"%s: overwriting P array with d_phi P\n",
		  program);
	  fprintf(stderr,"%s: sin(theta) factor in Legendre and DPTHETA factors\n",
		  program);
	}
	myrealloc_cp(&oneoverl,lmax1);
	// 1/sqrt(l(l+1)) factors, polar and toroidal harmonics are = 0 for l=0
	oneoverl[0]=1.0;
	for(l=1;l <= lmax;l++){
	  oneoverl[l]=1.0/sqrt((COMP_PRECISION)(l*(l+1)));
	}
	if(intmode != GAUSSIAN){
	  /*
	    we will have to multiply everything by sin(theta) for the block style 
	    integration with even spacing along y. nr_gauss_pts is simply nlat
	  */
	  for(j=0;j<nr_gauss_pts;j++){
	    for(l=0;l<=lmax;l++){
	      for(m=0;m<=l;m++){
		/* do not devide P by sin(theta) since we 
		   would have multiplied by sin(theta) 
		   for the integration */
		P(l,m,j)       *= ((COMP_PRECISION)m)*oneoverl[l];
		DPTHETA(l,m,j) *= oneoverl[l] * sinfac[j];
	      }
	    }
	  }
	}else{
	  /*
	    for Gauss integration, do not multiply with 
	    sin(theta) 
	  */
	  for(j=0;j < nr_gauss_pts;j++){
	    for(l=0;l<=lmax;l++){
	      for(m=0;m<=l;m++){
		P(l,m,j)       *= ((COMP_PRECISION)m)*oneoverl[l];
		P(l,m,j)       /= sinfac[j];
		DPTHETA(l,m,j) *= oneoverl[l];
	      }
	    }
	  }
	  free(absc);
	}
      }else{// no vectors, normal integration
	if(intmode != GAUSSIAN){
	  /* 
	     no derivatives involved, multiply Legendre 
	     function with the appropriate sin(theta) 
	     factor for summing up integration in y direction 
	     when y are equidistantly spaced. if we use Gauss integration,
	     we do not need this factor since the w_i are spaced
	     appropriately
	  */
	  if(verbose)
	    fprintf(stderr,"%s: including sin(theta) factor in Legendre factor\n",
		    program);
	  for(j=0;j<nlat;j++)
	    for(l=0;l<=lmax;l++)
	      for(m=0;m<=l;m++)
		P(l,m,j) *= sinfac[j];
	}
      }
      free(sinfac);
      /*

	BEGIN INTEGRATION HERE

      */
      switch(intmode){
	/*
	  select which integration mode should be used
	*/
      case TRAPEZOIDAL:{ 
	/* 
	   trapezoidal rule, we assume data is in the middle of cells and we
	   are interested in the integral on the whole surface. we can therefore
	   simply sum every contribution in longitude (last cell was eliminated above
	   if redundant second Greenwich value given). 
	   if range is -90 to 90, we have already divided the pole values by 2
	*/
	// spacing factors for integral, in case of y not identical to number of 
	// contributions
	normfac=TWO_PISQR/((COMP_PRECISION)(phisumlim * ysumn));
	sqrt2normfac=SQRT_TWO * normfac;

	if(verbose)
	  fprintf(stderr,"%s: Summing up block data, l_max=%i\n",
		  program,lmax);
#define AP atmp[0]
#define BP btmp[0]
#define AT atmp[1]
#define BT btmp[1]
	fprintf(stdout,"%i\n",lmax);
	if(!vectors){
	  /* loop through all m and l */
	  for(warned=l=0;l<=lmax;l++){
	    if(l <= lmax_max){// makes sense to expand
	      for(os1=m=0;m<=l;m++,os1+=nlon){
		AP=BP=0.0;
		/* do the integrations */
		// j would run from 0 to phisumlim but is 
		// used only in combination with os1+j,therefor use this 
		// as limit for j-loop
		os3=phisumlim+os1;
		for(j=os1,os2=0;j < os3;j++,os2+=nlat){ /* phi loop */
		  for(i=0;i < thetasumlim;i++){/* theta loop, sin(theta) 
						  was already mutliplied with P */
		    f[1][0] = P(l,m,i) * (COMP_PRECISION)(*(*func+os2+i));
		    if(m==0)
		      AP += f[1][0];
		    else{
		      AP += f[1][0] * cosarr[j];// this is really os1+j! 
		      BP += f[1][0] * sinarr[j];// same here
		    }
		  }
		}
		if(m==0){
		  AP *= normfac;
		}else{
		  AP *= sqrt2normfac;
		  BP *= sqrt2normfac;
		}
		fprintf(stdout,COEFF_DATA_FORMAT,AP,BP);
		fprintf(stdout,"\n");
	      }
	    }else{// we skip that since lmax too high
	      // warning message follows below
	      for(m=0;m<=l;m++){
		fprintf(stdout,COEFF_DATA_FORMAT,0.0,0.0);
		fprintf(stdout,"\n");
	      }
	    }
	    if(verbose)
	      fprintf(stderr,"%s: l=%i\r",program,l);
	  }
	  if(l > lmax_max+1){
	    fprintf(stderr,"\n%s: WARNING: zeroed out all coefficients above %i\n",
		    program,lmax_max);
	    fprintf(stderr,"%s: WARNING: since we only have %i points in latitude\n",
		    program,nlat);
	  }
	}else{
	  /* 
	     velocity expansion, that is vector harmonics  
	  */
	  if(verbose){
	    fprintf(stderr,"%s: using velocity expansion, spheroidal and toroidal output\n",
		    program);
	    fprintf(stderr,"%s: Dahlen and Tromp convention\n",
		    program);
	  }
#define UPHI(i,j)   (*(*func +                 (j)*nlat + (i)))
#define UTHETA(i,j) (*(*func + nlontimesnlat + (j)*nlat + (i)))
	  /* loop through all m and l */
	  for(l=0;l<=lmax;l++){
	    if(l<=lmax_max){// expand
	      for(os1=m=0;m<=l;m++,os1+=nlon){
		AP=BP=AT=BT=0.0;
		/* do the integrations */
		for(j=0;j < phisumlim;j++){ /* phi loop */	
		  for(i=0;i < thetasumlim;i++){/* theta loop  */
		    /* P and DPTHETA have already been 
		       multiplied by sin(theta) */
		    if(m == 0){
		      /* poloidal */
		      AP +=    DPTHETA(l,m,i) * 
			(COMP_PRECISION)UTHETA(i,j);
		      /* toroidal */
		      AT +=   -DPTHETA(l,m,i) * 
			(COMP_PRECISION)UPHI(i,j);
		    }else{
		      os3=os1+j;
		      cosphi =   *(cosarr+os3) * (COMP_PRECISION)UPHI(i,j);
		      costheta = *(cosarr+os3) * (COMP_PRECISION)UTHETA(i,j);
		      sinphi =   *(sinarr+os3) * (COMP_PRECISION)UPHI(i,j);
		      sintheta = *(sinarr+os3) * (COMP_PRECISION)UTHETA(i,j);
		      /* poloidal */
		      AP +=    DPTHETA(l,m,i) * costheta;
		      AP +=   -P(l,m,i)       * sinphi;
		      BP +=    P(l,m,i)       * cosphi;
		      BP +=    DPTHETA(l,m,i) * sintheta;
		      /* toroidal  */
		      AT +=   -DPTHETA(l,m,i) * cosphi;
		      AT +=   -P(l,m,i)       * sintheta;
		      BT +=    P(l,m,i)       * costheta;
		      BT +=   -DPTHETA(l,m,i) * sinphi;
		    }
		  }
		}
		if(m == 0){
		  AP *= normfac;
		  AT *= normfac;
		}else{
		  AP *= sqrt2normfac;
		  BP *= sqrt2normfac;
		  AT *= sqrt2normfac;
		  BT *= sqrt2normfac;
		}
#ifdef CHECK_FLOATING_EXCEPTION
		if(!myfinite(AP)||!myfinite(BP)||!myfinite(AT)||!myfinite(BT)){
		  fprintf(stderr,"l %i m %i AP %g BP %g AT %g BT %g\n",
			  l,m,AP,BP,AT,BT);
		  exit(-1);
		}
#endif
		fprintf(stdout,"%22.15e %22.15e %22.15e %22.15e\n",
			AP,BP,AT,BT);
	      }
	    }else{// we skip that since lmax to high
	      AP=BP=AT=BT=0.0;
	      for(m=0;m<=l;m++)
		fprintf(stdout,"%22.15e %22.15e %22.15e %22.15e\n",
			AP,BP,AT,BT);
	    }
	    if(verbose)fprintf(stderr,"%s: l=%i\r",program,l);
	  }
	  if(l>lmax_max+1){
	    fprintf(stderr,"\n%s: WARNING: zeroed out all coefficients above %i\n",
		    program,lmax_max);
	    fprintf(stderr,"%s: WARNING: since we only have %i points in latitude\n",
		    program,nlat);
	  }
	}
	if(verbose)fprintf(stderr,"\n");
	break;
      }
	/* 
	   use Gauss integration for theta direction 
	*/
      case GAUSSIAN:{ 
	if(verbose)
	  fprintf(stderr,"%s: Using Gaussian integration in theta...\n",program);
	if(nr_gauss_pts > nlat){
	  /*
	    change this, required because we want to overwrite the func
	    array
	  */
	  fprintf(stderr,"gauss interpolation assumes that nlat >= nr_gauss_pts\n");
	  fprintf(stderr,"now, nlat: %i and nr_gauss_pts: %i\n",
		  nlat, nr_gauss_pts);
	  exit(-1);
	}
	// normalization factors
	normfac=TWO_PI/((COMP_PRECISION)(phisumlim));
	sqrt2normfac=SQRT_TWO*normfac;
	/*
	  interpolate data to lie on nr_gauss_pts (l_max+1) 
	  Gauss points in latitude
	*/
	// arrays for spline interpolation
	mymalloc_cp(&tmp_func,nlat);
	mymalloc_cp(&sndder,nlat);
	for(k=0;// interpolate all values, might be two sets for vel
	    k < nlontimesnlat*((vectors)?(2):(1));
	    k += nlontimesnlat){
	  for(os1=j=0;j<phisumlim;j++,os1+=nlat){// loop through longitude values
	    // interpolate the y values on Gauss points
	    for(i=0;i<nlat;i++)// lat values to be interpolated
	      tmp_func[i]= *(*func+k+os1+i);
	    // linear interpolation
	    /*
	      for(i=0;i<nr_gauss_pts;i++)
	      *(func+k+j*nlat+i)=interpolate(y,tmp_func,absc[i],nlat);
	      */
	    /*
	      splines
	      obtain second derivatives, use natural condition for 
	      boundaries (again, numrec calling style)
	    */
	    spline(y-1,tmp_func-1,nlat,1e30,1e30,sndder-1);
	    for(i=0;i<nr_gauss_pts;i++)
	      /* interpolate tmp_func(y)
		 with sndder second derivatives
		 as func+j*nlat (absc)
		 func gets overwritten
	      */
	      splint(y-1,tmp_func-1,sndder-1,nlat,absc[i],(*func+k+os1+i));
	  }
	}
	free(tmp_func);free(sndder);
	//
	// start output
	//
	fprintf(stdout,"%i\n",lmax);
	if(!vectors){// simple field
	  for(l=0;l <= lmax;l++){
	    for(m=os1=0;m <= l;m++,os1+=nlon){
	      for(AP=BP=0.0,os2=j=0;j<phisumlim;j++,os2+=nlat){ 
		/* 
		   Gauss point integration in theta direction
		   data has been arranged to lie on Gauss points
		   weights are included in Legendre functions
		*/
		f[0][0]=0.0;
		for(i=0;i < nr_gauss_pts;i++)
		  f[0][0] += P(l,m,i)* *(*func+os2+i);
		if(m==0)
		  AP += f[0][0];
		else{
		  AP += f[0][0]* *(cosarr+os1+j);
		  BP += f[0][0]* *(sinarr+os1+j);
		}
	      }
	      if(m==0){
		AP *= normfac;
	      }else{
		AP *= sqrt2normfac;
		BP *= sqrt2normfac;
	      }
	      fprintf(stdout,COEFF_DATA_FORMAT,AP,BP);
	      fprintf(stdout,"\n");
	    }
	    if(verbose)
	      fprintf(stderr,"%s: l=%5i\r",program,l);
	  }
	  if(verbose)
	    fprintf(stderr,"\n");
	}else{// vector field
	  for(l=0;l<=lmax;l++){
	    for(os1=m=0;m<=l;m++,os1+=nlon){
	      AP=BP=AT=BT=0.0;
	      for(j=0;j<phisumlim;j++){ 
		for(i=0;i<nr_gauss_pts;i++){
		  if(m == 0){
		    AP +=    DPTHETA(l,m,i) * (COMP_PRECISION)UTHETA(i,j);
		    AT +=   -DPTHETA(l,m,i) * (COMP_PRECISION)UPHI(i,j);
		  }else{
		    os3=os1+j;
		    cosphi =   *(cosarr+os3) * (COMP_PRECISION)UPHI(i,j);
		    costheta = *(cosarr+os3) * (COMP_PRECISION)UTHETA(i,j);
		    sinphi =   *(sinarr+os3) * (COMP_PRECISION)UPHI(i,j);
		    sintheta = *(sinarr+os3) * (COMP_PRECISION)UTHETA(i,j);
		    AP +=    DPTHETA(l,m,i) * costheta;
		    AP +=   -P(l,m,i)       * sinphi;
		    BP +=    P(l,m,i)       * cosphi;
		    BP +=    DPTHETA(l,m,i) * sintheta;
		    AT +=   -DPTHETA(l,m,i) * cosphi;
		    AT +=   -P(l,m,i)       * sintheta;
		    BT +=    P(l,m,i)       * costheta;
		    BT +=   -DPTHETA(l,m,i) * sinphi;
		  }
		}
	      }
	      if(m == 0){
		AP *= normfac;AT *= normfac;
	      }else{
		AP *= sqrt2normfac;BP *= sqrt2normfac;
		AT *= sqrt2normfac;BT *= sqrt2normfac;
	      }
	      fprintf(stdout,"%22.15e %22.15e %22.15e %22.15e\n",
		      AP,BP,AT,BT);
	    }
	    if(verbose)
	      fprintf(stderr,"%s: l=%i\r",program,l);
	  }
	  if(verbose)
	    fprintf(stderr,"\n");
	}
	break;
      }
      default:{
	fprintf(stderr,"%s: block mode: integration mode %i not implemented.\n",program, intmode);
	phelp(program);exit(-1);
	break;
      }}
      break;
    }
    case POINTS_FOR_AB:
      /* use x/y locations to generate SHA of pointwise data  */
    case CONTOUR:{
      /* use data x/y/z to integrate along a contour 
	 using delta functions
      */

      if(verbose)fprintf(stderr,"%s: Initializing...\n",program);
      /* Legendre functions */
      /* check if space needed smaller than MEMLIMIT */
      memneeded=lmsize*nrp*sizeof(COMP_PRECISION)/ONE_MEGABYTE;
      if(memneeded > MEMLIMIT){
	lookup_p=0;
	if(verbose)
	  fprintf(stderr,"%s: using individual function calls for Assoc. Legendre pol.\n",
		  program);
      }else{
	if(verbose)fprintf(stderr,"%s: constructing Assoc. legendre functions...\n",program);
	myrealloc_cp(&p,lmsize*nrp);zero_cp(p,lmsize*nrp);
	plgndr2(cloc,nrp,p,lmax);
      }
      /* 
	 start output 
      */
      if(verbose){
	if(opmode==CONTOUR)
	  fprintf(stderr,"%s: Using Dirac function integration along scalar contours\n",program);
	else
	  fprintf(stderr,"%s: Using Dirac functions on x/y/z data\n",program);
      }
      fprintf(stdout,"%i\n",lmax);
      if(lookup_p){/* use look up table for the Legendre functions */
	for(l=0;l<=lmax;l++){
	  for(m=0;m<=l;m++){
	    ctmp[0] = ctmp[1] = 0.0;
	    if(nrp == 1){
	      f[0][0] = P(l,m,0) * (COMP_PRECISION)((*func)[0]);
	      if(m==0)
		ctmp[0] += f[0][0];
	      else{
		tmp = SQRT_TWO * f[0][0];
		tmp2 = (COMP_PRECISION)m*(COMP_PRECISION)cloc[0];
		ctmp[0] += tmp * cos(tmp2); 
		ctmp[1] += tmp * sin(tmp2);
	      }
	    }else{
	      /* 
		 first function value for a and b, 
		 f[first/second value][a or b] 
		 Legendre * function * { cos \atop sin} 
	      */
	      i=0;// if we just want the SHA of points
	      tmp= P(l,m,i)*(COMP_PRECISION)((*func)[i]);
	      tmp2 = (COMP_PRECISION)m*(COMP_PRECISION)(*(cloc));
	      f[0][0] = tmp * cos(tmp2);
	      f[0][1] = tmp * sin(tmp2);
	      if(opmode != CONTOUR){/* add first term to sum 
				       if we are not integrating */
		tmp=(m==0)?(1.0):(SQRT_TWO);
		for(coeff=0;coeff<2;coeff++)
		  ctmp[coeff] += tmp * f[0][coeff];
	      }
	      for(os1=2,i=1;i<nrp;i++,os1+=2){/* integrate over points 
						 along a contour */
		tmp= P(l,m,i)*(COMP_PRECISION)((*func)[i]);
		/* second function value for a and b */
		tmp2 = (COMP_PRECISION)m* 
		  (COMP_PRECISION)(*(cloc+os1));
		f[1][0] = tmp * cos(tmp2);
		f[1][1] = tmp * sin(tmp2);
		if(opmode == CONTOUR){/* we want to integrate 
					 over a contour */
		  /* 
		     this gives the scaling factor depending 
		     on distance between points, divides it by 
		     two for the trapezoidal rule 
		     and multiplies by sqrt(2) for m != 0
		     for the normalization of the coefficients
		  */
		  if(m==0)
		    dist= dist_rad(cloc,i-1,i) * 0.5;
		  else
		    dist= SQRT_TWO * dist_rad(cloc,i-1,i) * 0.5;
		  for(coeff=0;coeff<2;coeff++){/* we can't pull 
						  the averaging 
						  out of this loop 
						  since dist changes! */
		    /* sum contributions */
		    ctmp[coeff] += dist * (f[0][coeff]+f[1][coeff]); 
		    /* make the second function value the 
		       first for integration */
		    f[0][coeff]=f[1][coeff];
		  }
		}else{/* simply add up contributions of all 
			 points to this particular mode */
		  tmp=(m==0)?(1.0):(SQRT_TWO);
		  for(coeff=0;coeff<2;coeff++)
		    ctmp[coeff] += tmp * f[1][coeff];
		}
	      }
	    }
	    // write coefficients for this mode
	    fprintf(stdout,COEFF_DATA_FORMAT,ctmp[0],ctmp[1]);
	    fprintf(stdout,"\n");
	  }
	  if(verbose)fprintf(stderr,"%s: l=%i\r",program,l);
	}
      }else{
	/* 
	   construct Legendre function for every single 
	   point individually 
	*/
	for(l=0;l<=lmax;l++){
	  for(m=0;m<=l;m++){
	    ctmp[0]=ctmp[1]=0.0;
	    if(nrp==1){
	      f[0][0] = slgndr(l,m,cos((COMP_PRECISION)cloc[1])) * (COMP_PRECISION)((*func)[0]);
	      if(m==0)
		ctmp[0] += f[0][0];
	      else{
		tmp=SQRT_TWO * f[0][0];
		tmp2 = (COMP_PRECISION)m* (COMP_PRECISION)(*(cloc));
		ctmp[0] += tmp * cos(tmp2); 
		ctmp[1] += tmp * sin(tmp2);
	      }
	    } else {
	      tmp=slgndr(l,m,cos((COMP_PRECISION)cloc[1]));
	      tmp2 = (COMP_PRECISION)m*(COMP_PRECISION)(*(cloc)); 
	      f[0][0] = tmp*(COMP_PRECISION)((*func)[i])*cos(tmp2);
	      f[0][1] = tmp*(COMP_PRECISION)((*func)[i])*sin(tmp2);
	      if(opmode != CONTOUR){// add first term to sum if we are not integrating
		tmp=(m==0)?(1.0):(SQRT_TWO);
		for(coeff=0;coeff < 2;coeff++)
		  ctmp[coeff] += tmp * f[0][coeff];
	      }
	      for(os1=2,i=1;i<nrp;i++,os1+=2){/* integrate over points along a contour */
		tmp=slgndr(l,m,cos((COMP_PRECISION)cloc[i*2+1]));
		/* second function value for a and b */
		tmp2 = (COMP_PRECISION)m* (COMP_PRECISION)(*(cloc+os1));
		f[1][0] = tmp * (COMP_PRECISION)((*func)[i])*cos(tmp2);
		f[1][1] = tmp * (COMP_PRECISION)((*func)[i])*sin(tmp2);
		if(opmode == CONTOUR){
		  if(m==0)
		    dist= dist_rad(cloc,i-1,i)*0.5;
		  else
		    dist= SQRT_TWO * dist_rad(cloc,i-1,i) * 0.5;
		  for(coeff=0;coeff < 2;coeff++){
		    ctmp[coeff] += dist * (f[0][coeff]+f[1][coeff]); 
		    f[0][coeff]=f[1][coeff];
		  }
		}else{
		  tmp=(m==0)?(1.0):(SQRT_TWO);
		  for(coeff=0;coeff<2;coeff++)
		    ctmp[coeff] += tmp * f[1][coeff];
		}
	      }
	    }
	    fprintf(stdout,COEFF_DATA_FORMAT,ctmp[0],ctmp[1]);
	    fprintf(stdout,"\n");
	  }
	}
	if(verbose)fprintf(stderr,"\n");
      }
      break;
    }
      
      /* 
	 use point data to write A matrix for least 
	 square solution. matrix has A_00 A_10 A_11 B_11 A_20 A_21 B_21 A_22 B_22 ... as row i
	 at the location of data point i. colums are thus given by (lmax+1)(lmax+2)

	 A matrix is C style

      */
    case POINTS_FOR_AVEC_MATRIX_ASCII:
    case POINTS_FOR_AVEC_MATRIX:
    case POINTS_FOR_A_MATRIX_ASCII:
    case POINTS_FOR_A_MATRIX:{
      if((opmode == POINTS_FOR_AVEC_MATRIX_ASCII) || (opmode == POINTS_FOR_AVEC_MATRIX))
	vectors = 1;
      else
	vectors = 0;
      /* 
	 c style matrix, don't divide by sigma 
      */
      make_a(lmax,lmsize,nrp,&aflt,cloc,verbose,program,0,nsize,0.0,FALSE,FALSE,sigma,vectors);
      lmc = nsize[1];	/* number of non-zero paramteres */
      if(vectors)
	nrp2 = nrp * 2;
      else
	nrp2 = nrp;
      fprintf(stderr,"%s: (%i points, %i data (vector: %i) and %i different, non-zero coefficients)\n",
	      program,nrp,nrp2,vectors,lmc);
      if((opmode == POINTS_FOR_A_MATRIX_ASCII)||(opmode ==  POINTS_FOR_AVEC_MATRIX_ASCII)){
	for(i=0;i < nrp2;i++){
	  for(j=0;j < lmc;j++){	/* print C style */
	    fprintf(stdout,"%12.5e ",aflt[i*lmc+j]);
	  }
	  fprintf(stdout,"\n");
	}
      }else{
	fprintf(stderr,"%s: binary single precision FORTRAN style output will be %i rows by %i cols\n",
		program,nrp,lmc);
	
	for(j=0;j < lmc;j++)	/* print */
	  for(i=0,os1=j;i < nrp2;i++,os1+=lmc)
	    fwrite((aflt+os1),sizeof(float),(size_t)1,stdout);
      }
      free(aflt);
      break;
    }
    default:{
      fprintf(stderr,"%s: unknown opmode %i\n",program,opmode);
      exit(-1);
    }}
  }
  if(y)
    free(y);



  if(p)
    free(p);

  if(dptheta)
    free(dptheta);

  if(cosarr)
    free(cosarr);

  if(sinarr)
    free(sinarr);

  if(x)
    free(x);
}


/* 

   make an A matrix for least squares solution 

   pass a as NULL on init 

   
   asize[0] = m = number of rows    = number of data ( x 2 for vectors) + (damping ? number of parameters : 0)
   asize[1] = n = number of columns = number of parameters
   asize[2] = n * m

   damping mode: 0: no damping
                 1: norm damping constant set to damping
		 2: l/L scaled

   sigma is a weighting function
		 
*/


void make_a(int lmax, int lmsize, int nrp, float **a, DATA_PRECISION *cloc,int verbose, 
	    char *program,int damping_mode, int *asize, COMP_PRECISION damping,
	    int fortran_style,int use_sigma, 
	    float *sigma, int vectors)
{
  int lmax1,lmc;
  int i,j,k,klim,l,m,nrp2,ip;
  COMP_PRECISION tmp2,sinmphi,cosmphi,tmp_a,tmp_b,loc_sigma;
  COMP_PRECISION *p,*dptheta,*dp2,scale,*y,sig_l_factor,one_over_sin_theta;
  lmax1 = lmax + 1;
  
 
  /* 
     A(m,n) 

     m = number of rows    = number of data (x 2 for vectors ) + (damping ? number of parameters : 0)
     n = number of columns = number of parameters (lmc)

  */
  /* number of non-zero parameters */
  lmc = lmsize*2 - lmax1; 
  if(vectors){
    lmc--;			/* no L = 0 term */
    lmc *= 2;			/* two sets */
    klim = 2;
    nrp2 = nrp * 2;
  }else{
    klim = 1;
    nrp2 = nrp;
  }
  asize[0] = nrp2 + ((damping_mode > 0)?(lmc):(0));
  asize[1] = lmc;
  asize[2] = asize[1] * asize[0];
  if(verbose)
    fprintf(stderr,"%s: A matrix: %i data locations (%i values), %i free parameters, m: %i n: %i, dmode %i dfac %g\n",
	    program,nrp,nrp2,lmc,asize[0],asize[1],damping_mode,damping);
  if(verbose)fprintf(stderr,"%s: Initializing A matrix...\n",program);
  myrealloc_dp(a,asize[2]);
  for(i=0;i < asize[2];i++)
    *(*a+i) = 0.0;			/* init as zero */
  
#define LOC_AIJ(i,j) (*(*a + aij(i,j,asize[0],asize[1],fortran_style)))

  /* legendre functions */
  if(verbose)fprintf(stderr,"%s: constructing Assoc. legendre functions...lmsize %i nrp %i\n",
		     program,lmsize,nrp);
  mycalloc_cp(&p,lmsize * nrp);
  /* 

     get the legendre functions at the desired latitudes

  */
  plgndr2(cloc,nrp,p,lmax);
  if(vectors){
    /* get derivatives */
    mycalloc_cp(&y,nrp);
    for(i=0;i < nrp;i++)
      y[i] = cos(cloc[i*2+1]);
    /* 
       compute 
    */
    mycalloc_cp(&dptheta,lmsize * nrp);
    pdtheta_lgndr(y,nrp,p,dptheta,lmax,dp2,0);
    free(y);
  }else{
    dp2 = NULL;
  }

  /* 
     use look up table for the Legendre functions,
     assemble matrix in memory
  */

  for(ip=0;ip < nrp;ip++){		    /* loop through points */
    one_over_sin_theta = 1/sin(cloc[ip*2+1]); /* 1/sin(theta) */
    for(k=0;k < klim;k++){	/* for vectors, will loop; else only
				   executed once */
      i = ip + k  * nrp;	/* data counter, rows of matrix */
      if(use_sigma)
	loc_sigma = sigma[i];
      else
	loc_sigma = 1.0;
#ifdef DEBUG
      if(fabs(loc_sigma) < 1e-5){
	fprintf(stderr,"error, cannot use zero sigma weighting\n");
	exit(-1);
      }
#endif
      j = 0;			/* columns of matrix */
      for(l=0;l <= lmax;l++){
	if(vectors){
	  if(l > 0){
	    /* normalization factor */
	    sig_l_factor = (1/sqrt((COMP_PRECISION)(l*(l+1))));
	    sig_l_factor /= loc_sigma; /* use sigma here */

	    for(m=0;m <= l;m++){
	      tmp_a = (COMP_PRECISION)m * P(l,m,ip)*sig_l_factor*one_over_sin_theta;
	      tmp_b = DPTHETA(l,m,ip)*sig_l_factor;
	      if(m == 0){
		if(k == 0){
		  /* Blm phi */
		  LOC_AIJ(i,j) = (float)(tmp_a);j++;
		  /* Clm phi */
		  LOC_AIJ(i,j) = (float)(-tmp_b);j++;
		}else{
		  /* Blm theta */
		  LOC_AIJ(i,j) = (float)(tmp_b);j++;
		  /* Clm theta */
		  LOC_AIJ(i,j) = (float)(tmp_a);j++;
		}
	      }else{
		tmp2 = (COMP_PRECISION)m* (COMP_PRECISION)(*(cloc+ip*2)); /* m * phi */
		cosmphi = cos(tmp2);
		sinmphi = sin(tmp2);
		tmp_a = (COMP_PRECISION)m * P(l,m,ip)*sig_l_factor*one_over_sin_theta;
		tmp_b = DPTHETA(l,m,ip)*sig_l_factor;

		if(k == 0){
		  /* Blm phi */
		  LOC_AIJ(i,j) = (float)(tmp_a * cosmphi);j++;
		  LOC_AIJ(i,j) = (float)(tmp_a * sinmphi);j++;
		  /* Clm phi */
		  LOC_AIJ(i,j) = (float)(-tmp_b * cosmphi);j++;
		  LOC_AIJ(i,j) = (float)(-tmp_b * sinmphi);j++;
		}else{
		  /* Blm theta */
		  LOC_AIJ(i,j) = (float)(tmp_b * cosmphi);j++;
		  LOC_AIJ(i,j) = (float)(tmp_b * sinmphi);j++;
		  /* Clm theta */
		  LOC_AIJ(i,j) = (float)(tmp_a * cosmphi);j++;
		  LOC_AIJ(i,j) = (float)(tmp_a * sinmphi);j++;
		}
	      }
	    }
	  }
	}else{			/* scalar */
	  for(m=0;m <= l;m++){
	    if(m == 0){
	      LOC_AIJ(i,j) = (float)(P(l,m,ip))/loc_sigma;j++;
	    }else{
	      /* m times phi location of data i */
	      tmp2 = (COMP_PRECISION)m* (COMP_PRECISION)(*(cloc+ip*2)); /* m * phi */
	      LOC_AIJ(i,j) = ((float)(SQRT_TWO * P(l,m,ip) * cos(tmp2))/loc_sigma);j++;
	      LOC_AIJ(i,j) = ((float)(SQRT_TWO * P(l,m,ip) * sin(tmp2))/loc_sigma);j++;
	    }
	  }
	}
      }	/* l loop */
    } /* k loop for vectors */
  }   /* end point loop */

  if(verbose)fprintf(stderr,"%s: A matrix done.\n",program);

  free(p);
  if(vectors)
    free(dptheta);
    
  for(i=nrp2;i < asize[0];i++){
    j=0;
    for(l=0;l <= lmax;l++){
      switch(damping_mode){
      case 1:			/* norm damping */
	scale = damping;
	break;
      case 2:			/* roughness */
	scale = damping * (double)l/(double)lmax;
	break;
      default:
	fprintf(stderr,"damp mode %i undefined\n",damping_mode);
	exit(-1);
      }
      for(m=0;m <= l;m++){
	if(m == 0){
	  LOC_AIJ(i,j)   = ((i-nrp2) == j)?(scale):(0.0);j++;
	  if(vectors){
	    LOC_AIJ(i,j) = ((i-nrp2) == j)?(scale):(0.0);j++;
	  }
	}else{
	  LOC_AIJ(i,j)   = ((i-nrp2) == j)?(scale):(0.0);j++;
	  LOC_AIJ(i,j)   = ((i-nrp2) == j)?(scale):(0.0);j++;
	  if(vectors){
	    LOC_AIJ(i,j) = ((i-nrp2) == j)?(scale):(0.0);j++;
	    LOC_AIJ(i,j) = ((i-nrp2) == j)?(scale):(0.0);j++;
	  }
	}
      }
    }
  }
#undef LOC_AIJ

}
/* 
   get offset for A(M,N) matrix 

   SLOW!!!
*/
int aij(int i, int j, int m, int n, int fortran_style)
{
#ifdef DEBUG
  if((i >= m)||(j>=n)){fprintf(stderr,"aij: out of bounds i: %i m: %i j: %i n: %i\n",i,m,j,n);exit(-1);}
#endif
  if(fortran_style){
    return j * m + i;
  }else{			/* C style */
    return i * n + j;
  }
}

/* linear interpolation function */
COMP_PRECISION interpolate(COMP_PRECISION *x,DATA_PRECISION *y,COMP_PRECISION x0,int n)
{
  int j,k;
  COMP_PRECISION ifunc,a,b;
  /* 
     find the bordering data points 
  */
  j=k=0;
  while(k<n && x[k]<=x0)
    k++;
  if(k==n){
    k--;j=k-1;
  }else if(k==0){
    j=0;k=1;
  }else{
    j=k-1;
  }
  // linearly interpolate
  a=(x[k]-x0)/(x[k]-x[j]);
  b=1.0-a;
  ifunc =a*(COMP_PRECISION)y[j];// f[y[j]]
  ifunc+=b*(COMP_PRECISION)y[k];// f[y[k]]=f[y[j+1]]
  return ifunc;
}
/* 
   calculate the distance on the sphere based on my theta and phi 
   coordinates of the contour data points

   returns distance in radians 

*/
DATA_PRECISION dist_rad(DATA_PRECISION *cloc, int i, int j)
{
  COMP_PRECISION tmp1,tmp2,tmp3;
  int os1,os2;
  os1 = i * 2;os2 =j * 2;
  tmp1  = (cloc[os2+1]-cloc[os1+1])/2.0;
  tmp1  = sin(tmp1);
  tmp1  = tmp1 * tmp1;
  tmp2  = (cloc[os1]-cloc[os2])/2.0;
  tmp2  = sin(tmp2);
  tmp2  = tmp2 * tmp2;
  tmp3  = tmp1;
  tmp3 += sin(cloc[os1+1]) * sin(cloc[os2+1]) * tmp2;
  tmp3 = sqrt(tmp3);
  return (DATA_PRECISION)(2.0*asin(tmp3));
}
/* since we allow for grd filenames, have to do some funny
   checking for numerical opmodes
*/
void check_opmode(char *argument, int *op_mode, char *filename,
		  int *vectors, int *calculate_derivatives)
{
  char a[3],b[3],c[3],d[3],f[3],g[3],h[3],i[3],j[3];
  sprintf(a,"%1i",BLOCK);
  sprintf(b,"%1i",CONTOUR);
  sprintf(c,"%1i",ASCII_BLOCK);
  sprintf(d,"%1i",ASCII_BLOCK_HEADER);
  sprintf(f,"%1i",POINTS_FOR_A_MATRIX);
  sprintf(g,"%1i",POINTS_FOR_AB);
  sprintf(h,"%1i",POINTS_FOR_A_MATRIX_ASCII);
  sprintf(i,"%1i",POINTS_FOR_AVEC_MATRIX);
  sprintf(j,"%1i",POINTS_FOR_AVEC_MATRIX_ASCII);
 
  if(strcmp(argument,a) == 0)
    *op_mode=BLOCK;
  else if(strcmp(argument,b) == 0)
    *op_mode=CONTOUR;
  else if(strcmp(argument,c) == 0)
    *op_mode=ASCII_BLOCK;
  else if(strcmp(argument,d) == 0)
    *op_mode=ASCII_BLOCK_HEADER;
  else if(strcmp(argument,f) == 0)
    *op_mode=POINTS_FOR_A_MATRIX;
  else if(strcmp(argument,g) == 0)
    *op_mode=POINTS_FOR_AB;
  else if(strcmp(argument,h) == 0)
    *op_mode=POINTS_FOR_A_MATRIX_ASCII;
  else if(strcmp(argument,i) == 0)
    *op_mode=POINTS_FOR_AVEC_MATRIX;
  else if(strcmp(argument,j) == 0)
    *op_mode=POINTS_FOR_AVEC_MATRIX_ASCII;
  else{
    sprintf(filename,"%s.grd",argument);
    *op_mode=GMT_BLOCK;
    if(strcmp(argument,FST_VEC_FILE_START) == 0){
      sprintf((filename+STRING_LENGTH),SCD_VEC_FILE);
      *vectors=TRUE;
      *calculate_derivatives=TRUE;
    }
  }
}
/*
  
  scale, rotate, and check

*/

/* 
   version for GMT4, using array flip 
*/
void gmt2myconvention_rotate4(GMT_PRECISION *func,int nlon,int nlat,
			      GMT_PRECISION factor)
{
  int i,j,warned=0,os1,os2,os3,nm;
  GMT_PRECISION *array;
  mymalloc_dp(&array,nlon * nlat);
  // scale
  nm = nlon*nlat;
  for(i=0;i<nm;i++){
    if(!myfinite(func[i])){
      func[i]=0.0;
      if(!warned){
	fprintf(stderr,"WARNING: at least one NaN entry in the data has been replaced with zero\n");
	warned=1;
      }
    }
    array[i] =  func[i]*factor;
  }
  // rotate
  os3 = (nlat-1) * nlon;
  for(os1=j=0;j<nlon;j++,os1 += nlat){
    for(os2=os3,i=0;i<nlat;i++,os2 -= nlon)
      *(func + os1 + i)= *(array + os2 + j);
  }
  free(array);
}

#ifndef USE_GMT4
/* version for GMT > 4 */
void gmt2myconvention_rotate(GMT_PRECISION *func,int nlon,int nlat,
			     GMT_PRECISION factor, struct GMT_GRID *G)
{
  int i,j,warned=0,ij,ijm;
  //float mean = 0;
  for(j=0;j < nlon;j++){
    for(i=0;i < nlat;i++){
      ij = gmt_M_ijp (G->header, i, j); /* row, col */
      //ijm = j*nlat+i;
      ijm = j*nlat+nlat-1-i;		 /* convert from top down to
					    bottom up */
      func[ijm] =  G->data[ij] * factor; /* assign and scale */
      //mean +=  G->data[ij];
      if(!myfinite(func[ijm])){
	func[ijm] = 0.0;
	if(!warned){
	  fprintf(stderr,"WARNING: at least one NaN entry in the data has been replaced with zero\n");
	  warned=1;
	}
      }
    }
  }
  //fprintf(stderr,"converted %i by %i data, mean: %11g\n",nlon,nlat,mean/(float)nlon*nlat);
}
#endif

void phelp(char *name)
{
  fprintf(stderr,"\n%s [l_max, %i] [input mode, %1i] [integration mode, %i] [damping, %g]\n",
	  name,DEF_LMAX,DEF_OPMODE,DEF_INT_MODE,DEF_DAMPING);
  fprintf(stderr,"\t calculates spherical harmonic expansion for norm=1\n\t real theoretical physics spherical harmonic coefficients\n");
  fprintf(stderr,"\t as in Dahlen and Tromp (1998)\n\n");
  fprintf(stderr,"\t l_max:      maximum order of expansion\n");
  fprintf(stderr,"\t input mode: choose different forms of possible input\n\n");
  
  fprintf(stderr,"\t             name<.grd>: read from GMT grd-file\n");
  fprintf(stderr,"\t                         if name is %s, look for %s.grd and %s to\n",
	  FST_VEC_FILE_START,FST_VEC_FILE_START,SCD_VEC_FILE);
  fprintf(stderr,"\t                         expand vector field with phi theta components\n\n");
  fprintf(stderr,"\t             something like \"shana 31 mygrid\" where mygrid.grd is a netcdf/GMT grid file\n");
  fprintf(stderr,"\t             is likely what you want\n\n");

  fprintf(stderr,"\t             %i: read an ASCII file with points (lon, lat, d) for\n",CONTOUR);
  fprintf(stderr,"\t                 contour delta function integration (e.g. for slabs)\n");
  fprintf(stderr,"\t             %i: read an ASCII file with spotted points (lon, lat, z) for A matrix output\n",
	  POINTS_FOR_A_MATRIX);
  fprintf(stderr,"\t                 where A is m by n with n the spherical harmonics coefficients evaluated at m data locations (2m for vectors, u_phi and u_theta)\n");
  fprintf(stderr,"\t                 (A matrix will be in binary single precision, FORTRAN style), e.g. for later least squares inversion\n");
  fprintf(stderr,"\t             %i: same as %i, but ascii output, regular, C style\n",
	  POINTS_FOR_A_MATRIX_ASCII,POINTS_FOR_A_MATRIX);
  fprintf(stderr,"\t             %i: same as %i but output for poloidal/toroidal vector field expansion\n",
	  POINTS_FOR_AVEC_MATRIX,POINTS_FOR_A_MATRIX);
  fprintf(stderr,"\t             %i: same as %i but output for poloidal/toroidal vector field expansion\n\n",
	  POINTS_FOR_AVEC_MATRIX_ASCII,POINTS_FOR_A_MATRIX_ASCII);

  fprintf(stderr,"\t             %i: read an ASCII file with spotted points (lon, lat, z) for spotted SH analysis with delta functions\n",
	  POINTS_FOR_AB);
  fprintf(stderr,"\n\n\t             If the interpolation mode is %i, modes %i, %i, and %i will fit regular SH with least squares\n",
	  LEAST_SQUARES,CONTOUR,POINTS_FOR_A_MATRIX,POINTS_FOR_AB);
  fprintf(stderr,"\t             Quite likely, you will not want to use these modes else, as they only make sense for special cases.\n\n\n");



  fprintf(stderr,"\t             %i: read a binary block of n times m data values\n",BLOCK);
  fprintf(stderr,"\t                 (%i bytes/value) with dimensions\n",
	  (int)sizeof(DATA_PRECISION));
  fprintf(stderr,"\t                 as leading integers n and m (%i bytes/value))\n",(int)sizeof(int));
  fprintf(stderr,"\t                 the input is assumed to be dx/2<=x<=360-dx/2 and -90+dy/2<=y<=90-dy/2\n");
  fprintf(stderr,"\t                 Regular spacing in lon and lat is assumed\n");

  fprintf(stderr,"\t             %i: read an ASCII block with nlon data columns in nlat rows\n",ASCII_BLOCK);
  fprintf(stderr,"\t                 expects nlon=2*nlat, size determined from # of points\n");
  fprintf(stderr,"\t                 the input is assumed to be dx/2<=x<=360-dx/2 and -90+dy/2<=y<=90-dy/2\n");

  fprintf(stderr,"\t             %i: read an ascii block with header line (-%i for latitude flip)\n",
	  ASCII_BLOCK_HEADER,ASCII_BLOCK_HEADER);
  fprintf(stderr,"\t                 header has format:\n");
  fprintf(stderr,"\t                 nlon nlat lon_min lon_max lat_min lat_max\n\n");
  fprintf(stderr,"\t int. mode:  choose integration mode for block data (theta dir.)\n");
  fprintf(stderr,"\t             %i: sum up block data (\"trapezoidal\")\n",
	  TRAPEZOIDAL);
  fprintf(stderr,"\t             %i: Gaussian quadrature for latitudes\n",GAUSSIAN);
  fprintf(stderr,"\t             %i: least squares fit (works for input mode: %i, %i, %i (same), and GMT grid, interpolated)\n",
	  LEAST_SQUARES,CONTOUR,POINTS_FOR_A_MATRIX,POINTS_FOR_AB);
  fprintf(stderr,"\n\tdamping:     damping value for least-squares fit (%g), if negative use roughness, else normal\n",DEF_DAMPING);
  fprintf(stderr,"\n\t             lon and lat should be given in degrees\n");
  fprintf(stderr,"\t             and the code expects data in roughly 0-360,-90 to 90 format\n\n");
  fprintf(stderr,"EXAMPLE USAGE:\n\n");
  fprintf(stderr,"1) shana 100 topo\n\nexpand spatial data in netcdf topo.grd up to l_max = 100 by means of Gauss integration and print coefficnets to stdout\nuse abconvert to convert coefficients into other formats, compute power, etc.\n\n");
  fprintf(stderr,"2) shana 10 topo 3\n\nexpand spatial data in netcdf topo.grd up to l_max = 10 by means of least-squares inversion (requires LAPACK support compiled in)\n\n");
  fprintf(stderr,"3) shana 10 topo | shsyn 1 0 tmp.10\n\nexpand spatial data in topo.grd up to l_max = 10 and write a 1 deg resolution grid to tmp.10.grd\n\n");
  fprintf(stderr,"4) shana 10 1 3 -0.5 < tmp.lonlatz > tmp.ab\n\nread irregularly distributed lon lat value data in tmp.lonlatz and best-fit a l_max = 10 expansion with roughness damping of 0.5\n\n");
}


void minmax(DATA_PRECISION *min, DATA_PRECISION *max, DATA_PRECISION *avg, 
	    DATA_PRECISION *func, int size)
{
  int i;
  *min= *func;
  *max= *func;
  *avg= *func;
  for(i=1;i<size;i++){
    if(func[i] > *max)
      *max= func[i];
    if(*(func+i)< *min)
      *min= func[i];
    *avg += func[i];
  }
  if(i>0)
    *avg /= (DATA_PRECISION)i;
}
void field_message(DATA_PRECISION minphi,DATA_PRECISION dphi, DATA_PRECISION maxphi,
		   DATA_PRECISION mintheta,DATA_PRECISION dtheta, DATA_PRECISION maxtheta,
		   int nlon, int nlat,
		   DATA_PRECISION min,DATA_PRECISION max, DATA_PRECISION avg,
		   char *program)
{
  fprintf(stderr,"%s: input data field size: %7i by %7i\n",program,nlon,nlat);
  fprintf(stderr,"%s:             min:    dx/avg:      max\n",program);
  fprintf(stderr,"%s:    x:%10g:%10g:%10g\n",
	  program,PHI2LONGITUDE(minphi), RAD2DEG(dphi),PHI2LONGITUDE(maxphi));
  fprintf(stderr,"%s:    y:%10g:%10g:%10g\n",program,
	  THETA2LATITUDE(mintheta),RAD2DEG(dtheta),THETA2LATITUDE(maxtheta));
  fprintf(stderr,"%s:    z:%10g:%10g:%10g\n",
	  program,min,avg,max);
}

#ifdef USE_LAPACK

/* 

   solve A(m,n).x(n) = b(m) using LAPACK linear least squares 

   b needs to be dimensioned MAX(n,m) since it will be x on output

*/
//#define DEBUG_SOLVER
void solver_ab_lls(float *a, int m, int n, float *b,char *program)
{
  float *work;
  int info,lda,ldb,mn,lwork,nrhs=1,nb=256; /* nb is the blocksize */
  char trans='N';
#ifdef DEBUG_SOLVER
  FILE *out;
  int i,j;
#endif

  lda = m;
  ldb = MAX(m,n);
  mn = MIN(n,m);
  lwork =  2*(mn + MAX(mn, nrhs )*nb);

  mymalloc_dp(&work,lwork);
  fprintf(stderr,"%s: solver_ab_lls: m: %i n: %i (lwork: %i)\n",program,m,n,lwork);  
  /* LAPACK least squares */
  sgels_(&trans,&m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
  if(info != 0){
    fprintf(stderr,"%s: solver_ab_lls: LAPACK returned error code %i for m: %i n: %i\n",
	    program,info,m,n);
    exit(-1);
  }

  free(work);
}
#endif
