#include "shsyn.h"
//
// grd output part used by shsyn
//
// $Id: mygrdio.c,v 1.4 2009/04/16 00:55:29 becker Exp becker $
//
void grid_output(int out_mode,char *grdfilename,
		 GMT_PRECISION *gmtval,
		 int nlon,int nlat,
		 COMP_PRECISION xmin,
		 COMP_PRECISION xmax,
		 COMP_PRECISION ymin,
		 COMP_PRECISION ymax,
		 COMP_PRECISION dx,COMP_PRECISION dy,
		 int argc,char **argv,
		 int lmax,BOOLEAN verbose)
{

  switch(out_mode){
  case ONE_GRD:{
    my_gmt_write_grd(gmtval, verbose, argc,argv, 
		     grdfilename, 
		     nlon,nlat, xmin, xmax,
		     ymin, ymax, dx,  dy);
    break;
  }
  case BINARY_STDOUT:{ /* output to stdout as a binary block */
    if(verbose)fprintf(stderr,"%s: writing binary %i by %i grid with %i bytes/value to stdout\n",
		       argv[0],nlon,nlat,
		       (int)sizeof(GMT_PRECISION));
    if(fwrite(gmtval,sizeof(GMT_PRECISION),
	      nlon*nlat,stdout)!=nlon*nlat){
      fprintf(stderr,"%s: write error\n",argv[0]);
    }
    break;
  }
  default:{
    fprintf(stderr,"output mode mix up\n");exit(-1);
  }}
}

/* 
   write a GMT grd file using a GMT routine 
   which calls netCDF programs 

*/

void my_gmt_write_grd(float *phival, BOOLEAN verbose, 
		      int argc,char **argv, char *grdfilename, 
		      int nlon, int nlat,
		      COMP_PRECISION xmin, COMP_PRECISION xmax,
		      COMP_PRECISION ymin,COMP_PRECISION ymax,
		      COMP_PRECISION dx, COMP_PRECISION dy)
{
  struct GRD_HEADER grd;
  // initialize GMT grd file
  double wesn[4];
  GMT_grd_init (&grd, argc, argv, FALSE);
  wesn[0] = xmin;wesn[1] = xmax;wesn[2] = ymin;wesn[3] = ymax;
  grd.node_offset = 0;
  grd.x_min = xmin;	grd.x_max = xmax;
  grd.y_min = ymin;	grd.y_max = ymax;
  grd.x_inc = dx;	grd.y_inc = dy;

 
  grd.nx = nlon;//irint ((xmax-xmin)/dx)+ 1;
  grd.ny = nlat;//irint ((ymax-ymin)/dy)+ 1;
#ifndef USE_GMT3
  GMT_io.in_col_type[0] = GMT_io.out_col_type[0] = GMT_IS_LON;
  GMT_io.in_col_type[1] = GMT_io.out_col_type[1] = GMT_IS_LAT;
  GMT_set_xy_domain (wesn, &grd);
  GMT_RI_prepare (&grd);
  GMT_grd_RI_verify (&grd, 0);
  GMT_write_grd_info (grdfilename, &grd);
#endif
  if (GMT_write_grd (grdfilename, &grd, phival,
		     0.0, 0.0, 0.0, 0.0, 
		     GMT_pad, FALSE)) {
    fprintf (stderr, "%s: my_gmt_write_grd: error writing file %s\n", 
	     GMT_program, argv[1]);
    exit (EXIT_FAILURE);
  }
  if(verbose){
    fprintf(stderr,"%s: my_gmt_write_grd: written to %s\n",argv[0],grdfilename);
  }
}







