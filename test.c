#include "shansyn.h"

void main (int argc, char **argv)
{

  int i, j;
  //double w=0, e=360, s=-90, n=90, dx = 2, dy = 2;
  double w=0,e=360,s=-89,n=89,dx= 1.40625,dy= 1.4015748031496063;
  float *a;
  char *this=CNULL;
  struct GRD_HEADER grd;
  GMT_begin (argc, argv);
  GMT_grd_init (&grd, argc, argv, FALSE);
  grd.node_offset = 0;
  grd.x_min = w;	grd.x_max = e;
  grd.y_min = s;	grd.y_max = n;
  grd.x_inc = dx;	grd.y_inc = dy;
  GMT_grd_RI_verify (&grd, 0);
  grd.nx = irint ((e-w)/dx) + 1;
  grd.ny = irint ((n-s)/dy) + 1;
  a=(float *)malloc(sizeof(float)*grd.nx*grd.ny);
  for(i=0;i<grd.ny;i++)
    for(j=0;j<grd.nx;j++)
      a[i*grd.nx+j]=i*j;
  my_gmt_write_grd(a, TRUE, argc,argv, argv[1], grd.nx, grd.ny,
		   grd.x_min, grd.x_max,grd.y_min,grd.y_max,
		   grd.x_inc, grd.y_inc);
  /*
  if (GMT_write_grd (argv[1], &grd, a, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE)) {
    fprintf (stderr, "%s: Error writing file %s\n", GMT_program, argv[1]);
    exit (EXIT_FAILURE);
  }
  */

  GMT_free ((void *)a);
}
