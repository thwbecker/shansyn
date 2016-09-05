#include "shansyn.h"



COMP_PRECISION gammln(COMP_PRECISION xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

#include <stdio.h>
#define NR_END 1
#define FREE_ARG char*

void nrerror(error_text)
char error_text[];
/* Numerical Recipes standard error handler */
{


	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

COMP_PRECISION *vector(nl,nh)
long nh,nl;
/* allocate a COMP_PRECISION vector with subscript range v[nl..nh] */
{
	COMP_PRECISION *v;

	v=(COMP_PRECISION *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(COMP_PRECISION)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *cvector(nl,nh)
long nh,nl;
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;

	v=(unsigned char *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}

unsigned long *lvector(nl,nh)
long nh,nl;
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

COMP_PRECISION *dvector(nl,nh)
long nh,nl;
/* allocate a COMP_PRECISION vector with subscript range v[nl..nh] */
{
	COMP_PRECISION *v;

	v=(COMP_PRECISION *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(COMP_PRECISION)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

COMP_PRECISION **matrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a COMP_PRECISION matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	COMP_PRECISION **m;

	/* allocate pointers to rows */
	m=(COMP_PRECISION **) malloc((unsigned int)((nrow+NR_END)*sizeof(COMP_PRECISION*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(COMP_PRECISION *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(COMP_PRECISION)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

COMP_PRECISION **dmatrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a COMP_PRECISION matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	COMP_PRECISION **m;

	/* allocate pointers to rows */
	m=(COMP_PRECISION **) malloc((unsigned int)((nrow+NR_END)*sizeof(COMP_PRECISION*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(COMP_PRECISION *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(COMP_PRECISION)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((unsigned int)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

COMP_PRECISION **submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
COMP_PRECISION **a;
long newcl,newrl,oldch,oldcl,oldrh,oldrl;
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	COMP_PRECISION **m;

	/* allocate array of pointers to rows */
	m=(COMP_PRECISION **) malloc((unsigned int) ((nrow+NR_END)*sizeof(COMP_PRECISION*)));
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

COMP_PRECISION **convert_matrix(a,nrl,nrh,ncl,nch)
COMP_PRECISION *a;
long nch,ncl,nrh,nrl;
/* allocate a COMP_PRECISION matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	COMP_PRECISION **m;

	/* allocate pointers to rows */
	m=(COMP_PRECISION **) malloc((unsigned int) ((nrow+NR_END)*sizeof(COMP_PRECISION*)));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}


void free_vector(v,nl,nh)
COMP_PRECISION *v;
long nh,nl;
/* free a COMP_PRECISION vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(v,nl,nh)
int *v;
long nh,nl;
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(v,nl,nh)
long nh,nl;
unsigned char *v;
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(v,nl,nh)
long nh,nl;
unsigned long *v;
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(v,nl,nh)
COMP_PRECISION *v;
long nh,nl;
/* free a COMP_PRECISION vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(m,nrl,nrh,ncl,nch)
COMP_PRECISION **m;
long nch,ncl,nrh,nrl;
/* free a COMP_PRECISION matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(m,nrl,nrh,ncl,nch)
COMP_PRECISION **m;
long nch,ncl,nrh,nrl;
/* free a COMP_PRECISION matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
long nch,ncl,nrh,nrl;
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(b,nrl,nrh,ncl,nch)
COMP_PRECISION **b;
long nch,ncl,nrh,nrl;
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(b,nrl,nrh,ncl,nch)
COMP_PRECISION **b;
long nch,ncl,nrh,nrl;
/* free a matrix allocated by convert_matrix() */
{
  free((FREE_ARG) (b+nrl-NR_END));
}
