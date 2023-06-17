/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* 



*/
#include "shansyn.h"
#include <zlib.h>
int read_gz_line(gzFile *,char **);

int main(int argc, char *argv[] )
{
  int lmax,l,m,lmax_old;
  int nexp,lmsize,os,i,n,iexp;
  int cmode = 2;		/* 1L linear 2: pearson */
  float *ap=NULL,*bp=NULL,*at=NULL,*bt=NULL;
  gzFile in;
  char *line=NULL;
  COMP_PRECISION *x=NULL,*y=NULL,prob,corr,ppow,tpow;
  /*  */
  const int lmax_lim = 100, ide = 4, lmax_corr=50;
  /*  */
  if(argc<2){
    fprintf(stderr,"%s ptexpansion[.?.ab.gz]",argv[0]);
    exit(-1);
  }
  nexp = argc-1;
  for(iexp=0;iexp < nexp;iexp++){ /* read in all expansions */
    in = gzopen(argv[iexp+1],"r");
    if(!in){
      fprintf(stderr,"%s: cannot read %s\n",argv[0],argv[iexp+1]);
      exit(-1);
    }else{
      fprintf(stderr,"%s: reading %i PT expansion from %s\n",argv[0],iexp+1,argv[iexp+1]);
    }
    /* header */
    if(!read_gz_line(&in,&line)){
      fprintf(stderr,"%s: read error header lmax\n",argv[0]);
      exit(-1);
    }
    sscanf(line,"%i",&lmax);
    //fprintf(stderr,"%i\n",lmax);
    if(lmax>lmax_lim)
      lmax = lmax_lim;
    /*  */
    if(iexp > 0){
      if(lmax != lmax_old){
	fprintf(stderr,"%s: mismatch lmax %i %i \n",argv[0],lmax,lmax_old);
	exit(-1);
      }
    }else{
      lmax_old = lmax;
      /* 
	 
	 poloidal/toroidal expansion in Dahlen & Tromp format
	 
      */
      lmsize= (int)((((float)lmax)+1.0)*(((float)lmax)+2)/2.0);
      if(((ap=(float *)realloc(ap,lmsize*(argc-1)*sizeof(float))) ==NULL) ||
	 ((bp=(float *)realloc(bp,lmsize*(argc-1)*sizeof(float))) ==NULL) ||
	 ((at=(float *)realloc(at,lmsize*(argc-1)*sizeof(float))) ==NULL) ||
	 ((bt=(float *)realloc(bt,lmsize*(argc-1)*sizeof(float))) ==NULL)){
	fprintf(stderr,"%s: memerror, lmax=%i lmsize=%i iexp=%i \n",argv[0],lmax,lmsize,argc-1);
	exit(-1);
      }
    }
    /* input */
    for(l=0;l<=lmax;l++)
      for(m=0;m<=l;m++){ 
	os = POSLM(l, m) * nexp + iexp;
	if(!read_gz_line(&in,&line)){
	  fprintf(stderr,"%s: read error l %i m %i\n",argv[0],l,m);
	  exit(-1);
	}
	if(sscanf(line,"%f %f %f %f",(ap+os),(bp+os),(at+os),(bt+os))!=4){
	  fprintf(stderr,"%s: error with line %s\n",argv[0],line);
	  exit(-1);
	}
	//fprintf(stderr,"%22.15e %22.15e %22.15e %22.15e\n",ap[os],bp[os],at[os],bt[os]);
      }
    gzclose(in);
  }
  if(nexp < ide){
    fprintf(stderr,"%s: ide %i but only %i expansions\n",argv[0],ide,nexp);
    exit(-1);
  }
  for(iexp=ide;iexp < nexp-ide;iexp++){
    for(l=1;l<=lmax_corr;l++){
      n = 0;tpow=ppow=0.;
      for(m=0;m <= l;m++){
	for(i=iexp-ide;i<iexp+ide;i++){
	  os = POSLM(l, m) * nexp + i;
	  add_to_xy_float(  &x,&y,&n,ap[os],at[os]);
	  ppow += ap[os]*ap[os];
	  tpow += at[os]*at[os];
	  if(m != 0){
	    ppow += bp[os]*bp[os];
	    tpow += bt[os]*bt[os];
	    add_to_xy_float(&x,&y,&n,bp[os],bt[os]);
	  }
	}
      }
      ppow /= (COMP_PRECISION)((2*l)+1);
      tpow /= (COMP_PRECISION)((2*l)+1);
      corr = corr_sub(x,y,n,&prob,cmode);
      fprintf(stdout,"%i\t%i\t%12g\t%12e\t%12e\t%s\n",iexp+1,l,corr,ppow,tpow,argv[iexp+1]);
    }
  }
  free(x);free(y);free(ap);free(at);free(bp);free(bt);
}

/* 
   read one line of gzipped file 
*/
int read_gz_line(gzFile *in,char **line)
{
  int n;
  char tmp[1];
  
  n=0;
  while(gzread(*in,tmp,sizeof(char)) && (tmp[0]!='\n')){
    n++;
    if((*line = (char *)realloc(*line,sizeof(char)*(n+1)))==NULL){
      fprintf(stderr," read_gz_line: out of memory, n %i\n",n);
      exit(-1);
    }
    strncpy((*line+n-1), tmp, sizeof(char));
  }
  *line = (char *)realloc(*line,sizeof(char)*(n+1));
  *(*line+n)='\0';
  //fprintf(stderr,"%i: %s\n",n,*line);
  return n;
}
