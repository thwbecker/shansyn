/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: determine_coeff.h,v 1.4 2001/03/17 15:29:33 becker Exp $ */

// which to use?

//#define CHEBEV_FUNC chebev
//#define CHEBEV_DER_INT_FUNC chebev_der_int
#define CHEBEV_FUNC norm_chebev
#define CHEBEV_DER_INT_FUNC norm_chebev_der_int

COMP_PRECISION chebnorm(int );
COMP_PRECISION chebev(COMP_PRECISION ,COMP_PRECISION ,
		      COMP_PRECISION *,int ,COMP_PRECISION );
COMP_PRECISION norm_chebev(COMP_PRECISION ,COMP_PRECISION ,
		      COMP_PRECISION *,int ,COMP_PRECISION );

COMP_PRECISION cheb_norm_damping(int, int ,int );
COMP_PRECISION chebev_der_int(COMP_PRECISION,COMP_PRECISION,int);
COMP_PRECISION norm_chebev_der_int(COMP_PRECISION,COMP_PRECISION,int);

void determine_coeff(COMP_PRECISION *,int ,
		     COMP_PRECISION *,COMP_PRECISION *, 
		     int ,COMP_PRECISION , COMP_PRECISION ,
		     COMP_PRECISION ,COMP_PRECISION ,
		     COMP_PRECISION *,  COMP_PRECISION *, 
		     COMP_PRECISION *,
		     COMP_PRECISION (*)(COMP_PRECISION,
					COMP_PRECISION,
					COMP_PRECISION *,
					int,
					COMP_PRECISION),
		     COMP_PRECISION (*)(int,int,int));

void svd_solver(COMP_PRECISION *, 
		COMP_PRECISION *,
		COMP_PRECISION *,
		int ,int );
void svd_driver(COMP_PRECISION *,int ,int ,COMP_PRECISION *,
		COMP_PRECISION *);
COMP_PRECISION spline_base(COMP_PRECISION ,COMP_PRECISION ,
		      COMP_PRECISION *,int ,COMP_PRECISION );

COMP_PRECISION spline_norm_damping(int ,int ,int);

#define splhsetup splhsetup_
#define splh splh_
#define svdcmp svdcmp_
#define svbksb svbksb_


void svdcmp(COMP_PRECISION *,int *,int *,int *,
	    int *,COMP_PRECISION *,COMP_PRECISION *,
	    COMP_PRECISION *);
void svbksb(COMP_PRECISION *,COMP_PRECISION *,
	    COMP_PRECISION *,int *,int *,
	    int *,int *,COMP_PRECISION *,
	    COMP_PRECISION *,COMP_PRECISION *);

void splhsetup(COMP_PRECISION *,COMP_PRECISION *,
	       COMP_PRECISION *,COMP_PRECISION *, 
	       int *);
COMP_PRECISION splh(int *,COMP_PRECISION *,COMP_PRECISION *,
		    COMP_PRECISION *,int *,COMP_PRECISION *);

