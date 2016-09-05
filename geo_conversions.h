/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: geo_conversions.h,v 1.5 2001/06/16 21:37:45 becker Exp $ */
#define THETA2LATITUDE(x) ( (90.0 - (x)*ONEEIGHTYOVERPI) )
#define LATITUDE2THETA(x) ( (90.0 - (x))*PIOVERONEEIGHTY )
#define DEG2RAD(x) ( (x)*PIOVERONEEIGHTY )
#define PHI2LONGITUDE(x) ( (x) * ONEEIGHTYOVERPI )
#define LONGITUDE2PHI(x) ( (x) * PIOVERONEEIGHTY )
#define RAD2DEG(x) ( (x) * ONEEIGHTYOVERPI)

