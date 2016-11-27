/*****************************************************************************
 *                                                                           *
 *   random.c                                                                *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   written by Yoginho                                                      *
 *   modified by Vincent Noel                                                *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   functions for initializing and running erand48() random                 *
 *   number generator                                                        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (C) 2016 Vincent Noel                                         *
 *   the full GPL copyright notice can be found in lsa.c                     *
 *                                                                           *
 *****************************************************************************/

// #include <stdlib.h>

#include "random.h"

/* following to toggle between random number generators ********************
 * USE_... constants are defined in random.h or need to be specified using *
 * -D with the compiler; using any random number generator other than      *
 * erand48() is not supported and makes results hard to reproduce; it      *
 * should work fine though if you DO NOT RESTART RUNS FROM .state FILES    *
 * WHEN USING RANDOM NUMBER GENERATORS OTHER THAN ERAND48()                *
 ***************************************************************************/

/* #ifdef USE_DRAND48
  double drand48();
  long   lrand48();
  void   srand48();
#elif USE_ERAND48
  double erand48();
  #endif  */


/* STATIC VARIABLES ********************************************************/
/* an array needed by erand48 */

#ifdef USE_ERAND48
static unsigned short xsubj[3];        /* array used to initialize erand48 */
#endif


/*** RANDOM NUMBER FUNCTIONS ***********************************************/

/*** InitERand: initializes ERand by making xsubj static to random.c *******
 ***************************************************************************/

void InitERand(unsigned short xsubi[3])
{
    int i;
    for(i=0; i<3; i++)
	      xsubj[i] = xsubi[i];

    free(xsubi);  /* temporarily crashed for lorraine 1-19-03 */
}


/*** RandomReal: returns a random real number between 0 and 1 using the  **
 *               random number generator of choice                        *
 **************************************************************************/

double RandomReal(void)
{
    double i;

#ifdef USE_DRAND48
    i = drand48();
#elif  USE_ERAND48
    i = erand48(xsubj);
#else
    i = ((double)random()) / 2147483648.;
#endif

  return i;
}


/*** RandomInt: returns a random integer between 0 and max using the *******
 *              random number generator of choice                          *
 ***************************************************************************/

int RandomInt(int max)
{
    register int   i;

#ifdef USE_DRAND48
    i = ((int)(drand48() * max)) % max;
#elif  USE_ERAND48
    i = ((int)(erand48(xsubj) * max)) % max;
#else
    i = random() % max;
#endif

    return i;
}


/*** GetERandState: returns the xsubj array, which is used to initialize ***
 *                  erand48(); used for saving the erand state in a state  *
 *                  file                                                   *
 ***************************************************************************/

unsigned short * GetERandState(void)
{
    unsigned short * p;
    p = xsubj;
    return(p);
}
