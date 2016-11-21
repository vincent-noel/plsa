/*****************************************************************************
 *                                                                           *
 *   random.h                                                                *
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
#ifndef RANDOM_INCLUDED
#define RANDOM_INCLUDED

/*** A CONSTANT ************************************************************/

#define USE_ERAND48             1
       /* tells the code to use erand48 by default, use -D opt to override */
                /* alternative is: -DUSE_DRAND48 for using drand48 instead */

/*** FUNCTION PROTOTYPES ***************************************************/

/*** InitERand: initializes ERand by making xsubj static to lsa.c **********
 ***************************************************************************/

void InitERand(unsigned short *xsubi);

/*** RandomReal: returns a random real number between 0 and 1 using the  **
 *              random number generator of choice                          *
 **************************************************************************/

double RandomReal(void);

/*** RandomInt: returns a random integer between 0 and max using the *******
 *              random number generator of choice                          *
 ***************************************************************************/

int RandomInt(int max);

/*** GetERandState: returns the xsubj array, which is used to initialize **
 *                  erand48(); used for saving the erand state in a state *
 *                  file                                                  *
 **************************************************************************/

unsigned short *GetERandState(void);

#endif
