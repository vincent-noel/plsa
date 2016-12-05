/******************************************************************************
 *                                                                            *
 *   random.h                                                                 *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   written by Yoginho                                                       *
 *   modified by Vincent Noel                                                 *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   functions for initializing and running erand48() random                  *
 *   number generator                                                         *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   Copyright (C) 2016 Vincent Noel (vincent.noel@butantan.gov.br)           *
 *                                                                            *
 *   plsa is free software: you can redistribute it and/or modify             *
 *   it under the terms of the GNU General Public License as published by     *
 *   the Free Software Foundation, either version 3 of the License, or        *
 *   (at your option) any later version.                                      *
 *                                                                            *
 *   plsa is distributed in the hope that it will be useful,                  *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *   GNU General Public License for more details.                             *
 *                                                                            *
 *   You should have received a copy of the GNU General Public License        *
 *   along with plsa. If not, see <http://www.gnu.org/licenses/>.             *
 *                                                                            *
 ******************************************************************************/
 #ifndef _STDLIB_INCLUDED
 #include <stdlib.h>
 #endif
/*** A CONSTANT **************************************************************/

// #define USE_ERAND48             1
	   /* tells the code to use erand48 by default, use -D opt to override ***/
				/* alternative is: -DUSE_DRAND48 for using drand48 instead ***/

/*** FUNCTION PROTOTYPES *****************************************************/

/*** InitERand: initializes ERand by making xsubj static to lsa.c ************
 *****************************************************************************/

void InitERand(unsigned short *xsubi);

/*** RandomReal: returns a random real number between 0 and 1 using the  *****
 *              random number generator of choice                            *
 *****************************************************************************/

double RandomReal(void);

/*** RandomInt: returns a random integer between 0 and max using the *********
 *              random number generator of choice                            *
 *****************************************************************************/

int RandomInt(int max);

/*** GetERandState: returns the xsubj array, which is used to initialize *****
 *                  erand48(); used for saving the erand state in a state    *
 *                  file                                                     *
 *****************************************************************************/

unsigned short *GetERandState(void);
