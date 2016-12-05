/******************************************************************************
 *                                                                            *
 *   random.c                                                                 *
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


#include "random.h"


/* STATIC VARIABLES ********************************************************/
/* an array needed by erand48 */

static unsigned short int xsubj[3];        /* array used to initialize erand48 */


/*** RANDOM NUMBER FUNCTIONS ***********************************************/

/*** InitERand: initializes ERand by making xsubj static to random.c *******
 ***************************************************************************/

void InitERand(unsigned short int xsubi[3])
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

// #ifdef USE_DRAND48
//     i = drand48();
// #elif  USE_ERAND48
	i = erand48(xsubj);
// #else
//     i = ((double)random()) / 2147483648.;
// #endif

  return i;
}


/*** RandomInt: returns a random integer between 0 and max using the *******
 *              random number generator of choice                          *
 ***************************************************************************/

int RandomInt(int max)
{
	register int   i;

// #ifdef USE_DRAND48
	// i = ((int)(drand48() * max)) % max;
// #elif  USE_ERAND48
	i = ((int)(erand48(xsubj) * max)) % max;
// #else
	// i = random() % max;
// #endif

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
