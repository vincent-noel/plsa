/******************************************************************************
 *                                                                            *
 *   global.h                                                                 *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   written by JR, modified by Yoginho                                       *
 *   modified by Vincent Noel                                                 *
 *                                                                            *
 *   THIS HEADER CONTAINS *ONLY* STUFF THAT IS NEEDED IN BOTH                 *
 *   LSA.C AND COST FUNCTION-SPECIFIC FILES.                                  *
 *                                                                            *
 *   PLEASE THINK TWICE BEFORE PUTTING ANYTHING IN HERE!!!!                   *
 *                                                                            *
 *   global.h gets included automatically by other header files               *
 *   if needed and is not included explicitly in .c files.                    *
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

#ifndef GLOBAL_INCLUDED
#define GLOBAL_INCLUDED

#ifndef _FLOAT_INCLUDED
#include <float.h>
#endif


/*** A global in global for debugging **************************************/

int    debug;                                            /* debugging flag */

/*** Constants *************************************************************/

/* IMPORTANT NOTE: the following used to be only 80 to be backward compa-  */
/*                 tible with punch cards that were of that maximum length.*/
/*                 We have decided to abandon punch-card compatibility,    */
/*                 since few people are actually using them anymore...     */
/*                 It was a tough decision though.      JR & JJ, July 2001 */

#define MAX_RECORD            256   /* max. length of lines read from file */

/* The following defines the maximum float precision that is supported by  */
/* the code.                                                               */

#define MAX_PRECISION          16

/* the following constant as a score tells the annealer to reject a move,  */
/* no matter what. It had better not be a number that could actually be a  */
/* score.                                                                  */
#define FORBIDDEN_MOVE    DBL_MAX      /* the biggest possible score, ever */

#endif
