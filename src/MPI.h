/******************************************************************************
 *                                  				    					  									*
 *   MPI.h                                                      			  			*
 *                                                              			  			*
 ******************************************************************************
 *                                                              			  			*
 *   written by John Reinitz                                    			  			*
 *   modified by King-Wai Chu, Johannes Jaeger                  			 				*
 *   modified by Vincent Noel                                                 *
 *                                                              			  			*
 ******************************************************************************
 *                                                              			  			*
 *   IMPORTANT: IF YOU EVER CHANGE ANYTHING IN THIS FILE, LET ALL 			 			*
 *            YOUR FELLOW PROGRAMMERS KNOW WELL IN ADVANCE AND  			  			*
 *            CONSULT WITH THEM IF THEY AGREE ON YOUR CHANGES!! 			 				*
 *                                                              			  			*
 ******************************************************************************
 *                                                              			  			*
 *   MPI.h contains structs and constants that are specific to    			  		*
 *   parallel annealing code using MPI.                           			  		*
 *                                                              	          	*
 *   This includes prototypes for all the tuning functions below. 	 		  		*
 *                                                              			  			*
 *   There are two problem-specific functions declared below that 			  		*
 *   need to be defined in move(s).c.                             			  		*
 *                                                              			  			*
 ******************************************************************************                                                               			  *
 *   NOTE: this header only contains prototypes for functions used 			  		*
 *       for parallel annealing code only; all prototypes of    			  			*
 *       functions that include serial code need to go into sa.h	  		  		*
 *                                 											  										*
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
#ifndef MPI_INCLUDED
#define MPI_INCLUDED
#include "tuning.h"


/*** CONSTANTS ****************************************************************/

// #define MAX_MIX         10000      /* max number of mixes during a tuning run */
// 					   /* set this to a lower number if you run out of memory */
// #define GROUP_SIZE         10       /* group size for calculating upper bound */
//
// #define STOP_TUNE_CNT      20                              /* stop tune count */
// #define STOP_TUNE_CRIT   0.05                        /* tuning stop criterion */
//
// #define LSTAT_LENGTH        1       /* length of Lam msg array when annealing */
// #define LSTAT_LENGTH_TUNE  28          /* length of Lam msg array when tuning */



/*** PARALLEL GLOBALS ******************************************************/

int myid;                                  /* id of local node (processor) */
int nnodes;               /* number of nodes (processors), also known as P */
//
// int tuning;                           /* flag for switching on tuning mode */
// // int covar_index;    /* covariance sample index for tuning (in 'tau' units) */
// // int write_tune_stat;               /* how often to write tuning statistics */
// // int auto_stop_tune;       /* auto stop tune flag to stop tuning runs early */
// int write_llog;                        /* flag for writing local log files */



/*** FUNCTION PROTOTYPES ***************************************************/

/* move(s).c: functions for communicating move state for mixing */

/*** MakeStateMsg: function to prepare a message which is then passed ******
 *                 to other nodes via MPI telling the other nodes about    *
 *                 move stats; since we don't know about the move state    *
 *                 structs, but can assume that we'll only have to send    *
 *                 longs and doubles, we split the message in two arrays,  *
 *                 one for the longs and one for the doubles; then we re-  *
 *                 turn the arrays and their sizes for lsa.c to send them  *
 *                 to the dance partners                                   *
 *                                                                         *
 *                 note that all the arguments need to get passed by refe- *
 *                 rence since we need to allocate the arrays according to *
 *                 their problem-specific size in move(s).c                *
 ***************************************************************************/

void MakeStateMsg(long   **longbuf,   int *lsize,
		  double **doublebuf, int *dsize);

/*** AcceptMsg: communicates a message about move stats received via MPI ***
 *              to move(s).c; see the comment for MakeStateMsg for the ra- *
 *              tionale behind the two arrays that are passed              *
 ***************************************************************************/

void AcceptStateMsg(long *longbuf, double *doublebuf);

#endif
