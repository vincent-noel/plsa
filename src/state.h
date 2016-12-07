/******************************************************************************
 *                                                              					 	  *
 *   state.h                                                			 				*
 *                                                              						 	*
 ******************************************************************************
 *                                                               						  *
 *   written by JR, modified by Yoginho                         			 				*
 *   modified by Vincent Noel                                                	*
 *                                                              			 				*
 ******************************************************************************
 *                                                              			 				*
 *   state.c contains three functions that read, write and remove the    	*
 *	 state file for an annealing run. The frequency with which state are	 		*
 *   saved can be chosen by the command line option -b						 						*
 *   (for backup stepsize). The state file is very useful for the case		 		*
 *   when long annealing runs have to be interrupted or crash for some				*
 *   reason or another. The run can then be resumed by indicating the		 			*
 *	 state file as an additional argument to plsa.                      	 		*
 *                                                              			 				*
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

 #include "types.h"

 /* a function that writes the .state file (should live in state.c) */

 /*** StateWrite: collects Lam statistics, move state and the state of the **
  *               erand48 random number generator and writes all that into  *
  *               the state file, which can then be used to restore the run *
  *               in case it gets interrupted                               *
  ***************************************************************************/

void StateWrite(char * statefile, double energy);

void StateRead(char *statefile, Opts *options, MoveState *move_ptr,
		   LamState *lam_state, unsigned short *rand, double *delta);

void StateRm(void);
