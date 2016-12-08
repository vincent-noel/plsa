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
		   LamState *lam_state, unsigned short *rand, TuningSettings * t_settings,
		   double *delta);

void StateRm(void);



/* functions that are defined in lsa.c */

/*** RestoreMoves: restores move generator from state file *****************
 *           NOTE: InitMoves will be called before this function during    *
 *                 a restore                                               *
 ***************************************************************************/

void RestoreMoves(MoveState *MovePtr);

/* functions that communicate with other source files */

/*** MoveSave: returns a MoveState struct in which the current state of ****
 *             moves is saved; use for writing state file                  *
 ***************************************************************************/

MoveState *MoveSave(void);



/* functions that are defined in lsa.c */

/*** GetLamstats: returns Lam statistics in an array of doubles; used to ***
 *                store Lam statistics in a state file                     *
 ***************************************************************************/

LamState *GetLamstats(double energy);

/*** RestoreLamstats: restores static Lam statistics in lsa.c from an ******
 *                    array of doubles; used to restore runs from a state  *
 *                    file.                                                *
 ***************************************************************************/

double RestoreLamstats(LamState *stats);





/* functions that are defined in plsa.c */

/*** GetOptions: returns command line options to state.c ***************
 ***************************************************************************/

Opts *GetOptions(void);

/*** RestoreOptions: restores the values of the command line opt variables *
 ***************************************************************************/

void RestoreOptions(Opts *options);


/* Functions defined in tunning.c */
TuningSettings * GetTuningSettings();

void RestoreTuningSettings(TuningSettings * t_settings);



/*** GetTimes: returns a two-element array with the current wallclock and **
 *             user time to be saved in the state file                     *
 ***************************************************************************/
double *GetTimes(void);

/*** RestoreTimes: restores the wallclock and user times if -t is used *****
 ***************************************************************************/
void RestoreTimes(double *delta);
