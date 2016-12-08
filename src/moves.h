/******************************************************************************
 *                                                                            *
 *   moves.h                                                                  *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   written by JR, modified by Yoginho and Vincent Noel                      *
 *   modified by Vincent Noel                                                 *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   moves.h contains problem-specific stuff for the Lam annealer.            *
 *   There are functions for move generation and translation of               *
 *   model data structures into annealing parameter arrays. More-             *
 *   over there are a couple of functions that read annealing and             *
 *   tune parameters and the tweak table from the data file.                  *
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


/* this def needed for func. defs that refer to (* FILE) */
#ifndef _STDIO_INCLUDED
#include <stdio.h>
#endif

/* following for structures & consts used thruout */
#ifndef GLOBAL_INCLUDED
#include "global.h"
#endif

/* following for StopStyle enum style */
#ifndef SA_INCLUDED
#include "sa_shared.h"
#endif

/* ... and this for the Range struct */
#ifndef SCORE_INCLUDED
#include "types.h"
#endif



/*** CONSTANTS *************************************************************/

#define THETA_MIN  0.            /* minimum value of theta_bar (move size) */
#define THETA_MAX  20.            /* maximum value of theta_bar (move size) */
#define THETA_INIT 5.0      /* initial value for all theta_bar (move size) */

#define LOWBITS    0x330E                /* following two are for drand to */
#define BYTESIZE   8                                   /* erand conversion */


/*** FUNCTION PROTOTYPES ***************************************************/


/* move generation functions that are used in lsa.c (live in move(s).c) */

/* GenerateMove: evaluates the old energy, changes a parameter, then eval- *
 *               uates the new energy; returns the difference between old  *
 *               and new energy to the caller                              *
 ***************************************************************************/

double GenerateMove(void);

/*** AcceptMove: sets new energy as the old energy for the next step and ***
 *               keeps track of the number of successful moves             *
 ***************************************************************************/

void AcceptMove(void);

/*** RejectMove: simply resets the tweaked parameter to the pretweak value *
 ***************************************************************************/

void RejectMove(void);

/*** GetEnergy: returned the last computed value of the scoring function   *
 *    To avoid computing e = old_e + (new e - old e), and using directly   *
 *    e = new_e
 *    Avoid precision errors due floating number representation            *
 **************************************************************************/

double GetNewEnergy(void);
double GetOldEnergy(void);

/* plsa.c: I/O functions for miscellaneous stuff */

void PrintTimes(FILE *fp, double *delta);

/* moves.c: functions for move generation */

/* Initializing and restoring functions */

/*** InitMoves: initializes the following moves.c-specific stuff: **********
 *              - static annealing parameter struct (ap)                   *
 *              - static tweak struct (tweak) in translate.c               *
 *              - initializes random number generator in lsa.c             *
 *              - receives parameter list from Translate, stores nparams   *
 *              - initializes acc_tab for acceptance statistics            *
 *              - set mixing interval in lsa.c (parallel code only)        *
 *                                                                         *
 *              it then returns the initial temperature to the caller      *
 ***************************************************************************/

void InitMoves(SAType * state, PArrPtr * pl);
//
// /*** RestoreMoves: restores move generator from state file *****************
//  *           NOTE: InitMoves will be called before this function during    *
//  *                 a restore                                               *
//  ***************************************************************************/
//
// void RestoreMoves(MoveState *MovePtr);

/* a function for finalizing a run */

/*** GetFinalInfo: collects stop energy and final count for output to the **
 *                 data file                                               *
 ***************************************************************************/

AParms GetFinalInfo(void);
//
// /* functions that communicate with other source files */
//
// /*** MoveSave: returns a MoveState struct in which the current state of ****
//  *             moves is saved; use for writing state file                  *
//  ***************************************************************************/
//
// MoveState *MoveSave(void);

/* state.c */

/* funcs that read/write/erase the state file */

/*** StateRead: reads Lam statistics, move state and erand state from a ****
 *              state file and restores the annealer's state to the same   *
 *              state it was in before it got interrupted                  *
 *     CAUTION: InitMoves should be called before StateRead!               *
 ***************************************************************************/

// void StateRead(char *statefile, Opts *options, MoveState *move_ptr,
// 		   double *stats, unsigned short *rand, double *times);

/*** StateRm: removes the state file after the run has been completed ******
 ***************************************************************************/

// void StateRm(void);

void randomModelParameters(void);
