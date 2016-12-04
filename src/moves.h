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
#include "score.h"
#endif



/*** CONSTANTS *************************************************************/

#define THETA_MIN  0.            /* minimum value of theta_bar (move size) */
#define THETA_MAX  20.            /* maximum value of theta_bar (move size) */
#define THETA_INIT 5.0      /* initial value for all theta_bar (move size) */

#define LOWBITS    0x330E                /* following two are for drand to */
#define BYTESIZE   8                                   /* erand conversion */





/*** TONS OF STRUCTS *******************************************************/

/* The following are annealing parameters that are not specific to the Lam *
 * algorithm. In general they should be used in moves.c or plsa.c but    *
 * *not* in lsa.c. In the data file, the members of the struct labeled RO  *
 * are read from the $annealing_input section. They are used for initial   *
 * conditions of the annealer and do not change during a run. Members      *
 * labeled OUT are written to the $annealing_output section upon comple-   *
 * tion of a run.                                                          */

typedef struct {
  long   seed;                      /* seed for random number generator RO */
  double start_tempr;          /* the initial equilibration temperature RO */
  double gain;            /* gain for proportional control of move size RO */
  double stop_energy;                /* the final energy of the answer OUT */
  int    max_count;                      /* total number of iterations OUT */
  int    interval;       /* number of sweeps between updating theta_bar RO */
/*int    distribution;    1 - uniform; 2 - exp; 3 - normal; 4 - lorentz RO */
  // int    log_params;
} AParms;

/* Struct for acceptance statistics, which are used for keeping the accep- *
 * tance ratio for each parameter as close as possible to .44; there's one *
 * such struct for each parameter to be tweaked                            */

typedef struct {
  double acc_ratio;                      /* acceptance ratio for parameter */
  double theta_bar;              /* theta bar is proportional to move size */
  int    hits;       /* number of moves since last call to UpdateControl() */
  int    success;              /* number of these moves that were accepted */
} AccStats;





/* following struct contains copies of the static variables of moves.c     *
 * together with the values of parameters undergoing annealing             */

typedef struct {
  ParamList *pt;  /* Used during a save to point to annealed-on parameters */
  AccStats  *acc_tab_ptr;            /* points to current acceptance stats */
  double    *newval; /* points to array of annealed-on doubles for restore */
  double    old_energy;                     /* energy before the last move */
  int       nparams;                      /* # of parameters to be tweaked */
  int       index;      /* index of parameter to be tweaked during a sweep */
  int       nhits;                         /* number of moves already made */
  int       nsweeps;                         /* number of completed sweeps */
} MoveState;



/* Opts struct is for saving command line options in state.c */

typedef struct {
  StopStyle stop_flag;                                   /* stop criterion */
  int       time_flag;                             /* flag for timing code */
  int       log_flag;                                  /* log display flag */
  long      state_write;              /* frequency for writing state files */
  long      print_freq;         /* frequency for printing status to stdout */
  long      captions;                         /* opt for printing captions */
  int       quenchit;          /* flag for quenchit mode (T=0 immediately) */
  int       equil;              /* flag for equilibration mode (T = const) */
#ifdef MPI
  int       tuning;                                /* flag for tuning mode */
  int       covar_index; /* index for sample interval (=covar_index * tau) */
  int       write_tune_stat;    /* how many times to write tune statistics */
  int       auto_stop_tune;                      /* auto-stop tuning runs? */
#endif
} Opts;




/*** FUNCTION PROTOTYPES ***************************************************/

/* plsa.c: I/O functions for miscellaneous stuff */

void PrintTimes(FILE *fp, double *delta);


/* functions that communicate with state.c */

/*** GetOptions: returns command line options to state.c ***************
 ***************************************************************************/

Opts *GetOptions(void);

/*** RestoreOptions: restores the values of the command line opt variables *
 ***************************************************************************/

void RestoreOptions(Opts *options);



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

/*** RestoreMoves: restores move generator from state file *****************
 *           NOTE: InitMoves will be called before this function during    *
 *                 a restore                                               *
 ***************************************************************************/

void RestoreMoves(MoveState *MovePtr);

/*** RestoreProlix: restores prolix file after a run has been interrupted **
 ***************************************************************************/

void RestoreProlix(void);



/* a function for finalizing a run */

/*** GetFinalInfo: collects stop energy and final count for output to the **
 *                 data file                                               *
 ***************************************************************************/

AParms GetFinalInfo(void);



/* move generation functions that are used in moves.c, but not lsa.c */

/*** Move: tweaks the parameter-to-be-tweaked according to the current *****
 *         value of index; also calls UpdateControl if necessary           *
 ***************************************************************************/

int Move(void);

/*** UpdateControl: each interval number of steps, acceptance stats are ****
 *                  dated here; this function also prints prolix stuff, if *
 *                  required                                               *
 ***************************************************************************/

void UpdateControl(void);



/* functions that communicate with other source files */

/*** MoveSave: returns a MoveState struct in which the current state of ****
 *             moves is saved; use for writing state file                  *
 ***************************************************************************/

MoveState *MoveSave(void);





/* state.c */

/* funcs that read/write/erase the state file */

/*** StateRead: reads Lam statistics, move state and erand state from a ****
 *              state file and restores the annealer's state to the same   *
 *              state it was in before it got interrupted                  *
 *     CAUTION: InitMoves should be called before StateRead!               *
 ***************************************************************************/

void StateRead(char *statefile, Opts *options, MoveState *move_ptr,
		   double *stats, unsigned short *rand, double *times);

/*** StateRm: removes the state file after the run has been completed ******
 ***************************************************************************/

void StateRm(void);

void WriteParamsTrace(void);
void randomModelParameters(void);
