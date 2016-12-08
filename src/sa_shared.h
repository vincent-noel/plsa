/******************************************************************************
 *                                                                            *
 *   sa.h                                                                     *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   written by Jimmy Lam and Dan Greening                                    *
 *   modified by John Reinitz and Johannes Jaeger                             *
 *   modified by Vincent Noel                                                 *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   see ../doc/prog_man.ps for details of how this works                     *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   IMPORTANT: IF YOU EVER CHANGE ANYTHING IN THIS FILE, LET ALL             *
 *            YOUR FELLOW PROGRAMMERS KNOW WELL IN ADVANCE AND                *
 *            CONSULT WITH THEM IF THEY AGREE ON YOUR CHANGES!!               *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   The following came with lsa.c, obtained from Dan Greening,               *
 *   who got it from Jimmy Lam. It is probably copyright Jimmy Lam            *
 *   1988, with modifications by Greening. We (JR & JJ) have                  *
 *   mostly modified it by adding some comments and restructuring             *
 *   the code to make it more easily understandable. JR has added             *
 *   the criterion for continous search spaces and separated Lam              *
 *   statistics from move acceptance statistics. 							  *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   NOTE: this header only contains prototypes for functions used            *
 *       for serial annealing code; all prototypes that are spe-              *
 *       cific to parallel code need to go into MPI.h                         *
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

#ifndef SA_INCLUDED
#define SA_INCLUDED

/* this def needed for func. defs that refer to (* FILE) */
#ifndef _STDIO_INCLUDED
#include <stdio.h>
#endif

/* following for structures & consts used thruout */
#ifndef GLOBAL_INCLUDED
#include "global.h"
#endif


#include "types.h"

/*** CONSTANTS *************************************************************/

#define MIN_DELTA    -100.    /* minimum exponent for Metropolis criterion */
					/* provides a minimum probability for really bad moves */


/*** STRUCTS AND ENUMS *****************************************************/

/* These are the Lam parameters: *******************************************
 * (see King-Wai Chu's thesis for a detailed explanation)                  *
 *                                                                         *
 * lambda:              overall annealing accuracy                         *
 * lambda_mem_length_u: weighting factor for estimating mean energy        *
 * lambda_mem_length_v: weighting factor for estimating standard deviation *
 * initial_moves:       number of initial moves to randomize initial state *
 *                      of the system and gather initial statistics        *
 * tau:                 number of moves between updating statistics and    *
 *                      Lam parameters                                     *
 * freeze_count:        number of times system needs to be 'Frozen()' be-  *
 *                      fore we stop annealing                             *
 * update_S_skip:       number of iterations between changing the inverse  *
 *                      temperature; exists for increasing code efficiency *
 *                      and is probably obsolete now                       *
 * control:             OBSOLETE, used to be the control parameter for the *
 *                      move generation; acceptance statistics now live in *
 *                      move(s).c, since they are problem-specific         *
 * criterion:           defines the limits of what counts as 'Frozen()';   *
 *                      depends on the stop_flag (see below), i.e. it is   *
 *                      either a fixed energy (for discrete problems), a   *
 *                      limit to the change in absolute energy or the lim- *
 *                      it of proportional change in energy (good for con- *
 *                      tinous search spaces)                              *
 * mix_interval:        the mixing interval for parallel code in 'tau'     *
 *                      units, i.e. 100 means mix every 100 tau            *
 *                                                                         *
 ***************************************************************************
 *                                                                         *
 * NOTE: the stopping criterion is actually defined by the combination of  *
 *       'freeze_count' and 'criterion' below plus the stop_flag which     *
 *       defines the type of the criterion (absolute or proportional);     *
 *       usually when we talk about criterion (e.g. Chu, 2001) we mean     *
 *       only the value of criterion (kappa) though, since 'freeze_count'  *
 *       and the stop_flag should always be constant for specific problems *
 *                                                                         *
 ***************************************************************************/

/*** GLOBALS ***************************************************************/

// int         time_flag;                         /* flag for timing the code */
// int         log_flag;                 /* flag for displaying log to stdout */
//
// long        state_write;     /* frequency for writing state files (in tau) */
// long        print_freq;          /* frequency for printing to log (in tau) */
// long        captions;      /* option for printing freqency of log captions */

/* Special annealing modes: ************************************************
 *                                                                         *
 * the following are flags that should be set by cmd line opts for some    *
 * special annealing modes that we use for testing or gathering stats:     *
 *                                                                         *
 * - quenchit: means lower the temperature to zero *instantly* after fini- *
 *             shing initialize; this is not implemented for parallel code *
 *             since it doesn't take a long time to finish                 *
 *                                                                         *
 ***************************************************************************/

// int         quenchit;
// int         max_iter;
// int         max_seconds;
// int         start_time_seconds;


/*** FUNCTION PROTOTYPES ***************************************************/

/* problem independent functions that live in lsa.c ************************/

/* Initializing functions */

/*** Initialize: calls ParseCommandLine first; then does either initial ****
 *               randomization and collecting Lam stats or restores state  *
 *               of the annealer as saved in the state file                *
 ***************************************************************************/

// int InitializePLSA(char ** statefile);

/*** InitFilenames: initializes static file names needed in lsa.c **********
 ***************************************************************************/

void InitFilenames(void);

/*** InitialLoop: performs the two sets of initial moves: ******************
 *                   1. randomizing moves (not parallelized)               *
 *                   2. loop for initial collection of statistics          *
 ***************************************************************************/

double InitialLoop(SAType * state, double S_0);

/* Main loop and update functions */

/*** Loop: loops (making moves, updating stats etc.) until the system is ***
 *         considered frozen according to the stop criterion               *
 ***************************************************************************/

int Loop(SAType * state, char * statefile, StopStyle stop_flag);


void RestoreLog(SAType * state);



/* functions to write the .log */


void WriteLog(SAType * state);
/*** PrintLog: actually writes stuff to the .log file **********************
 ***************************************************************************/

void PrintLog(FILE *outptr, int init_moves);



#endif
