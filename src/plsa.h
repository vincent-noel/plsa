/******************************************************************************
 *                                                                            *
 *   plsa.h                                                                   *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   written by JR, modified by Yoginho                                       *
 *   modified by Vincent Noel                                                 *
 *   -D landscape option by Lorraine Greenwald Oct 2002                       *
 *   -g option by Yousong Wang, Feb 2002                                      *
 *   -a option by Marcel Wolf, Apr 2002                                       *
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
#include "types.h"


/* functions that communicate with state.c */

/*** GetOptions: returns command line options to state.c ***************
 ***************************************************************************/

Opts *GetOptions(void);

/*** RestoreOptions: restores the values of the command line opt variables *
 ***************************************************************************/

void RestoreOptions(Opts *options);





/*** InitialMove: initializes the following stuff: *************************
 *                - reads in Lam and other annealing parameters (passed to *
 *                  lsa.c through the state_ptr; the first three arguments *
 *                  are used to open the right data file etc.)             *
 *                - initializes the cost function, establishes link be-    *
 *                  tween cost function and annealer and passes the init-  *
 *                  tial energy to lsa.c by p_chisq)                       *
 *                - initializes move generation in move(s).c               *
 *                - sets initial energy by evaluating cost function for    *
 *                  the first time                                         *
 ***************************************************************************/

// void InitialMove(double *p_chisq);


/*** RestoreState: called when an interrupted run is restored; does the ****
 *                 following (see InitialMove for arguments, also see co-  *
 *                 mmunication functions below for how to restore the ran- *
 *                 dom number generator and Lam stats):                    *
 *                 - restores Lam and other annealing parameters           *
 *                 - reinitializes the cost function                       *
 *                 - restores move state in move(s).c                      *
 ***************************************************************************/

// void RestoreState(char *statefile, double *p_chisq);


/*** FinalMove: determines the final energy and move count and then prints *
 *              those to wherever they need to be printed to; also should  *
 *              do the cleaning up, i.e freeing stuff and such after a run *
 ***************************************************************************/

// double FinalMove();


/*** WriteTimes: writes the timing information to wherever it needs to be **
 *               written to at the end of a run                            *
 ***************************************************************************/

// void WriteTimes(double *times);
