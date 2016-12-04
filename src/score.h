/******************************************************************************
 *                                                                            *
 *   score.h                                                                  *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   written by JR, modified by Yoginho                                       *
 *   modified by Vincent Noel                                                 *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   This header file contains stuff that is needed for reading               *
 *   and defining facts and limits and for scoring functions. The             *
 *   functions declared here initialize or manipulate facts or                *
 *   data time tables, read and initialize limits and penalty (if             *
 *   needed and do the actual scoring.                                        *
 *                                                                            *
 *   Should be included in code that does scoring (e.g. printscore            *
 *   or annealing code).                                                      *
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

#ifndef SCORE_INCLUDED
#define SCORE_INCLUDED

/* following for structures & consts used thruout */
#ifndef GLOBAL_INCLUDED
#include "global.h"
#endif

#include "types.h"

/* FUNCTION PROTOTYPES *****************************************************/

/* Initialization Functions */

/*** InitScoring: sets the score function                             *******
 ***************************************************************************/

void        InitScoring     (SAType *) ;

/* Actual Scoring Functions */

/*** Score: as the name says, score runs the simulation, gets a solution   *
 *          and then compares it to the data using the Eval least squares  *
 *          function.                                                      *
 *   NOTE:  both InitZygote and InitScoring have to be called first!       *
 ***************************************************************************/

double      Score           (void);

// void        SaveBestScore   (double score);

#endif
