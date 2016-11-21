/*****************************************************************************
 *                                                                           *
 *   score.h                                                                 *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   written by JR, modified by Yoginho                                      *
 *   modified by Vincent Noel                                                *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   This header file contains stuff that is needed for reading              *
 *   and defining facts and limits and for scoring functions. The            *
 *   functions declared here initialize or manipulate facts or               *
 *   data time tables, read and initialize limits and penalty (if            *
 *   needed and do the actual scoring.                                       *
 *                                                                           *
 *   Should be included in code that does scoring (e.g. printscore           *
 *   or annealing code).                                                     *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (C) 2016 Vincent Noel                                         *
 *   the full GPL copyright notice can be found in lsa.c                     *
 *                                                                           *
 **************************************************************************** */

#ifndef SCORE_INCLUDED
#define SCORE_INCLUDED

/* following for structures & consts used thruout */
#ifndef GLOBAL_INCLUDED
#include "global.h"
#endif

#include "sa.h"

/* FUNCTION PROTOTYPES *****************************************************/

/* Initialization Functions */

/*** InitScoring: sets the score function                             *******
 ***************************************************************************/

void        InitScoring     (plsa_parameters * ) ;

/* Actual Scoring Functions */

/*** Score: as the name says, score runs the simulation, gets a solution   *
 *          and then compares it to the data using the Eval least squares  *
 *          function.                                                      *
 *   NOTE:  both InitZygote and InitScoring have to be called first!       *
 ***************************************************************************/

double      Score           (void);

void        SaveBestScore   (double score);

#endif
