/******************************************************************************
 *                                                                            *
 *   mixing.h                                                                 *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   Written by Vincent Noel                                                 *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   This file...											                  *
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
#include <stdio.h>
#include "types.h"

#define LSTAT_LENGTH        1       /* length of Lam msg array when annealing */

 /* tuning functions */

 /*** InitTuning: sets up/restores structs and variables for tuning runs ****
 ***************************************************************************/

void 		InitLocalFilenames();
void 		InitializeLocalParameters	(double S_0);
void 		InitializeLocalWeights		(double w_a, double w_b);

void 		InitializeMixing();

double * 	InitLocalStats				(double mean, double vari,
											int * success, int i,
											int initial_moves);

double *	UpdateLocalStats			(double mean, double vari, int i,
											int tau, int * success);

void 		ResetLocalStats				();
void 		AddLocalSuccess				();

/* lsa.c: parallel non-tuning funcs: update func for local Lam parameters */
/*** UpdateLParameter: update local parameters l_A, l_B, l_D and l_E and ***
*                    the local estimators for mean and standard deviation *
*                    for both upper and lower bounds of M_Opt for the     *
*                    current S                                            *
***************************************************************************/

void 		UpdateLParameter			(double S);


/* lsa.c: parallel non-tuning funcs: mixing functions */

/*** DoMix: does the mixing; sends move state and local Lam stats to the ***
*          dance partner(s)                                               *
***************************************************************************/

double 		DoMix						(double energy, double estimate_mean,
											double S, int tuning);


/*** MakeLamMsg: packages local Lam stats into send buffer *****************
***************************************************************************/

void 		MakeLamMsg					(double **sendbuf, double energy, int tuning);

/*** AcceptLamMsg: receives new energy and Lam stats upon mixing ***********
***************************************************************************/

double 		AcceptLamMsg				(double *recvbuf, int tuning);
void 		WriteLocalLog				(int count_steps, double S, double dS);
void 		PrintLocalLog				(FILE *outptr, int count_steps,
											double S, double dS);
