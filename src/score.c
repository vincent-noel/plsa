/******************************************************************************
 *                                                                            *
 *   score.c                                                                  *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   written by JR, modified by Yoginho                                       *
 *   modified by Vincent Noel                                                 *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   This file contains fuctions that are needed for reading and              *
 *   defining facts/limits and for scoring. The functions defined             *
 *   here initialize or manipulate facts or data time tables, read            *
 *   and initialize limits and penalty (if needed) and do the                 *
 *   actual scoring by least squares.                                         *
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

#include <math.h>
#include <stdio.h>

#include "config.h"    /* for bias and bcd reading functions & IGNORE */
#include "score.h"                                         /* obviously */

#ifdef MPI
#include "MPI.h"
#endif

double        (*scoreFunction)    ();
void          (*printFunction)    (char * path, int proc);

double        best_score;

void SaveBestScore(double score)
{


	char score_name[MAX_RECORD];
#ifdef MPI
	sprintf(score_name,"%s/best_score/score_%d", getLogDir(), myid);
#else
	sprintf(score_name,"%s/best_score/score_0", getLogDir());
#endif


	FILE * f_score = fopen(score_name,"w");
	fprintf(f_score, "%.5g\n", score);
	fclose(f_score);

}



/*** INITIALIZATION FUNCTIONS **********************************************/

/*** InitScoring: intializes a) facts-related structs and TTimes and *******
 *                           b) parameter ranges for the Score function.   *
 ***************************************************************************/

void InitScoring(SAType * tune)
{

	scoreFunction = tune->scoreFunction;
	printFunction = tune->printFunction;

	best_score = FORBIDDEN_MOVE;
}


/*** REAL SCORING CODE HERE ************************************************/


/*** Score: as the name says, score runs the simulation, gets a solution ***
 *          and then compares it to the data using the Eval least squares  *
 *          function                                                       *
 *   NOTE:  both InitZygote and InitScoring have to be called first!       *
 ***************************************************************************/
double Score(void)
{

	double     chisq   = 0;                    // summed squared differences

	chisq = scoreFunction();

	if (isnan(chisq) || isinf(chisq) || (chisq < 0))
		chisq = FORBIDDEN_MOVE;

	if (chisq < best_score)
	{
		if (logBestScore() > 0)
			SaveBestScore(chisq);

		if (logBestRes() > 0)
		{
			char log_dir[MAX_RECORD];
			sprintf(log_dir,"%s/best_res", getLogDir());

			if (printFunction != NULL)
#ifdef MPI
				printFunction(log_dir, myid);
#else
				printFunction(log_dir, 0);
#endif
		}
	}

	return chisq;
}
