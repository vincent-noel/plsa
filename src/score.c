/*****************************************************************************
 *                                                                           *
 *   score.c                                                                 *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   written by JR, modified by Yoginho                                      *
 *   modified by Vincent Noel                                                *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   This file contains fuctions that are needed for reading and             *
 *   defining facts/limits and for scoring. The functions defined            *
 *   here initialize or manipulate facts or data time tables, read           *
 *   and initialize limits and penalty (if needed) and do the                *
 *   actual scoring by least squares.                                        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (C) 2016 Vincent Noel                                         *
 *   the full GPL copyright notice can be found in lsa.c                     *
 *                                                                           *
 *****************************************************************************/

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
    sprintf(score_name,"%s/score/score_%d", getLogDir(), myid);
#else
    sprintf(score_name,"%s/score/score_0", getLogDir());
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
        if (logScore() > 0)
			SaveBestScore(chisq);

#ifdef MPI
        printFunction(getLogDir(), myid);
#else
        printFunction(getLogDir(), 0);
#endif
    }

    return chisq;
}
