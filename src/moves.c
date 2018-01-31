/******************************************************************************
 *                                                                            *
 *   moves.c                                                                  *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   written by JR, modified by Yoginho                                       *
 *   modified by Vincent Noel                                                 *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   This is the file that contains the functions that do the move            *
 *   generation for the Lam-schedule annealer. All this stuff is              *
 *   problem-dependent, but shouldn't have to worry about the ac-             *
 *   tual structure of the data in the model to be scored. That's             *
 *   why moves.c communicates with model-specific stuff via trans-            *
 *   late.c. Funcs in this file also shouldn't use any Lam-speci-             *
 *   fic stuff like SAType structs and such.                                  *
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

#ifdef ICC
#include <mathimf.h>
#else
#include <math.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "error.h"
#include "moves.h"
#include "random.h"
#include "score.h"
#include "distributions.h"               /* DistP.variables and prototypes */
#include "lsa.h"                /* *ONLY* for random number funcs and flags */
#include "config.h"
#include "state.h"

#ifdef MPI                      /* inludes parallel-specific stuff for MPI */
#include <mpi.h>                     /* this is the official MPI interface */
#include "MPI.h"  /* our own structs and such only needed by parallel code */
#endif



/*** STATIC VARIABLES ******************************************************/

static AParms    ap;                /* static copy of annealing parameters */

static ParamList *ptab;      /* array of pointers to parameters and ranges */
static AccStats  *acc_tab;   /* struct to accumulate acceptance statistics */

static int       nparams;            /* number of parameters to be tweaked */
static int       idx;              /* index in ptab of thing-to-be-tweaked */

static int       nhits;        /* number of moves since start of execution */
static int       nsweeps;     /* number of sweeps since start of execution */

static double    pretweak; /* used to restore param. value after rejection */

static double    old_energy; 		  /* two static vars used by Generate- */
static double    new_energy;                  /* Move to calculate delta_e */

#ifdef MPI
static long      *hits;                /* used for pooling number of moves */
static long      *success;         /* used for pooling number of successes */
static long      *tmp;             /* temp array for MPI_Allreduce sendbuf */
#endif

double GetNewEnergy(){return new_energy;}
double GetOldEnergy(){return old_energy;}



/*** FUNCTIONS *************************************************************/

/*** INITIALIZING AND RESTORING FUNCTIONS **********************************/

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


double InitMoves(SAType * state, PArrPtr * pl)
{
	int            i;                                  /* local loop counter */

	/* following is used to initialize erand48() */

	long           seedval;         /* contains random number generator seed */

	int            left;
	unsigned short left16;
	unsigned short middle16;
	unsigned short *xsubj;

	xsubj = (unsigned short *)calloc(3, sizeof(unsigned short));

	InitDistribution(state->dist_params);
	InitScoring(state);                   /* initializes facts and limits */
	
	old_energy = Score();
	if (old_energy == FORBIDDEN_MOVE)
		error("GenerateMove: 1st call gave forbidden move");
	/* read annealing paramters and parameters-to-be-tweaked */

	ap.max_count = 0;
	ap.seed = state->seed;
	ap.start_tempr = state->initial_temp;
	ap.gain = state->gain_for_jump_size_control;
	ap.interval = state->interval;
	ap.param_precision = state->param_precision;
	// ap.log_params = state->log_params;


	/* initialze the random number generator, now erand48() */
#ifdef MPI
	seedval  = ap.seed + myid;       /* each processor gets a different seed */
#else
	seedval  = ap.seed;
#endif

	xsubj[0] = LOWBITS;

	middle16 = (unsigned short)seedval;
	xsubj[1] = middle16;

	left     = seedval >> (BYTESIZE * sizeof(unsigned short));
	left16   = (unsigned short)left;
	xsubj[2] = left16;

	InitERand(xsubj);               /* makes the xsubj array static to lsa.c */

	/* Set up data structure for tweaking */

	nparams = pl->size;                       /* nparams is static to moves.c */
	ptab    = pl->array;            /* make parameter array static to moves.c */

	/* acc_tab is for statistics like acceptance ratio etc. */

	acc_tab = (AccStats *)calloc(nparams, sizeof(AccStats));

	for (i=0; i<nparams; i++)
	{
		acc_tab[i].acc_ratio = 0;
		acc_tab[i].theta_bar = THETA_INIT;
		acc_tab[i].hits      = 0;
		acc_tab[i].success   = 0;
	}

#ifdef MPI

	/* allocate static arrays for parallel code */

	hits    = (long *)calloc(nparams, sizeof(long));
	success = (long *)calloc(nparams, sizeof(long));
	tmp     = (long *)calloc(nparams, sizeof(long));

#endif

	/* Finally, return the initial energy. */
	return old_energy;
}



/*** RestoreMoves: restores move generator from state file *****************
 *           NOTE: InitMoves will be called before this function during    *
 *                 a restore                                               *
 ***************************************************************************/

void RestoreMoves(MoveState *MovePtr)
{
	int       i;                                       /* local loop counter */

	nparams    = MovePtr->nparams;                   /* restore static stuff */
	idx        = MovePtr->index;
	nhits      = MovePtr->nhits;
	nsweeps    = MovePtr->nsweeps;
	old_energy = MovePtr->old_energy;

	for(i=0; i < nparams; i++)                         /* restore parameters */
		*(ptab[i].param) = MovePtr->newval[i];

	free(MovePtr->newval);
	free(acc_tab);                     /* InitMoves has already been called. */

	acc_tab = MovePtr->acc_tab_ptr;             /* restore acceptance stats  */

	free(MovePtr);
}


/*** A FUNCTION FOR FINALIZING A RUN ***************************************/

/*** GetFinalInfo: collects stop energy and final count for output to the **
 *                 data file                                               *
 ***************************************************************************/

AParms * GetFinalInfo(void)
{
	free(acc_tab);

#ifdef MPI

	free(hits);
	free(success);
	free(tmp);

#endif

	ap.stop_energy = old_energy;
	ap.max_count   = nhits;

	return &ap;
}

/*** MOVE GENERATION *******************************************************/



/*** UpdateControl: each interval number of steps, acceptance stats are ****
 *                  dated here; this function also prints prolix stuff, if *
 *                  required.                                              *
 ***************************************************************************/

void UpdateControl(void)
{
	int        i;                                      /* local loop counter */
	double     x;                   /* temp variable to manipulate theta_bar */

#ifdef MPI

	for (i=0; i<nparams; i++)
	{
		hits[i]    = (long)acc_tab[i].hits;
		success[i] = (long)acc_tab[i].success;
	}

	for (i=0; i<nparams; i++)
	{
		tmp[i] = hits[i];
	}
	MPI_Allreduce(tmp, hits, nparams, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

	for(i=0; i<nparams; i++)
	{
		tmp[i] = success[i];
	}
	MPI_Allreduce(tmp, success, nparams, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

	for(i=0; i<nparams; i++)
	{
		acc_tab[i].hits    = (int)hits[i];
		acc_tab[i].success = (int)success[i];
	}

#endif

	for(i=0; i<nparams; i++)
	{

		acc_tab[i].acc_ratio =
			((double)acc_tab[i].success)/((double)acc_tab[i].hits);

		x  = log(acc_tab[i].theta_bar);
		x += ap.gain * (acc_tab[i].acc_ratio - 0.44);
		acc_tab[i].theta_bar = exp(x);


		if (acc_tab[i].theta_bar < THETA_MIN)
			acc_tab[i].theta_bar = THETA_MIN;
		if (acc_tab[i].theta_bar > THETA_MAX)
			acc_tab[i].theta_bar = THETA_MAX;

		acc_tab[i].hits    = 0;
		acc_tab[i].success = 0;

	}
}


/*** Move: tweaks the parameter-to-be-tweaked according to the current *****
 *         value of index; also calls UpdateControl if necessary           *
 ***************************************************************************/

int Move(void)
{
	double tweakee;
	double theta;
	// double sign;

	/* update counters */

	idx++;
	nhits++;

	idx     = idx % nparams;

#ifdef MPI
	nsweeps = (nhits / nparams) * nnodes;
#else
	nsweeps = (nhits / nparams);
#endif

	/* update statistics if interval passed & at least one sweep completed */

	if ( !(nsweeps % ap.interval) && !(idx) && (nsweeps) )
	{
		UpdateControl();                            /* see comments in moves.h */
	}

	tweakee  = *(ptab[idx].param);
	pretweak = tweakee;

	theta = generate_dev(acc_tab[idx].theta_bar);

	if (theta > THETA_MAX) theta = THETA_MAX;
	if (theta < -THETA_MAX) theta = -THETA_MAX;

	tweakee = tweakee * exp(theta);

	if (tweakee < ptab[idx].param_range.lower || tweakee > ptab[idx].param_range.upper)
		return -1;

	if (tweakee <= 0)
		return -1;

	if (isinf(tweakee) || isnan(tweakee))
		return -1;

    char buff[64];
    sprintf(buff, "%.*g", ptab[idx].significant_digits, tweakee);
    tweakee = atof(buff);

	*(ptab[idx].param) = tweakee;   /* original eqparms in score.c tweaked */

	return 0;

}





void WriteParamsTrace()
{
	FILE * trace_params;
	char trace_name[MAX_RECORD];
  #ifdef MPI
	sprintf(trace_name,"%s/trace/params/params_%d", getLogDir(), myid);
  #else
	sprintf(trace_name,"%s/trace/params/params_%d", getLogDir(), 0);
  #endif
	trace_params = fopen(trace_name, "a");

	fprintf(trace_params, "%g\t", new_energy);
	int i;
	for (i=0; i < nparams; i++)
		fprintf(trace_params, "%g\t", *(ptab[i].param));
	fprintf(trace_params, "\n");

	fclose(trace_params);
}
/* GenerateMove: wrapper for Move and Score; makes a move and reports ******
 *               difference of energies before and after the move          *
 ***************************************************************************/

double GenerateMove(void)
{

	// double     delta_e;           /* energy difference before and after move */

	/* for first call: check for valid parameters */
	// if (old_energy == -999.)
	// {
	// 	old_energy = Score();
	//
	// 	if (old_energy == FORBIDDEN_MOVE)
	// 		error("GenerateMove: 1st call gave forbidden move");
	// }

	/* make a move, score and return either FORBIDDEN_MOVE or delta_e */
	int res;
	res = Move();

	if (res == -1)
		return FORBIDDEN_MOVE;

	acc_tab[idx].hits++;

	new_energy = Score();

	if (logTraceParams() > 0)
		WriteParamsTrace();




	if (new_energy == FORBIDDEN_MOVE)
		return FORBIDDEN_MOVE;

	else
	{
		// printf("Generated Move : %.30f -> %.30f (%.30f)\n", old_energy, new_energy, (new_energy - old_energy));
		return (new_energy - old_energy);
	}
}



/*** AcceptMove: sets new energy as the old energy for the next step and ***
 *               keeps track of the number of successful moves             *
 ***************************************************************************/

void AcceptMove(void)
{
	old_energy = new_energy;
	acc_tab[idx].success++;
}



/*** RejectMove: simply resets the tweaked parameter to the pretweak value *
 ***************************************************************************/

void RejectMove(void)
{
	*(ptab[idx].param) = pretweak;
}

/*** MOVE GENERATION - PART 2: FUNCS NEEDED IN MOVES.C (BUT NOT LSA.C) *****/

/*** MoveSave: returns a MoveState struct in which the current state of ****
 *             moves is saved; use for writing state file                  *
 ***************************************************************************/

MoveState *MoveSave(void)
{
	MoveState    *move_stuff;

	move_stuff = (MoveState *)malloc(sizeof(MoveState));

	move_stuff->old_energy  = old_energy;     /* pretty straightforward, no? */
	move_stuff->pt          = ptab;
	move_stuff->acc_tab_ptr = acc_tab;
	move_stuff->newval      = NULL;
	move_stuff->index       = idx;
	move_stuff->nhits       = nhits;
	move_stuff->nparams     = nparams;
	move_stuff->nsweeps     = nsweeps;

	return move_stuff;
}





#ifdef MPI

/*** FUNCTIONS FOR SENDING MOVE STATES UPON MIXING (PARALLEL CODE ONLY) ****
 *   prototypes for these are in MPI.h; there's an extensive comment on    *
 *   how move state communication should be done at the beginning of lsa.c */

/*** MakeStateMsg: function to prepare a move state message which is then **
 *                 passed to other nodes via MPI; lsa.c doesn't know about *
 *                 the structs we use for acceptance statistics in         *
 *                 move(s).c, but we can safely assume that what we have   *
 *                 to send can be communicated as longs and doubles; thus, *
 *                 we split the move state message in two arrays, one for  *
 *                 the longs and one for the doubles; then we return the   *
 *                 arrays and their sizes to lsa.c                         *
 ***************************************************************************/

void MakeStateMsg(long **longbuf, int *lsize,
				  double **doublebuf, int *dsize)
{
	int    i, n_l, n_d;



	/* calculate buffer size, compare with move parameters below */

	*lsize = 2*nparams+3;
	*dsize = 3*nparams+1;

	/* allocate buffer */

	*longbuf   =   (long *)calloc(*lsize, sizeof(long));
	*doublebuf = (double *)calloc(*dsize, sizeof(double));

	/* pack longs into their buffer */

	(*doublebuf)[0]       = old_energy;

	(*longbuf)[0]         = (long)idx;
	(*longbuf)[1]         = (long)nhits;
	(*longbuf)[2]         = (long)nsweeps;

	for(i=0; i<nparams; i++)
	{

		n_l=i*2+3;
		n_d=i*3+1;

		(*longbuf)[n_l]     = (long)acc_tab[i].hits;
		(*longbuf)[n_l+1]   = (long)acc_tab[i].success;

		(*doublebuf)[n_d]   = *(ptab[i].param);
		(*doublebuf)[n_d+1] = acc_tab[i].acc_ratio;
		(*doublebuf)[n_d+2] = acc_tab[i].theta_bar;

	}
}



/*** AcceptMsg: gets the move state message from lsa.c and reinstalls acc- *
 *              eptance statistics into move.c; see the comment for Make-  *
 *              StateMsg above for the rationale behind the two arrays     *
 *              that are passed                                            *
 ***************************************************************************/

void AcceptStateMsg(long *longbuf, double *doublebuf)
{
	int i, n_l, n_d;

	idx        = (int)longbuf[0];
	nhits      = (int)longbuf[1];
	nsweeps    = (int)longbuf[2];

	old_energy = doublebuf[0];

	for (i=0; i<nparams; i++)
	{

		n_l=i*2+3;
		n_d=i*3+1;

		acc_tab[i].hits      = (int)longbuf[n_l];
		acc_tab[i].success   = (int)longbuf[n_l+1];

		*(ptab[i].param)     = doublebuf[n_d];
		acc_tab[i].acc_ratio = doublebuf[n_d+1];
		acc_tab[i].theta_bar = doublebuf[n_d+2];

	}

	free(longbuf);
	free(doublebuf);

}
#endif
