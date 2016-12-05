/******************************************************************************
 *                                                                            *
 *   lsa.c                                                                    *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   Credits:                                                                 *
 *                                                                            *
 *   Originally written by Jimmy Lam and Dan Greening                         *
 *   Adaptation for continuous problems and original implemen-                *
 *   tation of the parallel algorithm by John Reinitz                         *
 *   Tuning, quenchit mode by King-Wai Chu                  *
 *   Partly rewritten, extended & documented by Johannes Jaeger               *
 *   Repackaged as a non-problem-specific library by Vincent Noel             *
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
#include <string.h>
#include <sys/times.h>
#include <time.h>                                    /* this is for time() */
#include <unistd.h>          /* for command line option stuff and access() */
#include <sys/time.h>

#include "sa_shared.h"
#include "moves.h"
#include "error.h"
#include "random.h"
#include "config.h"
#include "state.h"
#ifdef MPI
#include <mpi.h>
#include "MPI.h"
#include "tuning.h"
#endif

/* STATIC VARIABLES ********************************************************/

int write_llog;                        /* flag for writing local log files */

/* log and input/output file related variables *****************************/

// static char   *statefile;                        /* name of the state file */
static char   *logfile;                    /* name of the global .log file */
int         	start_time_seconds;

/* some energy-related variables *******************************************/

// static double energy;                                    /* current energy */

static double S;                                 /* current inverse energy */
static double dS;                      /* delta S: change in S during move */
static double S_0;                      /* the initial inverse temperature */

static double exp_arg; /* the exponent of the Metropolis criterion (-dE/T) */

/* Lam stats stuff: estimators, stats and acceptance ratios ****************/

static double mean;          /* mean energy, collected from tau last steps */
static double vari;      /* energy variance, collected from tau last steps */

static double estimate_mean;              /* Lam estimator for mean energy */
static double estimate_sd;  /* Lam estimator for energy standard deviation */

static int    success;                       /* number of successful moves */

static double alpha;         /* the third term of the Lam schedule formula */
static double acc_ratio;    /* average acceptance ratio for all parameters */

/* Lam stats stuff: variables for calculating the estimators  **************/

/* mean estimator */

static double w_a;                       /* w_a is the weight for the mean */

static double usyy;         /* these parameters store intermediate results */
static double usxy;               /* for the updating formulas for A and B */
static double usxx;                       /* see Lam & Delosme, 1988b, p10 */
static double usx;
static double usy;
static double usum;

static double A;            /* A and B are the parameters for the rational */
static double B;                /* function for the estimation of the mean */

/* sd estimator */

static double w_b;         /* w_b is the weight for the standard deviation */

static double vsyy;         /* these parameters store intermediate results */
static double vsxy;               /* for the updating formulas for D and E */
static double vsxx;                       /* see Lam & Delosme, 1988b, p10 */
static double vsx;
static double vsy;
static double vsum;

static double D;            /* D and E are the parameters for the rational */
static double E;  /* function for the estimation of the standard deviation */

/* Lam stats stuff: variables related to tau *******************************/

// static double Tau;     /* double version of tau to calculate mean and vari */
static int    proc_tau;  /* proc_tau = tau                     in serial   */
/* proc_tau = tau / (# of processors) in parallel */
static long   count_tau;  /* how many times we did tau (or proc_tau) moves */

/* the actual number of moves for collecting initial statistics ************/

static int    proc_init;                        /* number of initial moves */

/* stuff used by Frozen ****************************************************/

static double old_mean;                    /* old mean as stored by Frozen */
static int    counter;                           /* counter used by Frozen */

/* skip tells the annealer how often it should update S ********************
 * implemented by Lam for increasing performance, i.e. currently we only   *
 * update the temperature every 10 moves (possibly obsolete now)           */

static int    skip = -1;

/* variables used for timing */

/* these variables are used to evaluate real time using time() */
static double      start;                     /* wallclock time before run */
static double      finish;                     /* wallclock time after run */

/* these structs are used to evaluate user time using times() or MPI_Wtime */
static struct tms *cpu_start;                      /* user time before run */
static struct tms *cpu_finish;                      /* user time after run */


/*** INITIALIZING FUNCTIONS ************************************************/


/*** InitFilenames: initializes static file names that depend on the out- **
 *                  put file name (i.e. this *must* be called after we     *
 *                  have called RestoreState())                            *
 ***************************************************************************/

void InitFilenames(void)
{
	/* allocate memory for static file names */

	logfile   = (char *)calloc(MAX_RECORD, sizeof(char));

#ifdef MPI
	if (myid == 0)
	{
#endif
		/* the global .log file: stores iterations, temperature, change in tempe-  *
		 * rature, global means and standard deviation, Lam estimators for mean    *
		 * and sd as well as global acceptance ratios over an annealing run        */

		sprintf(logfile, "plsa.log");

#ifdef MPI
	}

	InitLocalFilenames();
#endif
}


/*** InitializeWeights: initialize weights a and b for calculating Lam *****
 *                      estimators; these weights are computed from the    *
 *                      lambda memory length products                      *
 ***************************************************************************/

void InitializeWeights(SAType * state)
{
	/* w_a is the weight for the mean */

	w_a = state->lambda_mem_length_u / state->lambda;
	w_a = 1.0 - state->tau / w_a;
	if (w_a < 0.0)
		w_a = 0.0;

	/* w_b is the weight for the standard deviation */

	w_b = state->lambda_mem_length_v / state->lambda;
	w_b = 1.0 - state->tau / w_b;
	if (w_b < 0.0)
		w_b = 0.0;

#ifdef MPI
	/* local weights: there are two sets of local statistics, since we require *
	 * a local Lam estimator of the mean energy with a small weight for the    *
	 * calculation of cross-correlation of processors (for estimating the      *
	 * lower bound of M), whereas we need a local Lam estimator for the mean   *
	 * energy with a large weight for the calculation of the variance of local *
	 * means; these two different weights for local Lam stats have been deter- *
	 * mined experimentally by King-Wai Chu (although it's not mentioned in    *
	 * his thesis); it works for all practical purposes, but in the future, a  *
	 * better way of sampling local statistics will be needed (probably as     *
	 * part of a general theory of parallel Lam simulated annealing)           */

	if ( state->options->tuning && nnodes>1 )
		InitializeLocalWeights(w_a, w_b);

#endif

}

/*** InitializeParameter: initializes variables for Lam annealing: this is *
 *                        executed after doing the initial steps;          *
 *                        the local parameters only need to be set for tu- *
 *                        ning code                                        *
 ***************************************************************************/

void InitializeParameter(SAType * state)
{
	double   d;             /* d is used to store intermediate results below */

#ifdef MPI
	// double   l_d;               /* l_d is used to store intermediate results */
#endif


	/* 1. set global parameters (serial and parallel code) *********************/

	/* set estimators to stats collected during initializing phase of run */

	estimate_sd   = sqrt(vari);
	estimate_mean = mean;

	/* initialize A,B,D,E according to Lam & Delosme, 1988b, p10 */

	A = estimate_sd * estimate_sd / (estimate_mean * estimate_mean);
	B = (1.0 / estimate_mean) - (A * S_0);
	D = estimate_sd / estimate_mean;
	E = (1.0 / estimate_sd) - (D * S_0);

	/* initialize these intermediate variables for updating funcs for A,B,D,E */

	usum = vsum = 1.0;
	usxy = usxx = usx = 0.0;
	usy  = 1.0 / estimate_mean;
	usyy = usy * usy;
	vsxy = vsxx = vsx = 0.0;
	vsy  = 1.0 / estimate_sd;
	vsyy = vsy * vsy;

	/* set the initial temperature and the initial delta S */

	S  = S_0;
	dS = 0.5 / estimate_sd;            /* keep--may not need--based on s_0=0 */

	/* alpha is the third term of the main Lam schedule formula */

	d     = (1.0 - acc_ratio) / (2.0 - acc_ratio);
	alpha = 4.0 * acc_ratio * d * d;

#ifdef MPI
	/* 2. set local parameters *************************************************
	 * these are only needed when tuning...                                    *
	 * note that the two sets of local estimators used for tuning only differ  *
	 * in their weights, so they can actually be initialized the same way      */

	if ( state->options->tuning && nnodes>1 )
		InitializeLocalParameters(S_0);

#endif

	/* weights determine how the estimators are sampled for times before tau */
	InitializeWeights(state);

}





/*** InitialLoop: performs the two sets of initial moves: ******************
 *                   1. randomizing moves (not parallelized)               *
 *                   2. loop for initial collection of statistics          *
 ***************************************************************************/

double InitialLoop(SAType * state, double s0)
{
	int         i;                                     /* local loop counter */

	double 		energy;
	double      energy_change;               /* change of energy during move */

#ifdef MPI
	double * 		total;
	InitializeMixing();

	proc_tau  = state->tau / nnodes;               /* local copy of tau */
	proc_init = state->initial_moves / nnodes;    /* # of initial moves */
#else
	proc_tau  = state->tau;                       /* static copy to tau */
	proc_init = state->initial_moves;             /* # of initial moves */
#endif

	S_0 = s0;

	struct timeval tp;
	gettimeofday(&tp, NULL);
	start_time_seconds = (int) tp.tv_sec;
	/* randomize initial state; throw out results; DO NOT PARALLELIZE! */

	for (i=0; i<state->initial_moves; i++)
	{
		/* make a move: will either return the energy change or FORBIDDEN_MOVE */

		energy_change = GenerateMove();

		/* Metropolis stuff here; we usually want FORBIDDEN_MOVE to be very large  *
		 * that's why we want to prevent overflows here (hence the 'if')           */

		if ( energy_change != FORBIDDEN_MOVE )
			exp_arg = -S_0 * energy_change;

		/* MIN_DELTA provides a min. probability with which any move is accepted */

		if ( exp_arg <= MIN_DELTA )
			exp_arg = MIN_DELTA;



		/* below, we apply the Metropolis criterion to accept or reject a move */

		if ( energy_change == FORBIDDEN_MOVE )
			RejectMove();

		else if ( (energy_change <= 0.0) || (exp(exp_arg) > RandomReal()) ){
			energy = GetNewEnergy();
			AcceptMove();
		}
		else
			RejectMove();


	}   /* end randomize initial state */

	/* set all stats to zero, collection starts below */

	mean      = 0.0;
	vari      = 0.0;
	success   = 0;

	/* loop to collect initial statistics; this one is parallelized */

	for ( i=0; i<proc_init; i++ )
	{

		/* make a move: will either return the energy change or FORBIDDEN_MOVE */

		energy_change = GenerateMove();

		/* Metropolis stuff here; we usually want FORBIDDEN_MOVE to be very large  *
		 * that's why we want to prevent overflows here (hence the 'if')           */
		if ( energy_change != FORBIDDEN_MOVE )
			exp_arg = -S_0 * energy_change;

		/* MIN_DELTA provides a min. probability with which any move is accepted */

		if( exp_arg <= MIN_DELTA )
			exp_arg = MIN_DELTA;

		/* below, we apply the Metropolis criterion to accept or reject a move */

		if ( energy_change == FORBIDDEN_MOVE )
			RejectMove();

		else if ( (energy_change <= 0.0) || (exp(exp_arg) > RandomReal()) )
		{
			energy = GetNewEnergy();
			AcceptMove();
			success++;
		}
		else RejectMove();

		/* collect stats */

		mean += energy;
		vari += energy * energy;

	}

#ifdef MPI

	/* collect local stats if tuning */
	if ( state->options->tuning && nnodes>1 )
		InitLocalStatsTuning(mean, vari, success);

	total = InitLocalStats(mean, vari, &success, i, state->initial_moves);

	/* mean and variance are now summed over all nodes */
	mean = total[0];
	vari = total[1];
	free(total);

	/* local stats are calculated here (only if tuning) */
	if ( state->options->tuning && nnodes>1 )
		Init2LocalStatsTuning(proc_init);

#endif

	/* global stats are calculated here */
	mean     /= (double)state->initial_moves;
	vari      = vari / ((double)state->initial_moves) - mean * mean;
	acc_ratio = ((double)success) / ((double)state->initial_moves);

	/* initialize Lam parameters used for calculating Lam estimators */
	InitializeParameter(state);

	return energy;
}



/*** MAIN LOOP AND UPDATE FUNCTIONS ****************************************/
/*** UpdateS: update inverse temperature S at every Sskip step *************
 ***************************************************************************/

void UpdateS(SAType * state)
{
	register double    d;              /* used to store intermediate results */
	S += dS;                   /* here, inverse temperature is updated by dS */

	/* we need to update Lam parameters here, since S has changed; A, B, C and */
	/* D get updated in UpdateParameters                                       */

	estimate_mean = 1.0 / (A*S + B);    /* for temperature updating formulas */
	estimate_sd   = 1.0 / (D*S + E);        /* see Lam & Delosme, 1988b, p10 */

#ifdef MPI

	/* do the same for local Lam parameters (for both lower and upper bounds) */

	if ( state->options->tuning && nnodes>1 )
		UpdateLocalSTuning(S);

#endif

	d   = S * estimate_sd;             /* intermediate for the specific heat */

	/* following lines implement the main Lam schedule formula */

	dS  = state->lambda * alpha / (d*d * estimate_sd);
	dS *= state->update_S_skip;       /* ... we have to muliply by skip */

#ifdef MPI
	dS *= ((double) state->tau)/((double)proc_tau);
#endif

	/* reset skip */

	skip = state->update_S_skip;
}



/*** UpdateStats: updates mean, variance and acc_ratio after tau moves *****
 ***************************************************************************/

void UpdateStats(SAType * state, int i)
{

#ifdef MPI
	double * total;

	total = UpdateLocalStats(mean, vari, i, state->tau, &success);

	/* local stats are updated below */
	mean = total[0];      /* mean and variance are now summed over all nodes */
	vari = total[1];
	free(total);

	if ( state->options->tuning && nnodes>1 )
		UpdateLocalStatsTuning(proc_tau);

#else

	mean /= ((double) state->tau);                                  /* collect some statistics */
	vari /= ((double) state->tau);

#endif

	acc_ratio = ((double)success) / ((double) state->tau) ;         /* update acceptance ratio */
}



/*** UpdateParameter: update parameters A, B, D and E and the estimators ***
 *                    for mean and standard deviation for the current S    *
 ***************************************************************************/

void UpdateParameter(void)
{
	register double   d;    /* d is used to store intermediate results below */

	/* this part of the code updates the estimator for the mean */

	d     = 1.0 / mean;

	/* first: multiply all intermediate vars by weights */

	usyy *= w_a;
	usxy *= w_a;
	usy  *= w_a;
	usx  *= w_a;
	usxx *= w_a;
	usum *= w_a;

	/* then: update all intermediate vars */

	usyy += d*d;
	usxy += S*d;
	usy  += d;
	usx  += S;
	usxx += S*S;
	usum += 1.0;

	/* ... and use intermediate vars to update A and B ... */

	A = (usum * usxy - usx * usy) / (usum * usxx - usx * usx);
	B = (usy - A * usx)/usum;

	/* ... which are then used to update the estimator for the mean */

	estimate_mean = 1.0 / (A * S + B);


	/* this part of the code updates the estimator for the standard deviation */

	if (vari > 0.0)
	{

		d     = 1.0 / sqrt(vari);

		/* first: multiply all intermediate vars by weights */

		vsyy *= w_b;
		vsxy *= w_b;
		vsy  *= w_b;
		vsx  *= w_b;
		vsxx *= w_b;
		vsum *= w_b;

		/* then: update all intermediate vars */

		vsyy += d*d;
		vsxy += S*d;
		vsy  += d;
		vsx  += S;
		vsxx += S*S;
		vsum += 1.0;

		/* ... and use intermediate vars to update D and E ... */

		D = (vsum * vsxy - vsx * vsy) / (vsum * vsxx - vsx * vsx);
		E = (vsy - D * vsx) / vsum;
	}

	/* ... which are then used to update the estimator for the std dev */

	estimate_sd = 1.0 / (D*S + E);

	/* alpha corresponds to the third term in the main Lam schedule formula,   *
	 * which is a measure of how efficiently the space state is samples; this  *
	 * term is at a maximum for acc_ratio = 0.44, see Lam & Delosme, 1988b, p1 */

	d = (1.0 - acc_ratio) / (2.0 - acc_ratio);
	alpha = 4.0 * acc_ratio * d * d;

}

/*** Frozen: returns TRUE if frozen, FALSE otherwise; 'frozen' is defined **
 *           by 'freeze_count', 'stop_flag' and 'criterion' (see also sa.h *
 *           for a more extensive comment on this)                         *
 ***************************************************************************/

int Frozen(SAType * state, StopStyle stop_flag)
{
	double  delta;

	if(stop_flag == proportional_freeze)
		delta = (mean - old_mean)/mean;
	else if (stop_flag == absolute_freeze)
		delta = mean - old_mean;
	else if (stop_flag == absolute_energy)
		delta = mean;



	if (delta <= 0.)
		delta = -delta;

	if (delta <= state->criterion)
		counter++;

	else
	{
		counter  = 0;
		old_mean = mean;
	}

	return(counter >= state->freeze_count);
}


void WriteScoreTrace(double t_energy, int acceptance)
{
	FILE * trace_energy;
	char best_name[MAX_RECORD];

#ifdef MPI
	sprintf(best_name,"%s/trace/score/score_%d", getLogDir(), myid);
#else
	sprintf(best_name,"%s/trace/score/score_%d", getLogDir(), 0);
#endif

	trace_energy = fopen(best_name,"a");
	fprintf(trace_energy,"%g\t%d\n", t_energy, acceptance);
	fclose(trace_energy);
}

/*** Loop: loops (making moves, updating stats etc.) until the system is ***
 *         considered frozen according to the stop criterion               *
 ***************************************************************************/

int Loop(SAType * state, double energy, char * statefile, StopStyle stop_flag)
{
	int    i;                                          /* local loop counter */
	// double energy;
	double energy_change;                                   /* local Delta E */
	double d;                /* difference between energy and estimated mean */

	/* quenchit mode: set temperature to (approximately) zero immediately */

	if ( state->options->quenchit )
		S = DBL_MAX;

	/* loop till the end of the universe (or till the stop criterion applies) */

	while (1)
	{
		/* reset statistics */
		mean    = 0.0;
		vari    = 0.0;
		success = 0;

#ifdef MPI
		ResetLocalStats();
#endif

		/* do proc_tau moves here */

		for (i=0; i<proc_tau; i++)
		{
			/* make a move: will either return the energy change or FORBIDDEN_MOVE */
			energy_change = GenerateMove();

			/* Metropolis stuff here; we usually want FORBIDDEN_MOVE to be very large  *
			 * that's why we want to prevent overflows here (hence the 'if'); we also  *
			 * need to avoid overflows with quenchit (where S is (almost) infinite!)   */

			if ( !state->options->quenchit && (energy_change != FORBIDDEN_MOVE ) )
				exp_arg = -S * energy_change;

			/* MIN_DELTA provides a min. probability with which any move is accepted */

			if ( (exp_arg <= MIN_DELTA) )
				exp_arg = MIN_DELTA;


			/* below, we apply the Metropolis criterion to accept or reject a move; in *
			 * quenchit mode, only lower energies are accepted                         */
			int acceptance_result = 0; // Initialized to 0 aka rejected

			if (energy_change == FORBIDDEN_MOVE) {
				RejectMove();
			} else if ( (energy_change <= 0.0) ||
					  ( (!state->options->quenchit) && (exp(exp_arg) > RandomReal()) ) ) {
				AcceptMove();
				energy = GetNewEnergy();

				success++;
				acceptance_result = 1;
#ifdef MPI
				if ( state->options->tuning && nnodes>1 )
					AddLocalSuccess();
#endif
			} else
				RejectMove();

			if (logTraceScore() > 0)
				WriteScoreTrace(GetNewEnergy(), acceptance_result);

			/* update statistics */

			mean    += energy;
			d        = energy - estimate_mean;
			vari    += d * d;
#ifdef MPI

			/* if tuning: calculate the local mean and variance */
			if ( state->options->tuning && nnodes>1 )
				CalculateLocalStats(energy);
#endif

			/* update temperature every 'skip' steps; this was put here by Jimmy Lam,  *
			 * probably to save computation time on his old Spark; I guess it's obso-  *
			 * lete by now, but since it doesn't seem to do any harm and saves us some *
			 * time, we left it in here                                                */

			if ( !state->options->quenchit )
				if ( --skip <= 0 )
					UpdateS(state);


		}                 /* this is the end of the proc_tau loop */

		/* have done tau moves here: update the 'tau' counter */

		count_tau++;

		/* calculate mean, variance and acc_ratio for the last tau steps; i is     *
		 * passed as an argument for checking if all local moves add up to Tau     */

		UpdateStats(state, i);

		/* check if the stop criterion applies: annealing and tuning runs (that    *
		 * aren't stopped by the tuning stop criterion) leave the loop here; equi- *
		 * libration runs exit below */

		if (state->options->max_iter > 0 && count_tau >= state->options->max_iter)
			return 0;//FinalMove(state);

		if (state->options->max_seconds > 0)
		{
			struct timeval tp;
			gettimeofday(&tp, NULL);
			int duration = ((int) tp.tv_sec) - start_time_seconds;

			if (duration > state->options->max_seconds)
				return 0;

		}

		if ( Frozen(state, state->options->stop_flag) )
			return 0;

		/* update Lam stats: estimators for mean, sd and alpha from acc_ratio (we  *
		 * don't need this in quenchit mode since the temperature is fixed to 0)   */

		if ( !state->options->quenchit )
			UpdateParameter();

 #ifdef MPI

		/* tuning code: first update local Lam estimators */

		if ( state->options->tuning && nnodes>1 )
		{
			UpdateLParameter(S);

			if (UpdateTuning(logfile))
				return 0;//FinalMove(state);
		}

		/* at each mix_interval: do some mixing */

		if ( count_tau % state->mix_interval == 0 )
			DoMix(energy, estimate_mean, S, state->options->tuning);

#endif

		/* write the log every print_freq * tau (not proc_tau!) */

#ifdef MPI
		if ( (count_tau % (state->options->print_freq*nnodes) == 0) )
#else
		if ( (count_tau % state->options->print_freq == 0) )
#endif

			WriteLog(state);

		/* the state file gets written here every state_write * tau */

#ifdef MPI
		if ( (count_tau % (state->options->state_write*nnodes) == 0) && !state->options->tuning )
#else
		if ( (count_tau % state->options->state_write == 0) )
#endif
			StateWrite(statefile, energy);

	}                                /* this is the end of the while(1) loop */
	return -1;
}                          /* this is the end of Loop */


/*** GetLamstats: returns Lam statistics in an array of doubles; used to ***
 *                store Lam statistics in a state file                     *
 ***************************************************************************/

double *GetLamstats(double energy)
{
	double *stats;

	stats = (double *)calloc(31, sizeof(double));

	stats[0] = (double)counter;

	stats[1]  = old_mean;
	stats[2]  = energy;

	stats[3]  = mean;
	stats[4]  = vari;

	stats[5]  = estimate_mean;
	stats[6]  = estimate_sd;

	stats[7]  = S;
	stats[8]  = dS;
	stats[9]  = S_0;

	stats[10] = alpha;
	stats[11] = acc_ratio;

	stats[12] = w_b;
	stats[13] = vsyy;
	stats[14] = vsxy;
	stats[15] = vsxx;
	stats[16] = vsx;
	stats[17] = vsy;
	stats[18] = vsum;
	stats[19] = D;
	stats[20] = E;

	stats[21] = w_a;
	stats[22] = usyy;
	stats[23] = usxy;
	stats[24] = usxx;
	stats[25] = usx;
	stats[26] = usy;
	stats[27] = usum;
	stats[28] = A;
	stats[29] = B;

	stats[30] = (double)count_tau;

	return(stats);
}



/*** GetTimes: returns a two-element array with the current wallclock and **
 *             user time to be saved in the state file; for parallel code  *
 *             we average the times for all processes                      *
 ***************************************************************************/

double *GetTimes(void)
{
	double *delta;
	double clk_tck = (double)CLOCKS_PER_SEC; /* system-specific clock tick duration */
#ifdef MPI
	double temp;
#endif

	delta = (double *)calloc(2, sizeof(double));
	/* measure user time */
	times(cpu_finish);
	/* then wallclock time */
#ifdef MPI
	finish = MPI_Wtime();

	temp     = finish - start;
	MPI_Allreduce(&temp, &delta[0], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	delta[0] /= nnodes;

	temp     = (double)(cpu_finish->tms_utime - cpu_start->tms_utime)/clk_tck;
	MPI_Allreduce(&temp, &delta[1], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	delta[1] /= nnodes;

#else
	finish = time(NULL);

	delta[0] = finish - start;
	delta[1] = (cpu_finish->tms_utime - cpu_start->tms_utime)/clk_tck;
#endif

	return delta;
}


/*** FUNCTIONS TO RESTORE THINGS IN LSA.C AFTER A RESTART ******************/

/*** RestoreLamstats: restores static Lam statistics in lsa.c from an ******
 *                    array of doubles; used to restore runs from a state  *
 *                    file.                                                *
 ***************************************************************************/

double RestoreLamstats(double *stats)
{
	double energy;
	counter = (int)rint(stats[0]);

	old_mean      = stats[1];
	energy        = stats[2];
	mean          = stats[3];
	vari          = stats[4];

	estimate_mean = stats[5];
	estimate_sd   = stats[6];

	S             = stats[7];
	dS            = stats[8];
	S_0           = stats[9];
	alpha         = stats[10];
	acc_ratio     = stats[11];

	w_b           = stats[12];
	vsyy          = stats[13];
	vsxy          = stats[14];
	vsxx          = stats[15];
	vsx           = stats[16];
	vsy           = stats[17];
	vsum          = stats[18];
	D             = stats[19];
	E             = stats[20];

	w_a           = stats[21];
	usyy          = stats[22];
	usxy          = stats[23];
	usxx          = stats[24];
	usx           = stats[25];
	usy           = stats[26];
	usum          = stats[27];
	A             = stats[28];
	B             = stats[29];

	count_tau = (long)rint(stats[30]);

	free(stats);
	return energy;
}



/*** RestoreTimes: restores the wallclock and user times if -t is used *****
 ***************************************************************************/

void RestoreTimes(double *delta)
{
	double clk_tck = (double)CLOCKS_PER_SEC; /* system-specific clock tick duration */

#ifdef MPI
	start = MPI_Wtime();
#else
	start = time(NULL);
#endif
	start -= delta[0];

	times(cpu_start);
	cpu_start->tms_utime -= (delta[1] * clk_tck);
}



/*** RestoreLog: restores the .log (and the .llog files) upon restart ******
 ***************************************************************************/

void RestoreLog(SAType * state)
{
	char   *shell_cmd;                             /* used by 'system' below */
	char   *outfile;                           /* temporary output file name */

	FILE   *logptr;                                     /* .log file pointer */
	FILE   *outptr;

	char   *logline;                                /* array of read buffers */
	long   saved_count_tau;          /* count_tau as read from the .log file */
	long   max_saved_count;    /* last count_tau that was saved in .log file */
	long   i;                                                /* loop counter */


	/* this is the last line we've written into the .log file */

	max_saved_count = state->initial_moves+proc_init+count_tau*proc_tau;

	logline   = (char *)calloc(MAX_RECORD, sizeof(char));
	shell_cmd = (char *)calloc(MAX_RECORD, sizeof(char));

	/* get temporary file name for output file */

#ifdef MPI
	if ( myid == 0 )
	{
#endif
		outfile   = (char *)calloc(MAX_RECORD, sizeof(char));
		outfile   = strcpy(outfile,"logXXXXXX");      /* required by mkstemp() */
		if ( mkstemp(outfile) == -1 )         /* get unique name for temp file */
			error("RestoreLog: error creating temporary (log) file");
#ifdef MPI
	}

	if ( myid == 0 )
	{
#endif

		/* restore the global .log file */
		saved_count_tau = -1;                             /* reset the counter */

		/* open the log file and a temporary output file */
		logptr = fopen(logfile, "r");
		if ( !logptr )
			file_error("RestoreLog (at open log file for reading)");

		outptr = fopen(outfile, "w");
		if ( !outptr )
			file_error("RestoreLog (at open temp log file for writing)");

		/* read and write the first few title and caption lines */
		for (i=0; i<4; i++)
		{
			if ( NULL == fgets(logline, MAX_RECORD, logptr ) )
				error("RestoreLog: error reading log file captions");
			fprintf(outptr, "%s", logline);
		}

		/* read and write the actual log lines till we are at current time */
		while ( (saved_count_tau < max_saved_count) &&
				(NULL != fgets(logline, MAX_RECORD, logptr)) )
		{
			if ( 1 != sscanf(logline, "%ld", &saved_count_tau) )
				error("RestoreLog: error reading saved_count_tau (after %d)",
					  saved_count_tau);
			fprintf(outptr, "%s", logline);
		}

		fclose(logptr);
		fclose(outptr);

		/* rename tmpfile into new file */

		sprintf(shell_cmd, "cp -f %s %s", outfile, logfile);

		if ( -1 == system(shell_cmd) )
			error("RestoreLog: error renaming temp file %s", outfile);

		if ( remove(outfile) )
			warning("RestoreLog: temp file %s could not be deleted",
					outfile);

#ifdef MPI
	}

	/* now do the same stuff for the local .llog files */

	if ( state->options->tuning && nnodes>1 )
		RestoreLocalLogTuning(max_saved_count);

#endif

}



/*** WriteLog: writes things like mean and variation, Lam estimators, dS, **
 *             alpha and acceptance ratio to the log files and to stdout   *
 *             (if -l is chosen).                                          *
 ***************************************************************************/

void WriteLog(SAType * state)
{
	FILE   *logptr;                      /* file pointer for global log file */

#ifdef MPI

	if (myid == 0)
	{
#endif
		logptr = fopen(logfile, "a");   /* first write to the global .log file */
		if ( !logptr )
			file_error("WriteLog");
		PrintLog(logptr, state->initial_moves);
		fclose( logptr );

#ifdef MPI
	}

	if ( write_llog && state->options->tuning && nnodes>1 )
		WriteLocalLog((state->initial_moves+proc_init+count_tau*proc_tau),
						S, dS);

	if ( myid == 0 )
	{
#endif

		if ( state->options->log_flag )                          /* display log to the screen? */
		{
			PrintLog(stdout, state->initial_moves);
			fflush( stdout );
		}

#ifdef MPI
	}
#endif

}


/*** PrintLog: actually prints the log to wherever it needs to be printed **
 ***************************************************************************/

void PrintLog(FILE *outptr, int initial_moves)
{

	int t_secs_since_epoc = time(NULL);

	const char *format =
			"%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n";

	/* print data */
	fprintf(outptr, format,
			(initial_moves+proc_init+count_tau*proc_tau),
			1.0/S, dS/S,
			mean, sqrt(vari), estimate_mean, estimate_sd,
			acc_ratio, alpha,t_secs_since_epoc);
}
