/******************************************************************************
 *                                                                            *
 *   mixing.c                                                                 *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   Written by Vincent Noel                                                  *
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

 /*** TUNING CODE ***********************************************************/

 /*** DoTuning: calculates the cross-correlation (for lower bound) and the **
  *             variance of local means (for upper bound) for a sub_tune_   *
  *             interval; the results are added to arrays, which span the   *
  *             whole tuning_interval; when these values get written to the *
  *             'bound' files, they will be divided by the number of mixes  *
  *             we have done already to average them out (see also King-Wai *
  *             Chu's thesis, p. 63)                                        *
  ***************************************************************************/
#include <math.h>
#include <float.h>                                    /* for double limits */
#include <time.h>                                    /* this is for time() */
#include <string.h>

#include <mpi.h>


#include "MPI.h"
#include "tuning.h"
#include "mixing.h"

#include "error.h"
#include "random.h"
#include "state.h"


static int covar_index;    /* covariance sample index for tuning (in 'tau' units) */
static int write_tune_stat;               /* how often to write tuning statistics */
static int auto_stop_tune;       /* auto stop tune flag to stop tuning runs early */
/* parallel code: variables for local Lam stats for tuning */

/* local weights: there are two sets of local estimators for the mean ******
 * energy, since we require a local Lam estimator of the mean energy with  *
 * a small weight for the calculation of cross-correlation of processors   *
 * (for estimating the lower bound of M), whereas we need a local Lam esti-*
 * mator for the mean energy with a large weight for the calculation of    *
 * the variance of local means; these two different weights for local Lam  *
 * stats have been determined experimentally by King-Wai Chu (although     *
 * it's not mentioned in his thesis); it works for all practical purposes, *
 * but in the future, a better way of sampling local statistics will be    *
 * needed (probably as part of a general theory of parallel Lam simulated  *
 * annealing)                                                              *
 *                                                                         *
 * variables for estimators for the lower bound of M end with _l           *
 * variables for estimators for the upper bound of M end with _u           *
 *                                                                         *
 * note that the l_vari is calculated using the upper bound estimator for  *
 * the mean energy (l_estimate_mean_u); while tuning, the upper bound Lam  *
 * statistics are written to the .llog files                               *
 *                                                                         *
 * note that we don't need seperate estimators for the standard since the  *
 * only thing we do with the local sd estimators is writing them to local  *
 * .llog files (only the upper bound estimators get written there)         */
static int    l_success;               /* local number of successful moves */

static double l_mean;       /* local mean energy, from proc_tau last steps */
static double l_vari;   /* local energy variance, from proc_tau last steps */

static double l_estimate_mean_l;     /* estimator for mean for lower bound */
static double l_estimate_mean_u;     /* estimator for mean for upper bound */
static double l_estimate_sd;       /* Lam estimator for sd for upper bound */

static double l_acc_ratio;  /* average acceptance ratio for all parameters */


/* parallel code: variables for local sd estimator (see also comment on    *
 *                local weights above)                                     */

static double l_D;    /* l_D_l and l_E_l: parameters for the rational func */
static double l_E;     /* for the estimation of the sd for the lower bound */

/* parallel code: variables for calculating local mean estimator for lower *
 *                bound of M (see also comment on local weights above)     */

static double l_A_l;  /* l_A_l and l_B_l: parameters for the rational func */
static double l_B_l; /* for the estimation of the mean for the lower bound */

/* parallel code: variables for calculating local mean estimator for upper *
 *                bound of M (see also comment on local weights above)     */


static double l_A_u;        /* A and B are the parameters for the rational */
static double l_B_u;            /* function for the estimation of the mean */

/* vars used for tuning runs ***********************************************/

/* general variables that determine tune- and sub_tune_interval */

/* various counter variables and an array to save counts */

static int    count_sample = 0;        /* how many samples did we collect? */
static int    count_tune   = 0;      /* how many times did we do sub_tune? */
static int    count_mix    = 0;      /* counts the times we've been mixing */

static int    moves_tune = 0;   /* move counter: reset every tune_interval */

static int    *tau_count; /* # of tau's we've done for whole tune_interval */


/* a variable for the tuning stop criterion */

static int    stop_tune_count = STOP_TUNE_CNT;          /* stop tune count */

static char   *l_logfile;                  /* name of the local .llog file */

/* Files for tuning */
static char   *lbfile;               /* name of .lb file (for lower_bound) */
static char   *ubfile;               /* name of .ub file (for upper_bound) */
static char   *mbfile;                 /* name of .mb file (for mix_bound) */

/* vars used for tuning runs ***********************************************/

/* general variables that determine tune- and sub_tune_interval */

static int    covar_sample; /* # of moves over which we sample tuning stat */
static int    sample_size;                /* so many samples per processor */
static int    tune_interval;  /* total # of moves to be sampled for tuning */
static int    sub_tune_interval;     /* at end of this we write tune stats */


/* arrays for sampling stats for the lower bound */

static double *dev;         /* standard deviations for a sub_tune_interval */
static double *tot_dev;           /* gathers all local standard deviations */
static double *coll_dev;             /* collects local standard deviations */
static double *cross_correl;         /* cross-correlations for lower bound */

/* arrays for sampling stats for the upper bound */

static double *means;          /* used to save local means for upper bound */
static double *tot_means;                       /* gathers all local means */
static double *coll_means;                         /* collects local means */
static double *var_means;       /* variance of local means for lower bound */
static int    *midpoints;           /* midpoints of groups for upper bound */

/* an array needed for mixing in parallel code *****************************/

static int    *dance_partner;       /* stores dance partners for each node */

void InitLocalTuningFilenames()
{
	lbfile    = (char *)calloc(MAX_RECORD, sizeof(char));
	ubfile    = (char *)calloc(MAX_RECORD, sizeof(char));
	mbfile    = (char *)calloc(MAX_RECORD, sizeof(char));

	/* the following files are only needed for tuning */

	/* the lower_bound file: contains cross-correlations for a tune_interval   *
	 * averaged over all the past tune_intervals used to calculate the lower   *
	 * bound for the mixing_interval M                                         */

	sprintf(lbfile, "plsa.lb");

	/* the upper_bound file: for saving variance of local means for a tune_    *
	 * interval averaged over all the past tune_intervals; used to calculate   *
	 * the upper bound for the mixing_interval M when tuning                   */

	sprintf(ubfile, "plsa.ub");

	/* the mix_bound file: contains the history of both upper and lower bound  *
	 * for M over all past tune_intervals; used to check for the convergence   *
	 * of the estimate for the upper bound of M during a tuning run; also sto- *
	 * res the cross-corrlation and variance of local means for both bounds    */

	sprintf(mbfile, "plsa.mb");
}


TuningSettings * GetTuningSettings()
{
	TuningSettings * t_settings = (TuningSettings *) malloc(sizeof(TuningSettings));
	t_settings->covar_index = covar_index;
	t_settings->write_tune_stat = write_tune_stat;
	t_settings->auto_stop_tune = auto_stop_tune;
	return t_settings;
}

void RestoreTuningSettings(TuningSettings * t_settings)
{
	covar_index = t_settings->covar_index;
	write_tune_stat = t_settings->write_tune_stat;
	auto_stop_tune = t_settings->auto_stop_tune;
}


void RestoreTuningSettings(TuningSettings * tuning_settings);


// /*** DoTuning: calculates the cross-correlation (for lower bound) and the **
// *             variance of local means (for upper bound) for a sub_tune_   *
// *             interval; the results are added to arrays, which span the   *
// *             whole tuning_interval; when these values get written to the *
// *             'bound' files, they will be divided by the number of mixes  *
// *             we have done already to average them out (see also King-Wai *
// *             Chu's thesis, p. 63)                                        *
// ***************************************************************************/
void DoTuning(void)
{
	int    i,j,k,l;                                         /* loop counters */

	/* variables for calculating cross-correlation for lower bound */

	double var;                             /* variance over sample interval */
	double correl  = 0.;  /* temp variable for calculating cross-correlation */
	int    ncorrel = 0;   /* counter for cross-correlations we've calculated */

	/* variables for calculating variance of local means for upper bound */

	double avg_l_mean = 0.;  /* the average local mean across all processors */
	double var_l_mean = 0.;  /* temp variable for the var of the local means */



/* For the lower bound of M: ***********************************************
* here we calculate a intermediate product for the calculation of the     *
* cross-correlation between processors; what we end up having in each     *
* element of dev is the local dev / sqrt(var); we can then multi-         *
* ply this for two different processors to obtain the cross correlation   *
* below (according to King-Wai's thesis, pp. 55-56)                       *
***************************************************************************/

	for (j=0; j<sub_tune_interval; j++)
	{
		var = 0.;
		for (k=0; k<sample_size; k++)
			var += dev[j*sample_size+k]*dev[j*sample_size+k];
		var /= sample_size;

		for (k=0; k<sample_size; k++)
			if ( var > 0.0 )
				dev[j*sample_size+k] /= sqrt(var);
			else
				dev[j*sample_size+k] /= DBL_MIN;     /* not quite right, but close */
	}

/***************************************************************************
* distribute devs and means for this sub_tune_interval among all          *
* processors                                                              *
***************************************************************************/

	MPI_Allgather(dev, sample_size*sub_tune_interval, MPI_DOUBLE,
					tot_dev, sample_size*sub_tune_interval, MPI_DOUBLE,
					MPI_COMM_WORLD);

	MPI_Allgather(means, sample_size*sub_tune_interval, MPI_DOUBLE,
					tot_means, sample_size*sub_tune_interval, MPI_DOUBLE,
					MPI_COMM_WORLD);

/* the following stuff is done for each sample over the sub_tune_inteval */

	for (i=0; i<sub_tune_interval; i++)
	{

/***************************************************************************
* this loop takes all sample_size devs from all processors for the        *
* current sample within sub_tune_interval and rearranges them into the    *
* coll arrays, which we use below to calculate cross correlation and      *
* variances of the local means                                            *
***************************************************************************/

		for (j=0; j<nnodes; j++)
			for (k=0; k<sample_size; k++)
			{
				coll_dev[j*sample_size+k] =
					tot_dev[j*sub_tune_interval*sample_size + i*sample_size + k];
				coll_means[j*sample_size+k] =
					tot_means[j*sub_tune_interval*sample_size + i*sample_size + k];
			}


/* For the lower bound of M: ***********************************************
* calculate cross-correlation for those processors that had the same dan- *
* ce partner at the last mix; most of the work has already been done for  *
* us in a loop above, where we calculate the intermediate product that    *
* corresponds to dev / sqrt(var); here, we only need to multiply such     *
* intermediate products for processors that shared the same state after   *
* the previous mix; then save the results in the cross_correl array       *
***************************************************************************/

		correl  = 0.;
		ncorrel = 0;
		for (j=1; j<nnodes; j++)           /* loop through all combinations of */
			for (k=0; k<j; k++)                                    /* processors */
				if (dance_partner[j] == dance_partner[k])
					for (l=0; l<sample_size; l++)
					 {
						correl +=
							coll_dev[j*sample_size+l] *
							coll_dev[k*sample_size+l];
						ncorrel++;
					}

		if (ncorrel > 0)             /* only divide if anything was calculated */
			correl /= (double)ncorrel;

		cross_correl[count_tune*sub_tune_interval+i] += correl;

/* For the upper bound of M: ***********************************************
 * here we calculate the variance of the local mean for the current sample *
 * interval; avg_l_mean is the average local mean from all processors over *
 * the covar_sample; it is then used to calculate the variance of the      *
 * local means (var_l_mean) over the sample interval, which is then stored *
 * in the var_means array                                                  *
  ***************************************************************************/

		avg_l_mean = 0.;
		var_l_mean = 0.;

		for (j=0; j<covar_sample; j++)
			avg_l_mean += coll_means[j];
		avg_l_mean /= (double)covar_sample;

		for (j=0; j<covar_sample; j++)
		{
			if ( avg_l_mean == 0.0 )
				coll_means[j] = (coll_means[j] - DBL_MIN) / DBL_MIN;
			else
				coll_means[j] = (coll_means[j] - avg_l_mean) / avg_l_mean;

			var_l_mean += coll_means[j] * coll_means[j];
		}
		var_l_mean /= (double)covar_sample;

		var_means[count_tune*sub_tune_interval+i] += var_l_mean;

	}

	/* we're done with this sub_tune_interval: on to the next one */

	count_tune++;

	/* if we're at a tune_interval: reset the count_tune counter */

	if ( count_tune == write_tune_stat )
		count_tune = 0;

	/* restart the moves counter */

	moves_tune = 0;

}


/*** WriteTuning: every tune_interval (which basically corresponds to the **
*                mix_interval), we write the cross-correlations and vari- *
*                ances of the local mean to various files; the way this   *
*                is done is a little bizarre: the two .lb and .ub files   *
*                are averaged over all previous tune_intervals and then   *
*                overwritten every time, the .mb file on the other hand   *
*                contains the history of the upper and lower bounds over  *
*                past tune_intervals; here are the details in the words   *
*                of the Great Dr. Chu:                                    *
*                When tuning finishes, there is a file called <infile>.lb *
*                (Fig. 5.1 in thesis) and a file called <infile>.ub (Fig. *
*                5.2). The history of the lower and upper bound values    *
*                during tuning are in <infile>.mb (Fig. 5.3). The average *
*                of the last pair of lower and upper bound values is the  *
*                M from tuning.                                           *
***************************************************************************/

void WriteTuning(void)
{
	int    i,j;                                             /* loop counters */

	/* file pointers */

	FILE   *lbptr;                           /* pointer for lower_bound file */
	FILE   *ubptr;                           /* pointer for upper_bound file */
	FILE   *mbptr;                             /* pointer for mix_bound file */

	/* variables for lower bound */

	double long_avg;         /* average cross_correl over whole mix_interval */

	/* variables for upper bound */

	double min;                     /* these are for storing the minimum and */
	double max;      /* maximum variances of local means over a mix_interval */
	int    group_size = GROUP_SIZE;      /* group for evaluating upper bound */
	int    midpoint;                                  /* midpoint of a group */
	double sum;                               /* summed variances of a group */
	double avg;                                        /* average of a group */

	/* both upper and lower bound information goes into the mix_bound file     *
	   * that's why we keep it open during both writing lower and upper bounds   */

	if ( myid == 0 )
	{
		mbptr = fopen(mbfile, "a");
		if ( !mbptr )
			file_error("WriteTuning");

		 /* first write cross-correlation for lower bounds */

		lbptr = fopen(lbfile, "w");
		if ( !lbptr )
			file_error("WriteTuning");

		fprintf(lbptr, "# tau no.    cross_corr\n\n");

		for (i=0; i<tune_interval; i++)
			fprintf(lbptr, "   %6d   %11.8f\n",
					tau_count[i], cross_correl[i]/(double)count_mix);

		fclose(lbptr);


		/* lower bound: write sample-count for first cross-correlation which is    *
		* below the average cross-correlation over the whole tune_interval to     *
		* the mix-bound file (see King-Wai's thesis, p. 58)                       */

		long_avg = 0.;
		for (i=0; i<tune_interval; i++)
			long_avg += cross_correl[i];
		long_avg /= (double)tune_interval;

		for (i=0; i<tune_interval; i++)
			if (cross_correl[i] <= long_avg)
			{
				fprintf(mbptr, "       %4d   %6d   %11.8f",
						count_mix, tau_count[i], cross_correl[i]/(double)count_mix);
				break;
			}

		/* write variance of local means for upper bounds */

		ubptr = fopen(ubfile, "w");
		if ( !ubptr )
			file_error("WriteTuning");

		fprintf(ubptr, "# tau no.    var_l_mean\n\n");
	}

	/* here we evaluate min and max value of var_means over the whole interval */

	min = var_means[0]/(double)count_mix;
	max = -DBL_MIN;
	for (i=0; i<tune_interval; i++)
	 {
		if ( (var_means[i]/(double)count_mix) > max )
			max = var_means[i]/(double)count_mix;
		if ( myid == 0 )
		fprintf(ubptr, "   %6d   %11.8f\n",
				tau_count[i], var_means[i]/(double)count_mix);
	}

	if ( myid == 0 )
		fclose(ubptr);

 /* upper bounds: evaluate group of variances that differs from the first   *
  * group of variances by more than 7% (see King-Wai's thesis, p. 61)       */

	group_size /= covar_index;             /* set group size for upper bound */
	if (group_size < 1)
		group_size = 1;

	for (i=0; i<tune_interval-group_size; i++)
	{
		midpoint = i + group_size/2;
		sum = 0.;
		for (j=0; j<group_size; j++)
			sum += var_means[i+j];
		avg = sum / ((double)group_size * (double)count_mix);

		if ( fabs( (avg-min)/(max-min) ) >= 0.07 )
			break;
	}

	midpoints[count_mix-1] = midpoint;

	if ( myid == 0 )
	{
		fprintf(mbptr, "   %6d   %11.8f\n", midpoint, avg);
		fclose(mbptr);
	}

	/* reset the sample counter */

	count_sample = 0;

}



/*** StopTuning: is to tuning runs what Frozen() is to a normal annealing **
*               run: it basically checks if the tuning stop criterion     *
*               applies and returns true if that's the case               *
***************************************************************************/

int StopTuning(char * logfile)
{
	int    i;                                                /* loop counter */

	FILE   *logptr;                                  /* pointer to .log file */

	double tol_tune = STOP_TUNE_CRIT;		  /* tuning stop criterion */
	double avg;                             /* temp variable for the average */

	/* calculate average upper bound for the last 'stop_tune_count' mixes */

	avg = 0.;
	for (i=count_mix-stop_tune_count; i<count_mix; i++)
		avg += (double)midpoints[i];
	avg /= (double)stop_tune_count;

	/* check if any upper bound is further than tol_tune away from that avg */

	if ( avg == 0 )
		error("StopTuning: average midpoint was zero!!?!\n");
	else
		for (i=count_mix-stop_tune_count; i<count_mix; i++)
			if ( (fabs((double)midpoints[i] - avg) / avg) >= tol_tune)
				return 0;                           /* if yes, return 'false' */

	/* if we're done: write a message into the log file and quit */

	if (myid == 0)
	{

		logptr = fopen(logfile, "a");     /* write comment to global .log file */
		if ( !logptr )
			file_error("StopTuning");

		fprintf(logptr, "Tuning stops before the end of an annealing run.\n");
		fprintf(logptr, "Therefore, the score and iterations will not be\n");
		fprintf(logptr, "the true final score and iterations.\n");

		fclose(logptr);
	}

	return 1;

}


void UpdateLocalSTuning(double S)
{
	/* do the same for local Lam parameters (for both lower and upper bounds) */

	l_estimate_mean_l = 1.0 / (l_A_l*S + l_B_l);
	l_estimate_mean_u = 1.0 / (l_A_u*S + l_B_u);
	l_estimate_sd     = 1.0 / (l_D  *S + l_E);

}

void InitLocalStatsTuning(double mean, double vari, int success)
{
	l_mean    = mean;
	l_vari    = vari;
	l_success = success;
}

void Init2LocalStatsTuning(int proc_init)
{
	l_mean      /= (double)proc_init;
	l_vari       = l_vari / ((double)proc_init) - l_mean * l_mean;
	l_acc_ratio  = ((double)l_success) / ((double)proc_init);
}


void UpdateLocalStatsTuning(int proc_tau)
{
	l_mean      /=  (double)proc_tau;
	l_vari      /=  (double)proc_tau;
	l_acc_ratio  = ((double)l_success) / ((double)proc_tau);
}

void RestoreLocalLogTuning(int max_saved_count)
{
	int i;
	char   *shell_cmd;                             /* used by 'system' below */
	char   *logline;                                /* array of read buffers */
	char   *l_outfile;                         /* temporary output file name */
	FILE   *l_logptr;                            /* local .llog file pointer */
	FILE   *l_outptr;

	logline   = (char *)calloc(MAX_RECORD, sizeof(char));
	shell_cmd = (char *)calloc(MAX_RECORD, sizeof(char));

	l_outfile = (char *)calloc(MAX_RECORD, sizeof(char));
	l_outfile = strcpy(l_outfile,"llogXXXXXX");   /* required by mkstemp() */
	if ( mkstemp(l_outfile) == -1 )       /* get unique name for temp file */
		error("RestoreLog: error creating temporary (llog) file");



	int saved_count_tau = -1;                             /* reset the counter */

	/*open the log file and a temporary output file */

	l_logptr = fopen(l_logfile, "r");
	if ( !l_logptr )
		file_error("RestoreLog (at open llog file for reading)");

	l_outptr = fopen(l_outfile, "w");
	if ( !l_outptr )
		file_error("RestoreLog (at open temp llog file for writing)");

	/* read and write the first few title and caption lines */

	for (i=0; i<4; i++)
	{
		if ( NULL == fgets(logline, MAX_RECORD, l_logptr ) )
			error("RestoreLog: error reading llog file captions");
		fprintf(l_outptr, "%s", logline);
	}

	/* read and write the actual llog lines till we are at current time */

	while ( (saved_count_tau < max_saved_count) &&
			(NULL != fgets(logline, MAX_RECORD, l_logptr)) )
	{
		if ( 1 != sscanf(logline, "%d", (int *) &saved_count_tau) )
			error("RestoreLog: error reading saved_count_tau (after %d)",
				  saved_count_tau);
		fprintf(l_outptr, "%s", logline);
	}

	fclose(l_logptr);
	fclose(l_outptr);

	/* rename tmpfile into new file */

	sprintf(shell_cmd, "cp -f %s %s", l_outfile, l_logfile);

	if ( -1 == system(shell_cmd) )
		error("RestoreLog: error renaming temp file %s", l_outfile);

	if ( remove(l_outfile) )
		warning("RestoreLog: temp file %s could not be deleted",
				l_outfile);

	free(shell_cmd);
	free(logline);

}

void CalculateLocalStats(double energy)
{
	double d;                /* difference between energy and estimated mean */

	l_mean  += energy;
	d        = energy - l_estimate_mean_u;
	l_vari  += d * d;

	/* Collect samples for tuning: local estimated deviations (for calculating *
	 * cross-correlation between processors for the lower bound of M) and lo-  *
	 * cal estimated means (for calculating the variance of those local means  *
	 * between processors for the upper bound of M); here, we also need to     *
	 * keep track of the moves we've done in the current sub_tune_interval and *
	 * of the number of samples we've collected in the tune_interval           */

	dev[moves_tune]   = energy - l_estimate_mean_l;
	means[moves_tune] = l_estimate_mean_u;
	moves_tune++;
	if ( (moves_tune % sample_size) == 0 )
		count_sample++;
}

int UpdateTuning(char * logfile)
{
	/* Do tuning every sub_tune_interval, only after first mix */

	if ( (count_sample % sub_tune_interval) == 0 && count_sample > 0 )
	{
		if (count_mix > 0)
		{
			DoTuning();
		}
		else
			moves_tune = 0;
	}
	/* At the end of each tune interval: ***************************************
	 * 1. Root process writes tuning every tune_interval, only after first mix *
	 *    NOTE: this needs to get done BEFORE mixing, otherwise count_mix      *
	 *          isn't right and you'll get a division by zero in WriteTuning() *
	 * 2. If we've done more that 'stop_tune_count' mixes, we check if we can  *
	 *    stop tuning run (runs can be forced to continue using the -S option) *
	 ***************************************************************************/

	if ( (count_sample % tune_interval) == 0 && count_sample > 0 )
	{

		if ( count_mix > 0 )
			WriteTuning();

		if ( (count_mix >= stop_tune_count) && auto_stop_tune )
			if ( StopTuning(logfile) )
			{
				// FreeLocalVariables();
				free(dev);
				free(tot_dev);
				free(coll_dev);
				free(means);
				free(tot_means);
				free(coll_means);
				free(tau_count);
				free(cross_correl);
				free(var_means);
				free(midpoints);

				return 1;
			}
	}
	return 0;
}

/*** InitTuning: sets up/restores structs and variables for tuning runs ****
 ***************************************************************************/

void InitTuning(SAType * state)
{
	int     i;                                               /* loop counter */
	FILE    *mbptr;                            /* pointer for mix_bound file */

	/* Initializing variables */
	covar_index     = state->tuning_settings->covar_index;      /* covariance sample will be covar_index * tau */
	write_tune_stat = state->tuning_settings->write_tune_stat;         /* how many times do we write tuning stats? */
	auto_stop_tune  = state->tuning_settings->auto_stop_tune;               /* auto stop tuning runs? default: on */
	/* error check */

	if ( nnodes <= 1 )
		error("plsa: tuning does not make sense on one processor");

	if ( covar_index > state->mix_interval )
		error("plsa: you can't sample over more than the whole mix interval");

	/* set size of sample interval */

	covar_sample = covar_index * state->tau;

	/* set size of sample interval per processor */
	/* by the way: sample_size corresponds to covar_index * proc_tau */

	if ( (covar_sample % nnodes) != 0 )
		error("plsa: covar_sample (%d) not divisible by nnodes", covar_sample);

	sample_size = covar_sample / nnodes;

	/* tune_interval: how many covar_samples per mix_interval? */
	if ( (state->mix_interval % covar_index) != 0 )
		error("plsa: mix interval (%d) not divisible by covar_index",
			  state->mix_interval);

	tune_interval = state->mix_interval / covar_index;

	/* size of every tune interval in between writing tuning stats */

	if ( write_tune_stat > tune_interval )
		error("plsa: freq of writing tune stats (%d) must be smaller than %d",
			  write_tune_stat, tune_interval);

	if ( (tune_interval % write_tune_stat) != 0 )
		error("plsa: tune_interval (%d) not divisible by write_tune_stat (%d)",
			  tune_interval, write_tune_stat);

	sub_tune_interval = tune_interval / write_tune_stat;

	/* allocate memory for various tuning-specific arrays */

	dev = (double *)calloc(sub_tune_interval*sample_size, sizeof(double));
	tot_dev =
		(double *)calloc(sub_tune_interval*covar_sample, sizeof(double));
	coll_dev = (double *)calloc(covar_sample, sizeof(double));

	means = (double *)calloc(sub_tune_interval*sample_size, sizeof(double));
	tot_means =
		(double *)calloc(sub_tune_interval*covar_sample, sizeof(double));
	coll_means = (double *)calloc(covar_sample, sizeof(double));

	tau_count    =    (int *)calloc(tune_interval, sizeof(int));
	cross_correl = (double *)calloc(tune_interval, sizeof(double));
	var_means    = (double *)calloc(tune_interval, sizeof(double));

	midpoints    =    (int *)calloc(MAX_MIX, sizeof(int));

	/* initialize the arrays used to collect tuning information */

	for (i=0; i<tune_interval; i++)
	{
		tau_count[i]    = (i+1)*covar_index;
		cross_correl[i] = 0.;
		var_means[i]    = 0.;
	}

	/* initialize the mb file */

	if ( myid == 0 )
	{
		mbptr = fopen(mbfile, "w");
		if ( !mbptr )
			file_error("InitTuning (writing mb captions)");
		fprintf(mbptr, "# mix_count  M_lower    cross_corr");
		fprintf(mbptr, "  M_upper    var_l_mean\n\n");
		fclose(mbptr);
	}
}
