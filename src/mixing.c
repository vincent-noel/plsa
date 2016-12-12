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
#include "mixing.h"

#include "error.h"
#include "random.h"
#include "state.h"
#include "tuning.h"

// static int covar_index;    /* covariance sample index for tuning (in 'tau' units) */
// static int write_tune_stat;               /* how often to write tuning statistics */
// static int auto_stop_tune;       /* auto stop tune flag to stop tuning runs early */
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

static double l_alpha;       /* the third term of the Lam schedule formula */
static double l_acc_ratio;  /* average acceptance ratio for all parameters */


/* parallel code: variables for local sd estimator (see also comment on    *
 *                local weights above)                                     */

static double l_w_b;                         /* l_w_b: local weight for sd */

static double l_vsyy;       /* these parameters store intermediate results */
static double l_vsxy;         /* for the updating formulas for l_D and l_E */
static double l_vsxx;
static double l_vsx;
static double l_vsy;
static double l_vsum;

static double l_D;    /* l_D_l and l_E_l: parameters for the rational func */
static double l_E;     /* for the estimation of the sd for the lower bound */

/* parallel code: variables for calculating local mean estimator for lower *
 *                bound of M (see also comment on local weights above)     */

static double l_w_a_l; /* l_w_a_l: weight for the mean for the lower bound */

static double l_usyy_l;     /* these parameters store intermediate results */
static double l_usxy_l;   /* for the updating formulas for l_A_l and l_B_l */
static double l_usxx_l;
static double l_usx_l;
static double l_usy_l;
static double l_usum_l;

static double l_A_l;  /* l_A_l and l_B_l: parameters for the rational func */
static double l_B_l; /* for the estimation of the mean for the lower bound */

/* parallel code: variables for calculating local mean estimator for upper *
 *                bound of M (see also comment on local weights above)     */

static double l_w_a_u; /* l_w_a_u: weight for the mean for the upper bound */

static double l_usyy_u;     /* these parameters store intermediate results */
static double l_usxy_u;   /* for the updating formulas for l_A_u and l_B_u */
static double l_usxx_u;
static double l_usx_u;
static double l_usy_u;
static double l_usum_u;

static double l_A_u;        /* A and B are the parameters for the rational */
static double l_B_u;            /* function for the estimation of the mean */

/* vars used for tuning runs ***********************************************/

/* general variables that determine tune- and sub_tune_interval */

/* various counter variables and an array to save counts */

// static int    count_sample = 0;        /* how many samples did we collect? */
// static int    count_tune   = 0;      /* how many times did we do sub_tune? */
static int    count_mix    = 0;      /* counts the times we've been mixing */

// static int    moves_tune = 0;   /* move counter: reset every tune_interval */

static int    *tau_count; /* # of tau's we've done for whole tune_interval */


/* a variable for the tuning stop criterion */

static char   *l_logfile;                  /* name of the local .llog file */

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


void InitializeMixing()
{
	dance_partner = (int *)calloc(nnodes, sizeof(int));
}

void FreeLocalVariables()
{
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

}
void InitLocalFilenames()
{
	l_logfile = (char *)calloc(MAX_RECORD, sizeof(char));

	if (myid == 0)
		InitLocalTuningFilenames();

	/* the local .llog file: used to store iterations, temperature, tempera-   *
	 * ture change, local mean and stamdard deviation, local Lam estimators    *
	 * for mean and sd (for the upper bound) and local acceptance ratios;      *
	 * this is only needed when tuning, otherwise we just write one global log */

	if ( nnodes <= 10 )
		sprintf(l_logfile,   "plsa_%d.llog", myid);
	else if ( (nnodes > 10) && (nnodes <= 100) )
		sprintf(l_logfile, "plsa_%02d.llog", myid);
	else if ( (nnodes > 100) && (nnodes <= 1000) )
		sprintf(l_logfile, "plsa_%03d.llog", myid);
	else if ( (nnodes > 1000) && (nnodes <= 10000) )
		sprintf(l_logfile, "plsa_%04d.llog", myid);
	else if ( (nnodes > 10000) && (nnodes <= 100000) )
		sprintf(l_logfile, "plsa_%05d.llog", myid);
	else
		error("Initialize: can't open more than 100'000 llog files");
}

void InitializeLocalParameters(double S_0)
{
	/* 2. set local parameters *************************************************
	 * these are only needed when tuning...                                    *
	 * note that the two sets of local estimators used for tuning only differ  *
	 * in their weights, so they can actually be initialized the same way      */

	 double   l_d;               /* l_d is used to store intermediate results */

	/* set estimators to stats collected during initializing phase of run */

	l_estimate_sd                         = sqrt(l_vari);
	l_estimate_mean_l = l_estimate_mean_u = l_mean;

	/* initialize local A,B,D,E */

	l_A_l = l_A_u = l_estimate_sd * l_estimate_sd /
					(l_estimate_mean_u * l_estimate_mean_u);
	l_B_l = l_B_u = (1.0 / l_estimate_mean_u) - (l_A_u * S_0);

	l_D           = l_estimate_sd / l_estimate_mean_u;
	l_E           = (1.0 / l_estimate_sd) - (l_D * S_0);

	/* initialize these intermediate variables for updating funcs for A,B,D,E */

	l_usum_l = l_usum_u = l_vsum = 1.0;

	l_usxy_l = l_usxy_u = l_usxx_l = l_usxx_u = l_usx_l = l_usx_u = 0.0;
	l_usy_l  = l_usy_u  = 1.0 / l_estimate_mean_u;
	l_usyy_l = l_usyy_u = l_usy_u * l_usy_u;

	l_vsxy              = l_vsxx              = l_vsx             = 0.0;
	l_vsy               = 1.0 / l_estimate_sd;
	l_vsyy              = l_vsy * l_vsy;

	/* alpha is the third term of the main Lam schedule formula */

	l_d     = (1.0 - l_acc_ratio) / (2.0 - l_acc_ratio);
	l_alpha = 4.0 * l_acc_ratio * l_d * l_d;

}

void InitializeLocalWeights(double w_a, double w_b)
{
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

	/* l_w_a_{l|u} are the local weights for the mean */

	l_w_a_l = w_a / (double)nnodes;
	l_w_a_u = w_a;

	/* l_w_b is the local weight for the standard deviation */

	l_w_b   = w_b;
}

double * InitLocalStats(double mean, double vari, int * success,
							int i, int initial_moves)
{
	double      *total;        /* will hold mean & vari when pooling stats */
	double      tmptotal[2];                    /* temporary array for above */

	long        nodesuccess[2];   /* array for success and summed init moves */
	long        tmpsuccess[2];                  /* temporary array for above */

	/* parallel code: pool initial statistics from all nodes */
	total = (double *) malloc(sizeof(double)*2);
	total[0]       = mean;
	total[1]       = vari;
	nodesuccess[0] = (long)(*success);
	nodesuccess[1] = (long)i;         /* sum of i's must equal initial moves */

	/* the tmp arrays are used to send messages */

	tmptotal[0]    = total[0];
	tmptotal[1]    = total[1];
	tmpsuccess[0]  = nodesuccess[0];
	tmpsuccess[1]  = nodesuccess[1];

	/* stats from all nodes are summed up here */

	MPI_Allreduce(tmptotal, total, 2, MPI_DOUBLE, MPI_SUM,
				  MPI_COMM_WORLD);
	MPI_Allreduce(tmpsuccess, nodesuccess, 2, MPI_LONG, MPI_SUM,
				  MPI_COMM_WORLD);

	*success = (int)nodesuccess[0];         /* success is now global success! */

	/* sanity check: have we done the correct number of initial moves? */

	if ( nodesuccess[1] - initial_moves )
		error("InitialLoop: initial moves was %d?!\n", nodesuccess[1]);


	return total;
}


double * UpdateLocalStats(double mean, double vari, int i, int tau, int * success)
{
	double * total;             /* will hold mean & vari when pooling stats */
	double tmptotal[2];	                      /* temporary array for above */

	long   nodesuccess[2];             /* array for success and summed moves */
	long   tmpsuccess[2];                       /* temporary array for above */

	/* parallel code: pool statistics from all nodes */
	total = (double *) malloc(sizeof(double)*2);
	total[0] = mean;
	total[1] = vari;

	nodesuccess[0] = (long)(*success);
	nodesuccess[1] = (long)i;       /* sum of i's must equal total tau moves */

	tmptotal[0] = total[0];      /* the tmp arrays are used to send messages */
	tmptotal[1] = total[1];

	tmpsuccess[0] = nodesuccess[0];
	tmpsuccess[1] = nodesuccess[1];

	/* stats from all nodes are summed up here */

	MPI_Allreduce(tmptotal, total, 2, MPI_DOUBLE, MPI_SUM,
				  MPI_COMM_WORLD);
	MPI_Allreduce(tmpsuccess, nodesuccess, 2, MPI_LONG, MPI_SUM,
				  MPI_COMM_WORLD);

	total[0] /= ((double) tau);            /* need to divide mean and vari by Tau steps */
	total[1] /= ((double) tau);

	*success = (int)nodesuccess[0];         /* success is now global success! */

	/* sanity check: have we done the correct number of initial moves? */
	if (nodesuccess[1] - tau)
		error("Loop: total moves was %d ?!\n", nodesuccess[1]);

	return total;
}

void ResetLocalStats()
{
	l_mean     = 0.0;
	l_vari     = 0.0;
	l_success  = 0;
}

void AddLocalSuccess()
{
	l_success++;
}
/*** UpdateLParameter: update local parameters l_A_{l|u}, l_B_{l|u}, l_D   *
 *                     and l_E and the local estimators for mean and stan- *
 *                     dard deviation for both upper and lower bounds of M *
 *                     when tuning                                         *
 ***************************************************************************/

void UpdateLParameter(double S)
{
	register double d;      /* d is used to store intermediate results below */

	/* update the estimator for the local mean */

	d = 1.0 / l_mean;

	/* first: multiply all intermediate vars by weights */
	/* lower bound variables: */

	l_usyy_l *= l_w_a_l;
	l_usxy_l *= l_w_a_l;
	l_usy_l  *= l_w_a_l;
	l_usx_l  *= l_w_a_l;
	l_usxx_l *= l_w_a_l;
	l_usum_l *= l_w_a_l;

	/* upper bound variables: */

	l_usyy_u *= l_w_a_u;
	l_usxy_u *= l_w_a_u;
	l_usy_u  *= l_w_a_u;
	l_usx_u  *= l_w_a_u;
	l_usxx_u *= l_w_a_u;
	l_usum_u *= l_w_a_u;

	/* then: update all intermediate vars */
	/* lower bound variables: */

	l_usyy_l += d*d;
	l_usxy_l += S*d;
	l_usy_l  += d;
	l_usx_l  += S;
	l_usxx_l += S*S;
	l_usum_l += 1.0;

	/* upper bound variables: */

	l_usyy_u += d*d;
	l_usxy_u += S*d;
	l_usy_u  += d;
	l_usx_u  += S;
	l_usxx_u += S*S;
	l_usum_u += 1.0;

	/* ... and use intermediate vars to update l_A_{l|u} and l_B_{l|u} ... */

	/* lower bound variables: */

	l_A_l = (l_usum_l * l_usxy_l - l_usx_l * l_usy_l) /
			(l_usum_l * l_usxx_l - l_usx_l * l_usx_l);
	l_B_l = (l_usy_l - l_A_l * l_usx_l) / l_usum_l;

	/* upper bound variables: */

	l_A_u = (l_usum_u * l_usxy_u - l_usx_u * l_usy_u) /
			(l_usum_u * l_usxx_u - l_usx_u * l_usx_u);
	l_B_u = (l_usy_u - l_A_u * l_usx_u) / l_usum_u;

	/* ... which are then used to update the local estimators for the mean */

	l_estimate_mean_l = 1.0 / (l_A_l * S + l_B_l);
	l_estimate_mean_u = 1.0 / (l_A_u * S + l_B_u);

	/* update the local estimator for the standard deviation */

	if ( l_vari > 0.0 )
	{

		d = 1.0 / sqrt(l_vari);

		/* first: multiply all intermediate vars by weights */

		l_vsyy *= l_w_b;
		l_vsxy *= l_w_b;
		l_vsy  *= l_w_b;
		l_vsx  *= l_w_b;
		l_vsxx *= l_w_b;
		l_vsum *= l_w_b;

		/* then: update all intermediate vars */

		l_vsyy += d*d;
		l_vsxy += S*d;
		l_vsy  += d;
		l_vsx  += S;
		l_vsxx += S*S;
		l_vsum += 1.0;

		/* ... and use intermediate vars to update l_D and l_E ... */

		l_D = (l_vsum*l_vsxy - l_vsx*l_vsy) / (l_vsum*l_vsxx - l_vsx*l_vsx);
		l_E = (l_vsy - l_D*l_vsx) / l_vsum;
	}

	/* ... which are then used to update the local estimator for the std dev */

	l_estimate_sd = 1.0 / (l_D*S + l_E);

	/* alpha corresponds to the third term in the main Lam schedule formula,   *
	 * which is a measure of how efficiently the space state is samples; this  *
	 * term is at a maximum for acc_ratio = 0.44, see Lam & Delosme, 1988b, p1 */

	d = (1.0 - l_acc_ratio) / (2.0 - l_acc_ratio);
	l_alpha = 4.0 * l_acc_ratio * d * d;

}
/*** MIXING ****************************************************************/

/*** DoMix: does the mixing; sends move state and local Lam stats to the ***
 *          dance partner(s)                                               *
 ***************************************************************************/

double DoMix(double energy, double estimate_mean, double S, int tuning)
{
	int    i;                                                /* loop counter */

	/* variables needed for evaluating the dance partners; note that the dance *
	 * partner array is static to lsa.c, since it's also needed by tuning code */

	double prob;      /* probability of choosing a node's energy upon mixing */
	double norm;                      /* sum used to normalize probabilities */
	double theirprob;                    /* probability of the dance partner */
	double psum;         /* sum of probabilities for a certain dance partner */
	double *node_prob;                    /* all probabilities for all nodes */
	double l_energy = energy;
	/* buffers for sending/receiving Lam stats at mix */

	double *sendbuf;
	double *recvbuf;

	/* buffers for send/recv move state from/to move(s).c (incl. sizes) */

	int    lsize;
	int    dsize;

	long   *recv_longbuf;
	double *recv_doublebuf;

	long   *send_longbuf;
	double *send_doublebuf;

	/* MPI status and handle (request) arrays */

	MPI_Status    *status;       /* status array of the MPI_Waitall function */
	MPI_Request   *request;           /* handle array for receiving messages */

	/* allocate probability array and arrays for MPI_Waitall below */

	node_prob = (double *)calloc(nnodes, sizeof(double));

	request   = (MPI_Request *)calloc(3, sizeof(MPI_Request));
	status    =  (MPI_Status *)calloc(3, sizeof(MPI_Status));

	/* initialize probability & dance partner arrays */


	for (i=0; i<nnodes; i++)
	{
		dance_partner[i] = 0;       /* static to lsa.c since needed for tuning */
		node_prob[i]     = 0.;
	}


	/* update mix counter (used by tuning code only) */

	if ( tuning )
		count_mix++;              /* counts the number of mixings we have done */

	/* calculate probabilities for accepting local state upon mixing; we sub-  *
	 * tract estimate_mean (which is the same in all processes from the total  *
	 * energy to avoid overflows; see also Chu (2001, p.41) for details        */

	prob = exp((estimate_mean - l_energy)*S);

	if (prob < DBL_MIN) /* set some limits to prob to avoid under-/overflows */
		prob = DBL_MIN;
	else if (prob >= HUGE_VAL)
		prob = DBL_MAX/nnodes;


	/* sum up probabilities from all nodes for normalization, then normalize */

	MPI_Allreduce(&prob, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	prob /= norm;
	/* gather local probabilities from all nodes */

	MPI_Allgather(&prob, 1, MPI_DOUBLE, node_prob, 1, MPI_DOUBLE,
				  MPI_COMM_WORLD);

	/* theirprob determines the dance partner we choose */

	if (nnodes > 1)
		theirprob = RandomReal();
	else if (nnodes == 1)
		theirprob = 0.0;         /* guarantee same Markov proc. as serial case */
	else
		error("DoMix: you can't compute on %d nodes!", nnodes);

	/* sum up probabilities to determine dance partner */

	psum = 0.;
	for (i=0; i<nnodes; i++)
	{
		psum += node_prob[i];
		if (psum > theirprob)
			break;
	}
	/* pool dance partners */

	MPI_Allgather(&i, 1, MPI_INT, dance_partner, 1, MPI_INT, MPI_COMM_WORLD);

	/* get move state from move(s).c and collect local Lam stats for sending */

	MakeStateMsg(&send_longbuf, &lsize, &send_doublebuf, &dsize);
	MakeLamMsg(&sendbuf, l_energy, tuning);



	/* if I'm not dancing with myself: receive new state and Lam stats         */

	if (dance_partner[myid] != myid && dance_partner[myid] < nnodes)
	{

		/* allocate receive buffers for message; receive buffers need to be of the *
		 * same length as the sending buffers                                      */

		recv_longbuf   =   (long *)calloc(lsize,        sizeof(long));
		recv_doublebuf = (double *)calloc(dsize,        sizeof(double));
		if ( tuning && nnodes>1 )
			recvbuf      = (double *)calloc(LSTAT_LENGTH_TUNE, sizeof(double));
		else
			recvbuf      = (double *)calloc(LSTAT_LENGTH, sizeof(double));

		/* MPI_Irecv is a non-blocking receive which allows a process that is      *
		 * waiting to receive to send out its message while waiting i.e. we post   *
		 * the receive now, check for reception after send (using MPI_Waitall)     */

		MPI_Irecv(recv_doublebuf, dsize, MPI_DOUBLE, dance_partner[myid],
				  dance_partner[myid], MPI_COMM_WORLD, &request[0]);
		MPI_Irecv(recv_longbuf,   lsize, MPI_LONG,   dance_partner[myid],
				  dance_partner[myid], MPI_COMM_WORLD, &request[1]);

		/* the length of the Lam statistics message depends on if we're tuning or  *
		 * not (local stats are only sent when tuning)                             */

		if ( tuning && nnodes>1 )
			MPI_Irecv(recvbuf, LSTAT_LENGTH_TUNE, MPI_DOUBLE, dance_partner[myid],
					  dance_partner[myid], MPI_COMM_WORLD, &request[2]);
		else
			MPI_Irecv(recvbuf, LSTAT_LENGTH,      MPI_DOUBLE, dance_partner[myid],
					  dance_partner[myid], MPI_COMM_WORLD, &request[2]);

	}

	/* send messages to dance partners, if requested */

	for (i=0; i<nnodes; i++)
		if ( (dance_partner[i] == myid) && (i != myid) && (dance_partner[i] < nnodes))
		{
			MPI_Send(send_doublebuf, dsize, MPI_DOUBLE, i, myid, MPI_COMM_WORLD);
			MPI_Send(send_longbuf,   lsize, MPI_LONG,   i, myid, MPI_COMM_WORLD);
			if ( tuning && nnodes>1 )
				MPI_Send(sendbuf, LSTAT_LENGTH_TUNE, MPI_DOUBLE, i, myid,
						 MPI_COMM_WORLD);
			else
				MPI_Send(sendbuf, LSTAT_LENGTH,      MPI_DOUBLE, i, myid,
						 MPI_COMM_WORLD);
		}

	/* if I'm not dancing with myself, we need a new state; MPI_Waitall will   *
	 * collect the three messages we need to receive; then we install the move *
	 * state in move(s).c and the Lam stats in lsa.c                           */

	if (dance_partner[myid] != myid && dance_partner[myid] < nnodes)
	{

		MPI_Waitall(3, request, status);

		AcceptStateMsg(recv_longbuf, recv_doublebuf);
		l_energy = AcceptLamMsg(recvbuf, tuning);

	}


	/* clean up message buffers and MPI arrays ... */

	free(send_longbuf);
	free(send_doublebuf);
	free(sendbuf);

	free(request);
	free(status);

	/* ... and the probability array */

	free(node_prob);
	return l_energy;
}


/*** MakeLamMsg: packages local Lam stats into send buffer *****************
 ***************************************************************************/

void MakeLamMsg(double **sendbuf, double energy, int tuning)
{
	if ( tuning && nnodes>1 )
		*sendbuf = (double *)calloc(LSTAT_LENGTH_TUNE, sizeof(double));
	else
		*sendbuf = (double *)calloc(LSTAT_LENGTH, sizeof(double));

	(*sendbuf)[0]   = energy;

	if ( tuning && nnodes>1 )
	{

		(*sendbuf)[1]   = l_estimate_mean_l;
		(*sendbuf)[2]   = l_estimate_mean_u;
		(*sendbuf)[3]   = l_estimate_sd;

		(*sendbuf)[4]   = l_usyy_l;
		(*sendbuf)[5]   = l_usxy_l;
		(*sendbuf)[6]   = l_usy_l;
		(*sendbuf)[7]   = l_usx_l;
		(*sendbuf)[8]   = l_usxx_l;
		(*sendbuf)[9]   = l_usum_l;

		(*sendbuf)[10]  = l_usyy_u;
		(*sendbuf)[11]  = l_usxy_u;
		(*sendbuf)[12]  = l_usy_u;
		(*sendbuf)[13]  = l_usx_u;
		(*sendbuf)[14]  = l_usxx_u;
		(*sendbuf)[15]  = l_usum_u;

		(*sendbuf)[16]  = l_A_l;
		(*sendbuf)[17]  = l_B_l;

		(*sendbuf)[18]  = l_A_u;
		(*sendbuf)[19]  = l_B_u;

		(*sendbuf)[20]  = l_vsyy;
		(*sendbuf)[21]  = l_vsxy;
		(*sendbuf)[22]  = l_vsy;
		(*sendbuf)[23]  = l_vsx;
		(*sendbuf)[24]  = l_vsxx;
		(*sendbuf)[25]  = l_vsum;

		(*sendbuf)[26]  = l_D;
		(*sendbuf)[27]  = l_E;

	}
}



/*** AcceptLamMsg: receives new energy and Lam stats upon mixing ***********
 ***************************************************************************/

double AcceptLamMsg(double *recvbuf, int tuning)
{
	double energy = recvbuf[0];

	if ( tuning && nnodes>1 )
	{

		l_estimate_mean_l = recvbuf[1];
		l_estimate_mean_u = recvbuf[2];
		l_estimate_sd     = recvbuf[3];

		l_usyy_l = recvbuf[4];
		l_usxy_l = recvbuf[5];
		l_usy_l  = recvbuf[6];
		l_usx_l  = recvbuf[7];
		l_usxx_l = recvbuf[8];
		l_usum_l = recvbuf[9];

		l_usyy_u = recvbuf[10];
		l_usxy_u = recvbuf[11];
		l_usy_u  = recvbuf[12];
		l_usx_u  = recvbuf[13];
		l_usxx_u = recvbuf[14];
		l_usum_u = recvbuf[15];

		l_A_l = recvbuf[16];
		l_B_l = recvbuf[17];

		l_A_u = recvbuf[18];
		l_B_u = recvbuf[19];

		l_vsyy = recvbuf[20];
		l_vsxy = recvbuf[21];
		l_vsy  = recvbuf[22];
		l_vsx  = recvbuf[23];
		l_vsxx = recvbuf[24];
		l_vsum = recvbuf[25];

		l_D = recvbuf[26];
		l_E = recvbuf[27];

	}

	free(recvbuf);
	return energy;
}



void WriteLocalLog(int count_steps, double S, double dS)
{
	FILE   *l_logptr;             /* file pointer for local log file (.llog) */

	l_logptr = fopen(l_logfile, "a");   /* then do the same for .llog file */
	if ( !l_logptr )
		file_error("WriteLog");

	PrintLocalLog(l_logptr, S, count_steps, dS);
	fclose( l_logptr );
}


/*** PrintLog: actually prints the log to wherever it needs to be printed **
***************************************************************************/

void PrintLocalLog(FILE *outptr, int count_steps, double S, double dS)
{

	int t_secs_since_epoc = time(NULL);

	const char *format =
		"%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n";

	/* print data */

	fprintf(outptr, format,
		count_steps,
		1.0/S, dS/S,
		l_mean, sqrt(l_vari), l_estimate_mean_u, l_estimate_sd,
		l_acc_ratio, l_alpha,t_secs_since_epoc);

}
