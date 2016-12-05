/******************************************************************************
 *                                                                            *
 *   plsa.c                                                                   *
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
 *   Although main() is in lsa.c, this is the file that 'defines'             *
 *   the plsa program, since it contains most of its problem-                 *
 *   specific code (except for move generation -> moves.c, saving             *
 *   of intermediate state files -> state.c and communication             *
 *   with the specific cost function that is used -> translate.c).            *
 *                                                                            *
 *   After I've told you all that's NOT in this file, here's what             *
 *   the funcs below actually do: parsing plsa command line opts              *
 *   is one of its jobs; there are funcs that make the first and              *
 *   last moves and funcs that read and write Lam and Lam-indepen-            *
 *   dent annealing parameters to the problem-specific data file.             *
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

#include <float.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <sys/times.h>
#include <unistd.h>          /* for command line option stuff and access() */
#include <time.h>

#include "error.h"									  /* error handling funcs */
#include "config.h"								  /* for olddivstyle and such */
#include "random.h"										   /* for InitERand() */

// #include "distributions.h"               /* DistP.variables and prototypes */
#include "moves.h"						  /* problem-specific annealing funcs */
#include "sa_shared.h"				   /* problem-independent annealing funcs */
#include "score.h"								  /* for init and Score funcs */
#include "sa.h"
#include "plsa.h"
#include "state.h"
#ifdef MPI                 /* this inludes parallel-specific stuff for MPI */
#include <mpi.h>                     /* this is the official MPI interface */
#include "MPI.h"  /* our own structs and such only needed by parallel code */
#endif


/*** STATIC VARIABLES ******************************************************/

static char version[MAX_RECORD];                 /* version gets set below */

/* other static variables */

// static int    precision   = 16;                    /* precision for eqparms */
static PArrPtr 	 	plsa_params;
static DistParms		dist_params;
static SAType			state;
static Opts			options;
static StopStyle   	stop_flag;               /* type of stop criterion (see above) */
static SALogs			logs;
static int    			stateflag;                              /* state file or not? */
static char *	statefile;                        /* name of the state file */

static double energy;                                    /* current energy */

/* variables used for timing */

/* these variables are used to evaluate real time using time() */
static double      start;                     /* wallclock time before run */
// static double      finish;                     /* wallclock time after run */

/* these structs are used to evaluate user time using times() or MPI_Wtime */
static struct tms *cpu_start;                      /* user time before run */
static struct tms *cpu_finish;                      /* user time after run */

/*** FUNCTIONS *************************************************************/

/*** COMMAND LINE OPTS ARE ALL PARSED HERE *********************************/

/*** ParseCommandLine: well, parses the command line and returns an index **
 *                     to the 1st argument after the command line options  *
 ***************************************************************************/
#define VERS 0.1


/*** Initialize: calls ParseCommandLine first; then does either initial ****
 *               randomization and collecting Lam stats or restores state  *
 *               of the annealer as saved in the state file.               *
 ***************************************************************************/
void InitializePLSA()
{

#ifdef MPI
	// tuning          = 0;                    /* tuning mode is off by default */

	// write_llog      = 0; /* write local llog files when tuning; default: off */
#endif

	stateflag = 0;                              /* state file or not? */
#ifdef MPI
	int    flagsum;                       /* used for state file check below */
#endif


	/* allocate memory for static file names */
	statefile = (char *)calloc(MAX_RECORD, sizeof(char));

	/* state files: used for the case that a run terminates or crashes unex-   *
	 * pectedly; we can then restore the state of the run *precisely* as it    *
	 * was before the crash by restarting it from the state file               */
#ifdef MPI

	if ( nnodes>1 )
	{
		if ( nnodes <= 10 )
			sprintf(statefile,   "plsa_%d.state", myid);
		else if ( (nnodes > 10) && (nnodes <= 100) )
			sprintf(statefile, "plsa_%02d.state", myid);
		else if ( (nnodes > 100) && (nnodes <= 1000) )
			sprintf(statefile, "plsa_%03d.state", myid);
		else if ( (nnodes > 1000) && (nnodes <= 10000) )
			sprintf(statefile, "plsa_%04d.state", myid);
		else if ( (nnodes > 10000) && (nnodes <= 100000) )
			sprintf(statefile, "plsa_%05d.state", myid);
		else
			error("Initialize: can't open more than 100'000 state files");
	}
	else
		sprintf(statefile, "plsa.state");

#else
	sprintf(statefile, "plsa.state");
#endif

	/* check if a state file exists (access() is in unistd.h) */

	if ( 0 == access(statefile, F_OK) )
		stateflag = 1;

#ifdef MPI
	/* parallel code: make sure that all state files are present */

	MPI_Allreduce(&stateflag, &flagsum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	if ( (flagsum > 0) && (flagsum != nnodes) && (stateflag == 0) )
		error("Initialize: state file for process %d is missing");
#endif


}

void SetDefaultOptions()
{
	/* set the version string */

#ifdef MPI
	sprintf(version, "plsa version %f parallel", VERS);
#else
	sprintf(version, "plsa version %f serial", VERS);
#endif

	/* following part sets default values for command line options */

	options.captions        = 100000000;  /* default freq for writing captions (off) */
	options.print_freq      = 1;           /* default freq for writing to log file */
	options.state_write     = 100;        /* default freq for writing to state file */

	options.stop_flag       = absolute_freeze;             /* type of stop criterion */
	options.time_flag       = 0;                                  /* flag for timing */
	options.log_flag        = 0;              /* flag for writing logs to the screen */
	options.max_iter        = 0;
	options.max_seconds     = 0;
	options.quenchit = 0;          /* flag for quenchit mode (T=0 immediately) */
#ifdef MPI
	options.tuning = 0;
	options.covar_index     = 1;      /* covariance sample will be covar_index * tau */
	options.write_tune_stat = 1;         /* how many times do we write tuning stats? */
	options.auto_stop_tune  = 1;               /* auto stop tuning runs? default: on */
#endif


	state.seed = -6.60489e+08;
	state.initial_temp = 1000;

	state.lambda = 0.01;
	state.lambda_mem_length_u = 200;
	state.lambda_mem_length_v = 1000;

	state.initial_moves = 200;
	state.tau = 100;
	state.freeze_count = 100;

	state.update_S_skip = 1;
	state.control = 1;
	state.criterion = 0.01;
#ifdef MPI
	state.mix_interval = 10;
#endif
	state.gain_for_jump_size_control = 5;
	state.interval = 100;

	dist_params.distribution = 1;
	dist_params.q = 1;

	state.scoreFunction = NULL;
	state.printFunction = NULL;

	/* create time file name by appending .times to input file name */
	logs.dir = (char *)calloc(MAX_RECORD, sizeof(char));
	logs.dir = strcpy(logs.dir, "logs");
	logs.trace_score = 0;
	logs.trace_params = 0;
	logs.params = 0;
	logs.res = 0;
	logs.score = 0;
	logs.pid = 0;
	logs.best_score = 0;
	logs.best_res = 0;

	state.logs = &logs;
	state.options = &options;
	state.dist_params = &dist_params;
}

void PrintMyPid()
{
	FILE * t_file = fopen("pid","w");
	fprintf(t_file, "%d", getpid());
	fclose(t_file);
}


PArrPtr * InitPLSAParameters(int nb_dimensions)
{
	// plsa_params = (PArrPtr *) malloc(sizeof(PArrPtr));

	ParamList *p = (ParamList *) malloc(nb_dimensions * sizeof(ParamList));

	plsa_params.size  = nb_dimensions;
	plsa_params.array = p;

	return &plsa_params;
}


SAType * InitPLSA(int * nb_procs, int * my_id)
{

#ifdef MPI
	// MPI initialization steps

	int rc = MPI_Init(NULL, NULL); 	     /* initializes the MPI environment */
	if (rc != MPI_SUCCESS)
		printf (" > Error starting MPI program. \n");

	MPI_Comm_size(MPI_COMM_WORLD, nb_procs);        /* number of processors? */
	MPI_Comm_rank(MPI_COMM_WORLD, my_id);         /* ID of local processor? */

	nnodes = *nb_procs;
	myid = *my_id;
#endif

	SetDefaultOptions();


	InitializePLSA();

	return &state;
}



/*** THE FINAL MOVE FUNCTION ***********************************************/

/*** FinalMove: reads final energy and move count, then prints $version, ***
 *              $annealing_output and $eqparms sections and removes the    *
 *              state file                                                 *
 ***************************************************************************/

double FinalMove()
{
#ifdef MPI
	int      i;                                              /* loop counter */
#endif

	AParms   ap;
	double 	 final_score = DBL_MAX;
#ifdef MPI
	int      winner = 0;                           /* id of the winning node */
	double   minyet = DBL_MAX;         /* minimum score, used to find winner */

	double   *final_e;               /* array of final energies of all nodes */

	final_e = (double *)calloc(nnodes, sizeof(double));
#endif

	/* get the answer and some additional info */

	ap   = GetFinalInfo();              /* reads final energy and move count */


#ifdef MPI

	/* parallel code: find the node with the lowest score */

	for(i=0; i<nnodes; i++)                   /* initialize the energy array */
		final_e[i] = 0;
	/* collect the final scores from all nodes */
	MPI_Allgather(&ap.stop_energy, 1, MPI_DOUBLE, final_e, 1, MPI_DOUBLE,
				  MPI_COMM_WORLD);

	for(i=0; i<nnodes; i++)       /* evaluate the node with the lowest score */
	{
		if ( final_e[i] <= minyet )
		{
			minyet = final_e[i];
			winner = i;
		}
	}

	final_score = minyet;
#else
	final_score = ap.stop_energy;
#endif

#ifdef MPI

	/* First we share the final parameters values */
	double * res_params = malloc(sizeof(double)*plsa_params.size);
	if ( myid == winner )
	{
		int ii;
		for (ii=0; ii < plsa_params.size; ii++)
			res_params[ii] = *(plsa_params.array[ii].param);
	}

	MPI_Bcast(res_params, plsa_params.size, MPI_DOUBLE, winner, MPI_COMM_WORLD);

	if (myid != winner)
	{
		int ii;
		for (ii=0; ii < plsa_params.size; ii++)
			*(plsa_params.array[ii].param) = res_params[ii];

	}


	/* Then we write result logs */
	if (myid == winner)
	{
#endif

		if (logScore() > 0)
		{
			char score_final_name[MAX_RECORD];
			sprintf(score_final_name, "%s/score/score", getLogDir());

			FILE * score_final_f = fopen(score_final_name,"w");
			fprintf(score_final_f,"%g",ap.stop_energy);
			fclose(score_final_f);
		}

		if (logParams()> 0)
		{
			char parameters_output[MAX_RECORD];
			sprintf(parameters_output, "%s/params/output", getLogDir());

			FILE * parameters = fopen(parameters_output,"w");

			int ii;
			for (ii=0; ii < plsa_params.size; ii++)
				fprintf(parameters, "%s : %.16g\n",
							plsa_params.array[ii].name,
							*(plsa_params.array[ii].param));

			fclose(parameters);
		}

		if (logRes() > 0)
		{
			char res_dir[MAX_RECORD];
			sprintf(res_dir, "%s/res/", getLogDir());

			state.printFunction(res_dir, 0);
		}

#ifdef MPI
	}
#endif

	/* clean up the state file and free memory */

#ifdef MPI
	if ( ! options.tuning )
#endif
	StateRm();

	return final_score;
}





/*** WriteTimes: writes the time-structure to a .times file ****************
 ***************************************************************************/

void WriteTimes(double *times)
{
	char   *timefile;                             /* name of the .times file */
	FILE   *timeptr;                         /* file pointer for .times file */

	/* create time file name by appending .times to input file name */
	timefile = (char *)calloc(MAX_RECORD, sizeof(char));
	timefile = strcpy(timefile, "plsa.times");

	timeptr = fopen(timefile, "w");
	if ( !timeptr )
		file_error("main");

	PrintTimes(timeptr, times);                /* write times to .times file */

	fclose(timeptr);                                             /* clean up */
	free(timefile);
}



/*** PrintTimes: writes two (parallel: three) times sections ***************
 ***************************************************************************/

void PrintTimes(FILE *fp, double *times)
{
	fprintf(fp, "wallclock: %.3f\n", times[0]);
	fprintf(fp, "user:      %.3f\n", times[1]);
}




/*** FUNCTIONS THAT INITIALIZE/RESTORE THINGS ******************************/


/*** InitialMove: initializes the following stuff: *************************
 *                - various static things for filenames and such           *
 *                - Lam parameter struct for use in lsa.c (from tune sect) *
 *                - model and scoring funcs & solver stepsize              *
 *                - move generation in moves.c                             *
 *                - sets initial energy and checks validity of initial     *
 *                  parameters according to limit ranges                   *
 *                then it returns the initial temperature to the caller    *
 ***************************************************************************/

void InitialMove(double *p_chisq)
{

	char    *p;
	// double  i_temp;
	double  energy;

	/* initialize some Lam/Greening structures */
	// plsa_params = params;


	p = state.progname;     /* tune.progname contains program name */
	p = strcpy(p, version);

	state.debuglevel = 0;          /* following stuff not used now */
	p = state.tunename;
	p = strcpy(p, "The Other One");        /* Grateful Dead tune, what else? */

	InitScoring(&state);                   /* initializes facts and limits */
	InitMoves(&state, &plsa_params);     /* set initial temperature and *
											   *  initialize                 */
	// InitDistribution();   /* initialize distribution stuff */

	energy = -999;

	*p_chisq = energy;                                  /* set initial score */

	return;//(i_temp);                                   /* initial temperature */
}



/*** RestoreState: called when an interrupted run is restored; does the ****
 *                 following:                                              *
 *                 - stores various static things for filenames and such   *
 *                 - initializes Lam parameters for lsa.c                  *
 *                 - initializes model and scoring funcs & solver stepsize *
 *                 - initializes move generation in moves.c                *
 *                 - restores state at which previous run was interrupted  *
 *                   according to state file                               *
 *                                                                         *
 * Comment by JR:  RestoreState was originally going to be implemented with*
 * branches in InitialMove. That won't work because when when this func.   *
 * returns, control should go right to Loop(), skipping all the initiali-  *
 * zation stuff in Initialize. Hence most of the code in InitialMove is    *
 * just repeated here.                                                     *
 ***************************************************************************/

void RestoreState(char *statefile, double *p_chisq)
{
	char           *p;                                   /* temporary string */

	// Opts           *options;         /* used to restore command line options */
	MoveState      *move_ptr;                       /* used to restore moves */
	double         *stats;                      /* used to restore Lam stats */
	unsigned short *rand;                         /* used to restore ERand48 */
	double         delta[2];                        /* used to restore times */

	/* allocate memory for structures that will be returned (stats gets allo-  *
	 * cated in StateRead(), since we need to know if we're tuning or not)     */

	// options = (Opts *)malloc(sizeof(Opts));
	// options->inname    = (char *)calloc(MAX_RECORD, sizeof(char));
	// options->outname   = (char *)calloc(MAX_RECORD, sizeof(char));

	stats    = (double *)calloc(31, sizeof(double));
	move_ptr = (MoveState *)malloc(sizeof(MoveState));
	rand     = (unsigned short *)calloc(3, sizeof(unsigned short));

	StateRead(statefile, &options, move_ptr, stats, rand, delta);

	/* restore options in plsa.c (and some in lsa.c) */

	RestoreOptions(&options);

	/* initialize some Lam/Greening structures */

	p = state.progname;     /* tune.progname contains program name */
	p = strcpy(p, version);

	state.debuglevel = 0;          /* following stuff not used now */
	p = state.tunename;
	p = strcpy(p,"The Other One");         /* Grateful Dead tune, what else? */

	/* init initial cond., mutator and deriv */
	InitScoring(&state);                            /* init facts and limits */
	InitMoves(&state, &plsa_params);   /* set initial temperature and initialize */
	// InitDistribution();   /* initialize distribution stuff */

	RestoreMoves(move_ptr);
	energy = RestoreLamstats(stats);
	if ( options.time_flag )
		RestoreTimes(delta);
	InitERand(rand);

}



void StartPLSA()
{
	/* first get Lam parameters, initial temp and energy and initialize S_0 */
	/* if we restore a run from a state file: call RestoreState() */
	double S_0;

#ifdef MPI
	if (myid == 0)
	{
#endif
	if (logParams() > 0)
	{
		char parameters_input[MAX_RECORD];
		sprintf(parameters_input, "%s/params/input", getLogDir());
		FILE * parameters = fopen(parameters_input,"w");

		int ii;
		for (ii=0; ii < plsa_params.size; ii++)
		{
			fprintf(parameters, "%s : %.16g\n",
						plsa_params.array[ii].name,
						*(plsa_params.array[ii].param));
		}
		fclose(parameters);
	}
#ifdef MPI
	}
#endif


	if ( !stateflag )
	{
		InitialMove(&energy);
		S_0 = 1./state.initial_temp;
	}
	else
	   RestoreState(statefile, &energy);

	/* initialize those static file names that depend on the output file name */

	InitFilenames();

#ifdef MPI
	/* note that for parallel code both tau and init must be divided by nnodes */
	/* and we need to account for the case when tau isn't divisible by nnodes  */

	if ( (state.tau % nnodes) != 0 )
		error("plsa: the number of processors (%d) must be a divisor of tau",
			  nnodes);

	if ( (state.initial_moves % nnodes) != 0 )
		error("plsa: number of init moves must be divisible by nnodes (%d)",
			  nnodes);

#endif

	/* if we're not restarting: do the initial moves for randomizing and ga-   *
	 * thering initial statistics                                              */

	if ( !stateflag )
		energy = InitialLoop(&state, S_0);

	/* write first .log entry and write first statefile right after init; */


	if ( !stateflag )
	{
		WriteLog(&state);
#ifdef MPI
		if ( !options.tuning )
#endif
		StateWrite(statefile, energy);
	}
	else
		RestoreLog(&state);



#ifdef MPI
	/* if we are in tuning mode: initialize/restore tuning structs */
	if ( options.tuning )
		InitTuning(state.mix_interval, (double) state.tau);
#endif

}



PLSARes * runPLSA()
{
	double *delta;                            /* used to store elapsed times */
	double final_score;
	/* code for timing: wallclock and user times */

	// printf("calling initializelogs : %s\n", logs.dir);

	InitializeLogs(&logs);							/* build the log folders */

	if (logPid() > 0)
		PrintMyPid();

	cpu_start  = (struct tms *)malloc(sizeof(struct tms));      /* user time */
	cpu_finish = (struct tms *)malloc(sizeof(struct tms));

	times(cpu_start);

#ifdef MPI
	start = MPI_Wtime();       /* returns wallclock time on the calling node */
#else
	start = time(NULL);     /* returns wallclock time since EPOCH (1/1/1970) */
#endif

	/* initialize cost function and move state, do initial moves (or restore   */
	/* annealing state if restart                                              */

	StartPLSA();

	if (!Loop(&state, energy, statefile, stop_flag))
		final_score = FinalMove();

	/* code for timing */

	if ( options.time_flag )
	{
		delta = GetTimes();                  /* calculates times to be printed */
#ifdef MPI
		if ( myid == 0 )
#endif
			WriteTimes(delta);                            /* then write them out */
		free(delta);
	}

	PLSARes * res = (PLSARes *) malloc(sizeof(PLSARes));
	res->params = malloc(sizeof(double)*plsa_params.size);

	int ii;
	for (ii=0; ii < plsa_params.size; ii++)
		res->params[ii] = *plsa_params.array[ii].param;
	res->flag = 0;
	res->score = final_score;


	free(cpu_start);
	free(cpu_finish);

#ifdef MPI
	// terminates MPI execution environment
	MPI_Finalize();
#endif

	return res;

}



/*** FUNCTIONS THAT COMMUNICATE WITH SAVESTATE.C ***************************/

/*** GetOptions: returns command line options to state.c ***************
 *               for the detailed meaning of all these options see Parse-  *
 *               CommandLine() above); Opts struct defined in moves.h      *
 ***************************************************************************/

Opts *GetOptions(void)
{
// 	Opts       *options;
//
// 	options = (Opts *)malloc(sizeof(Opts));
// 	options->stop_flag   = stop_flag;
// 	options->log_flag    = log_flag;
// 	options->time_flag   = time_flag;
// 	options->state_write = state_write;
// 	options->print_freq  = print_freq;
// 	options->captions    = captions;
// 	// options->precision   = precision;
// 	options->quenchit    = quenchit;
//
// #ifdef MPI
// 	options->covar_index     = covar_index;
// 	options->write_tune_stat = write_tune_stat;
// 	options->auto_stop_tune  = auto_stop_tune;
// #endif

	return &options;
}



/*** RestoreOptions: restores the values of the command line opt variables *
 *                   from the Opts struct (used for restoring a run)       *
 ***************************************************************************/

void RestoreOptions(Opts *opts)
{
//
// 	/* all the other options */
// 	stop_flag   = options->stop_flag;
// 	log_flag    = options->log_flag;
// 	time_flag   = options->time_flag;
// 	state_write = options->state_write;
// 	print_freq  = options->print_freq;
// 	captions    = options->captions;
// 	// precision   = options->precision;
// 	quenchit    = options->quenchit;
//
// #ifdef MPI
// 	covar_index     = options->covar_index;
// 	write_tune_stat = options->write_tune_stat;
// 	auto_stop_tune  = options->auto_stop_tune;
// #endif
//
// 	free(options);
	options = *opts;
}
