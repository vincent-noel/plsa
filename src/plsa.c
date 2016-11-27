/*****************************************************************************
 *                                                                           *
 *   plsa.c                                                                *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   written by JR, modified by Yoginho                                      *
 *   modified by Vincent Noel                                                *
 *   -D landscape option by Lorraine Greenwald Oct 2002                      *
 *   -g option by Yousong Wang, Feb 2002                                     *
 *   -a option by Marcel Wolf, Apr 2002                                      *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Although main() is in lsa.c, this is the file that 'defines'            *
 *   the plsa program, since it contains most of its problem-              *
 *   specific code (except for move generation -> moves.c, saving            *
 *   of intermediate state files -> savestate.c and communication            *
 *   with the specific cost function that is used -> translate.c).           *
 *                                                                           *
 *   After I've told you all that's NOT in this file, here's what            *
 *   the funcs below actually do: parsing plsa command line opts           *
 *   is one of its jobs; there are funcs that make the first and             *
 *   last moves and funcs that read and write Lam and Lam-indepen-           *
 *   dent annealing parameters to the problem-specific data file.            *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (C) 2016 Vincent Noel                                         *
 *   the full GPL copyright notice can be found in lsa.c                     *
 *                                                                           *
 *****************************************************************************/

#include <float.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

#include "error.h"                                 /* error handling funcs */
#include "distributions.h"               /* DistP.variables and prototypes */
#include "config.h"                          /* for olddivstyle and such */
#include "moves.h"                     /* problem-specific annealing funcs */
#include "random.h"                                     /* for InitERand() */
#include "sa.h"                     /* problem-independent annealing funcs */
#include "score.h"                             /* for init and Score funcs */

#ifdef MPI                 /* this inludes parallel-specific stuff for MPI */
#include <mpi.h>                     /* this is the official MPI interface */
#include "MPI.h"  /* our own structs and such only needed by parallel code */
#endif


/*** STATIC VARIABLES ******************************************************/

static char version[MAX_RECORD];                 /* version gets set below */

/* other static variables */

// static int    precision   = 16;                    /* precision for eqparms */


/*** FUNCTIONS *************************************************************/

/*** COMMAND LINE OPTS ARE ALL PARSED HERE *********************************/

/*** ParseCommandLine: well, parses the command line and returns an index **
 *                     to the 1st argument after the command line options  *
 ***************************************************************************/
#define VERS 0.1

void ParseCommandLine()
{
    /* set the version string */

#ifdef MPI
    sprintf(version, "plsa version %f parallel", VERS);
#else
    sprintf(version, "plsa version %f serial", VERS);
#endif

    /* following part sets default values for command line options */

    captions        = 100000000;  /* default freq for writing captions (off) */
    print_freq      = 1;           /* default freq for writing to log file */
    state_write     = 100;        /* default freq for writing to state file */

    stop_flag       = absolute_freeze;             /* type of stop criterion */
    time_flag       = 0;                                  /* flag for timing */
    log_flag        = 0;              /* flag for writing logs to the screen */
    nofile_flag     = 0;       /* flog for not writing .log and .state files */
    max_iter        = 0;
    max_seconds     = 0;
    struct timeval tp;
    gettimeofday(&tp, NULL);
    start_time_seconds = (int) tp.tv_sec;
#ifdef MPI
    tuning          = 0;                    /* tuning mode is off by default */
    covar_index     = 1;      /* covariance sample will be covar_index * tau */
    write_tune_stat = 1;         /* how many times do we write tuning stats? */
    auto_stop_tune  = 1;               /* auto stop tuning runs? default: on */
    write_llog      = 0; /* write local llog files when tuning; default: off */
#endif


}




/*** FUNCTIONS THAT INITIALIZE/RESTORE THINGS ******************************/

/**********************************************************************
 * This routine reads the distribution parameters  from the input file*
 * and stores them in DistP.xx from distributions.h                   *
 * LG 03-02: need q for gen visiting distribution input file          *
 * LG 05-02: set factors only dependent on q for general visiting     *
 * distribution by calling qgt2_init or qlt2_init from distributions.c*
 **********************************************************************/

void InitDistribution(SAType * tune)
{

    DistP.distribution = tune->distribution;
    DistP.q = tune->q;


    if ((DistP.distribution > 11) || (DistP.distribution < 1))
    {
        error ("plsa: distribution must be int from 1 to 11 \n");
    }
    else if  ((DistP.distribution == 4)||(DistP.distribution == 3))
    {
        error ("plsa: PLEASE use 5 for Lorentz or 10 for normal distribution \n");
    }
    else if ((DistP.distribution == 6)||(DistP.distribution == 9))
    {
        error ("plsa: 6=poisson or 9=pareto distribution returns positive values--do not use for fly \n");
    }
    else if (DistP.distribution == 7 )
    {
        /* general distribution */
        if ((DistP.q >= 3.0) || (DistP.q <= 1.0))
        {
            error ("tsp_sa: q must be between 1 and 3 \n");
        }
        else if (DistP.q == 2.0)
        {
            DistP.distribution = 5;
            /* fly needs lorentz, tsp use abs lorentz(4)*/
            printf ("plsa: q=2 is lorentz--setting distribution to 5\n");
        }
        else if (DistP.q > 2.0)
        {
            qgt2_init();
        }
        else
        {
            qlt2_init();
        }
    } /* end of general distribution */


}  /* end InitDistribution */

/*** InitialMove: initializes the following stuff: *************************
 *                - various static things for filenames and such           *
 *                - Lam parameter struct for use in lsa.c (from tune sect) *
 *                - model and scoring funcs & solver stepsize              *
 *                - move generation in moves.c                             *
 *                - sets initial energy and checks validity of initial     *
 *                  parameters according to limit ranges                   *
 *                then it returns the initial temperature to the caller    *
 ***************************************************************************/

void InitialMove(SAType * state, double *p_chisq, PArrPtr * params)
{

    char    *p;

    // double  i_temp;
    double  energy;

    /* initialize some Lam/Greening structures */

	p = state->progname;     /* tune.progname contains program name */
    p = strcpy(p, version);

    state->debuglevel = 0;          /* following stuff not used now */
    p = state->tunename;
    p = strcpy(p, "The Other One");        /* Grateful Dead tune, what else? */
    state->tunefile = NULL;

    InitScoring(state);                   /* initializes facts and limits */
    InitMoves(state, params);     /* set initial temperature and *
                                               *  initialize                 */
    InitDistribution(state);   /* initialize distribution stuff */

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

void RestoreState(char *statefile, SAType * state, double *p_chisq,
                    PArrPtr * params)
{
    char           *p;                                   /* temporary string */

    Opts           *options;         /* used to restore command line options */
    MoveState      *move_ptr;                       /* used to restore moves */
    double         *stats;                      /* used to restore Lam stats */
    unsigned short *rand;                         /* used to restore ERand48 */
    double         delta[2];                        /* used to restore times */

    /* allocate memory for structures that will be returned (stats gets allo-  *
     * cated in StateRead(), since we need to know if we're tuning or not)     */

    options = (Opts *)malloc(sizeof(Opts));
    // options->inname    = (char *)calloc(MAX_RECORD, sizeof(char));
    // options->outname   = (char *)calloc(MAX_RECORD, sizeof(char));

    stats    = (double *)calloc(31, sizeof(double));
    move_ptr = (MoveState *)malloc(sizeof(MoveState));
    rand     = (unsigned short *)calloc(3, sizeof(unsigned short));

    StateRead(statefile, options, move_ptr, stats, rand, delta);

    /* restore options in plsa.c (and some in lsa.c) */

    RestoreOptions(options);

    /* initialize some Lam/Greening structures */

	p = state->progname;     /* tune.progname contains program name */
	p = strcpy(p, version);

	state->debuglevel = 0;          /* following stuff not used now */
	p = state->tunename;
	p = strcpy(p,"The Other One");         /* Grateful Dead tune, what else? */
	state->tunefile = NULL;

    /* init initial cond., mutator and deriv */
    InitScoring(state);                            /* init facts and limits */
    InitMoves(state, params);   /* set initial temperature and initialize */
    InitDistribution(state);   /* initialize distribution stuff */

    RestoreMoves(move_ptr);
    RestoreLamstats(stats);
    if ( time_flag )
        RestoreTimes(delta);
    InitERand(rand);

}





/*** THE FINAL MOVE FUNCTION ***********************************************/

/*** FinalMove: reads final energy and move count, then prints $version, ***
 *              $annealing_output and $eqparms sections and removes the    *
 *              state file                                                 *
 ***************************************************************************/

double FinalMove(void)
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
    /* write the answer */

    if ( myid == winner )
    {
#endif

        /* all the funcs below write a section at its appropriate position in the  */
        /* data file; to achieve this, they create a temporary file which is then  */
        /* renamed to the final output file name                                   */

        FILE * score_final_f;
        char score_final_name[MAX_RECORD];
        sprintf(score_final_name,"%s/../final_score", getLogDir());
        score_final_f = fopen(score_final_name,"w");
        fprintf(score_final_f,"%g",ap.stop_energy);
        fclose(score_final_f);



#ifdef MPI
    }
#endif

    /* clean up the state file and free memory */

    if ( !equil && !nofile_flag )
#ifdef MPI
        if ( ! tuning )
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
    // timefile = strcpy(timefile, outname);
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





/*** FUNCTIONS THAT COMMUNICATE WITH SAVESTATE.C ***************************/

/*** GetOptions: returns command line options to savestate.c ***************
 *               for the detailed meaning of all these options see Parse-  *
 *               CommandLine() above); Opts struct defined in moves.h      *
 ***************************************************************************/

Opts *GetOptions(void)
{
    Opts       *options;

    options = (Opts *)malloc(sizeof(Opts));
    options->stop_flag   = stop_flag;
    options->log_flag    = log_flag;
    options->time_flag   = time_flag;
    options->state_write = state_write;
    options->print_freq  = print_freq;
    options->captions    = captions;
    // options->precision   = precision;
    options->quenchit    = quenchit;

#ifdef MPI
    options->covar_index     = covar_index;
    options->write_tune_stat = write_tune_stat;
    options->auto_stop_tune  = auto_stop_tune;
#endif

    return options;
}



/*** RestoreOptions: restores the values of the command line opt variables *
 *                   from the Opts struct (used for restoring a run)       *
 ***************************************************************************/

void RestoreOptions(Opts *options)
{

    /* all the other options */
    stop_flag   = options->stop_flag;
    log_flag    = options->log_flag;
    time_flag   = options->time_flag;
    state_write = options->state_write;
    print_freq  = options->print_freq;
    captions    = options->captions;
    // precision   = options->precision;
    quenchit    = options->quenchit;

#ifdef MPI
    covar_index     = options->covar_index;
    write_tune_stat = options->write_tune_stat;
    auto_stop_tune  = options->auto_stop_tune;
#endif

    free(options);
}
