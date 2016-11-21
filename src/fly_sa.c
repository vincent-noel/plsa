/*****************************************************************************
 *                                                                           *
 *   fly_sa.c                                                                *
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
 *   the fly_sa program, since it contains most of its problem-              *
 *   specific code (except for move generation -> moves.c, saving            *
 *   of intermediate state files -> savestate.c and communication            *
 *   with the specific cost function that is used -> translate.c).           *
 *                                                                           *
 *   After I've told you all that's NOT in this file, here's what            *
 *   the funcs below actually do: parsing fly_sa command line opts           *
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

// static char   *inname;                           /* filename of input file */
// static char   *outname;                         /* filename of output file */
static int    precision   = 10;                    /* precision for eqparms */
// static int    landscape_flag = 0;        /* generate energy landscape data */


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
    sprintf(version, "fly_sa version %f parallel", VERS);
#else
    sprintf(version, "fly_sa version %f serial", VERS);
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

void InitDistribution(plsa_parameters * t_plsa_params)
{

    DistP.distribution = t_plsa_params->distribution;
    DistP.q = t_plsa_params->q;


    if ((DistP.distribution > 11) || (DistP.distribution < 1))
    {
        error ("fly_sa: distribution must be int from 1 to 11 \n");
    }
    else if  ((DistP.distribution == 4)||(DistP.distribution == 3))
    {
        error ("fly_sa: PLEASE use 5 for Lorentz or 10 for normal distribution \n");
    }
    else if ((DistP.distribution == 6)||(DistP.distribution == 9))
    {
        error ("fly_sa: 6=poisson or 9=pareto distribution returns positive values--do not use for fly \n");
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
            printf ("fly_sa: q=2 is lorentz--setting distribution to 5\n");
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

double InitialMove(NucStatePtr state_ptr, double *p_chisq,
                    plsa_parameters * settings, PArrPtr * params)
{

    char    *p;

    double  i_temp;
    double  energy;

    // if ( landscape_flag )                      /* 1 means gen landscape files */
    //     InitLandscape(landscape_flag);                       /*lives in lsa.c */

    /* initialize some Lam/Greening structures */

    p = state_ptr->tune.progname;     /* tune.progname contains program name */
    p = strcpy(p, version);

    state_ptr->tune.debuglevel = 0;          /* following stuff not used now */
    p = state_ptr->tune.tunename;
    p = strcpy(p, "The Other One");        /* Grateful Dead tune, what else? */
    state_ptr->tune.tunefile = NULL;

    InitScoring(settings);                   /* initializes facts and limits */
    i_temp  = InitMoves(settings, params);     /* set initial temperature and *
                                               *  initialize                 */
    InitDistribution(settings);   /* initialize distribution stuff */

    /* initialize Lam parameters (see sa.h for further detail) */

    state_ptr->tune.lambda              = settings->lambda;
    state_ptr->tune.lambda_mem_length_u = settings->lambda_mem_length_u;
    state_ptr->tune.lambda_mem_length_v = settings->lambda_mem_length_v;
    state_ptr->tune.control             = settings->control;
    state_ptr->tune.initial_moves       = settings->initial_moves;
    state_ptr->tune.tau                 = settings->tau;
    state_ptr->tune.freeze_count        = settings->freeze_count;
    state_ptr->tune.update_S_skip       = settings->update_S_skip;
    state_ptr->tune.criterion           = settings->criterion;
#ifdef MPI
    state_ptr->tune.mix_interval        = settings->mix_interval;
#endif
    state_ptr->tune.log_trace           = settings->log_trace;
    state_ptr->tune.log_params          = settings->log_params;

    energy = -999;

    *p_chisq = energy;                                  /* set initial score */

    return(i_temp);                                   /* initial temperature */
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

void RestoreState(char *statefile, NucStatePtr state_ptr, double *p_chisq,
                    plsa_parameters * settings, PArrPtr * params)
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

    /* restore options in fly_sa.c (and some in lsa.c) */

    RestoreOptions(options);

    /* initialize some Lam/Greening structures */

    p = state_ptr->tune.progname;     /* tune.progname contains program name */
    p = strcpy(p, version);

    state_ptr->tune.debuglevel = 0;          /* following stuff not used now */
    p = state_ptr->tune.tunename;
    p = strcpy(p,"The Other One");         /* Grateful Dead tune, what else? */
    state_ptr->tune.tunefile = NULL;

    /* init initial cond., mutator and deriv */
    InitScoring(settings);                            /* init facts and limits */
    InitMoves(settings, params);   /* set initial temperature and initialize */
    InitDistribution(settings);   /* initialize distribution stuff */

    /* initialize Lam parameters (see sa.h for further detail) */

    state_ptr->tune.lambda              = settings->lambda;
    state_ptr->tune.lambda_mem_length_u = settings->lambda_mem_length_u;
    state_ptr->tune.lambda_mem_length_v = settings->lambda_mem_length_v;
    state_ptr->tune.control             = settings->control;
    state_ptr->tune.initial_moves       = settings->initial_moves;
    state_ptr->tune.tau                 = settings->tau;
    state_ptr->tune.freeze_count        = settings->freeze_count;
    state_ptr->tune.update_S_skip       = settings->update_S_skip;
    state_ptr->tune.criterion           = settings->criterion;
#ifdef MPI
    state_ptr->tune.mix_interval        = settings->mix_interval;
#endif

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

void FinalMove(void)
{
#ifdef MPI
    int      i;                                              /* loop counter */
#endif

    AParms   ap;

    double   equil_var[2];         /* array for results of equilibration run */

#ifdef MPI
    int      winner = 0;                           /* id of the winning node */
    double   minyet = DBL_MAX;         /* minimum score, used to find winner */

    double   *final_e;               /* array of final energies of all nodes */

    final_e = (double *)calloc(nnodes, sizeof(double));
#endif

    /* get the answer and some additional info */

    ap   = GetFinalInfo();              /* reads final energy and move count */

    if ( equil )
        GetEquil(equil_var);            /* get the final equilibration results */

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


        // ExportVariablesOptim(get_optim_parameters());

#ifdef MPI
    }
#endif

    /* clean up the state file and free memory */

    if ( !equil && !nofile_flag )
#ifdef MPI
        if ( ! tuning )
#endif
            StateRm();


}






/*** FILE FUNCTIONS FOR READING AND WRITING MISCELLANEOUS STUFF ************/

/*** InitEquilibrate: reads the equilibrate section of the data file, ******
 *                    which is needed for equilibration runs; then puts    *
 *                    the parameters into a static struct in lsa.c         *
 ***************************************************************************/

 void InitEquilibrate()
{
    ChuParam   l_equil_param;                    /* local equil_param struct */

    // fp = FindSection(fp, "equilibrate");                /* find tune section */
    // if( !fp )
    //     error("ReadEquilibrate: could not locate equilibrate section");
    //
    // fscanf(fp,"%*s\n");                         /* advance past title line 1 */
    //
    // if ( 1 != fscanf(fp, "%lf\n", &(l_equil_param.end_T)) )
    //     error("ReadEquilibrate: error reading equilibration temperature");
    //
    // fscanf(fp,"%*s\n");                         /* advance past title line 2 */
    //
    // if ( 1 != fscanf(fp, "%d\n", &(l_equil_param.fix_T_skip)) )
    //     error("ReadEquilibrate: error reading fixed temperature skip");
    //
    // fscanf(fp,"%*s\n");                         /* advance past title line 3 */
    //
    // if ( 1 != fscanf(fp, "%d\n", &(l_equil_param.fix_T_step)) )
    //     error("ReadEquilibrate: error reading fixed temperature step");

    SetEquilibrate(l_equil_param);
}


/*** WriteEquil: writes the equilibrate_variance section to the data file **
 *               right after the $equilibrate section                      *
 ***************************************************************************/

void WriteEquil(char *filename, double *equil_var)
{
    // char   *temp;                                     /* temporary file name */
    // char   *record;                         /* record to be read and written */
    // char   *record_ptr;        /* pointer used to remember record for 'free' */
    // char   *saverec;                 /* used to save following section title */
    // char   *shell_cmd;                             /* used by 'system' below */
    //
    // FILE   *outfile;                                  /* name of output file */
    // FILE   *tmpfile;                               /* name of temporary file */
    //
    //
    // temp      = (char *)calloc(MAX_RECORD, sizeof(char));
    // record    = (char *)calloc(MAX_RECORD, sizeof(char));
    // saverec   = (char *)calloc(MAX_RECORD, sizeof(char));
    // shell_cmd = (char *)calloc(MAX_RECORD, sizeof(char));
    //
    // record_ptr = record;            /* this is to remember record for 'free' */
    //
    // /* open output and temporary file */
    //
    // outfile = fopen(filename, "r");              /* open outfile for reading */
    // if ( !outfile )                              /* sorry for the confusion! */
    //     error("WriteEquil: error opening %s", filename);
    //
    // temp = strcpy(temp,"equilXXXXXX");              /* required by mkstemp() */
    // if ( mkstemp(temp) == -1 )              /* get unique name for temp file */
    //     error("WriteEquil: error creating temporary file name");
    //
    // tmpfile = fopen(temp, "w");               /* ... and open it for writing */
    // if ( !tmpfile )
    //     error("WriteEquil: error opening temporary file %s", temp);
    //
    // if ( FindSection(outfile, "equilibrate_variance") )
    // {
    //     fclose(outfile);                   /* erase section if already present */
    //     KillSection(filename, "equilibrate_variance");
    //     outfile = fopen(filename, "r");
    // }
    // rewind(outfile);
    //
    // /* the follwoing three loops look for the appropriate file position to     */
    // /* write the equilibrate_variance section                                  */
    //
    // while ( strncmp(record=fgets(record, MAX_RECORD, outfile),
    //                 "$equilibrate", 12) )
    //     fputs(record, tmpfile);
    // fputs(record, tmpfile);
    //
    // while ( strncmp(record=fgets(record, MAX_RECORD, outfile), "$$", 2) )
    //     fputs(record, tmpfile);
    // fputs(record, tmpfile);
    //
    // do
    // {
    //     record = fgets(record, MAX_RECORD, outfile);
    //     if ( !record ) break;
    // }
    // while ( strncmp(record, "$", 1) );
    //
    // fputs("\n", tmpfile);
    //
    // if ( record )
    //     saverec = strcpy(saverec, record);
    //
    // /* now we write the eqparms section into the tmpfile */
    //
    // PrintEquil(tmpfile, equil_var, "equilibrate_variance");
    //
    // fprintf(tmpfile, "\n");
    //
    // /* ... and then write all the rest, if there is any */
    //
    // if ( record )
    //     fputs(saverec, tmpfile);
    //
    // while ( (record=fgets(record, MAX_RECORD, outfile)) )
    //     fputs(record, tmpfile);
    //
    // fclose(outfile);
    // fclose(tmpfile);
    //
    // /* rename tmpfile into new file */
    //
    // sprintf(shell_cmd, "cp -f %s %s", temp, filename);
    //
    // if ( -1 == system(shell_cmd) )
    //     error("WriteEquil: error renaming temp file %s", temp);
    //
    // if ( remove(temp) )
    //     warning("WriteEquil: temp file %s could not be deleted",
    //             temp);
    //
    // /* clean up */
    //
    // free(record_ptr);
    // free(saverec);
    // free(temp);
    // free(shell_cmd);
}


//
// /*** PrintEquil: writes an 'equilibrate_variance' section with 'title' *****
//  *               to the stream specified by fp                             *
//  ***************************************************************************/
//
// void PrintEquil(FILE *fp, double *equil_var, char *title)
// {
//     fprintf(fp, "$%s\n", title);
//     fprintf(fp, "variance:\n");
//     fprintf(fp, "%g\n", equil_var[0]);
//     fprintf(fp, "equilibrate_final_energy:\n");
//     fprintf(fp, "%g\n", equil_var[1]);
//     fprintf(fp, "$$\n");
// }
//
//
//


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

    // options->inname  = inname;
    // options->outname = outname;

    options->stop_flag   = stop_flag;
    // options->landscape_flag = landscape_flag;
    options->log_flag    = log_flag;
    options->time_flag   = time_flag;
    options->state_write = state_write;
    options->print_freq  = print_freq;
    options->captions    = captions;
    options->precision   = precision;
    options->quenchit    = quenchit;
    options->equil       = equil;

#ifdef MPI
    options->tuning          = tuning;
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

    /* restore input/output file names and the full command line string; note  *
     * that the output file name needs to be communicated to lsa.c, since we   *
     * need it there for setting up log and tuning file names                  */

    // inname  = options->inname;
    // outname = options->outname;
    // SetOutname(outname);

    /* all the other options */
    stop_flag   = options->stop_flag;

    // landscape_flag = options->landscape_flag;
    // if ( landscape_flag )
    //     error("RestoreOptions: cannot restore an equilibration (lanDscape) run");

    log_flag    = options->log_flag;
    time_flag   = options->time_flag;
    state_write = options->state_write;
    print_freq  = options->print_freq;
    captions    = options->captions;
    precision   = options->precision;
    quenchit    = options->quenchit;
    equil       = options->equil;
    if ( equil )
        error("RestoreOptions: cannot restore an equilibration run");

#ifdef MPI
    tuning          = options->tuning;
    if ( tuning )
        error("RestoreOptions: cannot restore a tuning run");
    covar_index     = options->covar_index;
    write_tune_stat = options->write_tune_stat;
    auto_stop_tune  = options->auto_stop_tune;
#endif

    free(options);
}
