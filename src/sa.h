/*****************************************************************************
 *                                                                           *
 *   sa.h                                                                    *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   written by Jimmy Lam and Dan Greening                                   *
 *   modified by John Reinitz and Johannes Jaeger                            *
 *   modified by Vincent Noel                                                *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   see ../doc/prog_man.ps for details of how this works                    *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   IMPORTANT: IF YOU EVER CHANGE ANYTHING IN THIS FILE, LET ALL            *
 *            YOUR FELLOW PROGRAMMERS KNOW WELL IN ADVANCE AND               *
 *            CONSULT WITH THEM IF THEY AGREE ON YOUR CHANGES!!              *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   The following came with lsa.c, obtained from Dan Greening,              *
 *   who got it from Jimmy Lam. It is probably copyright Jimmy Lam           *
 *   1988, with modifications by Greening. We (JR & JJ) have                 *
 *   mostly modified it by adding some comments and restructuring            *
 *   the code to make it more easily understandable. JR has added            *
 *   the criterion for continous search spaces and separated Lam             *
 *   statistics from move acceptance statistics. King-Wai Chu has            *
 *   written the equilibration code.                                         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   NOTE: this header only contains prototypes for functions used           *
 *       for serial annealing code; all prototypes that are spe-             *
 *       cific to parallel code need to go into MPI.h                        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (C) 2016 Vincent Noel                                         *
 *   the full GPL copyright notice can be found in lsa.c                     *
 *                                                                           *
 *****************************************************************************/

#ifndef SA_INCLUDED
#define SA_INCLUDED

/* this def needed for func. defs that refer to (* FILE) */
#ifndef _STDIO_INCLUDED
#include <stdio.h>
#endif

/* following for structures & consts used thruout */
#ifndef GLOBAL_INCLUDED
#include "global.h"
#endif




/*** CONSTANTS *************************************************************/

#define MIN_DELTA    -100.    /* minimum exponent for Metropolis criterion */
                    /* provides a minimum probability for really bad moves */




/*** STRUCTS AND ENUMS *****************************************************/

/* These are the Lam parameters: *******************************************
 * (see King-Wai Chu's thesis for a detailed explanation)                  *
 *                                                                         *
 * lambda:              overall annealing accuracy                         *
 * lambda_mem_length_u: weighting factor for estimating mean energy        *
 * lambda_mem_length_v: weighting factor for estimating standard deviation *
 * initial_moves:       number of initial moves to randomize initial state *
 *                      of the system and gather initial statistics        *
 * tau:                 number of moves between updating statistics and    *
 *                      Lam parameters                                     *
 * freeze_count:        number of times system needs to be 'Frozen()' be-  *
 *                      fore we stop annealing                             *
 * update_S_skip:       number of iterations between changing the inverse  *
 *                      temperature; exists for increasing code efficiency *
 *                      and is probably obsolete now                       *
 * control:             OBSOLETE, used to be the control parameter for the *
 *                      move generation; acceptance statistics now live in *
 *                      move(s).c, since they are problem-specific         *
 * criterion:           defines the limits of what counts as 'Frozen()';   *
 *                      depends on the stop_flag (see below), i.e. it is   *
 *                      either a fixed energy (for discrete problems), a   *
 *                      limit to the change in absolute energy or the lim- *
 *                      it of proportional change in energy (good for con- *
 *                      tinous search spaces)                              *
 * mix_interval:        the mixing interval for parallel code in 'tau'     *
 *                      units, i.e. 100 means mix every 100 tau            *
 *                                                                         *
 ***************************************************************************
 *                                                                         *
 * NOTE: the stopping criterion is actually defined by the combination of  *
 *       'freeze_count' and 'criterion' below plus the stop_flag which     *
 *       defines the type of the criterion (absolute or proportional);     *
 *       usually when we talk about criterion (e.g. Chu, 2001) we mean     *
 *       only the value of criterion (kappa) though, since 'freeze_count'  *
 *       and the stop_flag should always be constant for specific problems *
 *                                                                         *
 ***************************************************************************/

typedef struct
{
	long   seed;

	double initial_temp;                 /* initial temperature for annealer */

	double lambda;
	double lambda_mem_length_u;
	double lambda_mem_length_v;
	int    initial_moves;
	int    tau;
	int    freeze_count;
	int    update_S_skip;
	double control;
	double criterion;
#ifdef MPI
	int    mix_interval;
#endif

/* These were marked "Application program must set these." in the code     */
/* from Greening/Lam; only progname is used at the moment; we kept them    */
/* mainly for historical reasons (and to write funny things into tunename) */

	char   progname[128];                                 /* name of my prog */
	int    debuglevel;                                          /* who cares */
	char   tunename[128];                                         /* nothing */
	FILE   *tunefile;                                                /* NULL */


	double      gain_for_jump_size_control;
	double      interval;



	int         distribution;        /* move generation distribution type RO */
							/* 1 - exp; 2 - uni; 3 - absnor; 4 - abs lorentz */
							  /*  LG: 07-05-00 formerly dist_type in lj code */

	double      q;                 /* gen visiting distribution parameter RO */
											 /* 1=guassian; 2=lor; but 1<q<3 */
						  /*  LG: 03-02 need q and factors for GSA visit dist*/
						   /*  1   <  q < 2   uses qlt2_visit       */
						   /*  2   <  q < 2.6 uses qgt2_visit       */
						   /*  2.6 <= q < 3   uses binom_qgt2_visit */

	double      (*scoreFunction)   ();
	void        (*printFunction)   (char * path, int proc);

} SAType;


/* Flag for type of stopping criterion (added by JR) ***********************
 *                                                                         *
 * 'criterion' sets the limit within which the following must lie for      *
 * 'freeze_count' number of times:                                         *
 *                                                                         *
 * proportional freeze: (mean - old_mean)/mean                             *
 * absolute freeze:      mean - old_mean                                   *
 * absolute energy:      mean                                              *
 *                                                                         *
 ***************************************************************************/

typedef enum StopStyle {
  proportional_freeze,
  absolute_freeze,
  absolute_energy
} StopStyle;


typedef struct Range {
  double   lower;
  double   upper;
} Range;

typedef struct {
  double    *param;                /* pointers to parameters to be tweaked */
  Range     param_range;        /* pointers to corresponding range limits */
} ParamList;


typedef struct {
  int       size;                           /* size of the ParamList array */
  ParamList *array;            /* points to 1st element of ParamList array */
} PArrPtr;




/*** GLOBALS ***************************************************************/

// NucStatePtr state;                    /* global annealing parameter struct */
SAType		state;

StopStyle   stop_flag;               /* type of stop criterion (see above) */

int         time_flag;                         /* flag for timing the code */
int         log_flag;                 /* flag for displaying log to stdout */
int         nofile_flag;      /* flag for not writing .state or .log files */

long        state_write;     /* frequency for writing state files (in tau) */
long        print_freq;          /* frequency for printing to log (in tau) */
long        captions;      /* option for printing freqency of log captions */

/* Special annealing modes: ************************************************
 *                                                                         *
 * the following are flags that should be set by cmd line opts for some    *
 * special annealing modes that we use for testing or gathering stats:     *
 *                                                                         *
 * - bench:    does no loop, i.e. only runs through initial steps; this    *
 *             can be useful to do solver benchmarks in -B mode            *
 * - equil:    runs at a fixed temperature to collect equilibrium stats    *
 * - quenchit: means lower the temperature to zero *instantly* after fini- *
 *             shing initialize; this is not implemented for parallel code *
 *             since it doesn't take a long time to finish                 *
 *                                                                         *
 ***************************************************************************/

int         bench;
int         equil;
int         quenchit;
int         max_iter;
int         max_seconds;
int         start_time_seconds;




/*** FUNCTION PROTOTYPES ***************************************************/

/* problem independent functions that live in lsa.c ************************/

/* Initializing functions */

/*** Initialize: calls ParseCommandLine first; then does either initial ****
 *               randomization and collecting Lam stats or restores state  *
 *               of the annealer as saved in the state file                *
 ***************************************************************************/
#ifdef MPI
SAType * InitPLSA(int nb_procs, int my_id);
#else
SAType * InitPLSA();
#endif

PArrPtr * InitPLSAParameters(int nb_dimensions);
double runPLSA( PArrPtr * params);

SAType * InitializePLSA(void);
void StartPLSA(PArrPtr * params);

/*** InitFilenames: initializes static file names needed in lsa.c **********
 ***************************************************************************/

void InitFilenames(void);

/*** InitialLoop: performs the two sets of initial moves: ******************
 *                   1. randomizing moves (not parallelized)               *
 *                   2. loop for initial collection of statistics          *
 ***************************************************************************/

void InitialLoop(void);

/*** InitializeParameter: initializes variables used for calculating the ***
 *                        parameters A, B, D and E and the estimate_mean   *
 *                        plus estimate_sd (see Lam & Delosme, 1988a)      *
 ***************************************************************************/

void InitializeParameter(void);

/*** InitializeWeights: initialize weights a and b; these weights are ******
 *                      computed from the lambda memory length products    *
 ***************************************************************************/

void InitializeWeights(void);



/* Main loop and update functions */

/*** Loop: loops (making moves, updating stats etc.) until the system is ***
 *         considered frozen according to the stop criterion               *
 ***************************************************************************/

double Loop(void);

/*** UpdateS: update inverse temperature S at every Sskip step *************
 ***************************************************************************/

void UpdateS(void);

/*** UpdateStats: updates mean, variance and acc_ratio after tau moves *****
 *                it needs i to do sanity check in parallel code           *
 ***************************************************************************/

void UpdateStats(int i);

/*** UpdateParameter: update parameters A, B, D and E and the estimators ***
 *                    for mean and standard deviation for the current S    *
 ***************************************************************************/

void UpdateParameter(void);

/*** Frozen: returns TRUE if frozen, FALSE otherwise ***********************
 ***************************************************************************/

int Frozen(void);


/* functions which communicate with other source files, these are needed  *
 * for reading/writing .state files                                       */

/*** GetLamstats: returns Lam statistics in an array of doubles; used to ***
 *                store Lam statistics in a state file                     *
 ***************************************************************************/

double *GetLamstats(void);

/*** GetTimes: returns a two-element array with the current wallclock and **
 *             user time to be saved in the state file                     *
 ***************************************************************************/

double *GetTimes(void);

/* functions which restore things in lsa.c upon restart from state file */

/*** RestoreLamstats: restores static Lam statistics in lsa.c from an ******
 *                    array of doubles; used to restore runs from a state  *
 *                    file.                                                *
 ***************************************************************************/

void RestoreLamstats(double *stats);

/*** RestoreLog: restores .log and .prolix files after upon restart ********
 ***************************************************************************/

void RestoreLog(void);

/*** RestoreTimes: restores the wallclock and user times if -t is used *****
 ***************************************************************************/

void RestoreTimes(double *delta);




/* functions to write the .log */


void WriteLog(void);

/*** PrintLog: actually writes stuff to the .log file **********************
 ***************************************************************************/

void PrintLog(FILE *outptr, int local_flag);





/* PROBLEM SPECIFIC FUNCTIONS THAT NEED TO BE DEFINED OUTSIDE LSA.C ********/

/* miscellaneous functions that usually live in <problem>_sa.c */

/*** ParseCommandLine: well, parses the command line and returns an index **
 *                     to the 1st argument after the command line options  *
 ***************************************************************************/

void ParseCommandLine();

/*** InitialMove: initializes the following stuff: *************************
 *                - reads in Lam and other annealing parameters (passed to *
 *                  lsa.c through the state_ptr; the first three arguments *
 *                  are used to open the right data file etc.)             *
 *                - initializes the cost function, establishes link be-    *
 *                  tween cost function and annealer and passes the init-  *
 *                  tial energy to lsa.c by p_chisq)                       *
 *                - initializes move generation in move(s).c               *
 *                - sets initial energy by evaluating cost function for    *
 *                  the first time                                         *
 ***************************************************************************/

void InitialMove(SAType * state_ptr_vs, double *p_chisq, PArrPtr * params);

/*** RestoreState: called when an interrupted run is restored; does the ****
 *                 following (see InitialMove for arguments, also see co-  *
 *                 mmunication functions below for how to restore the ran- *
 *                 dom number generator and Lam stats):                    *
 *                 - restores Lam and other annealing parameters           *
 *                 - reinitializes the cost function                       *
 *                 - restores move state in move(s).c                      *
 ***************************************************************************/

void RestoreState(char *statefile, SAType * state_ptr_vs, double *p_chisq,
	                      PArrPtr * params);

/*** FinalMove: determines the final energy and move count and then prints *
 *              those to wherever they need to be printed to; also should  *
 *              do the cleaning up, i.e freeing stuff and such after a run *
 ***************************************************************************/

double FinalMove(void);

/*** WriteTimes: writes the timing information to wherever it needs to be **
 *               written to at the end of a run                            *
 ***************************************************************************/

void WriteTimes(double *times);



/* move generation functions that are used in lsa.c (live in move(s).c) */

/* GenerateMove: evaluates the old energy, changes a parameter, then eval- *
 *               uates the new energy; returns the difference between old  *
 *               and new energy to the caller                              *
 ***************************************************************************/

double GenerateMove(void);

/*** AcceptMove: sets new energy as the old energy for the next step and ***
 *               keeps track of the number of successful moves             *
 ***************************************************************************/

void AcceptMove(void);

/*** RejectMove: simply resets the tweaked parameter to the pretweak value *
 ***************************************************************************/

void RejectMove(void);

/*** GetEnergy: returned the last computed value of the scoring function   *
 *    To avoid computing e = old_e + (new e - old e), and using directly   *
 *    e = new_e
 *    Avoid precision errors due floating number representation            *
 **************************************************************************/

double GetNewEnergy(void);
double GetOldEnergy(void);
void WriteScoreTrace(double t_energy, int acceptance);



/* a function that writes the .state file (should live in savestate.c) */

/*** StateWrite: collects Lam statistics, move state and the state of the **
 *               erand48 random number generator and writes all that into  *
 *               the state file, which can then be used to restore the run *
 *               in case it gets interrupted                               *
 ***************************************************************************/

void StateWrite();

/* THE function that starts the optimization */

/*** StateWrite: collects Lam statistics, move state and the state of the **
 *               erand48 random number generator and writes all that into  *
 *               the state file, which can then be used to restore the run *
 *               in case it gets interrupted                               *
 ***************************************************************************/


#endif
