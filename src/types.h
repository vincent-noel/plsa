/******************************************************************************
 *                                                                            *
 *   types.h                                                                  *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   Written by Vincent Noel                                                  *
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
#ifndef PLSA_TYPES
#define PLSA_TYPES


/* The following are annealing parameters that are not specific to the Lam *
 * algorithm. In general they should be used in moves.c or plsa.c but    *
 * *not* in lsa.c. In the data file, the members of the struct labeled RO  *
 * are read from the $annealing_input section. They are used for initial   *
 * conditions of the annealer and do not change during a run. Members      *
 * labeled OUT are written to the $annealing_output section upon comple-   *
 * tion of a run.                                                          */

typedef struct
{
  long   seed;                      /* seed for random number generator RO */
  double start_tempr;          /* the initial equilibration temperature RO */
  double gain;            /* gain for proportional control of move size RO */
  double stop_energy;                /* the final energy of the answer OUT */
  int    max_count;                      /* total number of iterations OUT */
  int    interval;       /* number of sweeps between updating theta_bar RO */
/*int    distribution;    1 - uniform; 2 - exp; 3 - normal; 4 - lorentz RO */
  // int    log_params;
} AParms;

typedef enum StopStyle
{
	proportional_freeze,
	absolute_freeze,
	absolute_energy

} StopStyle;

/* Opts struct is for saving command line options in state.c */

typedef struct
{
	StopStyle stop_flag;                                   /* stop criterion */
	int       time_flag;                             /* flag for timing code */
	long      state_write;              /* frequency for writing state files */
	long      print_freq;         /* frequency for printing status to stdout */
	int       quenchit;          /* flag for quenchit mode (T=0 immediately) */

#ifdef MPI

	int       tuning;                                /* flag for tuning mode */
	int       covar_index; /* index for sample interval (=covar_index * tau) */
	int       write_tune_stat;    /* how many times to write tune statistics */
	int       auto_stop_tune;                      /* auto-stop tuning runs? */

#endif

	int max_iter;
	int max_seconds;

} Opts;



typedef struct
{
	int distribution;      /* move generation distribution type RO */
						 /* 1 - exp; 2 - uni; 3 - absnor; 4 - abs lorentz */
						 /*  LG: 07-05-00 formerly dist_type in lj code */
	double q;              /* gen visiting distribution parameter RO */
						 /* 1=guassian; 2=lor; but 1<q<3 */
						 /*  LG: 03-02 need q and factors for GSA visit dist*/
						 /*  1   <  q < 2   uses qlt2_visit       */
						 /*  2   <  q < 2.6 uses qgt2_visit       */
						 /*  2.6 <= q < 3   uses binom_qgt2_visit */


} DistParms;

typedef struct
{
	/* Added some log options here */
	char * 	dir;

	int    	trace_score;
	int 	trace_params;
	int    	params;
	int  	res;
	int 	score;
	int 	pid;
	int 	best_score;
	int 	best_res;

} SALogs;


typedef struct
{
	double   lower;
	double   upper;

} Range;


typedef struct
{
	double    *param;                /* pointers to parameters to be tweaked */
	Range     param_range;        /* pointers to corresponding range limits */
	char *	name;

} ParamList;


typedef struct
{
	int       size;                           /* size of the ParamList array */
	ParamList *array;            /* points to 1st element of ParamList array */

} PArrPtr;

typedef struct
{
	double acc_ratio;                      /* acceptance ratio for parameter */
	double theta_bar;              /* theta bar is proportional to move size */
	int    hits;       /* number of moves since last call to UpdateControl() */
	int    success;              /* number of these moves that were accepted */

} AccStats;

typedef struct
{
	int    counter;                           /* counter used by Frozen */
	double old_mean;                    /* old mean as stored by Frozen */
	double energy;                                    /* current energy */

	double mean;          /* mean energy, collected from tau last steps */
	double vari;      /* energy variance, collected from tau last steps */
	double estimate_mean;              /* Lam estimator for mean energy */
	double estimate_sd;  /* Lam estimator for energy standard deviation */

	double S;                                 /* current inverse energy */
	double dS;                      /* delta S: change in S during move */
	double S_0;                      /* the initial inverse temperature */

	double alpha;         /* the third term of the Lam schedule formula */
	double acc_ratio;    /* average acceptance ratio for all parameters */


	/* Lam stats stuff: variables for calculating the estimators  **************/

	/* mean estimator */

	double w_a;                       /* w_a is the weight for the mean */

	double usyy;         /* these parameters store intermediate results */
	double usxy;               /* for the updating formulas for A and B */
	double usxx;                       /* see Lam & Delosme, 1988b, p10 */
	double usx;
	double usy;
	double usum;

	double A;            /* A and B are the parameters for the rational */
	double B;                /* function for the estimation of the mean */

	/* sd estimator */

	double w_b;         /* w_b is the weight for the standard deviation */

	double vsyy;         /* these parameters store intermediate results */
	double vsxy;               /* for the updating formulas for D and E */
	double vsxx;                       /* see Lam & Delosme, 1988b, p10 */
	double vsx;
	double vsy;
	double vsum;

	double D;            /* D and E are the parameters for the rational */
	double E;  /* function for the estimation of the standard deviation */

	/* Lam stats stuff: variables related to tau *******************************/
	int    proc_tau;  /* proc_tau = tau                     in serial   */
	/* proc_tau = tau / (# of processors) in parallel */
	long   count_tau;  /* how many times we did tau (or proc_tau) moves */

	/* the actual number of moves for collecting initial statistics ************/

	int    proc_init;                        /* number of initial moves */


} LamState;

typedef struct
{
	int 	covar_index;
	int		write_tune_stat;
	int 	auto_stop_tune;

} TuningSettings;
/* following struct contains copies of the static variables of moves.c     *
 * together with the values of parameters undergoing annealing             */

typedef struct
{
	ParamList *pt;  /* Used during a save to point to annealed-on parameters */
	AccStats  *acc_tab_ptr;            /* points to current acceptance stats */
	double    *newval; /* points to array of annealed-on doubles for restore */
	double    old_energy;                     /* energy before the last move */

	int       nparams;                      /* # of parameters to be tweaked */
	int       index;      /* index of parameter to be tweaked during a sweep */
	int       nhits;                         /* number of moves already made */
	int       nsweeps;                         /* number of completed sweeps */

} MoveState;


typedef struct
{
	int 	flag;
	double 	score;
	double *params;

} PLSARes;


typedef struct
{
	/* PLSA, aka top level settings */
	long   				seed;

	/* LSA General settings (printing and stuff )*/
	// Opts * 				options;
	int 				print_freq;   /* default freq for writing to log file */
	int 				state_write;/* default freq for writing to state file */
	int 				time_flag;
	/* LSA settings */
	double 				initial_temp;     /* initial temperature for annealer */

	double 				lambda;
	double 				lambda_mem_length_u;
	double 				lambda_mem_length_v;
	int    				initial_moves;
	int    				tau;
	int    				freeze_count;
	int    				update_S_skip;
	double 				control;
	double 				criterion;

#ifdef MPI
/* Parallel settings */
	int    				mix_interval;
	int 				tuning;
	TuningSettings *	tuning_settings;
#endif


	int 				max_iter;
	int					max_seconds;
	int 				quenchit;
	StopStyle 			stop_flag;                          /* stop criterion */

	/* These were marked "Application program must set these." in the code    */
	/* from Greening/Lam; only progname is used at the moment; we kept them   */
	/* mainly for historical reasons (and to write funny things into tunename)*/

	char   				progname[128];                     /* name of my prog */
	int    				debuglevel;                              /* who cares */
	char   				tunename[128];                             /* nothing */


	double      		gain_for_jump_size_control;
	double      		interval;

	/* Score settings */
	double      		(*scoreFunction)   	();
	void        		(*printFunction)   	(char * path, int proc);

	/* Distribution settings */
	DistParms * 		dist_params;

	/* Logs settings */
	SALogs * 			logs;

} SAType ;

#endif
