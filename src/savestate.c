/*****************************************************************************
 *                                                              						 *
 *   savestate.c                                                						 *
 *                                                              						 *
 *****************************************************************************
 *                                                               						 *
 *   written by JR, modified by Yoginho                         						 *
 *   modified by Vincent Noel                                                *
 *                                                              						 *
 *****************************************************************************
 *                                                              						 *
 *   savestate.c contains three functions that read, write and remove the    *
 *	 state file for an annealing run. The frequency with which state are		 *
 *   saved can be chosen by the command line option -b											 *
 *   (for backup stepsize). The state file is very useful for the case			 *
 *   when long annealing runs have to be interrupted or crash for some			 *
 *   reason or another. The run can then be resumed by indicating the				 *
 *	 state file as an additional argument to fly_sa.                         *
 *                                                              						 *
 *****************************************************************************
 *                                                              						 *
 *   Copyright (C) 2016 Vincent Noel                            						 *
 *   the full GPL copyright notice can be found in lsa.c        						 *
 *                                                              						 *
 *****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>                                        /* getopt stuff */

#include "error.h"
#include "moves.h"
#include "random.h"
#include "sa.h"

#ifdef MPI
#include "MPI.h"                                              /* for myid */
#endif



/*** A STATIC VARIABLE *****************************************************/

static char           *filename;                     /* name of state file */





/*** FUNCTION DEFINITIONS **************************************************/

/*** StateRead: reads Lam statistics, move state and erand state from a ****
 *              state file and restores the annealer's state to the same   *
 *              state it was in before it got interrupted                  *
 *     CAUTION: InitMoves must be called before calling StateRead!         *
 ***************************************************************************/

void StateRead(char *statefile, Opts *options, MoveState *move_ptr,
	       double *stats, unsigned short *rand, double *delta)
{
  int            i;                                  /* local loop counter */
  FILE           *infile;                         /* pointer to state file */


/* make the state file name static to savestate.c */

  filename = (char *)calloc(MAX_RECORD, sizeof(char));
  filename = strcpy(filename, statefile);

/* open the state file and read it */

  infile = fopen(filename, "r");
  if( !infile )
    file_error("StateRead");

  // options->argv = fgets(options->argv, MAX_RECORD, infile);

  // fscanf(infile, "%s\n",    options->inname);
  // fscanf(infile, "%s\n",    options->outname);
 // fscanf(infile, "%s\n",    options->derivfunc);
//  fscanf(infile, "%s\n",    options->solver);
  fscanf(infile, "%d\n",  (int*) &(options->stop_flag));
  // fscanf(infile, "%d\n",  &(options->prolix_flag));
  // fscanf(infile, "%d\n",  &(options->landscape_flag));
  fscanf(infile, "%d\n",  &(options->log_flag));
  fscanf(infile, "%d\n",  &(options->time_flag));
  fscanf(infile, "%ld\n", &(options->state_write));
  fscanf(infile, "%ld\n", &(options->print_freq));
  fscanf(infile, "%ld\n", &(options->captions));
  fscanf(infile, "%d\n",  &(options->precision));
  // fscanf(infile, "%lg\n", &(options->stepsize));
  fscanf(infile, "%d\n",  &(options->quenchit));
  fscanf(infile, "%d\n",  &(options->equil));
#ifdef MPI
  fscanf(infile, "%d\n",  &(options->tuning));
  fscanf(infile, "%d\n",  &(options->covar_index));
  fscanf(infile, "%d\n",  &(options->write_tune_stat));
  fscanf(infile, "%d\n",  &(options->auto_stop_tune));
#endif

  if ( options->time_flag ) {
    fscanf(infile, "%lf\n", &(delta[0]));
    fscanf(infile, "%lf\n", &(delta[1]));
  }

  fscanf(infile, "%d\n",  &(move_ptr->nparams));
  fscanf(infile, "%d\n",  &(move_ptr->index));
  fscanf(infile, "%d\n",  &(move_ptr->nhits));
  fscanf(infile, "%d\n",  &(move_ptr->nsweeps));

  move_ptr->newval      =
    (double *)calloc(move_ptr->nparams, sizeof(double));
  move_ptr->pt          = NULL;

  move_ptr->acc_tab_ptr =
    (AccStats *)calloc(move_ptr->nparams, sizeof(AccStats));

  for(i=0; i < move_ptr->nparams; i++)
    fscanf(infile, "%lg\n", &(move_ptr->newval[i]));
  fscanf(infile,"%lg\n", &(move_ptr->old_energy));

  for(i=0; i < move_ptr->nparams; i++)
    fscanf( infile, "%lg %lg %d %d\n",
	   &(move_ptr->acc_tab_ptr[i].acc_ratio),
	   &(move_ptr->acc_tab_ptr[i].theta_bar),
	   &(move_ptr->acc_tab_ptr[i].hits),
	   &(move_ptr->acc_tab_ptr[i].success) );

  for(i=0; i<31; i++)
    fscanf(infile, "%lg\n", &(stats[i]));

  for(i=0; i<3; i++)
    fscanf(infile, "%hu\n", &(rand[i]));

  fclose(infile);

}



/*** StateWrite: collects command line options and arguemnts, Lam statis- **
 *               tics, move state and the state of the erand48 random num- *
 *               ber generator and writes all that into the state file,    *
 *               which can then be used to restore the run in case it gets *
 *               interrupted                                               *
 ***************************************************************************/

void StateWrite(char *statefile)
{
  int i;                                             /* local loop counter */
  FILE           *outfile;                           /* state file pointer */
  Opts           *options;                /* command line opts to be saved */
  MoveState      *move_status;                    /* MoveState to be saved */
  double         *lamsave;                        /* Lam stats to be saved */
  unsigned short *prand;                    /* erand48() state to be saved */
  double         *delta;            /* wallclock and user time to be saved */


/* if StateWrite() called for the first time: make filename static */

  if (filename == NULL) {
    filename = (char *)calloc(MAX_RECORD, sizeof(char));
    filename = strcpy(filename, statefile);
  }

/* collect the state and the options */

  options     = GetOptions();
  move_status = MoveSave();
  lamsave     = GetLamstats();
  prand       = GetERandState();
  if ( time_flag )
    delta     = GetTimes();

/* write the answer; now *fully* portable, no binary!!! */

  outfile = fopen(filename, "w");
  if ( !outfile )
    file_error("StateWrite");

  // fprintf(outfile, "%s\n",    options->inname);
  // fprintf(outfile, "%s\n",    options->outname);
  fprintf(outfile, "%d\n",    options->stop_flag);
  // fprintf(outfile, "%d\n",    options->prolix_flag);
  // fprintf(outfile, "%d\n",    options->landscape_flag);
  fprintf(outfile, "%d\n",    options->log_flag);
  fprintf(outfile, "%d\n",    options->time_flag);
  fprintf(outfile, "%ld\n",   options->state_write);
  fprintf(outfile, "%ld\n",   options->print_freq);
  fprintf(outfile, "%ld\n",   options->captions);
	fprintf(outfile, "%d\n",    options->precision);
  // fprintf(outfile, "%.16g\n", options->stepsize);
  fprintf(outfile, "%d\n",    options->quenchit);
  fprintf(outfile, "%d\n",    options->equil);
#ifdef MPI
  fprintf(outfile, "%d\n",    options->tuning);
  fprintf(outfile, "%d\n",    options->covar_index);
  fprintf(outfile, "%d\n",    options->write_tune_stat);
  fprintf(outfile, "%d\n",    options->auto_stop_tune);
#endif

  if ( time_flag ) {
    fprintf(outfile, "%.3f\n", delta[0]);
    fprintf(outfile, "%.3f\n", delta[1]);
  }

  fprintf(outfile, "%d\n",    move_status->nparams);
  fprintf(outfile, "%d\n",    move_status->index);
  fprintf(outfile, "%d\n",    move_status->nhits);
  fprintf(outfile, "%d\n",    move_status->nsweeps);

  for(i=0; i < move_status->nparams; i++)
    fprintf(outfile, "%.16g\n", *(move_status->pt[i].param));

  fprintf(outfile, "%.16g\n", move_status->old_energy);

  for(i=0; i < move_status->nparams; i++)
    fprintf(outfile, "%.16g %.16g %d %d\n",
	    move_status->acc_tab_ptr[i].acc_ratio,
	    move_status->acc_tab_ptr[i].theta_bar,
	    move_status->acc_tab_ptr[i].hits,
	    move_status->acc_tab_ptr[i].success);

  for(i=0; i < 31; i++)
    fprintf(outfile,"%.16g\n", lamsave[i]);

  for(i=0; i < 3; i++)
    fprintf(outfile,"%d\n", prand[i]);

  fclose(outfile);

  free(options);
  free(move_status);
  free(lamsave);
  if ( time_flag )
    free(delta);

}



/*** StateRm: removes the state file after the run has been completed; *****
 *            unless we're tuning in parallel, only the root node needs to *
 *            delete a state file                                          *
 ***************************************************************************/

void StateRm(void)
{
  if ( remove(filename) )
    warning("StateRm: could not delete %s", filename);

  free(filename);
}
