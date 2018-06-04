/******************************************************************************
 *																			  *
 *   state.c                                                			 	  *
 *                                                              			  *
 ******************************************************************************
 *                                                               			  *
 *   written by JR, modified by Yoginho                         			  *
 *   modified by Vincent Noel                                                 *
 *                                                              			  *
 ******************************************************************************
 *                                                              			 				*
 *   state.c contains three functions that read, write and remove the    	*
 *	 state file for an annealing run. The frequency with which state are	 		*
 *   saved can be chosen by the command line option -b						 						*
 *   (for backup stepsize). The state file is very useful for the case		 		*
 *   when long annealing runs have to be interrupted or crash for some				*
 *   reason or another. The run can then be resumed by indicating the		 			*
 *	 state file as an additional argument to plsa.                      	 		*
 *                                                              			 				*
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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>                                        /* getopt stuff */

#include "state.h"
#include "error.h"
#include "random.h"

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
		   LamState *lam_state, unsigned short *rand, TuningSettings * t_settings, double *delta)
{
  int            i;                                  /* local loop counter */
  FILE           *infile;                         /* pointer to state file */


/* make the state file name static to state.c */

  filename = (char *)calloc(MAX_RECORD, sizeof(char));
  filename = strcpy(filename, statefile);

/* open the state file and read it */

  infile = fopen(filename, "r");
  if( !infile )
	file_error("StateRead");

  fscanf(infile, "%d\n",  (int*) &(options->stop_flag));
  fscanf(infile, "%d\n",  &(options->time_flag));
  fscanf(infile, "%ld\n", &(options->state_write));
  fscanf(infile, "%ld\n", &(options->print_freq));
  fscanf(infile, "%d\n",  &(options->quenchit));

#ifdef MPI
  fscanf(infile, "%d\n",  &(options->tuning));
  // fscanf(infile, "%d\n",  &(options->covar_index));
  // fscanf(infile, "%d\n",  &(options->write_tune_stat));
  // fscanf(infile, "%d\n",  &(options->auto_stop_tune));
#endif

  if ( options->time_flag ) {
	fscanf(infile, "%lf\n", &(delta[0]));
	fscanf(infile, "%lf\n", &(delta[1]));
  }

  fscanf(infile, "%d\n",  &(options->max_iter));
  fscanf(infile, "%d\n",  &(options->max_seconds));

  fscanf(infile, "%d\n",  &(move_ptr->nparams));
  fscanf(infile, "%d\n",  &(move_ptr->index));
  fscanf(infile, "%d\n",  &(move_ptr->nhits));
  fscanf(infile, "%d\n",  &(move_ptr->nsweeps));

  // Here we put the values in a double array, and not in the ParamList array
  // TODO : just remove the ParamList, and always use the doubles ?
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

  // for(i=0; i<33; i++)
  // scanf(infile, "%lg\n", &(stats[i]));

	fscanf(infile, "%d\n", &(lam_state->counter));

	fscanf(infile, "%lg\n", &(lam_state->old_mean));
	fscanf(infile, "%lg\n", &(lam_state->energy));

	fscanf(infile, "%lg\n", &(lam_state->mean));
	fscanf(infile, "%lg\n", &(lam_state->vari));

	fscanf(infile, "%lg\n", &(lam_state->estimate_mean));
	fscanf(infile, "%lg\n", &(lam_state->estimate_sd));

	fscanf(infile, "%lg\n", &(lam_state->S));
	fscanf(infile, "%lg\n", &(lam_state->dS));
	fscanf(infile, "%lg\n", &(lam_state->S_0));

	fscanf(infile, "%lg\n", &(lam_state->alpha));
	fscanf(infile, "%lg\n", &(lam_state->acc_ratio));

	fscanf(infile, "%lg\n", &(lam_state->w_b));
	fscanf(infile, "%lg\n", &(lam_state->vsyy));
	fscanf(infile, "%lg\n", &(lam_state->vsxy));
	fscanf(infile, "%lg\n", &(lam_state->vsxx));
	fscanf(infile, "%lg\n", &(lam_state->vsx));
	fscanf(infile, "%lg\n", &(lam_state->vsy));
	fscanf(infile, "%lg\n", &(lam_state->vsum));
	fscanf(infile, "%lg\n", &(lam_state->D));
	fscanf(infile, "%lg\n", &(lam_state->E));

	fscanf(infile, "%lg\n", &(lam_state->w_a));
	fscanf(infile, "%lg\n", &(lam_state->usyy));
	fscanf(infile, "%lg\n", &(lam_state->usxy));
	fscanf(infile, "%lg\n", &(lam_state->usxx));
	fscanf(infile, "%lg\n", &(lam_state->usx));
	fscanf(infile, "%lg\n", &(lam_state->usy));
	fscanf(infile, "%lg\n", &(lam_state->usum));
	fscanf(infile, "%lg\n", &(lam_state->A));
	fscanf(infile, "%lg\n", &(lam_state->B));

	fscanf(infile, "%ld\n", &(lam_state->count_tau));

	fscanf(infile, "%d\n", &(lam_state->proc_init));
	fscanf(infile, "%d\n", &(lam_state->proc_tau));

	for(i=0; i<3; i++)
		fscanf(infile, "%hu\n", &(rand[i]));

#ifdef MPI
	fscanf(infile, "%d\n", &(t_settings->covar_index));
	fscanf(infile, "%d\n", &(t_settings->write_tune_stat));
	fscanf(infile, "%d\n", &(t_settings->auto_stop_tune));
#endif

	fclose(infile);

}



/*** StateWrite: collects command line options and arguemnts, Lam statis- **
 *               tics, move state and the state of the erand48 random num- *
 *               ber generator and writes all that into the state file,    *
 *               which can then be used to restore the run in case it gets *
 *               interrupted                                               *
 ***************************************************************************/

void StateWrite(char *statefile, double energy)
{
	int i;                                             /* local loop counter */
	FILE           *outfile;                           /* state file pointer */
	Opts           *options;                /* command line opts to be saved */
	MoveState      *move_status;                    /* MoveState to be saved */
#ifdef MPI
	TuningSettings *t_settings;
#endif
	LamState       *lam_state;                        /* Lam stats to be saved */
	unsigned short *prand;                    /* erand48() state to be saved */
	double         *delta;            /* wallclock and user time to be saved */


	/* if StateWrite() called for the first time: make filename static */

	if (filename == NULL)
	{
		filename = (char *)calloc(MAX_RECORD, sizeof(char));
		filename = strcpy(filename, statefile);
	}

	/* collect the state and the options */

	options     = GetOptions();
	move_status = MoveSave();
	lam_state     = GetLamstats(energy);
	prand       = GetERandState();
	if ( options->time_flag )
		delta     = GetTimes();

#ifdef MPI
	t_settings = GetTuningSettings();
#endif


	/* write the answer; now *fully* portable, no binary!!! */

	outfile = fopen(filename, "w");
	if ( !outfile )
		file_error("StateWrite");


	/* Start of the Opts struct */
	fprintf(outfile, "%d\n",    options->stop_flag);
	fprintf(outfile, "%d\n",    options->time_flag);
	fprintf(outfile, "%ld\n",   options->state_write);
	fprintf(outfile, "%ld\n",   options->print_freq);
	fprintf(outfile, "%d\n",    options->quenchit);
#ifdef MPI
	fprintf(outfile, "%d\n", 	options->tuning);
#endif

	if ( options->time_flag )
	{
		fprintf(outfile, "%.3f\n", delta[0]);
		fprintf(outfile, "%.3f\n", delta[1]);
	}

	fprintf(outfile, "%d\n", 	options->max_iter);
	fprintf(outfile, "%d\n", 	options->max_seconds);


	/* Start of the move_state struct */
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


	/* LAM stats*/

	fprintf(outfile, "%d\n", lam_state->counter);

	fprintf(outfile, "%.16g\n", lam_state->old_mean);
	fprintf(outfile, "%.16g\n", lam_state->energy);

	fprintf(outfile, "%.16g\n", lam_state->mean);
	fprintf(outfile, "%.16g\n", lam_state->vari);

	fprintf(outfile, "%.16g\n", lam_state->estimate_mean);
	fprintf(outfile, "%.16g\n", lam_state->estimate_sd);

	fprintf(outfile, "%.16g\n", lam_state->S);
	fprintf(outfile, "%.16g\n", lam_state->dS);
	fprintf(outfile, "%.16g\n", lam_state->S_0);

	fprintf(outfile, "%.16g\n", lam_state->alpha);
	fprintf(outfile, "%.16g\n", lam_state->acc_ratio);

	fprintf(outfile, "%.16g\n", lam_state->w_b);
	fprintf(outfile, "%.16g\n", lam_state->vsyy);
	fprintf(outfile, "%.16g\n", lam_state->vsxy);
	fprintf(outfile, "%.16g\n", lam_state->vsxx);
	fprintf(outfile, "%.16g\n", lam_state->vsx);
	fprintf(outfile, "%.16g\n", lam_state->vsy);
	fprintf(outfile, "%.16g\n", lam_state->vsum);
	fprintf(outfile, "%.16g\n", lam_state->D);
	fprintf(outfile, "%.16g\n", lam_state->E);

	fprintf(outfile, "%.16g\n", lam_state->w_a);
	fprintf(outfile, "%.16g\n", lam_state->usyy);
	fprintf(outfile, "%.16g\n", lam_state->usxy);
	fprintf(outfile, "%.16g\n", lam_state->usxx);
	fprintf(outfile, "%.16g\n", lam_state->usx);
	fprintf(outfile, "%.16g\n", lam_state->usy);
	fprintf(outfile, "%.16g\n", lam_state->usum);
	fprintf(outfile, "%.16g\n", lam_state->A);
	fprintf(outfile, "%.16g\n", lam_state->B);

	fprintf(outfile, "%ld\n", lam_state->count_tau);

	fprintf(outfile, "%d\n", lam_state->proc_init);
	fprintf(outfile, "%d\n", lam_state->proc_tau);


	/* The three values to initialize ERand I guess */
	for(i=0; i < 3; i++)
		fprintf(outfile,"%d\n", prand[i]);


#ifdef MPI
	fprintf(outfile, "%d\n", t_settings->covar_index);
	fprintf(outfile, "%d\n", t_settings->write_tune_stat);
	fprintf(outfile, "%d\n", t_settings->auto_stop_tune);
#endif

	fclose(outfile);
	free(move_status);
	free(lam_state);
	if ( options->time_flag )
		free(delta);
	free(options);


#ifdef MPI
	free(t_settings);
#endif
}



/*** StateRm: removes the state file after the run has been completed; *****
 *            unless we're tuning in parallel, only the root node needs to *
 *            delete a state file                                          *
 ***************************************************************************/

void StateRm(void)
{
  if ( remove(filename) )
	warning("StateRm: could not delete %s", filename);

  // free(filename);
}
