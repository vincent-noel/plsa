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

typedef enum StopStyle
{
	proportional_freeze,
	absolute_freeze,
	absolute_energy

} StopStyle;


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
	int 	flag;
	double 	score;
	double *params;

} PLSARes;


typedef struct
{
	long   				seed;
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
	int    				mix_interval;
	#endif

	/* These were marked "Application program must set these." in the code        */
	/* from Greening/Lam; only progname is used at the moment; we kept them       */
	/* mainly for historical reasons (and to write funny things into tunename)    */

	char   				progname[128];                     /* name of my prog */
	int    				debuglevel;                              /* who cares */
	char   				tunename[128];                             /* nothing */


	double      		gain_for_jump_size_control;
	double      		interval;



	int         distribution;        /* move generation distribution type RO */
			 /* 1 - exp; 2 - uni; 3 - absnor; 4 - abs lorentz */
			   /*  LG: 07-05-00 formerly dist_type in lj code */

	double      q;                 /* gen visiting distribution parameter RO */
							  /* 1=guassian; 2=lor; but 1<q<3 */
		   /*  LG: 03-02 need q and factors for GSA visit dist*/
			/*  1   <  q < 2   uses qlt2_visit       */
			/*  2   <  q < 2.6 uses qgt2_visit       */
			/*  2.6 <= q < 3   uses binom_qgt2_visit */

	double      		(*scoreFunction)   	();
	void        		(*printFunction)   	(char * path, int proc);

	SALogs * logs;

} SAType ;

#endif
