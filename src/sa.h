/******************************************************************************
 *                                                                            *
 *   sa.h                                                                     *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   written by Jimmy Lam and Dan Greening                                    *
 *   modified by John Reinitz and Johannes Jaeger                             *
 *   modified by Vincent Noel                                                 *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   see ../doc/prog_man.ps for details of how this works                     *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   IMPORTANT: IF YOU EVER CHANGE ANYTHING IN THIS FILE, LET ALL             *
 *            YOUR FELLOW PROGRAMMERS KNOW WELL IN ADVANCE AND                *
 *            CONSULT WITH THEM IF THEY AGREE ON YOUR CHANGES!!               *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   The following came with lsa.c, obtained from Dan Greening,               *
 *   who got it from Jimmy Lam. It is probably copyright Jimmy Lam            *
 *   1988, with modifications by Greening. We (JR & JJ) have                  *
 *   mostly modified it by adding some comments and restructuring             *
 *   the code to make it more easily understandable. JR has added             *
 *   the criterion for continous search spaces and separated Lam              *
 *   statistics from move acceptance statistics. King-Wai Chu has             *
 *   written the equilibration code.                                          *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   NOTE: this header only contains prototypes for functions used            *
 *       for serial annealing code; all prototypes that are spe-              *
 *       cific to parallel code need to go into MPI.h                         *
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

} SAType ;



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
  char *	name;
} ParamList;


typedef struct {
  int       size;                           /* size of the ParamList array */
  ParamList *array;            /* points to 1st element of ParamList array */
} PArrPtr;

typedef struct {
	int 	flag;
	double 	score;
	double *params;
} PLSARes;

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
PLSARes * runPLSA();
