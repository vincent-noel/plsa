/******************************************************************************
 *                                                                            *
 *   main.c                                                                   *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   written by Vincent Noel                                                  *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   This file contains:                                                      *
 *                                                                            *
 *    The first example of optimization, on a one-parameter funnel            *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   Copyright (C) 2016 Vincent Noel (vincent.noel@butantan.gov.br)           *
 *                                                                            *
 *   plsa is free software: you can redistribute it and/or modify     				*
 *   it under the terms of the GNU General Public License as published by     *
 *   the Free Software Foundation, either version 3 of the License, or        *
 *   (at your option) any later version.                                      *
 *                                                                            *
 *   plsa is distributed in the hope that it will be useful,          				*
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *   GNU General Public License for more details.                             *
 *                                                                            *
 *   You should have received a copy of the GNU General Public License        *
 *   along with plsa. If not, see <http://www.gnu.org/licenses/>.     				*
 *                                                                            *
 ******************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#ifdef MPI
	#include <mpi.h>
#endif

#include "../../src/sa.h"


//////////////////////////////////////////////////////////////////////////////
// Problem definition :
int nb_inputs = 6;
double input[6] = {0.50000, 0.80000, 1.10000, 1.40000, 1.70000, 2.00000};
double solution[6] = {0.15246, 0.16622, 0.96904, 2.54871, 2.71236, 2.72350};

double k;
double n;
double theta;
double basal;
long seed;

int nb_tests = 200;

// And the function is the distance to the solution
double 	score_function()
{
	int i;
	double score, t_score;
	score = 0;
	for (i=0; i < nb_inputs; i++)
	{
		t_score = k*pow(input[i], n)/(pow(input[i], n) + pow(theta,n)) + basal;
		score += fabs((t_score-solution[i])/solution[i]);
	}
	return score;
}

// And we print it for each new best score
void 	print_function()
{
	int i;
	double t_score;
	for (i=0; i < nb_inputs; i++)
	{
		t_score = k*pow(input[i], n)/(pow(input[i], n) + pow(theta,n)) + basal;
		printf("%g", t_score);
		if (i < (nb_inputs-1))
			printf(", ");
	}
	printf("\n");
}

//////////////////////////////////////////////////////////////////////////////
// Now the optimization

int 	main (int argc, char ** argv)
{

#ifdef MPI
	// MPI initialization steps
	int nnodes, myid;

	int rc = MPI_Init(NULL, NULL); 	     /* initializes the MPI environment */
	if (rc != MPI_SUCCESS)
		printf (" > Error starting MPI program. \n");

	MPI_Comm_size(MPI_COMM_WORLD, &nnodes);        /* number of processors? */
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);         /* ID of local processor? */


#endif


	int i, success;
	success = 0;
	for (i=0; i < nb_tests; i++)
	{

		k = 1;
		n = 1;
		theta = 1;
		basal = 1;

		srand ( time(NULL) );
		seed = rand();
		seed *= rand();


		// define the optimization settings
#ifdef MPI
		SAType * t_sa = InitPLSA(nnodes, myid);

#else
		SAType * t_sa = InitPLSA();

#endif
		t_sa->seed = seed;
		t_sa->scoreFunction = &score_function;
		t_sa->initial_temp = 1;
		t_sa->lambda = 0.0001;
		t_sa->initial_moves = 20000;
		t_sa->tau = 10000;
		t_sa->interval = 1000;
		t_sa->criterion = 1e-4;
		// define the optimization parameters
		PArrPtr * params = InitPLSAParameters(4);
		params->array[0] = (ParamList) { &k, (Range) {0,1e+16}, "k"};
		params->array[1] = (ParamList) { &n, (Range) {0,1e+16}, "n"};
		params->array[2] = (ParamList) { &theta, (Range) {0,1e+16}, "theta"};
		params->array[3] = (ParamList) { &basal, (Range) {0,1e+16}, "basal"};

		// run the optimization
		PLSARes * res = runPLSA();

		// print final parameter value and score
	#ifdef MPI
		if (myid == 0)
		{
	#endif

			printf("final score : %g\n", res->score);

	#ifdef MPI
		}
	#endif

		if (res->score < 1e-4)
			success++;


		free(params->array);
		free(res->params);
		// free(params);
		free(res);
	}
#ifdef MPI
	if (myid == 0)
	{
#endif
		printf("> Success = %2.0f%%\n", (double) success*100/nb_tests);
#ifdef MPI
	}

	// terminates MPI execution environment
	MPI_Finalize();
#endif

	return 0;
}
