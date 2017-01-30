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
#include <time.h>
#include <sys/time.h>
#ifdef MPI
	#include <mpi.h>
#endif

#include "../src/sa.h"
#include "../examples/sigmoid/problem.h"

long seed;
int nb_tests = 1000;


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
	FILE * f_res;

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
		params->array[0] = (ParamList) { &k, 1.0, (Range) {0,1e+16}, "k"};
		params->array[1] = (ParamList) { &n, 1.0, (Range) {0,1e+16}, "n"};
		params->array[2] = (ParamList) { &theta, 1.0, (Range) {0,1e+16}, "theta"};
		params->array[3] = (ParamList) { &basal, 1.0, (Range) {0,1e+16}, "basal"};

		// run the optimization
		PLSARes * res = runPLSA();

		if (res->score < 1e-3)
			success++;

#ifdef MPI
		if (myid == 0)
		{
#endif
			f_res = fopen("res", "a");

			fprintf(f_res, "%g, %ld, %lu\n", res->score, res->niters, res->duration);
			fclose(f_res);



#ifdef MPI
		}
#endif




		free(params->array);
		free(res->params);
		free(res);
	}

#ifdef MPI
	if (myid == 0)
	{
#endif




#ifdef MPI
	}

	// terminates MPI execution environment
	MPI_Finalize();
#endif

	if (success >= nb_tests*0.75)
		return 0;
	else
		return 1;
}
