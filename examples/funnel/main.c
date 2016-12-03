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

#include "../../src/sa.h"
#include <math.h>

#ifdef MPI
#include <mpi.h>
#endif


//////////////////////////////////////////////////////////////////////////////
// Problem definition :

// One parameter, one solution
double param = 1e-8;
double solution = 45651632.4203623;

// And the function is the distance to the solution
double 	score_function()
{
	return pow((solution-param),2);
}

// And we print it for each new best score
void 	print_function()
{
	// printf("New best score : %.5g (param = %10.7f)\n",
	// 		score_function(),
	// 		param
	// );
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


	// define the optimization settings
#ifdef MPI
	SAType * t_sa = InitPLSA(nnodes, myid);
#else
	SAType * t_sa = InitPLSA();
#endif
	t_sa->scoreFunction = &score_function;
	t_sa->printFunction = &print_function;


	// define the optimization parameters
	PArrPtr * params = InitPLSAParameters(1);
	params->array[0] = (ParamList) { &param, (Range) {0,1e+16}, "k"};


	// run the optimization
	PLSARes * res = runPLSA();


	// print final parameter value and score
#ifdef MPI
	if (myid == 0)
	{
#endif
		printf("final value : %10.7f\n", res->params[0]);
		printf("final score : %g\n", res->score);
#ifdef MPI
	}


	// terminates MPI execution environment
	MPI_Finalize();
#endif

	free(res->params);
	free(res);

	return 0;
}
