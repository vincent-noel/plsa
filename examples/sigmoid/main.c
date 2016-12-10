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


//////////////////////////////////////////////////////////////////////////////
// Problem definition :
int nb_inputs = 13;
double ras_input[13] = {0.27473296046875062, 0.75856928418982039,
							0.92019179219384561, 1.0010676507878311,
							1.0496138592497219, 1.0819866250315191,
							1.1051142559845317, 1.1359565780285732,
							1.1555867723155182, 1.1867693559230004,
							1.1976600204978742, 1.2125920395855168,
							1.2202407154431443};
double solution[13] = {0.15243705991527873, 0.15900202466660762,
						 0.24715087493852866, 0.43688580407507915,
						0.65242365783591083, 0.84597284362150726,
						1.0056920670435818, 1.2376189221670717,
						1.3902345246714087, 1.6285371043491492,
						1.7081269468499927, 1.812622261520743,
						1.8637588068489317};

double k = 1e-6;
double n = 1e-6;
double theta = 1e-6;
double ras_basal = 1e-6;



// And the function is the distance to the solution
double 	score_function()
{
	int i;
	double score, t_score;
	score = 0;
	for (i=0; i < nb_inputs; i++)
	{
		t_score = k*pow(ras_input[i], n)/(pow(ras_input[i], n) + pow(theta,n)) + ras_basal;
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
		t_score = k*pow(ras_input[i], n)/(pow(ras_input[i], n) + pow(theta,n)) + ras_basal;
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

	int nnodes, myid;

	// define the optimization settings
	SAType * t_sa = InitPLSA(&nnodes, &myid);
	t_sa->scoreFunction = &score_function;
	t_sa->lambda = 0.0001;
	t_sa->initial_moves = 2000;
	t_sa->tau = 1000;
	t_sa->criterion = 0;

	// define the optimization parameters
	PArrPtr * params = InitPLSAParameters(4);
	params->array[0] = (ParamList) { &k, (Range) {0,1e+16}, "k"};
	params->array[1] = (ParamList) { &n, (Range) {0,1e+16}, "n"};
	params->array[2] = (ParamList) { &theta, (Range) {0,1e+16}, "theta"};
	params->array[3] = (ParamList) { &ras_basal, (Range) {0,1e+16}, "ras_basal"};

	// run the optimization
	PLSARes * res = runPLSA();

	// print final parameter value and score
#ifdef MPI
	if (myid == 0)
	{
#endif
		printf("final score : %g\n", res->score);
		// printf("k : %g\n", k);
		// printf("n : %g\n", n);
		// printf("theta : %g\n", theta);
		// printf("ras_basal : %g\n", ras_basal);
		//
		// int i;
		// for (i=0; i < nb_inputs; i++)
		// {
		// 	printf("%g", solution[i]);
		// 	if (i < (nb_inputs-1))
		// 		printf(", ");
		// }
		// printf("\n");
		// print_function();
#ifdef MPI
	}
#endif

	free(res->params);
	free(res);
	free(params->array);

	return 0;
}
