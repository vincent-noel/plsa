/******************************************************************************
 *                                                                            *
 *   problem.c                                                                *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   written by Vincent Noel                                                  *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   This file contains:                                                      *
 *                                                                            *
 *    The second example of optimization, on a 4 parameters sigmoid           *
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
#include <stdio.h>
#include <math.h>
#include "problem.h"

//////////////////////////////////////////////////////////////////////////////
// Problem definition :
int nb_inputs = 6;
double input[6] = {0.50000, 0.80000, 1.10000, 1.40000, 1.70000, 2.00000};
double solution[6] = {0.15246, 0.16622, 0.96904, 2.54871, 2.71236, 2.72350};


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
