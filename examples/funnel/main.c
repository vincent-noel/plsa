/*****************************************************************************
 *                                                                           *
 *   main.c                                                                  *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   written by Vincent Noel                                                 *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   This file contains:                                                     *
 *                                                                           *
 *    The first example of optimization, on a one-parameter funnel           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (C) 2016 Vincent Noel                                         *
 *   the full GPL copyright notice can be found in src/lsa.c                 *
 *                                                                           *
 *****************************************************************************/


#include <stdlib.h>
#include <plsa/sa.h>
#include <plsa/config.h>
#include <math.h>

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
	//printf("New best score : %10.7f\n", param);
}


//////////////////////////////////////////////////////////////////////////////
// Now the optimization

// The list of parameters, and the settings variable
PArrPtr params;
plsa_parameters * settings;

int		main (char * argv, int argc)
{

	setLogDir("logs");

	settings = malloc(sizeof(plsa_parameters));
	settings->seed = -6.60489e+08;
	settings->initial_temperature = 1;
	settings->gain_for_jump_size_control = 5;
	settings->interval = 100;
	settings->lambda = 0.001;
	settings->lambda_mem_length_u = 200;
	settings->lambda_mem_length_v = 1000;
	settings->control = 1;
	settings->initial_moves = 200;
	settings->tau = 100;
	settings->freeze_count = 100;
	settings->update_S_skip = 1;
	settings->criterion = 0.001;
	settings->mix_interval = 10;
	settings->distribution = 1;
	settings->q = 1;
	settings->log_trace = 0;
	settings->log_params = 0;
	settings->scoreFunction = &score_function;
	settings->printFunction = &print_function;


	int t_size = 1;
	ParamList *p = (ParamList *)malloc(t_size * sizeof(ParamList));

	p[0] = (ParamList) { &param, (Range) {0,1e+16}};

	params.size  = t_size;
	params.array = p;

  	run(settings, &params);

	printf("final value : %10.7f\n", param);

	return 0;
}
