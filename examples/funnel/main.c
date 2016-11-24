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
 *   the full GPL copyright notice can be found in lsa.c                     *
 *                                                                           *
 *****************************************************************************/


#include <stdlib.h>
#include <plsa/sa.h>
#include <plsa/config.h>
#include <math.h>


double param = 1e-8;
PArrPtr params;
plsa_parameters * settings;


double function()
{

	return pow((42.4203623-param),2);
}


void print_function()
{
	;
}

int main (char * argv, int argc)
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
	settings->freeze_count = 1000;
	settings->update_S_skip = 1;
	settings->criterion = 0.0001;
	settings->mix_interval = 10;
	settings->distribution = 1;
	settings->q = 1;
	settings->log_trace = 0;
	settings->log_params = 0;
	settings->scoreFunction = &function;
	settings->printFunction = &print_function;

	ParamList *p;
	int t_size = 1;

	p = (ParamList *)malloc(t_size * sizeof(ParamList));

	p[0].param       = &param;
	p[0].param_range = (Range){0, 100};


	params.size  = t_size;
	params.array = p;

  	run(settings, &params);

	printf("final value : %16g\n", param);

	return 0;
}
