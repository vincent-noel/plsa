/*****************************************************************************
 *                                                                           *
 *   config.h                                                                *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   written by JR, modified by Yoginho                                      *
 *   modified by Vincent Noel                                                *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (C) 2016 Vincent Noel                                         *
 *   the full GPL copyright notice can be found in lsa.c                     *
 *                                                                           *
 *****************************************************************************/

/* this def needed for func. defs that refer to (* FILE) *********************/
#ifndef _STDIO_INCLUDED
#include <stdio.h>
#endif

/* following for structures & consts used thruout ****************************/

#ifndef GLOBAL_INCLUDED
#include "global.h"
#endif


///*** GLOBALS ** ************************************************************/

typedef struct
{
	/* Added some log options here */
	char * 	dir;

	int    	trace;
	int    	params;
	int  	res;
	int 	score;
	int 	pid;

} SALogs;

/* Following are functions that are needed by many other reading funcs */

char *	getLogDir();

int		logTrace();
int 	logParams();
int 	logRes();
int 	logScore();
int 	logPid();

void 	InitLogs();
void 	BuildLogs();
