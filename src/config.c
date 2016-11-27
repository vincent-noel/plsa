/*****************************************************************************
 *                                                                           *
 *   config.c                                                                *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   written by JR, modified by Yoginho                                      *
 *   modified by Vincent Noel                                                *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   This file contains:                                                     *
 *                                                                           *
 *    some small I/O functions (FindSection(), KillSection())                *
 *    that are used throughout in the fly code                               *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (C) 2016 Vincent Noel                                         *
 *   the full GPL copyright notice can be found in lsa.c                     *
 *                                                                           *
 *****************************************************************************/

#include <sys/stat.h>
#include <sys/types.h>

#include "config.h"

SALogs logs;


char *      getLogDir()
{
    return logs.dir;
}

int logTrace()
{
	return logs.trace;
}

int logParams()
{
	return logs.params;
}

int logRes()
{
	return logs.res;
}

int logScore()
{
	return logs.score;
}

int logPid()
{
	return logs.pid;
}

void InitLogs()
{
	logs.dir = "logs";
	logs.trace = 0;
	logs.params = 0;
	logs.res = 0;
	logs.score = 0;
	logs.pid = 0;
}

void BuildLogs()
{
    mkdir(logs.dir, 0777);

    // char res[250];
    // sprintf(res, "%s/res", logs.dir);
    // mkdir(res, 0777);
	//
    // char score[200];
    // sprintf(score, "%s/score", logs.dir);
    // mkdir(score, 0777);

}
