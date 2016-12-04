/******************************************************************************
 *                                                                            *
 *   config.c                                                                 *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   written by JR, modified by Yoginho                                       *
 *   modified by Vincent Noel                                                 *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   This file contains:                                                      *
 *                                                                            *
 *    some small I/O functions (FindSection(), KillSection())                 *
 *    that are used throughout in the fly code                                *
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
#include <sys/stat.h>
#include <sys/types.h>

#include "config.h"

SALogs logs;
SALogs * getLogSettings()
{
	return &logs;
}

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

int logBestScore()
{
	return logs.best_score;
}

int logBestRes()
{
	return logs.best_res;
}

void InitLogs()
{
	logs.dir = "logs";
	logs.trace = 0;
	logs.params = 0;
	logs.res = 0;
	logs.score = 0;
	logs.pid = 0;
	logs.best_score = 0;
	logs.best_res = 0;
}

void BuildLogs()
{
	mkdir(logs.dir, 0777);

	if (logScore() > 0)
	{
		char score[200];
		sprintf(score, "%s/score", logs.dir);
		mkdir(score, 0777);
	}

	if (logRes() > 0)
	{
		char res[250];
		sprintf(res, "%s/res", logs.dir);
		mkdir(res, 0777);
	}

	if (logParams() > 0)
	{
		char res[250];
		sprintf(res, "%s/params", logs.dir);
		mkdir(res, 0777);

	}

	if (logBestScore() > 0)
	{
		char score[200];
		sprintf(score, "%s/best_score", logs.dir);
		mkdir(score, 0777);
	}

	if (logBestRes() > 0)
	{
		char res[250];
		sprintf(res, "%s/best_res", logs.dir);
		mkdir(res, 0777);

	}

}
