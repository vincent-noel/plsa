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

static SALogs * log_settings;

char *      getLogDir()
{
	return log_settings->dir;
}

int logTraceScore()
{
	return log_settings->trace_score;
}

int logTraceParams()
{
	return log_settings->trace_params;
}

int logParams()
{
	return log_settings->params;
}

int logRes()
{
	return log_settings->res;
}

int logScore()
{
	return log_settings->score;
}

int logPid()
{
	return log_settings->pid;
}

int logBestScore()
{
	return log_settings->best_score;
}

int logBestRes()
{
	return log_settings->best_res;
}

void InitializeLogs(SALogs * t_logs)
{
	log_settings = t_logs;
	
	mkdir(log_settings->dir, 0777);

	if (logScore() > 0)
	{
		char score[200];
		sprintf(score, "%s/score", log_settings->dir);
		mkdir(score, 0777);
	}

	if (logRes() > 0)
	{
		char res[250];
		sprintf(res, "%s/res", log_settings->dir);
		mkdir(res, 0777);
	}

	if (logParams() > 0)
	{
		char res[250];
		sprintf(res, "%s/params", log_settings->dir);
		mkdir(res, 0777);

	}

	if (logTraceScore() > 0)
	{
		char score[200];
		sprintf(score, "%s/trace/score", log_settings->dir);
		mkdir(score, 0777);
	}

	if (logBestRes() > 0)
	{
		char res[250];
		sprintf(res, "%s/trace/res", log_settings->dir);
		mkdir(res, 0777);

	}

	if (logBestScore() > 0)
	{
		char score[200];
		sprintf(score, "%s/best_score", log_settings->dir);
		mkdir(score, 0777);
	}

	if (logBestRes() > 0)
	{
		char res[250];
		sprintf(res, "%s/best_res", log_settings->dir);
		mkdir(res, 0777);

	}

}
