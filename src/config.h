/******************************************************************************
 *                                                                            *
 *   config.h                                                                 *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   written by JR, modified by Yoginho                                       *
 *   modified by Vincent Noel                                                 *
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

/* this def needed for func. defs that refer to (* FILE) **********************/
#ifndef _STDIO_INCLUDED
#include <stdio.h>
#endif

/* following for structures & consts used thruout *****************************/

#ifndef GLOBAL_INCLUDED
#include "global.h"
#endif

#include "types.h"
/* GLOBALS ********************************************************************/


/* Following are functions that are needed by many other reading funcs */
SALogs * getLogSettings();

char *	getLogDir();

int		logTrace();
int 	logParams();
int 	logRes();
int 	logScore();
int 	logPid();

int 	logBestScore();
int 	logBestRes();

void 	InitLogs();
void 	BuildLogs();
