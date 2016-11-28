/******************************************************************************
 *                                                                            *
 *   error.h                                                                  *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   written by JR, modified by Yoginho                                       *
 *   modified by Vincent Noel                                                 *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   contains error and warning functions                                     *
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

/* following for structures & consts used thruout *****************************/
#ifndef GLOBAL_INCLUDED
#include "global.h"
#endif


/*** A CONSTANT ***************************************************************/

#define MAX_ARGS     25   /* max number of error() and warning() arguments ****/




/*** FUNCTION PROTOTYPES *****************************************************/

/*** The following two routines print error messages with value. 'error' *****
 *   then exits, while 'warning' returns to the calling function.            *
 *                                                                           *
 *   Both functions can only handle %c, %d, %f, %g, %s (char, string, int    *
 *   (includes long integers too) & double); it does not yet work for        *
 *   floats or long doubles or the more esoteric format specifiers of        *
 *   printf.                                                                 *
 *                                                                           *
 *   No % sign means just print msg and exit.                                *
 *                                                                           *
 *****************************************************************************/

void error  (const char *format, ... );
void warning(const char *format, ... );

/*** file_error prints an file handling error using perror(); it only ********
 *   needs the name of the calling function as an argument                   *
 *****************************************************************************/

void file_error(const char *call_name);

/*** PrintMsg: prints a string to stderr, then quits with exit status 0; *****
 *             useful for printing help and usage messages                   *
 *****************************************************************************/

void PrintMsg(const char *msg, int exit_status);
