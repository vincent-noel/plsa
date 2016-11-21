/*****************************************************************************
 *                                                                           *
 *   error.h                                                                 *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   written by JR, modified by Yoginho                                      *
 *   modified by Vincent Noel                                                *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   contains error and warning functions                                    *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (C) 2016 Vincent Noel                                         *
 *   the full GPL copyright notice can be found in lsa.c                     *
 *                                                                           *
 *****************************************************************************/


#ifndef ERROR_INCLUDED
#define ERROR_INCLUDED

/* following for structures & consts used thruout */
#ifndef GLOBAL_INCLUDED
#include "global.h"
#endif


/*** A CONSTANT ************************************************************/

#define MAX_ARGS     25   /* max number of error() and warning() arguments */




/*** FUNCTION PROTOTYPES ***************************************************/

/*** The following two routines print error messages with value. 'error' ***
 *   then exits, while 'warning' returns to the calling function.          *
 *                                                                         *
 *   Both functions can only handle %c, %d, %f, %g, %s (char, string, int  *
 *   (includes long integers too) & double); it does not yet work for      *
 *   floats or long doubles or the more esoteric format specifiers of      *
 *   printf.                                                               *
 *                                                                         *
 *   No % sign means just print msg and exit.                              *
 *                                                                         *
 ***************************************************************************/

void error  (const char *format, ... );
void warning(const char *format, ... );

/*** file_error prints an file handling error using perror(); it only ******
 *   needs the name of the calling function as an argument                 *
 ***************************************************************************/

void file_error(const char *call_name);

/*** PrintMsg: prints a string to stderr, then quits with exit status 0; ***
 *             useful for printing help and usage messages                 *
 ***************************************************************************/

void PrintMsg(const char *msg, int exit_status);

#endif
