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

#ifndef CONFIG_INCLUDED
#define CONFIG_INCLUDED

/* this def needed for func. defs that refer to (* FILE) *******************/
#ifndef _STDIO_INCLUDED
#include <stdio.h>
#endif

/* following for structures & consts used thruout **************************/

#ifndef GLOBAL_INCLUDED
#include "global.h"
#endif


///*** GLOBALS ** ************************************************************/

int    debug;                                            /* debugging flag */

/* Following are functions that are needed by many other reading funcs */

char *      getLogDir();
void        setLogDir(char * dir);

#endif
