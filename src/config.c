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


#include "config.h"


char   *log_dir;

char *      getLogDir()
{
    return log_dir;
}

void      setLogDir(char * dir)
{
    log_dir = dir;
}
