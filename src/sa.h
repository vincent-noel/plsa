/******************************************************************************
 *                                                                            *
 *   sa.h                                                                     *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   written by Jimmy Lam and Dan Greening                                    *
 *   modified by John Reinitz and Johannes Jaeger                             *
 *   modified by Vincent Noel                                                 *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   This is the header imported by optimizations, which contains just		  *
 *   What is needed. The rest of the code must stay private. The functions	  *
 *   definitions are in lsa.c (for now).									  *
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
#include "types.h"


/*** FUNCTION PROTOTYPES ******************************************************/


/*** InitPLSA: Initialize the settings variable, are returns it               */
#ifdef MPI
SAType * InitPLSA(int nb_procs, int my_id);
#else
SAType * InitPLSA();
#endif

/*** InitPLSAParameters: Initialize the parameters variable, are returns it   */
PArrPtr * InitPLSAParameters(int nb_dimensions);

/*** runPLSA: Perform the optimization, and returns the result                */
PLSARes * runPLSA();
