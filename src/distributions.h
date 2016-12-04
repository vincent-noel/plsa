/******************************************************************************
 *                                                                            *
 *   distributions.h                                                          *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   created 5-02 for function prototypes                                     *
 *   found in distributions.c                                                 *
 *   used in tsp_sa.c (for qgt2_init and qlt2_init), and move.c               *
 *           need pi for distributions                                        *
 *           used in input file to set factors, in distributions              *
 *   modified by Vincent Noel                                                 *
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
typedef struct my_distrib
{
	int distribution;      /* move generation distribution type RO */
						 /* 1 - exp; 2 - uni; 3 - absnor; 4 - abs lorentz */
						 /*  LG: 07-05-00 formerly dist_type in lj code */
	double q;              /* gen visiting distribution parameter RO */
						 /* 1=guassian; 2=lor; but 1<q<3 */
						 /*  LG: 03-02 need q and factors for GSA visit dist*/
						 /*  1   <  q < 2   uses qlt2_visit       */
						 /*  2   <  q < 2.6 uses qgt2_visit       */
						 /*  2.6 <= q < 3   uses binom_qgt2_visit */

	/*****variables that depend on q ********/
	double gam1;           /* gammaln of 1/(q-1)       */
	double gam2;           /* gammaln of 1/((q-1)-0.5) */
	double fact2;          /* exp(gam1-gam2)           */
	double alpha;          /* sqrt((q-1)/pi) * exp(gam1-gam2) NOT dependent on theta_bar*/
	double alpha2;         /* (sqrt(2.)/ 2. )* sqrt((q-1)/pi) * exp(gam1-gam2) NOT dependent on theta_bar*/
	double c;              /*  1/(q-1)                  */
	double rejects;        /*  number of recursive calls; useful for efficiency analysis */
	double trunc;          /*  range for x is [-trunc, trunc]             */
						 /*  2    <  q < 2.6   trunc is OK at 2 million */
						 /*  2.6  <= q < 2.85  trunc is 99 and 16 zeros*/
						 /*  2.85 <= q < 3     trunc is 24 nines */




} DistParms;
 DistParms DistP;   /* variables for distributions */

/* prototype functions for distributions */

double generate_dev(double theta_bar, int distribution, double q);
