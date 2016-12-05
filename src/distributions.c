/******************************************************************************
 *                                                                            *
 *   distributions.c                                                          *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *   Created 5-02 by: Lorraine Greenwald                                      *
 *   used with deviates.c to generate only deviates                           *
 *   used with the SA applications for move generation                        *
 *   uses RandomReal  which uses and modifies xsubj                           *
 *   elsewhere i.e. move.c or deviates.c                                      *
 *   takes distribution type q value theta_bar                                *
 *   and output file on command line                                          *
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

#ifdef ICC
#include <mathimf.h>
#else
#include <math.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include "random.h"   /* for random number generation */

#include "distributions.h"   /* problem independent distibutions need q_init prototype */

static double pi=3.14159265;     /* pi used for several distributions */
static DistParms * params_dist;   /* variables for distributions */

void InitDistribution(DistParms * dist)
{
	params_dist = dist;
}
/*******************************************
* gasdev returns a normally distributed
* deviate with zero mean and unit variance
* use RandomReal() rather than ran1
* taken from Numerical Recipes in C p.289
*******************************************/
double gasdev(void)
{
	/* begin gasdev */
	static double gset;
	double fac,rsq,v1,v2;

	static int iset=0;
	if (iset == 0)
	{
		do
		{
			v1 = 2.0*RandomReal()-1.0;
			v2 = 2.0*RandomReal()-1.0;
			rsq = v1*v1 + v2*v2;
		}
		while (rsq >= 1.0 || rsq == 0.0);

		fac=sqrt(-2.0*log(rsq)/rsq);
		gset = v1 * fac;
		iset=1;
		return v2*fac;
	}
	else
	{
		iset = 0;
		return gset;
	}
}  /* end gasdev */

/*********************************************
* gammln returns the ln(gamma(xx)) for xx > 0
* uses ln because without it too many overflows
* taken from Numerical Recipes in C p.214
**********************************************/
double gammln(double xx)
{
	/* begin gammln */
	double x,y,tmp,ser;
	static double cof[6] = {76.18009172947146, -86.50532032941677,
							24.01409824083091, -1.231739572450155,
							0.1208650973866179e-2, -0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;

	for (j=0; j<=5; j++)
		ser += cof[j]/++y;

	return -tmp+log(2.5066282746310005*ser/x);
}  /* end gammln */

/******************************************************
* poidev returns as a float an integer value
* that is a random deviate from a Poisson distribution
* with mean = xm.  Uses RandomReal() rather than ran1
* uses ln because without it too many overflows
* taken from Numerical Recipes in C p.294
*******************************************************/
double poidev(double xm)
{
	/* begin poidev */
	static double sq,alxm,g,oldm=(-1.0);
	double em,t,y;

	if (xm < 12.0)
	{
		if (xm != oldm)
		{
			oldm = xm;
			g = exp(-xm);
		}

		em = -1;
		t = 1.0;

		do
		{
			++em;
			t *= RandomReal();
		}
		while (t > g);
	}
	else
	{
		if (xm != oldm)
		{
			oldm = xm;
			sq = sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-gammln(xm+1.0);
		}

		do
		{
			do
			{
				y = tan(pi*RandomReal());
				em = sq*y+xm;
			}
			while (em < 0.0);

			em = floor(em);
			t = 0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		}
		while (RandomReal() > t);
	}
	return em;
}   /* end poidev */



/**************************************************************
* generate_dev generates basic deviates for any distribution  *
* The distribution type and theta_bar are input as parameters *
* deviate is the value returned.                              *
***************************************************************/
double generate_dev( double theta_bar)
{
	/* begin generate_dev*/

	double xi;  /* uniform random variable for deviates */
	double theta; /* resulting deviate */

/******************************************************************************
 *																			  *
 *   1  = 	exponential        exp(-x/mu)    ORIGINAL lam's WAY				  *
 *																			  *
 *   2  = 	controlled uniform        interval=[0, theta_bar]				  *
 *																			  *
 *   3  = 	normal absolute values w/mean theta_bar							  *
 *				exp[-(x*x)/(sigma*sigma)/2]/[sigma*sqrt(2pi)])				  *
 *																			  *
 * 	4  = 	lorentz use absolute value (sigma)/[(x*x)+(sigma*sigma/4)]/pi	  *
 *																			  *
 * 	5  = 	lorentz2 returns negative values								  *
 *																			  *
 *	6  = 	poisson deviates are all positive, and can get very large		  *
 *																			  *
 * 	7  = 	generalized visiting distribution stariolo and tsallis '96		  *
 *																			  *
 * 	8  = 	standard normal returns negative values							  *
 *																			  *
 * 	9  =  	pareto deviates are all positive, and can get very large		  *
 *   				pareto(x)= {(a*b)**a}/(x**(a+1))						  *
 *					a=value and b=value (a=theta_bar)						  *
 *      		deviate = b / [xi **(1/a)]									  *
 *																			  *
 * 	10 = 	normal returns negative values w/mean theta_bar					  *
*******************************************************************************/


	if (params_dist->distribution == 1)
	{
		/*****exponential distribution, exp(-x/mu) */
		xi= RandomReal();
		theta = (-1) * theta_bar * log(xi); /* exp dist */

		double sign  = RandomReal() - 0.5;
		if (sign <= 0)
				theta = -theta;

	}

	else if ( params_dist->distribution == 2 )
	{
		/* controlled uniform distribution*/
		xi= RandomReal();  /* uniform dist interval =[0, theta_bar] */
		theta = theta_bar * xi; /* controlled uniform */
	}

	else if ( params_dist->distribution == 3 )
	{
		/* normal distribution */
		theta = fabs (theta_bar * gasdev()); /* absolute values of normal dist */
	}

	else if ( params_dist->distribution == 4 )
	{
		/* lorentzian distribution:(sigma) / [(x*x)+(sigma*sigma/4)]/pi */
		/* tan of 0.5*pi is not good */
		do xi= RandomReal();
		while (xi == 0.5);

		/* king's extra /2  theta = fabs (theta_bar/2 * tan(xi*pi) );  */
		theta = fabs (theta_bar * tan(xi*pi) ); /* corrected 5-02 lorentz dist */
	}

	else if ( params_dist->distribution == 5 )
	{
		/* lorentz2 distribution   */
		/* tan of 0.5*pi is not good */
		do xi= RandomReal();
		while (xi == 0.5);  /* do not want tan (0) */

		/* true lorentzian not abs value or round robin index modulo */
		/*  theta = (theta_bar/2 * tan(xi*pi) ); king's extra /2 */
		theta = (theta_bar * tan(xi*pi) ); /* corrected 5-02 lorentz dist */

	}

	else if (params_dist->distribution == 6 )
	{
		/* begin poisson ---needs work still */
		/* poisson deviates are all positive, and can get very large */
		theta = poidev(theta_bar);
	}

	/*******************************************************************/
	else if (params_dist->distribution == 8)
	{
		/* gasdev */
		theta = gasdev(); /* just try standard normal */
	}

	else if (params_dist->distribution == 9)
	{
		/* pareto */
		/* pareto deviates are all positive, and can get very large */
		/*  pareto(x)= {(a*b)**a}/(x**(a+1))  */

		xi= RandomReal(); /*uniform*/
		theta = 1/(pow(xi,(1/theta_bar))); /* pareto with b=1, theta_bar=a */
		/* deviate = b / [xi **(1/a)] */
	}

	else if (params_dist->distribution == 10)
	{
		/* normal returning negative values  */
		theta = theta_bar * gasdev(); /****  normal dist */
	}

	else
	{
		/* unknown distribution, exit */
		printf ("unknown distribution type=%d\n",params_dist->distribution);
		exit (1);
	}

	return theta;

}
