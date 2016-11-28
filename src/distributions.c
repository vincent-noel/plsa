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
/************************************************************
* print_qgt2_visit created for calculation debugging
* is OK to keep around. 1-23-03  LG
************************************************************/

void print_qgt2_visit( double theta_bar)

{  /* begin print_qgt2 */
 /* theta_bar was stariolo's temperature (value 1.)*/


         double fact1,fact3,ov_fact3,fact4,fact5;
	 double a,b,kappa,binom_kappa,area,binom_area;
	 double alpha,krappa,binom_krappa; /* debug 1-13-03 */
	 double star_alpha; /* debug 1-13-03 */
	 double overstar_alpha; /* debug 1-13-03 */
	 double binom_fact4 ,binom_deviate; /* debug 1-13-03 */
	 double deviate,deviate_squared,r1,f,g;
	 double tempf,tempg,tempgc;  /* complicated parameters */


 /*************** begin stariolo rejection method code *****/
  fact1=exp(log(theta_bar)/(3.-DistP.q));     /* dependent on theta_bar */
  fact3=exp(2.*log(theta_bar)/(3.-DistP.q));   /* dependent on theta_bar */
printf ("factor1 = %16.40f  factor3= %16.40f \n ",fact1, fact3 ); /* debug 1-13-03 */
  fact3=fact1*fact1;   /* dependent on theta_bar */
printf ("factor1 = %16.40f  quick fact3 = %16.40f fact1 same as \n s ",fact1,fact3); /* debug 1-13-03 */
printf ("sqrt factor3 = %16.40f \n ", sqrt(fact3) ); /* debug 1-13-03 */
  ov_fact3=1./fact3;                       /* dependent on theta_bar */
printf ("1./factor3 = %16.40f \n ", ov_fact3 ); /* debug 1-13-03 */

  fact5=sqrt(2./fact3);   /* debug 1-13-03 */
printf ("factor5= %16.40f \n",fact5 ); /*debug 1-13-03*/

  a = DistP.alpha / fact1; /* dependent on theta_bar  use alpha 1-17-03 */
printf ("a = %lf  \n",a); /* debug 1-13-03 */
  a = sqrt((DistP.q-1.)/ pi)* DistP.fact2/ fact1; /* dependent on theta_bar */
  b = (DistP.q-1.)/fact3;  /* dependent on theta_bar */
printf ("a = %16.40f b= %16.40f  \n",a, b); /* debug 1-13-03 */
printf ("DistP.trunc %16.40f   \n",DistP.trunc); /* debug 1-13-03 */

  krappa=log(DistP.trunc+sqrt(DistP.trunc*DistP.trunc + fact3/2.)); /* debug 1-13-03 */
  /*  krappa=log(trunc+sqrt(trunc*trunc + 1./(2.*b*DistP.c)));  debug 1-13-03 */
printf ("       krappa = %16.40f \n", krappa ); /* debug 1-13-03 */

  binom_krappa=log(2*DistP.trunc + fact3/(4.*DistP.trunc)); /* debug 1-13-03 */
printf (" binom_krappa = %16.40f \n", binom_krappa ); /* debug 1-13-03 */

/* lucas magic
	 char calcmd[80];
	 char buf[80];
	 double kappa2;
	 FILE *fptr;

	 sprintf(calcmd,"calc 'ln(sqrt(%f**2+%f)-%f)'",trunc,fact3/2.,trunc);
	 printf(calcmd);
	 fptr = popen(calcmd,"r");
	 fgets(buf,79,fptr);
	 printf(buf);
	 kappa2 = atof(buf);
	 fclose(fptr);
	 system(calcmd);
printf ("      kappa2 = %16.40f \n", kappa2 );
*/
  binom_kappa=log(fact3/(4.*DistP.trunc));   /* debug 1-13-03 */
printf ("binom_kappa  = %16.40f   \n", binom_kappa ); /* debug 1-13-03 */

   kappa=binom_kappa;

printf ("  DistP.alpha  = %16.8f \n", DistP.alpha ); /* debug 1-13-03 */
  alpha = (sqrt(2.)* DistP.alpha / 2. );  /* debug 1-13-03 */
printf ("          alpha  = %16.40f should be same \n DistP.alpha2. = %16.40f \n", alpha, DistP.alpha2 ); /* debug 1-13-03 */
  star_alpha = (a/fact5);  /* debug 1-13-03 */
  /*  star_alpha = (a/sqrt(2.*b*DistP.c));   debug 1-13-03 */
  overstar_alpha = (fact5/a);  /* debug 1-13-03 */
printf ("    star_alpha  = %16.40f should be same \n   DistP.alpha2. = %16.40f \n", star_alpha, DistP.alpha2 ); /* debug 1-13-03 */
printf (" overstar_alpha  = %16.40f should be same \n 1/   DistP.alpha2. = %16.40f \n", overstar_alpha, 1./DistP.alpha2 ); /* debug 1-13-03 */

/*  area=(a/sqrt(2./fact3))*(log(trunc+sqrt(trunc*trunc + fact3/2.)) -kappa); */
/*  area=(a/sqrt(2./fact3))*(krappa -kappa); */
  area=(a/fact5)*(krappa -kappa);
printf ("      area    = %16.8f \n", area ); /* debug 1-13-03 */
  binom_area=(alpha)*(krappa -kappa);
printf ("alpha_area    = %16.8f \n", binom_area ); /* debug 1-13-03 */
  binom_area=(a/fact5)*(binom_krappa -binom_kappa);
printf ("binom_area    = %16.8f \n", binom_area ); /* debug 1-13-03 */
  binom_area=(alpha)*(binom_krappa -binom_kappa);
printf ("alpha_area    = %16.8f \n", binom_area ); /* debug 1-13-03 */

/* may want to calc 2 times integral from 0 to trunc rather than -trunc to trunc */

  /*  deviate is a deviate with distribution given by g (q-distribution) */
	   r1=area* RandomReal();
	   fact4=r1*fact5/a+binom_kappa+log(fact5);
printf ("      factor4= %16.30f \n",fact4 ); /*debug 1-13-03*/

	   deviate=(exp(fact4)-exp(-fact4))/(2.*fact5);
printf ("      deviate = %16.8f  \n", deviate ); /* debug 1-13-03 */


	   binom_fact4=r1/alpha+binom_kappa+log(fact5);
printf ("alpha_factor4= %16.30f \n",binom_fact4 ); /* debug 1-13-03 */
	   binom_deviate=(exp(binom_fact4)-exp(-binom_fact4))/(2.*fact5);
printf ("binom_deviate= %16.8f \n", binom_deviate ); /* debug 1-13-03 */

	   deviate_squared=deviate*deviate;
          /***********************************************
           * new calcs for tempf(x) and tempg(x)         *
           ***********************************************/
	   tempf=1.+2.*deviate_squared*ov_fact3 ;/* multiplication is quicker 1-19-03 */
	   tempg=1.+deviate_squared*(DistP.q-1.)*ov_fact3;  /*multiplication is quicker 1-19-03 */
printf ("newer tempf = %16.30f newer tempg= %16.30f \n",tempf, tempg ); /* debug 1-13-03 */

	   tempf=1.+2.*deviate_squared/fact3 ;
	   tempg=1.+deviate_squared*(DistP.q-1.)/fact3;
printf ("   newtempf = %16.30f    newtempg= %16.30f \n",tempf, tempg ); /* debug 1-13-03 */

          /*****************************
           * calculates f(x) and g(x)  *
           *****************************/
	   tempf=1.+2.*b*DistP.c*deviate*deviate;
	   tempg=1.+b*deviate*deviate;
printf ("      tempf = %16.30f       tempg= %16.30f \n",tempf, tempg ); /* debug 1-13-03 */
	   tempgc=pow(tempg,DistP.c);
	   f = a/sqrt(tempf);
	   g = a/(tempgc);
printf ("f = %lf g= %lf \n",f, g ); /* debug 1-13-03 */

     return;  /* return */
} /*end print_qgt2*/

/************************************************************
* binom_qgt2_visit generates random deviates for 2.6 <= q < 3
* for Tsallis GSA q-distribution using the REJECTION METHOD.
* when q is 2.6 or bigger, the tails of the distribution are
* significant, so it becomes necessary to use a larger
* range for x (the variable trunc).  When q >2.85, an even
* larger range for x is needed. This value is set in qgt2_init.
* These large values computationally challenge the computer,
* so the first two terms of the binomial approximation are used
* in the calculation of kappa and krappa.
* This function may be called recursively because a
* deviate value must be returned.  Due to theoretic limit
* the comparison function is the limit of genvis at q=3.
* for values of q between 1 and 2, use qlt2_visit
* for values of q=2 use lorentz, q=1 use normal
* for values of 2< q <2.6 use qgt2_visit
* for values of 2.6<= q < 3 use binom_qgt2_visit (this routine)
* trunc must be set to a larger number for q > 2.8 (see qgt2_init)
*********************************************************/
double binom_qgt2_visit( double theta_bar)

{  /* begin binom_qgt2_visit */
 /* theta_bar was stariolo's temperature (value 1.)*/
        double fact1,fact3,ov_fact3,fact4,fact5;
	double kappa,krappa,area;
	double deviate,deviate_squared,r1,r2,f,g;
	double a      ;   /* NOT dependent on theta_bar should be in init*/
	double tempf,tempg,tempgc;  /* complicated parameters */

  /*************** begin stariolo rejection method code *****/

  fact1=exp(log(theta_bar)/(3.-DistP.q));     /* dependent on theta_bar */
  a = DistP.alpha / fact1; /* dependent on theta_bar  use alpha 1-17-03 */

  fact3=fact1*fact1;                                  /* dependent on theta_bar */
  ov_fact3=1./fact3 ;      /* multiplication is quicker; dependent on theta_bar */

  fact5=sqrt(2.*ov_fact3);  /* multiplication is quicker dependent on theta_bar */


  /*************** calculations made using binomial expansion  *****/
  krappa=log(2.*DistP.trunc + fact3/(4.*DistP.trunc));
  kappa=log(fact3/(4.*DistP.trunc));
  area=(DistP.alpha2)*(krappa -kappa);
  /*************** calculations made using binomial expansion  *****/

  /******************************************************************************************/
  /*  original code needs more efficient calculations above to keep float exceptions at bay */
  /*   NOTE: sqrt(2.*b*DistP.c);    =sqrt(2./fact3);   =sqrt(2)/fact1;                      */
  /*         DistP.alpha = sqrt((DistP.q-1.)/ pi)* DistP.fact2                              */
  /*         fact3=fact1*fact1                                                              */
  /*         ov_fact3 replaces some / which are slow                                        */
  /*         krappa replaces some exp and log which are slow                                */
  /*         deviate_squared=deviate*deviate                                                */
  /*  fact1=exp(log(theta_bar)/(3.-DistP.q));                                               */
  /*  fact3=exp(2.*log(theta_bar)/(3.-DistP.q));                                            */
  /*  a = sqrt((DistP.q-1.)/ pi)* DistP.fact2/ fact1;                                       */
  /*  b = (DistP.q-1.)/fact3;                                                               */
  /*  kappa=log(-trunc+sqrt(trunc*trunc + 1./(2.*b*DistP.c)));                              */
  /*  area=(a/sqrt(2.*b*DistP.c))*(log(trunc+sqrt(trunc*trunc +1./(2.*b*DistP.c))) -kappa); */
  /*  fact5=sqrt(2.*b*DistP.c);    =sqrt(2./fact3);                                         */
  /*  fact4=r1*fact5/a+kappa+log(fact5); fact5/a = DistP.alpha2 1-17-03                     */
  /*  tempf=1.+2.*b*DistP.c*deviate*deviate;  replaced b*c with /fact3                      */
  /*  tempg=1.+b*deviate*deviate;             replaced b with (DistP.q-1.)/fact3            */
  /******************************************************************************************/

  /*  deviate is a deviate with distribution given by g (q-distribution) */
	   r1=area* RandomReal();
	   fact4=(r1/DistP.alpha2)+kappa+log(fact5);
	   deviate=(exp(fact4)-exp(-fact4))/(2.*fact5);
	   deviate_squared=deviate*deviate;

          /*****************************
           * calculates f(x) and g(x)  *
           *****************************/
	   tempf=1.+2.*deviate_squared*ov_fact3 ;/* multiplication is quicker 1-19-03 */
	   tempg=1.+deviate_squared*(DistP.q-1.)*ov_fact3;  /*multiplication is quicker 1-19-03 */
	   tempgc=pow(tempg,DistP.c);
	   f = a/sqrt(tempf);
	   g = a/(tempgc);

           /************************************************************
            * generates a second uniform random number between 0 and f *
            ***********************************************************/

	    r2 = f*RandomReal();

	    /********************************************************************************
	     * recursive rejection code because we need to ALWAYS have a deviate             *
	     * reject when deviate greater than or equal to desired probability distribution *
	     ********************************************************************************/
	    if(r2 >= g)
	      {	/* when reject, try again */

		deviate= binom_qgt2_visit(theta_bar);  /* recursive call*/
		/* DistP.rejects = DistP.rejects +1;     may be used for debugging */

	      } /* endif */

	    return deviate;  /* now always have a deviate to return */

} /*end  binom_qgt2_visit*/


/********************************************************
* qgt2_visit generates random deviates for q between 2 and 2.6
* for Tsallis GSA q-distribution using the REJECTION METHOD.
* This function may be called recursively because a
* deviate value must be returned.  Due to theoretic limit
* the comparison function is the limit of genvis at q=3.
* for values of q between 1 and 2, use qlt2_visit
* for values of q=2 use lorentz, q=1 use normal
* for values of 2< q <2.6 use qgt2_visit (this routine)
* for values of 2.6<= q < 3 use binom_qgt2_visit
********************************************************/
double qgt2_visit( double theta_bar)

{  /* begin qgt2_visit */
 /* theta_bar was stariolo's temperature (value 1.)*/


        double fact1,fact3,ov_fact3,fact4,fact5;
	double a,krappa,kappa,area;
	double deviate,deviate_squared,r1,r2,f,g;
	double tempf,tempg,tempgc;  /* complicated parameters */


  /*************** begin stariolo rejection method code *****/

  fact1=exp(log(theta_bar)/(3.-DistP.q));     /* dependent on theta_bar */

  fact3=fact1*fact1;
  ov_fact3=1./fact3 ;   /*  multiplication is quicker ; dependent on theta_bar */


  a = DistP.alpha / fact1;         /* dependent on theta_bar  */

  kappa=log(-DistP.trunc+sqrt(DistP.trunc*DistP.trunc + fact3/2.));

  krappa=log(DistP.trunc+sqrt(DistP.trunc*DistP.trunc + fact3/2.));


  area=(DistP.alpha2)*(krappa -kappa); /*  more efficient calculation of area */


  fact5=sqrt(2.*ov_fact3);    /* multiplication is quicker 1-19-03 */


  /******************************************************************************************/
  /*  original code needs more efficient calculations above to keep float exceptions at bay */
  /*   NOTE: sqrt(2.*b*DistP.c);    =sqrt(2./fact3);   =sqrt(2)/fact1;                      */
  /*         DistP.alpha = sqrt((DistP.q-1.)/ pi)* DistP.fact2                              */
  /*         fact3=fact1*fact1                                                              */
  /*         ov_fact3 replaces some / which are slow                                        */
  /*         krappa replaces some exp and log which are slow                                */
  /*         deviate_squared=deviate*deviate                                                */
  /*                                                                                        */
  /*  fact1=exp(log(theta_bar)/(3.-DistP.q));                        dependent on theta_bar */
  /*  fact3=exp(2.*log(theta_bar)/(3.-DistP.q));                     dependent on theta_bar */
  /*  a = sqrt((DistP.q-1.)/ pi)* DistP.fact2/ fact1;                dependent on theta_bar */
  /*  b = (DistP.q-1.)/fact3;                                        dependent on theta_bar */
  /*  kappa=log(-DistP.trunc+sqrt(DistP.trunc*DistP.trunc + 1./(2.*b*DistP.c)));            */
  /*  area=(a/sqrt(2.*b*DistP.c))*                                                          */
  /*            (log(DistP.trunc+sqrt(DistP.trunc*DistP.trunc +1./(2.*b*DistP.c))) -kappa); */
  /*  fact5=sqrt(2.*b*DistP.c);   NOTE: sqrt(2.*b*DistP.c)=sqrt(2./fact3)=sqrt(2.)/fact1    */
  /*                                    1./fact5  = fact1* (sqrt(2.)/2)                     */
  /*  fact4=r1*fact5/a+kappa+log(fact5);                                   inefficient way  */
  /*  deviate=(exp(fact4)-exp(-fact4))/(2.*fact5);                                          */
  /*  tempf=1.+2.*b*DistP.c*deviate*deviate;                                                */
  /*  tempg=1.+b*deviate*deviate;                                                           */
  /******************************************************************************************/




  /*  deviate is a deviate with distribution given by g (q-distribution) */
	   r1=area* RandomReal();

	   fact4=(r1/DistP.alpha2)+kappa+log(fact5);      /* fact4 efficient way */
	   deviate=(exp(fact4)-exp(-fact4))/(2.*fact5);

	   deviate_squared=deviate*deviate;

          /*****************************
           * calculates f(x) and g(x)  *
           *****************************/
	   tempf=1.+2.* deviate_squared * ov_fact3 ;
	   tempg=1.+deviate_squared*(DistP.q-1.)*ov_fact3;

	   tempgc=pow(tempg,DistP.c);

	   f = a/sqrt(tempf);
	   g = a/(tempgc);

           /************************************************************
            * generates a second uniform random number between 0 and f *
            ***********************************************************/

	    r2 = f*RandomReal();

	    /*****************************
	     * recursive rejection code
	     *****************************/
            if(r2 >= g)
	      {	/* when we need to reject, must try again */
		deviate=qgt2_visit(theta_bar); /*recursive call*/
		/* DistP.rejects = DistP.rejects +1;   may be used for debugging */
	      } /* endif */
	     /*********************************
	     * now  we ALWAYS have a deviate  *
	     **********************************/
	    return deviate;  /* now always have a deviate to return */

} /*end  qgt2_visit*/

/********************************************************
* qlt2_visit generates random deviates for values of
* q between 1 and 2,  distributed with Tsallis
* GSA q-distribution using the REJECTION METHOD.
* This function may be called recursively because a deviate
* value must be returned.
* This uses a lorentzian bounding function.
* for values of q between 1 and 2, use qlt2_visit (this routine)
* for q between 2 and 3 use qgt2_visit
* for values of q=2 use lorentz, q=1 use normal
* for values of 2< q <2.6 use qgt2_visit
* for values of 2.6<= q < 3 use binom_qgt2_visit
********************************************************/
double qlt2_visit( double theta_bar)

{  /* begin qlt2_visit */

        double fact1,fact3,ov_fact3;
	double a;
	double deviate,deviate_squared,r1,r2,f,g;
	double tempf,tempg,tempgc;  /* complicated parameters */

  /*************** begin stariolo rejection method code *****/

  fact1=exp(log(theta_bar)/(3.-DistP.q));     /* dependent on theta_bar */

  fact3=fact1*fact1;
  ov_fact3=1./fact3 ;   /*  multiplication is quicker 1-19-03; dependent on theta_bar */

  a = DistP.alpha / fact1;        /*  dependent on theta_bar  use alpha 1-17-03 */

  /*  deviate is a deviate with distribution given by f Lorentzian */
	   r1=RandomReal();

           deviate=(tan(pi*r1))*fact1;
	   deviate_squared=deviate*deviate;

          /*****************************
           * calculates f(x) and g(x)  *
           *****************************/

	   tempf=1.+deviate_squared*ov_fact3 ;
	   tempg=1.+deviate_squared*(DistP.q-1.)*ov_fact3;

	   tempgc=pow(tempg,DistP.c);
	   f = a/tempf;
	   g = a/(tempgc);

  /******************************************************************************************/
  /*  original code needs more efficient calculations above to keep float exceptions at bay */
  /*   NOTE: sqrt(2.*b*DistP.c);    =sqrt(2./fact3);   =sqrt(2)/fact1;                      */
  /*         DistP.alpha = sqrt((DistP.q-1.)/ pi)* DistP.fact2                              */
  /*         fact3=fact1*fact1                                                              */
  /*         ov_fact3 replaces some / which are slow                                        */
  /*         krappa replaces some exp and log which are slow                                */
  /*         deviate_squared=deviate*deviate                                                */
  /*                                                                                        */
  /*  fact1=exp(log(theta_bar)/(3.-DistP.q));                        dependent on theta_bar */
  /*  fact3=exp(2.*log(theta_bar)/(3.-DistP.q));                     dependent on theta_bar */
  /*  a = sqrt((DistP.q-1.)/ pi)* DistP.fact2/ fact1;                dependent on theta_bar */
  /*  b = (DistP.q-1.)/fact3;                                        dependent on theta_bar */
  /*  d = sqrt(b*DistP.c);                    remember sqrt(b*c)= sqrt(1./fact3) = 1./fact1 */
  /*  deviate=(tan(pi*r1))/d;                                                  1./d = fact1 */
  /*  tempf=1.+b*DistP.c*deviate*deviate;                                                   */
  /*  tempg=1.+b*deviate*deviate;                                                           */
  /******************************************************************************************/

           /************************************************************
            * generates a second uniform random number between 0 and f *
            ***********************************************************/

	    r2 = f*RandomReal();

	    /*****************************
	     * recursive rejection code
	     *****************************/
            if(r2 >= g)
	    {	/* when we need to reject, must try again */
	      deviate=qlt2_visit(theta_bar); /*recursive call*/
            } /* endif */
	     /*********************************
	     * now  we ALWAYS have a deviate  *
	     **********************************/

	    return deviate;  /* now always have a deviate to return */
} /*end  qlt2_visit*/


/**********************************************************
* qgt2_init initializes the factors used for the  general *
* visiting distribution for q between 2 and 3, but are    *
* not dependent on theta_bar.  These get calculated  only *
* once during initialization of the SA run and remain     *
* during the entire run.                                  *
***********************************************************/

void qgt2_init() /* this routine in distributions.c */

{  /* begin qgt2_init */
  /* these factors are in distributions.h and are static for entire run   */
  DistP.gam1=gammln(1./( DistP.q-1.));
  DistP.gam2=gammln(1./( DistP.q-1.)-0.5);
  DistP.fact2=exp( DistP.gam1 - DistP.gam2);
  DistP.c = 1./( DistP.q-1.);
  DistP.rejects = 0;
  DistP.alpha = sqrt((DistP.q-1.)/ pi)* DistP.fact2; /* NOT dependent on theta_bar */
  DistP.alpha2 = (sqrt(2.)/ 2.)* DistP.alpha; /* NOT dependent on theta_bar */
  /* printf ("DistP.alpha   = %16.8f \n", DistP.alpha );  debug 1-13-03 */
  /* printf ("DistP.alpha2  = %16.8f \n", DistP.alpha2 );  debug 1-13-03 */
  /* printf ("q = %lf factor2= %lf gam1= %lf gam2= %lf \n",DistP.q, DistP.fact2, DistP.gam1, DistP.gam2 ); */

 /*******************************************************************/
 /* trunc is depedent on q, but NOT dependent on theta_bar          */
 /* q= 2.1 10**3   	  q= 2.2 10**3  	     q= 2.3 10**4   */
 /* q= 2.4 10**5  	  q= 2.5 10**6	             q= 2.6 10**9   */
 /* q= 2.6 10**9        for q>2.67 use 10**12                       */
 /* q= 2.7 10**12	for q>2.77 use 10**19                 	    */
 /* q= 2.8 10**19	for q>2.85 use 10**24	  q= 2.9 9*10**38   */
 /*******************************************************************/

  if (DistP.q<2.25)
    { DistP.trunc = 1000.;           /*  2 <  q <2.25 trunc is  10**3 */
    }
  else    /*  q>=2.25  */
    {if (DistP.q<2.35)
      { DistP.trunc = 10000.;       /* 2.25 < q <2.35 trunc is 10**4  */
      }
    else    /*  q>=2.35  */
      {if (DistP.q<2.45)
 	{ DistP.trunc = 100000.;       /* 2.35 < q <2.45 trunc is 10**5  */
 	}
      else    /*  q>=2.45 */
 	{if (DistP.q<2.6)
	  { DistP.trunc = 1000000.;   /* 2.45 < q <2.6 trunc is 10**6  */
	  }
	else   /*  q>=2.6 */
	  {if (DistP.q<2.67)
	    { DistP.trunc = 1000000000.;    /* 2.6 <= q <2.67 trunc is 10**9  use binomial */
	    }
	  else /*  q>=2.67 */
	    {if (DistP.q<2.77)
	      {DistP.trunc = 1000000000000.;  /* 2.67 <= q <2.77 trunc is 10**12 use binomial */
	      }
	    else /*  q>=2.77 */
	      {if (DistP.q<2.85)
		{DistP.trunc = 10000000000000000000.;  /* 2.77 <= q <2.85 trunc is 10**19 use binomial */
		}
	      else /*  q>=2.85 */
		{if (DistP.q<2.9)
		  {DistP.trunc = 1000000000000000000000000.;  /* 2.85 <= q <2.9 trunc is 10**24 use binomial */
		  }
		else
		  {DistP.trunc =99999999999999999999999999999999999999. ;
		  }   /* 2.9 <= q <3 trunc is 38 nines  */

		} /*  q>=2.85 */

	      } /*  q>=2.77 */

	    } /*  q>=2.67 */

	  } /*  q>=2.6 */

	} /*  q>=2.45 */

      } /*  q>=2.35  */

    } /*  q>=2.25  */


/*   printf ("DistP.trunc   = %40.10f \n", DistP.trunc );  debug 1-13-03 */

}   /* end qgt2_init */

/**********************************************************
* qlt2_init initializes the factors used for the  general *
* visiting distribution for q between 1 and 2, but are    *
* not dependent on theta_bar.  These get calculated  only *
* once during initialization of the SA run and remain     *
* during the entire run.                                  *
***********************************************************/

void qlt2_init() /* this routine in distributions.c */

{  /* begin qlt2_init */
  /* these factors are in distributions.h and are static for entire run */
  DistP.gam1=gammln(1./( DistP.q-1.));
  DistP.gam2=gammln(1./( DistP.q-1.)-0.5);
  DistP.fact2=exp( DistP.gam1 - DistP.gam2);
  DistP.c = 1./( DistP.q-1.);
  DistP.alpha = sqrt((DistP.q-1.)/ pi)* DistP.fact2; /* NOT dependent on theta_bar */
  DistP.alpha2 = (sqrt(2.)/ 2.)* DistP.alpha; /* NOT dependent on theta_bar */

  /* printf ("DistP.alpha   = %16.8f \n", DistP.alpha );  debug 1-13-03 */
  /* printf ("DistP.alpha2  = %16.8f \n", DistP.alpha2 );  debug 1-13-03 */
  /* printf ("q = %lf factor2= %lf gam1= %lf gam2= %lf \n",DistP.q, DistP.fact2, DistP.gam1, DistP.gam2 ); */

 }   /* end qlt2_init */

/*******************************************
* gasdev returns a normally distributed
* deviate with zero mean and unit variance
* use RandomReal() rather than ran1
* taken from Numerical Recipes in C p.289
*******************************************/
 double gasdev(void)
{ /* begin gasdev */
     static int iset=0;
     static double gset;
     double fac,rsq,v1,v2;
      if (iset == 0) {
       do {
         v1 = 2.0*RandomReal()-1.0;
         v2 = 2.0*RandomReal()-1.0;
         rsq = v1*v1 + v2*v2;
       } while (rsq >= 1.0 || rsq == 0.0);
       fac=sqrt(-2.0*log(rsq)/rsq);
       gset = v1 * fac;
       iset=1;
       return v2*fac;
     } else {
       iset = 0;
       return gset;
     }
}  /* end gasdev */

/******************************************************
* poidev returns as a float an integer value
* that is a random deviate from a Poisson distribution
* with mean = xm.  Uses RandomReal() rather than ran1
* uses ln because without it too many overflows
* taken from Numerical Recipes in C p.294
*******************************************************/
 double poidev(double xm)
{  /* begin poidev */
   static double sq,alxm,g,oldm=(-1.0);
   double em,t,y;
    if (xm < 12.0) {
     if (xm != oldm) {
         oldm = xm;
         g = exp(-xm);
     }
     em = -1;
     t = 1.0;
     do {
        ++em;
        t *= RandomReal();
     } while (t > g);
   } else {
      if (xm != oldm) {
         oldm = xm;
         sq = sqrt(2.0*xm);
         alxm=log(xm);
         g=xm*alxm-gammln(xm+1.0);
      }
      do {
         do {
            y = tan(pi*RandomReal());
            em = sq*y+xm;
         } while (em < 0.0);
         em = floor(em);
         t = 0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
      } while (RandomReal() > t);
   }
   return em;
}   /* end poidev */

/*********************************************
* gammln returns the ln(gamma(xx)) for xx > 0
* uses ln because without it too many overflows
* taken from Numerical Recipes in C p.214
**********************************************/
double gammln(double xx)
{  /* begin gammln */
   double x,y,tmp,ser;
   static double cof[6] = {76.18009172947146, -86.50532032941677,
                           24.01409824083091, -1.231739572450155,
                           0.1208650973866179e-2, -0.5395239384953e-5};
   int j;
   y=x=xx;
   tmp=x+5.5;
   tmp -= (x+0.5)*log(tmp);
   ser=1.000000000190015;
   for (j=0; j<=5; j++) ser += cof[j]/++y;
   return -tmp+log(2.5066282746310005*ser/x);
}  /* end gammln */


/**************************************************************
* generate_dev generates basic deviates for any distribution  *
* The distribution type and theta_bar are input as parameters *
* deviate is the value returned.                              *
***************************************************************/
double generate_dev( double theta_bar, int distribution, double q)
{  /* begin generate_dev*/

  double xi;  /* uniform random variable for deviates */
  double theta; /* resulting deviate */

/***************************************************************************
****************************************************************************
* 1  = exponential        exp(-x/mu)    ORIGINAL lam's WAY
* 2  = controlled uniform        interval=[0, theta_bar]
* 3  = normal absolute values w/mean theta_bar
       exp[-(x*x)/(sigma*sigma)/2]/[sigma*sqrt(2pi)])
* 4  = lorentz use absolute value (sigma)/[(x*x)+(sigma*sigma/4)]/pi
* 5  = lorentz2 returns negative values
* 6  = poisson
      poisson deviates are all positive, and can get very large
* 7  = generalized visiting distribution stariolo and tsallis '96
* 8  = standard normal returns negative values
* 9  =  pareto deviates are all positive, and can get very large
       pareto(x)= {(a*b)**a}/(x**(a+1)) a=value and b=value (a=theta_bar)
*      deviate = b / [xi **(1/a)]
* 10 = normal returns negative values w/mean theta_bar
*****************************************************************************
*****************************************************************************/


 if (distribution == 1)
/*****************************************************************/
    {  /*****exponential distribution, exp(-x/mu) */
      xi= RandomReal();
      theta = (-1) * theta_bar * log(xi); /* exp dist */

    }   /* end exponential distribution, exp(-x/mu) */
/*****************************************************************/
/********************************************************************/
   else if ( distribution == 2 )
/*******************************************************************/
  { /* controlled uniform distribution*/
      xi= RandomReal();  /* uniform dist interval =[0, theta_bar] */
      theta = theta_bar * xi; /* controlled uniform */
    } /* end controlled uniform distribution  */
/*****************************************************************/
/*****************************************************************/
    else if ( distribution == 3 )
/*****************************************************************/
   { /* normal distribution */
    theta = fabs (theta_bar * gasdev()); /* absolute values of normal dist */
    }  /*end normal distribution */
/*****************************************************************/
   else if ( distribution == 4 )
/*****************************************************************/
   { /* lorentzian distribution:(sigma) / [(x*x)+(sigma*sigma/4)]/pi */
     /* tan of 0.5*pi is not good */
     do
         xi= RandomReal();
     while (xi == 0.5);
     /* king's extra /2  theta = fabs (theta_bar/2 * tan(xi*pi) );  */
     theta = fabs (theta_bar * tan(xi*pi) ); /* corrected 5-02 lorentz dist */
    }  /* end lorentzian distribution:(sigma)/[(x*x)+(sigma*sigma/4)]/pi */
/*****************************************************************/
   else if ( distribution == 5 )
/***********************************************/
   { /* lorentz2 distribution   */
     /* tan of 0.5*pi is not good */
     do
        xi= RandomReal();
     while (xi == 0.5);  /* do not want tan (0) */
   /* true lorentzian not abs value or round robin index modulo */
     /*  theta = (theta_bar/2 * tan(xi*pi) ); king's extra /2 */
     theta = (theta_bar * tan(xi*pi) ); /* corrected 5-02 lorentz dist */

     } /* end lorentz 2 */
/********************************************************/
/************************************************************/
   else if (distribution == 6 )
/************************************************************/
     {  /* begin poisson ---needs work still */
  /* poisson deviates are all positive, and can get very large */
        theta = poidev(theta_bar);
} /* end poisson */
/************************************************************/
/*******************************************************************/
   else if (distribution == 7 )
   {   /* generalized visiting distribution */
     if (q <  2.0) /* when q = 2 should use lorentz */
       {theta = qlt2_visit(theta_bar); /* gsa q between 1 and 2 */
       }
     else /* less than or equal to 2.6 */
       if (q <  2.6) /* do not use binomial */
	 {theta = qgt2_visit(theta_bar); /* gsa q between 2 and 2.6 */}
       else /* q = 2.6 or more, OK to use binomial */
	 {theta = binom_qgt2_visit(theta_bar); /* gsa 2.6 <= q < 3 */ }
   } /* end  generalized visiting distribution */

/*******************************************************************/
 else if (distribution == 8)
    {   /* gasdev */
      theta = gasdev(); /* just try standard normal */
    }   /* end gasdev returns pos and neg values*/
/*******************************************************************/
/*******************************************************************/

 else if (distribution == 9)
    {   /* pareto */
  /* pareto deviates are all positive, and can get very large */
  /*  pareto(x)= {(a*b)**a}/(x**(a+1))  */

       xi= RandomReal(); /*uniform*/
      theta = 1/(pow(xi,(1/theta_bar))); /* pareto with b=1, theta_bar=a */
      /* deviate = b / [xi **(1/a)] */
   }   /* end pareto */
/*******************************************************************/

 else if (distribution == 10)
    {   /* normal returning negative values  */
      theta = theta_bar * gasdev(); /****  normal dist */
    }   /* end normal returning negative values */

   else
/*******************************************************************/
      { /* unknown distribution, exit */
        printf ("unknown distribution type=%d\n",distribution);
        exit (1);
      } /* unknown distribution */

 return theta;

/*******************************************************************/
 }  /* end generate_dev*/
