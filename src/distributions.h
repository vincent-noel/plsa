/*****************************************************************************
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
*   Copyright (C) 2016 Vincent Noel                                          *
*   the full GPL copyright notice can be found in lsa.c                      *
*                                                                            *
******************************************************************************/
typedef struct my_distrib {
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

  /*****variables that depend on theta_bar ********/



} DistParms;
 DistParms DistP;   /* variables for distributions */

/* prototype functions for distributions */

double gasdev(void);  /* normal distribution */
double gammln(double xx);  /* gamma ln distribution */
double poidev(double xm);/* poisson distribution */

double binom_qgt2_visit(double theta_bar); /* 2.6 <= q < 3 */
double qgt2_visit(double theta_bar);       /* 2 < q < 2.6 */
double qlt2_visit(double theta_bar);       /* q between 1 and 2 */

void print_qgt2_visit( double theta_bar); /* debug 1-13-03 */

double generate_dev(double theta_bar, int distribution, double q);

/**********************************************************
* qgt2_init initializes the factors used for the  general *
* visiting distribution for q between 2 and 3, but are    *
* not dependent on theta_bar.  These get calculated  only *
* once during initialization of the SA run and remain     *
* during the entire run.                                  *
***********************************************************/

void qgt2_init(void); /* this routine in distributions.c */

/**********************************************************
* qlt2_init initializes the factors used for the  general *
* visiting distribution for q between 1 and 2, but are    *
* not dependent on theta_bar.  These get calculated  only *
* once during initialization of the SA run and remain     *
* during the entire run.                                  *
***********************************************************/
void qlt2_init(void); /* this routine in distributions.c */
