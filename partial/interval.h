/*
interval.h

Header for time interval structures

Darren Wilkinson
d.j.wilkinson@ncl.ac.uk
http://www.staff.ncl.ac.uk/d.j.wilkinson/

Last updated: 27/2/2002
*/

#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

/* definition of the "interval" structure */
typedef struct
{
  int n;                   /* total # of reactions */
  int n_max;               /* max allowed total without a realloc */
  double len;              /* length of time interval */
  int n_types;             /* number of reaction types */
  gsl_vector_int *nv;      /* numbers of each type */
  gsl_vector_int *nv_temp; /* temp vec used during simulation */
  double *time;            /* array of reaction times */
  int *type;               /* array of corresponding types */
} interval;


/* function prototypes */

/* external interface */
interval * interval_alloc(int n_types,double length);
void interval_set(gsl_rng *r,interval *in,gsl_vector_int *nv);
void interval_cpy(interval *dest,interval *src);
void interval_print(interval *in);

/* private functions */
int vector_int_sum(gsl_vector_int *nv);
double gen_min_uni(gsl_rng *r,int n,double a);
int gen_type(gsl_rng *r,gsl_vector_int *nv);







/* eof */

