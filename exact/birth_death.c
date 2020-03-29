/*
birth_death.c (v7)

Code for exact MCMC inference in the birth-death model

Based on an inhomogeneous Poisson process proposal
with rate changing linearly from the initial rate to the
final rate for each time interval

Births updated using a difference between two Poissons.  
Mean depends on current state, so PMF (which involves a 
Bessel function) is part of the acceptance probability


Darren Wilkinson
d.j.wilkinson@ncl.ac.uk
http://www.staff.ncl.ac.uk/d.j.wilkinson/

Last updated: 8/5/2002
*/

#include "interval.h"

/* function prototypes */
void read_data(void);
void init_in(gsl_vector_int *type_vec);
void output_header(void);
void output_state(long it);
void update_c(double alpha,double beta);
void update_in(interval *in_prop,gsl_vector_int *type_vec);
int accept_lp(double lp);
double log_rnd(interval *invl,double xm,int x0,int x1);
double log_pi(int x,double mean);
double int_in(interval *invl,int x0);
double log_sum_in(interval *invl,int x0,int x1);
double log_pi_poisd(int y,double theta);
void transform_in(interval *invl,int x0,int x1);


/* global variables */
gsl_vector *c;           /* parameter vector */
interval **in;           /* array of pointers to interval structures */
int n_in;                /* number of intervals */
gsl_vector_int *data;    /* data vector */
gsl_rng *rng;            /* random number stream */

/* main function */
int main(int argc,char *argv[])
{
  /* main variables */
  long it,it_max;
  int nobs,i;
  double alpha,beta;
  interval *in_prop;
  gsl_vector_int *type_vec;
  
  /* process command line */
  if (argc != 3) {
    fprintf(stderr,"Usage: %s <iters> <num obs>\n",argv[0]);
    exit(1);
  }
  it_max=atoi(argv[1]);
  nobs=atoi(argv[2]);
  n_in=nobs-1;

  /* prior */
  alpha=1.0;
  beta=0.01;

  /* inits */
  rng=gsl_rng_alloc(gsl_rng_mt19937);
  c=gsl_vector_alloc(2);
  for (i=0;i<(c->size);i++) {
    gsl_vector_set(c,i,1.0);
  }
  in_prop=interval_alloc(c->size,1.0);
  type_vec=gsl_vector_int_calloc(c->size);
  data=gsl_vector_int_calloc(nobs);
  read_data();
  init_in(type_vec);  

  /* main MCMC loop */
  output_header();
  for (it=0;it<it_max;it++) {
    update_c(alpha,beta);
    update_in(in_prop,type_vec);
    output_state(it);
  }
  /* end */
  return(0);
}


/* helper functions */


void update_c(double alpha,double beta)
{
  int i,births,deaths;
  double lambda,mu,in_x;
  births=0;
  deaths=0;
  in_x=0.0;
  for (i=0;i<n_in;i++) {
    births+=gsl_vector_int_get(in[i]->nv,0);
    deaths+=gsl_vector_int_get(in[i]->nv,1);
    in_x+=int_in(in[i],gsl_vector_int_get(data,i));
  }
#ifdef DEBUG
  printf("B: %d, D: %d, Int: %f\n",births,deaths,in_x);
#endif
  lambda=gsl_ran_gamma(rng,alpha+births,1.0/(beta+in_x));
  mu=gsl_ran_gamma(rng,alpha+deaths,1.0/(beta+in_x));
  gsl_vector_set(c,0,lambda);
  gsl_vector_set(c,1,mu);
}

void update_in(interval *in_prop,gsl_vector_int *type_vec)
{
  int i,births,deaths,old_births,old_deaths,x0,x1,y;
  double a,xm,theta;
  for (i=0;i<n_in;i++) {
    x0=gsl_vector_int_get(data,i);
    x1=gsl_vector_int_get(data,i+1);
    old_births=gsl_vector_int_get(in[i]->nv,0);
    theta=(double) old_births*old_births/200.0 + 1.0;
    y=gsl_ran_poisson(rng,theta)-gsl_ran_poisson(rng,theta);
    births=old_births+y;
    deaths=births+x0-x1;
    if ( (births >= 0) && (deaths >= 0) ) {
      xm=0.5*(x0+x1);
      /* first construct proposed interval */
      gsl_vector_int_set(type_vec,0,births);
      gsl_vector_int_set(type_vec,1,deaths);
      interval_set(rng,in_prop,type_vec);
      transform_in(in_prop,x0,x1);   /* make inhomogeneous */
      /* now evaluate log acceptance probability */
      old_deaths=gsl_vector_int_get(in[i]->nv,1);
      a=log_rnd(in_prop,xm,x0,x1) - log_rnd(in[i],xm,x0,x1)
	-log_pi_poisd(y,theta)
	+log_pi(births,xm*gsl_vector_get(c,0))
	+log_pi(deaths,xm*gsl_vector_get(c,1))
	+log_pi_poisd(y,(double) births*births/200.0 + 1.0)
	-log_pi(old_births,xm*gsl_vector_get(c,0))
	-log_pi(old_deaths,xm*gsl_vector_get(c,1));
#ifdef DEBUG
      printf("i: %d, ob: %d, od: %d, b: %d, d: %d, a: %f. ",i,old_births,old_deaths,births,deaths,a);
      printf("lrnd: %f, olrnd: %f. ",log_rnd(in_prop,xm,x0,x1),log_rnd(in[i],xm,x0,x1));
#endif
      /* test and update */
      if (accept_lp(a)) {
	interval_cpy(in[i],in_prop);
#ifdef DEBUG
	printf("Accept. ");
#endif
      }
    }
  }
#ifdef DEBUG
  printf("\n");
#endif
}

void transform_in(interval *invl,int x0,int x1)
{
  int j;
  double t,tp;
  if (x0 != x1) {
    for (j=0;j<(invl->n);j++) {
      tp=invl->time[j];
      t=( sqrt( x0*x0 + (x1*x1 - x0*x0)*tp ) - x0 )/(x1-x0);
      invl->time[j] = t;
    }
  }
}

double log_rnd(interval *invl,double xm,int x0,int x1)
{
  double lambda,mu,in_x,log_sum;
  lambda=gsl_vector_get(c,0);
  mu=gsl_vector_get(c,1);
  in_x=int_in(invl,x0);
  log_sum=log_sum_in(invl,x0,x1);
  return( (lambda+mu)*(xm-in_x) + log_sum );
}

double int_in(interval *invl,int x0)
{
  int i,type;
  double sum;
  sum=x0;
  for (i=0;i<(invl->n);i++) {
    type=invl->type[i];
    if (type == 0) {
      sum+=(1.0-(invl->time[i]));
    } else {
      sum-=(1.0-(invl->time[i]));
    }
  }
  return(sum);
}

double int_in_old(interval *invl,int x0)
{
  int i,x,type;
  double sum;
  x=x0;
  sum=x*(invl->time[0]);
  for (i=0;i<(invl->n)-1;i++) {
    type=invl->type[i];
    if (type == 0) {x++;} else {x--;}
    sum+=(invl->time[i+1] - invl->time[i])*x;
  }
  type=invl->type[invl->n];
  if (type == 0) {x++;} else {x--;}
  sum+=(1.0 - invl->time[(invl->n)-1])*x;
  return(sum);
}

double log_sum_in(interval *invl,int x0,int x1)
{
  int i,x,type;
  double sum,t;
  sum=0.0;
  t=0.0;
  x=x0;
  for (i=0;i<(invl->n);i++) {
    /* first update sum */
    sum += log( x );
    sum -= log( (1-t)*x0 + t*x1 ) ;
    /* now update state */
    t=invl->time[i];
    type=invl->type[i];
    if (type == 0) {x++;} else {x--;}
  }
  return(sum);
}

double log_pi(int x,double mean)
{
  return( -mean + x*log(mean) - gsl_sf_lngamma((double) x+1) );
}

double log_pi_poisd(int y,double theta)
{
  return( log(gsl_sf_bessel_In_scaled(y,2.0*theta)) );
}

int accept_lp(double lp)
{
  return(log(gsl_ran_flat(rng,0,1)) < lp);
}

void output_state(long it)
{
  int i;
  printf("%ld ",it);
  printf("%f ",gsl_vector_get(c,0));
  printf("%f ",gsl_vector_get(c,1));
  printf("%f ",gsl_vector_get(c,0)+gsl_vector_get(c,1));
  printf("%f ",gsl_vector_get(c,0)-gsl_vector_get(c,1));
  for (i=0;i<n_in;i++) {
    printf("%d ",gsl_vector_int_get(in[i]->nv,0));
    printf("%d ",gsl_vector_int_get(in[i]->nv,1));
  }
  printf("\n");
#ifdef DEBUG2
  for (i=0;i<n_in;i++) {
    interval_print(in[i]);
  }
#endif
}

void output_header(void)
{
  int i;
  printf("Iter lambda mu sum diff ");
  for (i=0;i<n_in;i++) {
    printf("b[%d] ",i);
    printf("d[%d] ",i);
  }
  printf("\n");
}

void init_in(gsl_vector_int *type_vec)
{
  int i,x,nx;
  in=malloc(n_in*sizeof(interval *));
  for (i=0;i<n_in;i++) {
    in[i]=interval_alloc(2,1.0);
    x=gsl_vector_int_get(data,i);
    nx=gsl_vector_int_get(data,i+1);
    if (nx>x) {
      gsl_vector_int_set(type_vec,0,11*(nx-x)+15);
      gsl_vector_int_set(type_vec,1,10*(nx-x)+15);
    }
    else {
      gsl_vector_int_set(type_vec,0,10*(x-nx)+15);
      gsl_vector_int_set(type_vec,1,11*(x-nx)+15);
    }
    interval_set(rng,in[i],type_vec);
  }
}

void read_data(void)
{
  FILE *s;
  int i,temp;
  s=fopen("bd.dat","r");
  if (s==NULL) {
    perror("error opening bd.dat");
    exit(EXIT_FAILURE);
  }
  for (i=0;i<(data->size);i++) {
    fscanf(s,"%d",&temp);
    gsl_vector_int_set(data,i,temp);
  }
  fclose(s);
}


/* eof */


