/*
lotka_volterra.c (v2)

Code for approximate MCMC inference in the lotka-volterra model

Based on a Poisson process approximation; rates are the average
rates for the interval

Births updated using a difference between two Poissons.  
Mean depends on current state, so PMF (which involves a 
Bessel function) is part of the acceptance probability


Darren Wilkinson
d.j.wilkinson@ncl.ac.uk
http://www.staff.ncl.ac.uk/d.j.wilkinson/

Last updated: 7/9/2002
*/

#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sf.h>


/* function prototypes */
gsl_matrix_int * read_data(char *filename);
void init_in(void);
void output_header(void);
void output_state(long it);
void update_c(double alpha,double beta);
void update_in(void);
int accept_lp(double lp);
double log_pi(int x,double mean);
double log_pi_poisd(int y,double theta);

/* macros */
#define PARAM 3
#define DATA 2

/* global variables */
gsl_vector *c;           /* parameter vector */
gsl_matrix_int *r;       /* reaction matrix */
int n_in;                /* number of intervals */
gsl_matrix_int *data;    /* data matrix */
gsl_rng *rng;            /* random number stream */

/* main function */
int main(int argc,char *argv[])
{
  /* main variables */
  long it,it_max;
  int nobs,i;
  double alpha,beta;
  
  /* process command line */
  if (argc != 3) {
    fprintf(stderr,"Usage: %s <iters> <num obs>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  it_max=atoi(argv[1]);
  nobs=atoi(argv[2]);
  n_in=nobs-1;

  /* prior */
  alpha=1.0;
  beta=0.01;

  /* inits */
  rng=gsl_rng_alloc(gsl_rng_mt19937);
  c=gsl_vector_alloc(PARAM);
  for (i=0;i<(c->size);i++) {
    gsl_vector_set(c,i,1.0);
  }
  data=read_data("lv.dat");
  init_in();  

  /* main MCMC loop */
  output_header();
  for (it=0;it<it_max;it++) {
    update_c(alpha,beta);
    update_in();
    output_state(it);
  }
  /* end */
  return(0);
}


/* helper functions */


void update_c(double alpha,double beta)
{
  int i,r0,r1,r2;
  double c0,c1,c2,in0,in1,in2;
  r0=0; r1=0; r2=0;
  in0=0.0; in1=0.0; in2=0.0;
  for (i=0;i<n_in;i++) {
    r0+=gsl_matrix_int_get(r,i,0);
    r1+=gsl_matrix_int_get(r,i,1);
    r2+=gsl_matrix_int_get(r,i,2);
    in0+=0.5*(gsl_matrix_int_get(data,i,0)
	      + gsl_matrix_int_get(data,i+1,0) );
    in1+=0.5*( gsl_matrix_int_get(data,i,0)*gsl_matrix_int_get(data,i,1)
	       + gsl_matrix_int_get(data,i+1,0)*gsl_matrix_int_get(data,i+1,1)
	       );
    in2+=0.5*(gsl_matrix_int_get(data,i,1)
	      + gsl_matrix_int_get(data,i+1,1) );
  }
  c0=gsl_ran_gamma(rng,alpha+r0,1.0/(beta+in0));
  c1=gsl_ran_gamma(rng,alpha+r1,1.0/(beta+in1));
  c2=gsl_ran_gamma(rng,alpha+r2,1.0/(beta+in2));
  gsl_vector_set(c,0,c0);
  gsl_vector_set(c,1,c1);
  gsl_vector_set(c,2,c2);
}

void update_in(void)
{
  int i,r0,r1,r2,r0o,r1o,r2o,y0,y1,y0n,y1n,y;
  double a,theta,lam0,lam1,lam2,lam0n,lam1n,lam2n,c0,c1,c2;
  for (i=0;i<n_in;i++) {
    y0=gsl_matrix_int_get(data,i,0);
    y1=gsl_matrix_int_get(data,i,1);
    y0n=gsl_matrix_int_get(data,i+1,0);
    y1n=gsl_matrix_int_get(data,i+1,1);
    r0o=gsl_matrix_int_get(r,i,0);
    r1o=gsl_matrix_int_get(r,i,1);
    r2o=gsl_matrix_int_get(r,i,2);
    theta=(double) r0o*r0o/200.0 + 1.0;
    y=gsl_ran_poisson(rng,theta)-gsl_ran_poisson(rng,theta);
    r0=r0o+y;
    r1=y0-y0n+r0;
    r2=y0-y0n+y1-y1n+r0;
    if ( (r0 >= 0) && (r1 >= 0) && (r2 >= 0) ) {
      /* first construct proposed interval */
      c0=gsl_vector_get(c,0);
      c1=gsl_vector_get(c,1);
      c2=gsl_vector_get(c,2);
      lam0=c0*y0;
      lam1=c1*y0*y1;
      lam2=c2*y1;
      lam0n=c0*y0n;
      lam1n=c1*y0n*y1n;
      lam2n=c2*y1n;
      /* now evaluate log acceptance probability */
      a=-log_pi_poisd(y,theta)
	+log_pi(r0,c0*(y0+y0n)/2.0)
	+log_pi(r1,c1*(y0*y1+y0n*y1n)/2.0)
	+log_pi(r2,c2*(y1+y1n)/2.0)
	+log_pi_poisd(y,(double) r0*r0/200.0 + 1.0)
	-log_pi(r0o,c0*(y0+y0n)/2.0)
	-log_pi(r1o,c1*(y0*y1+y0n*y1n)/2.0)
	-log_pi(r2o,c2*(y1+y1n)/2.0);
#ifdef DEBUG
      printf("[a:%f ",a);
      printf("] ");
#endif
      /* test and update */
      if (accept_lp(a)) {
	gsl_matrix_int_set(r,i,0,r0);
	gsl_matrix_int_set(r,i,1,r1);
	gsl_matrix_int_set(r,i,2,r2);
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
  int i,j;
  printf("%ld ",it);
  for (i=0;i<(c->size);i++) {
    printf("%f ",gsl_vector_get(c,i));
  }
  for (i=0;i<n_in;i++) {
    for (j=0;j<1;j++) {
      printf("%d ",gsl_matrix_int_get(r,i,j));
    }
  }
  printf("\n");
}

void output_header(void)
{
  int i,j;
  printf("Iter ");
  for (i=0;i<(c->size);i++) {
    printf("c%d ",i);
  }
  for (i=0;i<n_in;i++) {
    for (j=0;j<1;j++) {
      printf("r[%d,%d] ",i,j);
    }
  }
  printf("\n");
}

void init_in(void)
{
  /* re-write this in a more generic way! */
  int i,r0,r1,r2,y0,y0n,y1,y1n;
  r=gsl_matrix_int_calloc(n_in+1,PARAM);
  for (i=0;i<n_in;i++) {
    y0=gsl_matrix_int_get(data,i,0);
    y0n=gsl_matrix_int_get(data,i+1,0);
    y1=gsl_matrix_int_get(data,i,1);
    y1n=gsl_matrix_int_get(data,i+1,1);
    r0=20; /* do something better here! */
    r0=5;
    if (y0n>y0) { r0+=y0n-y0; }
    if (y1n>y1) { r0+=y1n-y1; }
    r1=y0-y0n+r0;
    r2=y0-y0n+y1-y1n+r0;
    if ( (r0<0)||(r1<0)||(r2<0) ) {
      perror("Initialisation failed");
      exit(EXIT_FAILURE);
    }
    gsl_matrix_int_set(r,i,0,r0);
    gsl_matrix_int_set(r,i,1,r1);
    gsl_matrix_int_set(r,i,2,r2);
  }
}

gsl_matrix_int * read_data(char *filename)
{
  FILE *s;
  int i,j,tmp;
  double temp;
  data=gsl_matrix_int_calloc(n_in + 1,DATA);
  s=fopen(filename,"r");
  if (s==NULL) {
    perror("error opening data file");
    exit(EXIT_FAILURE);
  }
  for (i=0;i<=n_in;i++) {
    for (j=0;j<DATA;j++) {
      fscanf(s,"%lf",&temp);
      tmp=(int) temp;
      gsl_matrix_int_set(data,i,j,tmp);
    }
  }
  fclose(s);
  return(data);
}


/* eof */


