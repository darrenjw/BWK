/*
lotka_volterra.c (v3)

Code for exact MCMC inference in the lotka-volterra model

Based on an inhomogeneous Poisson process proposal
with rate changing linearly from the initial rate to the
final rate for each time interval

Births updated using a difference between two Poissons.  
Mean depends on current state, so PMF (which involves a 
Bessel function) is part of the acceptance probability

FIXED BUG IN log_sum_in!!!

ADDING code for sticky intervals...

Darren Wilkinson
d.j.wilkinson@ncl.ac.uk
http://www.staff.ncl.ac.uk/d.j.wilkinson/

Last updated: 23/4/2007
*/

#include "interval.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_sort_vector.h>

/* function prototypes */
gsl_matrix_int * read_data(char *filename);
void init_in(gsl_vector_int *type_vec);
void output_header(void);
void output_state(long it);
void update_c(double alpha,double beta);
void update_in(interval *in_prop,gsl_vector_int *type_vec);
int accept_lp(double lp);
double log_rnd(interval *invl,int y0,int y1,int y0n,int y1n);
double log_pi(int x,double mean);
double int_in0(interval *invl,int x0);
double int_in1(interval *invl,int y0,int y1);
double int_in2(interval *invl,int x0);
double int_in(interval *invl,int y0,int y1,double c0,double c1,double c2);
double log_sum_in(interval *invl,int y0,int y1,int y0n,int y1n,double c0,double c1,double c2);
double log_pi_poisd(int y,double theta);
void transform_in(interval *invl,double,double,double,double,double,double);
int interval_ok(interval *invl,int y0,int y1);

/* macros */
#define PARAM 3
#define DATA 2

/* global variables */
gsl_vector *c;           /* parameter vector */
interval **in;           /* array of pointers to interval structures */
int n_in;                /* number of intervals */
gsl_matrix_int *data;    /* data vector */
gsl_vector_int *sticky;  /* count since last acceptance */
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
  gsl_rng_set(rng,1234); /* non-default seed for -e2- */
  c=gsl_vector_alloc(PARAM);
  for (i=0;i<(c->size);i++) {
    gsl_vector_set(c,i,1.0);
  }

  gsl_vector_set(c,0,0.25);
  gsl_vector_set(c,1,0.0025);
  gsl_vector_set(c,2,0.125);

  in_prop=interval_alloc(c->size,1.0);
  type_vec=gsl_vector_int_calloc(c->size);
  data=read_data("lv.dat");
  init_in(type_vec);  

  sticky=gsl_vector_int_calloc(n_in);

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
  int i,r0,r1,r2;
  double c0,c1,c2,in0,in1,in2;
  r0=0; r1=0; r2=0;
  in0=0.0; in1=0.0; in2=0.0;
  for (i=0;i<n_in;i++) {
    r0+=gsl_vector_int_get(in[i]->nv,0);
    r1+=gsl_vector_int_get(in[i]->nv,1);
    r2+=gsl_vector_int_get(in[i]->nv,2);
    in0+=int_in0(in[i],gsl_matrix_int_get(data,i,0));
    in1+=int_in1(in[i],gsl_matrix_int_get(data,i,0)
		 ,gsl_matrix_int_get(data,i,1));
    in2+=int_in2(in[i],gsl_matrix_int_get(data,i,1));
  }
  c0=gsl_ran_gamma(rng,alpha+r0,1.0/(beta+in0));
  c1=gsl_ran_gamma(rng,alpha+r1,1.0/(beta+in1));
  c2=gsl_ran_gamma(rng,alpha+r2,1.0/(beta+in2));
  gsl_vector_set(c,0,c0);
  gsl_vector_set(c,1,c1);
  gsl_vector_set(c,2,c2);
}

void update_in(interval *in_prop,gsl_vector_int *type_vec)
{
  int i,r0,r1,r2,r0o,r1o,r2o,y0,y1,y0n,y1n,y;
  double a,theta,lam0,lam1,lam2,lam0n,lam1n,lam2n,c0,c1,c2;
  for (i=0;i<n_in;i++) {
    gsl_vector_int_set(sticky,i,gsl_vector_int_get(sticky,i)+1);
    y0=gsl_matrix_int_get(data,i,0);
    y1=gsl_matrix_int_get(data,i,1);
    y0n=gsl_matrix_int_get(data,i+1,0);
    y1n=gsl_matrix_int_get(data,i+1,1);
    r0o=gsl_vector_int_get(in[i]->nv,0);
    r1o=gsl_vector_int_get(in[i]->nv,1);
    r2o=gsl_vector_int_get(in[i]->nv,2);
    theta=(double) r0o*r0o/200.0 + 1.0;
    y=gsl_ran_poisson(rng,theta)-gsl_ran_poisson(rng,theta);
    r0=r0o+y;
    r1=y0-y0n+r0;
    r2=y0-y0n+y1-y1n+r0;
    if ( (r0 >= 0) && (r1 >= 0) && (r2 >= 0) ) {
      /* first construct proposed interval */
      gsl_vector_int_set(type_vec,0,r0);
      gsl_vector_int_set(type_vec,1,r1);
      gsl_vector_int_set(type_vec,2,r2);
      c0=gsl_vector_get(c,0);
      c1=gsl_vector_get(c,1);
      c2=gsl_vector_get(c,2);
      lam0=c0*y0;
      lam1=c1*y0*y1;
      lam2=c2*y1;
      lam0n=c0*y0n;
      lam1n=c1*y0n*y1n;
      lam2n=c2*y1n;
      interval_set(rng,in_prop,type_vec);
      /* now make it inhomogeneous */
      transform_in(in_prop,lam0,lam1,lam2,lam0n,lam1n,lam2n);
      /* if interval is bad, don't bother going any further */
      if (interval_ok(in_prop,y0,y1)) {
	/* now evaluate log acceptance probability */
	a=log_rnd(in_prop,y0,y1,y0n,y1n) - log_rnd(in[i],y0,y1,y0n,y1n)
	  -log_pi_poisd(y,theta)
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
	if ( accept_lp(a) || (gsl_vector_int_get(sticky,i)>1000) ) {
	  interval_cpy(in[i],in_prop);
	  gsl_vector_int_set(sticky,i,0);
#ifdef DEBUG
	  printf("Accept. ");
#endif
	}
      }
    }
  }
#ifdef DEBUG
  printf("\n");
#endif
}

void transform_in(interval *invl,double l0,double l1,double l2,double l0n,double l1n,double l2n)
{
  int j,type;
  double t,tp;
  gsl_vector_view timeView;
  gsl_vector_int_view  typeView;
  gsl_permutation *perm;
  for (j=0;j<(invl->n);j++) {
    tp=invl->time[j];
    t=tp;
    type=invl->type[j];
    if ((type == 0)&&(l0 != l0n)) {
      t=( sqrt(l0*l0 + (l0n*l0n - l0*l0)*tp) - l0)/(l0n-l0);
    }
    else if ((type == 1)&&(l1 != l1n)) {
      t=( sqrt(l1*l1 + (l1n*l1n - l1*l1)*tp) - l1)/(l1n-l1);
    }
    else if ((type == 2)&&(l2 != l2n)) {
      t=( sqrt(l2*l2 + (l2n*l2n - l2*l2)*tp) - l2)/(l2n-l2);
    }
    invl->time[j] = t;
  }
  /* now re-sort times */
  if (invl->n > 0) {
    perm=gsl_permutation_alloc(invl->n);
    if (!perm) {
      perror("failed to alloc perm");
      exit(EXIT_FAILURE);
    }
    timeView=gsl_vector_view_array(invl->time,invl->n);
    typeView=gsl_vector_int_view_array(invl->type,invl->n);
    gsl_sort_vector_index(perm,&(timeView.vector));
    gsl_permute_vector(perm,&(timeView.vector));
    gsl_permute_vector_int(perm,&(typeView.vector));
    gsl_permutation_free(perm);
  }
}

int interval_ok(interval *invl,int y0,int y1)
{
  int j,ok,type;
  ok=1;
  for (j=0;j<(invl->n);j++) {
    if ((y0<=0)||(y1<=0)) { ok=0; }
    type=(invl->type)[j];
    if (type==0) {y0++;}
    else if (type==1) {y0--; y1++; }
    else {y1--;}
  }
  return(ok);
}

double log_rnd(interval *invl,int y0,int y1,int y0n,int y1n)
{
  double c0,c1,c2,in_x,log_sum,a0,a1;
  c0=gsl_vector_get(c,0);
  c1=gsl_vector_get(c,1);
  c2=gsl_vector_get(c,2);
  in_x=int_in(invl,y0,y1,c0,c1,c2);
  log_sum=log_sum_in(invl,y0,y1,y0n,y1n,c0,c1,c2);
  a0=c0*y0 + c1*y0*y1 + c2*y1;
  a1=c0*y0n + c1*y0n*y1n + c2*y1n;
#ifdef DEBUG
printf("{inx:%f, ls:%f, a0:%f, a1:%f}",in_x,log_sum,a0,a1);
#endif
  return( (a0+a1)/2.0 - in_x + log_sum );
}

double int_in0(interval *invl,int x0)
{
  int i,type;
  double sum;
  sum=x0;
  for (i=0;i<(invl->n);i++) {
    type=invl->type[i];
    if (type == 0) {
      sum+=(1.0-(invl->time[i]));
    }
    if (type == 1) {
      sum-=(1.0-(invl->time[i]));
    }
  }
  return(sum);
}

double int_in1(interval *invl,int y0,int y1)
{
  int i,type;
  double sum;
  sum=y0*y1;
  for (i=0;i<(invl->n);i++) {
    type=invl->type[i];
    if (type == 0) {
      sum+=y1*(1.0-(invl->time[i]));
      y0++;
    }
    if (type == 1) {
      sum+=(y0-y1-1)*(1.0-(invl->time[i]));
      y0--; y1++;
    }
    if (type == 2) {
      sum-=y0*(1.0-(invl->time[i]));
      y1--;
    }
  }
  return(sum);
}

double int_in2(interval *invl,int x0)
{
  int i,type;
  double sum;
  sum=x0;
  for (i=0;i<(invl->n);i++) {
    type=invl->type[i];
    if (type == 1) {
      sum+=(1.0-(invl->time[i]));
    }
    if (type == 2) {
      sum-=(1.0-(invl->time[i]));
    }
  }
  return(sum);
}

double int_in(interval *invl,int y0,int y1,double c0,double c1,double c2)
{
  return( c0*int_in0(invl,y0) + c1*int_in1(invl,y0,y1) + c2*int_in2(invl,y1) );
}

double log_sum_in(interval *invl,int y0,int y1,int y0n,int y1n,double c0,double c1,double c2)
{
  int i,type,y0c,y1c;
  double sum,t,a,a0,a1;
  sum=0.0;
  t=0.0;
  y0c=y0; y1c=y1;
  for (i=0;i<(invl->n);i++) {
    type=invl->type[i];
    /* first update sum */
    if (type==0) {
      a=c0*y0c;
      a0=c0*y0;
      a1=c0*y0n;
    } 
    else if (type == 1) {
      a=c1*y0c*y1c;
      a0=c1*y0*y1;
      a1=c1*y0n*y1n;
    }
    else {
      a=c2*y1c;
      a0=c2*y1;
      a1=c2*y1n;
    }
    if (a<=0) {
      perror("Bad interval");
      fprintf(stderr,"a=%f, type=%d, c0=%f,c1=%f,c2=%f, y0c=%d,y1c=%d, i=%d, y0=%d,y1=%d,y0n=%d,y1n=%d\n",a,type,c0,c1,c2,y0c,y1c,i,y0,y1,y0n,y1n);
      exit(EXIT_FAILURE);
    }
    sum += log( a );
    t=invl->time[i];
    sum -= log( (1-t)*a0 + t*a1 ) ;
    /* now update state */
    if (type == 0)
      {y0c++;} 
    else if (type == 1)
      {y0c--; y1c++;}
    else if (type == 2)
      {y1c--;}
    else {
      perror("illegal type");
      exit(EXIT_FAILURE);
    }
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
  int i,j;
  printf("%ld ",it);
  for (i=0;i<(c->size);i++) {
    printf("%f ",gsl_vector_get(c,i));
  }
  for (i=0;i<n_in;i++) {
    for (j=0;j<1;j++) {
      printf("%d ",gsl_vector_int_get(in[i]->nv,j));
    }
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

void init_in(gsl_vector_int *type_vec)
{
  /* re-write this in a more generic way! */
  int i,r0,r1,r2,y0,y0n,y1,y1n;
  in=malloc(n_in*sizeof(interval *));
  for (i=0;i<n_in;i++) {
    in[i]=interval_alloc(PARAM,1.0);
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
    gsl_vector_int_set(type_vec,0,r0);
    gsl_vector_int_set(type_vec,1,r1);
    gsl_vector_int_set(type_vec,2,r2);
    while (1) {
      interval_set(rng,in[i],type_vec);
      if (interval_ok(in[i],y0,y1)) { break; }
    }
  }
}

gsl_matrix_int * read_data(char *filename)
{
  FILE *s;
  int i,j;
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
      gsl_matrix_int_set(data,i,j,(int) temp);
    }
  }
  fclose(s);
  return(data);
}


/* eof */


