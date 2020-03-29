/*
test-code.c

Code to properly test the interval code and inhomogenous PP proposal
code...

*/

#include "interval.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_sort_vector.h>


/* global variables */

long stats_n;
double stats_sumx,stats_sumxx;
int test_fails;


/* function prototypes */

int dtest(char *name, double actual, double expected, double tol);
void stats_clear();
void stats_add(double x);
double stats_mean();
double stats_var();
int mean_test(char *name,double smean,double tmean,double tvar,double n);
int poisson_test(char *name,gsl_rng *r,double mean,long n);
void poisson_tests(gsl_rng *r,long n);
void exp_tests(gsl_rng *r,long n);
int exp_test(char *name,gsl_rng *r,double mean,long n);
void first_event_test(gsl_rng *r,long n);
void last_event_test(gsl_rng *r,long n);
void half_poisson_test(gsl_rng *r,long n);
void rescale_poisson_test(gsl_rng *r,long n);
void transform_in(interval *invl,double l0,double l1,double l2,double l0n,double l1n,double l2n);
void forwards_sim(gsl_rng *r,long n);

/* main function */

int main(int argc,char *argv[])
{
  /* main variables */
  gsl_rng *r;
  /* initialisations */
  test_fails=0;
  r=gsl_rng_alloc(gsl_rng_mt19937);
  /* tests */
  dtest("pass",1.001,1.0,0.01);
  dtest("fail",1.001,1.0,0.00001);
  poisson_tests(r,100000);
  exp_tests(r,100000);
  first_event_test(r,100000);
  last_event_test(r,100000);
  half_poisson_test(r,100000);
  rescale_poisson_test(r,100000);
  forwards_sim(r,10000);
  /* summary */
  printf("%d fail(s) in total\n",test_fails);
  return(EXIT_SUCCESS);
}


/* helper functions */


void forwards_sim(gsl_rng *r,long n)
{
  long i;
  int y0,y1;
  double u,t,a0,a1,a2,a,c0,c1,c2;
  stats_clear();
  c0=0.5; c1=0.0025; c2=0.3;
  for (i=0;i<n;i++) {
    y0=0;
    while (y0 != 137) {
      y0=97; y1=78;
      t=0;
      while (t<=1) {
	a0=c0*y0; a1=c1*y0*y1; a2=c2*y1;
	a=a0+a1+a2;
	t+=gsl_ran_exponential(r,1.0/a);
	if (t>1.0)
	  break;
	u=gsl_rng_uniform(r);
	if (u<a0/a) {
	  y0++;
	} else if (u<(a0+a1)/a) {
	  y0--; y1++;
	} else {
	  y1--;
	}
      }
    }
    stats_add(y1);
  }
  printf("mean of y1 is %f\n",stats_mean());
  printf("var of y1 is %f\n",stats_var());
}

void rescale_poisson_test(gsl_rng *r,long n)
{
  long i;
  int j,count;
  interval *in;
  gsl_vector_int *types;
  int n0,n1,n2;
  stats_clear();
  in=interval_alloc(3,1.0);
  types=gsl_vector_int_calloc(3);
  for (i=0;i<n;i++) {
    n0=gsl_ran_poisson(r,30);
    n1=gsl_ran_poisson(r,40);
    n2=gsl_ran_poisson(r,10);
    gsl_vector_int_set(types,0,n0);
    gsl_vector_int_set(types,1,n1);
    gsl_vector_int_set(types,2,n2);
    interval_set(r,in,types);
    transform_in(in,20,50,10,40,30,10);
    count=0;
    for (j=0;j<n;j++) {
      if (in->time[j] > 0.5) break;
      if (in->type[j] == 1) count++;
    }
    stats_add(count);
  }
  mean_test("rescale poisson",stats_mean(),22.5,22.5,n);
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


void half_poisson_test(gsl_rng *r,long n)
{
  long i;
  int j,count;
  interval *in;
  gsl_vector_int *types;
  int n0,n1;
  stats_clear();
  in=interval_alloc(2,1.0);
  types=gsl_vector_int_calloc(2);
  for (i=0;i<n;i++) {
    n0=gsl_ran_poisson(r,10);
    n1=gsl_ran_poisson(r,15);
    gsl_vector_int_set(types,0,n0);
    gsl_vector_int_set(types,1,n1);
    interval_set(r,in,types);
    count=0;
    for (j=0;j<n;j++) {
      if (in->time[j] > 0.5) break;
      if (in->type[j] == 0) count++;
    }
    stats_add(count);
  }
  mean_test("half poisson",stats_mean(),5,5,n);
}

void last_event_test(gsl_rng *r,long n)
{
  long i;
  interval *in;
  gsl_vector_int *types;
  int n0,n1;
  stats_clear();
  in=interval_alloc(2,1.0);
  types=gsl_vector_int_calloc(2);
  for (i=0;i<n;i++) {
    n0=gsl_ran_poisson(r,10);
    n1=gsl_ran_poisson(r,15);
    gsl_vector_int_set(types,0,n0);
    gsl_vector_int_set(types,1,n1);
    interval_set(r,in,types);
    stats_add(1-(in->time[(in->n)-1]));
  }
  mean_test("last event",stats_mean(),1.0/25,1.0/225,n);
}

void first_event_test(gsl_rng *r,long n)
{
  long i;
  interval *in;
  gsl_vector_int *types;
  int n0,n1;
  stats_clear();
  in=interval_alloc(2,1.0);
  types=gsl_vector_int_calloc(2);
  for (i=0;i<n;i++) {
    n0=gsl_ran_poisson(r,10);
    n1=gsl_ran_poisson(r,15);
    gsl_vector_int_set(types,0,n0);
    gsl_vector_int_set(types,1,n1);
    interval_set(r,in,types);
    stats_add(in->time[0]);
  }
  mean_test("first event",stats_mean(),1.0/25,1.0/225,n);
}

int dtest(char *name, double actual, double expected, double tol)
{
  if ( fabs(actual-expected)>tol ) {
    printf("test %s FAILED: actual=%f, expected=%f, tol=%f\n",
	   name,actual,expected,tol);
    fprintf(stderr,"test %s FAILED\n",name);
    test_fails++;
    return(1);
  } else {
    printf("test %s passed: actual=%f, expected=%f, tol=%f\n",
	   name,actual,expected,tol);
    return(0);
  }
}

void stats_clear()
{
  stats_n=0;
  stats_sumx=0.0;
  stats_sumxx=0.0;
}

void stats_add(double x)
{
  stats_n++;
  stats_sumx+=x;
  stats_sumxx+=x*x;
}

double stats_mean()
{
  return(stats_sumx/stats_n);
}

double stats_var()
{
  return((stats_sumxx-stats_sumx*stats_sumx/stats_n)/(stats_n-1));
}

int mean_test(char *name,double smean,double tmean,double tvar,double n)
{
  return(dtest(name,smean,tmean,3*sqrt(tvar/n)));
}

int exp_test(char *name,gsl_rng *r,double mean,long n)
{
  long i;
  stats_clear();
  for (i=0;i<n;i++) {
    stats_add(gsl_ran_exponential(r,mean));
  }
  return(mean_test(name,stats_mean(),mean,mean*mean,n));
}

int poisson_test(char *name,gsl_rng *r,double mean,long n)
{
  long i;
  stats_clear();
  for (i=0;i<n;i++) {
    stats_add(gsl_ran_poisson(r,mean));
  }
  return(mean_test(name,stats_mean(),mean,mean,n));
}

void poisson_tests(gsl_rng *r,long n)
{
  poisson_test("poisson 2.3",r,2.3,n);
  poisson_test("poisson 5.5",r,5.5,n);
  poisson_test("poisson 23.1",r,23.1,n);
}

void exp_tests(gsl_rng *r,long n)
{
  exp_test("exp 1.4",r,1.4,n);
  exp_test("exp 0.1",r,0.1,n);
  exp_test("exp 5.3",r,5.3,n);
}


/* eof */

