/*
test-mcmc.c

Code to do a vanilla MCMC on one end interval...
*/

/* includes */
#include "interval.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_sort_vector.h>

/* global variables */
double c0,c1,c2;

/* helper function prototypes */
void interval_mcmc(gsl_rng *r,long n);
void transform_in(interval *invl,double l0,double l1,double l2,double l0n,double l1n,double l2n); 
int interval_ok(interval *invl,int y0,int y1);
double log_rnd(interval *invl,int y0,int y1,int y0n,int y1n);
double int_in0(interval *invl,int x0);
double int_in1(interval *invl,int y0,int y1);
double int_in2(interval *invl,int x0);
double int_in(interval *invl,int y0,int y1);
double log_sum_in(interval *invl,int y0,int y1,int y0n,int y1n);
double log_pi(int x,double mean);


/* main function */
int main(int argc,char *argv[])
{
  gsl_rng *r;
  long iters;
  /* process command line */
  if (argc!=2) {
    fprintf(stderr,"Usage: %s <NumIters>\n",argv[0]);
    return(EXIT_FAILURE);
  }
  iters=atoi(argv[1]);
  /* initialise */
  r=gsl_rng_alloc(gsl_rng_mt19937);
  printf("Iter r0 r1 r2 y1\n");
  /* main loop */
  interval_mcmc(r,iters);
  return(EXIT_SUCCESS);
}


/* helper functions */
void interval_mcmc(gsl_rng *r,long n)
{
  long i;
  int y1,y1o,r0,r1,r2,r0o,r1o,r2o,e0,e2;
  double a,a0,a1,a2,a0n,a1n,a2n;
  interval *in,*inprop;
  gsl_vector_int *typevec;
  c0=0.5; c1=0.0025; c2=0.3;
  a0=c0*97; a1=c1*97*78; a2=c2*78;
  a0n=c0*137;
  in=interval_alloc(3,1.0);
  inprop=interval_alloc(3,1.0);
  typevec=gsl_vector_int_alloc(3);
  gsl_vector_int_set(typevec,0,53);
  gsl_vector_int_set(typevec,1,13);
  gsl_vector_int_set(typevec,2,10);
  interval_set(r,in,typevec);
  for (i=0;i<n;i++) {
    r0o=gsl_vector_int_get(in->nv,0);
    r1o=gsl_vector_int_get(in->nv,1);
    r2o=gsl_vector_int_get(in->nv,2);
    y1o=78+r1o-r2o;
    e0=gsl_ran_poisson(r,3)-gsl_ran_poisson(r,3);
    e2=gsl_ran_poisson(r,2)-gsl_ran_poisson(r,2);
    r0=r0o+e0;
    r2=r2o+e2;
    r1=r0-40;
    y1=78+r1-r2;
    if ((r0>=0)&&(r1>=0)&&(r2>=0)&&(y1>=0)) {
      /* construct proposal */
      gsl_vector_int_set(typevec,0,r0);
      gsl_vector_int_set(typevec,1,r1);
      gsl_vector_int_set(typevec,2,r2);
      interval_set(r,inprop,typevec);
      a1n=c1*137*y1; a2n=c2*y1;
      transform_in(inprop,a0,a1,a2,a0n,a1n,a2n);
      /* now evaluate acceptance probability */
      if (interval_ok(inprop,97,78)) {
	/* now evaluate log acceptance probability */
	a=log_rnd(inprop,97,78,137,y1) 
	  -log_rnd(in,97,78,137,y1o)
	  +log_pi(r0,c0*(97+137)/2.0)
	  +log_pi(r1,c1*(97*78+137*y1)/2.0)
	  +log_pi(r2,c2*(78+y1)/2.0)
	  -log_pi(r0o,c0*(97+137)/2.0)
	  -log_pi(r1o,c1*(97*78+137*y1o)/2.0)
	  -log_pi(r2o,c2*(78+y1o)/2.0);
	/* accept! */
	if (log(gsl_rng_uniform(r))<a)
	  interval_cpy(in,inprop);
      }
    }
    /* output current state */
    printf("%ld ",i);
    printf("%d ",gsl_vector_int_get(in->nv,0));
    printf("%d ",gsl_vector_int_get(in->nv,1));
    printf("%d ",gsl_vector_int_get(in->nv,2));
    printf("%d\n",78+gsl_vector_int_get(in->nv,1)
	   - gsl_vector_int_get(in->nv,2));
  }
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
  double in_x,log_sum,a0,a1;
  in_x=int_in(invl,y0,y1);
  log_sum=log_sum_in(invl,y0,y1,y0n,y1n);
  a0=c0*y0 + c1*y0*y1 + c2*y1;
  a1=c0*y0n + c1*y0n*y1n + c2*y1n;
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

double int_in(interval *invl,int y0,int y1)
{
  return( c0*int_in0(invl,y0) + c1*int_in1(invl,y0,y1) + c2*int_in2(invl,y1) );
}

double log_sum_in(interval *invl,int y0,int y1,int y0n,int y1n)
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




/* eof */

