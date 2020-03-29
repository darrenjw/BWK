/*
interval.c

Code for time interval structures

Darren Wilkinson
d.j.wilkinson@ncl.ac.uk
http://www.staff.ncl.ac.uk/d.j.wilkinson/

Last updated: 27/2/2002
*/

#include "interval.h"


/* functions available externally */

interval * interval_alloc(int n_types,double length)
{
  interval *in;
  in=malloc(sizeof(interval));
  (in->n)=0;
  (in->n_max)=10;
  (in->len)=length;
  (in->n_types)=n_types;
  (in->nv)=gsl_vector_int_calloc(n_types);
  (in->nv_temp)=gsl_vector_int_calloc(n_types);
  (in->time)=malloc(10*sizeof(double));
  (in->type)=malloc(10*sizeof(int));
  return(in);
}

void interval_set(gsl_rng *r,interval *in,gsl_vector_int *nv)
{
  int n,i,type;
  double t,nt;
  t=0;
  n=vector_int_sum(nv);
  (in->n)=n;
  gsl_vector_int_memcpy(in->nv,nv);
  gsl_vector_int_memcpy(in->nv_temp,nv);
  /* first realloc if necessary */
  if (n > (in->n_max)) {
    (in->n_max)=2*n;
    (in->time)=realloc((in->time),(in->n_max)*sizeof(double));
    (in->type)=realloc((in->type),(in->n_max)*sizeof(int));
  }
  /* now fill in times and types */
  for (i=0;i<n;i++) {
    type=gen_type(r,in->nv_temp);
    nt=gen_min_uni(r,n-i,(in->len)-t);
    t+=nt;
    (in->time)[i] = t;
    (in->type)[i] = type;
    gsl_vector_int_set(in->nv_temp,type,
		       gsl_vector_int_get(in->nv_temp,type)-1);
  }
}

void interval_cpy(interval *dest,interval *src)
{
  int i;
  /* first realloc if necessary */
  if ( (src->n) > (dest->n_max) ) {
    (dest->n_max)=2*(src->n);
    (dest->time)=realloc((dest->time),(dest->n_max)*sizeof(double));
    (dest->type)=realloc((dest->type),(dest->n_max)*sizeof(int));
  }
  /* now copy across */
  (dest->n)=(src->n);
  (dest->len)=(src->len);
  (dest->n_types)=(src->n_types); /* unnecessary! */
  gsl_vector_int_memcpy(dest->nv,src->nv);
  gsl_vector_int_memcpy(dest->nv_temp,src->nv_temp);
  for (i=0;i<(dest->n);i++) {
    dest->time[i] = src->time[i];
    dest->type[i] = src->type[i];
  }
}

void interval_print(interval *in)
{
  int i;
  printf("n: %d\n",in->n);
  printf("n_max: %d\n",in->n_max);
  printf("len: %f\n",in->len);
  printf("n_types: %d\n",in->n_types);
  printf("nv: ");
  for (i=0;i<(in->n_types);i++) {
    printf("%d ",gsl_vector_int_get(in->nv,i));
  }
  printf("\n");
  printf("nv_temp: ");
  for (i=0;i<(in->n_types);i++) {
    printf("%d ",gsl_vector_int_get(in->nv_temp,i));
  }
  printf("\n");
  printf("Time Type\n");
  for (i=0;i<(in->n);i++) {
    printf("%f %d\n",in->time[i],in->type[i]);
  }
}


/* private functions */

int vector_int_sum(gsl_vector_int *nv)
{
  int sum,i;
  sum=0;
  for (i=0;i<(nv->size);i++) {
    sum+=gsl_vector_int_get(nv,i);
  }
  return(sum);
}

double gen_min_uni(gsl_rng *r,int n,double a)
{
  double u;
  u=gsl_ran_flat(r,0,1);
  return( a*(1-pow(u,1.0/n)) );
}

int gen_type(gsl_rng *r,gsl_vector_int *v)
{
  int ind,runningSum,sum;
  double u;
  sum=vector_int_sum(v);
  runningSum=0;
  u=sum*gsl_ran_flat(r,0,1);
  for (ind=0;ind < (v->size);ind++) {
    runningSum+=gsl_vector_int_get(v,ind);
    if (runningSum > u) break;
  }
  if (ind == (v->size)) ind--;
  return(ind);
}



/* eof */

