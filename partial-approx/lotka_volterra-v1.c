/*
lotka_volterra.c (v1)

Code for approximate MCMC inference in the lotka-volterra model

Partially observed version - only prey observed

This version uses true predator numbers for initialisation
purposes

Darren Wilkinson
d.j.wilkinson@ncl.ac.uk
http://www.staff.ncl.ac.uk/d.j.wilkinson/

Last updated: 14/5/2004
*/

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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
#define GMIG gsl_matrix_int_get
#define GMIS gsl_matrix_int_set
#define GVG gsl_vector_get
#define GVS gsl_vector_set

/* global variables */
gsl_vector *c;           /* parameter vector */
gsl_matrix_int *r;       /* reaction matrix */
int n_in;                /* number of intervals */
gsl_matrix_int *data;    /* data vector */
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
  
  /* fix params */
  gsl_vector_set(c,0,0.5);
  gsl_vector_set(c,1,0.0025);
  gsl_vector_set(c,2,0.3);

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
  return(EXIT_SUCCESS);
}


/* helper functions */

void update_c(double alpha,double beta)
{
  int i,r0,r1,r2;
  double c0,c1,c2,in0,in1,in2;
  r0=0; r1=0; r2=0;
  in0=0.0; in1=0.0; in2=0.0;
  for (i=0;i<n_in;i++) {
    r0+=GMIG(r,i,0);
    r1+=GMIG(r,i,1);
    r2+=GMIG(r,i,2);
    in0+=0.5*(GMIG(data,i,0)+GMIG(data,i+1,0));
    in1+=0.5*(GMIG(data,i,0)*GMIG(data,i,1)+GMIG(data,i+1,0)*GMIG(data,i+1,1));
    in2+=0.5*(GMIG(data,i,1)+GMIG(data,i+1,1));
  }
  c0=gsl_ran_gamma(rng,alpha+r0,1.0/(beta+in0));
  c1=gsl_ran_gamma(rng,alpha+r1,1.0/(beta+in1));
  c2=gsl_ran_gamma(rng,alpha+r2,1.0/(beta+in2));
  GVS(c,0,c0);
  GVS(c,1,c1);
  GVS(c,2,c2);
}

void update_in(void)
{
  int i,r0,r1,r2,r0n,r1n,r2n,r0o,r1o,r2o,r0no,r1no,r2no;
  int y0,y1,y1s,y0n,y1n,y1ns,y0nn,y1nn,e0,e1,e2,en;
  double a,theta0,theta1,theta2,thetan;
  double c0,c1,c2;
  c0=GVG(c,0); c1=GVG(c,1); c2=GVG(c,2);
  /* update first interval */
  /*
  y0=GMIG(data,0,0); y1=GMIG(data,0,1);
  y0n=GMIG(data,1,0); y1n=GMIG(data,1,1);
  r0o=GMIG(r,0,0); r1o=GMIG(r,0,1); r2o=GMIG(r,0,2);
  theta0=(double) r0o*r0o/200.0 + 1.0;
  e0=gsl_ran_poisson(rng,theta0)-gsl_ran_poisson(rng,theta0);
  theta1=(double) y1*y1/200.0 + 1.0;
  e1=gsl_ran_poisson(rng,theta1)-gsl_ran_poisson(rng,theta1);
  y1s=y1+e1;
  r0=r0o+e0;
  r1=y0-y0n+r0; r2=y0-y0n+y1s-y1n+r0;
  if ( (r0>=0) && (r1>=0) && (r2>=0) && (y1s>=0) ) {
    a=-log_pi_poisd(e0,theta0)
      -log_pi_poisd(e1,theta1)
      +log_pi(r0,c0*(y0+y0n)/2.0)
      +log_pi(r1,c1*(y0*y1s+y0n*y1n)/2.0)
      +log_pi(r2,c2*(y1s+y1n)/2.0)
      +log_pi_poisd(e0,(double) r0*r0/200.0 + 1.0)
      +log_pi_poisd(e1,(double) y1s*y1s/200.0 + 1.0)
      -log_pi(r0o,c0*(y0+y0n)/2.0)
      -log_pi(r1o,c1*(y0*y1s+y0n*y1n)/2.0)
      -log_pi(r2o,c2*(y1s+y1n)/2.0);
    if (accept_lp(a)) {
      GMIS(r,0,0,r0); GMIS(r,0,1,r1); GMIS(r,0,2,r2);
      GMIS(data,0,1,y1s);
    }
  }
  */
  /* update interval pairs */
  for (i=0;i<n_in-1;i++) {
    y0=GMIG(data,i,0); y1=GMIG(data,i,1);
    y0n=GMIG(data,i+1,0); y1n=GMIG(data,i+1,1);
    r0o=GMIG(r,i,0); r1o=GMIG(r,i,1); r2o=GMIG(r,i,2);
    y0nn=GMIG(data,i+2,0); y1nn=GMIG(data,i+2,1);
    r0no=GMIG(r,i+1,0); r1no=GMIG(r,i+1,1); r2no=GMIG(r,i+1,2);
    theta0=(double) r0o*r0o/200.0 + 1.0;
    e0=gsl_ran_poisson(rng,theta0)-gsl_ran_poisson(rng,theta0);
    r0=r0o+e0;
    theta2=(double) r2o*r2o/200.0 + 1.0;
    e2=gsl_ran_poisson(rng,theta2)-gsl_ran_poisson(rng,theta2);
    r2=r2o+e2;
    r1=y0-y0n+r0;
    y1ns=y1+r1-r2;
    thetan=(double) r0no*r0no/200.0 + 1.0;
    en=gsl_ran_poisson(rng,thetan)-gsl_ran_poisson(rng,thetan);
    r0n=r0no+en;
    r1n=y0n-y0nn+r0n;
    r2n=y0n-y0nn+y1ns-y1nn+r0n;
    if ( (r0 >= 0) && (r1 >= 0) && (r2 >= 0)
	 && (r0n >= 0) && (r1n >= 0) && (r2n >= 0) && (y1ns >= 0) ) {
      /* now evaluate log acceptance probability */
      a=-log_pi_poisd(e0,theta0)
	-log_pi_poisd(e2,theta2)
	+log_pi(r0,c0*(y0+y0n)/2.0)
	+log_pi(r1,c1*(y0*y1+y0n*y1n)/2.0)
	+log_pi(r2,c2*(y1+y1n)/2.0)
	+log_pi_poisd(e0,(double) r0*r0/200.0 + 1.0)
	+log_pi_poisd(e2,(double) r2*r2/200.0 + 1.0)
	-log_pi(r0o,c0*(y0+y0n)/2.0)
	-log_pi(r1o,c1*(y0*y1+y0n*y1n)/2.0)
	-log_pi(r2o,c2*(y1+y1n)/2.0)
	-log_pi_poisd(en,thetan)
	+log_pi(r0n,c0*(y0n+y0nn)/2.0)
	+log_pi(r1n,c1*(y0n*y1ns+y0nn*y1nn)/2.0)
	+log_pi(r2n,c2*(y1ns+y1nn)/2.0)
	+log_pi_poisd(en,(double) r0n*r0n/200.0 + 1.0)
	-log_pi(r0no,c0*(y0n+y0nn)/2.0)
	-log_pi(r1no,c1*(y0n*y1n+y0nn*y1nn)/2.0)
	-log_pi(r2no,c2*(y1n+y1nn)/2.0);
#ifdef DEBUG
      printf("[a:%f ",a);
      printf("] ");
#endif
      /* test and update */
      if (accept_lp(a)) {
	GMIS(r,i,0,r0); GMIS(r,i,1,r1); GMIS(r,i,2,r2);
	GMIS(r,i+1,0,r0n); GMIS(r,i+1,1,r1n); GMIS(r,i+1,2,r2n);
	GMIS(data,i+1,1,y1ns);
#ifdef DEBUG
	printf("Accept. ");
#endif
      }
    }
  }
#ifdef DEBUG
  printf("\n");
#endif
  /* update last interval */
  
  y0=GMIG(data,n_in-1,0); y1=GMIG(data,n_in-1,1);
  y0n=GMIG(data,n_in,0); y1n=GMIG(data,n_in,1);
  r0o=GMIG(r,n_in-1,0); r1o=GMIG(r,n_in-1,1); r2o=GMIG(r,n_in-1,2);
  theta0=(double) r0o*r0o/200.0 + 1.0;
  e0=gsl_ran_poisson(rng,theta0)-gsl_ran_poisson(rng,theta0);
  r0=r0o+e0;
  theta2=(double) r2o*r2o/200.0 + 1.0;
  e2=gsl_ran_poisson(rng,theta2)-gsl_ran_poisson(rng,theta2);
  r2=r2o+e2; r1=y0-y0n+r0;
  y1ns=y1+r1-r2;
#ifdef DEBUG
  printf("last: r0=%d, r1=%d, r2=%d\n",r0,r1,r2);
  printf("last: y0=%d, y1=%d, y0n=%d, y1ns=%d\n",y0,y1,y0n,y1ns);
#endif
  if ( (r0>=0) && (r1>=0) && (r2>=0) && (y1ns>=0) ) {
    a=-log_pi_poisd(e0,theta0)
      -log_pi_poisd(e2,theta2)
      +log_pi(r0,c0*(y0+y0n)/2.0)
      +log_pi(r1,c1*(y0*y1+y0n*y1n)/2.0)
      +log_pi(r2,c2*(y1+y1n)/2.0)
      +log_pi_poisd(e0,(double) r0*r0/200.0 + 1.0)
      /* +log_pi_poisd(e1,(double) y1*y1/200.0 + 1.0) */
      +log_pi_poisd(e2,(double) r2*r2/200.0 + 1.0)
      -log_pi(r0o,c0*(y0+y0n)/2.0)
      -log_pi(r1o,c1*(y0*y1+y0n*y1ns)/2.0)
      -log_pi(r2o,c2*(y1+y1ns)/2.0);
    if (accept_lp(a)) {
      GMIS(r,n_in-1,0,r0); GMIS(r,n_in-1,1,r1); GMIS(r,n_in-1,2,r2);
      GMIS(data,n_in,1,y1ns);
    }
  }
  
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
    printf("%f ",GVG(c,i));
  }
  for (i=0;i<n_in;i++) {
    printf("%d ",GMIG(data,i,1));
    for (j=0;j<PARAM;j++) {
      printf("%d ",GMIG(r,i,j));
    }
  }
  printf("%d ",GMIG(data,n_in,1));
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
    printf("x[%d] ",i);
    for (j=0;j<PARAM;j++) {
      printf("r[%d,%d] ",i,j);
    }
  }
  printf("x[%d] ",n_in);
  printf("\n");
}

void init_in(void)
{
  /* re-write this in a more generic way! */
  int i,r0,r1,r2,y0,y0n,y1,y1n;
  r=gsl_matrix_int_calloc(n_in,PARAM);
  for (i=0;i<n_in;i++) {
    y0=GMIG(data,i,0); y0n=GMIG(data,i+1,0);
    y1=GMIG(data,i,1); y1n=GMIG(data,i+1,1);
    r0=20; /* do something better here! */
    r0=5;
    if (y0n>y0) { r0+=y0n-y0; }
    if (y1n>y1) { r0+=y1n-y1; }
    r1=y0-y0n+r0; r2=y0-y0n+y1-y1n+r0;
    if ( (r0<0)||(r1<0)||(r2<0) ) {
      perror("Initialisation failed");
      exit(EXIT_FAILURE);
    }
    GMIS(r,i,0,r0); GMIS(r,i,1,r1); GMIS(r,i,2,r2);
  }
#ifdef DEBUG
  for (i=0;i<n_in;i++) {
    printf("i=%d, r0=%d, r1=%d, r2=%d\n",
	   i,GMIG(r,i,0),GMIG(r,i,1),GMIG(r,i,2));
  }
#endif
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
      GMIS(data,i,j,(int) temp);
    }
  }
  fclose(s);
  return(data);
}


/* eof */
