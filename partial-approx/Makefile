# Makefile

CC=gcc

#CFLAGS=-g -Wall -ansi -pedantic
#CFLAGS=-g -Wall -ansi -pedantic -DDEBUG
CFLAGS=-O2 -ansi
LDFLAGS=-lgsl -lgslcblas -lm

OBS=50

lv.tab: lotka_volterra
	nice ./lotka_volterra 100000 $(OBS) | ./mcmc.py -t 5 > lv.tab

lv-big.tab: lotka_volterra
	nice ./lotka_volterra 10000000 $(OBS) | ./mcmc.py -t 500 > lv-big.tab

lotka_volterra: lotka_volterra.o
	gcc $(CFLAGS) -o lotka_volterra lotka_volterra.o $(LDFLAGS)

clean:
	rm -f *.tab *.sim *.ind *.o *.out *~ *.bak core interval_test lotka_volterra

print:
	a2ps interval.h interval.c interval_test.c Makefile







# eof


