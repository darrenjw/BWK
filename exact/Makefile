# Makefile

CFLAGS=-g -O2 -Wall -ansi -pedantic
#CFLAGS=-g -Wall -ansi -pedantic -DDEBUG
#CFLAGS=-O2 -ansi
LDFLAGS=-lgsl -lgslcblas -lm

OBS=50

lv.tab: lotka_volterra
	date
	nice ./lotka_volterra 50000 $(OBS)| mcmc.py -t 2 > lv.tab
	date

lv-big.tab: lotka_volterra
	nice ./lotka_volterra 500000 $(OBS) | mcmc.py -t 20 > lv-big.tab

lotka_volterra: lotka_volterra.o interval.o
	gcc $(CFLAGS) -o lotka_volterra interval.o lotka_volterra.o $(LDFLAGS)

clean:
	rm -f *.tab *.sim *.ind *.o *.out *~ *.bak core interval_test lotka_volterra

print:
	a2ps interval.h interval.c interval_test.c Makefile







# eof


