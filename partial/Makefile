# Makefile

CC=gcc

#CFLAGS=-g -Wall -ansi -pedantic
#CFLAGS=-g -Wall -ansi -pedantic -DDEBUG
CFLAGS=-O2 -ansi -Wall -pedantic
#CFLAGS=-ansi -Wall -pedantic
LDFLAGS=-lgsl -lgslcblas 

OBS=50
#OBS=3
#OBS=5

lv.tab: lotka_volterra
	nice ./lotka_volterra 100000 $(OBS) | mcmc.py -t 2 > lv.tab

lv-big.tab: lotka_volterra
	nice ./lotka_volterra 1000000 $(OBS) | mcmc.py -t 20 > lv-big.tab

test.tab: test-mcmc
	nice ./test-mcmc 1000000 | mcmc.py -t 20 > test.tab

test: test-code
	nice ./test-code

test-code: test-code.o interval.o
	gcc $(CFLAGS) $(LDFLAGS) -o test-code interval.o test-code.o

test-mcmc: test-mcmc.o interval.o
	gcc $(CFLAGS) $(LDFLAGS) -o test-mcmc interval.o test-mcmc.o

lotka_volterra: lotka_volterra.o interval.o
	gcc $(CFLAGS) $(LDFLAGS) -o lotka_volterra interval.o lotka_volterra.o

lotka_volterra-v4: lotka_volterra-v4.o interval.o
	gcc $(CFLAGS) $(LDFLAGS) -o lotka_volterra-v4 interval.o lotka_volterra-v4.o

clean:
	rm -f *.tab *.sim *.ind *.o *.out *~ *.bak core test-code test-mcmc lotka_volterra

print:
	a2ps interval.h interval.c interval_test.c Makefile







# eof


