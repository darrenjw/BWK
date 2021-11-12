# Makefile


toy.tab: a.out
	./job_script in0


a.out: prey.f
	cd ../randlib.f/src ; make
	f77 -O2 prey.f ../randlib.f/src/randlib.a



clean:
	rm a.out

# eof
