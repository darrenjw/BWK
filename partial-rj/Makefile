# Makefile


FORCE:
	make toy.txt


ALL:
	make toy.txt
	make rjb40-e.txt
	make rjb200-e.txt
	make rjb200-fine-e.txt
	make rjb40-epif.txt
	make rjb40-epif-2.txt
	make rjb40-epi.txt
	make rjb40-epi-2.txt
	make rjb40-e-2.txt
	make rjb40-e-3.txt


toy.tab: a.out in0
	./job_script in0

rjb40-e.tab: a.out in1
	./job_script in1

rjb200-e.tab: a.out in2
	./job_script in2

rjb200-fine-e.tab: a.out in3
	./job_script in3

rjb40-epif.tab: a.out in4
	./job_script in4

rjb40-epif-2.tab: a.out in5
	./job_script in5

rjb40-epi.tab: a.out in6
	./job_script in6

rjb40-epi-2.tab: a.out in7
	./job_script in7

rjb40-e-2.tab: a.out in8
	./job_script in8

rjb40-e-3.tab: a.out in9
	./job_script in9



a.out: prey.f
	cd ../randlib.f/src ; make
	gfortran -O2 prey.f ../randlib.f/src/randlib.a

%.txt: %.tab
	./summarise.R $< > $@ 2>&1
	cp Rplots.pdf $@.pdf


clean:
	rm a.out

# eof
