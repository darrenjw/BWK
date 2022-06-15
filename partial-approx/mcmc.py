#!/usr/bin/python2
# mcmc.py

# (C) 2005 Darren Wilkinson
# http://www.staff.ncl.ac.uk/d.j.wilkinson/


import sys,getopt

def mcmc(inS,outS,burn,thin,chop):
    count=0
    nextline=inS.readline()
    while (nextline):
        line=nextline
        nextline=inS.readline()
        if ((count%thin==0) and (count==0 or count>burn)):
            if ((not chop) or nextline):
                outS.write(line)
        count=count+1


if (__name__=='__main__'):
    try:
        opt,args=getopt.getopt(sys.argv[1:],'b:t:ch')
    except:
        sys.stderr.write('Do "'+sys.argv[0]+' -h" for a usage summary\n')
        sys.exit(1)
    burn=0
    thin=1
    chop=None
    for o,a in opt:
        if (o=='-b'):
            burn=eval(a)
        if (o=='-t'):
            thin=eval(a)
        if (o=='-c'):
            chop="true"
        if (o=='-h'):
            print 'Usage: '+sys.argv[0]+' [-b <burn>] [-t <thin>] [-c]'
            print 'Defaults: -b 0 -t 1'
            print '-c chops last line off the input'
            print 'Filters tabular standard input to standard output'
            print 'eg. "'+sys.argv[0]+' -b 1000 -t 10 < in.tab > out.tab"'
            print 'or "'+sys.argv[0]+' -c < running.tab > snapshot.tab"'
            sys.exit(1)
    mcmc(sys.stdin,sys.stdout,burn,thin,chop)
    

# eof

