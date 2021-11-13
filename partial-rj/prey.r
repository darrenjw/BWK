plottab<-function(filename="mcmc.tab",lm=400,rows=4)
{
        op<-par(mfrow=c(rows,3))
        tab<-read.table(filename,header=TRUE,check.names=FALSE)
        names<-attr(tab,"names")
        count<-0
        for (i in names[2:length(names)]) {
                plot(ts(tab[[i]]),
                        main=paste("Trace plot for",i),
                        ylab="Value",xlab="Iteration",col=4)
                v<-var(tab[[i]])
                if (v>1e-10) {
                        acf(tab[[i]],lag.max=lm,ci=0,
                                main=paste("ACF plot for",i),col=2)
                        d<-density(tab[[i]])
                        plot(d,
                                main=paste("Density for",i),
                                xlab="Value",col=4,lwd=3)
                } else {
                        plot(0,0)
                        plot(0,0)
                }
                print(paste("Mean for",i,"is",mean(tab[[i]]),"SD
is",sqrt(v)))
                count<-count+1
                if (count %% rows == 0) {
                        readline("Press return to continue... ")
                }
        }
        par(op)
        NULL
}

zscore<-function(m1,v1,n1,m2,v2,n2)
{
    (m1-m2)/sqrt(v1/n1+v2/n2)
}

cmp<-function(rjb,i,djw,j)
{
    z<-zscore(mean(rjb[,i]),var(rjb[,i]),length(rjb[,i]),mean(djw[,j]),var(djw[,j]),length(djw[,j]))
    print(paste("zscore for",names(rjb)[i],names(djw)[j],"is",z))
	z
}

compare<-function(rjb,rjblist,djw,djwlist)
{
	rjblen<-length(rjblist)
	sumz<-0
	for (i in 1:rjblen) {
	   z<-cmp(rjb,rjblist[[i]],djw,djwlist[[i]])
	sumz<-sumz+z*z
	}
	print(paste("Overall: chisquare",sumz,"on",rjblen,"d.f. - p=",1-pchisq(sumz,rjblen)))
}

library(stats)
mah<-function(rjb,djw)
{
	print(cor(rjb[,2:4]))
	ch<-mahalanobis(colMeans(as.matrix(djw[,2:4])),colMeans(as.matrix(rjb[,2:4])),2*var(as.matrix(rjb[,2:4]))/dim(as.matrix(rjb[,2:4]))[1])
	print(paste("Mahalanobis chisquare",ch,"on 3 d.f. - p=",1-pchisq(ch,3)))
}


#compare(rjb,rjblist,djw,djwlist)
#mah(rjb,djw)

# 40 intervals, Table 1
rjb<-read.table("rjb40-e.tab",header=TRUE,check.names=FALSE)
djw<-read.table("fred.tab",header=TRUE,check.names=FALSE)
#djw<-read.table("lv-41-e-5e6.tab",header=TRUE,check.names=FALSE)
rjblist<-c(2,3,4)
djwlist<-c(2,3,4)
compare(rjb,rjblist,djw,djwlist)
mah(rjb,djw)
#djwa<-read.table("lv-41-a-5e7.tab",header=TRUE,check.names=FALSE)
#compare(djw,djwlist,djwa,djwlist)

rjb3<-read.table("rjb40-e-3.tab",header=TRUE,check.names=FALSE)

compare(rjb3,rjblist,rjb,rjblist)
mah(rjb,rjb3)

compare(rjb3,rjblist,djw,djwlist)
mah(rjb,djw)

rjb<-read.table("rjb40-e.tab",header=TRUE,check.names=FALSE)
djw<-read.table("lv-41-e-5e5.tab",header=TRUE,check.names=FALSE)
rjblist<-c(2,3,4)
djwlist<-c(2,3,4)
compare(rjb,rjblist,djw,djwlist)
mah(rjb,djw)
#djwa<-read.table("lv-41-a-5e7.tab",header=TRUE,check.names=FALSE)
#compare(djw,djwlist,djwa,djwlist)

rjb3<-read.table("rjb40-e-3.tab",header=TRUE,check.names=FALSE)

compare(rjb3,rjblist,rjb,rjblist)
mah(rjb,rjb3)

compare(rjb3,rjblist,djw,djwlist)
mah(rjb,djw)

rjb<-read.table("rjb40-e-3.tab",header=TRUE,check.names=FALSE)
djw<-read.table("lv-41-e2-5e6.tab",header=TRUE,check.names=FALSE)
rjblist<-c(2,3,4)
djwlist<-c(2,3,4)
compare(rjb,rjblist,djw,djwlist)
mah(rjb,djw)
#djwa<-read.table("lv-41-a-5e7.tab",header=TRUE,check.names=FALSE)
#compare(djw,djwlist,djwa,djwlist)

rjb3<-read.table("rjb40-e-3.tab",header=TRUE,check.names=FALSE)

compare(rjb3,rjblist,rjb,rjblist)
mah(rjb,rjb3)

compare(rjb3,rjblist,djw,djwlist)
mah(rjb,djw)

rjb<-read.table("rjb40-epif-2.tab",header=TRUE,check.names=FALSE)
djw<-read.table("lv-41-epsif-5e6.tab",header=TRUE,check.names=FALSE)
rjblist<-c(2,3,4)
djwlist<-c(2,3,4)
compare(rjb,rjblist,djw,djwlist)
mah(rjb,djw)

rjb<-read.table("rjb40-epi-2.tab",header=TRUE,check.names=FALSE)
djw<-read.table("lv-41-epsi-5e6.tab",header=TRUE,check.names=FALSE)
rjblist<-c(2,3,4)
djwlist<-c(2,3,4)
compare(rjb,rjblist,djw,djwlist)
mah(rjb,djw)


rjb<-read.table("rjb200-e.tab",header=TRUE,check.names=FALSE)
djw<-read.table("lv-201-e-5e5.tab",header=TRUE,check.names=FALSE)
rjblist<-c(2,3,4)
djwlist<-c(2,3,4)
compare(rjb,rjblist,djw,djwlist)
mah(rjb,djw)

rjb<-read.table("rjb200-e.tab",header=TRUE,check.names=FALSE)
djw<-read.table("lv-201-es-5e6.tab",header=TRUE,check.names=FALSE)
rjblist<-c(2,3,4)
djwlist<-c(2,3,4)
compare(rjb,rjblist,djw,djwlist)
mah(rjb,djw)

rjb<-read.table("rjb200-fine-e.tab",header=TRUE,check.names=FALSE)
djw<-read.table("lv-fine-201-e-5e6.tab",header=TRUE,check.names=FALSE)
rjblist<-c(2,3,4)
djwlist<-c(2,3,4)
compare(rjb,rjblist,djw,djwlist)
mah(rjb,djw)


# 5 intervals
rjb<-read.table("rjb5-e.tab",header=TRUE,check.names=FALSE)
djw<-read.table("lv-6-e-5e6.tab",header=TRUE,check.names=FALSE)
rjblist<-c(2,3,4)
djwlist<-c(2,3,4)
compare(rjb,rjblist,djw,djwlist)
mah(rjb,djw)

rjb<-read.table("rjb5-e-3.tab",header=TRUE,check.names=FALSE)
djw<-read.table("lv-6-e-5e7-50k.tab",header=TRUE,check.names=FALSE)
rjblist<-c(2,3,4)
djwlist<-c(2,3,4)
compare(rjb,rjblist,djw,djwlist)
mah(rjb,djw)

# partials

rjb<-read.table("rjb40-epif.tab",header=TRUE,check.names=FALSE)
djw<-read.table("lv-41-epsif-5e5.tab",header=TRUE,check.names=FALSE)
rjblist<-c(2,3,4)
djwlist<-c(2,3,4)
compare(rjb,rjblist,djw,djwlist)
mah(rjb,djw)


















# 3 intervals, both ends fixed
rjb<-read.table("fred7.tab",header=TRUE,check.names=FALSE)
djw<-read.table("lv-4-epik.tab",header=TRUE,check.names=FALSE)
forward<-read.table("forward3.tab",header=TRUE,check.names=FALSE)
rjblist<-c(10,11,12,14,15,16,17,19,20,21,22)
djwlist<-c(6,7,8,9,10,11,12,13,14,15,16)
forlist<-c(2,3,4,5,6,7,8,9,10,11,12)
compare(rjb,rjblist,djw,djwlist)
compare(rjb,rjblist,forward,forlist)
compare(forward,forlist,djw,djwlist)

djw2<-read.table("lv-4-epifk.tab",header=TRUE,check.names=FALSE)
compare(rjb,rjblist,djw2,djwlist)
compare(forward,forlist,djw2,djwlist)


# 2 intervals, both ends fixed
rjb<-read.table("fred5.tab",header=TRUE,check.names=FALSE)
djw<-read.table("lv-3-epifk.tab",header=TRUE,check.names=FALSE)
forward<-read.table("forward.tab",header=TRUE,check.names=FALSE)
rjblist<-c(10,11,12,14,15,16,17)
djwlist<-c(6,7,8,9,10,11,12,13,14,15,16)
forlist<-c(2,3,4,5,6,7,8)
compare(rjb,rjblist,djw,djwlist)
compare(rjb,rjblist,forward,forlist)
compare(forward,forlist,djw,djwlist)

# 2 intervals, initial end fixed
rjb<-read.table("fred4.tab",header=TRUE,check.names=FALSE)
djw<-read.table("lv-3-epik.tab",header=TRUE,check.names=FALSE)
forward<-read.table("forward2.tab",header=TRUE,check.names=FALSE)
rjblist<-c(10,11,12,14,15,16,17,19)
djwlist<-c(6:13)
forlist<-c(2:9)
compare(rjb,rjblist,djw,djwlist)
compare(rjb,rjblist,forward,forlist)
compare(forward,forlist,djw,djwlist)

rjb<-read.table("fred4.tab",header=TRUE,check.names=FALSE)
djw<-read.table("lv-3-epik-new.tab",header=TRUE,check.names=FALSE)
#ag<-read.table("ag-3-epik.tab",header=TRUE,check.names=FALSE)
rjblist<-c(10,11,12,14,15,16,17,19)
djwlist<-c(6:13)
#aglist<-c(6:13)
compare(rjb,rjblist,djw,djwlist)
#compare(rjb,rjblist,ag,aglist)
#compare(ag,aglist,djw,djwlist)



forward<-read.table("forward2.tab",header=TRUE,check.names=FALSE)
djwf<-read.table("lv-3-epik-forwards.tab",header=TRUE,check.names=FALSE)
forlist<-c(2:8)
djwflist<-c(2:8)
compare(forward,forlist,djwf,djwflist)

djw<-read.table("lv-3-epik.tab",header=TRUE,check.names=FALSE)
djwf<-read.table("lv-3-epik-forwards.tab",header=TRUE,check.names=FALSE)
djwlist<-c(6:12)
djwflist<-c(2:8)
compare(djw,djwlist,djwf,djwflist)

# 4 intervals, both ends fixed, variable rate constants
rjb<-read.table("fred6.tab",header=TRUE,check.names=FALSE)
djw<-read.table("lv-5-epif.tab",header=TRUE,check.names=FALSE)
rjblist<-c(2,3,4,10:12,14:17,19:22,24:26)
djwlist<-c(2,3,4,6:20)
compare(rjb,rjblist,djw,djwlist)

rjb<-read.table("fred6.tab",header=TRUE,check.names=FALSE)
djw<-read.table("lv-5-epif-djw-c.tab",header=TRUE,check.names=FALSE)
ag<-read.table("lv-5-epif-ag-c.tab",header=TRUE,check.names=FALSE)
rjblist<-c(2,3,4,10:12,14:17,19:22,24:26)
djwlist<-c(2,3,4,6:20)
aglist<-c(2,3,4,6:20)
compare(rjb,rjblist,djw,djwlist)
compare(rjb,rjblist,ag,aglist)

rjb<-read.table("fred6.tab",header=TRUE,check.names=FALSE)
djw<-read.table("lv-5-epif-new.tab",header=TRUE,check.names=FALSE)
rjb<-read.table("rub.tab",header=TRUE,check.names=FALSE)
djw<-read.table("rubd.tab",header=TRUE,check.names=FALSE)
rjblist<-c(2,3,4,10:12,14:17,19:22,24:26)
djwlist<-c(2,3,4,6:20)
compare(rjb,rjblist,djw,djwlist)

rjb<-read.table("fred6.tab",header=TRUE,check.names=FALSE)
djw<-read.table("lv-5-epif-big.tab",header=TRUE,check.names=FALSE)
rjblist<-c(2,3,4,10:12,14:17,19:22,24:26)
djwlist<-c(2,3,4,6:20)
compare(rjb,rjblist,djw,djwlist)

rjb<-read.table("jobs.tab",header=TRUE,check.names=FALSE)
djw<-read.table("lv-5-epif-big.tab",header=TRUE,check.names=FALSE)
rjblist<-c(2,3,4,10:12,14:17,19:22,24:26)
djwlist<-c(2,3,4,6:20)
compare(rjb,rjblist,djw,djwlist)


# 1 interval - the 2nd one - initial end fixed
rjb<-read.table("fred8.tab",header=TRUE,check.names=FALSE)
djw<-read.table("lv-3-ep-endonly-new.tab",header=TRUE,check.names=FALSE)
forward<-read.table("for2.tab",header=TRUE,check.names=FALSE)
rjblist<-c(15,16,17,19)
djwlist<-c(10:13)
forlist<-c(2:5)
compare(rjb,rjblist,djw,djwlist)
compare(rjb,rjblist,forward,forlist)
compare(forward,forlist,djw,djwlist)

forward<-read.table("for2.tab",header=TRUE,check.names=FALSE)
djwf<-read.table("lv-3-ep-endonly-forwards.tab",header=TRUE,check.names=FALSE)
forlist<-c(2,5)
djwflist<-c(2,3)
compare(forward,forlist,djwf,djwflist)

djw<-read.table("lv-3-ep-endonly.tab",header=TRUE,check.names=FALSE)
djwf<-read.table("lv-3-ep-endonly-forwards.tab",header=TRUE,check.names=FALSE)
djwlist<-c(10,13)
djwflist<-c(2,3)
compare(djw,djwlist,djwf,djwflist)
