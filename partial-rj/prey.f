c
c Initialisation method changed to depend on initial parameter values
c
c Works on numbers of reactions - not triples!
c
      program prey
      external gengamm
      integer rmax,nmax
      parameter (rmax=10000, nmax=251)
      integer i,j,k,x(2,0:nmax),rtype(nmax,rmax),r(3,nmax),niter,nthin,
     &     seed1,seed2,mmaxr(3),npts,sumr(3),nrun,nburn,countm2,counta,
     &     countd,countm,countm22,counta2,counta22,countd2,countd22,
     &     cm2,ca,cd,cm,cm22,ca2,ca22,cd2,cd22,fixedpred,counta3,
     &     countd3,countm3,ca3,cd3,cm3
      double precision times(nmax,rmax),ath(3),bth(3),gengamm,
     &     theta(3),int(3,nmax),lpopprod(nmax),sumint(3)
      character(len=50) fin,fout
      data counta,countd,countm,counta2,counta22,countd2,countd22,
     & countm2,countm22,counta3,countd3,countm3/12*0/
      data ca,cd,cm,ca2,ca22,cd2,cd22,cm2,cm22,ca3,cd3,cm3/12*0/
      data ath,bth/3*1d0,3*0.01d0/
c  read data and initialise missing data
      print *, 'Enter name for data file'
      read '(A50)', fin
      open(unit=7,file=fin)
      print *,'Enter number of intervals (npts)'
      read(5,*) npts
      print *, 'Enter name for output file'
      read '(A50)', fout
      open (unit=8,file=fout)
      if (npts.lt.9) then
         write(8,'(a57,10(a4,4(i1,a5),i1,a1))')
     &        'inter theta[1] theta[2] theta[3] maxr[1] maxr[2] maxr[3]'
     &        ,(' x1[',i,'] x2[',i,'] r1[',i,'] r2[',i,'] r3[',i,']',
     &        i=1,min0(9,npts)),' x1[',npts+1,'] x2[',npts+1,']    '
      else
         write(8,'(a57,9(a4,4(i1,a5),i1,a1),30(a4,4(i2,a5),i2,a1))')
     &        'inter theta[1] theta[2] theta[3] maxr[1] maxr[2] maxr[3]' 
     &        ,(' x1[',i,'] x2[',i,'] r1[',i,'] r2[',i,'] r3[',i,']',
     &        i=1,min0(19,npts)),' x1[',min0(19,npts)+1,'] x2[',
     &        min0(19,npts)+1,']     '
      end if
      do i=0,npts
         read(7,*) x(1,i),x(2,i)
      end do
      close(7)
      print *,'Enter data type for predators: 1=fixed, 2=fixed at ends',
     &     ', 3=fixed start, 4=fixed start on interval npts'
      read(5,*) fixedpred
      print *,'Enter initial values for theta[1-3]'
      read(5,*) (theta(k),k=1,3)
      call initialise(theta,npts,x,rtype,times,r,sumr,int,lpopprod,
     &     sumint)
c
      print *,'Enter seeds'
      read(5,*) seed1,seed2
      call setall(seed1,seed2)
c   repeatable seeds
c      call g05cbf(seed1)
c   non-repeatable seeds
c      call g05ccf
c
      print *,'Enter burn-in, run length and nthin'
      read(5,*) nburn,nrun,nthin
      if (fixedpred.eq.1) then
         do niter=1,nburn
            do k=1,3
               theta(k)=gengamm(bth(k)+sumint(k),ath(k)+sumr(k))
            end do
            call allfixed(npts,theta,rtype,times,r,sumr,x,int,
     &           sumint,lpopprod,ca,cd,cm)
         end do
         do niter=1,nrun
            do j=1,nthin
               do k=1,3
                  theta(k)=gengamm(bth(k)+sumint(k),ath(k)+sumr(k))
               end do
               call allfixed(npts,theta,rtype,times,r,sumr,x,int,
     &              sumint,lpopprod,counta,countd,countm)
            end do
            do k=1,3
               mmaxr(k)=-1
               do i=1,npts
                  mmaxr(k)=max0(mmaxr(k),r(k,i))
               end do
            end do
            write(8,'(i5,3f12.7,3i6,110i5)') niter,(theta(k),k=1,3),
     &        (mmaxr(k),k=1,3),(x(1,i-1),x(2,i-1),r(1,i),r(2,i),r(3,i),
     &        i=1,min0(19,npts)),x(1,min0(19,npts)),x(2,min0(19,npts))
         end do
      end if
      if (fixedpred.eq.2) then
         do niter=1,nburn
            do k=1,3
               theta(k)=gengamm(bth(k)+sumint(k),ath(k)+sumr(k))
            end do
            call fixedmidprey(npts,theta,rtype,times,r,sumr,x,int,
     &           sumint,lpopprod,ca2,cd2,cm2,ca22,cd22,cm22)
         end do
         do niter=1,nrun
            do j=1,nthin
               do k=1,3
                  theta(k)=gengamm(bth(k)+sumint(k),ath(k)+sumr(k))
               end do
               call fixedmidprey(npts,theta,rtype,times,r,sumr,x,int,
     &              sumint,lpopprod,counta2,countd2,countm2,counta22,
     &              countd22,countm22)
            end do
            do k=1,3
               mmaxr(k)=-1
               do i=1,npts
                  mmaxr(k)=max0(mmaxr(k),r(k,i))
               end do
            end do
            write(8,'(i5,3f12.7,3i6,110i5)') niter,(theta(k),k=1,3),
     &        (mmaxr(k),k=1,3),(x(1,i-1),x(2,i-1),r(1,i),r(2,i),r(3,i),
     &        i=1,min0(19,npts)),x(1,min0(19,npts)),x(2,min0(19,npts))
         end do
      end if
      if (fixedpred.eq.3) then
         do niter=1,nburn
            do k=1,3
               theta(k)=gengamm(bth(k)+sumint(k),ath(k)+sumr(k))
            end do
            call fixedmidprey(npts,theta,rtype,times,r,sumr,x,int,
     &           sumint,lpopprod,ca2,cd2,cm2,ca22,cd22,cm22)
            call rkmoves(npts,3,theta,rtype,times,r,sumr,x,int,
     &           sumint,lpopprod,ca3,cd3,cm3)
         end do
         do niter=1,nrun
            do j=1,nthin
               do k=1,3
                  theta(k)=gengamm(bth(k)+sumint(k),ath(k)+sumr(k))
               end do
               call fixedmidprey(npts,theta,rtype,times,r,sumr,x,int,
     &              sumint,lpopprod,counta2,countd2,countm2,counta22,
     &              countd22,countm22)
               call rkmoves(npts,3,theta,rtype,times,r,sumr,x,int,
     &              sumint,lpopprod,counta3,countd3,countm3)
            end do
            do k=1,3
               mmaxr(k)=-1
               do i=1,npts
                  mmaxr(k)=max0(mmaxr(k),r(k,i))
               end do
            end do
            write(8,'(i5,3f12.7,3i6,110i5)') niter,(theta(k),k=1,3),
     &        (mmaxr(k),k=1,3),(x(1,i-1),x(2,i-1),r(1,i),r(2,i),r(3,i),
     &        i=1,min0(19,npts)),x(1,min0(19,npts)),x(2,min0(19,npts))
         end do
      end if
      if (fixedpred.eq.4) then
         do niter=1,nburn
            do k=1,3
               theta(k)=gengamm(bth(k)+sumint(k),ath(k)+sumr(k))
            end do
            call rkmoves(npts,3,theta,rtype,times,r,sumr,x,int,
     &           sumint,lpopprod,ca3,cd3,cm3)
         end do
         do niter=1,nrun
            do j=1,nthin
               do k=1,3
                  theta(k)=gengamm(bth(k)+sumint(k),ath(k)+sumr(k))
               end do
               call allfixed(npts,theta,rtype,times,r,sumr,x,int,
     &              sumint,lpopprod,counta,countd,countm)
               call rkmoves(npts,3,theta,rtype,times,r,sumr,x,int,
     &              sumint,lpopprod,counta3,countd3,countm3)
            end do
            do k=1,3
               mmaxr(k)=-1
               do i=1,npts
                  mmaxr(k)=max0(mmaxr(k),r(k,i))
               end do
            end do
            write(8,'(i5,3f12.7,3i6,110i5)') niter,(theta(k),k=1,3),
     &        (mmaxr(k),k=1,3),(x(1,i-1),x(2,i-1),r(1,i),r(2,i),r(3,i),
     &        i=1,min0(19,npts)),x(1,min0(19,npts)),x(2,min0(19,npts))
         end do
      end if
      close(8)
      print '(a9,3f10.5)','singles: ',float(counta)/(nrun*nthin*npts),
     &     float(countd)/(nrun*nthin*npts),
     &     float(countm)/(nrun*nthin*npts)
      print '(a9,6f10.5)','doubles: ',
     &     float(counta2)/(nrun*nthin*(npts-1)),
     &     float(counta22)/(nrun*nthin*(npts-1)),
     &     float(countd2)/(nrun*nthin*(npts-1)),
     &     float(countd22)/(nrun*nthin*(npts-1)),
     &     float(countm2)/(nrun*nthin*(npts-1)),
     &     float(countm22)/(nrun*nthin*(npts-1))
      print '(a9,3f10.5)','end: ',float(counta3)/(nrun*nthin),
     &     float(countd3)/(nrun*nthin),float(countm3)/(nrun*nthin)
      end program prey
      
      subroutine rkmoves(i,k,theta,rtype,times,r,sumr,x,integral,
     &     sumint,lpopprod,ca3,cd3,cm3)
      external genunf
      integer rmax,nmax
      parameter (rmax=10000, nmax=251)
      integer rtype(nmax,rmax),r(3,nmax),i,x(2,0:nmax),sumr(3),k,
     &     ca3,cd3,cm3
      double precision times(nmax,rmax),theta(3),genunf,u,
     &     lpopprod(nmax),integral(3,nmax),sumint(3)
      u=genunf(0d0,1d0)
      if (u.lt.0.3d0) then
         call addsingle(i,k,theta,rtype,times,r,x,integral,lpopprod,
     &        sumr,sumint,ca3)
      else if (u.lt.0.6d0) then
         call delsingle(i,k,theta,rtype,times,r,x,integral,lpopprod,
     &        sumr,sumint,cd3)
      else
         call shiftsingle(i,k,theta,rtype,times,r,x,integral,lpopprod,
     &        sumint,cm3)
      end if
      end subroutine rkmoves
      
      subroutine fixedmidprey(npts,theta,rtype,times,r,sumr,x,integral,
     &     sumint,lpopprod,ca2,cd2,cm2,ca22,cd22,cm22)
      external genunf
      integer rmax,nmax
      parameter (rmax=10000, nmax=251)
      integer rtype(nmax,rmax),r(3,nmax),i,x(2,0:nmax),sumr(3),npts,
     &     ca2,cd2,cm2,ca22,cd22,cm22
      double precision times(nmax,rmax),theta(3),genunf,u,
     &     lpopprod(nmax),integral(3,nmax),sumint(3)
      do i=2,npts
         u=genunf(0d0,1d0)
         if (u.lt.0.3d0) then
            call add2move(i,theta,rtype,times,r,x,integral,lpopprod,
     &           sumr,sumint,ca2,ca22)
         else if (u.lt.0.6d0) then
            call del2move(i,theta,rtype,times,r,x,integral,lpopprod,
     &           sumr,sumint,cd2,cd22)
         else
            call shift2move(i,theta,rtype,times,r,x,integral,lpopprod,
     &           sumint,cm2,cm22)
         end if
      end do
      end subroutine fixedmidprey


      subroutine allfixed(npts,theta,rtype,times,r,sumr,x,integral,
     &     sumint,lpopprod,counta,countd,countm)
      external genunf
      integer rmax,nmax
      parameter (rmax=10000, nmax=251)
      integer rtype(nmax,rmax),r(3,nmax),i,x(2,0:nmax),sumr(3),npts,
     &     counta,countd,countm
      double precision times(nmax,rmax),theta(3),genunf,u,
     &     lpopprod(nmax),integral(3,nmax),sumint(3)
      do i=1,npts
         u=genunf(0d0,1d0)
         if (u.lt.0.3d0) then
            call addmove(i,theta,rtype,times,r,sumr,x,integral,sumint,
     &           lpopprod,counta)
         else if (u.lt.0.6d0) then
            call delmove(i,theta,rtype,times,r,sumr,x,integral,sumint,
     &           lpopprod,countd)
         else
            call shiftmove(i,theta,rtype,times,r,x,integral,sumint,
     &           lpopprod,countm)
         end if
      end do
      end subroutine allfixed


      subroutine add2move(i,theta,rtype,times,r,x,integral,
     &     lpopprod,sumr,sumint,count,count2)
      external genunf
      integer rmax,nmax
      parameter (rmax=10000, nmax=251)
      integer rtype(nmax,rmax),rtypes1(rmax),numi,i,j,k,numim1,count,
     &     count2,r(3,nmax),rtypes2(rmax),x(2,0:nmax),ps1(2,0:rmax),
     &     add1(3),add2(3),sumr(3),ps2(2,0:rmax)
      double precision times(nmax,rmax),timess1(rmax),theta(3),jtlr2,
     &     sum2,genunf,jtlr1,lratio,u,lpopprod(nmax),int(3),intm(3),
     &     integral(3,nmax),sum2m,timess2(rmax),newtime,sumint(3)
      logical feasible
c     print *,'---------------------- In add2move'
      numim1=r(1,i-1)+r(2,i-1)+r(3,i-1)
      numi=r(1,i)+r(2,i)+r(3,i)
      do j=1,numim1
         timess1(j)=times(i-1,j)
         rtypes1(j)=rtype(i-1,j)
      end do
      do j=1,numi
         timess2(j)=times(i,j)
         rtypes2(j)=rtype(i,j)
      end do
c     insert reaction types
      do k=1,3
         add1(k)=0
         add2(k)=0
         newtime=genunf(0d0,2d0)
         if (newtime.lt.1d0) then
            call insert_event(k,newtime,numim1,timess1,rtypes1)
            add1(k)=1
         else
            call insert_event(k,newtime-1,numi,timess2,rtypes2)
            add2(k)=1
         end if
      end do
c  calculate popn sizes in first interval
      ps1(1,0)=x(1,i-2)
      ps1(2,0)=x(2,i-2)
      call calpop(numim1,rtypes1,timess1,ps1,sum2m,intm,feasible)
      if (.not.feasible) then
c         print *,'Not feasible in 1st interval'
         goto 1
      end if
      if (ps1(1,numim1).ne.x(1,i-1)) then
c         print *,'add2move: prey mismatch at end of 1st interval',
c     &        ps1(1,numim1),x(1,i-1)
         goto 1
      end if
c  calculate popn sizes in second interval
      ps2(1,0)=x(1,i-1)
      ps2(2,0)=ps1(2,numim1)
      call calpop(numi,rtypes2,timess2,ps2,sum2,int,feasible)
      if (.not.feasible) then
c         print *,'Not feasible in 2nd interval'
         goto 1
      end if
      if (ps2(1,numi).ne.x(1,i)) then
         print *,'add2move: prey mismatch at end of 2nd interval',
     &        ps2(1,numi),x(1,i)
         stop
      end if
      if (ps2(2,numi).ne.x(2,i)) then
         print *,'add2move: predator mismatch at end of 2nd interval',
     &        ps2(2,numi),x(2,i)
         stop
      end if
c
      jtlr1=sum2m-lpopprod(i-1)
      jtlr2=sum2-lpopprod(i)
      do k=1,3
         jtlr1=jtlr1-theta(k)*(intm(k)-integral(k,i-1))
         jtlr2=jtlr2-theta(k)*(int(k)-integral(k,i))
      end do
c
      lratio=jtlr1+jtlr2+dlog(8d0)
      do k=1,3
         lratio=lratio+dlog(theta(k))-dlog(r(k,i-1)+r(k,i)+1d0)
      end do
      u=genunf(0d0,1d0)
      count2=count2+1
      if (dlog(u).lt.lratio) then
         count=count+1
         do j=1,numim1
            times(i-1,j)=timess1(j)
            rtype(i-1,j)=rtypes1(j)
         end do
         x(2,i-1)=ps2(2,0)
         do j=1,numi
            times(i,j)=timess2(j)
            rtype(i,j)=rtypes2(j)
         end do
         lpopprod(i-1)=sum2m
         lpopprod(i)=sum2
         do k=1,3
            r(k,i-1)=r(k,i-1)+add1(k)
            r(k,i)=r(k,i)+add2(k)
            sumr(k)=sumr(k)+1
            sumint(k)=sumint(k)+intm(k)-integral(k,i-1)
     &           +int(k)-integral(k,i)
            integral(k,i-1)=intm(k)
            integral(k,i)=int(k)
         end do
      end if
 1    end subroutine add2move
      
      subroutine del2move(i,theta,rtype,times,r,x,integral,lpopprod,
     &     sumr,sumint,count,count2)
      external genunf,ignuin
      integer rmax,nmax
      parameter (rmax=10000, nmax=251)
      integer rtype(nmax,rmax),rtypes1(rmax),numi,i,numim1,r(3,nmax),
     &     count,count2,j,rtypes2(rmax),x(2,0:nmax),ps1(2,0:rmax),
     &     ps2(2,0:rmax),ignuin,k,kill,add1(3),add2(3),sumr(3)
      double precision times(nmax,rmax),timess1(rmax),theta(3),
     &     jtlr2,sum2,genunf,jtlr1,lratio,u,lpopprod(nmax),int(3),
     &     intm(3),integral(3,nmax),sum2m,timess2(rmax),sumint(3)
      logical feasible
c     print *,'---------------------- In del2move'
      if ((r(1,i-1)+r(1,i).eq.0).or.(r(2,i-1)+r(2,i).eq.0)
     &     .or.(r(3,i-1)+r(3,i).eq.0)) goto 1
      numim1=r(1,i-1)+r(2,i-1)+r(3,i-1)
      numi=r(1,i)+r(2,i)+r(3,i)
      do j=1,numim1
         timess1(j)=times(i-1,j)
         rtypes1(j)=rtype(i-1,j)
      end do
      do j=1,numi
         timess2(j)=times(i,j)
         rtypes2(j)=rtype(i,j)
      end do
c     delete reaction types
      do k=1,3
         kill=ignuin(1,r(k,i-1)+r(k,i))
         if (kill.le.r(k,i-1)) then
            call delete_event(k,numim1,timess1,rtypes1)
            add1(k)=-1
            add2(k)=0
         else
            call delete_event(k,numi,timess2,rtypes2)
            add1(k)=0
            add2(k)=-1
         end if
      end do
c  calculate popn sizes in first interval
      ps1(1,0)=x(1,i-2)
      ps1(2,0)=x(2,i-2)
      call calpop(numim1,rtypes1,timess1,ps1,sum2m,intm,feasible)
      if (.not.feasible) then
c         print *,'Not feasible in 1st interval'
         goto 1
      end if
      if (ps1(1,numim1).ne.x(1,i-1)) then
c         print *,'del2move: prey mismatch at end of 1st interval',
c     &        ps1(1,numim1),x(1,i-1)
         goto 1
      end if
c  calculate popn sizes in second interval
      ps2(1,0)=x(1,i-1)
      ps2(2,0)=ps1(2,numim1)
      call calpop(numi,rtypes2,timess2,ps2,sum2,int,feasible)
      if (.not.feasible) then
c         print *,'Not feasible in 2nd interval'
         goto 1
      end if
      if (ps2(1,numi).ne.x(1,i)) then
         print *,'del2move: prey mismatch at end of 2nd interval',
     &        ps2(1,numi),x(1,i)
         stop
      end if
      if (ps2(2,numi).ne.x(2,i)) then
         print *,'del2move: predator mismatch at end of 2nd interval',
     &        ps2(2,numi),x(2,i)
         stop
      end if
c
      jtlr1=sum2m-lpopprod(i-1)
      jtlr2=sum2-lpopprod(i)
      do k=1,3
         jtlr1=jtlr1-theta(k)*(intm(k)-integral(k,i-1))
         jtlr2=jtlr2-theta(k)*(int(k)-integral(k,i))
      end do
c
      lratio=jtlr1+jtlr2-dlog(8d0)
      do k=1,3
         lratio=lratio-dlog(theta(k))+dlog(r(k,i-1)+r(k,i)+0d0)
      end do
      u=genunf(0d0,1d0)
      count2=count2+1
      if (dlog(u).lt.lratio) then
         count=count+1
         do j=1,numim1
            times(i-1,j)=timess1(j)
            rtype(i-1,j)=rtypes1(j)
         end do
         x(2,i-1)=ps2(2,0)
         do j=1,numi
            times(i,j)=timess2(j)
            rtype(i,j)=rtypes2(j)
         end do
         lpopprod(i-1)=sum2m
         lpopprod(i)=sum2
         do k=1,3
            r(k,i-1)=r(k,i-1)+add1(k)
            r(k,i)=r(k,i)+add2(k)
            sumr(k)=sumr(k)-1
            sumint(k)=sumint(k)+intm(k)-integral(k,i-1)
     &           +int(k)-integral(k,i)
            integral(k,i-1)=intm(k)
            integral(k,i)=int(k)
         end do
      end if
 1    end subroutine del2move
      
      subroutine shift2move(i,theta,rtype,times,r,x,integral,lpopprod,
     &     sumint,count,count2)
      external genunf,ignuin
      integer rmax,nmax
      parameter (rmax=10000, nmax=251)
      integer rtype(nmax,rmax),rtypes1(rmax),numi,i,j,k,numim1,
     &     count,count2,rtypes2(rmax),x(2,0:nmax),ps1(2,0:rmax),
     &     ps2(2,0:rmax),movertype,ignuin,kill,add1(3),add2(3),r(3,nmax)
      double precision times(nmax,rmax),timess1(rmax),theta(3),jtlr2,
     &     sum2,genunf,jtlr1,lratio,u,lpopprod(nmax),int(3),intm(3),
     &     integral(3,nmax),sum2m,timess2(rmax),newtime,sumint(3)
      logical feasible
c     print *,'---------------------- In shift2move'
      numim1=r(1,i-1)+r(2,i-1)+r(3,i-1)
      numi=r(1,i)+r(2,i)+r(3,i)
      if (numi+numim1.eq.0) then
         print *,'HELLO'
         goto 1
      end if
      do k=1,3
         add1(k)=0
         add2(k)=0
      end do
c     delete and insert reaction type
      do j=1,numim1
         timess1(j)=times(i-1,j)
         rtypes1(j)=rtype(i-1,j)
      end do
      do j=1,numi
         timess2(j)=times(i,j)
         rtypes2(j)=rtype(i,j)
      end do
      kill=ignuin(1,numi+numim1)
      if (kill.le.numim1) then
         movertype=rtype(i-1,kill)
         call delete_event(movertype,numim1,timess1,rtypes1)
         add1(movertype)=add1(movertype)-1
      else
         movertype=rtype(i,kill-numim1)
         call delete_event(movertype,numi,timess2,rtypes2)
         add2(movertype)=add2(movertype)-1
      end if
      newtime=genunf(0d0,2d0)
      if (newtime.lt.1d0) then
         call insert_event(movertype,newtime,numim1,timess1,rtypes1)
         add1(movertype)=add1(movertype)+1
      else
         call insert_event(movertype,newtime-1,numi,timess2,rtypes2)
         add2(movertype)=add2(movertype)+1
      end if
c  calculate popn sizes in first interval
      ps1(1,0)=x(1,i-2)
      ps1(2,0)=x(2,i-2)
      call calpop(numim1,rtypes1,timess1,ps1,sum2m,intm,feasible)
      if (.not.feasible) then
c         print *,'Not feasible in 1st interval'
         goto 1
      end if
      if (ps1(1,numim1).ne.x(1,i-1)) then
c         print *,'shift2move: prey mismatch at end of 1st interval',
c     &        ps1(1,numim1),x(1,i-1)
         goto 1
      end if
c  calculate popn sizes in second interval
      ps2(1,0)=x(1,i-1)
      ps2(2,0)=ps1(2,numim1)
      call calpop(numi,rtypes2,timess2,ps2,sum2,int,feasible)
      if (.not.feasible) then
c         print *,'Not feasible in 2nd interval'
         goto 1
      end if
      if (ps2(1,numi).ne.x(1,i)) then
         print *,'shift2move: prey mismatch at end of 2nd interval',
     &        ps2(1,numi),x(1,i)
         stop
      end if
      if (ps2(2,numi).ne.x(2,i)) then
         print *,'shift2move: predator mismatch at end of 2nd interval',
     &        ps2(2,numi),x(2,i)
         stop
      end if
c
      jtlr1=sum2m-lpopprod(i-1)
      jtlr2=sum2-lpopprod(i)
      do k=1,3
         jtlr1=jtlr1-theta(k)*(intm(k)-integral(k,i-1))
         jtlr2=jtlr2-theta(k)*(int(k)-integral(k,i))
      end do
c
      lratio=jtlr1+jtlr2
      u=genunf(0d0,1d0)
      count2=count2+1
      if (dlog(u).lt.lratio) then
         count=count+1
         do j=1,numim1
            times(i-1,j)=timess1(j)
            rtype(i-1,j)=rtypes1(j)
         end do
         x(2,i-1)=ps2(2,0)
         do j=1,numi
            times(i,j)=timess2(j)
            rtype(i,j)=rtypes2(j)
         end do
         lpopprod(i-1)=sum2m
         lpopprod(i)=sum2
         do k=1,3
            r(k,i-1)=r(k,i-1)+add1(k)
            r(k,i)=r(k,i)+add2(k)
            sumint(k)=sumint(k)+intm(k)-integral(k,i-1)
     &           +int(k)-integral(k,i)
            integral(k,i-1)=intm(k)
            integral(k,i)=int(k)
         end do
      end if
 1    end subroutine shift2move
      
      subroutine shiftmove(i,theta,rtype,times,r,x,integral,sumint,
     &     lpopprod,count)
      external genunf,ignuin
      integer rmax,nmax
      parameter (rmax=10000, nmax=251)
      integer rtype(nmax,rmax),rtypes(rmax),numi,i,j,k,r(3,nmax),
     &     x(2,0:nmax),ps(2,0:rmax),movertype,ignuin,count
      double precision times(nmax,rmax),timess(rmax),theta(3),newtime,
     &     sum2,genunf,jtlr,lratio,u,lpopprod(nmax),int(3),sumint(3),
     &     integral(3,nmax)
      logical feasible
c      print *,'In shiftmove'
      numi=r(1,i)+r(2,i)+r(3,i)
      if (numi.eq.0) goto 1
      if (numi.eq.1) then
         timess(1)=genunf(0d0,1d0)
         rtypes(1)=rtype(i,1)
      else 
         movertype=rtype(i,ignuin(1,numi))
c delete and insert reaction type
         do j=1,numi
            timess(j)=times(i,j)
            rtypes(j)=rtype(i,j)
         end do
         call delete_event(movertype,numi,timess,rtypes)
         newtime=genunf(0d0,1d0)
         call insert_event(movertype,newtime,numi,timess,rtypes)
      end if
c calculate popn sizes
      ps(1,0)=x(1,i-1)
      ps(2,0)=x(2,i-1)
      call calpop(numi,rtypes,timess,ps,sum2,int,feasible)
      if (.not.feasible) goto 1
c
      jtlr=sum2-lpopprod(i)
      do k=1,3
         jtlr=jtlr-theta(k)*(int(k)-integral(k,i))
      end do
      lratio=jtlr
      u=genunf(0d0,1d0)
      if (dlog(u).lt.lratio) then
         count=count+1
         do j=1,numi
            times(i,j)=timess(j)
            rtype(i,j)=rtypes(j)
         end do
         lpopprod(i)=sum2
         do k=1,3
            sumint(k)=sumint(k)+int(k)-integral(k,i)
            integral(k,i)=int(k)
         end do
      end if
 1    end subroutine shiftmove

      subroutine delmove(i,theta,rtype,times,r,sumr,x,integral,sumint,
     &     lpopprod,count)
      external genunf
      integer rmax,nmax
      parameter (rmax=10000, nmax=251)
      integer rtype(nmax,rmax),rtypes(rmax),i,j,numi,k,r(3,nmax),
     &     x(2,0:nmax),ps(2,0:rmax),count,sumr(3)
      double precision times(nmax,rmax),timess(rmax),theta(3),sum2,jtlr,
     &     genunf,lratio,u,lpopprod(nmax),int(3),integral(3,nmax),
     &     sumint(3)
      logical feasible
c      print *,'In delmove'
      if ((r(1,i).eq.0).or.(r(2,i).eq.0).or.(r(3,i).eq.0)) goto 1
      numi=r(1,i)+r(2,i)+r(3,i)
      do j=1,numi
         timess(j)=times(i,j)
         rtypes(j)=rtype(i,j)
      end do
c delete reaction types
      do k=1,3
         call delete_event(k,numi,timess,rtypes)
      end do
c calculate popn sizes
      ps(1,0)=x(1,i-1)
      ps(2,0)=x(2,i-1)
      call calpop(numi,rtypes,timess,ps,sum2,int,feasible)
      if (.not.feasible) goto 1
c
      jtlr=sum2-lpopprod(i)
      do k=1,3
         jtlr=jtlr-theta(k)*(int(k)-integral(k,i))-dlog(theta(k))
      end do
c
      lratio=jtlr
      do k=1,3
         lratio=lratio+dlog(r(k,i)+0d0)
      end do
c
      u=genunf(0d0,1d0)
      if (dlog(u).lt.lratio) then
         count=count+1
         do j=1,numi
            times(i,j)=timess(j)
            rtype(i,j)=rtypes(j)
         end do
         lpopprod(i)=sum2
         do k=1,3
            r(k,i)=r(k,i)-1
            sumr(k)=sumr(k)-1
            sumint(k)=sumint(k)+int(k)-integral(k,i)
            integral(k,i)=int(k)
         end do
      end if
 1    end subroutine delmove
      
      subroutine addmove(i,theta,rtype,times,r,sumr,x,integral,sumint,
     &     lpopprod,count)
      external genunf
      integer rmax,nmax
      parameter (rmax=10000, nmax=251)
      integer rtype(nmax,rmax),r(3,nmax),rtypes(rmax),i,j,k,numi,
     &     x(2,0:nmax),ps(2,0:rmax),count,sumr(3)
      double precision times(nmax,rmax),timess(rmax),theta(3),
     &     newtime,sum2,genunf,jtlr,lratio,u,lpopprod(nmax),int(3),
     &     integral(3,nmax),sumint(3)
      logical feasible
c      print *,'In addmove'
      numi=r(1,i)+r(2,i)+r(3,i)
      do j=1,numi
         timess(j)=times(i,j)
         rtypes(j)=rtype(i,j)
      end do
c insert reaction types
      do k=1,3
         newtime=genunf(0d0,1d0)
         call insert_event(k,newtime,numi,timess,rtypes)
      end do
c calculate popn sizes
      ps(1,0)=x(1,i-1)
      ps(2,0)=x(2,i-1)
      call calpop(numi,rtypes,timess,ps,sum2,int,feasible)
      if (.not.feasible) goto 1
c
      jtlr=sum2-lpopprod(i)
      do k=1,3
         jtlr=jtlr-theta(k)*(int(k)-integral(k,i))+dlog(theta(k))
      end do
c
      lratio=jtlr
      do k=1,3
         lratio=lratio-dlog(r(k,i)+1d0)
      end do
c 
      u=genunf(0d0,1d0)
      if (dlog(u).lt.lratio) then
         count=count+1
         do j=1,numi
            times(i,j)=timess(j)
            rtype(i,j)=rtypes(j)
         end do
         lpopprod(i)=sum2
         do k=1,3
            r(k,i)=r(k,i)+1
            sumr(k)=sumr(k)+1
            sumint(k)=sumint(k)+int(k)-integral(k,i)
            integral(k,i)=int(k)
         end do
      end if
 1    end subroutine addmove
      
      subroutine addsingle(i,k,theta,rtype,times,r,x,integral,lpopprod,
     &     sumr,sumint,count)
      external genunf
      integer rmax,nmax
      parameter (rmax=10000, nmax=251)
      integer rtype(nmax,rmax),r(3,nmax),rtypes(rmax),i,j,k,numi,
     &     x(2,0:nmax),ps(2,0:rmax),count,sumr(3),kk
      double precision times(nmax,rmax),timess(rmax),theta(3),
     &     newtime,sum2,genunf,jtlr,lratio,u,lpopprod(nmax),int(3),
     &     integral(3,nmax),sumint(3)
      logical feasible
c      print *,'In add3move'
      numi=r(1,i)+r(2,i)+r(3,i)
      do j=1,numi
         timess(j)=times(i,j)
         rtypes(j)=rtype(i,j)
      end do
c insert reaction type k
      newtime=genunf(0d0,1d0)
      call insert_event(k,newtime,numi,timess,rtypes)
c calculate popn sizes
      ps(1,0)=x(1,i-1)
      ps(2,0)=x(2,i-1)
      call calpop(numi,rtypes,timess,ps,sum2,int,feasible)
      if (.not.feasible) goto 1
c
      jtlr=sum2-lpopprod(i)
      do kk=1,3
         jtlr=jtlr-theta(kk)*(int(kk)-integral(kk,i))
      end do
c
      lratio=jtlr+dlog(theta(k))-dlog(r(k,i)+1d0)
c 
      u=genunf(0d0,1d0)
      if (dlog(u).lt.lratio) then
         count=count+1
         do j=1,numi
            times(i,j)=timess(j)
            rtype(i,j)=rtypes(j)
         end do
         lpopprod(i)=sum2
         r(k,i)=r(k,i)+1
         sumr(k)=sumr(k)+1
         x(2,i)=ps(2,numi)
         do kk=1,3
            sumint(kk)=sumint(kk)+int(kk)-integral(kk,i)
            integral(kk,i)=int(kk)
         end do
      end if
 1    end subroutine addsingle
      
      subroutine delsingle(i,k,theta,rtype,times,r,x,integral,lpopprod,
     &     sumr,sumint,count)
      external genunf
      integer rmax,nmax
      parameter (rmax=10000, nmax=251)
      integer rtype(nmax,rmax),rtypes(rmax),i,j,numi,k,r(3,nmax),
     &     x(2,0:nmax),ps(2,0:rmax),count,sumr(3),kk
      double precision times(nmax,rmax),timess(rmax),theta(3),sum2,jtlr,
     &     genunf,lratio,u,lpopprod(nmax),int(3),integral(3,nmax),
     &     sumint(3)
      logical feasible
c      print *,'In del3move'
      if (r(k,i).eq.0) goto 1
      numi=r(1,i)+r(2,i)+r(3,i)
      do j=1,numi
         timess(j)=times(i,j)
         rtypes(j)=rtype(i,j)
      end do
c delete reaction type k
      call delete_event(k,numi,timess,rtypes)
c calculate popn sizes
      ps(1,0)=x(1,i-1)
      ps(2,0)=x(2,i-1)
      call calpop(numi,rtypes,timess,ps,sum2,int,feasible)
      if (.not.feasible) goto 1
c
      jtlr=sum2-lpopprod(i)
      do kk=1,3
         jtlr=jtlr-theta(kk)*(int(kk)-integral(kk,i))
      end do
c
      lratio=jtlr-dlog(theta(k))+dlog(r(k,i)+0d0)
c
      u=genunf(0d0,1d0)
      if (dlog(u).lt.lratio) then
         count=count+1
         do j=1,numi
            times(i,j)=timess(j)
            rtype(i,j)=rtypes(j)
         end do
         lpopprod(i)=sum2
         r(k,i)=r(k,i)-1
         sumr(k)=sumr(k)-1
         x(2,i)=ps(2,numi)
         do kk=1,3
            sumint(kk)=sumint(kk)+int(kk)-integral(kk,i)
            integral(kk,i)=int(kk)
         end do
      end if
 1    end subroutine delsingle
      
      subroutine shiftsingle(i,k,theta,rtype,times,r,x,integral,
     &     lpopprod,sumint,count)
      external genunf
      integer rmax,nmax
      parameter (rmax=10000, nmax=251)
      integer rtype(nmax,rmax),rtypes(rmax),numi,i,j,k,r(3,nmax),
     &     x(2,0:nmax),ps(2,0:rmax),count,kk
      double precision times(nmax,rmax),timess(rmax),theta(3),newtime,
     &     sum2,genunf,jtlr,lratio,u,lpopprod(nmax),int(3),sumint(3),
     &     integral(3,nmax)
      logical feasible
c      print *,'In shift3move'
      numi=r(1,i)+r(2,i)+r(3,i)
      if (r(k,i).eq.0) goto 1
c delete and insert reaction type k
      do j=1,numi
         timess(j)=times(i,j)
         rtypes(j)=rtype(i,j)
      end do
      call delete_event(k,numi,timess,rtypes)
      newtime=genunf(0d0,1d0)
      call insert_event(k,newtime,numi,timess,rtypes)
c calculate popn sizes
      ps(1,0)=x(1,i-1)
      ps(2,0)=x(2,i-1)
      call calpop(numi,rtypes,timess,ps,sum2,int,feasible)
      if (.not.feasible) goto 1
c
      jtlr=sum2-lpopprod(i)
      do kk=1,3
         jtlr=jtlr-theta(kk)*(int(kk)-integral(kk,i))
      end do
      lratio=jtlr
      u=genunf(0d0,1d0)
      if (dlog(u).lt.lratio) then
         count=count+1
         do j=1,numi
            times(i,j)=timess(j)
            rtype(i,j)=rtypes(j)
         end do
         lpopprod(i)=sum2
         do kk=1,3
            sumint(kk)=sumint(kk)+int(kk)-integral(kk,i)
            integral(kk,i)=int(kk)
         end do
      end if
 1    end subroutine shiftsingle

      subroutine delete_event(deletertype,num,times,rtype)
      external ignuin
      integer rmax
      parameter (rmax=10000)
      integer rtype(rmax),num,j,ignuin,js,deletertype,n,label(rmax)
      double precision times(rmax)
      n=0
      do j=1,num
         if (rtype(j).eq.deletertype) then
            n=n+1
            label(n)=j
         end if
      end do
      if (n.eq.0) then
         print *,'delete_event: no reactions of type',deletertype
         stop
      end if
      js=label(ignuin(1,n))
      if (js.lt.num) then
         do j=js,num-1
            times(j)=times(j+1)
            rtype(j)=rtype(j+1)
         end do
      end if
      num=num-1
 1    end subroutine delete_event
      
      subroutine insert_event(insertrtype,newtime,num,times,rtype)
      integer rmax
      parameter (rmax=10000)
      integer rtype(rmax),num,j,js,insertrtype
      double precision times(rmax),newtime
      if (num.eq.0) then
         times(1)=newtime
         rtype(1)=insertrtype
      else
         if (newtime.lt.times(1)) then
            do j=num+1,2,-1
               times(j)=times(j-1)
               rtype(j)=rtype(j-1)
            end do
            times(1)=newtime
            rtype(1)=insertrtype
         else if (newtime.gt.times(num)) then
            times(num+1)=newtime
            rtype(num+1)=insertrtype
         else
            do j=1,num
               if (newtime.lt.times(j)) exit
            end do
            js=j-1
            do j=num+1,js+2,-1
               times(j)=times(j-1)
               rtype(j)=rtype(j-1)
            end do
            times(js+1)=newtime
            rtype(js+1)=insertrtype
         end if
      end if
      num=num+1
      end subroutine insert_event
      
      subroutine calpop(numi,rtypes,timess,ps,sum2,integral,feasible)
      external pop_rate
      integer rmax
      parameter (rmax=10000)
      integer rtypes(rmax),j,k,numi,ps(2,0:rmax),pop_rate
      double precision sum2,timess(rmax),integral(3)
      logical feasible
      feasible=.true.
      sum2=0d0
      if (numi.eq.0) then
         do k=1,3
            integral(k)=pop_rate(k,ps(1,0),ps(2,0))
         end do
      else
         do j=1,numi
            if (rtypes(j).eq.1) then
               ps(1,j)=ps(1,j-1)+1
               ps(2,j)=ps(2,j-1)
            end if
            if (rtypes(j).eq.2) then
               ps(1,j)=ps(1,j-1)-1
               if (ps(1,j).eq.0) then
                  feasible=.false.
                  goto 1
               end if
               ps(2,j)=ps(2,j-1)+1
            end if
            if (rtypes(j).eq.3) then
               ps(1,j)=ps(1,j-1)
               ps(2,j)=ps(2,j-1)-1
               if (ps(2,j).eq.0) then
                  feasible=.false.
                  goto 1
               end if
            end if
         end do
         do j=1,numi
            sum2=sum2+dlog(pop_rate(rtypes(j),ps(1,j-1),ps(2,j-1))+0d0)
         end do
         do k=1,3
            integral(k)=pop_rate(k,ps(1,0),ps(2,0))*timess(1) 
     &           +pop_rate(k,ps(1,numi),ps(2,numi))*(1-timess(numi))
            do j=2,numi
               integral(k)=integral(k)
     &          +pop_rate(k,ps(1,j-1),ps(2,j-1))*(timess(j)-timess(j-1))
            end do
         end do
      end if
 1    end subroutine calpop
      
      integer function pop_rate(rtype,pop1,pop2)
      integer rtype,pop1,pop2
      if (rtype.eq.1) then 
         pop_rate=pop1
      else if (rtype.eq.2) then 
         pop_rate=pop1*pop2
      else
         pop_rate=pop2
      end if
      end function pop_rate
      
      subroutine initconfig(x1start,x1end,x2start,x2end,r1,r2,r3)
      integer x1start,x1end,x2start,x2end,r1,r2,r3,diff1,diff2
      diff1=x1end-x1start
      diff2=x2end-x2start
      r1=0
      r2=-diff1
      r3=-diff1-diff2
      if ((r2.ge.0).and.(r3.ge.0)) goto 1
      r1=diff1
      r2=0
      r3=-diff2
      if ((r1.ge.0).and.(r3.ge.0)) goto 1
      r1=diff1+diff2
      r2=diff2
      r3=0
 1    end subroutine initconfig
      
c      subroutine initconfig(x1start,x1end,x2start,x2end,r1,r2,r3)
c      integer x1start,x1end,x2start,x2end,r1,r2,r3,diff
c      diff=x1end-x1start
c      r1=max0(diff,0)
c      r2=-min0(diff,0)
c      r3=x2start+r2-x2end
c      if (r3.lt.0) then
c         r1=r1-r3
c         r2=r2-r3
c         r3=0
c      end if
c      end subroutine initconfig
      
      subroutine initialise(theta,npts,x,rtype,times,r,sumr,integral,
     &     lpopprod,sumint)
      external pop_rate,genpoi
      integer rmax,nmax
      parameter (rmax=10000, nmax=251)
      integer x(2,0:nmax),r1,r2,r3,sumr(3),genpoi,rtype(nmax,rmax),num,
     &     i,j,k,nex,pop_rate,rtypes(rmax),ps(2,0:rmax),npts,r(3,nmax)
      double precision theta(3),times(nmax,rmax),incr,int(3),
     &     integral(3,nmax),lpopprod(nmax),rate,timess(rmax),sumint(3)
      logical feasible
      do i=1,npts
         call initconfig(x(1,i-1),x(1,i),x(2,i-1),x(2,i),r1,r2,r3)
         num=r1+r2+r3
         rate=0d0
         do k=1,3
            rate=rate+theta(k)*(pop_rate(k,x(1,i-1),x(2,i-1))
     &           +pop_rate(k,x(1,i),x(2,i)))/2
         end do
         nex=max0((genpoi(rate)-num)/3,0)
         incr=1d0/(num+3*nex+1)
         do j=1,num+3*nex
            times(i,j)=j*incr
         end do
         do j=1,r1
            rtype(i,j)=1
         end do
         do j=r1+1,r1+r2
            rtype(i,j)=2
         end do
         do j=r1+r2+1,r1+r2+r3
            rtype(i,j)=3
         end do
         if (nex.gt.0) then
            do j=num+1,num+3*nex,3
               rtype(i,j)=1
               rtype(i,j+1)=2
               rtype(i,j+2)=3
            end do
         end if
         num=num+3*nex
         do j=1,num
            rtypes(j)=rtype(i,j)
            timess(j)=times(i,j)
         end do
         ps(1,0)=x(1,i-1)
         ps(2,0)=x(2,i-1)
         call calpop(num,rtypes,timess,ps,lpopprod(i),int,feasible)
         if (.not.feasible) then
            print *,'Initialisation not feasible!'
            stop
         end if
         if (ps(1,num).ne.x(1,i)) then
            print *,'prey mismatch in initialisation'
            stop
         end if
         r(1,i)=r1+nex
         r(2,i)=r2+nex
         r(3,i)=r3+nex
         do k=1,3
            integral(k,i)=int(k)
         end do
      end do
      do k=1,3
         sumr(k)=0
         sumint(k)=0d0
         do i=1,npts
            sumr(k)=sumr(k)+r(k,i)
            sumint(k)=sumint(k)+integral(k,i)
         end do
      end do
      end subroutine initialise


c     Extensive DJW edits below...
      
      
      integer function genpoi(a)
      external ignpoi
      double precision a
      real ra
      ra=a
      genpoi=ignpoi(ra)
      end function genpoi

      double precision function gengamm(a,r)
      double precision a,r,s
      real ra,rs
      ra=a
      s=1d0/r
      rs=s
      gengamm=gengam(ra,rs)
      end function gengamm

      
c      double precision function genunf(a,b)
c      external g05caf
c      double precision a,b,g05caf,dum
c      genunf=a+(b-a)*g05caf(dum)
c      end function genunf

c      integer function ignuin(a,b)
c      external g05dyf
c      integer a,b,g05dyf
c      ignuin=g05dyf(a,b)
c      end function ignuin
