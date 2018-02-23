      module pkfit
! fac, amult, icyc, chi2, chiold, npars, amat
      integer, parameter ::  npars=8 , nbkpars=2, npkpars=5
      character(len=10), dimension(npars) :: names
      data names /'Intensity ', 'Const Bk  ','Slope Bk  ','Position  ', &
     &'Width     ','Eta       ','S/L       ','D/L       '  /
      integer, parameter :: scal=1, bk1=2, bk2=3, pos=4, wid=5, eta=6,  &
     &  s_l=7, d_l=8
      real, dimension(npars) :: bvec, d, s, e, pars, minpars
      real, dimension(npars,npars) :: amat
      real, dimension(:),allocatable :: fit
      real,allocatable :: data(:,:)
      real, parameter :: facst=10., ast=10.0
      real :: ycalc, ypeak, yback, ydiff, wydiff
      real :: fac, amult, chi2, chiold, chimin
      logical :: asy=.false. ! whether or not to use asymetry
      logical :: gau=.false. ! Gaussian only - eta is fixed
      logical posonly, poswd
      integer :: icyc, ndata ! number of data points
      contains
      subroutine backfunction(x,phi,pars,n,dpars)
      real, intent(in) :: x     ! two theta value
      real, intent(out):: phi   ! function value
      integer, intent(in):: n   ! Number of parameters
      real, intent(in),  dimension(n) :: pars  ! parameters
      real, intent(out), dimension(n) :: dpars 
      integer :: i                       ! derivs of phi wrt parameters
      phi=0.0; dpars=0.0 ! set out args to zero
      do i=1,n
        phi=phi+pars(i)*x**(i-1)                   ! polynomial 
        dpars(i)=x**(i-1)     ! derivatives
      enddo
      return
      end subroutine backfunction                                            
      subroutine peakfunction(x,phi,pars,n,dpars)
      real, intent(in) :: x     ! two theta value
      real, intent(out):: phi   ! function value
      integer, intent(in):: n   ! Number of parameters
      real, intent(in),  dimension(n) :: pars  ! parameters
      real, intent(out), dimension(n) :: dpars 
      real(kind=4):: Eta , Gamma , S_L , D_L , TwoTH 
      real(kind=4):: TwoTH0 , dPRdT, dPRdG, dPRdE , dPRdS , dPRdD
! also uses asy module variable
      logical Use_asym
      real(kind=4)::Profval
      external Profval
      TwoTH0=pars(1);Gamma=pars(2);Eta=pars(3)
      S_L=pars(4);D_L=pars(5); TwoTH=x
      if(asy.and.(S_L.gt.0.0 .or. D_L.gt.0.0))Use_Asym=.TRUE.
! Optimise to cut at 50 fwhm tails
      if(abs(TwoTH0-TwoTH).lt.50.0*Gamma)then
      phi = Profval(Eta,Gamma,S_L,D_L,TwoTH,TwoTH0 ,                    &
     &         dPRdT, dPRdG, dPRdE , dPRdS , dPRdD , Use_Asym)
      dpars(1)=dPRdT;dpars(2)=dPRdG;dpars(3)=dPRdE
      dpars(4)=dPRdS;dpars(5)=dPRdD
      else
       phi=0.; dpars(1:5)=0.
      endif
      return
      end subroutine peakfunction                                            
      subroutine calpat
      integer :: i, j, k
      real :: wt
      do i=1,ndata
        d=0.0 ! set derivs to zero
        call backfunction(data(1,i),yback,pars(bk1:bk2),nbkpars,d(bk1:bk2))
        call peakfunction(data(1,i),ypeak,pars(pos:d_l),npkpars,d(pos:d_l))
        ycalc = yback + pars(scal) * ypeak ! calc and diff
        d(pos:d_l)=d(pos:d_l) * pars(scal) ! correct peaks derivs for scale
        d(scal)  = ypeak                   ! deriv w.r.t to scal
        ydiff = data(2,i) - ycalc          ! obs - calc
        fit(i)=ycalc                       ! store calc
        if(gau)d(eta)=0.
        if(posonly)d(pos:d_l)=0.
        if(poswd)d(s_l:d_l)=0.
        wt=1./data(3,i)**2 
! make wt a function of difference like WIFD ??
        do j=1,npars                            ! sum in lsq contribuitions
          do k=1,npars
            amat(j,k)=amat(j,k)+d(j)*d(k)*wt
          enddo
          bvec(j)=bvec(j)+ydiff*d(j)*wt 
        enddo
        chi2=chi2+ydiff*ydiff*wt
      enddo               
      chi2=chi2/real(ndata-npars) ! make reduced chi2
      return; end subroutine calpat                                       
      ! Apply Marquardt damping (effectively observations of no change)
      subroutine marq
! fac, amult, icyc, chi2, chiold, npars, amat
      integer::i
      if(icyc.eq.0)fac=0                ! icyc=0 is last cycle flag
      if(icyc.gt.2) then
       if(chimin.gt.chi2)fac=fac/amult
       if(chimin.le.chi2 .and.icyc.gt.3)then
         fac=fac*amult**2        ! go back two lambda steps and make smaller
         amult=1.+amult/2.    ! reducing steps on lambda >= 1
       endif
      endif
      if(icyc.gt.1)then
       do i=1,npars; amat(i,i)=amat(i,i)*(1.0+fac) ; enddo
      endif 
      end subroutine marq     
! Scale matrix and take out non variable parameters
      subroutine matscal
      integer::i
      do i=1,npars
       if(amat(i,i).gt.0.) then
        d(i)=sqrt(amat(i,i))
        amat(:,i)=amat(:,i)/d(i); amat(i,:)=amat(i,:)/d(i)
        bvec(i)=bvec(i)/d(i)
       else
        d(i)=-1.0; amat(:,i)=0. ; amat(i,:)=0. ; amat(i,i)=1.0
       endif
      enddo
      end subroutine matscal

      subroutine lsqmagic()
      integer :: i, j
      allocate(fit(ndata))
      call estimatepars(pars(bk1:bk2),pars(pos:d_l),pars(scal))
      fac=facst; amult=ast; chiold=-1.; icyc=0 ; poswd=.true.; posonly=.true.
      write(*,*)fac
1     icyc=icyc+1
      amat=0.0; bvec=0.0; chi2=0.
      call calpat
      call marq
      call matscal
      call invert        ! overwrites amat with inverse
      s=0.
      do i=1,npars
        do j=1,npars; s(i)=s(i)+amat(i,j)*bvec(j) ; enddo
        if(d(i).gt.0)then ; s(i)=s(i)/d(i) ; else ; s(i)=0. ;endif
        if(d(i).lt.0.)s(i)=0.
      enddo
      do i=1,npars
        if(d(i).gt.0)then
! Scale back to original matrix
         amat(:,i)=amat(:,i)/d(i); amat(i,:)=amat(i,:)/d(i)      
        else
         amat(:,i)=0.; amat(i,:)=0.      
        endif
      enddo
      bvec=s ; e=0. ; s=0.
      do i=1,npars 
        if(d(i).gt.0.) then
         e(i)=sqrt(chi2*amat(i,i)) 
         s(i)=bvec(i)/e(i)
        endif
      enddo
      write(*,1000)icyc,chi2,fac,amult,maxval(s)
1000  format('Cycle ',i5,' Chi2 ',F15.6,3E15.6)
      if(icyc.gt.0) then ! not last cycle so update pars
         if(icyc.eq.1.and.posonly)then 
           chimin=chi2 ; minpars=pars
         endif
         chiold=chi2
         if(chi2.lt.chimin)then
           chimin=chi2
           minpars=pars ! save best parameters found
         endif
         call updatepars
         if(maxval(s).lt.0.01 .and. icyc.gt.5)then ! converged
           if(posonly)then
             posonly=.false. ; fac=facst ; amult =ast ; icyc=1
             write(*,*)'width now varying'
           else if(poswd) then
             poswd=.false. ; fac=facst ; amult = ast ; icyc=1
             write(*,*)'all parameters now varying'
           else 
             write(*,*)'Going to minpars'
             pars=minpars ! use best set found
             icyc=-1
           endif
         endif
         if(icyc.gt.100) icyc=-1 ! run out of patience
         goto 1
      else
         write(*,1001)chi2,icyc
1001     format('Perhaps converged somewhere with Chi2 = ',F15.8,1x,i5)
         write(21,*)
         write(21,'(a)')'% linecolor=3'
         do i=1,ndata
           write(21,*)data(1,i),fit(i)
         enddo
       write(21,*)
       write(21,'(a)')'% linecolor=5'
       do i=1,ndata
         write(21,*)data(1,i),data(2,i)-fit(i)
       enddo
       write(21,'(a)')'$ END'
       close(21)
       do i=1,npars
        j=i
        write(*,1002)names(i),pars(i),e(i),s(i)
       enddo
1002   format(a10,G15.8,' +/- ',G15.8,' with sh/esd =',G15.8)
      endif
      return
      end subroutine lsqmagic                                                
      subroutine estimatepars(bkpars,pkpars,scal)
      integer, parameter :: nbkpars=2, npkpars=5, npars=8 ! 1+5+2
      real, dimension(nbkpars), intent(out) :: bkpars
      real, dimension(npkpars), intent(out) :: pkpars
      real, intent(out) :: scal
      integer :: i, mxdata(1)
      real:: y,y1,y2
! background is linear a+bx+... so make a=0 and b=minimum value
      bkpars(1)=minval(data(2,1:ndata))
      bkpars(2:nbkpars)=0.0
! peak pars are tth0, width, eta, s/l and h/l
! tth0 - weighted average of data
      s=0.; sy=0. ; syt=0.; syt2=0.
      mxdata=maxloc(data(2,:))
      pkpars(1)=data(1,mxdata(1))      ! centre is max of scan
! Width from (mxdata back(1))/2  and walk array up and down to find 1/2
      y=(data(2,mxdata(1))+bkpars(1))/2.0
! walk array up
      do i=mxdata(1),ndata
        y1=data(1,i)
        if(data(2,i).lt.y) exit ! next loop
      enddo
      do i=mxdata(1),1,-1
        y2=data(1,i)
        if(data(2,i).lt.y)exit
      enddo      
      pkpars(2)=abs(y1-y2)
      write(*,*)'tth guess=' , pkpars(1),pkpars(2)
      pkpars(3)=0.9 ! guess a default eta ?
      if(gau)pkpars(3)=0.0
      if(asy)then
       pkpars(4)=0.001
       pkpars(5)=0.002
      else
       pkpars(4)=0.0
       pkpars(5)=0.0
      endif
      scal=1.
      return
      end subroutine estimatepars                                            
      subroutine invert
      !use imsl ! to get lsgrr
      real, dimension(npars,npars) :: ainv
      integer :: i
      ainv=0. ; i=1
      do i=1,5;call erset(i,1,0);enddo ! don't give up!
      call lsgrr(npars,npars,amat,npars,1.0e-5,i,ainv,npars)
!      print *,"Returned from LSGRR"
      amat=ainv
      return 
      end subroutine invert                                                  
      subroutine updatepars
       real :: step
       pars=pars+bvec
       if(pars(s_l).lt.0.)pars(s_l)=1.0e-6
       if(pars(d_l).lt.0.)pars(d_l)=1.0e-6      
       step=abs(data(1,1)-data(1,2))
       if(pars(wid).lt.step)pars(wid)=step
       if(pars(eta).lt.-0.5)pars(eta)=-0.5
       if(pars(eta).gt.2.0)pars(eta)=2.0 ! hard limits
       if(chi2.gt.chimin)then
         pars=(minpars+pars)/2.
       endif
      end subroutine updatepars                
       end module pkfit                       
      program fitit
      use pkfit
      real x,y,esd,start_time,end_time, lowtth, hightth
      integer i,j
      character(len=256) :: line
! read from stdin to a scratch file
      call cpu_time(start_time)
      i=0
      asy=.false.; lowtth=-360.0 ; hightth=360.0
      call getarg(1,line)
      read(line,*,end=20,err=20)lowtth 
      call getarg(2,line)
      read(line,*,end=20,err=20)hightth
      call getarg(3,line)
      if(line(1:2).eq.'+a') asy=.true.
      if(line(1:2).eq.'+g') gau=.true.
20    open(unit=20,status='SCRATCH',form='UNFORMATTED')
      open(unit=21,status='UNKNOWN',form='FORMATTED',                  &
     & ACCESS='SEQUENTIAL',FILE='temp.mtv')
!      open(unit=19,file='pk.dat')
!      do j=1,5
!        read(*,'(a256)')line
!        write(21,'(a)')line(1:len_trim(line))
!      enddo
1     read(*,*,end=10,err=100)x,y,esd
      if(x.gt.lowtth .and. x.lt.hightth) then
        write(20)x,y,esd
        write(21,*)x,y
        i=i+1
      endif
      goto 1
! Copy from scratch to allocatable data array
10    allocate(data(3,i),stat=j)
      ndata=i
      if(j.ne.0)stop 'Memory allocation error'
      rewind(20)
      do j=1,i
        read(20)data(1,j),data(2,j),data(3,j)
      enddo
      close(20)
      write(*,*)'No of datapoints=',ndata
! Data is read in, now do the fitting
      call lsqmagic()
      deallocate(data)
      call cpu_time(end_time)
      write(*,1000)end_time-start_time
1000  format('Time taken was ',F15.3,'/s')
      goto 101
100   STOP 'Data format problem'
101   end program fitit                                                      