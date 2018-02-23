      program sifit
      implicit none
      integer, parameter :: npks=20
      integer, dimension(3,npks) :: hkls
      integer :: n, i, j, ncyc, nobs
      data hkls / 1,1,1, 2,2,0, 3,1,1, 4,0,0, 3,3,1,       &
     &            4,2,2, 5,1,1, 4,4,0, 5,3,1, 6,2,0,       &
     &            5,3,3, 4,4,4, 5,5,1, 6,4,2, 5,5,3,       &
     &            8,0,0, 7,3,3, 8,2,2, 5,5,5, 8,4,0        /
      real(kind=8) :: alambda, dalam, zero, dzero ! lambda, deriv
      real(kind=8), dimension(2,npks) :: peaks ! obs, esd
      real(kind=8), dimension(2,npks) :: cpks  ! calc peaks, d, tth
      real(kind=8), dimension(npks)   :: d     ! differences obs-calc
      real(kind=8), dimension(2,2)    :: a, ai ! LSQ matrix & inverse
      real(kind=8), dimension(2)      :: b, x  ! vectors in A.x=b
      real(kind=8), parameter :: alatt=5.4311946d0  ! si lattice par
      real(kind=8) :: PI, RAD ! convert deg to rad by multiplying
      real(kind=8) :: t, chi2 ! temporary
!
      peaks=0.0d0 ! initialise peaks array to zero
      cpks=0.0d0
      PI=2.0d0*datan2(1.0d0,0.0d0)
      RAD=PI/180.0d0
!      write(*,*)'RAD = ',RAD, ' PI= ',PI
      i=1
!      write(*,'(a)')'  h   k   l      obs              esd             d'   
1     read(*,*,end=2,err=2)peaks(1,i),peaks(2,i) ! fill in obs and esd
      n=hkls(1,i)*hkls(1,i)+hkls(2,i)*hkls(2,i)+hkls(3,i)*hkls(3,i)
      cpks(1,i)=1.0d0/dsqrt(real(n,8)/alatt/alatt) ! in d-spacing
!      write(*,1000)(hkls(j,i),j=1,3),(peaks(j,i),j=1,2),cpks(1,i)
1000  format(3(i3,1x),3F16.8)
      i=i+1
      if(i.le.npks)then
        goto 1
      else
        write(*,*)'Only the first ',npks,' peaks can be used'
      endif
2     nobs=i-1
! Estimate lambda from highest hkl peak with zero as 0.0
      zero=0.0d0; ncyc=0
      alambda=2.0d0*cpks(1,nobs)*dsin(RAD*peaks(1,nobs)/2.0d0)
!      write(*,*)'Initial guess lambda = ',alambda,' zero ', zero
3     a=0.0d0 ; x=0.0d0; b=0.0d0; chi2=0.0d0
! Fill in calc positions with this lambda and zero
      do i=1,nobs
        t=alambda/2.0d0/cpks(1,i)
        cpks(2,i)=2.0d0*dasin(t)/RAD - zero 
        dalam=(1.0d0/cpks(1,i))*(1.0d0/dsqrt(1.0d0-t*t))/RAD
        dzero=-1.0d0
        dalam=dalam/peaks(2,i)
        dzero=dzero/peaks(2,i)
        a(1,1)=dalam*dalam+a(1,1)                          ! fill in lsq
        a(2,2)=dzero*dzero+a(2,2)
        a(2,1)=dzero*dalam+a(2,1)
        a(1,2)=dzero*dalam+a(1,2)
        d(i)=(peaks(1,i)-cpks(2,i))/(peaks(2,i))   ! wtd difference
        b(1)=b(1)+d(i)*dalam                        ! rhs
        b(2)=b(2)+d(i)*dzero
        chi2=chi2+d(i)*d(i)
      enddo
! invert matrix, and apply shifts 
      t=1.0d0/(a(1,1)*a(2,2)-a(1,2)*a(2,1)) 
             ! determinant  - should catch singularity here
!      write(*,*)'Determinant = ',t
      ai(1,1)=a(2,2)*t
      ai(2,1)=-a(1,2)*t
      ai(1,2)=-a(2,1)*t
      ai(2,2)=a(1,1)*t
! fill in shifts in x
      x(1)=ai(1,1)*b(1)+ai(1,2)*b(2)
      x(2)=ai(2,1)*b(1)+ai(2,2)*b(2)
      alambda=alambda+x(1)
      zero=zero+x(2)
! more cycles? say 5 for now - seems to converge in 1 (is it linear?)
!      write(*,*)alambda,zero
!      write(*,*)'Cycle ',ncyc,' Chi2 ',chi2
      ncyc=ncyc+1
      if(ncyc.lt.6)goto 3
      write(*,*)' h   k   l   obs           esd         calc        '//      &
     & 'diff/esd   diff'
      do i=1,nobs
       write(*,2000)(hkls(j,i),j=1,3),peaks(1,i),peaks(2,i),cpks(2,i),d(i),  &
     & peaks(1,i)-cpks(2,i)
2000  format(3(i3,1x),3(F12.7,1x),1F8.2,1x,F12.7)
      enddo
      write(*,'(a)')'Lambda = 2 * d * sin[(two_theta + zero)/2] '
      write(*,2001)'a(Si)  ',alatt
2001  format(a7,F16.8,a,F16.8)
      write(*,2001)'Chi2   ',chi2   ,'   reduced Chi2', chi2/real(nobs-2,8) 
!
      write(*,'(a)')'+------------------------------------------------------+'
      write(*,2002)'Lambda',alambda,'   Zero',zero
      write(*,2002)'   +/-',dsqrt(chi2*ai(1,1)),'    +/-', dsqrt(chi2*ai(2,2))
2002  format('|',a7,F16.8,8x,a7,F16.8,' |')
      write(*,'(a)')'+------------------------------------------------------+'
      write(*,2003)ai(1,2)/dsqrt(ai(1,1)*ai(2,2))
2003  format('Correlation coefficient of wavelength and zero = ',F16.8)
      end program sifit                                                      