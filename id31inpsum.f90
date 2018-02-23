      program id31inpsum
      implicit none
      character(len=256)::line
      integer :: ier, iarg, ioutunit, jjj
      integer, external :: iargc
      integer(kind=4) :: i, ndata, iscan, n, j, k, n3s, n6s
      real(kind=8) :: x,y,esd, diff, c2, wd, tthlow, tthhigh, step
      real(kind=8),allocatable,dimension(:,:) :: data
      integer(kind=4),external :: jbin
! Read epf file in
      i=0; ier=0
      call getarg(1,line) ! name of first file
      open(unit=20,file=line(1:len_trim(line)),status='OLD',iostat=ier)
      if(ier.ne.0)then
       write(*,'(a)')'Error opening '//line(1:len_trim(line))
       write(*,'(a)')'Usage: id31inpsum file1 file2 file3 ...'
       stop
      endif
      write(*,'(a)')'Reading'//line(1:len_trim(line))
      read(20,*,end=10,err=11)x0,y0,esd0
      read(20,*,end=10,err=11)x,y,esd      
      step=x-x0
      nbins=360.0d0/step
 10   allocate(data(3,nbins),stat=ier)
      if(ier.ne.0)stop 'Memory allocation problem'
      i=i+1
      goto 1
 11   stop 'Data format problem'
!      

      ndata=i
      rewind(21)
      do i=1,ndata
         read(21)data(1,i),data(2,i),data(3,i)
      enddo
      close(21)
      close(20)
      jjj=1; ioutunit=19
! write a file in plotmtv format with the 9 channels and sum
      open(unit=ioutunit,file='scans.mtv',status='UNKNOWN',              &
     & form='FORMATTED', access='SEQUENTIAL',iostat=ier) 
      if(ier.ne.0)stop 'error opening diagnostic file'
      write(ioutunit,'(a)')'$ DATA = CURVE2D'
      write(ioutunit,'(a)')'% xlabel = "Two Theta"'
      write(ioutunit,'(a)')'% ylabel = "Cts/Monitor"'
      write(ioutunit,'(a)')'% toplabel= "Scans and total plot"'
      write(ioutunit,'(a,i2,a)')'% linelabel = "Total" '      
      write(ioutunit,'(a,i2)')'% linecolor = ',jjj      
      write(ioutunit,'(a)')'% linetype=1 markertype=0'      
      do j=1,ndata
        write(ioutunit,'(2F15.8)')data(1,j),data(2,j)
      enddo
      write(ioutunit,*)
! Epf file is now in array data
      step=(data(1,ndata)-data(1,1))/real(ndata-1,8)
      tthlow=data(1,1)-step/2.0d0
      tthhigh=data(1,ndata)+step/2.0d0
      write(*,*)'Opened summed file ',line(1:len_trim(line)),' npts ',    &
     & ndata
      write(*,1000)'Step ',step,' tthlow ',tthlow,' tthhigh ',tthhigh
 1000 format(3(a,F15.8))
! Now read in the series of scans
      iarg=iargc() ! total number of cmdline arguments
      if(iarg.le.1)then
       stop 'need to supply a list of scans to check'
      endif
!                    1234567890123456789012345678901234567890
      write(*,999)
 999  format('      Scan      npts   3-sigma  6-sigma     chi**2')
      do i=2,iarg
       call getarg(i,line)
       line=adjustl(line)
       n=len_trim(line)
       write(line,'(a)')line(1:n)//'.inp'
       open(unit=20,file=line(1:len_trim(line)),status='OLD',iostat=ier)
       if(ier.ne.0)then
        write(*,'(a)') 'Error opening '//line(1:len_trim(line))
        stop
       endif
!       write(*,*)'Opened '//line(1:len_trim(line))
       write(ioutunit,*)
       write(ioutunit,'(a,a,a)')'% linelabel = "',                          &
     &   line(1:len_trim(line)),'"'      
       jjj=jjj+1
       write(ioutunit,'(a,i3)')'% linecolor = ',jjj      
       write(ioutunit,'(a)')'% linetype=1 markertype=0'
       c2=0.0d0; n3s=0; n6s=0; n=0
! read and check the .inp against the .epf now
 2     read(20,*,end=21,err=22)x,y,esd
         write(ioutunit,'(2F15.8)')x,y
         k=jbin(x,tthlow,tthhigh,step)
! test binning !       write(*,*)x,data(1,k)
         if(k.gt.0)then
          n=n+1
          diff=y-data(2,k)
          wd=diff*diff/(esd*esd+data(3,k)*data(3,k))
          c2=c2+wd
          if(wd.gt.9.0d0)n3s=n3s+1
          if(wd.gt.27.0d0)n6s=n6s+1
         endif
       goto 2 
 22    write(*,'(a)')'Data format problem in'//line(1:len_trim(line))
       stop
 21    close(20)
       write(*,1001)line(1:len_trim(line)),n,n3s,n6s,c2/real(n,8)
 1001  format(a10,3i10,F16.8)
      enddo
      write(ioutunit,'(a)')'$ END'
      close(ioutunit)
      end program id31check

      integer(kind=4) function jbin(tth,tthlow,tthhigh,step)
      implicit none
      real(kind=8),intent(in)::tth, tthlow, tthhigh, step
      real(kind=8)::x
      if(tth.ge.tthlow .and. tth.le.tthhigh)then
       x=(tth-tthlow)/step+0.5d0
       jbin=nint(x)
      else
       jbin=-1
      endif
      return
      end function jbin                                                     