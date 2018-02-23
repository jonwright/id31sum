      module specfiles
      integer(kind=4) :: iunit,iscan,ispecerr,ncolumns
      integer(kind=4),parameter :: LINELENGTH=512, WORDLENGTH=100
      integer(kind=4),parameter :: NWORDS=60, NLINES=50
      integer(kind=4),dimension(NLINES) :: headerwords
      real(kind=8) :: scanstart, scanend
      real(kind=8),dimension(NWORDS,NLINES) :: wordvalues
      real(kind=8),dimension(NWORDS) :: Q
      character(len=6) :: lnfmt
      character(len=LINELENGTH) :: filnam, line
      character(len=LINELENGTH) :: specsfilename, specsdate
      character(len=WORDLENGTH) :: scantype
      character(len=WORDLENGTH),dimension(NWORDS,NLINES) :: words
      character(len=WORDLENGTH),dimension(NWORDS) :: columnlabels
      integer(kind=4) :: i0, i9, id, isp, im
      character(len=WORDLENGTH) :: FIRSTDET,LASTDET,TWOTTH
      character(len=WORDLENGTH) :: MONITORCOL="Monitor"
      logical topofscan
      data lnfmt /'(a512)'/
      data iunit /15/                                                        
      contains
      subroutine getfile
      logical :: od
! See if there is already a file open, if there is then close it
      inquire(unit=iunit,opened=od)
      if(od)close(iunit)
! filnam must already be filled in, or we'll try to open garbage
1     open(unit=iunit,file=filnam,form='formatted',                     &
     & access='sequential', status='old',iostat=ispecerr)
      if(ispecerr.ne.0) then
      write(*,'(3a,i6)')'Problem with file ',filnam(1:len_trim(filnam)),&
     & ' iostat=',ispecerr
        write(*,'(a)')'Please supply a valid SPEC file name'
        read(*,lnfmt)filnam
        goto 1
      endif      
      topofscan=.false.
! Supply defaults in case dumb specfile has no #L
      ncolumns=14 ; columnlabels(1)='2_theta';  columnlabels(2)='Epoch' 
      columnlabels(3)='Seconds';   columnlabels(4)='MA0'     
      columnlabels(5)='MA1';       columnlabels(6)='MA2'     
      columnlabels(7)='MA3';       columnlabels(8)='MA4'     
      columnlabels(9)='MA5';       columnlabels(10)='MA6'     
      columnlabels(11)='MA7';      columnlabels(12)='MA8'     
      columnlabels(13)='Monitor';  columnlabels(14)='Fluo det'
      i0=ichar('0');i9=ichar('9')
      id=ichar('.');isp=ichar(' ');im=ichar('-')
      return
      end subroutine getfile                                                  
      subroutine findscan(n)
      integer(kind=4),intent(in)::n
      character(len=7) :: rd ! enough space for yes, no and unknown
      inquire(unit=iunit,read=rd)        ! Can we read the file?
      if(rd(1:3).ne.'YES') call getfile  ! If not go open it 
      if(iscan.eq.n .and. topofscan) return
      if(iscan.gt.n) rewind(n) ! 
2     call readheader ! in case top of file     
      if(iscan.ne.n)then ! work through file to find #
1      read(iunit,lnfmt,err=10,end=10)line
       if(line(1:1).eq.'#') goto 2
       go to 1
      else
! Scan found - check scantype elsewhere - not a job for specfile module
       return
      endif
10    write(*,*)
      write(*,'(a,i5,a)')'Scan ',n,' not found in file '//              &
     & filnam(1:len_trim(filnam))
      ispecerr=-2
      return
      end subroutine findscan                                                
      subroutine getdata(a,n)
! A is an array, dimensioned by the number of data items to be expected
      integer(kind=4),intent(in) :: n
      real(kind=8),dimension(n),intent(out) :: a
      integer(kind=4), parameter :: one = 1
      if(topofscan)then; topofscan=.false. ; else
        read(iunit,lnfmt,err=10,end=10)line                  ! goes via line
        topofscan=.false.
      endif
      call rdnums(one,n,a)                                   ! reads from line
      if(line(1:1).eq.'#')then; ispecerr=-1 ;goto 100 ; endif
      go to 100    
10    ispecerr=-1 ! End of file
100   return ; end subroutine getdata                                        
      subroutine readheader
      character(len=1) :: letter
      character(len=128) :: motor
      integer(kind=4) :: i,j
      real(kind=8) :: a
      integer(kind=4),parameter :: four=4
      if(line(1:1).eq.'#') goto 2 ! if already on header, don't skip
1     read(iunit,lnfmt,end=100)line
2     if(line(1:1).eq.'#') then; letter=line(2:2); topofscan=.true.
       select case((letter))
        case('O')
         read(line(3:3),'(i1)')i
         call split(line(4:LINELENGTH),words(:,i+1),WORDLENGTH,NWORDS,j)
         headerwords(i+1)=j
        case('P')
         read(line(3:3),'(i1)')i
         call rdnums(four,NWORDS,wordvalues(:,i+1))
        case('Q')
! Copy the Q from the header into our specfile module
         read(line(3:LINELENGTH),*,end=1, err=1) Q
        case('N')
         read(line(3:len_trim(line)),*,end=1)ncolumns
        case('L') ! signals end of header, bug out here !
         call split(line(3:LINELENGTH),columnlabels,WORDLENGTH,NWORDS,i)
         if(ncolumns.ne.i)                                              &
     &   write(*,'(a,i5)')'error reading header for scan',iscan 
        case('S')
         read(line(3:len_trim(line)),*,end=1)iscan,scantype,motor,      &
     &                 scanstart, scanend
        case('F')
         specsfilename=line(3:len_trim(line)) ! get original filename
        case('D')
         read(line(3:len_trim(line)),'(a256)',end=1)specsdate
        case('C'); continue ! comments
        case('G'); continue ! no idea - always zero
        case('E'); continue ! epoch - we don't care what it was for now.
        case default
       end select
      else ! Not a # line, so assume end of header
       read(line,*,err=1,end=1)a ! catch blank lines in header
       topofscan=.true.              !!! Must be able to read a number !
       return                        !!! escapes from routine here !!!!!
      endif
      goto 1
100   ispecerr=-1; return;  end subroutine readheader                        
      real(kind=8) function getheadervalue(string)
      character(len=*),intent(in) :: string
      integer(kind=4) :: i, j, k
      k=len(string); getheadervalue=0.0
      do i=1,NWORDS; do j=1,NLINES
        if(string(1:k).eq.words(i,j)(1:k))then
          getheadervalue=wordvalues(i,j); return  
        endif
      enddo; enddo
      return; end function getheadervalue                                    
      integer(kind=4) function whichcolumn(string)
! Interprets the #L line information
      character(len=*),intent(in) :: string
      integer(kind=4) :: i, j
      j=len(string) ;  whichcolumn=-1
      do i=1,NWORDS
        if(string(1:j).eq.columnlabels(i)(1:j))then
          whichcolumn=i ;  return
        endif
      enddo
      return; end function whichcolumn                                       
      subroutine rdnums(ic,n,values)
! Placed here in a subroutine in case of formatting or error handling problems
      integer(kind=4),intent(in) :: n, ic
      real(kind=8),dimension(n),intent(inout) :: values
      if(line(ic:ic).eq.'#')then; ispecerr=-1; return; endif
      read(line(ic:len_trim(line)),*,err=10,end=20)values(1:ncolumns)
      goto 100
10    write(*,*)'Error reading line, looking for ',ncolumns,' values'
      write(*,*)line(ic:len_trim(line))
      write(*,*)values
20    continue
100   return
      end subroutine rdnums  
      subroutine split(instring,outstrings,lenout,n,i)
      integer(kind=4),intent(in) :: lenout,n
      character(len=*),intent(in) :: instring
      character(len=lenout),dimension(n),intent(out) :: outstrings
      integer(kind=4),intent(out) :: i
      integer(kind=4) :: j,k,l
      j=1; k=1
      do i=1,len_trim(instring) ! hope len > len_trim or array overstep
      if(instring(i:i+1).ne.'  ')then
       outstrings(j)(k:k)=instring(i:i)
       if(k.ne.1 .and. k.lt.lenout) k=k+1
       if(instring(i:i).ne.' '.and.k.eq.1)k=k+1
      else 
       if(k.gt.1)then ; do l=k,lenout
         outstrings(j)(l:l)=' ' ! blank pad end of string
       enddo ; j=j+1; k=1;  if(j.gt.n)exit ; endif
      endif
      enddo    
! Blank pad last string if necessary
      if(k.gt.1)then
       do l=k,lenout;outstrings(j)(l:l)=' ';enddo
      endif
      i=j; return;  end subroutine split                                     
      end module specfiles                                                   
      program c2xye
      use specfiles
      integer(kind=4) :: i,j
      integer(kind=4) :: nscan, ncol
      character(LEN=WORDLENGTH) :: label
      real(kind=8),allocatable :: data(:), prevdata(:)
      real :: time1, time2
      call getarg(1,filnam)
      call getfile
      call getarg(2,line)
      read(line,*,err=10,end=10)nscan
      call findscan(nscan)
      call getarg(3,label)
      read(label,*,err=11,end=11)ncol
      goto 20
! allow labels as well as numbers for column headings
11    ncol=whichcolumn(label(1:len_trim(label)))
      if(ncol.lt.1)then
       write(0,'(a,a)')'Couldn''t find column labelled ',                 &
     &        label(1:len_trim(label))
       write(0,'(a,i5)')'Number of columns = ', ncolumns
       do j=1,ncolumns
        write(0,'(i5,1x,a)')j,columnlabels(j)
       enddo
       goto 10
      endif
      goto 20
10    write(0,*)'Probs with command line'
      goto 100
20    continue

      if(.not.allocated(data))allocate(data(ncolumns))
      if(.not.allocated(prevdata))allocate(prevdata(ncolumns))
      write(*,'(a)')'$ DATA=CURVE2D'
      write(*,'(a)')'% linetype=0 markertype=2'
      write(*,'(3a)') '% xlabel="',columnlabels(1),'"'
      write(*,'(3a)') '% ylabel="',columnlabels(ncol),'"'
      write(*,'(a)') '% title="A fit"'
1     prevdata=data                  ! Start of loop through reading data
      call getdata(data,ncolumns)
      if(ispecerr.eq.0)then
         write(*,*)data(1),data(ncol),sqrt(data(ncol)+1.) ! default weight
         goto 1
      endif
      deallocate(data); deallocate(prevdata)
      call cpu_time(time2)
100   end program c2xye                                                   