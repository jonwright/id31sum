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
      module rebin
      integer(kind=4), parameter :: NCHAN = 9  ! ID31 will have nine channels
      integer(kind=4) :: NCHANNEL = 9  ! ID31 will have nine channels
      real(kind=8), dimension(NCHAN) :: offset, mult, multerr                
      logical :: tempres 
      integer(kind=4)logexdet(NCHAN)
      integer(kind=4) :: iexrc ! region count
      integer(kind=4),allocatable, dimension(:) :: iexarray
      real(kind=8),allocatable, dimension(:,:) :: exarray                               
      data logexdet /9*0/ ! whole array initialised to false                
      data iexrc /0/    ! no excluded regions by default                  
      real(kind=8), allocatable ::  ascan(:,:)
      real(kind=8) :: step, tthlow, tthhigh, aminstep
! if requested by the user
      real(kind=8) :: user_step = 0.003d0 
      real(kind=8) :: user_tthlow = -30.0d0
      real(kind=8) :: user_tthhigh = 160.0d0
      integer(kind=4) :: npts
      data tthlow, tthhigh, step, aminstep /-30.0d0, 160.0d0, 0.003d0,  &
     & 0.0002d0 /                                                            
      real(kind=8) :: sumtotal(nchan), sumtotalmon, minmon=1.0d0
      real(kind=8) :: winhigh,winlow,winhighread,winlowread
      integer(kind=4) :: wincol
      real(kind=8) :: minrenormsig=5.0
      real(kind=8) :: randomstart = 0, randomval=0
      real(kind=8) :: randomend = 0, randomvalend=0
      character(len=90) :: wincnt
      logical :: winlog=.false. 
      logical :: userandomstart = .false.
      logical :: rstchan( nchan ) = .false.
      logical :: userandomend = .false.
      logical :: renchan( nchan ) = .false.
      logical :: wavelength_set = .false.
      real(kind=8) :: wavelength=0.0d0
      ! T = Two theta
      ! Q = 4 * pi * sin( two_theta * pi / 360.0 ) / wavelength
      ! R = [q squared] = Q * Q
      character(len=1) :: units = 'T'
      real(kind=8) :: four_pi, pi_over_360
      contains
      integer(kind=4) function ibin(tth)
      real(kind=8), intent(in) :: tth
      real(kind=8)::x
      if(tth.ge.tthlow .and. tth.le.tthhigh) then
        x=(tth-tthlow)/step+0.5d0
        ibin=int(x) 
      else ;  ibin=-1 ; endif       ! Out of range for this ascan array
      return ; end function ibin                                             
      real(kind=8) function tthhb(n)
      integer(kind=4), intent(in) :: n
      integer(kind=4)nplusone
      nplusone=n+1
      tthhb=tthlb(nplusone)
      return; end function tthhb                                             
      real(kind=8) function tthlb(n)
      integer(kind=4), intent(in) :: n
      tthlb = (real(n,8)-0.5d0)*step + tthlow
      return; end function tthlb                                      
      real(kind=8) function bincen(n)
      integer(kind=4), intent(in) :: n
      bincen = real(n,8)*step + tthlow
      return; end function bincen                                     
      subroutine checkrebinpars
      real(kind=8) :: x
      if ((units.ne.'T').and. (.not. wavelength_set) ) return
      if(step.lt.aminstep .and. (units.eq.'T'))then
       write(*,'(a,G12.4)')'Step is a bit small, resetting to ',aminstep
       step=aminstep
       user_step = step
      endif
      x=tthlow/step 
      tthlow=real(int(x),8)*step
      x=(tthhigh-tthlow)/step
      npts=int(x) 
      tthhigh=tthhb(npts)
      if (units .eq. 'T') then
      write(*,'(3(G12.5,a),G12.5)')tthlow,' < tth < ',tthhigh,            &
     &' step=',step,' npts=',npts
      endif
      if (units.eq.'Q') then
      write(*,'(3(G12.5,a),G12.5)')tthlow,' < 2pi/d < ',tthhigh,          &
     &' step=',step,' npts=',npts
      endif
      if (units.eq.'R') then
      write(*,'(3(G12.5,a),G12.5)')tthlow,' < Q^2 < ',tthhigh,            &
     &' step=',step,' npts=',npts
      endif
      return
      end subroutine checkrebinpars                                          
      subroutine initialiserebin
      integer(kind=4) :: ierr, i
      real(kind=8) :: junk
!      real :: time1, time2
! constants to machine precision
      four_pi = 8.0d0*asin(1.0d0)
      pi_over_360 = asin(1.0d0)/180.0d0
! default values
 ! Germanium numbers
      offset(1) =  7.88643510d0
      offset(2) =  5.91013889d0
      offset(3) =  3.89405184d0
      offset(4) =  1.97036951d0
      offset(5) =  0.00000000d0
      offset(6) = -2.12832415d0
      offset(7) = -4.03585040d0
      offset(8) = -6.00076222d0
      offset(9) = -8.03007953d0

 ! Silicons latest numbers
      offset(1) = 8.0634515d0
      offset(2) = 5.8865159d0
      offset(3) = 3.9594961d0
      offset(4) = 2.0986688d0
      offset(5) = 0.0000000d0
      offset(6) = -1.9488783d0
      offset(7) = -3.9966086d0
      offset(8) = -6.0474594d0
      offset(9) = -8.0536348d0

 ! Some new numbers for silicon - there seems to be some drift?
      offset(1) = 8.05624406d0
      offset(2) = 5.88332391d0
      offset(3) = 3.95657893d0
      offset(4) = 2.09530059d0
      offset(5) = 0.00000000d0
      offset(6) = -1.94883199d0
      offset(7) = -3.99982480d0
      offset(8) = -6.04725698d0
      offset(9) = -8.05585677d0

! some more new numbers for silicon - there is drift for sure

      offset(1) =  8.02703050d0
      offset(2) =  5.88348041d0
      offset(3) =  3.95722668d0
      offset(4) =  2.09585757d0
      offset(5) =  0.00000000d0
      offset(6) = -1.94681946d0
      offset(7) = -3.99878112d0
      offset(8) = -6.04566287d0
      offset(9) = -8.05515342d0


                                                                            
      open(unit=16,file='temp.res',form='FORMATTED',                    &
     & access='SEQUENTIAL', status='OLD', iostat=ierr)
      if(ierr.eq.0)then
         do i=1,nchannel
! Should clarify policy on whether to read/use these?
          read(16,*,err=10,end=12)offset(i),mult(i),multerr(i),junk
         enddo
! Report temp.res found and read in         
         write(*,'(a)')'temp.res file found and read in'
         tempres=.true.
        goto 11
      endif
10    write(*,'(a)')'Not able to read temp.res file,'                 
      write(*,'(a)')'You need this file for the detector calibration'
      write(*,'(a)')'Please copy this file to the current directory'
      write(*,'(a)')' eg:  " cp ~/temp.res ." '
      stop
12    write(*,'(a)')'Reached end of temp.res file early?'
11    close(16)
! two theta low, two theta high, step and that npts is correct
      if ((units.eq.'T') .or. wavelength_set) call unit_lims
      return ; end subroutine initialiserebin

      subroutine rebinallocate
      implicit none
      integer ierr
      if( allocated(ascan) ) deallocate(ascan) ! can reset
      allocate(ascan(nchannel*2,npts),stat=ierr)
      if(ierr.ne.0)stop 'Memory allocation error'
!      call cpu_time(time1)
      ascan=0.0d0   ! Clear any junk  
!      call cpu_time(time2)
!      write(*,*)'Time taken to zero array =',time2-time1,'/s'
      return ; end subroutine rebinallocate      
      
      subroutine bin(tth1,tth2,cts,ichan)
! tthh & tthl are the high & low two thetas from the SPEC file
! cts is the number of counts and ichan is the column to use in ascan
      real(kind=8), intent(in) :: tth1, tth2, cts
      integer(kind=4), intent(in) :: ichan
      integer(kind=4) :: ibl, ibh, j
      real(kind=8) :: frac, tthh, tthl
!      tthl=min(tth1,tth2)                       ! Ensures ascending data
!      tthh=max(tth1,tth2)    ! strange bug on Fe3O4 exp??
      if(tth1 .ge. tth2)then
        tthl=tth2; tthh=tth1
      else
        tthl=tth1; tthh=tth2
      endif
      ibl = ibin(tthl)                              ! Bin of lower point
      ibh = ibin(tthh)                              ! Bin of upper point
      if ((ibl.le.0).or.(ibh.lt.ibl).or.(ibh.gt.npts))return     ! error   
      if ( ibh .eq. ibl) then                           ! All in one bin
         ascan(ichan,ibl) = ascan(ichan,ibl) + cts
         return
      elseif ( ibh .eq. ibl+1) then               ! Spread over two bins
         frac=(tthhb(ibl)-tthl)/(tthh-tthl)          ! Fraction in bin 1
         ascan(ichan,ibl)=ascan(ichan,ibl) + cts*frac
         ascan(ichan,ibh)=ascan(ichan,ibh) + cts*(1.0d0-frac)
         return
      else
         frac=(tthhb(ibl)-tthl)/(tthh-tthl)          ! Fraction in bin 1
         ascan(ichan,ibl)=ascan(ichan,ibl) + cts*frac
         frac=(tthh-tthlb(ibh))/(tthh-tthl)       ! Fraction in last bin 
         ascan(ichan,ibh)=ascan(ichan,ibh) + cts*frac 
         frac= step/(tthh-tthl)                ! Fraction in middle bins
         do j = ibl+1,ibh-1,1
           ascan(ichan,j)=ascan(ichan,j) + cts*frac 
         enddo
         return
      endif
      end subroutine bin                                                     
      subroutine processline(tth1,tth2,cts,n,ctsmon)
      real(kind=8), intent(in) :: tth1, tth2, ctsmon
      integer(kind=4), intent(in) :: n
      real(kind=8), intent(in), dimension(n) :: cts
      integer(kind=4) :: i, ich
      real(kind=8) :: tthh, tthl
      do i = 1,n
        if(logexdet(i).eq.1)cycle ! skip if excluded completely
        if(userandomstart.and.rstchan(i).and.(tth2.lt.randomval))cycle
        if(userandomend.and.renchan(i).and.(tth2.gt.randomvalend))cycle
        tthl = tth1 - offset(i)
        tthh = tth2 - offset(i)
        if(logexdet(i).eq.2)then
! Is this an excluded region for this channel
          if(led(i,tthl,tthh))cycle  
        endif
        tthl = convert_unit_function( tthl )
        tthh = convert_unit_function( tthh )
        ich = 2*i-1
        call bin(tthl,tthh,cts(i),ich)
        ich = 2*i
        call bin(tthl,tthh,ctsmon,ich)
        sumtotal(i)=sumtotal(i)+cts(i)
      enddo
      sumtotalmon=sumtotalmon+ctsmon
      return
      end subroutine processline                                             
      logical function led(ic,xl,xh)
      use specfiles
      integer(kind=4),intent(in) :: ic
      real(kind=8),intent(in) :: xl, xh
      integer(kind=4)i
      if(.not.allocated(iexarray).or. .not. allocated(exarray)) then
! oh dear, function should never have been called
        write(*,'(a)')'Bug in program, bailing out, sorry!'
        write(*,'(a)')'Please mail a bug report to wright@esrf.fr'
        stop
      endif
      led=.false.
      do i=1,iexrc
        if((iexarray(i)+1).eq.ic)then !  found the right line
          if(int(exarray(i,3)).lt.1)then
           if(xl.gt.exarray(i,1) .and.xl.lt.exarray(i,2))led=.true.
           if(xh.gt.exarray(i,1) .and.xh.lt.exarray(i,2))led=.true.
          else
           if(xl.gt.exarray(i,1) .and.xl.lt.exarray(i,2) .and. &
     & iscan.eq.int(exarray(i,3)))led=.true.
           if(xh.gt.exarray(i,1) .and.xh.lt.exarray(i,2) .and. &
     & iscan.eq.int(exarray(i,3)))led=.true.
          endif
        endif
      enddo
      return
      end function led                                                  
      logical function window(ar,n)
      use specfiles ! 
      integer(kind=4),intent(in) :: n
      real(kind=8),dimension(n),intent(in) :: ar
      if(winlog)then 
        wincol=whichcolumn(wincnt(1:len_trim(wincnt)))
        if(wincol.lt.0)then
           write(*,*)'PROBLEMS WITH YOUR WINDOW - cannot find column', &
     & wincnt(1:len_trim(wincnt))
           write(*,*)'No further windowing will be attempted'
           winlog=.false. 
           window=.true.
        else
          if( (ar(wincol).gt.winlow) .and. (ar(wincol).lt.winhigh)) then
            window=.true.
          else
            window=.false.
          endif
        endif
      else
        window=.true.
      endif
      end function window                                                  

      real(kind=8) function convert_unit_function( x )
      real(kind=8), intent(in) :: x
      real(kind=8) :: qt
      select case((units))
        case('T') ! Two theta, do nothing
           convert_unit_function = x
           return
        case('Q') ! Q - constants from initialise_rebin
           convert_unit_function = four_pi*sin(x*pi_over_360)/wavelength
           return
        case('R') ! Q^2 
           qt = four_pi * sin( x * pi_over_360 )/wavelength
! use a signed quantity
!        write(*,*)qt,qt*qt,SIGN( qt * qt , qt ),x
           convert_unit_function = sign( qt * qt , qt )
           return
        case default
           convert_unit_function = x
           return
      end select
      return
      end function convert_unit_function
!
!
! To be called once per scan - picks up wavelength from command line
      subroutine convert_unit_setupQ( Q )
! from rebin module :::: logical :: wavelength_set = .false.
      real(kind=8), dimension(60), intent(in) :: Q
      if ( (units .eq. 'T') .or. wavelength_set) then
         return
      endif 
      if( Q(4) .gt. 1.0E-15 ) then
          wavelength = Q(4)
! Assume always the same for all scans:
          write(*,*)
          write(*,*)'Got wavelength',wavelength,'from spec file Q[3]'
          wavelength_set = .true.
      endif
! Check for errors
! let them have Angstroms as meters (1e-11) 
      if ((wavelength .lt. 1.0E-15) .and. (units .ne. 'T') ) then
        write(*,*)'Your wavelength is a bit small ',wavelength 
        write(*,*)'Giving up, try putting wvln=1.234 on command line'
        stop
      endif
!      write(*,*)'calling setupQ from convert_unit_setupQ'
      call unit_lims()
      return
      end subroutine convert_unit_setupQ

      subroutine unit_lims()
      real(kind=8) :: stmp
! apply the limits according to if the user changed them
      tthlow = convert_unit_function( user_tthlow )
      tthhigh = convert_unit_function( user_tthhigh )
! set the step size to be right at 30 degrees twotheta
      step = convert_unit_function( 30.0d0+user_step ) 
      step = step - convert_unit_function( 30.0d0 )
! round it to be some printable represention 
      write(*,*)
! round to something nicely printable
      if(units.ne.'T')then
        stmp = int(log10(step)-2.0d0)*1.0d0 
        stmp = exp( log(10.0d0) * stmp )
        step = int(step/stmp)*stmp
      endif
!      write(*,*)'from',tthlow,'to',tthhigh,'step',step
      call checkrebinpars
      call rebinallocate
      return
      end subroutine unit_lims



      real(kind=4) function rlcg( )
! This returns a float in the range 0 to 1 which is a 
! deterministic sequence. Will be the same over many runs.
      integer(kind=4) :: a = 69069, c = 1327217885 
      integer(kind=4), save :: x = 1
      real(kind=4), save :: m = 2147483647
      x = x*a + c
      rlcg = (real(x)/m+1.0)/2.0
!       write(*,'(I12,1X,F12.10)')x, rlcg
      end function rlcg

      end module rebin                                                        
      module summation
      integer(kind=4):: isc=0, i6s=0 !  which scale factor, 6sig pts
      integer(kind=4)::nzap=0, nsuperzap=0 ! for zapping
      real(kind=8),allocatable :: sumdata(:,:),hist(:,:)           
      real(kind=8) :: alp=0.5d0  ! Bayesian zero counts fudge
      real(kind=8) :: scalinp=1.0d5  ! Scale factor for .inp files
      real(kind=8) :: zap=6.0d0 ! level for outlier elimination (median filter?)
      real(kind=8) :: superzaplevel=1.0d0 ! level for zinger elimination
      logical :: renorm=.false.  ! to update temp.res file efficiencies
      logical :: zapping=.false.
      logical :: superzap=.false.
      logical ::  medianofchannels=.false.
      contains
      subroutine calibsum
      use specfiles
      use rebin
      integer(kind=4) :: m,n,ierr
      real(kind=8),dimension(NCHAN) :: tmult,tmulterr
      call effic(n,m,tmult,tmulterr)
      call reporteffic(n,m,tmult,tmulterr)
      ierr=0
      if(.not.allocated(sumdata))allocate(sumdata(3,npts),stat=ierr)
! 3 ! cts, mon, e(mon)
      if(ierr.ne.0)then
       write(*,'(a)') 'Memory allocation error in calibsum'
       stop
      endif
      if(.not.allocated(hist))allocate(hist(2,npts),stat=ierr) 
! final signal and esd
      if(ierr.ne.0)then
       write(*,'(a)') 'Memory allocation error in calibsum'
       stop
      endif
      call sumthem      ! combine detectors in sumdata array
      if(medianofchannels)call medianchannels
      call normerr      ! determine error bars and fill in hist
      call checkdets    ! check ascan matches hist
1     if(zapping)then
       write(*,*)
       write(*,'(A,F8.5,A)')"After zapping a sigma level ",zap," ... "
       call sumthem      ! combine detectors in sumdata array
       call normerr      ! determine error bars and fill in hist
       call checkdets    ! check ascan matches hist
      endif
      if(nzap.gt.0)goto 1
2     if(superzap .and. nsuperzap.gt.0)then
       write(*,*)
       write(*,'(A,F8.5,A)')"After superzapping at level ", superzaplevel,"..."
       call sumthem      ! combine detectors in sumdata array
       call normerr      ! determine error bars and fill in hist
       call checkdets    ! check ascan matches hist
      endif
      if(nsuperzap.gt.0)goto 2
      return 
      end subroutine calibsum                                                
      subroutine effic(n,m,tmult,tmulterr)
      use rebin ! mult, multerr, nchan, ascan, tempres, npts
      integer(kind=4),intent(inout)::n,m
      real(kind=8),intent(inout),dimension(NCHAN)::tmult,tmulterr
      integer(kind=4) :: i,j,k,nex
      real(kind=8), dimension(4,NCHAN) :: signal
      real(kind=8) :: sumsig,sumsige
      signal=0.0d0; n=0
      write(*,'(a)')'Determining detector efficiencies'
      nex=0; do j=1,nchannel; if(logexdet(j).eq.1)nex=nex+1; enddo
      do i=1,npts
        k=0
        do j=1,nchannel         ! must have more than 1 monitor in bin to use
          if(ascan(2*j,i).gt.1.0d0)then
             if(ascan(2*j-1,i).gt.minrenormsig)k=k+1
          endif 
        enddo
        if(k.eq.(nchannel-nex))then
          n=n+1
          do j=1,nchannel    ! signal is a big bin for all overlapping points
            if(logexdet(j).eq.1)cycle
! If a channel is excluded for the whole tth range in a file then this will fail?
            signal(1,j)=signal(1,j)+ascan(2*j-1,i) ! counts
            signal(2,j)=signal(2,j)+ascan(2*j,i) ! monitors
          enddo
        endif
      enddo
      if(n.gt.0)then
       do j=1,nchannel          ! Normalise Signal and get the error bar on it
         if(logexdet(j).eq.1)cycle
         signal(3,j)=signal(1,j)/signal(2,j) 
         signal(4,j)=signal(3,j) *                                      &
     &      sqrt(1.0d0/signal(1,j)+1.0d0/signal(2,j))
       enddo     
       sumsig=sum(signal(3,:)) ! sum of all signals in overlap region
       sumsige=0.0d0
       do i=1,nchannel          ! Error in sum of all signals
        sumsige=sumsige+signal(4,i)*signal(4,i)
       enddo
       sumsige=sqrt(sumsige)
       do j=1,nchannel          ! Finally get the channel efficiencies
        if(logexdet(j).eq.1)then
          tmult(j)=1.0d0; tmulterr(j)=0.0d0
        else
          tmult(j)=(real((nchannel-nex),8)*signal(3,j)/sumsig)
          tmulterr(j)=tmult(j)*sqrt((signal(4,j)/signal(3,j))**2 +      &
     &     (sumsige/sumsig)**2)
        endif
       enddo       
      if(.not.tempres)then    ! Copy these to rebin module if no temp.res file
       mult=tmult
       multerr=tmulterr
       m=1
      else ! tempres     
       m=2                ! Flag the need to print the temporary (unused vals)
      endif ! tempres
      else ! n.gt.0
        m=1               ! no overlap so only one to print
        if(.not.tempres)then
          mult=1.0d0   ! Make sure defaults are always 1.0d0
          multerr=0.0d0
        endif
      endif   ! if n.gt.0
      return
      end subroutine effic                                                   
      subroutine sumthem
      use rebin ! for ascan and mult
      integer(kind=4) :: i,j     
      sumdata=0.0d0 ! the big array where the sum of channels goes
      do i=1,npts 
       do j=1,nchannel
! counts
        sumdata(1,i)=sumdata(1,i)+ascan(2*j-1,i)
! normalised mon      = m1*e1 + m2*e2 + ....
        sumdata(2,i)=sumdata(2,i)+ascan(2*j,i)*mult(j)
! emon**2  = 
        sumdata(3,i)=sumdata(3,i)+                                       &
!             c         (     e**2     +         c   *  <e>**2     )   
     &     ascan(2*j,i) * ( mult(j)**2 + ascan(2*j,i)*multerr(j)**2)
       enddo
! emon=sqrt(emon**2)
       sumdata(3,i)=dsqrt(sumdata(3,i))
      enddo
      return
      end subroutine sumthem                                                 
      subroutine normerr
! alp, hist & sumdata available as member of summation
      use rebin ! npts
      real(kind=8):: msq, s
      integer(kind=4) :: i, n, n3
      n=0; s=0.0d0; n3=0
      do i=1,npts
        if(sumdata(2,i).gt.minmon)then ! needs at least 1. mon count 
          hist(1,i)=sumdata(1,i)/sumdata(2,i) ! correct
          msq=(sumdata(2,i)*sumdata(2,i)) ! always gt zero
          hist(2,i)=sqrt(                                               &
     &   (sumdata(1,i)+alp)/msq + (sumdata(1,i)*sumdata(3,i)/msq)**2)  
          n=n+1
          s=s+hist(1,i)**2/hist(2,i)**2
          if(hist(1,i)/hist(2,i).lt.3.0d0)n3=n3+1
      else
          hist(1,i)= 0.0d0  ! unobserved regions are filled with nonsense
          hist(2,i)=-1.0d0
      endif
      enddo
      write(*,1000)100.0d0*sqrt(real(n,8)/s),n3,n
1000  format('R_exp = ',f7.3,' with ',i7,                               &
     & ' pts having I/<I> less than 3, from ',i7,' pts obs')
      end subroutine normerr                                                 
      integer(kind=4) function ctchan(i)
      use rebin
      integer(kind=4),intent(in) :: i
      integer(kind=4) :: j
      ctchan=0
      do j=1,nchannel
        if(ascan(2*j,i).gt.minmon)ctchan=ctchan+1
      enddo; return; end function ctchan                                  
      subroutine checkdets
      use rebin ! for ascan and mult arrays
      integer(kind=4) :: i, j, ipt3s, ipt6s, iapts, myctchan
      real(kind=8) :: c,y,ey2,wd2 ! chi2,sig,error^2 and wtd diff^2
      real(kind=8),dimension(nchan)::sz ! superzap array
      ipt3s=0 ; ipt6s=0; c=0.0d0; iapts=0; nzap=0; nsuperzap=0;
      do i = 1, npts
        myctchan=ctchan(i)
        if(myctchan.ge.2)then
        if(superzap)sz=-1.0
        do j = 1, nchannel
          if(ascan(2*j,i).gt.1.0d0)then ! need one monitor count to bother
            y=ascan(2*j-1,i)/ascan(2*j,i)/mult(j)   ! OK
! assumes zero error in the efficiency    
            ey2= ((ascan(2*j-1,i)+alp)/ascan(2*j,i)**2                  &
     &           +  ascan(2*j-1,i)**2/ascan(2*j,i)**3)/mult(j)**2
            wd2=(y-hist(1,i))**2/(ey2+hist(2,i)**2)                          
            if(superzap)sz(j)=y ! 
            c=c+wd2
            iapts=iapts+1
            if(wd2.gt.9.0d0) ipt3s=ipt3s+1  ! 3s^2=9.0
            if(wd2.gt.36.0d0)ipt6s=ipt6s+1  ! 6s^2=36.0
            if(wd2.gt.zap*zap .and. zapping)then
               ascan(2*j,i)=0.0 ! set mon and det to zero for zapped points
               ascan(2*j-1,i)=0.0
               nzap=nzap+1
            endif
           endif
         enddo
         ! now check with superzap
         if(superzap .and. myctchan.ge.3)call superzapem(sz,ascan(:,i),nchannel)
         endif
      enddo
      if(iapts.gt.0) then
      c=c/real(iapts,8)
      write(*,1000)c,iapts
1000  format('Reduced chi**2 for channel merge =',F9.4,' from ',        &
     & i10,' pairs of pts')
      write(*,1001)100.0d0*real(ipt3s)/real(iapts),                     &
     &100.0d0*real(ipt6s)/real(iapts)
1001  format(f5.2,'% differ by >3 sigma, ',f7.4,                        &
     & '% by >6 sigma (ideally 0.04% and 0.0000%)' ) 
      i6s=ipt6s ! module variable
      if(zapping)then
       write(*,'(A,I8,A)')'Zapping ',nzap,' points'
      endif
      if(superzap)write(*,'(A,I8,A)')'Superzapping ',nsuperzap,' points'
      else
       write(*,*)'No channel overlap, whatsoever, was found'
       i6s=0
      endif
      return
      end subroutine checkdets                                             
      subroutine reporteffic(n,m,tmult,tmulterr)
      use rebin ! mult, multerr, NCHAN
      integer(kind=4),intent(in) :: n,m
      real(kind=8),intent(in),dimension(nchan)::tmult,tmulterr
      integer(kind=4) :: i
      logical :: trex
      if(n.gt.0)then
       write(*,'(a,i9,a)')'Channel efficiences found from ',n,          &
     &  ' points where all detectors overlap'
       if(m.eq.1)then
        write(*,'(a)')'Det     Offset      Effic    <Effic>'
        do i=1,NCHANNEL
! i-1 instead of i to start at channel zero
         write(*,'(i3,3(1X,F10.7))')i-1,offset(i),mult(i),multerr(i)
        enddo
       endif ! m.eq.1
       if(m.eq.2)then
        write(*,'(a)')'Efficiencies from temp.res file,'//              &
     &     ' the values found now are compared'
        if(.not.renorm)then ; write(*,'(a)')                            &
     &     'Det     Offset      Effic    <Effic> current unused values'
        else ; write(*,'(a)')                                           &
     &     'Det     Offset    Old Eff    Old <E>    New Eff    New <E>'
        endif
        do i=1,NCHANNEL
! i-1 instead of i to start at channel zero
         write(*,'(i3,5(1X,F10.7))')i-1,offset(i),mult(i),multerr(i),     &
     &   tmult(i),tmulterr(i)
        enddo 
        if(renorm)then ; mult=tmult ; multerr=tmulterr ; endif
       endif ! m.eq.2 
       inquire(file='temp.res',exist=trex)
       if(.not.trex .or. renorm)call tempreswrite
       if(renorm) renorm=.false. ! for sumall - can only renorm once
      else ! n.gt.0
       write(*,'(a,$)')'No detector overlap found, efficiencies'
       if(tempres)then
        write(*,'(a)')' from temp.res file'
       else
        write(*,'(a)')' are probably wrong'
        call tempreswrite
       endif
       write(*,'(a)')'Det     Offset      Effic    <Effic>'
       do i=1,NCHANNEL
! i-1 instead of i to start at channel zero
        write(*,'(i3,3(1X,F10.7))')i-1,offset(i),mult(i),multerr(i)
       enddo
      endif
      return
      end subroutine reporteffic                                             
    subroutine tempreswrite
      use rebin ! offset, nchan, mult
      integer(kind=4) :: ierr, i
      open(unit=16,file='temp.res',status='UNKNOWN',access='SEQUENTIAL',&
     & form='FORMATTED',iostat=ierr)
      if(ierr.ne.0)then
         write(*,'(a)')'Couldn''t open temp.res to write!!!'
         return !! bugs out
      endif
      write(*,'(a)')'Created temp.res file'
      do i=1,nchannel
        write(16,1000)offset(i),mult(i),multerr(i),0.0d0
      enddo
1000  format(4(f11.8,','))
      close(16)
      tempres=.true. ! Should exist and be readable now
      return; end subroutine tempreswrite                                    
      subroutine scaltot
! uses hist and sumdata arrays in summation module
      real(kind=8) :: factor, sum1, sum2
      sum1=sum(sumdata(1,:)) ! sum up real counts
      sum2=sum(hist(1,:))    ! sum up normed counts
      factor=sum1/sum2       ! determine multiplier
      hist=hist*factor       ! rescale hist array
      return
      end subroutine scaltot                                               
      subroutine scalpk
! uses hist and sumdata arrays in summation module
      real(kind=8) :: factor
      integer(kind=4), dimension(3) :: i
      i=maxloc(sumdata,dim=2)  ! get index of highest number of counts
      factor=sumdata(1,i(1))/hist(1,i(1)) ! get scaling factor
      hist=hist*factor
      return
      end subroutine scalpk                                                
      subroutine superzapem(sz,as,n)
      integer(kind=4), intent(in) :: n
      real(kind=8),dimension(n),intent(in) :: sz   ! signal
      real(kind=8),dimension(2*n),intent(inout) :: as ! ascan array
      integer(kind=4) :: i, ib, nc
      real(kind=8) :: s, ss, mean1, mean2, esd1, esd2, pc
      ! make mean and esd from sz
      nc=0;s=0.0;ss=0.0
      ib=1; pc=0.0
      do i = 1,n
        if(sz(i).ge.0.0)then 
           s=s+sz(i)
           ss=ss+sz(i)*sz(i)
           nc=nc+1
           if(sz(i).gt.sz(ib))ib=i
        endif
      enddo
      mean1=s/nc
      esd1=(ss-mean1*mean1)/nc
      ! now without ib
      s=0.0;ss=0.0;nc=0
      do i = 1,n
        if(sz(i).ge.0.0 .and. i.ne.ib)then 
           s=s+sz(i)
           ss=ss+sz(i)*sz(i)
           nc=nc+1
        endif
      enddo
      mean2=s/nc
      esd2=(ss-mean2*mean2)/nc
      if(esd1.gt.0.0) pc=(esd1-esd2)/esd1 ! should use mean info too.
      if(pc.gt.superzaplevel)then ! zap the ib channel
        as(2*ib)=0.0
        as(2*ib-1)=0.0
        nsuperzap=nsuperzap+1
      endif
      return
      end subroutine superzapem                                            

      ! netlib code

! *DECK DSORT
      SUBROUTINE DISORT (DX, DY, N, KFLAG)
!C***BEGIN PROLOGUE  DSORT
!C***PURPOSE  Sort an array and optionally make the same interchanges in
!C            an auxiliary array.  The array may be sorted in increasing
!C            or decreasing order.  A slightly modified QUICKSORT
!C            algorithm is used.
!C***LIBRARY   SLATEC
!C***CATEGORY  N6A2B
!C***TYPE      DOUBLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
!C***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
!C***AUTHOR  Jones, R. E., (SNLA)
!C           Wisniewski, J. A., (SNLA)
!C***DESCRIPTION
!C
!C   DSORT sorts array DX and optionally makes the same interchanges in
!C   array DY.  The array DX may be sorted in increasing order or
!C   decreasing order.  A slightly modified quicksort algorithm is used.
!C
!C   Description of Parameters
!C      DX - array of values to be sorted   (usually abscissas)
!C      DY - array to be (optionally) carried along
!C      N  - number of values in array DX to be sorted
!C      KFLAG - control parameter
!C            =  2  means sort DX in increasing order and carry DY along.
!C            =  1  means sort DX in increasing order (ignoring DY)
!C            = -1  means sort DX in decreasing order (ignoring DY)
!C            = -2  means sort DX in decreasing order and carry DY along.
!C
!C***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
!C                 for sorting with minimal storage, Communications of
!C                 the ACM, 12, 3 (1969), pp. 185-187.
!C***ROUTINES CALLED  XERMSG
!C***REVISION HISTORY  (YYMMDD)
!C   761101  DATE WRITTEN
!C   761118  Modified to use the Singleton quicksort algorithm.  (JAW)
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890831  Modified array declarations.  (WRB)
!C   891009  Removed unreferenced statement labels.  (WRB)
!C   891024  Changed category.  (WRB)
!C   891024  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!C   901012  Declared all variables; changed X,Y to DX,DY; changed
!C           code to parallel SSORT. (M. McClain)
!C   920501  Reformatted the REFERENCES section.  (DWL, WRB)
!C   920519  Clarified error messages.  (DWL)
!C   920801  Declarations section rebuilt and code restructured to use
!C           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
!C***END PROLOGUE  DSORT
!C     .. Scalar Arguments ..
      INTEGER(kind=4) KFLAG
      INTEGER(kind=4) N
!C     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
      INTEGER(kind=4) DY(*)
!C     .. Local Scalars ..
      DOUBLE PRECISION R, T, TT, TTY, TY
      INTEGER I, IJ, J, K, KK, L, M, NN
!C     .. Local Arrays ..
      INTEGER IL(21), IU(21)
!C     .. External Subroutines ..
!C      EXTERNAL XERMSG
!C     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
!C***FIRST EXECUTABLE STATEMENT  DSORT
      NN = N
      IF (NN .LT. 1) THEN
!         CALL XERMSG ('SLATEC', 'DSORT',
!     +      'The number of values to be sorted is not positive.', 1, 1)
         RETURN
      ENDIF
!C
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
!         CALL XERMSG ('SLATEC', 'DSORT',
!     +      'The sort control parameter, K, is not 2, 1, -1, or -2.', 2,
!     +      1)
         RETURN
      ENDIF
!C
!C     Alter array DX to get decreasing order if needed
!C
      IF (KFLAG .LE. -1) THEN
         DO 10 I=1,NN
            DX(I) = -DX(I)
   10    CONTINUE
      ENDIF
!C
      IF (KK .EQ. 2) GO TO 100
!C
!C     Sort DX only
!C
      M = 1
      I = 1
      J = NN
      R = 0.375D0
!C
   20 IF (I .EQ. J) GO TO 60
      IF (R .LE. 0.5898437D0) THEN
         R = R+3.90625D-2
      ELSE
         R = R-0.21875D0
      ENDIF
!C
   30 K = I
!C
!C     Select a central element of the array and save it in location T
!C
      IJ = I + INT((J-I)*R)
      T = DX(IJ)
!C
!C     If first element of array is greater than T, interchange with T
!C
      IF (DX(I) .GT. T) THEN
         DX(IJ) = DX(I)
         DX(I) = T
         T = DX(IJ)
      ENDIF
      L = J
!C
!C     If last element of array is less than than T, interchange with T
!C
      IF (DX(J) .LT. T) THEN
         DX(IJ) = DX(J)
         DX(J) = T
         T = DX(IJ)
!C
!C        If first element of array is greater than T, interchange with T
!C
         IF (DX(I) .GT. T) THEN
            DX(IJ) = DX(I)
            DX(I) = T
            T = DX(IJ)
         ENDIF
      ENDIF
!C
!C     Find an element in the second half of the array which is smaller
!C     than T
!C
   40 L = L-1
      IF (DX(L) .GT. T) GO TO 40
!C
!C     Find an element in the first half of the array which is greater
!C     than T
!C
   50 K = K+1
      IF (DX(K) .LT. T) GO TO 50
!C
!C     Interchange these elements
!C
      IF (K .LE. L) THEN
         TT = DX(L)
         DX(L) = DX(K)
         DX(K) = TT
         GO TO 40
      ENDIF
!C
!C     Save upper and lower subscripts of the array yet to be sorted
!C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 70
!C
!C     Begin again on another portion of the unsorted array
!C
   60 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
!C
   70 IF (J-I .GE. 1) GO TO 30
      IF (I .EQ. 1) GO TO 20
      I = I-1
!C
   80 I = I+1
      IF (I .EQ. J) GO TO 60
      T = DX(I+1)
      IF (DX(I) .LE. T) GO TO 80
      K = I
!C
   90 DX(K+1) = DX(K)
      K = K-1
      IF (T .LT. DX(K)) GO TO 90
      DX(K+1) = T
      GO TO 80
!C
!C     Sort DX and carry DY along
!C
  100 M = 1
      I = 1
      J = NN
      R = 0.375D0
!C
  110 IF (I .EQ. J) GO TO 150
      IF (R .LE. 0.5898437D0) THEN
         R = R+3.90625D-2
      ELSE
         R = R-0.21875D0
      ENDIF
!C
  120 K = I
!C
!C     Select a central element of the array and save it in location T
!C
      IJ = I + INT((J-I)*R)
      T = DX(IJ)
      TY = DY(IJ)
!C
!C     If first element of array is greater than T, interchange with T
!C
      IF (DX(I) .GT. T) THEN
         DX(IJ) = DX(I)
         DX(I) = T
         T = DX(IJ)
         DY(IJ) = DY(I)
         DY(I) = TY
         TY = DY(IJ)
      ENDIF
      L = J
!C
!C     If last element of array is less than T, interchange with T
!C
      IF (DX(J) .LT. T) THEN
         DX(IJ) = DX(J)
         DX(J) = T
         T = DX(IJ)
         DY(IJ) = DY(J)
         DY(J) = TY
         TY = DY(IJ)
!C
!C        If first element of array is greater than T, interchange with T
!C
         IF (DX(I) .GT. T) THEN
            DX(IJ) = DX(I)
            DX(I) = T
            T = DX(IJ)
            DY(IJ) = DY(I)
            DY(I) = TY
            TY = DY(IJ)
         ENDIF
      ENDIF
!C
!C     Find an element in the second half of the array which is smaller
!C     than T
!C
  130 L = L-1
      IF (DX(L) .GT. T) GO TO 130
!C
!C     Find an element in the first half of the array which is greater
!C     than T
!C
  140 K = K+1
      IF (DX(K) .LT. T) GO TO 140
!C
!C     Interchange these elements
!C
      IF (K .LE. L) THEN
         TT = DX(L)
         DX(L) = DX(K)
         DX(K) = TT
         TTY = DY(L)
         DY(L) = DY(K)
         DY(K) = TTY
         GO TO 130
      ENDIF
!C
!C     Save upper and lower subscripts of the array yet to be sorted
!C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 160
!C
!C     Begin again on another portion of the unsorted array
!C
  150 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
!C
  160 IF (J-I .GE. 1) GO TO 120
      IF (I .EQ. 1) GO TO 110
      I = I-1
!C
  170 I = I+1
      IF (I .EQ. J) GO TO 150
      T = DX(I+1)
      TY = DY(I+1)
      IF (DX(I) .LE. T) GO TO 170
      K = I
!C
  180 DX(K+1) = DX(K)
      DY(K+1) = DY(K)
      K = K-1
      IF (T .LT. DX(K)) GO TO 180
      DX(K+1) = T
      DY(K+1) = TY
      GO TO 170
!C
!C     Clean up
!C
  190 IF (KFLAG .LE. -1) THEN
         DO 200 I=1,NN
            DX(I) = -DX(I)
  200    CONTINUE
      ENDIF
      RETURN
      END      SUBROUTINE DISORT 

      ! netlib code 
      subroutine medianchannels
      use rebin
      integer(kind=4) :: i,j,k
      real(kind=8), dimension(nchan) :: signal ! local copy
      integer(kind=4), dimension(nchan) :: iactive ! local copy
      sumdata=0.0d0 ! the big array where the sum of channels goes
      do i=1,npts 
        k=0
        do j=1,nchannel
          if(ascan(2*j,i).gt.0)then
            k=k+1
            signal(k)=ascan(2*j-1,i)/(ascan(2*j,i)*mult(j))
            iactive(k)=j ! which channel was active
          endif  ! signal is computed
        enddo
        if(k.eq.0)then ! nothing recorded at all
           sumdata(1,i)=0.
           sumdata(2,i)=0.
           sumdata(3,i)=0.
        endif
        if(k.eq.1)then ! only one channel active
           j=iactive(1)
           sumdata(1,i)=ascan(2*j-1,i)
           sumdata(2,i)=ascan(2*j,i)*mult(j)
           sumdata(3,i)= &
     &     ascan(2*j,i) * ( mult(j)**2 + ascan(2*j,i)*multerr(j)**2)
           sumdata(3,i)=dsqrt(sumdata(3,i))
        endif
        if(k.gt.1)then ! more channels active
           CALL DISORT(SIGNAL,IACTIVE,K,2) ! 2 means increasing order
           j=iactive((k+1)/2) ! median channel
           sumdata(1,i)=ascan(2*j-1,i)*k
           sumdata(2,i)=ascan(2*j,i)*mult(j)*k
           sumdata(3,i)= k* &
     &     ascan(2*j,i) * ( mult(j)**2 + ascan(2*j,i)*multerr(j)**2)
           sumdata(3,i)=dsqrt(sumdata(3,i))
        endif
      enddo
      return
      end subroutine medianchannels                                          
      end module summation                                                   
      module outputfiles
!      use summation  ! to get the final dataset
       character(len=1024) :: outfile
      integer(kind=4) :: ioutunit=11, ilow, ihigh, itot, logfile=21
      logical :: gsas=.false. , spf=.false. , pds=.false., diag=.true.
      logical :: epf=.false. , xye=.true.
      character(len=1024) :: title                                           
      logical :: bcm=.false.
      contains
      subroutine outputformats(is)
      integer(kind=4),intent(in)::is
      if(epf) call epfout(is)
      if(xye) call xyeout(is)
      if(gsas) call gsasout(is)
      if(spf)call spfout(is)
      if(pds)call pdsout(is)
      return; end subroutine outputformats                                   
      subroutine gsasout(is)
      use rebin ! ibin, bincen
      use specfiles
      use summation
      integer(kind=4),intent(in)::is
      integer(kind=4) :: low, ipts, nrec, i, k, n
      real(kind=8), dimension(2,5) :: vals 
      real(kind=8) :: top ! for scaling to format
      character(len=10) :: str
      if(units.ne.'T')then
        write(*,*)'Can only write GSAS format for two theta data'
        stop
      endif
      call filext('.gsa',is)
      open(file=outfile,unit=ioutunit,status='UNKNOWN')
      write(title,1002)       specsfilename(1:len_trim(specsfilename)), &
     & specsdate(1:len_trim(specsdate))
1002  format(a,1X,a)      
      if(is.gt.0)then
       write(str,'(i10)')is; str=adjustl(str)
       title=title(1:len_trim(title))//' Scan '//str(1:len_trim(str))
      endif                                                                      
      write(ioutunit,'(a80)')title
      write(ioutunit,1015)
1015  format("Instrument parameter      id31.prm                  ",    &
     & "                            ")
      ilow=getfirstpoint(hist,2,npts,2)
      ihigh=getlastpoint(hist,2,npts,2)
      low=ibin(0.0d0)
      if(low.gt.ilow) then ! the data goes below zero !
       low=low+1           ! one point after zero
      else
       low=ilow         ! value from module
      endif
      ipts=ihigh-low    ! how many points for GSAS
      nrec=(ipts-mod(ipts,5))/5  ! how many records
      if(mod(ipts,5).ne.0)nrec=nrec+1
      write(ioutunit,1000)ipts,nrec,bincen(low)*100.0d0,step*100.0d0
!             1234567   +16     123456    20         123456789012 = 61(+19=80)
1000  format('BANK 1 ',2(I7,1X),'CONST ',2(F9.3,1X),' 0.0 0.0 ESD'      &
!        1234567890123456789
     & ,'                   ')
! scale factor to fit into F8.1 format? Max must be less than 999999
      do i = low,ihigh
        if(hist(1,i).gt.top)top=hist(1,i)
      enddo
      if (top.gt.999999.9d0)then
        top=999999.9d0/top
        write(*,'(a)')'Rescaled your data for GSAS to avoid overflows'
        write(*,*)'   The data were multiplied by',top
      else
        top=1.0d0
      endif
      n = ipts-mod(ipts,5)
      do i = low,n+low-1,5
       do k=1,5
         vals(1,k)=hist(1,k+i-1)*top
         vals(2,k)=hist(2,k+i-1)*top
         if (vals(2,k).lt.0.)then
           vals(1,k)=0.0
           vals(2,k)=999999.9
         endif
       enddo 
       write(ioutunit,'(10F8.1)')(vals(1,k),vals(2,k),k=1,5)
      enddo
      if(mod(ipts,5).ne.0) then
        do k=low+n,ihigh
          if (hist(2,k).gt.0.)then
           vals(1,k+1-low-n)=hist(1,k)*top
           vals(2,k+1-low-n)=hist(2,k)*top
          else
           vals(1,k+1-low-n)=0.
           vals(1,k+1-low-n)=999999.9
          endif
        enddo
        write(ioutunit,'(10F8.1)')(vals(1,k),vals(2,k),k=1,ihigh-low-n)
      endif
      close(ioutunit)
      write(*,'(a)')'Wrote '//outfile(1:len_trim(outfile))
      return
      end subroutine gsasout                                                 
      subroutine spfout(is)
      use rebin
      use summation
      use specfiles
      integer(kind=4),intent(in)::is
      integer(kind=4)::low,ipts
      character(len=18) :: time
      character(len=10) :: str
      integer(kind=4) :: n, i, k
      character(len=3) :: months(12)
      data months /'JAN','FEB','MAR','APR','MAY','JUN',                 &
     &             'JUL','AUG','SEP','OCT','NOV','DEC'/
      if(units.ne.'T')then
        write(*,*)'Can only write SPF format for two theta data'
        stop
      endif
      call filext('.spf',is)
      open(file=outfile,unit=ioutunit,status='UNKNOWN')
      ilow=getfirstpoint(hist,2,npts,2)
      ihigh=getlastpoint(hist,2,npts,2)
      low=ibin(0.0d0)
      if(low.gt.ilow) then ! the data goes below zero !
       low=low+1           ! one point after zero
      else
       low=ilow         ! value from module
      endif
      ipts=ihigh-low    ! how many points
      n = ipts-mod(ipts,10) 
      write(title,1002)specsdate(1:len_trim(specsdate)),                &
     &  specsfilename(1:len_trim(specsfilename))
1002  format(a,1X,a)     
      if(is.gt.0)then
       write(str,'(i10)')is; str=adjustl(str)
       title=title(1:len_trim(title))//' Scan '//str(1:len_trim(str))
      endif                                                                  
      write(ioutunit,'(a80)')title
      write(ioutunit,2)ipts,bincen(low),bincen(ihigh),step
 2    format("       1       1       1"/1i8,"       1",                 &
     & 3f10.3,"  0.000000  0.000000")
!                               CCYYMMDD      HHMMSS.SSS
       call date_and_time(date=time(1:8),time=time(9:18))
! 123456789012345678
! hh:mm:ss dd-mm-yy 
       read(time(5:6),'(i2)')i
       time=time(9:10)//':'//time(11:12)//':'//time(13:14)//' '//       & 
     & time(7:8)//'-'//months(i)//'-'//time(3:4)
       write(ioutunit,3)time(1:18)
 3     format("  100000.0       0.0     ",1a18,                         &
     & "      0      0 SYNCHROTRON  ID31 ESRF"/8("     0.000"))
       do i = low,n+low-1,10
       write(ioutunit,4)bincen(i),(nint(hist(1,k)),k=i,i+9)
       write(ioutunit,5)(hist(2,k),k=i,i+9)
 4     format(1f8.3,10i7)
 5     format(8x,10f7.2)
       enddo
       if(mod(ipts,10).ne.0)then 
       write(ioutunit,4)bincen(n+low),(nint(hist(1,k)),k=n+low,ipts)
       write(ioutunit,5)(hist(2,k),k=n+low,ipts)
       endif
       write(ioutunit,6)
 6     format('       0 HKL reflection(s)'/                             &
      & '       0 excluded region(s)')                                  
       close(ioutunit)
       write(*,'(a)')'Wrote '//outfile(1:len_trim(outfile))
       end subroutine spfout                                                 
       subroutine pdsout(is)
      use rebin
      use summation
      use specfiles
      integer(kind=4),intent(in)::is
      integer (kind=4) :: low, ipts, i, k, n
      character(len=10):: str
      if(units.ne.'T')then
        write(*,*)'Can only write PDS format for two theta data'
        stop
      endif
      call filext('.pds',is)
      open(file=outfile,unit=ioutunit,status='UNKNOWN')    
      write(title,1002)specsdate(1:len_trim(specsdate)),                &
     &  specsfilename(1:len_trim(specsfilename))
1002  format(a,1X,a)      
      if(is.gt.0)then
       write(str,'(i10)')is; str=adjustl(str)
       title=title(1:len_trim(title))//' Scan '//str(1:len_trim(str))
      endif                                                                  
      write(ioutunit,'(1a80)')title      
      write(ioutunit,7)
      ilow=getfirstpoint(hist,2,npts,2)
      ihigh=getlastpoint(hist,2,npts,2)
      low=ibin(0.0d0)
      if(low.gt.ilow) then ! the data goes below zero !
       low=low+1           ! one point after zero
      else
       low=ilow         ! value from module
      endif
      ipts=ihigh-low
7     format("   Start     Step      End     Monitor")
      write(ioutunit,8)bincen(low),step,bincen(ihigh)
8     format(3f9.3,'    10000')
      n = ipts-mod(ipts,10) 
      do 800 i = low,n+low-1,10
        write(ioutunit,9)(hist(2,k),k=i,i+9)
        write(ioutunit,10)(nint(hist(1,k)),k=i,i+9)
9     format(10f8.4)
10    format(10i8)
800   continue
      if(mod(ipts,10).ne.0)then
       write(ioutunit,9)(hist(2,k),k=n+1,ipts)
       write(ioutunit,10)(nint(hist(1,k)),k=n+1,ipts)
      endif
      write(ioutunit,11)
11    format('   -1000'/'  -10000')
      close(ioutunit)
      write(*,'(a)')'Wrote '//outfile(1:len_trim(outfile))
      return; end subroutine pdsout                                         
      subroutine epfout(n)
      use rebin
      use summation
      integer(kind=4) :: i, ilow, ihigh, n
      call filext('.epf',n) ! epfs never refer to a scan, always a sum
      open(unit=ioutunit,status='UNKNOWN',file=outfile)
      ilow=getfirstpoint(hist,2,npts,2)
      ihigh=getlastpoint(hist,2,npts,2)
      if(ilow.eq.-1 .or. ihigh.eq.-1)then
       ! No points with valid data found !
       return
      endif
      do i=ilow,ihigh     
        write(ioutunit,*)bincen(i),hist(1,i),hist(2,i)
      enddo
      close(ioutunit)
      write(*,'(a)')'Wrote '//outfile(1:len_trim(outfile))      
      return
      end subroutine epfout                                               
      subroutine xyeout(n)
      use rebin
      use summation
      integer(kind=4) :: i, ilow, ihigh, n
      character(len=4) :: extn
      if(units .eq. 'T') extn = '.xye'
      if(units .eq. 'Q') extn = '.qye'
      if(units .eq. 'R') extn = '.q2'
      call filext(extn,n) ! epfs never refer to a scan, always a sum
      open(unit=ioutunit,status='UNKNOWN',file=outfile)
      ilow=getfirstpoint(hist,2,npts,2)
      ihigh=getlastpoint(hist,2,npts,2)
      if(ilow.eq.-1 .or. ihigh.eq.-1)then
       ! No points with valid data found !
       return
      endif
      do i=ilow,ihigh     
        write(ioutunit,'(F12.6,2(1X,G14.8))')bincen(i),hist(1,i),hist(2,i)
      enddo
      close(ioutunit)
      write(*,'(a)')'Wrote '//outfile(1:len_trim(outfile))      
      return
      end subroutine xyeout                                               
      subroutine outputinp(n)
      use rebin
      use summation
      integer(kind=4),intent(in)::n
      integer(kind=4) :: i, ilow, ihigh
      character(len=15):: string
      write(string,'(i10)')n
      string=adjustl(string)
      if(units.eq.'T') string=string(1:len_trim(string))//'.inp'
      if(units.eq.'Q') string=string(1:len_trim(string))//'.inq'
      if(units.eq.'R') string=string(1:len_trim(string))//'.inq2'
!      write(*,*)'String was ',string
      open(unit=ioutunit,status='UNKNOWN',file=string)
      ilow=getfirstpoint(hist,2,npts,2)
      ihigh=getlastpoint(hist,2,npts,2)
      if(ilow.eq.-1 .or. ihigh.eq.-1)then
       ! No points with valid data found !
       return
      endif
      do i=ilow,ihigh     
        write(ioutunit,*)bincen(i),hist(1,i),hist(2,i)
      enddo
      close(ioutunit)
      write(*,'(a)')'Wrote '//string(1:len_trim(string)) 
      return
      end subroutine outputinp                                               
      subroutine outputdiagnostic(wd)
      use rebin ! ascan
      use summation
      real(kind=8),intent(in)::wd
      integer(kind=4) :: i, ihigh, ilow, j
      real(kind=8) :: y, ey2, wd2
! write a file in plotmtv format with the 9 channels and sum
      open(unit=ioutunit,file='diag.mtv',status='UNKNOWN',              &
     & form='FORMATTED', access='SEQUENTIAL',iostat=i) 
      if(i.ne.0)then
       write(*,'(a)')'error opening diagnostic file'
      endif
      write(ioutunit,'(a)')'$ DATA = CURVE2D'
      write(ioutunit,'(a)')'% xlabel = "Two Theta"'
      write(ioutunit,'(a)')'% ylabel = "Cts/Monitor"'
      write(ioutunit,'(a)')'% toplabel= "Diagnostic plot"'
      do j=1,nchannel
! j-1 to fix the zeroth channel being first
       write(ioutunit,'(a,i1,a)')'% linelabel = " MA',j-1,' "'      
       write(ioutunit,'(a,i2)')'% linecolor = ',j      
       write(ioutunit,'(a)')'% linetype=1 markertype=0'      
       ilow=getfirstpoint(ascan,2*nchannel,npts,2*j)
       ihigh=getlastpoint(ascan,2*nchannel,npts,2*j)
       do i=ilow,ihigh
        if(ascan(2*j,i).gt.1.0d0)                                       &
        write(ioutunit,'(2F15.8)')bincen(i),                            &
     & ascan(2*j-1,i)/ascan(2*j,i)/mult(j)
       enddo
       write(ioutunit,*)
      enddo
      write(ioutunit,'(a,i2,a)')'% linelabel = "Total"'      
      write(ioutunit,'(a,i2)')'% linecolor = ',nchannel+1      
      write(ioutunit,'(a)')'% linetype=1 markertype=0'
      ilow=getfirstpoint(hist,2,npts,2)
      ihigh=getlastpoint(hist,2,npts,2)
      do i=ilow,ihigh
        write(ioutunit,'(2F15.8)')bincen(i),hist(1,i)
      enddo
      write(ioutunit,*)
      write(ioutunit,'(a,i2,a)')'% linelabel = "outliers"'      
      write(ioutunit,'(a,i2)')'% markercolor = ',nchannel+2      
      write(ioutunit,'(a)')'% linetype=0 markertype=2'
      do i=ilow,ihigh
        if(ctchan(i).ge.2)then
        do j=1,nchannel
          if(ascan(2*j,i).gt.1.0d0)then
            y=ascan(2*j-1,i)/ascan(2*j,i)/mult(j)   ! OK
! assumes zero error in the efficiency    
            ey2= ((ascan(2*j-1,i)+alp)/ascan(2*j,i)**2                  &
     &           +  ascan(2*j-1,i)**2/ascan(2*j,i)**3)/mult(j)**2
            wd2=(y-hist(1,i))**2/(ey2+hist(2,i)**2)                          
            if(wd2.gt.wd)write(ioutunit,'(2F15.8)')bincen(i),y
          endif
        enddo
        endif
      enddo
      write(ioutunit,*)
      write(ioutunit,'(a)')'$ END'
      write(*,'(a,G12.5)')'Wrote diag.mtv file, outliers at ',sqrt(wd)
      close(ioutunit)
      return
      end subroutine outputdiagnostic                                        
      subroutine outputdiagnosticw32(wd)
      use rebin ! ascan
      use summation
      real(kind=8),intent(in)::wd
      integer(kind=4) :: i, j, k
      real(kind=8) :: y, ey2, wd2
! write a file in plotmtv format with the 9 channels and sum
      open(unit=ioutunit,file='diag.mtv',status='UNKNOWN',              &
     & form='FORMATTED', access='SEQUENTIAL',iostat=i) 
      if(i.ne.0)then
       write(*,'(a)') 'error opening diagnostic file'
       stop
      endif
      write(ioutunit,'(a)')'#@set[0].point.style none'
      write(ioutunit,'(a)')'#@set[1].point.style none'
      write(ioutunit,'(a)')'#@set[2].point.style none'
      write(ioutunit,'(a)')'#@set[3].point.style none'
      write(ioutunit,'(a)')'#@set[4].point.style none'
      write(ioutunit,'(a)')'#@set[5].point.style none'
      write(ioutunit,'(a)')'#@set[6].point.style none'
      write(ioutunit,'(a)')'#@set[7].point.style none'
      write(ioutunit,'(a)')'#@set[8].point.style none'
      write(ioutunit,'(a)')'#@set[9].point.style none'
      write(ioutunit,'(a)')'#@set[10].point.style circle'
      write(ioutunit,'(a)')'#@set[0].line.style solid'
      write(ioutunit,'(a)')'#@set[1].line.style solid'
      write(ioutunit,'(a)')'#@set[2].line.style solid'
      write(ioutunit,'(a)')'#@set[3].line.style solid'
      write(ioutunit,'(a)')'#@set[4].line.style solid'
      write(ioutunit,'(a)')'#@set[5].line.style solid'
      write(ioutunit,'(a)')'#@set[6].line.style solid'
      write(ioutunit,'(a)')'#@set[7].line.style solid'
      write(ioutunit,'(a)')'#@set[8].line.style solid'
      write(ioutunit,'(a)')'#@set[9].line.style solid'
      write(ioutunit,'(a)')'#@set[10].line.style none'
      write(ioutunit,'(a)')'#@set[0].line.color custom 255 0   0'
      write(ioutunit,'(a)')'#@set[1].line.color custom 215 30  0'
      write(ioutunit,'(a)')'#@set[2].line.color custom 180 60  0'
      write(ioutunit,'(a)')'#@set[3].line.color custom 150 90  0'
      write(ioutunit,'(a)')'#@set[4].line.color custom 120 120 0'
      write(ioutunit,'(a)')'#@set[5].line.color custom 90  150 0'
      write(ioutunit,'(a)')'#@set[6].line.color custom 60  180 0'
      write(ioutunit,'(a)')'#@set[7].line.color custom 30  215 0'
      write(ioutunit,'(a)')'#@set[8].line.color custom 0   25  0'
      write(ioutunit,'(a)')'#@set[9].line.color custom 0   0   0'
      write(ioutunit,'(a)')'#@set[10].point.color custom 0 0 0'
      do i=1,npts
      if(ctchan(i).lt.1)cycle
      write(ioutunit,'(11F15.8)',advance='no')bincen(i),                &
     &(ascan(2*j-1,i)/(1.0d0+ascan(2*j,i)*mult(j)), j=1,nchannel),      &
     &hist(1,i)
      k=0
      if(ctchan(i).ge.2)then  
        do j=1,nchannel
          if(ascan(2*j,i).gt.1.0d0)then
            y=ascan(2*j-1,i)/ascan(2*j,i)/mult(j)   ! OK
! assumes zero error in the efficiency    
            ey2= ((ascan(2*j-1,i)+alp)/ascan(2*j,i)**2                  &
     &           +  ascan(2*j-1,i)**2/ascan(2*j,i)**3)/mult(j)**2
            wd2=(y-hist(1,i))**2/(ey2+hist(2,i)**2)                          
            if(wd2.gt.wd .and. k.eq.0)then
              write(ioutunit,'(F15.8)')                                 &
     &            ascan(2*j-1,i)/(ascan(2*j,i)*mult(j))   
              k=1
            endif
          endif
        enddo
      endif
      if(k.eq.0)write(ioutunit,'(F15.8)')-0.001d0
      enddo
      write(ioutunit,*)
      write(*,'(a,G12.5)')'Wrote diag.mtv file, outliers at ',sqrt(wd)
      close(ioutunit)
      return; end subroutine outputdiagnosticw32                             
      integer(kind=4) function getfirstpoint(array,m,n,l)
      integer(kind=4), intent(in) :: l, m, n
      real(kind=8), intent(in), dimension(m,n) :: array
      integer(kind=4) :: i
      getfirstpoint=1
      do i=1,n
       if(array(l,i).gt.0.0d0)then
        getfirstpoint=i
        return
       endif
      enddo
      return
      end function getfirstpoint                                             
      integer(kind=4) function getlastpoint(array,m,n,l)
      integer(kind=4), intent(in) :: m,n,l
      real(kind=8), intent(in), dimension(m,n) :: array
      integer(kind=4) :: i
      getlastpoint=n
      do i=n,1,-1
       if(array(l,i).gt.0.0d0)then
        getlastpoint=i
        return
       endif
      enddo
      return
      end function getlastpoint                                              
      subroutine dumpscan(is)
      use rebin
      integer(kind=4) :: i, j
      integer(kind=4),intent(in)::is
      real(kind=8),dimension(9) :: sig
      call filext('.dum',is)
      open(unit=ioutunit,file=outfile,status='UNKNOWN',                 &
     & form='FORMATTED',access='SEQUENTIAL',iostat=i)
      if(i.ne.0) stop 'error opening output file'
      do i=1,npts
        do j=1,9
          if(ascan(2*j,i).gt.0)then
            sig(j)=ascan(2*j-1,i)/ascan(2*j,i)
          else
             sig(j)=0.0
          endif
        enddo
        write(ioutunit,'(10(F12.8,1X))')bincen(i),(sig(j),j=1,9)
      enddo
      close(ioutunit)
      return  
      end subroutine dumpscan                                                
      subroutine filext(extn,is)
      use specfiles
      character(len=4),intent(in) :: extn
      integer(kind=4),intent(in) :: is
      integer(kind=4) :: n, i
      character(len=15):: string
      if(is.gt.0)then
       write(string,'(i10)')is ! scan number into string
       string=adjustl(string)  ! shift to left end of string
      endif
      n=len_trim(filnam)
      if(filnam(n-3:n).eq.'.dat')n=n-4 ! overwrite .dat if exists
      if(is.gt.0) then
       write(outfile,'(a)')filnam(1:n)//'_'//string(1:len_trim(string)) &
     & //extn(1:4)
      else
       write(outfile,'(a)')filnam(1:n)//extn(1:4)
      endif
! strip any / or \ from the start of the filename so that all output
! is in the working directory.
      n=len_trim(outfile)
      do i=n,1,-1
        if(outfile(i:i).eq.'\' .or. outfile(i:i).eq.'/')then
        outfile=outfile(i+1:n); exit; endif; enddo
      end subroutine filext                                                  
      subroutine openlogfile
      integer(kind=4) :: i
      call filext('.log',0) ! logfiles for all scans, not one each
      open(unit=logfile,file=outfile,status='UNKNOWN',                  &
     & form='FORMATTED',access='SEQUENTIAL',iostat=i)
      if(i.ne.0)stop 'could not open logfile'
      return; end subroutine openlogfile                                     
      subroutine wlogfile(s)
      character(len=*) :: s
      write(logfile,'(a)')s
      return; end subroutine wlogfile                                        
      subroutine bcmfile(n)
      use rebin ! for ascan and mult
      use summation
      use specfiles
      implicit none
      integer(kind=4) :: n,ma0
      character(len=WORDLENGTH) :: c
      integer(kind=4) :: ilow,ihigh,i,j     
      call filext('.bcm',n) 
      open(unit=ioutunit,status='UNKNOWN',file=outfile)
      ilow=getfirstpoint(hist,2,npts,2) ! range of real points
      ihigh=getlastpoint(hist,2,npts,2)
      c=columnlabels(1)
      ma0=whichcolumn(FIRSTDET) ! FIXME deal with absent dets
      write(ioutunit,'(a,2x,a,$)')'#',c(1:len_trim(c))
      do i=0,nchannel-1
        if(logexdet(i+1).eq.1)cycle ! skip if excluded completely
        c=columnlabels(ma0+i)
        write(ioutunit,'(2X,a,2x,a,$)')c(1:len_trim(c)),                      &
     & c(1:len_trim(c))//'_mon'                                  
      enddo
      write(ioutunit,*) ! end of line
      do i=ilow,ihigh
       write(ioutunit,'(G14.8,2X,$)')bincen(i)
       do j=1,nchannel
        if(logexdet(j).eq.1)cycle ! skip if excluded completely         
         write(ioutunit,'(2(G14.8,2X),$)')ascan(2*j-1,i),ascan(2*j,i)
       enddo
       write(ioutunit,*) ! end of line
      enddo
      close(ioutunit)
      write(*,'(a)')'Wrote '//outfile(1:len_trim(outfile))
      return
      end subroutine bcmfile                                                  
      subroutine rstset(string)
      use rebin
      character(len=*), intent(in) :: string
      integer itok, iprev, itot, ichan, i
      itot = len_trim(string)
      itok = index( string(5:itot), ",")
      if (itok.eq.0) then
         read(string(5:itot),*,err=10) randomstart
         do i = 1, nchan
           rstchan(i) = .true.
         enddo
         goto 100
      else
         read(string(5:5+itok-2),*,err=10) randomstart
         iprev = 5+itok
         do
            itok = index( string(iprev:itot), ",")
            if (itok.eq.0) then
               if(iprev.le.itot)then
                  read( string(iprev:itot), *, err=10, end=10) ichan
                  call rstadd(ichan)
               endif
               exit ! loop
            endif
            read( string(iprev:iprev+itok-1), *, err=10, end=10) ichan
            call rstadd(ichan)         
            iprev = iprev+itok
         enddo
      endif

      goto 100
10    write(*,*)itok,iprev,itot,string
      STOP 'Could not understand your rst=x.xxx request'
100   userandomstart = .true.
      write(*,'(A,1X,f7.5,1X,A,$)')'Applying random start of',randomstart, &
     & 'to channels:'
      do i = 1,nchan
         if (rstchan(i)) write(*,'(1X,I2,$)') i-1
      enddo
      write(*,*)
      return; end subroutine rstset                                       
      
      subroutine rstadd(ichan)
      use rebin
      integer ichan
      if ((ichan .lt. 0).or. (ichan .ge. NCHAN)) then
        write(*,*) 'rst channel out of range',ichan
        stop
      endif
! the +1 makes it go from zero
      rstchan(ichan+1) = .true.
      return
      end subroutine rstadd

      subroutine renset(string)
      use rebin
      character(len=*), intent(in) :: string
      integer itok, iprev, itot, ichan, i
      itot = len_trim(string)
      itok = index( string(5:itot), ",")
      if (itok.eq.0) then
         read(string(5:itot),*,err=10) randomend
         do i = 1, nchan
           renchan(i) = .true.
         enddo
         goto 100
      else
         read(string(5:5+itok-2),*,err=10) randomend
         iprev = 5+itok
         do
            itok = index( string(iprev:itot), ",")
            if (itok.eq.0) then
               if(iprev.le.itot)then
                  read( string(iprev:itot), *, err=10, end=10) ichan
                  call renadd(ichan)
               endif
               exit ! loop
            endif
            read( string(iprev:iprev+itok-1), *, err=10, end=10) ichan
            call renadd(ichan)         
            iprev = iprev+itok
         enddo
      endif

      goto 100
10    write(*,*)itok,iprev,itot,string
      STOP 'Could not understand your ren=x.xxx request'
100   userandomend = .true.
      write(*,'(A,1X,f7.5,1X,A,$)')'Applying random end of',randomend, &
     & 'to channels:'
      do i = 1,nchan
         if (renchan(i)) write(*,'(1X,I2,$)') i-1
      enddo
      write(*,*)
      return; end subroutine renset                                       
      
      subroutine renadd(ichan)
      use rebin
      integer ichan
      if ((ichan .lt. 0).or. (ichan .ge. NCHAN)) then
        write(*,*) 'ren channel out of range',ichan
        stop
      endif
! the +1 makes it go from zero
      renchan(ichan+1) = .true.
      return
      end subroutine renadd

      end module outputfiles                                                 
      subroutine mma(arr,n,mode)
      use specfiles
      use outputfiles
      implicit none
      integer(kind=4),intent(in) :: mode,n
      real(kind=8),dimension(n),intent(in) :: arr
      real(kind=8),allocatable, dimension(:,:),save :: mmasum
      real(kind=8) :: rn
      integer(kind=4) :: i
      integer(kind=4),save :: nt, nwas
      character(len=80) :: s
      if(mode.eq.2 .and. nt.ne.0)then ! write to log file
        rn=real(nt,8) 
        do i=1,nwas
        mmasum(i,1)=mmasum(i,1)/rn
! The absolute value is to avoid problems with rounding errors when there
! is zero variance in the data. eg: blower records a constant when the
! serial line is unplugged so we can end up trying to take a square root 
! of a very tiny negative number, instead of a tiny positive number.
        mmasum(i,2)=sqrt(abs(rn*(mmasum(i,2)/rn-mmasum(i,1)**2)/(rn-1.0d0)))
        write(s,1000)iscan,columnlabels(i),mmasum(i,1),mmasum(i,2),nt
1000    format('Scan ',i5,' Ctr ',a10,' Avg =',F10.4,' +/- ',           &
     &     F10.4,' npts = ',i10)
        call wlogfile(s)
        write(s,1001)iscan,columnlabels(i),mmasum(i,4),mmasum(i,3)
1001    format('Scan ',i5,' Ctr ',a10,' T max = ',F10.4,' T min = ',F10.4)
        call wlogfile(s)
        enddo
        nt=0
        deallocate(mmasum)
      endif
      if(mode.eq.1)then
        if(.not. allocated(mmasum))then 
          allocate(mmasum(ncolumns,4))
          mmasum=0.0d0
          nt=0
          nwas=ncolumns
        endif
      mmasum(:,1)=mmasum(:,1)+arr                     ! for avg
      mmasum(:,2)=mmasum(:,2)+arr*arr                 ! for std dev
      if(nt.eq.0)mmasum(:,3)=arr
      if(nt.eq.0)mmasum(:,4)=arr
      nt=nt+1                             ! no of points
      if(nwas.eq.ncolumns)then
       do i=1,ncolumns
        if(mmasum(i,3).gt.arr(i))mmasum(i,3)=arr(i)   ! min vals
        if(mmasum(i,4).lt.arr(i))mmasum(i,4)=arr(i)   ! max vals
       enddo
      endif
      endif
      return; end subroutine mma                                           
      subroutine logmotors
      use specfiles
      use outputfiles
      integer(kind=4) :: i, j
      character(len=80):: s      
      do j=1,NLINES 
        if(headerwords(j).gt.0)then
          do i=1,headerwords(j) 
            write(s,1000)iscan,words(i,j),wordvalues(i,j)
1000        format('Scan ',i5,1X,a20,1X,f16.8)
            call wlogfile(s)
          enddo
        endif ! headerwords(j.gt.0)
      enddo
      write(s,1001)iscan, ncolumns
1001  format('Scan ',i5,' number of columns = ',i5)
      call wlogfile(s)
      do j=1,ncolumns
        write(s,1002)iscan,j,columnlabels(j)
1002    format('Scan ',i5,' column ',i3,1x,a30)
        call wlogfile(s)
      enddo
      return; end subroutine logmotors                                       
      module pointfilter
      logical :: filterlogical=.false.
      real(kind=8), dimension(10,3) :: lt
      real(kind=8) :: tth, tthnew
      integer(kind=4) :: nc=0, np=1
      contains
      subroutine filterinit()
      lt(:,1)=0.0
      lt(:,2)=0.0
      lt(:,3)=0.0
      tth=-9999.0
      np=1
      return; end subroutine filterinit
      subroutine pf(a,itth,MA0,MA8,M)
      integer(kind=4), intent(in) :: itth, MA0, MA8, M
      real(kind=8), intent(inout), dimension(:) :: a
      integer :: i, j
      real(kind=8) :: tthnew
      tthnew=a(itth)
      a(itth)=tth
      tth=tthnew ! swap for last point
      lt(1:9,np)=a(MA0:MA8)
      lt(10,np)=a(M)
      do i = MA0, MA8
         j = i-MA0+1
         a(i)=lt(j,middle(lt(j,:)))
      enddo
      a(M)=lt(10,middle(lt(10,:)))
      np=np+1
      if(np.eq.4)np=1 ! Rolling buffer
      end subroutine pf
      integer function middle(x)
      real(kind=8), dimension(3), intent(in) :: x
      integer,dimension(1) :: i,j
      i=maxloc(x)
      j=minloc(x)
      middle=6-i(1)-j(1) ! 1+2+3 == 6
      if(middle.gt.3)middle=3
      end function middle  
      end module pointfilter                                           
      subroutine processscan(n)
      use specfiles     
      use rebin
      use outputfiles    
      use pointfilter
      integer(kind=4),intent(in) :: n
      real(kind=8) :: tth1, tth2, tthf, x
      real(kind=8), allocatable :: chantot(:), a(:) ! ncolumns
      integer(kind=4) :: mon, ma0, ma8, itth, i
      integer(kind=4) :: igoodpoints, ibadpoints
      character(len=80) :: s
      call findscan(n)   ! Positions the file on that scan 
      call logmotors  ! dumps all the starting motor positions to the log file
      if(ispecerr.ne.0) return
      write(*,'(i5,$)')iscan
      if(scantype.ne.'turboscan' .and. scantype.ne.'hookscan'               &
     &   .and. scantype.ne.'cscan' .and. scantype.ne.'zapline')then
        write(*,'(a)')' is not a turboscan or a hookscan, ignoring it'
        ispecerr=-1
        return
      endif
      tth1=getheadervalue('2_theta')
      tthf=tth1
      mon=-1;ma0=-1;itth=-1;itth=-1    ! FIXME
      mon=whichcolumn(MONITORCOL)
      if(mon.lt.1)goto 2
      ma0=whichcolumn(FIRSTDET) ! FIXME deal with absent dets
      if(ma0.lt.1)goto 2
      ma8=whichcolumn(LASTDET) ! FIXME deal with absent dets.      
      itth=whichcolumn(TWOTTH)
!      write(*,*)"columns",mon,ma0,ma8,itth
      if(mon.lt.1)goto 2
! Allocate a and chantot if necessary (cannot do this till after findscan)
      if( allocated(a) ) deallocate(a) ! resets
      if( allocated(chantot) ) deallocate(chantot)
      allocate(a(ncolumns))
      allocate(chantot(ncolumns))
      chantot=0.0d0 ; igoodpoints=0 ; ibadpoints=0; a=0.0d0
! Set up any unit conversions which are requested
!   Send the current Q line from the specfile module to the rebin module
      call convert_unit_setupQ( Q )
! Skip the first line of data
      call getdata(a,ncolumns)
      if(wincnt(1:6) .eq. 'Epoch') then
         write(*,*)'Scan Epoch starts at',a(wincol)
         winlow = winlowread+a(wincol)
         winhigh = winhighread+a(wincol)
         write(*,*)'Adjusting your Epoch window to',winlow,winhigh
      endif
      tth1=a(itth)
      if(filterlogical)then  ! initialise
        call filterinit()
        do i=1,3
          call getdata(a,ncolumns)
          if(ispecerr.ne.0)goto 2
          call pf(a,itth,ma0,ma8,mon)
          tth1=a(itth)
        enddo
      endif
      if( userandomstart ) then
        randomval = randomstart * rlcg() + tth1
        write(*,*)'rst using',randomval,'from',tth1
      endif
      if( userandomend ) then
        randomvalend = scanend - randomend * rlcg() 
        write(*,*)'ren using',randomvalend,'from',scanend
      endif
! Main loop
1     call getdata(a,ncolumns)
      if(filterlogical) call pf(a,itth,ma0,ma8,mon)
      do i=0,nchannel-1
        if(a(ma0+i).lt.0)then
         a(mon)=0.0d0 ! set monitor to zero if negative
! why not set all channels to zero for negative counts??
         write(*,*)'negative counts found for line'
         write(*,*)line(1:len_trim(line))
        endif
      enddo
      if(ispecerr.eq.0)then
        tth2=a(itth)  
        if(a(mon).gt.minmon .and. window(a,ncolumns)) then
         igoodpoints=igoodpoints+1
         chantot=chantot+a ! sum on totals
         call processline(tth1,tth2,a(ma0:ma8),nchannel,a(mon))
         call mma(a,ncolumns,1)
        else
         ibadpoints=ibadpoints+1
        endif
        tth1=tth2
        goto 1            ! loops here
      endif
2     continue
      write(s,1000)iscan,chantot(mon)
1000  format('Scan ',i5, ' Total monitor ',F25.0)
      call wlogfile(s)
      write(s,1001)iscan,igoodpoints,igoodpoints+ibadpoints
1001  format('Scan ',i5, ' used ',i20,' of ',i20,' points' )
      call wlogfile(s)
      do i=ma0,ncolumns
       write(s,1002)iscan,columnlabels(i),chantot(i)
1002   format('Scan ',i5,1x,a20,' total ',F25.0)      
       call wlogfile(s)
      enddo
      write(s,1003)iscan,chantot( whichcolumn('Seconds') )
1003  format('Scan ',i5,1x,'Total time          ',F25.0,' seconds')
      call wlogfile(s)
! Average step size = total range / total points
      x=abs((tthf-tth1)/real((igoodpoints+ibadpoints),8))
      write(s,1004)iscan,x
1004  format('Scan ',i5,1x,'Average stepsize',E25.15)
      call wlogfile(s)
! no of bins = step/(avg step)
      write(s,1005)iscan,step/x
1005  format('Scan ',i5,1x,'Average points per bin',E25.10)
      call wlogfile(s)                                                     
      write(s,1006)iscan,tthf,tth1
1006  format('Scan ',i5,1x,'low tth',E25.10,1x,'last tth',E25.10)
      call wlogfile(s)
      write(s,1007)iscan,igoodpoints,ibadpoints
1007  format('Scan ',i5,1x,'gd',i15,1x,'bad',i15)
      call wlogfile(s)
      call mma(a,ncolumns,2)                                           
      return
      end subroutine processscan                                             
      module useroptions
      integer(kind=4) :: ifirstscan, ilastscan, icurrscan=0,nexcld
      integer(kind=4),allocatable :: exscanlist(:),exdetlist(:)
      real(kind=8)::wd=36.0d0
      logical :: snbl=.false.
      contains
      subroutine getcmdline
      use rebin
      use specfiles
      integer(kind=4) :: iarg, i
      integer(kind=4),external :: iargc
      character(len=256) :: string ! massive string in case of silly user
      step=0.001d0; ifirstscan=0; ilastscan=0
      iarg=iargc()
      if(iarg.eq.0)then ; call helpmsg ; endif
      call getarg(1,filnam) ! get's filename
      if(iarg.eq.1) then
        write(*,'(a)')'No stepsize, assuming 0.001, and all scans'
        goto 100 ! return
      endif
      if(iarg.ge.2)then
        call getarg(2,string)   ! get the stepsize
        read(string,*,err=10,end=10)step
        user_step = step
      endif
      if(iarg.ge.3)then      
        call getarg(3,string)   ! get the first scan to bin
        read(string,*,err=10,end=10)ifirstscan
        if(ifirstscan.eq.0)write(*,*)                                   &
     &'Zero for first scan is going to cause problems... sorry' ! FIXME?
      endif   
      if(iarg.ge.4)then
        call getarg(4,string)   ! get last scan to bin
        read(string,*,err=10,end=10)ilastscan
      endif
      if(iarg.ge.5)then
        do i=5,iarg
          call getarg(i,string); call option(string)
        enddo
      endif
      goto 100
10    write(*,*)'Problems interpreting your command line'
      call helpmsg
100   return
      end subroutine getcmdline                                              
      subroutine option(string)
      use outputfiles ! logicals for which files to write
      character(len=*),intent(in) :: string
      if(string(1:1).ne.' ') then
       select case (string(1:3))
        case('mon') ; call setmonitorcol(string)
        case('wvl') ; call setwavelength(string)
        case('uni') ; call setunits(string)
        case('ed=') ; call exdet(string)
        case('es=') ; call exscan(string)
        case('is=') ; call incscan(string)
        case('ef=') ; call exfile(string)
        case('low') ; call lowtth(string)
        case('hig') ; call hightth(string)
        case('sca') ; call scale(string)
        case('ste') ; call minstepset(string)
        case('wd=') ; call wdset(string)
        case('alp') ; call alpset(string)
        case('ren') ; call renormset(string)
        case('mm=') ; call minmonset(string)
        case('mr=') ; call minrenormset(string)
        case('win') ; call windowset(string)
        case('zap') ; call zapset(string)
        case('sup') ; call superzapset(string)
        case('med') ; call medianchannelset(string)
        case('3pf') ; call filterset(string) 
        case('bcm') ; call bcmset(string)
        case('snb') ; call snblset(string)      
        case('rst') ; call rstset(string)
        case('rnd') ; call renset(string)
        case('nod') ; if(string(1:6).eq.'nodiag')then 
          diag=.false. ; else ;
          write(*,*)'Sorry, I did not understand the command line'
          write(*,*)string(1:len_trim(string)) ; endif
        case('gsa') ; gsas=.true.
        case('spf') ; spf=.true.
        case('pds') ; pds=.true.
        case('epf') ; epf=.true.
        case default
          write(*,*)'Sorry, I did not understand the command line'
          write(*,*)string(1:len_trim(string))
          call helpmsg
        end select
      endif
      return; end subroutine option                                          
      subroutine exdet(string)
! Read a comma separated list of detectors to ignore for the
! final sum
      use rebin
      character(len=*), intent(in) :: string
      integer(kind=4) :: i, j
      i=ncommas(string)+1
      allocate(exdetlist(i))
      read(string(4:len(string)),*,err=10)exdetlist
      do j=1,i
! Added a plus one here to go from channel zero
       logexdet(exdetlist(j)+1)=1
       write(*,'(a,i2)')'Intending to exclude channel ',exdetlist(j)
      enddo
      goto 100
10    STOP 'Could not understand your list of excluded detectors'
100   return; end subroutine exdet                                        
      subroutine exscan(string)
! Read a comma separated list of scans to skip over when summing
! Place these in an array which only nextscan knows about
      character(len=*),intent(in) :: string
      nexcld=ncommas(string)+1
      allocate(exscanlist(nexcld))
      read(string(4:len(string)),*)exscanlist
      return;  end subroutine exscan                                       
      subroutine incscan(string)
! Read a comma separated list of scans to skip over when summing
! Place these in an array which only nextscan knows about
      character(len=*),intent(in) :: string
      integer(kind=4),allocatable :: temp(:)
      integer(kind=4) :: i, j, nincld
      nincld=ncommas(string)+1
      allocate(temp(nincld))
      read(string(4:len(string)),*)temp
! total = ilastscan - ifirstscan + 1   ... eg 2 -> 10 is 2,3,4,5,6,7,8,9,10
!                                                        1 2 3 4 5 6 7 8 9
      nexcld=ilastscan - ifirstscan + 1 - nincld
      allocate(exscanlist(nexcld)) 
      i=1
      
      do j=ifirstscan,ilastscan
         if (.not. inlist(j,temp,nincld)) then
           exscanlist(i)=j
           i=i+1
         endif
      enddo
      if (i-1 .ne. nexcld) then
         write(*,*)' problem - debug is= please'
      endif
      deallocate(temp)
      return;  end subroutine incscan                                      
      subroutine exfile(string)
! Given a filename in string read in the file and set up for 
! excluding regions of channels
      use rebin
      character(len=*), intent(in) :: string
      integer(kind=4)i,ier,ns
      real(kind=8)xl,xh
      open(unit=29,file=string(4:len_trim(string)),                     &
     & status='OLD',iostat=ier)
      if(.not.(ier.eq.0))then
        write(*,'(a,i5)')'Error opening your file '//                   &
     & string(4:len_trim(string))//' iostat=',ier
      write(*,'(a)')                                                    &
     &'want exclusion file, lines must have channel lowtth hightth scan'
      write(*,'(a)')'Make scan number less than zero for all'
        stop
      endif
      write(*,'(a)')                                                    &
     &'Reading exclusion file, lines must have channel lowtth hightth scan'
      write(*,'(a)')'Make scan number less than zero for all'
1     read(29,*,end=2,err=2)i,xl,xh,ns
      if((i.ge.0).and.(i.le.NCHAN)) then
         if(logexdet(i+1).eq.0) then
            logexdet(i+1)=2 ! set that this has excluded regions
         else
            write(*,'(a,i3,a)')'Channel',i,' already excluded'
         endif
      else
         write(*,'(a,i3,a)')'Error in your exfile, channel',i,          &
     &      'not allowed'  
         stop ! Operator error - give up
      endif
      iexrc=iexrc+1 ! increment number of excluded regions
      goto 1
2     allocate(exarray(iexrc,3)) ! holds low tth, high tth pairs
      allocate(iexarray(iexrc))  ! holds channel number
      rewind(29)
      do i=1,iexrc
        read(29,*)iexarray(i),exarray(i,1:2),exarray(i,3)
        write(*,'(a,i2,a,f10.6,a,f10.6)')'Excluding channel ',          &
     & iexarray(i),' from ', exarray(i,1),' to ',exarray(i,2)
        if(exarray(i,3).lt.1) then
          write(*,*)"in all scans"
        else
          write(*,'(a,i4)')'in scan',int(exarray(i,3))
        endif
      enddo
      close(29)
      return             
      end subroutine exfile

      integer(kind=4) function ncommas(string)
      character(len=*),intent(in) :: string
      integer(kind=4) :: i
      ncommas=0
      do i=1,len(string)
        if(string(i:i).eq.',')ncommas=ncommas+1
      enddo
      return; end function ncommas                                         
      subroutine lowtth(s)
      use rebin
      character(len=*),intent(in) :: s
      if(s(1:7).eq.'lowtth=')read(s(8:len_trim(s)),*)tthlow
      user_tthlow = tthlow
      return ; end subroutine lowtth
      subroutine hightth(s)
      use rebin
      character(len=*),intent(in) :: s
      if(s(1:8).eq.'hightth=')read(s(9:len_trim(s)),*)tthhigh
      user_tthhigh = tthhigh
      return; end subroutine hightth                                       
      subroutine renormset(s)
      use summation
      character(len=*)::s
      if(s(1:6).eq.'renorm') renorm=.true.
      return; end subroutine renormset                                    
      subroutine wdset(s)
      character(len=*),intent(in)::s
      if(s(1:3).eq.'wd=')read(s(4:len_trim(s)),*)wd
      wd=wd*wd
      return ; end subroutine wdset                                        
      subroutine minstepset(s)
      use rebin
      character(len=*),intent(in):: s
      if(s(1:5).eq.'step=')read(s(6:len_trim(s)),*)aminstep
      user_step = aminstep
      return; end subroutine minstepset                                    
      subroutine minmonset(s)
      use rebin
      character(len=*),intent(in) :: s
      if(s(1:3).eq.'mm=')read(s(4:len_trim(s)),*)minmon
      return; end subroutine minmonset                                     
      subroutine minrenormset(s)
      use rebin
      character(len=*),intent(in) :: s
      if(s(1:3).eq.'mr=')read(s(4:len_trim(s)),*)minrenormsig
      write(*,1)"Only using where all channels have more than ",minrenormsig, &
     & " counts for renorm"
1     format(a,f8.2,a) 
      return; end subroutine minrenormset                                  
      subroutine windowset(s)
      use rebin
      character(len=*),intent(in) :: s
      character(len=256) :: mys
      integer i
      if(s(1:7).eq.'window=')then
         do i=8,len_trim(s)
           if(s(i:i).eq.',')then
              mys(i:i)=' '
           else
              mys(i:i)=s(i:i)
           endif
         enddo
         read(mys(8:len_trim(s)),*)wincnt,winlow,winhigh
      endif
      write(*,*)'windowing counter ',wincnt(1:len_trim(wincnt)),' from ', &
     & winlow, ' to ', winhigh
      winhighread=winhigh
      winlowread=winlow
      winlog=.true.
      return
      end subroutine windowset                                             
      subroutine scale(string)
      use summation
      character(len=*),intent(in) :: string
      if(string(1:6).eq.'scalpk' )then ; isc=1 ; return ;endif
      if(string(1:7).eq.'scaltot')then ; isc=2 ; return ;endif
! Replace scalinp string with scalmon string
      if(string(1:8).eq.'scalmon=')then 
       read(string(9:len_trim(string)),*)scalinp;isc=0
       write(*,'(a,E12.6,a)')'Units will be counts per ',scalinp,       &
     &  ' monitor ct' ; return ; endif
      if(string(1:7).eq.'scalmon')then ; isc=0 ; return ;endif
      write(*,'(a)')'Didn''t understand your optional argument:'
      write(*,'(a)')string(1:len_trim(string))
      return
      end subroutine scale                                                 
      integer(kind=4) function nextscan()      
! Returns to next scan to process or -1 for all done
1     if(ifirstscan.eq.0 .and. ilastscan.eq.0) then          ! all scans
         icurrscan=icurrscan+1
         nextscan=icurrscan
         return  ! exit here !
      elseif (ilastscan.eq.0) then                       ! one scan only
         if(icurrscan.eq.ifirstscan)then
           icurrscan=-1 ! stop
         else
           icurrscan=ifirstscan
           nextscan=ifirstscan
           return ! exit here !
         endif           
      else
        if(icurrscan.eq.0)then       ! needs to deal with ifirstscan=0
           icurrscan=ifirstscan
        else
           icurrscan=icurrscan+1
        endif
      endif
      if(allocated(exscanlist))then ! check if this is an excluded scan
        if(inlist(icurrscan,exscanlist,nexcld)) goto 1
      endif   
      if(icurrscan.le.ilastscan)then
         nextscan=icurrscan
      else
         nextscan=-1
      endif
      return
      end function nextscan                                                  
      logical function inlist(i,iarr,n)
      integer(kind=4),intent(in) :: i, n
      integer(kind=4),dimension(n),intent(in)::iarr
      integer(kind=4) :: j
      inlist=.FALSE.
      do j=1,n
        if(iarr(j).eq.i) inlist=.TRUE.
      enddo
      return; end function inlist                                            
      subroutine alpset(s)
      use summation
      character(len=*), intent(in) :: s
      if(s(1:4).eq.'alp=')read(s(5:len_trim(s)),*)alp
      return; end subroutine alpset                                        
      subroutine zapset(s)
      use summation
      character(len=*), intent(in) :: s
      if(s(1:4).eq.'zap=')then
        read(s(5:len_trim(s)),*)zap
        zapping=.true.
        write(*,*)"Aiming to zap points which off by ",zap," sigma"
      else
        write(*,*)"did not understand your option ",s
      endif
      return; end subroutine zapset                                        
      subroutine superzapset(s)
      use summation
      character(len=*),intent(in) :: s
      if(s(1:9).eq.'superzap=')then
        read(s(10:len_trim(s)),*,err=1)superzaplevel
        superzap=.true.
        write(*,'(A,F8.5)')"Superzapping at level ",superzaplevel
      else 
        write(*,*)'Do you mean to say "superzap=0.9" ? '
        write(*,*)'So say it then, how else can I read the number!!'
        goto 1
      endif
      return
1     write(*,*)'Could not figure out what you meant by'
      write(*,*)s(1:len_trim(s))
      end subroutine superzapset                                  
      subroutine medianchannelset(s)
      use summation
      character(len=*),intent(in) :: s
      if(s(1:len_trim(s)).eq.'medianofchannels')then
        medianofchannels=.TRUE.
        write(*,'(A)')"Signal and statistics reflect the median", &
     & " channel * number_of_active_channels"
      else
        goto 1
      endif
      return
1     write(*,*)'Could not figure out what you meant by'
      write(*,*)s(1:len_trim(s))
      write(*,*)'Nearest guess is medianofchannels ??'
      write(*,*)'IGNORING:',s(1:len_trim(s))
      end subroutine medianchannelset                                  
      subroutine filterset(s)
      use pointfilter
      character(len=*), intent(in) :: s
      if(s(1:3).eq.'3pf')filterlogical=.true.
      write(*,*)"Aiming to apply a 3 point median filter"
      return; end subroutine filterset                                    
      subroutine bcmset(s)
      use outputfiles
      character(len=*)::s
      if(s(1:6).eq.'bcm') bcm=.true.
      return; end subroutine bcmset                                       
      subroutine snblset(s)
      use rebin
      character(len=*)::s
      if(s(1:4).eq.'snbl') snbl=.true.
      nchannel=6
      logexdet(7)=1
      logexdet(8)=1
      logexdet(9)=1
      write(*,'(a,i2)')"Exclude last 3 channels, as they don't exist"
      return; end subroutine snblset                                      

      subroutine setmonitorcol(s)
      use specfiles
      character(len=*), intent(in) :: s
      if(s(1:4).eq."mon=")then
         MONITORCOL=s(5:len_trim(s))
         write(*,*)"Using monitor counts from column ",MONITORCOL
      else
         write(*,*)"I do not understand your argument",s(1:len_trim(s))
         write(*,*)"Perhaps you mean mon=Monitor ?"
      endif
      end subroutine setmonitorcol
      


! These will end up in the useroptions module
!        case('wvl')  call setwavelength(string)
!
      subroutine setwavelength(s)
      use rebin
      character(len=*), intent(in) :: s
      if(s(1:5).eq.'wvln=')read(s(6:len_trim(s)),*,err=1)wavelength
      wavelength_set = .true.
      write(*,*)'Using wavelength of ',wavelength
      return
1     write(*,*)'Error reading wavelength',s
      stop
      end subroutine setwavelength
!
!        case('uni')  call setunits(string)
      subroutine setunits(s)
      use rebin
      character(len=*), intent(in) :: s
!                      12345678
      if( s(1:8) .eq. 'units=Q2' ) then
          units = 'R'
          write(*,*)'Binning into Q^2 = 16.pi^2.sin^2(theta)/wavelength^2'
          return
      endif
!                      1234567
      if( s(1:7) .eq. 'units=Q' ) then
          units = 'Q'
          write(*,*)'Binning into Q = 4.pi.sin(theta)/wavelength'
          return
      endif

      write(*,*)'Did not understand your option',s
      write(*,*)'Use nothing for twotheta, units=Q or units=Q2'
      stop 
      end subroutine setunits

end module useroptions                                            
      subroutine tidyup
! Free all allocated memory and close all opened files
      use specfiles    ! iunit
      use rebin        ! ascan
      use summation    ! hist, sumdata
      use outputfiles  ! ioutunit
      use useroptions  ! exdetlist, exscanlist
      if(allocated(ascan     ))deallocate(ascan)
      if(allocated(hist      ))deallocate(hist )
      if(allocated(sumdata   ))deallocate(sumdata)
      if(allocated(exdetlist ))deallocate(exdetlist)
      if(allocated(exscanlist))deallocate(exscanlist)
! do an inquire on these to see if they are actually open
      close(iunit)    ! specfile
      close(ioutunit) ! output data
      close(logfile)  ! logfile
      close(16) ! temp.res unit 
      return
      end subroutine tidyup                                                  
      subroutine helpmsg
! Writes a helpful message to stdout
      character(len=80)::name      
      call getarg(0,name)
      write(*,'(a)')name(1:len_trim(name))//' version 10-05-2011'
      write(*,*)
      write(*,1000)name(1:len_trim(name))
                                                                        
1000  format('example: ',a,' file.dat 0.01 1 10' /                    &
     &'will process the scans 1 to 10  from file.dat with binsize 0.01' &
     &//'Optional arguments:'                                           &
     &/'   ed=n1,n2 to exclude detectors n1 and n2 (default=none)'      &
     &/'   es=m1,m2 to exclude scans m1 and m2 (default=none)'          &
     &/'   is=m1,m2 to include only scans m1 and m2 (default=none)'     &
     &/'   lowtth=xx.xx to set min two theta to use (default=-30.0)'    &
     &/'   hightth=xx.xx to set max two theta to use (default=160.0)'   &
     &/'   step=x.xxxx to force a stepsize (if very small)'             &
     &/'   wd=xx to set a limit for outliers on diagnostic file',       &
     & ', diag.mtv'                                                     &
     &/'   wvln=xx.xx to set the wavelength used for unit conversion',  &
     &/'   units=Q bins into Q = 4*pi*sin(theta)/wavelength',           &
     &/'   units=Q2 bins into Q^2 = 16*pi^2*sin^2(theta)/wavelength^2', &
     &/'   alp=xx for esd=sqrt(cts+alp), default is 0.5'                &
     &/'   mm=xx sets the minimum monitor counts threshold (default=5)' &
     &/'   mr=xx for minimum counts needed in all channels to use '     &
     &/'      points for determining detector efficiencies (default=1)' &
     &/'   zap=xx for esd level in filtering operation                ' &
     &/'   superzap=xx for zinger elimination, 0.0 < xx < 1.0         ' &
     &/'          where xx is the fractional reduction is esd for     ' &
     &/'          removing the highest point'                           &
     &/'   medianofchannels to get median_channel*n_active_channels   ' &
     &/'   3pf for a 3 point median filter (you need to be desparate) ' &
     &/'   rst=x.xx to randomly offset scan start points by rnd*x.xx  ' &
     &/'   rst=x.xx,1,2,3 to do that for channels 1,2,3 only          ' &
     &/'   rnd=x.xx to randomly offset scan end points by rnd*x.xx    ' &
     &/'   rnd=x.xx,1,2,3 to do that for channels 1,2,3 only          ' &
     &/'   renorm to use current effics instead of values from temp.res'&
     &/'   nodiag to prevent creation of diagnostic file, diag.mtv'     &
     &/'   gsas to output a .gsa file for gsas'                         &
     &/'   spf to output a .spf file for profil'                        &
     &/'   pds to output a .pds file for a pds file'                    &
     &/'   epf to output a .epf file for a epf file'                    &
     &/'   bcm to output a .bcm file (binned counts,monitor)'           &
     &/'   window=counter,low,high to reject bad points'                &
     &/'       eg: window=lake,4,6 for temperatures between 4 and 6 K'  &
     &/'   snbl  if your data are in the Swiss Norwegian format'        &
     &/'   *** PATTERN SCALING ********************************************'& 
     &/'     scalpk to scale highest peak to correct number of obs counts'  &
     &/'     scaltot to make total number of counts in scan correct   '     &
     &/'     scalmon=xxx for counts per xxx moniter counts (default=100000)'&
     &/'   ****************************************************************') 
     stop
     end subroutine helpmsg                                                  
      program id31sum
      use specfiles
      use rebin
      use useroptions
      use summation
      use outputfiles
      integer(kind=4)::i
      real start_time, end_time
      call cpu_time(start_time)
      isc=1 ! scalpk
      call getcmdline                                 ! Get users options
! fill out names of monitor, tth, first and last columns
      if(snbl)then
       MONITORCOL="Mon";FIRSTDET="Det1";LASTDET="Det6";TWOTTH="TwoTheta"
      else
       FIRSTDET="MA0";LASTDET="MA8";TWOTTH="2_theta"      
      endif     
      call initialiserebin          ! Allocate space and check parameters
      call getfile                                  ! Opens the spec file
      call openlogfile            ! Get the file for recording T etc open
      write(*,'(a,$)')'Processing scan '
1     i=nextscan()                       ! Gets next scan to be processed
      if(i.gt.0) then
        call processscan(i)                       ! sums into ascan array
        if(ispecerr.eq.-1)then
           ispecerr=0
           goto 1                           ! next only if current was OK
        endif   
      endif
      write(*,*)            ! after non advancing io list of scans binned
      call calibsum
      if(bcm)call bcmfile(0)
      !if(diag) call outputdiagnosticw32(wd) ! for windows and prestoplot
      if(diag)call outputdiagnostic(wd)!If problems in combining chans
      hist=hist*scalinp ! apply scale factor after diagnostic plot
      select case(isc)                             ! rescale if necessary
      case(0) ; continue                     ! do nothing
      case(1) ; call scalpk        ! scale to peak height
      case(2) ; call scaltot       ! scale to total counts
      case default ; continue               ! do nothing
      end select
!     call outputxye                           ! Writes out sumdata array
      call outputformats(0)                ! Any other (GSAS etc formats)
      call tidyup                                ! Frees allocated memory
      call cpu_time(end_time)
      write(*,'(a,f6.2,a)')'Time taken was ',end_time-start_time,'/s'
      end program id31sum                                                    