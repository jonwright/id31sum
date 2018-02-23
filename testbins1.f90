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
      subroutine report
      use rebin
      integer(kind=4) :: i, j, k
1000  format('low ',f12.5,' high ',f8.5,' step ',f6.4, ' npts ',i5)
1001  format('bin ', i5,' low ',f12.5,' cen ',f12.5,' high ',f12.5)
      write(*,1000)tthlow,tthhigh,step,npts
      do i=1,4
        j=npts+2-i
        k=i-2
        write(*,1001)j,tthlb(j),bincen(j),tthhb(j)
        write(*,1001)k,tthlb(k),bincen(k),tthhb(k)
      enddo
      i=ibin(0.0d0)
      write(*,1001)i,tthlb(i),bincen(i),tthhb(i)      
      end subroutine report                                                  
      program testbins1 
      use rebin ! pull in the module defined above
      integer(kind=4) :: i, j, n
      real(kind=8) :: tth, steplocal
      character *80 string
      step=0.1d0 ; tthlow=-10.023487d0 ; tthhigh=60.0129384d0 
      call getarg(1,string)
      read(string,*,err=1, end=1)step
      call getarg(2,string)
      read(string,*,err=1, end=1)tthlow
      call getarg(3,string)
      read(string,*,err=1, end=1)tthhigh
1     write(*,'(a)')'Calling initialise rebin'
      call initialiserebin
      call report
      write(*,'(a)')'Calling initialise rebin'
      call initialiserebin
      call report
! Now for some numbers to bin
      write(*,'(a)')'Two Theta, IB, LOW, HIGH, CENTRE'
      steplocal=0.0242349876
      n=nint((tthhigh-tthlow)/steplocal)-1   ! random hard to bin numbers?
      do i=1,n
        tth=(tthhigh-tthlow)*float(i)/float(n) + tthlow
        j=ibin(tth)
        write(*,'(f15.9,I9,3f12.4)')tth,j,tthlb(j),tthhb(j),bincen(j)
      enddo
      end program testbins1                                                  