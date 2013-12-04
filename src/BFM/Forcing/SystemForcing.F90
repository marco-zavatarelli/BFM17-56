#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bfm
!
! !INTERFACE:
   module SystemForcing
!
! !DESCRIPTION: 
! System Forcing I/O definition for the BFM. 
!
! !USE:
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY:NO_BOXES_XY
#endif

   use global_mem,  only:RLEN,ZERO,bfm_lwp,LOGUNIT
   use constants, only: SEC_PER_DAY
   use api_bfm, ONLY: SRFindices
   use time
   use netcdf

   implicit none
   PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
  
! Possible convention to differentiate input field:
! 0 = No input / Set Model constant value
! 1 = Analytic
! 2 = 1D timeseries
! 3 = 2D Field
! 4 = Coupling with external model

   type, public :: ForcingName
      integer           :: init
      character(LEN=80) :: filename
      logical           :: filetype
      character(LEN=40) :: varname
      character(LEN=40) :: RefTime
      character(LEN=40) :: cltype
      logical           :: tinterp 
   end type ForcingName
!
   type, public :: ForcingField
      integer           :: init
      integer           :: lun
      logical           :: filetype
      integer           :: varID
      integer           :: nrec 
      character(LEN=40) :: cltype
      logical           :: tinterp 
      real(RLEN)                          :: tnow
      real(RLEN),allocatable,dimension(:) :: fnow
      real(RLEN)                          :: tbef
      real(RLEN),allocatable,dimension(:) :: fbef
      integer                             :: nbef
      real(RLEN)                          :: taft
      real(RLEN),allocatable,dimension(:) :: faft
      integer                             :: naft
   end type ForcingField

!
! !PUBLIC DATA MEMBERS:
   public FieldInit, FieldRead, FieldClose
!
! !REVISION HISTORY:
!  Author(s): Marcello 
!
! ! PRIVATE MEMBERS:
!

   contains
!-------------------------------------------------------------------------!
!BOC

   subroutine FieldInit(FName, FData)

   use netcdf_bfm, only: check_err

   integer :: yy, mm, dd, hh, nn, jday, LUN, icon
   integer :: iy, im, id, iyend, imend, idend
   integer :: jh, jn 
   integer :: RecordDimID, nRecords
   character(len = nf90_max_name) :: RecordDimName
   character(len = 50) :: InpDate, c1
   logical :: timedone
   real(RLEN) :: jday0, jdayaft, jdaybef, ValAft, ValBef
   type(ForcingName),  intent(IN)     :: FName
   type(ForcingField), intent(INOUT)  :: FData
   
   iyend=0
   imend=0
   idend=0
   timedone=.FALSE. 
   icon = 0
   hh = 0
   nn = 0
   ! additional check
   if (FName%init == 0 .OR. FName%init == 4 ) then 
      allocate (FData%fnow(NO_BOXES_XY))
      write(LOGUNIT,*) 'FieldInit Warning ',trim(FName%varname), &
                       ': Data will not be read, ONLY %fnow memory structure is allocated !'
      return
   endif
   if (FName%init < 0 .OR. FName%init > 4 ) then
      LEVEL1  'FieldInit Error ',trim(FName%varname), &
                       ': initialization flag %init is  invalid! Allowed 0-4.'
      STOP
   endif 
   ! Temporary check to avoid unavailable data structure selection 
   if (FName%init == 1 .OR. FName%init == 3 ) then
      LEVEL1  'FieldInit Error ',trim(FName%varname), &
                       ': initialization flag %init with option 1 and 3 is not yet implemented in the code.'
      STOP
   endif
   ! Initialize FData structure 
   allocate (FData%fbef(NO_BOXES_XY), FData%fnow(NO_BOXES_XY), FData%faft(NO_BOXES_XY))
   FData%init     = Fname%init 
   FData%filetype = Fname%filetype
   FData%cltype   = Fname%cltype
   FData%tbef     = bfmtime%time0      
   FData%tnow     = bfmtime%time0      
   FData%taft     = bfmtime%time0      
   FData%tinterp  = FName%tinterp 
   FData%nbef     = 0
   FData%naft     = 0

   read (FName%RefTime,'(I4,a1,I2,a1,I2,1x,I2,a1,I2)',ERR=903) yy,c1,mm,c1,dd,hh,c1,nn
   call julian_day(yy,mm,dd,hh,nn,jday0)

   ! ACCESS EXTERNAL INPUT & FIND TIMELINE BOUNDARIES
   IF (FName%filetype) THEN 
      ! NETCDF FILE 
      call check_err(NF90_OPEN(trim(FName%filename),NF90_NOWRITE,FData%lun), FName%filename)
      ! Get unlimited dimension name and length (e.g. time and # of records)
      call check_err(NF90_INQUIRE(FData%lun , unlimitedDimId = RecordDimID), FName%filename) 
      call check_err(NF90_INQUIRE_DIMENSION(FData%lun, RecordDimID, name = RecordDimName, len = nRecords), FName%filename) 
      FData%nrec = nRecords
      ! Get variable ID
      call check_err(nf90_inq_varid(FData%lun, trim(FName%varname), FData%varID), FName%filename)
      ! Setup centered time for inputs
      select case (FData%cltype)
         case('yearly')
           iyend = FData%nrec - 1
           mm = 0
           call halftime(yy,mm,dd,hh) 
         case('monthly')
           iyend = (FData%nrec / 12 )
           imend = 1
           call halftime(yy,mm,dd,hh)
         case('daily')
           iyend = (FData%nrec / 365 ) + 1
           imend = 1
           idend = 1
           hh = 12 
         case default 
           write(LOGUNIT,*) 'FieldInit: Unrecognized time format for file : ', FName%filename
           write(LOGUNIT,*) 'FieldInit supports: yearly , monthly , daily'
           stop 
      end select

      ! Find the before and after input steps 
      jdaybef = bfmtime%time0
      do iy  = yy , yy + iyend
        if ( imend .NE. 0 ) imend = 12 - mm 
        do im = mm , mm + imend 
           if ( idend .NE. 0 ) idend = eomdays(iy,im) - dd 
           do id = dd , dd + idend
               icon = icon + 1 
               call julian_day(iy,im,id,hh,nn,jdayaft)
               if ( jdayaft > bfmtime%time0 .and. .not. timedone) then
                  timedone=.TRUE.
                  FData%taft = jdayaft
                  FData%naft = icon   
                  FData%nbef = icon-1   
                  FData%tbef = jdaybef
               endif
               jdaybef = jdayaft
           enddo 
           if (FData%cltype == 'daily') dd = 1
        enddo
        if (FData%cltype == 'monthly') mm = 1
      enddo
      ! Set uniform value : Backward in time 
      if (FData%tbef .EQ. bfmtime%time0 .AND. jday0 > bfmtime%time0) then
         write(LOGUNIT,*) 'FieldInit: Backward use in time of the last input value as a constant.'
         FData%nbef = FData%naft
      endif
      ! forward in time
      if (.not. timedone) then
         write(LOGUNIT,*) 'FieldInit: Forward use in time of the last input value as a constant.'
         FData%tbef = jdaybef
         FData%nbef = icon         
         FData%taft = bfmtime%timeEnd
         FData%naft = icon
      endif
      ! Load initial data
      call FieldGet(FData)
   ELSE 
   ! SEQUENTIAL FILE  (Text file with : time , data)
   ! it works for timeseries and 1D data in the form time, value 1...#total
      LUN=GetLun()
      FData%lun = LUN
      ! Open file and read Header
      open(FData%lun,file=trim(FName%filename),ERR=902) 
      read(FData%lun,*) InpDate
      jdayaft = 0
      jdaybef = bfmtime%time0
      do 
         ValBef  = ValAft
         read(FData%lun,*,END=901,ERR=902) InpDate,ValAft
         write(*,*) InpDate,ValAft
         read (InpDate,'(I4,a1,I2,a1,I2,1x,I2,a1,I2)') yy,c1,mm,c1,dd,hh,c1,nn
         call julian_day(yy,mm,dd,hh,nn,jdayaft)
         if ( jdayaft > bfmtime%time0 ) then
            timedone=.TRUE.
            FData%taft = jdayaft
            FData%tbef = jdaybef
            FData%faft = ValAft
            FData%fbef = ValBef
            EXIT
         endif 
         jdaybef = jdayaft
         ValBef = ValAft
      enddo
      ! Set uniform value : Backward in time 
      if (FData%tbef .EQ. bfmtime%time0 .AND. jday0 > bfmtime%time0 ) then
         write(LOGUNIT,*) 'FieldInit: Backward use in time of the last input value as a constant.'
         FData%fbef = FData%faft
      endif 
      ! forward in time
901   if (.not. timedone) then 
         write(LOGUNIT,*) 'FieldInit: Forward use in time of the last input value as a constant.'
         FData%fbef = ValBef
         FData%faft = ValBef
         FData%taft = bfmtime%timeEnd 
      endif
   ENDIF

   call FieldRead(FData)

   return

902 LEVEL1 'FieldInit: Error opening sequential input file:', FName%varname
   write(LOGUNIT,*) 'FieldInit: Error opening sequential input file:', FName%varname
   stop 

903 LEVEL1 'FieldInit: Error reading namelist RefTime for variable: ', FName%varname
   write(LOGUNIT,*) 'FieldInit: Error reading namelist RefTime for variable: ', FName%varname
   stop
   end subroutine FieldInit
!-------------------------------------------------------------------------!
!-------------------------------------------------------------------------!
   subroutine FieldRead(FData)

   use netcdf_bfm, only:check_err

   integer        :: yy, mm, dd, hh, nn
   real(RLEN)     :: ValAft, jday
   real(RLEN),allocatable,dimension(:) :: diff
   logical        :: readok
   character(len = 50) :: InpDate, c1
   type(ForcingField), intent(INOUT)   :: FData
     
   ! additional check 
   if (FData%init == 0 .OR. FData%init == 4 ) then
      if (bfmtime%stepnow ==  bfmtime%step0) &
         write(LOGUNIT,*) 'FieldRead Warning: Data will not be read because filed %init is either 0 or 4!'
      return
   endif
   ! Temporary check to avoid unavailable data structure selection 
   if (FData%init == 1 .OR. FData%init == 3 ) then
      LEVEL1  'FieldRead Error : initialization flag %init ', &
              'with option 1 and 3 is not yet implemented in the code.'
      STOP
   endif

   allocate(diff((NO_BOXES_XY)))
 
   ! Set actual time 
   FData%tnow =( ( (bfmtime%stepnow - bfmtime%step0) * bfmtime%timestep ) &
                 / SEC_PER_DAY ) + bfmtime%time0 

   ! Read in Forcing file
   IF ( (FData%tnow > FData%taft) )  THEN 
      ! NETCDF 
      IF (FData%filetype) THEN 
         if (FData%naft == FData%nrec) then
            write(LOGUNIT,*) 'FieldRead: Forward use in time of the last input value as a constant.'
            FData%tbef = FData%taft
            FData%nbef = FData%naft
            FData%taft = bfmtime%timeEnd
            FData%fbef = FData%faft
         else
            ! Update records to read
            FData%nbef = FData%naft
            FData%naft = FData%naft + 1
            call FieldGet(FData)
            ! update time 
            FData%tbef = FData%taft
            call calendar_date(FData%taft,yy,mm,dd,hh,nn)
            select case (FData%cltype)
               case('yearly') 
                  yy = yy + 1
               case('monthly') 
                  mm = mm + 1
                  if (mm == 13) then
                     mm = 1
                     yy = yy + 1
                  endif
               case('daily') 
                  dd = dd + 1 
                  if (dd > eomdays(yy,mm)) then
                     dd = 1 
                     mm = mm + 1
                     if (mm == 13) then
                        mm = 1
                        yy = yy + 1
                     endif
                  endif
            end select
            call julian_day(yy,mm,dd,hh,nn,FData%taft)
         endif
      ELSE
      ! SEQUENTIAL
         readok = .FALSE. 
         select case (FData%init)
            case(2) ! Timeseries
               read(FData%lun,*,END=905,ERR=906) Inpdate, ValAft
               readok = .TRUE.
               read (InpDate,'(I4,a1,I2,a1,I2,1x,I2,a1,I2)') yy,c1,mm,c1,dd,hh,c1,nn
               call julian_day(yy,mm,dd,hh,nn,jday)
               FData%tbef = FData%taft 
               FData%taft = jday
               FData%fbef = FData%faft 
               FData%faft = ValAft 
905            if (.not. readok) then
                  write(LOGUNIT,*) 'FieldRead: Forward use in time of the last input value as a constant.'
                  FData%tbef = FData%taft      
                  FData%fbef = FData%faft
                  FData%taft = bfmtime%timeEnd
               endif
            case(0,1,3,4) 
               Stop ' FieldRead: Wrong type of input. Only timeseries for sequential file.'
         end select
      ENDIF
   ENDIF

   ! Linear time interpolation
   if (FData%tinterp) then
      diff = ( FData%faft - FData%fbef ) / ( FData%taft - FData%tbef )
      FData%fnow = FData%fbef + (FData%tnow - FData%tbef) * diff
   else
      FData%fnow = FData%fbef
   endif   

   deallocate(diff)
   return

906 write(LOGUNIT,*) ' FieldRead: Wrong type of input. Only timeseries for sequential file.'
    stop
   end subroutine FieldRead
!-------------------------------------------------------------------------!
!-------------------------------------------------------------------------!
   subroutine FieldGet(FData)

   use netcdf_bfm, only:check_err

   real(RLEN)     :: Val0D
   type(ForcingField), intent(INOUT)   :: FData

   ! read in records
      select case (FData%init)
         case(1) ! Analytic

         case(2) ! Timeseries
            call check_err( nf90_get_var(FData%lun, FData%varID, Val0D, start = (/FData%nbef/)), 'FieldGet') 
            FData%fbef = Val0D
            call check_err( nf90_get_var(FData%lun, FData%varID, Val0D, start = (/FData%naft/)), 'FieldGet')
            FData%faft = Val0D 
         case(3) ! 1D fields

         case(4) ! 2D fields

      end select
      return

   end subroutine FieldGet
!-------------------------------------------------------------------------!
!-------------------------------------------------------------------------!
 subroutine FieldClose(FName, FData)
   ! This routine close the forcing file

   use netcdf_bfm, only: check_err

   implicit none
   type(ForcingName),  intent(IN) :: FName
   type(ForcingField), intent(IN) :: FData

   IF (FData%filetype) THEN 
      call check_err(NF90_CLOSE(FData%lun))
   ELSE
      close(FData%lun)
   ENDIF

   return
 end subroutine 
!-------------------------------------------------------------------------!
!-------------------------------------------------------------------------!
 integer function GetLun ()
! adapted from Lionel, Shepherd, Clodius, Page, Drummond.
! and others as posed at comp.lang.fortran on 1997-09-01
      implicit none
      logical :: exs, opn
      integer :: i
      getlun = -1  ! returned if no units are available.
      i = 1000
 L1:  do 
        inquire (unit=i,exist=exs,opened=opn)
          if (exs .and. .not. opn) then
            getlun = i
            exit L1
          end if
        i = i + 1
      end do L1
      return
      stop "There are no free Fortran logical units available."
 end function getlun
!-------------------------------------------------------------------------!
!-------------------------------------------------------------------------!
 subroutine halftime(Year, Month, Day, Hour)
 ! This routine compute the month or year central date in mm-dd-hh
     implicit none
     integer, intent(IN):: Year
     integer, intent(INOUT):: Month, Day, Hour
     real(RLEN)  :: ref, res, now, good
     if (Month == 0 ) then  ! year case
        ref = FLOAT(yeardays(Year)) / 2
  L1:   do 
          Month = Month + 1 
          now = now + FLOAT(eomdays(Year,Month))        
          res = ref - now
          if (res < 0)  exit L1
        enddo L1
        good = res + FLOAT(eomdays(Year,Month-1))
     else ! Month case
        good = FLOAT(eomdays(Year,Month)) / 2
     endif 
     Day = FLOOR(good) 
     Hour = (good - FLOAT(Day)) * 24
     ! correct Day value when Floor of # < 1 gives 0 (like 0.5)
     if (Day == 0) Day = 1
     return
 end subroutine 
!-------------------------------------------------------------------------!
!EOC
!-------------------------------------------------------------------------!
! BFM Module
!-------------------------------------------------------------------------!
   end module SystemForcing  
