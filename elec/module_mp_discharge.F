!
#define FLOAT_PRECISION MPI_DOUBLE_PRECISION

      Module  module_mp_discharge 
#ifdef DM_PARALLEL
      USE module_dm, ONLY : &
         local_communicator, mytask, ntasks, ntasks_x, ntasks_y,mytask_x, mytask_y
      USE module_commasmpi
#endif

      CONTAINS

      subroutine light1d(iter,loccur,   &
     &   nor,nx,ny,nz      &
     &  ,dxx,dyy,lightrad,ibrkd,elght1,elgt1,dzz,sctot &
     &  ,scw,scs,sci,scr,sch,schl &
     &  ,qc,qs,qi,qr,qg,qh &
     &  ,qnc,qns,qni,qnr,qng,qnh &
     &  ,elec,scion,sciona        &
     &  ,nbw,nbe,nbs,nbn   &
     &  ,iunit, db,light,lightdens,etot,lightdis)
     
      implicit none

#ifdef DM_PARALLEL
      INCLUDE 'mpif.h'
#endif
      
      integer iunit,tag2,iter,loccur
      real elec(0:nx+1,nz,0:ny+1) 
      real scion(0:nx+1,nz,0:ny+1,2)  ! 2=+, 1=- 
      real dzz(0:nx+1,nz,0:ny+1),dxx,dyy   
      real db(0:nx+1,nz,0:ny+1)
      integer, parameter    :: ng1 = 1
      integer nbw,nbe,nbs,nbn
      
      double precision chg,chgfracn,chgfracp,chgneg,chgpos ! total net negative/positive charge in storm (not including corona)
      real  :: chgmax
      
      integer ix,jy,kz,kz1,i,j
      integer nor,nx,ny,nz
      real pi
      integer ietotx,ietoty,ietotz

!  ibrkd:    ibrkd is they of break down; 1=based on emax; 2=based
!            on emax(z)

      real  elght1,lightrad,elgt1

      integer nlgt1,ibrkd
!      integer :: ibrkd    = 2
!      integer :: lightrad = 6000.0
      real    :: chgthr   = 0.1e-9  ! threshold charge density
!      real    :: elght1   = 120000. ! threshold e-field to init lightning
!      real    :: elgt1    = 0.3     ! fractional reduction in charge in excess of chgthr
      real    :: eperao=8.8592e-12
      real    :: ec = 1.602e-19 ! fundamental unit of charge 
      real    :: rho00 = 1.225
      real    :: qmin = 1.e-12 !1.e-5
!
      real    :: rhos = 100.
      real    :: rhog = 400.
      real    :: rhoh = 900.
      real    :: rhor = 1000.
      real    :: rhoi = 900.
      real    :: rhow = 1000.

      real elgt,temp
      double precision etot
      real sctot(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real sciona(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real sch(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real,optional :: schl(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real scw(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real scr(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real sci(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real scs(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)

      real qg(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real,optional :: qh(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real qr(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real qc(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real qi(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real qs(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      
      real qng(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real,optional :: qnh(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real qnr(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real qnc(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real qni(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real qns(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)

      real frachl
      real frach
      real fracr
      real fracw
      real fraci
      real fracs
      real fractot 
      double precision dv,deninv, tmp

      real :: light(-nbw:nx+1+nbe,-nbs:ny+1+nbn)
      real :: lightdens(-nbw:nx+1+nbe,-nbs:ny+1+nbn)
      integer :: lightdis(-nbw:nx+1+nbe,-nbs:ny+1+nbn)
         
       real,allocatable,save :: ebrkdp(:,:,:),ebrkd(:,:,:),zlev(:,:,:)
 
      real, allocatable :: scdisch(:,:,:) 
      real,   allocatable  :: xy_init(:,:)
      integer, allocatable :: ij_reduce(:,:)
      integer, allocatable,save  :: ijmask(:,:)
      integer, save :: nxm, nym
      integer  :: nr, ninit
      real  :: x,y,r


!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES

      integer, parameter :: ntot = 50
      real  mpitotin(ntot), mpitotout(ntot)
      double precision  mpitotindp(ntot), mpitotoutdp(ntot)
      real totneg1,totpos1,totneg2,totpos2,totneg3,totpos3,chglost
      real totnegbef,totnegaft,totposbef,totposaft

      integer :: westward_tag, eastward_tag
      integer :: northward_tag, southward_tag
      INTEGER :: mpi_error_code = 1      
      INTEGER,DIMENSION(:),ALLOCATABLE :: w_proc, e_proc, s_proc, n_proc
      INTEGER :: memory_status,process
      integer :: my_rank = 0
! #####################################################################
!   Begin Executable Code
! #####################################################################
       ALLOCATE(w_proc(0 : ntasks - 1),STAT=memory_status)
       ALLOCATE(e_proc(0 : ntasks - 1),STAT=memory_status)
       ALLOCATE(s_proc(0 : ntasks - 1),STAT=memory_status)
       ALLOCATE(n_proc(0 : ntasks - 1),STAT=memory_status)

      allocate( scdisch(-nbw:nx+1+nbe,-nbs:ny+1+nbn,nz) )

       scdisch(:,:,:)=0. 
       nlgt1 =0
       etot =0.
       loccur = 0

!       write(0,*) ,'DIS: elght1,elgt1,lightrad=',elght1,elgt1,lightrad

#ifdef DM_PARALLEL
 
              w_proc(mytask) = MPI_UNDEFINED 
              e_proc(mytask) = MPI_UNDEFINED
              s_proc(mytask) = MPI_UNDEFINED
              n_proc(mytask) = MPI_UNDEFINED

              my_rank = mytask

#endif

              if (mytask_x.gt.0.and.mytask_x.lt.ntasks_x-1) THEN   ! in the "middles"

              e_proc(mytask) = mytask+1
              w_proc(mytask) = mytask-1

              ENDIF

              if (mytask_y.gt.0.and.mytask_y.lt.ntasks_y-1) THEN   ! in the "middles"

              s_proc(mytask) = mytask-ntasks_x
              n_proc(mytask) = mytask+ntasks_x

              ENDIF

              if (mytask_x.eq.0) THEN   ! W boundary
      
              e_proc(mytask) = mytask+1

              ENDIF

              if (mytask_x.eq.ntasks_x-1) THEN ! E boundary
              
              w_proc(mytask) = mytask-1
 
              ENDIF

              if (mytask_y.eq.0) THEN  ! S boundary

              n_proc(mytask) = mytask+ntasks_x

              ENDIF

              if (mytask_y.eq.ntasks_y-1) THEN  ! N boundary

              s_proc(mytask) = mytask-ntasks_x 
     
              ENDIF

!        write(0,*) '+============BEFORE DISCHARGE============+'
!        write(0,*) 'SCTOT MIN', MINVAL(sctot(:,:,:))
!        write(0,*) 'SCTOT MAX', MAXVAL(sctot(:,:,:))
!        write(0,*) 'ELECMAG MAX',MAXVAL(elec(1:nx,1:nz,1:ny)),MAXLOC(elec(1:nx,1:nz,1:ny))
!        write(0,*) 'SCW MIN', MINVAL(scw(:,:,:))
!        write(0,*) 'SCW MAX', MAXVAL(scw(:,:,:))
!        write(0,*) 'SCR MIN', MINVAL(scr(:,:,:))
!        write(0,*) 'SCR MAX', MAXVAL(scr(:,:,:))
!        write(0,*) 'SCH MIN', MINVAL(sch(:,:,:))
!        write(0,*) 'SCH MAX', MAXVAL(sch(:,:,:))
!        write(0,*) 'SCI MIN', MINVAL(sci(:,:,:))
!        write(0,*) 'SCI MAX', MAXVAL(sci(:,:,:))
!        write(0,*) 'SCS MIN', MINVAL(scs(:,:,:))
!        write(0,*) 'SCS MAX', MAXVAL(scs(:,:,:))

!
!  compute ebrkz (breakdown field as a function of z)
!  function computed by Straka using Marshall and Rust plot
!  provided by MacGorman
!
!       zlev=0. ; ebrkd= 1.e8 ; ebrkdp= 1.e8 ! dont want zero ecrit otherwise light might be out of control
       fracw = 0. ; fracs = 0. ; frachl=0. ; fraci =0. ; frach = 0. ; fracr = 0.       

      IF (.not.allocated(ijmask)) THEN
       
!      create disks mask for cylinders 

      nxm = Int(abs(lightrad)/dxx) + 1
      nym = nxm

      allocate( ijmask(-nxm:nxm,-nym:nym) )

      ijmask(:,:) = 0

       IF ( lightrad <= 0.0 ) THEN
         ijmask(:,:) = 1
       ELSE
       DO j = -nym,nym
        DO i = -nxm,nxm
!   assume dx=dy here
          r = dxx*Sqrt( Float(i**2 + j**2) )
          IF ( r <= Abs(lightrad) ) ijmask(i,j) = 1
        ENDDO
       ENDDO
       ENDIF

       ENDIF

      IF ( .not. allocated( ebrkd ) ) THEN ! {
       allocate( ebrkd(-nbw:nx+1+nbe,-nbs:ny+1+nbn,nz) )
       allocate( ebrkdp(-nbw:nx+1+nbe,-nbs:ny+1+nbn,nz) )
       allocate( zlev(nz,-nbw:nx+1+nbe,-nbs:ny+1+nbn) )

        ebrkd(:,:,:) = 1.e8
        ebrkdp(:,:,:) = 1.e8

       IF (ibrkd .le. 1) ebrkd(:,:,:) = elght1

       IF (ibrkd .ge. 2) THEN

      DO ix=-nbw, nx+1+nbe
      DO jy=-nbs, ny+1+nbn
      do kz = 1,nz
      
      zlev(kz,ix,jy) = sum(dzz(ix,1:kz,jy)) - 0.5*dzz(ix,kz,jy) ! dz8w is actually for w points?
!      ebrkd(ix,kz,jy) = 197.81
!     >          - 1.9606e-2*zlev(ix,kz,jy)
!     >          + 5.7664e-7*zlev(ix,kz,jy)**2
      ebrkd(ix,jy,kz) = 1.208*167.0*Exp(-zlev(kz,ix,jy)/8.4e3)
      ebrkd(ix,jy,kz) = ebrkd(ix,jy,kz) * 1.e3
      kz1 = Min(kz,nz-1)
      ebrkdp(ix,jy,kz) = 284.e3*db(ix,kz,jy)/rho00

      IF (ibrkd .eq. 2) THEN
        ebrkd(ix,jy,kz)= Min( ebrkd(ix,jy,kz), 125.0e3 )
      ELSEIF ( ibrkd .eq. 3 ) THEN
        ebrkd(ix,jy,kz)= Min( ebrkd(ix,jy,kz), 100.0e3 )
      ELSEIF ( ibrkd .eq. 4 ) THEN
        ebrkd(ix,jy,kz)= Max( 50.e3, Min( ebrkdp(ix,jy,kz), 180.0e3 ) )
!        ebrkd(ix,jy,kz)= Max( 65.e3, Min( ebrkdp(ix,jy,kz), 180.0e3 ) )
      ELSEIF ( ibrkd .eq. 5 ) THEN
        ebrkd(ix,jy,kz)= Min( ebrkdp(ix,jy,kz), 100.0e3 )
      ENDIF
    
!      if (ix.eq.1.and.jy.eq.1) write(0,*) 'erbd', ebrkd(ix,jy,kz),kz,zlev(ix,jy,kz),dzz(ix,kz,jy)
 
!      ecrit(ix,kz,jy)=ebrkd(ix,jy,kz)

 
      IF ( ebrkd(ix,jy,kz) .ne. 0.0 ) THEN
        dv = ebrkdp(ix,jy,kz)/ebrkd(ix,jy,kz)
      ELSE
        dv = 0.0
      ENDIF
      
      end do
      end do
      end do
      
      ENDIF ! }
!
       END IF 

       
      pi = 4.*atan(1.)
!
!  This scheme is the simplist of all.  If |E| => elgtt1 
!  (400 kV/m; REF) then the space charge is decreased by some percent
!  elgtp1 (30%:  Rawlins, 1981, QJRMS) and Ziegler and MacGorman (1994,JAS).
!  

      IF ( lightrad < 0 ) THEN 
        nr = 0
      ELSE
        nr = NInt( lightrad/dxx )
      ENDIF
      
      allocate( xy_init(-nr+1:nx+nr,-nr+1:ny+nr))
      allocate( ij_reduce(-nr+1:nx+nr,-nr+1:ny+nr))
      
      xy_init(:,:) = 0.0 ! zero out in order not to add lightning...even though these are deallocated.

      IF ( lightrad < 0 ) THEN 
        ij_reduce(:,:) = 1
      ELSE
        ij_reduce(:,:) = 0
      ENDIF
!
!  find largest etot and location, but first zero etot and location
!  indices first
!
      etot = 0.0
      ietotx = 0
      ietoty = 0
      ietotz = 0
      ninit = 0
!
      DO ix=-nbw, nx+1+nbe
      DO jy=-nbs, ny+1+nbn
      do kz = 1,nz
!
      if ( elec(ix,kz,jy) .gt. etot ) then
      ietotx = ix
      ietoty = jy
      ietotz = kz
      etot = elec(ix,kz,jy)
      end if
     
!      write(*,*) 'elec(ix,kz,jy),ebrkd(ix,jy,kz)',elec(ix,kz,jy),ebrkd(ix,jy,kz)

!      if ( elec(ix,kz,jy) .gt. elght1 ) then
!      if ( elec(ix,kz,jy,4) .gt. ebrkd(ix,jy,kz)   ) then
      if ( elec(ix,kz,jy) .gt. ebrkd(ix,jy,kz) ) then
        xy_init(ix,jy) = 1.0
        ninit = ninit + 1
      ENDIF
   
!        light(ix,jy) = light(ix,jy) + xy_init(ix,jy) ! set to zero ar every iterations..added up outaide in driver
!
      end do
      end do
      end do
      
!      write(0,*) 'ninit = ',ninit
      
#ifdef DM_PARALLEL
! find global integrated rate max
       mpitotindp(1)  = etot
       mpitotindp(2)  = ninit

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_MAX, local_communicator, mpi_error_code)

       etot = mpitotoutdp(1)
       ninit = Nint( mpitotoutdp(2) )

       

      IF ( ntasks .gt. 1 .and. nr > 0 .and. ninit > 0 ) THEN ! because cylinder might encompass several tiles..

        IF ( ntasks_x > 1 ) THEN
!        write(0,*) 'send_west, ntask_x = ',ntasks_x
        westward_tag = 201
        CALL sendrecv_westward(nx,ny,1,nr,nr,0,nr,1,  &
     &        w_proc(mytask),e_proc(mytask),westward_tag,xy_init)

        eastward_tag = 202
        CALL sendrecv_eastward(nx,ny,1,nr,nr,0,nr,1,  &
     &        w_proc(mytask),e_proc(mytask),eastward_tag,xy_init)

        ENDIF
        
        southward_tag = 203
        
        IF ( ntasks_y > 1 ) THEN
!        write(0,*) 'send_north, ntask_y = ',ntasks_y
        CALL sendrecv_southward(nx,ny,1,nr,nr,0,nr,1,  &
     &        n_proc(mytask),s_proc(mytask),southward_tag,xy_init)

        northward_tag = 204
        CALL sendrecv_northward(nx,ny,1,nr,nr,0,nr,1,  &
     &        n_proc(mytask),s_proc(mytask),northward_tag,xy_init)
        
        ENDIF
      ENDIF
#endif

       IF ( my_rank == 0 ) write(0,*) 'ETOT,ninit', ETOT,ninit
       
       IF ( ninit > 0 ) loccur = 1

!      IF ( mytask == 0 ) THEN
!        write(iunit,90) 'LOCAL IETOTX,IETOTY,IETOTZ=: ',ietotx,ietoty,ietotz
!   90  format(A,3(i4,1x))
!      ENDIF
!
!  decrease space charge by percentage elgtp1
!
      chgpos = 0
      chgneg = 0

      IF ( ninit > 0 ) THEN ! BIG LOOP ALL THE WAY TO THE END
!       IF ( etot .gt. elght1 ) THEN ! BIG LOOP ALL THE WAY TO THE END
!      IF ( etot .gt. ebrkd(ix,jy,kz)  ) THEN ! BIG LOOP ALL THE WAY TO THE END

      nlgt1 = nlgt1 + 1

      IF ( nr > 0 ) THEN ! create square region around each init point
       DO jy = -nr+1,ny+nr
        DO ix = -nr+1,nx+nr
          IF ( xy_init(ix,jy) > 0.5 ) THEN ! means that if elec(i,k,j) > elght1
            DO j = Max(-nr+1,jy-nr),Min(ny+nr, jy+nr)
             DO i = Max(-nr+1,ix-nr),Min(nx+nr, ix+nr)
! cuboids  
!               ij_reduce(i,j) = 1
! cylinders     oiverlapping ones are taken into account
               ij_reduce(i,j) = Max( ij_reduce(i,j) , ijmask( i-ix, j-jy))
             ENDDO
            ENDDO
          ENDIF
        ENDDO
       ENDDO
      ELSE
        ij_reduce(:,:) = 1 ! nr < 0 never attained - just used for whole domain
      ENDIF

!  space charge reduction routine      
!  finds total pos/neg space charge within all cylinders/cuboids for
!  sctot >  chgthr
!  compute fraction of +/- charge to be removed upon discharge based on
!  magnitude of the latter. Also decrease of charge on each hydrometeor
!  specie
!  is proportional to their fraction relative to the total surface area.
!  update of sctot for recomputation of E.

       totneg1 = 0.0
       totpos1 = 0.0
      
      DO jy = 1,ny
      DO ix = 1,nx


      IF ( ij_reduce(ix,jy) > 0 ) THEN !{

! check if any point in the column actually has charge density > chgthr

      chgmax = 0.0

      DO kz = 1,nz-1
       chgmax = Max( chgmax, Abs(sctot(ix,kz,jy) ) )
      ENDDO
      
        IF ( chgmax <= chgthr ) THEN
! then no point has chgthr in this column, turn it off as a discharge column
          ij_reduce(ix,jy) = 0

        ELSE
! the column does have chgthr, so count up the available charge

        
        do kz = 1,nz-1
      
         dv = dxx*dyy*dzz(ix,kz,jy)
!
! t1 is the non-ion charge
        IF ( sctot(ix,kz,jy) .gt. chgthr ) THEN
          chgpos = chgpos + (sctot(ix,kz,jy) - chgthr)*dv 
        ELSEIF ( sctot(ix,kz,jy) .lt. -chgthr ) THEN
          chgneg = chgneg + (sctot(ix,kz,jy) + chgthr)*dv
        ENDIF
        
         totneg1 = totneg1 + Min(0.0, sctot(ix,kz,jy) )*dv
         totpos1 = totpos1 + Max(0.0, sctot(ix,kz,jy) )*dv
      
       end do
      
      ENDIF ! chgmax
      ENDIF ! } ( ij_reduce(ix,jy) > 0 ) 
      
      ENDDO
      ENDDO

!   lightning count here per tile - Note that light is real here.


      DO jy = 1,ny
      DO ix = 1,nx

!      light(ix,jy)=FLOAT(ij_reduce(ix,jy)) / ( pi*(lightrad*1.e-03)**2 )
       lightdens(ix,jy)=lightdens(ix,jy)+FLOAT(ij_reduce(ix,jy))
       light(ix,jy)=light(ix,jy)+xy_init(ix,jy)

      ENDDO
      ENDDO


! find global integrated rate sum within cuboids
       mpitotindp(1)  = chgpos
       mpitotindp(2)  = chgneg
       mpitotindp(3)  = totpos1
       mpitotindp(4)  = totneg1

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 4, MPI_DOUBLE_PRECISION, MPI_SUM, local_communicator, mpi_error_code)

       
       chgpos  = mpitotoutdp(1)  
       chgneg  = mpitotoutdp(2)  
       totpos1 = mpitotoutdp(3)  
       totneg1 = mpitotoutdp(4)  

!       write(0,*) 'CHG SUM' , chgpos,chgneg
!#endif

!     otherwise X=Y=0 --> disch=0. so if enough + charge alone for elecmag > ecrit will loop infinitely because scdich =0 and no changes will occur on sctot and elec.
      if (chgpos.lt.1.) chgpos=1.
      if (chgneg.gt.-1.) chgneg=-1.
     

!      what this does: if amount of +/- charge are comparable, then just take 30% of max . IF much more of one polarity, will take just the min of (chgpos, Abs(chgneg))


      IF ( .false. ) THEN
      ! original scheme to set reduction fraction
      chg = elgt1*Max(chgpos, Abs(chgneg))
!      chg = (1.0 - elgt1)*Max(chgpos, Abs(chgneg))
      IF ( chg .gt. chgpos ) chg = chgpos
      IF ( chg .gt. Abs(chgneg) ) chg = Abs(chgneg)

      if (abs(chgneg).gt.0.) then
!      x = chg/(Abs(chgneg)) 
      x = (Abs(chgneg)-chg)/(Abs(chgneg)) 
      else
      x=0.
      endif

      if (abs(chgpos).gt.0.) then 
!      y = chg/(Abs(chgpos))
      y = (chgpos-chg)/(Abs(chgpos))
      else
      y=0.
      endif
     
     ELSE
     ! newer version that tries to offset net charge within the discharge regions
      chg = elgt1*Max(totpos1, Abs(totneg1))
!      chgfracn = elgt1*Abs(totneg1)
!      chgfracp = elgt1*totpos1
! try using average:
      chgfracn = elgt1*0.5*(Abs(chgneg) + Abs(totneg1))
      chgfracp = elgt1*0.5*(chgpos+totpos1)
      
!      IF ( Sign(1.0,totneg1+totpos1) /= Sign(1.0,chgneg + totneg1 + chgpos + totpos1) ) THEN
      IF ( Sign(1.0d0,dble(totneg1+totpos1)) /= Sign(1.0d0,dble(chgneg + totneg1 + chgpos + totpos1)) ) THEN
       ! check whether the sum of the averages maintains the same sign as the total charge. If not, then treat as IC flash with conserved charge
        chgfracn = Min(chgfracn,chgfracp)
        chgfracp = chgfracn
      ENDIF
      ! need to test these checks:
    !  IF ( chgfracp .gt. 0.5*chgpos ) chgfracp = 0.5*chgpos
    !  IF ( chgfracn .gt. 0.5*Abs(chgneg) ) chgfracn = 0.5*Abs(chgneg)

      if (abs(totneg1).gt.0.) then
      x = (chgneg+chgfracn)/chgneg
      else
      x=0.
      endif

      if (abs(totpos1).gt.0.) then 
      y = (chgpos-chgfracp)/(Abs(chgpos))
      else
      y=0.
      endif
      
!        IF ( x > 0.9 .or. y > 0.9 ) THEN
!          x = x * 0.9/(Max(x,y))
!          y = y * 0.9/(Max(x,y))
!        ENDIF

        IF ( x < 0.1 .or. y < 0.1 ) THEN ! use old method
         write(0,*) 'small X,Y: ', x,y,chg,chgpos,chgneg,totpos1,totneg1
      chg = elgt1*Max(chgpos, Abs(chgneg))
!      chg = (1.0 - elgt1)*Max(chgpos, Abs(chgneg))
      IF ( chg .gt. chgpos ) chg = chgpos
      IF ( chg .gt. Abs(chgneg) ) chg = Abs(chgneg)

      if (abs(chgneg).gt.0.) then
!      x = chg/(Abs(chgneg)) 
      x = (Abs(chgneg)-chg)/(Abs(chgneg)) 
      else
      x=0.
      endif

      if (abs(chgpos).gt.0.) then 
!      y = chg/(Abs(chgpos))
      y = (chgpos-chg)/(Abs(chgpos))
      else
      y=0.
      endif
        ENDIF
      
      ENDIF

      write(0,*) 'X,Y', x,y,chg,chgpos,chgneg,totpos1,totneg1
      write(0,*) 'chgfracn,chgfracp = ',chgfracn,chgfracp

       totneg2 = 0.0
       totpos2 = 0.0
       totneg3 = 0.0
       totpos3 = 0.0
       chglost = 0.0
       totnegbef = 0.0
       totnegaft = 0.0
       totposbef = 0.0
       totposaft = 0.0
       

      DO jy = 1,ny
      do kz = 1,nz
      DO ix = 1,nx

         dv = dxx*dyy*dzz(ix,kz,jy)

         totnegbef = totnegbef + Min(0.0, sctot(ix,kz,jy) )*dv
         totposbef = totposbef + Max(0.0, sctot(ix,kz,jy) )*dv

!         totnegaft = totnegaft + Min(0.0, sch(ix,kz,jy) )*dv*db(ix,kz,jy)
!         totposaft = totposaft + Max(0.0, sch(ix,kz,jy) )*dv*db(ix,kz,jy)

      IF ( ij_reduce(ix,jy) > 0 ) THEN


!      NEED THIS IN C/kg and have C/m3 --> divide by RHO

!      chg = chg*(1.0 - elgt1)
!      IF ( t1(ix,jy,kz) .lt. 0.0 ) THEN ! local chg is neg. so add pos. ions
      IF ( sctot(ix,kz,jy) .lt. -chgthr ) THEN ! local chg is neg. so add pos. ions
!        an(ix,jy,kz,lscpi) = an(ix,jy,kz,lscpi) - t1(ix,jy,kz)*x
!        scdisch(ix,jy,kz) = (-(sctot(ix,kz,jy)+chgthr)*y )/ db(ix,kz,jy)!x ! OPPOSITE TO DECREASE IT !
        scdisch(ix,jy,kz) = (-(sctot(ix,kz,jy)+chgthr)*(1.0 - x) )/ db(ix,kz,jy)!x ! OPPOSITE TO DECREASE IT !
        totneg2 = totneg2 + scdisch(ix,jy,kz)*dv*db(ix,kz,jy)
!      ELSEIF ( t1(ix,jy,kz) .gt. 0.0 ) THEN ! local chg is pos. so add neg. ions
      ELSEIF ( sctot(ix,kz,jy) .gt. chgthr ) THEN ! local chg is pos. so add neg. ions
!        an(ix,jy,kz,lscni) = an(ix,jy,kz,lscni) + t1(ix,jy,kz)*y
!        scdisch(ix,jy,kz) = ( (sctot(ix,kz,jy) - chgthr)*x )/db(ix,kz,jy) !y
        scdisch(ix,jy,kz) = ( (sctot(ix,kz,jy) - chgthr)*(1.0 - y) )/db(ix,kz,jy) !y
        totpos2 = totpos2 - scdisch(ix,jy,kz)*dv*db(ix,kz,jy)
      ENDIF
      
       if (abs(scdisch(ix,jy,kz)).gt.1.e-3) write(0,*) 'scdisch,sctot', scdisch(ix,jy,kz),sctot(ix,kz,jy)     


       IF ( scdisch(ix,jy,kz) /= 0.0 ) THEN

        IF ( Abs(scdisch(ix,jy,kz)) > 0.1*chgthr ) lightdis(ix,jy) = lightdis(ix,jy) + 1
       
       IF ( .true. ) THEN ! put the lightning charge into sciona, which is then distributed to hydrometeors in the ionattach subroutine
       
        
      if (sctot(ix,kz,jy) .gt.  chgthr) sciona(ix,kz,jy)=sciona(ix,kz,jy)-scdisch(ix,jy,kz)
      if (sctot(ix,kz,jy) .lt. -chgthr) sciona(ix,kz,jy)=sciona(ix,kz,jy)+scdisch(ix,jy,kz)
      
        IF ( present( schl ) ) THEN
          sctot(ix,kz,jy) = db(ix,kz,jy)*(schl(ix,kz,jy)+sch(ix,kz,jy)+scr(ix,kz,jy)+sci(ix,kz,jy)+scs(ix,kz,jy)+scw(ix,kz,jy)+sciona(ix,kz,jy))
        ELSE
          sctot(ix,kz,jy) = db(ix,kz,jy)*(sch(ix,kz,jy)+scr(ix,kz,jy)+sci(ix,kz,jy)+scs(ix,kz,jy)+scw(ix,kz,jy)+sciona(ix,kz,jy))
        ENDIF
       ELSE

!     compute area fraction Ax/sum(Ax). where Ax=(Qx/rhox)**2/3 * Nx 
!   but qnr etc. is N/rhoair, and need to convert scdisch to C/kg
!   So have Ax=(Qx*rhoair/(Nx*rhoair*rhox))**2/3 * (Nx*rhoair) = (Qx/(Nx*rhox))**2/3 * (Nx*rhoair)
!  or Ax = rhoair*(Qx/(rhox))**2/3 * (Nx)**(1/3), and the factor of rhoair cancels out in the ratio Ax/Atot
    


      fracs = 0.0
      fraci = 0.0
      frach = 0.0
      fracw = 0.0
      fracr = 0.0
      frachl = 0.0
      if ( qs(ix,kz,jy).gt.qmin .and. qns(ix,kz,jy) > 1.e-9 ) fracs=((qs(ix,kz,jy)/rhos)**(2./3.) ) * qns(ix,kz,jy)**(1./3.)
      if ( qi(ix,kz,jy).gt.qmin .and. qni(ix,kz,jy) > 1.e-9 ) fraci=((qi(ix,kz,jy)/rhoi)**(2./3.) ) * qni(ix,kz,jy)**(1./3.)
      if ( qg(ix,kz,jy).gt.qmin .and. qng(ix,kz,jy) > 1.e-9 ) frach=((qg(ix,kz,jy)/rhog)**(2./3.) ) * qng(ix,kz,jy)**(1./3.)
      if ( qc(ix,kz,jy).gt.qmin .and. qnc(ix,kz,jy) > 1.e-9 ) fracw=((qc(ix,kz,jy)/rhow)**(2./3.) ) * qnc(ix,kz,jy)**(1./3.)
      if ( qr(ix,kz,jy).gt.qmin .and. qnr(ix,kz,jy) > 1.e-9 ) fracr=((qr(ix,kz,jy)/rhor)**(2./3.) ) * qnr(ix,kz,jy)**(1./3.)

      deninv = 1. ! /db(ix,kz,jy)
      fractot=fracs+fraci+frach+fracr+fracw


      if (present(qh)) then
        IF ( qh(ix,kz,jy).gt.qmin .and. qns(ix,kz,jy) > 1.e-9 ) frachl=((qh(ix,kz,jy)/rhoh)**2./3. ) * qnh(ix,kz,jy)**(1./3.)
        fractot = fractot + frachl
      endif

     
      IF ( fractot > 0.0 ) THEN
!'     reduce/redistribute scdisch on scx based on area fraction 
   
      if (qs(ix,kz,jy).gt.qmin) then ! avoid division by zero if fractot = 0 because qx=0
      if (sctot(ix,kz,jy) .gt.  chgthr) scs(ix,kz,jy)=scs(ix,kz,jy)-scdisch(ix,jy,kz)*(fracs/fractot)
      if (sctot(ix,kz,jy) .lt. -chgthr) scs(ix,kz,jy)=scs(ix,kz,jy)+scdisch(ix,jy,kz)*(fracs/fractot)
      endif

      if (qi(ix,kz,jy).gt.qmin) then
      if (sctot(ix,kz,jy) .gt.  chgthr) sci(ix,kz,jy)=sci(ix,kz,jy)-scdisch(ix,jy,kz)*(fraci/fractot)
      if (sctot(ix,kz,jy) .lt. -chgthr) sci(ix,kz,jy)=sci(ix,kz,jy)+scdisch(ix,jy,kz)*(fraci/fractot)
      endif

      if (qr(ix,kz,jy).gt.qmin) then
      if (sctot(ix,kz,jy) .gt.  chgthr) scr(ix,kz,jy)=scr(ix,kz,jy)-scdisch(ix,jy,kz)*(fracr/fractot)
      if (sctot(ix,kz,jy) .lt. -chgthr) scr(ix,kz,jy)=scr(ix,kz,jy)+scdisch(ix,jy,kz)*(fracr/fractot)
      endif

      if (qg(ix,kz,jy).gt.qmin) then
       tmp = sch(ix,kz,jy)
      if (sctot(ix,kz,jy) .gt.  chgthr) sch(ix,kz,jy)=sch(ix,kz,jy)-scdisch(ix,jy,kz)*(frach/fractot)
      if (sctot(ix,kz,jy) .lt. -chgthr) sch(ix,kz,jy)=sch(ix,kz,jy)+scdisch(ix,jy,kz)*(frach/fractot)
!            IF ( .not. ( sch(ix,kz,jy) > -1.e-6 .and. sch(ix,kz,jy) < 1.e-6 ) ) THEN
!              write(0,*) 'discharge: bad sch! ix,kz,jy, sch = ',ix,kz,jy,tmp,sch(ix,kz,jy),qg(ix,kz,jy),qng(ix,kz,jy)
!              write(0,*) 'scdisch(ix,jy,kz),frach,fractot,rhog: ',scdisch(ix,jy,kz),frach,fractot,rhog
!              write(0,*) 'parts= ',(qg(ix,kz,jy)/rhog)**(2./3.),qng(ix,kz,jy)**(1./3.), (qg(ix,kz,jy)/rhog)**(2./3.)*qng(ix,kz,jy)**(1./3.)
!              write(0,*) 'fracs,fraci,frach,fracw,fracr,frachl: ',fracs,fraci,frach,fracw,fracr,frachl
!              STOP
!            ENDIF
      endif

      if (qc(ix,kz,jy).gt.qmin) then
      if (sctot(ix,kz,jy) .gt.  chgthr) scw(ix,kz,jy)=scw(ix,kz,jy)-scdisch(ix,jy,kz)*(fracw/fractot)
      if (sctot(ix,kz,jy) .lt. -chgthr) scw(ix,kz,jy)=scw(ix,kz,jy)+scdisch(ix,jy,kz)*(fracw/fractot)
      endif
   
      if (present(qh)) then
      if (qh(ix,kz,jy).gt.qmin) then
      if (sctot(ix,kz,jy) .gt.  chgthr) schl(ix,kz,jy)=schl(ix,kz,jy)-scdisch(ix,jy,kz)*(frachl/fractot)
      if (sctot(ix,kz,jy) .lt. -chgthr) schl(ix,kz,jy)=schl(ix,kz,jy)+scdisch(ix,jy,kz)*(frachl/fractot)
      endif
      endif
      
       sctot(ix,kz,jy) = db(ix,kz,jy)*(sch(ix,kz,jy)+scr(ix,kz,jy)+sci(ix,kz,jy)+scs(ix,kz,jy)+scw(ix,kz,jy)+sciona(ix,kz,jy))
       if (present(qh)) sctot(ix,kz,jy)=sctot(ix,kz,jy)+db(ix,kz,jy)*schl(ix,kz,jy)

      ELSE
      ! have scdisch /= 0, but fractot=0
       chglost = chglost + scdisch(ix,jy,kz)*dv*db(ix,kz,jy)
      

      ENDIF ! fractot > 0.0

      ENDIF ! true/false
      
      ENDIF ! IF scdisch /= 0

      ENDIF ! ij_reduce

!     UPDATE sctot     


         totneg3 = totneg3 + Min(0.0, sctot(ix,kz,jy) )*dv
         totpos3 = totpos3 + Max(0.0, sctot(ix,kz,jy) )*dv

!         totnegaft = totnegaft + Min(0.0, sctot(ix,kz,jy) )*dv
!         totposaft = totposaft + Max(0.0, sctot(ix,kz,jy) )*dv
     
      end do
      end do
      end do

#ifdef DM_PARALLEL
       mpitotindp(1)  = totpos3
       mpitotindp(2)  = totneg3
       mpitotindp(3)  = totpos2
       mpitotindp(4)  = totneg2
       mpitotindp(5)  = chglost
       mpitotindp(6)  = totposbef
       mpitotindp(7)  = totnegbef
       mpitotindp(8)  = totposaft
       mpitotindp(9)  = totnegaft

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 9, MPI_DOUBLE_PRECISION, MPI_SUM, local_communicator, mpi_error_code)


      IF ( mytask == 0 ) THEN ! only print from task 0

       totpos3 = mpitotoutdp(1)
       totneg3 = mpitotoutdp(2)
       totpos2 = mpitotoutdp(3)
       totneg2 = mpitotoutdp(4)
       chglost = mpitotoutdp(5)
       totposbef = mpitotoutdp(6)
       totnegbef = mpitotoutdp(7)
       totposaft = mpitotoutdp(8)
       totnegaft = mpitotoutdp(9)

      ENDIF
#endif

!  Reserved for future use
      IF ( mytask == 0 ) THEN ! only print from task 0
        write(0,'(a,3(2x,1pe12.5))') 'pos/neg1: ',totpos1,totneg1,totpos1+totneg1
        write(0,'(a,3(2x,1pe12.5))') 'pos/neg2 scdisch: ',totpos2,totneg2,totpos2+totneg2
        write(0,'(a,3(2x,1pe12.5))') 'pos/negbef: ',totposbef,totnegbef,totposbef+totnegbef
        write(0,'(a,3(2x,1pe12.5))') 'pos/neg3: ',totpos3,totneg3,totpos3+totneg3
!        write(iunit,'(a,3(2x,1pe12.5))') 'pos/negaft: ',totposaft,totnegaft,totposaft+totnegaft
        write(0,'(a,1(2x,1pe12.5))') 'chglost: ',chglost
      ENDIF

!        write(0,*) '+============AFTER DISCHARGE============+'
!        write(0,*) 'SCTOT MIN', MINVAL(sctot(:,:,:))
!        write(0,*) 'SCTOT MAX', MAXVAL(sctot(:,:,:))
!        write(0,*) 'ELEC MAG MAX',MAXVAL(elec(1:nx,1:nz,1:ny)),MAXLOC(elec(1:nx,1:nz,1:ny))
!        write(0,*) 'SCW MIN', MINVAL(scw(:,:,:))
!        write(0,*) 'SCW MAX', MAXVAL(scw(:,:,:))
!        write(0,*) 'SCR MIN', MINVAL(scr(:,:,:))
!        write(0,*) 'SCR MAX', MAXVAL(scr(:,:,:))
!        write(0,*) 'SCH MIN', MINVAL(sch(:,:,:))
!        write(0,*) 'SCH MAX', MAXVAL(sch(:,:,:))
!        write(0,*) 'SCI MIN', MINVAL(sci(:,:,:))
!        write(0,*) 'SCI MAX', MAXVAL(sci(:,:,:))
!        write(0,*) 'SCS MIN', MINVAL(scs(:,:,:))
!        write(0,*) 'SCS MAX', MAXVAL(scs(:,:,:))
!        if (iter .gt.0) write(0,*) 'SCDIS MAX', MAXVAL(scdisch(:,:,:))
!        if (iter .gt.0) write(0,*) 'FRACTOT MAX', MAXVAL(fractot(:,:,:))
      
!
!  write out that there is  a lightning discharge
!
!!      IF ( mytask == 0 ) THEN
      write(iunit,*) '______________________________'
      write(iunit,*) 'LIGHTNING DISCHARGE:  SCHEME 1'
!!      write(iunit,*) 'TIME STEP IS: ',time_real
      write(iunit,*) 'TOTAL NUMBER OF DISCHARGES IS: ',nlgt1
      write(iunit,*) 'chgpos,chgneg : ',chgpos,chgneg
!      write(iunit,*) '______________________________'
!!      ENDIF
!

      ENDIF ! etot .gt elght1 !!!!
!
      nlgt1 =0      

      IF ( allocated(xy_init) ) deallocate( xy_init   )
      IF ( allocated(ij_reduce) ) deallocate( ij_reduce )
      deallocate( scdisch )
      deallocate( w_proc,e_proc,n_proc,s_proc )

      write(iunit,*) 'END OF LIGHTNING DISCHARGE'

      end subroutine light1d

      end module module_mp_discharge
