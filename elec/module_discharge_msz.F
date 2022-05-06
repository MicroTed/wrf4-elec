!
#define FLOAT_PRECISION MPI_DOUBLE_PRECISION

      Module  module_discharge_msz
#ifdef DM_PARALLEL
#define MPI 1
#endif
      
      private ran0mpi
      
      CONTAINS

      subroutine lightmsz(iliter,lgtstp,loccur,dtp,nstep,ilight,elgtfestop,cnoxtot,   &
     &   nor,nx,ny,nz,ids,ide,jds,jde,kds,kde,its,jts,      &
     &   ims,ime,jms,jme,kms,kme,temg    &
     &  ,dxx,dyy,lightrad,elght1,elgt1,dzz,sctot,pot &
     &  ,scw,scs,sci,scr,sch,schl &
     &  ,qc,qs,qi,qr,qg,qh &
     &  ,qnc,qns,qni,qnr,qng,qnh &
     &  ,elec        &
     &  ,sciona            &
     &  ,nbw,nbe,nbs,nbn,nbw0,nbe0,nbs0,nbn0    &
     &  ,iunit0,db,pn,light,lightdens,etot,lightdis  &
     &  ,flshfedic,flshfedicp,flshfedicn,flshfedcg,flshfedcgp,flshfedcgn &
     &  ,flshi,flshn,flshp,lightmask,lightinit,initmask &
     &  ,eint,tgrnd,zgrnd &
     &  ,f_clnox,ixbeg,jybeg,flashcount,clnox,kcldtop)
     
#ifdef DM_PARALLEL
      USE module_dm, ONLY : &
         local_communicator, mytask, ntasks, ntasks_x, ntasks_y,mytask_x, mytask_y
      USE module_commasmpi ! , ONLY : w_proc,e_proc,s_proc,n_proc,mpi_setup_flag
#endif
      use module_boxmgsetup, only: igslg0, jgslg0, kgslg0

      implicit none

#ifdef DM_PARALLEL
      INCLUDE 'mpif.h'
      INTEGER :: my_mpi_status(MPI_Status_size)
#else
      integer :: mytask = 0
#endif
      
      integer, intent(IN) :: ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,its,jts
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in)::  temg
      integer iunit,iunit0
      real elec(0:nx+1,nz,0:ny+1,4)
!      real scion(0:nx+1,nz,0:ny+1,2)  ! 2=+, 1=- 
      real dzz(0:nx+1,nz,0:ny+1),dxx,dyy   ! dzz is dz8w
      real db(0:nx+1,nz,0:ny+1)  ! air density (kg m**-3)
      real pn(0:nx+1,nz,0:ny+1)  ! air pressure
      integer, parameter    :: ng1 = 1
      integer nbw,nbe,nbs,nbn,nbw0,nbe0,nbs0,nbn0
      
      double precision chg,chgneg,chgpos,chgtot,chgpnts,cnoxtot ! total net negative/positive charge in storm (not including corona)
      real  :: chgmax
      
      integer ix,jy,kz,kz1,i,j,k,count,counttot,n,m
      integer nor,nx,ny,nz
      real pi
      integer ietotx,ietoty,ietotz
      integer :: ixbeg,jybeg
      integer :: flashcount

!  ibrkd:    ibrkd is they of break down; 1=based on emax; 2=based
!            on emax(z)

      real  elght1,lightrad,elgt1

      integer nlgt1
      integer :: ibrkd    = 4
!      integer :: lightrad = 6000.0
      real    :: chgthr   = 0.1e-9  ! threshold charge density
!      real    :: elght1   = 120000. ! threshold e-field to init lightning
!      real    :: elgt1    = 0.3     ! fractional reduction in charge in excess of chgthr
      real,parameter    :: eperao=8.8592e-12
      real,parameter    :: ec = 1.602e-19 ! fundamental unit of charge 
      real,parameter    :: rho00 = 1.225
      real    :: qmin = 1.e-13 !1.e-5
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
!      real t0(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real pot(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
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
      
      integer kcldtop(-nbw:nx+1+nbe,-nbs:ny+1+nbn)

      real frachl
      real frach
      real fracr
      real fracw
      real fraci
      real fracs
      real fractot
      double precision dv, deninv, tmp
      integer ncldtop
      real zcldtoptot

!      real t2(-nbw:nx+1+nbe,-nbs:ny+1+nbn,nz)
!      real t4(-nbw:nx+1+nbe,-nbs:ny+1+nbn,nz)
!      real, allocatable, save :: t0(:,:,:) ! t2(0:nx+1,0:ny+1,0:nz+1)
      real, allocatable, save :: t2(:,:,:) ! t2(0:nx+1,0:ny+1,0:nz+1)
      real, allocatable, save :: t4(:,:,:) ! t4(0:nx+1,0:ny+1,0:nz+1)
      ! t3 -> scdisch (is i,j,k)

!      real,save :: zlev(nz,-nbw:nx+1+nbe,-nbs:ny+1+nbn)
      real :: z1d2(nz,4)
      real,allocatable,save :: ebrkd(:,:,:),zlev(:,:,:)
      real,allocatable,save :: cghgt(:,:)
!      real :: ecrit(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real :: light(-nbw:nx+1+nbe,-nbs:ny+1+nbn)
      real :: lightdis(-nbw:nx+1+nbe,-nbs:ny+1+nbn)
      integer :: lightinit(-nbw:nx+1+nbe,-nbs:ny+1+nbn)
      integer :: lightmask(-nbw:nx+1+nbe,-nbs:ny+1+nbn)
      integer :: initmask(-nbw:nx+1+nbe,-nbs:ny+1+nbn,nz)
      real :: lightdens(-nbw:nx+1+nbe,-nbs:ny+1+nbn)
      REAL, DIMENSION( ims:ime , jms:jme ) , INTENT(INOUT)   :: flshfedic,flshfedicn,flshfedicp
      REAL, DIMENSION( ims:ime , jms:jme ) , INTENT(INOUT)   :: flshfedcg,flshfedcgn,flshfedcgp

      real,allocatable,save :: ebrkdp(:,:,:) !(-nbw:nx+1+nbe,-nbs:ny+1+nbn,nz)

      real flshi(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real flshn(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real flshp(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real, optional :: clnox(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn) ! moles/kg of LNOx
      real nox(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn),lnox(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real noxtmp
      logical :: f_clnox  ! whether the clnox array really exists or not
      
      real, allocatable :: scdisch(:,:,:) 
      real, allocatable :: scdisch0(:,:,:) 
      real,   allocatable  :: xy_init(:,:)
      integer, allocatable :: ij_reduce(:,:)
      real, allocatable, save :: flsh_map(:,:,:)
      integer, allocatable,save  :: ijmask(:,:)
      integer, save :: nxm, nym
      integer  :: nr, nrflsh, ninit,ninit1
      real  :: x,y,r

      integer nttim,nenum,ninfo
      integer llx,lly,llz,lex,ley,lez
      integer lemag,lpot,lstop,lemag2,lcharge
      integer, parameter :: itest = 0
      integer, parameter :: ndebug = 1
      integer lgtinit
      integer ieout,ndiscg
      integer ndisic,ndisca
      real    zsfc
      integer iliter,lgtstp,ilghtest
      integer iopen,ifail
      data iopen /0/
      integer ilgt
      integer ixmx,jymx,kzmx,ixmn,jymn,kzmn
      integer iemaxx,iemaxy,iemaxz
      integer lgtgo
      integer lightsuccess
      integer islgt,jslgt,kslgt
      integer islgt_local,jslgt_local,kslgt_local
      integer ixst,jyst,kzst
      integer ixst0,jyst0,kzst0
      real    t0mx,t0mn,scramt,cndamt,times
      integer in,it,ic,jc,kc,iep,je,ke     ! local tile-relative indices
      integer :: icdom,jcdom,iepdom,jedom  ! domain-relative indices
      integer  idebug,idtsign, isign1,isign2
      real    eps,facx,facy,facz,facxe,facye,facze
      parameter (eps = 1.0e-7)
      real    telec,tcharge,dtcomp
      integer i1,j1,k1,ichneu,idist,ipos,ineg
      real :: rdist
      real :: e1,e2,tem
      integer ichpos,ichneg,sizetest
      real    scthend1,scthend2,scthendp,scthendn
      real    scthp,scthn,potthrp, potthrn
      integer isendn,jsendn,ksendn,isendp,jsendp,ksendp
      integer iic,ica,icg,icgp,icgn,icgpn,i2,j2,k2,icgpflg,icgnflg
      integer :: iredis = 1
      integer, save :: numic = 0, numca = 0, numcg = 0, numlgt = 0, numcgn = 0, numcgp = 0
      integer :: iicp = 0
      integer iicn,itype, icountp,icountn
      real    rndnum,cons1,cons2,scnetpi,scnetni,volzone
      integer irepeat,nscnetp,nscnetn,nscnet,nscnetnp
      real    voln,volp,scnetp,scnetn,scnet,scnetnp,scnetnn
      integer nscnetnn,ia
      real :: efracinitwire = 1.0
      real :: overvolt = 1.0
      integer :: nliter2
      integer,save :: iseed    = 789, imkseed = 0
      real :: ran
      
      character*40 flgtchn
      character*6 rstim,htail
      character*2 rnl
      character*11 yesno
      
      real zgrnd,tgrnd,eint   !  threshold height to declare channel end at ground
!      parameter (zgrnd = 2000.0)
      integer ichend(2)   ! tells whether channel end at
                         ! 1 = ground
                         ! 2 = cloud
                         ! 3 = air
!
      integer lc,lr,li,lip,lir,ls,lgl,lgm,lgh,lf,lh,lhl,lschl
      integer lscw,lscr,lsci,lscip,lscir,lscs
      integer lscgl,lscgm,lscgh,lscf,lsch
!      integer nxelec,nyelec,nzelec,neelec
      integer icgyn
      integer ipend,inend
      real qci  ! interpolated value of cloud content at end of channel
      real chgmin
!      
!  vars for dipole calculation
!
      double precision qpx,qpy,qpz    ! dipole moment vector components
      double precision qpxn,qpyn,qpzn ! negative drho centroid
      double precision qpxp, qpyp, qpzp  ! positive drho centroid
      double precision qpxm, qpym, qpzm  ! dipole midpoint location
      double precision qpxcn,qpycn,qpzcn
      double precision qpxcp,qpycp,qpzcp  ! neg and pos centroid locations
      double precision qnetn,qnetp
      integer npntn,npntp
      double precision fracn,fracp,cgfrac,chgrtp,chgrtn,chgrtp2,chgrtn2
      double precision totneg1,totpos1,totneg2,totpos2
      double precision totneg0,totpos0

      real multp,multn
      real sctarg
      integer npsurfp,npsurfn  ! number points inside pot. surface
      integer nppntp,nppntn   ! number of potential flash points
      double precision :: total_nox
!      
!  vars for lightning channel 'cone'
!
!    rbase,theta define initial radius and angle of cone (read from inlight2)
!
!        channel
!           ^        edge of cone
!  \    |   |   |    /
!   \   |   |   |   /
!    \  |   |   |  /
!     \ |   |   |t/    (t=theta)
!      \|___|___|/
!           |<->|
!           rbase
      
      integer imx,jmy,kmz
      real rtest,chanlen,rchan,tanth,dl
      !  rtest is the distance from a gridpoint to the channel test point
      !  chanlen is the length along the channel of the test point
      !  rchan = rbase + chanlen*Tan(theta)
      !  tanth = Tan(theta)
      integer nchadd(2)  !  number of points added along channel
      integer isignt
      integer ichcnt
      real chcnt
      real delx,dely,delz  !  relative coords of channel point from corner
      
      real scdrop, sckeepp, sckeepn
      real, save :: scdropnet = 0.0
      integer idrop
      
! variables from 'inlight2'
      real    :: elgtthx  = 100.0e3          ! elgtthx (e-max before discharge)
      real    :: elgtthn  =  90.0e3          ! elgtthn (end of lightning 'bolt' at e<elgtthn)
      real    :: elgtthncg=  15.0e3          ! elgtthncg (end of lightning 'bolt' at e<elgtthncg to check for CG flash)
      real    :: elgtdel  = 10.0e3           ! elgtdel (discharge possible [random] where emax-elgtthn<e) 
      real    :: elgtfdel =  0.9             ! elgtfdel e magnitude fraction of delta for [f(z)] breakdown
      real    :: scth     = 0.11e-9 ! 0.08e-9          ! scth (space charge threshold for discharge region +/- scth)
      real    :: fscth    = 0.10             ! fscth (fraction of scth removed in iteration part of scheme)
      real    :: scthmn   = -1.0e-12         ! scthmn (minimum value [negative space charge] for work)
      real    :: scthmx   = 1.0e-12          ! scthmx (maximum value [positive space charge] for work)
      real    :: fprcnt   = 0.70             ! fprcnt (percent of charge removed in discharge volume)
      integer :: niter    = 2                ! niter (number of iterations on space charge threshold)
      integer :: icountth0 = 04               ! icountth (number of points required in discharge region)
      integer :: icountth
      integer :: nretry   = 06               ! nretry (number of channel retries)
      real    :: rbase    = 250.0            ! rbase (meters)
      real    :: theta    = 22.0             ! theta (degrees)
      integer :: ibal     = 1                ! 1=allow unbalanced IC flash when it reduces the net charge in a local region around the flash
                                             ! 2=strict charge conservation for IC flashes
      real    :: cgfr     = 1.0
      real    :: elgtfestop                  ! fraction of Emag for stopping channel
      real    :: elgtfestopcg = 0.15         ! fraction of Emag for stopping channel for CG
      
      real, parameter :: delq0 = 1.e-9       ! reference charge density change to scale NOx production
      
      integer :: idownward
      integer :: ifoundic

      integer, save :: itermax = 0 ! maximum iterations for flash search
      integer, save :: maxdist = 0 ! maximum rectilinear distance for flash search
      integer, save :: ibalsize = 0 ! size of 1D array for tracking volume/charge by iteration
      integer :: itern, iterp ! which iteration of neg/pos channel volume to stop at
      real, allocatable, save :: chgavailiter(:),chgavailitersum(:) ! array to sum up available lightning charge in each wildfire iteration
      real, allocatable, save :: voliter(:),volitersum(:)
      real :: chanzn(nz),chanzp(nz), channtot,chanptot

      integer iretry   ! counter for nretry
      integer kzbrkd,iread
      integer loccur,nliter
    !  integer :: itermax ! maximum iterations for flash search
      integer nstep,nxl,nyl,nzl,ngt,na
      integer :: nxend, nyend
      integer,parameter :: id1=1,jd1=1,kd1=1,istag=1,jstag=1,kstag=1
      integer :: ilight
      integer iex,iey,iez,iemag,ipot,ie
      integer, save :: klgtmin = 1 ! lowest level to extend lightning: can set higher than 1 if explicit ions are used (to avoid corona layers)
      real, parameter :: zlgtmin = 600.
      integer imapx,imapy,imapz,ig13,ig23,igsr

!       integer ninit  ! number of possible initiation points
       integer, save :: npick = 0
!       parameter(npick = 10000)
       integer, allocatable, save :: ipick(:,:) ! ipick(npick,3)
       integer ichoose
       integer idum

      integer :: ichooserank = 0
      integer :: ichoosenumber
      integer :: ichoosenextrank

      integer ioffset, joffset, koffset
      integer ishift, jshift, kshift
      integer NLxchoose, NLychoose, NLzchoose
      
      integer, parameter :: isa = 0
      
! variable of lnox calculation
      integer itt
      real amin,bmin,cmin
      real alnox,blnox  ! parameters from Wang et al. 1998 (JGR)
      real avogadro     ! Avogadro's number
      real air_molmass  ! molecular mass of air
      parameter(alnox=0.34e21,blnox=1.30e16,avogadro=6.0221415e23,air_molmass=28.96e-3)

!
!  nttim:  maximum number of intervals in a lightning channel
!  nenum:  number of ends to a channel from a common location
!          (a positive and a negative end)
!  ninfo:  number of parameters in the trje array 
!
      parameter  (nttim=4000,nenum=2,ninfo=12)
!
!  llx:    x location of trajectory of lightning channel in trje
!  lly:    y location of trajectory of lightning channel in trje
!  llz:    z location of trajectory of lightning channel in trje
!  lex:    x component of e-field along lightning channel in trje
!  ley:    y component of e-field along lightning channel in trje
!  lez:    z component of e-field along lightning channel in trje
!  lemag:  magnitude of e along lightning channel in traj
!  lpot:   potential along lightning channel in traj
!  lstop:  stop indictator along lightning channel in traj
!  lemag2: magnitude from lex,y,z of e along lightning channel in traj
!  lcharge:  charge along lightning channel in trje
!  ndebug:  ndebug parameter (1=on)
!  itest:   itest    parameter (1=on)
!  ilgtinit: iligtinit is a parameter that must be = 1 for lightning
!            discharges to be calculated in the code
!  ibrkd:    ibrkd is they of break down; 1=based on emax; 2=based
!            on emax(z)
!  ieout:    ieout is interval to output lightning information
!
      parameter  (llx=1,lly=2,llz=3,lex=4,ley=5,lez=6)
      parameter  (lemag=7,lpot=8,lstop=9,lemag2=10)
      parameter  (lcharge=11)
      integer, parameter :: ltemg = 12
      parameter  (lgtinit=1)
!      ,ibrkd=2)
      parameter  (ieout=1)
!
      integer :: find_indexlgt

!
!  ndis** = 1 (redistribute charge based on sfc area)
!  ndis** = 2 (redistribute charge based just removing a fraction)
!    cg = cloud-ground
!    ic = intracloud
!    ca = cloud-air
!
      parameter  (ndiscg=3)
      parameter  (ndisic=3)
      parameter  (ndisca=1)
!
      integer    nttimend(nenum),idxlgt(nenum,nttim,3),nttimendcg(nenum)
!
!
      real dtelec(nenum),dt,dtp,dx,dy
      real dxi,dyi,dzi
      real dxi2,dyi2,dzi2
!
      integer nemaxx(1000),nemaxy(1000),nemaxz(1000)
      integer kelec(1000)
!      real emax, zlev(1000), ebrkd(1000), ebrkdp(1000)
      real emax, emaxkz(1000),emaxlightkz(1000),exmaxkz(1000),eymaxkz(1000),ezmaxkz(1000)
      real emag, ez, zinit, pot1
!      save zlev,ebrkd,ebrkdp
      real estop(1000)
      real estopcg(1000)
      real trje(nenum,0:nttim,ninfo)


!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES

      integer, parameter :: ntot = 50
      integer :: ichoose_send_buffer(ntot)
      integer :: ichoose_recv_buffer(ntot)
      real    :: rchoose_send_buffer(ntot)
      integer :: mpitotinint(ntot), mpitotoutint(ntot)
      real  mpitotin(ntot), mpitotout(ntot)
      double precision  mpitotindp(ntot), mpitotoutdp(ntot)

      integer :: westward_tag, eastward_tag
      integer :: northward_tag, southward_tag
      INTEGER :: mpi_error_code = 1      
!      INTEGER,DIMENSION(:),ALLOCATABLE :: w_proc, e_proc, s_proc, n_proc
      INTEGER :: memory_status,process
      integer :: my_rank = 0

      integer :: ixb, jyb, kzb
      integer :: ixe, jye, kze
      integer :: nz1
      integer :: imxmn
      integer :: irank1,irank2
   
      integer :: listsize

      integer :: mxmn_size,max_rank,min_rank
      integer :: loc_send_size,loc_recv_size
      integer :: exyz_send_size,exyz_recv_size

      integer,allocatable,save :: mpiloc(:,:), mpiloc0(:,:,:)
      integer,allocatable,save :: loc_send_buffer(:,:),loc_recv_buffer(:,:)

      real,allocatable,save :: mpiexyz(:,:), mpiexyz0(:,:,:)
      real,allocatable,save :: exyz_send_buffer(:,:),exyz_recv_buffer(:,:)

      real,allocatable,save :: max_send_buffer(:,:),max_recv_buffer(:,:)
      real,allocatable,save :: min_send_buffer(:,:),min_recv_buffer(:,:)

      real,allocatable,save :: mpimxmn(:,:)

      integer :: init_send_size, init_recv_size, ninittmp, ninittot
      
      integer,allocatable,save :: init_send_buffer(:)
      integer,allocatable,save :: init_recv_buffer(:,:)
      integer,allocatable,save :: mpi_init0(:,:)
      
      
      integer :: itile,myproci,nproci
      integer :: jtile,myprocj,nprocj
      integer :: number_of_processes,mpi_abort_error_code
      
      integer :: my_comm


! #####################################################################
!   Begin Executable Code
! #####################################################################
       dx = dxx
       dy = dyy
       tanth = Tan(theta*4.*atan(1.)/180.0)
       
       nxend = ide-ids+1
       nyend = jde-jds+1

       icountth = Min(icountth0, Max(2, Nint(icountth0*1000./dx)))



       IF ( .not. allocated( chgavailiter ) ) THEN
         itermax = Max( nxend, nyend, kde - kds)
         itermax = Min( itermax-1, Int(0.5*100000./dx) ) ! limit lightning half-width to 100 km
         IF ( ibal == 1 ) THEN
           maxdist = itermax
         ELSEIF ( ibal == 2 ) THEN
           maxdist = 3*itermax
         ELSE
           maxdist = itermax
         ENDIF
         allocate( chgavailiter(-maxdist:maxdist) )
         allocate( chgavailitersum(-maxdist:maxdist) )
         allocate( voliter(-maxdist:maxdist) )
         allocate( volitersum(-maxdist:maxdist) )

       ENDIF
       
       IF ( npick == 0 ) THEN
         npick = nxend*nyend*(kde - kds)/4
         allocate ( ipick(npick,3) )

       ENDIF


       iunit = iunit0
       
       IF ( mpi_setup_flag < 0 ) THEN
         CALL TASK_PROC_MAP()
       ENDIF

!      allocate( scdisch(-nbw:nx+1+nbe,-nbs:ny+1+nbn,nz) )
!      allocate( scdisch0(-nbw:nx+1+nbe,-nbs:ny+1+nbn,nz) )
       allocate( scdisch(0:nx+1,0:ny+1,0:nz+1) )
       allocate( scdisch0(0:nx+1,0:ny+1,0:nz+1) )

       scdisch(:,:,:)=0. 
       scdisch0(:,:,:)=0. 
       nlgt1 =0
       etot =0.
       lightsuccess = 0

!       write(0,*) ,'DIS: elght1,elgt1,lightrad=',elght1,elgt1,lightrad

#ifdef DM_PARALLEL 

              my_rank = mytask
              nproci = ntasks_x
              nprocj = ntasks_y
              myproci = mytask_x+1
              myprocj = mytask_y+1
              number_of_processes = ntasks

      IF ( .not. allocated(init_recv_buffer) ) THEN
      allocate(init_send_buffer(2))
      allocate(init_recv_buffer(number_of_processes,2))

      allocate(mpi_init0(2,0:number_of_processes-1))
      
      ENDIF
              
#endif

              itile = nx
              jtile = ny

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
!       zlev=0. ! dont want zero ecrit otherwise light might become crazy
       ilgt = 0

      IF (.not.allocated(ijmask)) THEN
       
!      create disks mask for cylinder instead of cuboids

      IF ( ilight == 1 ) THEN
        nxm = Int(abs(lightrad)/dxx) + 1
      ELSE
        nxm = Int(abs(Max(20000.,lightrad))/dxx) + 1
      ENDIF
      nym = nxm

      allocate( ijmask(-nxm:nxm,-nym:nym) )
      
      ijmask(:,:) = 0
      ix = 0

       IF ( lightrad <= 0.0 ) THEN
         ijmask(:,:) = 1
       ELSE
       ix = 0
       DO j = -nym,nym
        DO i = -nxm,nxm
!   assume dx=dy here
          r = dxx*Sqrt( Float(i**2 + j**2) )
          IF ( r <= Abs(lightrad) ) THEN
            ijmask(i,j) = 1
            ix = ix + 1
          ENDIF
        ENDDO
       ENDDO
       ENDIF
       IF ( ilight >= 2 ) THEN
         ijmask(-1:1,-1:1) = 2
       ENDIF
       
!      write(0,*) 'allocated ijmask. nxm = ',nxm,lightrad,ix

       ENDIF

      IF ( .not. allocated( ebrkd ) ) THEN ! {
       IF ( ndebug > 2 ) THEN
         write(0,*) 'msz allocate: nx,ny,nz,nbw,nbe,nbe,nbn = ',nx,ny,nz,nbw,nbe,nbe,nbn
       ENDIF
       allocate( ebrkd(-nbw:nx+1+nbe,-nbs:ny+1+nbn,nz) )
       allocate( ebrkdp(-nbw:nx+1+nbe,-nbs:ny+1+nbn,nz) )
       allocate( zlev(nz,-nbw:nx+1+nbe,-nbs:ny+1+nbn) )
       allocate( cghgt(-nbw:nx+1+nbe,-nbs:ny+1+nbn) )
!       allocate( t0(0:nx+1,0:ny+1,0:nz+1) )
       allocate( t2(0:nx+1,0:ny+1,0:nz+1) )
       allocate( t4(0:nx+1,0:ny+1,0:nz+1) )

        ebrkd(:,:,:) = 1.e8
        ebrkdp(:,:,:) = 1.e8
        
       IF ( zgrnd > 0.0 ) THEN 
         cghgt(:,:) = zgrnd
       ELSE
         cghgt(:,:) = 0.0 ! 
       ENDIF

       IF (ibrkd .ge. 2) THEN

      DO ix=-nbw, nx+1+nbe
      DO jy=-nbs, ny+1+nbn
      do kz = 1,nz
      
      zlev(kz,ix,jy) = sum(dzz(ix,1:kz,jy)) - 0.5*dzz(ix,kz,jy) ! dz8w is actually for w points?

!       IF ( zlev(kz,ix,jy) < zlgtmin ) THEN
!         klgtmin = Max(kz,klgtmin)
!       ENDIF

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

      IF ( ilight >= 2 ) THEN
        nr = 1
      ELSE !  ( ilight == 1 ) 
        IF ( lightrad < 0 ) THEN 
          nr = 0
        ELSE
          nr = NInt( lightrad/dxx )
        ENDIF
      ENDIF

      ! nrflsh = Max(2, NInt(lightextendmsz/dx))
      nrflsh = Max(2, NInt(6000./dx))
      
      allocate( xy_init(-nr+1:nx+nr,-nr+1:ny+nr))
      allocate( ij_reduce(-nr+1:nx+nr,-nr+1:ny+nr))
      IF ( .not. allocated( flsh_map ) ) THEN
        allocate( flsh_map(-nrflsh+1:nx+nrflsh,-nrflsh+1:ny+nrflsh,2) )
      ENDIF
      
      xy_init(:,:) = 0.0 ! zero out in order not to add lightning...even though these are deallocated.
      flsh_map(:,:,:) = 0.0

      IF ( lightrad < 0 ) THEN 
        ij_reduce(:,:) = 1
      ELSE
        ij_reduce(:,:) = 0
      ENDIF

      ixb = 1
      ixe = nx
      
      jyb = 1
      jye = ny


      IF ( ilight == 1 ) THEN
!
!  LIGHTNING DISCHARGE:  SCHEME 1.....
!
!  This scheme is the simplist of all.  If |E| => elgtt1 
!  (400 kV/m; REF) then the space charge is decreased by some percent
!  elgtp1 (30%:  Rawlins, 1981, QJRMS).
!  
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
      if ( elec(ix,kz,jy,4) .gt. etot ) then
      ietotx = ix
      ietoty = jy
      ietotz = kz
      etot = elec(ix,kz,jy,4)
      end if

      if ( elec(ix,kz,jy,4) .gt. elght1 ) then
!      if ( elec(ix,kz,jy,4) .gt. ebrkd(ix,jy,kz)   ) then
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
     &        w_proc,e_proc,westward_tag,xy_init)

        eastward_tag = 202
        CALL sendrecv_eastward(nx,ny,1,nr,nr,0,nr,1,  &
     &        w_proc,e_proc,eastward_tag,xy_init)

        ENDIF
        
        southward_tag = 203
        
        IF ( ntasks_y > 1 ) THEN
!        write(0,*) 'send_north, ntask_y = ',ntasks_y
        CALL sendrecv_southward(nx,ny,1,nr,nr,0,nr,1,  &
     &        n_proc,s_proc,southward_tag,xy_init)

        northward_tag = 204
        CALL sendrecv_northward(nx,ny,1,nr,nr,0,nr,1,  &
     &        n_proc,s_proc,northward_tag,xy_init)
        
        ENDIF
      ENDIF
#endif

       write(6,*) 'ETOT,ninit', ETOT,ninit
       loccur = 0
       IF ( ninit > 0 ) loccur = 1

!      IF ( mytask == 0 ) THEN
!        write(iunit,90) 'LOCAL IETOTX,IETOTY,IETOTZ=: ',ietotx,ietoty,ietotz
   90  format(A,3(3I,x))
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
          IF ( xy_init(ix,jy) > 0.5 ) THEN ! means that if elec(i,k,j,4) > elght1
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

        lightdis(ix,jy) = lightdis(ix,jy) + 1
        
        do kz = 1,nz-1
      
         dv = dxx*dyy*dzz(ix,kz,jy)
!
! t1 is the non-ion charge
        IF ( sctot(ix,kz,jy) .gt. chgthr ) THEN
          chgpos = chgpos + (sctot(ix,kz,jy) - chgthr)*dv 
        ELSEIF ( sctot(ix,kz,jy) .lt. -chgthr ) THEN
          chgneg = chgneg + (sctot(ix,kz,jy) + chgthr)*dv
        ENDIF
      
       end do
      
      ENDIF ! chgmax
      ENDIF ! } ( ij_reduce(ix,jy) > 0 ) 
      
      ENDDO
      ENDDO

!   lightning count here per tile - Note that light is real here.


      DO jy = 1,ny
      DO ix = 1,nx

!      light(ix,jy)=FLOAT(ij_reduce(ix,jy)) / ( pi*(lightrad*1.e-03)**2 )
       lightdens(ix,jy)=lightdens(ix,jy)+FLOAT(ij_reduce(ix,jy)) / ( pi*(lightrad*1.e-03)**2 )
       light(ix,jy)=light(ix,jy)+xy_init(ix,jy)

      ENDDO
      ENDDO


! find global integrated rate sum within cuboids
       mpitotindp(1)  = chgpos
       mpitotindp(2)  = chgneg

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, local_communicator, mpi_error_code)

       
       chgpos = mpitotoutdp(1)  
       chgneg = mpitotoutdp(2)  

       write(6,*) 'CHG SUM' , chgpos,chgneg
!#endif

!     otherwise X=Y=0 --> disch=0. so if enough + charge alone for elecmag > ecrit will loop infinitely because scdich =0 and no changes will occur on sctot and elec.
      if (chgpos.lt.1.) chgpos=1.
      if (chgneg.gt.-1.) chgneg=-1.
     

!      what this does: if amount of +/- charge are comparable, then just take 30% of max . IF much more of one polarity, will take just the min of (chgpos, Abs(chgneg))


      chg = elgt1*Max(chgpos, Abs(chgneg))
!      chg = (1.0 - elgt1)*Max(chgpos, Abs(chgneg))
      IF ( chg .gt. chgpos ) chg = chgpos
      IF ( chg .gt. Abs(chgneg) ) chg = Abs(chgneg)
     

      if (abs(chgneg).gt.0.) then
      x = chg/(Abs(chgneg)) 
      else
      x=0.
      endif

      if (abs(chgpos).gt.0.) then 
      y = chg/(Abs(chgpos))
      else
      y=0.
      endif
      

      write(6,*) 'X,Y', x,y,chg! ,chgpos,chgneg

       totneg1 = 0.0
       totpos1 = 0.0
       totneg2 = 0.0
       totpos2 = 0.0

      DO jy = 1,jye ! ny
      do kz = 1,nz
      DO ix = 1,ixe ! nx

      IF ( ij_reduce(ix,jy) > 0 ) THEN

!      NEED THIS IN C/kg and have C/m3 --> divide by RHO

!      chg = chg*(1.0 - elgt1)
!      IF ( t1(ix,jy,kz) .lt. 0.0 ) THEN ! local chg is neg. so add pos. ions
      IF ( sctot(ix,kz,jy) .lt. -chgthr ) THEN ! local chg is neg. so add pos. ions
!        an(ix,jy,kz,lscpi) = an(ix,jy,kz,lscpi) - t1(ix,jy,kz)*x
        scdisch(ix,jy,kz) = (-(sctot(ix,kz,jy)+chgthr)*y )/ db(ix,kz,jy)!x ! OPPOSITE TO DECREASE IT !
!      ELSEIF ( t1(ix,jy,kz) .gt. 0.0 ) THEN ! local chg is pos. so add neg. ions
      ELSEIF ( sctot(ix,kz,jy) .gt. chgthr ) THEN ! local chg is pos. so add neg. ions
!        an(ix,jy,kz,lscni) = an(ix,jy,kz,lscni) + t1(ix,jy,kz)*y
        scdisch(ix,jy,kz) = ( (sctot(ix,kz,jy) - chgthr)*x )/db(ix,kz,jy) !y
      ENDIF
      
       if (abs(scdisch(ix,jy,kz)).gt.1.e-3) write(6,*) 'scdisch,sctot', scdisch(ix,jy,kz),sctot(ix,kz,jy)     


      ENDIF

!     compute area fraction Ax/sum(Ax). where Ax=(Qx*rhoair/(Nx*rhox))**2/3 * Nx 
!   but qnr etc. is N/rhoair, and need to convert scdisch to C/kg
!   So have Ax=(Qx*rhoair/(Nx*rhoair*rhox))**2/3 * (Nx*rhoair) = (Qx/(Nx*rhox))**2/3 * (Nx*rhoair)
!  or Ax = rhoair*(Qx/(rhox))**2/3 * (Nx)**(1/3), and the factor of rhoair cancels out in the ratio Ax/Atot
    
      if(qnr(ix,kz,jy).lt.0.) qnr(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      if(qr(ix,kz,jy).lt.0.) qr(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      if(qng(ix,kz,jy).lt.0.) qng(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      if(qg(ix,kz,jy).lt.0.) qg(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      if(qni(ix,kz,jy).lt.0.) qni(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      if(qi(ix,kz,jy).lt.0.) qi(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      if(qns(ix,kz,jy).lt.0.) qns(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      if(qs(ix,kz,jy).lt.0.) qs(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      if(qnc(ix,kz,jy).lt.0.) qnc(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      if(qc(ix,kz,jy).lt.0.) qc(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      if (present(qh)) then
      if(qnh(ix,kz,jy).lt.0.) qnh(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      if(qh(ix,kz,jy).lt.0.) qh(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      endif

      fracs = 0.0
      fraci = 0.0
      frach = 0.0
      fracw = 0.0
      fracr = 0.0
      if (qs(ix,kz,jy).gt.qmin) fracs=((qs(ix,kz,jy)/rhos)**(2./3.) ) * qns(ix,kz,jy)**(1./3.)
      if (qi(ix,kz,jy).gt.qmin) fraci=((qi(ix,kz,jy)/rhoi)**(2./3.) ) * qni(ix,kz,jy)**(1./3.)
      if (qh(ix,kz,jy).gt.qmin) frach=((qg(ix,kz,jy)/rhog)**(2./3.) ) * qng(ix,kz,jy)**(1./3.)
      if (qc(ix,kz,jy).gt.qmin) fracw=((qc(ix,kz,jy)/rhow)**(2./3.) ) * qnc(ix,kz,jy)**(1./3.)
      if (qr(ix,kz,jy).gt.qmin) fracr=((qr(ix,kz,jy)/rhor)**(2./3.) ) * qnr(ix,kz,jy)**(1./3.)

      deninv = 1. ! /db(ix,kz,jy)
      fractot=fracs+fraci+frach &
     &             +fracr+fracw

      if (present(qh)) frachl=((qh(ix,kz,jy)/rhoh)**2./3. ) * qnh(ix,kz,jy)**(1./3.)

      if (present(qh))  fractot=fractot+frachl
     
!'     reduce/redistribute scdisch on scx based on area fraction 
   
      IF ( Abs(sctot(ix,kz,jy)) > chgthr) THEN
      IF ( fractot > 0. ) THEN
      
      if (qs(ix,kz,jy).gt.qmin) then ! avoid division by zero if fractot = 0 because qx=0
      if (sctot(ix,kz,jy) .gt.  chgthr) scs(ix,kz,jy)=scs(ix,kz,jy)-scdisch(ix,jy,kz)*(fracs/fractot)*deninv
      if (sctot(ix,kz,jy) .lt. -chgthr) scs(ix,kz,jy)=scs(ix,kz,jy)+scdisch(ix,jy,kz)*(fracs/fractot)*deninv
      endif

      if (qi(ix,kz,jy).gt.qmin) then
      if (sctot(ix,kz,jy) .gt.  chgthr) sci(ix,kz,jy)=sci(ix,kz,jy)-scdisch(ix,jy,kz)*(fraci/fractot)*deninv
      if (sctot(ix,kz,jy) .lt. -chgthr) sci(ix,kz,jy)=sci(ix,kz,jy)+scdisch(ix,jy,kz)*(fraci/fractot)*deninv
      endif

      if (qr(ix,kz,jy).gt.qmin) then
      if (sctot(ix,kz,jy) .gt.  chgthr) scr(ix,kz,jy)=scr(ix,kz,jy)-scdisch(ix,jy,kz)*(fracr/fractot)*deninv
      if (sctot(ix,kz,jy) .lt. -chgthr) scr(ix,kz,jy)=scr(ix,kz,jy)+scdisch(ix,jy,kz)*(fracr/fractot)*deninv
      endif

      if (qg(ix,kz,jy).gt.qmin) then
      if (sctot(ix,kz,jy) .gt.  chgthr) sch(ix,kz,jy)=sch(ix,kz,jy)-scdisch(ix,jy,kz)*(frach/fractot)*deninv
      if (sctot(ix,kz,jy) .lt. -chgthr) sch(ix,kz,jy)=sch(ix,kz,jy)+scdisch(ix,jy,kz)*(frach/fractot)*deninv
      endif

      if (qc(ix,kz,jy).gt.qmin) then
      if (sctot(ix,kz,jy) .gt.  chgthr) scw(ix,kz,jy)=scw(ix,kz,jy)-scdisch(ix,jy,kz)*(fracw/fractot)*deninv
      if (sctot(ix,kz,jy) .lt. -chgthr) scw(ix,kz,jy)=scw(ix,kz,jy)+scdisch(ix,jy,kz)*(fracw/fractot)*deninv
      endif
   
      if (present(qh)) then
      if (qh(ix,kz,jy).gt.qmin) then
      if (sctot(ix,kz,jy) .gt.  chgthr) schl(ix,kz,jy)=schl(ix,kz,jy)-scdisch(ix,jy,kz)*(frachl/fractot)*deninv
      if (sctot(ix,kz,jy) .lt. -chgthr) schl(ix,kz,jy)=schl(ix,kz,jy)+scdisch(ix,jy,kz)*(frachl/fractot)*deninv
      endif
      endif
      
      ELSE

      if (sctot(ix,kz,jy) > 0.0) then
        sciona(ix,kz,jy) = sciona(ix,kz,jy) - scdisch(ix,jy,kz)*deninv
!        chgrtn = chgrtn + scdisch(ix,jy,kz)*dv
      elseif (sctot(ix,kz,jy) <= 0.0) THEN 
        sciona(ix,kz,jy) = sciona(ix,kz,jy) + scdisch(ix,jy,kz)*deninv
!        chgrtp = chgrtp + scdisch(ix,jy,kz)*dv
      endif
      
      ENDIF

!     UPDATE sctot     

       totneg1 = totneg1 + Min(0.0, sctot(ix,kz,jy)*dv)
       totpos1 = totpos1 + Max(0.0, sctot(ix,kz,jy)*dv)
       

       sctot(ix,kz,jy) = db(ix,kz,jy)*(sch(ix,kz,jy)+scr(ix,kz,jy)+sci(ix,kz,jy)+scs(ix,kz,jy)+scw(ix,kz,jy)+sciona(ix,kz,jy))
       if (present(qh)) sctot(ix,kz,jy)=sctot(ix,kz,jy)+db(ix,kz,jy)*schl(ix,kz,jy)

       totneg2 = totneg2 + Min(0.0, sctot(ix,kz,jy)*dv)
       totpos2 = totpos2 + Max(0.0, sctot(ix,kz,jy)*dv)
     
!      if (sctot(ix,kz,jy).lt.0.) write(0,*) 'IN DISCHARGE', sctot(ix,kz,jy),sch(ix,kz,jy),scr(ix,kz,jy),sci(ix,kz,jy),scs(ix,kz,jy),scw(ix,kz,jy),ix,jy,kz
!      if (sctot(ix,kz,jy).lt.0.) write(0,*) 'IN DISCHARGE', scw(ix,kz,jy),scdisch(ix,jy,kz),fracw,fractot,qc(ix,kz,jy),rhow,qnc(ix,kz,jy),fracs,fraci,frach,fracr,qr(ix,kz,jy),rhor,qnr(ix,kz,jy)

      ENDIF ! Abs(sctot(ix,kz,jy)) > chgthr)
      
      end do
      end do
      end do

#ifdef MPI
       mpitotindp(1)  = totpos1
       mpitotindp(2)  = totneg1
       mpitotindp(3)  = totpos2
       mpitotindp(4)  = totneg2

      CALL MPI_Reduce(mpitotindp, mpitotoutdp, 4, MPI_DOUBLE_PRECISION, MPI_SUM, 0, local_communicator, mpi_error_code)


      IF ( mytask == 0 ) THEN ! only print from task 0

       totpos1 = mpitotoutdp(1)  
       totneg1 = mpitotoutdp(2)  
       totpos2 = mpitotoutdp(3)  
       totneg2 = mpitotoutdp(4)  

      ENDIF
#endif

      IF ( mytask == 0 ) THEN ! only print from task 0
        write(iunit,'(a,2(2x,1pe12.5))') 'pos/neg1: ',totpos1,totneg1,totpos1+totneg1
        write(iunit,'(a,2(2x,1pe12.5))') 'pos/neg2: ',totpos2,totneg2,totpos2+totneg2
      ENDIF


!        write(6,*) '+============AFTER DISCHARGE============+'
!        write(6,*) 'SCTOT MIN', MINVAL(sctot(:,:,:))
!        write(6,*) 'SCTOT MAX', MAXVAL(sctot(:,:,:))
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
!      IF ( mytask == 0 ) THEN
      write(iunit,*) '______________________________'
      write(iunit,*) 'LIGHTNING DISCHARGE:  SCHEME 1'
!      write(iunit,*) 'TIME STEP IS: ',time_real
      write(iunit,*) 'TOTAL NUMBER OF DISCHARGES IS: ',nlgt1
      write(iunit,*) 'chgpos,chgneg : ',chgpos,chgneg
      write(iunit,*) '______________________________'
!      ENDIF
!

      ENDIF ! etot .gt elght1 !!!!
!
      nlgt1 =0      

!      write(iunit,*) 'END OF LIGHTNING DISCHARGE'

!    END OF LIGHTNING PARAMETERIZATION VERSION 1

     ELSEIF ( ilight >= 2 ) THEN
     
      if ( iliter .eq. 1 ) then
      lgtstp = 0
      end if
      
      scthn = -scth
      scthp =  scth

#ifdef MPI
      my_comm = local_communicator
      nz1 = nz-kd1
      
      IF ( .not. allocated( mpimxmn ) ) THEN
      allocate(mpimxmn(nz1,5))
      allocate(max_send_buffer(2,nz1))
      allocate(max_recv_buffer(2,nz1))
      allocate(min_send_buffer(2,nz1))
      allocate(min_recv_buffer(2,nz1))

      allocate(mpiloc(nz1,2))
      allocate(mpiloc0(nz1,2,0:ntasks-1))
      allocate(loc_send_buffer(nz1,2))
      allocate(loc_recv_buffer(nz1*ntasks,2))

!      real,allocatable :: mpiexyz(:,:), mpiexyz0(:,:,:)
!      real,allocatable :: exyz_send_buffer(:,:),exyz_recv_buffer(:,:)

      allocate(mpiexyz(nz1,3))
      allocate(mpiexyz0(nz1,3,0:ntasks-1))
      allocate(exyz_send_buffer(nz1,3))
      allocate(exyz_recv_buffer(nz1*ntasks,3))
      ENDIF

#endif
     
!
!  find largest emag, but first zero emag
!
      if ( ndebug > 2 .and. my_rank == 0 ) write(iunit,*) 'emax'
      emax = 0.0
      iemaxx = 1
      iemaxy = 1
      iemaxz = 1
!
!  find largest emag as a f(z) and locations, 
!  but first zero emag and locations
!
      ixb = 1
      ixe = nx
      
!      IF ( myproci == 1 ) ixb = 2
!      IF ( myproci == nproci ) ixe = nx - 1
      
      jyb = 1
      jye = ny
      
!      IF ( myprocj == 1 ) jyb = 2
!      IF ( myprocj == nprocj ) jye = ny - 1
      IF ( ny .le. 2 ) THEN
        jyb = 1
        jye = 1
      ENDIF

      ninit1 = 0
      
      do kz = 1,nz-kd1
!
      emaxkz(kz) = 0
      emaxlightkz(kz) = -1.
      nemaxx(kz) = 1
      nemaxy(kz) = 1
      nemaxz(kz) = kz
!
      do jy = jyb,jye ! Min(ny-1,2),Max(1,ny-2*jd1)
      do ix = ixb,ixe ! 2,nx-2*id1
!
      if ( elec(ix,kz,jy,4) .gt. emaxkz(kz) ) then
      nemaxx(kz) = ix
      nemaxy(kz) = jy
      nemaxz(kz) = kz
      emaxkz(kz) = elec(ix,kz,jy,4)
      exmaxkz(kz) = elec(ix,kz,jy,1)
      eymaxkz(kz) = elec(ix,kz,jy,2)
      ezmaxkz(kz) = elec(ix,kz,jy,3)
      IF ( emaxkz(kz) >= ebrkd(ix,jy,kz) ) THEN
        emaxlightkz(kz) = emaxkz(kz)
        ninit1 = ninit1 + 1
      ENDIF
      end if
!
      if ( elec(ix,kz,jy,4) .gt. emax ) then
      iemaxx = ix
      iemaxy = jy
      iemaxz = kz
      emax = elec(ix,kz,jy,4)
      end if
!
      end do
      end do


#ifdef MPI
      mpimxmn(kz,1) = exmaxkz(kz) ! elec(iex)%flt3d(nemaxx(kz),nemaxy(kz),kz)
      mpimxmn(kz,2) = eymaxkz(kz) ! elec(iey)%flt3d(nemaxx(kz),nemaxy(kz),kz)
      mpimxmn(kz,3) = ezmaxkz(kz) ! elec(iez)%flt3d(nemaxx(kz),nemaxy(kz),kz)
      mpimxmn(kz,4) = emaxkz(kz)
!      mpimxmn(kz,5) = emaxlightkz(kz)
      mpiloc(kz,1) = nemaxx(kz)
      mpiloc(kz,2) = nemaxy(kz)
!      mpiloc(kz,3) = nemaxz(kz)

#endif
      end do
!

#ifdef MPI

      IF ( ntasks > 1 ) THEN
      DO kz = 1,nz1

       exyz_send_buffer(kz,1) = emaxlightkz(kz)
!       IF ( emaxlightkz(kz) > 0.0 ) THEN
!         write(iunit,*) 'kz, emaxlightkz,ezmaxkz = ',kz,emaxlightkz(kz),ezmaxkz(kz)
!       ENDIF

      ENDDO
      mxmn_size = nz1

! find out the max emag at each level, and which processor it is on
      CALL MPI_Allreduce(exyz_send_buffer, exyz_recv_buffer, mxmn_size, MPI_REAL,  &
     &                  MPI_MAX, local_communicator, mpi_error_code)

      IF (mpi_error_code /= MPI_SUCCESS) THEN

        WRITE(0,*) my_rank, "mudlight: about to CALL MPI_Abort after Allreduce attempt"

        CALL MPI_Abort(local_communicator, mpi_error_code,mpi_abort_error_code)
      ENDIF

! put global maximum emax into mpimxmn
      DO kz = 1,nz1

       emaxlightkz(kz) = exyz_recv_buffer(kz,1)
!       write(iunit,*) 'kz, emaxlightkz,ezmaxkz = ',kz,emaxlightkz(kz),ezmaxkz(kz)

      ENDDO



      DO kz = 1,nz1

       max_send_buffer(1,kz) = mpimxmn(kz,4)
!       max_send_buffer(3,kz) = mpimxmn(kz,5)
       max_send_buffer(2,kz) = real(my_rank)

      ENDDO
    
      mxmn_size = nz1

! find out the max emag at each level, and which processor it is on
      CALL MPI_Allreduce(max_send_buffer, max_recv_buffer, mxmn_size, MPI_2REAL,  &
     &                  MPI_MAXLOC, local_communicator, mpi_error_code)

      IF (mpi_error_code /= MPI_SUCCESS) THEN

        WRITE(0,*) my_rank, "mudlight: about to CALL MPI_Abort after Allreduce attempt"

        CALL MPI_Abort(local_communicator, mpi_error_code,mpi_abort_error_code)
      ENDIF

! put global maximum emax into mpimxmn
      emax = 0.0
      DO kz = 1,nz1

       mpimxmn(kz,4) = max_recv_buffer(1,kz)
       emaxkz(kz) = max_recv_buffer(1,kz)

      ENDDO

! send e-field components of emag max
      DO kz = 1,nz1
       exyz_send_buffer(kz,1) = mpimxmn(kz,1)
       exyz_send_buffer(kz,2) = mpimxmn(kz,2)
       exyz_send_buffer(kz,3) = mpimxmn(kz,3)
      ENDDO

      exyz_send_size = nz1*3
      exyz_recv_size = nz1*3

      CALL MPI_Gather(exyz_send_buffer,exyz_send_size,MPI_REAL,  &
     &               exyz_recv_buffer,exyz_recv_size,MPI_REAL,   &
     &               0,local_communicator,mpi_error_code)

      IF (mpi_error_code /= MPI_SUCCESS) THEN

        WRITE(0,*) my_rank, "mudlight: about to CALL MPI_Abort after Gather attempt"

        CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,mpi_abort_error_code)

       END IF !! (debug_mpi)

       IF (my_rank == 0) THEN

        mpiexyz0 = reshape(exyz_recv_buffer,(/nz1,3,number_of_processes/))

           DO kz = 1,nz1
       
            irank1 = int(max_recv_buffer(2,kz))

            exmaxkz(kz) = mpiexyz0(kz, 1, irank1)
            eymaxkz(kz) = mpiexyz0(kz, 2, irank1)
            ezmaxkz(kz) = mpiexyz0(kz, 3, irank1)

          ENDDO
        ENDIF


! send locations (x,y) of emag max
      DO kz = 1,nz1
       loc_send_buffer(kz,1) = ixbeg-1+mpiloc(kz,1)
       loc_send_buffer(kz,2) = jybeg-1+mpiloc(kz,2)
      ENDDO

      loc_send_size = nz1*2
      loc_recv_size = nz1*2

      CALL MPI_Gather(loc_send_buffer,loc_send_size,MPI_INTEGER, &
     &               loc_recv_buffer,loc_recv_size,MPI_INTEGER,  &
     &               0,local_communicator,mpi_error_code)

      IF (mpi_error_code /= MPI_SUCCESS) THEN

        WRITE(0,*) my_rank, "mudlight: about to CALL MPI_Abort after Gather attempt"

        CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,mpi_abort_error_code)

       END IF !! (debug_mpi)

       IF (my_rank == 0) THEN

        mpiloc0 = reshape(loc_recv_buffer,(/nz1,2,number_of_processes/))

!          DO itype = 1,2 ! type is really the ix or jy index
           DO kz = 1,nz1
       
            irank1 = int(max_recv_buffer(2,kz))

            mpiloc(kz,1) = mpiloc0(kz,1,irank1)
            mpiloc(kz,2) = mpiloc0(kz,2,irank1)
            nemaxx(kz) = mpiloc0(kz,1,irank1)
            nemaxy(kz) = mpiloc0(kz,2,irank1)

          ENDDO
!         ENDDO

        DO kz = nz1,1,-1
!         write(iunit,'(a,3(i4,1x),2(g15.3,1x))') 
!     :  'IEMAXX, IEMAXY, IEMAXZ, EMAX =: ',mpimxmn(kz,4),
!     :  iemaxy, iemaxz, emax, elec(iez)%flt3d(iemaxx,iemaxy,iemaxz)

!         write(iunit,'(a,3(i4,1x),3(g15.5,1x))') 
!     :  'gIEMAXX, gIEMAXY, gIEMAXZ, gEMAX =: ',
!     :  mpiloc(kz,1) , mpiloc(kz,2) , kz, exmaxkz(kz),eymaxkz(kz),ezmaxkz(kz) !, elec(iez)%flt3d(iemaxx,iemaxy,iemaxz)
          
        ENDDO


      ENDIF !! my_rank==0
      
      ENDIF ! ntasks

#endif

      DO kz = 1,nz-1
      if ( emaxkz(kz) .gt. emax ) then
      iemaxx = nemaxx(kz)
      iemaxy = nemaxy(kz)
      iemaxz = kz
      emax = emaxkz(kz) 
      end if
      ENDDO

      IF ( my_rank == 0 ) THEN
      write(iunit,*) 'Location and value of max(|E|)'
!      write(iunit,'(a,3(i4,1x),2(g15.3,1x))') 
!     :  'IEMAXX, IEMAXY, IEMAXZ, EMAX =: ',iemaxx,
!     :  iemaxy, iemaxz, emax, elec(iez)%flt3d(iemaxx,iemaxy,iemaxz)

      write(iunit,'(a,3(i4,1x),2(g15.3,1x))')  &
     &  'IEMAXX, IEMAXY, IEMAXZ, EMAX =: ',iemaxx, &
     &  iemaxy, iemaxz, emax, emaxkz(iemaxz) ! elec(iez)%flt3d(iemaxx,iemaxy,iemaxz)
!
      write(iunit,'(a)') 'Z, NEMAXX, NEMAXY, NEMAXZ, EMAXKZ, EZ, EH, yesno =: '
     
       ENDIF
     
!
!  check if lightning should occur (if emax .lt. elgtthx )
!  there are two methods....One based on a typical value
!  for the whole atmosphere....another based on a breakdown
!  efield that is a function of z
!
!  test if maximum breakdown approach is going to allow lightning
!  note that here lgtstp is the important parameter
!
      if ( ibrkd .eq. 0 ) then
      lgtstp = 0
      if ( emax .lt. elgtthx ) then
      IF ( my_rank == 0 ) write(iunit,*) 'NO LIGHTNING AT STEP',nstep
      yesno = 'NO'
      lgtstp = 1
      go to 2999
      else
      IF ( my_rank == 0 ) write(iunit,*) 'START LIGHTNING AT STEP',nstep
      lgtstp = 0
      end if
      end if
!
!  test if variable breakdown approach is going to allow lightning
!  note that here lgtgo is the important parameter
!
      if ( ibrkd .ge. 2 ) then
      
      lgtstp = 0
      lgtgo = 0
      do kz = nz-kd1,1,-1
      if ( emaxlightkz(kz) <= 0. ) then
!      write(iunit,'(a,i4,1x,f7.0)') 
!     :       'NO LIGHTNING AT STEP AND Z OF',nstep,zlev(kz)
      yesno = 'NO LIGHT'
      lgtstp = 1
      else
      kzbrkd = kz
!      write(iunit,'(a,i4,1x,f7.0,1x,i4)') 
!     >   'START LIGHTNING AT STEP AND Z OF', nstep,zlev(kz),kzbrkd
      yesno = 'START LIGHT'
      lgtgo = 1
      end if
      IF ( my_rank == 0 .and. (yesno == 'START LIGHT' .or. ndebug > 1 ) ) THEN
      write(iunit,'(3(i4,1x),1x,3(1x,1pe12.4),1x,a,1x,0pf7.0)') &
     &   nemaxx(kz),nemaxy(kz), nemaxz(kz),                     &
     &   emaxkz(kz), ezmaxkz(kz),                               &
     &  Sqrt(exmaxkz(kz)**2 + eymaxkz(kz)**2 ),yesno,           &
     &   zlev(kz,1,1)
!      write(iunit,'(3i4,2x,3(1x,1pe12.4),1x,a,1x,0pf7.0)') 
!     >   nemaxx(kz),nemaxy(kz), nemaxz(kz),
!     >   emaxkz(kz), elec(iez)%flt3d(nemaxx(kz),nemaxy(kz),nemaxz(kz)),
!     >  Sqrt(elec(iex)%flt3d(nemaxx(kz),nemaxy(kz),nemaxz(kz))**2 +
!     >       elec(iey)%flt3d(nemaxx(kz),nemaxy(kz),nemaxz(kz))**2 ),yesno,
!     >   zlev(kz)
      
      ENDIF
      end do
      if ( lgtgo .eq. 0 ) go to 2999 
      IF (lgtgo .eq. 1) lgtstp = 0
      
      end if
      
!      IF ( isa .eq. 1 ) THEN
!        nliter2 = nliter
!      ELSE
        nliter2 = 1
!      ENDIF
      
      IF ( ny .le. 2 ) THEN
!        lgtstp = 1
!        RETURN
      ENDIF
     
      
! make list of possible initiation points
      ninit = 0
      DO kz = 2,nz-1
      DO jy = jyb,jye
      DO ix = ixb,ixe
       IF ( elec(ix,kz,jy,4) >= ebrkd(ix,jy,kz)*elgtfdel ) THEN
       IF ( lightmask(ix,jy) == 0 .and. initmask(ix,jy,kz) == 0 ) THEN
        ninit = ninit + 1
        lightinit(ix,jy) = 1
        IF (ninit .gt. npick) THEN
          write(6,'(a)') 'Ee-gad! Over-ran the ipick array! ABORT!'
          STOP
        END IF
        ipick(ninit,1) = ix
        ipick(ninit,2) = jy
        ipick(ninit,3) = kz
       ENDIF
       END IF ! 
      ENDDO  ! i loop
      ENDDO  ! j loop
      ENDDO  ! k loop

#ifdef MPI
! send e-field components of emag max
      nz1 = 2
        init_send_buffer(1) = my_rank
        init_send_buffer(2) = ninit
!        init_send_buffer(3) = ninit1
        
      IF ( ndebug > 1 ) write(iunit,*) 'my_rank,ninit = ',my_rank,ninit,ninit1

      init_send_size = nz1
      init_recv_size = nz1

      CALL MPI_Gather(init_send_buffer,init_send_size,MPI_INTEGER,  &
     &               init_recv_buffer,init_recv_size,MPI_INTEGER,   &
     &               0,local_communicator,mpi_error_code)

       IF (mpi_error_code /= MPI_SUCCESS) THEN

        WRITE(0,*) my_rank, "discharge_msz: about to CALL MPI_Abort after Gather attempt"

        CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,mpi_abort_error_code)

       END IF !! (debug_mpi)

       IF (my_rank == 0) THEN

        mpi_init0 = reshape(init_recv_buffer,(/2,number_of_processes/))

           ninittot = 0
           DO kz = 0,number_of_processes-1
       
            ninittot = ninittot +  mpi_init0(2, kz)

          ENDDO
          
       IF ( ndebug > 1)  write(iunit,'(1x,i6,a,i6)') ninittot, ' possible starting points' ! on rank =',my_rank
          

       ENDIF

       mpitotinint(1)  = ninit

      CALL MPI_AllReduce(mpitotinint, mpitotoutint, 1, MPI_INTEGER, MPI_SUM, local_communicator, mpi_error_code)

       
!      IF ( my_rank == 0 ) THEN
       ninittot = mpitotoutint(1)  
          

#else

       ninittot = ninit

#endif

      IF ( my_rank == 0 .and. ndebug >= 1 ) write(iunit,'(1x,i6,a)') ninittot, ' possible starting points'
!
         IF ( ninittot == 0 ) THEN
           lgtstp = 1
           loccur = 0
           GOTO 2999
         ENDIF
!
! put continue here to make a new channel if previous one turned out 
! to be a dud.  Maybe implement a list of starting points....
!
      iretry = 0
      chgmin = (1.0 - fscth*(niter-1))*scth  !niter=2

 1101 CONTINUE
       t2(:,:,:) = 0.
       t4(:,:,:) = 0.
       IF ( f_clnox ) THEN
       lnox(:,:,:) = 0. 
       nox(:,:,:) = 0. 
       ENDIF
        
      iretry = iretry+1
!      IF (iretry .gt. nretry) THEN
!         write(6,'(a)') 'Too many channel tries...abort.'
!           lgtstp = 1
!           loccur = 0
!         go to 2999
!      END IF
            
      if ( ndebug >= 1 .and. my_rank == 0 ) write(*,*)  'find lightning initiation pt. Try # ',iretry 
      IF (ninittot .eq. 0) THEN
        IF ( my_rank == 0 ) write(6, '(a)') ' No more starting points, try next timestep'
        lgtstp = 1
!        loccur = 0 ! do not set this 
        GOTO 2999
      END IF

      ichoose = int(real(ninittot)*ran0mpi(iseed)) + 1
      IF ( my_rank == 0 ) THEN
      write (6,'(a,i4)') 'ichoose = ',ichoose
!      write(6,'(a,3(1x,i3))') 'Picked point ',islgt,jslgt,kslgt
      ENDIF

      ichooserank = 0
#ifdef MPI
! find which processor has the initiation point
      IF ( my_rank == 0 ) THEN
       ninittmp = 0
           ichooserank = -1
           DO m = 0,number_of_processes-1
       
!            write(iunit,*) 'm = ',m
            
            IF ( ichoose .le. ninittmp + mpi_init0(2, m) ) THEN
               IF ( ndebug > 1 ) write(iunit,*) 'ichooserank = ',mpi_init0(1,m)
               ichooserank = mpi_init0(1,m)
               ichoosenumber = ichoose - ninittmp
               EXIT
!            ELSE
!              write(iunit,*) 'init point not on processor ',mpi_init0(1,m)
            ENDIF
            
            ninittmp = ninittmp + mpi_init0(2, m)

          ENDDO
          
          IF ( ichooserank == -1 ) THEN
            write(iunit,*) 'Problem choosing init point! ichooserank = ',ichooserank
          ENDIF

          ichoose_send_buffer(1) = ichooserank
          ichoose_send_buffer(2) = ichoosenumber
          
      ELSE
          ichoose_send_buffer(1) = -1
          ichoose_send_buffer(2) = -1
      ENDIF ! (my_rank == 0 )

       
       CALL MPI_Bcast(ichoose_send_buffer, 2, MPI_INTEGER, 0, local_communicator, mpi_error_code)
       ichooserank = ichoose_send_buffer(1)
       ichoose = ichoose_send_buffer(2)
       
       IF ( ndebug > 1 .and. my_rank == 0 )  write(iunit,*) 'my_rank: ichooserank,ichoose = ',my_rank,ichooserank,ichoose
#endif

       
      IF ( my_rank == ichooserank ) THEN
      
      i1 = ipick(ichoose,1) ! + Igslg - 1
      j1 = ipick(ichoose,2) ! + Jgslg - 1
      k1 = ipick(ichoose,3)
      ixst = i1
      jyst = j1
      kzst = k1
      ixst0 = i1
      jyst0 = j1
      kzst0 = k1
      initmask(i1,j1,k1) = 1 ! mark this point as having been used for this time step
!      write(iunit,'(a,3(1x,i3))') 'Picked starting point ',i1,j1,k1
! update the points array (move up the points lower on the list)     
      DO idum=ichoose,ninit-1
        ipick(idum,1) = ipick(idum+1,1)
        ipick(idum,2) = ipick(idum+1,2)
        ipick(idum,3) = ipick(idum+1,3)
      END DO

          ichoose_send_buffer(1) = i1
          ichoose_send_buffer(2) = j1
          ichoose_send_buffer(3) = k1
          ichoose_send_buffer(4) = igslg0
          ichoose_send_buffer(5) = jgslg0
          ichoose_send_buffer(6) = kgslg0
          ichoose_send_buffer(7) = Nx
          ichoose_send_buffer(8) = Ny
          ichoose_send_buffer(9) = Nz
      
          mpitotin(1) = elec(i1,k1,j1,4) ! emag
          mpitotin(2) = elec(i1,k1,j1,3) ! ez
          mpitotin(3) = pot(i1,k1,j1) ! ez
          mpitotin(4) = zlev(k1,i1,j1)
          
!          write(iunit,*) 'EMAG,k-1,k,k+1: ',elec(iemag)%flt3d(i1,j1,k1-1),elec(iemag)%flt3d(i1,j1,k1),elec(iemag)%flt3d(i1,j1,k1+1)
!          write(iunit,*) 'EZ,k-1,k,k+1: ',elec(iez)%flt3d(i1,j1,k1-1),elec(iez)%flt3d(i1,j1,k1),elec(iez)%flt3d(i1,j1,k1+1)
!          write(iunit,*) 'x,y,z = ',gxt(i1,1),gyt(j1,1),gzt(k1,1)
      
       ninit = ninit - 1
      ENDIF ! my_rank

#ifdef MPI
       CALL MPI_Bcast(ichoose_send_buffer, 9, MPI_INTEGER, ichooserank, local_communicator, mpi_error_code)
#endif
       i1 = ichoose_send_buffer(1)
       j1 = ichoose_send_buffer(2)
       k1 = ichoose_send_buffer(3)
       ixst = i1 + ichoose_send_buffer(4)
       jyst = j1 + ichoose_send_buffer(5)
       kzst = k1 + ichoose_send_buffer(6)
       ioffset = ichoose_send_buffer(4)
       joffset = ichoose_send_buffer(5)
       koffset = ichoose_send_buffer(6)
       NLxchoose = ichoose_send_buffer(7)
       NLychoose = ichoose_send_buffer(8)
       NLzchoose = ichoose_send_buffer(9)

#ifdef MPI
       CALL MPI_Bcast(mpitotin, 4, MPI_REAL, ichooserank, local_communicator, mpi_error_code)
#endif
        emag = mpitotin(1)
        ez = mpitotin(2)
        pot1 = mpitotin(3)
        zinit = mpitotin(4)

      IF ( my_rank == 0 ) THEN
      write(iunit,'(a,3(1x,i3))') 'Picked starting point ',i1+ioffset,j1+joffset,k1+koffset
      ENDIF


      islgt = ixst
      jslgt = jyst
      kslgt = kzst

      islgt_local = i1
      jslgt_local = j1
      kslgt_local = k1

! one fewer init point now      
      ninittot = ninittot - 1
! find neighboring point with largest E and make it the first bond.

    IF ( my_rank == 0 ) THEN
      ! write(6,*) 'LIGHTNING INITIATION LOCATION AND |E|'
      write(6,*) 'NSTEP,ILITER = ', nstep,iliter
      write(iunit,*) 'ISLGT, JSLGT, KSLGT : ',islgt,jslgt, kslgt
      write(iunit,*) 'emag,pot,ez,zc = ',emag,pot1,ez,zinit ! z1d2(kslgt_local,1)
    ENDIF



!
! For now, the lightning channel will stay vertical to avoid making a streamline across an MPI patch boundary
!
!
!  draw a streamline along the efield lines until the efield threshold
!  is found and store the potential value for the + and - regions
!
      if ( ndebug > 1 .and. my_rank == 0 )write(*,*) 'generate +/- channels'
      dt = 4.0 ! dtp
      times = dtp*(nstep-1)
!
!  Initialize lstop--indicates trajectory end.
!
      trje(1:2,0,:) = 0
      trje(1,1,lstop) = 0.
      trje(2,1,lstop) = 0.
!
!  initiate position of trajectories
!
      if ( ndebug > 2 ) write(0,*) 'check 2: islgt_local,jslgt_local,dx,dy = ',  &
     &       islgt_local,jslgt_local,kslgt,dx,dy,zlev(kslgt,islgt_local,jslgt_local)
      IF ( my_rank == ichooserank ) THEN !{
      trje(1,1,llx) = (islgt-1)*dx
      trje(2,1,llx) = (islgt-1)*dx
!      if ( ndebug > 2 ) write(0,*) 'check 2a'
      trje(1,1,lly) = (jslgt-1)*dy
      trje(2,1,lly) = (jslgt-1)*dy
!      if ( ndebug > 2 ) write(0,*) 'check 2b'
      trje(1,1,llz) = zlev(kslgt,islgt_local,jslgt_local) ! z1d2(kslgt,1) ! (kslgt-1)*dz
      trje(2,1,llz) = zlev(kslgt,islgt_local,jslgt_local) ! z1d2(kslgt,1) ! (kslgt-1)*dz
!      trje(1,1,llz) = (kslgt-1)*dz
!      trje(2,1,llz) = (kslgt-1)*dz
      
      if ( ndebug > 2 ) THEN 
        write(0,*) 'check 2c',elgtfestop,elgtfestopcg
        write(0,*) 'ebrkd: ',ebrkd(islgt_local,jslgt_local,kslgt)
      ENDIF
        
      DO kz=1,nz
!        estop(kz) =  Min( (elgtthn)*ebrkd(islgt,jslgt,kz)/elgtthx,elgtthn)
!        if ( ndebug > 2 ) write(0,*) 'check 2d kz=',kz
        estop(kz)   =  elgtfestop*ebrkd(islgt_local,jslgt_local,kz)
        estopcg(kz) =  elgtfestopcg*ebrkd(islgt_local,jslgt_local,kz)
      ENDDO
      
      ENDIF !}
!
!  Set nttimend
!
      if ( ndebug > 2 ) write(0,*) 'check 3'
      nttimend(1) = 0
      nttimend(2) = 0
      nttimendcg(1) = 0
      nttimendcg(2) = 0
!
! set end indicators to zero for safety      
      ichend(1) = 0
      ichend(2) = 0
!      
! zero out the number of points added in 'channel cone'
!
      nchadd(1) = 0
      nchadd(2) = 0

! chosen tile with starting point now draws a trajectory (vertical for now in MPI)
      IF ( my_rank == ichooserank ) THEN !{

       ! shift starting height to w/Ez point with largest efield:
       e1 = elec(i1,k1,j1,3)   ! w-point below the scalar point
       e2 = elec(i1,k1+1,j1,3)
       IF ( Abs(e1) > Abs(e2) ) THEN
         trje(1,1,llz) = zlev(kslgt,islgt_local,jslgt_local) - 0.5*dzz(islgt_local,kslgt,jslgt_local)
         trje(2,1,llz) = zlev(kslgt,islgt_local,jslgt_local) - 0.5*dzz(islgt_local,kslgt,jslgt_local)

       ELSE
!         trje(1,1,llz) =  z1d2(kslgt+1,2) ! (kslgt-1)*dz
!         trje(2,1,llz) =  z1d2(kslgt+1,2) ! (kslgt-1)*dz
         trje(1,1,llz) = zlev(kslgt,islgt_local,jslgt_local) + 0.5*dzz(islgt_local,kslgt+1,jslgt_local)
         trje(2,1,llz) = zlev(kslgt,islgt_local,jslgt_local) + 0.5*dzz(islgt_local,kslgt+1,jslgt_local)
       ENDIF
         
!
!
!  Trajectory number counter
!
      do in = 1,nenum
      
      ifoundic = 0
      idownward = 0
!
!  Iteration counter on drawing channel
!  (Number of iterations to find decenet channel and volume)
!
      do it = 1,nttim-1
!
!
!  actual interpolation (mlint1 does tri-linear interpolation)
!
!
!     subroutine mlint1
!    > (ndebug,ac,n1b,n1e,n2b,n2e,n3b,n3e,n4b,n4e,
!    >  nx,ny,nz,aint,facx,facy,facz,ic,jc,kc,lc)
!
!
!  find lower left corner of grid cell that trajectory
!  end point is currently in for interpolation of scalars
!
      if ( ndebug > 2 ) write(0,*) 'get ic, jc, kc, '

! subtract off igslg0/jgslg0 to translate to local coordinates

      icdom = ifix((trje(in,it,llx)+eps)/dx)  + 1
      ic = icdom - igslg0
      jcdom = ifix((trje(in,it,lly)+eps)/dy) + 1
      jc = jcdom  - jgslg0
      iepdom = ifix((trje(in,it,llx)+eps)/dx + 0.5) + 1  
      iep = iepdom - igslg0
      jedom = ifix((trje(in,it,lly)+eps)/dy + 0.5) + 1
      je = jedom - jgslg0
!      kc = ifix((trje(in,it,llz)+eps)/dz) + 1
      z1d2(1:nz,1) = zlev(1:nz,ic,jc)
      z1d2(1:nz,4) = dzz(ic,1:nz,jc)
!      kc = find_indexlgt(trje(in,it,llz),zlev(1:nz,ic,jc),nz)
      kc = find_indexlgt(trje(in,it,llz),z1d2(1,1),nz)

      if ( ndebug > 2 ) THEN
        write(0,*) 'ic, jc, kc =',ic,jc,kc,iep,je
        write(0,*) 'trje(llz),z1d2(1) = ',trje(in,it,llz),z1d2(1,1)
      ENDIF

!      z1d2(1,1) = zlev(1,ic,jc)
      z1d2(1,2) = 0.0
      z1d2(1,3) = dzz(ic,1,jc)
      DO kz = 2,nz-1
!       z1d2(kz,1) = zlev(kz,ic,jc)
       z1d2(kz,2) = z1d2(kz,1) + 0.5*dzz(ic,kz,jc) ! 0.5*(z1d2(kz-1,1) + z1d2(kz,1)) ! set scalar-point as average of w points
       z1d2(kz-1,4) = dzz(ic,kz,jc) ! z1d2(kz,2) - z1d2(kz-1,2)
       z1d2(kz,3) = 0.5*(dzz(ic,kz-1,jc) + dzz(ic,kz,jc) )
!       write(0,*) 'kz,z1d2 = ',kz,z1d2(kz,1),z1d2(kz,2),z1d2(kz,3),z1d2(kz,4),z1d2(kz,2) - z1d2(kz-1,2)
      ENDDO
      z1d2(nz,4) = z1d2(nz-1,4)
      z1d2(nz,2) = z1d2(nz-1,2) + z1d2(nz,4)
      
      
      ke = find_indexlgt(trje(in,it,llz),z1d2(1,2),nz)
      facx = (trje(in,it,llx)-(icdom-1)*dx)/dx
      facy = (trje(in,it,lly)-(jcdom-1)*dy)/dy
!      facz = (trje(in,it,llz)-(kc-1)*dz)/dz
      facz = (trje(in,it,llz)- z1d2(kc,1))/z1d2(kc,3) ! *z1d2(kc,3)
      facxe = (trje(in,it,llx)-(iepdom-1+0.5)*dx)/dx
      facye = (trje(in,it,lly)-(jedom-1+0.5)*dy)/dy
      facze = (trje(in,it,llz)- z1d2(ke,2))/(z1d2(ke,4)) ! dzz(kc) ! *z1d2(kc,3)
!
      idxlgt(in,it,1) = icdom
      idxlgt(in,it,2) = jcdom
      idxlgt(in,it,3) = kc
!     
!   print if in test mode.
!
      if ( ndebug >= 2 ) then
      write(*,*) '==================================================='
      write(*,*) '==================================================='
      write(*,*) 'FINDING SCALARS'
      write(*,*) 'IN=:',in,' FACX,Y,Z:',  facx, facy, facz
      write(*,*) 'IN=:',in,' IC,JC,KC:',  IC,JC,KC
      write(*,*)  
      end if
!
      if ( ic .ge. 1 .and. ic .le. ixe .and.    & 
     &     jc .ge. 1 .and. jc .le. jye .and.    & 
     &     kc .ge. 1 .and. kc .le. nz-2*kd1 ) then
!
!  Interpolate
!  efield stuff
!
      idebug = 0
!      if ( ndebug > 0 ) idebug = 1

      IF ( ntasks == 0 ) THEN
      ie = 1 ! ex
      call mlint2    & 
     & (idebug, elec, -nbw, nx+1+nbe, 1, nz, -nbs, ny+1+nbn, 1, 4,    & 
     &  nx, nz, ny, telec, facxe, facz,facy, iep,kc,jc, ie)
      trje(in,it,(lex)) = telec              

      ie = 2 ! ey
      call mlint2    & 
     & (idebug, elec,  -nbw, nx+1+nbe, 1, nz, -nbs, ny+1+nbn, 1, 4,    & 
     &  nx, nz, ny, telec, facx, facz,facye, ic,kc,je, ie)
      trje(in,it,(ley)) = telec              

      ELSE
       trje(in,it,(lex)) = 0.0
       trje(in,it,(ley)) = 0.0
      ENDIF
      
      ie = 3 ! ez
      call mlint2    & 
     & (idebug, elec,  -nbw, nx+1+nbe, 1, nz, -nbs, ny+1+nbn, 1, 4,    & 
     &  nx, nz, ny, telec, facx, facze,facy, ic,ke,jc, ie)
      trje(in,it,(lez)) = telec              

      ie = 4
      call mlint2    & 
     & (idebug, elec,  -nbw, nx+1+nbe, 1, nz, -nbs, ny+1+nbn, 1, 4,    & 
     &  nx, nz, ny, telec, facx, facz,facy, ic,kc,jc, ie)
      trje(in,it,(lemag)) = telec              

      call mlint2    & 
     & (idebug, pot,  -nbw, nx+1+nbe, 1, nz, -nbs, ny+1+nbn,     & 
     &  1, 1, nx, nz, ny, telec, facx, facz, facy,ic,kc,jc, 1)
      trje(in,it,lpot) = telec              
      

      call mlint2    & 
!     & (idebug, temg,  nbw0, nx+1+nbe0, kms, kme, nbs0, ny+1+nbn0,     & 
     & (idebug, temg,  ims, ime, kms, kme, jms, jme,     & 
     &  1, 1, nx, nz, ny, tem, facx, facz, facy,ic,kc,jc, 1)
      trje(in,it,ltemg) = tem
      if ( ndebug > 2 ) THEN
        write(0,*) 'trje(ltemg),tgrd = ',trje(in,it,ltemg), tgrnd
      ENDIF
!      write(0,*) 'temint: ic,kc,jc,facz = ',ic,kc,jc,facz
!      write(0,*) 'temint: nbw',nbw0,nx+1+nbe0, 1, nz,nbs0,ny+1+nbn0
!      write(0,*) 'temint: ims,ime = ',ims,ime,kms,kme,jms,jme
!      write(0,*) 'temint: tem,temg = ',tem,temg(ic,kc,jc),temg(ic,kc-1,jc),temg(ic,kc+1,jc)

!      do ie = 1,neelec
!      call mlint1    & 
!     & (idebug, elec,  1, nxelec, 1, nyelec, 1, nzelec, 1, neelec,    & 
!     &  nx, ny, nz, telec, facx, facy, facz,ic,jc,kc, ie)
!      trje(in,it,(lex-1+ie)) = telec              
!      end do
!
!  Interpolate
!  total charge sctot(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
!
      call mlint2    & 
     & (idebug, sctot,  -nbw, nx+1+nbe, 1, nz, -nbs, ny+1+nbn,     & 
     &  1, 1, nx, nz, ny, tcharge, facx, facz, facy,ic,kc,jc, 1)
      trje(in,it,lcharge) = tcharge              
!
!  Compute
!  lemag from lex,ley,lez
!
      trje(in,it,lemag2) =     & 
     &  ((trje(in,it,lex)**2)+(trje(in,it,ley)**2)    & 
     &  +(trje(in,it,lez)**2))**(0.5)
!
!  Set "pseudo" time step:  They are very small increments...With
!  a breakdown of 100 kV/m they are < 0.1 of a grid interval
!
      if ( in .eq. 1 ) idtsign = 1  ! toward neg. charge (pos. leader)
      if ( in .eq. 2 ) idtsign = -1 ! toward pos. charge (neg. leader)
!
      dtelec(in) = idtsign*dt*1.0e-7
!
      if ( it .gt. 1 ) then  
      dtcomp = (.1)*min(dx,dzz(ic,kc,jc))    & 
     &       / max(trje(in,it,lemag2),elgtthn)
      dtelec(in) = idtsign*dtcomp
      if ( ndebug .ge. 2 ) write(6,*) 'dtelec',dtelec(in),min(dx,dzz(ic,kc,jc)),    & 
     &  trje(in,it,lemag2),dtelec(in)*trje(in,it,lemag2)
      end if
!
!  Integrate forward in time (along streamline)
!
      IF ( .false. .and. ntasks == 1 ) THEN ! only do horizontal if single processor -- for now.
!      IF ( ntasks == 1 ) THEN ! only do horizontal if single processor -- for now.
      trje(in,it+1,llx) =     & 
     &  trje(in,it,llx) + dtelec(in) * trje(in,it,lex)
      trje(in,it+1,lly) =     & 
     &  trje(in,it,lly) + dtelec(in) * trje(in,it,ley)
      ELSEIF ( it >= 1 ) THEN
      trje(in,it+1,llx) = trje(in,it,llx)
      trje(in,it+1,lly) = trje(in,it,lly)
      ENDIF

      IF ( it >= 1 ) THEN  
        trje(in,it+1,llz) = trje(in,it,llz) + dtelec(in) * trje(in,it,lez)
     
     IF ( it == 1 .and. dtelec(in) * trje(in,it,lez) < 0.0 ) idownward = 1
!
!  Print debug stuff
!
      if ( ndebug .ge. 2 ) then
      write(*,*) 'TRAJ'
      write(*,*) 'in, it, times'
      write(*,*) 'trje(in,it,llx),trje(in,it,lly),trje(in,it,llz)'
      write(*,*) 'trje(in,it,lex),trje(in,it,ley),trje(in,it,lez)'
      write(*,*) 'trje(in,it,lpot),trje(in,it,lemag),trje(in,it,lemag2)'
      write(*,*) 'trje(in,it,lcharge)'
      write(*,*) 'elec(ic+1,jc,kc,ipot),elec(ic+1,jc+1,kc,ipot)'
      write(*,*) 'elec(ic+1,jc,kc+1,ipot),elec(ic+1,jc+1,kc+1,ipot)'
      write(*,*) 'elec(ic,jc,kc,ipot),elec(ic,jc,kc,ipot)'
      write(*,*) in, it, times
      write(*,*) trje(in,it,llx),trje(in,it,lly),trje(in,it,llz)
      write(*,*) trje(in,it,lex),trje(in,it,ley),trje(in,it,lez)
      write(*,*) trje(in,it,lpot),trje(in,it,lemag),trje(in,it,lemag2)
      write(*,*) trje(in,it,lcharge)
      write(*,*) pot(ic+1,kc,jc),pot(ic+1,kc,jc+1)
      write(*,*) pot(ic+1,kc+1,jc),pot(ic+1,kc+1,jc+1)
      write(*,*) pot(ic,kc,jc),pot(ic,kc,jc)
      write(*,*) '==================================================='
      write(*,*) '==================================================='
      end if

      ENDIF

      nttimendcg(in) = it
!
! check for ground flash .... 
      if (  idownward == 1 .and. &
          (trje(in,it,llz) .lt. zgrnd .or. (zgrnd < 0. .and. trje(in,it,ltemg) > tgrnd) &
            .or. kc <= 2 ) &
         .and. &
       ( ( trje(in,1,lez) >  10.e3 .and. trje(in,1,lpot) < -50.e6 ) .or.   &
         ( trje(in,1,lez) < -10.e3 .and. trje(in,1,lpot) >  50.e6 ) ) ) then
      IF ( ifoundic == 0 ) nttimend(in) = it
      nttimendcg(in) = it
      ichend(in) = 1 
      go to 1999
      end if

      if ( trje(in,it,lemag2) .lt. estop(kc) .and. ifoundic == 0 ) then
!      trje(in,it,lstop) = 10.
        IF ( ( in .eq. 1 .and. trje(in,it,lcharge) .lt. -chgmin ) .or.    & 
     &     ( in .eq. 2 .and. trje(in,it,lcharge) .gt.  chgmin ) .or.    & 
     &     trje(in,it,lemag2) .lt. 0.1*estop(kc) ) THEN
      
        ifoundic = 1
        nttimend(in) = Max(1,it-1) ! Max(1,it)
        ichend(in) = 2  ! set end type as 'cloud' for now -- check for air later
          IF ( idownward == 0 ) THEN
            go to 1999
          ENDIF
        ENDIF
      ELSEIF ( trje(in,it,lemag2) .lt. estopcg(kc) .and. it >= 1 ) THEN
    ! get here if we hit the estopcg going downward, but didn't reach zgrnd
    ! leave nttimend and ichend at the values when estop was reached
          go to 1999
      end if
 
      
!
!  Else ic,jc, kc at boundary
!  set lstop > 1 and nttimend(in) = it-1
!
      else 
!
!      trje(in,it,lstop) = 10.
      nttimend(in) = max(it-1,1)
      go to 1999
!
!  end if on computing trajectory
!
      end if
!
!  end do on it
!
      end do  ! it
!
1999  continue
      it = nttimend(in)

!
!  Bound computed values....set it+1 values = it (for plotting purposes)
!
      trje(in,it+1,llz) = max(0.0,trje(in,it,llz))
      trje(in,it+1,llz) = min( z1d2(nz-1,1),trje(in,it,llz))
      trje(in,it+1,lly) = max(0.0,trje(in,it,lly))
      trje(in,it+1,lly) = min((ny-2*jstag)*dy,trje(in,it,lly))
      trje(in,it+1,llx) = max(0.0,trje(in,it,llx))
      trje(in,it+1,llx) = min((nx-2*istag)*dx,trje(in,it,llx))
      trje(in,it+1,lex) = trje(in,it,lex)
      trje(in,it+1,ley) = trje(in,it,ley)
      trje(in,it+1,lez) = trje(in,it,lez)
      trje(in,it+1,lemag) = trje(in,it,lemag)
      trje(in,it+1,lemag2) = trje(in,it,lemag2)
      trje(in,it+1,lpot) = trje(in,it,lpot)
      trje(in,it+1,lcharge) = trje(in,it,lcharge)
      idxlgt(in,it+1,1) = idxlgt(in,it,1)
      idxlgt(in,it+1,2) = idxlgt(in,it,2)
      idxlgt(in,it+1,3) = idxlgt(in,it,3)
!     
      IF ( ndebug > 1 ) THEN
      write(*,*) '==================================================='
      write(*,*) '==================================================='
      write(*,*) 'TRAJ END'
      write(*,*) in, it, times
      write(*,*) trje(in,it,llx),trje(in,it,lly),trje(in,it,llz)
      write(*,*) trje(in,it,lex),trje(in,it,ley),trje(in,it,lez)
      write(*,*) trje(in,it,lpot),trje(in,it,lemag),trje(in,it,lemag2)
      write(*,*) trje(in,it,lcharge)
      write(*,*) '==================================================='
      write(*,*) '==================================================='
      ENDIF
!
      if ( ndebug > 1 .and. my_rank == 0 ) write(*,*) 'nttimend',nttimend(1),nttimend(2)
!
!  dump lightning x, i coordinates
!
      if ( ilgt .gt. 0 ) then
!     write(rnl,911) int(iliter)
!911  format(i2.2)
!     write(rstime,912) int(times)
!912  format(i6.6)
!     flgtchn = 'hlgtchn'//rnl//rstime
!     open(unit=9,file=flgtchn,form='formatted')
!     do in = 1,nenum
!     do it = 1,nttimend(in)
!     write(9,*) trje(in,it,llx),trje(in,it,lly),trje(in,it,llz)
!     end do
!     end do
!     close(9)
      end if
!
      end do  ! in = 1, nenum

      ENDIF ! } my_rank == ichooserank
!
#ifdef MPI
        CALL MPI_Bcast(nttimend, nenum, MPI_INTEGER, ichooserank, my_comm, mpi_error_code)
        CALL MPI_Bcast(nttimendcg, nenum, MPI_INTEGER, ichooserank, my_comm, mpi_error_code)
        CALL MPI_Bcast(ichend, nenum, MPI_INTEGER, ichooserank, my_comm, mpi_error_code)
     !   call mpi_bcast( nttimend

       ! send trje info to all processors
        n = nenum*nttim*ninfo
        CALL MPI_Bcast(trje, n, MPI_REAL, ichooserank, my_comm, mpi_error_code)

! idxlgt(nenum,nttim,3)       
        n = nenum*nttim*3
        CALL MPI_Bcast(idxlgt, n, MPI_INTEGER, ichooserank, my_comm, mpi_error_code)
      
!      IF ( ichooserank /= 0 ) THEN
!       ! communicate trje from ichooserank to zero
!       IF ( my_rank == 0 .or. my_rank == ichooserank ) THEN
!        ! trje(nenum,nttim,ninfo)
!         n = nenum*nttim*ninfo
!        
!         westward_tag = 10
!         IF ( my_rank == ichooserank ) THEN
!           Call MPI_Send(trje, n, MPI_REAL, 0, westward_tag, local_communicator, mpi_error_code)
!         ENDIF

         IF ( my_rank == 0 .and. ndebug > 1 ) THEN
!           Call MPI_Recv(trje, n, MPI_REAL, ichooserank, westward_tag, local_communicator, my_mpi_status, mpi_error_code)
      DO in = 1,nenum
       it = nttimend(in)
!      write(iunit,*) '==================================================='
      write(iunit,*) '==================================================='
      write(iunit,*) 'TRAJ END'
      write(iunit,*) in, it, times
      write(iunit,*) trje(in,it,llx),trje(in,it,lly),trje(in,it,llz)
      write(iunit,*) trje(in,it,lex),trje(in,it,ley),trje(in,it,lez)
      write(iunit,*) trje(in,it,lpot),trje(in,it,lemag),trje(in,it,lemag2)
      write(iunit,*) trje(in,it,lcharge),trje(in,it,ltemg)
       tmp = trje(in,it,llz)
       it = nttimendcg(in)
       IF ( trje(in,it,llz) < tmp .or. ichend(in) == 1 ) THEN
      write(iunit,*) 'last downward point: ', it
      write(iunit,*) trje(in,it,llx),trje(in,it,lly),trje(in,it,llz)
      write(iunit,*) trje(in,it,lex),trje(in,it,ley),trje(in,it,lez)
      write(iunit,*) trje(in,it,lpot),trje(in,it,lemag),trje(in,it,lemag2)
      write(iunit,*) trje(in,it,lcharge),trje(in,it,ltemg)
      write(iunit,*) 'starting potential = ',trje(in,1,lpot)
      
      write(iunit,*) 'zgrnd,ez,pot = ',zgrnd,trje(in,1,lez),trje(in,1,lpot)
        ENDIF
      write(iunit,*) '==================================================='
!      write(iunit,*) '==================================================='
        ENDDO
      if ( ndebug > 1 .and. my_rank == 0 ) write(iunit,*) 'nttimend',nttimend(1),nttimend(2)
        ENDIF
        
!       ENDIF
      
!      ENDIF
#endif


!  Determine which end of channel is at pos and neg charge.
!  Also, permit no air-to-air end point discharges or other un-physical
!  discharges.  Initiate isign1 & isign2 = 999
!
      isign1 = 999
      isign2 = 999
      
      chgmin = (1.0-(fscth)*(niter-1))*scth
!
!  Set flags for isign1:  if end of channel is in very small charge
!  then isign1 = 999
!
      IF (ichend(1) .gt. 1) THEN
!      if ( trje(1,nttimend(1),lcharge) .gt. (1.0e-14) ) then
      if ( trje(1,nttimend(1),lcharge) .gt. chgmin*0.01 ) then
      isign1 =  1
      end if
!      if ( trje(1,nttimend(1),lcharge) .lt. (-1.0e-14) ) then
      if ( trje(1,nttimend(1),lcharge) .lt. -chgmin*0.01 ) then
      isign1 = -1
      end if
      if ( abs(trje(1,nttimend(1),lcharge)) .lt. chgmin*0.01 ) then
      isign1 = 999
      end if
      
        IF ( trje(1,2,lpot) <  trje(1,1,lpot) .and. isign1 /= -1 ) THEN
          isign1 = -1
          trje(1,nttimend(1),lcharge) = -chgmin*0.1
        ELSEIF ( trje(1,2,lpot) >  trje(1,1,lpot) .and. isign1 /= 1 ) THEN
          isign1 = 1
          trje(1,nttimend(1),lcharge) = chgmin*0.1
        ENDIF
      
      END IF  !  (ichend(1) .ne. 1) 
!
!  Set flags for isign2:  if end of channel is in very small charge
!  then isign2 = 999
!
      IF (ichend(2) .gt. 1) THEN
      if ( trje(2,nttimend(2),lcharge) .gt. chgmin*0.01 ) then
      isign2 =  1
      end if
      if ( trje(2,nttimend(2),lcharge) .lt. -chgmin*0.01 ) then
      isign2 = -1
      end if
      if ( abs(trje(2,nttimend(2),lcharge)) .lt. chgmin*0.01 ) then
      isign2 = 999
      end if

        IF ( trje(2,2,lpot) >  trje(2,1,lpot)  .and. isign2 /= 1  ) THEN
          isign2 = 1
          trje(2,nttimend(2),lcharge) = chgmin*0.1
        ELSEIF ( trje(2,2,lpot) <   trje(2,1,lpot)  .and. isign2 /= -1  ) THEN
          isign2 = -1
          trje(2,nttimend(2),lcharge) = -chgmin*0.1
        ENDIF


      END IF  !  (ichend(2) .ne. 1) 
!

      IF ( my_rank == ichooserank ) THEN !{

!  Now calculate lnox using by chanlen and dv
        
     IF ( f_clnox ) THEN
     DO in = 1,2 
      chanlen = 0.0
      itt=2
        DO it=itt,nttimend(in)
            ic = ifix((trje(in,it,llx)+eps)/dx) + 1 - igslg0
            jc = ifix((trje(in,it,lly)+eps)/dy) + 1 - jgslg0  ! set lower-left corner
            kc = idxlgt(in,it,3) ! find_indexlgt(trje(in,it,llz),z1d2(1,1),nz)
         dl = Sqrt(    &
     &      (trje(in,it,llx) - trje(in,it-1,llx))**2 +    &
     &      (trje(in,it,lly) - trje(in,it-1,lly))**2 +    &
     &      (trje(in,it,llz) - trje(in,it-1,llz))**2 )
         chanlen = chanlen + dl

!     IF (chanlen .ge. dzz(ic,kc,jc)) THEN
     IF (chanlen .ge. dl .and. dl > 0.0) THEN
!        amin=4000.
!          DO ix=ic-1,ic+1
!            IF (Abs(trje(in,it,llx)-(ix-1)*dx).lt.amin) THEN
!              amin=Abs(trje(in,it,llx)-(ix-1)*dx)
!                IF (amin.lt.1.) THEN
!                  ic=ix
!                END IF
!             END IF
!           END DO
           
        cmin=500.
!         DO kz=kc-1,kc+1
!            IF (Abs(trje(in,it,llz)-(kz-1)*dzz(kz)).lt.cmin) THEN
!             cmin=Abs(trje(in,it,llz)-(kz-1)*dzz(kz))
!!                IF (cmin.lt.1.) THEN
!                   kc=kz
!                END IF
!             END IF
!           END DO
           
!        bmin=4000.
!          DO jy=jc-1,jc+1
!            IF (Abs(trje(in,it,lly)-(jy-1)*dy).lt.bmin) THEN
!              bmin=Abs(trje(in,it,lly)-(jy-1)*dy)
!                IF (bmin.lt.1.) THEN
!                   jc=jy
!                END IF
!             END IF
!           END DO
      ! Here, chanlen is the channel length, pn = pressure, db=air density, dv=dx*dy*dzz (cell volume)
              dv=dx*dy*dzz(ic,kc,jc)
              IF ( f_clnox ) THEN
                noxtmp = dl*(alnox+blnox*pn(ic,kc,jc))/avogadro
                nox(ic,kc,jc) = nox(ic,kc,jc) + noxtmp
                lnox(ic,kc,jc)=lnox(ic,kc,jc)+air_molmass*noxtmp/(db(ic,kc,jc)*dv)    
              ENDIF
      END IF
          
        END DO
        
          itt = it  
          chanlen = 0.0
        
      END DO
            
      ENDIF ! lnox > 1
            
! Now look for points within a widening cone along the channel in 
!  each direction that ends above ground

     DO in = 1,2
       IF (in .eq. 1) isignt = isign1
       IF (in .eq. 2) isignt = isign2
        nchadd(in) = 0
      IF (ichend(in) .gt. 1 .and. isignt .ne. 999) THEN
        chanlen = 0.0
        chcnt = 1.0
        DO it=2,nttimend(in)
         chanlen = chanlen + Sqrt(    &
     &      (trje(in,it,llx) - trje(in,it-1,llx))**2 +    &
     &      (trje(in,it,lly) - trje(in,it-1,lly))**2 +    &
     &      (trje(in,it,llz) - trje(in,it-1,llz) )**2 )
!          IF (chanlen .gt. 500.0*chcnt) THEN
          IF (chanlen .gt. z1d2(kc,4)*chcnt) THEN
            chcnt = chcnt + 1.0  ! increase to next check length
            rchan = (rbase + chanlen*tanth)  ! radius of sphere to look into         
            ic = ifix((trje(in,it,llx)+eps)/dx) + 1 - igslg0
            jc = ifix((trje(in,it,lly)+eps)/dy) + 1 - jgslg0  ! set lower-left corner
!            kc = ifix((trje(in,it,llz)+eps)* z1d2(kz,3)) + 1
            kc = idxlgt(in,it,3) ! find_indexlgt(trje(in,it,llz),z1d2(1,1),nz)
            imx = Int(rchan/dx) 
            jmy = Int(rchan/dy)   ! how many grid lengths can radius reach
            kmz = Int(rchan/dzz(ic,kc,jc))
            delx = trje(in,it,llx) - Real(ic - 1)*dx ! relative coords of
            dely = trje(in,it,lly) - Real(jc - 1)*dy ! channel point from
            delz = trje(in,it,llz) - z1d2(kc,1) ! Real(kc - 1)*dz ! lowerleft corner
            rchan = (rchan)**2  ! use r^2 to check and save doing square roots
            
            DO kz = Max(1,kc-kmz),Min(nz,kc+kmz+1)
            DO jy = Max(jyb,jc-jmy),Min(jye,jc+jmy+1)
            DO ix = Max(ixb,ic-imx),Min(ixe,ic+imx+1)
               rtest = (dx*Real(ix-ic) - delx)**2 +    &
     &                 (dy*Real(jy-jc) - dely)**2 +    &
     &                 (dzz(ic,kc,jc) - delz)**2
               IF (rtest .le. rchan .and.     &
     &             Abs(t2(ix,jy,kz)) .lt. 0.5)  THEN
                    
               IF ( Int(Sign(1.1,sctot(ix,kz,jy))) .eq. isignt .and.    &
     &               Abs(sctot(ix,kz,jy)) .gt. scth ) THEN
                    t2(ix,jy,kz) = -2.1*Real(isignt) 
                    t4(ix,jy,kz) = -1.01*Real(isignt)
                    nchadd(in) = nchadd(in) + 1
                    dv=dx*dy*dzz(ix,kz,jy)
                    IF ( f_clnox ) THEN
                      noxtmp = Sqrt(rtest)*(alnox+blnox*pn(ix,kz,jy))/avogadro   !here rtest is chanlen
                      nox(ix,kz,jy) = nox(ix,kz,jy) + noxtmp
                      lnox(ix,kz,jy) = lnox(ix,kz,jy)+air_molmass*noxtmp/(db(ix,kz,jy)*dv)
                    ENDIF
               END IF
               END IF        
            END DO
            END DO
            END DO
          END IF ! (chanlen .gt. 500.0*chcnt)
        END DO  ! it
     IF ( my_rank == 0 ) write(6,'(2(a,i3))') 'added ',nchadd(in),    &
     &       ' points along channel end ',in
        
      END IF    ! (ichend(in) .gt. 1 .and. isignt .ne. 999) 
      END DO    ! in

      

      IF (ichend(1) .gt. 1) THEN
! decide between cloud/precip vs. CG
      in = 1
!      it = nttimend(in)
!      i1 = idxlgt(in,it,1)
!      j1 = idxlgt(in,it,2)
!      k1 = idxlgt(in,it,3)
!      facx = (trje(in,it,llx)-(i1-1)*dx)/dx
!      facy = (trje(in,it,lly)-(j1-1)*dy)/dy
!      facz = (trje(in,it,llz)-z1d2(k1,1))/dzz(i1,k1,j1)

! leaving out the check for cloudy/precip vs. clear air for now
!       IF ( k1 > 5 ) THEN 
!        ichend(in) = 2 ! IC
!       ENDIF
        
      ENDIF

! leaving out the check for cloudy/precip vs. clear air for now for end 2


      ENDIF ! } my_rank == ichooserank



!
!
!  so, now have isign = 999 for a flash at ground or in clear air
!  otherwise, isign shows sign of charge in the cloud/precip at the
!  ends of the channel.  Also have isign = 999 if space charge is
!  below the absolute minimum we can go to.
!
!  Channel will fail for ground-ground or air-air or ground-air
!
      ifail = 0
      
      IF ( my_rank == 0 ) write(iunit,'(a,2(i3))') 'ichend(1),ichend(2)',ichend(1),ichend(2)
      IF ( ichend(1) .eq. 0 .or. ichend(2) .eq. 0 ) THEN
        IF ( my_rank == 0 ) THEN
        write(iunit,'(a)') 'FAILED LIGHTNING: one end not determined'
        ENDIF
        ifail = 1 !GOTO 1101
      END IF 
        IF ( my_rank == ichooserank .or. my_rank == 0 ) THEN
      if ( ndebug > 1 .and. my_rank == 0 ) write(iunit,*) 'isign1,isign2',isign1,isign2
      ENDIF
!
      if ( isign1 .eq. isign2 .and. ifail == 0 ) then
!        IF ( my_rank == ichooserank .or. my_rank == 0 ) THEN
      IF ( my_rank == 0 ) THEN
      write(iunit,*) 'Same sign of charge for both ends of lightning'
      write(iunit,*) 'FAILED LIGHTNING END SIGN PROBLEM-try again'
      write(iunit,*) 'isign1,trje(1,nttimend(1),lcharge)'
      write(iunit,*) isign1,trje(1,nttimend(1),lcharge)
      write(iunit,*) 'isign2,trje(2,nttimend(2),lcharge)'
      write(iunit,*) isign2,trje(2,nttimend(2),lcharge)
      ENDIF
!      go to 2999
      lgtstp = 1
      ifail = 1 ! GOTO 1101
      end if
!      
! check for ends in cloud but not enough charge
!  
! Change this to check instead for added points along channel....
! if at least two points were added, then allow flash.
      IF (ichend(1).eq.2 .and. isign1.eq.999 .and.  (nchadd(1) .lt. icountth .and. ifail == 0)  ) THEN
       IF ( my_rank == 0 ) write(iunit,'(a)') 'FAILED LIGHTNING at end1: in cloud but < chgmin'
      ifail = 1 ! GOTO 1101
      end if
      IF (ichend(2).eq.2 .and. isign2.eq.999 .and. (nchadd(2) .lt. icountth .and. ifail == 0) ) THEN
       IF ( my_rank == 0 )  write(iunit,'(a)') 'FAILED LIGHTNING at end2: in cloud but < chgmin'
      ifail = 1 ! GOTO 1101
      end if

      
!      IF ( my_rank /= ichooserank ) ifail = 0
!      
!#ifdef MPI
!       ichoose_send_buffer(1) = ifail
!       ichoose_send_buffer(2) = lgtstp
!       CALL MPI_Bcast(ichoose_send_buffer, 2, MPI_INTEGER, ichooserank, my_comm, mpi_error_code)
!       ifail = ichoose_send_buffer(1)
!       lgtstp = ichoose_send_buffer(2)
!#endif
       
       IF ( ifail > 0 ) GO TO 1101

!
!  Set indices to be used in rest of code to indicate which channel 
!  ends in positive charge and which in negative charge.
!
      ichneu = 0
      ichpos = 0
      ichneg = 0
!      if ( idxlgt(1,nttimend(1),3) .gt. 1 ) then
      IF ( ichend(1) .eq. 2 ) THEN
      if ( trje(1,nttimend(1),lcharge) .eq. 0.0 ) then
      ichneu=1
      end if
      if ( trje(1,nttimend(1),lcharge) .gt. 0.0 ) then
      ichpos=1
      ichneg=2
      end if
      if ( trje(1,nttimend(1),lcharge) .lt. 0.0 ) then
      ichneg=1
      ichpos=2
      end if
      end if
!
!      if ( idxlgt(2,nttimend(2),3) .gt. 1 ) then
      IF ( ichend(2) .eq. 2 ) THEN
      if ( trje(2,nttimend(2),lcharge) .eq. 0.0 ) then
      ichneu=2
      end if
      if ( trje(2,nttimend(2),lcharge) .gt. 0.0 ) then
      ichpos=2
      ichneg=1
      end if
      if ( trje(2,nttimend(2),lcharge) .lt. 0.0 ) then
      ichpos=1
      ichneg=2
      end if
      end if
!
      if ( ndebug > 1 .and. my_rank == 0 ) then
      write(iunit,*) 'ichneu,ichpos,ichneg',ichneu,ichpos,ichneg
      end if
!
!  Setting end points of lighting and other assorted goodies...
!
!      iretry = 0
!      iterst = 1
!
! 88887 continue

!
!  Iteration to ensure that space charge and potential values at end of
!  each channel is adequate and within tolerance chosen by user.
!  After niter iterations an attempt to discharge is aborted and is 
!  retried later.
!
!      do 89999 iter = iterst,niter
!
!  scthend# are temporary values for the charge at the end of 
!  channel 1 and 2.
!
      scthend1 =  1.0e10
      scthend2 = -1.0e10
      IF (ichend(ichpos).eq.2) scthend1 = trje(1,nttimend(1),lcharge)
      IF (ichend(ichneg).eq.2) scthend2 = trje(2,nttimend(2),lcharge)
!
!  scthendp and scthendn are the values of charge at the ends of 
!  the channels. 
!  
      scthendp =  1.0e10
      scthendn = -1.0e10
      IF (ichend(ichpos).eq.2) scthendp=trje(ichpos,nttimend(ichpos),lcharge)
      IF (ichend(ichneg).eq.2) scthendn=trje(ichneg,nttimend(ichneg),lcharge)
!
!  scthp and scthn are working values (minimums) of the space charge 
!  at the end of the channel required to outline the area of the 
!  discharge.
!
!      IF (ichend(ichpos).eq.2) scthp = scth ! Min(scth,scthendp,Abs(scthendn))
!      IF (ichend(ichneg).eq.2) scthn = -scth !-Min(scth,scthendp,Abs(scthendn))
      IF (ichend(ichpos).eq.2) scthp = Min(scth,scthendp,Abs(scthendn))
      IF (ichend(ichneg).eq.2) scthn = -Min(scth,scthendp,Abs(scthendn))
!
!  potthrp and potthrn are the values of potential at the ends of 
!  the channels.
!
      potthrp = 1.0e20
      potthrn = -1.0e20
      IF (ichend(ichpos).eq.2) potthrp=max(trje(1,nttimend(1),lpot),trje(2,nttimend(2),lpot))
      IF (ichend(ichneg).eq.2) potthrn=min(trje(1,nttimend(1),lpot),trje(2,nttimend(2),lpot))
!
!  isendn, jsendn, ksendn are the indices of the lower left hand corner
!  of the grid cell where the end of the negative channel is found.
!
!      isendn = idxlgt(ichneg,nttimend(ichneg),1) - igslg0 + 1
!      jsendn = idxlgt(ichneg,nttimend(ichneg),2) - jgslg0 + 1
      isendn = idxlgt(ichneg,nttimend(ichneg),1)
      jsendn = idxlgt(ichneg,nttimend(ichneg),2)
      ksendn = idxlgt(ichneg,nttimend(ichneg),3)
!
!  isendp, jsendp, ksendp are the indices of the lower left hand corner
!  of the grid cell where the end of the positive channel is found.
!
!      isendp = idxlgt(ichpos,nttimend(ichpos),1) - igslg0 + 1
!      jsendp = idxlgt(ichpos,nttimend(ichpos),2) - jgslg0 + 1
      isendp = idxlgt(ichpos,nttimend(ichpos),1)
      jsendp = idxlgt(ichpos,nttimend(ichpos),2)
      ksendp = idxlgt(ichpos,nttimend(ichpos),3)
!
! MPI: need to communicate the starting points for pos/neg branching regions

      if ( ndebug > 1 .and. my_rank == 0 ) then
      print*,'scthend1',trje(1,nttimend(1),lcharge)
      print*,'scthend2',trje(2,nttimend(2),lcharge)
      print*,'sendn',isendn,jsendn,ksendn
      
      print*,'sendp',isendp,jsendp,ksendp
      end if
      if ( ndebug > 1 .and. my_rank == 0 ) then
      write(iunit,*) 'scthp',scthp
      write(iunit,*) 'scthn',scthn
      write(iunit,*) 'scthendp',scthendp
      write(iunit,*) 'scthendn',scthendn
      write(iunit,*) 'potthrp',potthrp
      write(iunit,*) 'potthrn',potthrn
      write(iunit,*) 'ichend p/n = ',ichend(ichpos),ichend(ichneg)
      end if

      IF ( my_rank == ichooserank ) THEN !{
      ENDIF ! } my_rank == ichooserank
      

#ifdef MPI
    IF ( ntasks > 1 ) THEN
      IF ( my_rank == ichooserank ) THEN

          rchoose_send_buffer(1) = isendn
          rchoose_send_buffer(2) = jsendn
          rchoose_send_buffer(3) = ksendn
          rchoose_send_buffer(4) = isendp
          rchoose_send_buffer(5) = jsendp
          rchoose_send_buffer(6) = ksendp
          rchoose_send_buffer(7) = scthp
          rchoose_send_buffer(8) = scthn
          rchoose_send_buffer(9) = potthrp
          rchoose_send_buffer(10) = potthrn
          rchoose_send_buffer(11) = ichpos
          rchoose_send_buffer(12) = ichneg
          rchoose_send_buffer(13) = ichend(1)
          rchoose_send_buffer(14) = ichend(2)

      ENDIF ! } my_rank == ichooserank


       k = 14
       
       CALL MPI_Bcast(rchoose_send_buffer, k, MPI_REAL, ichooserank, local_communicator, mpi_error_code)

       IF ( my_rank /= ichooserank ) THEN
       isendn = Nint( rchoose_send_buffer(1) )
       jsendn = Nint( rchoose_send_buffer(2) )
       ksendn = Nint( rchoose_send_buffer(3) )
       isendp = Nint( rchoose_send_buffer(4) )
       jsendp = Nint( rchoose_send_buffer(5) )
       ksendp = Nint( rchoose_send_buffer(6) )
       scthp   = rchoose_send_buffer(7)
       scthn   = rchoose_send_buffer(8)
       potthrp = rchoose_send_buffer(9)
       potthrn = rchoose_send_buffer(10)
       ichpos    = Nint( rchoose_send_buffer(11) )
       ichneg    = Nint( rchoose_send_buffer(12) )
       ichend(1) = Nint( rchoose_send_buffer(13) )
       ichend(2) = Nint( rchoose_send_buffer(14) )
       
       ENDIF

    ENDIF
#endif
!

       
!       IF ( lgtstp .gt. 0 ) GOTO 2999

!  Count number of ICs and CGs and CGPs and CGNs and CAs
!
      ilgt = 0
      iic  = 0
      ica  = 0
      icg  = 0
      icgp = 0
      icgn = 0
      icgpn = 0
      icgyn = 0

!      i1=idxlgt(1,nttimend(1),1)
!      j1=idxlgt(1,nttimend(1),2)
!      k1=idxlgt(1,nttimend(1),3)
!      i2=idxlgt(2,nttimend(2),1)
!      j2=idxlgt(2,nttimend(2),2)
!      k2=idxlgt(2,nttimend(2),3)
      
      IF ( ichend(1) == 1 .or. ichend(2) == 1 ) THEN
        icgyn = 1
        IF ( ichend(ichpos) .eq. 2 ) icgpn = 1
        IF ( ichend(ichneg) .eq. 2 ) icgpn = -1
      ENDIF

!
! Is discharge an ic
!
      if ( icgyn .eq. 0 ) then

      iic = 1
      ilgt = 1
      iredis = 1 ! ndisic
      numic = numic + 1
      numlgt = numlgt + 1
      IF ( my_rank == 0 ) THEN
      write(iunit,*) '      IC DISCHARGE'
      write(iunit,*) '         TOTAL NUMBER OF ICs',numic
      write(iunit,*) '         TOTAL NUMBER OF LGT',numlgt
      ENDIF
      IF ( isa .lt. 1 ) THEN
!      write(24,'(i5,a,f13.3,3x,i4)') nstep,' IC ', dtp*(nstep - 1),
!     : iliter
      ENDIF

      end if
!

      if ( icgyn .eq. 1 ) then
!
      icg = 1
      ilgt = 1
      iredis = 1 ! ndiscg
      numcg = numcg + 1
      numlgt = numlgt + 1
      IF ( my_rank == 0 ) THEN
      write(iunit,*) '      CG DISCHARGE'
      write(iunit,*) '         TOTAL NUMBER OF CGs',numcg
      write(iunit,*) '         TOTAL NUMBER OF LGT',numlgt
      ENDIF
!      icgyn = 1
!
!
      if ( icgpn .eq. 1 ) then
      iicp = 1
      ilgt = 1
      loccur = 3
      numcgp = numcgp + 1
      IF ( my_rank == 0 ) THEN
      write(iunit,*) '      DISCHARGE IS POSITIVE'
      write(iunit,*) '         TOTAL NUMBER OF CGPs',numcgp
      write(iunit,*) '         TOTAL NUMBER OF LGT',numlgt
      ENDIF
      IF ( isa .lt. 1 ) THEN
!        write(24,'(i5,a,f13.3,3x,i4)') nstep,' CG POS ', 
!     :  dtp*(nstep - 1), iliter
      ENDIF
!      end if
!
      elseif ( icgpn .eq. -1 ) then
      iicn = 1
      ilgt = 1
      loccur = 2
      numcgn = numcgn + 1
      IF ( my_rank == 0 ) THEN
      write(iunit,*) '      DISCHARGE IS NEGATIVE'
      write(iunit,*) '         TOTAL NUMBER OF CGNs',numcgn
      write(iunit,*) '         TOTAL NUMBER OF LGT',numlgt
      ENDIF
      IF ( isa .lt. 1 ) THEN
!       write(24,'(i5,a,f13.3,3x,i4)') nstep,' CG NEG ', dtp*(nstep - 1),
!     : iliter
      ENDIF
      end if
!
      end if  ! ( icgyn.eq.1 ) 


!
!
!  flag the regions using t2 to state where there is potential for
!  lightning discharge:  sc>scth and pot>potthr 
!  also find scdel
!
!  Set t2 arrays...the first guess array...
!
!    number of points in each potential surface
      npsurfp = 0
      npsurfn = 0
      nppntp = 0
      nppntn = 0
      
      do kz = 1,nz-kd1
      do jy = jyb,jye
      do ix = ixb,ixe
! already zeroed and may have added points already...
!       t2(ix,jy,kz) = 0.0
      ! use rectilinear distance for potential change along channel
      idist = Abs((ixbeg + ix) - islgt) + Abs( (jybeg + jy ) - jslgt)
      rdist = Sqrt(float( ((ixbeg + ix) - islgt)**2 + ( (jybeg + jy ) - jslgt)**2 ))
      if ( ichend(ichpos) .eq. 2 ) then
!        if ( pot(ix,kz,jy) .gt. potthrp  + eint*(dx*idist + Abs(zlev(kslgt,ix,jy) - zlev(kz,ix,jy)))  ) then
        if ( pot(ix,kz,jy) .gt. potthrp + eint*Sqrt( (dx*rdist)**2 + Abs(zlev(kslgt,ix,jy) - zlev(kz,ix,jy))**2) ) then
!        if ( pot(ix,kz,jy) .gt. potthrp ) then ! + eint*Sqrt( (dx*rdist)**2 + Abs(zlev(kslgt,ix,jy) - zlev(kz,ix,jy))**2) ) then
         npsurfp = npsurfp + 1
          if ( sctot(ix,kz,jy) .gt. scthp ) then
!            t2(ix,jy,kz) = Max( 1.01, t2(ix,jy,kz) )
            t2(ix,jy,kz) = Min( -1.01, t2(ix,jy,kz) )
            nppntp = nppntp + 1
          end if
        end if
      end if
      if ( ichend(ichneg) .eq. 2) then
!        if ( pot(ix,kz,jy) .lt. potthrn - eint*(dx*idist + Abs(zlev(kslgt,ix,jy) - zlev(kz,ix,jy))) ) then
        if ( pot(ix,kz,jy) .lt. potthrn - eint*Sqrt( (dx*rdist)**2 + Abs(zlev(kslgt,ix,jy) - zlev(kz,ix,jy))**2)  ) then
!        if ( pot(ix,kz,jy) .lt. potthrn ) then ! - eint*Sqrt( (dx*rdist)**2 + Abs(zlev(kslgt,ix,jy) - zlev(kz,ix,jy))**2)  ) then
         npsurfn = npsurfn + 1
          if ( sctot(ix,kz,jy) .lt. scthn ) then
!            t2(ix,jy,kz) = Min( -1.01, t2(ix,jy,kz) )
            t2(ix,jy,kz) = Max( 1.01, t2(ix,jy,kz) )
            nppntn = nppntn + 1
          end if
        end if
      end if

      end do
      end do
      end do
      
#ifdef MPI
! find global integrated rate max
         mpitotinint(1)  = npsurfp
         mpitotinint(2)  = npsurfn
         mpitotinint(3)  = nppntp
         mpitotinint(4)  = nppntn


        CALL MPI_AllReduce(mpitotinint, mpitotoutint, 4, MPI_INTEGER, MPI_SUM, local_communicator, mpi_error_code)

        npsurfp = mpitotoutint(1)
        npsurfn = mpitotoutint(2)
        nppntp = mpitotoutint(3)
        nppntn = mpitotoutint(4)
       
#endif

      IF ( my_rank == 0 ) THEN
      write(iunit,*) 'potthrp,potthrn = ',potthrp,potthrn
      write(iunit,*) 'pnts in pos/neg pot. surfs ',npsurfp,npsurfn
      write(iunit,*) 'pot pnts in pos/neg regions',nppntp,nppntn
      ENDIF
!
!  search for positive continuous discharge region
!
      itype = -1
      cons1 = 0.0
      cons2 = 0.0
      icountp = 0
      icgnflg = 0
      if ( ichend(ichpos) .eq. 2 ) then
!      subroutine searchl    &
!     &  (itype,icount,is,js,ks,id1,jd1,kd1    &
!     &  ,nx,ny,nz,nxl,nyl,nzl,nor    &
!     &  ,cons1,cons2,tthr,tflg)

      IF ( mytask == ichooserank ) THEN
! set starting point
        t4(isendp-igslg0,jsendp-jgslg0,ksendp) = itype*maxdist
      ENDIF
      
!      write(0,*) 'call searchl for pos end'
      call searchl   &
     &  (itype,itermax,maxdist,icountp,isendp,jsendp,ksendp,id1,jd1,kd1   &
     &  ,nx,ny,nz,nx,ny,nz,nxend,nyend,1   &
     &  ,cons1,cons2,t2,t4,icgnflg,klgtmin,ibal)
      end if
!
!  search for negative continuous discharge region = positive leader breakdown region
!
      itype = 1
      cons1 = 0
      cons2 = 0
      icountn = 0
      icgpflg = 0
      if ( ichend(ichneg) .eq. 2 ) then
      IF ( mytask == ichooserank ) THEN 
! set starting point
        t4(isendn-igslg0,jsendn-jgslg0,ksendn) = itype*maxdist
      ENDIF
!      write(0,*) 'call searchl for neg end'
      call searchl   &
     &  (itype,itermax,maxdist,icountn,isendn,jsendn,ksendn,id1,jd1,kd1   &
     &  ,nx,ny,nz,nx,ny,nz,nxend,nyend,1   &
     &  ,cons1,cons2,t2,t4,icgpflg,klgtmin,ibal)
      end if
!      write(0,*) 'done searchl'
! calculate lnox in positive and negative continuous discharge regions
      IF ( f_clnox ) THEN
      do kz = 1,nz-kd1
      do jy = jyb,jye
      do ix = ixb,ixe
      IF (abs(t4(ix,jy,kz)).lt.abs(itype)*maxdist.and.abs(t4(ix,jy,kz)).gt.1.01) THEN
       dv=dx*dy*dzz(ix,kz,jy)
!       nox(ix,kz,jy)=min(dx,dzz(ix,kz,jy))*(alnox+blnox*pn(ix,kz,jy))/(db(ix,kz,jy)*avogadro)
!        nox(ix,kz,jy)=dzz(ix,kz,jy)*(alnox+blnox*pn(ix,kz,jy))/(db(ix,kz,jy)*avogadro)
         noxtmp = dx*(alnox+blnox*pn(ix,kz,jy))/avogadro
         nox(ix,kz,jy) = nox(ix,kz,jy) + noxtmp
       lnox(ix,kz,jy) = lnox(ix,kz,jy)+air_molmass*noxtmp/(db(ix,kz,jy)*dv)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDIF
      
!      write(0,*) 'icgpflg,icgnflg= ',icgpflg,icgnflg
      
      IF ( icgpflg /= 0 .or. icgnflg /= 0 ) THEN
!        write(iunit,*) 'CG? icgpflg,icgnflg = ',icgpflg,icgnflg
        IF ( icgpflg == 0 .or. icgnflg == 0 ) THEN
           icgpflg = Min( 1, icgpflg )
           icgnflg = Min( 1, icgnflg )
        ELSE ! both ends at ground, so treat as IC
         icgpflg = 0
         icgnflg = 0
        ENDIF
      ENDIF


!
      if ( ndebug > 1 .and. my_rank == 0 ) then
      write(iunit,*) 'ICNP: icountn,icountp=',icountn,icountp
      end if
!
      sizetest = 1
      
      if ( icountn .lt. icountth .and. ichend(ichneg) .eq. 2 )  then
      IF ( my_rank == 0 ) THEN
      write(iunit,*) 'icountn = ',icountn
      write(iunit,*) 'RETRY:IRETRY:ICOUNTN:STEP ',iretry,icountn,nstep
!      iterst = iterst+1
!      iretry = iretry+1
!      if ( iterst .gt. niter ) then
      write(iunit,*) 'Cant find enough points..try a new starting point'
      ENDIF
      
       sizetest = 0
      
      
!      go to 2999
!      end if
!      go to 1101
      end if
!
      if ( icountp .lt. icountth .and. ichend(ichpos) .eq. 2 ) then
      IF ( my_rank == 0 ) THEN
      write(*,*) 'icountp = ',icountp
!      iterst = iterst+1
      write(*,*) 'RETRY:IRETRY:ICOUNTP:STEP ',iretry,icountp,nstep
!      iretry = iretry+1
!      if ( iterst .gt. niter ) then
      write(*,*) 'Cant find enough points..try a new starting point'
      ENDIF

       sizetest = 0

!      go to 2999
!      end if
!      go to 1101
      end if
      

      IF ( sizetest == 0 ) THEN
        numlgt = numlgt - 1
        IF ( icgyn .eq. 0 ) THEN
          numic = numic - 1
        ELSE
          numcg = numcg - 1
          IF ( icgpn .eq. 1 ) THEN
            numcgp = numcgp - 1
          ELSEIF ( icgpn .eq. -1 ) THEN
            numcgn = numcgn - 1
          ENDIF
        ENDIF
       
        go to 1101
        
      ENDIF


!
!  define t3 (scdisch) using t4 to define regions where actual discharge 
!  will occur...
!
      irepeat = 0
59999 continue
      irepeat = irepeat +1
      scnetpi = 0.0
      scnetni = 0.0
!
      xy_init(:,:) = 0
      flsh_map(:,:,:) = 0.0
      
      do kz = 1,nz-kd1
      do jy = jyb,jye
      do ix = ixb,ixe
      
      scdisch(ix,jy,kz) = 0.0
      scdisch0(ix,jy,kz) = 0.0
      
       IF ( t4(ix,jy,kz) /= 0.0 ) THEN
         xy_init(ix,jy) = 1. ! using xy_init to designate the lightning footprint
         IF ( t4(ix,jy,kz) > 0.5 ) THEN
           flsh_map(ix,jy,1) = 1
         ELSE
           flsh_map(ix,jy,2) = 1
         ENDIF
       ENDIF

       IF ( t4(ix,jy,kz) /= 0.0 .and. lightmask(ix,jy) >= 2 ) THEN
         t4(ix,jy,kz) = 0.
       ENDIF

      if ( ichend(ichneg) .eq. 2) then
       if ( t4(ix,jy,kz) .gt. 0.5 ) then
        if ( sctot(ix,kz,jy) .lt. scthn ) then
!         scdisch0(ix,jy,kz) = -fprcnt*(Min(0.0,sctot(ix,kz,jy)+scth))
         scdisch0(ix,jy,kz) = -fprcnt*(Min(0.0,sctot(ix,kz,jy)))
!         IF ( t2(ix,jy,kz) .gt. 1.5 )  scdisch(ix,jy,kz) = -fprcnt*sctot(ix,kz,jy)
         scnetpi = scnetpi + sctot(ix,kz,jy)
        end if
       end if
      end if
      
      if (ichend(ichpos) .eq. 2 ) then
       if ( t4(ix,jy,kz) .lt. -0.5 ) then
        if ( sctot(ix,kz,jy) .gt. scthp ) then
!         scdisch0(ix,jy,kz) = -fprcnt*(Max(0.0,sctot(ix,kz,jy)-scth))
         scdisch0(ix,jy,kz) = -fprcnt*(Max(0.0,sctot(ix,kz,jy)))
!         IF ( t2(ix,jy,kz) .lt. -1.5 ) scdisch(ix,jy,kz) = -fprcnt*sctot(ix,kz,jy)
         scnetni = scnetni + sctot(ix,kz,jy)
        end if
       end if
      end if
      
      end do
      end do
      end do
!
#ifdef MPI
       mpitotindp(1) = scnetpi
       mpitotindp(2) = scnetni

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
       scnetpi = mpitotoutdp(1)
       scnetni = mpitotoutdp(2)
#endif

        IF ( my_rank == 0 ) THEN
        write(iunit,*) 'scnetpi,scnetni = ',scnetpi,scnetni
        ENDIF

#ifdef MPI


      IF ( number_of_processes .gt. 1 .and. nr > 0 .and. ninittot > 0 ) THEN
        
        westward_tag = 201
        CALL sendrecv_westward(nx,ny,1,nr,nr,0,nr,1,                &
     &        w_proc,e_proc,westward_tag,xy_init)

        eastward_tag = 202
        CALL sendrecv_eastward(nx,ny,1,nr,nr,0,nr,1,                 &
     &        w_proc,e_proc,eastward_tag,xy_init)

        southward_tag = 203
        CALL sendrecv_southward(nx,ny,1,nr,nr,0,nr,1,                 &
     &        n_proc,s_proc,southward_tag,xy_init)

        northward_tag = 204
        CALL sendrecv_northward(nx,ny,1,nr,nr,0,nr,1,                 &
     &        n_proc,s_proc,northward_tag,xy_init)

        westward_tag = 201
        CALL sendrecv_westward(nx,ny,2,nrflsh,nrflsh,0,nrflsh,1,                &
     &        w_proc,e_proc,westward_tag,flsh_map)

        eastward_tag = 202
        CALL sendrecv_eastward(nx,ny,2,nrflsh,nrflsh,0,nrflsh,1,                 &
     &        w_proc,e_proc,eastward_tag,flsh_map)

        southward_tag = 203
        CALL sendrecv_southward(nx,ny,2,nrflsh,nrflsh,0,nrflsh,1,                 &
     &        n_proc,s_proc,southward_tag,flsh_map)

        northward_tag = 204
        CALL sendrecv_northward(nx,ny,2,nrflsh,nrflsh,0,nrflsh,1,                 &
     &        n_proc,s_proc,northward_tag,flsh_map)
    
      ENDIF
#endif

       ! for IC flash, compare footprints of positive and negative discharge regions
       ! to limit extent of one beyond the other
       IF ( ichend(ichneg) .eq. 2 .and. ichend(ichpos) .eq. 2 ) THEN
       DO jy = 1,ny
        DO ix = 1,nx
          IF ( flsh_map(ix,jy,1) > 0.5 .or. flsh_map(ix,jy,2) > 0.5 ) THEN
            
            i1 = 0
            i2 = 0
            DO j = Max(-nrflsh+1,jy-nrflsh),Min(ny+nrflsh, jy+nrflsh)
             DO i = Max(-nrflsh+1,ix-nrflsh),Min(nx+nrflsh, ix+nrflsh)
               i1 = i1 + Nint(flsh_map(i,j,1))
               i2 = i2 + Nint(flsh_map(i,j,2))
             ENDDO
            ENDDO
            
            IF ( flsh_map(ix,jy,1) > 0.5 .and. i2 == 0 ) THEN
              flsh_map(ix,jy,1) = 0.0
              t4(ix,jy,1:nz-1) = 0.0
              scdisch(ix,jy,1:nz-1) = 0.0
              scdisch0(ix,jy,1:nz-1) = 0.0
            ENDIF
            IF ( flsh_map(ix,jy,2) > 0.5 .and. i1 == 0 ) THEN
              flsh_map(ix,jy,2) = 0.0
              t4(ix,jy,1:nz-1) = 0.0
              scdisch(ix,jy,1:nz-1) = 0.0
              scdisch0(ix,jy,1:nz-1) = 0.0
            ENDIF
            
          ENDIF
        ENDDO
       ENDDO
       ENDIF


!      ij_reduce(:,:) = 0

!  Project lightning footprint onto lightmask
!      IF ( nr > 0 ) THEN ! create circular region around each flash point 
       count = 0 
       tmp = -1
       DO jy = -nr+1,ny+nr
        DO ix = -nr+1,nx+nr
          IF ( xy_init(ix,jy) > 0.5 ) THEN
            count = count + 1

             lightmask(ix,jy) = 2 ! set to 2 in the lightning footprint + 1 in the boundary

!            DO j = Max(jy-1,1),Min(jy+1, ny)
!              DO i = Max(1,ix-1),Min(nx, ix+1)
!               lightmask(i,j) = 2 ! set to 2 in the lightning footprint + 1 in the boundary
!              ENDDO
!            ENDDO
            
            DO j = Max(-nr+1,jy-nym),Min(ny+nr, jy+nym)
             DO i = Max(-nr+1,ix-nxm),Min(nx+nr, ix+nxm)
               lightmask(i,j) = Max( lightmask(i,j) , ijmask( i-ix, j-jy ) ) ! set to 1 around the perifery to disallow nearby initiations
!               IF ( tmp < 0 ) THEN
!                write(0,*) 'mask: i,j,i-ix,j-jy,ijmask: ',i,j,i-ix,j-jy,ijmask( i-ix, j-jy )
!               ENDIF
             ENDDO
            ENDDO
            tmp = 1
          ENDIF
        ENDDO
       ENDDO
       
       counttot = count
#if MPI
      mpitotinint(1) = count
      CALL MPI_AllReduce(mpitotinint, mpitotoutint, 1, MPI_INTEGER, MPI_SUM, local_communicator, mpi_error_code)
      
      counttot = mpitotoutint(1)
      
#endif
       IF ( my_rank == 0 ) write(iunit,*) 'xy_init count,tot = ',count,counttot
       
       
!      ELSE
!        ij_reduce(:,:) = 1
!      ENDIF

      if ( ndebug .ge. 1 .and. my_rank == 0 ) then
!      write(*,*) 'scnetni,scnetpi',scnetni,scnetpi 
      end if
!
!  accumululate scnet and nscnet
!
      if ( ndebug > 2 .and. my_rank == 0 ) write(*,*) 'set scnet,nscnet'
!      volzone = dx*dy*dzz(
      voln = 0. ! icountn ! *volzone
      volp = 0. ! icountp ! *volzone
      scnetp = 0.0
      nscnetp = 0
      scnetn = 0.0
      nscnetn = 0
      scnet = 0.0
      nscnet = 0
      do kz = 1,nz-kd1
      do jy = jyb,jye
      do ix = ixb,ixe
       dv = dx*dy*dzz(ix,kz,jy)
       volzone = dv ! /gt(1,1,kz,imapz)
      if ( ichend(ichpos) .eq. 2 ) then
      if ( t4(ix,jy,kz) .gt. 0.5 ) then
      if ( sctot(ix,kz,jy) .lt. scthn ) then
      scnetp = scnetp + scdisch0(ix,jy,kz)*dv
      nscnetp = nscnetp + 1
      scnet = scnet + scdisch0(ix,jy,kz)*dv
      nscnet = nscnet + 1
      volp = volp + dv
      end if
      end if
      end if
      if ( ichend(ichneg) .eq. 2 ) then
      if ( t4(ix,jy,kz) .lt. -0.5 ) then
      if ( sctot(ix,kz,jy) .gt. scthp ) then
      scnetn = scnetn + scdisch0(ix,jy,kz)*dv
      nscnetn = nscnetn + 1
      scnet = scnet + scdisch0(ix,jy,kz)*dv
      nscnet = nscnet + 1
      voln = voln + dv
      end if
      end if
      end if
      end do
      end do
      end do
!
!      if ( ndebug .ge. 1 ) print*,'comp scnetp',scnetp,nscnetp
!      if ( ndebug .ge. 1 ) print*,'comp scnetn',scnetn,nscnetn
!      if ( ndebug .ge. 1 ) print*,'comp scnet',scnet,nscnet
!      if ( ndebug .ge. 1 ) print*,'comp volp',volp
!      if ( ndebug .ge. 1 ) print*,'comp voln',voln
!      if ( ndebug .ge. 1 ) print*,'COULOMB SCNETP',nstep,scnetp ! *volzone
!      if ( ndebug .ge. 1 ) print*,'COULOMB SCNETN',nstep,scnetn ! *volzone
!      if ( ndebug .ge. 1 ) print*,'COULOMB SCNET',nstep,scnet ! *volzone

!      lgtstp = 1
!      return

!
!
!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  If scnetn or scnetp is too small to accomodate the charge in
!  voln or volp, respectively, then make a volume sufficient
!  to accomodate the charge so breakdown will not occur.
!  after that go to 59999 and repeat this  procedure above.
!  only one attempt will be required.
!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
      if ( irepeat .gt. 1 ) then
      write(iunit,*) 'failed lightning owing to failure to redistribute charge'
      lgtstp =1
      GOTO 2999 ! return
      end if

!
!  accumululate scnet and nscnet
!
      scnetpi = 0.0
      scnetni = 0.0
!       chgrtp = 0.0
!       chgrtn = 0.0
       scdrop = 0.0
       sckeepp = 0.0
       sckeepn = 0.0
       idrop = 0

#ifndef MPI
! C$DOACROSS LOCAL(kz,jy,ix)
! c$omp  PARALLEL DO DEFAULT(SHARED), 
! c$omp+ PRIVATE(kz,jy,ix)
#endif
      do kz = 1,nz-kd1
      do jy = jyb,jye
      do ix = ixb,ixe
       dv = dx*dy*dzz(ix,kz,jy)
!      t6(ix,jy,kz) = 0.0
!      t7(ix,jy,kz) = 0.0
      if ( abs(t4(ix,jy,kz)) .gt. 0.5 ) then
       
     
!
!  If no mass is around, then zero the lightning charge,
!    except when the ambient charge (sctot) is not zero, in which
!    case we apply t3/scdisch, but only up to the amount in sctot.
!
!      IF ( t6(ix,jy,kz) .gt. 1.0e-12 ) THEN
!        t7(ix,jy,kz) = 1/t6(ix,jy,kz)
!        chgrtn = chgrtn + Min( 0.0, scdisch(ix,jy,kz) )
!        chgrtp = chgrtp + Max( 0.0, scdisch(ix,jy,kz) )
!      ELSE
       IF ( sctot(ix,kz,jy) .lt. 0.0 .and. scdisch0(ix,jy,kz) .gt. 0.0) THEN
        chg = Min ( Abs(sctot(ix,kz,jy)), scdisch0(ix,jy,kz) )
        scdrop = scdrop + (scdisch0(ix,jy,kz) - chg)*dv ! /gt(ix,jy,kz,imapz)
        scdisch0(ix,jy,kz) = chg
        sckeepp = sckeepp + chg*dv ! /gt(ix,jy,kz,imapz)
!        t7(ix,jy,kz) = 1.0
       ELSEIF ( sctot(ix,kz,jy) .gt. 0.0 .and. scdisch0(ix,jy,kz) .lt. 0.0) THEN
         chg = Max ( -(sctot(ix,kz,jy)), scdisch0(ix,jy,kz) )
        scdrop = scdrop + (scdisch0(ix,jy,kz) - chg)*dv ! /gt(ix,jy,kz,imapz)
        scdisch0(ix,jy,kz) = chg
        sckeepn = sckeepn + chg*dv ! /gt(ix,jy,kz,imapz)
!        t7(ix,jy,kz) = 1.0
       ELSE
        scdrop = scdrop + scdisch0(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
!        scdisch0(ix,jy,kz) = 0.0
!        t4(ix,jy,kz) = 0.0
       ENDIF
        idrop = idrop + 1
!        scnetpi = scnetpi + Max(0.0, sctot(ix,kz,jy) )*dv ! /gt(ix,jy,kz,imapz)
!        scnetni = scnetni + Min(0.0, sctot(ix,kz,jy) )*dv ! /gt(ix,jy,kz,imapz)
!      ENDIF
      end if
      end do
      end do
      end do
      
       scdropnet = scdropnet + scdrop

!       ENDIF !   ( isa .lt. 1 ) 

      
!      scalp = 1.0
!      scaln = 1.0
!      IF ( iic .eq. 1 ) THEN
!        IF ( Abs(chgrtn) .gt. chgrtp ) THEN
!          scaln = chgrtp/Abs(chgrtn)
!        ELESEIF ( Abs(chgrtn) .lt. chgrtp ) THEN
!          scalp = Abs(chgrtn)/chgrtp
!        ENDIF
!      ENDIF



      scnetpi = 0.0d0
      scnetni = 0.0d0
      if ( ndebug > 2 ) write(iunit,*) 'set scnet,nscnet'
      scnetp = 0.0d0
      nscnetp = 0
      scnetn = 0.0d0
      nscnetn = 0
      scnet = 0.0d0
      nscnet = 0
      voln = 0.0d0
      volp = 0.0d0
      chgtot = 0.0d0
      chgpnts = 0

!      ixe = nx
!      IF ( myproci == nproci ) ixe = nx-1
      
!      jye = ny
!      IF ( myprocj == nprocj ) jye = ny-1

      chgavailiter(:) = 0.0
      voliter(:) = 0.0

      do kz = 1,nz-kd1
      do jy = 1,jye
      do ix = 1,ixe
       dv = dx*dy*dzz(ix,kz,jy)
       volzone = dv ! /gt(1,1,kz,imapz)
        if ( t4(ix,jy,kz) .gt. 0.5 ) then
          scnetpi = scnetpi + sctot(ix,kz,jy)*dv ! /gt(1,1,kz,imapz)
          scnetp = scnetp + scdisch0(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
          nscnetp = nscnetp + 1
          scnet = scnet + scdisch0(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
          nscnet = nscnet + 1
          volp = volp + volzone
          n = Nint( t4(ix,jy,kz) )
          chgavailiter(n) = chgavailiter(n) + scdisch0(ix,jy,kz)*dv
          voliter(n) = voliter(n) + volzone
        end if
        if ( t4(ix,jy,kz) .lt. -0.5 ) then
          scnetni = scnetni + sctot(ix,kz,jy)*dv ! /gt(1,1,kz,imapz)
          scnetn = scnetn + scdisch0(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
          nscnetn = nscnetn + 1
          scnet = scnet + scdisch0(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
          nscnet = nscnet + 1
          voln = voln + volzone
          n = Nint( t4(ix,jy,kz) )
          chgavailiter(n) = chgavailiter(n) + scdisch0(ix,jy,kz)*dv
          voliter(n) = voliter(n) + volzone
        end if
        
        IF ( ij_reduce(ix,jy) == 1 ) THEN ! count up charge in all points in the flash columns
          chgtot = chgtot + sctot(ix,kz,jy)*dv
          IF ( kz == 1 ) chgpnts = chgpnts + 1
        ENDIF
        
      end do
      end do
      end do

#ifdef MPI

       n = 2*maxdist+1
       CALL MPI_AllReduce(chgavailiter(-maxdist), chgavailitersum(-maxdist), n, MPI_REAL, MPI_SUM, my_comm, mpi_error_code)
       
       chgavailiter(-maxdist:maxdist) = chgavailitersum(-maxdist:maxdist)
       chgavailitersum(:) = 0.0

       CALL MPI_AllReduce(voliter(-maxdist), chgavailitersum(-maxdist), n, MPI_REAL, MPI_SUM, my_comm, mpi_error_code)
       voliter(-maxdist:maxdist) = chgavailitersum(-maxdist:maxdist)
       chgavailitersum(:) = 0.0


! find global integrated rate max
       mpitotindp(1)  = scnetpi
       mpitotindp(2)  = scnetp
       mpitotindp(3)  = nscnetp
       mpitotindp(4)  = volp

       mpitotindp(5)  = scnetni
       mpitotindp(6)  = scnetn
       mpitotindp(7)  = nscnetn
       mpitotindp(8)  = voln

       mpitotindp(9)  = scnet
       mpitotindp(10) = nscnet

       mpitotindp(11)  = idrop
       mpitotindp(12)  = scdrop
       mpitotindp(13)  = sckeepp
       mpitotindp(14)  = sckeepn
       mpitotindp(15)  = scdropnet
       mpitotindp(16) = chgtot
       mpitotindp(17) = chgpnts

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 17, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
       scnetpi = mpitotoutdp(1)
       scnetp  = mpitotoutdp(2)
       nscnetp = Nint( mpitotoutdp(3) )
       volp    = mpitotoutdp(4)

       scnetni = mpitotoutdp(5)
       scnetn  = mpitotoutdp(6)
       nscnetn = Nint( mpitotoutdp(7) )
       voln    = mpitotoutdp(8)

       scnet   = mpitotoutdp(9)
       nscnet  = Nint( mpitotoutdp(10) )

       idrop     = Nint( mpitotoutdp(11) )
       scdrop    = mpitotoutdp(12)
       sckeepp   = mpitotoutdp(13)
       sckeepn   = mpitotoutdp(14)
       scdropnet = mpitotoutdp(15)
       chgtot    = mpitotoutdp(16)
       chgpnts   = mpitotoutdp(17)

#endif


        chgavailitersum(maxdist) = chgavailiter(maxdist)
        chgavailitersum(-maxdist) = chgavailiter(-maxdist)
        volitersum( maxdist) = voliter( maxdist)
        volitersum(-maxdist) = voliter(-maxdist)
          IF ( my_rank == 0 ) THEN
!            write(iunit,*) 'n, chgavailitersum(n), chgavailitersum(-n)'
          ENDIF
        DO n = maxdist-1,1,-1
          chgavailitersum( n) = chgavailitersum( n+1) + chgavailiter( n)
          chgavailitersum(-n) = chgavailitersum(-n-1) + chgavailiter(-n)

          volitersum( n) = volitersum( n+1) + voliter( n)
          volitersum(-n) = volitersum(-n-1) + voliter(-n)

          IF ( my_rank == 0 ) THEN
!            write(iunit,*) n, chgavailitersum(n), chgavailitersum(-n),volitersum(n),voliter(n)
          ENDIF
        ENDDO

!
!      voln = nscnetn*volzone
!      volp = nscnetp*volzone

      IF ( my_rank == 0 .and. ndebug > 1 ) THEN
      write(iunit,'(a,i6,2x,1pe12.5)')    &
     &  'Points, net charge dropped outside cloud: ',    &
     &   idrop,scdrop
      write(iunit,'(a,2(1x,1pe12.5))')    &
     &  'Net kept pos/neg lightning charge = ',    &
     &   sckeepp,sckeepn
      write(iunit,'(a,1pe12.5)')    &
     &  'Cumulative net charge dropped outside cloud: ',    &
     &   scdropnet
!      write(iunit,'(a,2(1x,1pe12.5))')    &
!     &  'Ambient net pos/neg charge outside cloud: ',    &
!     &   scnetpi, scnetni
      ENDIF
      
      IF ( my_rank == 0 ) THEN
      if ( ndebug .ge. 1 ) then
      write(iunit,'(a,2(1x,1pe12.5))') 'Ambient charge scnetni,scnetpi ',    &
     &    scnetni,scnetpi
      end if
!
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,1pe12.5,1x,i7)')     &
     &      'comp scnetp',scnetp,nscnetp
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,1pe12.5,1x,i7)')     &
     &     'comp scnetn',scnetn,nscnetn
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,1pe12.5,1x,i7)')     &
     &     'comp scnet',scnet,nscnet
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,1pe12.5)')     &
     &     'comp volp',volp
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,1pe12.5)')     &
     &     'comp voln',voln
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,i6,3x,g20.9)')     &
     &     'COULOMB SCNETP',nstep,scnetp
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,i6,3x,g20.9)')     &
     &     'COULOMB SCNETN',nstep,scnetn
      if ( ndebug .ge. 1 ) write(iunit,'(a,1x,i6,3(3x,g20.9))')     &
     &     'COULOMB SCNET',nstep, scnet, chgtot, chgpnts
      ENDIF
!      IF ( isa .lt. 1 ) THEN

!      write(24,'(a,i6,3x,g20.9,2(1x,i6))') 'COULOMB SCNETP',nstep,
!     >                       scnetp*volzone, nscnetp,pchnl
!      write(24,'(a,i6,3x,g20.9,2(1x,i6))') 'COULOMB SCNETN',nstep,
!     >                       scnetn*volzone, nscnetn,nchnl
!      write(24,'(a,i6,3x,g20.9)') 'COULOMB SCNET',nstep,
!     >                            scnet*volzone
      
!      ENDIF
!
!
!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  If scnetn or scnetp is to0 small to accomodate the charge in
!  voln or volp, respectively, then make a volume sufficient
!  to accomodate the charge so breakdown will not occur.
!  after that go to 59999 and repeat this  procedure above.
!  only one attempt will be required.
!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
      irepeat = 0
      if ( irepeat .gt. 1 ) then
      write(iunit,*) 'failed lightning owing to failure to redistribute charge'
      lgtstp =1
      GOTO 9898
      end if
!
!
!  compute scdelp (t3) [modified from scdel or t3]
!
!
      
      fracn = 1.0d0
      fracp = 1.0d0
      multp = 0.0d0
      multn = 0.0d0

!
! if ibal = 0 or 3 then want to use average maybe...
!
      if ( iic .eq. 1 ) then  
      if ( nscnet .gt. 0 ) then
!        sctarg = Abs( scnetn - scnet/Float(nscnet)*Float(nscnetn) )
        sctarg = Min( Abs(scnetn),Abs(scnetp) )
!        
! for ibal = 0 or 3 set sctarg to minimum of average or 2.0*the smaller.
!
        IF ( ibal .eq. 0 .or. ibal .ge. 3 ) THEN

         IF ( ibal .eq. 0 ) THEN
           sctarg = Min( 0.5d0*(Abs(scnetn) + Abs(scnetp)), 1.5d0*sctarg )
         ELSE
           sctarg = Min( 0.5d0*(Abs(scnetn) + Abs(scnetp)), 1.1d0*sctarg )
         ENDIF
         
         
          fracn = sctarg/Abs(scnetn)
          multn = 0.0
          fracp = sctarg/Abs(scnetp)
          multp = 0.0

        ELSEIF ( ibal == 1 ) THEN

         tmp = Min(  Abs(scnetn), Abs(scnetp) )
         iterp = 1
         itern = -1
          IF ( Abs(scnetn) .gt. Abs(scnetp) ) THEN
            ! find index itern that just exceeds scnetp
            DO n = -maxdist,-1
     !         IF ( my_rank == 0 ) write(iunit,*) 'n:',n,chgavailitersum(n),tmp
              IF ( Abs(chgavailitersum(n)) > tmp ) THEN
                itern = n
                EXIT
              ENDIF
            ENDDO
            scnetn = chgavailitersum(itern)
            voln = volitersum(itern)
            IF ( my_rank == 0 ) THEN
              write(iunit,*) 'New scnetn = ',scnetn
              write(iunit,*) 'New voln = ',voln
            ENDIF
            WHERE ( (Nint(t4) > itern ) .and. t4 < -0.5 )
              scdisch0 = 0.0
            END WHERE
          ELSE
            ! find index iterp that just exceeds scnetn
            DO n = maxdist,1,-1
       !       IF ( my_rank == 0 ) write(iunit,*) 'p:',n,chgavailitersum(n),tmp
              IF ( Abs(chgavailitersum(n)) > tmp ) THEN
                iterp = n
                EXIT
              ENDIF
            ENDDO

            scnetp = chgavailitersum(iterp)
            volp = volitersum(iterp)
            IF ( my_rank == 0 ) THEN
              write(iunit,*) 'New scnetp = ',scnetp
              write(iunit,*) 'New volp = ',volp
            ENDIF
            WHERE ( (NInt(t4) < iterp ) .and. t4 > 0.5 )
              scdisch0 = 0.0
            END WHERE
          
          ENDIF
          
          IF ( my_rank == 0 ) THEN
            write(iunit,*) 'itern,iterp = ',itern,iterp
          ENDIF
          
        if ( Abs(scnetn) .gt. Abs(scnetp) ) then
          IF ( chgtot < 0 ) THEN
          ! allow slightly unbalanced flash to account for net charge, but not more than the net charge
          sctarg = Min( Abs(scnetn),Min( 1.25*Abs(scnetp), Abs(scnetp) - chgtot)  )
          ELSE
          sctarg = Min( Abs(scnetn),Abs(scnetp) )
          ENDIF
          fracn = sctarg/Abs(scnetn)
          multn = 0.0d0
          fracp = 1.0d0
          multp = 0.0d0
        elseif ( Abs(scnetn) .lt. Abs(scnetp) ) then
          IF ( chgtot >  0 ) THEN
          sctarg = Min( Min( 1.25*Abs(scnetn), Abs(scnetn) + chgtot ) ,Abs(scnetp) )
          ELSE
          sctarg = Min( Abs(scnetn),Abs(scnetp) )
          ENDIF
          fracn = 1.0d0
          multn = 0.0d0
          fracp = sctarg/Abs(scnetp)
          multp = 0.0d0
        end if
          
        ELSE ! ibal = 2
!        write(iunit,*)  'sctarg,scnetn,scnetp = ',sctarg,scnetn,scnetp

         IF ( .true. ) THEN
         tmp = Min(  Abs(scnetn), Abs(scnetp) )
         iterp = 1
         itern = -1
          IF ( Abs(scnetn) .gt. Abs(scnetp) ) THEN
            ! find index itern that just exceeds scnetp
            DO n = -maxdist,-1
       !       IF ( my_rank == 0 ) write(iunit,*) 'n:',n,chgavailitersum(n),tmp
              IF ( Abs(chgavailitersum(n)) > tmp ) THEN
                itern = n
                EXIT
              ENDIF
            ENDDO
            scnetn = chgavailitersum(itern)
            voln = volitersum(itern)
            IF ( my_rank == 0 ) THEN
              write(iunit,*) 'New scnetn = ',scnetn
              write(iunit,*) 'New voln = ',voln
            ENDIF
            WHERE ( (t4 > itern + 0.5) .and. t4 < -0.5 )
              scdisch0 = 0.0
            END WHERE
          ELSE
            ! find index iterp that just exceeds scnetn
            DO n = maxdist,1,-1
        !      IF ( my_rank == 0 ) write(iunit,*) 'p:',n,chgavailitersum(n),tmp
              IF ( Abs(chgavailitersum(n)) > tmp ) THEN
                iterp = n
                EXIT
              ENDIF
            ENDDO

            scnetp = chgavailitersum(iterp)
            volp = volitersum(iterp)
            IF ( my_rank == 0 ) THEN
              write(iunit,*) 'New scnetp = ',scnetp
              write(iunit,*) 'New volp = ',volp
            ENDIF
            WHERE ( (t4 < iterp - 0.5) .and. t4 > 0.5 )
              scdisch0 = 0.0
            END WHERE
          
          ENDIF
          
          IF ( my_rank == 0 ) THEN
            write(iunit,*) 'itern,iterp = ',itern,iterp
          ENDIF
          
         ENDIF


        if ( Abs(scnetn) .gt. Abs(scnetp) ) then
          fracn = sctarg/Abs(scnetn)
          multn = 0.0d0
          fracp = 1.0d0
          multp = 0.0d0
        elseif ( Abs(scnetn) .lt. Abs(scnetp) ) then
          fracn = 1.0d0
          multn = 0.0d0
          fracp = sctarg/Abs(scnetp)
          multp = 0.0d0
        end if
        ENDIF
      end if
       IF ( my_rank == 0 ) THEN
         write(iunit,'(a,2(1x,f12.5))') 'fracn,multn',fracn,multn
         write(iunit,'(a,2(1x,f12.5))') 'fracp,multp',fracp,multp
       ENDIF
      end if

      if ( iic .eq. 1 ) then
      scnetnp = 0.0d0
      nscnetnp = 0
      scnetnn = 0.0d0
      nscnetnn = 0
      nscnetn = 0
      scnetn = 0.0d0
        if ( nscnet .gt. 0 ) then
          scnet = scnet/nscnet
        else 
          scnet = 0. 
        end if
      if ( ndebug .ge. 1 .and. my_rank == 0 ) write(iunit,'(a,1x,1pe12.5,1x,i7)')    &
     &      'scnet=',scnet,nscnet
      if ( ndebug > 2 ) write(iunit,'(a)') 'set scdelp--t3'

 !     ixe = nx
 !     IF ( myproci == nproci ) ixe = nx-1
 !     
 !     jye = ny
 !     IF ( myprocj == nprocj ) jye = ny-1

      channtot = 0.0
      chanptot = 0.0
      chanzn(:) = 0.0
      chanzp(:) = 0.0
      
   !   write(6,*) 'fracn, multn, scnetnn = ',fracn,multn,scnetnn
   !   write(6,*) 'fracp, multp, scnetnp = ',fracp,multp,scnetnp
      
      do kz = 1,nz-1
      do jy = 1,jye
      do ix = 1,ixe
       dv = dxx*dyy*dzz(ix,kz,jy)
      if ( t4(ix,jy,kz) .gt. 0.5 .and. scdisch0(ix,jy,kz) /= 0.0 ) then
      scdisch(ix,jy,kz) = fracp*scdisch0(ix,jy,kz) - multp*scnet
      scnetnp = scnetnp + scdisch(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
      nscnetnp = nscnetnp + 1
      scnetn = scnetn + scdisch(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
      nscnetn = nscnetn + 1
      chanzp(kz) = chanzp(kz) + dx
      chanptot = chanptot + dx
      end if
      if ( t4(ix,jy,kz) .lt. -0.5 .and. scdisch0(ix,jy,kz) /= 0.0 ) then
      scdisch(ix,jy,kz) = fracn*scdisch0(ix,jy,kz) - multn*scnet
      scnetnn = scnetnn + scdisch(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
      nscnetnn = nscnetnn + 1
      scnetn = scnetn + scdisch(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
      nscnetn = nscnetn + 1
      chanzn(kz) = chanzn(kz) + dx
      channtot = channtot + dx
      end if
      end do
      end do
      end do

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = scnetnp
       mpitotindp(2)  = nscnetnp
       mpitotindp(3)  = scnetnn
       mpitotindp(4)  = nscnetnn

       mpitotindp(5)  = scnetn
       mpitotindp(6)  = nscnetn
       mpitotindp(7)  = channtot
       mpitotindp(8)  = chanptot

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 8, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
       scnetnp  = mpitotoutdp(1)
       nscnetnp = Nint( mpitotoutdp(2) )
       scnetnn  = mpitotoutdp(3) 
       nscnetnn = Nint( mpitotoutdp(4) )

       scnetn = mpitotoutdp(5)
       nscnetn = Nint( mpitotoutdp(6) )
       channtot = mpitotoutdp(7)
       chanptot = mpitotoutdp(8)
#endif


      IF ( my_rank == 0 ) THEN

      if ( nscnetn .gt. 0 ) then
      write(iunit,'(a,1(1x,g20.9),1x,i7)') 'ADJ COULOMBS SCNETP',    &
     &  scnetnp, nscnetnp
      write(iunit,'(a,1(1x,g20.9),1x,i7)') 'ADJ COULOMBS SCNETN',    &
     &  scnetnn, nscnetnn
      write(iunit,'(a,1(1x,g20.9),1x,i7)') 'ADJ COULOMBS SCNET',    &
     &  scnetn, nscnetn
      write(iunit,'(a,3(1x,g20.9))') 'CHANNEL LENGTH N,P,TOT',    &
     &  channtot*1.e-3, chanptot*1.e-3, (channtot + chanptot)*1.e-3
     
      ENDIF

!        IF ( isa .lt. 1 ) THEN

!          write(24,'(a,2(1x,g20.9))') 'scnetn, scnetp',
!     :    scnetnn*(volzone),scnetnp*(volzone)
      
!        ENDIF
      
      end if
      
      
!
! calculate dipole moment of the neutralizing charge
!  (change in dipole moment)
       qpx = 0.0d0
       qpy = 0.0d0
       qpz = 0.0d0
       qpxn = 0.0d0
       qpyn = 0.0d0
       qpzn = 0.0d0
       qpxp = 0.0d0
       qpyp = 0.0d0
       qpzp = 0.0d0
       qpxm = 0.0d0
       qpym = 0.0d0
       qpzm = 0.0d0
       qpxcn = 0.0d0
       qpycn = 0.0d0
       qpzcn = 0.0d0
       qpxcp = 0.0d0
       qpycp = 0.0d0
       qpzcp = 0.0d0
       qnetn = 0.0d0
       qnetp = 0.0d0
       npntn = 0
       npntp = 0
      do kz = 1,nz-1
      do jy = 1,jye !ny-1
      do ix = 1,ixe ! nx-1
       dv = dxx*dyy*dzz(ix,kz,jy)
       IF (t4(ix,jy,kz) .lt. -1.0e-18) THEN
        qpxn = qpxn + (its+ix-0.5)*dx*scdisch(ix,jy,kz)*dv  ! fix "ix" and "jy" for MPI
        qpyn = qpyn + (jts+jy-0.5)*dy*scdisch(ix,jy,kz)*dv
        qpzn = qpzn +     &
     &         z1d2(kz,1)*scdisch(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
        qnetn = qnetn + scdisch(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
        npntn = npntn + 1
       ELSEIF (t4(ix,jy,kz) .gt. 1.0e-18) THEN
        qpxp = qpxp + (its+ix-0.5)*dx*scdisch(ix,jy,kz)*dv
        qpyp = qpyp + (jts+jy-0.5)*dy*scdisch(ix,jy,kz)*dv
        qpzp = qpzp +     &
     &         z1d2(kz,1)*scdisch(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
        qnetp = qnetp + scdisch(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
        npntp = npntp + 1
       END IF
      end do
      end do
      end do

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = qpxn
       mpitotindp(2)  = qpyn
       mpitotindp(3)  = qpzn

       mpitotindp(4)  = qpxp
       mpitotindp(5)  = qpyp
       mpitotindp(6)  = qpzp

       mpitotindp(7)  = qnetn
       mpitotindp(8)  = npntn

       mpitotindp(9)  = qnetp
       mpitotindp(10) = npntp

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 10, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
       qpxn = mpitotoutdp(1)
       qpyn = mpitotoutdp(2)
       qpzn = Nint( mpitotoutdp(3) )

       qpxp = mpitotoutdp(4)
       qpyp = mpitotoutdp(5)
       qpzp = mpitotoutdp(6)

       qnetn = mpitotoutdp(7)
       npntn = Nint( mpitotoutdp(8) )

       qnetp = mpitotoutdp(9)
       npntp = Nint( mpitotoutdp(10) )
#endif

! dipole moment components (C km)
        qpx = (qpxn + qpxp)*0.001d0
        qpy = (qpyn + qpyp)*0.001d0
        qpz = (qpzn + qpzp)*0.001d0
! location of midpoint of dipole (x,y,z in km)        
       IF ( qnetp .ne. 0.0 ) THEN
       qpxcn = -qpxn*0.001/(qnetp)
       qpycn = -qpyn*0.001/(qnetp)
       qpzcn = -qpzn*0.001/(qnetp)
       qpxcp = qpxp*0.001/(qnetp)
       qpycp = qpyp*0.001/(qnetp)
       qpzcp = qpzp*0.001/(qnetp)
        qpxm = (-qpxn + qpxp)*0.001/(2.0*qnetp)
        qpym = (-qpyn + qpyp)*0.001/(2.0*qnetp)
        qpzm = (-qpzn + qpzp)*0.001/(2.0*qnetp)
        ENDIF


      IF ( my_rank == 0 ) THEN
      
      write(iunit,'(a,2(1x,i7))') 'npntn, npntp ',npntn,npntp
      
      write(iunit,'(a,2(1x,1pe12.5))')     &
     &    'qnetn, qnetp',qnetn,qnetp

      IF ( isa .lt. 1 ) THEN

!      write(24,'(a,2(1x,1pe12.5))') 
!     :    'qnetn, qnetp',qnetn*(volzone),qnetp*(volzone)
      
      ENDIF
      
        write(iunit,'(a,4(1x,1pe12.5))') 'dipole moment (px,py,pz) ',    &
     &       qpx,qpy,qpz,Sqrt( qpx**2 + qpy**2 + qpz**2 )

      IF ( isa .lt. 1 ) THEN
!        write(24,'(a,3(1x,1pe12.5))') 'dipole moment (px,py,pz) ',
!     :       qpx,qpy,qpz,Sqrt( qpx**2 + qpy**2 + qpz**2 )
      ENDIF

        write(iunit,'(a,3(1x,f12.5))') 'dipole midpoint at (x,y,z) ',    &
     &    qpxm,qpym,qpzm   
      IF ( isa .lt. 1 ) THEN
!        write(24,'(a,3(1x,f12.5))') 'dipole midpoint at (x,y,z) ',
!     :    qpxm,qpym,qpzm   
      ENDIF

        write(iunit,'(a,3(1x,f12.5))') 'pos. drho center (x,y,z) ',    &
     &    qpxcp,qpycp,qpzcp   
        write(iunit,'(a,3(1x,f12.5))') 'neg. drho center (x,y,z) ',    &
     &    qpxcn,qpycn,qpzcn   

      IF ( isa .lt. 1 ) THEN
!        write(24,'(a,3(1x,f12.5))') 'pos. drho center (x,y,z) ',
!     :    qpxcp,qpycp,qpzcp   
!        write(24,'(a,3(1x,f12.5))') 'neg. drho center (x,y,z) ',
!     :    qpxcn,qpycn,qpzcn   
      ENDIF
      
      ENDIF ! my_rank == 0
        
      ELSE  ! cg or ca flash: calculate charge center
      
      IF ( ibal .eq. 2 .and. cgfr .ge. 1.0 ) THEN
        cgfrac = 1.0
        fracp = 1.0
        fracn = 1.0
        multp = 0.0
        multn = 0.0
      ELSE
        cgfrac = cgfr
      ENDIF
      
      IF ( (ibal .eq. 5)  .or.     &
     &    ( ibal .eq. 2 .and. cgfr .lt. 1.0 ) ) THEN
       IF ( my_rank == 0 ) THEN
        write(iunit,*) 'scnetp,scnetn',scnetp,scnetn
       ENDIF
      IF ( icgpn .eq. 1 ) THEN
! DISCHARGE IS POSITIVE, set scnetp <= cgfrac*scnetn
      fracn = 1.0
      multp = 0.0
      multn = 0.0
       IF ( Abs(scnetp) .gt. cgfrac*(Abs(scnetn))) THEN
          fracp = cgfrac*Abs(scnetn/scnetp)
       ELSE
          fracp = 1.0
       ENDIF
        
      ELSEIF ( icgpn .eq. -1 ) THEN
! DISCHARGE IS NEGATIVE, set scnetn <= cgfrac*scnetp
      fracp = 1.0
      multp = 0.0
      multn = 0.0
       IF ( Abs(scnetn) .gt. cgfrac*(Abs(scnetp))) THEN
          fracn = cgfrac*Abs(scnetp/scnetn)
       ELSE
          fracn = 1.0
       ENDIF
      
      ENDIF ! icgpn
        
        IF ( my_rank == 0 ) THEN
          write(iunit,'(a,(1x,f12.5))') 'fracn',fracn
          write(iunit,'(a,(1x,f12.5))') 'fracp',fracp
        ENDIF
      
      ENDIF ! (ibal .eq. 5 .or. (ibal .eq. 2 .and. cgfr .lt. 1.0) )

      scnetnp = 0.0
      nscnetnp = 0
      scnetnn = 0.0
      nscnetnn = 0
      scnetn = 0.0
      nscnetn = 0
      if ( ndebug .ge. 1 .and. my_rank == 0 ) write(iunit,'(a,1x,1pe12.5,1x,i7)')    &
     &      'scnet=',scnet,nscnet
      if ( ndebug > 2 ) write(iunit,'(a)') 'set scdelp--t3'

!      ixe = nx
!      IF ( myproci == nproci ) ixe = nx-1
      
!      jye = ny
!      IF ( myprocj == nprocj ) jye = ny-1

      do kz = 1,nz-1
       do jy = 1,jye
        do ix = 1,ixe
        dv = dxx*dyy*dzz(ix,kz,jy)
        if ( t4(ix,jy,kz) .gt. 0.5 ) then
          scdisch(ix,jy,kz) = fracp*scdisch0(ix,jy,kz)
          scnetnp = scnetnp + scdisch(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
          nscnetnp = nscnetnp + 1
          scnetn = scnetn + scdisch(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
          nscnetn = nscnetn + 1
        end if
        if ( t4(ix,jy,kz) .lt. -0.5 ) then
          scdisch(ix,jy,kz) = fracn*scdisch0(ix,jy,kz)
          scnetnn = scnetnn + scdisch(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
          nscnetnn = nscnetnn + 1
          scnetn = scnetn + scdisch(ix,jy,kz)*dv ! /gt(1,1,kz,imapz)
          nscnetn = nscnetn + 1
        end if
        end do
       end do
      end do

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = scnetnp
       mpitotindp(2)  = nscnetnp
       mpitotindp(3)  = scnetnn
       mpitotindp(4)  = nscnetnn

       mpitotindp(5)  = scnetn
       mpitotindp(6)  = nscnetn

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 6, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
       scnetnp  = mpitotoutdp(1)
       nscnetnp = Nint( mpitotoutdp(2) )
       scnetnn  = mpitotoutdp(3) 
       nscnetnn = Nint( mpitotoutdp(4) )

       scnetn = mpitotoutdp(5)
       nscnetn = Nint( mpitotoutdp(6) )
#endif

      if ( nscnetn .gt. 0 .and. my_rank == 0 ) then
      write(iunit,'(a,1(1x,g20.9),1x,i7)') 'ADJ COULOMBS SCNETP',    &
     &  scnetnp, nscnetnp
      write(iunit,'(a,1(1x,g20.9),1x,i7)') 'ADJ COULOMBS SCNETN',    &
     &  scnetnn, nscnetnn
      write(iunit,'(a,1(1x,g20.9),1x,i7)') 'ADJ COULOMBS SCNET ',    &
     &  scnetn, nscnetn
      ENDIF ! nscnetn
        
      
       qpx = 0.0
       qpy = 0.0
       qpz = 0.0
       qpxn = 0.0
       qpyn = 0.0
       qpzn = 0.0
       qpxcn = 0.0
       qpycn = 0.0
       qpzcn = 0.0
       qpxcp = 0.0
       qpycp = 0.0
       qpzcp = 0.0
       qnetn = 0.0
!       qnetp = 0.0
       npntn = 0
!       npntp = 0

!      ixe = nx
!      IF ( myproci == nproci ) ixe = nx-1
      
!      jye = ny
!      IF ( myprocj == nprocj ) jye = ny-1

      do kz = 1,nz-1
       do jy = 1,jye
        do ix = 1,ixe

       IF ((t4(ix,jy,kz)) .ne. 0.0) THEN

         dv = dxx*dyy*dzz(ix,kz,jy)

       IF ( (icgpn .eq. -1 .or. ibal .eq. 5 .or. ibal .eq. 2)    &
     &      .and. scdisch(ix,jy,kz) .gt. 0.0 ) THEN
        qpxn = qpxn + dx*(its+ix-0.5)*scdisch(ix,jy,kz)*dv
        qpyn = qpyn + dy*(jts+jy-0.5)*scdisch(ix,jy,kz)*dv
!        qpzn = qpzn + (float(kz)-0.5)*scdisch(ix,jy,kz)

        qpzn = qpzn +     &
     &         z1d2(kz,1)*scdisch(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)

        qnetn = qnetn + scdisch(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
        npntn = npntn + 1
       ELSEIF ( (icgpn .eq. 1 .or. ibal .eq. 5 .or. ibal .eq. 2)    &
     &      .and. scdisch(ix,jy,kz) .lt. 0.0 ) THEN
        qpxn = qpxn + dx*(its+ix-0.5)*scdisch(ix,jy,kz)*dv
        qpyn = qpyn + dy*(jts+jy-0.5)*scdisch(ix,jy,kz)*dv
!        qpzn = qpzn + (float(kz)-0.5)*scdisch(ix,jy,kz)

        qpzn = qpzn +     &
     &         z1d2(kz,1)*scdisch(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)

        qnetn = qnetn + scdisch(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
        npntn = npntn + 1
       ENDIF

       ENDIF ! t4
      end do
      end do
      end do


#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = qpxn
       mpitotindp(2)  = qpyn
       mpitotindp(3)  = qpzn

!       mpitotindp(4)  = qpxp
!       mpitotindp(5)  = qpyp
!       mpitotindp(6)  = qpzp

       mpitotindp(7)  = qnetn
       mpitotindp(8)  = npntn

!       mpitotindp(9)  = qnetp
!       mpitotindp(10) = npntp

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 8, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
       qpxn = mpitotoutdp(1)
       qpyn = mpitotoutdp(2)
       qpzn = Nint( mpitotoutdp(3) )

!       qpxp = mpitotoutdp(4)
!       qpyp = mpitotoutdp(5)
!       qpzp = mpitotoutdp(6)

       qnetn = mpitotoutdp(7)
       npntn = Nint( mpitotoutdp(8) )

!       qnetp = mpitotoutdp(9)
!       npntp = Nint( mpitotoutdp(10) )
#endif



! dipole moment components (C km)
!        qpx = (qpxn + qpxp)*volzone*dx*0.001
!        qpy = (qpyn + qpyp)*volzone*dy*0.001
!        qpz = (qpzn + qpzp)*volzone*dz*0.001
! location of midpoint of dipole (x,y,z in km)
       IF ( qnetn .ne. 0.0 ) THEN
         qpxcn = qpxn*0.001/(qnetn)
         qpycn = qpyn*0.001/(qnetn)
         qpzcn = qpzn*0.001/(qnetn)
       ELSE
         qpxcn = 0.0
         qpycn = 0.0
         qpzcn = 0.0
       ENDIF
       
       IF ( my_rank == 0 ) THEN
       write(iunit,'(a,1pe12.5)')     &
     &   'Charge center of ',qnetn ! *volzone
       write(iunit,'(a,3(2x,f11.4))')     &
     &   'at x,y,z (km)', qpxcn,qpycn,qpzcn
       ENDIF

      IF ( isa .lt. 1 ) THEN
!       write(24,'(a,1pe12.5)') 
!     :   'Charge center of ',qnetn*volzone
!       write(24,'(a,3(2x,f11.4))') 
!     :   'at x,y,z (km)', qpxcn,qpycn,qpzcn
      ENDIF

!       write(24,'(7i6)') nstep, ilgt,iic,ica,icg,icgp,icgn
!       write(24,'(a,1pe12.5,a,3(2x,0pf11.4))') 
!     :   'Charge center of ',qnetn*volzone,' C is at x,y,z (km)',
!     :    qpxcn,qpycn,qpzcn
      
      end if  ! ( iic .eq. 1 ) 
!
!  sum total surface area exposed by hydrometeor species
!
      if ( ndebug .ge. 10 ) write(iunit,*) 'sum total particle area'
!
!  now redistribute the charge....
!
!      if ( ndebug .ge. 2 ) write(iunit,*)'redistribute charge',lscb,lsceq
      if ( nscnet .gt. 0 ) then

      IF ( isa .lt. 1 ) THEN
!
!  redistribute charge based on surface area for WRF version of code
!  SURFACE AREA ratios are calculated here
!
      if ( iredis .eq. 1 ) then
       
!       write(0,*) 'iredis=1 is not available anymore!  Stopping!'
!       STOP
! C$DOACROSS LOCAL(ia,kz,jy,ix,chgrtp,chgrtn )
! c$omp  PARALLEL DO DEFAULT(SHARED), 
! c$omp+ PRIVATE(ia,kz,jy,ix,chgrtp,chgrtn)
!      do ia = lscb,lsceq
       chgrtp = 0.0
       chgrtn = 0.0
       chgrtp2 = 0.0
       chgrtn2 = 0.0
       totneg0 = 0.0
       totpos0 = 0.0
       totneg1 = 0.0
       totpos1 = 0.0
       totneg2 = 0.0
       totpos2 = 0.0
       scnetnp = 0.0
       scnetnn = 0.0
       qnetp   = 0.0
       qnetn   = 0.0
       
      do kz = 1,nz-kd1
      do jy = jyb,jye
      do ix = ixb,ixe

       dv = dxx*dyy*dzz(ix,kz,jy)
       deninv = 1./db(ix,kz,jy)

       tmp = db(ix,kz,jy)*(sch(ix,kz,jy)+scr(ix,kz,jy)+sci(ix,kz,jy)+scs(ix,kz,jy)+scw(ix,kz,jy)+sciona(ix,kz,jy))
       if (present(qh)) tmp = tmp+db(ix,kz,jy)*schl(ix,kz,jy)

       totneg0 = totneg0 + Min(0.0, tmp*dv)
       totpos0 = totpos0 + Max(0.0, tmp*dv)


      scnetnp = scnetnp + Max(0.0, scdisch(ix,jy,kz) )*dv ! /gt(1,1,kz,imapz)
      scnetnn = scnetnn + Min(0.0, scdisch(ix,jy,kz) )*dv ! /gt(1,1,kz,imapz)

       IF ( scdisch(ix,jy,kz) /= 0.0 ) THEN
       

!     compute area fraction Ax/sum(Ax). where Ax=(Qx*rhoair/(Nx*rhox))**2/3 * Nx 
!   but qnr etc. is N/rhoair, and need to convert scdisch to C/kg
!   So have Ax=(Qx*rhoair/(Nx*rhoair*rhox))**2/3 * (Nx*rhoair) = (Qx/(Nx*rhox))**2/3 * (Nx*rhoair)
!  or Ax = rhoair*(Qx/(rhox))**2/3 * (Nx)**(1/3), and the factor of rhoair cancels out in the ratio Ax/Atot
    

      if(qnr(ix,kz,jy).lt.0.) qnr(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      if(qr(ix,kz,jy).lt.0.) qr(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      if(qng(ix,kz,jy).lt.0.) qng(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      if(qg(ix,kz,jy).lt.0.) qg(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      if(qni(ix,kz,jy).lt.0.) qni(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      if(qi(ix,kz,jy).lt.0.) qi(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      if(qns(ix,kz,jy).lt.0.) qns(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      if(qs(ix,kz,jy).lt.0.) qs(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      if(qnc(ix,kz,jy).lt.0.) qnc(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      if(qc(ix,kz,jy).lt.0.) qc(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      if (present(qh)) then
      if(qnh(ix,kz,jy).lt.0.) qnh(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      if(qh(ix,kz,jy).lt.0.) qh(ix,kz,jy)=0. ! if neg numbers get SQRT3 (-X) = NaN??...FORTRAN's favorite number.
      endif

      fracs=((qs(ix,kz,jy)/rhos)**(2./3.) ) * qns(ix,kz,jy)**(1./3.)
      fraci=((qi(ix,kz,jy)/rhoi)**(2./3.) ) * qni(ix,kz,jy)**(1./3.)
      frach=((qg(ix,kz,jy)/rhog)**(2./3.) ) * qng(ix,kz,jy)**(1./3.)
      fracw=((qc(ix,kz,jy)/rhow)**(2./3.) ) * qnc(ix,kz,jy)**(1./3.)
      fracr=((qr(ix,kz,jy)/rhor)**(2./3.) ) * qnr(ix,kz,jy)**(1./3.)

      fractot=fracs+fraci+frach &
     &             +fracr+fracw

      if (present(qh)) THEN
        frachl=((qh(ix,kz,jy)/rhoh)**2./3. ) * qnh(ix,kz,jy)**(1./3.)
        fractot=fractot+frachl
      endif
     
      IF ( fractot > 0.0 ) THEN
      
!'     reduce/redistribute scdisch on scx based on area fraction 
   
!      if (qs(ix,kz,jy).gt.qmin) then ! avoid division by zero if fractot = 0 because qx=0
      IF ( fracs > 0.0 ) THEN
      if (sctot(ix,kz,jy) >=  scthp) THEN
        scs(ix,kz,jy) = scs(ix,kz,jy) + scdisch(ix,jy,kz)*(fracs/fractot)*deninv
        chgrtn = chgrtn + scdisch(ix,jy,kz)*(fracs/fractot)*dv
      elseif (sctot(ix,kz,jy) <= -scthn) THEN 
        scs(ix,kz,jy) = scs(ix,kz,jy) + scdisch(ix,jy,kz)*(fracs/fractot)*deninv
        chgrtp = chgrtp + scdisch(ix,jy,kz)*(fracs/fractot)*dv
      endif
      endif

!      if (qi(ix,kz,jy).gt.qmin) then
      IF ( fraci > 0.0 ) THEN
      if (sctot(ix,kz,jy) >=  scthp) then
        sci(ix,kz,jy) = sci(ix,kz,jy) + scdisch(ix,jy,kz)*(fraci/fractot)*deninv
        chgrtn = chgrtn + scdisch(ix,jy,kz)*(fraci/fractot)*dv
      elseif (sctot(ix,kz,jy) <= -scthn) THEN 
         sci(ix,kz,jy) = sci(ix,kz,jy) + scdisch(ix,jy,kz)*(fraci/fractot)*deninv
        chgrtp = chgrtp + scdisch(ix,jy,kz)*(fraci/fractot)*dv
      endif
      endif

!      if (qr(ix,kz,jy).gt.qmin) then
      IF ( fracr > 0.0 ) THEN
      if (sctot(ix,kz,jy) >=  scthp) then
        scr(ix,kz,jy) = scr(ix,kz,jy) + scdisch(ix,jy,kz)*(fracr/fractot)*deninv
        chgrtn = chgrtn + scdisch(ix,jy,kz)*(fracr/fractot)*dv
      elseif (sctot(ix,kz,jy) <= -scthn) THEN 
        scr(ix,kz,jy) = scr(ix,kz,jy) + scdisch(ix,jy,kz)*(fracr/fractot)*deninv
        chgrtp = chgrtp + scdisch(ix,jy,kz)*(fracr/fractot)*dv
      endif
      endif

!      if (qg(ix,kz,jy).gt.qmin) then
      IF ( frach > 0.0 ) THEN
      if (sctot(ix,kz,jy) >=  scthp) then
        sch(ix,kz,jy) = sch(ix,kz,jy) + scdisch(ix,jy,kz)*(frach/fractot)*deninv
        chgrtn = chgrtn + scdisch(ix,jy,kz)*(frach/fractot)*dv
      elseif (sctot(ix,kz,jy) <= -scthn) THEN 
        sch(ix,kz,jy) = sch(ix,kz,jy) + scdisch(ix,jy,kz)*(frach/fractot)*deninv
        chgrtp = chgrtp + scdisch(ix,jy,kz)*(frach/fractot)*dv
      endif
      endif

!      if (qc(ix,kz,jy).gt.qmin) then
      IF ( fracw > 0.0 ) THEN
      if (sctot(ix,kz,jy) >=  scthp) then
        scw(ix,kz,jy) = scw(ix,kz,jy) + scdisch(ix,jy,kz)*(fracw/fractot)*deninv
        chgrtn = chgrtn + scdisch(ix,jy,kz)*(fracw/fractot)*dv
      elseif (sctot(ix,kz,jy) <= -scthn) THEN 
        scw(ix,kz,jy) = scw(ix,kz,jy) + scdisch(ix,jy,kz)*(fracw/fractot)*deninv
        chgrtp = chgrtp + scdisch(ix,jy,kz)*(fracw/fractot)*dv
      endif
      endif
   
      if (present(qh)) then
!      if (qh(ix,kz,jy).gt.qmin) then
      IF ( frachl > 0.0 ) THEN
      if (sctot(ix,kz,jy) >=  scthp) then
        schl(ix,kz,jy) = schl(ix,kz,jy) + scdisch(ix,jy,kz)*(frachl/fractot)*deninv
        chgrtn = chgrtn + scdisch(ix,jy,kz)*(frachl/fractot)*dv
      elseif (sctot(ix,kz,jy) <= -scthn) THEN 
        schl(ix,kz,jy) = schl(ix,kz,jy) + scdisch(ix,jy,kz)*(frachl/fractot)*deninv
        chgrtp = chgrtp + scdisch(ix,jy,kz)*(frachl/fractot)*dv
      endif
      endif
      endif

      qnetp = qnetp + Max(0.0, scdisch(ix,jy,kz) )*dv ! /gt(1,1,kz,imapz)
      qnetn = qnetn + Min(0.0, scdisch(ix,jy,kz) )*dv ! /gt(1,1,kz,imapz)
      
      ELSE
      
      if (sctot(ix,kz,jy) > 0.0) then
        sciona(ix,kz,jy) = sciona(ix,kz,jy) + scdisch(ix,jy,kz)*deninv
        chgrtn = chgrtn + scdisch(ix,jy,kz)*dv
      elseif (sctot(ix,kz,jy) <= 0.0) THEN 
        sciona(ix,kz,jy) = sciona(ix,kz,jy) + scdisch(ix,jy,kz)*deninv
        chgrtp = chgrtp + scdisch(ix,jy,kz)*dv
      endif

      qnetp = qnetp + Max(0.0, scdisch(ix,jy,kz) )*dv ! /gt(1,1,kz,imapz)
      qnetn = qnetn + Min(0.0, scdisch(ix,jy,kz) )*dv ! /gt(1,1,kz,imapz)
      
      ENDIF ! fractot > 0

      IF ( scdisch(ix,jy,kz) .lt. 0.0 ) THEN
!        an(ix,jy,kz,lscni) = an(ix,jy,kz,lscni)-scdisch(ix,jy,kz)/ec
        chgrtn2 = chgrtn2 + scdisch(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
      ELSEIF ( scdisch(ix,jy,kz) .gt. 0.0 ) THEN
!        an(ix,jy,kz,lscpi) = an(ix,jy,kz,lscpi)+scdisch(ix,jy,kz)/ec
        chgrtp2 = chgrtp2 + scdisch(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
      ENDIF
      
      ENDIF ! IF scdisch /= 0

!     UPDATE sctot     

       totneg1 = totneg1 + Min(0.0, sctot(ix,kz,jy)*dv)
       totpos1 = totpos1 + Max(0.0, sctot(ix,kz,jy)*dv)
       
       sctot(ix,kz,jy) = db(ix,kz,jy)*(sch(ix,kz,jy)+scr(ix,kz,jy)+sci(ix,kz,jy)+scs(ix,kz,jy)+scw(ix,kz,jy)+sciona(ix,kz,jy))
       if (present(qh)) sctot(ix,kz,jy) = sctot(ix,kz,jy)+db(ix,kz,jy)*schl(ix,kz,jy)

       totneg2 = totneg2 + Min(0.0, sctot(ix,kz,jy)*dv)
       totpos2 = totpos2 + Max(0.0, sctot(ix,kz,jy)*dv)
     
!      if (sctot(ix,kz,jy).lt.0.) write(0,*) 'IN DISCHARGE', sctot(ix,kz,jy),sch(ix,kz,jy),scr(ix,kz,jy),sci(ix,kz,jy),scs(ix,kz,jy),scw(ix,kz,jy),ix,jy,kz
!      if (sctot(ix,kz,jy).lt.0.) write(0,*) 'IN DISCHARGE', scw(ix,kz,jy),scdisch(ix,jy,kz),fracw,fractot,qc(ix,kz,jy),rhow,qnc(ix,kz,jy),fracs,fraci,frach,fracr,qr(ix,kz,jy),rhor,qnr(ix,kz,jy)
      end do ! ix
      end do ! jy
      end do ! kz
!       chgrt(1,ia-lscb+1) = chgrtp
!       chgrt(2,ia-lscb+1) = chgrtn
!      end do ! ia
      
      
#ifdef MPI
       mpitotindp(1)  = totpos0
       mpitotindp(2)  = totneg0
       mpitotindp(3)  = totpos1
       mpitotindp(4)  = totneg1
       mpitotindp(5)  = totpos2
       mpitotindp(6)  = totneg2
       mpitotindp(7)  = chgrtp
       mpitotindp(8)  = chgrtn
       mpitotindp(9)  = chgrtp2
       mpitotindp(10) = chgrtn2
       mpitotindp(11) = scnetnp
       mpitotindp(12) = scnetnn
       mpitotindp(13) = qnetp
       mpitotindp(14) = qnetn

      CALL MPI_Reduce(mpitotindp, mpitotoutdp, 14, MPI_DOUBLE_PRECISION, MPI_SUM, 0, local_communicator, mpi_error_code)


      IF ( mytask == 0 ) THEN ! only print from task 0

       totpos0 = mpitotoutdp(1) 
       totneg0 = mpitotoutdp(2) 
       totpos1 = mpitotoutdp(3) 
       totneg1 = mpitotoutdp(4) 
       totpos2 = mpitotoutdp(5) 
       totneg2 = mpitotoutdp(6) 
       chgrtp  = mpitotoutdp(7) 
       chgrtn  = mpitotoutdp(8) 
       chgrtp2 = mpitotoutdp(9) 
       chgrtn2 = mpitotoutdp(10)
       scnetnp = mpitotoutdp(11)
       scnetnn = mpitotoutdp(12)
       qnetp   = mpitotoutdp(13)
       qnetn   = mpitotoutdp(14)

      ENDIF
#endif

      IF ( mytask == 0 ) THEN ! only print from task 0

      write(iunit,'(a)') 'lightning charge breakdown by category: pos,neg'
!      DO ia = lscb,lsceq
!       chgrtp = chgrtp + chgrt(1,ia-lscb+1)
!       chgrtn = chgrtn + chgrt(2,ia-lscb+1)
!       write(iunit,'(i3,2(2x,1pe12.5))' ) ia,    &
!     &       chgrt(1,ia-lscb+1),chgrt(2,ia-lscb+1)
!      END DO ! ia
      
      write(iunit,'(a,2(2x,1pe12.5))') 'check totals: ',chgrtp,chgrtn
!      write(iunit,'(a,2(2x,1pe12.5))') 'check totals2: ',chgrtp2,chgrtn2
!      write(iunit,'(a,2(2x,1pe12.5))') 'check totals3: ',scnetnp,scnetnn
!      write(iunit,'(a,2(2x,1pe12.5))') 'check totals4: ',qnetp,qnetn
      write(iunit,'(a,3(2x,1pe12.5))') 'pos/neg0: ',totpos0,totneg0,totpos0+totneg0
!      write(iunit,'(a,3(2x,1pe12.5))') 'pos/neg1: ',totpos1,totneg1,totpos1+totneg1
      write(iunit,'(a,3(2x,1pe12.5))') 'pos/neg2: ',totpos2,totneg2,totpos2+totneg2
      
      ENDIF
      
      end if
!
!  redistribute based on fraction
!
      if ( iredis .eq. 2 ) then
!      do ia = lscb,lsceq
      do kz = 1,nz-kd1
      do jy = jyb,jye
      do ix = ixb,ixe
      if ( abs(t4(ix,jy,kz)) .gt. 0.5 ) then
!      an(ix,jy,kz,ia) = elgt1*an(ix,jy,kz,ia)
      end if
      end do
      end do
      end do
!      end do
      end if
!
!  redistribute to ions
!
      IF ( iredis .eq. 3 ) THEN

      chgrtp = 0.0
      chgrtn = 0.0

!      ixe = nx
!      IF ( myproci == nproci ) ixe = nx-1
      
!      jye = ny
!      IF ( myprocj == nprocj ) jye = ny-1

      do kz = 1,nz-kd1
      do jy = 1,jye
      do ix = 1,ixe
       dv = dxx*dyy*dzz(ix,kz,jy)
      IF ( scdisch(ix,jy,kz) .lt. 0.0 ) THEN
!        an(ix,jy,kz,lscni) = an(ix,jy,kz,lscni)-scdisch(ix,jy,kz)/ec
        chgrtn = chgrtn + scdisch(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
      ELSEIF ( scdisch(ix,jy,kz) .gt. 0.0 ) THEN
!        an(ix,jy,kz,lscpi) = an(ix,jy,kz,lscpi)+scdisch(ix,jy,kz)/ec
        chgrtp = chgrtp + scdisch(ix,jy,kz)*dv ! /gt(ix,jy,kz,imapz)
      ENDIF
      end do
      end do
      end do

#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = chgrtn
       mpitotindp(2)  = chgrtp

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
       chgrtn = mpitotoutdp(1)
       chgrtp  = mpitotoutdp(2)
#endif

      IF ( my_rank == 0 ) THEN
       write(iunit,'(a,2(2x,1pe12.5))') 'check totals: ',    &
     &  chgrtp,chgrtn
       ENDIF
      END IF
!
!  end nscnet
!
      if ( ndebug > 2 .and. my_rank == 0 ) write(iunit,*) 'redistribute charge done'
      end if
      
      ENDIF  ! ( isa .lt. 1 ) 
!      ENDDO ! ilit

!
!2999  continue
!
!  zero approprioate parts of the ab array
!
!      IF ( isa .ne. 1 ) THEN
!      if ( iliter .eq. nliter .or. lgtstp .eq. 1 ) then
! C$DOACROSS LOCAL(ia,kz,jy,ix)
! c$omp  PARALLEL DO DEFAULT(SHARED), 
! c$omp+ PRIVATE(ia,kz,jy,ix)
!      do ia = lscb,lsceq
!      do kz = 1,nz
!      do jy = 1,ny
!      do ix = 1,nx
!      ab(ix,jy,kz,ia) = 0.0
!      end do
!      end do
!      end do
!      end do
!      end if
!      
!      ENDIF !  ( isa .ne. 1 )
      

!    determine the clnox      
      total_nox = 0.0d0
      IF ( f_clnox ) THEN
      do kz = 1,nz-kd1
      do jy = jyb,jye
      do ix = ixb,ixe     
!         tmp = 0.1*Min(1.0, Max(0.0,Abs(scdisch0(ix,jy,kz)) - scth )/delq0) ! only make lnox for delta-q exceeding scth
         tmp = 0.1*Min(1.0, Max(0.0,Abs(scdisch(ix,jy,kz)) )/delq0) ! only make lnox for delta-q /= 0 and scale by delq0
         clnox(ix,kz,jy)=clnox(ix,kz,jy)+ tmp*lnox(ix,kz,jy)
         dv=dx*dy*dzz(ix,kz,jy)
         total_nox = total_nox + tmp*nox(ix,kz,jy)
!         IF ( kz == 2 .and. tmp > 0.0 ) THEN
!           write(0,*) 'scdisch1: ',ix,jy,kz, scdisch0(ix,jy,kz)
!           write(0,*) 'scdisch2: ',scdisch(ix,jy,kz),tmp,lnox(ix,kz,jy)
!         ENDIF
      end do
      end do
      end do

       write(6,*) 'total_nox1 = ',total_nox



#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = total_nox

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
       total_nox = mpitotoutdp(1)
#endif

       cnoxtot = total_nox
       write(6,*) 'total_nox = ',total_nox

      ENDIF ! f_clnox

!
!
!  put the lightning position in the t4 array (holding the lightning
!  location)...
!
!
49999  continue
!
!
     
!      IF ( iliter > 1 ) STOP

       t2(:,:,:)  = 0.0
       
       DO kz = 1,nz-kd1
       DO jy = 1,jye
       DO ix = 1,ixe
       
!        IF (Abs(t2(ix,jy,kz)) .lt. 0.5 .and.     &
!     &      Abs(t4(ix,jy,kz)) .gt. 0.5 ) THEN
!            t4(ix,jy,kz) = Sign( 0.6,t4(ix,jy,kz) )
!        ELSE
!            t4(ix,jy,kz) = t2(ix,jy,kz)
!        END IF

         IF ( Abs(t4(ix,jy,kz)) > 0.5 .and. scdisch(ix,jy,kz) /= 0.0 ) THEN
           
           IF ( Abs(scdisch(ix,jy,kz)) > scth ) t2(ix,jy,kz) =  Sign(1.0,t4(ix,jy,kz))  ! flip signs for consistency with DBM branched lightning
           IF ( scdisch(ix,jy,kz) < 0.0 ) flshn(ix,kz,jy) = flshn(ix,kz,jy) + Min(0.0, t2(ix,jy,kz))
           IF ( scdisch(ix,jy,kz) > 0.0 ) flshp(ix,kz,jy) = flshp(ix,kz,jy) + Max(0.0, t2(ix,jy,kz))
           lightdis(ix,jy) = lightdis(ix,jy) + 1
         ELSE
           t2(ix,jy,kz)  = 0.0
         ENDIF
         IF ( kz == 1 .and.  t2(ix,jy,kz) /= 0.0 ) THEN
!          write(iunit,*) 'got to kz=1? STOP. ix,jy,kz,t2 = ',ix,jy,kz,t2(ix,jy,kz),scdisch(ix,jy,kz),t4(ix,jy,kz) 
!          STOP
         ENDIF
         
!!         IF ( t2(ix,jy,kz) /= 0.0 ) THEN
!          dv=dx*dy*dzz(ix,kz,jy)
!          clnox(ix,kz,jy)=clnox(ix,kz,jy)+((1.5e21*nx/6.02e23)*273*8.314)/(dv*(1013.25-zlev(kz,ix,jy)/9.)*100)
!          clnox(ix,kz,jy)=clnox(ix,kz,jy)+((nx*1.5e21/6.02e23)*273*8.314)/(dv*(((44.3-zlev(kz,ix,jy)*1e-3)/4.95)**5.26))
!!         clnox(ix,kz,jy)=clnox(ix,kz,jy)+((dxx*1.5e21/6.02e23)*273*8.314)/(dv*pn(ix,kz,jy))
!          clnox(ix,kz,jy)=clnox(ix,jy,kz)+((nx*1.5e21/6.02e23)*273*8.314)/(dv*pn(ix,kz,jy))

!!         ENDIF
         
       END DO
       END DO
       END DO


       ncldtop = 0
       zcldtoptot = 0.0
       DO jy = 1,jye
       DO ix = 1,ixe
         ipos = 0
         ineg = 0
         IF ( ANY(t2(ix,jy,1:nz-1) /= 0.0 ) ) THEN
           IF ( ANY(t2(ix,jy,1:nz-1) > 0.0 ) ) ipos = 1
           IF ( ANY(t2(ix,jy,1:nz-1) < 0.0 ) ) ineg = 1
           lightdens(ix,jy) = lightdens(ix,jy) + 1  ! flash extent density
           IF ( icgyn == 0 ) THEN 
             flshfedic(ix,jy) = flshfedic(ix,jy) + 1
             flshfedicp(ix,jy) = flshfedicp(ix,jy) + ipos
             flshfedicn(ix,jy) = flshfedicn(ix,jy) + ineg
           ELSE
             flshfedcg(ix,jy) = flshfedcg(ix,jy) + 1
             flshfedcgp(ix,jy) = flshfedcgp(ix,jy) + ipos
             flshfedcgn(ix,jy) = flshfedcgn(ix,jy) + ineg
           ENDIF
           
           IF ( kcldtop(ix,jy) > 1 ) THEN
             kz = kcldtop(ix,jy)
             ncldtop = ncldtop + 1
             zcldtoptot = zcldtoptot + zlev(kz,ix,jy)
           ENDIF
           
         ENDIF
       END DO
       END DO
       
#ifdef MPI
! find global integrated rate max
       mpitotindp(1)  = zcldtoptot
       mpitotindp(2)  = ncldtop

      CALL MPI_AllReduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, mpi_error_code)

       
       zcldtoptot = mpitotoutdp(1)
       ncldtop    = Nint(mpitotoutdp(2))
#endif
       
       IF ( my_rank == 0 ) THEN
         write(6,*) 'Avg. cloud top, area ',zcldtoptot/ncldtop*1.e-3, ncldtop*dx*dy*1.e-3**2
       ENDIF
       
       lightsuccess = 1
       
       IF ( my_rank == ichooserank ) THEN
           light(islgt_local,jslgt_local) = light(islgt_local,jslgt_local) + 1
           flashcount = flashcount + 1
            flshi(islgt_local,kslgt_local,jslgt_local) = flshi(islgt_local,kslgt_local,jslgt_local) + 1
       ENDIF

      IF  ( ilgt .ge. 1 ) THEN
        loccur = Max(1, loccur)
      ELSE
        loccur = 0
      ENDIF
     
     ENDIF ! ilight == 2

 2999 CONTINUE

      IF ( allocated(xy_init) ) deallocate( xy_init   )
      IF ( allocated(ij_reduce) ) deallocate( ij_reduce )

      IF ( ilight == 2 .and. lightsuccess == 0 ) loccur = 0
      
      if ( ndebug > 1 .and. my_rank == 0 ) write(iunit,*) 'done, ilgt,loccur,lgtstp = ',ilgt,loccur,lgtstp
!
!
!  END OF LIGHTNING PARAMETERIZATION VERSION 2
!
!
!
 9898 continue

      deallocate( scdisch )
      deallocate( scdisch0 )

      end subroutine lightmsz


!
!
!     ##################################################################
!     ##################################################################
!     ######                                                      ######
!     ######                  FUNCTION ran0mpi                    ######
!     ######                                                      ######
!     ##################################################################
!     ##################################################################
!

      REAL FUNCTION ran0mpi(idum)

!  ##################################################################
!  PURPOSE: to generate a random number.  To continue a sequence,
!           do not change idum between calls to ran0mpi
! Minimal" random number generator of Park and Miller. Returns a
! uniform random deviate between 0.0 and 1.0. Set or reset idum to
! any integer value (except the unlikely value MASK) to initialize
! the sequence; idum must not be altered between calls for
! successive deviates in a sequence.
!
!  AUTHOR:  lifted straight out of Numerical Recipes in Fortran 77
!
!  ##################################################################
      
      implicit none
      
      INTEGER idum,IA,IM,IQ,IR,MASK
      REAL AM
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,MASK=123459876)
     
      INTEGER k
      
      idum=ieor(idum,MASK) 
!      XORing with MASK allows use of zero and other simple
!      bit patterns for idum. 
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k 
! Compute idum=mod(IA*idum,IM) without overflows by
!      Schrage's method. 
      if (idum.lt.0) idum=idum+IM
      ran0mpi=AM*idum                !   Convert idum to afloating result.
      idum=ieor(idum,MASK)        !    Unmask before return.
      return
      END FUNCTION ran0mpi

!  ##################################################################




      end module module_discharge_msz
      
      
 !
!
! #####################################################################
! #####################################################################

!
      subroutine searchl    &
     &  (itype,itermax,maxdist,icount,is,js,ks,id1,jd1,kd1    &
     &  ,nx,ny,nz,nxl,nyl,nzl,nxend,nyend,nor    &
     &  ,cons1,cons2,tthr,tflg,icgflg,klgtmin,ibal)
! 
#ifdef MPI
      USE module_dm, ONLY : &
         local_communicator, mytask, ntasks, ntasks_x, ntasks_y,mytask_x, mytask_y
      USE module_commasmpi ! , ONLY : w_proc,e_proc,s_proc,n_proc,mpi_setup_flag
#endif
      use module_boxmgsetup, only: igslg0, jgslg0, kgslg0

      implicit none

#ifdef MPI
      INCLUDE 'mpif.h'
#else
      integer :: mytask = 0
#endif
!
      integer :: klgtmin
      integer ndebug
      parameter (ndebug=0)
      integer n,icheck,myicheck,itype,icount,icount2,numflg
      integer :: ichecksum, icountsum
      integer, intent(IN) :: itermax,maxdist,nxend,nyend,ibal
      integer icgflg
      INTEGER :: mpi_error_code = 1      
      integer i0,j0,k0
      integer i1,j1,k1
      integer i2,j2,k2
      integer is,js,ks
      integer ix,jy,kz
      integer id1,jd1,kd1
      integer nx,ny,nz,nmax
      integer nxl,nyl,nzl,nor
      integer nxs0,nys0,nzs0
      integer nxs1,nys1,nzs1
      integer nxs2,nys2,nzs2
      integer nxe0,nye0,nze0
      integer nxe1,nye1,nze1
      integer nxe2,nye2,nze2
      integer :: iunit = 6
 
      integer mynxs1,mynys1
      integer mynxe1,mynye1
      integer mynxe2,mynye2
      integer mynxs2,mynys2
!
      real cons1,cons2,xsval
      real tthr(-nor+1:nxl+nor,-nor+1:nyl+nor,-nor+1:nzl+nor) ! t2
      real tflg(-nor+1:nxl+nor,-nor+1:nyl+nor,-nor+1:nzl+nor) ! t4
!
      integer westward_tag,eastward_tag,southward_tag,northward_tag

!  -----------------------------------

      nmax = itermax ! max(nx,nz,ny) + 2
#ifdef MPI
! communication for tflg
      IF ( ntasks > 1 ) THEN 
        IF ( ntasks_x > 1 ) THEN
        westward_tag = 201
        CALL sendrecv_westward(nx,ny,nz,nor,nor,nor,1,1,  &
     &        w_proc,e_proc,westward_tag,tflg)
        eastward_tag = 204
        CALL sendrecv_eastward(nx,ny,nz,nor,nor,nor,1,1,  &
     &        w_proc,e_proc,eastward_tag,tflg)
        ENDIF

        IF ( ntasks_y > 1 ) THEN
        southward_tag = 207
        CALL sendrecv_southward(nx,ny,nz,nor,nor,nor,1,1, &
     &        n_proc,s_proc,southward_tag,tflg)
        northward_tag = 210
        CALL sendrecv_northward(nx,ny,nz,nor,nor,nor,1,1,  &
     &        n_proc,s_proc,northward_tag,tflg)
        ENDIF


        IF ( ntasks_x > 1 ) THEN
        westward_tag = 201
        CALL sendrecv_westward(nx,ny,nz,nor,nor,nor,1,1,  &
     &        w_proc,e_proc,westward_tag,tthr)
        eastward_tag = 204
        CALL sendrecv_eastward(nx,ny,nz,nor,nor,nor,1,1,  &
     &        w_proc,e_proc,eastward_tag,tthr)
        ENDIF

        IF ( ntasks_y > 1 ) THEN
        southward_tag = 207
        CALL sendrecv_southward(nx,ny,nz,nor,nor,nor,1,1, &
     &        n_proc,s_proc,southward_tag,tthr)
        northward_tag = 210
        CALL sendrecv_northward(nx,ny,nz,nor,nor,nor,1,1,  &
     &        n_proc,s_proc,northward_tag,tthr)
        ENDIF
      
      ENDIF
#endif

      icgflg = 0
!
      nzs0=ks
      nze0=ks
      nys0=js
      nye0=js
      nxs0=is
      nxe0=is
!      do k0 = nzs0,nze0
!      do j0 = nys0,nye0
!      do i0 = nxs0,nxe0
!      if ( itype .eq. 1 ) then
!      tflg(i0,j0,k0) =nmax
!      end if
!      if ( itype .eq. -1 ) then
!      tflg(i0,j0,k0) =-nmax
!      end if
!      end do
!      end do
!      end do
!
      do n = 1,nmax
      if ( ndebug .ge. 2 ) write(iunit,*) 'n,is,js,ks',n,is,js,ks
      icheck = 0
      numflg = 0
      xsval = (nmax-n)
      nxs1 = max(1,is-n)
      nxe1 = min(nxend,is+n)
      nys1 = max(1,js-n)
      nye1 = min(nyend,js+n)
      nzs1 = max(klgtmin,ks-n)
      nze1 = min(nz,ks+n)
      
! get limits relative to this tile:
      mynxs1 = nxs1 - igslg0
      mynxe1 = nxe1 - igslg0
      mynys1 = nys1 - jgslg0
      mynye1 = nye1 - jgslg0
      
!      write(6,*) 'mynxs1 etc = ',mynxs1,mynxe1,mynys1,mynye1
      
      IF ( ( mynxs1 <= nx .and. mynxe1 >= 1 ) .and. &
     &     ( mynys1 <= ny .and. mynye1 >= 1 ) ) THEN
      
      nxs1 = Max(mynxs1,1)
      nxe1 = Min(mynxe1,nx)
      nys1 = Max(mynys1,1)
      nye1 = Min(mynye1,ny)

!      write(6,*) 'nxs1 etc = ',nxs1,nxe1,nys1,nye1
      
      do k1 = nzs1,nze1
      do j1 = nys1,nye1
      do i1 = nxs1,nxe1
      nxs2 = max(0,i1-1)
      nxe2 = min(nx+1,i1+1)
      nys2 = max(0,j1-1)
      nye2 = min(ny+1,j1+1)
      nzs2 = max(klgtmin,k1-1)
      nze2 = min(nz,k1+1)
      do k2 = nzs2,nze2
      do j2 = nys2,nye2
      do i2 = nxs2,nxe2
!
!  itype .eq. 1
!
      if ( itype .eq. 1 ) then
        if ( tthr(i1,j1,k1) .gt. 0.5 ) then
          numflg = numflg + 1
          if ( tflg(i2,j2,k2) .gt. 0.5 ) then 
            if(tflg(i1,j1,k1).lt.0.5) then
              icheck = icheck+1
              tflg(i1,j1,k1) = xsval
     
            end if
          end if
        end if
      end if
!
!  itype .eq. -1
!
      if ( itype .eq. -1 ) then
        if ( tthr(i1,j1,k1) .lt. -0.5 ) then
         numflg = numflg + 1
         if ( tflg(i2,j2,k2) .lt. -0.5 ) then 
           if(tflg(i1,j1,k1).gt.-0.5) then
             icheck = icheck+1
             tflg(i1,j1,k1) = -xsval
       
           end if
         end if
        end if
      end if
!
      end do
      end do
      end do
!
      end do
      end do
      end do
      
      ENDIF ! limits within tile
!

#ifdef MPI
! communication for tflg
      IF ( ntasks > 1 ) THEN 
        IF ( ntasks_x > 1 ) THEN
!        write(0,*) 'send_west, ntask_x = ',ntasks_x
!        CALL sendrecv_westward(nx,ny,nz,1,1,0,1,1,  &
!        CALL sendrecv_westward(nxslm,nyslm,nzslm,ghosti1,0,0,ghosti1,1,  &
        westward_tag = 201
        CALL sendrecv_westward(nx,ny,nz,nor,nor,nor,1,1,  &
     &        w_proc,e_proc,westward_tag,tflg)
        eastward_tag = 204
        CALL sendrecv_eastward(nx,ny,nz,nor,nor,nor,1,1,  &
     &        w_proc,e_proc,eastward_tag,tflg)
        ENDIF

        IF ( ntasks_y > 1 ) THEN
        southward_tag = 207
        CALL sendrecv_southward(nx,ny,nz,nor,nor,nor,1,1, &
     &        n_proc,s_proc,southward_tag,tflg)
        northward_tag = 210
        CALL sendrecv_northward(nx,ny,nz,nor,nor,nor,1,1,  &
     &        n_proc,s_proc,northward_tag,tflg)
        ENDIF
      
      ENDIF
#endif
#ifdef MPI
!
!
! MPI: reduce icheck as SUM
!
      myicheck = icheck
      IF ( ntasks > 1 ) THEN
        CALL MPI_AllReduce(icheck, ichecksum, 1, MPI_INTEGER, MPI_SUM, local_communicator, mpi_error_code)
        icheck = ichecksum
      ENDIF
#endif

      if ( ndebug .ge. 2 ) then
      write(iunit,*) 'search iteration at n,icheck',n,icheck
      write(iunit,*) 'icount,',icount
      write(iunit,*)  nzs1,nze1
      write(iunit,*)  nys1,nye1
      write(iunit,*)  nxs1,nxe1

      IF ( nzs1 == 1 ) THEN
       write(*,*) 'got to kz=1? STOP. is,js,ks = ',is,js,ks
!       STOP
      ENDIF
      end if


      if ( icheck .eq. 0 ) EXIT ! leave loop when no points are found
!
      end do
777   continue
!
      if ( ndebug .ge. 3 ) then
      write(iunit,*) 'done'
      do jy = 1,ny
      write(iunit,*) 'jy',jy
      do kz = nz,1,-1
      write(iunit,101) (tflg(ix,jy,kz),ix=is-10,is+10)
101   format(1x,22f3.0)
      end do
      end do
      end if
!
      icount = 0
!
      do kz = 1,nz
      do jy = 1,ny
      do ix = 1,nx
      if ( itype .eq. 1 ) then
      if ( tflg(ix,jy,kz) .gt. 0.5 ) then
!      tflg(ix,jy,kz) = 1.01
      icount = icount + 1
      end if
      end if
      if ( itype .eq. -1 ) then
      if ( tflg(ix,jy,kz) .lt. -0.5 ) then
!      tflg(ix,jy,kz) = -1.01
      icount = icount + 1
      end if
      end if
      IF ( kz == 1 .and. tflg(ix,jy,kz) /= 0.0 .and. Int(Sign(1., tflg(ix,jy,kz))) == itype ) icgflg = icgflg + 1
      end do
      end do
      end do

!        IF ( Any(tflg(:,:,1)  /= 0 ) ) icgflg = 1
#ifdef MPI
!
! MPI: reduce icount as SUM
!
      icount2 = icount
      IF ( ntasks > 1 ) THEN
        CALL MPI_AllReduce(icount, icountsum, 1, MPI_INTEGER, MPI_SUM, local_communicator, mpi_error_code)
        icount = icountsum
      ENDIF
#endif

      if ( ndebug .ge. 1 ) then
      write(iunit,*) 'icount',icount,icount2,icheck
      end if
!
      return
      end
!
!
!###########################################################################
!
!     ##################################################################
!     ######                                                      ######
!     ######             INTEGER FUNCTION FIND_INDEX              ######
!     ######                                                      ######
!     ##################################################################
!
!     PURPOSE:
!
!     This function returns the array index (here, the value returned by
!     find_index is designated as i) such that x is between xa(i) and xa(i+1).
!     If x is less than xa(1), then i=0 is returned.  If x is greater than
!     xa(n), then i=n is returned.  It is assumed that the values of
!     xa increase monotonically with increasing i.
!
!############################################################################
!
!     Author:  David Dowell (based on "locate" algorithm in Numerical Recipes)
!
!     Creation Date:  17 November 2004
!
!############################################################################

      integer function find_indexlgt(x, xa, n)

      implicit none

      integer n                          ! array size
      real xa(n)                         ! array of locations
      real x                             ! location of interest
      integer il, im, iu                 ! lower and upper limits, and midpoint


      il = 0
      iu = n+1
      
      IF ( x .gt. xa(n) ) THEN
        il = n+2
      ELSEIF ( x .lt. xa(1) ) THEN
        il = -1
      ELSE

 10   if ((iu-il).gt.1) then
        im=(il+iu)/2
        if (x.ge.xa(im)) then
          il=im
        else
          iu=im
        endif
        go to 10
      endif
      
      ENDIF

!      if ( (il.eq.0) .and. (n.gt.1) .and.
!           ( (xa(1)-x) .lt. 0.001*(xa(2)-xa(1)) ) ) then
!        il = 1
!      endif
!      if ( (il.eq.n) .and. (n.gt.1) .and.
!           ( (x-xa(n)) .lt. 0.001*(xa(n)-xa(n-1)) ) ) then
!        il = n-1
!      endif

      find_indexlgt = il

      return
      end
