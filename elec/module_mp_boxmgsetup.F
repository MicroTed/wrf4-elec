#ifdef WRF_ELEC

#ifndef RKIND
#define RKIND 4
#endif

#if RKIND == 4
#define FLOAT_PRECISION MPI_REAL
#else
#define FLOAT_PRECISION MPI_DOUBLE_PRECISION
#endif

#ifndef RELEASE
! multigrid solver for the potential poisson eq. (Dave Moulton)
#define BOXMG 
#endif

       MODULE module_mp_boxmgsetup
#ifdef DM_PARALLEL
      USE module_dm, ONLY : &
                  local_communicator, mytask, ntasks, ntasks_x, ntasks_y
#endif
#include        "BMG_constants.h"
      INCLUDE   'BMG_workspace.h'
      INCLUDE   'BMG_parameters.h'
!
!  Stuff for BOXMG solver (used for MPI)
!
! ------------------------------------------------
!      Multigrid/Workspace Memory Allocation:
! ------------------------------------------------

      ! Workspace pointers and logicals
      !
      INTEGER   BMG_pWORK(NBMG_pWORK)
      LOGICAL   BMG_InWORK(NBMG_InWORK)

      !
      ! Workspace pointer shift variables
      !
      INTEGER   pSI, pSR

! -------------------------------------------------
!      Multigrid:  Variables
! -------------------------------------------------

      !
      ! BoxMG Cycle and I/O Parameters
      !
      INTEGER   BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND    BMG_rPARMS(NBMG_rPARMS)
      LOGICAL   BMG_IOFLAG(NBMG_IOFLAG)

      !
      ! Workspace: Plane Relaxation
      !
      INTEGER, allocatable, dimension(:) ::  BMG_iWORK_PL
      REAL*RKIND,  allocatable, dimension(:) ::  BMG_rWORK_PL

      !
      ! Workspace: Generic
      !
      INTEGER, allocatable, dimension(:)  :: BMG_iWORK
      REAL*RKIND,  allocatable, dimension(:)  :: BMG_rWORK

      !
      ! Workspace: Coarse-grid Solve
      !
      INTEGER, allocatable, dimension(:)  ::   BMG_iWORK_CS
      REAL*RKIND,  allocatable, dimension(:)  ::   BMG_rWORK_CS

      !
      ! Solution, Source/RHS, and Stencil
      !
      REAL*RKIND, allocatable, dimension(:) :: Q, QF, SO

      !
      ! Miscellaneous
      !
      INTEGER   NCbmg, NCBW, NCI, NCU, NF, NOG, NOGc, NSO, NSOR
      REAL*RKIND    TOL, TOL_SAVE

!      INTEGER NGx, NGy, NGz
!      INTEGER   iGs, jGs, kGs, NGx, NGy, NGz, NLx, NLy, NLz
      INTEGER   NGx, NGy, NGz, NLx, NLy, NLz

      INTEGER, allocatable, dimension(:) ::  pMSG, pMSGSO

! -------------------------------------------------
!      MPI/MSG:  Variables
! -------------------------------------------------

      INTEGER   NBMG_MSG_iGRIDm, NBMG_MSG_iGRID

      INTEGER, allocatable, dimension(:) :: BMG_MSG_iGRID

      INTEGER  BMG_MSG_iGRIDdum
      INTEGER  BMG_MSG_pGRID(NBMG_MSG_pGRID)

      INTEGER   NMSGi, NMSGr

      integer :: nlxdg, nlydg, nlzdg, ngxdg, ngydg, ngzdg
      real*8  dt11,dt2
      
      integer, parameter :: maxbzone = 0

!     local arrays
      integer :: nbw,nbe,nbs,nbn
      integer :: nbwmg(1),nbemg(1),nbsmg(1),nbnmg(1)
      integer :: nnxs(1),nnys(1),nnzs(1)
      integer :: nxslm,nyslm,nzslm    !  dimensions for the SLM lightning

!     Took the necessary from ELEC_MODULE:

!   BoxMG control parameters (for MPI)

       integer ::  bmg_cycletype = 1 ! 0=FMG, 1=n-cycle
       integer ::  bmg_maxiter = 2
       integer ::  bmg_sep     = 0  ! whether to use separable mod or not: 0=BMG default, 1=simplified 1-D coeffs for fine grid
       integer ::  bmg_cgsolve   = 1  ! 0 = LU direct solver; 1 = boxmg_serial solve
       real    ::  bmg_tol     = 1.d-3

       integer :: igsdg, jgsdg, kgsdg

       integer :: NOGdg
       integer :: NOGcdg
       integer :: NFdg
       integer :: NCbmgdg
       integer :: NCIdg
       integer :: NSOdg
       integer :: NSORdg
       integer :: NCBWdg
       integer :: NCUdg
       integer :: NMSGidg
       integer :: NMSGrdg

       integer :: nlxlg, nlylg, nlzlg, ngxlg, ngylg, ngzlg
       integer :: igslg, jgslg, kgslg
       integer :: igslg0, jgslg0, kgslg0

       integer :: NOGlg
       integer :: NOGclg
       integer :: NFlg
       integer :: NCbmglg
       integer :: NCIlg
       integer :: NSOlg
       integer :: NSORlg
       integer :: NCBWlg
       integer :: NCUlg
       integer :: NMSGilg
       integer :: NMSGrlg

! BoxMG arrays for dynamics grid (dg):

       integer, parameter :: ip_RESdg = 3
       integer, parameter :: ip_SORdg = 5
       integer, parameter :: ip_CIdg  = 6
       integer, parameter :: ip_iGdg  = 10
       integer, parameter :: ip_MSGdg = 14
       integer, parameter :: ip_MSG_BUFdg = 16
       integer, parameter :: id_BMG3_SETUPdg = 52
       integer, parameter :: id_BMG3_MAX_ITERSdg = 61


      integer, parameter :: NBMG_pWORKdg = 18
      integer, parameter :: NBMG_InWORKdg = 7

      INTEGER   BMG_pWORKdg(NBMG_pWORKdg)
      LOGICAL   BMG_InWORKdg(NBMG_InWORKdg)

      !
      ! Workspace pointer shift variables
      !
      INTEGER   pSIdg, pSRdg


      integer, parameter :: NBMG_iPARMSdg = 72
      integer, parameter :: NBMG_rPARMSdg = 35
      integer, parameter :: NBMG_IOFLAGdg = 44

      INTEGER   BMG_iPARMSdg(NBMG_iPARMSdg)
      REAL      BMG_rPARMSdg(NBMG_rPARMSdg)
      LOGICAL   BMG_IOFLAGdg(NBMG_IOFLAGdg)

      INTEGER   NFmdg, NOGmdg, NSOmdg
      INTEGER   NBMG_iWORKdg, NBMG_rWORKdg
      INTEGER   NBMG_iWORK_PLdg, NBMG_rWORK_PLdg
      INTEGER   NBMG_iWORK_CSdg, NBMG_rWORK_CSdg



      REAL, allocatable, dimension(:) :: Qdg, QFdg, SOdg

      INTEGER, allocatable, dimension(:)  :: BMG_iWORKdg
      REAL, allocatable, dimension(:)  :: BMG_rWORKdg

      INTEGER, allocatable, dimension(:) ::  BMG_iWORK_PLdg
      REAL,  allocatable, dimension(:) ::  BMG_rWORK_PLdg

      INTEGER, allocatable, dimension(:)  ::   BMG_iWORK_CSdg
      REAL,  allocatable, dimension(:)  ::   BMG_rWORK_CSdg

      INTEGER, allocatable, dimension(:) ::  pMSGdg, pMSGSOdg


      INTEGER   NBMG_MSG_iGRIDmdg, NBMG_MSG_iGRIDdg

      INTEGER, allocatable, dimension(:) :: BMG_MSG_iGRIDdg
!      INTEGER, allocatable, dimension(:) :: BMG_MSG_iGRIDdum

      INTEGER, parameter :: NBMG_MSG_pGRIDdg = 8
      INTEGER  BMG_MSG_pGRIDdg(NBMG_MSG_pGRIDdg)

      INTEGER :: myproci,myprocj,myprock

      INTEGER :: nor


       CONTAINS

#if defined (BOXMG)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine mgsetupdg(nx,ny,nz,nor,   &
     &   nnxs,nnys,nnzs,  &
     &   nbw,nbe,nbs,nbn,  &
     &  iunit0)

!      USE ELEC_MODULE
!      USE COMMASMPI_MODULE
!      USE GRID_MODULE ! because use pointers
!      USE INDEX_MODULE ! only used for neelec and ipot (elec array indices ) - not used here.

#ifdef DM_PARALLEL
      USE module_dm, ONLY : &
                  local_communicator, mytask, ntasks, ntasks_x, ntasks_y, mytask_x, mytask_y
#endif
      IMPLICIT NONE

! ------------------------------------------------
!     Includes
! ------------------------------------------------

      INCLUDE   'mpif.h'

#include        "BMG_constants.h"
      INCLUDE   'BMG_workspace.h'
      INCLUDE   'BMG_parameters.h'

! ------------------------------------------------
!      Multigrid/Workspace Memory Allocation:
! ------------------------------------------------

      integer nx,ny,nz,nor

      integer nnxs(1),nnys(1),nnzs(1)
      integer nbw,nbe,nbs,nbn
      real dx,dy
      integer iunit0,iunit
      integer iestag
      logical lcalce   ! whether to calculate electric field components when iestag=1


      integer bcx, bcy

      integer i

!     Took the necessary from COMMASMPI_MODULE:
      INTEGER :: mpi_error_code

      INTEGER :: nprock = 1

!     HOW TO SET THOSE ?     


      !
      ! Workspace pointers and logicals

      INTEGER  BMG_MSG_iGRIDdum

      REAL*8    T1, T2, T_SOLVE_LOOP
      
      integer, allocatable :: nx_tiles(:), ny_tiles(:)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      iunit=0
      
      myproci = mytask_x+1
      myprocj = mytask_y+1
      myprock = 1

      IF ( NBMG_iPARMSdg .ne. NBMG_iPARMS .or.     &
     &     NBMG_rPARMSdg .ne. NBMG_rPARMS .or.      &
     &     NBMG_IOFLAGdg .ne. NBMG_IOFLAG .or.      &
     &      NBMG_pWORKdg  .ne. NBMG_pWORK  .or.    &
     &       NBMG_InWORKdg .ne. NBMG_InWORK .or.   &
     &        NBMG_MSG_pGRIDdg .ne. NBMG_MSG_pGRID) THEN

        write(0,*) 'Problem with BoxMG parameter constants!!'
!'
        STOP
      ENDIF

      BMG_iPARMSdg(id_BMG3_POINTERS) = BMG_USE_pointers

!      allocate(BMG_MSG_iGRIDdum(1:NBMG_MSG_iGRIDmdg))
! #ifdef __apple__
!      BMG_iPARMSdg(id_BMG3_CG_COMM) = BMG_CG_ALLGATHER
! #endif
      CALL WRF_SETUP_ProcGrid_file(     &
     &                 BMG_MSG_iGRIDdum, BMG_MSG_pGRIDdg,  &
     &                 BMG_iPARMSdg, NBMG_MSG_iGRIDmdg, NBMG_MSG_iGRIDdg, &
     &                 nbw,nbe,nbs,nbn,nnxs(1),nnys(1),nnzs(1),nx,ny,nz, &
     &                 mytask, ntasks, local_communicator, mpi_error_code &
     &                 )

      IF ( .not. allocated( BMG_MSG_iGRIDdg )) allocate(BMG_MSG_iGRIDdg(1:NBMG_MSG_iGRIDdg))
      
!      BMG_MSG_iGRIDdg(1:NBMG_MSG_iGRIDdg) = 0

      BMG_iPARMSdg(id_BMG3_POINTERS) = BMG_NO_pointers
      NBMG_MSG_iGRIDmdg = NBMG_MSG_iGRIDdg ! such that BMG_MSG_iGRIDdg(NBMG_MSG_iGRIDmdg)


      CALL WRF_SETUP_ProcGrid_file(  &
     &                 BMG_MSG_iGRIDdg, BMG_MSG_pGRIDdg,  &
     &                 BMG_iPARMSdg, NBMG_MSG_iGRIDmdg, NBMG_MSG_iGRIDdg,  &
     &                 nbw,nbe,nbs,nbn,nnxs(1),nnys(1),nnzs(1),nx,ny,nz,  &
     &                 mytask, ntasks, local_communicator, mpi_error_code  &
     &                 )


      NGxdg = BMG_MSG_iGRIDdg(id_BMG_MSG_NGx)
      NGydg = BMG_MSG_iGRIDdg(id_BMG_MSG_NGy)
      NGzdg = BMG_MSG_iGRIDdg(id_BMG_MSG_NGz)

      NLxdg = BMG_MSG_iGRIDdg(id_BMG_MSG_NLx)
      NLydg = BMG_MSG_iGRIDdg(id_BMG_MSG_NLy)
      NLzdg = BMG_MSG_iGRIDdg(id_BMG_MSG_NLz)

      iGsdg = BMG_MSG_iGRIDdg(id_BMG_MSG_iGs)
      jGsdg = BMG_MSG_iGRIDdg(id_BMG_MSG_jGs)
      kGsdg = BMG_MSG_iGRIDdg(id_BMG_MSG_kGs)

      iGslg0 = iGsdg - 1
      jGslg0 = jGsdg - 1
      kGslg0 = kGsdg - 1

      IF ( MyProcI /= BMG_MSG_iGRIDdg(id_BMG_MSG_MyProcI) .or.  &
     &      MyProcJ /= BMG_MSG_iGRIDdg(id_BMG_MSG_MyProcJ) .or.  &
     &      MyProcK /= BMG_MSG_iGRIDdg(id_BMG_MSG_MyProcK) ) THEN
           write(0,*) 'problem with myproci,j, or k! ', mytask
           write(0,*) mytask, MyProcI, BMG_MSG_iGRIDdg(id_BMG_MSG_MyProcI)
           write(0,*) mytask, MyProcJ, BMG_MSG_iGRIDdg(id_BMG_MSG_MyProcJ)
           write(0,*) mytask, MyProcK, BMG_MSG_iGRIDdg(id_BMG_MSG_MyProcK)
       ENDIF

!      IF ( mytask.EQ.0 ) THEN
        write(iunit,*) 'Global grid dims are (NGxdg,NGydg,NGzdg) = ',  &
     &       NGxdg,NGydg,NGzdg
        write(iunit,*) 'My local dims are NLxdg,NLydg,NLzdg = ',NLxdg,NLydg,NLzdg
        write(iunit,*) 'id_BMG_MSG_NLx,id_BMG_MSG_NLy,id_BMG_MSG_NLz = ',id_BMG_MSG_NLx,id_BMG_MSG_NLy,id_BMG_MSG_NLz
        write(iunit,*) 'iGsdg,jGsdg,kGsdg = ',iGsdg,jGsdg,kGsdg
        write(iunit,*) 'nx,ny,nz,nnxs,nnys = ',nx,ny,nz,nnxs(1),nnys(1)
!      ENDIF

      CALL EX_SETUP_BMG_parms(  &
     &               BMG_iPARMSdg, BMG_rPARMSdg, BMG_IOFLAGdg,  &
     &               mytask, ntasks, local_communicator, mpi_error_code  &
     &               )

      BMG_IOFLAGdg(iBMG3_BUG_RES_RELAX)    =  .FALSE.
      BMG_IOFLAGdg(iBMG3_BUG_RES_CG_SOLVE) =  .FALSE.
      BMG_IOFLAGdg(iBMG3_OUT_TIME_SETUP)   =  .FALSE.
      BMG_IOFLAGdg(iBMG3_OUT_TIME_CYCLING) =  .FALSE.
      BMG_IOFLAGdg(iBMG3_OUT_TIME_TOTAL)   =  .FALSE.

! separable fine GS speed-up doesn't seem to work with stretching for some reason,
! so make sure to use default if stretching is active.

! ------------------------------
!      Workspace Allocation:
! ------------------------------

      pSIdg=1
      pSRdg=1

      !
      ! Start with nothing stored in the workspace array
      !
      DO i=1, NBMG_InWORKdg
         BMG_InWORKdg(i) = .FALSE.
      ENDDO

      BMG_InWORKdg(i_InWORK_SO)  = .FALSE.    ! use a separate array for SO
      BMG_InWORKdg(i_InWORK_U)   = .FALSE.    ! use a separate array for Q
      BMG_InWORKdg(i_InWORK_Q)   = .FALSE.    ! use a separate array for QF
      BMG_InWORKdg(i_InWORK_RES) = .TRUE.     ! store RES in work array

      BMG_iPARMSdg(id_BMG3_POINTERS) = BMG_USE_pointers

      CALL BMG3_SymStd_SETUP_PtrWork(  &
     &          NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,  &
     &          iGsdg, jGsdg, kGsdg,  &
     &          BMG_iPARMSdg, BMG_rPARMSdg, BMG_IOFLAGdg,  &
     &          mytask, ntasks, ntasks_x, ntasks_y, NProcK,  &
     &          local_communicator,  &
     &          NOGmdg, NFmdg, NSOmdg,  &
     &          NBMG_iWORKdg, NBMG_rWORKdg,  &
     &          NBMG_iWORK_PLdg, NBMG_rWORK_PLdg,  &
     &          NBMG_iWORK_CSdg, NBMG_rWORK_CSdg,  &
     &          BMG_pWORKdg, BMG_InWORKdg, pSRdg, pSIdg  &
     &          )

      NOGdg  = BMG_iPARMSdg(id_BMG3_DIM_NOG)
      NFdg   = BMG_iPARMSdg(id_BMG3_DIM_NF)
      NCbmgdg   = BMG_iPARMSdg(id_BMG3_DIM_NC)
      NCIdg  = BMG_iPARMSdg(id_BMG3_DIM_NCI)
      NSOdg  = BMG_iPARMSdg(id_BMG3_DIM_NSO)
      NSORdg = BMG_iPARMSdg(id_BMG3_DIM_NSOR)
      NCBWdg = BMG_iPARMSdg(id_BMG3_DIM_NCBW)
      NCUdg  = BMG_iPARMSdg(id_BMG3_DIM_NCU)

      NMSGidg = BMG_iPARMSdg(id_BMG3_DIM_NMSGi)
      NMSGrdg = BMG_iPARMSdg(id_BMG3_DIM_NMSGr)

! ------------------------------
!     Report space requirements
! ------------------------------

!      IF (mytask .eq. 0) THEN
         WRITE(iunit,*)
         WRITE(iunit,*) 'Memory allocation in MB'
         WRITE(iunit,*) '======================='
         WRITE(iunit,*) 'size of real = ',RKIND
         WRITE(iunit,1210) 'SO           : ',NSOdg*float(RKIND)/(1024*1024)
         WRITE(iunit,1210) 'Q            : ',NFdg*float(RKIND)/(1024*1024)
         WRITE(iunit,1210) 'QF           : ',NFdg*float(RKIND)/(1024*1024)
         WRITE(iunit,1210) 'BMG_rWORK    : ',NBMG_rWORKdg*float(RKIND)/(1024*1024)
         WRITE(iunit,1210) 'BMG_iWORK    : ',NBMG_iWORKdg*4./(1024*1024)
         WRITE(iunit,1210) 'pMSG         : ',NBMG_pMSG*NOGdg*4./(1024*1024)
         WRITE(iunit,1210) 'pMSGSO       : ',NBMG_pMSG*NOGdg*4./(1024*1024)
         WRITE(iunit,1210) 'pLS          : ',NBMG_pLS*NOGdg*4./(1024*1024)
         WRITE(iunit,1210) 'BMG_iWORK_PL : ',NBMG_iWORK_PLdg*4./(1024*1024)
         WRITE(iunit,1210) 'BMG_rWORK_PL : ',NBMG_rWORK_PLdg*float(RKIND)/(1024*1024)
         WRITE(iunit,1210) 'BMG_iWORK_CS : ',NBMG_iWORK_CSdg*4./(1024*1024)
         WRITE(iunit,1210) 'BMG_rWORK_CS : ',NBMG_rWORK_CSdg*float(RKIND)/(1024*1024)
         WRITE(iunit,*)
!       ENDIF

 1210 FORMAT(4X,A,1X,F18.6,1X,'MBytes')

! ------------------------------
!     Allocate the space
! ------------------------------

       IF ( .not. allocated( SOdg ) ) THEN
       allocate(SOdg(1:NSOdg))
       allocate(Qdg(1:NFdg))
       allocate(QFdg(1:NFdg))
       allocate(BMG_rWORKdg(1:NBMG_rWORKdg))
       allocate(BMG_iWORKdg(1:NBMG_iWORKdg))
       allocate(BMG_iWORK_PLdg(1:NBMG_iWORK_PLdg))
       allocate(BMG_rWORK_PLdg(1:NBMG_rWORK_PLdg))
       allocate(BMG_iWORK_CSdg(1:NBMG_iWORK_CSdg))
       allocate(BMG_rWORK_CSdg(1:NBMG_rWORK_CSdg))
       allocate(pMSGdg(1:NBMG_pMSG*NOGdg))
       allocate(pMSGSOdg(1:NBMG_pMSG*NOGdg))
       ENDIF

         ! ------------------------------
         !     Zero the Stencil:
         ! ------------------------------

         DO i=1, NSOdg
            SOdg(i)=rZERO
         END DO

         DO i=1, NFdg
            Qdg(i)=rZERO
            QFdg(i)=rZERO
         END DO

         ! ------------------------------
         !     Zero the workspace:
         ! ------------------------------

         DO i=1, NBMG_rWORKdg
            BMG_rWORKdg(i)=rZERO
         END DO

         DO i=1, NBMG_iWORKdg
            BMG_iWORKdg(i)=iZERO
         END DO

         DO i=1, NBMG_rWORK_PLdg
            BMG_rWORK_PLdg(i)=rZERO
         END DO

         DO i=1, NBMG_iWORK_PLdg
            BMG_iWORK_PLdg(i)=iZERO
         END DO

         DO i=1, NBMG_rWORK_CSdg
            BMG_rWORK_CSdg(i)=rZERO
         END DO

         DO i=1, NBMG_iWORK_CSdg
            BMG_iWORK_CSdg(i)=iZERO
         END DO


      RETURN
      END subroutine mgsetupdg

! ==========================================================================
! ==========================================================================
! ==========================================================================

!      SUBROUTINE PUTF( SO, QF, QFW,Q,QW,sctot,nx,ny,nz,nor,  &
      SUBROUTINE PUTF( SO, QFW, Q, sctot,potfair,  &
     &                 nx,ny,nz,nor,                      &
     &                 nbw,nbe,nbs,nbn,                   &
     &                 NLx, NLy, NLz, NGx, NGy, NGz,      &
     &                 iGs, jGs, kGs, hx, hy, hz, isetup  &
     &               )
! ==========================================================================

!      USE ELEC_MODULE
!      USE COMMASMPI_MODULE
#ifdef DM_PARALLEL
      USE module_dm, ONLY : &
                  local_communicator, mytask, ntasks, ntasks_x, ntasks_y, mytask_x, mytask_y
#endif
      IMPLICIT NONE

      INTEGER :: nprock = 1

! ----------------
!     Includes
!
#include       "BMG_constants.h"
      INCLUDE  'BMG_stencils.h'

! ----------------
!     Parmeters
!
      REAL*RKIND   rHALF
      PARAMETER ( rHALF = 0.5D0 )

! ----------------
!     Arguments
!
      INTEGER  iGs, jGs, kGs, NLx, NLy, NLz, NGx, NGy, NGz, isetup
      REAL*RKIND   hx, hy,                   & 
     &         hz(0:NLx+1,NLz,0:NLy+1),      &
     &         Q(0:NLx+1,0:NLy+1,0:NLz+1),   &
     &         QFW(-nbw:nx+1+nbe,-nbs:ny+1+nbn,0:NLz+1),   &
     &         SO(0:NLx+2,0:NLy+2,0:NLz+2,4)


      integer nx,ny,nz,nor
      integer nbw,nbe,nbs,nbn
      integer, parameter    :: ng1 = 1

      real sctot(0:nx+1,nz,0:ny+1)
      real potfair(nz)

! ----------------
!     Local
!
      INTEGER  i,  iGf, j,  jGf, k,  kGf,  &
     &         iBEG, iEND, jBEG, jEND, kBEG, kEND

      integer ibegw,iendw,jbegw,jendw, ix,jy,kz
      REAL*RKIND   cm

      real eperao
      parameter (eperao  = 8.8592e-12 )

! ==========================================================================

! ----------------------------------------
!     Initalize constants
! ----------------------------------------

!         hz(0,:,:)=hz(1,:,:)
!         hz(NLx+1,:,:)=hz(NLx,:,:)
!         hz(:,:,0)=hz(:,:,1)
!         hz(:,:,NLy+1)=hz(:,:,NLy)

      !  Global indexing
      !
      iGf=iGs+NLx-1
      jGf=jGs+NLy-1
      kGf=kGs+NLz-1

!      write(*,*) 'igf,jgf,kgf = ',igf,jgf,kgf
!      write(*,*) 'NGx,NGy,NGz = ',NGx,NGy,NGz

      !  Loop Boundaries: iGs, jGs, kGs
      !
      IF ( iGs.EQ.1 ) THEN
         iBEG = 2
         ibegw = 1
      ELSE
         iBEG = 1
         ibegw = 1
      ENDIF

      IF ( jGs.EQ.1 ) THEN
         jBEG = 2
         jbegw = 1
      ELSE
         jBEG = 1
         jbegw = 1
      ENDIF

      IF ( kGs.EQ.1 ) THEN
         kBEG = 2
      ELSE
         kBEG = 1
      ENDIF

      !
      !  Loop Boundaries: iGf, jGf, kGf
      !
      IF ( iGf.EQ.NGx ) THEN
         iEND = NLx
         iendw = nx-1
      ELSE
         iEND = NLx+1
         iendw = nx
      ENDIF

      IF ( jGf.EQ.NGy ) THEN
         jEND = NLy
         jendw = ny-1
      ELSE
         jEND = NLy+1
         jendw = ny
      ENDIF

      IF ( kGf.EQ.NGz ) THEN
         kEND = NLz
      ELSE
         kEND = NLz+1
      ENDIF


      IF ( isetup == 1 ) THEN

! ----------------------------------------
!     Zero solution - actually initial guess !
! ----------------------------------------

      DO k=0, NLz+1
         DO j=0, NLy+1
            DO i=0, NLx+1

               Q(i,j,k)   = rZERO
!               Q(i,j,k) = potfair(i,k,j) ! 0 ! BMG_rand(r)

            END DO
         END DO
      END DO

! ----------------------------------------
!     Zero Stencil
! ----------------------------------------

      DO k=0, NLz+2
         DO j=0, NLy+2
            DO i=0, NLx+2

               SO(i,j,k,kp)  = rZERO
               SO(i,j,k,kpw) = rZERO
               SO(i,j,k,kps) = rZERO
               SO(i,j,k,kb)  = rZERO

            END DO
         END DO
      END DO

      ENDIF

!     SO       User defined real array that contains the coefficient
!              matrix. See above description for format.  If SETUP is
!              specified then coarse-grid operators are created and returned.
!
!     QFW       The user defined array that contains the right hand side.
!              It must be dimensioned to at least NFMAX.
!
!     Q        the user defined array that contains the solution vector.
!              It must be dimensioned to at least NFMAX.
!

! ----------------------------------------
!     Compute the source (right-hand-side) multiplied by the grid point volume
!    Note that this is actually -rhs, so multiply by -1
!    For electric poisson eq., then qf = (rho/epsilon)*(grid volume)
! ----------------------------------------

!      write(*,*) 'putf: rank, i1 = ',mytask,nbw,nbe,i1-1,nx+nbe-(-nbw+1)

      QFW=0.

      DO k=1,NLz 
         DO j=-nbs, ny+1+nbn
            DO i=-nbw, nx+1+nbe
               ix = Max( 1, Min(i,nx) )
               jy = Max( 1, Min(j,ny) )
               kz = Max( 1, Min(k,NLz-1) )
            IF ( k .ge. 1 .and. k .le. nz-1 .and.   &
     &           j .ge. 1 .and. j .le. ny  .and.    &
     &           i .ge. 1 .and. i .le. nx ) THEN
               QFW(i,j,k)    = (sctot(ix,kz,jy)/eperao)*hx*hy*hz(ix,kz,jy) ! actual potential* jacobian elements
            ELSE
               QFW(i,j,k)    = 0.
            ENDIF
!            IF ( .not. ( QFW(i,j,k) > -1.d30 .and. QFW(i,j,k) < 1.d30  ) ) THEN
!             write(0,*) 'bad value in QFW 1 = ',i,j,k,QFW(i,j,k)
!             write(0,*) 'ix,jy,kz,sctot,hz = ',ix,jy,kz,sctot(ix,kz,jy),hz(ix,kz,jy),hx,hy
!             STOP
!            ENDIF
            
            ENDDO
         ENDDO
      ENDDO

      IF ( isetup == 1 ) THEN  

! ----------------------------------------
!     Stencil coefficients
! ----------------------------------------
      !
      !  South
      !
      DO k=1, kEND
         DO j=jBEG,jEND
            DO i=1, iEND
               ix = Max( 1, Min(i,nx) )
               jy = Max( 1, Min(j,ny) )
               kz = Max( 1, Min(k,NLz-1) )
               SO(i,j,k,kps) = hx*hz(ix,kz,jy)/hy
            ENDDO
         ENDDO
      ENDDO

      !
      !  West
      !
      DO k=1, kEND
         DO j=1, jEND
            DO i=iBEG, iEND
               ix = Max( 1, Min(i,nx) )
               jy = Max( 1, Min(j,ny) )
               kz = Max( 1, Min(k,NLz-1) )
               SO(i,j,k,kpw) =  hy*hz(ix,kz,jy)/hx
            ENDDO
         ENDDO
      ENDDO

      !
      !  Bottom
      !
      DO k=kBEG, kEND
         DO j=1, jEND
            DO i=1, iEND
               ix = Max( 1, Min(i,nx) )
               jy = Max( 1, Min(j,ny) )
               kz = Max( 1, Min(k,NLz-1) )
               SO(i,j,k,kb) = hx*hy/hz(ix,kz,jy)
            ENDDO
         ENDDO
      ENDDO

! ----------------------------------------
!     Initial Guess
! ----------------------------------------

      DO k=1, NLz
         DO j=1, NLy
            DO i=1, NLz
!               Q(i,j,k) = potfair(i,k,j) ! 0 ! BMG_rand(r)
            ENDDO
         ENDDO
      ENDDO

! ---------------------------------------
!     Boundary Conditions
! ---------------------------------------

      !  Bottom Boundary
      !
      IF ( kGs.EQ.1 ) THEN

           DO k=1,1
            DO j=0, nly+1 ! jEND
              DO i=0, nlx+1 ! iEND
               Q(i,j,k-1) = 0.  ! set potential at one point LOWER
               ix = Max( 1, Min(i,nx) )
               jy = Max( 1, Min(j,ny) )
               kz = Max( 1, Min(k,NLz-1) )
               SO(i,j,k,kb) =  hx*hy/hz(ix,kz,jy)
              ENDDO
            ENDDO
           ENDDO

      ENDIF


      !
      !  Top boundary
      !
      IF ( kGf.EQ.NGz ) THEN

         DO k= NLz+1,NLz+1 
!           DO j=jBEG-1, jEND
!              DO i=iBEG-1, iEND
            DO j=0, nly+1
               DO i=0, nlx+1
                  Q(i,j,k) =  potfair(nlz-1)  ! potfair(i,k,j) 
               ix = Max( 1, Min(i,nx) )
               jy = Max( 1, Min(j,ny) )
               kz = Max( 1, Min(k,NLz-1) )
                 SO(i,j,k,kb) = SO(i,j,k,kb) + hx*hy/hz(ix,kz,jy)
              ENDDO
           ENDDO
         ENDDO

      ENDIF ! endif boundary stuff isetup=1


      DO k=1, NLz
         DO j=1, NLy
            DO i=1, NLx
               cm=rONE
               IF ( ( i.EQ.2 .AND. iGs.EQ.1 )  &
     &              .OR. ( i.EQ.NLx .AND. iGf.EQ.NGx ) ) THEN
                  cm = rHALF*cm
               ENDIF
               IF ( ( j.EQ.2 .AND. jGs.EQ.1 )  &
     &              .OR. ( j.EQ.NLy .AND. jGf.EQ.NGy ) ) THEN
                  cm = rHALF*cm
               ENDIF
               IF ( ( k.EQ.2 .AND. kGs.EQ.1 )  &
     &              .OR. ( k.EQ.NLz .AND. kGf.EQ.NGz ) ) THEN
                  cm = rHALF*cm
               ENDIF
               !
               SO(i,j,k,kp) = cm*SO(i,j,k,kp)  &
     &                      + SO(i,j,k,kps)+SO(i  ,j+1,k  ,kps)  &
     &                      + SO(i,j,k,kpw)+SO(i+1,j  ,k  ,kpw)  &
     &                      + SO(i,j,k,kb )+SO(i  ,j  ,k+1,kb )
               !
            ENDDO
         ENDDO
      ENDDO

      ENDIF
      
!      write(0,*) 'putf: pot2 (Q) : ',maxval(q),q(1,1,1)

!      DO k=0,NLz+1
!         DO j=-nbs,ny+1+nbn
!            DO i=-nbw,nx+1+nbe

!            IF ( .not. ( QFW(i,j,k) > -1.d30 .and. QFW(i,j,k) < 1.d30  ) ) THEN
!             write(0,*) 'bad value in QFW 2 = ',i,j,k,QFW(i,j,k)
!             STOP
!            ENDIF
!            ENDDO
!         ENDDO
!      ENDDO

! ==========================================================================

      RETURN
      END SUBROUTINE PUTF
! ==========================================================================
!            rho = qf(i,j,k)*eperao/(hx*hy*hz(i,k,j))
! ==========================================================================
! ==========================================================================

      SUBROUTINE PUTPHI( pot2,potwrf,elec,nx,ny,nz,nor,  &
     &               nbw,nbe,nbs,nbn,iunit,  &
     &                 NLx, NLy, NLz, NGx, NGy, NGz,  &
     &                 iGs, jGs, kGs,dx,dy,dz  &
     &               )
! ==========================================================================

      IMPLICIT NONE
!
! ----------------
!     Arguments

      INTEGER  iGs, jGs, kGs, NLx, NLy, NLz, NGx, NGy, NGz, kg
!     Local
      REAL*RKIND   &
     &         electmp(-nbw:nx+1+nbe,-nbs:ny+1+nbn,nz,4),  &  
     &         dztmp(-nbw:nx+1+nbe,-nbs:ny+1+nbn,nz)  

!!    solution that is  passed in - pot2=Qdg 
      real*RKIND pot2(-nbw:nx+1+nbe,-nbs:ny+1+nbn,0:nz+1)   ! must be same size of Qdg=Q. 


!      FINAL SOL on the WRF IKJ GRID

      real elec(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn,4)
      real potwrf(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)
      real dz(-nbw:nx+1+nbe,nz,-nbs:ny+1+nbn)

      integer nx,ny,nz,nor
      integer nbw,nbe,nbs,nbn
      integer iunit
      real dx,dy!, dz(-nbw:nx+1+nbe,-nbs:ny+1+nbn,0:NLz+1)

!     Local
!
      INTEGER  i, j, k,l, ix,jy,kz

! ==========================================================================


! calculate E 

          electmp(:,:,:,:)=0
!         might want to need setting up those from -nbs etc
!          dz(0,:,:)=dz(1,:,:)
!          dz(NLx+1,:,:)=dz(NLx,:,:)
!          dz(:,:,0)=dz(:,:,1)
!          dz(:,:,NLy+1)=dz(:,:,NLy)

!          do l=1,4
!          electmp(0,:,:,l)=electmp(3,:,:,l)
!          electmp(NLx+1,:,:,l)=electmp(NLx-3,:,:,l)
!          electmp(:,0,:,l)=electmp(:,3,:,l)
!          electmp(:,NLy+1,:,l)=electmp(:,NLy-3,:,l)
!          enddo

      do k = 1,nz
      DO j=-nbs, ny+1+nbn
      DO i=-nbw, nx+1+nbe
               ix = Max( 1, Min(i,nx) )
               jy = Max( 1, Min(j,ny) )
               kz = Max( 1, Min(k,NLz-1) )
          dztmp(i,j,k)   = dz(ix,kz,jy)
      end do
      end do
      end do
!      write(0,*) 'efieldstag',NLx,NLy,NLz,nx,ny,nz
!      write(0,*) 'pot2 max = ', Maxval( pot2 ), pot2(1,1,1)

        call efieldstag(  &
     &   nx,ny,nz,nor  &
     &  ,electmp,dx,dy,dztmp   &
     &  ,nbw,nbe+1,nbs,nbn+1 &  
     &  ,pot2 &
     &  ,iunit)

!      write(0,*) 'pot2 max2 = ', Maxval( pot2 ), pot2(1,1,1)

      do k = 1,nz
      DO j=-nbs, ny+1+nbn
      DO i=-nbw, nx+1+nbe
        do l = 1,4
            elec(i,k,j,l) = electmp(i,j,k,l)
        end do
        potwrf(i,k,j)   = pot2(i,j,k)
      end do
      end do
      end do


      RETURN
      END SUBROUTINE PUTPHI

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  START OF EFIELDSTAG: Calculates efield on staggered C-grid
!
      subroutine efieldstag(  &
     &  nx,ny,nz,nor &
     &  ,elec,dx,dy,dz &
     &  ,nbw,nbe,nbs,nbn & !  ,nnz &
     &  ,phi0  &
     &  ,iunit)  
     
#ifdef DM_PARALLEL
      USE module_dm, ONLY : &
                  local_communicator, mytask, ntasks, ntasks_x, ntasks_y
#endif
!
      implicit none
!
      integer nx,ny,nz,nor

      integer    nbw,nbe,nbs,nbn ! ,nnz

      real phi0(-nbw:nx+nbe,-nbs:ny+nbn,0:nz+1)
      real elec(-nbw:nx+nbe,-nbs:ny+nbn,nz,4) ! 3 components + mag

      integer ix,jy,kz,n,i,j,k
      integer iunit
      integer ibg  ! =1 for phi=0 at k=1, -1 for phi=0 at k=0
      
      integer imx,jmx,kmx, imx2, jmx2, kmx2
!
      real emag
      real emax, emax2
      real emin

!      real dz(-nbw:nx+nbe,-nbs:ny+nbn,0:nnz-1)
      real dx,dy,dz(0:nx+1,0:ny+1,1:nz)

!
!  compute electric field - three components

     elec(:,:,:,:)=0. 
!
!  x-component
!
!      if ( nx .gt. 2 ) then
      if ( nx .gt. 1 ) then
!
!      elec(:,ny,:,1) = 0.0
!      elec(:,:,nz,1) = 0.0
      elec(:,ny+1,:,1) = 0.0
      elec(:,:,nz,1) = 0.0

      do kz = 1,nz-1
       do jy = 1,ny
        do ix = 1,nx
       
         elec(ix,jy,kz,1) = -(phi0(ix,jy,kz) - phi0(ix-1,jy,kz))/dx

        end do
       end do
      end do
!
      end if
!
!  y-component
!
!      if ( ny .gt. 2 ) then
      if ( ny .gt. 1 ) then
!
!      elec(nx,:,:,2) = 0.0
      elec(nx+1,:,:,2) = 0.0
      elec(:,:,nz,2) = 0.0

      do kz = 1,nz-1
       do jy = 1,ny
        do ix = 1,nx
      
         elec(ix,jy,kz,2) = -(phi0(ix,jy,kz) - phi0(ix,jy-1,kz))/dy
     
        end do
       end do
      end do
!
      ELSE
       elec(:,:,:,2) = 0.0
      end if
!
!  z-component
!
!      if ( nz .gt. 2 ) then
      if ( nz .gt. 1 ) then
!

!      elec(nx,:,:,3) = 0.0
!      elec(:,ny,:,3) = 0.0
      elec(nx+1,:,:,3) = 0.0
      elec(:,ny+1,:,3) = 0.0

!      do kz = 2,nz
!       do jy = 1,ny-1
!       do ix = 1,nx-1
       do kz = 1,nz
       do jy = 1,ny
       do ix = 1,nx
      
         elec(ix,jy,kz,3) = -(phi0(ix,jy,kz) - phi0(ix,jy,kz-1))/dz(ix,jy,kz)
      
        end do
       end do
      end do
!
      endif

!
! 1st order difference
!      elec()%flt3d(ix,jy,01,iez) = -elec()%flt3d(ix,jy,01,ipot)*dzi*
!     :   gt(ix,jy,01,imapz)

! 2nd order difference
!      elec()%flt3d(ix,jy,01,iez) = -elec()%flt3d(ix,jy,02,ipot)*dzi2*
!     :   gt(ix,jy,01,imapz)


! 2nd order difference -- not working right?
!      elec()%flt3d(ix,jy,01,iez) = 
!     :    -0.5*gt(ix,jy,2,imapz)/dz* 
!     :          (  -3.0*elec()%flt3d(ix,jy,1,ipot)  
!     :              + 4*elec()%flt3d(ix,jy,2,ipot) 
!     :              -   elec()%flt3d(ix,jy,3,ipot)   )

!      elec()%flt3d(ix,jy,nz-kstag,iez) = 
!     > -(elec()%flt3d(ix,jy,nz-kstag,ipot)-elec()%flt3d(ix,jy,nz-2*kstag,ipot))*dzi*
!     :   gt(ix,jy,nz-kstag,imapz)

!
!  e-magnitude
!

      do kz = 1,nz !nz-1
      do jy = 1,ny !ny-1
      do ix = 1,nx !nx-1

        elec(ix,jy,kz,4) = &
!     & Sqrt( 0.25*(( elec(ix,jy,kz,1) + elec(ix+1,jy ,kz  ,1) )**2 &
!     &     +       ( elec(ix,jy,kz,2) + elec(ix ,jy+1,kz  ,2) )**2 &
!!     &     +       ( elec(ix,jy,kz,3) + elec(ix ,jy ,kz+1,3) )**2 )) 
!     &     +       ( elec(ix,jy,kz,3) + elec(ix ,jy ,kz,3) )**2 )) 
     & Sqrt(  ( elec(ix,jy,kz,1) )**2 &
     &     +  ( elec(ix,jy,kz,2) )**2 &
     &     +  ( elec(ix,jy,kz,3) )**2 ) 

      end do
      end do
      end do
!

!      IF ( mytask == 0 ) THEN
!        write(iunit,*) 'EFIELDSTAG--MAG-MAX=',MAXVAL(elec(1:nx,1:ny,1:nz,4)),MAXLOC(elec(1:nx,1:ny,1:nz,4))
!      ENDIF
!
!
      return
      end subroutine efieldstag
!
!
!  END OF EFIELD
!
!

      END MODULE module_mp_boxmgsetup

! ###############################################################
! ###############################################################
! ###############################################################

      SUBROUTINE WRF_SETUP_ProcGrid_file(  &
     &                BMG_MSG_iGRID, BMG_MSG_pGRID,  &
     &                BMG_iPARMS, NBMG_MSG_iGRIDm, NBMG_MSG_iGRID,  &
     &                nbw,nbe,nbs,nbn,nnxs,nnys,nnzs,NLx,NLy,NLz,  &
     &                MPI_MyProc, NProc, BMG_MPI_COMM, MPI_IERROR  &
     &                )

! ==========================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG3_SymStd_SETUP_ProcGrid_file reads and communicates minimal
!     data from a file that describes the processor grid and creates
!     related information for use locally.
!
! $license_flag$

!      USE COMMASMPI_MODULE
#ifdef DM_PARALLEL
      USE module_dm, ONLY : &
                  local_communicator, mytask, ntasks, ntasks_x, ntasks_y, mytask_x, mytask_y
#endif

      IMPLICIT NONE

! -----------------------------
!     Includes
!
      INCLUDE 'mpif.h'
      INCLUDE   'BMG_workspace.h'
      INCLUDE   'BMG_parameters.h'

! -----------------------------
!     Argument Declarations
!
      INTEGER  BMG_iPARMS(NBMG_iPARMS)

      INTEGER  NBMG_MSG_iGRIDm, NBMG_MSG_iGRID
      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRIDm),  &
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)

      INTEGER   BMG_MPI_COMM, MPI_MyProc, MPI_IERROR! , NProc

      integer  nnxs,nnys,nnzs
      integer  nbw,nbe,nbs,nbn

! -----------------------------
!     Local Declarations
!
      INTEGER  CE, GRIDDIST, i, iGs, j, jGs, k, kGs,  &
!     &         MyProc, & !  MyProcI, MyProcJ, MyProcK, &  ! ditto ??
     &         NGx, NGy, NGz, NLx, NLy, NLz,  &
     &         pSI, & 
     &         p_iGs, p_jGs, p_kGs, p_NLx, p_NLy, p_NLz,  &
     &         p_PC, p_PG, p_ProcCoord, p_ProcGrid

      INTEGER :: NProcK,NProcJ,NProcI,NProc,MyProcI,MyProcJ,MyProcK,MyProc

      LOGICAL  exist

      integer, allocatable :: nx_tiles(:), ny_tiles(:)
      integer, allocatable :: nx_tilesin(:), ny_tilesin(:)

      INTEGER :: mpi_error_code = 1      

! ==========================================================================

       NProcK=1
       NProcI=ntasks_x
       NProcJ=ntasks_y
       NProc=ntasks
       MYProcI=mytask_x+1
       MYProcJ=mytask_y+1
       MYProcK=NProcK
       
       allocate(nx_tilesin(NProcI))
       allocate(ny_tilesin(NProcJ))
       allocate(nx_tiles(NProcI))
       allocate(ny_tiles(NProcJ))
       
       nx_tilesin(:) = 0
       ny_tilesin(:) = 0
       
       nx_tilesin(mytask_x+1) = nlx
       ny_tilesin(mytask_y+1) = nly
       
#ifdef DM_PARALLEL
        write(0,*) 'NProcJ = ',NProcJ
        DO j = 1,ntasks_y
          write(0,*) 'j,ny_tiles = ',j,ny_tilesin(j)
        ENDDO

        CALL MPI_AllReduce(nx_tilesin, nx_tiles, NProcI, MPI_INTEGER, MPI_MAX, local_communicator, mpi_error_code)
        CALL MPI_AllReduce(ny_tilesin, ny_tiles, NProcJ, MPI_INTEGER, MPI_MAX, local_communicator, mpi_error_code)

        write(0,*) 'NProcJ = ',NProcJ
        DO j = 1,ntasks_y
          write(0,*) 'j,ny_tiles = ',j,ny_tiles(j)
        ENDDO

        write(0,*) 'NProcI = ',NProcI
        DO j = 1,ntasks_x
          write(0,*) 'j,nx_tiles = ',j,nx_tiles(j)
        ENDDO
#endif
      
! -------------------------------------------
!     MPI -> MSG indexing:
! -------------------------------------------

      MyProc = MPI_MyProc+1

      GRIDDIST = BMG_MSG_GRIDDIST_CUSTOM

! -------------------------------------------
!     Open the data file, read the header
! -------------------------------------------

      !
      !  Master:
      !
      IF ( MyProc .EQ. 1 ) THEN

!      write(0,*) 'INFO', NProcK,NProcJ,NProcI,NProc,MyProcI,MyProcJ,MyProcK,MyProc
         !
         ! Eliminate trailing blanks
         !

      ENDIF

! --------------------------------------------
!     Process and distribute the grid data
! --------------------------------------------

      IF ( GRIDDIST.EQ.BMG_MSG_GRIDDIST_AUTOMATIC ) THEN

         !
         !  This section would take NGx, NGy, NGz and
         !  distribute them as evenly as possible on
         !  the ntasks_x x ntasks_y x NProcK processor grid.
         !

!
      ELSEIF ( GRIDDIST.EQ.BMG_MSG_GRIDDIST_CUSTOM  &
     &        .OR. GRIDDIST.EQ.BMG_MSG_GRIDDIST_UNIFORM ) THEN
         !
         !  Master:
         !
!         IF ( MyProc.EQ.1 ) THEN
            !
!             Read the number of processors in i and j
! ALL THIS WAS ALSO COMMENTED OUT ??
            !
!            READ(10,*) ntasks_x, ntasks_y, NProcK
            !
!         ENDIF

         !
!         All:
         !
!         CALL MPI_Bcast( ntasks_x, 1, MPI_INTEGER,  &
!     &                   0, BMG_MPI_COMM, MPI_IERROR )
!
!         CALL MPI_BCast( ntasks_y, 1, MPI_INTEGER,  &
!     &                   0, BMG_MPI_COMM, MPI_IERROR )
!
!!         CALL MPI_BCast( NProcK, 1, MPI_INTEGER,  &
!         CALL MPI_BCast( 1, 1, MPI_INTEGER,  &
!     &                   0, BMG_MPI_COMM, MPI_IERROR )
!
!         IF ( ntasks_x*ntasks_y*1 .NE. NProc ) THEN
!            WRITE(*,*) 'Procgrid Error : ',   &
!     &                 ntasks_x,' x ', ntasks_y, 'x 1 not eq to ', NProc
!            STOP
!         ENDIF
!   UP TO HERE
         !
         ! Compute storage requirements and pointers
         !
         pSI = NBMG_MSG_iGRID_data + 1
         !
         BMG_MSG_pGRID(ip_BMG_MSG_NLx_Grid)  = pSI
         pSI = pSI + NProcI
         !
         BMG_MSG_pGRID(ip_BMG_MSG_NLy_Grid)  = pSI
         pSI = pSI + NProcJ
         !
         BMG_MSG_pGRID(ip_BMG_MSG_NLz_Grid)  = pSI
         pSI = pSI + NProcK
         !
         BMG_MSG_pGRID(ip_BMG_MSG_iGs_Grid)  = pSI
         pSI = pSI + NProcI
         !
         BMG_MSG_pGRID(ip_BMG_MSG_jGs_Grid)  = pSI
         pSI = pSI + NProcJ
         !
         BMG_MSG_pGRID(ip_BMG_MSG_kGs_Grid)  = pSI
         pSI = pSI + NProcK
         !
         BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)  = pSI
         pSI = pSI + NProcI*NProcJ*NProcK ! or just  + nstasks ?
         !
         BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord) = pSI
         pSI = pSI + 3*NProcI*NProcJ*NProcK
         !
         NBMG_MSG_iGRID = pSI-1

         IF ( BMG_iPARMS(id_BMG3_POINTERS).EQ.BMG_USE_pointers ) THEN

            NBMG_MSG_iGRIDm = NBMG_MSG_iGRID

!            CLOSE(10)
            RETURN

         ELSEIF ( BMG_iPARMS(id_BMG3_POINTERS).EQ.BMG_NO_pointers ) THEN

            IF ( NBMG_MSG_iGRID.GT.NBMG_MSG_iGRIDm ) THEN

               IF (MyProc.EQ.1) THEN
                  WRITE(*,500) 'Space for BMG_MSG_iGRID is too small!'
                  WRITE(*,510) 'HAVE: NBMG_MSG_iGRIDm = ',  &
     &                         NBMG_MSG_iGRIDm
                  WRITE(*,520) 'NEED: NBMG_MSG_iGRIDm = ',  &
     &                         NBMG_MSG_iGRID
               ENDIF

               CALL MPI_FINALIZE(MPI_IERROR)
               STOP

            ENDIF

         ELSE

            IF (MyProc.EQ.1) THEN
               WRITE(*,500) 'Value of BMG_iPARMS(id_BMG3_pointers) is '  &
     &                      //'out of range ...'
               WRITE(*,510) 'HAVE: BMG_iPARMS(id_BMG3_POINTERS) = ',  &
     &                      BMG_iPARMS(id_BMG3_POINTERS)
               WRITE(*,505) 'NEED: BMG_iPARMS(id_BMG3_POINTERS) = '  &
     &                      //'0 or 1'
            ENDIF

            CALL MPI_FINALIZE(MPI_IERROR)
            STOP

         ENDIF

         !
         !  Make local copies of pointers
         !
         p_NLx = BMG_MSG_pGRID(ip_BMG_MSG_NLx_Grid)
         p_NLy = BMG_MSG_pGRID(ip_BMG_MSG_NLy_Grid)
         p_NLz = BMG_MSG_pGRID(ip_BMG_MSG_NLz_Grid)
         p_iGs = BMG_MSG_pGRID(ip_BMG_MSG_iGs_Grid)
         p_jGs = BMG_MSG_pGRID(ip_BMG_MSG_jGs_Grid)
         p_kGs = BMG_MSG_pGRID(ip_BMG_MSG_kGs_Grid)
         !
         p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)
         p_ProcCoord = BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord)

         
         
         
         !
         !  Master:  ! 
         !
         IF ( MyProc.EQ.1 ) THEN
            !
            IF ( GRIDDIST.EQ.BMG_MSG_GRIDDIST_CUSTOM ) THEN
               !
               ! Read arrays of fine-grid local dimensions
               !
               IF ( NProcI == 1 ) THEN
                 BMG_MSG_iGRID(p_NLx) = nnxs
               ELSE
               DO i=1, NProcI
                  IF ( i == 1 ) THEN
                  BMG_MSG_iGRID(p_NLx+i-1) = nx_tiles(i) + nbw 
                  ELSEIF ( i == NProcI ) THEN  
                  BMG_MSG_iGRID(p_NLx+i-1) = nx_tiles(i) + nbw
                  ELSE
                  BMG_MSG_iGRID(p_NLx+i-1) = nx_tiles(i)
                  ENDIF
               
               ENDDO
               ENDIF
               !

               IF ( NProcJ == 1 ) THEN
                 BMG_MSG_iGRID(p_NLy) = nnys
               ELSE

               DO j=1, NProcJ
                 IF ( j == 1 ) THEN
                  BMG_MSG_iGRID(p_NLy+j-1) = ny_tiles(j) + nbs

                 ELSEIF ( j == NProcJ ) THEN
                  BMG_MSG_iGRID(p_NLy+j-1) = ny_tiles(j) + nbs

                 ELSE
                  BMG_MSG_iGRID(p_NLy+j-1) = ny_tiles(j)
                 ENDIF
               ENDDO
               ENDIF
               !


               IF ( NProcK == 1 ) THEN
                 BMG_MSG_iGRID(p_NLz) = nnzs !- 1
               ELSE

               DO k=1, NProcK
                 IF ( k == 1 ) THEN
                  BMG_MSG_iGRID(p_NLz+k-1) = NLz+1 ! does there need to be an extra point below ground??

                 ELSEIF ( k == NProcK ) THEN
                  BMG_MSG_iGRID(p_NLz+k-1) = NLz +nbs

                 ELSE
                  BMG_MSG_iGRID(p_NLz+k-1) = NLz
                 ENDIF
               ENDDO
               ENDIF
! WHY THIS COMMENTED ?
!               DO k=1, NProcK
!                  BMG_MSG_iGRID(p_NLz+k-1) = nnzs-1
!               ENDDO
               !
               ! Close the grid file.
               !
!               CLOSE (10)
            ELSE
               !
!              why was this COMMENTED ?
               ! Read fine-grid local dimensions
               !
!               READ(10,*) NLx, NLy, NLz
               !
               ! Close the grid file.
               !
!               CLOSE (10)
               !
               ! Fill arrays of fine-grid local dimensions
               !
               DO i=1, NProcI
                  BMG_MSG_iGRID(p_NLx+i-1) = NLx
               ENDDO
               !
               DO j=1, NProcJ
                  BMG_MSG_iGRID(p_NLy+j-1) = NLy
               ENDDO
               !
               DO k=1, NProcK
                  BMG_MSG_iGRID(p_NLz+k-1) = NLz
               ENDDO
               !
            ENDIF

         ENDIF

         !
         !  All:
         !
         CALL MPI_Bcast( BMG_MSG_iGRID(p_NLx), NProcI, MPI_INTEGER,  &
     &                   0, BMG_MPI_COMM, MPI_IERROR )

         CALL MPI_BCast( BMG_MSG_iGRID(p_NLy), NProcJ, MPI_INTEGER,  &
     &                   0, BMG_MPI_COMM, MPI_IERROR )

         CALL MPI_BCast( BMG_MSG_iGRID(p_NLz), NProcK, MPI_INTEGER,  &
     &                   0, BMG_MPI_COMM, MPI_IERROR )


!              why is this COMMENTED ?
         !
         ! Setup processor grid:
         !
         !   MyProc  range is [1,ntasks_x*ntasks_y*NProcK]
         !   MyProcI range is [1,ntasks_x]
         !   MyProcJ range is [1,ntasks_y]
         !   MyProcJ range is [1,NProcK]
         !
!         MyProcK = (MyProc-1)/(ntasks_x*ntasks_y) + 1
!         MyProcJ = (MyProc-(MyProcK-1)*ntasks_x*ntasks_y-1)/ntasks_x + 1
!         MyProcI = MyProc
!     &           - (MyProcK-1)*ntasks_x*ntasks_y - (MyProcJ-1)*ntasks_x

         p_PG = p_ProcGrid
         IF ( NProcK > 1 ) THEN
            write(0,*) 'STOP! Cannot handle NProcK > 1 !!'
!            call WRFmpi_shutdown()
            STOP
         ENDIF
         DO k=1, NProcK
            DO j=1, NProcJ
               DO i=1, NProcI
!                  BMG_MSG_iGRID(p_PG) = procmap(i,j) + 1
                  BMG_MSG_iGRID(p_PG) =  i + (j-1)*NProcI  &
     &                                    + (k-1)*NProcJ*NProcI
!                  IF (  BMG_MSG_iGRID(p_PG) .ne. procmap(i,j)+1 ) THEN
!                    write(0,*) 'BMG_MSG_iGRID, procmap: ',BMG_MSG_iGRID(p_PG), procmap(i,j)+1
!                  ENDIF

                  p_PG = p_PG + 1
               ENDDO
            ENDDO
         ENDDO

         !
         ! Store coordinates of each processor
         !
         p_PC = p_ProcCoord
         DO k=1, NProcK
            DO j=1, NProcJ
               DO i=1, NProcI
                  BMG_MSG_iGRID(p_PC)   = i
                  BMG_MSG_iGRID(p_PC+1) = j
                  BMG_MSG_iGRID(p_PC+2) = k
                  p_PC = p_PC + 3
               ENDDO
            ENDDO
         ENDDO

         !
         !  Global problem size
         !
         NGx = 0
         DO i=1, NProcI
            NGx = NGx + BMG_MSG_iGRID(p_NLx+i-1)
         ENDDO

         NGy = 0
         DO j=1, NProcJ
            NGy = NGy + BMG_MSG_iGRID(p_NLy+j-1)
         ENDDO

         NGz = 0
         DO k=1, NProcK
            NGz = NGz + BMG_MSG_iGRID(p_NLz+k-1)
         ENDDO

         !
         ! Global index (i,j) of the lower left corner
         !
         iGs = 1
         DO i=1, NProcI
            BMG_MSG_iGRID(p_iGs+i-1) = iGs
            iGs = iGs + BMG_MSG_iGRID(p_NLx+i-1)
         ENDDO

         jGs = 1
         DO j=1, NProcJ
            BMG_MSG_iGRID(p_jGs+j-1) = jGs
            jGs = jGs + BMG_MSG_iGRID(p_NLy+j-1)
         ENDDO

         kGs = 1
         DO k=1, NProcK
            BMG_MSG_iGRID(p_kGs+k-1) = kGs
            kGs = kGs + BMG_MSG_iGRID(p_NLz+k-1)
         ENDDO

         !
         !  Sanity check
         !
         IF ( NGx.NE.iGs-1 .OR. NGy.NE.jGs-1 .OR. NGz.NE.kGs-1 ) THEN

            IF ( MyProc.EQ.1 ) THEN
               WRITE(*,*) ' THIS IS STRANGE: NGx = ', NGx
               WRITE(*,*) '                  iGs = ', iGs
               WRITE(*,*) '                  NGy = ', NGy
               WRITE(*,*) '                  jGs = ', jGs
               WRITE(*,*) '                  NGz = ', NGz
               WRITE(*,*) '                  kGs = ', kGs
            ENDIF

            CALL MPI_FINALIZE(MPI_IERROR)
            STOP

         ENDIF

      ENDIF
      
       deallocate(nx_tilesin)
       deallocate(ny_tilesin)
       deallocate(nx_tiles)
       deallocate(ny_tiles)



! --------------------------------------------
!     Packup the local grid data
! --------------------------------------------

      !
      !  Store: NGx, NGy, NGz, NLx, NLy, NLz, iGs, jGs, kGs
      !
      BMG_MSG_iGRID(id_BMG_MSG_NGx) = NGx
      BMG_MSG_iGRID(id_BMG_MSG_NGy) = NGy
      BMG_MSG_iGRID(id_BMG_MSG_NGz) = NGz
      BMG_MSG_iGRID(id_BMG_MSG_NLx) = BMG_MSG_iGRID(p_NLx+MyProcI-1)
      BMG_MSG_iGRID(id_BMG_MSG_NLy) = BMG_MSG_iGRID(p_NLy+MyProcJ-1)
      BMG_MSG_iGRID(id_BMG_MSG_NLz) = BMG_MSG_iGRID(p_NLz+MyProcK-1)
      BMG_MSG_iGRID(id_BMG_MSG_iGs) = BMG_MSG_iGRID(p_iGs+MyProcI-1)
      BMG_MSG_iGRID(id_BMG_MSG_jGs) = BMG_MSG_iGRID(p_jGs+MyProcJ-1)
      BMG_MSG_iGRID(id_BMG_MSG_kGs) = BMG_MSG_iGRID(p_kGs+MyProcK-1)
      !
      !  Store: Nproc, ntasks_x, ntasks_y, NProcK
      !
      BMG_MSG_iGRID(id_BMG_MSG_NProc)   = NProc
      BMG_MSG_iGRID(id_BMG_MSG_NProcI)  = NProcI
      BMG_MSG_iGRID(id_BMG_MSG_NProcJ)  = NProcJ
      BMG_MSG_iGRID(id_BMG_MSG_NProcK)  = NProcK
      !
      !  Store: MyProc, MyProcI, MyProcJ, MyProcK
      !
      BMG_MSG_iGRID(id_BMG_MSG_MyProc)  = MyProc
      BMG_MSG_iGRID(id_BMG_MSG_MyProcI) = MyProcI
      BMG_MSG_iGRID(id_BMG_MSG_MyProcJ) = MyProcJ
      BMG_MSG_iGRID(id_BMG_MSG_MyProcK) = MyProcK 
      !
      !  Store/set communicators:
      !
      BMG_MSG_iGRID(id_BMG_MSG_COMM)    = BMG_MPI_COMM
      BMG_MSG_iGRID(id_BMG_MSG_COMM_xy) = 0
      BMG_MSG_iGRID(id_BMG_MSG_COMM_yz) = 0
      BMG_MSG_iGRID(id_BMG_MSG_COMM_xz) = 0

! ==========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SymStd_SETUP_ProcGrid_file.f',  &
     &        //,5X,A,/)
 505  FORMAT (5X,A,/)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

! ===========================================

      RETURN
      END subroutine WRF_SETUP_ProcGrid_file

! ==========================================================================
! ==========================================================================
! ==========================================================================
! ==========================================================================

      SUBROUTINE EX_SETUP_BMG_parms(  &
     &              BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,  &
     &              MPI_MyProc, NProc, EX_MPI_COMM, EX_MPI_IERROR  &
     &              )

! ==========================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     EX_SETUP_BMG_parms on the master ( MPI_MyProc=0 ) first sets the
!     default cycle parameters, then reads customization parameters
!     from a data file, and then broadcasts the arrays.
!
! =======================================================================
! $license_flag$
! =======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   INPUT/OUTPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   OUTPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   LOCAL:
!  --------------------
!
!
! ==========================================================================

      USE module_mp_boxmgsetup, Only: bmg_cycletype, bmg_maxiter, bmg_sep, bmg_cgsolve, bmg_tol


      IMPLICIT NONE
!   BoxMG control parameters (for MPI)

!       integer ::  bmg_cycletype = 1 ! 0=FMG, 1=n-cycle
!       integer ::  bmg_maxiter = 2
!       integer ::  bmg_sep     = 0  ! whether to use separable mod or not: 0=BMG default, 1=simplified 1-D coeffs for fine grid
!       integer ::  bmg_cgsolve   = 1  ! 0 = LU direct solver; 1 = boxmg_serial solve
!       real    ::  bmg_tol     = 1.d-3

! -----------------------------
!     Includes
! -----------------------------

      INCLUDE 'mpif.h'

#include      "BMG_constants.h"
      INCLUDE 'BMG_parameters.h'

! ----------------------------
!     Argument Declarations
! ----------------------------

      !
      ! BoxMG Cycle and I/O Parameters
      !
      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      INTEGER   EX_MPI_COMM, EX_MPI_IERROR, MPI_MyProc, NProc

! ----------------------------
!     Local Declarations
! ----------------------------

      INTEGER   CE, i, IO_in, MyProc, READ_STATUS
      LOGICAL   exist

! ==========================================================================

! -------------------------------------------
!     MPI indexing:
! -------------------------------------------

      MyProc = MPI_MyProc


! -------------------------------------------
!     Set default parameter values
! -------------------------------------------

      IF ( MyProc.EQ.0 ) THEN

         CALL BMG3_SymStd_SETUP_parms(  &
     &                    BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG  &
     &                    )

      ENDIF

! -------------------------------------------
!     Read parameter data
! -------------------------------------------

      IF ( MyProc.EQ.0 ) THEN



            !
            !  Open data file
            !

            !
            !  Read the parameters
            !
             BMG_iPARMS(id_BMG3_SETUP) = 0
             BMG_iPARMS(id_BMG3_STENCIL) = 1
             BMG_iPARMS(id_BMG3_SEPARABLE) = bmg_sep
!            write(0,*) 'BMG_iPARMS(id_BMG3_SEPARABLE) = ',
!     :        BMG_iPARMS(id_BMG3_SEPARABLE),id_BMG3_SEPARABLE
             BMG_iPARMS(id_BMG3_RELAX)  = 1
             BMG_iPARMS(id_BMG3_NRELAX_DOWN) = 1
             BMG_iPARMS(id_BMG3_NRELAX_UP)  = 1
             BMG_iPARMS(id_BMG3_RELAX_SYM )  = 0

             BMG_iPARMS(id_BMG3_CYCLE_CLASS) = bmg_cycletype
             BMG_iPARMS(id_BMG3_NCYCLE_TYPE) = 1  ! 1=v-cycles, 2=w-cycles

             BMG_iPARMS(id_BMG3_CG_SOLVER) = bmg_cgsolve
             BMG_iPARMS(id_BMG3_CG_MIN_DIM) = 3   ! min dimension for coarse grid

             BMG_iPARMS(id_BMG2_CYCLE_CLASS) = 1 ! for 2-d solver?
             BMG_iPARMS(id_BMG2_NCYCLE_TYPE) = 1 ! for 2-d solver?

             BMG_iPARMS(id_BMG2_CG_SOLVER) = 0
             BMG_iPARMS(id_BMG2_CG_MIN_DIM) = 3

             BMG_iPARMS(id_BMG3_MAX_ITERS) = bmg_maxiter
             BMG_rPARMS(id_BMG3_STOP_TOL) = bmg_tol


            !
            !  Override default stopping criteria
            !
!            BMG_iPARMS(id_BMG3_STOP_TEST) = BMG_STOP_REL_RES_L2

            !
            !  Override debbugging flags
            !
!            BMG_IOFLAG(iBMG3_BUG_RES_RELAX)    = .TRUE.
!            BMG_IOFLAG(iBMG3_BUG_RES_CG_SOLVE) = .TRUE.
!            BMG_IOFLAG(iBMG3_BUG_STENCIL_CG) = .TRUE.
!             BMG_IOFLAG(iBMG3_BUG_RESTRICT)  = .true.

!            BMG_IOFLAG(iBMG3_OUT_STOP_ERROR) = .true.
            !
            !  Set the Coarse Grid Comm Pattern
            !
            BMG_iPARMS(id_BMG3_CG_COMM) = BMG_CG_GATHER_SCATTER

            BMG_iPARMS(id_BMG3_SYNC_INITIAL_GUESS)  &
     &           = BMG_SYNC_INITIAL_GUESS

      ENDIF

! -------------------------------------------
!     Broadcast read status
! -------------------------------------------

!      CALL MPI_Bcast( READ_STATUS, 1, MPI_INTEGER,
!     &                0, EX_MPI_COMM, EX_MPI_IERROR )


! -------------------------------------------
!     Broadcast parameters arrays
! -------------------------------------------

      CALL MPI_Bcast( BMG_iPARMS, NBMG_iPARMS, MPI_INTEGER,  &
     &                0, EX_MPI_COMM, EX_MPI_IERROR )

      CALL MPI_Bcast( BMG_rPARMS, NBMG_rPARMS, FLOAT_PRECISION,  &
     &                0, EX_MPI_COMM, EX_MPI_IERROR )

      CALL MPI_Bcast( BMG_IOFLAG, NBMG_IOFLAG, MPI_LOGICAL,  &
     &                0, EX_MPI_COMM, EX_MPI_IERROR )

! -------------------------------------------
!     Override some I/O on the master
! -------------------------------------------

!         BMG_IOFLAG(iBMG3_OUT_ITERATIONS)   = .TRUE.

      IF ( MyProc.EQ.0 ) THEN

         !
         !  Override I/O Parameters
         !
         BMG_IOFLAG(iBMG3_OUT_WSPACE_POINT) = .FALSE.
         BMG_IOFLAG(iBMG3_OUT_WSPACE_SIZE)  = .FALSE.

         BMG_IOFLAG(iBMG3_OUT_ITERATIONS)   = .TRUE.
         BMG_IOFLAG(iBMG3_OUT_TIME_SETUP)   = .TRUE.
         BMG_IOFLAG(iBMG3_OUT_TIME_CYCLING) = .TRUE.
         BMG_IOFLAG(iBMG3_OUT_TIME_TOTAL)   = .TRUE.

         BMG_IOFLAG(iBMG3_BUG_PARAMETERS)   = .FALSE.
         !
         !  Debugging parameters?
         !
         IF ( BMG_IOFLAG(iBMG3_BUG_PARAMETERS) ) THEN
            CALL BMG3_SymStd_DUMP_parms(  &
     &                       BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG  &
     &                       )
         ENDIF

      ENDIF

! ==========================================================================

 500  FORMAT (/,'FATAL ERROR: EX_SETUP_BMG_parms.f',  &
     &        //,5X,A,/)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

! ===========================================


      RETURN
      END SUBROUTINE EX_SETUP_BMG_parms




!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE WRF_SymStd_ncycle(  &
     &                KC, KCF, NOG,  &
     &                IFD, IU, ID, IVW, IRELAX, IRELAX_SYM, ISEP,  &
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,  &
     &                Q, QF, RES, NF, NC,  &
     &                SO, NSO, SOR, NSOR, CI, NCI, IGRD, NOGm,  &
     &                BMG_iWORK_PL, NBMG_iWORK_PL,  &
     &                BMG_rWORK_PL, NBMG_rWORK_PL,  &
     &                BMG_iWORK_CS, NBMG_iWORK_CS,  &
     &                BMG_rWORK_CS, NBMG_rWORK_CS,  &
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,  &
     &                MSG_Buffer, NMSGr,  &
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID  &
     &                ,Qmask)

! ==========================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!   BMG3_SymStd_ncycle performs a single multigrid n-cycle.
!
! =======================================================================
! $license_flag$
! =======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   INPUT/OUTPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   OUTPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   LOCAL:
!  --------------------
!
!
! ==========================================================================

      IMPLICIT NONE

! -----------------------------
!     Includes
!
      INCLUDE 'mpif.h'

#include    "BMG_constants.h"
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_parameters.h'

! ----------------------------
!     Argument Declarations
!
      INTEGER  KC, KCF, NOG, IFD, IU, ID, IVW, IRELAX, IRELAX_SYM, ISEP,  &
     &         NF, NC, NSO, NSOR, NCI, NOGm, NMSGi, NMSGr,  &
     &         NBMG_MSG_iGRID, NBMG_iWORK_CS, NBMG_iWORK_PL,  &
     &         NBMG_rWORK_CS, NBMG_rWORK_PL

      INTEGER  BMG_iPARMS(NBMG_iPARMS),    &
     &         BMG_iWORK_CS(NBMG_iWORK_CS),  &
     &         BMG_iWORK_PL(NBMG_iWORK_PL),  &
     &         iWorkMSG(NMSGi), pMSG(NBMG_pMSG,NOGm),  &
     &         pMSGSO(NBMG_pMSG,NOGm), IGRD(NOGm,NBMG_pIGRD),  &
     &         BMG_MSG_iGRID(NBMG_MSG_iGRID),  &
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS),  &
     &         BMG_rWORK_CS(NBMG_rWORK_CS),  &
     &         BMG_rWORK_PL(NBMG_rWORK_PL),  &
     &         CI(NCI), Q(NF), QF(NF), RES(NF),  &
     &         SO(NSO), SOR(NSOR), MSG_Buffer(NMSGr)
      REAL*RKIND :: Qmask(NF)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

! ----------------------------
!     Local Declarations
!
      INTEGER  MyProc, MyProcI, MyProcJ, MyProcK, NProc, ntasks_x, ntasks_y,  &
     &         NProcK, MPICOMM, MPI_IERROR, NStncl

      INTEGER  iGs, jGs, kGs, NLx, NLy, NLz, NGx, NGy, NGz

      INTEGER  p_ProcGrid, p_ProcCoord

      INTEGER  i, j, k, p_ABD, p_BBD, p_CI, p_SO, p_SOR, p_U,  &
     &         NCBW, NCU
      REAL*8   RES_L2
      REAL*8  TT1, TT2

! ==========================================================================

! --------------------------------------------------------------------
! >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
! --------------------------------------------------------------------
! --------------------------------------------
!     MPI Parallel Data Decomposition:
! --------------------------------------------
      !
      ! MPI Communicators
      !
      MPICOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM)
      !
      ! Processor Grid:
      !
      NProc  = BMG_MSG_iGRID(id_BMG_MSG_NProc)
      ntasks_x = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      ntasks_y = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)
      NProcK = 1 !BMG_MSG_iGRID(id_BMG_MSG_NProcK)

      MyProc  = BMG_MSG_iGRID(id_BMG_MSG_MyProc)
      MyProcI = BMG_MSG_iGRID(id_BMG_MSG_MyProcI)
      MyProcJ = BMG_MSG_iGRID(id_BMG_MSG_MyProcJ)
      MyProcK = BMG_MSG_iGRID(id_BMG_MSG_MyProcK)

      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)
      p_ProcCoord = BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord)

! --------------------------------------------------------------------
! >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
! --------------------------------------------------------------------

      k = KCF   ! set the current grid index

      IF (KCF.EQ.KC) THEN
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN
            WRITE(*,*) '*** FATAL ERROR: BMG3_SymStd_ncycle.f ****'
            WRITE(*,*) 'The current finest grid = ', KCF
            WRITE(*,*) 'The coarest grid        = ', KC
            WRITE(*,*) 'ERROR; minimum number of grids is 2 !'
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,18)
         RETURN

      ENDIF

      ! set the n-cycle counter to zero on all grids
      DO i = KC, KCF
         IGRD(i,idL_BMG_IVW)=iZERO
      END DO


! -----------------------------------------------------------
!     Begin the n-cycle
!     (this is also the n-cycle recursive loop boundary)
! -----------------------------------------------------------

 140  CONTINUE

         IF ( K.EQ.KC .AND. IGRD(k,idL_BMG_IVW).EQ.iZERO ) THEN
            !
            ! Solve on the coarsest grid
            !
            CALL BMG3_SymStd_GET_pointers(  &
     &                KC, IGRD, NOGm,  &
     &                p_U, p_SO, p_SOR, p_CI,  &
     &                NLx, NLy, NLz, NGx, NGy, NGz  &
     &                )
            !
            iGs = IGRD(KC,idL_BMG_Icoord)
            jGs = IGRD(KC,idL_BMG_Jcoord)
            kGs = IGRD(KC,idL_BMG_Kcoord)
            !
            IF ( BMG_iPARMS(id_BMG3_CG_SOLVER).EQ.BMG_CG_SOLVE_LU ) THEN
               !
               p_ABD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CSO)
               p_BBD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CU)
               !
               NCBW  = BMG_iWORK_CS(id_BMG_iWORK_CS_NCBW)
               NCU   = BMG_iWORK_CS(id_BMG_iWORK_CS_NCU)
               !
               TT1 = MPI_Wtime()
               !
               CALL BMG3_SymStd_SOLVE_cg_LU(  &
     &                   q(p_U), qf(p_U), NLx, NLy, NLz,  &
     &                   BMG_rWORK_CS(p_ABD), BMG_rWORK_CS(p_BBD),  &
     &                   NCBW, NCU, NOGm,  &
     &                   ntasks_x, ntasks_y, NProcK, NProc, MyProc,  &
     &                   BMG_MSG_iGRID(p_ProcGrid),  &
     &                   BMG_MSG_iGRID(p_ProcCoord),  &
     &                   iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),  &
     &                   MSG_Buffer, NMSGr, MPICOMM, BMG_iPARMS  &
     &                   )
               !
               TT2 = MPI_Wtime()
               !
            ELSE
               !
               TT1 = MPI_Wtime()
               !
               CALL BMG3_SymStd_SOLVE_cg_boxmg(   &
     &                   Q(p_U), QF(p_U), NOG,  &
     &                   BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,  &
     &                   NLx, NLy, NLz, NGx, NGy, NGz,  &
     &                   iGs, jGs, kGs,  &
     &                   BMG_iWORK_CS, NBMG_iWORK_CS,  &
     &                   BMG_rWORK_CS, NBMG_rWORK_CS,  &
     &                   ntasks_x, ntasks_y, NProcK, NProc, MyProc,  &
     &                   BMG_MSG_iGRID(p_ProcGrid),  &
     &                   BMG_MSG_iGRID(p_ProcCoord),  &
     &                   iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),  &
     &                   MSG_Buffer, NMSGr, MPICOMM  &
     &                   )

               !
               TT2 = MPI_Wtime()
               !
            ENDIF

            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            BMG_rPARMS(id_BMG3_TIME_SOLVE_CG)   &
     &          = BMG_rPARMS(id_BMG3_TIME_SOLVE_CG) + TT2 - TT1

            !
            IGRD(k,idL_BMG_IVW) = IVW
            IF( BMG_IOFLAG(iBMG3_BUG_RES_CG_SOLVE)   &
     &          .AND. BMG_iPARMS(id_BMG3_CG_SOLVER).EQ.BMG_CG_SOLVE_LU  &
     &        ) THEN
               NStncl = 14   ! always the same on coarse grid if NOG > 1
               CALL BMG3_SymStd_residual(   &
     &                   KC, NOG, IFD,  &
     &                   q(p_U), qf(p_U), so(p_SO), RES(p_U),  &
     &                   NLx, NLy, NLz, NStncl,  &
     &                   iWorkMSG, NMSGi, pMSG,  &
     &                   MSG_Buffer, NMSGr, MPICOMM  &
     &                   )
               CALL BMG3_SymStd_UTILS_norm_l2(  &
     &                   RES(p_U), NLx, NLy, NLz, RES_L2, MPICOMM  &
     &                   )
               !
               IF ( MyProc.EQ.1 ) THEN
                  WRITE (*,100) k, RES_L2
               ENDIF
               !
            ENDIF
            !
            ! continue the n-cycle
            !
            GOTO 140
            !
         ELSE IF ( IGRD(k,idL_BMG_IVW).LT.IVW ) THEN
            !
            ! Move down
            !
            CALL WRF_SymStd_updown(   &
     &                k, KCF, NOG, BMG_DOWN,   &
     &                IFD, IU, ID, IRELAX, IRELAX_SYM, ISEP,  &
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,  &
     &                Q, QF, RES, NF, NC,  &
     &                SO, NSO, SOR, NSOR, CI, NCI,  &
     &                IGRD, NOGm,  &
     &                BMG_iWORK_PL, NBMG_iWORK_PL,  &
     &                BMG_rWORK_PL, NBMG_rWORK_PL,  &
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,  &
     &                MSG_Buffer, NMSGr,  &
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID  &
     &                ,Qmask)
            !
            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            IGRD(k,idL_BMG_IVW) = IGRD(k,idL_BMG_IVW) + 1
            k=k-1
            !
            ! Continue the n-cycle
            !
            GOTO 140
            !
         ELSE IF ( IGRD(k,idL_BMG_IVW).EQ.IVW ) THEN
            !
            ! Move up
            !
            IGRD(k,idL_BMG_IVW) = iZERO   ! reset n-cycle counter to zero
            k=k+1
            !
            CALL WRF_SymStd_updown(   &
     &                k, KCF, NOG, BMG_UP,  &
     &                IFD, IU, ID, IRELAX, IRELAX_SYM, ISEP,  &
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,  &
     &                Q, QF, RES, NF, NC,  &
     &                SO, NSO, SOR, NSOR, CI, NCI,  &
     &                IGRD, NOGm,  &
     &                BMG_iWORK_PL, NBMG_iWORK_PL,  &
     &                BMG_rWORK_PL, NBMG_rWORK_PL,  &
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,  &
     &                MSG_Buffer, NMSGr,  &
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID  &
     &                ,Qmask)
            !
            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            !
            ! Are we on the finest grid?
            !
            IF (K.EQ.KCF) THEN
               ! n-cycle is done
               RETURN
            ELSE
               ! continue the n-cycle
               GOTO 140
            ENDIF
            !
         ENDIF

! ==========================================================================

 100     FORMAT (1X,'(3D)LEVEL ',I2,' RESIDUAL NORM = ',1P,E12.5)

! ===========================================

         RETURN
         END SUBROUTINE WRF_SymStd_ncycle

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE WRF_SymStd_updown (  &
     &                K, KF, NOG, UPDOWN, &
     &                IFD, IU, ID, IRELAX, IRELAX_SYM, ISEP, &
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, &
     &                Q, QF, RES, NF, NC,  &
     &                SO, NSO, SOR, NSOR, CI, NCI, &
     &                IGRD, NOGm,  &
     &                BMG_iWORK_PL, NBMG_iWORK_PL,  &
     &                BMG_rWORK_PL, NBMG_rWORK_PL,  &
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,  &
     &                MSG_Buffer, NMSGr,  &
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID &
     &                ,Qmask)


! ==========================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG3_SymStd_updown performs all the necessary tasks associated
!     with moving "down" to a coarser grid or moving "up" to a finer
!     grid.  It is necessary because f77 does not support recursion.
!
! =======================================================================
! $license_flag$
! =======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   INPUT/OUTPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   OUTPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   LOCAL:
!  --------------------
!
!
! ==========================================================================

      IMPLICIT NONE

! -----------------------------
!     Includes
! -----------------------------

      INCLUDE 'mpif.h'

#include    "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_parameters.h'
      INCLUDE 'BMG_workspace.h'

! ----------------------------
!     Argument Declarations
! ----------------------------

      !
      !  BoxMG parameters
      !
      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      !
      !  BoxMG Cycle Parameters
      !
      INTEGER  k, KC, KCF, KF, IFD, IU, ID, IVW,   &
     &         IRELAX, IRELAX_SYM, UPDOWN, ISEP

      !
      !  BoxMG dimensions
      !
      INTEGER  NF, NC, NSO, NSOR, NCI, NOGm, NOG,   &
     &         NBMG_iWORK_CS, NBMG_iWORK_PL,  &
     &         NBMG_rWORK_CS, NBMG_rWORK_PL

      !
      !  BoxMG internal pointers
      !
      INTEGER  IGRD(NOGm,NBMG_pIGRD)

      !
      !  BoxMG data arrays
      !
      REAL*RKIND   CI(NCI), Q(NF), QF(NF), RES(NF), SO(NSO), SOR(NSOR)
      REAL*RKIND :: Qmask(NF)

      !
      !  BoxMG plane relaxation workspace
      !
      INTEGER BMG_iWORK_PL(NBMG_iWORK_PL)
      REAL*RKIND  BMG_rWORK_PL(NBMG_rWORK_PL)

      !
      !  MSG workspace dimensions
      !
      INTEGER  NMSGi, NMSGr
      !
      !  MSG workspace and pointers
      !
      INTEGER  iWorkMSG(NMSGi),   &
     &         pMSG(NBMG_pMSG,NOG),  &
     &         pMSGSO(NBMG_pMSG,NOG)
      REAL*RKIND   MSG_Buffer(NMSGr)

      !
      !  BMG_MSG grid stuff
      !
      INTEGER  NBMG_MSG_iGRID
      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRID),   &
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)

! ----------------------------
!     Local Declarations
! ----------------------------

      INTEGER  I, J, NLx, NLy, NLz, NLx_c, NLy_c, NLz_c,  &
     &         NGx, NGy, NGz, NGx_c, NGy_c, NGz_c,  &
     &         p_CI, p_CIC, p_SO, p_SOC, p_SOR, p_SORC, p_U, p_UC
      INTEGER  MPICOMM, XYPLANECOMM, XZPLANECOMM, YZPLANECOMM,  &
     &         iGs, jGs, kGs, NSORv, NStncl
      INTEGER  MyProc, MyProcI, MyProcJ, MyProcK, MPI_IERROR,  &
     &         NProc, ntasks_x, ntasks_y, NProcK, p_ProcGrid, p_ProcCoord
      REAL*8   RES_L2
      REAL*8   TT1, TT2

! ==========================================================================

! --------------------------------------------------------------------
! >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
! --------------------------------------------------------------------
! --------------------------------------------
!     MPI Parallel Data Decomposition:
! --------------------------------------------
      !
      ! MPI Communicators
      !
      MPICOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM)
      !
      ! MPI Communicatorsn for plane relaxation
      !
      XYPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_xy)
      XZPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_xz)
      YZPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_yz)

      !
      ! Processor Grid:
      !
      NProc  = BMG_MSG_iGRID(id_BMG_MSG_NProc)
      ntasks_x = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      ntasks_y = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)
      NProcK = 1! BMG_MSG_iGRID(id_BMG_MSG_NProcK)

      MyProc  = BMG_MSG_iGRID(id_BMG_MSG_MyProc)
      MyProcI = BMG_MSG_iGRID(id_BMG_MSG_MyProcI)
      MyProcJ = BMG_MSG_iGRID(id_BMG_MSG_MyProcJ)
      MyProcK = BMG_MSG_iGRID(id_BMG_MSG_MyProcK)

      !
      ! Processor Coordinates
      !
      iGs = IGRD(k,idL_BMG_Icoord)
      jGs = IGRD(k,idL_BMG_Jcoord)
      kGs = IGRD(k,idL_BMG_Kcoord)

      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)
      p_ProcCoord = BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord)

! --------------------------------------------------------------------
! >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
! --------------------------------------------------------------------

      !
      ! Set NStncl for the finest grid
      !
      IF ( K.LT.NOG .OR. IFD.NE.BMG_STENCIL_7pt ) THEN
         NStncl=14
      ELSE
         NStncl=4
      ENDIF

      !
      !  Number of temporary vectors in SOR
      !
      IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
         NSORv = 2
      ELSE
         NSORv = 2  ! the same for now but this may change
      ENDIF

      IF (UPDOWN.EQ.BMG_DOWN) THEN

         ! Get pointers for grids k and k-1
         CALL BMG3_SymStd_GET_pointers(    &
     &               K, IGRD, NOGm,  &
     &               p_U, p_SO, p_SOR, p_CI,  &
     &               NLx, NLy, NLz, NGx, NGy, NGz   &
     &               )
         CALL BMG3_SymStd_GET_pointers(   &
     &               K-1, IGRD, NOGm,  &
     &               p_UC, p_SOC, p_SORC, p_CIC,  &
     &               NLx_c, NLy_c, NLz_c, NGx_c, NGy_c, NGz_c  &
     &               )

         ! Relaxation
         DO 10 i=1, ID
            !
            TT1 = MPI_Wtime()

            CALL WRF_SymStd_relax(   &
     &                k, KF, NOG, UPDOWN,   &
     &                IFD, NStncl, IRELAX, iRELAX_SYM, RES_L2, ISEP,  &
     &                Q, QF, RES, NF,  &
     &                SO, NSO, SOR, NSOR, CI, NCI,  &
     &                IGRD, NOGm,  &
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,  &
     &                BMG_iWORK_PL, NBMG_iWORK_PL,  &
     &                BMG_rWORK_PL, NBMG_rWORK_PL,  &
     &                iGs, jGs, kGs,  &
     &                iWorkMSG, NMSGi, pMSG, MSG_Buffer, NMSGr,  &
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,  &
     &                MPICOMM, XYPLANECOMM, XZPLANECOMM, YZPLANECOMM  &
     &                ,Qmask)

            TT2 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_relax) =   &
     &           BMG_rPARMS(id_BMG3_TIME_relax) + TT2 - TT1


            !
            IF( BMG_IOFLAG(iBMG3_BUG_RES_RELAX) .AND. MyProc.EQ.1 ) THEN
               WRITE(*,100) K, RES_L2
            ENDIF
            !
 10      CONTINUE

         ! restrict the residual
         TT1 = MPI_Wtime()

         ! calculate the residual
         CALL BMG3_SymStd_residual(   &
     &             K, NOG, IFD,  &
     &             Q(p_U), QF(p_U), SO(p_SO), RES(p_U),  &
     &             NLx, NLy, NLz, NStncl,  &
     &             iWorkMSG, NMSGi, pMSG,  &
     &             MSG_Buffer, NMSGr, MPICOMM  &
     &             )

         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_PCG_PRECON) =  &
     &        BMG_rPARMS(id_BMG3_TIME_PCG_PRECON) + TT2 - TT1

!$$$         CALL BMG3_SymStd_DUMP_vector(  &
!$$$     &             BMG_IOFLAG, RES(p_U), K, NOG,
!$$$     &             NLx, NLy, NLz, NGx, NGy, NGz,
!$$$     &             IGRD(K,idL_BMG_ICoord),
!$$$     &             IGRD(K,idL_BMG_JCoord),
!$$$     &             IGRD(K,idL_BMG_KCoord),
!$$$     &             'DOWN-RES-post-relax', .FALSE.,
!$$$     &             BMG_MSG_iGRID(p_ProcGrid),
!$$$     &             ntasks_x, ntasks_y, NProcK, NProc,
!$$$     &             MyProcI, MyProcJ, MyProcK, MPICOMM  &
!$$$     &             )

!         ! restrict the residual
!         TT1 = MPI_Wtime()

         CALL BMG3_SymStd_restrict(  &
     &             K, K-1,  &
     &             RES(p_U), QF(p_UC), CI(p_CIC),  &
     &             NLx, NLy, NLz, NLx_c, NLy_c, NLz_c,  &
     &             iGs, jGs, kGs   &
     &             )

         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_restrict) =   &
     &        BMG_rPARMS(id_BMG3_TIME_restrict) + TT2 - TT1


!$$$         CALL BMG3_SymStd_DUMP_vector(  &
!$$$     &             BMG_IOFLAG, QF(p_UC), K-1, NOG,
!$$$     &             NLx_c, NLy_c, NLz_c, NGx_c, NGy, NGz,
!$$$     &             IGRD(K-1,idL_BMG_ICoord),
!$$$     &             IGRD(K-1,idL_BMG_JCoord),
!$$$     &             IGRD(K-1,idL_BMG_KCoord),
!$$$     &             'DOWN-QF-post-restrict', .FALSE.,
!$$$     &             BMG_MSG_iGRID(p_ProcGrid),
!$$$     &             ntasks_x, ntasks_y, NProcK, NProc,
!$$$     &             MyProcI, MyProcJ, MyProcK, MPICOMM  &
!$$$     &             )

         ! Zero the initial guess
         CALL BMG3_SymStd_UTILS_rV_zero(   &
     &               k-1, q(p_UC), NLx_c, NLy_c, NLz_c  &
     &               )


      ELSE IF (UPDOWN.EQ.BMG_UP) THEN

         ! Get pointers for grids k and k-1
         CALL BMG3_SymStd_GET_pointers(  &
     &             k, IGRD, NOGm,  &
     &             p_U, p_SO, p_SOR, p_CI,  &
     &             NLx, NLy, NLz, NGx, NGy, NGz  &
     &             )
         CALL BMG3_SymStd_GET_pointers(  &
     &             k-1, IGRD, NOGm,  &
     &             p_UC, p_SOC, p_SORC, p_CIC,  &
     &             NLx_c, NLy_c, NLz_c, NGx_c, NGy_c, NGz_c  &
     &             )


         ! Interpolate and Correct
         TT1 = MPI_Wtime()

!        write(0,*) 'SymStd_updown: k, NOG = ',k,NOG
            IF ( k < NOG ) THEN
            CALL BMG3_SymStd_interp_add(  &
     &             k-1, k, NOG,  &
     &             Q(p_U) ,Q(p_UC), RES(p_U),  &
     &             SO(p_SO), NStncl, CI(p_CIC),  &
     &             NLx_c, NLy_c, NLz_c,  &
     &             NLx, NLy, NLz,  &
     &             iGs, jGs, kGs,  &
     &             iWorkMSG, NMSGi, pMSG,  &
     &             MSG_Buffer, NMSGr, MPICOMM  &
     &             )

            ELSE ! special version for fine grid (checks mask array)
!        write(0,*) 'SymStd_updown: call WRF_SymStd_interp_add = ',k,NOG
         CALL WRF_SymStd_interp_add(  &
     &             k-1, k, NOG,  &
     &             Q(p_U) ,Q(p_UC), RES(p_U),  &
     &             SO(p_SO), NStncl, CI(p_CIC),  &
     &             NLx_c, NLy_c, NLz_c,  &
     &             NLx, NLy, NLz,  &
     &             iGs, jGs, kGs,  &
     &             iWorkMSG, NMSGi, pMSG, Qmask,   &
     &             MSG_Buffer, NMSGr, MPICOMM   &
     &             )
         ENDIF

         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_interp_add) =   &
     &        BMG_rPARMS(id_BMG3_TIME_interp_add) + TT2 - TT1



         IF ( BMG_IOFLAG(iBMG3_BUG_RES_INTERP) ) THEN
            CALL BMG3_SymStd_residual(  &
     &                K, NOG, IFD,  &
     &                Q(p_U), QF(p_U), SO(p_SO), RES(p_U),  &
     &                NLx, NLy, NLz, NStncl,  &
     &                iWorkMSG, NMSGi, pMSG,  &
     &                MSG_Buffer, NMSGr, MPICOMM   &
     &                )
            CALL BMG3_SymStd_UTILS_norm_l2(   &
     &                RES(p_U), NLx, NLy, NLz, RES_L2, MPICOMM  &
     &                )
            WRITE(*,110) K, RES_L2
         ENDIF

         ! Relaxation
         DO 20 i=1, IU
            !
            TT1 = MPI_Wtime()

            CALL WRF_SymStd_relax(  &
     &                k, KF, NOG, UPDOWN,  &
     &                IFD, NStncl, IRELAX, iRELAX_SYM, RES_L2, ISEP,  &
     &                Q, QF, RES, NF,  &
     &                SO, NSO, SOR, NSOR, CI, NCI,  &
     &                IGRD, NOGm,  &
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,  &
     &                BMG_iWORK_PL, NBMG_iWORK_PL,  &
     &                BMG_rWORK_PL, NBMG_rWORK_PL,  &
     &                iGs, jGs, kGs,  &
     &                iWorkMSG, NMSGi, pMSG, MSG_Buffer, NMSGr,  &
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,  &
     &                MPICOMM, XYPLANECOMM, XZPLANECOMM, YZPLANECOMM  &
     &                ,Qmask)
            TT2 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_relax) =   &
     &           BMG_rPARMS(id_BMG3_TIME_relax) + TT2 - TT1

            !
            IF( BMG_IOFLAG(iBMG3_BUG_RES_RELAX) .AND. MyProc.EQ.1 ) THEN
               WRITE(*,100) K, RES_L2
            ENDIF
            !
 20      CONTINUE

      ELSE
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN
            WRITE(*,*) '*** UPDOWN out of range: '
            WRITE(*,*) 'HAVE: UPDOWN = ', UPDOWN
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,19)
         RETURN

      ENDIF

! ==========================================================================

 100  FORMAT (1X,'(3D)LEVEL ',I2,' RESIDUAL NORM = ',1P,E12.5)
 110  FORMAT (1X,'(3D)LEVEL ',I2,  &
     &           ' AFTER INTERPOLATION RES_L2 = ',1P,E12.5 )

! ===========================================

      RETURN
      END SUBROUTINE WRF_SymStd_updown

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE WRF_SymStd_relax(   &
     &                KG, KF, NOG, UPDOWN,   &
     &                IFD, NStncl, IRELAX, iRELAX_SYM, RES_L2, ISEP,  &
     &                Q, QF, RES, NF,  &
     &                SO, NSO, SOR, NSOR, CI, NCI,   &
     &                IGRD, NOGm,  &
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,  &
     &                BMG_iWORK_PL, NBMG_iWORK_PL,  &
     &                BMG_rWORK_PL, NBMG_rWORK_PL,  &
     &                iGs, jGs, kGs,  &
     &                iWorkMSG, NMSGi, pMSG, MSG_Buffer, NMSGr,  &
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,  &
     &                MPICOMM, XYPLANECOMM, XZPLANECOMM, YZPLANECOMM  &
     &                ,Qmask)

! ==========================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG3_SymStd_relax performs one relaxation step (either colored
!     point Gauss Seidel, or alternating plane red-black Gauss-Seidel)
!
! =======================================================================
! $license_flag$
! =======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   INPUT/OUTPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   OUTPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   LOCAL:
!  --------------------
!
!
! ==========================================================================

      IMPLICIT NONE

! ------------------------------------------------
!     Includes
!
#include    "BMG_constants.h"
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

! ------------------------------------------------
!     Argument Declarations
!
      INTEGER  NCI, NSO, NSOR, NSORv, NStncl, kg, kf, NBMG_iWORK_PL,   &
     &         NBMG_rWORK_PL, NOG, NF, NOGm, NMSGi, NMSGr

      INTEGER  NBMG_MSG_iGRID
      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRID),   &
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)

      INTEGER  BMG_iPARMS(NBMG_iPARMS), BMG_iWORK_PL(NBMG_iWORK_PL),   &
     &         IFD, IGRD(NOGm,24), IRELAX_SYM, UPDOWN, IRELAX, ISEP
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS), BMG_rWORK_PL(NBMG_rWORK_PL),   &
     &         CI(NCI), q(NF), qf(NF), RES(NF),   &
     &         SO(NSO), SOR(NSOR), MSG_Buffer(NMSGr)

      REAL*8   RES_L2

      REAL*RKIND :: Qmask(NF)

      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      INTEGER  iGs, jGs, kGs,   &
     &         MPICOMM, iWorkMSG(NMSGi), pMSG(NBMG_pMSG,NOGm)
      INTEGER  XYPLANECOMM, XZPLANECOMM, YZPLANECOMM
      INTEGER  MyProc, MyProcI, MyProcJ, MyProcK, MPI_IERROR,  &
     &         NProc, ntasks_x, ntasks_y, NProcK


! ----------------------------
!     Local Declarations
!
      INTEGER  NLx, NLy, NLz, NGx, NGy, NGz,  &
     &         NLx_f, NLy_f, NLz_f, NGx_f, NGy_f, NGz_f,  &
     &         p_CI, p_SO, p_SOR, p_U, p_ProcGrid

! ==========================================================================

! --------------------------------------------------------------------
! >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
! --------------------------------------------------------------------

      !
      ! MPI Communicators
      !
      MPICOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM)
      !
      ! MPI Communicatorsn for plane relaxation
      !
      XYPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_xy)
      XZPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_xz)
      YZPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_yz)

      !
      ! Processor Grid:
      !
      NProc  = BMG_MSG_iGRID(id_BMG_MSG_NProc)
      ntasks_x = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      ntasks_y = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)
      NProcK = 1 !BMG_MSG_iGRID(id_BMG_MSG_NProcK)

      MyProc  = BMG_MSG_iGRID(id_BMG_MSG_MyProc)
      MyProcI = BMG_MSG_iGRID(id_BMG_MSG_MyProcI)
      MyProcJ = BMG_MSG_iGRID(id_BMG_MSG_MyProcJ)
      MyProcK = BMG_MSG_iGRID(id_BMG_MSG_MyProcK)

      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)

! --------------------------------------------------------------------
! >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
! --------------------------------------------------------------------

      !
      !  Number of temporary vectors in SOR
      !
      IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
         NSORv = 2
      ELSE
         NSORv = 2  ! the same for now but this may change
      ENDIF

      !
      ! Collect fine-grid dimensions
      ! ( need to check: KF or NOG )
      !
      CALL BMG3_SymStd_GET_pointers(   &
     &            KF, IGRD, NOGm,  &
     &            p_U, p_SO, p_SOR, p_CI,  &
     &            NLx_f, NLy_f, NLz_f, NGx_f, NGy_f, NGz_f  &
     &            )

      CALL BMG3_SymStd_GET_pointers(  &
     &            KG, IGRD, NOGm,  &
     &            p_U, p_SO, p_SOR, p_CI,   &
     &            NLx, NLy, NLz, NGx, NGy, NGz  &
     &            )


      IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
         !
         ! Gauss-Seidel relaxation
         !
         CALL WRF_SymStd_relax_GS(  &
     &             KG, SO(p_SO), QF(p_U), Q(p_U), RES(p_U), SOR(p_SOR),  &
     &             NLx, NLy, NLz, NGx, NGy, NGz,  &
     &             RES_L2, BMG_IOFLAG, NOG, NOGm,  &
     &             IFD, NStncl, NSORv, IRELAX, IRELAX_SYM, UPDOWN,  &
     &             iGs, jGs, kGs,  &
     &             MyProcI, MyProcJ, MyProcK, MyProc,  &
     &             BMG_MSG_iGRID(p_ProcGrid),  &
     &             ntasks_x, ntasks_y, NProcK, NProc,  &
     &             iWorkMSG, NMSGi, pMSG, MSG_Buffer, NMSGr, MPICOMM  &
     &            ,ISEP,Qmask)
         !
      ELSEIF ( IRELAX.EQ.BMG_GS_RB_planes_xy_yz_xz ) THEN

         WRITE(*,*) 'Fatal Error: To busy to code planes !!!'
         CALL MPI_FINALIZE(MPI_IERROR)
         STOP

         !
      ENDIF

! ==========================================================================

      RETURN
      END SUBROUTINE WRF_SymStd_relax


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      SUBROUTINE WRF_SymStd_relax_GS(  &
     &                KG, SO, QF, Q, RES, SOR,  &
     &                NLx, NLy, NLz, NGx, NGy, NGz,  &
     &                RES_L2, BMG_IOFLAG, NOG, NOGm,  &
     &                IFD, NStncl, NSORv, IRELAX, IRELAX_SYM, UPDOWN,  &
     &                iGs, jGs, kGs,  &
     &                MyProcI, MyProcJ, MyProcK, MyProc,  &
     &                ProcGrid, ntasks_x, ntasks_y, NProcK, NProc,  &
     &                iWork, NMSGi, pMSG, MSG_Buffer, NMSGr, MPICOMM  &
     &               ,ISEP,Qmask)

! ==========================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG3_SymStd_relax_GS performs one sweep of Gauss Seidel (with the
!     correct ordering depending on whether we are on the way down, or
!     up in the symmetric cycling case)
!
! =======================================================================
! $license_flag$
! =======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   INPUT/OUTPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   OUTPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   LOCAL:
!  --------------------
!
!
! ==========================================================================

      IMPLICIT NONE

! ------------------------------------------------
!     Includes
!
      INCLUDE 'mpif.h'
      INCLUDE 'MSG_f90.h'

#include    "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_parameters.h'
      INCLUDE 'BMG_workspace.h'

! ------------------------------------------------
!     Argument Declarations
!
      INTEGER NGx, NGy, NGz, NLx, NLy, NLz

      INTEGER IFD, IRELAX, IRELAX_SYM, KG, NOG, NOGm,  &
     &        NSORv, NStncl, UPDOWN, NMSGi, NMSGr

      INTEGER iGs, jGs, kGs
      INTEGER MyProc, MyProcI, MyProcJ, MyProcK,  &
     &        NProc, ntasks_x, ntasks_y, NProcK
      INTEGER ProcGrid(ntasks_x,ntasks_y,NProcK)

      INTEGER MPICOMM, iWork(NMSGi), pMSG(NBMG_pMSG,NOGm)
      REAL*RKIND  Q(NLx,NLy,NLz), QF(NLx,NLy,NLz),  &
     &        RES(NLx,NLy,NLz),  &
     &        SO(NLx+1,NLy+1,NLz+1,NStncl), SOR(NLx,NLy,NLz,NSORv)

      REAL*8 RES_L2

      REAL*RKIND  MSG_Buffer(NMSGr)

      REAL*RKIND :: Qmask(NLx,Nly,Nlz)

      LOGICAL BMG_IOFLAG(NBMG_IOFLAG)

! ------------------------------------------------
!     Local Declarations
!
      INTEGER i, i1, j ,j1, k, k1

      INTEGER ie, iw, js, jn, kbot, kt

      INTEGER pts, ibeg, iend, ptsred, ptsblk
      INTEGER ptstart, ptend, ptstride
      INTEGER ptrn, ierror
      INTEGER relaxtest
      LOGICAL comhalf

      REAL*RKIND SOx(NLx+1),SOy(NLy+1),SOz(NLz+1)

      INTEGER ISEP

! ==========================================================================

      j1 = NLy-1
      i1 = NLx-1
      k1 = NLz-1

!      write(0,*) 'mproc, mpi,mpj,mpk, igs,jgs,kgs = ',myproc,myproci,
!     :   myprocj,myprock,igs,jgs,kgs
!      write(0,*) 'mproc, isep, BMG_SEPARABLE = ',myproc,isep,
!     :   BMG_SEPARABLE, BMG_NONSEPARABLE


!      write(0,*) 'RELAX_GS: kg,nog,ifd = ',KG,NOG,ifd

      comhalf = .false.


      IF ( KG.LT.NOG .OR. ifd.NE.1 ) THEN

         !
         !   27-point relaxation (8-color)
         !

!         IF ((UPDOWN.eq.BMG_UP).or.(IRELAX_SYM.EQ.BMG_RELAX_NONSYM))
!     &              THEN
!            ptstart = 1
!            ptend   = 8
!            ptstride = 1
!         ELSE
!            ptstart = 8
!            ptend = 1
!            ptstride = -1
!         ENDIF
!
!         DO pts = ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<
!
!            DO k=2+mod(mod((pts-1)/4,2)+mod(kGs+1,2),2),K1,2
!               !
!               DO j=2+mod(mod((pts-1)/2,2)+mod(jGs+1,2),2),J1,2
!                  !
!                  DO i=2+mod(mod(pts-1,2)+mod(iGs+1,2),2),I1,2
                     !

         IF ( (UPDOWN.eq.BMG_UP)   &
     &        .OR. (IRELAX_SYM.EQ.BMG_RELAX_NONSYM)) THEN
            ptstart  = 0
            ptend    = 1
            ptstride = 1
         ELSE
            ptstart  = 1
            ptend    = 0
            ptstride = -1
         ENDIF


        relaxtest = 2

        IF ( relaxtest .eq. 0 ) THEN ! default

         DO pts = ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<

!            DO k=2+mod(mod((pts-1)/4,2)+mod(kGs+1,2),2),K1,2
!               !
!               DO j=2+mod(mod((pts-1)/2,2)+mod(jGs+1,2),2),J1,2
!                  !
!                  DO i=2+mod(mod(pts-1,2)+mod(iGs+1,2),2),I1,2
             DO k = 2,k1

               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                     !
                     Q(i,j,k) = ( QF(i,j,k)  &
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)  &
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)  &
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)  &
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)  &
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)  &
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)  &
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)  &
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)  &
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)  &
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)  &
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)  &
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)  &
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)  &
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)  &
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)  &
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)  &
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)  &
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)  &
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)  &
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)  &
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)  &
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)  &
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)  &
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)  &
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)  &
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)  &
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO
            !

            IF ( comhalf ) THEN

            ptrn = 1

            CALL MSG_tbdx_send(Q, MSG_Buffer,  &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Proc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),  &
     &           iWork(pMSG(ipL_MSG_Index,KG)),  &
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,  &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Proc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),  &
     &           iWork(pMSG(ipL_MSG_Index,KG)),  &
     &           ptrn, ierror)
            !

            ENDIF
          ENDDO

        ELSEIF ( relaxtest .eq. 1 ) THEN ! non-redblack GS

             DO k = 2,k1

               DO j=2, j1
                  !
                  !
                  DO i = 2,i1
                     !
                     Q(i,j,k) = ( QF(i,j,k)  &
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)  &
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)  &
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)  &
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)  &
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)  &
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)  &
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)  &
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)  &
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)  &
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)  &
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)  &
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)  &
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)  &
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)  &
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)  &
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)  &
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)  &
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)  &
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)  &
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)  &
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)  &
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)  &
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)  &
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)  &
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)  &
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)  &
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO

        ELSEIF ( relaxtest .eq. 2 ) THEN ! layered redblack GS

!         DO pts=ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<

! do 'red' points on first plane
            ptsred = ptstart
            ptsblk = ptend
            DO  k=2, 2
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+ptsred,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                     Q(i,j,k) = ( QF(i,j,k)  &
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)  &
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)  &
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)  &
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)  &
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)  &
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)  &
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)  &
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)  &
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)  &
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)  &
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)  &
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)  &
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)  &
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)  &
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)  &
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)  &
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)  &
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)  &
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)  &
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)  &
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)  &
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)  &
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)  &
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)  &
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)  &
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)  &
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO

! do 'red' points on k plane and 'black' points on k-1 plane

            DO  k=3, k1

               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+ptsred,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                     Q(i,j,k) = ( QF(i,j,k)  &
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)  &
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)  &
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)  &
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)  &
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)  &
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)  &
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)  &
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)  &
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)  &
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)  &
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)  &
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)  &
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)  &
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)  &
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)  &
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)  &
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)  &
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)  &
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)  &
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)  &
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)  &
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)  &
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)  &
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)  &
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)  &
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)  &
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO

            IF ( comhalf ) THEN

            ptrn = 1

            CALL MSG_tbdx_send(Q, MSG_Buffer,  &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Proc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),  &
     &           iWork(pMSG(ipL_MSG_Index,KG)),  &
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,  &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Proc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),  &
     &           iWork(pMSG(ipL_MSG_Index,KG)),  &
     &           ptrn, ierror)
            !

            ENDIF

! black points
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-2)+ptsblk,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                     Q(i,j,k-1) = ( QF(i,j,k-1)  &
     &                    + SO(i,j,k-1,kpw)*Q(i-1,j,k-1)  &
     &                    + SO(i,j+1,k-1,kpnw)*Q(i-1,j+1,k-1)  &
     &                    + SO(i,j+1,k-1,kps)*Q(i,j+1,k-1)  &
     &                    + SO(i+1,j+1,k-1,kpsw)*Q(i+1,j+1,k-1)  &
     &                    + SO(i+1,j,k-1,kpw)*Q(i+1,j,k-1)  &
     &                    + SO(i+1,j,k-1,kpnw)*Q(i+1,j-1,k-1)  &
     &                    + SO(i,j,k-1,kps)*Q(i,j-1,k-1)  &
     &                    + SO(i,j,k-1,kpsw)*Q(i-1,j-1,k-1)  &
     &                    + SO(i,j,k-1,kb)*Q(i,j,k-1-1)  &
     &                    + SO(i,j,k-1,kbw)*Q(i-1,j,k-1-1)  &
     &                    + SO(i,j+1,k-1,kbnw)*Q(i-1,j+1,k-1-1)  &
     &                    + SO(i,j+1,k-1,kbn)*Q(i,j+1,k-1-1)  &
     &                    + SO(i+1,j+1,k-1,kbne)*Q(i+1,j+1,k-1-1)  &
     &                    + SO(i+1,j,k-1,kbe)*Q(i+1,j,k-1-1)  &
     &                    + SO(i+1,j,k-1,kbse)*Q(i+1,j-1,k-1-1)  &
     &                    + SO(i,j,k-1,kbs)*Q(i,j-1,k-1-1)  &
     &                    + SO(i,j,k-1,kbsw)*Q(i-1,j-1,k-1-1)  &
     &                    + SO(i,j,k-1+1,kb)*Q(i,j,k-1+1)  &
     &                    + SO(i,j,k-1+1,kbe)*Q(i-1,j,k-1+1)  &
     &                    + SO(i,j+1,k-1+1,kbse)*Q(i-1,j+1,k-1+1)  &
     &                    + SO(i,j+1,k-1+1,kbs)*Q(i,j+1,k-1+1)  &
     &                    + SO(i+1,j+1,k-1+1,kbsw)*Q(i+1,j+1,k-1+1)  &
     &                    + SO(i+1,j,k-1+1,kbw)*Q(i+1,j,k-1+1)  &
     &                    + SO(i+1,j,k-1+1,kbnw)*Q(i+1,j-1,k-1+1)  &
     &                    + SO(i,j,k-1+1,kbn)*Q(i,j-1,k-1+1)  &
     &                    + SO(i,j,k-1+1,kbne)*Q(i-1,j-1,k-1+1)  &
     &                    )*SOR(i,j,k-1,msor)
                     !
                  ENDDO
                  !
               ENDDO

               !
            ENDDO

! do 'black' points on last plane
            DO  k=k1, k1
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+ptsblk,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                     Q(i,j,k) = ( QF(i,j,k)  &
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)  &
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)  &
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)  &
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)  &
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)  &
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)  &
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)  &
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)  &
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)  &
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)  &
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)  &
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)  &
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)  &
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)  &
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)  &
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)  &
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)  &
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)  &
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)  &
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)  &
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)  &
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)  &
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)  &
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)  &
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)  &
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)  &
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO

       ENDIF ! relaxtest
       !
!         ENDDO   ! >>> END: loop over colors <<<<<<<<<<<<<<<<<<<<<<<<<<

        IF ( relaxtest .ne. 0 ) THEN

            ptrn = 1
            !
            CALL MSG_tbdx_send(Q, MSG_Buffer,  &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)), &
     &           iWork(pMSG(ipL_MSG_Proc,KG)),   &
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),  &
     &           iWork(pMSG(ipL_MSG_Index,KG)),  &
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,  &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Proc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),  &
     &           iWork(pMSG(ipL_MSG_Index,KG)),  &
     &           ptrn, ierror)

        ENDIF
            !
!         ENDDO   ! >>> END: loop over colors <<<<<<<<<<<<<<<<<<<<<<<<<<


         !
      ELSE
         !
         !  7 point relaxation (2-color)
         !
         IF ( (UPDOWN.eq.BMG_UP) &
     &        .OR. (IRELAX_SYM.EQ.BMG_RELAX_NONSYM)) THEN
            ptstart  = 0
            ptend    = 1
            ptstride = 1
         ELSE
            ptstart  = 1
            ptend    = 0
            ptstride = -1
         ENDIF

      IF ( ISEP .eq. BMG_SEPARABLE ) THEN


      DO  k=1, k1+1
        SOz(k) = SO(3,3,k,kb)
        SOx(k) = SO(3,3,k,kpw)
        SOy(k) = SO(3,3,k,kps)
!       IF ( myproc == 1 ) write(0,*) 'k, SOz,x,y = ',k,SOz(k),
!     &  SO(3,3,k,kpw),SO(3,3,k,kps)
!        IF ( myproc == 1 ) write(0,*) 'j, SOy = ',
      ENDDO
!        SOz(2) = SOz(3)
!        SOz(k1+1) = SOz(k1)


      DO j=1, j1+1
!        SOy(j) = SO(3,j,3,kps)
!        IF ( myproc == 1 ) write(0,*) 'j, SOy,x,z = ',j,SOy(j),
!     &   SO(3,j,3,kpw),SO(3,j,3,kb)
!       IF ( myproc == 1 ) write(0,*) 'k, SOz = ',j,SO(3,j,3,kb)
!        IF ( myproc == 1 ) write(0,*) 'i, SOx = ',j,SO(3,j,3,kpw)
!        write(0,*) 'j, SOy, myproc = ',j,SOy(j),myproc
      ENDDO
!        SOy(2) = SOy(3)
!        SOy(j1+1) = SOy(j1)


      DO i=1,i1+1
!        SOx(i) = SO(i,3,3,kpw)
!        IF ( myproc == 1 ) write(0,*) 'i, SOx,y,z = ',i,SOx(i),
!     &   SO(i,3,3,kps),SO(i,3,3,kb)
      ENDDO
!        SOx(2) = SOx(3)
!        SOx(i1+1) = SOx(i1)

      ENDIF

         DO pts=ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<

          IF ( ISEP .eq. BMG_SEPARABLE ) THEN

          IF ( .true. ) THEN


            DO  k=2, k1
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !
                     Q(i,j,k) = ( QF(i,j,k) &
     &                        + SOx(k)*Q(i-1,j,k) &
     &                        + SOy(k)*Q(i,j+1,k) &
     &                        + SOx(k)*Q(i+1,j,k) &
     &                        + SOy(k)*Q(i,j-1,k) &
     &                        + SOz(k)*Q(i,j,k-1) &
     &                        + SOz(k+1)*Q(i,j,k+1) &
     &                        )*SOR(i,j,k,msor)
                  END DO
                  !
               END DO
               !
            END DO

          ELSE ! old try

          kbot = 0
          kt = 0
          js = 0
          jn = 0
          iw = 0
          ie = 0

          IF (  MyProcK .eq. 1 )     kbot = 1
          IF (  MyProcK .eq. NProcK) kt = 1
          IF (  MyProcJ .eq. 1)      js = 1
          IF (  MyProcJ .eq. ntasks_y) jn = 1
          IF (  MyProcI .eq. 1)      iw = 1
          IF (  MyProcI .eq. ntasks_x) ie = 1

            DO  k=2+kbot, k1-kt
               DO j=2+js, j1-jn
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg

                  IF ( ibeg .eq. 2 ) ibeg = ibeg + 2*iw
                  IF ( iend .eq. i1 ) iend = iend - 2*ie

                  !
                  DO i=ibeg,iend,2
                  !
!                  IF ( .not. (( (k.eq.2      .and. MyProcK.eq.1 ) .or.
!     &                        ( k .eq. k1    .and. MyProcK .eq. NProcK) .or.
!     &                        ( j .eq. 2     .and. MyProcJ .eq. 1)  .or.
!     &                        ( j .eq. j1    .and. MyProcJ .eq. ntasks_y) .or.
!     &                        ( i .eq. ibeg  .and. MyProcI .eq. 1)  .or.
!     &                        ( i .eq. iend  .and. MyProcI .eq. ntasks_x) ) ) THEN
                     Q(i,j,k) = ( QF(i,j,k) &
     &                        + SOx(i)*Q(i-1,j,k) &
     &                        + SOy(j+1)*Q(i,j+1,k) &
     &                        + SOx(i+1)*Q(i+1,j,k) &
     &                        + SOy(j)*Q(i,j-1,k) &
     &                        + SOz(k)*Q(i,j,k-1) &
     &                        + SOz(k+1)*Q(i,j,k+1) &
     &                        )/(SOx(i) + SOx(i+1) + SOy(j) + &
     &                           SOy(j+1) + SOz(k) + SOz(k+1))

                  END DO
                  !
               END DO
               !
            END DO

! edges:

           IF ( kbot .eq. 1 ) THEN
            k=2

               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2

                     Q(i,j,k) = ( QF(i,j,k) &
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k) &
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k) &
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k) &
     &                        + SO(i,j,k,kps)*Q(i,j-1,k) &
     &                        + SO(i,j,k,kb)*Q(i,j,k-1) &
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1) &
     &                        )*SOR(i,j,k,msor)


                 END DO
              END DO

           ENDIF

           IF ( kt .eq. 1 ) THEN
            k=k1

               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2

                     Q(i,j,k) = ( QF(i,j,k) &
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k) &
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k) &
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k) &
     &                        + SO(i,j,k,kps)*Q(i,j-1,k) &
     &                        + SO(i,j,k,kb)*Q(i,j,k-1) &
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1) &
     &                        )*SOR(i,j,k,msor)


                 END DO
              END DO

           ENDIF

           IF ( js .eq. 1 ) THEN

            DO  k=2+kbot, k1-kt
              j = 2
              ibeg = mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
              iend = 2*((i1-ibeg)/2)+ibeg

               DO i=ibeg,iend,2

                     Q(i,j,k) = ( QF(i,j,k) &
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k) &
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k) &
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k) &
     &                        + SO(i,j,k,kps)*Q(i,j-1,k) &
     &                        + SO(i,j,k,kb)*Q(i,j,k-1) &
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1) &
     &                        )*SOR(i,j,k,msor)

                END DO
             END DO
          ENDIF


           IF ( jn .eq. 1 ) THEN

            DO  k=2+kbot, k1-kt
              j = j1
              ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
              iend=2*((i1-ibeg)/2)+ibeg
                  DO i=ibeg,iend,2


                     Q(i,j,k) = ( QF(i,j,k) &
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k) &
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k) &
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k) &
     &                        + SO(i,j,k,kps)*Q(i,j-1,k) &
     &                        + SO(i,j,k,kb)*Q(i,j,k-1) &
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1) &
     &                        )*SOR(i,j,k,msor)

                END DO
             END DO
          ENDIF


           IF ( ie .eq. 1 ) THEN
            DO  k=2+kbot, k1-kt
               DO j=2+js, j1-jn

              ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
              iend=2*((i1-ibeg)/2)+ibeg
              i = ibeg
               IF ( ibeg .eq. 2 ) THEN

                     Q(i,j,k) = ( QF(i,j,k) &
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k) &
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k) &
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k) &
     &                        + SO(i,j,k,kps)*Q(i,j-1,k) &
     &                        + SO(i,j,k,kb)*Q(i,j,k-1) &
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1) &
     &                        )*SOR(i,j,k,msor)

                ENDIF
                END DO
             END DO
          ENDIF

           IF ( iw .eq. 1 ) THEN
            DO  k=2+kbot, k1-kt
               DO j=2+js, j1-jn

              ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
              iend=2*((i1-ibeg)/2)+ibeg

              i = iend

               IF ( i .eq. i1 ) THEN

                     Q(i,j,k) = ( QF(i,j,k) &
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)   &
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)  &
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k) &
     &                        + SO(i,j,k,kps)*Q(i,j-1,k) &
     &                        + SO(i,j,k,kb)*Q(i,j,k-1) &
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1) &
     &                        )*SOR(i,j,k,msor)

                ENDIF
                END DO
             END DO
          ENDIF

         ENDIF


         ELSE  ! Non separable

            DO  k=2, k1
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                  IF ( Qmask(i,j,k) .gt. 0.9e30 ) THEN
                     Q(i,j,k) = ( QF(i,j,k)   &
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k) &
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)  &
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)  &
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)  &
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)  &
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)  &
     &                        )*SOR(i,j,k,msor)
                  ELSE
                      Q(i,j,k) = Qmask(i,j,k)
!                      QF(i,j,k) =  Qmask(i,j,k,1)*SO(i,j,k,kp) -
                      QF(i,j,k) =  Qmask(i,j,k)/SOR(i,j,k,msor) - &
     &                        ( SO(i,j,k,kpw)*Q(i-1,j,k) &
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k) &
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k) &
     &                        + SO(i,j,k,kps)*Q(i,j-1,k) &
     &                        + SO(i,j,k,kb)*Q(i,j,k-1) &
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1) &
     &                        )

                  ENDIF
                  END DO
                  !
               END DO
               !
            END DO

          ENDIF
            !

            IF ( comhalf ) THEN

            ptrn = 1

            CALL MSG_tbdx_send(Q, MSG_Buffer, &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)), &
     &           iWork(pMSG(ipL_MSG_Proc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),  &
     &           iWork(pMSG(ipL_MSG_Index,KG)),  &
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,  &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Proc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),  &
     &           iWork(pMSG(ipL_MSG_Index,KG)),  &
     &           ptrn, ierror)
            !

            ENDIF


            !
         END DO   ! >>> END: loop over colors <<<<<<<<<<<<<<<<<<<<<<<<<<
            !
            IF ( .not. comhalf ) THEN

            ptrn = 1

            CALL MSG_tbdx_send(Q, MSG_Buffer,   &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Proc,KG)),   &
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),   &
     &           iWork(pMSG(ipL_MSG_Index,KG)),  &
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,  &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)), &
     &           iWork(pMSG(ipL_MSG_Proc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),  &
     &           iWork(pMSG(ipL_MSG_Index,KG)),  &
     &           ptrn, ierror)
            !

            ENDIF
            !
!         END DO   ! >>> END: loop over colors <<<<<<<<<<<<<<<<<<<<<<<<<<

      ENDIF

      ptrn = 1

      CALL MSG_tbdx_close(Q, MSG_Buffer,   &
     &     iWork(pMSG(ipL_MSG_NumAdjProc,KG)), &
     &     iWork(pMSG(ipL_MSG_Proc,KG)),   &
     &     iWork(pMSG(ipL_MSG_Ipr,KG)),   &
     &     iWork(pMSG(ipL_MSG_Index,KG)),  &
     &     ptrn, ierror)

      IF( BMG_IOFLAG(iBMG3_BUG_RES_RELAX) ) THEN
         !
         CALL BMG3_SymStd_residual(    &
     &             KG, NOG, IFD,      &
     &             Q, QF, SO, RES, NLx, NLy, NLz, NStncl, &
     &             iWork, NMSGi, pMSG,   &
     &             MSG_Buffer, NMSGr, MPICOMM &
     &             )
         CALL BMG3_SymStd_UTILS_norm_l2(   &
     &             RES, NLx, NLy, NLz, RES_L2, MPICOMM   &
     &             )
         !
      ENDIF

! ==========================================================================

      RETURN
      END SUBROUTINE WRF_SymStd_relax_GS

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE WRF_SymStd_interp_add(   &
     &                KCG, KFG, NOG,    &
     &                Q ,QC, RES, SO, NStncl, CI,  &
     &                IIC, JJC, KKC,   &
     &                IIF, JJF, KKF,  &
     &                iGs, jGs, kGs,  &
     &                iWork, NMSGi, pMSG, Qmask, &
     &                BUFFER, NMSGr, MPICOMM    &
     &                )

! ==========================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!   WRF_SymStd_interp_add.f interpolates Q from the coarse mesh, KCG, to
!   the fine mesh, KFG, and adds the result to Q on fine mesh.
!
!  Special version for finest mesh: checks mask array and does not update Q
!  in block around lightning flash
!
! =======================================================================
! $license_flag$
! =======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   INPUT/OUTPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   OUTPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   LOCAL:
!  --------------------
!
!
! ==========================================================================

      IMPLICIT NONE

! -----------------------------
!     Includes
!
      INCLUDE 'mpif.h'
      INCLUDE 'MSG_f90.h'

      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'

! ----------------------------
!     Argument Declarations
!
      INTEGER IIC, IIF, JJC, JJF, KCG, KFG, KKC, KKF
      INTEGER iGs, jGs, kGs, NMSGi, NMSGr, NOG, NStncl

      INTEGER iWork(NMSGi), pMSG(NBMG_pMSG,NOG), MPICOMM
      REAL*RKIND  CI(IIC,JJC,KKC,26), Q(IIF,JJF,KKF), QC(IIC,JJC,KKC)
      REAL*RKIND  RES(IIF,JJF,KKF), SO(IIF+1,JJF+1,KKF+1,NStncl)
      REAL*RKIND  BUFFER(NMSGr)

      REAL*RKIND :: Qmask(IIF,JJF,KKF)

! ----------------------------
!     Local Declarations
!
      INTEGER IC, I, IICF, IICF1, IIC1,  &
     &        JC, J, JJCF, JJCF1, JJC1,  &
     &        KC, K, KKCF, KKCF1, KKC1
      INTEGER ISTART, JSTART, KSTART, ICSTART, JCSTART, KCSTART
      INTEGER ierror, ptrn
      REAL*RKIND  A, AQ

      INTEGER imn,imx, jmn,jmx, kmn,kmx
      LOGICAL flash
      REAL    maskmn, maskmx

! ==========================================================================

! -------------------------------------------------
!     Useful index bounds:
! -------------------------------------------------

      iic1=iic-1
      jjc1=jjc-1
      kkc1=kkc-1

      iicf=(iif-2)/2+3
      jjcf=(jjf-2)/2+3
      kkcf=(kkf-2)/2+3

      IF ((mod(IIF,2).eq.0).and.(mod(iGs,2).eq.0)) THEN
         IICF = IICF - 1
      ENDIF

      IF ((mod(JJF,2).eq.0).and.(mod(jGs,2).eq.0)) THEN
         JJCF = JJCF - 1
      ENDIF

      IF ((mod(KKF,2).eq.0).and.(mod(kGs,2).eq.0)) THEN
         KKCF = KKCF - 1
      ENDIF

      IICF1=IICF-1
      JJCF1=JJCF-1
      KKCF1=KKCF-1

! --------------------------------------------------
!   interpolate answers from coarse to fine mesh
!   and add to answers on fine mesh.
! --------------------------------------------------

!     odd values of iGs, jGs, or kGs indicate the
!     the case similar to the serial version


      IF (mod(iGs,2).eq.1) THEN
         ISTART = 2
         ICSTART = 2
      ELSE
         ISTART = 1
         ICSTART = 1
      ENDIF

      IF (mod(jGs,2).eq.1) THEN
         JSTART = 2
         JCSTART = 2
      ELSE
         JSTART = 1
         JCSTART = 1
      ENDIF

      IF (mod(kGs,2).eq.1) THEN
         KSTART = 2
         KCSTART = 2
      ELSE
         KSTART = 1
         KCSTART = 1
      ENDIF

! --------------------------------------------------
!     NB: division is only possible in the interior
! --------------------------------------------------

      imn = iif
      imx = 0
      jmn = jjf
      jmx = 0
      kmn = kkf
      kmx = 0
      flash = .false.

      maskmn = 1.e31
      maskmx = 0

      do k=2, kkf-1
         do j=2, jjf-1
            do i=2, iif-1
               RES(i,j,k)=RES(i,j,k)/so(i,j,k,kp)

               maskmn = Min( maskmn, Qmask(i,j,k) )
               maskmx = Max( maskmx, Qmask(i,j,k) )

               IF ( Qmask(i,j,k) .lt. 0.9e30 ) THEN
                 flash = .true.
                 imn = Min( imn, i-1 )
                 jmn = Min( jmn, j-1 )
                 kmn = Min( kmn, k-1 )

                 imx = Max( imx, i+1 )
                 jmx = Max( jmx, j+1 )
                 kmx = Max( kmx, k+1)
!           write(0,*) 'WRF_SymStd_interp_add: ijk,qmask= ',i,j,k,Qmask(i,j,k)
!            write(0,*) 'WRF_SymStd_interp_add: imn: ',imn,imx
               ENDIF
            enddo
         enddo
      enddo

      IF ( imn .le. imx ) THEN
!      write(0,*) 'flash detected in WRF_SymStd_interp_add'
      flash = .true.
      imn = Max(1, imn-3)
      jmn = Max(1, jmn-3)
      kmn = Max(1, kmn-3)

      imx = Max(1, imx+3)
      jmx = Max(1, jmx+3)
      kmx = Max(1, kmx+3)

      ELSE
!      write(0,*) 'no flash detected in WRF_SymStd_interp_add, mx,mn=',maskmn,maskmx
!      write(0,*) 'WRF_SymStd_interp_add: iif,jjf,kkf = ',iif,jjf,kkf
      ENDIF


      IF ( .not. flash ) THEN

      k=KSTART-2
      j=JSTART
      i=ISTART

      do kc=KCSTART,kkc1

         i=ISTART
         j=JSTART
         k=k+2

         ic=ICSTART
         jc=JCSTART

         q(i,j,k)=q(i,j,k)+qc(ic,jc,kc)

         do ic=ICSTART+1,iicf1
            i=i+2
            q(i,j,k) = q(i,j,k) + qc(ic,jc,kc)

            a = ci(ic,jc,kc,lxyr)*qc(ic,jc,kc) &
     &           + ci(ic,jc,kc,lxyl)*qc(ic-1,jc,kc)

            q(i-1,j,k) = q(i-1,j,k) + a + RES(i-1,j,k)

         enddo

         do jc=JCSTART+1,jjcf1

            i=ISTART
            j=j+2

            ic=ICSTART

            q(i,j,k) = q(i,j,k) + qc(ic,jc,kc)

            aq = ci(2,jc,kc,lxya)*qc(2,jc,kc) &
     &         + ci(2,jc,kc,lxyb)*qc(2,jc-1,kc)

            q(i,j-1,k) = q(i,j-1,k) + aq + RES(i,j-1,k)

            do ic=ICSTART+1,iicf1

               i=i+2

               q(i,j,k) = q(i,j,k) + qc(ic,jc,kc)

               a = ci(ic,jc,kc,lxyr)*qc(ic,jc,kc) &
     &              + ci(ic,jc,kc,lxyl)*qc(ic-1,jc,kc)

               q(i-1,j,k) = q(i-1,j,k) + a + RES(i-1,j,k)

               aq = ci(ic,jc,kc,lxya)*qc(ic,jc,kc) &
     &              + ci(ic,jc,kc,lxyb)*qc(ic,jc-1,kc)

               q(i,j-1,k) = q(i,j-1,k) + aq + RES(i,j-1,k)

               a = ci(ic,jc,kc,lxysw)*qc(ic-1,jc-1,kc) &
     &              + ci(ic,jc,kc,lxynw)*qc(ic-1,jc,kc) &
     &              + ci(ic,jc,kc,lxyne)*qc(ic,jc,kc) &
     &              + ci(ic,jc,kc,lxyse)*qc(ic,jc-1,kc)

               q(i-1,j-1,k) = q(i-1,j-1,k) + a + RES(i-1,j-1,k)

            enddo
         enddo
      enddo


      k=KSTART-1
      do kc=KCSTART+1,kkcf1

         i=ISTART-2
         j=JSTART
         k=k+2

         jc=JCSTART

         do ic=ICSTART,iic1

            i=i+2

            q(i,j,k) = q(i,j,k) + ci(ic,jc,kc,lxza)*qc(ic,jc,kc) &
     &           + ci(ic,jc,kc,lxzb)*qc(ic,jc,kc-1) + RES(i,j,k)

         enddo

         j=JSTART
         do jc=JCSTART+1,jjcf1

            i=ISTART-2
            j=j+2

            do ic=ICSTART,iic1

               i=i+2

               q(i,j,k) = q(i,j,k) &
     &                  + ci(ic,jc,kc,lxza)*qc(ic,jc,kc) &
     &                  + ci(ic,jc,kc,lxzb)*qc(ic,jc,kc-1) &
     &                  + RES(i,j,k)

               q(i,j-1,k) = q(i,j-1,k) &
     &                    + ci(ic,jc,kc,lyznw)*qc(ic,jc,kc) &
     &                    + ci(ic,jc,kc,lyzne)*qc(ic,jc-1,kc) &
     &                    + ci(ic,jc,kc,lyzsw)*qc(ic,jc,kc-1) &
     &                    + ci(ic,jc,kc,lyzse)*qc(ic,jc-1,kc-1) &
     &                    + RES(i,j-1,k)

            enddo
         enddo
      enddo


      k=KSTART-1
      do kc=KCSTART+1,kkcf1
         k=k+2

         j=JSTART
         jc=JCSTART
         i=ISTART-1

         do ic=ICSTART+1,iicf1
            i=i+2

            q(i,j,k) = q(i,j,k) &
     &               + ci(ic,jc,kc,lxznw)*qc(ic-1,jc,kc) &
     &               + ci(ic,jc,kc,lxzne)*qc(ic,jc,kc) &
     &               + ci(ic,jc,kc,lxzsw)*qc(ic-1,jc,kc-1) &
     &               + ci(ic,jc,kc,lxzse)*qc(ic,jc,kc-1) &
     &               + RES(i,j,k)

         enddo

         j=JSTART
         do jc=JCSTART+1,jjcf1
            j=j+2

            i=ISTART-1
            do ic=ICSTART+1,iicf1
               i=i+2

               q(i,j,k) = q(i,j,k)  &
     &                  + ci(ic,jc,kc,lxznw)*qc(ic-1,jc,kc) &
     &                  + ci(ic,jc,kc,lxzne)*qc(ic,jc,kc) &
     &                  + ci(ic,jc,kc,lxzsw)*qc(ic-1,jc,kc-1) &
     &                  + ci(ic,jc,kc,lxzse)*qc(ic,jc,kc-1) &
     &                  + RES(i,j,k)

               q(i,j-1,k) = q(i,j-1,k) &
     &                    + ci(ic,jc,kc,ltnw)*qc(ic-1,jc,kc) &
     &                    + ci(ic,jc,kc,ltne)*qc(ic,jc,kc) &
     &                    + ci(ic,jc,kc,ltsw)*qc(ic-1,jc-1,kc) &
     &                    + ci(ic,jc,kc,ltse)*qc(ic,jc-1,kc) &
     &                    + ci(ic,jc,kc,lbnw)*qc(ic-1,jc,kc-1) &
     &                    + ci(ic,jc,kc,lbne)*qc(ic,jc,kc-1) &
     &                    + ci(ic,jc,kc,lbsw)*qc(ic-1,jc-1,kc-1) &
     &                    + ci(ic,jc,kc,lbse)*qc(ic,jc-1,kc-1) &
     &                    + RES(i,j-1,k)
            enddo

         enddo

      enddo

      ELSE ! there are channel points in the tile
      k=KSTART-2
      j=JSTART
      i=ISTART

      do kc=KCSTART,kkc1

         i=ISTART
         j=JSTART
         k=k+2

         ic=ICSTART
         jc=JCSTART

         IF ( .not. ( i .ge. imn .and. i .le. imx .and. &
     &                j .ge. jmn .and. j .le. jmx .and. &
     &                k .ge. kmn .and. k .le. kmx ) ) THEN

         q(i,j,k)=q(i,j,k)+qc(ic,jc,kc)

         ENDIF

         do ic=ICSTART+1,iicf1

            i=i+2

         IF ( .not. ( i .ge. imn .and. i .le. imx .and.  &
     &                j .ge. jmn .and. j .le. jmx .and.   &
     &                k .ge. kmn .and. k .le. kmx ) ) THEN

            q(i,j,k) = q(i,j,k) + qc(ic,jc,kc)

            a = ci(ic,jc,kc,lxyr)*qc(ic,jc,kc)  &
     &           + ci(ic,jc,kc,lxyl)*qc(ic-1,jc,kc)

            q(i-1,j,k) = q(i-1,j,k) + a + RES(i-1,j,k)

         ENDIF
         enddo

         do jc=JCSTART+1,jjcf1

            i=ISTART
            j=j+2

            ic=ICSTART

            IF ( .not. ( i .ge. imn .and. i .le. imx .and. &
     &                   j .ge. jmn .and. j .le. jmx .and. &
     &                   k .ge. kmn .and. k .le. kmx ) ) THEN

            q(i,j,k) = q(i,j,k) + qc(ic,jc,kc)

            aq = ci(2,jc,kc,lxya)*qc(2,jc,kc)   &
     &         + ci(2,jc,kc,lxyb)*qc(2,jc-1,kc)

            q(i,j-1,k) = q(i,j-1,k) + aq + RES(i,j-1,k)

            ENDIF

            do ic=ICSTART+1,iicf1

               i=i+2

             IF ( .not. ( i .ge. imn .and. i .le. imx .and. &
     &                    j .ge. jmn .and. j .le. jmx .and. &
     &                    k .ge. kmn .and. k .le. kmx ) ) THEN
               q(i,j,k) = q(i,j,k) + qc(ic,jc,kc)
               a = ci(ic,jc,kc,lxyr)*qc(ic,jc,kc) &

     &              + ci(ic,jc,kc,lxyl)*qc(ic-1,jc,kc)

               q(i-1,j,k) = q(i-1,j,k) + a + RES(i-1,j,k)

               aq = ci(ic,jc,kc,lxya)*qc(ic,jc,kc)  &
     &              + ci(ic,jc,kc,lxyb)*qc(ic,jc-1,kc)

               q(i,j-1,k) = q(i,j-1,k) + aq + RES(i,j-1,k)

               a = ci(ic,jc,kc,lxysw)*qc(ic-1,jc-1,kc)  &
     &              + ci(ic,jc,kc,lxynw)*qc(ic-1,jc,kc)  &
     &              + ci(ic,jc,kc,lxyne)*qc(ic,jc,kc)  &
     &              + ci(ic,jc,kc,lxyse)*qc(ic,jc-1,kc)

               q(i-1,j-1,k) = q(i-1,j-1,k) + a + RES(i-1,j-1,k)

               ENDIF

            enddo
         enddo
      enddo


      k=KSTART-1
      do kc=KCSTART+1,kkcf1

         i=ISTART-2
         j=JSTART
         k=k+2

         jc=JCSTART

         do ic=ICSTART,iic1

            i=i+2

           IF ( .not. ( i .ge. imn .and. i .le. imx .and. &
     &                  j .ge. jmn .and. j .le. jmx .and.  &
     &                  k .ge. kmn .and. k .le. kmx ) ) THEN

            q(i,j,k) = q(i,j,k) + ci(ic,jc,kc,lxza)*qc(ic,jc,kc)   &
     &           + ci(ic,jc,kc,lxzb)*qc(ic,jc,kc-1) + RES(i,j,k)

            ENDIF
         enddo

         j=JSTART
         do jc=JCSTART+1,jjcf1

            i=ISTART-2
            j=j+2

            do ic=ICSTART,iic1

               i=i+2

            IF ( .not. ( i .ge. imn .and. i .le. imx .and. &
     &                   j .ge. jmn .and. j .le. jmx .and.  &
     &                   k .ge. kmn .and. k .le. kmx ) ) THEN
               q(i,j,k) = q(i,j,k)     &
     &                  + ci(ic,jc,kc,lxza)*qc(ic,jc,kc)  &
     &                  + ci(ic,jc,kc,lxzb)*qc(ic,jc,kc-1)  &
     &                  + RES(i,j,k)

               q(i,j-1,k) = q(i,j-1,k)     &
     &                    + ci(ic,jc,kc,lyznw)*qc(ic,jc,kc)    &
     &                    + ci(ic,jc,kc,lyzne)*qc(ic,jc-1,kc)   &
     &                    + ci(ic,jc,kc,lyzsw)*qc(ic,jc,kc-1)    &
     &                    + ci(ic,jc,kc,lyzse)*qc(ic,jc-1,kc-1)   &
     &                    + RES(i,j-1,k)

              ENDIF
            enddo
         enddo
      enddo


      k=KSTART-1
      do kc=KCSTART+1,kkcf1
         k=k+2

         j=JSTART
         jc=JCSTART
         i=ISTART-1

         do ic=ICSTART+1,iicf1
            i=i+2

         IF ( .not. ( i .ge. imn .and. i .le. imx .and.  &
     &                j .ge. jmn .and. j .le. jmx .and. &
     &                k .ge. kmn .and. k .le. kmx ) ) THEN
            q(i,j,k) = q(i,j,k)  &
     &               + ci(ic,jc,kc,lxznw)*qc(ic-1,jc,kc)  &
     &               + ci(ic,jc,kc,lxzne)*qc(ic,jc,kc)  &
     &               + ci(ic,jc,kc,lxzsw)*qc(ic-1,jc,kc-1)  &
     &               + ci(ic,jc,kc,lxzse)*qc(ic,jc,kc-1)  &
     &               + RES(i,j,k)

          ENDIF
         enddo

         j=JSTART
         do jc=JCSTART+1,jjcf1
            j=j+2

            i=ISTART-1
            do ic=ICSTART+1,iicf1
               i=i+2

         IF ( .not. ( i .ge. imn .and. i .le. imx .and.  &
     &                j .ge. jmn .and. j .le. jmx .and.  &
     &                k .ge. kmn .and. k .le. kmx ) ) THEN
               q(i,j,k) = q(i,j,k)   &
     &                  + ci(ic,jc,kc,lxznw)*qc(ic-1,jc,kc)  &
     &                  + ci(ic,jc,kc,lxzne)*qc(ic,jc,kc)  &
     &                  + ci(ic,jc,kc,lxzsw)*qc(ic-1,jc,kc-1)  &
     &                  + ci(ic,jc,kc,lxzse)*qc(ic,jc,kc-1)  &
     &                  + RES(i,j,k)

               q(i,j-1,k) = q(i,j-1,k)  &
     &                    + ci(ic,jc,kc,ltnw)*qc(ic-1,jc,kc)  &
     &                    + ci(ic,jc,kc,ltne)*qc(ic,jc,kc)  &
     &                    + ci(ic,jc,kc,ltsw)*qc(ic-1,jc-1,kc)  &
     &                    + ci(ic,jc,kc,ltse)*qc(ic,jc-1,kc)  &
     &                    + ci(ic,jc,kc,lbnw)*qc(ic-1,jc,kc-1)  &
     &                    + ci(ic,jc,kc,lbne)*qc(ic,jc,kc-1)  &
     &                    + ci(ic,jc,kc,lbsw)*qc(ic-1,jc-1,kc-1)  &
     &                    + ci(ic,jc,kc,lbse)*qc(ic,jc-1,kc-1)  &
     &                    + RES(i,j-1,k)
            ENDIF

            enddo

         enddo

      enddo

      ENDIF





      ptrn = 1

      call MSG_tbdx_send(Q, buffer,  &
     &     iWork(pMSG(ipL_MSG_NumAdjProc,KFG)),  &
     &     iWork(pMSG(ipL_MSG_Proc,KFG)),  &
     &     iWork(pMSG(ipL_MSG_Ipr,KFG)),  &
     &     iWork(pMSG(ipL_MSG_Index,KFG)),  &
     &     ptrn, ierror)

      call MSG_tbdx_receive(Q, buffer,  &
     &     iWork(pMSG(ipL_MSG_NumAdjProc,KFG)),  &
     &     iWork(pMSG(ipL_MSG_Proc,KFG)),  &
     &     iWork(pMSG(ipL_MSG_Ipr,KFG)),  &
     &     iWork(pMSG(ipL_MSG_Index,KFG)),  &
     &     ptrn, ierror)

      call MSG_tbdx_close(Q, buffer,  &
     &     iWork(pMSG(ipL_MSG_NumAdjProc,KFG)),  &
     &     iWork(pMSG(ipL_MSG_Proc,KFG)),  &
     &     iWork(pMSG(ipL_MSG_Ipr,KFG)),  &
     &     iWork(pMSG(ipL_MSG_Index,KFG)),  &
     &     ptrn, ierror)

! ==========================================================================

      RETURN
      END SUBROUTINE WRF_SymStd_interp_add



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE WRF_SymStd_SOLVE_boxmgdg(  &
     &                NLx, NLy, NLz, NGx, NGy, NGz, iGs, jGs, kGs,  &
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,  &
     &                Q, QF, RES, NFm, NCm,  &
     &                SO, NSOm,  &
     &                SOR, NSORm,  &
     &                CI, NCIm,  &
     &                IGRD, NOGm, NOG,  &
     &                BMG_iWORK_PL, NBMG_iWORK_PLm,  &
     &                BMG_rWORK_PL, NBMG_rWORK_PLm,  &
     &                BMG_iWORK_CS, NBMG_iWORK_CSm,  &
     &                BMG_rWORK_CS, NBMG_rWORK_CSm,  &
     &                iWorkMSG, NMSGim,  &
     &                pMSG, pMSGSO,  &
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,  &
     &                NProc, MSG_Buffer, NMSGrm,  &
     &                MPICOMM)
!     &                ,hx,hy,hz)

! ==========================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG3_SymStd_SOLVE_boxmg is a black box matrix solver. It takes a
!     matrix defined by the user from a given (fine) grid and constructs
!     coarser grids and their associated coefficient matrices. BMG3D
!     then performs the multigrid algorithm and returns the solution
!     vec- tor. The matrix is set by a user written subroutine. The
!     difference stencil at (i,j,k) spans three planes. The k-plane has
!     difference stencil
!
!                        xnw xn xne
!                        xw  xo xe
!                        xsw xs xse.
!
!     The (k+1)-plane has difference stencil
!
!                        xtnw xtn xtne
!                        xtw  xt  xte
!                        xtsw xts xtse.
!
!     The(k-1)-plane has difference stencil
!
!                        xbnw xbn xbne
!                        xbw  xb  xbe
!                        xbsw xbs xbse.
!
!     Here
!           xo   = SO(i,j,k,ko)
!           xw   = - SO(i,j,k,kw)
!           xs   = - SO(i,j,k,ks)
!           xb   = - SO(i,j,k,kb)
!           xsw  = - SO(i,j,k,ksw)
!           xnw  = - SO(i,j+1,k,knw)
!           xbw  = - SO(i,j,k,kbw)
!           xbnw = - SO(i,j+1,k,kbnw)
!           xbn  = - SO(i,j+1,k,kbn)
!           xbne = - SO(i+1,j+1,k,kbne)
!           xbe  = - SO(i+1,j,k,kbe)
!           xbse = - SO(i+1,j,k,kbse)
!           xbs  = - SO(i,j,k,kbs)
!           xbsw = - SO(i,j,k,kbsw)
!
!     The remaining coefficents:
!
!     xn, xne,xe,xse, xt, xtw, xtnw, xtn, xtne, xte, xtse, xts, and xtsw
!
!     are specified by symmetry. In the case of a seven point operator,
!     only xo, xw, xs, and xb need to be set. The difference scheme is
!     assumed to be positive definite.
!
!     The fictitious points
!
!         (1,j,k), j=1,jj, k=1,kk, (ii,j,k), j=1,jj, k=1,kk,
!         (i,1,k), i=1,ii, k=1,kk, (i,jj,k), i=1,ii, k=1,kk,
!         (i,j,1), i=1,ii, j=1,jj, (i,j,kk), i=1,ii, j=1,jj,
!
!     are assumed for ease of programming. However, it is also assumed
!     that the user has set any difference coefficient referring
!     to these points to zero. For example, if the above difference
!     template is centered at (2,2), then xb, xbw, xbnw, xbn, xbne, xbe,
!     xbse, xbs xnw, xw, xsw, xs, xse, xtnw, xtw, xtsw, xts, xtse
!     would be assumed zero.
!
! =======================================================================
! $license_flag$
! =======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!     NXM      x-dimension of the grid, excluding fictitious points
!     NYM      y-dimension of the grid, excluding fictitious points
!     NZM      z-dimension of the grid, excluding fictitious points
!
!     TOL      Convergence tolerance.
!
!     ISTOP    Indicates the stopping criteria. At present there
!              is a relative residual test
!
!                ISTOP = BMG_STOP_REL_RES_L2
!
!                => iteration has converged if the l2-norm of the
!                   current residual, divided by the l2-norm of the
!                   initial residual is less than TOL.
!
!              and an absolute residual test
!
!                ISTOP = BMG_STOP_REL_ABS_L2
!
!                => iteration has converged if the l2-norm of the
!                  current residual is less than TOL.
!
!              NOTE: if the selected test fails and the number of
!              iterations is greater than ABS(ISTRT) then the iteration
!              has failed => set TOL=-RES_L2 and return.
!
!
!     IFD      Indicator for difference scheme. IFD=1 means a seven
!              point operator. I.e., in the template given above, xnw
!              xne, xse, xsw, xbw, xbnw, xbn, xbne, xbe, xbse, xbs, xbsw,
!              xtw, xtnw, xtn, xtne, xte, xtse, xts, xtsw.
!              are assumed zero. IFD.ne.1. means a 27 point
!              operator. (Note that 27 point operators will be
!              generated on the coarser grids even in the case of a seven
!              point operator.
!
!     IU       Number of relaxation sweeps to be performed
!              on a coarse grid before interpolation to a fine grid.
!
!     ID       Number of relaxation sweeps to be performed on a fine
!              grid before the problem is transferred to a coarse grid.
!
!     IM       Unused
!
!     ISTRT    Indicator for number of multigrid cycles to be performed
!              if ISTRT.gt.0 the algorithm will begin on the coarsest grid
!              with a full multigrid cycle (however,
!              without cubic interpolation). It will continue
!              cycling until ISTRT cycles have been performed or until
!              the error criterion is satisfied. If ISTRT.lt.0 the
!              algorithm will begin on the finest grid and will perform
!              -ISTRT cycles unless the error criterion is satisfied
!              first.
!
!     ISETUP   Indicator for whether or not to skip initial setup,
!              initial guess for Q, computation of pointers, and
!              generation of coefficients on the coarser grids.
!
!                ISETUP=0  => perform setup, including,
!                             - compute workspace pointers
!                             - compute coarse-grid operators
!                             - compute interpolation transfer operators
!
!              If ISETUP.NE.0 then we have the following cases:
!
!                ISETUP=1  => skip computation of workspace pointers
!                             - compute coarse-grid operators
!                             - compute interpolation operators
!                ISETUP=2  => skip all setup
!
!              It is important to note the interaction that ISETUP
!              has with ISTRT.  Specifically, with ISETUP.NE.0 we have
!
!                 ISTRT > 0 => FMG so residual is coarsened appropriately
!                 ISTRT < 0 => n-cycle so nothing extra to do here.
!
!              An example of usage would be a constant coefficient
!              time-dependent problem where ISETUP=0 on the first time step
!              and ISETUP=2 on subsequent time steps.
!
!     BMG_IOFLAG   Logical array of I/O and debugging flags. Admittedly
!                  this is not well documented yet, but look in the include
!                  file "BMG_parameters.h" for hints (the names are verbose.)
!
!     IRELAX   Indicator for the relaxation algorithm:
!
!              IRELAX=1  => Colored Gauss-Seidel relaxation:
!                        -  red-black on the finest grid for 5-point stencils
!                        -  four color otherwise
!
!              IRELAX=5  => alternating red-black plane Gauss-Seidel relaxation
!                        -  ordering is xy planes, yz planes, xz planes
!                        -  for symmetric relaxation this is reversed
!                           as necessary when moving "up" or "down" levels
!
!     IRELAX_SYM  Symmetric cycles require that the relaxation ordering
!                 be reversed between the restriction and interpolation
!                 stage of the cycle.   This is important when BMG3D
!                 is being used as a preconditioner.
!
!                    IRELAX_SYM = BMG_RELAX_NONSYM
!                    => use the stated relaxation ordering (nonsymmetric)
!
!                    IRELAX_SYM = BMG_RELAX_SYM
!                    => make the cycle symmetric.
!
!                 as defined in BMG_parameters.h
!
!     ISEP     Indicator for separable equation such that the SO array can be
!              reduced to 3 1-D arrays in the interior, saving memory bandwidth
!
!               ISEP = BMG_NONSEPARABLE (0): Always use 3d SO arrays
!               ISEP = BMG_SEPARABLE (1): only use 3d SO arrays on boundaries.
!
!
!     IVW      Indicator for type of cycle. IVW=1 means v-cycles will
!              be performed. IVW=2 means w-cycles will be performed, etc.
!
!     MCYCL    Currently disabled!
!
!     NFm      Maximum storage for a vector on all grids. This should be
!              computed by BMG3_SymStd_SETUP_space and is checked in boxmg.
!
!     NCm      Maximum storage for a vector on all coarse grids.  This should
!              be computed in BMG3_SymStd_SETUP_space and is checked in boxmg.
!
!     NSOm     See INPUT/OUTPUT
!
!     NSORm    See INPUT/OUTPUT
!
!     NCIm     See INPUT/OUTPUT
!
!     NOGm     Maximum number of grids that can be supported by the
!              pointer array IGRD. It should be set in calling program.
!
! =======================================================================
!  --------------------
!   INPUT/OUTPUT:
!  --------------------
!
!     SO       User defined real array that contains the coefficient
!              matrix. See above description for format.  If SETUP is
!              specified then coarse-grid operators are created and returned.
!
!     NSOm     Dimension of SO, the calling program should have computed
!              this with a call to BMG2_SymStd_SETUP_space.  Clearly,
!
!              IFD .EQ. 1  =>    NSOm = 4*NFm + 10*NCm
!              IFD .NE. 1  =>    NSOm = 14*NFm
!
!              BMG3_SymStd_SOLVE_boxmg verifies that NSOm is large enough.
!
!     CI       CI is a real array which should be dimensioned
!              to (NCIm) in calling program. It is used to contain
!              the interpolation coefficients.
!
!     NCIm     Dimension of CI, set in calling program.  NCIm = 26*NCm is
!              enough.  BOXMG checks to see if NCIm is large enough.
!
!     QF       The user defined array that contains the right hand side.
!              It must be dimensioned to at least NFMAX.
!
!     Q        the user defined array that contains the solution vector.
!              It must be dimensioned to at least NFMAX.
!
! =======================================================================
!  --------------------
!   OUTPUT:
!  --------------------
!
!
!
! =======================================================================
!  --------------------
!   WORK ARRAYS:
!  --------------------
!
!     ABD      User declared two dimensional real array, which
!              is used to store the coefficient matrix for the coarsest
!              grid. It is then used by the LAPACK routines. It
!              should be dimensioned to (NCBWm,NCUm).
!
!     BBD      User declared real array of dimension NCUm for use in the
!              the LAPACK routines
!
!     NCBWm    Maximum first subscript of ABD, which needs to be > or =
!              (x-dimension+1)*(y-dimension+2) on coarsest grid.
!
!     NCUm     Maximum second subscript of ABD, which needs to be > or =
!              (x-dimension+1)*(y-dimension+1)*(z-dimension+1) on
!              coarsest grid.
!
!     IGRD     A work array. IGRD should be dimensioned to
!              IGRD(NOGm,NBMG_pIGRD) in the calling program.
!
!     SOR      SOR is a real array, which should be dimensioned to
!              (NSORm) in the calling program. It is used to store
!              residuals, reciprocals of SO(.,1) if IRELAX=1, and lu
!              decompositions if IRELAX.gt.1.
!
!     NSORm    Dimension of SOR, set in calling program. NSORm=2*NFMAX
!              is enough for IRELAX=1 or 2. For IRELAX=3 or 4, 4*NFMAX
!              is enough.
!
!
! =======================================================================
!  --------------------
!   REFERENCES:
!  --------------------
!
!     Dendy, J. E. Jr., "Two Multigrid Methods for Three Dimensional
!     Problems with Discountinuos and Anisotropic Coefficients", SIAM
!     Journal of Scientific and Satatistical Computing, Vol. 8, No. 2,
!     September 1987
!
! ==========================================================================
!  --------------------
!   Error Codes:
!  --------------------
!
!     -------------------------------------------------------------------
!     CODE    DESCRIPTION                           ORIGIN
!     -------------------------------------------------------------------
!      1     initial residual has become zero     SOLVE_boxmg
!      2     NOG = 0                              SETUP_parts
!      3     NOG < 0                              SETUP_parts
!      4     NStncl out of range                  SETUP_cg_LU
!      5     Cholesky decomposition failed, look
!            in BMG_iPARMS(id_BMG3_Ext_Err_Code)
!            for the return code of the LAPACK
!            routine                              SETUP_cg_LU
!      6    NOG = 1                               SETUP_relax
!      7    NOG = 0                               SETUP_relax
!      8    NOG < 0                               SETUP_relax
!      9    fatal setup error                     SETUP_PtrGrid
!     10    IRELAX out of range                   SETUP_PtrGrid
!     11    min coarse grid dim < 3 ( NXYZc )     SETUP_space
!     12    min coarse grid dim < 3 ( NXYc )      SETUP_space
!     13    IRELAX out of range                   SETUP_space
!     14    computed number of grids too small    SETUP_space
!     15    computed number of grids < 1          SETUP_space
!     16    fatal error (check output)            SETUP_PtrWork
!     17    memory allocation mode unspecified    SETUP_PtrWork
!     18    KCF = KC                              ncycle
!     19    updown out of range                   updown
!     20    coarse grid solve failed, look in
!           BMG_iPARMS(id_BMG3_Ext_Err_Code)
!           for the return code of the LAPACK
!           routine                               SOLVE_cg
!     21    inconsistent # of stencil entries     COPY_SO_xy
!     22    inconsistent # of stencil entries     COPY_SO_yz
!     23    inconsistent # of stencil entries     COPY_SO_xz
!     24    invalid cg_operator type              SETUP_parts
!     25    invalid cg_operator construction      SETUP_parts
!
! ==========================================================================

      IMPLICIT NONE

! ------------------------------------------------
!     Includes
!
      INCLUDE 'geom_param_fort90.h'
      INCLUDE 'mpi_param_fort90.h'
      INCLUDE 'MSG_f90.h'

#include    "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

! ------------------------------------------------
!     Argument Declarations
!
      INTEGER  NLx, NLy, NLz, NGx, NGy, NGz, iGs, jGs, kGs,  &
     &         NFm, NCm, NSOm, NSORm, NCIm, NCBWm, NCUm, NOGm, NOG,  &
     &         NMSGim, NMSGrm, NProc, MPICOMM, NBMG_MSG_iGRID,  &
     &         NBMG_iWORK_CSm, NBMG_iWORK_PLm,  &
     &         NBMG_rWORK_CSm, NBMG_rWORK_PLm

      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRID),  &
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)

      INTEGER  iWorkMSG(NMSGim),  &
     &         pMSG(NBMG_pMSG,NOGm),  &
     &         pMSGSO(NBMG_pMSG,NOGm)
      REAL*RKIND   MSG_Buffer(NMSGrm)

      INTEGER  BMG_iWORK_PL(NBMG_iWORK_PLm)
      REAL*RKIND   BMG_rWORK_PL(NBMG_rWORK_PLm)

      INTEGER  BMG_iWORK_CS(NBMG_iWORK_CSm)
      REAL*RKIND   BMG_rWORK_CS(NBMG_rWORK_CSm)

      INTEGER  IGRD(NOGm,NBMG_pIGRD)
      REAL*RKIND CI(NCIm),Q(NFm),QF(NFm),RES(NFm),SO(NSOm),SOR(NSORm)


! ------------------------------------------------
!     Local Declarations
!
      INTEGER i, iBC, id, ifd, iil, iilc, iilf, im, irelax, isep,  &
     &        IRELAX_SYM, ISETUP, ISTOP, istrt, istrt2, iu, ivw,  &
     &        jjl, jjlc, jjlf, k, kc, kfmg, kf, kkl, kklc, kklf,  &
     &        mcyc, mcycl, ncyc, NStncl, NMSGi, NMSGr, p_ABD, p_BBD,  &
     &        p_CI, p_CIC, p_SO, p_SOC, p_SOR, p_SORC, p_U, p_UC,  &
     &        IIG, IIGC, IIGF, JJG, JJGC, JJGF, KKG, KKGC, KKGF,  &
     &        NC, NCBW, NCI, NCU, NF, NSO, NSOR, NXYZc,  &
     &        NBMG_iWORK_CS, NBMG_iWORK_PL, NBMG_rWORK_CS,  &
     &        NBMG_rWORK_PL


      INTEGER IMdum

      INTEGER  iGs_c, iGs_f, jGs_c, jGs_f, kGs_c, kGs_f

      INTEGER  MyProc, MyProcI, MyProcJ, MyProcK,  &
     &         ntasks_x, ntasks_y, NProcK,  &
     &         p_NLx, p_NLy, p_NLz,  &
     &         p_NLx_kg, p_NLy_kg, p_NLz_kg,  &
     &         p_ProcCoord, p_ProcGrid, pSI_MSG

      INTEGER  BMG_MSG_MyProc, MPI_IERR

      REAL*8  REL_RES_L2, RES_L2, RES_L2_0
      REAL*8 T, T1, T2, TOL
      REAL*8  TT1, TT2
      LOGICAL NCYCLE_FLAG, INIT_CHECK, FINE_ONLY

      logical, parameter :: debug = .false.

! ==========================================================================

!      do i=1,min(NFm,NSOm)
!      write(0,*) 'BEGIN Q,QF,SO', Q(i),QF(i),SO(i),i,NFm,NSOm
!      enddo
      

      IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
         RETURN
      END IF

! -------------------
!     Zero Times:
! -------------------

      t=dZERO
      t1=dZERO
      t2=dZERO

      IF ( debug )  write(0,*) 'myProc, start WRF_SymStd_SOLVE_boxmgdg =',myProc 
!'
      IF ( BMG_IOFLAG(iBMG3_BUG_PARAMETERS) ) THEN
         CALL BMG3_SymStd_DUMP_parms(  &
     &                    BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG  &
     &                    )
      ENDIF

! --------------------------------------------------------------------
! >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
! --------------------------------------------------------------------
! --------------------------------------------
!     MPI Parallel Data Decomposition:
! --------------------------------------------

      !
      ! Processor Grid:
      !
      ntasks_x = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      ntasks_y = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)
      NProcK = 1! BMG_MSG_iGRID(id_BMG_MSG_NProcK)

      MyProc  = BMG_MSG_iGRID(id_BMG_MSG_MyProc)
      MyProcI = BMG_MSG_iGRID(id_BMG_MSG_MyProcI)
      MyProcJ = BMG_MSG_iGRID(id_BMG_MSG_MyProcJ)
      MyProcK = BMG_MSG_iGRID(id_BMG_MSG_MyProcK)

      p_NLx       = BMG_MSG_pGRID(ip_BMG_MSG_NLx_Grid)
      p_NLy       = BMG_MSG_pGRID(ip_BMG_MSG_NLy_Grid)
      p_NLz       = BMG_MSG_pGRID(ip_BMG_MSG_NLz_Grid)
      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)
      p_ProcCoord = BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord)

! ----------------------------------------
!     MSG workspace pointers:
! ----------------------------------------

      p_NLx_kg = 1
      p_NLy_kg = p_NLx_kg + ntasks_x*NOGm
      p_NLz_kg = p_NLy_kg + ntasks_y*NOGm

      pSI_MSG  = p_NLz_kg + NProcK*NOGm  !!NB: initalization of pSI_MSG

! ----------------------------------------
!     Multigrid Cycle parameters:
! ----------------------------------------

      IFD    = BMG_iPARMS(id_BMG3_STENCIL)
      IBC    = BMG_iPARMS(id_BMG3_BC)
      ISETUP = BMG_iPARMS(id_BMG3_SETUP)

      IRELAX = BMG_iPARMS(id_BMG3_RELAX)
      IRELAX_SYM = BMG_iPARMS(id_BMG3_RELAX_SYM )

      ISEP = BMG_iPARMS(id_BMG3_SEPARABLE)

      ID = BMG_iPARMS(id_BMG3_NRELAX_DOWN)
      IU = BMG_iPARMS(id_BMG3_NRELAX_UP)
      IM = BMG_iPARMS(id_BMG3_NRELAX_FG)
      IMdum = 0

      IF (BMG_iPARMS(id_BMG3_CYCLE_CLASS).EQ.BMG_N_CYCLE) THEN
         ISTRT = - BMG_iPARMS(id_BMG3_MAX_ITERS)
      ELSEIF (BMG_iPARMS(id_BMG3_CYCLE_CLASS).EQ.BMG_FMG_CYCLE) THEN
         ISTRT = BMG_iPARMS(id_BMG3_MAX_ITERS)
      ENDIF

      IVW   = BMG_iPARMS(id_BMG3_NCYCLE_TYPE)
      MCYCL = BMG_iPARMS(id_BMG3_FMG_NNCYCLE)

      ISTOP = BMG_iPARMS(id_BMG3_STOP_TEST)
      TOL   = BMG_rPARMS(id_BMG3_STOP_TOL)

      IF ( BMG_IOFLAG(iBMG3_BUG_PARAMETERS) ) THEN
         write(0,*) 'myProc, nproci,nprocj,nprock =',myProc, ntasks_x,ntasks_y, 1

      ENDIF

! ================================================

      IF (BMG_iPARMS(id_BMG2_CYCLE_CLASS).EQ.BMG_N_CYCLE) THEN
         ISTRT2 = - BMG_iPARMS(id_BMG2_MAX_ITERS)
      ELSEIF (BMG_iPARMS(id_BMG2_CYCLE_CLASS).EQ.BMG_FMG_CYCLE) THEN
         ISTRT2 = BMG_iPARMS(id_BMG2_MAX_ITERS)
      ENDIF

! ================================================

! --------------------------------------------------------------------
! >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
! --------------------------------------------------------------------

! -------------------------------------------------------------------
!     Enable MSG communications checking to make sure all initialized
! -------------------------------------------------------------------

      IF ( debug )  write(0,*) 'myProc, call MPI_INITIALIZED =',myProc

      CALL MPI_INITIALIZED(INIT_CHECK, MPI_IERR)
      IF( .NOT. INIT_CHECK ) THEN
         PRINT *,' MPI NOT INITIALIZED!'
         RETURN
      ENDIF

      IF( ISETUP.EQ.BMG_SETUP_ptrs_opers .OR.  &
     &    ISETUP.EQ.BMG_SETUP_only ) THEN

         CALL MSG_set_comm_parent(MPICOMM)
         CALL MSG_enable(BMG_MSG_MyProc,NProc)
         BMG_MSG_MyProc = MSG_myproc()
         NProc  = MSG_nproc()

      ELSE

         IF( MPICOMM .EQ. MSG_COMM_PARENT) THEN
            BMG_MSG_MyProc = MSG_myproc()
            NProc  = MSG_nproc()
         ELSE
            IF( MyProc.EQ.1 ) PRINT *,' MSG NOT INITIALIZED!'
         ENDIF

      ENDIF

! ---------------------------------
!     Initialize timers
! ---------------------------------

      T  = dZERO
      T1 = dZERO
      T2 = dZERO

! ----------------------------------
!     Sanity checks!
! ----------------------------------

      !!
      !!   We need to add way more some sanity checks here  !!!
      !!   Maybe a subroutine to check all parameter values??
      !!

      !
      ! Processor indexing
      !
      IF ( BMG_MSG_MyProc.NE.BMG_MSG_iGRID(id_BMG_MSG_MyProc) ) THEN
         WRITE(*,*) '**** ERROR: BMG3_SymStd_SOLVE_boxmg.f ... '
         WRITE(*,*) ' MyProc =         ', MyProc
         WRITE(*,*) ' BMG_MSG_MyProc = ', BMG_MSG_MyProc
         STOP
      ENDIF


! --------------------------------------------------------------------
! >>>>>>>>>>>>>>>>>>>> BEGIN:  POINTER SETUP <<<<<<<<<<<<<<<<<<<<<<<<<
! --------------------------------------------------------------------

      IF ( ISETUP.EQ.BMG_SETUP_ptrs_opers .OR.  &
     &     ISETUP.EQ.BMG_SETUP_only .OR.  &
     &     ISETUP.EQ.BMG_SETUP_ptrs ) THEN

         TT1 = MPI_Wtime()

      IF ( debug )  write(0,*) 'myProc, call BMG3_SymStd_SETUP_PtrGrid =',myProc
!'
         CALL BMG3_SymStd_SETUP_PtrGrid(  &
     &             NLx, NLy, NLz, NGx, NGy, NGz,  &
     &             iGs, jGs, kGs,  &
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,  &
     &             MyProc, NProc, ntasks_x, ntasks_y, NProcK,  &
     &             MPICOMM, BMG_iWORK_PL, BMG_iWORK_CS,  &
     &             NOGm, NFm, NCm, NSOm, NSORm, NCIm,  &
     &             NBMG_iWORK_PLm, NBMG_rWORK_PLm,  &
     &             NBMG_iWORK_CSm, NBMG_rWORK_CSm, NMSGim, NMSGrm,  &
     &             NOG, NF, NC, NSO, NSOR, NCI, IGRD,  &
     &             NBMG_iWORK_PL, NBMG_rWORK_PL,  &
     &             NBMG_iWORK_CS, NBMG_rWORK_CS, NMSGi, NMSGr &
     &             )
         IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF


         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_SETUP_PTR_GRID) =  &
     &        BMG_rPARMS(id_BMG3_TIME_SETUP_PTR_GRID) + TT2 - TT1

         ! Output
         IF( BMG_IOFLAG(iBMG3_OUT_WSPACE_SIZE) ) THEN
            WRITE (*,260)  &
     &            'Storage for a vector on all grids = ', NF
         ENDIF

         ! Output
         IF( BMG_IOFLAG(iBMG3_OUT_WSPACE_SIZE) ) THEN
            WRITE (*,270)  &
     &            'Storage for a vector on all coarse grids = ', NC
         ENDIF

         !
         ! Setup of the MSG arrays (pSI_MSG initialized above).
         !
         TT1 = MPI_Wtime()

      IF ( debug )  write(0,*) 'myProc, call BMG3_SymStd_SETUP_MSG =',myProc

         CALL BMG3_SymStd_SETUP_MSG(  &
     &             pMSG, pMSGSO, iWorkMSG, NMSGi, pSI_MSG,  &
     &             IGRD, NOG, NOGm, NProc, MyProc,  &
     &             iWorkMSG(p_NLx_kg), iWorkMSG(p_NLy_kg),  &
     &             iWorkMSG(p_NLz_kg),  &
     &             BMG_MSG_iGRID(p_NLx), BMG_MSG_iGRID(p_NLy),  &
     &             BMG_MSG_iGRID(p_NLz),  &
     &             BMG_MSG_iGRID(p_ProcGrid),  &
     &             ntasks_x, ntasks_y, NProcK,  &
     &             MPICOMM  &
     &             )

         TT2 = MPI_Wtime()

         BMG_rPARMS(id_BMG3_TIME_SETUP_MSG) =  &
     &        BMG_rPARMS(id_BMG3_TIME_SETUP_MSG) + TT2 - TT1

      ELSE

         !
         !  Maximum dimensions are the actual dimensions.
         !

         NF   = NFm
         NC   = NCm
         NSO  = NSOm
         NSOR = NSORm
         NCI  = NCIm

         NBMG_iWORK_CS = NBMG_iWORK_CSm
         NBMG_rWORK_CS = NBMG_rWORK_CSm

         NMSGi = NMSGim
         NMSGr = NMSGrm

      ENDIF

      !
      ! Create local pointers for LU solve
      !
      IF ( BMG_iPARMS(id_BMG3_CG_SOLVER).EQ. BMG_CG_SOLVE_LU ) THEN

         p_ABD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CSO)
         p_BBD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CU)

         NCBW  = BMG_iWORK_CS(id_BMG_iWORK_CS_NCBW)
         NCU   = BMG_iWORK_CS(id_BMG_iWORK_CS_NCU)

      ENDIF

      !
      !  No components just setup pointers, thanks!
      !
      IF ( ISETUP.EQ.BMG_SETUP_ptrs ) RETURN

! --------------------------------------------------------------------
! >>>>>>>>>>>>>>>>>>> END:  POINTER SETUP <<<<<<<<<<<<<<<<<<<<<<<<<<<<
! --------------------------------------------------------------------

! --------------------------------------------------------------------
! >>>>>>>>>>>>>>>>>>> BEGIN:  COMPONENT SETUP <<<<<<<<<<<<<<<<<<<<<<<<
! --------------------------------------------------------------------


 89   CONTINUE  ! >>>>>>>>> BEGIN: Component Setup

      IF ( ISETUP.EQ.BMG_SETUP_ptrs_opers .OR.  &
     &     ISETUP.EQ.BMG_SETUP_opers .OR.  &
     &     ISETUP.EQ.BMG_SETUP_only ) THEN

         !  Start the timer
         t1 = MPI_Wtime()

         !
         ! Construct coarse-grid and interpolation operators.  In
         ! addition, setup any necessary components for relaxation.
         !
         TT1 = MPI_Wtime()

      IF ( debug )  write(0,*) 'myProc, call BMG3_SymStd_SETUP_parts =',myProc

         CALL BMG3_SymStd_SETUP_parts(  &
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,  &
     &             SO, NSOm, SOR, NSORm, CI, NCIm, IGRD, NOGm, NOG,  &
     &             BMG_iWORK_PL, NBMG_iWORK_PL,  &
     &             BMG_rWORK_PL, NBMG_rWORK_PL,  &
     &             BMG_iWORK_CS, NBMG_iWORK_CS,  &
     &             BMG_rWORK_CS, NBMG_rWORK_CS,  &
     &             iWorkMSG, NMSGi, pSI_MSG, pMSG, pMSGSO,  &
     &             MSG_Buffer, NMSGr,  &
     &             MyProcI, MyProcJ, MyProcK,  &
     &             ntasks_x, ntasks_y, NProcK, NProc,  &
     &             BMG_MSG_iGRID(p_ProcGrid),  &
     &             BMG_MSG_iGRID(p_ProcCoord),  &
     &             iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),  &
     &             iWorkMSG(p_NLx_kg), iWorkMSG(p_NLy_kg),  &
     &             iWorkMSG(p_NLz_kg),  &
     &             MPICOMM  &
     &             )

      IF ( debug ) THEN
        write(0,*) 'myProc, done BMG3_SymStd_SETUP_parts =',myProc
        write(0,*) 'ISETUP = ',ISETUP,BMG_SETUP_only
        write(0,*) 'BMG_iPARMS(id_BMG3_Err_Code) = ',BMG_iPARMS(id_BMG3_Err_Code)
        write(0,*) 'BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)',BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)
      ENDIF
      IF (debug) CALL MPI_Barrier(mpicomm, i)

         IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF

         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_SETUP_PARTS) =  &
     &        BMG_rPARMS(id_BMG3_TIME_SETUP_PARTS) + TT2 - TT1

         ! Compute the setup time
         t2 = MPI_Wtime()
         t=t+t2-t1
         BMG_rPARMS(id_BMG3_TIME_SETUP_TOTAL) =  &
     &        BMG_rPARMS(id_BMG3_TIME_SETUP_TOTAL) + t2 - t1


         ! Output the setup time
         IF( BMG_IOFLAG(iBMG3_OUT_TIME_SETUP) )  THEN
            WRITE (*,220) '(3D) SETUP TIME =', T
         ENDIF

      ENDIF

      IF ( debug )  write(0,*) 'myProc, check BMG_SETUP_only =',myProc

      !
      !  No solve just the setup, thanks!
      !
      IF (ISETUP.EQ.BMG_SETUP_only) RETURN

! --------------------------------------------------------------------
! >>>>>>>>>>>>>>>>>>>> END:  COMPONENT SETUP <<<<<<<<<<<<<<<<<<<<<<<<<
! --------------------------------------------------------------------

! --------------------------------------------------------------------
! >>>>>>>>>>>>>>>>>> BEGIN: MULTIGRID CYCLING <<<<<<<<<<<<<<<<<<<<<<<
! --------------------------------------------------------------------

 95   CONTINUE

      ! Start the timer
      t1 = MPI_Wtime()

! ----------------------------------------------------------------
!     Initial residual:
! ----------------------------------------------------------------

      !
      ! Get pointers for the finest grid (GRID LEVEL = NOG)
      !
      IF (debug) CALL MPI_Barrier(mpicomm, i)
      IF ( debug )  write(0,*) 'myProc, call BMG3_SymStd_GET_pointers1 =',myProc

      CALL BMG3_SymStd_GET_pointers(  &
     &          NOG, IGRD, NOGm,  &
     &          p_U, p_SO, p_SOR, p_CI,  &
     &          IIL, JJL, KKL, IIG, JJG, KKG  &
     &          )

!      FIXME: The arguments passed in should be correct, but this
!             is a reasonable place to double check
!
!      iGs = IGRD(NOG,idL_BMG_Icoord)
!      jGs = IGRD(NOG,idL_BMG_Jcoord)
!      kGs = IGRD(NOG,idL_BMG_Kcoord)

      !
      ! Set NStncl for the finest grid
      !
      IF ( IFD.NE.1 ) THEN
         NStncl=14
      ELSE
         NStncl=4
      ENDIF

!      do i=1,p_U
!       write(*,*) 'BEFORE Q,QF,RES', Q(i), QF(i), RES(i)
!      enddo
!
!      do i=1,p_SO
!       write(*,*) 'BEFORE SO',SO(i)
!      enddo
      !
      ! if requested, sync the initial guess
      !
      IF ( BMG_iPARMS(id_BMG3_SYNC_INITIAL_GUESS)  &
     &     .eq. BMG_SYNC_INITIAL_GUESS ) THEN

      IF (debug) CALL MPI_Barrier(mpicomm, i)
      IF ( debug )  write(0,*) 'myProc, call BMG3_SymStd_UTILS_update_ghosts =',myProc
! '
         CALL BMG3_SymStd_UTILS_update_ghosts(  &
     &             NOG, Q(p_U), IIL, JJL, KKL, iWorkMSG, pMSG,  &
     &             MSG_Buffer, MPICOMM &
     &             )
!      do i=1,p_U
!       write(*,*) 'Q', Q(i)
!      enddo

      END IF

      IF (debug) CALL MPI_Barrier(mpicomm, i)
      IF ( debug )  write(0,*) 'myProc, call BMG3_SymStd_residual =',myProc

      CALL BMG3_SymStd_residual(  &
     &          NOG, NOGm, IFD, &
     &          Q(p_U), QF(p_U), SO(p_SO), RES(p_U),  &
     &          IIL, JJL, KKL, NStncl,  &
     &          iWorkMSG, NMSGi, pMSG,  &
     &          MSG_Buffer, NMSGr, MPICOMM  &
     &          )

      IF ( debug )  write(0,*) 'myProc, call BMG3_SymStd_UTILS_norm_l2 =',myProc
! '
      CALL BMG3_SymStd_UTILS_norm_l2(  &
     &           RES(p_U), IIL, JJL, KKL, RES_L2_0, MPICOMM  &
     &           )

!      write(*,*) 'NOG, NOGm, IFD', NOG, NOGm, IFD

!      do i=1,p_U
!       write(*,*) 'Q,QF,RES,IIL,JJL,KKL', Q(i), QF(i), RES(i),IIL, JJL, KKL
!      enddo

!     do i=1,p_SO
!       write(*,*) 'SO',SO(i)
!      enddo

      !
      !  Output
      !
      IF ( BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
         IF ( IM == 0 ) THEN
           WRITE (*,405) '*** THE INITIAL RESIDUAL DG (l2-NORM) = ', RES_L2_0
         ELSE
           WRITE (*,405) '*** THE INIT RES DG-FG = ', RES_L2_0
         ENDIF
      ENDIF

      IF ( RES_L2_0 .EQ. dZERO ) THEN
         IF ( BMG_IOFLAG(iBMG3_WARN_ZERO_RESIDUAL) ) THEN
            WRITE(*,505) '*** WARNING: bmg3d.f '
            WRITE(*,505) '    Zero initial residual !!!! '
            WRITE(*,505) '    D''ooh! Initial guess is the solution!'
         ENDIF
         RETURN
      ENDIF

! ----------------------------------------------------------------
!     Direct solve in the degenerate case of one grid!
! ----------------------------------------------------------------

      IF ( NOG.EQ.1 ) THEN

         !
         ! Solve on the coarsest grid
         !
         IF ( MyProc.EQ.1 ) THEN
            WRITE(*,*) '*** Performing direct solve on the fine grid!'
         ENDIF
         !

         IF ( BMG_iPARMS(id_BMG3_CG_SOLVER).EQ. BMG_CG_SOLVE_LU ) THEN
            !
            TT1 = MPI_Wtime()
            !
            CALL BMG3_SymStd_SOLVE_cg_LU(  &
     &                q(p_U), qf(p_U), IIL, JJL, KKL,  &
     &                BMG_rWORK_CS(p_ABD), BMG_rWORK_CS(p_BBD),  &
     &                NCBW, NCU, NOGm,  &
     &                ntasks_x, ntasks_y, NProcK, NProc, MyProc,  &
     &                BMG_MSG_iGRID(p_ProcGrid),  &
     &                BMG_MSG_iGRID(p_ProcCoord),  &
     &                iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),  &
     &                MSG_Buffer, NMSGr, MPICOMM, BMG_iPARMS  &
     &                )
            !
            TT2 = MPI_Wtime()
            !
         ELSE
            !
            TT1 = MPI_Wtime()
            !

            CALL BMG3_SymStd_SOLVE_cg_boxmg(  &
     &                Q(p_U), QF(p_U), NOG,  &
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,  &
     &                IIL, JJL, KKL, IIG, JJG, KKG,  &
     &                iGs, jGs, kGs,  &
     &                BMG_iWORK_CS, NBMG_iWORK_CS,  &
     &                BMG_rWORK_CS, NBMG_rWORK_CS,  &
     &                ntasks_x, ntasks_y, NProcK, NProc, MyProc,  &
     &                BMG_MSG_iGRID(p_ProcGrid),  &
     &                BMG_MSG_iGRID(p_ProcCoord),  &
     &                iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),  &
     &                MSG_Buffer, NMSGr, MPICOMM  &
     &                )
            !
            TT2 = MPI_Wtime()
            !
         ENDIF

         IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF

         BMG_rPARMS(id_BMG3_TIME_SOLVE_CG)  &
     &       = BMG_rPARMS(id_BMG3_TIME_SOLVE_CG) + TT2 - TT1

         !
         ! Compute the final residual
         !
         CALL BMG3_SymStd_residual(  &
     &             NOG, NOG, IFD,  &
     &             Q(p_U), QF(p_U), SO(p_SO), RES(p_U),  &
     &             IIL, JJL, KKL, NStncl,  &
     &             iWorkMSG, NMSGi, pMSG,  &
     &             MSG_Buffer, NMSGr, MPICOMM  &
     &             )
         !
         CALL BMG3_SymStd_UTILS_norm_l2(  &
     &             RES(p_U), IIL, JJL, KKL, RES_L2, MPICOMM  &
     &             )
         !
         IF( BMG_IOFLAG(iBMG3_BUG_RES_CG_SOLVE) .AND. MyProc.EQ.1 ) THEN
            WRITE (*,405) '*** THE FINAL RESIDUAL (l2-NORM) = ', RES_L2
         ENDIF
         !
         TOL = RES_L2
         ! Jump to final I/O
         GOTO 150
         !
      ENDIF



! --------------------------------
!     Cycle Preliminaries:
! --------------------------------

      KF = NOG                  ! finest grid index
      KC = 1                    ! coarsest grid index
      KFMG = 1                  ! FMG current finest grid
      MCYC = ABS(ISTRT)         ! Maximum number of multigrid cycles
      NCYC = 1                  ! multigrid cycle counter
      NCYCLE_FLAG = .TRUE.      ! default to do n-cycles
      FINE_ONLY = .FALSE.       ! run relaxation on finest grid only

      IF ( IM .gt. 0 ) FINE_ONLY = .TRUE.

! ----------------------------------------------------------------
!     Perform multigrid F-cycle:
!     (if necessary, then continue with n-cycles)
! ----------------------------------------------------------------

      IF ( ISTRT.GT.0 .and. .not. FINE_ONLY ) THEN

         ! Coarsen the right hand side
         DO k=KF-1, 1, -1

            CALL BMG3_SymStd_GET_pointers( k+1, IGRD, NOGm,  &
     &                p_U, p_SO, p_SOR, p_CI,  &
     &                IILF, JJLF, KKLF, IIGF, JJGF, KKGF  &
     &                )
            iGs_f = IGRD(k+1,idL_BMG_Icoord)
            jGs_f = IGRD(k+1,idL_BMG_Jcoord)
            kGs_f = IGRD(k+1,idL_BMG_Kcoord)

            CALL BMG3_SymStd_GET_pointers( k, IGRD, NOGm,  &
     &                p_UC, p_SOC, p_SORC, p_CIC,  &
     &                IILC, JJLC, KKLC, IIGC, JJGC, KKGC  &
     &                )
            iGs_c = IGRD(k,idL_BMG_Icoord)
            jGs_c = IGRD(k,idL_BMG_Jcoord)
            kGs_c = IGRD(k,idL_BMG_Kcoord)

            ! restrict the right hand side

            TT2 = MPI_Wtime()

            CALL BMG3_SymStd_restrict(  &
     &                K+1, K,  &
     &                QF(p_U), QF(p_UC), CI(p_CIC),  &
     &                IILF, JJLF, KKLF, IILC, JJLC, KKLC,  &
     &                iGs_f, jGs_f, kGs_f  &
     &                )

            TT1 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_restrict) =  &
     &           BMG_rPARMS(id_BMG3_TIME_restrict) + TT2 - TT1

            ! Copy right hand side into RES(i,j,k) for interpolation
            DO i=0, IILF*JJLF*KKLF-1
               RES(p_U+i)=QF(p_U+i)
            END DO
         END DO


         ! FMG current finest grid index is set to the coarsest grid
         KFMG=KC
         ! FMG: solve exactly on the coarsest grid
         CALL BMG3_SymStd_GET_pointers(  &
     &             KFMG, IGRD, NOGm,  &
     &             p_U, p_SO, p_SOR, p_CI,  &
     &             IIL, JJL, KKL, IIG, JJG, KKG  &
     &             )
         iGs_c = IGRD(KFMG,idL_BMG_Icoord)
         jGs_c = IGRD(KFMG,idL_BMG_Jcoord)
         kGs_c = IGRD(KFMG,idL_BMG_Kcoord)

         IF ( BMG_iPARMS(id_BMG3_CG_SOLVER).EQ.BMG_CG_SOLVE_LU ) THEN
            !
            TT1 = MPI_Wtime()
            !
            CALL BMG3_SymStd_SOLVE_cg_LU(  &
     &                q(p_U), qf(p_U), IIL, JJL, KKL,  &
     &                BMG_rWORK_CS(p_ABD), BMG_rWORK_CS(p_BBD),  &
     &                NCBW, NCU, NOGm,  &
     &                ntasks_x, ntasks_y, NProcK, NProc, MyProc,  &
     &                BMG_MSG_iGRID(p_ProcGrid),  &
     &                BMG_MSG_iGRID(p_ProcCoord),  &
     &                iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),  &
     &                MSG_Buffer, NMSGr, MPICOMM, BMG_iPARMS  &
     &                )
            !
            TT2 = MPI_Wtime()
            !
         ELSE
            !
            TT1 = MPI_Wtime()
            !

            CALL BMG3_SymStd_SOLVE_cg_boxmg(  &
     &                Q(p_U), QF(p_U), NOG,  &
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,  &
     &                IIL, JJL, KKL, IIG, JJG, KKG,  &
     &                iGs_c, jGs_c, kGs_c,  &
     &                BMG_iWORK_CS, NBMG_iWORK_CS,  &
     &                BMG_rWORK_CS, NBMG_rWORK_CS,  &
     &                ntasks_x, ntasks_y, NProcK, NProc, MyProc,  &
     &                BMG_MSG_iGRID(p_ProcGrid),  &
     &                BMG_MSG_iGRID(p_ProcCoord),  &
     &                iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),  &
     &                MSG_Buffer, NMSGr, MPICOMM  &
     &                )
            !
            TT2 = MPI_Wtime()
            !
         ENDIF

         IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF

         BMG_rPARMS(id_BMG3_TIME_SOLVE_CG)  &
     &       = BMG_rPARMS(id_BMG3_TIME_SOLVE_CG) + TT2 - TT1

         RES_L2=dZERO
         IF ( BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
            WRITE (*,230) KFMG, RES_L2
         ENDIF

 100     CONTINUE    ! >>>>>>>> LOOP BOUNDARY: multigrid F-cycle


            KFMG=KFMG+1
            !
            ! FMG: Interpolate
            !
            CALL BMG3_SymStd_GET_pointers(  &
     &                KFMG, IGRD, NOGm,  &
     &                p_U, p_SO, p_SOR, p_CI, IIL, JJL, KKL,  &
     &                IIG, JJG, KKG  &
     &                )
            iGs_f = IGRD(KFMG,idL_BMG_Icoord)
            jGs_f = IGRD(KFMG,idL_BMG_Jcoord)
            kGs_f = IGRD(KFMG,idL_BMG_Kcoord)

            CALL BMG3_SymStd_GET_pointers(  &
     &                KFMG-1, IGRD, NOGm,  &
     &                p_UC, p_SOC, p_SORC, p_CIC, IILC, JJLC, KKLC,  &
     &                IIGC, JJGC, KKGC  &
     &                )
            iGs_c = IGRD(KFMG-1,idL_BMG_Icoord)
            jGs_c = IGRD(KFMG-1,idL_BMG_Jcoord)
            kGs_c = IGRD(KFMG-1,idL_BMG_Kcoord)

            !
            ! Set NStncl for grid KFMG
            !
            IF ( KFMG.LT.NOG .OR. IFD.NE.BMG_STENCIL_7pt ) THEN
               NStncl=14
            ELSE
               NStncl=4
            ENDIF
            !


            TT1 = MPI_Wtime()

            CALL BMG3_SymStd_interp_add(  &
     &                KFMG-1, KFMG, NOG,  &
     &                Q(p_U) ,Q(p_UC), RES(p_U),  &
     &                SO(p_SO), NStncl, CI(p_CIC),  &
     &                IILC, JJLC, KKLC,  &
     &                IIL, JJL, KKL,  &
     &                iGs_f, jGs_f, kGs_f,  &
     &                iWorkMSG, NMSGi, pMSG,  &
     &                MSG_Buffer, NMSGr, MPICOMM  &
     &                )


            TT2 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_interp_add) =  &
     &           BMG_rPARMS(id_BMG3_TIME_interp_add) + TT2 - TT1
            !
            ! FMG; Perform an n-cycle
            !
            IF ( BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
               WRITE(*,410) '*** Performing n-cycle:  Coarsest =', KC,  &
     &                      'Finest =', KFMG
            ENDIF


            CALL WRF_SymStd_ncycledg(  &
     &                KC, KFMG, KF,  &
     &                IFD, IU, ID, IMdum, IVW, IRELAX, IRELAX_SYM, ISEP,  &
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,  &
     &                Q, QF, RES, NFm, NCm,  &
     &                SO, NSOm, SOR, NSORm, CI, NCIm, IGRD, NOGm,  &
     &                BMG_iWORK_PL, NBMG_iWORK_PL,  &
     &                BMG_rWORK_PL, NBMG_rWORK_PL,  &
     &                BMG_iWORK_CS, NBMG_iWORK_CS,  &
     &                BMG_rWORK_CS, NBMG_rWORK_CS,  &
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,  &
     &                MSG_Buffer, NMSGr,  &
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID  &
     &                ) ! ,Qmask)


            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            IF ( KFMG.EQ.KF ) THEN
               !
               !  Compute the final residual
               !
               CALL BMG3_SymStd_residual(  &
     &                   KF, KF, IFD,  &
     &                   Q(p_U), QF(p_U), SO(p_SO), RES(p_U),  &
     &                   IIL, JJL, KKL, NStncl,  &
     &                   iWorkMSG, NMSGi, pMSG,  &
     &                   MSG_Buffer, NMSGr, MPICOMM  &
     &                   )
               CALL BMG3_SymStd_UTILS_norm_l2(  &
     &                   RES(p_U), IIL, JJL, KKL, RES_L2, MPICOMM  &
     &                   )

               ! Note this is still vulnerable to very small RES_L2_0
               IF ( RES_L2_0.EQ.dZERO ) THEN
                  IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR).AND.MyProc.EQ.1)  &
     &                 THEN
                    WRITE(*,505)  &
     &                    '*** FATAL ERROR: BMG3_SymStd_SOLVE_boxmg.f '
                    WRITE(*,505)'    Initial residual has become zero!'
                  END IF

                  CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,1)
                  RETURN

               ELSE
                  REL_RES_L2 = RES_L2/RES_L2_0
               ENDIF

               !
               !  Output
               !
               IF ( ( ISTOP.EQ.BMG_STOP_REL_RES_L2 )  &
     &            .AND. BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
                  WRITE (*,400) '*** ITERATION =', NCYC,  &
     &                          '*** RELATIVE RESIDUAL = ', REL_RES_L2
               ELSE IF ( ( ISTOP.EQ.BMG_STOP_ABS_RES_L2 )  &
     &                 .AND. BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
                  WRITE (*,400) '*** ITERATION =', NCYC,  &
     &                          '*** ABSOLUTE RESIDUAL = ', RES_L2
               ENDIF
               !
               !  Check Convergence
               !
               IF ( ( ISTOP.EQ.BMG_STOP_REL_RES_L2 )  &
     &            .AND. ( REL_RES_L2.LT.TOL )     ) THEN
                  ! FMG: converged
                  !      - set TOL and skip n-cycles
                  TOL = REL_RES_L2
                  NCYCLE_FLAG=.FALSE.
               ELSE IF ( ( ISTOP.EQ.BMG_STOP_ABS_RES_L2 )  &
     &                 .AND. ( RES_L2.LT.TOL )         ) THEN
                  ! FMG: converged
                  !      - set TOL and skip n-cycles
                  TOL = RES_L2
                  NCYCLE_FLAG=.FALSE.
               ELSE IF ( NCYC.GE.MCYC ) THEN
                  ! FMG: cycle limit reached
                  !      - set TOL and skip n-cycles
                  IF ( ISTOP.EQ.BMG_STOP_REL_RES_L2 ) THEN
                     TOL = -REL_RES_L2
                  ELSE
                     TOL = -RES_L2
                  ENDIF
                  NCYCLE_FLAG=.FALSE.
               ELSE
                  ! FMG: didn't converge
                  !      - increase counter and continue with n-cycles
                  NCYC=NCYC+1
               ENDIF
               !
            ELSE
               ! FMG: Continue the F-cycle
               GOTO 100
            ENDIF

 120     CONTINUE    ! <<<<<<<< LOOP BOUNDARY: multigrid F-cycle

      ENDIF  ! ISTRT > 0

! ----------------------------------------------------------------
!     Perform IM relaxations on fine grid
! ----------------------------------------------------------------

      IF ( FINE_ONLY ) THEN
        NCYCLE_FLAG = .FALSE.

         CALL BMG3_SymStd_GET_pointers(  &
     &             KF, IGRD, NOGm, p_U, p_SO,  &
     &             p_SOR, p_CI, IIL, JJL, KKL,  &
     &             IIG, JJG, KKG  &
     &             )

            CALL WRF_SymStd_ncycledg(  &
     &                KC, KF, KF,  &
     &                IFD, IU, ID, IM, IVW, IRELAX, IRELAX_SYM, ISEP,  &
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,  &
     &                Q, QF, RES, NFm, NCm,  &
     &                SO, NSOm, SOR, NSORm, CI, NCIm, IGRD, NOGm,  &
     &                BMG_iWORK_PL, NBMG_iWORK_PL,  &
     &                BMG_rWORK_PL, NBMG_rWORK_PL,  &
     &                BMG_iWORK_CS, NBMG_iWORK_CS,  &
     &                BMG_rWORK_CS, NBMG_rWORK_CS,  &
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,  &
     &                MSG_Buffer, NMSGr,  &
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID  &
     &                ) ! ,Qmask)

            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            !
            !  Compute the final residual
            !
            CALL BMG3_SymStd_residual(  &
     &                KF, KF, IFD,  &
     &                Q(p_U), QF(p_U), SO(p_SO), RES(p_U),  &
     &                IIL, JJL, KKL, NStncl,  &
     &                iWorkMSG, NMSGi, pMSG,  &
     &                MSG_Buffer, NMSGr, MPICOMM  &
     &                )
            CALL BMG3_SymStd_UTILS_norm_l2(  &
     &                RES(p_U), IIL, JJL, KKL, RES_L2, MPICOMM  &
     &                )

            ! Note this is still vulnerable to very small RES_L2_0
            IF ( RES_L2_0.EQ.dZERO ) THEN
               WRITE(*,505) '*** FATAL ERROR: bmg3d.f '
               WRITE(*,505) '    Initial residual has become zero!!!'
               STOP
            ELSE
               REL_RES_L2 = RES_L2/RES_L2_0
            ENDIF

            !
            !  Output
            !
            IF ( ( ISTOP.EQ.BMG_STOP_REL_RES_L2 )  &
     &           .AND. BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
               WRITE (*,400) '*** FG_ONLY ITERATION =', NCYC,  &
     &                       '*** FG_ONLY RELATIVE RESIDUAL = ', REL_RES_L2
            ELSE IF ( ( ISTOP.EQ.BMG_STOP_ABS_RES_L2 )  &
     &              .AND. BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
               WRITE (*,400) '*** FG_ONLY ITERATION =', NCYC,  &
     &                      '*** FG_ONLY ABSOLUTE RESIDUAL = ', RES_L2
            ENDIF

      ENDIF

! ----------------------------------------------------------------
!     Perform multigrid n-cycle(s):
! ----------------------------------------------------------------

      IF ( NCYCLE_FLAG ) THEN

      IF ( debug )  write(0,*) 'myProc, call BMG3_SymStd_GET_pointers2 =',myProc

         CALL BMG3_SymStd_GET_pointers(  &
     &             KF, IGRD, NOGm, p_U, p_SO,  &
     &             p_SOR, p_CI, IIL, JJL, KKL,  &
     &             IIG, JJG, KKG  &
     &             )


 130     CONTINUE   ! >>>>>>>> LOOP BOUNDARY: multigrid n-cycles


      IF ( debug )  write(0,*) 'myProc, call WRF_SymStd_ncycledg =',myProc

            !
            ! Perform a multigrid n-cycle
            !
            CALL WRF_SymStd_ncycledg(  &
     &                KC, KF, KF,  &
     &                IFD, IU, ID, IMdum, IVW, IRELAX, IRELAX_SYM, ISEP,  &
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,  &
     &                Q, QF, RES, NFm, NCm,  &
     &                SO, NSOm, SOR, NSORm, CI, NCIm, IGRD, NOGm,  &
     &                BMG_iWORK_PL, NBMG_iWORK_PL,  &
     &                BMG_rWORK_PL, NBMG_rWORK_PL,  &
     &                BMG_iWORK_CS, NBMG_iWORK_CS,  &
     &                BMG_rWORK_CS, NBMG_rWORK_CS,  &
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,  &
     &                MSG_Buffer, NMSGr,  &
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID  &
     &                ) ! ,Qmask)

            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            !
            !  Compute the final residual
            !
            CALL BMG3_SymStd_residual(  &
     &                KF, KF, IFD,  &
     &                Q(p_U), QF(p_U), SO(p_SO), RES(p_U),  &
     &                IIL, JJL, KKL, NStncl,  &
     &                iWorkMSG, NMSGi, pMSG,  &
     &                MSG_Buffer, NMSGr, MPICOMM  &
     &                )
            CALL BMG3_SymStd_UTILS_norm_l2(  &
     &                RES(p_U), IIL, JJL, KKL, RES_L2, MPICOMM  &
     &                )

            ! Note this is still vulnerable to very small RES_L2_0
            IF ( RES_L2_0.EQ.dZERO ) THEN
               WRITE(*,505) '*** FATAL ERROR: bmg3d.f '
               WRITE(*,505) '    Initial residual has become zero!!!'
               STOP
            ELSE
               REL_RES_L2 = RES_L2/RES_L2_0
            ENDIF

            !
            !  Output
            !
            IF ( ( ISTOP.EQ.BMG_STOP_REL_RES_L2 )    &
     &           .AND. BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
               WRITE (*,400) '*** ITERATION =', NCYC,  &
     &                       '*** RELATIVE RESIDUAL = ', REL_RES_L2
            ELSE IF ( ( ISTOP.EQ.BMG_STOP_ABS_RES_L2 )   &
     &              .AND. BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
               WRITE (*,400) '*** ITERATION =', NCYC,  &
     &                      '*** ABSOLUTE RESIDUAL = ', RES_L2
            ENDIF

            !
            ! Check Convergence
            !
            IF ( ( ISTOP.EQ.BMG_STOP_REL_RES_L2 )  &
     &         .AND. ( REL_RES_L2.LT.TOL )     ) THEN
               !
               ! n-cycles have converged in the RELATIVE RESIDUAL
               ! - set TOL and return
               !
               TOL = REL_RES_L2
            ELSE IF ( ( ISTOP.EQ.BMG_STOP_ABS_RES_L2 )  &
     &              .AND. ( RES_L2.LT.TOL )         ) THEN
               !
               ! n-cycles have converged in the ABSOLUTE RESIDUAL
               ! - set TOL and return
               !
               TOL = RES_L2
            ELSE IF ( NCYC.GE.MCYC ) THEN
               !
               ! n-cycles failed to converge, cycle limit reached
               ! - set TOL and return
               !
               IF ( ISTOP.EQ.BMG_STOP_REL_RES_L2 ) THEN
                  TOL = -REL_RES_L2
               ELSE
                  TOL = -RES_L2
               ENDIF
            ELSE
               !
               ! increase cycle count and continue n-cycling
               !
               NCYC=NCYC+1
               GOTO 130
            ENDIF

 140     CONTINUE    ! <<<<<<<< LOOP BOUNDARY: multigrid n-cycles

      ENDIF

 150  CONTINUE       ! <<<<<<<< FINAL I/O
      !
      !
      BMG_rPARMS(id_BMG3_STOP_TOL) = TOL
      !
      !  Compute the solve time and update the total time
      !
      t2 = MPI_Wtime()
      T=T+T2-T1
      BMG_rPARMS(id_BMG3_TIME_SOLVE_total) =  &
     &     BMG_rPARMS(id_BMG3_TIME_SOLVE_total) + t2 - t1


      ! Output the multigrid cycling time
      IF ( BMG_IOFLAG(iBMG3_OUT_TIME_CYCLING) ) THEN
         WRITE(*,240) '(3D) MULTIGRID CYCLING TIME =', T2-T1

      ENDIF

      IF ( BMG_IOFLAG(iBMG3_OUT_TIME_TOTAL) ) THEN
         WRITE(*,240) '(3D) TOTAL TIME = ', T
      ENDIF

! -------------------------------------
!     Output:
! -------------------------------------

! ==========================================================================

 220  FORMAT (/,1X,A,1X,F9.3,/)
 230  FORMAT ('(3D)LEVEL',I2,' RESIDUAL NORM= ',1P,E10.3)
 240  FORMAT (/,1X,A,1X,F12.3,/)
 260  FORMAT (/,/,2X,A,1X,I5)
 270  FORMAT (2X,A,1X,I5,/)

 400  FORMAT (1X,A,1X,I2,4X,A,1X,1P,E16.9)
 405  FORMAT (/,1X,A,1X,1P,E16.9,/)
 410  FORMAT (1X,A,1X,I2,4X,A,1X,I2)

 500  FORMAT (/,2X,A)
 505  FORMAT (5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)

! ===========================================

      RETURN
      END SUBROUTINE WRF_SymStd_SOLVE_boxmgdg

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE WRF_SymStd_ncycledg(  &
     &                KC, KCF, NOG,  &
     &                IFD, IU, ID, IM, IVW, IRELAX, IRELAX_SYM, ISEP,  &
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,  &
     &                Q, QF, RES, NF, NC,  &
     &                SO, NSO, SOR, NSOR, CI, NCI, IGRD, NOGm,  &
     &                BMG_iWORK_PL, NBMG_iWORK_PL,  &
     &                BMG_rWORK_PL, NBMG_rWORK_PL,  &
     &                BMG_iWORK_CS, NBMG_iWORK_CS,  &
     &                BMG_rWORK_CS, NBMG_rWORK_CS,  &
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,  &
     &                MSG_Buffer, NMSGr,  &
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID  &
     &                ) !,Qmask)

! ==========================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!   BMG3_SymStd_ncycle performs a single multigrid n-cycle.
!
! =======================================================================
! $license_flag$
! =======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   INPUT/OUTPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   OUTPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   LOCAL:
!  --------------------
!
!
! ==========================================================================

      IMPLICIT NONE

! -----------------------------
!     Includes
!
      INCLUDE 'mpif.h'

#include    "BMG_constants.h"
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_parameters.h'

! ----------------------------
!     Argument Declarations
!
      INTEGER  KC, KCF, NOG, &
     &         IFD, IU, ID, IM, IVW, IRELAX, IRELAX_SYM, ISEP, &
     &         NF, NC, NSO, NSOR, NCI, NOGm, NMSGi, NMSGr, &
     &         NBMG_MSG_iGRID, NBMG_iWORK_CS, NBMG_iWORK_PL, &
     &         NBMG_rWORK_CS, NBMG_rWORK_PL

      INTEGER  BMG_iPARMS(NBMG_iPARMS), &
     &         BMG_iWORK_CS(NBMG_iWORK_CS), &
     &         BMG_iWORK_PL(NBMG_iWORK_PL), &
     &         iWorkMSG(NMSGi), pMSG(NBMG_pMSG,NOGm), &
     &         pMSGSO(NBMG_pMSG,NOGm), IGRD(NOGm,NBMG_pIGRD), &
     &         BMG_MSG_iGRID(NBMG_MSG_iGRID), &
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS), &
     &         BMG_rWORK_CS(NBMG_rWORK_CS), &
     &         BMG_rWORK_PL(NBMG_rWORK_PL), &
     &         CI(NCI), Q(NF), QF(NF), RES(NF), &
     &         SO(NSO), SOR(NSOR), MSG_Buffer(NMSGr)
!      REAL*RKIND :: Qmask(NF)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

! ----------------------------
!     Local Declarations
!
      INTEGER  MyProc, MyProcI, MyProcJ, MyProcK, NProc, ntasks_x, ntasks_y, &
     &         NProcK, MPICOMM, MPI_IERROR, NStncl

      INTEGER  iGs, jGs, kGs, NLx, NLy, NLz, NGx, NGy, NGz

      INTEGER  p_ProcGrid, p_ProcCoord

      INTEGER  i, j, k, p_ABD, p_BBD, p_CI, p_SO, p_SOR, p_U, &
     &         NCBW, NCU
      REAL*8   RES_L2
      REAL*8  TT1, TT2

! ==========================================================================

! --------------------------------------------------------------------
! >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
! --------------------------------------------------------------------
! --------------------------------------------
!     MPI Parallel Data Decomposition:
! --------------------------------------------
      !
      ! MPI Communicators
      !
      MPICOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM)
      !
      ! Processor Grid:
      !
      NProc  = BMG_MSG_iGRID(id_BMG_MSG_NProc)
      ntasks_x = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      ntasks_y = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)
      NProcK = 1! BMG_MSG_iGRID(id_BMG_MSG_NProcK)

      MyProc  = BMG_MSG_iGRID(id_BMG_MSG_MyProc)
      MyProcI = BMG_MSG_iGRID(id_BMG_MSG_MyProcI)
      MyProcJ = BMG_MSG_iGRID(id_BMG_MSG_MyProcJ)
      MyProcK = BMG_MSG_iGRID(id_BMG_MSG_MyProcK)

      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)
      p_ProcCoord = BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord)

! --------------------------------------------------------------------
! >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
! --------------------------------------------------------------------

      k = KCF   ! set the current grid index

      IF (KCF.EQ.KC) THEN
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN
            WRITE(*,*) '*** FATAL ERROR: BMG3_SymStd_ncycle.f ****'
            WRITE(*,*) 'The current finest grid = ', KCF
            WRITE(*,*) 'The coarest grid        = ', KC
            WRITE(*,*) 'ERROR; minimum number of grids is 2 !'
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,18)
         RETURN

      ENDIF

      ! set the n-cycle counter to zero on all grids
      DO i = KC, KCF
         IGRD(i,idL_BMG_IVW)=iZERO
      END DO


! -----------------------------------------------------------
!     If IM > 0 then run fine grid relax only and return.
! -----------------------------------------------------------
      IF ( IM .gt. 0 ) THEN

            CALL WRF_SymStd_updowndg( &
     &                KCF, KCF, NOG, BMG_UP, &
     &                IFD, IU, ID, IM, IRELAX, IRELAX_SYM, ISEP, &
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, &
     &                Q, QF, RES, NF, NC, &
     &                SO, NSO, SOR, NSOR, CI, NCI, &
     &                IGRD, NOGm, &
     &                BMG_iWORK_PL, NBMG_iWORK_PL, &
     &                BMG_rWORK_PL, NBMG_rWORK_PL, &
     &                iWorkMSG, NMSGi, pMSG, pMSGSO, &
     &                MSG_Buffer, NMSGr, &
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID &
     &                ) ! ,Qmask)
            !
            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF




         RETURN
      ENDIF


! -----------------------------------------------------------
!     Begin the n-cycle
!     (this is also the n-cycle recursive loop boundary)
! -----------------------------------------------------------

 140  CONTINUE

         IF ( K.EQ.KC .AND. IGRD(k,idL_BMG_IVW).EQ.iZERO ) THEN
            !
            ! Solve on the coarsest grid
            !
            CALL BMG3_SymStd_GET_pointers( &
     &                KC, IGRD, NOGm, &
     &                p_U, p_SO, p_SOR, p_CI, &
     &                NLx, NLy, NLz, NGx, NGy, NGz &
     &                )
            !
            iGs = IGRD(KC,idL_BMG_Icoord)
            jGs = IGRD(KC,idL_BMG_Jcoord)
            kGs = IGRD(KC,idL_BMG_Kcoord)
            !
            IF ( BMG_iPARMS(id_BMG3_CG_SOLVER).EQ.BMG_CG_SOLVE_LU ) THEN
               !
               p_ABD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CSO)
               p_BBD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CU)
               !
               NCBW  = BMG_iWORK_CS(id_BMG_iWORK_CS_NCBW)
               NCU   = BMG_iWORK_CS(id_BMG_iWORK_CS_NCU)
               !
               TT1 = MPI_Wtime()
               !
               CALL BMG3_SymStd_SOLVE_cg_LU( &
     &                   q(p_U), qf(p_U), NLx, NLy, NLz, &
     &                   BMG_rWORK_CS(p_ABD), BMG_rWORK_CS(p_BBD), &
     &                   NCBW, NCU, NOGm, &
     &                   ntasks_x, ntasks_y, NProcK, NProc, MyProc, &
     &                   BMG_MSG_iGRID(p_ProcGrid), &
     &                   BMG_MSG_iGRID(p_ProcCoord), &
     &                   iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)), &
     &                   MSG_Buffer, NMSGr, MPICOMM, BMG_iPARMS &
     &                   )
               !
               TT2 = MPI_Wtime()
               !
            ELSE
               !
               TT1 = MPI_Wtime()
               !
               CALL BMG3_SymStd_SOLVE_cg_boxmg( &
     &                   Q(p_U), QF(p_U), NOG, &
     &                   BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, &
     &                   NLx, NLy, NLz, NGx, NGy, NGz, &
     &                   iGs, jGs, kGs, &
     &                   BMG_iWORK_CS, NBMG_iWORK_CS, &
     &                   BMG_rWORK_CS, NBMG_rWORK_CS, &
     &                   ntasks_x, ntasks_y, NProcK, NProc, MyProc, &
     &                   BMG_MSG_iGRID(p_ProcGrid), &
     &                   BMG_MSG_iGRID(p_ProcCoord), &
     &                   iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)), &
     &                   MSG_Buffer, NMSGr, MPICOMM &
     &                   )

               !
               TT2 = MPI_Wtime()
               !
            ENDIF

            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            BMG_rPARMS(id_BMG3_TIME_SOLVE_CG) &
     &          = BMG_rPARMS(id_BMG3_TIME_SOLVE_CG) + TT2 - TT1

            !
            IGRD(k,idL_BMG_IVW) = IVW
            IF( BMG_IOFLAG(iBMG3_BUG_RES_CG_SOLVE) &
     &          .AND. BMG_iPARMS(id_BMG3_CG_SOLVER).EQ.BMG_CG_SOLVE_LU &
     &        ) THEN
               NStncl = 14   ! always the same on coarse grid if NOG > 1
               CALL BMG3_SymStd_residual( &
     &                   KC, NOG, IFD, &
     &                   q(p_U), qf(p_U), so(p_SO), RES(p_U), &
     &                   NLx, NLy, NLz, NStncl, &
     &                   iWorkMSG, NMSGi, pMSG, &
     &                   MSG_Buffer, NMSGr, MPICOMM &
     &                   )
               CALL BMG3_SymStd_UTILS_norm_l2( &
     &                   RES(p_U), NLx, NLy, NLz, RES_L2, MPICOMM &
     &                   )
               !
               IF ( MyProc.EQ.1 ) THEN
                  WRITE (*,100) k, RES_L2
               ENDIF
               !
            ENDIF
            !
            ! continue the n-cycle
            !
            GOTO 140
            !
         ELSE IF ( IGRD(k,idL_BMG_IVW).LT.IVW ) THEN
            !
            ! Move down
            !
            CALL WRF_SymStd_updowndg( &
     &                k, KCF, NOG, BMG_DOWN, &
     &                IFD, IU, ID, IM, IRELAX, IRELAX_SYM, ISEP, &
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, &
     &                Q, QF, RES, NF, NC, &
     &                SO, NSO, SOR, NSOR, CI, NCI, &
     &                IGRD, NOGm, &
     &                BMG_iWORK_PL, NBMG_iWORK_PL, &
     &                BMG_rWORK_PL, NBMG_rWORK_PL, &
     &                iWorkMSG, NMSGi, pMSG, pMSGSO, &
     &                MSG_Buffer, NMSGr, &
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID &
     &                ) ! ,Qmask)
            !
            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            IGRD(k,idL_BMG_IVW) = IGRD(k,idL_BMG_IVW) + 1
            k=k-1
            !
            ! Continue the n-cycle
            !
            GOTO 140
            !
         ELSE IF ( IGRD(k,idL_BMG_IVW).EQ.IVW ) THEN
            !
            ! Move up
            !
            IGRD(k,idL_BMG_IVW) = iZERO   ! reset n-cycle counter to zero
            k=k+1
            !
            CALL WRF_SymStd_updowndg( &
     &                k, KCF, NOG, BMG_UP, &
     &                IFD, IU, ID, IM, IRELAX, IRELAX_SYM, ISEP, &
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, &
     &                Q, QF, RES, NF, NC, &
     &                SO, NSO, SOR, NSOR, CI, NCI, &
     &                IGRD, NOGm, &
     &                BMG_iWORK_PL, NBMG_iWORK_PL, &
     &                BMG_rWORK_PL, NBMG_rWORK_PL, &
     &                iWorkMSG, NMSGi, pMSG, pMSGSO, &
     &                MSG_Buffer, NMSGr, &
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID &
     &                ) ! ,Qmask)
            !
            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            !
            ! Are we on the finest grid?
            !
            IF (K.EQ.KCF) THEN
               ! n-cycle is done
               RETURN
            ELSE
               ! continue the n-cycle
               GOTO 140
            ENDIF
            !
         ENDIF

! ==========================================================================

 100     FORMAT (1X,'(3D)LEVEL ',I2,' RESIDUAL NORM = ',1P,E12.5)

! ===========================================

         RETURN
         END SUBROUTINE WRF_SymStd_ncycledg

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE WRF_SymStd_updowndg ( &
     &                K, KF, NOG, UPDOWN, &
     &                IFD, IU, ID, IM, IRELAX, IRELAX_SYM, ISEP, &
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, &
     &                Q, QF, RES, NF, NC, &
     &                SO, NSO, SOR, NSOR, CI, NCI, &
     &                IGRD, NOGm, &
     &                BMG_iWORK_PL, NBMG_iWORK_PL, &
     &                BMG_rWORK_PL, NBMG_rWORK_PL, &
     &                iWorkMSG, NMSGi, pMSG, pMSGSO, &
     &                MSG_Buffer, NMSGr, &
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID &
     &                ) ! ,Qmask)


! ==========================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG3_SymStd_updown performs all the necessary tasks associated
!     with moving "down" to a coarser grid or moving "up" to a finer
!     grid.  It is necessary because f77 does not support recursion.
!
! =======================================================================
! $license_flag$
! =======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   INPUT/OUTPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   OUTPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   LOCAL:
!  --------------------
!
!
! ==========================================================================

      IMPLICIT NONE

! -----------------------------
!     Includes
! -----------------------------

      INCLUDE 'mpif.h'

#include    "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_parameters.h'
      INCLUDE 'BMG_workspace.h'

! ----------------------------
!     Argument Declarations
! ----------------------------

      !
      !  BoxMG parameters
      !
      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      !
      !  BoxMG Cycle Parameters
      !
      INTEGER  k, KC, KCF, KF, IFD, IU, ID, IM, IVW, &
     &         IRELAX, IRELAX_SYM, UPDOWN, ISEP

      !
      !  BoxMG dimensions
      !
      INTEGER  NF, NC, NSO, NSOR, NCI, NOGm, NOG, &
     &         NBMG_iWORK_CS, NBMG_iWORK_PL, &
     &         NBMG_rWORK_CS, NBMG_rWORK_PL

      !
      !  BoxMG internal pointers
      !
      INTEGER  IGRD(NOGm,NBMG_pIGRD)

      !
      !  BoxMG data arrays
      !
      REAL*RKIND   CI(NCI), Q(NF), QF(NF), RES(NF), SO(NSO), SOR(NSOR)
!      REAL*RKIND :: Qmask(NF)

      !
      !  BoxMG plane relaxation workspace
      !
      INTEGER BMG_iWORK_PL(NBMG_iWORK_PL)
      REAL*RKIND  BMG_rWORK_PL(NBMG_rWORK_PL)

      !
      !  MSG workspace dimensions
      !
      INTEGER  NMSGi, NMSGr
      !
      !  MSG workspace and pointers
      !
      INTEGER  iWorkMSG(NMSGi), &
     &         pMSG(NBMG_pMSG,NOG), &
     &         pMSGSO(NBMG_pMSG,NOG)
      REAL*RKIND   MSG_Buffer(NMSGr)

      !
      !  BMG_MSG grid stuff
      !
      INTEGER  NBMG_MSG_iGRID
      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRID), &
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)

! ----------------------------
!     Local Declarations
! ----------------------------

      INTEGER  I, J, NLx, NLy, NLz, NLx_c, NLy_c, NLz_c, &
     &         NGx, NGy, NGz, NGx_c, NGy_c, NGz_c, &
     &         p_CI, p_CIC, p_SO, p_SOC, p_SOR, p_SORC, p_U, p_UC
      INTEGER  MPICOMM, XYPLANECOMM, XZPLANECOMM, YZPLANECOMM,  &
     &         iGs, jGs, kGs, NSORv, NStncl
      INTEGER  MyProc, MyProcI, MyProcJ, MyProcK, MPI_IERROR,  &
     &         NProc, ntasks_x, ntasks_y, NProcK, p_ProcGrid, p_ProcCoord
      INTEGER  NPASS
      REAL*8   RES_L2
      REAL*8   TT1, TT2

! ==========================================================================

! --------------------------------------------------------------------
! >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
! --------------------------------------------------------------------
! --------------------------------------------
!     MPI Parallel Data Decomposition:
! --------------------------------------------
      !
      ! MPI Communicators
      !
      MPICOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM)
      !
      ! MPI Communicatorsn for plane relaxation
      !
      XYPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_xy)
      XZPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_xz)
      YZPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_yz)

      !
      ! Processor Grid:
      !
      NProc  = BMG_MSG_iGRID(id_BMG_MSG_NProc)
      ntasks_x = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      ntasks_y = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)
      NProcK = 1! BMG_MSG_iGRID(id_BMG_MSG_NProcK)

      MyProc  = BMG_MSG_iGRID(id_BMG_MSG_MyProc)
      MyProcI = BMG_MSG_iGRID(id_BMG_MSG_MyProcI)
      MyProcJ = BMG_MSG_iGRID(id_BMG_MSG_MyProcJ)
      MyProcK = BMG_MSG_iGRID(id_BMG_MSG_MyProcK)

      !
      ! Processor Coordinates
      !
      iGs = IGRD(k,idL_BMG_Icoord)
      jGs = IGRD(k,idL_BMG_Jcoord)
      kGs = IGRD(k,idL_BMG_Kcoord)

      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)
      p_ProcCoord = BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord)

! --------------------------------------------------------------------
! >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
! --------------------------------------------------------------------

      !
      ! Set NStncl for the finest grid
      !
      IF ( K.LT.NOG .OR. IFD.NE.BMG_STENCIL_7pt ) THEN
         NStncl=14
      ELSE
         NStncl=4
      ENDIF

      !
      !  Number of temporary vectors in SOR
      !
      IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
         NSORv = 2
      ELSE
         NSORv = 2  ! the same for now but this may change
      ENDIF

      IF ( IM .gt. 0 ) THEN
         DO i=1, IM
            !
            TT1 = MPI_Wtime()

            CALL WRF_SymStd_relaxdg(  &
     &                KF, KF, NOG, UPDOWN,  &
     &                IFD, NStncl, IRELAX, iRELAX_SYM, RES_L2, ISEP,  &
     &                Q, QF, RES, NF,  &
     &                SO, NSO, SOR, NSOR, CI, NCI,  &
     &                IGRD, NOGm,  &
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,  &
     &                BMG_iWORK_PL, NBMG_iWORK_PL,  &
     &                BMG_rWORK_PL, NBMG_rWORK_PL,  &
     &                iGs, jGs, kGs,  &
     &                iWorkMSG, NMSGi, pMSG, MSG_Buffer, NMSGr,  &
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,  &
     &                MPICOMM, XYPLANECOMM, XZPLANECOMM, YZPLANECOMM  &
     &                ) ! ,Qmask)

            TT2 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_relax) = &
     &           BMG_rPARMS(id_BMG3_TIME_relax) + TT2 - TT1


            !
            IF( BMG_IOFLAG(iBMG3_BUG_RES_RELAX) .AND. MyProc.EQ.1 ) THEN
               WRITE(*,100) K, RES_L2
            ENDIF
            !
          ENDDO

        RETURN

      ENDIF


      IF (UPDOWN.EQ.BMG_DOWN) THEN

         ! Get pointers for grids k and k-1
         CALL BMG3_SymStd_GET_pointers(  &
     &               K, IGRD, NOGm,  &
     &               p_U, p_SO, p_SOR, p_CI,  &
     &               NLx, NLy, NLz, NGx, NGy, NGz  &
     &               )
         CALL BMG3_SymStd_GET_pointers(  &
     &               K-1, IGRD, NOGm,  &
     &               p_UC, p_SOC, p_SORC, p_CIC,  &
     &               NLx_c, NLy_c, NLz_c, NGx_c, NGy_c, NGz_c  &
     &               )

         ! Relaxation

         DO 10 i=1, ID
            !
            TT1 = MPI_Wtime()

            CALL WRF_SymStd_relaxdg(  &
     &                k, KF, NOG, UPDOWN,  &
     &                IFD, NStncl, IRELAX, iRELAX_SYM, RES_L2, ISEP,  &
     &                Q, QF, RES, NF,  &
     &                SO, NSO, SOR, NSOR, CI, NCI,  &
     &                IGRD, NOGm,  &
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,  &
     &                BMG_iWORK_PL, NBMG_iWORK_PL,  &
     &                BMG_rWORK_PL, NBMG_rWORK_PL,  &
     &                iGs, jGs, kGs,  &
     &                iWorkMSG, NMSGi, pMSG, MSG_Buffer, NMSGr,  &
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,  &
     &                MPICOMM, XYPLANECOMM, XZPLANECOMM, YZPLANECOMM  &
     &                ) ! ,Qmask)

            TT2 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_relax) =  &
     &           BMG_rPARMS(id_BMG3_TIME_relax) + TT2 - TT1


            !
            IF( BMG_IOFLAG(iBMG3_BUG_RES_RELAX) .AND. MyProc.EQ.1 ) THEN
               WRITE(*,100) K, RES_L2
            ENDIF
            !
 10      CONTINUE

         ! restrict the residual
         TT1 = MPI_Wtime()

         ! calculate the residual
         CALL BMG3_SymStd_residual(  &
     &             K, NOG, IFD,  &
     &             Q(p_U), QF(p_U), SO(p_SO), RES(p_U),  &
     &             NLx, NLy, NLz, NStncl,  &
     &             iWorkMSG, NMSGi, pMSG,  &
     &             MSG_Buffer, NMSGr, MPICOMM  &
     &             )

         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_PCG_PRECON) =  &
     &        BMG_rPARMS(id_BMG3_TIME_PCG_PRECON) + TT2 - TT1

!$$$         CALL BMG3_SymStd_DUMP_vector(  &
!$$$     &             BMG_IOFLAG, RES(p_U), K, NOG,  &
!$$$     &             NLx, NLy, NLz, NGx, NGy, NGz,
!$$$     &             IGRD(K,idL_BMG_ICoord),
!$$$     &             IGRD(K,idL_BMG_JCoord),
!$$$     &             IGRD(K,idL_BMG_KCoord),
!$$$     &             'DOWN-RES-post-relax', .FALSE.,
!$$$     &             BMG_MSG_iGRID(p_ProcGrid),
!$$$     &             ntasks_x, ntasks_y, NProcK, NProc,
!$$$     &             MyProcI, MyProcJ, MyProcK, MPICOMM
!$$$     &             )

!         ! restrict the residual
!         TT1 = MPI_Wtime()

         CALL BMG3_SymStd_restrict(  &
     &             K, K-1,  &
     &             RES(p_U), QF(p_UC), CI(p_CIC),  &
     &             NLx, NLy, NLz, NLx_c, NLy_c, NLz_c,  &
     &             iGs, jGs, kGs  &
     &             )

         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_restrict) =  &
     &        BMG_rPARMS(id_BMG3_TIME_restrict) + TT2 - TT1


!$$$         CALL BMG3_SymStd_DUMP_vector(  &
!$$$     &             BMG_IOFLAG, QF(p_UC), K-1, NOG,  &
!$$$     &             NLx_c, NLy_c, NLz_c, NGx_c, NGy, NGz,
!$$$     &             IGRD(K-1,idL_BMG_ICoord),
!$$$     &             IGRD(K-1,idL_BMG_JCoord),
!$$$     &             IGRD(K-1,idL_BMG_KCoord),
!$$$     &             'DOWN-QF-post-restrict', .FALSE.,
!$$$     &             BMG_MSG_iGRID(p_ProcGrid),
!$$$     &             ntasks_x, ntasks_y, NProcK, NProc,
!$$$     &             MyProcI, MyProcJ, MyProcK, MPICOMM
!$$$     &             )

         ! Zero the initial guess
         CALL BMG3_SymStd_UTILS_rV_zero(  &
     &               k-1, q(p_UC), NLx_c, NLy_c, NLz_c  &
     &               )


      ELSE IF (UPDOWN.EQ.BMG_UP) THEN

         ! Get pointers for grids k and k-1
         CALL BMG3_SymStd_GET_pointers(  &
     &             k, IGRD, NOGm,  &
     &             p_U, p_SO, p_SOR, p_CI,  &
     &             NLx, NLy, NLz, NGx, NGy, NGz  &
     &             )
         CALL BMG3_SymStd_GET_pointers(  &
     &             k-1, IGRD, NOGm,  &
     &             p_UC, p_SOC, p_SORC, p_CIC,  &
     &             NLx_c, NLy_c, NLz_c, NGx_c, NGy_c, NGz_c  &
     &             )


         ! Interpolate and Correct
         TT1 = MPI_Wtime()

            CALL BMG3_SymStd_interp_add(  &
     &             k-1, k, NOG,  &
     &             Q(p_U) ,Q(p_UC), RES(p_U),  &
     &             SO(p_SO), NStncl, CI(p_CIC),  &
     &             NLx_c, NLy_c, NLz_c,  &
     &             NLx, NLy, NLz,  &
     &             iGs, jGs, kGs,  &
     &             iWorkMSG, NMSGi, pMSG,  &
     &             MSG_Buffer, NMSGr, MPICOMM  &
     &             )


         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_interp_add) =  &
     &        BMG_rPARMS(id_BMG3_TIME_interp_add) + TT2 - TT1



         IF ( BMG_IOFLAG(iBMG3_BUG_RES_INTERP) ) THEN
            CALL BMG3_SymStd_residual(  &
     &                K, NOG, IFD,  &
     &                Q(p_U), QF(p_U), SO(p_SO), RES(p_U),  &
     &                NLx, NLy, NLz, NStncl,  &
     &                iWorkMSG, NMSGi, pMSG,  &
     &                MSG_Buffer, NMSGr, MPICOMM  &
     &                )
            CALL BMG3_SymStd_UTILS_norm_l2(  &
     &                RES(p_U), NLx, NLy, NLz, RES_L2, MPICOMM  &
     &                )
            WRITE(*,110) K, RES_L2
         ENDIF

         ! Relaxation
         DO 20 i=1, IU
            !
            TT1 = MPI_Wtime()

            CALL WRF_SymStd_relaxdg(  &
     &                k, KF, NOG, UPDOWN,  &
     &                IFD, NStncl, IRELAX, iRELAX_SYM, RES_L2, ISEP,  &
     &                Q, QF, RES, NF,  &
     &                SO, NSO, SOR, NSOR, CI, NCI,  &
     &                IGRD, NOGm,  &
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,  &
     &                BMG_iWORK_PL, NBMG_iWORK_PL,  &
     &                BMG_rWORK_PL, NBMG_rWORK_PL,  &
     &                iGs, jGs, kGs,  &
     &                iWorkMSG, NMSGi, pMSG, MSG_Buffer, NMSGr,  &
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,  &
     &                MPICOMM, XYPLANECOMM, XZPLANECOMM, YZPLANECOMM  &
     &                ) ! ,Qmask)
            TT2 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_relax) =  &
     &           BMG_rPARMS(id_BMG3_TIME_relax) + TT2 - TT1

            !
            IF( BMG_IOFLAG(iBMG3_BUG_RES_RELAX) .AND. MyProc.EQ.1 ) THEN
               WRITE(*,100) K, RES_L2
            ENDIF
            !
 20      CONTINUE

      ELSE
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN
            WRITE(*,*) '*** UPDOWN out of range: '
            WRITE(*,*) 'HAVE: UPDOWN = ', UPDOWN
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,19)
         RETURN

      ENDIF

! ==========================================================================

 100  FORMAT (1X,'(3D)LEVEL ',I2,' RESIDUAL NORM = ',1P,E12.5)
 110  FORMAT (1X,'(3D)LEVEL ',I2,  &
     &           ' AFTER INTERPOLATION RES_L2 = ',1P,E12.5 )

! ===========================================

      RETURN
      END SUBROUTINE WRF_SymStd_updowndg

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE WRF_SymStd_relaxdg(  &
     &                KG, KF, NOG, UPDOWN,  &
     &                IFD, NStncl, IRELAX, iRELAX_SYM, RES_L2, ISEP,  &
     &                Q, QF, RES, NF,  &
     &                SO, NSO, SOR, NSOR, CI, NCI,  &
     &                IGRD, NOGm,  &
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,  &
     &                BMG_iWORK_PL, NBMG_iWORK_PL,  &
     &                BMG_rWORK_PL, NBMG_rWORK_PL,  &
     &                iGs, jGs, kGs,  &
     &                iWorkMSG, NMSGi, pMSG, MSG_Buffer, NMSGr,  &
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,  &
     &                MPICOMM, XYPLANECOMM, XZPLANECOMM, YZPLANECOMM  &
     &                ) ! ,Qmask)

! ==========================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG3_SymStd_relax performs one relaxation step (either colored
!     point Gauss Seidel, or alternating plane red-black Gauss-Seidel)
!
! =======================================================================
! $license_flag$
! =======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   INPUT/OUTPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   OUTPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   LOCAL:
!  --------------------
!
!
! ==========================================================================

      IMPLICIT NONE

! ------------------------------------------------
!     Includes
!
#include    "BMG_constants.h"
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

! ------------------------------------------------
!     Argument Declarations
!
      INTEGER  NCI, NSO, NSOR, NSORv, NStncl, kg, kf, NBMG_iWORK_PL,  &
     &         NBMG_rWORK_PL, NOG, NF, NOGm, NMSGi, NMSGr

      INTEGER  NBMG_MSG_iGRID
      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRID),  &
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)

      INTEGER  BMG_iPARMS(NBMG_iPARMS), BMG_iWORK_PL(NBMG_iWORK_PL),  &
     &         IFD, IGRD(NOGm,24), IRELAX_SYM, UPDOWN, IRELAX, ISEP
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS), BMG_rWORK_PL(NBMG_rWORK_PL),  &
     &         CI(NCI), q(NF), qf(NF), RES(NF),  &
     &         SO(NSO), SOR(NSOR), MSG_Buffer(NMSGr)

      REAL*8  RES_L2

!      REAL*RKIND :: Qmask(NF)

      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      INTEGER  iGs, jGs, kGs,  &
     &         MPICOMM, iWorkMSG(NMSGi), pMSG(NBMG_pMSG,NOGm)
      INTEGER  XYPLANECOMM, XZPLANECOMM, YZPLANECOMM
      INTEGER  MyProc, MyProcI, MyProcJ, MyProcK, MPI_IERROR,  &
     &         NProc, ntasks_x, ntasks_y, NProcK


! ----------------------------
!     Local Declarations
!
      INTEGER  NLx, NLy, NLz, NGx, NGy, NGz,  &
     &         NLx_f, NLy_f, NLz_f, NGx_f, NGy_f, NGz_f,  &
     &         p_CI, p_SO, p_SOR, p_U, p_ProcGrid

! ==========================================================================

! --------------------------------------------------------------------
! >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
! --------------------------------------------------------------------

      !
      ! MPI Communicators
      !
      MPICOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM)
      !
      ! MPI Communicatorsn for plane relaxation
      !
      XYPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_xy)
      XZPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_xz)
      YZPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_yz)

      !
      ! Processor Grid:
      !
      NProc  = BMG_MSG_iGRID(id_BMG_MSG_NProc)
      ntasks_x = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      ntasks_y = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)
      NProcK = 1! BMG_MSG_iGRID(id_BMG_MSG_NProcK)

      MyProc  = BMG_MSG_iGRID(id_BMG_MSG_MyProc)
      MyProcI = BMG_MSG_iGRID(id_BMG_MSG_MyProcI)
      MyProcJ = BMG_MSG_iGRID(id_BMG_MSG_MyProcJ)
      MyProcK = BMG_MSG_iGRID(id_BMG_MSG_MyProcK)

      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)

! --------------------------------------------------------------------
! >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
! --------------------------------------------------------------------

      !
      !  Number of temporary vectors in SOR
      !
      IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
         NSORv = 2
      ELSE
         NSORv = 2  ! the same for now but this may change
      ENDIF

      !
      ! Collect fine-grid dimensions
      ! ( need to check: KF or NOG )
      !
      CALL BMG3_SymStd_GET_pointers(  &
     &            KF, IGRD, NOGm,  &
     &            p_U, p_SO, p_SOR, p_CI,  &
     &            NLx_f, NLy_f, NLz_f, NGx_f, NGy_f, NGz_f  &
     &            )

      CALL BMG3_SymStd_GET_pointers(  &
     &            KG, IGRD, NOGm,  &
     &            p_U, p_SO, p_SOR, p_CI,  &
     &            NLx, NLy, NLz, NGx, NGy, NGz  &
     &            )


      IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
         !
         ! Gauss-Seidel relaxation
         !
         CALL WRF_SymStd_relax_GSdg(  &
     &             KG, SO(p_SO), QF(p_U), Q(p_U), RES(p_U), SOR(p_SOR),  &
     &             NLx, NLy, NLz, NGx, NGy, NGz,  &
     &             RES_L2, BMG_IOFLAG, NOG, NOGm,  &
     &             IFD, NStncl, NSORv, IRELAX, IRELAX_SYM, UPDOWN,  &
     &             iGs, jGs, kGs,  &
     &             MyProcI, MyProcJ, MyProcK, MyProc,  &
     &             BMG_MSG_iGRID(p_ProcGrid),  &
     &             ntasks_x, ntasks_y, NProcK, NProc,  &
     &             iWorkMSG, NMSGi, pMSG, MSG_Buffer, NMSGr, MPICOMM  &
     &            ,ISEP ) ! ,Qmask)
         !
      ELSEIF ( IRELAX.EQ.BMG_GS_RB_planes_xy_yz_xz ) THEN

         WRITE(*,*) 'Fatal Error: To busy to code planes !!!'
         CALL MPI_FINALIZE(MPI_IERROR)
         STOP

         !
      ENDIF

! ==========================================================================

      RETURN
      END SUBROUTINE WRF_SymStd_relaxdg


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      SUBROUTINE WRF_SymStd_relax_GSdg(  &
     &                KG, SO, QF, Q, RES, SOR,  &
     &                NLx, NLy, NLz, NGx, NGy, NGz,  &
     &                RES_L2, BMG_IOFLAG, NOG, NOGm,  &
     &                IFD, NStncl, NSORv, IRELAX, IRELAX_SYM, UPDOWN,  &
     &                iGs, jGs, kGs,  &
     &                MyProcI, MyProcJ, MyProcK, MyProc,  &
     &                ProcGrid, ntasks_x, ntasks_y, NProcK, NProc,  &
     &                iWork, NMSGi, pMSG, MSG_Buffer, NMSGr, MPICOMM  &
     &               ,ISEP ) ! ,Qmask)

! ==========================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG3_SymStd_relax_GS performs one sweep of Gauss Seidel (with the
!     correct ordering depending on whether we are on the way down, or
!     up in the symmetric cycling case)
!
! =======================================================================
! $license_flag$
! =======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   INPUT/OUTPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   OUTPUT:
!  --------------------
!
!
! =======================================================================
!  --------------------
!   LOCAL:
!  --------------------
!
!
! ==========================================================================

!      USE CPUTIME_MODULE

      IMPLICIT NONE

! ------------------------------------------------
!     Includes
!
      INCLUDE 'mpif.h'
      INCLUDE 'MSG_f90.h'

#include    "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_parameters.h'
      INCLUDE 'BMG_workspace.h'

! ------------------------------------------------
!     Argument Declarations
!
      INTEGER NGx, NGy, NGz, NLx, NLy, NLz

      INTEGER IFD, IRELAX, IRELAX_SYM, KG, NOG, NOGm,  &
     &        NSORv, NStncl, UPDOWN, NMSGi, NMSGr

      INTEGER iGs, jGs, kGs
      INTEGER MyProc, MyProcI, MyProcJ, MyProcK,  &
     &        NProc, ntasks_x, ntasks_y, NProcK
      INTEGER ProcGrid(ntasks_x,ntasks_y,NProcK)

      INTEGER MPICOMM, iWork(NMSGi), pMSG(NBMG_pMSG,NOGm)
      REAL*RKIND  Q(NLx,NLy,NLz), QF(NLx,NLy,NLz),  &
     &        RES(NLx,NLy,NLz),  &
     &        SO(NLx+1,NLy+1,NLz+1,NStncl), SOR(NLx,NLy,NLz,NSORv)
      REAL*RKIND  MSG_Buffer(NMSGr)

      REAL*8 RES_L2

!      REAL*RKIND :: Qmask(NLx,Nly,Nlz)

      LOGICAL BMG_IOFLAG(NBMG_IOFLAG)

! ------------------------------------------------
!     Local Declarations
!
      INTEGER i, i1, j ,j1, k, k1

      INTEGER ie, iw, js, jn, kbot, kt

      INTEGER pts, ibeg, iend, ptsred, ptsblk
      INTEGER ptstart, ptend, ptstride
      INTEGER ptrn, ierror
      INTEGER relaxtest
      LOGICAL comhalf

      REAL*RKIND SOx(NLx+1),SOy(NLy+1),SOz(NLz+1)

      INTEGER ISEP

! ==========================================================================

      j1 = NLy-1
      i1 = NLx-1
      k1 = NLz-1

!      write(0,*) 'mproc, mpi,mpj,mpk, igs,jgs,kgs = ',myproc,myproci,
!     :   myprocj,myprock,igs,jgs,kgs
!      write(0,*) 'mproc, isep, BMG_SEPARABLE = ',myproc,isep,
!     :   BMG_SEPARABLE, BMG_NONSEPARABLE


!      write(0,*) 'RELAX_GS: kg,nog,ifd = ',KG,NOG,ifd

      comhalf = .true.

!      CALL cld_cpu('BOXMG-GS')


      IF ( KG.LT.NOG .OR. ifd.NE.1 ) THEN

         !
         !   27-point relaxation (8-color)
         !

!         IF ((UPDOWN.eq.BMG_UP).or.(IRELAX_SYM.EQ.BMG_RELAX_NONSYM))
!     &              THEN
!            ptstart = 1
!            ptend   = 8
!            ptstride = 1
!         ELSE
!            ptstart = 8
!            ptend = 1
!            ptstride = -1
!         ENDIF
!
!         DO pts = ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<
!
!            DO k=2+mod(mod((pts-1)/4,2)+mod(kGs+1,2),2),K1,2
!               !
!               DO j=2+mod(mod((pts-1)/2,2)+mod(jGs+1,2),2),J1,2
!                  !
!                  DO i=2+mod(mod(pts-1,2)+mod(iGs+1,2),2),I1,2
                     !

         IF ( (UPDOWN.eq.BMG_UP)  &
     &        .OR. (IRELAX_SYM.EQ.BMG_RELAX_NONSYM)) THEN
            ptstart  = 0
            ptend    = 1
            ptstride = 1
         ELSE
            ptstart  = 1
            ptend    = 0
            ptstride = -1
         ENDIF


        relaxtest = 0

        IF ( relaxtest .eq. 0 ) THEN ! default

         DO pts = ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<

!            DO k=2+mod(mod((pts-1)/4,2)+mod(kGs+1,2),2),K1,2
!               !
!               DO j=2+mod(mod((pts-1)/2,2)+mod(jGs+1,2),2),J1,2
!                  !
!                  DO i=2+mod(mod(pts-1,2)+mod(iGs+1,2),2),I1,2
             DO k = 2,k1

               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                     !
                     Q(i,j,k) = ( QF(i,j,k)  &
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)  &
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)  &
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)  &
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)  &
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)  &
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)  &
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)  &
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)  &
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)  &
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)  &
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)  &
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)  &
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)  &
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)  &
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)  &
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)  &
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)  &
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)  &
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)  &
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)  &
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)  &
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)  &
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)  &
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)  &
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)  &
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)  &
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO
            !

            IF ( comhalf ) THEN

            ptrn = 1

            CALL MSG_tbdx_send(Q, MSG_Buffer,  &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Proc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),  &
     &           iWork(pMSG(ipL_MSG_Index,KG)),  &
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,  &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Proc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),  &
     &           iWork(pMSG(ipL_MSG_Index,KG)),  &
     &           ptrn, ierror)
            !

            ENDIF
          ENDDO

        ELSEIF ( relaxtest .eq. 1 ) THEN ! non-redblack GS

             DO k = 2,k1

               DO j=2, j1
                  !
                  !
                  DO i = 2,i1
                     !
                     Q(i,j,k) = ( QF(i,j,k)  &
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)  &
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)  &
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)  &
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)  &
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)  &
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)  &
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)  &
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)  &
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)  &
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)  &
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)  &
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)  &
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)  &
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)  &
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)  &
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)  &
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)  &
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)  &
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)  &
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)  &
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)  &
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)  &
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)  &
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)  &
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)  &
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)  &
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO

        ELSEIF ( relaxtest .eq. 2 ) THEN ! layered redblack GS

!         DO pts=ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<

! do 'red' points on first plane
            ptsred = ptstart
            ptsblk = ptend
            DO  k=2, 2
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+ptsred,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                     Q(i,j,k) = ( QF(i,j,k)  &
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)  &
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)  &
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)  &
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)  &
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)  &
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)  &
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)  &
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)  &
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)  &
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)  &
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)  &
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)  &
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)  &
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)  &
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)  &
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)  &
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)  &
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)  &
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)  &
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)  &
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)  &
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)  &
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)  &
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)  &
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)  &
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)  &
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO

! do 'red' points on k plane and 'black' points on k-1 plane

            DO  k=3, k1

               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+ptsred,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                     Q(i,j,k) = ( QF(i,j,k)  &
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)  &
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)  &
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)  &
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)  &
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)  &
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)  &
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)  &
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)  &
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)  &
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)  &
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)  &
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)  &
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)  &
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)  &
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)  &
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)  &
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)  &
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)  &
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)  &
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)  &
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)  &
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)  &
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)  &
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)  &
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)  &
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)  &
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO

            IF ( comhalf ) THEN

            ptrn = 1

            CALL MSG_tbdx_send(Q, MSG_Buffer,  &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Proc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),  &
     &           iWork(pMSG(ipL_MSG_Index,KG)),  &
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,  &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Proc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),  &
     &           iWork(pMSG(ipL_MSG_Index,KG)),  &
     &           ptrn, ierror)
            !

            ENDIF

! black points
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-2)+ptsblk,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                     Q(i,j,k-1) = ( QF(i,j,k-1)  &
     &                    + SO(i,j,k-1,kpw)*Q(i-1,j,k-1)  &
     &                    + SO(i,j+1,k-1,kpnw)*Q(i-1,j+1,k-1)  &
     &                    + SO(i,j+1,k-1,kps)*Q(i,j+1,k-1)  &
     &                    + SO(i+1,j+1,k-1,kpsw)*Q(i+1,j+1,k-1)  &
     &                    + SO(i+1,j,k-1,kpw)*Q(i+1,j,k-1)  &
     &                    + SO(i+1,j,k-1,kpnw)*Q(i+1,j-1,k-1)  &
     &                    + SO(i,j,k-1,kps)*Q(i,j-1,k-1)  &
     &                    + SO(i,j,k-1,kpsw)*Q(i-1,j-1,k-1)  &
     &                    + SO(i,j,k-1,kb)*Q(i,j,k-1-1)  &
     &                    + SO(i,j,k-1,kbw)*Q(i-1,j,k-1-1)  &
     &                    + SO(i,j+1,k-1,kbnw)*Q(i-1,j+1,k-1-1)  &
     &                    + SO(i,j+1,k-1,kbn)*Q(i,j+1,k-1-1)  &
     &                    + SO(i+1,j+1,k-1,kbne)*Q(i+1,j+1,k-1-1)  &
     &                    + SO(i+1,j,k-1,kbe)*Q(i+1,j,k-1-1)  &
     &                    + SO(i+1,j,k-1,kbse)*Q(i+1,j-1,k-1-1)  &
     &                    + SO(i,j,k-1,kbs)*Q(i,j-1,k-1-1)  &
     &                    + SO(i,j,k-1,kbsw)*Q(i-1,j-1,k-1-1)  &
     &                    + SO(i,j,k-1+1,kb)*Q(i,j,k-1+1)  &
     &                    + SO(i,j,k-1+1,kbe)*Q(i-1,j,k-1+1)  &
     &                    + SO(i,j+1,k-1+1,kbse)*Q(i-1,j+1,k-1+1)  &
     &                    + SO(i,j+1,k-1+1,kbs)*Q(i,j+1,k-1+1)  &
     &                    + SO(i+1,j+1,k-1+1,kbsw)*Q(i+1,j+1,k-1+1)  &
     &                    + SO(i+1,j,k-1+1,kbw)*Q(i+1,j,k-1+1)  &
     &                    + SO(i+1,j,k-1+1,kbnw)*Q(i+1,j-1,k-1+1)  &
     &                    + SO(i,j,k-1+1,kbn)*Q(i,j-1,k-1+1)  &
     &                    + SO(i,j,k-1+1,kbne)*Q(i-1,j-1,k-1+1)  &
     &                    )*SOR(i,j,k-1,msor)
                     !
                  ENDDO
                  !
               ENDDO

               !
            ENDDO

! do 'black' points on last plane
            DO  k=k1, k1
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+ptsblk,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                     Q(i,j,k) = ( QF(i,j,k)  &
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)  &
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)  &
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)  &
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)  &
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)  &
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)  &
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)  &
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)  &
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)  &
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)  &
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)  &
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)  &
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)  &
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)  &
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)  &
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)  &
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)  &
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)  &
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)  &
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)  &
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)  &
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)  &
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)  &
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)  &
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)  &
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)  &
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO

       ENDIF ! relaxtest
       !
!         ENDDO   ! >>> END: loop over colors <<<<<<<<<<<<<<<<<<<<<<<<<<

        IF ( relaxtest .ne. 0 ) THEN

            ptrn = 1
            !
            CALL MSG_tbdx_send(Q, MSG_Buffer,  &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Proc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),  &
     &           iWork(pMSG(ipL_MSG_Index,KG)),  &
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,  &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Proc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),  &
     &           iWork(pMSG(ipL_MSG_Index,KG)),  &
     &           ptrn, ierror)

        ENDIF
            !
!         ENDDO   ! >>> END: loop over colors <<<<<<<<<<<<<<<<<<<<<<<<<<


         !
      ELSE

!      CALL cld_cpu('BOXMG-GSfine')
         !
         !  7 point relaxation (2-color)
         !
         IF ( (UPDOWN.eq.BMG_UP) &
     &        .OR. (IRELAX_SYM.EQ.BMG_RELAX_NONSYM)) THEN
            ptstart  = 0
            ptend    = 1
            ptstride = 1
         ELSE
            ptstart  = 1
            ptend    = 0
            ptstride = -1
         ENDIF

      IF ( ISEP .eq. BMG_SEPARABLE ) THEN


      DO  k=1, k1+1
        SOz(k) = SO(3,3,k,kb)
        SOx(k) = SO(3,3,k,kpw)
        SOy(k) = SO(3,3,k,kps)
!       IF ( myproc == 1 ) write(0,*) 'k, SOz,x,y = ',k,SOz(k),
!     &  SO(3,3,k,kpw),SO(3,3,k,kps)
!        IF ( myproc == 1 ) write(0,*) 'j, SOy = ',
      ENDDO
!        SOz(2) = SOz(3)
!        SOz(k1+1) = SOz(k1)


      DO j=1, j1+1
!        SOy(j) = SO(3,j,3,kps)
!        IF ( myproc == 1 ) write(0,*) 'j, SOy,x,z = ',j,SOy(j),
!     &   SO(3,j,3,kpw),SO(3,j,3,kb)
!       IF ( myproc == 1 ) write(0,*) 'k, SOz = ',j,SO(3,j,3,kb)
!        IF ( myproc == 1 ) write(0,*) 'i, SOx = ',j,SO(3,j,3,kpw)
!        write(0,*) 'j, SOy, myproc = ',j,SOy(j),myproc
      ENDDO
!        SOy(2) = SOy(3)
!        SOy(j1+1) = SOy(j1)


      DO i=1,i1+1
!        SOx(i) = SO(i,3,3,kpw)
!        IF ( myproc == 1 ) write(0,*) 'i, SOx,y,z = ',i,SOx(i),
!     &   SO(i,3,3,kps),SO(i,3,3,kb)
      ENDDO
!        SOx(2) = SOx(3)
!        SOx(i1+1) = SOx(i1)

      ENDIF

         DO pts=ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<

          IF ( ISEP .eq. BMG_SEPARABLE ) THEN

          IF ( .true. ) THEN


            DO  k=2, k1
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !
                     Q(i,j,k) = ( QF(i,j,k)  &
     &                        + SOx(k)*Q(i-1,j,k)  &
     &                        + SOy(k)*Q(i,j+1,k)  &
     &                        + SOx(k)*Q(i+1,j,k)  &
     &                        + SOy(k)*Q(i,j-1,k)  &
     &                        + SOz(k)*Q(i,j,k-1)  &
     &                        + SOz(k+1)*Q(i,j,k+1)  &
     &                        )*SOR(i,j,k,msor)
                  END DO
                  !
               END DO
               !
            END DO

          ELSE ! old try

          kbot = 0
          kt = 0
          js = 0
          jn = 0
          iw = 0
          ie = 0

          IF (  MyProcK .eq. 1 )     kbot = 1
          IF (  MyProcK .eq. NProcK) kt = 1
          IF (  MyProcJ .eq. 1)      js = 1
          IF (  MyProcJ .eq. ntasks_y) jn = 1
          IF (  MyProcI .eq. 1)      iw = 1
          IF (  MyProcI .eq. ntasks_x) ie = 1

            DO  k=2+kbot, k1-kt
               DO j=2+js, j1-jn
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg

                  IF ( ibeg .eq. 2 ) ibeg = ibeg + 2*iw
                  IF ( iend .eq. i1 ) iend = iend - 2*ie

                  !
                  DO i=ibeg,iend,2
                  !
!                  IF ( .not. (( (k.eq.2      .and. MyProcK.eq.1 ) .or.
!     &                        ( k .eq. k1    .and. MyProcK .eq. NProcK) .or.
!     &                        ( j .eq. 2     .and. MyProcJ .eq. 1)  .or.
!     &                        ( j .eq. j1    .and. MyProcJ .eq. ntasks_y) .or.
!     &                        ( i .eq. ibeg  .and. MyProcI .eq. 1)  .or.
!     &                        ( i .eq. iend  .and. MyProcI .eq. ntasks_x) ) ) THEN
                     Q(i,j,k) = ( QF(i,j,k)  &
     &                        + SOx(i)*Q(i-1,j,k)  &
     &                        + SOy(j+1)*Q(i,j+1,k)  &
     &                        + SOx(i+1)*Q(i+1,j,k)  &
     &                        + SOy(j)*Q(i,j-1,k)  &
     &                        + SOz(k)*Q(i,j,k-1)  &
     &                        + SOz(k+1)*Q(i,j,k+1)  &
     &                        )/(SOx(i) + SOx(i+1) + SOy(j) +  &
     &                           SOy(j+1) + SOz(k) + SOz(k+1))

                  END DO
                  !
               END DO
               !
            END DO

! edges:

           IF ( kbot .eq. 1 ) THEN
            k=2

               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2

                     Q(i,j,k) = ( QF(i,j,k)  &
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)  &
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)  &
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)  &
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)  &
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)  &
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)  &
     &                        )*SOR(i,j,k,msor)


                 END DO
              END DO

           ENDIF

           IF ( kt .eq. 1 ) THEN
            k=k1

               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2

                     Q(i,j,k) = ( QF(i,j,k)  &
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)  &
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)  &
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)  &
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)  &
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)  &
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)  &
     &                        )*SOR(i,j,k,msor)


                 END DO
              END DO

           ENDIF

           IF ( js .eq. 1 ) THEN

            DO  k=2+kbot, k1-kt
              j = 2
              ibeg = mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
              iend = 2*((i1-ibeg)/2)+ibeg

               DO i=ibeg,iend,2

                     Q(i,j,k) = ( QF(i,j,k)  &
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)  &
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)  &
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)  &
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)  &
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)  &
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)  &
     &                        )*SOR(i,j,k,msor)

                END DO
             END DO
          ENDIF


           IF ( jn .eq. 1 ) THEN

            DO  k=2+kbot, k1-kt
              j = j1
              ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
              iend=2*((i1-ibeg)/2)+ibeg
                  DO i=ibeg,iend,2


                     Q(i,j,k) = ( QF(i,j,k)  &
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)  &
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)  &
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)   &
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)  &
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)  &
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)  &
     &                        )*SOR(i,j,k,msor)

                END DO
             END DO
          ENDIF


           IF ( ie .eq. 1 ) THEN
            DO  k=2+kbot, k1-kt
               DO j=2+js, j1-jn

              ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
              iend=2*((i1-ibeg)/2)+ibeg
              i = ibeg
               IF ( ibeg .eq. 2 ) THEN

                     Q(i,j,k) = ( QF(i,j,k)  &
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)  &
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)  &
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)  &
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)  &
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)  &
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)  &
     &                        )*SOR(i,j,k,msor)

                ENDIF
                END DO
             END DO
          ENDIF

           IF ( iw .eq. 1 ) THEN
            DO  k=2+kbot, k1-kt
               DO j=2+js, j1-jn

              ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
              iend=2*((i1-ibeg)/2)+ibeg

              i = iend

               IF ( i .eq. i1 ) THEN

                     Q(i,j,k) = ( QF(i,j,k)  &
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)  &
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)  &
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)  &
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)  &
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)  &
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)  &
     &                        )*SOR(i,j,k,msor)

                ENDIF
                END DO
             END DO
          ENDIF

         ENDIF


         ELSE  ! Non separable

            DO  k=2, k1
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                     Q(i,j,k) = ( QF(i,j,k)  &
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)  &
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)  &
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)  &
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)  &
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)  &
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)  &
     &                        )*SOR(i,j,k,msor)
                  END DO
                  !
               END DO
               !
            END DO

          ENDIF
            !

            IF ( comhalf ) THEN

            ptrn = 1

            CALL MSG_tbdx_send(Q, MSG_Buffer,  &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Proc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),  &
     &           iWork(pMSG(ipL_MSG_Index,KG)),  &
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,  &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Proc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),  &
     &           iWork(pMSG(ipL_MSG_Index,KG)),  &
     &           ptrn, ierror)
            !

            ENDIF


            !
         END DO   ! >>> END: loop over colors <<<<<<<<<<<<<<<<<<<<<<<<<<
            !
            IF ( .not. comhalf ) THEN

            ptrn = 1

            CALL MSG_tbdx_send(Q, MSG_Buffer,  &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Proc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),  &
     &           iWork(pMSG(ipL_MSG_Index,KG)),  &
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,  &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Proc,KG)),  &
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),  &
     &           iWork(pMSG(ipL_MSG_Index,KG)),  &
     &           ptrn, ierror)
            !

            ENDIF
            !
!         END DO   ! >>> END: loop over colors <<<<<<<<<<<<<<<<<<<<<<<<<<

!      CALL cld_cpu('BOXMG-GSfine')

      ENDIF

      ptrn = 1

      CALL MSG_tbdx_close(Q, MSG_Buffer,  &
     &     iWork(pMSG(ipL_MSG_NumAdjProc,KG)),  &
     &     iWork(pMSG(ipL_MSG_Proc,KG)),  &
     &     iWork(pMSG(ipL_MSG_Ipr,KG)),  &
     &     iWork(pMSG(ipL_MSG_Index,KG)),  &
     &     ptrn, ierror)

      IF( BMG_IOFLAG(iBMG3_BUG_RES_RELAX) ) THEN
         !
         CALL BMG3_SymStd_residual(  &
     &             KG, NOG, IFD,  &
     &             Q, QF, SO, RES, NLx, NLy, NLz, NStncl,  &
     &             iWork, NMSGi, pMSG,  &
     &             MSG_Buffer, NMSGr, MPICOMM &
     &             )
         CALL BMG3_SymStd_UTILS_norm_l2(  &
     &             RES, NLx, NLy, NLz, RES_L2, MPICOMM  &
     &             )
         !


      ENDIF

!      CALL cld_cpu('BOXMG-GS')

! ==========================================================================

      RETURN
      END SUBROUTINE WRF_SymStd_relax_GSdg
! ==========================================================================
! ==========================================================================
! ==========================================================================
! ==========================================================================

      

#endif
!      END MODULE module_mp_boxmgsetup


!-------------------------------------------------------------------------
!
!  subroutine mlint1 linearly interpolates, in three dimensions,
!  for any variable 'ac'
!
      subroutine mlint2                             &
       (idebug,ac,n1b,n1e,n2b,n2e,n3b,n3e,n4b,n4e,  &
        nx,ny,nz,aint,facx,facy,facz,ic,jc,kc,lc)
!
      integer  idebug
      integer  ic
      integer  kc
      integer  jc
      integer  lc
!
      real     acx1
      real     acx2
      real     acx3
      real     acx4
      real     acy1
      real     acy2
      real     aint
      real     facx
      real     facy
      real     facz
      real     ac(n1b:n1e,n2b:n2e,n3b:n3e,n4b:n4e)
      
      
      acx1 = 0.0
      acx2 = 0.0
      acx3 = 0.0
      acx4 = 0.0
      acy1 = 0.0
      acy2 = 0.0
      
!
!  x-direction interpolation
!
      if( nx .gt. 1 ) then
!
      acx1 = facx*(ac(ic+1,jc,kc,lc)-ac(ic,jc,kc,lc)) +  &
             ac(ic,jc,kc,lc)
      acx2 = facx*(ac(ic+1,jc,kc+1,lc)-ac(ic,jc,kc+1,lc)) + &
             ac(ic,jc,kc+1,lc)
!
      else if( nx .le. 1 .and. ny .gt. 1 ) then
!
      acy1 = facy*(ac(ic,jc+1,kc,lc)-ac(ic,jc,kc,lc)) +       &
             ac(ic,jc,kc,lc)
      acy2 = facy*(ac(ic,jc+1,kc+1,lc)-ac(ic,jc+1,kc+1,lc)) + &
            ac(ic,jc+1,kc+1,lc)
!
      
      else  ! 1-d vertical interpolation

      acy1 = ac(ic,jc,kc,lc)
      acy2 = ac(ic,jc,kc+1,lc)
      
      end if
!
!  y-direction interpolation of data interpolated in x-direction
!  for bottem of cube if in 3-d; eg x:y:z trajectory
!
      if( nx .gt. 1 .and. ny .gt.1 ) then
!
      acx3 = facx*(ac(ic+1,jc+1,kc,lc)-ac(ic,jc+1,kc,lc)) +     &
             ac(ic,jc+1,kc,lc)
      acx4 = facx*(ac(ic+1,jc+1,kc+1,lc)-ac(ic,jc+1,kc+1,lc)) + &
             ac(ic,jc+1,kc+1,lc)
!
! y-direction interpolation of x interpolated data
!
      acy1 = facy*(acx3-acx1) + acx1
      acy2 = facy*(acx4-acx2) + acx2
!
!  z-direction interpolation of data interpolated in x:y-directions
!
      aint = facz*(acy2-acy1) + acy1
!
!  z-direction interpolation of data interpolated in y-direction for
!  2-d y:z trajectory
!
      else if ( nx .le. 1 ) then
      aint = facz*(acy2-acy1) + acy1
!
!  z-direction interpolation of data interpolated in x-direction for
!  2-d x:z trajectory
!
      else if ( ny .le. 1 ) then
      aint = facz*(acx2-acx1) + acx1
      end if
!
      if ( idebug .ge. 2 ) then
      IF ( nx .gt. 1 ) print*,'acx1,acx2: ',acx1, acx2, facx
      print*,'acx3,acx4: ',acx3, acx4, facx
      print*,'acy1,acy2: ',acy1, acy2, facy
      print*,'aint: ', aint, facz
      end if
!
      return
      end
!
      
#endif
