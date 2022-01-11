!
#define FLOAT_PRECISION MPI_DOUBLE_PRECISION

      Module  module_mp_screen
#ifdef DM_PARALLEL 
      USE module_dm, ONLY : &
         local_communicator, mytask, ntasks, ntasks_x, ntasks_y,mytask_x, mytask_y
      USE module_commasmpi
#endif

      CONTAINS

      subroutine screen(iscreen,nx,ny,nz,dtp,dx,dy,ng   &
     &    ,ids,ide,jds,jde,kds,kde,its,ite,jts,jte,ims,ime,kms,kme,jms,jme    &
     &    ,dz,pot,qc,qi,scw,sci,sctot,dn,iunit)
!
!
!  5.13.2002  Added vertical grid stretch corrections
!
!  5/15/2001  Added reporting of the amount of charge per second added by
!             the parameterization.

!  7/19/2001   version sam.screen.1a.f
!             Changed the way the electric field in- and outside the cloud
!             boundary is calculated.  Instead of averaging the second-order
!             field components between the inside and outside points, now take
!             first-order differences without using the boundary point at all.
!
!  7/19/2001  version sam.screen.1b.f
!             Test using the boundary point for the calculation of the field
!             inside, but not the outside point.
!
!  7/21/2001  Changed nsmth1 and nsmth2 to 1 instead of 4
!             Also now distribute (smoothed) screening charge to any point with
!             cloud mixing ratio .gt. qctol (instead of just the points identified
!             as being boundary points) which allows interior points to get
!             some of the charge, too.  May try a value smaller than qctol in the
!             future so that some charge goes to the point just outside the
!             arbitrary cloud boundary.
!
!  7/24/2001 version sam.screen.1c.f
!            Change electric field calculations from 1st order (2 point) to
!            2nd-order forward (or backward) differences.
!
!  8/22/2001 test of reduced smoothing by 50%
!
!  9/13/2001 sam.screen.1d.f:
!          further reduced smoothing to 1/32 (instead of previous 1/16)
!

      implicit none

#ifdef DM_PARALLEL
      INCLUDE 'mpif.h'

      integer, parameter :: ntot = 50
      double precision  mpitotindp(ntot), mpitotoutdp(ntot)
      INTEGER :: mpi_error_code = 1      
#endif

      integer iunit,iscreen
      real dccld,dccldn, ecnrml, eanrml
      real screenn,screenp,screent
      integer ix,jy,kz,kk
      integer, INTENT(IN) :: its,ite,jts,jte
      integer, INTENT(IN) :: ims,ime,jms,jme,kms,kme
      integer, INTENT(IN) :: ids,ide,jds,jde,kds,kde
      integer nx,ny,nz,ng
      integer westward_tag,eastward_tag,southward_tag,northward_tag


      logical, parameter :: firstord = .true.
      
      integer :: id1 =1
      integer :: jd1 =1
      integer :: id2 =1
      integer :: jd2 =1
      integer :: kd1 =1
      integer :: istag =1
      integer :: jstag =1
      integer :: kstag =1
      logical :: ndebug = .false.
      real :: eperao = 8.8592e-12
      real :: epascl = 5.0
      real :: qctol = 1.e-5 !5.e-5
      real :: qctol2 = 1.e-6 !5.e-5
      real :: fin1 = 1.5
      real :: fin2 = 1.5*1.5
      real, parameter :: sc_thresh = 0.25e-9 ! 1.e-9 ! 0.2e-10
!     

!     NOTE NX,NY and NZ here ARE = ime-ims+1,jme-jms+1,kme-kms+1

!      real  dslay(1:nx,1:ny,1:nz)
      real  dslay(-ng+1:nx+ng,-ng+1:ny+ng,nz)

      real dtp,dx,dy
      real zlayer(-ng+1:nx+ng,-ng+1:ny+ng,nz)
      real dxi2,dyi2,dzi2(-ng+1:nx+ng,-ng+1:ny+ng,nz)
      real dxi,dyi
!
      real  bconda(21)
      real  bzcond(21)
!
      real  conda(-ng+1:nx+ng,-ng+1:ny+ng,nz), condc(-ng+1:nx+ng,-ng+1:ny+ng,nz)
!
!  polar conductivity
!
      data bconda/5.5e-14,6.e-14,7.e-14,7.8e-14,8.8e-14, &
      1.e-13,1.1e-13,1.3e-13,1.5e-13,1.7e-13,1.9e-13,    &
      2.1e-13,2.3e-13,2.6e-13,2.8e-13,3.1e-13,3.4e-13,   &
      3.8e-13,4.2e-13,4.5e-13,4.9e-13/

!     Potential gradients outside/inside cloud..

      real  ascl(-ng+1:2000)
      real  dccx(-ng+1:2000), dcpx(-ng+1:2000), dapx(-ng+1:2000)
      real  dccy(-ng+1:2000), dcpy(-ng+1:2000), dapy(-ng+1:2000)
      real  dccz(-ng+1:2000), dcpz(-ng+1:2000), dapz(-ng+1:2000)

!     local temp arrays...

!      real, allocatable, dimension(:,:,:): conda,condc,dzi2,dslay,zlayer 
!      real, allocatable, dimension(:,:,:): dz,pot,scw,sci,qc,qi,t0,t1,t4,t5,t6,t7,t8,t9 

      real t0(-ng+1:nx+ng,-ng+1:ny+ng,nz)
      real t1(-ng+1:nx+ng,-ng+1:ny+ng,nz)
      real t4(-ng+1:nx+ng,-ng+1:ny+ng,nz)
      real t5(-ng+1:nx+ng,-ng+1:ny+ng,nz)
      real t6(-ng+1:nx+ng,-ng+1:ny+ng,nz)
      real t7(-ng+1:nx+ng,-ng+1:ny+ng,nz)
      real t8(-ng+1:nx+ng,-ng+1:ny+ng,nz)
      real t9(-ng+1:nx+ng,-ng+1:ny+ng,nz)
      real dz2(-ng+1:nx+ng,-ng+1:ny+ng,nz)
!
      real pot(-ng+1:nx+ng,nz,-ng+1:ny+ng)
      real dz(-ng+1:nx+ng,nz,-ng+1:ny+ng)
      real scw(-ng+1:nx+ng,nz,-ng+1:ny+ng)
      real sci(-ng+1:nx+ng,nz,-ng+1:ny+ng)
      real sctot(-ng+1:nx+ng,nz,-ng+1:ny+ng)
      real dn(-ng+1:nx+ng,nz,-ng+1:ny+ng)
      real qc(-ng+1:nx+ng,nz,-ng+1:ny+ng)
      real qi(-ng+1:nx+ng,nz,-ng+1:ny+ng)

      integer kcldtop(-ng+1:nx+ng,-ng+1:ny+ng)

      integer, save :: imap2=0
!      integer  ghost
      integer  ghosti1,ghostj1,ghosti2,ghostj2
      integer  ghosti,ghostj
      integer km5,km4,kp5,kp4
      double precision dv
      integer :: nz1, kz1, kzb,kze
      
! ########################################################
!
!
!        ghost=min(ime-ite,jme-jte,its-abs(ims),jts-abs(jms))
!
!        ALLOCATE(dz(its-ghost:ite+ghost;kts:kte;jts-ghost:jte+ghost))
!        ALLOCATE(pot(its-ghost:ite+ghost;kts:kte;jts-ghost:jte+ghost))
!        ALLOCATE(sci(its-ghost:ite+ghost;kts:kte;jts-ghost:jte+ghost))
!        ALLOCATE(scw(its-ghost:ite+ghost;kts:kte;jts-ghost:jte+ghost))
!        ALLOCATE(qc(its-ghost:ite+ghost;kts:kte;jts-ghost:jte+ghost))
!        ALLOCATE(qi(its-ghost:ite+ghost;kts:kte;jts-ghost:jte+ghost))
!
!        ALLOCATE(t0(its-ghost:ite+ghost;jts-ghost:jte+ghost;kts:kte))
!        ALLOCATE(t1(its-ghost:ite+ghost;jts-ghost:jte+ghost;kts:kte))
!        ALLOCATE(t4(its-ghost:ite+ghost;jts-ghost:jte+ghost;kts:kte))
!        ALLOCATE(t5(its-ghost:ite+ghost;jts-ghost:jte+ghost;kts:kte))
!        ALLOCATE(t6(its-ghost:ite+ghost;jts-ghost:jte+ghost;kts:kte))
!        ALLOCATE(t7(its-ghost:ite+ghost;jts-ghost:jte+ghost;kts:kte))
!        ALLOCATE(t8(its-ghost:ite+ghost;jts-ghost:jte+ghost;kts:kte))
!        ALLOCATE(t9(its-ghost:ite+ghost;jts-ghost:jte+ghost;kts:kte))
!        ALLOCATE(conda(its-ghost:ite+ghost;jts-ghost:jte+ghost;kts:kte))
!        ALLOCATE(condc(its-ghost:ite+ghost;jts-ghost:jte+ghost;kts:kte))
!        ALLOCATE(dzi2(its-ghost:ite+ghost;jts-ghost:jte+ghost;kts:kte))
!        ALLOCATE(dslay(its-ghost:ite+ghost;jts-ghost:jte+ghost;kts:kte))
!        ALLOCATE(zlayer(its-ghost:ite+ghost;jts-ghost:jte+ghost;kts:kte))

        t0(:,:,:)=0.0 ; t1(:,:,:)=0.0
        t4(:,:,:)=0.0 ; t5(:,:,:)=0.0 ; t6(:,:,:)=0.0 ; t7(:,:,:)=0.0
        t8(:,:,:)=0.0 ; t9(:,:,:)=0.0 ; dslay(:,:,:) = 0.0  
        ecnrml=0.0 ; eanrml =0.0 ;dslay(:,:,:)=0.0 ; zlayer(:,:,:)=0.0
        condc(:,:,:)=1.0e-15 ; conda(:,:,:)=1.0e-14 ; dzi2(:,:,:) = 0.0 

         write(iunit,*) 'kds,kde',kds,kde

!         dz(0,:,:)=dz(1,:,:)
!         dz(nx+1,:,:)=dz(nx,:,:)
!         dz(:,:,0)=dz(:,:,1)
!         dz(:,:,ny+1)=dz(:,:,ny)

        ghosti1=abs(its-ims)
        ghosti2=ime-ite
        ghostj1=abs(jts-jms)
        ghostj2=jme-jte

!        ghosti=3 ! must always be  < =min(ghosti1,ghosti2,ghostj1,ghostj2) 
!        ghostj=3 

        IF (ng.gt.ghosti1.or.ng.gt.ghosti2.or.ng.gt.ghostj1.or.ng.gt.ghostj2) THEN 

        WRITE(iunit,*) 'NG > ghosti1-2,j1-2', ng,ghosti1,ghosti2,ghostj1,ghostj2 
        STOP

        ENDIF

!     NEED TO PLACE ALL WRF ARRAYS IN XYZ BEFORE CALLING THE MPI SEND_RCV STUFF...

      kcldtop(:,:) = 0
      
      do kz = 1,nz
      do jy = 1,ny
      do ix = 1,nx

      t0(ix,jy,kz) = max(qc(ix,kz,jy),0.0) + max(qi(ix,kz,jy),0.0)

!     put pot in t1() because commas mpi uses X,Y,Z ordering

       t1(ix,jy,kz) = pot(ix,kz,jy)
       dz2(ix,jy,kz) = dz(ix,kz,jy)     

! find cloud top
      IF ( iscreen >= 2 ) THEN
         IF ( t0(ix,jy,kz) > qctol ) THEN
           kcldtop(ix,jy) = Max(kz,kcldtop(ix,jy))
         ENDIF
      ENDIF
     
!      if (ndebug) write(iunit,*) 'pot(ix,kz,jy)',pot(ix,kz,jy),ix,kz,jy
      if (ndebug) then
      if(dz(ix,kz,jy).le.0) write(iunit,*) 'dz(ix,kz,jy)',dz(ix,kz,jy),ix,kz,jy
      endif

      end do
      end do
      end do

#ifdef DM_PARALLEL


      IF ( ntasks .gt. 0) THEN 

        IF ( imap2 == 0) THEN 
        imap2=1
        CALL TASK_PROC_MAP()
        ENDIF

! Passing-sendrcv data for gradient computation near tile boundaries 

!       ARRAYS MUST BE OF ORDER XYZ and NOT XZY as in ****  WRF...

!-------T1---------

        IF ( ntasks_x > 1 ) THEN
!        write(0,*) 'send_west, ntask_x = ',ntasks_x
!        CALL sendrecv_westward(nx,ny,nz,1,1,0,1,1,  &
!        CALL sendrecv_westward(nxslm,nyslm,nzslm,ghosti1,0,0,ghosti1,1,  &
        westward_tag = 201
        CALL sendrecv_westward(nx,ny,nz,ng,ng,0,ng,1,  &
     &        w_proc,e_proc,westward_tag,t1)
        eastward_tag = 204
        CALL sendrecv_eastward(nx,ny,nz,ng,ng,0,ng,1,  &
     &        w_proc,e_proc,eastward_tag,t1)
        ENDIF

        IF ( ntasks_y > 1 ) THEN
        southward_tag = 207
        CALL sendrecv_southward(nx,ny,nz,ng,ng,0,ng,1, &
     &        n_proc,s_proc,southward_tag,t1)
        northward_tag = 210
        CALL sendrecv_northward(nx,ny,nz,ng,ng,0,ng,1,  &
     &        n_proc,s_proc,northward_tag,t1)
        ENDIF

!-------T0---------

        IF ( ntasks_x > 1 ) THEN
        westward_tag = 202
        CALL sendrecv_westward(nx,ny,nz,ng,ng,0,ng,1,  &
     &        w_proc,e_proc,westward_tag,t0)
        eastward_tag = 205
        CALL sendrecv_eastward(nx,ny,nz,ng,ng,0,ng,1,  &
     &        w_proc,e_proc,eastward_tag,t0)
        ENDIF

        IF ( ntasks_y > 1 ) THEN
        southward_tag = 208
        CALL sendrecv_southward(nx,ny,nz,ng,ng,0,ng,1, &
     &        n_proc,s_proc,southward_tag,t0)
        northward_tag = 211
        CALL sendrecv_northward(nx,ny,nz,ng,ng,0,ng,1,  &
     &        n_proc,s_proc,northward_tag,t0)
        ENDIF

!-------DZ2---------ARRAY OF PAIN...

        IF ( ntasks_x > 1 ) THEN
        westward_tag = 203
!    should only be done once for dz...the array of PAIN !!!
        CALL sendrecv_westward(nx,ny,nz,ng,ng,0,ng,1,  &
     &        w_proc,e_proc,westward_tag,dz2)

        eastward_tag = 206
        CALL sendrecv_eastward(nx,ny,nz,ng,ng,0,ng,1,  &
     &        w_proc,e_proc,eastward_tag,dz2)
        ENDIF
        
        
        IF ( ntasks_y > 1 ) THEN
!        write(0,*) 'send_north, ntask_y = ',ntasks_y
        southward_tag = 209
        CALL sendrecv_southward(nx,ny,nz,ng,ng,0,ng,1, &
     &        n_proc,s_proc,southward_tag,dz2)

        northward_tag = 212
        CALL sendrecv_northward(nx,ny,nz,ng,ng,0,ng,1,  &
     &        n_proc,s_proc,northward_tag,dz2)
        ENDIF

      ENDIF ! ntasts > 0

#endif


!  variables for gradients calculations
!
      dxi2 = 0.5/dx
      dyi2 = 0.5/dy
      dxi = 1.0/dx
      dyi = 1.0/dy

      do kz = 1,     nz
      do jy = 1,ny
      do ix = 1,nx

!      IF (dz2(ix,jy,kz).gt.1.0) THEN ! because of mysterious small dz values in ghost zones...
        dzi2(ix,jy,kz) = 0.5/dz2(ix,jy,kz)
         if (ndebug) then
          IF (dz2(ix,jy,kz).le.1.or.dzi2(ix,jy,kz).gt.20000.or.dzi2(ix,jy,kz).lt.0.) THEN
          write(iunit,*) 'dzi2',dzi2(ix,jy,kz),dz2(ix,jy,kz),ix,jy,kz
          ENDIF
         endif
       dslay(ix,jy,kz) = (dx*dy*dz2(ix,jy,kz))**(0.333333)
!      ELSE
!      dz2(ix,jy,kz)=500. ! hack because stupid fortran puts in ZEROS as usual.
!      dzi2(ix,jy,kz)= 0.5/dz2(ix,jy,kz)   
!      dslay(ix,jy,kz)= (dx*dy*dz2(ix,jy,kz))**(0.333333)
!      ENDIF

      enddo
      enddo
      ENDDO

! polar conductivity interpolation
!
!      write(iunit,*) 'POLAR CONDUCTIVITY (OBS) AIR AND Z'
      do  kk = 1,21
      bzcond(kk) = (kk-1)*500.
!      write(iunit,*)  bconda(kk),bzcond(kk)
      enddo

!  interpolations to get condc
!
!      write(iunit,*) 'POLAR CONDUCTIVITY: CLD, AIR'
!
      do kz = 1,     nz
      do jy = -ng+1,ny+ng
      do ix = -ng+1,nx+ng
!     zlayer = (kz-kd1)*dz + 0.5*dz
      zlayer(ix,jy,kz) = sum(dz(ix,1:kz,jy)) - 0.5*dz(ix,kz,jy)
       if (ndebug) then
        if (zlayer(ix,jy,kz).gt.20000.or.zlayer(ix,jy,kz).lt.0.) write(iunit,*) 'zlayer', zlayer(ix,jy,kz),dz(ix,kz,jy)
       endif
      enddo
      enddo
      enddo

      do kz = 1, nz-kd1
      do jy = -ng+1,ny+ng
      do ix = -ng+1,nx+ng
!
!  interpolations for z less than 21*dz
!
      do  kk = 1,20
      if ( zlayer(ix,jy,kz) .ge. bzcond(kk) .and. zlayer(ix,jy,kz) .lt. bzcond(kk+kd1) ) then
      conda(ix,jy,kz) = bconda(kk) + (bconda(kk+kd1)-bconda(kk)) &
      * (zlayer(ix,jy,kz)-bzcond(kk)) / (bzcond(kk+kd1)-bzcond(kk))
      end if
      enddo
!
!  extrapolate for z higher than 21*dz
!
      if ( zlayer(ix,jy,kz) .ge. bzcond(21) ) then
      conda(ix,jy,kz) = bconda(20) + (bconda(21)-bconda(20))     &
      * (zlayer(ix,jy,kz)-bzcond(20)) / (bzcond(21)-bzcond(20))
      end if
!
      condc(ix,jy,kz) = (0.1)*conda(ix,jy,kz)
!
!      if (jy.eq.int(ny/2).and.ix.eq.int(nx/2) ) write(iunit,*) kz,condc(int(nx/2),int(ny/2),kz),conda(int(nx/2),int(ny/2),kz)
!
      enddo
      enddo
      enddo

      id1=ng+1 ;jd1=ng+1 ; jd2=ng+1 ; id2=ng+1
!
!  compute gradients components of qi+qc ! bewareof tiles vs domain borders
!
      IF (its.gt.ids+1) id1=0 ! is zero if W boundary is internal
      IF (ite.lt.ide-1) id2=0 
      IF (jts.gt.jds+1) jd1=0 
      IF (jte.lt.jde-1) jd2=0 

      do kz = 1+kd1,nz-kd1*1!2
!      do jy = -ng+1+jd1,ny+ng-jd2
!      do ix = -ng+1+id1,nx+ng-id2
      do jy = -ng+1+jd1+jstag,ny+ng-jd2-jstag
      do ix = -ng+1+id1+istag,nx+ng-id2-istag

      t6(ix,jy,kz)=(t0(ix,jy,kz+kstag)-t0(ix,jy,kz-kstag))*dzi2(ix,jy,kz)
      t5(ix,jy,kz)=(t0(ix,jy+jstag,kz)-t0(ix,jy-jstag,kz))*dyi2
      t4(ix,jy,kz)=(t0(ix+istag,jy,kz)-t0(ix-istag,jy,kz))*dxi2

!      if (ndebug) then
!      write(iunit,*) 'dzi2(ix,jy,kz),dxi2,dyi2',dzi2(ix,jy,kz),dxi2,dyi2
!      if (t4(ix,jy,kz).gt.0.) write(iunit,*) 't4(ix,jy,kz)',t4(ix,jy,kz),ix,jy,kz
!      if (t5(ix,jy,kz).gt.0.) write(iunit,*) 't5(ix,jy,kz)',t5(ix,jy,kz),ix,jy,kz
!      if (t6(ix,jy,kz).gt.0.) write(iunit,*) 't6(ix,jy,kz)',t6(ix,jy,kz),ix,jy,kz
!      endif

      end do
      end do
      end do
!
!  now screening layer computations
!

      IF ( iscreen == 2 ) THEN
        nz1 = 1
        kzb = 1
        kze = 1
      ELSE
        nz1 = nz
        kzb = 1+kd1
        kze = nz-kd1
      ENDIF



      do kz1 = kzb,kze

!      do jy = -ng+1+jd1,ny+ng-jd2
!      do ix = -ng+1+id1,nx+ng-id2
      do jy = -ng+1+jd1+jstag,ny+ng-jd2-jstag
      do ix = -ng+1+id1+istag,nx+ng-id2-istag

      IF ( iscreen == 2 ) THEN
        kz = kcldtop(ix,jy)
      ELSE
        kz = kz1
      ENDIF

        km5 = Max(1,kz-5)
        km4 = Max(1,kz-4)
        kp5 = Min(nz-kd1,kz+5)
        kp4 = Min(nz-kd1,kz+4)

!  zero
!
      ascl(ix)=0.0
      dccx(ix)=0.0
      dcpx(ix)=0.0
      dapx(ix)=0.0
      dccy(ix)=0.0
      dcpy(ix)=0.0
      dapy(ix)=0.0
      dccz(ix)=0.0
      dcpz(ix)=0.0
      dapz(ix)=0.0
!
!  start from point of concern:  do only if ac(lc) + ac(li) > qctol
!
      if ( t0(ix,jy,kz) .gt. qctol ) then
      ascl(ix)=6.
!
!
!  x-direction
!
      IF ( iscreen == 1 .and. (ix .gt. 3 .or. id1 == 0 ) .and. ( ix .lt. nx-3*istag .or. id2 == 0 ) ) THEN

       if ( t0(ix+istag,jy,kz) .lt. qctol .and.     &
          t0(ix+istag+1,jy,kz) .lt. qctol2 .and.   &
          t0(ix+istag+2,jy,kz) .lt. qctol2 .and.   &
          t0(ix-istag,jy,kz) .gt. qctol*fin1 .and.     &
          t0(ix-istag-1,jy,kz) .gt. qctol*fin2 ) then  

! cloud points at ix, ix-1, and ix-2, clear air points at ix+1 and ix+2

         ascl(ix)=ascl(ix)-1.0

! dccx is the average cloud mixing ratio gradient in the x-direction
! between the 'cloudy' point and the 'clear' point

         dccx(ix)=0.5*(t4(ix,jy,kz)+t4(ix-istag,jy,kz))

! dcpx is the electric field on the cloudy side
        IF ( firstord ) THEN
         dcpx(ix) = -dxi*(t1(ix,jy,kz)-t1(ix-istag,jy,kz))
         dapx(ix) = -dxi*(t1(ix+istag+1,jy,kz)-t1(ix+istag,jy,kz))
        ELSE
         dcpx(ix) =-dxi2* (  3.0*t1(ix,jy,kz) -     &
                     4*t1(ix-istag,jy,kz) +      &
                       t1(ix-istag-1,jy,kz) )
! dapx is the electric field on the clear side
         dapx(ix) =-dxi2* (  -3.0*t1(ix+istag,jy,kz) +  &
                      4*t1(ix+istag+1,jy,kz) -       &
                        t1(ix+istag+2,jy,kz) )
        ENDIF
        end if

        if ( t0(ix-istag,jy,kz) .lt. qctol .and.   &
          t0(ix-istag-1,jy,kz) .lt. qctol2 .and.  &
          t0(ix-istag-2,jy,kz) .lt. qctol2 .and.  &
          t0(ix+istag,jy,kz) .gt. qctol*fin1 .and.    &
          t0(ix+istag+1,jy,kz) .gt. qctol*fin2 ) then

        ascl(ix)=ascl(ix)-1.0

        dccx(ix)=0.5*(t4(ix+istag,jy,kz)+t4(ix,jy,kz))

        IF ( firstord ) THEN
         dcpx(ix) = -dxi*(t1(ix+istag,jy,kz)-t1(ix,jy,kz))
         dapx(ix) = -dxi*(t1(ix-istag,jy,kz)-t1(ix-istag-1,jy,kz))
        ELSE
        dcpx(ix) =-dxi2* (  -3.0*t1(ix,jy,kz) +    &
                      4*t1(ix+istag,jy,kz) -     &
                        t1(ix+istag+1,jy,kz) )

        dapx(ix) =-dxi2* (  3.0*t1(ix-istag,jy,kz) - &
                     4*t1(ix-istag-1,jy,kz) +      &
                       t1(ix-istag-2,jy,kz) )
        ENDIF
        
        end if

! check that electric field has same sign on both sides
        IF ( Sign(1.,dcpx(ix)) /= Sign(1.,dapx(ix) ) ) THEN
          dcpx(ix) = 0.0
          dapx(ix) = 0.0
        ENDIF
      ENDIF
!
!  y-direction
!
      IF ( iscreen == 1 .and. ( jy .gt. 3 .or. jd1 == 0 ) .and. ( jy .lt. ny-3*jstag .or. jd2 == 0 ) ) THEN

        if ( t0(ix,jy+jstag,kz) .lt. qctol .and.    &
          t0(ix,jy+jstag+1,kz) .lt. qctol2 .and.   &
          t0(ix,jy+jstag+2,kz) .lt. qctol2 .and.   &
          t0(ix,jy-jstag,kz) .gt. qctol*fin1 .and.     &
          t0(ix,jy-jstag-1,kz) .gt. qctol*fin2 ) then

        ascl(ix)=ascl(ix)-1.0

        dccy(ix)=0.5*(t5(ix,jy,kz)+t5(ix,jy-jstag,kz))

        IF ( firstord ) THEN
         dcpy(ix) = -dyi*(t1(ix,jy,kz)-t1(ix,jy-jstag,kz))
         dapy(ix) = -dyi*(t1(ix,jy+jstag+1,kz)-t1(ix,jy+jstag,kz))
        ELSE

        dcpy(ix) =-dyi2* (  3.0*t1(ix,jy,kz) -   &
                     4*t1(ix,jy-jstag,kz) +    &
                       t1(ix,jy-jstag-1,kz) )

        dapy(ix) =-dyi2* (  -3.0*t1(ix,jy+jstag,kz) +  &
                      4*t1(ix,jy+jstag+1,kz) -       &
                        t1(ix,jy+jstag+2,kz) )
        ENDIF
        
        end if

        if ( t0(ix,jy-jstag,kz) .lt. qctol .and.   &
          t0(ix,jy-jstag-1,kz) .lt. qctol2 .and. &
          t0(ix,jy-jstag-2,kz) .lt. qctol2 .and. &
          t0(ix,jy+jstag,kz) .gt. qctol*fin1 .and.   &
          t0(ix,jy+jstag+1,kz) .gt. qctol*fin2 ) then

        ascl(ix)=ascl(ix)-1.0

        dccy(ix)=0.5*(t5(ix,jy+jstag,kz)+t5(ix,jy,kz))

        IF ( firstord ) THEN
         dcpy(ix) = -dyi*(t1(ix,jy+jstag,kz)-t1(ix,jy,kz))
         dapy(ix) = -dyi*(t1(ix,jy-jstag,kz)-t1(ix,jy-jstag-1,kz))
        ELSE
        dcpy(ix) =-dyi2* (  -3.0*t1(ix,jy,kz) +  &
                      4*t1(ix,jy+jstag,kz) -   &
                        t1(ix,jy+jstag+1,kz) )

        dapy(ix) =-dyi2* (  3.0*t1(ix,jy-jstag,kz) -  &
                     4*t1(ix,jy-jstag-1,kz) +       &
                       t1(ix,jy-jstag-2,kz) )
        ENDIF
        
        end if

! check that electric field has same sign on both sides
        IF ( Sign(1.,dcpy(ix)) /= Sign(1.,dapy(ix) ) ) THEN
          dcpy(ix) = 0.0
          dapy(ix) = 0.0
        ENDIF

       ENDIF
!
!  z-direction  ! not MPI-ed
!
      IF ( kz .gt. 3 .and. kz .lt. nz-3*kstag .or. kcldtop(ix,jy) > 1  ) THEN

        if ( t0(ix,jy,kz+kstag) .lt. qctol .and.     &
     &       t0(ix,jy,kz+kstag+1) .lt. qctol2 .and.   &
     &       t0(ix,jy,kz+kstag+2) .lt. qctol2 .and.   &
     &       t0(ix,jy,kp4) .lt. qctol2 .and.   &
     &       t0(ix,jy,kp5) .lt. qctol2 .and.   &
     &       t0(ix,jy,kz-kstag) .gt. qctol*fin1 .and.     &
     &       t0(ix,jy,kz-kstag-1) .gt. qctol*fin2 ) then

        ascl(ix)=ascl(ix)-1.0

        dccz(ix)=0.5*(t6(ix,jy,kz)+t6(ix,jy,kz-kstag))

        IF ( firstord ) THEN
         dcpz(ix) = -2.0*dzi2(ix,jy,kz-1)*(t1(ix,jy,kz)-t1(ix,jy,kz-kstag))
         dapz(ix) = -2.0*dzi2(ix,jy,kz-1)*(t1(ix,jy,kz+kstag+1)-t1(ix,jy,kz+kstag))
        ELSE
        dcpz(ix) = -dzi2(ix,jy,kz-1)* (  3.0*t1(ix,jy,kz) -    &
                     4*t1(ix,jy,kz-kstag) +           &
                       t1(ix,jy,kz-kstag-1) )

        dapz(ix) = -dzi2(ix,jy,kz+1)* (  -3.0*t1(ix,jy,kz+kstag) +   &
                      4*t1(ix,jy,kz+kstag+1) -               &
                        t1(ix,jy,kz+kstag+2) )
        ENDIF
        
        if (ndebug) write(iunit,*) 'dapz1',ix,jy,kz,dzi2(ix,jy,kz+1),t1(ix,jy,kz+kstag),t1(ix,jy,kz+kstag+1),t1(ix,jy,kz+kstag+2) 
        if (ndebug) write(iunit,*) 'dapz1-dz',ix,jy,kz,dz2(ix,jy,kz+1)

        end if
        if ( t0(ix,jy,kz-kstag) .lt. qctol .and.     &
     &       t0(ix,jy,kz-kstag-1) .lt. qctol2 .and.   &
     &       t0(ix,jy,kz-kstag-2) .lt. qctol2 .and.   &
     &       t0(ix,jy,km4) .lt. qctol2 .and.   &
     &       t0(ix,jy,km5) .lt. qctol2 .and.   &
     &       t0(ix,jy,kz+kstag) .gt. qctol*fin1 .and.     &
     &       t0(ix,jy,kz+kstag+1) .gt. qctol*fin2 ) then

        ascl(ix)=ascl(ix)-1.0

        dccz(ix)=0.5*(t6(ix,jy,kz+kstag)+t6(ix,jy,kz))

        IF ( firstord ) THEN
         dcpz(ix) = -2.0*dzi2(ix,jy,kz-1)*(t1(ix,jy,kz+kstag)-t1(ix,jy,kz))
         dapz(ix) = -2.0*dzi2(ix,jy,kz-1)*(t1(ix,jy,kz-kstag)-t1(ix,jy,kz-kstag-1))
        ELSE

        dcpz(ix) = -dzi2(ix,jy,kz+1)* (  -3.0*t1(ix,jy,kz) +  &
                      4*t1(ix,jy,kz+kstag) -          &
                        t1(ix,jy,kz+kstag+1) )

        dapz(ix) = -dzi2(ix,jy,kz-1)* (  3.0*t1(ix,jy,kz-kstag) -   &
                     4*t1(ix,jy,kz-kstag-1) +              &
                       t1(ix,jy,kz-kstag-2) )
        ENDIF
      if (ndebug) write(iunit,*) 'dapz',ix,jy,kz,dzi2(ix,jy,kz-1),t1(ix,jy,kz-kstag),t1(ix,jy,kz-kstag-1),t1(ix,jy,kz-kstag-2) 
  
        end if
        
! check that electric field has same sign on both sides
        IF ( Sign(1.,dcpz(ix)) /= Sign(1.,dapz(ix) ) ) THEN
          dcpz(ix) = 0.0
          dapz(ix) = 0.0
        ENDIF

      ENDIF
!
!  t7 is the count of points with no cloud

!
      t7(ix,jy,kz) = ascl(ix)

!
!  now screening layer calculation
!
      IF( ascl(ix) .le. epascl ) THEN
!
      dccld  = dccx(ix) + dccy(ix) + dccz(ix)
      dccldn = (dccx(ix)**2 + dccy(ix)**2 + dccz(ix)**2)**(0.50)

        if (ndebug) write(iunit,*) dccldn,dccx(ix),dapx(ix),dccy(ix),dapy(ix),dccz(ix),dapz(ix) 
!
!  only where there is a gradient in qc...
!
       if ( dccldn .ge. 1.e-10 ) then
!
        ecnrml = - ((1.0)/(dccldn))* (dccx(ix)*dcpx(ix)+dccy(ix)*dcpy(ix)+dccz(ix)*dcpz(ix))
        eanrml = - ((1.0)/(dccldn))* (dccx(ix)*dapx(ix)+dccy(ix)*dapy(ix)+dccz(ix)*dapz(ix))

        if (ndebug) write(iunit,*) dccldn,dccx(ix),dapx(ix),dccy(ix),dapy(ix),dccz(ix),dapz(ix) 

        t8(ix,jy,kz) = -(2.0)*eperao*(1.0/dslay(ix,jy,kz))                                    &
        *(eanrml*conda(ix,jy,kz)-ecnrml*condc(ix,jy,kz))/(conda(ix,jy,kz)+condc(ix,jy,kz))    & 
        *((1.0)-exp(-dtp*(conda(ix,jy,kz)+condc(ix,jy,kz))/(2.0*eperao)))

        if (ndebug) write(iunit,*) ix,jy,kz,t8(ix,jy,kz),ecnrml,eanrml,dslay(ix,jy,kz)
!
!  end grad in qc
!
       end if
!
!  end if ascl => epascl
!
      ENDIF
!
!  end if qc() > qctol
!
      end if
!
      end do ! ix = -ng+1+id1,ny+ng-id2*1 
      end do ! jy = -ng+1+jd1,ny+ng-jd2*1 
      end do ! kz =  1+kd1 ,nz-kd1*2 
!
!    now smooth t8
!
!      call smooth(nx,ny,nz,ng,its,ite,jts,jte,ids,ide,jds,jde,   &
!     &           ng+1,ng+1,ng+1,ng+1,1,t8,t9)

!      t9(:,:,:)=0.0

      screenn = 0.0
      screenp = 0.0
!
!  add del rho to scw, sci
!
      do kz = 1+kstag,  nz-1*kstag
      do jy = 1,ny ! -ng+1+jd1,  ny+ng-jd2
      do ix = 1,nx ! -ng+1+id1,  nx+ng-id2
      if ( t0(ix,jy,kz) .ge. qctol .and. t8(ix,jy,kz) /= 0.0 ) then
        dv = dx*dy*dz(ix,kz,jy)
!      if ( t7(ix,jy,kz) .le. epascl ) then
      IF ( t8(ix,jy,kz) > 0.0 ) THEN
        IF ( sctot(ix,kz,jy) > sc_thresh ) THEN
          t8(ix,jy,kz) = 0.0
        ELSE
          IF ( sctot(ix,kz,jy) + t8(ix,jy,kz) > sc_thresh ) THEN
          ! want sctot(ix,kz,jy) + t8new(ix,jy,kz) = 1.e-9
          ! excess pos charge is sctot(ix,kz,jy) + t8(ix,jy,kz) - 1.e-9, subtract this from t8:
          ! t8(ix,jy,kz) =  t8(ix,jy,kz) - (sctot(ix,kz,jy) + t8(ix,jy,kz) - 1.e-9)
            t8(ix,jy,kz) =  - (sctot(ix,kz,jy) - sc_thresh)
          ENDIF
        ENDIF
      ELSE !  ( t8(ix,jy,kz) < 0.0 )
        IF ( sctot(ix,kz,jy) < -sc_thresh ) THEN
         t8(ix,jy,kz) = 0.0 
        ELSE
          IF ( sctot(ix,kz,jy) + t8(ix,jy,kz) < -sc_thresh ) THEN
          ! want sctot(ix,kz,jy) + t8new(ix,jy,kz) = -1.e-9
          ! excess neg charge is sctot(ix,kz,jy) + t8(ix,jy,kz) + 1.e-9, subtract this from t8:
          ! t8(ix,jy,kz) =  t8(ix,jy,kz) - (sctot(ix,kz,jy) + t8(ix,jy,kz) + 1.e-9)
            t8(ix,jy,kz) =  - (sctot(ix,kz,jy) + sc_thresh)
          ENDIF
      
        ENDIF
      ENDIF
      scw(ix,kz,jy) = scw(ix,kz,jy)+t8(ix,jy,kz)*max(qc(ix,kz,jy),0.0)/(dn(ix,kz,jy)*t0(ix,jy,kz))
      sci(ix,kz,jy) = sci(ix,kz,jy)+t8(ix,jy,kz)*max(qi(ix,kz,jy),0.0)/(dn(ix,kz,jy)*t0(ix,jy,kz))
      sctot(ix,kz,jy) = sctot(ix,kz,jy)+t8(ix,jy,kz)
      
        IF ( t8(ix,jy,kz) .gt. 0.) THEN
           screenp = screenp + t8(ix,jy,kz)*dv
        ELSE
           screenn = screenn + t8(ix,jy,kz)*dv
        ENDIF
      end if
      end do
      end do
      end do

! fix: need to sum screenp and screenn in MPI

#ifdef DM_PARALLEL
       mpitotindp(1)  = screenp
       mpitotindp(2)  = screenn

      CALL MPI_Reduce(mpitotindp, mpitotoutdp, 2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, local_communicator, mpi_error_code)


      IF ( mytask == 0 ) THEN ! only print from task 0

       screenp = mpitotoutdp(1)  
       screenn = mpitotoutdp(2)  

#endif


!      if (ndebug) then
      write(iunit,'(a,4(1pe12.5,2x))') 'Screening layer charge (rates) (pos,neg,tot)', &
           screenp/dtp, screenn/dtp, (screenp+screenn)/dtp
      write(iunit,'(a,4(1pe12.5,2x))') 'Screening layer charges (pos,neg,tot)', &
           screenp, screenn, (screenp+screenn)

#ifdef DM_PARALLEL
      ENDIF
#endif
!!'
!      do kz = 1+kstag,  nz-2*kstag
!      do jy = -ng+1+jstag,  ny+ng-1*jstag
!      do ix = -ng+1+istag,  nx+ng-1*istag
!      ix=INT(nx/2)
!      jy=INT(ny/2)

!       IF(abs(screenp).gt.0.or.abs(screenn).gt.0) THEN
!       write(iunit,*) 'ix,jy,kz',ix,jy,kz
!       write(iunit,'(a,3(1pe12.5,2x))') 'Screening layer charge (pos,neg,tot)',screenp,screenn,dx,dy,dz2(ix,jy,kz), &
!           screenp*dx*dy*dz2(ix,jy,kz), screenn*dx*dy*dz2(ix,jy,kz), (screenp+screenn)*dx*dy*dz2(ix,jy,kz)
!!'
!       ENDIF

!      enddo
!      enddo
!      enddo

!      endif

!      DEALLOCATE(dz,pot,scw,sci,qc,qi,t0,t1,t4,t5,t6,t7,t8,t9) 
!      DEALLOCATE(conda,condc,dzi2,dslay,zlayer) 

!  end of screening layer computation
!
      return
      end subroutine screen

!!!!!!!!!!! LOCAL SUBROUTINES !!!!!!!!!!!!
!
      subroutine smooth(nx,ny,nz,ng,its,ite,jts,jte,ids,ide,jds,jde,   &
     &                   id1,id2,jd1,jd2,kd1,a,b)
!
      implicit none
!
      integer   ia,na,ng
      integer   ix,jy,kz
      integer   id1,jd1,id2,jd2
      integer   kd1,its,ite,jts,jte,ids,ide,jds,jde
      integer   nx,ny,nz
!
      real      a(-ng+1:nx+ng,-ng+1:ny+ng,nz)
      real      b(-ng+1:nx+ng,-ng+1:ny+ng,nz)
!
! zero array
!
      b(:,:,:) = a(:,:,:)
!
!  smooth
!
      IF (its.gt.ids+1) id1=0 
      IF (ite.lt.ide-1) id2=0 
      IF (jts.gt.jds+1) jd1=0 
      IF (jte.lt.jde-1) jd2=0 
!      do kz = 1+kd1,nz-1*kd1
!      do jy = -ng+1+jd1,ny+ng-1*jd2
!      do ix = -ng+1+id1,nx+ng-1*id2
!      do kz = 1+kd1,nz-kd1*1!2
!      do jy = -ng+1+jd1,ny+ng-jd2
!      do ix = -ng+1+id1,nx+ng-id2
      do kz = 1+kd1,nz-kd1
      do jy = -ng+1+jd1+1,ny+ng-jd2-1 ! last 1 here is jstag - too lazy to pass it in
      do ix = -ng+1+id1+1,nx+ng-id2-1 ! last 1 here is istag
      a(ix,jy,kz) = b(ix,jy,kz)                                    &
       +1.0/32.0*(b(ix+1,jy,kz)-2.*b(ix,jy,kz)+b(ix-1,jy,kz))  &
       +1.0/32.0*(b(ix,jy+1,kz)-2.*b(ix,jy,kz)+b(ix,jy-1,kz))  &
       +1.0/32.0*(b(ix,jy,kz+1)-2.*b(ix,jy,kz)+b(ix,jy,kz-1))
      end do
      end do
      end do

      return

      end subroutine smooth
!
!
!
      end module module_mp_screen
