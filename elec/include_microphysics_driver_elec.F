       IF ( (mp_physics == NSSL_2MOM .or. mp_physics == NSSL_2MOMCCN ) .and. elec_physics > 0 ) THEN

     if (ipelectmp.gt.0)  then ! 

     elec(:,:,:,:)=0.
     elecmag(:,:,:)=0.

!    add up flashes up to history dump 

     if (MOD(NINT(curr_secs),history_interval*60).eq.0) then
!      write(0,*) 'clear arrays'
      light(:,:)=0.
      lightdis(:,:)=0.
      lightdens(:,:)=0.
     endif

       IF ( num_tiles > 1 ) THEN
         write(0,*) 'CANNOT RUN WITH num_tiles > 1 !'
         STOP
       ENDIF

  93  format(A,9(I3,2x))

       nxdg = ime-ims+1
       nydg = jme-jms+1
       nzdg = kme-kms+1

       nxslm = ite-its+1
       nyslm = jte-jts+1
       nzslm = kte-kts+1 ! for now. need to make this deeper with extended domain

       nnxs(1) = nxslm
       nnys(1) = nyslm
       nnzs(1) = nzslm

!       write(iunit,*) 'BOXMG: nnxs,nnys,nnzs = ',nnxs(1),nnys(1),nnzs(1)

         nbw = 0
         nbe = 0
         nbs = 0
         nbn = 0


! Compute fair weather elec field used as first guess:
! fair weather field from Gish 1944 / Mansell et al 2005


       do k=kts,kte
       do j=jts-1,jte+1
       do i=its-1,ite+1

       height2(i,k,j)=sum( dz8w(i,1:k,j) )

      ezfair(i,k,j) =   ezfairo *                           &
         &      ( efb1*exp(-efa1*height2(i,k,j))             &
         &        +efb2*exp(-efa2*height2(i,k,j))            &
         &        +efb3*exp(-efa3*height2(i,k,j)) )


      enddo
      enddo

      potfair(k) = ezfairo *                            &
         &       ( -efa1*efb1*exp(-efa1*height2(its,k,jts))      &
         &         -efa2*efb2*exp(-efa2*height2(its,k,jts))      &
         &         -efa3*efb3*exp(-efa3*height2(its,k,jts)) )

      enddo

!        potfair(its-1,:,:)=potfair(its,:,:)
!        potfair(ite+1,:,:)=potfair(ite,:,:)
!        potfair(:,:,jts-1)=potfair(:,:,jts)
!        potfair(:,:,jte+1)=potfair(:,:,jte)

        potfair(:) = 0.0 ! set to zero for now, since we do not have an ion field

        ezfair(its-1,:,:)=ezfair(its,:,:)
        ezfair(ite+1,:,:)=ezfair(ite,:,:)
        ezfair(:,:,jts-1)=ezfair(:,:,jts)
        ezfair(:,:,jte+1)=ezfair(:,:,jte)


      IF ( iboxmgsetup == 0 ) THEN ! do not make more than one call to boxmgsetup otherwise will allocate arrays more than once 

      iboxmgsetup = 1
      CALL mgsetupdg(nxslm,nyslm,nzslm,nor,    &
     &   nnxs,nnys,nnzs,              &
     &   nbw,nbe,nbs,nbn,    &
     &   iunit )

      ENDIF
        
        iter=0

      timeboxmg = 0.0
      
     lgtstp = 0
     loccur = 1

     DO WHILE (loccur > 0 .and. .not. lgtstp == 1 )
 105   iter=iter+1

  94  format(A,9(I3,2x))
       write(iunit,*) 'lightning iter', iter
#if DM_PARALLEL
!         dt1 = MPI_Wtime()
#endif

!        SCTOT used here to compute RHS and initial guess...get SOs Qdg and QFgs

       if ( isetup == 1 ) then       
       BMG_rPARMSdg(id_BMG3_STOP_TOL) =  1.e-10 
       BMG_iPARMSdg(id_BMG3_MAX_ITERSdg) = 12 
       BMG_iPARMS(id_BMG3_CYCLE_CLASS) = 0 ! 0 = FMG for first cycle
       else
       BMG_rPARMSdg(id_BMG3_STOP_TOL) = 1.e-4 
       BMG_iPARMSdg(id_BMG3_MAX_ITERSdg) = 3 
       BMG_iPARMS(id_BMG3_CYCLE_CLASS) = 1 ! 1 = ncycles (skips FMG initial solve)
       endif


      CALL PUTF( SOdg, QFdg, Qdg,                                                 &
     &               sctot(its-1:ite+1,kts:kte,jts-1:jte+1),                      &
     &               potfair(kts:kte),                    &
     &               nxslm,nyslm,nzslm,nor,                                       &
     &               nbw,nbe,nbs,nbn,                                             &
     &               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,                    &
     &               iGsdg, jGsdg, kGsdg, dx, dy, dz8w(its-1:ite+1,kts:kte,jts-1:jte+1), isetup )  
     
       isetup = 0

!      write(0,*) 'sctot max2 = ', Maxval( sctot )
!      write(0,*) 'QFdg max = ', Maxval( QFdg )
!      write(0,*) 'Qdg max = ', Maxval( Qdg )

         CALL BMG3_SymStd_UTILS_zero_times(BMG_rPARMSdg)
         CALL MPI_Barrier(local_communicator, mpi_error_code)

! ==========================================================================
!     >>>>>>>>>>>>>>>>     END: WORKSPACE SETUP   <<<<<<<<<<<<<<<<<<<<<<<<<<
! ==========================================================================

!       i =   BMG_iPARMSdg(id_BMG3_MAX_ITERSdg)

!       if (iter.eq.1) write(0,*) 'PUTF DONE'

!      get POT in the form of Q
         if (.true.) then

         CALL WRF_SymStd_SOLVE_boxmgdg(                                     &
     &             NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg, iGsdg, jGsdg, kGsdg, &
     &             BMG_iPARMSdg, BMG_rPARMSdg, BMG_IOFLAGdg,                      &
     &             Qdg, QFdg, BMG_rWORKdg(BMG_pWORKdg(ip_RESdg)), NFdg, NCbmgdg,  &
     &             SOdg, NSOdg,                                                   &
     &             BMG_rWORKdg(BMG_pWORKdg(ip_SORdg)), NSORdg,                    &
     &             BMG_rWORKdg(BMG_pWORKdg(ip_CIdg)), NCIdg,                      &
     &             BMG_iWORKdg(BMG_pWORKdg(ip_iGdg)), NOGdg, NOGcdg,               &
     &             BMG_iWORK_PLdg, NBMG_iWORK_PLdg,                          &
     &             BMG_rWORK_PLdg, NBMG_rWORK_PLdg,                          &
     &             BMG_iWORK_CSdg, NBMG_iWORK_CSdg,                          &
     &             BMG_rWORK_CSdg, NBMG_rWORK_CSdg,                          &
     &             BMG_iWORKdg(BMG_pWORKdg(ip_MSGdg)), NMSGidg,                  &
     &             pMSGdg, pMSGSOdg,                                         &
     &             BMG_MSG_iGRIDdg, NBMG_MSG_iGRIDdg, BMG_MSG_pGRIDdg,         &
     &             ntasks, BMG_rWORKdg(BMG_pWORKdg(ip_MSG_BUFdg)), NMSGrdg, &
     &             local_communicator                                          &
     &             )

!       BMG_iPARMSdg(id_BMG3_MAX_ITERSdg) = i
!       BMG_rPARMSdg(id_BMG3_STOP_TOL) = bmg_tol

         ENDIF

         CALL MPI_Barrier(local_communicator, mpi_error_code)

!c Turn off future setup:
         BMG_iPARMSdg(id_BMG3_SETUPdg) = 2
         BMG_iPARMS(id_BMG3_CYCLE_CLASS) = bmg_cycletype


! calculate E field and flip ijk back to ikj - spew out new values for pot and elec are IKJ


         pot(:,:,:) = 0.0
         CALL PUTPHI( Qdg, pot(its-1:ite+1,kts:kte,jts-1:jte+1),elec(its-1:ite+1,kts:kte,jts-1:jte+1,1:4), &
     &               nxslm,nyslm,nzslm,nor,     &
     &               nbw,nbe,nbs,nbn,iunit,         &
     &               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,   &
     &               iGsdg, jGsdg, kGsdg,dx, dy, dz8w(its-1:ite+1,kts:kte,jts-1:jte+1) )


!      relax top elec field to ezfair
!     recompute elecmag for 4 top points (without weighted average)
!     would need a deeper domain for the potential solve (avoid e.g., mirror of charge effect)-

       if (irelax) then

       do j=jts-1,jte+1
       do k=kte-6,kte
       do i=its-1,ite+1

       if (abs(elec(i,k,j,4)).gt.abs(ezfair(i,k,j))) then
       
       do l=1,3
       elec(i,k,j,l)=( elec(i,k,j,l) + (k-(kte-6))*ezfair(i,k,j) ) / (iter*(1+(k-(kte-6)))) 
       enddo
       endif

        elec(i,k,j,4) = &                                                                                        
     & Sqrt( elec(i,k,j,1)**2 &                                                
     &     + elec(i,k,j,2)**2 &                                                
     &     + elec(i,k,j,3)**2 )

       enddo
       enddo
       enddo
     
       endif

       elecmag(:,:,:)=elec(:,:,:,4)
       loccur = 0
       lgtstp = 0
       
#if DM_PARALLEL
!         dt4 = MPI_Wtime()
!         timeboxmg = timeboxmg + (dt4 - dt1)
#endif
       
       
!       IF ( (mp_physics == NSSL_2MOM .or. mp_physics == NSSL_2MOMCCN ) .and. elec_physics > 0 ) THEN
       IF (idischarge.eq.1) THEN

!    Ziegler and MacGorman (1994,JAS) simple column discharge within cylinders
       
      CALL  light1d(iter,loccur,   &
      &   nor,nxslm,nyslm,nzslm      &
      &  ,dx,dy,lightradtmp,ibrkdtmp,ecrittmp,disfrac,dz8w(its-1:ite+1,kts:kte,jts-1:jte+1),sctot(its-1:ite+1,kts:kte,jts-1:jte+1) &
      &  ,scw(its-1:ite+1,kts:kte,jts-1:jte+1),scs(its-1:ite+1,kts:kte,jts-1:jte+1),sci(its-1:ite+1,kts:kte,jts-1:jte+1)  &
      &  ,scr(its-1:ite+1,kts:kte,jts-1:jte+1),sch(its-1:ite+1,kts:kte,jts-1:jte+1),schl(its-1:ite+1,kts:kte,jts-1:jte+1) &
      &  ,qc_curr(its-1:ite+1,kts:kte,jts-1:jte+1),qs_curr(its-1:ite+1,kts:kte,jts-1:jte+1)                               &
      &  ,qi_curr(its-1:ite+1,kts:kte,jts-1:jte+1),qr_curr(its-1:ite+1,kts:kte,jts-1:jte+1)                               &
      &  ,qg_curr(its-1:ite+1,kts:kte,jts-1:jte+1),qh_curr(its-1:ite+1,kts:kte,jts-1:jte+1)                               &
      &  ,qnc_curr(its-1:ite+1,kts:kte,jts-1:jte+1),qns_curr(its-1:ite+1,kts:kte,jts-1:jte+1)                             &
      &  ,qni_curr(its-1:ite+1,kts:kte,jts-1:jte+1),qnr_curr(its-1:ite+1,kts:kte,jts-1:jte+1)                             &
      &  ,qng_curr(its-1:ite+1,kts:kte,jts-1:jte+1),qnh_curr(its-1:ite+1,kts:kte,jts-1:jte+1)                             &
      &  ,elec(its-1:ite+1,kts:kte,jts-1:jte+1,4),scion(its-1:ite+1,kts:kte,jts-1:jte+1,1:2)                              &
      &  ,sciona(its-1:ite+1,kts:kte,jts-1:jte+1) &
      &  ,nbw,nbe,nbs,nbn                                                                                                 &
      &  ,iunit, rho(its-1:ite+1,kts:kte,jts-1:jte+1)                           &
      &  ,light(its-1:ite+1,jts-1:jte+1),lightdens(its-1:ite+1,jts-1:jte+1),etot,lightdis(its-1:ite+1,jts-1:jte+1))

       ELSEIF (idischarge .eq. 2) THEN
        
        
        
       ENDIF
!       ENDIF


     write(iunit,*) 'discharge done', iter

!       REDUCE charge everywhere by 5 %  to reduce E if iterations do not
!       succeed in doing so. 
    
       if (iter.ge.maxiter) then

       do k=kts,kte
       do i=its-1,ite+1
       do j=jts-1,jte+1
       scw(i,k,j)=0.95*scw(i,k,j)
       scr(i,k,j)=0.95*scr(i,k,j)
       sch(i,k,j)=0.95*sch(i,k,j)
       sci(i,k,j)=0.95*sci(i,k,j)
       scs(i,k,j)=0.95*scs(i,k,j)
       if (present(qh_curr)) schl(i,k,j)=0.95*schl(i,k,j)
       enddo
       enddo
       enddo

       EXIT

       endif
 
       if ( loccur > 0 .and. .not. lgtstp == 1 ) write(iunit,*) 'loccur,lgtstp = ',loccur,lgtstp

       ENDDO ! while



!       CALL SCREENING LAYER HERE -- Ziegler et al. (1991) 
!       NOTE: Using potential from the end of the previous time step, which is probably OK. 
!           Could move to a point after the potential solve and call only on first iteration 
!           of lightning (iter == 1)
!       Also should update sctot after the screening charge application.

      IF ( iscreen .ge. 1  ) THEN 

      CALL screen(iscreen,nxslm,nyslm,nzslm,dt,dx,dy,ng &
       &  ,ids,ide,jds,jde,kds,kde,its,ite,jts,jte,ims,ime,kms,kme,jms,jme                 &
       &  ,dz8w(its-3:ite+3,kts:kte,jts-3:jte+3),pot(its-3:ite+3,kts:kte,jts-3:jte+3)      &
       &  ,qc_curr(its-3:ite+3,kts:kte,jts-3:jte+3)                                        &
       &  ,qi_curr(its-3:ite+3,kts:kte,jts-3:jte+3),scw(its-3:ite+3,kts:kte,jts-3:jte+3)   &
       &  ,sci(its-3:ite+3,kts:kte,jts-3:jte+3),sctot(its-3:ite+3,kts:kte,jts-3:jte+3)     &
       &  ,rho(its-3:ite+3,kts:kte,jts-3:jte+3)                                            &
       & ,iunit)

      ENDIF

      IF ( idischarge == 1 ) THEN
        write (0,*) 'MAX LIGHTNING FLASH RATE in this tile THIS TIME STEP', MAXVAL(light)
      ENDIF
!        write(0,*) 'times for boxmg = ',timeboxmg
!        write(iunit,*) 'times for boxmg = ',timeboxmg


!      save final values

       do j=jts,jte
       do k=kts,kte
       do i=its,ite
       elecz(i,k,j)=elec(i,k,j,3)  ! for induc charging=dezcomp
        IF ( f_elecx ) THEN ! if horizontal components are being stored
         elecx(i,k,j)=elec(i,k,j,1)
         elecy(i,k,j)=elec(i,k,j,2)
        ENDIF
       enddo
       enddo
       enddo


     
      ELSE   ! ipelectmp.gt.0 - computes charging but not E and discharge- 
     
     elec(:,:,:,:)=0.
     elecz(:,:,:)=0.
     light(:,:)=0.
     lightdens(:,:)=0.
     lightdis(:,:)=0.
     elecmag(:,:,:)=0.

      ENDIF  !  ipelectmp.gt.1

    ENDIF
