       IF ( mp_physics == NSSL_2MOM .and. elec_physics > 0 ) THEN

!     its = ips
!     ite = ipe
!     jts = jps
!     jte = jpe
     
        flashtype = 'Time Step'

     if (ipelectmp.gt.0)  then ! 

     elec(:,:,:,:)=0.
     elecmag(:,:,:)=0.

!    add up flashes up to history dump 

!     if (MOD(NINT(curr_secs),history_interval*60).eq.0) then
      IF ( clear_history ) THEN
        clear_history = .false.
!      write(0,*) 'clear arrays'
      light(:,:) = 0.
      lightdis(:,:) = 0.
      lightdens(:,:) = 0.
!      ionsource(:,:,:) = 0.0
      induc(:,:,:) = 0.0
      noninduc(:,:,:) = 0.0
       IF ( idischarge == 1 ) THEN
         lightfod(:,:) = 0.
       ELSEIF ( idischarge >= 2 ) THEN
         flshfedic(:,:) = 0.
         flshfedicp(:,:) = 0.
         flshfedicn(:,:) = 0.
         flshfedcg(:,:) = 0.
         flshfedcgp(:,:) = 0.
         flshfedcgn(:,:) = 0.
         flshi(:,:,:) = 0.
         flshn(:,:,:) = 0.
         flshp(:,:,:) = 0.
       ENDIF
      
     endif
     
     ! this sets the clear flag for the next time step if a history dump is occurring now
     IF ( diag_flag_hist ) THEN
       clear_history = .true.
     ENDIF

       IF ( num_tiles > 1 ) THEN
         write(0,*) 'CANNOT RUN WITH num_tiles > 1 !'
         STOP
       ENDIF

  93  format(A,9(I3,2x))

!       nxdg = ime-ims+1
!       nydg = jme-jms+1
!       nzdg = kme-kms+1

        itssav=its
        itesav=ite
        jtssav=jts
        jtesav=jte

       IF (channel) THEN
         its = max(i_start(1),ids)
         ite = min(i_end(num_tiles),ide-1)
       ELSE
         its = ips ! max(i_start(1),ids+sz)
         ite = min(i_end(num_tiles),ide-1-sz)
         IF ( ime > ide ) ite = ide - 1 
       ENDIF
       jts = jps ! max(j_start(1),jds+sz)
       jte = min(j_end(1),jde-1-sz)
       IF ( jme > jde ) jte = jde - 1

!       write(*,*) 'its,ite,jts,jte = ',its,ite,jts,jte
!       write(*,*) 'ips,ipe,jps,jpe = ',ips,ipe,jps,jpe
!       write(*,*) 'kts,kte, kps,kpe, kms,kme = ',kts,kte, kps,kpe, kms,kme 
       nxslm = ite-its+1 
       nyslm = jte-jts+1
!       nxslm = ipe-ips+1 ! should use "patch" size (ipe) not "tile" size (ite), BUT need to update calls...
!       nyslm = jpe-jps+1
!       IF ( ipe == ide ) nxslm = nxslm - 1
!       IF ( jpe == jde ) nyslm = nyslm - 1

       nzslm = kte-kts+1 ! for now. need to make this deeper with extended domain

       IF ( ims0 == 0 ) THEN
       write(0,*) 'ims,ime,its,ite = ',ims,ime,its,ite
       write(0,*) 'jms,jme,jts,jte = ',jms,jme,jts,jte
       write(0,*) 'ips,ipe,jps,jpe = ',ips,ipe,jps,jpe
       write(0,*) 'kts,kte,kms,kme = ',kts,kte,kms,kme
!       write(0,*) 'nxdg,nydg = ',nxdg,nydg,nzdg
       write(0,*) 'nxslm,nyslm = ',nxslm,nyslm,nzslm
       write(0,*) 'ng = ',its-ims,ime-ite,jts-jms,jme-jte
       ims0 = ims - its + 1
       ime0 = ime - ims + ims0
       jms0 = jms - jts + 1
       jme0 = jme - jms + jms0
       write(0,*) 'ims0,ime0,jms0,jme0 = ',ims0,ime0,jms0,jme0
       write(0,*) 'nbw,nbe,nbs,nbn = ',nbw,nbe,nbs,nbn
       ENDIF

       IF ( mpi_setup_flag < 0 ) THEN
         CALL TASK_PROC_MAP()

         IF ( irealmpitype == mpi_real ) THEN
           write(0,*) 'irealmpitype = ',irealmpitype,' = mpi_real'
         ELSEIF ( irealmpitype == mpi_double_precision ) THEN
           write(0,*) 'irealmpitype = ',mpi_double_precision, ' = mpi_double_precision'
         ELSE
           write(0,*) 'irealmpitype not correctly defined'
           WRITE( wrf_err_message , * ) 'WRF-ELEC: irealmpitype not correctly defined. Stop. '
           CALL wrf_error_fatal ( wrf_err_message )
         ENDIF

       ENDIF

       ngx1 = -ims0 ! ips - ims - 1
       ngx2 = ime - ite
       ngy1 = -jms0 ! jps - jms - 1
       ngy2 = jme - jte
       
!       write(0,*) 'ngx1/2,ngy1/2 = ',ngx1,ngx2,ngy1,ngy2
!       write(0,*) 'dims x/y = ', (nxslm+ngx2 - (-ngx1))+1, (nyslm+ngy2+ngy1+1), ime-ims+1,jme-jms+1
!       write(0,*) 'dims of light1 = ', size(light,1), size(light,2)
       
       ! while its,ite are being used, need to check for north/east boundary
!       IF ( ipe == ide ) ngx2 = ngx2 + 1
!       IF ( jpe == jde ) ngy2 = ngy2 + 1
!       write(*,*) 'arrays x: ',-ngx1,nxslm,ngx2, nxslm+ngx2 + ngx1 + 1, ime-ims+1
!       write(*,*) 'arrays y: ',-ngy1,nyslm,ngy2, nyslm+ngy2 + ngy1 + 1, jme-jms+1
!       write(*,*) 'arrays z: ',nzslm, kme-kms+1
!       write(*,*) 'ims:ime , jms:jme = ',ims,ime , jms,jme
!       write(*,*) 'ips,ipe,jps,jpe = ',ips,ipe,jps,jpe

       nnxs(1) = nxslm
       nnys(1) = nyslm
       nnzs(1) = nzslm

!       write(iunit,*) 'BOXMG: nnxs,nnys,nnzs = ',nnxs(1),nnys(1),nnzs(1)

         nbw = 0
         nbe = 0
         nbs = 0
         nbn = 0
         
         ! Future test of extended lateral boundary -- do nnxs,nnys need to be updated, too?
         ! Not necessarily restricted to WRF array sizes since Boxmg arrays are separately allocated
         ! IF ( ips == 1 ) nbw = ngx1
         ! IF ( jps == 1 ) nbs = ngy1
         ! IF ( ipe == ide ) nbe = ngx2
         ! IF ( jpe == jde ) nbn = ngy2


! Compute fair weather elec field used as first guess:
! fair weather field from Gish 1944 / Mansell et al 2005


!      do k=kts,kte
!      do j=jts-1,jte+1
!      do i=its-1,ite+1
!
!      height2(i,k,j)=sum( dz8w(i,1:k,j) )
!
!     ezfair(i,k,j) =   ezfairo *                           &
!        &      ( efb1*exp(-efa1*height2(i,k,j))             &
!        &        +efb2*exp(-efa2*height2(i,k,j))            &
!        &        +efb3*exp(-efa3*height2(i,k,j)) )
!
!
!     enddo
!     enddo
!
!     potfair(k) = ezfairo *                            &
!        &       ( -efa1*efb1*exp(-efa1*height2(its,k,jts))      &
!        &         -efa2*efb2*exp(-efa2*height2(its,k,jts))      &
!        &         -efa3*efb3*exp(-efa3*height2(its,k,jts)) )
!
!     enddo

!        potfair(its-1,:,:)=potfair(its,:,:)
!        potfair(ite+1,:,:)=potfair(ite,:,:)
!        potfair(:,:,jts-1)=potfair(:,:,jts)
!        potfair(:,:,jte+1)=potfair(:,:,jte)

        potfair(:) = 0.0 ! set to zero for now, since we do not have an ion field

!       ezfair(its-1,:,:)=ezfair(its,:,:)
!       ezfair(ite+1,:,:)=ezfair(ite,:,:)
!       ezfair(:,:,jts-1)=ezfair(:,:,jts)
!       ezfair(:,:,jte+1)=ezfair(:,:,jte)


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
      flashcount = 0
      IF ( idischarge == 2 ) THEN
        initmask(:,:,:) = 0
      ENDIF

!         allocate( iontmp(ims:ime, kms:kme, jms:jme ) )

!         iontmp(:,:,:) = sciona(:,:,:)

           nt_sub = Max(1, Nint(dt/0.75) ) ! approx 1s sub time step

     dt_sub = dt/nt_sub

     DO it_sub = 1,nt_sub
      IF ( wrf_dm_on_monitor() ) THEN
       write(iunit,*) 'START ITERSUB ', it_sub
      ENDIF
      lgtstp = 0
      loccur = 1
      iter=0

         
        IF ( mp_physics == NSSL_2MOM ) THEN

        CALL nssl_2mom_sed_driver(                          &
                     ITIMESTEP=itimestep,                &
                     TH=th,                              &
                     QV=qv_curr,                         &
                     QC=qc_curr,                         &
                     QR=qr_curr,                         &
                     QI=qi_curr,                         &
                     QS=qs_curr,                         &
                     QH=qg_curr,                         &
                     QHL=qh_curr,                        &
                     cn=qnn_curr,                        &
                     CCW=qndrop_curr,                    &
                     CRW=qnr_curr,                       &
                     CCI=qni_curr,                       &
                     CSW=qns_curr,                       &
                     CHW=qng_curr,                       &
                     CHL=qnh_curr,                       &
                     SCR=scr,                            &
                     SCW=scw,                            &
                     SCI=sci,                            &
                     SCS=scs,                            &
                     SCH=sch,                            &
                     SCHL=schl,                          &
                     SCTOT=sctot,                        &
                     NONINDUC = noninduc,                &
                     rscghis_2d = rscghis_2d,            &
                     INDUC = induc,                      & 
                     ELECZ = elecz,                      & 
                     SCIONA = sciona, f_sciona=f_sciona, & 
                     scion = cion, f_scion=f_cnion,      &
                     VHW=qvolg_curr, f_vhw=F_QVOLG,      &
                     VHL=qvolh_curr, f_vhl=F_QVOLH,      &
                     ZRW=qzr_curr,  f_zrw = f_qzr,       &
                     ZHW=qzg_curr,  f_zhw = f_qzg,       &
                     ZHL=qzh_curr,  f_zhl = f_qzh,       &
                     PII=pi_phy,                         &
                     P=p,                                &
                     W=w,                                &
                     DZ=dz8w,                            &
                     dx=dx, dy=dy,                       &
                     DTP=dt_sub,                         &
                     DN=rho,                             &
                     RAINNC   = RAINNC,                  &
                     RAINNCV  = RAINNCV,                 &
                     SNOWNC   = SNOWNC,                  &
                     SNOWNCV  = SNOWNCV,                 &
                     HAILNC   = HAILNC,                  &
                     HAILNCV  = HAILNCV,                 &
                     GRPLNC   = GRAUPELNC,               &
                     GRPLNCV  = GRAUPELNCV,              &
                     SR=SR,                              &
                     diagflag = diagflag,                &
                     ke_diag = ke_diag,                  &
                     ipelectmp = ipelectmp,              &
                     isedonly_in = 0,                    &
                     lastloop = ( it_sub == nt_sub ),    &
                     dbz      = refl_10cm,               &
                     nssl_progn=nssl_progn,              &
                     re_cloud=re_cloud,                  &
                     re_ice=re_ice,                      &
                     re_snow=re_snow,                    &
                     has_reqc=has_reqc,                  & ! ala G. Thompson
                     has_reqi=has_reqi,                  & ! ala G. Thompson
                     has_reqs=has_reqs,                  & ! ala G. Thompson
                     hail_maxk1=hail_maxk1,              &
                     hail_max2d=hail_max2d,              &
                     nwp_diagnostics=config_flags%nwp_diagnostics, &
                  IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde, &
                  IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme, &
                  ITS=itssav,ITE=itesav, JTS=jtssav,JTE=jtesav, KTS=kts,KTE=kte  &
                             )

!       CALL nssl_2mom_sed_driver(                       &
!                    ITIMESTEP=itimestep,                &
!                    TH=th,                              &
!                    QV=qv_curr,                         &
!                    QC=qc_curr,                         &
!                    QR=qr_curr,                         &
!                    QI=qi_curr,                         &
!                    QS=qs_curr,                         &
!                    QH=qg_curr,                         &
!                    QHL=qh_curr,                        &
!!                     cn=qnn_curr,                        &
!                    CCW=qndrop_curr,                    &
!                    CRW=qnr_curr,                       &
!                    CCI=qni_curr,                       &
!                    CSW=qns_curr,                       &
!                    CHW=qng_curr,                       &
!                    CHL=qnh_curr,                       &
!                    SCR=scr,                            &
!                    SCW=scw,                            &
!                    SCI=sci,                            &
!                    SCS=scs,                            &
!                    SCH=sch,                            &
!                    SCHL=schl,                          &
!                    SCTOT=sctot,                        &
!                    NONINDUC = noninduc,                &
!                    rscghis_2d = rscghis_2d,            &
!                    INDUC = induc,                      &
!                    ELECZ = elecz,                      &
!                    SCIONA = sciona, f_sciona=f_sciona, & 
!                    scion = cion, f_scion=f_cnion,      &
!                    VHW=qvolg_curr, f_vhw=F_QVOLG,      &
!                    VHL=qvolh_curr, f_vhl=F_QVOLH,      &
!                    ZRW=qzr_curr,  f_zrw = f_qzr,       &
!                    ZHW=qzg_curr,  f_zhw = f_qzg,       &
!                    ZHL=qzh_curr,  f_zhl = f_qzh,       &
!                    cn=qnn_curr,  f_cn=f_qnn,           &
!                    PII=pi_phy,                         &
!                    P=p,                                &
!                    W=w,                                &
!                    DZ=dz8w,                            &
!                    dx=dx, dy=dy,                       &
!                    DTP=dt_sub,                         &
!                    DN=rho,                             &
!                    RAINNC   = RAINNC,                  &
!                    RAINNCV  = RAINNCV,                 &
!                    SNOWNC   = SNOWNC,                  &
!                    SNOWNCV  = SNOWNCV,                 &
!                    HAILNC   = HAILNC,                  &
!                    HAILNCV  = HAILNCV,                 &
!                    GRPLNC   = GRAUPELNC,               &
!                    GRPLNCV  = GRAUPELNCV,              &
!                    SR=SR,                              &
!                    diagflag = diagflag,                &
!                    ke_diag = ke_diag,                  &
!                    ipelectmp = ipelectmp,              &
!                    lastloop = ( it_sub == nt_sub ),    &
!                    dbz      = refl_10cm,               &
!                    nssl_progn=nssl_progn,              &
!                    re_cloud=re_cloud,                  &
!                    re_ice=re_ice,                      &
!                    re_snow=re_snow,                    &
!                    has_reqc=has_reqc,                  & ! ala G. Thompson
!                    has_reqi=has_reqi,                  & ! ala G. Thompson
!                    has_reqs=has_reqs,                  & ! ala G. Thompson
!                 IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde, &
!                 IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme, &
!                 ITS=itssav,ITE=itesav, JTS=jtssav,JTE=jtesav, KTS=kts,KTE=kte  &
!                                                      )

           ELSEIF (  mp_physics == FAST_KHAIN_LYNN_SHPUND ) THEN


           ENDIF
!
! Start of lightning loop
!
     DO WHILE (loccur > 0 .and. .not. lgtstp == 1 )
 105   iter=iter+1

  94  format(A,9(I3,2x))
      IF ( wrf_dm_on_monitor() ) THEN
      ! write(iunit,*) 'lightning iter', iter
      ENDIF
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


! calculate E field and flip ijk back to ikj - new values for pot and elec are IKJ


         pot(:,:,:) = 0.0
         CALL PUTPHI( Qdg, pot(its-1:ite+1,kts:kte,jts-1:jte+1),elec(its-1:ite+1,kts:kte,jts-1:jte+1,1:4), &
     &               nxslm,nyslm,nzslm,nor,     &
     &               nbw,nbe,nbs,nbn,iunit,         &
     &               NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,   &
     &               iGsdg, jGsdg, kGsdg,dx, dy, dz8w(its-1:ite+1,kts:kte,jts-1:jte+1) )

         energynew = 0.0d0
         
         DO j = jts,jte
         DO k = kts,kte-1
         DO i = its,ite 
             energynew = energynew + 0.5*pot(i,k,j)*sctot(i,k,j)*dz8w(i,k,j)*dx*dy
         ENDDO
         ENDDO
         ENDDO
         
#if ( defined(DM_PARALLEL)  &&   ! defined(STUBMPI) )
         energynew = wrf_dm_sum_real8(energynew)
#endif

       IF ( wrf_dm_on_monitor() ) THEN
        IF ( energyold .gt. 0.0 ) THEN
          energychange = 100.0*(energynew - energyold)/energyold
        ELSE
          energychange = 0.0
        ENDIF
       write(6,'(a,3(1x,1pe13.5),1x,0pf7.2,a)')   &
     &       'Energy info: old,new,difference: ',     &
     &                 energyold, energynew, energynew - energyold,energychange,' %'
       
        IF ( one_flash_at_a_time ) THEN
         write(6,'(a,a,1pe13.5)') trim(flashtype),' energy: ',energynew - energyold
        ENDIF
       ENDIF

        IF ( energyold .gt. 0.0 ) THEN
          energy_diff = -(energynew - energyold)
        ELSE
          energy_diff = 0.
        ENDIF

         energyold = energynew


!      relax top elec field to ezfair
!     recompute elecmag for 4 top points (without weighted average)
!     would need a deeper domain for the potential solve (avoid e.g., mirror of charge effect)-

!      if (irelax) then
!
!      do j=jts-1,jte+1
!      do k=kte-6,kte
!      do i=its-1,ite+1
!
!      if (abs(elec(i,k,j,4)).gt.abs(ezfair(i,k,j))) then
!      
!      do l=1,3
!      elec(i,k,j,l)=( elec(i,k,j,l) + (k-(kte-6))*ezfair(i,k,j) ) / (iter*(1+(k-(kte-6)))) 
!      enddo
!      endif
!
!       elec(i,k,j,4) = &                                                                                        
!    & Sqrt( elec(i,k,j,1)**2 &                                                
!    &     + elec(i,k,j,2)**2 &                                                
!    &     + elec(i,k,j,3)**2 )
!
!      enddo
!      enddo
!      enddo
!    
!      endif

       elecmag(:,:,:)=elec(:,:,:,4)
       loccur = 0
       lgtstp = 0
       
#if DM_PARALLEL
!         dt4 = MPI_Wtime()
!         timeboxmg = timeboxmg + (dt4 - dt1)
#endif
       
       
!       IF ( (mp_physics == NSSL_2MOM .or. mp_physics == NSSL_2MOMCCN ) .and. elec_physics > 0 ) THEN
       IF ( idischarge == 1 ) THEN

!    Ziegler and MacGorman (1994,JAS) simple column discharge within cylinders

!        write(0,*) 'call light1d: ngx1+its:its+nxslm+ngx1,jts+ngy1:nyslm+ngy1+jts = ', &
!                ngx1+its,its+nxslm+ngx1,jts+ngy1,nyslm+ngy1+jts
      CALL  light1d(iter,loccur,   &
      &   nor,nxslm,nyslm,nzslm,ngx1,ngx2,ngy1,ngy2,      &
      &   dx,dy,lightradtmp,ibrkdtmp,ecrittmp,disfrac,dz8w,sctot, &
      &   elecmag, sciona,                               &
      &   nbw,nbe,nbs,nbn,                                        &
      &   iunit, rho,                                             &
      &   light,lightdens,lightdis,etot )
     
        lightfod(:,:) = (dx**2/lightradtmp**2)*light(:,:)
     
       ELSEIF ( idischarge == 2 ) THEN

         allocate( temg(ims:ime, kms:kme, jms:jme ) )
         
         temg(:,:,:) = 0.0
       
         DO j = jts,jte
           DO k = kts,kte
             DO i = its,ite
               temg(i,k,j) = th(i,k,j)*pi_phy(i,k,j)
!               IF ( j == jts + (jte-jts)/2 .and. i == its + (ite-its)/2 ) THEN
!                 write(0,*) 'i,j,k,temg = ',i,j,k,temg(i,k,j)
!               ENDIF
             ENDDO
           ENDDO
         ENDDO

         IF ( .not. allocated( kcldtop ) ) THEN
           allocate( kcldtop(its-1:ite+1,jts-1:jte+1) ) ! (ims:ime, jms:jme ) )
         ENDIF
       
         kcldtop(:,:) = 1
         DO j = jts,jte
           DO k = kts,kte
             DO i = its,ite
               tmp = qi_curr(i,k,j) + qc_curr(i,k,j)
               IF ( tmp >= 1.0e-5 ) THEN
                 kcldtop(i,j) = Max(k, kcldtop(i,j))
               ENDIF
             ENDDO
           ENDDO
         ENDDO

         cnoxtot = 0.0

         lightinit(:,:) = 0
         lightmask(:,:) = 0
         
         inittot = 1
         count = 0
         
         
         DO WHILE ( inittot > 0 .and. lgtstp == 0 )
         
!         write(0,*) 'inittot,lgtstp = ',inittot,lgtstp
         
!         write(0,*) 'call msz'
         CALL  lightmsz(iter,lgtstp,loccur,dt,itimestep,idischarge,nssl_elgtfestop,cnoxtot,   &
      &   nor,nxslm,nyslm,nzslm,ids,ide,jds,jde,kds,kde,its,jts, &
      &   ims0,ime0,jms0,jme0,kms,kme,temg     &
      &  ,dx,dy,lightradtmp,ecrittmp,disfrac,dz8w(its-1:ite+1,kts:kte,jts-1:jte+1),sctot(its-1:ite+1,kts:kte,jts-1:jte+1) &
      &  ,pot(its-1:ite+1,kts:kte,jts-1:jte+1) &
      &  ,elec(its-1:ite+1,kts:kte,jts-1:jte+1,1:4) & ! ,scion(its-1:ite+1,kts:kte,jts-1:jte+1,1:2)                            & ! could also pass elecmag...
      &  ,sciona(its-1:ite+1,kts:kte,jts-1:jte+1)                                                                         &
      &  ,nbw,nbe,nbs,nbn ,nbw0,nbe0,nbs0,nbn0                                                                                                 &
      &  ,iunit, rho(its-1:ite+1,kts:kte,jts-1:jte+1),p(its-1:ite+1,kts:kte,jts-1:jte+1)                                  &
      &  ,light(its-1:ite+1,jts-1:jte+1),lightdens(its-1:ite+1,jts-1:jte+1),etot,lightdis(its-1:ite+1,jts-1:jte+1)        &
      &  ,flshfedic,flshfedicp,flshfedicn,flshfedcg,flshfedcgp,flshfedcgn                                                 &
      &  ,flshi(its-1:ite+1,kts:kte,jts-1:jte+1) &
      &  ,flshn(its-1:ite+1,kts:kte,jts-1:jte+1),flshp(its-1:ite+1,kts:kte,jts-1:jte+1)                                   &
      &  ,lightmask(its-1:ite+1,jts-1:jte+1),lightinit(its-1:ite+1,jts-1:jte+1),initmask(its-1:ite+1,jts-1:jte+1,kts:kte) &
      &  ,nssl_einternal,nssl_tgrnd,nssl_zgrnd &
      &  ,f_clnox,its,jts,flashcount,clnox(its-1:ite+1,kts:kte,jts-1:jte+1),kcldtop )
       
!       write(0,*) 'done msz'
      CALL MPI_Barrier(local_communicator, mpi_error_code)

       ! check lightinit and lightmask for hits in lightinit that are not covered by lightmask
       init = 0
       init2 = 0
       maskcount = 0
       maskcounttot = 0
        outer: DO  i=its-1,ite+1
         DO j=jts-1,jte+1
           ! count how many possible initiation points are *remaining* (not used and not in footprint of previous flash)
           IF ( lightmask(i,j) == 0 .and. lightinit(i,j) >= 1 ) THEN 
            init = init + 1
           ENDIF
           IF ( lightmask(i,j) > 0 ) THEN
             maskcount = maskcount + 1
             IF ( lightinit(i,j) >= 1 ) init2 = init2 + 1
           ENDIF
         ENDDO
       ENDDO outer
       
       inittot = init
       inittot2 = init2
       maskcounttot = maskcount
       

! #if DM_PARALLEL
#if ( defined(DM_PARALLEL)  &&   ! defined(STUBMPI) )
      impitotin(1) = init
      impitotin(2) = maskcount
      impitotin(3) = init2
      n = 3
      CALL MPI_AllReduce(impitotin, impitotout, n, MPI_INTEGER, MPI_SUM, local_communicator, ierr)
!      call wrf_dm_sum_integers(impitotin,impitotout)
      
      inittot = impitotout(1)
      maskcounttot = impitotout(2)
      inittot2 = impitotout(3)
      
#endif
!     write(iunit,*) 'count,init,inittot,maskcount,maskcounttot',count,init,inittot,maskcount,maskcounttot
!     write(0,*) 'count,init,inittot,maskcount,maskcounttot',count,init,inittot,maskcount,maskcounttot

       
!       write(0,*) 'micro_driver: init,inittot,count,lgtstp,loccur = ',init,inittot,count,lgtstp,loccur
!       write(0,*) 'micro_driver: maskcount,maskcounttot = ',maskcount,maskcounttot

       
        cnoxall = cnoxall + cnoxtot
        flashtype = 'No Flash'
        IF ( loccur == 1 ) THEN
           numic = numic + 1
           cnoxic = cnoxic + cnoxtot
           flashtype = 'IC flash'
        ENDIF
        IF ( loccur == 2 ) THEN
           numcgn = numcgn + 1
           cnoxcg = cnoxcg + cnoxtot
           flashtype = 'Neg. CG flash'
        ENDIF
        IF ( loccur == 3 ) THEN
           numcgp = numcgp + 1
           cnoxcg = cnoxcg + cnoxtot
           flashtype = 'Pos. CG flash'
        ENDIF
         
        count = count + 1

       IF ( wrf_dm_on_monitor() ) THEN
         write(iunit,*) 'msz done: iter,count,loccur,lgtstp,remaining inits', iter,count,loccur,lgtstp,inittot
       ENDIF

       IF ( maskcounttot > 0 .and. inittot >= 0 .and. inittot2 > 0 ) THEN
         ! recalculate electric field
         loccur = 1
         lgtstp = 0
         
       ENDIF
       
       ! temporary hack to force one flash per iter
       IF ( one_flash_at_a_time ) THEN
         inittot = 0
       ENDIF
       

        IF ( count >= maxiter ) EXIT
         
       ENDDO ! while (init > 0 )
        
       deallocate( temg ) 
        
        
       ENDIF
!       ENDIF

     IF ( wrf_dm_on_monitor() ) THEN
 !      write(iunit,*) 'discharge done  iter ', iter
 !      write(0,*) 'discharge done  iter ', iter
     ENDIF

!       REDUCE charge everywhere by 5 %  to reduce E if iterations do not
!       succeed in doing so. 
    
       if (iter.ge.maxiter) then

       do k=kts,kte
       do i=its-1,ite+1
       do j=jts-1,jte+1
!       scw(i,k,j)=0.95*scw(i,k,j)
!       scr(i,k,j)=0.95*scr(i,k,j)
!       sch(i,k,j)=0.95*sch(i,k,j)
!       sci(i,k,j)=0.95*sci(i,k,j)
!       scs(i,k,j)=0.95*scs(i,k,j)
!       if (f_schl) schl(i,k,j)=0.95*schl(i,k,j)
        IF ( f_sciona ) sciona(i,k,j) = sciona(i,k,j) - 0.05*sctot(i,k,j)
        sctot(i,k,j)  = sctot(i,k,j)  - 0.05*sctot(i,k,j)
       enddo
       enddo
       enddo

       EXIT

       endif
 
       IF ( wrf_dm_on_monitor() ) THEN
!         if ( loccur > 0 .and. .not. lgtstp == 1 ) write(iunit,*) 'loccur,lgtstp = ',loccur,lgtstp
       ENDIF

       ENDDO ! while

       ENDDO ! substep

! #if DM_PARALLEL
#if ( defined(DM_PARALLEL)  &&   ! defined(STUBMPI) )
      n = 1
      CALL MPI_AllReduce(flashcount, impitotout, n, MPI_INTEGER, MPI_SUM, local_communicator, ierr)
      
      flashcount = impitotout(1)
      
#endif


!       CALL SCREENING LAYER HERE -- Ziegler et al. (1991) 
!       NOTE: Using potential from the end of the previous time step, which is probably OK. 
!           Could move to a point after the potential solve and call only on first iteration 
!           of lightning (iter == 1)
!       Also should update sctot after the screening charge application.

      IF ( iscreen .ge. 1 .and. f_sciona ) THEN 
      
!      write(0,*) 'call screen'

      CALL screen(iscreen,nxslm,nyslm,nzslm,dt,dx,dy,ng &
       &  ,ids,ide,jds,jde,kds,kde,its,ite,jts,jte,ims,ime,kms,kme,jms,jme                 &
       &  ,dz8w(its-3:ite+3,kts:kte,jts-3:jte+3),pot(its-3:ite+3,kts:kte,jts-3:jte+3)      &
       &  ,qc_curr(its-3:ite+3,kts:kte,jts-3:jte+3)                                        &
       &  ,qi_curr(its-3:ite+3,kts:kte,jts-3:jte+3),scw(its-3:ite+3,kts:kte,jts-3:jte+3)   &
       &  ,sci(its-3:ite+3,kts:kte,jts-3:jte+3),sctot(its-3:ite+3,kts:kte,jts-3:jte+3)     &
       &  ,rho(its-3:ite+3,kts:kte,jts-3:jte+3)                                            &
       &  ,sciona(its-3:ite+3,kts:kte,jts-3:jte+3)                                         &
       & ,iunit)
       
!       write(0,*) 'done screen'

      ENDIF

      IF ( idischarge > 1 .and. wrf_dm_on_monitor() ) THEN
!        write (0,*) 'MAX LIGHTNING FLASH RATE in this tile THIS TIME STEP', flashcount ! MAXVAL(light)
!        write (iunit,*) 'MAX LIGHTNING FLASH RATE in this tile THIS TIME STEP', flashcount ! MAXVAL(light)
        write (0,*) 'Number of flashes in this time step:', flashcount ! MAXVAL(light)
        write (iunit,*) 'Number of flashes in this time step:', flashcount ! MAXVAL(light)
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

!        ionsource(ims:ime,kms:kme,jms:jme) =  ionsource(ims:ime,kms:kme,jms:jme) +  &
!                    ( sciona(ims:ime,kms:kme,jms:jme) - iontmp(ims:ime,kms:kme,jms:jme) )

!      deallocate( iontmp )
    ENDIF
