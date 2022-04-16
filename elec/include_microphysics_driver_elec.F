       IF ( (mp_physics == NSSL_2MOM .or. mp_physics == NSSL_2MOMG .or.  &
             mp_physics == NSSL_2MOMCCN ) .and. elec_physics > 0 ) THEN

!     its = ips
!     ite = ipe
!     jts = jps
!     jte = jpe
     
     if (ipelectmp.gt.0)  then ! 

     elec(:,:,:,:)=0.
     elecmag(:,:,:)=0.

!    add up flashes up to history dump 

     if (MOD(NINT(curr_secs),history_interval*60).eq.0) then
!      write(0,*) 'clear arrays'
      light(:,:) = 0.
      lightdis(:,:) = 0.
      lightdens(:,:) = 0.
       IF ( idischarge == 1 ) THEN
         lightfod(:,:) = 0.
       ELSEIF ( idischarge == 2 ) THEN
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

       IF ( num_tiles > 1 ) THEN
         write(0,*) 'CANNOT RUN WITH num_tiles > 1 !'
         STOP
       ENDIF

  93  format(A,9(I3,2x))

!       nxdg = ime-ims+1
!       nydg = jme-jms+1
!       nzdg = kme-kms+1

       IF (channel) THEN
         its = max(i_start(1),ids)
         ite = min(i_end(1),ide-1)
       ELSE
         its = max(i_start(1),ids+sz)
         ite = min(i_end(1),ide-1-sz)
       ENDIF
       jts = max(j_start(1),jds+sz)
       jte = min(j_end(1),jde-1-sz)

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
!       write(0,*) 'nxdg,nydg = ',nxdg,nydg,nzdg
       write(0,*) 'nxslm,nyslm = ',nxslm,nyslm,nzslm
       write(0,*) 'ng = ',its-ims,ime-ite,jts-jms,jme-jte
       ims0 = ims - its + 1
       ime0 = ime - ims + ims0
       jms0 = jms - jts + 1
       jme0 = jme - jms + jms0
       write(0,*) 'ims0,ime0,jms0,jme0 = ',ims0,ime0,jms0,jme0
       ENDIF


       ngx1 = ips - ims - 1
       ngx2 = ime - ipe
       ngy1 = jps - jms - 1
       ngy2 = jme - jpe
       
       ! while its,ite are being used, need to check for north/east boundary
       IF ( ipe == ide ) ngx2 = ngx2 + 1
       IF ( jpe == jde ) ngy2 = ngy2 + 1
!       write(*,*) 'arrays x: ',-ngx1,nxslm,ngx2, nxslm+ngx2 + ngx1 + 1, ime-ims+1
!       write(*,*) 'arrays y: ',-ngy1,nyslm,ngy2, nyslm+ngy2 + ngy1 + 1, jme-jms+1
!       write(*,*) 'arrays z: ',nzslm, kme-kms+1

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

      CALL  light1d(iter,loccur,   &
      &   nor,nxslm,nyslm,nzslm,ngx1,ngx2,ngy1,ngy2      &
      &  ,dx,dy,lightradtmp,ibrkdtmp,ecrittmp,disfrac,dz8w,sctot &
      &  ,elecmag                               &
      &  ,sciona &
      &  ,nbw,nbe,nbs,nbn                                                                                                 &
      &  ,iunit, rho                           &
      &  ,light,lightdens,etot,lightdis )
     
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
      &  ,scw(its-1:ite+1,kts:kte,jts-1:jte+1),scs(its-1:ite+1,kts:kte,jts-1:jte+1),sci(its-1:ite+1,kts:kte,jts-1:jte+1)  &
      &  ,scr(its-1:ite+1,kts:kte,jts-1:jte+1),sch(its-1:ite+1,kts:kte,jts-1:jte+1),schl(its-1:ite+1,kts:kte,jts-1:jte+1) &
      &  ,qc_curr(its-1:ite+1,kts:kte,jts-1:jte+1),qs_curr(its-1:ite+1,kts:kte,jts-1:jte+1)                               &
      &  ,qi_curr(its-1:ite+1,kts:kte,jts-1:jte+1),qr_curr(its-1:ite+1,kts:kte,jts-1:jte+1)                               &
      &  ,qg_curr(its-1:ite+1,kts:kte,jts-1:jte+1),qh_curr(its-1:ite+1,kts:kte,jts-1:jte+1)                               &
      &  ,qnc_curr(its-1:ite+1,kts:kte,jts-1:jte+1),qns_curr(its-1:ite+1,kts:kte,jts-1:jte+1)                             &
      &  ,qni_curr(its-1:ite+1,kts:kte,jts-1:jte+1),qnr_curr(its-1:ite+1,kts:kte,jts-1:jte+1)                             &
      &  ,qng_curr(its-1:ite+1,kts:kte,jts-1:jte+1),qnh_curr(its-1:ite+1,kts:kte,jts-1:jte+1)                             &
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
       maskcount = 0
       maskcounttot = 0
        outer: DO  i=its-1,ite+1
         DO j=jts-1,jte+1
           ! count how many possible initiation points are *remaining* (not used and not in footprint of previous flash)
           IF ( lightmask(i,j) == 0 .and. lightinit(i,j) >= 1 ) THEN 
            init = init + 1
           ENDIF
           IF ( lightmask(i,j) > 0 ) maskcount = maskcount + 1
         ENDDO
       ENDDO outer
       
       inittot = init
       maskcounttot = maskcount


! #if DM_PARALLEL
#if ( defined(DM_PARALLEL)  &&   ! defined(STUBMPI) )
      impitotin(1) = init
      impitotin(2) = maskcount
      n = 2
      CALL MPI_AllReduce(impitotin, impitotout, n, MPI_INTEGER, MPI_SUM, local_communicator, ierr)
!      call wrf_dm_sum_integers(impitotin,impitotout)
      
      inittot = impitotout(1)
      maskcounttot = impitotout(2)
      
#endif
!     write(iunit,*) 'count,init,inittot,maskcount,maskcounttot',count,init,inittot,maskcount,maskcounttot
!     write(0,*) 'count,init,inittot,maskcount,maskcounttot',count,init,inittot,maskcount,maskcounttot

       
!       write(0,*) 'micro_driver: init,inittot,count,lgtstp,loccur = ',init,inittot,count,lgtstp,loccur
!       write(0,*) 'micro_driver: maskcount,maskcounttot = ',maskcount,maskcounttot

       
        cnoxall = cnoxall + cnoxtot
        IF ( loccur == 1 ) THEN
           numic = numic + 1
           cnoxic = cnoxic + cnoxtot
        ENDIF
        IF ( loccur == 2 ) THEN
           numcgn = numcgn + 1
           cnoxcg = cnoxcg + cnoxtot
        ENDIF
        IF ( loccur == 3 ) THEN
           numcgp = numcgp + 1
           cnoxcg = cnoxcg + cnoxtot
        ENDIF
         
        count = count + 1

       IF ( wrf_dm_on_monitor() ) THEN
         write(iunit,*) 'msz done: iter,count,loccur,lgtstp,remaining inits', iter,count,loccur,lgtstp,inittot
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
       scw(i,k,j)=0.95*scw(i,k,j)
       scr(i,k,j)=0.95*scr(i,k,j)
       sch(i,k,j)=0.95*sch(i,k,j)
       sci(i,k,j)=0.95*sci(i,k,j)
       scs(i,k,j)=0.95*scs(i,k,j)
       if (f_schl) schl(i,k,j)=0.95*schl(i,k,j)
       enddo
       enddo
       enddo

       EXIT

       endif
 
       IF ( wrf_dm_on_monitor() ) THEN
!         if ( loccur > 0 .and. .not. lgtstp == 1 ) write(iunit,*) 'loccur,lgtstp = ',loccur,lgtstp
       ENDIF

       ENDDO ! while


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

      IF ( iscreen .ge. 1  ) THEN 
      
!      write(0,*) 'call screen'

      CALL screen(iscreen,nxslm,nyslm,nzslm,dt,dx,dy,ng &
       &  ,ids,ide,jds,jde,kds,kde,its,ite,jts,jte,ims,ime,kms,kme,jms,jme                 &
       &  ,dz8w(its-3:ite+3,kts:kte,jts-3:jte+3),pot(its-3:ite+3,kts:kte,jts-3:jte+3)      &
       &  ,qc_curr(its-3:ite+3,kts:kte,jts-3:jte+3)                                        &
       &  ,qi_curr(its-3:ite+3,kts:kte,jts-3:jte+3),scw(its-3:ite+3,kts:kte,jts-3:jte+3)   &
       &  ,sci(its-3:ite+3,kts:kte,jts-3:jte+3),sctot(its-3:ite+3,kts:kte,jts-3:jte+3)     &
       &  ,rho(its-3:ite+3,kts:kte,jts-3:jte+3)                                            &
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

    ENDIF
