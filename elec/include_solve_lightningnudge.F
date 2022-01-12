
! LIGHTNING QV NUDGING, Fierro et al. (2012, MWR)


!   gridlight3d(:,:,:)=0 ! Set to zero to avoid using lightning data when assimilation is turned off  
   tmpTimeInterval = domain_get_current_time ( grid ) - domain_get_start_time ( grid )
   curr_secs = real_time(tmpTimeInterval)

!    IF (LIGHT_ASSIM_FLAG.AND.grid%id.eq.2) THEN  ! for nested domain 
!    IF (LIGHT_ASSIM_FLAG) THEN    ! assumes one parent domain

      IF (grid%nudge_lightning.eq.1) THEN    

      IF ( curr_secs .EQ. grid%nudge_light_times) light_assim_int_old=0
 
      IF ( curr_secs .GE. grid%nudge_light_times .AND. curr_secs .LT. grid%nudge_light_timee) THEN  ! LIGHT ASSIM interval

       light_assim_int = INT( curr_secs / grid%nudge_light_int ) + 1 

      IF (light_assim_int .ne. light_assim_int_old) THEN

        light_assim_int_old=light_assim_int

!       IF ( light_assim_int .LT. 10 ) THEN
!         WRITE ( char_assim_int, '(i1.1)' ) light_assim_int
!       ELSEIF ( light_assim_int .LT. 100 ) THEN
!         WRITE ( char_assim_int, '(i2.2)' ) light_assim_int
!       ELSEIF ( light_assim_int .LT. 1000 ) THEN
!         WRITE ( char_assim_int, '(i3.3)' ) light_assim_int
!       ELSEIF ( light_assim_int .LT. 10000 ) THEN
!         WRITE ( char_assim_int, '(i4.4)' ) light_assim_int
!       ELSE
!         WRITE ( char_assim_int, '(i5.5)' ) light_assim_int
         WRITE ( char_assim_int, '(i5.5)' ) light_assim_int-(INT( grid%nudge_light_times/grid%nudge_light_int ))
!       ENDIF

       WRITE ( dom, '(i1.1)' ) grid%id

       print *,'filename= ',' light.d0',trim(dom),'.out.',trim(char_assim_int),'.txt'

       OPEN ( UNIT=14, FILE=trim(grid%path_to_files)//'light.d0'//trim(dom)//'.out.'//trim(char_assim_int)//'.txt',  &
         FORM='formatted', STATUS='old' )
 
       DO i = ids, ide
       DO j = jds, jde
         READ ( 14, * ) ltg_val
      
         IF ( i .GE. ims .AND. i .LE. ime .AND. j .GE. jms .AND. j .LE. jme ) THEN
           grid%ltg_dat(i,j) = ltg_val
!           DO k = kms, kme
!             gridlight3d(i,k,j) =  grid%ltg_dat(i,j)
!           ENDDO
         ENDIF
       ENDDO
       ENDDO

       CALL wrf_debug ( 0 , ' DONE reading lightning data ' )
       print *,' light.d0',trim(dom),'.out.',trim(char_assim_int),'.txt'
       print *,'curr_secs ',curr_secs

       CLOSE ( 14 )

      ENDIF  ! light_assim_int .ne. light_assim_int_old

   ENDIF ! curr_secs if-test

!   ELSE
!
!   gridlight3d(:,:,:)=0 ! Set to zero to avoid using lightning data when assimilation is turned off  
!
   ENDIF  ! LIGHT NUDGING 


