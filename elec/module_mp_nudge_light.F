       MODULE module_mp_nudge_light

!      lightning assimilation subroutine - Fierro et al. (2012, MWR)

#ifdef DM_PARALLEL
         USE module_dm, ONLY : &
            local_communicator, mytask, ntasks, ntasks_x, ntasks_y,mytask_x, mytask_y
#endif

        CONTAINS

         SUBROUTINE nudge_light(nx,ny,nz,qg,qv,gridlight,temp,press,dz,dx,dy,dt)

          implicit none
#ifdef DM_PARALLEL
        INCLUDE 'mpif.h'
#endif

         real ::  RH= 0.95 ! minimum RH for nudging: larger RH -> larger grid volume activated for nudging
         real ::  QG_TRESH=0.003 ! threshold nudging on graupel mass to further constraint nudging volume.
         real ::  t_bot=273.15 !top of nudging layer
         real ::  t_up= 253.15 ! bottom of nudging layer
         real ::  MAX_RH= 1.02 
         real ::  BB= 0.25 ! controls the slope of TANH curve
         real ::  CC= 0.25 ! controls spacing of QG curves in eq 1 o
         real ::  LVo=2.501E+06 ! Latent heat of vaporization - assume const
         real ::  Rv=461.5
         real ::  To=273.15
         real ::  ep1=0.6222
         real ::  ep2=611.73 ! in Pa
         real ::  const=2.2 ! do not change

         real AA ! controls how much % is increased above lowest threshold (set by RH) 

          integer i,j,k,nx,ny,nz
          real  gridlight(nx,ny)
          real qg(nx,nz,ny)
          real qv(nx,nz,ny)
          real temp(nx,nz,ny) ! in deg Kelvin
          real press(nx,nz,ny) ! in Pa
          real qvsat ! in kg/kg
          real esw ! in Pa
          real angle
          real nudge_coeff

          real, dimension(nx,nz,ny), intent(in):: dz
          real, intent(IN) :: dx,dy,dt
          
          double precision :: dv, total_vapor, tmp, total_vapor_out
          integer mpi_error_code

         if (maxval(gridlight).gt.0.) then
         write(0,*) 'MAXVAL GRIDLIGHT',maxval(gridlight)
!         write(0,*) 'MAXVAL ARRAYS',maxval(temp),maxval(press),maxval(qg),maxval(qv)
!         write(0,*) 'dx,dy,dz',dx,dy,maxval(dz)
         endif

         AA = MAX_RH - RH
        
        nudge_coeff = 1/dt ! direct insertion as in Fierro et al. (2012).

        total_vapor = 0.
        
        do i=1,nx
           do k=1,nz
              do j=1,ny
              

              dv = dx*dy*dz(i,k,j)

!          for QVS: Clausius Clapeyron Equation (constant Latent Heat of
!          Vaporization Lv) - Taken from Rogers and Yau (1989)

           qvsat=(ep1/press(i,k,j))*(ep2*exp( (LVo/Rv)*((1./To)-(1./(temp(i,k,j)) ) ) ) )

           angle=(1000.*qg(i,k,j))**const

          if (gridlight(i,j).ge.1.) then
          if (qg(i,k,j).lt.QG_TRESH) then
          if (qv(i,k,j).lt.RH*qvsat) then 
            if (temp(i,k,j).lt.t_bot.and.temp(i,k,j).gt.t_up) then 

            tmp = RH*qvsat+AA*qvsat*tanh(BB*gridlight(i,j))*(1.-tanh(CC*angle))
            
            total_vapor = total_vapor + Min(1.,nudge_coeff*dt)*(tmp-qv(i,k,j))*dv

            qv(i,k,j) = qv(i,k,j) + Min(1.,nudge_coeff*dt)*(tmp -qv(i,k,j))
 
            endif
          endif
          endif
          endif

              enddo
            enddo
          enddo

#ifdef DM_PARALLEL
      CALL MPI_Reduce(total_vapor, total_vapor_out, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, local_communicator, mpi_error_code)

      
      
      IF ( mytask == 0 ) THEN
         total_vapor = total_vapor_out
#endif
       write(0,*) 'Total vapor mass added: ',total_vapor
       
#ifdef DM_PARALLEL
       ENDIF
#endif
       
       
        RETURN
        END SUBROUTINE nudge_light
        END MODULE module_mp_nudge_light
