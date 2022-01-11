   Module module_commasmpi
   
      INTEGER :: w_proc = -1, e_proc = -1, s_proc = -1, n_proc = -1
      INTEGER :: mpi_setup_flag = -1
   
   CONTAINS

!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE TASKS-PROCS MAP  <<<<<<<<<<<<<<<<<<<<<  !
!
!-------------------------------------------------------------------------------


      SUBROUTINE TASK_PROC_MAP()

#ifdef DM_PARALLEL
        USE module_dm, ONLY : &
         local_communicator, mytask, ntasks, ntasks_x, ntasks_y,mytask_x, mytask_y

          IMPLICIT NONE

          INTEGER :: memory_status

          INCLUDE "mpif.h"

       IF (  mpi_setup_flag < 0  ) THEN
         mpi_setup_flag = 1
       ELSE
         RETURN
       ENDIF

              w_proc = MPI_UNDEFINED 
              e_proc = MPI_UNDEFINED
              s_proc = MPI_UNDEFINED
              n_proc = MPI_UNDEFINED

! 3x3 grid rsl.000x where x=ntasks


!       NTASKS                 mytasks_x           mytasks_y
! ------------------   ------------------   -------------------
! | 6  |  7  |  8  |   | 0  |  1  |  2  |   |  2  |  2  |  2  |
! ------------------   ------------------   -------------------
! | 3  |  4  |  5  |   | 0  |  1  |  2  |   |  1  |  1  |  1  |
! ------------------   ------------------   -------------------
! | 0  |  1  |  2  |   | 0  |  1  |  2  |   |  0  |  0  |  0  |
! ------------------   ------------------   -------------------

!mytask, ntasks, ntasks_x, ntasks_y,mytask_x, mytask_y  0   9   3   3   0   0
!mytask, ntasks, ntasks_x, ntasks_y,mytask_x, mytask_y  1   9   3   3   1   0
!mytask, ntasks, ntasks_x, ntasks_y,mytask_x, mytask_y  2   9   3   3   2   0
!mytask, ntasks, ntasks_x, ntasks_y,mytask_x, mytask_y  3   9   3   3   0   1
!mytask, ntasks, ntasks_x, ntasks_y,mytask_x, mytask_y  4   9   3   3   1   1
!mytask, ntasks, ntasks_x, ntasks_y,mytask_x, mytask_y  5   9   3   3   2   1
!mytask, ntasks, ntasks_x, ntasks_y,mytask_x, mytask_y  6   9   3   3   0   2
!mytask, ntasks, ntasks_x, ntasks_y,mytask_x, mytask_y  7   9   3   3   1   2
!mytask, ntasks, ntasks_x, ntasks_y,mytask_x, mytask_y  8   9   3   3   2   2                      

              if (mytask_x.gt.0.and.mytask_x.lt.ntasks_x-1) THEN   ! in the "middles"

              e_proc = mytask+1
              w_proc = mytask-1

              ENDIF

              if (mytask_y.gt.0.and.mytask_y.lt.ntasks_y-1) THEN   ! in the "middles"

              s_proc = mytask-ntasks_x
              n_proc = mytask+ntasks_x

              ENDIF

              if (mytask_x.eq.0) THEN   ! W boundary
      
              e_proc = mytask+1

              ENDIF

              if (mytask_x.eq.ntasks_x-1) THEN ! E boundary
              
              w_proc = mytask-1
 
              ENDIF

              if (mytask_y.eq.0) THEN  ! S boundary

              n_proc = mytask+ntasks_x

              ENDIF

              if (mytask_y.eq.ntasks_y-1) THEN  ! N boundary

              s_proc = mytask-ntasks_x 
     
              ENDIF


#endif
      END SUBROUTINE TASK_PROC_MAP



!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE SENDRECV_WESTWARD  <<<<<<<<<<<<<<<<<<<<<  !
!
!-------------------------------------------------------------------------------
      SUBROUTINE sendrecv_westward (                                          &
     &               nxt, nyt, nzt, ngx, ngy, ngz, ng, na,                    &
     &               west_neighbor_process, east_neighbor_process,            &
     &               westward_tag, variable)

        USE module_dm, ONLY : &
         local_communicator, mytask, ntasks, ntasks_x, ntasks_y,mytask_x, mytask_y

          IMPLICIT NONE

          INCLUDE "mpif.h"

          INTEGER,INTENT(IN) :: nxt, nyt, nzt, ngx, ngy, ngz, ng, na
          INTEGER,INTENT(IN) :: west_neighbor_process
          INTEGER,INTENT(IN) :: east_neighbor_process
          INTEGER,INTENT(IN) :: westward_tag

          REAL,                                                               &
     &     DIMENSION(-ngx+1:nxt+ngx,-ngy+1:nyt+ngy,-ngz+1:nzt+ngz,na),              &
     &     INTENT(INOUT) :: variable

          INTEGER,PARAMETER :: memory_success = 0

          REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: westward_send_buffer
          REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: westward_recv_buffer

          INTEGER :: westward_send_size, westward_recv_size
          INTEGER :: n, i, j, k, ghost_zone
          INTEGER :: memory_status
          INTEGER :: mpi_error_code, mpi_abort_error_code

          LOGICAL :: has_western_neighbor, has_eastern_neighbor
          INTEGER :: offset

         LOGICAL :: debug_mpi = .false.
         INTEGER :: myproci,nproci
         INTEGER :: my_rank       
         LOGICAL :: xperiodic
         INTEGER :: mpi_status(MPI_Status_size)

        xperiodic = .false.

         my_rank = mytask
!         number_of_processes=ntasks 
         myproci = mytask_x+1 
         nproci  = ntasks_x  
!         my_com=local_communicator     
!         MPI_COMM_WORLD=local_communicator     
!         myprocj = mytask_y+1 
!         nprocj  = ntasks_y  


          offset = 0
          
          IF ( xperiodic .and. myproci == nproci ) THEN
              offset = 1
          ENDIF


          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Enter sendrecv_westward"
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            "nxt=", nxt, ",nyt=", nyt, ",nzt=", nzt, ",ng=", ng
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            "na=", na
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            "west_neighbor_process=", west_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            "east_neighbor_process=", east_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            "westward_tag=", westward_tag
          END IF !! (debug_mpi)

!........ Send process's west side to east side halo of the neighbor on
!........ the west and recv into process's east side halo the west side
!........ of the neighbor on the east.
!........ west process <-- process <-- east process 

          has_western_neighbor =                                              &
     &        (west_neighbor_process /= MPI_UNDEFINED) .AND.                  &
     &        (west_neighbor_process /= my_rank)
          has_eastern_neighbor =                                              &
     &        (east_neighbor_process /= MPI_UNDEFINED) .AND.                  &
     &        (east_neighbor_process /= my_rank)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            "has_western_neighbor=", has_western_neighbor,              &
     &            "has_eastern_neighbor=", has_eastern_neighbor
          END IF !! (debug_mpi)

          ALLOCATE(westward_send_buffer(                                      &
     &                 -ng  + 1 : 0,                                           &
     &                 -ngy + 1 : nyt + ngy,                                    &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "westward_send_buffer of bounds ",                          &
     &            -ng + 1,":",0,",",                                          &
     &            -ngy + 1,":",nyt + ngy, ",",                                  &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(local_communicator, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(westward_recv_buffer(                                      &
     &                 -ng + 1 : 0,                                           &
     &                 -ngy + 1 : nyt + ngy,                                    &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "westward_recv_buffer of bounds ",                          &
     &            -ng + 1,":",0,",",                                          &
     &            -ngy + 1,":",nyt + ngy,",",                                   &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(local_communicator, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          westward_send_buffer = 0.0
          westward_recv_buffer = 0.0

          IF (has_western_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "about to prepare westward_send_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO j = -ngy+1,nyt+ngy

                      DO ghost_zone = -ng+1, 0

                          westward_send_buffer(ghost_zone,j,k,n) =            &
     &                        variable(-ghost_zone+1,j,k,n)

                      END DO !! ghost_zone

                  END DO !! j

              END DO !! k

              END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "done preparing westward_send_buffer"
              END IF !! (debug_mpi)

          END IF !! (has_western_neighbor)

          westward_send_size =                                                &
     &        (ng)           *                                                &
     &        (nyt + ngy * 2) *                                                &
     &        (nzt + ngz * 2) *                                                &
     &        (na)
          westward_recv_size =                                                &
     &        (ng)           *                                                &
     &        (nyt + ngy * 2) *                                                &
     &        (nzt + ngz * 2) *                                                &
     &        (na)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            " westward_send_size=", westward_send_size
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            " westward_recv_size=", westward_recv_size
          END IF !! (debug_mpi)

          IF (has_western_neighbor .or. has_eastern_neighbor) THEN

              IF (.NOT. has_western_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_westward: ",             &
     &                    "about to call MPI_Recv"
                  END IF !! (debug_mpi)

                  CALL MPI_Recv(                                              &
     &                     westward_recv_buffer,                              &
     &                     westward_recv_size,                                &
     &                     MPI_REAL,                                          &
     &                     east_neighbor_process,                             &
     &                     westward_tag,                                      &
     &                     local_communicator,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_westward: ",             &
     &                    "done calling MPI_Recv, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_westward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(local_communicator, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_westward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE IF (.NOT. has_eastern_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_westward: ",             &
     &                    "about to call MPI_Send"
                  END IF !! (debug_mpi)

                  CALL MPI_Send(                                              &
     &                     westward_send_buffer,                              &
     &                     westward_send_size,                                &
     &                     MPI_REAL,                                          &
     &                     west_neighbor_process,                             &
     &                     westward_tag,                                      &
     &                     local_communicator,                                    &
     &                     mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_westward: ",             &
     &                    "done calling MPI_Send, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_westward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(local_communicator, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_westward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE   !! (.NOT. has_eastern_neighbor)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_westward: ",             &
     &                    "about to call MPI_Sendrecv"
                  END IF !! (debug_mpi)


                  CALL MPI_Sendrecv(                                          &
     &                     westward_send_buffer,                              &
     &                     westward_send_size,                                &
     &                     MPI_REAL,                                          &
     &                     west_neighbor_process,                             &
     &                     westward_tag,                                      &
     &                     westward_recv_buffer,                              &
     &                     westward_recv_size,                                &
     &                     MPI_REAL,                                          &
     &                     east_neighbor_process,                             &
     &                     westward_tag,                                      &
     &                     local_communicator,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_westward: ",             &
     &                    "done calling  MPI_Sendrecv, ",                     &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_westward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(local_communicator, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_westward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              END IF !! (.NOT. has_eastern_neighbor)...ELSE

          END IF !! (has_western_neighbor .or. has_eastern_neighbor)


          IF (has_eastern_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "about to extract westward_recv_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO j = -ngy+1,nyt+ngy

                      DO ghost_zone = -ng+1,0

                          variable(nxt-ghost_zone+1-offset,j,k,n) =                  &
     &                        westward_recv_buffer(ghost_zone,j,k,n)

                      END DO !! ghost_zone

                   END DO !! j

               END DO !! k

               END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "done extracting westward_recv_buffer"
              END IF !! (debug_mpi)

          END IF

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Exit sendrecv_westward"
          END IF !! (debug_mpi)
          
          DEALLOCATE( westward_send_buffer )
          DEALLOCATE( westward_recv_buffer )
          

      END SUBROUTINE sendrecv_westward

!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE SENDRECV_EASTWARD  <<<<<<<<<<<<<<<<<<<<<  !
!
!-------------------------------------------------------------------------------

      SUBROUTINE sendrecv_eastward (                                          &
     &               nxt, nyt, nzt, ngx, ngy, ngz, ng, na,                                   &
     &               west_neighbor_process, east_neighbor_process,            &
     &               eastward_tag, variable)

      USE module_dm, ONLY : &
         local_communicator, mytask, ntasks, ntasks_x, ntasks_y,mytask_x, mytask_y

          IMPLICIT NONE

          INCLUDE 'mpif.h'

          INTEGER,INTENT(IN) :: nxt, nyt, nzt, ngx, ngy, ngz, ng, na
          INTEGER,INTENT(IN) :: west_neighbor_process
          INTEGER,INTENT(IN) :: east_neighbor_process
          INTEGER,INTENT(IN) :: eastward_tag

          REAL,                                                               &
     &     DIMENSION(-ngx+1:nxt+ngx,-ngy+1:nyt+ngy,-ngz+1:nzt+ngz,na),              &
     &     INTENT(INOUT) :: variable

          INTEGER,PARAMETER :: memory_success = 0

          REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: eastward_send_buffer
          REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: eastward_recv_buffer

          INTEGER :: eastward_send_size, eastward_recv_size
          INTEGER :: n, i, j, k, ghost_zone
          INTEGER :: memory_status
          INTEGER :: mpi_error_code, mpi_abort_error_code

          LOGICAL :: has_western_neighbor, has_eastern_neighbor
          INTEGER :: offset
         LOGICAL :: debug_mpi=.false.
        INTEGER :: myproci,nproci
        INTEGER :: my_rank       
        LOGICAL :: xperiodic
        INTEGER :: mpi_status(MPI_Status_size)

        xperiodic = .false. 

         my_rank = mytask
!         number_of_processes=ntasks 
         myproci = mytask_x+1 
         nproci  = ntasks_x  
!         my_com=local_communicator     
!         myprocj = mytask_y+1 
!         nprocj  = ntasks_y  
          
          offset = 0
          
          IF ( xperiodic .and. myproci == nproci ) THEN
              offset = 1
          ENDIF
          

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Enter sendrecv_eastward"
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            "nxt=", nxt, ",nyt=", nyt, ",nzt=", nzt, ",ng=", ng
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            "na=", na
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            "west_neighbor_process=", west_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            "east_neighbor_process=", east_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            "eastward_tag=", eastward_tag
          END IF !! (debug_mpi)

!........ Send process's east side to west side halo of the neighbor on
!........ the east and recv into process's west side halo the east side
!........ of the neighbor on the west.
!........ west process --> process --> east process 

          has_western_neighbor =                                              &
     &        (west_neighbor_process /= MPI_UNDEFINED) .AND.                  &
     &        (west_neighbor_process /= my_rank)
          has_eastern_neighbor =                                              &
     &        (east_neighbor_process /= MPI_UNDEFINED) .AND.                  &
     &        (east_neighbor_process /= my_rank)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            "has_western_neighbor=", has_western_neighbor,              &
     &            "has_eastern_neighbor=", has_eastern_neighbor
          END IF !! (debug_mpi)

          ALLOCATE(eastward_send_buffer(                                      &
     &                 -ng + 1 : 0,                                           &
     &                 -ngy + 1 : nyt + ngy,                                    &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "eastward_send_buffer of bounds ",                          &
     &            -ng + 1,":",0,",",                                          &
     &            -ngy + 1,":",nyt + ngy,",",                                   &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(local_communicator, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(eastward_recv_buffer(                                      &
     &                 -ng + 1 : 0,                                           &
     &                 -ngy + 1 : nyt + ngy,                                    &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "eastward_recv_buffer of bounds ",                          &
     &            -ng + 1,":",0,",",                                          &
     &            -ngy + 1,":",nyt + ngy, ",",                                  &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(local_communicator, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          eastward_send_buffer = 0.0
          eastward_recv_buffer = 0.0

          IF (has_eastern_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "about to prepare eastward_send_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO j = -ngy+1,nyt+ngy

                      DO ghost_zone = -ng+1,0

                          eastward_send_buffer(ghost_zone,j,k,n) =            &
                              variable(nxt+ghost_zone-offset,j,k,n)

                      END DO !! ghost_zone

                  END DO !! j

              END DO !! k

              END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "done preparing eastward_send_buffer"
              END IF !! (debug_mpi)

          END IF !! (has_western_neighbor)

          eastward_send_size =                                                &
     &        (ng)           *                                                &
     &        (nyt + ngy * 2) *                                                &
     &        (nzt + ngz * 2) *                                                &
     &        (na)
          eastward_recv_size =                                                &
     &        (ng)           *                                                &
     &        (nyt + ngy * 2) *                                                &
     &        (nzt + ngz * 2) *                                                &
     &        (na)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            " eastward_send_size=", eastward_send_size
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            " eastward_recv_size=", eastward_recv_size
          END IF !! (debug_mpi)

          IF (has_western_neighbor .or. has_eastern_neighbor) THEN

              IF (.NOT. has_eastern_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_eastward: ",             &
     &                    "about to call MPI_Recv"
                  END IF !! (debug_mpi)

                  CALL MPI_Recv(                                              &
     &                     eastward_recv_buffer,                              &
     &                     eastward_recv_size,                                &
     &                     MPI_REAL,                                          &
     &                     west_neighbor_process,                             &
     &                     eastward_tag,                                      &
     &                     local_communicator,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_eastward: ",             &
     &                    "done calling MPI_Recv, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_eastward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(local_communicator, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_eastward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE IF (.NOT. has_western_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_eastward: ",             &
     &                    "about to call MPI_Send"
                  END IF !! (debug_mpi)

                  CALL MPI_Send(                                              &
     &                     eastward_send_buffer,                              &
     &                     eastward_send_size,                                &
     &                     MPI_REAL,                                          &
     &                     east_neighbor_process,                             &
     &                     eastward_tag,                                      &
     &                     local_communicator,                                    &
     &                     mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_eastward: ",             &
     &                    "done calling MPI_Send, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_eastward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(local_communicator, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_eastward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE   !! (.NOT. has_eastern_neighbor)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_eastward: ",             &
     &                    "about to call MPI_Sendrecv"
                  END IF !! (debug_mpi)

                  CALL MPI_Sendrecv(                                          &
     &                     eastward_send_buffer,                              &
     &                     eastward_send_size,                                &
     &                     MPI_REAL,                                          &
     &                     east_neighbor_process,                             &
     &                     eastward_tag,                                      &
     &                     eastward_recv_buffer,                              &
     &                     eastward_recv_size,                                &
     &                     MPI_REAL,                                          &
     &                     west_neighbor_process,                             &
     &                     eastward_tag,                                      &
     &                     local_communicator,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_eastward: ",             &
     &                    "done calling MPI_Sendrecv, ",                      &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_eastward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(local_communicator, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_eastward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              END IF !! (.NOT. has_eastern_neighbor)...ELSE

          END IF !! (has_western_neighbor .or. has_eastern_neighbor)


          IF (has_western_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "about to extract eastward_recv_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO j = -ngy+1,nyt+ngy

                      DO ghost_zone = -ng+1,0

                          variable(ghost_zone,j,k,n) =                        &
                              eastward_recv_buffer(ghost_zone,j,k,n)

                      END DO !! ghost_zone

                   END DO !! j

               END DO !! k

               END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "done extracting eastward_recv_buffer"
              END IF !! (debug_mpi)

          END IF

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Exit sendrecv_eastward"
          END IF !! (debug_mpi)

          DEALLOCATE( eastward_send_buffer )
          DEALLOCATE( eastward_recv_buffer )
          

      END SUBROUTINE sendrecv_eastward

      
!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE SENDRECV_S_WRD  <<<<<<<<<<<<<<<<<<<<<<<  !
!
!-------------------------------------------------------------------------------
      SUBROUTINE sendrecv_southward (                                         &
     &               nxt, nyt, nzt, ngx, ngy, ngz, ng, na,                                   &
     &               north_neighbor_process, south_neighbor_process,          &
     &               southward_tag, variable, voffset)

      USE module_dm, ONLY : &
         local_communicator, mytask, ntasks, ntasks_x, ntasks_y,mytask_x, mytask_y

          IMPLICIT NONE

          INCLUDE 'mpif.h'

          INTEGER,INTENT(IN) :: nxt, nyt, nzt, ngx, ngy, ngz, ng, na
          INTEGER,INTENT(IN) :: north_neighbor_process
          INTEGER,INTENT(IN) :: south_neighbor_process
          INTEGER,INTENT(IN) :: southward_tag

          REAL,                                                               &
     &     DIMENSION(-ngx+1:nxt+ngx,-ngy+1:nyt+ngy,-ngz+1:nzt+ngz,na),              &
     &     INTENT(INOUT) :: variable

          INTEGER,OPTIONAL,INTENT(IN) :: voffset(na)
          
          INTEGER,PARAMETER :: memory_success = 0

          REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: southward_send_buffer
          REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: southward_recv_buffer

          INTEGER :: southward_send_size, southward_recv_size
          INTEGER :: n, i, j, k, ghost_zone
          INTEGER :: memory_status
          INTEGER :: mpi_error_code, mpi_abort_error_code

          LOGICAL :: has_northern_neighbor, has_southern_neighbor
          INTEGER :: offset(na),isoffset(na)
         LOGICAL :: debug_mpi=.false.
        INTEGER :: myprocj,nprocj
        INTEGER :: my_rank       
        LOGICAL :: yperiodic
        INTEGER :: mpi_status(MPI_Status_size)

        yperiodic = .false. 

         my_rank = mytask
!         number_of_processes=ntasks 
!         myproci = mytask_x+1 
!         nproci  = ntasks_x  
!         my_com=local_communicator     
         myprocj = mytask_y+1 
         nprocj  = ntasks_y  
          
          offset(:) = 0
          isoffset(:) = 0
          
          IF ( yperiodic .and. myprocj == nprocj ) THEN
            IF ( present( voffset ) ) THEN
              offset(:) = voffset(:)
            ELSE
              offset(:) = 1
            ENDIF
          ENDIF

          IF ( yperiodic .and. myprocj == 1 ) THEN
            IF ( present( voffset ) ) THEN
              isoffset(:) = voffset(:)
            ENDIF
          ENDIF

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Enter sendrecv_southward"
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            "nxt=", nxt, ",nyt=", nyt, ",nzt=", nzt, ",ng=", ng
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            "na=", na
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            "north_neighbor_process=", north_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            "south_neighbor_process=", south_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            "southward_tag=", southward_tag
          END IF !! (debug_mpi)

!.......  Send process's south side to north side halo of the neighbor on the south
!.......  and recv into process's north side halo the south side of the neighbor 
!.......  on the north.
!.......  north process --> process --> south process 

          has_northern_neighbor =                                             &
     &        (north_neighbor_process /= MPI_UNDEFINED) .AND.                 &
     &        (north_neighbor_process /= my_rank)
          has_southern_neighbor =                                             &
     &        (south_neighbor_process /= MPI_UNDEFINED) .AND.                 &
     &        (south_neighbor_process /= my_rank)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            "has_northern_neighbor=", has_northern_neighbor,            &
     &            "has_southern_neighbor=", has_southern_neighbor
          END IF !! (debug_mpi)

          ALLOCATE(southward_send_buffer(                                     &
     &                 -ngx + 1 : nxt + ngx,                                    &
     &                 -ng + 1 : 0,                                           &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "southward_send_buffer of bounds ",                         &
     &            -ngx + 1,":",nxt + ngx, ",",                                  &
     &            -ng + 1,":",0,",",                                          &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(local_communicator, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(southward_recv_buffer(                                     &
     &                 -ngx + 1 : nxt + ngx,                                    &
     &                 -ng + 1 : 0,                                           &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "southward_recv_buffer of bounds ",                         &
     &            -ngx + 1,":",nxt + ngx,",",                                   &
     &            -ng + 1,":",0,",",                                          &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(local_communicator, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          southward_send_buffer = 0.0
          southward_recv_buffer = 0.0

          IF (has_southern_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "about to prepare southward_send_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO i = -ngx+1,nxt+ngx

                      DO ghost_zone = -ng+1,0

                          southward_send_buffer(i,ghost_zone,k,n) =           &
                              variable(i,-ghost_zone+1+isoffset(na),k,n)

                      END DO !! ghost_zone

                  END DO !! i

              END DO !! k

              END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "done preparing southward_send_buffer"
              END IF !! (debug_mpi)

          END IF !! (has_northern_neighbor)

          southward_send_size =                                               &
     &        (nxt + ngx * 2) *                             &
     &        (ng)           *                             &
     &        (nzt + ngz * 2) *                                                &
     &        (na)
          southward_recv_size =                                               &
     &        (nxt + ngx * 2) *                             &
     &        (ng)           *                             &
     &        (nzt + ngz * 2) *                                                &
     &        (na)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            " southward_send_size=", southward_send_size
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            " southward_recv_size=", southward_recv_size
          END IF !! (debug_mpi)

          IF (has_northern_neighbor .or. has_southern_neighbor) THEN

              IF (.NOT. has_southern_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_southward: ",            &
     &                    "about to call MPI_Recv"
                  END IF !! (debug_mpi)

                  CALL MPI_Recv(                                              &
     &                     southward_recv_buffer,                             &
     &                     southward_recv_size,                               &
     &                     MPI_REAL,                                          &
     &                     north_neighbor_process,                            &
     &                     southward_tag,                                     &
     &                     local_communicator,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_southward: ",            &
     &                    "done calling MPI_Recv, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_southward: ",        &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(local_communicator, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_southward: ",        &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE IF (.NOT. has_northern_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_southward: ",            &
     &                    "about to call MPI_Send"
                  END IF !! (debug_mpi)

                  CALL MPI_Send(                                              &
     &                     southward_send_buffer,                             &
     &                     southward_send_size,                               &
     &                     MPI_REAL,                                          &
     &                     south_neighbor_process,                            &
     &                     southward_tag,                                     &
     &                     local_communicator,                                    &
     &                     mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_southward: ",            &
     &                    "done calling MPI_Send, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_southward: ",        &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(local_communicator, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_southward: ",        &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE   !! (.NOT. has_southern_neighbor)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_southward: ",            &
     &                    "about to call MPI_Sendrecv"
                  END IF !! (debug_mpi)

                  CALL MPI_Sendrecv(                                          &
     &                     southward_send_buffer,                             &
     &                     southward_send_size,                               &
     &                     MPI_REAL,                                          &
     &                     south_neighbor_process,                            &
     &                     southward_tag,                                     &
     &                     southward_recv_buffer,                             &
     &                     southward_recv_size,                               &
     &                     MPI_REAL,                                          &
     &                     north_neighbor_process,                            &
     &                     southward_tag,                                     &
     &                     local_communicator,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_southward: ",            &
     &                    "done calling  MPI_Sendrecv, ",                     &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_southward: ",        &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(local_communicator, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_southward: ",        &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              END IF !! (.NOT. has_southern_neighbor)...ELSE

          END IF !! (has_northern_neighbor .or. has_southern_neighbor)


          IF (has_northern_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "about to extract southward_recv_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO i = -ngx+1,nxt+ngx

                      DO ghost_zone = -ng+1,0

                          variable(i,nyt-ghost_zone+1-offset(na),k,n) =           &
                              southward_recv_buffer(i,ghost_zone,k,n)

                      END DO !! ghost_zone

                   END DO !! i

               END DO !! k

               END DO !! n


              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "done extracting southward_recv_buffer"
              END IF !! (debug_mpi)

          END IF

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Exit sendrecv_southward"
          END IF !! (debug_mpi)

          DEALLOCATE( southward_send_buffer )
          DEALLOCATE( southward_recv_buffer )
          

      END SUBROUTINE sendrecv_southward

!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE SENDRECV_N_WRD  <<<<<<<<<<<<<<<<<<<<<<<  !
!
!-------------------------------------------------------------------------------
      SUBROUTINE sendrecv_northward (                                         &
     &               nxt, nyt, nzt, ngx, ngy, ngz, ng, na,                                   &
     &               north_neighbor_process, south_neighbor_process,          &
     &               northward_tag, variable, voffset )

      USE module_dm, ONLY : &
         local_communicator, mytask, ntasks, ntasks_x, ntasks_y,mytask_x, mytask_y

          IMPLICIT NONE

          INCLUDE 'mpif.h'

          INTEGER,INTENT(IN) :: nxt, nyt, nzt, ngx, ngy, ngz, ng, na
          INTEGER,INTENT(IN) :: north_neighbor_process
          INTEGER,INTENT(IN) :: south_neighbor_process
          INTEGER,INTENT(IN) :: northward_tag

          REAL,                                                               &
     &     DIMENSION(-ngx+1:nxt+ngx,-ngy+1:nyt+ngy,-ngz+1:nzt+ngz,na),              &
     &     INTENT(INOUT) :: variable

          INTEGER,OPTIONAL,INTENT(IN) :: voffset(na)

          INTEGER,PARAMETER :: memory_success = 0

          REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: northward_send_buffer
          REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: northward_recv_buffer

          INTEGER :: northward_send_size, northward_recv_size
          INTEGER :: n, i, j, k, ghost_zone
          INTEGER :: memory_status
          INTEGER :: mpi_error_code, mpi_abort_error_code

          LOGICAL :: has_northern_neighbor, has_southern_neighbor

          INTEGER :: offset(na)
         LOGICAL :: debug_mpi=.false.
        INTEGER :: myprocj,nprocj
        INTEGER :: my_rank       
        LOGICAL :: yperiodic
        INTEGER :: mpi_status(MPI_Status_size)

        yperiodic = .false. 

         my_rank = mytask
!         number_of_processes=ntasks 
!         myproci = mytask_x+1 
!         nproci  = ntasks_x  
!         my_com=local_communicator     
         myprocj = mytask_y+1 
         nprocj  = ntasks_y  
          
          offset(:) = 0
          
          IF ( yperiodic .and. myprocj == nprocj ) THEN
            IF ( present( voffset ) ) THEN
              offset(:) = voffset(:)
            ELSE
              offset(:) = 1
            ENDIF
          ENDIF
          
          

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Enter sendrecv_northward"
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            "nxt=", nxt, ",nyt=", nyt, ",nzt=", nzt, ",ng=", ng
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            "na=", na
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            "south_neighbor_process=", south_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            "north_neighbor_process=", north_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            "northward_tag=", northward_tag
          END IF !! (debug_mpi)

!.......  Send process's north side to south side halo of the neighbor on the north
!.......  and recv into process's south side halo the north side of the neighbor 
!.......  on the south.
!.......  north process <-- process <-- south process 

          has_northern_neighbor =                                             &
     &        (north_neighbor_process /= MPI_UNDEFINED) .AND.                 &
     &        (north_neighbor_process /= my_rank)
          has_southern_neighbor =                                             &
     &        (south_neighbor_process /= MPI_UNDEFINED) .AND.                 &
     &        (south_neighbor_process /= my_rank)


          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            "has_northern_neighbor=", has_northern_neighbor,            &
     &            "has_southern_neighbor=", has_southern_neighbor
          END IF !! (debug_mpi)

          ALLOCATE(northward_send_buffer(                                     &
     &                 -ngx + 1 : nxt + ngx,                                    &
     &                 -ng + 1 : 0,                                           &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "northward_send_buffer of bounds ",                         &
     &            -ngx + 1,":",nxt + ngx,",",                                   &
     &            -ng + 1,":",0,",",                                          &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(local_communicator, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(northward_recv_buffer(                                     &
     &                 -ngx + 1 : nxt + ngx,                                    &
     &                 -ng + 1 : 0,                                           &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "northward_recv_buffer of bounds ",                         &
     &            -ngx + 1,":",nxt + ngx,",",                                   &
     &            -ng + 1,":",0,",",                                          &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(local_communicator, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          northward_send_buffer = 0.0
          northward_recv_buffer = 0.0

          IF (has_northern_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "about to prepare northward_send_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO i = -ngx+1,nxt+ngx

                      DO ghost_zone = -ng+1,0

                          northward_send_buffer(i,ghost_zone,k,n) =           &
                             variable(i,nyt+ghost_zone-offset(na),k,n)

                      END DO !! ghost_zone

                  END DO !! i

              END DO !! k

              END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "done preparing northward_send_buffer"
              END IF !! (debug_mpi)

          END IF !! (has_southern_neighbor)

          northward_send_size =                                               &
     &        (nxt + ngx * 2) *                             &
     &        (ng)           *                             &
     &        (nzt + ngz * 2) *                                                &
     &        (na)
          northward_recv_size =                                               &
     &        (nxt + ngx * 2) *                             &
     &        (ng)           *                             &
     &        (nzt + ngx * 2) *                                                &
     &        (na)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            " northward_send_size=", northward_send_size
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            " northward_recv_size=", northward_recv_size
          END IF !! (debug_mpi)

          IF (has_northern_neighbor .or. has_southern_neighbor) THEN

              IF (.NOT. has_northern_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_northward: ",            &
     &                    "about to call MPI_Recv"
                  END IF !! (debug_mpi)

                  CALL MPI_Recv(                                              &
     &                     northward_recv_buffer,                             &
     &                     northward_recv_size,                               &
     &                     MPI_REAL,                                          &
     &                     south_neighbor_process,                            &
     &                     northward_tag,                                     &
     &                     local_communicator,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_northward: ",            &
     &                    "done calling MPI_Recv, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_northward: ",        &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(local_communicator, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_northward: ",        &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE IF (.NOT. has_southern_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_northward: ",            &
     &                    "about to call MPI_Send"
                  END IF !! (debug_mpi)

                  CALL MPI_Send(                                              &
     &                     northward_send_buffer,                             &
     &                     northward_send_size,                               &
     &                     MPI_REAL,                                          &
     &                     north_neighbor_process,                            &
     &                     northward_tag,                                     &
     &                     local_communicator,                                    &
     &                     mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_northward: ",            &
     &                    "done calling MPI_Send, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_northward: ",        &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(local_communicator, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_northward: ",        &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE   !! (.NOT. has_northern_neighbor)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_northward: ",            &
     &                    "about to call MPI_Sendrecv"
                  END IF !! (debug_mpi)

                  CALL MPI_Sendrecv(                                          &
     &                     northward_send_buffer,                             &
     &                     northward_send_size,                               &
     &                     MPI_REAL,                                          &
     &                     north_neighbor_process,                            &
     &                     northward_tag,                                     &
     &                     northward_recv_buffer,                             &
     &                     northward_recv_size,                               &
     &                     MPI_REAL,                                          &
     &                     south_neighbor_process,                            &
     &                     northward_tag,                                     &
     &                     local_communicator,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_northward: ",            &
     &                    "done calling  MPI_Sendrecv, ",                     &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_northward: ",        &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(local_communicator, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_northward: ",        &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              END IF !! (.NOT. has_northern_neighbor)...ELSE

          END IF !! (has_southern_neighbor .or. has_northern_neighbor)


          IF (has_southern_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "about to extract northward_recv_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO i = -ngx+1,nxt+ngx

                      DO ghost_zone = -ng+1,0

                          variable(i,ghost_zone,k,n) =                        &
                              northward_recv_buffer(i,ghost_zone,k,n)

                      END DO !! ghost_zone

                   END DO !! i

               END DO !! k

               END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "done extracting northward_recv_buffer"
              END IF !! (debug_mpi)

          END IF

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Exit sendrecv_northward"
          END IF !! (debug_mpi)

          DEALLOCATE( northward_send_buffer )
          DEALLOCATE( northward_recv_buffer )
        

      END SUBROUTINE sendrecv_northward
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      END MODULE module_commasmpi
