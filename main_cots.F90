
!!!=== ver 2015/06/30   Copyright (c) 2013-2015 Takashi NAKAMURA  =====

#include "cppdefs.h"


      PROGRAM main_cots
! **********************************************************************
! *                                                                    *
! *   Main program of cots                                             *
! *                                                                    *
! **********************************************************************
!
      USE netcdf
      USE mod_netcdf
      USE mod_cots
      USE mod_p_coral
      
      implicit none
      

      real(8), parameter :: dt = 1.0d-1      ! time step (days)
! -------------------------------------------------------------------------
      integer, parameter :: Syear  = 2000   ! Starting year
      integer, parameter :: Smonth = 6      ! Starting month
      integer, parameter :: Sday   = 1      ! Starting day
     ! NetCDF file     
      character(len=*), parameter :: GRID_FILE = "input/Yaeyama2_grd_v9b.nc"
      character(len=*), parameter :: INI_FILE  = "input/cots_ini.nc"
      character(len=*), parameter :: IN_FILE1  = "input/cots_larvae_his.nc"
      character(len=*), parameter :: IN_FILE2  = "input/coral_larvae_his.nc"
      character(len=*), parameter :: OUT_FILE  = "D:/COTS_model/output/cots_his.nc"
! -------------------------------------------------------------------------
      integer, parameter :: Nburst  = 20
      integer, parameter :: N  = 1
      integer, parameter :: ng  = 1
      

      real(8), allocatable :: h(:,:)              ! depth (meter)
      real(8), allocatable :: p_sand(:,:)         ! sand coverage (0-1)
      real(8), allocatable :: cots_sum(:,:)       ! density of COTS on all life stages (indiv. m-2)
      
      real(8), allocatable :: LARcots(:,:,:)        ! COTS larval recruitment rate (individual m-3)
      real(8), allocatable :: LARcoral(:,:,:)       ! Coral larval recruitment rate (individual m-3)

      real(8), allocatable :: LARcots_time(:)     ! 
      real(8), allocatable :: LARcoral_time(:)    ! 
     
      real(8) :: cots_time(1)            ! COTS model time

      real(8) :: dx  = 300.d0  ! X Grid size (m)
      real(8) :: dy  = 300.d0  ! Y Grid size (m)   
!      real(8) :: dz(Im,Jm,N)  = 1.d0   
!      real(8) :: C(2,0:Im,0:Jm) = 0.d0
!      real(8) :: dC_dt(2,0:Im,0:Jm) = 0.d0

      integer :: Im, Jm
      integer :: i,j,k,id, Nid
      integer :: istep, iburst, iprint
      integer :: n_larcots, n_larcoral
      integer :: IntPrint
      real(8) :: time    !(day)
      
      character(33) :: TIME_ATT  = "days since 1992-01-01 00:00:00"
      
      integer :: N_xi_rho, N_eta_rho
      integer :: N_LARcots_time, N_LARcoral_time
      integer :: ncid, var_id
      integer :: start1D(1), count1D(1)
      integer :: start3D(3), count3D(3)
      integer :: start4D(4), count4D(4)
      character(4) :: YYYY
      character(2) :: MM
      character(2) :: DD

      
      write (YYYY, "(I4.4)") Syear
      write (MM, "(I2.2)") Smonth
      write (DD, "(I2.2)") Sday
      
!---- Modify time-unit description ---------------------------------
      
      TIME_ATT(12:15)=YYYY
      TIME_ATT(17:18)=MM
      TIME_ATT(20:21)=DD

!---- Read ROMS grid netCDF file --------------------------------

      write(*,*) "OPEN: ", GRID_FILE
      
      ! Open NetCDF grid file
      call check( nf90_open(GRID_FILE, nf90_nowrite, ncid) )
      ! Get dimension data
      call get_dimension(ncid, 'xi_rho',  N_xi_rho)
      call get_dimension(ncid, 'eta_rho', N_eta_rho)
      
      allocate(h(N_xi_rho, N_eta_rho))
      allocate(p_sand (N_xi_rho, N_eta_rho))
      allocate(cots_sum(N_xi_rho, N_eta_rho))
      
      Im = N_xi_rho -1
      Jm = N_eta_rho-1
      
!----- Initialization -------------------------
      CALL initialize_cots(1, 1, 0, Im, 0, Jm)
      CALL initialize_p_coral(1, 1, 0, Im, 0, Jm)
      
      ! Get variable id
      call check( nf90_inq_varid(ncid, 'h', var_id) ) 
      call check( nf90_get_var(ncid, var_id, h) )
      call check( nf90_inq_varid(ncid, 'p_coral', var_id) ) 
      call check( nf90_get_var(ncid, var_id, P_CORAL(ng)%cover(:,:)) )
      call check( nf90_inq_varid(ncid, 'p_sand', var_id) ) 
      call check( nf90_get_var(ncid, var_id, p_sand ) )
      ! Close NetCDF file
      call check( nf90_close(ncid) )
      
      write(*,*) "CLOSE: ", GRID_FILE
      
      
      cots_time(1)=0.0d0
      IntPrint = int(30.0d0/dt)
      
!---- Read COTS initial condition netCDF file --------------------------------

      write(*,*) "OPEN: ", INI_FILE
      
      ! Open NetCDF grid file
      call check( nf90_open(INI_FILE, nf90_nowrite, ncid) )
      
      ! Get variable id
      call check( nf90_inq_varid(ncid, 'cots1', var_id) ) 
      call check( nf90_get_var(ncid, var_id, COTS(ng)%dens(1,:,:)) )
      call check( nf90_inq_varid(ncid, 'cots2', var_id) ) 
      call check( nf90_get_var(ncid, var_id, COTS(ng)%dens(2,:,:)) )
      call check( nf90_inq_varid(ncid, 'cots3', var_id) ) 
      call check( nf90_get_var(ncid, var_id, COTS(ng)%dens(3,:,:)) )
      call check( nf90_inq_varid(ncid, 'cots4', var_id) ) 
      call check( nf90_get_var(ncid, var_id, COTS(ng)%dens(4,:,:)) )
      call check( nf90_inq_varid(ncid, 'cots5', var_id) ) 
      call check( nf90_get_var(ncid, var_id, COTS(ng)%dens(5,:,:)) )
      ! Close NetCDF file
      call check( nf90_close(ncid) )
      
      write(*,*) "CLOSE: ", INI_FILE
      
!---- Read COTS larvalrecruitment rate netCDF file --------------------------------

      write(*,*) "OPEN: ", IN_FILE1
      
      ! Open NetCDF grid file
      call check( nf90_open(IN_FILE1, nf90_nowrite, ncid) )
      call get_dimension(ncid, 'larval_time',  N_LARcots_time)
      
      allocate(LARcots_time(N_LARcots_time))
      allocate(LARcots(N_xi_rho, N_eta_rho, N_LARcots_time+1))
      
      ! Get variable id
      call check( nf90_inq_varid(ncid, 'larval_time', var_id) ) 
      call check( nf90_get_var(ncid, var_id, LARcots_time(:)) )
      call check( nf90_inq_varid(ncid, 'larvae', var_id) ) 
      call check( nf90_get_var(ncid, var_id, LARcots(:,:,1:N_LARcots_time)) )
      ! Close NetCDF file
      call check( nf90_close(ncid) )
      
      write(*,*) "CLOSE: ", IN_FILE1
      
      LARcots(:,:,N_LARcots_time+1)=0.0d0

!---- Read Coral larvalrecruitment rate netCDF file --------------------------------

      write(*,*) "OPEN: ", IN_FILE2
      
      ! Open NetCDF grid file
      call check( nf90_open(IN_FILE2, nf90_nowrite, ncid) )
      call get_dimension(ncid, 'larval_time',  N_LARcoral_time)
      
      allocate(LARcoral_time(N_LARcoral_time))
      allocate(LARcoral(N_xi_rho, N_eta_rho, N_LARcoral_time+1))
      
      ! Get variable id
      call check( nf90_inq_varid(ncid, 'larval_time', var_id) ) 
      call check( nf90_get_var(ncid, var_id, LARcoral_time(:)) )
      call check( nf90_inq_varid(ncid, 'larvae', var_id) ) 
      call check( nf90_get_var(ncid, var_id, LARcoral(:,:,1:N_LARcoral_time)) )
      ! Close NetCDF file
      call check( nf90_close(ncid) )
      
      write(*,*) "CLOSE: ", IN_FILE2

      LARcoral(:,:,N_LARcoral_time+1)=0.0d0

!---- Create the COTS_his NetCDF file --------------------------------
          
      call createNetCDFcots_his(           &
!          input parameters
     &        OUT_FILE                     &
     &      , TIME_ATT                     &
     &      , N_xi_rho, N_eta_rho, 1       &
     &      )
     
!---- Write COTS initial condition --------------------------------

      write(*,*) '                                        '
      write(*,*) '****************************************'
      write(*,*) 'time (days): ', cots_time(1)
      write(*,*) '                                        '

      start1D = (/ 1 /)
      count1D = (/ 1 /)

      call writeNetCDF_1d(           &
!          input parameters
     &        'cots_time'            &
     &      , OUT_FILE               &
     &      , 1                      &
     &      , cots_time              &
     &      , start1D, count1D       &
     &      )
     
      start3D = (/ 1,  1,  1 /)
      count3D = (/ N_xi_rho, N_eta_rho, 1 /)
          
      call writeNetCDF_3d(                      &
!          input parameters
     &        'cots1'                           &
     &      , OUT_FILE                          &
     &      , N_xi_rho, N_eta_rho, 1            &
     &      , COTS(ng)%dens(1,:,:)              &
     &      , start3D, count3D                  &
     &      )
      call writeNetCDF_3d(                      &
!          input parameters
     &        'cots2'                           &
     &      , OUT_FILE                          &
     &      , N_xi_rho, N_eta_rho, 1            &
     &      , COTS(ng)%dens(2,:,:)              &
     &      , start3D, count3D                  &
     &      )
      call writeNetCDF_3d(                      &
!          input parameters
     &        'cots3'                           &
     &      , OUT_FILE                          &
     &      , N_xi_rho, N_eta_rho, 1            &
     &      , COTS(ng)%dens(3,:,:)              &
     &      , start3D, count3D                  &
     &      )
      call writeNetCDF_3d(                      &
!          input parameters
     &        'cots4'                           &
     &      , OUT_FILE                          &
     &      , N_xi_rho, N_eta_rho, 1            &
     &      , COTS(ng)%dens(4,:,:)              &
     &      , start3D, count3D                  &
     &      )
      call writeNetCDF_3d(                      &
!          input parameters
     &        'cots5'                           &
     &      , OUT_FILE                          &
     &      , N_xi_rho, N_eta_rho, 1            &
     &      , COTS(ng)%dens(5,:,:)              &
     &      , start3D, count3D                  &
     &      )
      
      cots_sum(:,:)=0.0d0
      do k=1,COTS(ng)%Nstg
        cots_sum(:,:)=cots_sum(:,:)+COTS(ng)%dens(k,:,:)
      end do
      
      call writeNetCDF_3d(                      &
!          input parameters
     &        'cots_sum'                        &
     &      , OUT_FILE                          &
     &      , N_xi_rho, N_eta_rho, 1            &
     &      , cots_sum(:,:)                     &
     &      , start3D, count3D                  &
     &      )
     
      call writeNetCDF_3d(                      &
!          input parameters
     &        'p_coral'                         &
     &      , OUT_FILE                          &
     &      , N_xi_rho, N_eta_rho, 1            &
     &      , P_CORAL(ng)%cover(:,:)            &
     &      , start3D, count3D                  &
     &      )
     
      write(*,*) '****************************************'
      write(*,*) '                                        '

!-----------------------------------------------------------

      istep=0
      iprint=1

!----- Main loop -------------------------------------------

      do istep=1, int(365.0d0/dt) * 15 +1      ! 15 years
        
        do i=1,N_LARcots_time
          if(int(cots_time(1))==int(LARcots_time(i))) then
            n_larcots=i
            exit
          else
            n_larcots=N_LARcots_time+1
          end if
        end do
        
        do i=1,N_LARcoral_time
          if(int(cots_time(1))==int(LARcoral_time(i))) then
            n_larcoral=i
            exit
          else
            n_larcoral=N_LARcoral_time+1
          end if
        end do
        
!        write(*,*) n_larcots, n_larcoral
!        write(*,*) int(cots_time(1)),int(LARcots_time(1))

!-----------------------------------------------------------
!    COTS model
!-----------------------------------------------------------

        CALL cots_behavior                  &
!          input parameters
     &            (ng                       &   ! ng: nested grid number; n: coral compartment; i,j: position
     &            ,0, Im, 0, Jm             &
     &            ,dt                       &   ! Time step (day)
     &            ,dx                       &   ! dx: x grid size (m)
     &            ,dy                       &   ! dy: y grid size (m)
     &            ,h                        &   ! h(LBi:UBi,LBj:UBj) : water depth (m)
     &            ,P_CORAL(ng)%cover        &   ! sand(LBi:UBi,LBj:UBj) : sand coverage (0-1)
     &            ,p_sand                   &   ! sand(LBi:UBi,LBj:UBj) : sand coverage (0-1)
     &            ,LARcots(:,:,n_larcots)   &   ! COTS larval recruitment rate (individual m-2 d-1)
     &            ,cots_time(1)             &   ! COTS model time (days since initialization)
     &             )


!-----------------------------------------------------------
!    Coral population dynamics model
!-----------------------------------------------------------
        
        CALL p_coral_dynamics                & 
!          input parameters
     &            (ng                        &   ! ng: nested grid number; n: coral compartment; i,j: position
     &            ,0, Im, 0, Jm              &
     &            ,dt                        &   ! Time step (day)
     &            ,dx                        &   ! dx: x grid size (m)
     &            ,dy                        &   ! dy: y grid size (m)
     &            ,h                         &   ! h(LBi:UBi,LBj:UBj) : water depth (m)
     &            ,COTS(ng)%Nstg             &   ! sand(LBi:UBi,LBj:UBj) : sand coverage (0-1)
     &            ,COTS(ng)%pred             &   ! COTS_pred(Nstage, LBi:UBi,LBj:UBj) : sand coverage (0-1)
     &            ,p_sand                    &   ! sand(LBi:UBi,LBj:UBj) : sand coverage (0-1)
     &            ,LARcoral(:,:,n_larcoral)  &   ! Coral larval recruitment rate (individual m-2 d-1)
     &             )
        
!-----------------------------------------------------------

        cots_time(1) = cots_time(1) + dt
        
        write(*,*) 'time (days): ', cots_time(1)
        
! ***** Print section *****
        
        if ( mod(istep, IntPrint) == 0 ) then
          iprint = iprint + 1
          
          write(*,*) '                                        '
          write(*,*) '****************************************'
          write(*,*) 'time (days): ', cots_time(1)
          write(*,*) '                                        '

          
!---- Write COTS initial condition --------------------------------

          start1D = (/ iprint /)
          count1D = (/ 1 /)

          call writeNetCDF_1d(           &
!              input parameters
     &            'cots_time'            &
     &          , OUT_FILE               &
     &          , 1                      &
     &          , cots_time              &
     &          , start1D, count1D       &
     &          )
     
          start3D = (/ 1,  1,  iprint /)
          count3D = (/ N_xi_rho, N_eta_rho, 1 /)
              
          call writeNetCDF_3d(                      &
!              input parameters
     &            'cots1'                           &
     &          , OUT_FILE                          &
     &          , N_xi_rho, N_eta_rho, 1            &
     &          , COTS(ng)%dens(1,:,:)              &
     &          , start3D, count3D                  &
     &          )
          call writeNetCDF_3d(                      &
!              input parameters
     &            'cots2'                           &
     &          , OUT_FILE                          &
     &          , N_xi_rho, N_eta_rho, 1            &
     &          , COTS(ng)%dens(2,:,:)              &
     &          , start3D, count3D                  &
     &          )
          call writeNetCDF_3d(                      &
!              input parameters
     &            'cots3'                           &
     &          , OUT_FILE                          &
     &          , N_xi_rho, N_eta_rho, 1            &
     &          , COTS(ng)%dens(3,:,:)              &
     &          , start3D, count3D                  &
     &          )
          call writeNetCDF_3d(                      &
!              input parameters
     &            'cots4'                           &
     &          , OUT_FILE                          &
     &          , N_xi_rho, N_eta_rho, 1            &
     &          , COTS(ng)%dens(4,:,:)              &
     &          , start3D, count3D                  &
     &          )
          call writeNetCDF_3d(                      &
!              input parameters
     &            'cots5'                           &
     &          , OUT_FILE                          &
     &          , N_xi_rho, N_eta_rho, 1            &
     &          , COTS(ng)%dens(5,:,:)              &
     &          , start3D, count3D                  &
     &          )

          cots_sum(:,:)=0.0d0
          do k=1,COTS(ng)%Nstg
            cots_sum(:,:)=cots_sum(:,:)+COTS(ng)%dens(k,:,:)
          end do
          
          call writeNetCDF_3d(                      &
!              input parameters
     &            'cots_sum'                        &
     &          , OUT_FILE                          &
     &          , N_xi_rho, N_eta_rho, 1            &
     &          , cots_sum(:,:)                     &
     &          , start3D, count3D                  &
     &          )
     
          call writeNetCDF_3d(                      &
!              input parameters
     &            'p_coral'                         &
     &          , OUT_FILE                          &
     &          , N_xi_rho, N_eta_rho, 1            &
     &          , P_CORAL(ng)%cover(:,:)            &
     &          , start3D, count3D                  &
     &          )
     
          write(*,*) '****************************************'
          write(*,*) '                                        '
        end if
     
      enddo
      
!----- End loop --------------------------------------

      write(*,*) 'FINISH!!'

      RETURN

      END PROGRAM main_cots
!----------------------------------------------------------------------!

!     End of main program

!-----------------------------------------------------------------------

