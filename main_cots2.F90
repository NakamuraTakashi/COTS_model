
!!!=== ver 2016/03/10   Copyright (c) 2013-2016 Takashi NAKAMURA  =====

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
      USE mod_cots_larvae
      USE mod_p_coral
      
      implicit none
      
! -------------------------------------------------------------------------

      real(8), parameter :: dt = 1.0d0/12.0d0      ! time step (days)
      
      integer, parameter :: Syear  = 2000   ! Starting year
      integer, parameter :: Smonth = 6      ! Starting month
      integer, parameter :: Sday   = 1      ! Starting day
     ! NetCDF file     
      character(len=*), parameter :: GRID_FILE = &
     &   "D:/ROMS/Yaeyama/Data/Yaeyama2_grd_v9.nc"
      character(len=*), parameter :: INI_FILE  = "input/cots_ini_v3.nc"
      character(len=*), parameter :: HIS_FILE  = &
     &   "D:/ROMS/Yaeyama/Y2_1g_eco_130607_7/ocean_his_Yaeyama2_1g_eco_130607.nc"
      character(len=*), parameter :: FLT_FILE  = &
     &   "D:/ROMS/Yaeyama/Y2_1g_eco_off_130607/ocean_flt_Yaeyama2_1g_eco_off_2.nc"
      character(len=*), parameter :: OUT_FILE  = "D:/COTS_model/output/cots_his.nc"
      character(38) :: FLT_OUT_FILE  = "D:/COTS_model/output/cots_larvae_01.nc"

      real(8), parameter :: dx  = 300.0d0  ! X Grid size (m)
      real(8), parameter :: dy  = 300.0d0  ! Y Grid size (m) 
      
      integer, parameter :: Nburst  = 20
      integer, parameter :: N  = 1
      integer, parameter :: ng  = 1
      
! -------------------------------------------------------------------------

      real(8), allocatable :: h(:,:)              ! depth (meter)
      real(8), allocatable :: p_sand(:,:)         ! sand coverage (0-1)
      real(8), allocatable :: cots_sum(:,:)       ! density of COTS on all life stages (indiv. m-2)
      
      real(8), allocatable :: cots_rec(:,:)        ! COTS larval recruitment rate (individual m-3)
      real(8), allocatable :: coral_rec(:,:)       ! Coral larval recruitment rate (individual m-3)

!      real(8), allocatable :: LARcots_time(:)     ! 
!      real(8), allocatable :: LARcoral_time(:)    ! 
      real(8), allocatable :: his_time(:)     ! 
      real(8), allocatable :: flt_time(:)     ! 
      real(8), allocatable :: phy1(:,:)     ! 
      real(8), allocatable :: phy2(:,:)     ! 
      real(8), allocatable :: phy_tot(:,:)     ! 
     
      real(8) :: cots_time(1)            ! COTS model time

      integer :: Im, Jm
      integer :: i,j,k,id, Nid
      integer :: istep, iburst, iprint, iflt, iINFILE
      integer :: n_larcots, n_larcoral
      real(8) :: IntPrint
      real(8) :: time    !(day)
      
      character(33) :: TIME_ATT  = "days since 1992-01-01 00:00:00"
      character(37) :: TIME_ATT2 = "seconds since 1992-01-01 00:00:00"
      
      integer :: N_xi_rho, N_eta_rho
      integer :: N_flt
      integer :: N_LARcots_time, N_LARcoral_time
      integer :: N_his_time, N_flt_time
      integer :: ncid, var_id
      integer :: ncid_his, ncid_flt
      integer :: var_id_phy1, var_id_phy2
      integer :: var_id_xgrd, var_id_ygrd, var_id_zgrd
      integer :: start1D(1), count1D(1)
      integer :: start2D(2), count2D(2)
      integer :: start3D(3), count3D(3)
      integer :: start4D(4), count4D(4)
      character(4) :: YYYY
      character(2) :: MM
      character(2) :: DD
      
      integer :: first_step =1
      real(8) :: flt_start_day
      real(8) :: nanVal
      character(2) :: NN

      
      write (YYYY, "(I4.4)") Syear
      write (MM, "(I2.2)") Smonth
      write (DD, "(I2.2)") Sday
      
!---- Modify time-unit description ---------------------------------
      
      TIME_ATT(12:15)=YYYY
      TIME_ATT(17:18)=MM
      TIME_ATT(20:21)=DD
      
      TIME_ATT2(15:18)=YYYY
      TIME_ATT2(20:21)=MM
      TIME_ATT2(23:24)=DD

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
      allocate(coral_rec(N_xi_rho, N_eta_rho))
      allocate(cots_rec(N_xi_rho, N_eta_rho))
      
      allocate(phy1(N_xi_rho, N_eta_rho))
      allocate(phy2(N_xi_rho, N_eta_rho))
      allocate(phy_tot(N_xi_rho, N_eta_rho))
      
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
      IntPrint = 365.25d0/12.0d0
      
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
     
      COTS(ng)%dens(1,:,:)=COTS(ng)%dens(1,:,:) * 1.0d0 !!! 100.0d0  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
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


      istep=0
      iprint=1
      iflt = 1
      iINFILE = 1
      
      flt_start_day = 1.0d0

!****** Main loop ******************************************************

      do istep=1, int(365.25d0/dt) * 17 +1      ! 17 years
      
        cots_rec(:,:) = 0.0d0
       
!-----------------------------------------------------------
!    COTS larval behavior
!-----------------------------------------------------------
        if (flt_start_day <= cots_time(1) ) then
          if (first_step == 1 ) then
            !  ---- Read ROMS ecosystem output NetCDF file --------------------------------

            write(*,*) "OPEN: ", HIS_FILE
            
            ! Open NetCDF grid file
            call check( nf90_open(HIS_FILE, nf90_nowrite, ncid_his) )
            call get_dimension(ncid_his, 'ocean_time',  N_his_time)
            
            N_his_time = 1359  !!!! TEST mode **********************************************************
            
            allocate(his_time(N_his_time))
            
            ! Get variable id
            call check( nf90_inq_varid(ncid_his, 'ocean_time', var_id) ) 
            call check( nf90_get_var(ncid_his, var_id, his_time(:)) )
            call check( nf90_inq_varid(ncid_his, 'phytoplankton1', var_id_phy1) )
            call check( nf90_inq_varid(ncid_his, 'phytoplankton2', var_id_phy2) )
            call check( nf90_get_att(ncid_his,var_id_phy1,'_FillValue',nanVal) )    ! to get _FillValue

!  ---      - Read ROMS float output NetCDF file --------------------------------

            write(*,*) "OPEN: ", FLT_FILE
            
            ! Open NetCDF grid file
            call check( nf90_open(FLT_FILE, nf90_nowrite, ncid_flt) )
            call get_dimension(ncid_flt, 'ocean_time',  N_flt_time)
            call get_dimension(ncid_flt, 'drifter',  N_flt)  !!!!!!!!!!!!!!!!!!!!!
            
            N_flt_time = 1359  !!!! TEST mode **********************************************************
            
            allocate(flt_time(N_flt_time))
            
            CALL initialize_cots_larvae(1, 1, N_flt)
            
            ! Get variable id
            call check( nf90_inq_varid(ncid_flt, 'ocean_time', var_id) ) 
            call check( nf90_get_var(ncid_flt, var_id, flt_time(:)) )
            call check( nf90_inq_varid(ncid_flt, 'Xgrid', var_id_xgrd) ) 
            call check( nf90_inq_varid(ncid_flt, 'Ygrid', var_id_ygrd) ) 
            call check( nf90_inq_varid(ncid_flt, 'Zgrid', var_id_zgrd) ) 
            
!----       Create the COTS_flt NetCDF file --------------------------------
                
            call createNetCDFcots_flt(           &
!                input parameters
     &              FLT_OUT_FILE                 &
     &            , TIME_ATT2                    &
     &            , N_flt, N_flt_time            &
     &            )
     

          end if
          
          do i=1,4
           start4D = (/ 1,        1,         1,  iflt /)
           count4D = (/ N_xi_rho, N_eta_rho, 1,  1     /)
          
            call check( nf90_get_var(ncid_his, var_id_phy1, phy1, start=start4D, count=count4D) )
            call check( nf90_get_var(ncid_his, var_id_phy2, phy2, start=start4D, count=count4D) )
            
            start2D = (/ 1,     iflt  /)
            count2D = (/ N_flt, 1     /)
            call check( nf90_get_var(ncid_flt, var_id_xgrd, COTS_LAR(ng)%Xgrid, start=start2D, count=count2D ) )
            call check( nf90_get_var(ncid_flt, var_id_ygrd, COTS_LAR(ng)%Ygrid, start=start2D, count=count2D ) )
            call check( nf90_get_var(ncid_flt, var_id_zgrd, COTS_LAR(ng)%Zgrid, start=start2D, count=count2D ) )
            

!            phy_tot(:,:) = phy1(:,:)+phy2(:,:)   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Test 
            phy_tot(:,:) = 0.45d0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Test for low phy condition

            CALL cots_larvae               &
!              input parameters
     &                (ng                    &   ! ng: nested grid number; n: coral compartment; i,j: position
     &                ,0, Im, 0, Jm          &
     &                ,N_flt, N_flt_time     &
     &                ,dt*0.25d0             &   ! Time step (day)
     &                ,dx                    &   ! dx: x grid size (m)
     &                ,dy                    &   ! dy: y grid size (m)
     &                ,COTS(ng)%dens(5,:,:)  &   ! COTSdens(LBi:UBi,LBj:UBj) : adult COTS density (indiv./m2)
     &                ,phy_tot               &   ! PHY(LBi:UBi,LBj:UBj): Phytoplankton density at surface (umolC/L)
     &                ,nanVal                &   ! Nan value
     &                ,h                     &   ! h(LBi:UBi,LBj:UBj) : water depth (m)
     &                ,P_CORAL(ng)%cover     &   ! coral(LBi:UBi,LBj:UBj) : coral coverage (0-1)
     &                ,p_sand                &   ! sand(LBi:UBi,LBj:UBj) : sand coverage (0-1)
     &                ,time                  &   ! COTS model time (days since initialization)
!              input & output parameters
     &                ,first_step            &   ! Status 1: first_step
!              output parameters
     &                ,cots_rec              &   ! COTS larval recruitment rate (individual m-2 d-1)
     &                 )
     
! ***** Print section *****
        
            write(*,*) '                                        '
            write(*,*) '****************************************'
            write(*,*) 'time (days): ', cots_time(1)
            write(*,*) '                                        '

          
!---- Write COTS larval condition --------------------------------

            start1D = (/ iflt /)
            count1D = (/ 1 /)

            call writeNetCDF_1d(           &
!                input parameters
     &              'ocean_time'           &
     &            , FLT_OUT_FILE           &
     &            , 1                      &
     &            , flt_time(iflt)         &
     &            , start1D, count1D       &
     &            )
     
            start2D = (/ 1,     iflt /)
            count2D = (/ N_flt, 1    /)
                
            call writeNetCDF_2d(                      &
!                input parameters
     &              'num_COTS_larvae'                 &
     &            , FLT_OUT_FILE                      &
     &            , N_flt, N_flt_time                 &
     &            , COTS_LAR(ng)%num                  &
     &            , start2D, count2D                  &
     &            )
            call writeNetCDF_2d(                      &
!                input parameters
     &              'mort_COTS_larvae'                &
     &            , FLT_OUT_FILE                      &
     &            , N_flt, N_flt_time                 &
     &            , COTS_LAR(ng)%mort                 &
     &            , start2D, count2D                  &
     &            )
     
            call writeNetCDF_2d(                      &
!                input parameters
     &              'PHY_COTS_larvae'                 &
     &            , FLT_OUT_FILE                      &
     &            , N_flt, N_flt_time                 &
     &            , COTS_LAR(ng)%phy                  &
     &            , start2D, count2D                  &
     &            )
            call writeNetCDF_2d(                      &
!                input parameters
     &              'lipid_COTS_larvae'               &
     &            , FLT_OUT_FILE                      &
     &            , N_flt, N_flt_time                 &
     &            , COTS_LAR(ng)%lipid                &
     &            , start2D, count2D                  &
     &            )
            call writeNetCDF_2d(                      &
!                input parameters
     &              'age_COTS_larvae'               &
     &            , FLT_OUT_FILE                      &
     &            , N_flt, N_flt_time                 &
     &            , COTS_LAR(ng)%age                &
     &            , start2D, count2D                  &
     &            )
            call writeIntNetCDF_2d(                   &
!                input parameters
     &              'status_COTS_larvae'              &
     &            , FLT_OUT_FILE                      &
     &            , N_flt, N_flt_time                 &
     &            , COTS_LAR(ng)%status               &
     &            , start2D, count2D                  &
     &            )
     
     
            write(*,*) '****************************************'
            write(*,*) '                                        '
            
            iflt = iflt+1
            
            if (iflt == N_flt_time+1) then
              call check( nf90_close(ncid_his) )
              write(*,*) "CLOSE: ", HIS_FILE
              call check( nf90_close(ncid_flt) )
              write(*,*) "CLOSE: ", FLT_FILE
              CALL finalize_cots_larvae(1)
              deallocate( his_time )
              deallocate( flt_time )
              
              first_step = 1  !!! Reset first-step flag.
              flt_start_day = 1.0d0 + 365.25d0*dble(iINFILE)  !!!!!!!!! —vC³
              iINFILE = iINFILE + 1
              iflt = 1
!        ---- Modify time-unit description ---------------------------------
              write (NN, "(I2.2)") iINFILE
              FLT_OUT_FILE(34:35)=NN
              
              exit
            end if
      
          end do
          
        end if


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
     &            ,cots_rec                 &   ! COTS larval recruitment rate (individual m-2 d-1)
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
     &            ,coral_rec                 &   ! Coral larval recruitment rate (individual m-2 d-1)
     &             )
        
!-----------------------------------------------------------

        cots_time(1) = cots_time(1) + dt
        
        write(*,*) 'time (days): ', cots_time(1)
        
! ***** Print section *****
        
        if ( dble(istep)*dt >= IntPrint*iprint ) then
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
      ! Close NetCDF file

      write(*,*) 'FINISH!!'

      RETURN

      END PROGRAM main_cots
!----------------------------------------------------------------------!

!     End of main program

!-----------------------------------------------------------------------

