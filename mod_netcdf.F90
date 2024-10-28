
!!!=== ver 2015/10/30   Copyright (c) 2013-2015 Takashi NAKAMURA  =====

#include "cppdefs.h"

!!!**** MODULE OF NETCDF READ & WRITE ****************

  module mod_netcdf
  
    USE netcdf
    
    implicit none

  CONTAINS


!**** create cost_his NetCDF file **********************************************

      SUBROUTINE createNetCDFcots_his(   &
!        input parameters
     &      OUT_FILE                     &
     &    , TIME_ATT                     &
     &    , Im, Jm, Nt                   &
     &)
                               
!    input parameters
      character(len=*),  intent( in) :: OUT_FILE
      character(len=*),  intent( in) :: TIME_ATT
      integer, intent( in) :: Im, Jm, Nt
      
      integer :: ncid,var_id
      integer :: xi_rho_dimid, eta_rho_dimid
      integer :: time_dimid
      integer :: dim3Dids(3)
      

      write(*,*) "CREATE: ", OUT_FILE

      call check( nf90_create(OUT_FILE, nf90_clobber, ncid) )

      call check( nf90_def_dim(ncid, 'xi_rho', Im, xi_rho_dimid) )
      call check( nf90_def_dim(ncid, 'eta_rho', Jm, eta_rho_dimid) )
      call check( nf90_def_dim(ncid, 'cots_time', NF90_UNLIMITED, time_dimid) )
      
    ! Define the netCDF variables.

      call check( nf90_def_var(ncid, 'cots_time', NF90_DOUBLE, time_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'time since initialization') )
      call check( nf90_put_att(ncid, var_id, 'units',     TIME_ATT ) )

      dim3Dids = (/ xi_rho_dimid, eta_rho_dimid, time_dimid /)

      call check( nf90_def_var(ncid, 'cots1', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name'              &
     &    , 'density of COTS on life stage 1: 0.5-5.5 month (corallin algae eater)') )
      call check( nf90_put_att(ncid, var_id, 'units',     'individuals meter-2') )
      call check( nf90_put_att(ncid, var_id, 'time',      'cots_time') )

      call check( nf90_def_var(ncid, 'cots2', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name'              &
     &    , 'density of COTS on life stage 2: 5.5 month - 1 year (coral eater)') )
      call check( nf90_put_att(ncid, var_id, 'units',     'individuals meter-2') )
      call check( nf90_put_att(ncid, var_id, 'time',      'cots_time') )

      call check( nf90_def_var(ncid, 'cots3', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name'              &
     &    , 'density of COTS on life stage 3: 1-2 years (coral eater)') )
      call check( nf90_put_att(ncid, var_id, 'units',     'individuals meter-2') )
      call check( nf90_put_att(ncid, var_id, 'time',      'cots_time') )
      
      call check( nf90_def_var(ncid, 'cots4', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name'              &
     &    , 'density of COTS on life stage 4: 2-3 years (coral eater)') )
      call check( nf90_put_att(ncid, var_id, 'units',     'individuals meter-2') )
      call check( nf90_put_att(ncid, var_id, 'time',      'cots_time') )
      
      call check( nf90_def_var(ncid, 'cots5', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name'              &
     &    , 'density of COTS on life stage 5: 3- years (coral eater, adult stage)') )
      call check( nf90_put_att(ncid, var_id, 'units',     'individuals meter-2') )
      call check( nf90_put_att(ncid, var_id, 'time',      'cots_time') )
      
      call check( nf90_def_var(ncid, 'cots_sum', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name'              &
     &    , 'density of COTS on all life stages') )
      call check( nf90_put_att(ncid, var_id, 'units',     'individuals meter-2') )
      call check( nf90_put_att(ncid, var_id, 'time',      'cots_time') )

      call check( nf90_def_var(ncid, 'p_coral', NF90_DOUBLE, dim3Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'coral coverage') )
      call check( nf90_put_att(ncid, var_id, 'units',     '0 to 1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'cots_time') )

  ! End define mode.
      call check( nf90_enddef(ncid) )
      call check( nf90_close(ncid) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE createNetCDFcots_his


!**** create cost_flt NetCDF file **********************************************

      SUBROUTINE createNetCDFcots_flt(   &
!        input parameters
     &      OUT_FILE                     &
     &    , TIME_ATT                     &
     &    , Im, Jm                       &
     &)
                               
!    input parameters
      character(len=*),  intent( in) :: OUT_FILE
      character(len=*),  intent( in) :: TIME_ATT
      integer, intent( in) :: Im, Jm
      
      integer :: ncid,var_id
      integer :: id_flt_dimid
      integer :: time_dimid
      integer :: dim2Dids(2)
      
      write(*,*) "CREATE: ", OUT_FILE

      call check( nf90_create(OUT_FILE, nf90_clobber, ncid) )

      call check( nf90_def_dim(ncid, 'id_flt', Im, id_flt_dimid) )
      call check( nf90_def_dim(ncid, 'ocean_time', NF90_UNLIMITED, time_dimid) )
      
    ! Define the netCDF variables.

      call check( nf90_def_var(ncid, 'ocean_time', NF90_DOUBLE, time_dimid, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name', 'time since initialization') )
      call check( nf90_put_att(ncid, var_id, 'units',     TIME_ATT ) )

      dim2Dids = (/ id_flt_dimid, time_dimid /)

      call check( nf90_def_var(ncid, 'num_COTS_larvae', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name'              &
     &    , 'Number of COTS larvae per particle') )
      call check( nf90_put_att(ncid, var_id, 'units',     'individuals particle-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )

      call check( nf90_def_var(ncid, 'mort_COTS_larvae', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name'              &
     &    , 'Mortality of COTS larvae') )
      call check( nf90_put_att(ncid, var_id, 'units',     'day-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )
      
      call check( nf90_def_var(ncid, 'PHY_COTS_larvae', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name'              &
     &    , 'Phytoplankton density around COTS larvae') )
      call check( nf90_put_att(ncid, var_id, 'units',     'umolC L-1') )
      call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )

      call check( nf90_def_var(ncid, 'lipid_COTS_larvae', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name'              &
     &    , 'Lipid of COTS larvae') )
      call check( nf90_put_att(ncid, var_id, 'units',     '???') )
      call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )

      call check( nf90_def_var(ncid, 'age_COTS_larvae', NF90_DOUBLE, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name'              &
     &    , 'Age of COTS larvae') )
      call check( nf90_put_att(ncid, var_id, 'units',     'days') )
      call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )

      call check( nf90_def_var(ncid, 'status_COTS_larvae', NF90_INT, dim2Dids, var_id) )
      call check( nf90_put_att(ncid, var_id, 'long_name'              &
     &    , 'Status of COTS larvae') )
      call check( nf90_put_att(ncid, var_id, 'units',     '??') )  !!!0 not started, 1 peragic stage1, 2 peragic stage2, 3 peragic stage3, 4 settled, 5 out
      call check( nf90_put_att(ncid, var_id, 'time',      'ocean_time') )


  ! End define mode.
      call check( nf90_enddef(ncid) )
      call check( nf90_close(ncid) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE createNetCDFcots_flt
      
      
!**** writeNetCDF_1d **********************************************
      
      SUBROUTINE writeNetCDF_1d(   &
!        input parameters
     &      NCNAME                 &
     &    , OUT_FILE               &
     &    , Im                     &
     &    , data                   &
     &    , start1D, count1D       &
     &)
                               
!    input parameters
      character(len=*), intent( in) :: NCNAME
      character(len=*), intent( in) :: OUT_FILE
      integer, intent( in) :: Im 
      real(8), intent( in) :: data(Im )
      integer, intent( in) :: start1D(1), count1D(1)
      
      integer :: ncid,var_id
      
! --- Write NetCDF file ------------------------
      
      write(*,*) "WRITE ", NCNAME," to ", OUT_FILE
      call check( nf90_open(OUT_FILE, NF90_WRITE, ncid) )
      call check( nf90_inq_varid(ncid, NCNAME, var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start1D, count = count1D) )
      call check( nf90_close(ncid) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE writeNetCDF_1d
      
!**** writeNetCDF_2d **********************************************
      
      SUBROUTINE writeNetCDF_2d(   &
!        input parameters
     &      NCNAME                 &
     &    , OUT_FILE               &
     &    , Im, Jm                 &
     &    , data                   &
     &    , start2D, count2D       &
     &)
                               
!    input parameters
      character(len=*), intent( in) :: NCNAME
      character(len=*), intent( in) :: OUT_FILE
      integer, intent( in) :: Im, Jm
      real(8), intent( in) :: data(Im, Jm)
      integer, intent( in) :: start2D(2), count2D(2)
      
      integer :: ncid,var_id
      
! --- Write NetCDF file ------------------------
      
      write(*,*) "WRITE ", NCNAME," to ", OUT_FILE
      call check( nf90_open(OUT_FILE, NF90_WRITE, ncid) )
      call check( nf90_inq_varid(ncid, NCNAME, var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start2D, count = count2D) )
      call check( nf90_close(ncid) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE writeNetCDF_2d
      
!**** writeIntNetCDF_2d **********************************************

      SUBROUTINE writeIntNetCDF_2d(   &
!        input parameters
     &      NCNAME                 &
     &    , OUT_FILE               &
     &    , Im, Jm                 &
     &    , data                   &
     &    , start2D, count2D       &
     &)
                               
!    input parameters
      character(len=*), intent( in) :: NCNAME
      character(len=*), intent( in) :: OUT_FILE
      integer, intent( in) :: Im, Jm
      integer, intent( in) :: data(Im, Jm)
      integer, intent( in) :: start2D(2), count2D(2)
      
      integer :: ncid,var_id
      
! --- Write NetCDF file ------------------------
      
      write(*,*) "WRITE ", NCNAME," to ", OUT_FILE
      call check( nf90_open(OUT_FILE, NF90_WRITE, ncid) )
      call check( nf90_inq_varid(ncid, NCNAME, var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start2D, count = count2D) )
      call check( nf90_close(ncid) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE writeIntNetCDF_2d
      
!**** writeNetCDF_3d **********************************************
      
      SUBROUTINE writeNetCDF_3d(   &
!        input parameters
     &      NCNAME                 &
     &    , OUT_FILE               &
     &    , Im, Jm, Nt             &
     &    , data                   &
     &    , start3D, count3D       &
     &)
                               
!    input parameters
      character(len=*), intent( in) :: NCNAME
      character(len=*), intent( in) :: OUT_FILE
      integer, intent( in) :: Im, Jm, Nt 
      real(8), intent( in) :: data(Im, Jm, Nt )
      integer, intent( in) :: start3D(3), count3D(3)
      
      integer :: ncid,var_id
      
! --- Write NetCDF file ------------------------
      
      write(*,*) "WRITE ", NCNAME," to ", OUT_FILE
      call check( nf90_open(OUT_FILE, NF90_WRITE, ncid) )
      call check( nf90_inq_varid(ncid, NCNAME, var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start3D, count = count3D) )
      call check( nf90_close(ncid) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE writeNetCDF_3d
      
!**** writeNetCDF_4d **********************************************
      
      SUBROUTINE writeNetCDF_4d(   &
!        input parameters
     &      NCNAME                 &
     &    , OUT_FILE               &
     &    , Im, Jm, Nz, Nt         &
     &    , data                   &
     &    , start4D, count4D       &
     &)
                               
!    input parameters
      character(len=*), intent( in) :: NCNAME
      character(len=*), intent( in) :: OUT_FILE
      integer, intent( in) :: Im, Jm, Nz, Nt
      real(8), intent( in) :: data(Im, Jm, Nz, Nt)
      integer, intent( in) :: start4D(4), count4D(4)
      
      integer :: ncid,var_id
      
! --- Write NetCDF file ------------------------
      
      write(*,*) "WRITE ", NCNAME," to ", OUT_FILE
      call check( nf90_open(OUT_FILE, NF90_WRITE, ncid) )
      call check( nf90_inq_varid(ncid, NCNAME, var_id) )
      call check( nf90_put_var(ncid, var_id, data, start = start4D, count = count4D) )
      call check( nf90_close(ncid) )
      write(*,*) '*** SUCCESS'

      END SUBROUTINE writeNetCDF_4d
      
!**** readNetCDF_3d **********************************************
      
      SUBROUTINE readNetCDF_3d(    &
!        input parameters
     &      NC_FILE                &
     &    , NCNAME                 &
     &    , Im, Jm, Nt             &
     &    , start3D, count3D       &
!        output parameters
     &    , data                   &
     &)
                               
!    input parameters
      character(len=*), intent( in) :: NC_FILE
      character(len=*), intent( in) :: NCNAME
      integer, intent( in) :: Im, Jm, Nt
      integer, intent( in) :: start3D(3), count3D(3)
      real(8), intent(out) :: data(Im, Jm, Nt)
      
      integer :: ncid,var_id
      integer :: err_flag
      real(8) :: sf, off
      
! --- Read NetCDF file ------------------------
      
      do
        write(*,*) "OPEN: ", NC_FILE
        call check( nf90_open(NC_FILE, nf90_nowrite, ncid) )
      
        write(*,*) 'DOWNLOAD ', NCNAME
        
!       Get variable id
        call check2( nf90_inq_varid(ncid, NCNAME, var_id), err_flag ) ! Water Temperature (degC)
        if(err_flag == 1) then
          write(*,*) '*** FAILED 1: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_var(ncid, var_id, data, start=start3D, count=count3D), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** DOWNLOAD FAILED: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_att(ncid, var_id, 'scale_factor', sf), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** FAILED 2: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_att(ncid, var_id, 'add_offset', off), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** FAILED 3: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check( nf90_close(ncid) )
        write(*,*) "CLOSE: ", NC_FILE
        exit
      end do

      data(:,:,:)=data(:,:,:)*sf+off
      write(*,*) '*** SUCCESS'


      END SUBROUTINE readNetCDF_3d
      
!**** readNetCDF_4d **********************************************
      
      SUBROUTINE readNetCDF_4d(    &
!        input parameters
     &      NC_FILE                &
     &    , NCNAME                 &
     &    , Im, Jm, Nz, Nt         &
     &    , start4D, count4D       &
!        output parameters
     &    , data                   &
     &)
                               
!    input parameters
      character(len=*), intent( in) :: NC_FILE
      character(len=*), intent( in) :: NCNAME
      integer, intent( in) :: Im, Jm, Nz, Nt
      integer, intent( in) :: start4D(4), count4D(4)
      real(8), intent(out) :: data(Im, Jm, Nz, Nt)
      
      integer :: ncid,var_id
      integer :: err_flag
      real(8) :: sf, off
      
! --- Read NetCDF file ------------------------
      
      do
        write(*,*) "OPEN: ", NC_FILE
        call check( nf90_open(NC_FILE, nf90_nowrite, ncid) )
      
        write(*,*) 'DOWNLOAD ', NCNAME
        
!       Get variable id
        call check2( nf90_inq_varid(ncid, NCNAME, var_id), err_flag ) ! Water Temperature (degC)
        if(err_flag == 1) then
          write(*,*) '*** FAILED 1: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_var(ncid, var_id, data, start=start4D, count=count4D), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** DOWNLOAD FAILED: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_att(ncid, var_id, 'scale_factor', sf), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** FAILED 2: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check2( nf90_get_att(ncid, var_id, 'add_offset', off), err_flag  )
        if(err_flag == 1) then
          write(*,*) '*** FAILED 3: Retry!'
          call check( nf90_close(ncid) )
          write(*,*) "CLOSE: ", NC_FILE
          cycle
        end if
        
        call check( nf90_close(ncid) )
        write(*,*) "CLOSE: ", NC_FILE
        exit
      end do
      
      data(:,:,:,:)=data(:,:,:,:)*sf+off
      write(*,*) '*** SUCCESS'


      END SUBROUTINE readNetCDF_4d
      
!**** NetCDF utility **********************************************
      
      SUBROUTINE get_dimension(ncid, name, dim)
      
      integer,           intent( in) :: ncid
      character(len=*),  intent( in) :: name
      integer,           intent(out) :: dim

      integer :: dimid
      call check( nf90_inq_dimid(ncid, name, dimid) )
      call check( nf90_inquire_dimension(ncid, dimid, len=dim) )
      
      END SUBROUTINE get_dimension
! -------------------------------------------------------------------------
      SUBROUTINE  get_dimension2(ncid, name, dim, dims)
      integer,           intent( in) :: ncid
      character(len=*),  intent( in) :: name
      integer,           intent(out) :: dim
      real, allocatable, intent(out) :: dims(:)

      integer :: varid, dimid
      integer :: err

      call check( nf90_inq_dimid(ncid, name, dimid) )
      call check( nf90_inquire_dimension(ncid, dimid, len=dim) )
      allocate(dims(dim), stat=err)
      if (err /= 0) print *, name, ": Allocation request denied"
      call check( nf90_inq_varid(ncid, name, varid) )
      call check( nf90_get_var(ncid, varid, dims) )
      END SUBROUTINE get_dimension2

! -------------------------------------------------------------------------

      SUBROUTINE check(status)
      
      integer, intent(in) :: status

      if (status /= nf90_noerr) then 
          print *, trim(nf90_strerror(status))
          stop "Stopped"
      end if
      
      END SUBROUTINE check
      
! -------------------------------------------------------------------------
      
      SUBROUTINE check2(status, err_flag)
      
      integer, intent( in) :: status
      integer, intent(out) :: err_flag
      
      err_flag = 0

      if (status /= nf90_noerr) then 
          print *, trim(nf90_strerror(status))
          err_flag = 1
!          stop "Stopped"
      end if
      
      END SUBROUTINE check2

  END MODULE mod_netcdf

