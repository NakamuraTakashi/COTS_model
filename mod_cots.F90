
!!!=== ver 2016/03/09   Copyright (c) 2013-2016 Takashi NAKAMURA  =====

#include "cppdefs.h"


!!!**** MODULE OF COTS (Crown-of-thorns starfish) MODEL ****************

  module mod_cots

    implicit none

    integer, public, parameter :: Nstage = 5    !! Number of COTS life stage
    
    TYPE T_COTS
    
      real(8), pointer :: dens(:,:,:)       ! COTS density (Individual m-2)
      real(8), pointer :: pred(:,:,:)       ! COTS predetion rate (m2 d-1)
      real(8), pointer :: egg(:,:)          ! COTS egg production rate (m2 d-1)
      integer, pointer :: Nstg              ! Number of COTS life stage

    END TYPE T_COTS

    TYPE (T_COTS), allocatable :: COTS(:)

!-----------------------------------------------------------------------

  CONTAINS


!!! **********************************************************************
!!!  set initial conditions for COTS (Crown-of-thorns starfish) model
!!! **********************************************************************

    SUBROUTINE initialize_cots(ng, Ngrids, LBi, UBi, LBj, UBj)

      implicit none
! input parameters
      integer, intent(in) :: ng, Ngrids, LBi, UBi, LBj, UBj

      integer i,j,n

      IF (ng.eq.1) allocate ( COTS(Ngrids) )

      allocate( COTS(ng)%dens(Nstage,LBi:UBi,LBj:UBj)     )
      allocate( COTS(ng)%pred(Nstage,LBi:UBi,LBj:UBj)     )
      allocate( COTS(ng)%egg (LBi:UBi,LBj:UBj)     )
      allocate( COTS(ng)%Nstg     )


!----------set data -----------------------
      COTS(ng)%Nstg = Nstage
!  Set initial conditions

      do j=LBj,UBj
        do i=LBi,UBi
          COTS(ng)%dens(1,i,j)=0.0d0 ! Life stage 1: 0.5-5.5 month (corallin algae eater)
          COTS(ng)%dens(2,i,j)=0.0d0 ! Life stage 2: 5.5 month - 1 year (coral eater)
          COTS(ng)%dens(3,i,j)=0.0d0 ! Life stage 3: 1-2 years (coral eater)
          COTS(ng)%dens(4,i,j)=0.0d0 ! Life stage 4: 2-3 years (coral eater)
          COTS(ng)%dens(5,i,j)=0.0d0 ! Life stage 5: 3- years (coral eater, adult stage)
        enddo
      enddo
      
      do j=LBj,UBj
        do i=LBi,UBi
          do n=1,Nstage
            COTS(ng)%pred(n,i,j)=0.d0 ! COTS initial predetion rate
          enddo
        enddo
      enddo

      RETURN
    END SUBROUTINE initialize_cots
      

!!! **********************************************************************
!!!  Main program of COTS (Crown-of-thorns starfish) model 
!!! **********************************************************************

    SUBROUTINE cots_behavior             &
!          input parameters
     &            (ng                     &   ! ng: nested grid number; n: coral compartment; i,j: position
     &            ,LBi, UBi, LBj, UBj     &
     &            ,dt             &   ! Time step (day)
     &            ,dx             &   ! dx: x grid size (m)
     &            ,dy             &   ! dy: y grid size (m)
     &            ,h              &   ! h(LBi:UBi,LBj:UBj) : water depth (m)
     &            ,coral          &   ! coral(LBi:UBi,LBj:UBj) : coral coverage (0-1)
     &            ,sand           &   ! sand(LBi:UBi,LBj:UBj) : sand coverage (0-1)
     &            ,LARcots        &   ! COTS larval recruitment rate (individual m-2 d-1)
     &            ,time           &   ! COTS model time (days since initialization)
     &             )
!
      
      implicit none
      
! input parameters
      integer, intent(in) :: ng
      integer, intent(in) :: LBi, UBi, LBj, UBj
      real(8), intent(in) :: dt
      real(8), intent(in) :: dx      
      real(8), intent(in) :: dy      
      real(8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(8), intent(in) :: coral(LBi:UBi,LBj:UBj)
      real(8), intent(in) :: sand(LBi:UBi,LBj:UBj)
      real(8), intent(in) :: LARcots(LBi:UBi,LBj:UBj)
      real(8), intent(in) :: time

!!!------------Set parameters  ----------------------------------

!--- Other variables ----------------------------------------------
      real(8) :: mask_rho(LBi:UBi,LBj:UBj)      ! mask for rho points (m)   
      real(8) :: mask_u(LBi+1:UBi,LBj:UBj)      ! mask for u points (m)   
      real(8) :: mask_v(LBi:UBi,LBj+1:UBj)      ! mask for v points (m)   
      real(8) :: Fx(Nstage,LBi+1:UBi,LBj:UBj)      ! COTS flux at u points (d-1)   
      real(8) :: Fy(Nstage,LBi:UBi,LBj+1:UBj)      ! COTS flux at v points (d-1)   
      real(8) :: Umax(Nstage)
      real(8) :: Diff(Nstage)
      real(8) :: Mort(Nstage)
      real(8) :: Recr(Nstage)
      integer :: i, j, n
      integer :: flag150 = 0
      integer :: flag365 = 0
      real(8) :: u
!  Output
      real(8), save :: t_days = 0.d0  !day after spawning

      
! COTS maximum moving (m d-1)
      Umax(1)=1.0d1
      Umax(2)=1.0d1
      Umax(3)=5.0d1
      Umax(4)=1.0d2
      Umax(5)=3.0d2

! COTS random movement coefficients (m2 d-1)
      Diff(1)=1.0d2
      Diff(2)=1.0d2
      Diff(3)=1.0d2
      Diff(4)=1.0d2
      Diff(5)=1.0d2
      
      
! ---- COTS mask out -----------------------------------------------
!$omp parallel private(i,j,n)
!$omp do
      do i=LBi,UBi
        do j=LBj,UBj
          if(h(i,j) > 0.0d0) then
            mask_rho(i,j)=1.0d0
          else
            mask_rho(i,j)=0.0d0
          end if
        end do
      end do
!$omp end do
!$omp do
      do i=LBi+1,UBi
        do j=LBj,UBj
          mask_u(i,j)=mask_rho(i-1,j)*mask_rho(i,j)
        end do
      end do
!$omp end do
!$omp do
      do i=LBi,UBi
        do j=LBj+1,UBj
          mask_v(i,j)=mask_rho(i,j-1)*mask_rho(i,j)
        end do
      end do
!$omp end do

! ---- COTS flux calculation -----------------------------------------------

!$omp do private(u)

      do i=LBi+1,UBi
        do j=LBj,UBj
          do n=1,Nstage
!          Coral fastidious behavior
            u = Umax(n) * (coral(i,j)-coral(i-1,j))! * coral(i,j)
!          Coral fastidious behavior
!            u = u + Umax(n) * (h(i,j)-h(i-1,j))/dx
!          Rock fastidious behavior
            u = 0.5d0*(u+ABS(u)) * ((1.0d0-sand(i,j))+0.1d0*sand(i,j))         & ! u >= 0
     &         +0.5d0*(u-ABS(u)) * ((1.0d0-sand(i-1,j))+0.1d0*sand(i-1,j))       ! u <  0

            Fx(n,i,j) =                                                &
!      Advective flux
     &        ( 0.5d0*(u+ABS(u)) * COTS(ng)%dens(n,i-1,j)              & ! u >= 0
     &         +0.5d0*(u-ABS(u)) * COTS(ng)%dens(n,i,j)                & ! u <  0
     &        )/dx                                                     &
!      Diffusive flux
     &       +Diff(n)*( COTS(ng)%dens(n,i-1,j)                         &
     &                 -COTS(ng)%dens(n,i,j))/dx/dx

            Fx(n,i,j) = Fx(n,i,j) * mask_u(i,j)

          end do
        end do
      end do
!$omp end do

!$omp do private(u)
      do i=LBi,UBi
        do j=LBj+1,UBj
          do n=1,Nstage
!          Coral fastidious behavior
            u = Umax(n) * (coral(i,j)-coral(i,j-1))! * coral(i,j)
!          Coral fastidious behavior
!            u = u + Umax(n) * (h(i,j)-h(i,j-1))/dy
!          rock fastidious behavior
            u = 0.5d0*(u+ABS(u)) * ((1.0d0-sand(i,j))+0.1d0*sand(i,j))         & ! u >= 0
     &         +0.5d0*(u-ABS(u)) * ((1.0d0-sand(i,j-1))+0.1d0*sand(i,j-1))       ! u <  0

            Fy(n,i,j) =                                                &
!      Advective flux
     &        ( 0.5d0*(u+ABS(u)) * COTS(ng)%dens(n,i,j-1)              & ! u >= 0
     &         +0.5d0*(u-ABS(u)) * COTS(ng)%dens(n,i,j)                & ! u <  0
     &        )/dy                                                     &
!      Diffusive flux
     &       +Diff(n)*( COTS(ng)%dens(n,i,j-1)                         &
     &                 -COTS(ng)%dens(n,i,j))/dy/dy

            Fy(n,i,j) = Fy(n,i,j) * mask_v(i,j)

          end do
        end do
      end do
!$omp end do

! ---- COTS dynamics -----------------------------------------------

!$omp do private(Mort,Recr)
      do i=LBi+1,UBi-1
        do j=LBj+1,UBj-1

!      COTS mortality rate (Individual d-1)
          if(coral(i,j) <= 1.0d-2) then
            Mort(1)=5.0d-2
            Mort(2)=3.0d-2  !4.d-2 * (1.d0 - coral(i,j))
            Mort(3)=2.0d-2  !3.d-2 * (1.d0 - coral(i,j))
            Mort(4)=1.0d-2  !2.d-2 * (1.d0 - coral(i,j))
            Mort(5)=1.0d-2  !1.d-2 * (1.d0 - coral(i,j))
          else
            Mort(1)=1.0d-4!1.0d-3  !1.0d-2 * exp(-1.0d2*coral(i,j))  !0.d0  !1.d-2
            Mort(2)=2.0d-3  !4.d-2 * (1.d0 - coral(i,j))
            Mort(3)=1.0d-3  !3.d-2 * (1.d0 - coral(i,j))
            Mort(4)=1.0d-3  !2.d-2 * (1.d0 - coral(i,j))
            Mort(5)=1.0d-3  !1.d-2 * (1.d0 - coral(i,j))
          end if
          
!      COTS recruitment rate (Individual m-2 d-1)
          Recr(:)=0.0d0
          Recr(1)=LARcots(i,j)

          do n=1,Nstage
      
            COTS(ng)%dens(n,i,j)=COTS(ng)%dens(n,i,j)+(                 &
          ! Advection term (upstream difference scheme)
          !   x axis
     &          Fx(n,i,j)-Fx(n,i+1,j)                                   &
          !   y axis
     &        + Fy(n,i,j)-Fy(n,i,j+1)                                   &
          ! Reaction term
     &        + Recr(n)                                                 &
     &        - Mort(n) * COTS(ng)%dens(n,i,j)                          &
     &       ) *dt

            COTS(ng)%dens(n,i,j)=COTS(ng)%dens(n,i,j) * mask_rho(i,j)
            
            COTS(ng)%dens(n,i,j)=max(COTS(ng)%dens(n,i,j), 0.0d0)  !!! For Error handring
     
          end do

!--- Coral predation rates --------------------------------------

          COTS(ng)%pred(1,i,j)=COTS(ng)%dens(1,i,j)*coral(i,j) * 1.0d-5
          COTS(ng)%pred(2,i,j)=COTS(ng)%dens(2,i,j)*coral(i,j) * 1.0d-4
          COTS(ng)%pred(3,i,j)=COTS(ng)%dens(3,i,j)*coral(i,j) * 5.0d-3
          COTS(ng)%pred(4,i,j)=COTS(ng)%dens(4,i,j)*coral(i,j) * 5.0d-2
          COTS(ng)%pred(5,i,j)=COTS(ng)%dens(5,i,j)*coral(i,j) * 1.0d-1
        
        end do
      end do
!$omp end do
!$omp end parallel

!--- Change each life stage --------------------------------------
      
      t_days = time - aint(time/365.25)*365.25
      
      if ((int(t_days)==150.0d0) .and. (flag150==0)) then
        COTS(ng)%dens(2,:,:)=COTS(ng)%dens(1,:,:)
        COTS(ng)%dens(1,:,:)=0.d0
        flag150 = 1
        flag365 = 0
      end if
      if ((int(t_days)==0.0d0) .and. (flag365 == 0)) then
        COTS(ng)%dens(5,:,:)=COTS(ng)%dens(5,:,:)+COTS(ng)%dens(4,:,:)
        COTS(ng)%dens(4,:,:)=COTS(ng)%dens(3,:,:)
        COTS(ng)%dens(3,:,:)=COTS(ng)%dens(2,:,:)
        COTS(ng)%dens(2,:,:)=0.d0
        flag365 = 1
        flag150 = 0
      end if

!---------------------------------------------------------------------
      
      
      RETURN

     END SUBROUTINE cots_behavior

  END MODULE mod_cots

