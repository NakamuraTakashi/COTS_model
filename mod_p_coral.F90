
!!!=== ver 2016/03/09   Copyright (c) 2013-2016 Takashi NAKAMURA  =====

#include "cppdefs.h"


!!!**** MODULE OF CORAL POPULATION DYNAMICS MODEL ****************

  module mod_p_coral

    implicit none

    TYPE T_P_CORAL
    
      real(8), pointer :: cover(:,:)      ! coral coverage (0-1)

    END TYPE T_P_CORAL

    TYPE (T_P_CORAL), allocatable :: P_CORAL(:)

!-----------------------------------------------------------------------

  CONTAINS


!!! **********************************************************************
!!!  set initial conditions for coral population dynamics model
!!! **********************************************************************

    SUBROUTINE initialize_p_coral(ng, Ngrids, LBi, UBi, LBj, UBj)

      implicit none
! input parameters
      integer, intent(in) :: ng, Ngrids, LBi, UBi, LBj, UBj

      integer i,j,n

      IF (ng.eq.1) allocate ( P_CORAL(Ngrids) )

      allocate( P_CORAL(ng)%cover(LBi:UBi,LBj:UBj)     )


!----------set data -----------------------

!  Set initial conditions

      do j=LBj,UBj
        do i=LBi,UBi
          P_CORAL(ng)%cover(i,j)=0.0d0
        enddo
      enddo
      
      RETURN
    END SUBROUTINE initialize_p_coral
      

!!! **********************************************************************
!!!  Main program of coral population dynamics model
!!! **********************************************************************

    SUBROUTINE p_coral_dynamics           &
!          input parameters
     &            (ng                     &   ! ng: nested grid number; n: coral compartment; i,j: position
     &            ,LBi, UBi, LBj, UBj     &
     &            ,dt             &   ! Time step (day)
     &            ,dx             &   ! dx: x grid size (m)
     &            ,dy             &   ! dy: y grid size (m)
     &            ,h              &   ! h(LBi:UBi,LBj:UBj) : water depth (m)
     &            ,Nstage         &   ! sand(LBi:UBi,LBj:UBj) : sand coverage (0-1)
     &            ,COTS_pred      &   ! COTS_pred(Nstage, LBi:UBi,LBj:UBj) : sand coverage (0-1)
     &            ,sand           &   ! sand(LBi:UBi,LBj:UBj) : sand coverage (0-1)
     &            ,LARcoral       &   ! Coral larval recruitment rate (individual m-3 d-1)
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
      integer, intent(in) :: Nstage
      real(8), intent(in) :: COTS_pred(Nstage, LBi:UBi,LBj:UBj)
      real(8), intent(in) :: sand(LBi:UBi,LBj:UBj)
      real(8), intent(in) :: LARcoral(LBi:UBi,LBj:UBj)

!!!------------Set parameters  ----------------------------------
      real(8), parameter :: Gc = 1.0d-3      ! Growt rate (d-1)

      real(8) :: Growth
      real(8) :: Mort
      real(8) :: Recr
      real(8) :: PFD
      integer :: i, j, n
      real(8) :: u
!  Output
      real(8), save :: t_day = 0.d0  !day after spawning

! ---- Coral population dynamics -----------------------------------------------

!$omp parallel private(i,j,n)
!$omp do private(PFD, Mort, Growth, Recr)
      do i=LBi+1,UBi-1
        do j=LBj+1,UBj-1

          PFD = 1000.0d0*exp(-0.12*h(i,j))

!      Coral mortality rate (d-1)
          Mort = 1.0d-4*(1.0d0-tanh(PFD/275.0d0))*P_CORAL(ng)%cover(i,j)
          do n=1,Nstage
            Mort = Mort + COTS_pred(n,i,j)
          end do
          
!      Coral growth rate (d-1)
          
          Growth = Gc*tanh(PFD/275.0d0)*( 1.0d0-P_CORAL(ng)%cover(i,j)/(1.0d0-sand(i,j)) )*P_CORAL(ng)%cover(i,j)
          
          Recr = 1.0d-4*LARcoral(i,j)
          
          P_CORAL(ng)%cover(i,j)=P_CORAL(ng)%cover(i,j)+(   &
          ! Reaction term
     &        + Growth                                      &
     &        - Mort                                        &
     &       ) *dt
          P_CORAL(ng)%cover(i,j) = max (P_CORAL(ng)%cover(i,j), 0.0d0)    !!! For Error handling
        end do
      end do
!$omp end do
!$omp end parallel

      
!---------------------------------------------------------------------
      
      
      RETURN

     END SUBROUTINE p_coral_dynamics

  END MODULE mod_p_coral

