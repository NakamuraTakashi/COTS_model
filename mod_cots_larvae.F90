
!!!=== ver 2016/03/10   Copyright (c) 2015-2016 Takashi NAKAMURA  =====

#include "cppdefs.h"


!!!**** MODULE OF COTS (Crown-of-thorns starfish) MODEL ****************

  module mod_cots_larvae

    implicit none

    TYPE T_COTS_LAR
    
      real(8), pointer :: Xgrid(:)      ! X position of COTS larval (index base)
      real(8), pointer :: Ygrid(:)      ! Y position of COTS larval (index base)
      real(8), pointer :: Zgrid(:)      ! Z position of COTS larval (index base)
      real(8), pointer :: num(:)        ! Number of COTS larval (Individuals/particle)
      real(8), pointer :: mort(:)       ! Mortrality of COTS larval (d-1)
      real(8), pointer :: phy(:)        ! Phytoplankton density around COTS larvae (umolC/L)
      real(8), pointer :: age(:)        ! Age of COTS larvae (days)
      real(8), pointer :: lipid(:)      ! Lipid of COTS larvae (umolC/L)
      integer, pointer :: status(:)     ! Status of COTS larvae. 
                                        !    0: Not started
                                        !    1: Peragic stage1; 
                                        !    2: Peragic stage2; 
                                        !    3: Peragic stage3; 
                                        !    4: settled
                                        !    5: out

    END TYPE T_COTS_LAR

    TYPE (T_COTS_LAR), allocatable :: COTS_LAR(:)

!-----------------------------------------------------------------------

  CONTAINS


!!! **********************************************************************
!!!  set initial conditions for COTS (Crown-of-thorns starfish) model
!!! **********************************************************************

    SUBROUTINE initialize_cots_larvae(ng, Ngrids, N_flt)

      implicit none
! input parameters
      integer, intent(in) :: ng, Ngrids, N_flt

      integer i,j,n

      IF (ng.eq.1) allocate ( COTS_LAR(Ngrids) )

      allocate( COTS_LAR(ng)%Xgrid( N_flt ) )
      allocate( COTS_LAR(ng)%Ygrid( N_flt ) )
      allocate( COTS_LAR(ng)%Zgrid( N_flt ) )
      allocate( COTS_LAR(ng)%num ( N_flt ) )
      allocate( COTS_LAR(ng)%mort( N_flt ) )
      allocate( COTS_LAR(ng)%phy ( N_flt ) )
      allocate( COTS_LAR(ng)%age ( N_flt ) )
      allocate( COTS_LAR(ng)%lipid( N_flt ) )
      allocate( COTS_LAR(ng)%status( N_flt ) )

!----------set data -----------------------
!  Set initial conditions

      COTS_LAR(ng)%Xgrid(:)=0.0d0 ! 
      COTS_LAR(ng)%Ygrid(:)=0.0d0 ! 
      COTS_LAR(ng)%Zgrid(:)=0.0d0 ! 
      COTS_LAR(ng)%num (:)=0.0d0  ! 
      COTS_LAR(ng)%mort(:)=0.0d0  ! 
      COTS_LAR(ng)%phy (:)=0.0d0  ! 
      COTS_LAR(ng)%age (:)=0.0d0  ! 
      COTS_LAR(ng)%lipid(:)=0.0d0  ! 
      COTS_LAR(ng)%status(:)=0    ! 

      RETURN
    END SUBROUTINE initialize_cots_larvae
      
    SUBROUTINE finalize_cots_larvae(ng)

      implicit none
! input parameters
      integer, intent(in) :: ng

      IF (ng.eq.1) deallocate ( COTS_LAR )

!      deallocate( COTS_LAR(ng)%Xgrid )
!      deallocate( COTS_LAR(ng)%Ygrid )
!      deallocate( COTS_LAR(ng)%Zgrid )
!      deallocate( COTS_LAR(ng)%num  )
!      deallocate( COTS_LAR(ng)%mort )
!      deallocate( COTS_LAR(ng)%phy  )
!      deallocate( COTS_LAR(ng)%age  )
!      deallocate( COTS_LAR(ng)%lipid )
!      deallocate( COTS_LAR(ng)%status  )

      RETURN
    END SUBROUTINE finalize_cots_larvae

!!! **********************************************************************
!!!  Main program of COTS (Crown-of-thorns starfish) larvae model 
!!! **********************************************************************

    SUBROUTINE cots_larvae               &
!          input parameters
     &            (ng                    &   ! ng: nested grid number; n: coral compartment; i,j: position
     &            ,LBi, UBi, LBj, UBj    &
     &            ,N_flt, N_flt_time     &
     &            ,dt                    &   ! Time step (day)
     &            ,dx                    &   ! dx: x grid size (m)
     &            ,dy                    &   ! dy: y grid size (m)
     &            ,COTSdens              &   ! COTSdens(LBi:UBi,LBj:UBj) : adult COTS density (indiv./m2)
     &            ,PHY                   &   ! PHY(LBi:UBi,LBj:UBj): Phytoplankton density at surface (umolC/L)
     &            ,nanVal                &   ! Nan value
     &            ,h                     &   ! h(LBi:UBi,LBj:UBj) : water depth (m)
     &            ,coral                 &   ! coral(LBi:UBi,LBj:UBj) : coral coverage (0-1)
     &            ,sand                  &   ! sand(LBi:UBi,LBj:UBj) : sand coverage (0-1)
     &            ,time                  &   ! COTS model time (days since initialization)
!          input & output parameters
     &            ,first_step            &   ! Status 1: first_step
     &            ,cots_rec              &   ! COTS larval recruitment rate (individual m-2 d-1)
     &             )
!
      
      implicit none
      
! input parameters
      integer, intent(in) :: ng
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: N_flt, N_flt_time
      real(8), intent(in) :: dt
      real(8), intent(in) :: dx      
      real(8), intent(in) :: dy      
      real(8), intent(in) :: COTSdens(LBi:UBi,LBj:UBj)
      real(8), intent(in) :: PHY(LBi:UBi,LBj:UBj)
      real(8), intent(in) :: nanVal
      real(8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(8), intent(in) :: coral(LBi:UBi,LBj:UBj)
      real(8), intent(in) :: sand(LBi:UBi,LBj:UBj)
      real(8), intent(in) :: time
      integer, intent(inout) :: first_step
      real(8), intent(inout) :: cots_rec(LBi:UBi,LBj:UBj)

!!!------------Set parameters  ----------------------------------
      real(8), parameter :: iniLipid = 1.0d0
!      real(8), parameter :: num_egg = 5.0d7/30.0d0  ! Number of egg per one COTS
      real(8), parameter :: num_egg = 5.0d7/30.0d0  ! Number of egg per one COTS
      real(8), parameter :: p_ferit = 0.3           ! Probability of fertilization
      real(8), parameter :: phy2chla = 0.4d0      ! Chl-a/C: 0.24d0 gChl-a/molC (Faure et al. 2006) !!!
      
      real(8), parameter :: Vp   = 0.159d0         ! Phytoplankton grazing rate (/chl-a)
      real(8), parameter :: Kp   = 4.55d0           ! Phytoplankton grazing rate (/chl-a)
      real(8), parameter :: resp = 0.0563d0         ! Respiration rate
      real(8), parameter :: Kmor = 4.82d0           ! 
      real(8), parameter :: mor0 = 0.791d0          ! 

!--- Other variables ----------------------------------------------
      integer :: i, j, n
      integer :: xi, eta
      real(8) :: chla           ! 
      real(8), save :: t_days = 0.d0  !day after spawning

      
     
      
! ---- Initial Status laveling -----------------------------------------------

      if (first_step == 1) then
        first_step =0
        do i=1,N_flt
          if(COTS_LAR(ng)%Xgrid(i) < nanVal) then
            COTS_LAR(ng)%status(i)=1
          else
            COTS_LAR(ng)%status(i)=0
          end if
        end do
      end if
      
      
! ---- Calculate number of COTS larvae per one perticle -----------------------------------------------
      
!      cots_rec(:,:) = 0.0d0

!$omp parallel private(i,j)
!$omp do private(xi,eta)
      do i=1,N_flt
      
! ----- Status check -----

        if (COTS_LAR(ng)%status(i)==0) then
        
          if (COTS_LAR(ng)%Xgrid(i) < nanVal) then
!   ----- Initiallize the properties -----
            COTS_LAR(ng)%status(i)=1
            
            xi  = nint(COTS_LAR(ng)%Xgrid(i))
            eta = nint(COTS_LAR(ng)%Ygrid(i))
            
            COTS_LAR(ng)%num(i) = COTSdens(xi,eta) * dx * dy * num_egg * p_ferit
            COTS_LAR(ng)%num(i) = max(COTS_LAR(ng)%num(i), 0.0d0)    !!! Error handling
            COTS_LAR(ng)%age(i) = 0.0d0
            COTS_LAR(ng)%lipid(i) = iniLipid
            
          else
            cycle
          end if
          
          
        else if (COTS_LAR(ng)%status(i)==1) then
          if (COTS_LAR(ng)%age(i)>2.0d0) then
            COTS_LAR(ng)%status(i)=2
          end if
          
        else if (COTS_LAR(ng)%status(i)==2) then
          if (COTS_LAR(ng)%age(i)>5.5d0) then
            COTS_LAR(ng)%status(i)=3
          end if
          
          
        else if (COTS_LAR(ng)%status(i)==3) then
          
          xi  = nint(COTS_LAR(ng)%Xgrid(i))
          eta = nint(COTS_LAR(ng)%Ygrid(i))
          
          if ( COTS_LAR(ng)%age(i) >= 19.0d0 .and.              &
          & h(xi,eta) <= 30.0d0 .and. sand(xi,eta) <= 0.95d0) then  ! settlement
            COTS_LAR(ng)%status(i)=4
            
            cots_rec(xi,eta) = cots_rec(xi,eta) + COTS_LAR(ng)%num(i)/dx/dy
            
            IF(cots_rec(xi,eta)*0.0d0 /= 0.0d0) THEN  !!!---------Error Handling: Check NAN
              write(50,*) xi, eta, COTS_LAR(ng)%num(i)
            END IF

            cycle
          end if
            
        else if (COTS_LAR(ng)%status(i)==4) then
          cycle
        else if (COTS_LAR(ng)%status(i)>=1     .and.            &
     &           COTS_LAR(ng)%Xgrid(i) < nanVal) then
          COTS_LAR(ng)%status(i)=5
          cycle
        end if
        
        xi  = nint(COTS_LAR(ng)%Xgrid(i))
        eta = nint(COTS_LAR(ng)%Ygrid(i))
        COTS_LAR(ng)%phy(i) = PHY(xi,eta)
!        chla = 0.7d0 * phy2chla  !!!!!!!!!!! Test mode ******************************************************
        chla = COTS_LAR(ng)%phy(i) * phy2chla
        
        COTS_LAR(ng)%lipid(i) = COTS_LAR(ng)%lipid(i)       &
     &          +(Vp*chla/(chla+Kp)                         &
     &          - resp                                      &
     &           )*dt
        COTS_LAR(ng)%lipid(i) = max(COTS_LAR(ng)%lipid(i), 0.0d0)    !!! Error handling
!        COTS_LAR(ng)%lipid(i) = min(COTS_LAR(ng)%lipid(i), 1.0d0)  !!! Error handling
        
        COTS_LAR(ng)%mort(i) = mor0 * exp( - Kmor * COTS_LAR(ng)%lipid(i) )
        
        COTS_LAR(ng)%num(i) = COTS_LAR(ng)%num(i)           &
     &          -(COTS_LAR(ng)%mort(i)*COTS_LAR(ng)%num(i)  &
     &           )*dt
        COTS_LAR(ng)%num(i) = max (COTS_LAR(ng)%num(i), 0.0d0)    !!! For Error handling
        
        IF(COTS_LAR(ng)%num(i)*0.0d0 /= 0.0d0) THEN  !!!---------Error Handling: Check NAN
          COTS_LAR(ng)%num(i) = 0.0d0
          write(50,*) COTS_LAR(ng)%num(i), COTS_LAR(ng)%phy(i), COTS_LAR(ng)%mort(i), COTS_LAR(ng)%lipid(i)
        END IF

     
        COTS_LAR(ng)%age(i)=COTS_LAR(ng)%age(i)+dt
        
      end do
!$omp end do
!$omp end parallel

      
      RETURN

     END SUBROUTINE cots_larvae

  END MODULE mod_cots_larvae

