!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE NUMERICAL_FLUXES

  USE GLOBAL
  USE GRID
  USE PDE
  use implicit_electrons_utilities

  implicit none


  CONTAINS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine COMPUTE_other_WAVES(ul,ur,SL,SR,S_star)

  implicit none

  REAL(KIND=8), DIMENSION(4), intent(in) :: UL, UR
  real(kind=8), intent(in) :: SL,SR
  real(kind=8), intent(inout) :: S_star
  real(kind=8)  :: rhol,uxl,uyl,El,Pl
  real(kind=8)  :: rhor,uxr,uyr,Er,Pr

  call cons2primvar(ul,rhol,uxl,uyl,El,Pl)
  call cons2primvar(ur,rhor,uxr,uyr,Er,Pr)

  S_star = ( Pr - Pl + rhol*uxl*(SL - uxl) - rhor*uxr*(SR - uxr)) / &
            (rhol*(Sl - uxl) - rhor*(Sr - uxr))

  end subroutine COMPUTE_other_WAVES

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine cons2primvar(sol,rho,ux,uy,E,P)

    implicit NONE

    REAL(KIND=8), DIMENSION(4), intent(in) :: sol
    real(kind=8), intent(inout) :: rho,ux,uy,E,P
    real(kind=8) :: rhoE,rhoux,rhouy

    rho    = sol(1)
    rhoux  = sol(2)
    rhouy  = sol(3)
    ux     = rhoux/(rho + 1.0e-35)
    uy     = rhouy/(rho + 1.0e-35)
    rhoE   = sol(4)
    E = rhoE/rho
    P      = (GAS_GAMMA_1 - 1)*(rhoE - rhoux**2/rho/2.0 - rhouy**2/rho/2.0)

    end subroutine cons2primvar

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine COMPUTE_NUMERICAL_FLUXES(sol,Flx)

  implicit none

  REAL(KIND=8), DIMENSION(N_EQ, N_cells), INTENT(IN)  :: sol
  REAL(KIND=8), DIMENSION(N_EQ, N_int), INTENT(INOUT) :: Flx
  REAL(KIND=8), DIMENSION(N_EQ):: num_flux

  INTEGER :: i_INT
  REAL(KIND=8), DIMENSION(N_EQ) :: UL, UR, FL, FR

  REAL(KIND=8) :: L_AVE        ! For Lax-Friedrichs
  INTEGER :: IC_LEFT, IC_RIGHT ! For Lax-Friedrichs

  REAL(KIND=8) :: rho, ux, uy, P, T, a, sL, sR, dummy ! For Lax-Friedrichs

  CALL UPDATE_CFL(sol)

  ! Loop on interfaces
  DO I_INT = 1, N_int

    ! Compute interface left and right states
    CALL COMPUTE_INT_LEFT_RIGHT_STATES(sol, I_INT, UL, UR)

    ! Compute PDE flux vectors from the left and right states
    call apply_numerical_scheme(UL,UR,num_flux)

    if(NaN_BOOL)THEN
      print*," NaN found at ",I_int
      STOP
    end if

    flx(:,i_int) = num_flux

  end do

  end subroutine COMPUTE_NUMERICAL_FLUXES

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine apply_numerical_scheme(ul,ur,num_flux)


    implicit NONE

    REAL(KIND=8) :: Smin_L1, Smin_R1, Smax_L1, Smax_R1, Smin1, Smax1, S_star1
    REAL(KIND=8) :: Smin_L2, Smin_R2, Smax_L2, Smax_R2, Smin2, Smax2 , S_star2
    REAL(KIND=8), DIMENSION(N_EQ) :: FL, FR, U_starL, U_starR
    REAL(KIND=8), DIMENSION(N_EQ), intent(in) :: UL, UR
    REAL(KIND=8), DIMENSION(N_EQ), INTENT(INOUT) :: num_flux
    real(kind=8) :: Smin_L, Smax_L,Smin_r, Smax_r
    REAL(KIND=8) :: Smin, Smax, S_star

    CALL COMPUTE_PDE_FLUXES(UL, FL) ! Left flux
    CALL COMPUTE_PDE_FLUXES(UR, FR) ! Right flux


    if (HLL_BOOL) then

        ! Multi fluid
        if (SYS_Euler_MF_BOOL) then

            CALL COMPUTE_MAX_MIN_WAVESPEEDS_MF(UL, Smin_L1, Smax_L1, Smin_L2, Smax_L2)
            CALL COMPUTE_MAX_MIN_WAVESPEEDS_MF(UR, Smin_R1, Smax_R1, Smin_R2, Smax_R2)

            Smin1 = MIN(Smin_L1, Smin_R1)
            Smax1 = MAX(Smax_L1, Smax_R1)

            smin2 = MIN(Smin_L2, Smin_R2)
            Smax2 = MAX(Smax_L2, Smax_R2)



            ! COMPUTE flux for FIRST FLUID
            IF (Smin1 .GT. 0.0) THEN
              num_flux(1:4) = FL(1:4)
            ELSE IF (Smax1 .LT. 0.0) THEN
              num_flux(1:4) = FR(1:4)
            ELSE
              num_flux(1:4) = (Smax1*FL(1:4) - Smin1*FR(1:4) + Smin1*Smax1*(UR(1:4) - UL(1:4)))/(Smax1-Smin1 + 1.e-15)
            END IF

            ! COMPUTE flux for SECOND FLUID
            IF (Smin2 .GT. 0.0) THEN
              num_flux(5:8) = FL(5:8)
            ELSE IF (Smax2 .LT. 0.0) THEN
              num_flux(5:8) = FR(5:8)
            ELSE
              num_flux(5:8) = (Smax2*FL(5:8) - Smin2*FR(5:8) + Smin2*Smax2*(UR(5:8) - UL(5:8)))/(Smax2-Smin2 + 1.e-15)
            END IF

        ! Single Fluid
        ELSE

              CALL COMPUTE_MAX_MIN_WAVESPEEDS(UL, Smin_L, Smax_L)
              CALL COMPUTE_MAX_MIN_WAVESPEEDS(UR, Smin_R, Smax_R)

              Smin = MIN(Smin_L, Smin_R)
              Smax = MAX(Smax_L, Smax_R)



              ! COMPUTE flux now
              IF (Smin .GT. 0.0) THEN
                num_flux(:) = FL
              ELSE IF (Smax .LT. 0.0) THEN
                num_flux(:) = FR
              ELSE
                num_flux(:) = (Smax*FL - Smin*FR + Smin*Smax*(UR - UL))/(Smax-Smin + 1.e-15)
              END IF

        end if

    else if(HLLC_BOOL) then

      ! Multi fluid
      if (SYS_Euler_MF_BOOL) then

        CALL COMPUTE_MAX_MIN_WAVESPEEDS_MF(UL, Smin_L1, Smax_L1, Smin_L2, Smax_L2)
        CALL COMPUTE_MAX_MIN_WAVESPEEDS_MF(UR, Smin_R1, Smax_R1, Smin_R2, Smax_R2)

        Smin1 = MIN(Smin_L1, Smin_R1)
        Smax1 = MAX(Smax_L1, Smax_R1)

        smin2 = MIN(Smin_L2, Smin_R2)
        Smax2 = MAX(Smax_L2, Smax_R2)

        CALL COMPUTE_other_WAVES(UL(1:4),UR(1:4), Smin1, Smax1,S_star1)
        CALL COMPUTE_other_WAVES(UL(5:8),UR(5:8), Smin2, Smax2,S_star2)

        call compute_Ustar_left_right_MF(UL,Smin1,S_star1,Smin2,S_star2,U_starL)
        call compute_Ustar_left_right_MF(UR,Smax1,S_star1,Smax2,S_star2,U_starR)

        ! COMPUTE flux for FIRST FLUID
        IF (Smin1 .Ge. 0.0) THEN
          num_flux(1:4) = FL(1:4)
        else if((Smin1 .le. 0.0) .and. (S_star1 .ge. 0.0)) then
          num_flux(1:4) = FL(1:4) + Smin1 * (U_starL(1:4) - UL(1:4))
        else if((Smax1 .ge. 0.0) .and. (S_star1 .le. 0.0)) then
          num_flux(1:4) = FR(1:4) + Smax1 * (U_starR(1:4) - UR(1:4))
        else if((Smax1 .le. 0.0)) then
          num_flux(1:4) = FR(1:4)
        else
          print*, " FATAL ERROR! EXPLICIT METHOD COULD NOT FIND RIGHT HLLC FLUX "
          NaN_BOOL = .True.
        end if


        ! COMPUTE flux for SECOND FLUID
        IF (Smin2 .Ge. 0.0) THEN
          num_flux(5:8) = FL(5:8)
        else if((Smin2 .le. 0.0) .and. (S_star2 .ge. 0.0)) then
          num_flux(5:8) = FL(5:8) + Smin2 * (U_starL(5:8) - UL(5:8))
        else if((Smax2 .ge. 0.0) .and. (S_star2 .le. 0.0)) then
          num_flux(5:8) = FR(5:8) + Smax2 * (U_starR(5:8) - UR(5:8))
        else if((Smax2 .le. 0.0)) then
          num_flux(5:8) = FR(5:8)
        else
          print*, " FATAL ERROR! EXPLICIT METHOD COULD NOT FIND RIGHT HLLC FLUX "
          NaN_BOOL = .True.
        end if

      ELSE

        CALL COMPUTE_MAX_MIN_WAVESPEEDS(UL, Smin_L, Smax_L)
        CALL COMPUTE_MAX_MIN_WAVESPEEDS(UR, Smin_R, Smax_R)

        Smin = MIN(Smin_L, Smin_R)
        Smax = MAX(Smax_L, Smax_R)

        CALL COMPUTE_other_WAVES(UL,UR, Smin, Smax,S_star)


        call compute_Ustar_left_right(UL,Smin,S_star,U_starL)
        call compute_Ustar_left_right(UR,Smax,S_star,U_starR)

        ! COMPUTE flux for FIRST FLUID
        IF (Smin .Ge. 0.0) THEN
          num_flux(1:4) = FL(1:4)
        else if((Smin .le. 0.0) .and. (S_star .ge. 0.0)) then
          num_flux(1:4) = FL(1:4) + Smin * (U_starL(1:4) - UL(1:4))
        else if((Smax .ge. 0.0) .and. (S_star .le. 0.0)) then
          num_flux(1:4) = FR(1:4) + Smax * (U_starR(1:4) - UR(1:4))
        else if((Smax .le. 0.0)) then
          num_flux(1:4) = FR(1:4)
        else
          print*, " FATAL ERROR! IMPLICIT METHOD COULD NOT FIND RIGHT HLLC FLUX "
          NaN_BOOL = .True.
        end if

      end if

    end if


    end subroutine apply_numerical_scheme

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  IMPLICIT STUFF  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE COMPUTE_JACOBIAN(SOL,F_N)

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ, N_cells), INTENT(IN)  :: sol
  REAL(KIND=8), DIMENSION(N_EQ, N_int), INTENT(IN) :: F_N
  INTEGER :: i_INT
  REAL(KIND=8), DIMENSION(N_EQ) :: UL, UR

  ! Fix first interface assembling
  call first_interface(sol)

  ! Loop on other interfaces
  DO I_INT = 2, N_int-1

    ! Compute interface left and right states
    CALL COMPUTE_INT_LEFT_RIGHT_STATES(sol, I_INT, UL, UR)

    call assemble_tensor_flux(I_INT,UL,UR,F_N(:,I_INT-1:I_INT+1))

  END do

  call last_interface(sol,F_N(:,N_int-1:N_int))

  END SUBROUTINE COMPUTE_JACOBIAN

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine first_interface(sol)

  implicit none

  REAL(KIND=8), DIMENSION(N_EQ, N_cells), INTENT(IN)  :: sol
  REAL(KIND=8), DIMENSION(N_EQ) :: UL, UR
  INTeGER :: cons_var,i_INT
  REAL(KIND=8), DIMENSION(N_EQ) :: pert_r
  REAL(KIND=8), DIMENSION(N_EQ) :: PERTURBED_CONS_r,num_flx_l,num_flx_r

  i_int = 1

  CALL COMPUTE_INT_LEFT_RIGHT_STATES(sol, I_INT, UL, UR)

  call COMPUTE_PERTURBED_CONS(UR,PERTURBED_CONS_R)

  do CONS_VAR=1,4

    call apply_perturbation(ur,perturbed_cons_R,pert_R,cons_var)

    call apply_numerical_scheme(ul,pert_r,num_flx_r)

    INTER_FLUX(:,cons_var) = num_flx_r

  end do

  end subroutine first_interface


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine last_interface(sol,num_flx)

  implicit none

  REAL(KIND=8), DIMENSION(N_EQ, N_cells), INTENT(IN)  :: sol
  REAL(KIND=8), DIMENSION(N_EQ,2),intent(in) :: num_flx
  REAL(KIND=8), DIMENSION(N_EQ) :: UL, UR
  INTeGER :: cons_var,i_INT,random_i
  REAL(KIND=8), DIMENSION(N_EQ) :: pert_l
  REAL(KIND=8), DIMENSION(N_EQ) :: PERTURBED_CONS_l,num_flx_l,num_flx_r
  REAL(KIND=8), DIMENSION(4,4) :: Tensor_Flux_Center

  i_int = N_int

  CALL COMPUTE_INT_LEFT_RIGHT_STATES(sol, I_INT, UL, UR)

  call COMPUTE_PERTURBED_CONS(UL,PERTURBED_CONS_L)


  do CONS_VAR=1,4

    call apply_perturbation(ul,perturbed_cons_l,pert_l,cons_var)

    call apply_numerical_scheme(pert_l,ur,num_flx_l)

      do random_i = 1,4

      Tensor_Flux_Center(random_i,cons_var) = ((num_flx_l(random_i) - INTER_FLUX(random_i,cons_var))- &
                                            (num_flx(random_i,2) - num_flx(random_i,1))) / &
                                          ((ul(cons_var) + perturbed_cons_l(cons_var)) - ul(cons_var) + 1d-35)
      end do

  end do

  call Insert_in_Final_Tensor(Tensor_Flux_Center,i_int - 1,i_int - 1)


  end subroutine last_interface


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine assemble_tensor_flux(inter_index,ul,ur,num_flx)

  implicit NONE

  REAL(KIND=8), DIMENSION(N_EQ,3),intent(in) :: num_flx
  REAL(KIND=8), DIMENSION(N_EQ),intent(in) :: UL, UR
  integer, intent(in) :: inter_index
  REAL(KIND=8), DIMENSION(N_EQ) :: pert_r,pert_l
  REAL(KIND=8), DIMENSION(N_EQ) :: PERTURBED_CONS_L,PERTURBED_CONS_r,num_flx_l,num_flx_r
  integer :: cons_var,col_index,row_index,random_i
  REAL(KIND=8), DIMENSION(4,4) ::Tensor_Flux_Down,Tensor_Flux_Right,Tensor_Flux_Center ! Fluxes rows - Perturbed conservative cols


  ! GC *------1-----|(1)-------2--------|(2)-------3---------|(3)-----4----|(N_int)-----N_cells-----*GC

  ! Subroutine that returns conservative variables perturbed
  call COMPUTE_PERTURBED_CONS(UL,PERTURBED_CONS_L)
  call COMPUTE_PERTURBED_CONS(UR,PERTURBED_CONS_R)

  !-------------(UL)----------|(inter_index)----------------(UR)---------
  ! Note that subroutine still does not know if the interface is the RHS or
  ! the LHS of the cell_index


  ! Before computing the Jacobian I need to retrieve the numerical scheme for fluxes:
  ! indeed perturbations concern electrons conservative variables only

    do CONS_VAR=1,4

      call apply_perturbation(ul,perturbed_cons_l,pert_l,cons_var)
      call apply_perturbation(ur,perturbed_cons_R,pert_R,cons_var)

      ! Computing fluxes as function of perturbed state according to the
      ! cell (L - Left; R - Right ) the perturbation has been computed
      call apply_numerical_scheme(pert_l,ur,num_flx_l)
      call apply_numerical_scheme(ul,pert_r,num_flx_r)

      do random_i = 1,4

        ! Indeed we are interested in electrons population perturbation only
        Tensor_Flux_Right(random_i,cons_var) = ((num_flx_r(random_i) - num_flx(random_i,1))- &
                                              (num_flx(random_i,2) - num_flx(random_i,1))) / &
                                            ((ur(cons_var) + perturbed_cons_r(cons_var)) - ur(cons_var) + 1d-35)

        Tensor_Flux_Down(random_i,cons_var) = ((num_flx(random_i,3) - num_flx_l(random_i))- &
                                              (num_flx(random_i,3) - num_flx(random_i,2))) / &
                                            ((ul(cons_var) + perturbed_cons_l(cons_var)) - ul(cons_var) + 1d-35)

        Tensor_Flux_Center(random_i,cons_var) = ((num_flx_l(random_i) - INTER_FLUX(random_i,cons_var))- &
                                              (num_flx(random_i,2) - num_flx(random_i,1))) / &
                                            ((ul(cons_var) + perturbed_cons_l(cons_var)) - ul(cons_var) + 1d-35)


      end do

      INTER_FLUX(:,cons_var) = num_flx_r

   end do


   call Insert_in_Final_Tensor(Tensor_Flux_Right,inter_index - 1,inter_index)
   call Insert_in_Final_Tensor(Tensor_Flux_Down,inter_index,inter_index-1)
   call Insert_in_Final_Tensor(Tensor_Flux_Center,inter_index - 1,inter_index - 1)



end subroutine assemble_tensor_flux


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE COMPUTE_INT_LEFT_RIGHT_STATES(sol, I_INT, UL, UR)

    ! This subroutine computes the left state UL and the right state UR.
    ! If first-order is requested, it just takes the value in the left and right cells.
    ! If second-order is requested, it follows the MUSCL approach.

    IMPLICIT NONE

    INTEGER,                                INTENT(IN)  :: I_INT
    REAL(KIND=8), DIMENSION(N_EQ, N_cells), INTENT(IN)  :: sol
    REAL(KIND=8), DIMENSION(N_EQ),          INTENT(OUT) :: UL, UR

    INTEGER :: IC_LEFT, IC_RIGHT, I_EQ
    REAL(KIND=8) :: TH_LEFT, TH_RIGHT

    ! I will always need these
    CALL LEFT_CELL_FROM_INT(I_INT, IC_LEFT)   ! Ask to the grid
    CALL RIGHT_CELL_FROM_INT(I_INT, IC_RIGHT) ! Ask to the grid


    ! =========== FIRST ORDER SCHEME - NO RECONSTRUCTION =============
    IF (RECONSTRUCTION_ORDER .EQ. 0) THEN

      UL = sol(:,IC_LEFT)  ! Left state, order 1
      UR = sol(:,IC_RIGHT) ! Right state, order 1

    ! =========== SECOND ORDER SCHEME - LINEAR RECONSTRUCTION ========
    ELSE IF (RECONSTRUCTION_ORDER .EQ. 1) THEN

      ! +++++++ First and last interfaces don't have enough cells. Just use first order.
      IF ((I_INT .EQ. 1) .OR. (I_INT .EQ. N_INT) ) THEN

        UL = sol(:,IC_LEFT)  ! Left state, order 1
        UR = sol(:,IC_RIGHT) ! Right state, order 1

        RETURN

      ! +++++++ Internal interfaces
      ELSE

        ! Loop on equations and do reconstruction
        DO I_EQ = 1, N_EQ

          ! ---- Left state
          !
          TH_LEFT =   (sol(I_EQ, IC_LEFT + 1) - sol(I_EQ, IC_LEFT))/(x_cc(IC_LEFT+1) - x_cc(IC_LEFT)) &
                    * (x_cc(IC_LEFT) - x_cc(IC_LEFT-1))/(sol(I_EQ, IC_LEFT) - sol(I_EQ, IC_LEFT-1) + 1.0e-15)

          UL(I_EQ) = sol(I_EQ, IC_LEFT) + (x_int(I_INT) - x_cc(IC_LEFT))*(sol(I_EQ, IC_LEFT) - sol(I_EQ, IC_LEFT-1))&
                                       /(x_cc(IC_LEFT) - x_cc(IC_LEFT - 1))*LIMITER_FUN(TH_LEFT)


          ! ---- Right state
          TH_RIGHT =  (sol(I_EQ, IC_RIGHT + 1) - sol(I_EQ, IC_RIGHT))/(x_cc(IC_RIGHT+1) - x_cc(IC_RIGHT)) &
                    * (x_cc(IC_RIGHT) - x_cc(IC_RIGHT-1))/(sol(I_EQ, IC_RIGHT) - sol(I_EQ,IC_RIGHT-1) + 1.0e-15)

          UR(I_EQ) = sol(I_EQ, IC_RIGHT) - (x_cc(IC_RIGHT) - x_int(I_INT))*(sol(I_EQ, IC_RIGHT+1) - sol(I_EQ, IC_RIGHT))&
                                        /(x_cc(IC_RIGHT+1) - x_cc(IC_RIGHT))*LIMITER_FUN(1/(TH_RIGHT+1.0e-15))


        END DO

      END IF ! Internal or border interfaces

    END IF ! End second order

  END SUBROUTINE COMPUTE_INT_LEFT_RIGHT_STATES

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION LIMITER_FUN(THETA)

    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN)  :: THETA
    REAL(KIND=8) :: LIMITER_FUN

    LIMITER_FUN = (THETA**2 + THETA)/(THETA**2 + 1) ! Van Albada
    ! LIMITER_FUN = (THETA + ABS(THETA))/(1 + ABS(THETA)) ! Van Leer

    RETURN

  END FUNCTION LIMITER_FUN

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine UPDATE_CFL(sol)

  implicit none

  REAL(KIND=8), DIMENSION(N_EQ, N_cells), INTENT(IN)  :: sol
  INTEGER :: IC
  REAL(KIND=8) :: Smin1, Smax1 ! For HLL
  REAL(KIND=8) :: Smin2, Smax2 ! For HLL
  real(kind = 8):: eig_max_ions,eig_max_electrons,cell_max_ions,cell_max_electrons
  REAL(KIND=8) :: Smin, Smax,eig_max,cell_max,dx_ions,dx_electrons


  eig_max_electrons = 0.0d0
  eig_max_ions = 0.0d0
  CFL_MAX = 0.0d0
  dx_electrons=1d6
  dx_ions=1d6

  ! Multi Fluid
  if (SYS_Euler_MF_BOOL) then
  ! Loop on cells and find the maximum CFL (approximated maybe)
    DO IC = 1, N_cells

      CALL COMPUTE_MAX_MIN_WAVESPEEDS_MF(sol(:,IC), sMIN1, sMAX1, sMIN2, sMAX2)

      if(MAX(ABS(sMIN1), ABS(sMAX1))/L_cells(IC) .gt. eig_max_electrons/dx_electrons)THEN

        eig_max_electrons = MAX(ABS(sMIN1), ABS(sMAX1) )
        dx_electrons = L_cells(IC)

      end if

      if(MAX(ABS(sMIN2), ABS(sMAX2))/L_cells(IC) .gt. eig_max_ions/dx_ions)THEN

        eig_max_ions = MAX(ABS(sMIN2), ABS(sMAX2) )
        dx_ions = L_cells(IC)

      end if


    END DO

    if(ELECTRONS_IMPLCIT_BOOL)THEN

      ! dt decided by maximum CFL desired for electrons implicit physical iteration
      dt =  CFL_PHYS_MAX * dx_electrons / eig_max_electrons

      CFL_max_ions = dt/dx_ions * eig_max_ions

      if(CFL_max_ions .gt. CFL_INP_IONS)THEN

        dt = CFL_INP_IONS * dx_ions / eig_max_ions

      end if

    else

      ! dt decided by electrons dynamic in general if simulation is explicit
      dt =  CFL_INP_electrons * dx_electrons / eig_max_electrons

    end if

    CFL_max_electrons = dt/dx_electrons * eig_max_electrons
    CFL_max_ions = dt/dx_ions * eig_max_ions


    ! Single Fluid
    ELSE if (SYS_EULER_BOOL .and. FLUID_ELECTRONS)then

      if(ELECTRONS_IMPLCIT_BOOL)THEN


              ! Loop on cells and find the maximum CFL (approximated maybe)
              DO IC = 1, N_cells

                CALL COMPUTE_MAX_MIN_WAVESPEEDS(sol(:,IC), sMIN, sMAX)

                eig_max_electrons = MAX(MAX( ABS(sMIN), ABS(sMAX) ),eig_max_electrons)

                ! Local time-stepping for pseudo-iteration
                d_tau_vect(IC) = NUM_CFL_EVOLUTION_FACTOR * L_cells(IC) / MAX( ABS(sMIN), ABS(sMAX) )


              END DO

              ! Global time-stepping for pseudo-iteration
              !d_tau = NUM_CFL_EVOLUTION_FACTOR * L_cells(20) / eig_max_electrons


      else

            ! Loop on cells and find the maximum CFL (approximated maybe)
            DO IC = 1, N_cells

              CALL COMPUTE_MAX_MIN_WAVESPEEDS(sol(:,IC), sMIN, sMAX)

              eig_max_electrons = MAX(MAX(ABS(sMIN), ABS(sMAX) ),eig_max_electrons)

            END DO

            ! dt decided by electrons dynamic in case of explicit simulation
            dt =  CFL_INP_electrons * L_cells(20) / eig_max_electrons

            CFL_max_electrons = CFL_INP_electrons

      end if

else if (SYS_EULER_BOOL .and. FLUID_IONS)THEN

            if(ELECTRONS_IMPLCIT_BOOL)THEN

            ELSE

              ! Loop on cells and find the maximum CFL (approximated maybe)
              DO IC = 1, N_cells

                CALL COMPUTE_MAX_MIN_WAVESPEEDS(sol(:,IC), sMIN, sMAX)

                eig_max_ions = MAX(MAX(ABS(sMIN), ABS(sMAX) ),eig_max_ions)

              END DO

              ! dt decided by electrons dynamic in case of explicit simulation
              dt =  CFL_INP_ions * L_cells(20) / eig_max_ions

              CFL_max_ions = CFL_INP_IONS

            end if

end if



  end subroutine UPDATE_CFL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE NUMERICAL_FLUXES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
