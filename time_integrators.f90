
MODULE TIME_INTEGRATORS

  USE GLOBAL
  USE PDE
  USE NUMERICAL_FLUXES
  use poisson_solver
  use implicit_electrons_utilities
  !use omp_lib

  implicit none

  integer :: activate
  real(kind=8) :: save_res
  logical :: BOOL_KRYLOV

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  CONTAINS




  !MANGOTURB

  SUBROUTINE save_E_TO_FILE(E_x_ID)

  IMPLICIT NONE

  REAL(KIND=8), INTENT(IN) :: E_x_ID
  CHARACTER(LEN=512) :: file_name
  INTEGER :: IC

  OPEN(12399, FILE='./write_E/E.dat', position ="append", FORM="formatted")

  WRITE(12399, *) E_x_ID

  CLOSE(12399)

  END SUBROUTINE save_E_TO_FILE


  SUBROUTINE COMPUTE_TIMESTEP(sol)
  ! Computes timestep by calling the requested integrator

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ, N_cells), INTENT(INOUT) :: sol

  ! Perform only one of the following
  IF (FORWARD_EULER_BOOL) THEN
    CALL FORWARD_EULER(sol)
  ELSE IF (MIDPOINT_EULER_BOOL) THEN
    CALL MIDPOINT_EULER(sol)
  ELSE IF (ELECTRONS_IMPLCIT_BOOL) THEN
    call IMPLICIT_ELECTRONS(SOL)
  ELSE
    WRITE(*,*) "ATTENTION! TIME INTEGRATOR BOOL NOT RECOGNIZED! ABORTING!"
    STOP
  END IF

  END SUBROUTINE COMPUTE_TIMESTEP

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE FORWARD_EULER(sol)
  ! The subroutine updates the solution "sol" using the Forward Euler (Explicit Euler)
  ! scheme. The solution "sol" is flagged as INOUT, since its value is updated from the previous value.

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ, N_cells), INTENT(INOUT) :: sol
  REAL(KIND=8), DIMENSION(n_cells) :: E_SC
  INTEGER :: C_ID, I_INT, int_left,int_right
  REAL(KIND=8), DIMENSION(N_EQ) :: PDE_SRC,PDE_SRC_EM,PDE_SRC_COLL
  REAL(KIND=8), DIMENSION(n_cells) :: Ey_field,b_field
  REAL(KIND=8), DIMENSION(N_EQ) :: dt_vect
  REAL(KIND=8), DIMENSION(N_EQ, N_int) :: F_n ! Fluxes at timestep n
  real(kind=8), dimension(3) :: coll_vect
  real(kind=8) :: p

  PDE_SRC = 0.0d0
  PDE_SRC_EM = 0.0d0
  PDE_SRC_COLL = 0.0d0

  E_SC = 0.0d0


  ! Compute fluxes at interfaces and save them into F_int
  CALL COMPUTE_NUMERICAL_FLUXES(sol, F_n)


  ! MANGOTURB
  ! Update self-consistent Electromagnetic field
  call compute_SC_EM_field(sol,E_SC(2:N_cells-1),Ey_field,B_field)


    ! Update solution (internal cells only, no ghost cells)
    DO C_ID = 2, N_cells-1

      CALL COMPUTE_PDE_SOURCES(sol(:,C_ID), E_SC(C_ID), Ey_field(C_ID), B_field(C_ID), PDE_SRC_EM, &
          PDE_SRC_COLL,C_ID)! Compute sources in the cell

      coll_elastic_vect_el(C_ID) = nu_implicit_el
      coll_elastic_vect_io(C_ID) = nu_implicit_io
      nu_ions_vect(1,C_ID) = nu_implicit_io
      nu_ions_vect(2,C_ID) = nu_implicit_el_iz
      nu_ions_vect(3,C_ID) = nu_implicit_el_ex


      PDE_SRC = PDE_SRC_EM + PDE_SRC_COLL

      INT_LEFT  = C_ID - 1
      INT_RIGHT = C_ID
      sol(:,C_ID) = sol(:,C_ID) + dt/(L_cells(C_ID))*(F_n(:, INT_LEFT) - F_n(:, INT_RIGHT)) &
                                + dt*PDE_SRC

      call check_Temperature(sol(:,C_ID))

      call check_non_physical_points(sol)

      call update_residual(F_n(:, INT_LEFT),F_n(:, INT_RIGHT),PDE_SRC,L_cells(C_ID))

   END DO

  call compute_residual()

  ! Enforcing boundary conditions
  CALL IMPOSE_BC(SOL)


  END SUBROUTINE FORWARD_EULER

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE COMPUTE_DUAL_TIMESTEP(sol)

  ! Implicit method applied to pseudo-time in resolving electrons, dt is the ions one.
  ! Ions are updated explicitly

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ, N_cells), INTENT(INOUT) :: sol
  real(kind=8), dimension(n_eq,N_cells):: sol_auxiliary,sol_2
  INTeGER :: I,index
  i = 0

  NUM_CFL_EVOLUTION_FACTOR = CFL_NUM_EL

  CFL_max_electrons = 0.0d0
  CFL_max_ions = 0.0d0

  ! At each physical (ions) timestep, goal is to reach convergence (electrons residual) solving
  ! M times a linear system, in which the electrons iterating variable Q_el_star changes at every
  ! iteration

  BOOL_KRYLOV = .False.

  ! Saving solution
  sol_auxiliary = sol(1:8,:)

  ! Inlcuding ions density at the n step
  call set_ions_primitive(sol(5,:),n_background)

  call set_Q_ns(sol_auxiliary)

  call set_MF_CFL(sol(1:8,:))

  ! Call Poisson to try and predict Electric Field at timestep n+1
  !call Poisson_prediction(sol_auxiliary,sol_n_minus1)

  ! In the pseudo-time iterations ions value are fixed
  CALL SET_PSEUDO_ITERATION

  res_rho_el = 1000000000000000000000.0d0

  do while (log10(res_rho_el) .gt. toll)

    res_rho_el = 0.0d0

    I = I + 1

    !call Poisson_prediction(sol_auxiliary,sol_n_minus1)

    !sol_n_minus1 = sol_auxiliary

    ! Pseudo-Time iteration
    call IMPLICIT_ELECTRONS(sol_auxiliary(1:4,:))

    NUM_CFL_EVOLUTION_FACTOR = evolution(I,log10(res_rho_el))



    !-------------------------- DEBUG ---------------------------------- !

    IF (DEBUG_BOOL) THEN

      if(mod(I,10) .eq. 0)then

        print*," Residual: ",log10(res_rho_el)
        print*," Minimum d_tau :",minval(d_tau_vect)
        print*," Maximum d_tau :",maxval(d_tau_vect)
        print*," CFL : ",NUM_CFL_EVOLUTION_FACTOR
      end if

    end if
      !------------------------------------------------------------------!


  end do


  !-------------------------- DEBUG ---------------------------------- !
  IF (DEBUG_BOOL) THEN

      PRINT*, "-------- CONVERGENCE UPDATE AT END OF SUB-ITERATIONS --------------"

          WRITE(*,*) "Density residuals :"
          WRITE(*,*) "Delta_Q :",log10(res_rho_el)
          WRITE(*,*) "N of subiterations :",I
          write(*,*) "Max numerical electrons CFL : ", NUM_CFL_EVOLUTION_FACTOR

     PRINT*, "--------------------------------------------"

  end if
  !------------------------------------------------------------------!

  call SET_PHYSICAL_ITERATION

  call set_electrons_primitive(sol(1:4,:))

  sol_n_minus1 = sol

  ! Iterate ions explicitly
  call MIDPOINT_EULER(sol(5:8,:))

  ! In the pseudo-time iterations ions value are fixed
  CALL SET_PSEUDO_ITERATION

  ! Reintroducing electrons solutions
  sol(1:4,:) = sol_auxiliary(1:4,:)

  END SUBROUTINE COMPUTE_DUAL_TIMESTEP


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_MF_CFL(sol)

  implicit NONE

  REAL(KIND=8), DIMENSION(8, N_cells), INTENT(INOUT) :: sol

  SYS_EULER_BOOL = .False.
  SYS_EULER_MF_BOOL = .True.
  N_EQ = 8

  ELECTRONS_IMPLCIT_BOOL = .True.

  call UPDATE_CFL(sol)

  end subroutine set_MF_CFL

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function evolution(iii,res1) result(out)

    implicit none

    real(kind=8) :: out
    integer, intent(in) :: iii
    real(kind=8), intent(in) :: res1

    if(iii .eq. 1)then

      out = CFL_NUM_EL
      save_res = res1

    else if (iii .gt. 1) THEN

      if((res1 .le. save_res) .and. (.not. BOOL_KRYLOV))then

        BOOL_KRYLOV = .True.

      end if

      if(.not. BOOL_KRYLOV)THEN

        out = CFL_NUM_EL
        save_res = res1

      else if (BOOL_KRYLOV) THEN

        out = alfa_1 * ( 1/(abs(res1 - save_res)))**alfa_2
        save_res = res1

      end if

    end if

  end function evolution

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE SET_PHYSICAL_ITERATION

  IMPLICIT NONE

  SYS_EULER_BOOL = .True.
  FLUID_ELECTRONS = .False.
  FLUID_IONS = .True.
  N_EQ = 4
  GAS_M=GAS_M_2
  GAS_Q=GAS_Q_2
  GAS_GAMMA=GAS_GAMMA_2




  END SUBROUTINE SET_PHYSICAL_ITERATION

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE SET_INITIAL_ITERATION

  IMPLICIT none

  !Non c'è bisogno dal momento in cui sono gia input alla prima iterazione
  !SYS_EULER_MF_BOOL = .True.
  !SYS_EULER_BOOL = .False.
  !N_EQ = 8


  ELECTRONS_IMPLCIT_BOOL = .False.

  END SUBROUTINE SET_INITIAL_ITERATION

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE SET_PSEUDO_ITERATION

  IMPLICIT NONE

  if (mod(activate,2) .eq. 0) then

    SYS_EULER_MF_BOOL = .False.
    SYS_EULER_BOOL = .True.
    FLUID_ELECTRONS = .True.
    FLUID_IONS = .False.
    N_EQ = 4
    GAS_M=GAS_M_1
    GAS_Q=GAS_Q_1
    GAS_GAMMA=GAS_GAMMA_1

  else

    SYS_EULER_MF_BOOL = .True.
    SYS_EULER_BOOL = .False.
    N_EQ = 8

  end if



  activate = activate + 1

  ELECTRONS_IMPLCIT_BOOL = .True.

  end SUBROUTINE SET_PSEUDO_ITERATION

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_Q_ns(sol_el)

  implicit none

  REAL(KIND=8), DIMENSION(N_eq, N_cells), intent(inout) :: sol_el
  integer :: iii,row_index,col_index
  REAL(KIND=8), DIMENSION(N_eq,N_cells) :: auxiliary_sol
  REAL(KIND=8), DIMENSION(4,N_cells) :: aux_small_1,aux_small_2
  real(kind=8) :: auxiliary_dt

  if (dual_iteration .eq. 0) then

    ! dt for first 2 explicit steps for electrons
    call SET_INITIAL_ITERATION

    auxiliary_sol = sol_el

    call FORWARD_EULER(auxiliary_sol)

    ! Saving n-1 iterations from explicit step
    sol_n_minus1(1:8,:) = sol_el(1:8,:)

    aux_small_1 = sol_el(1:4,:)

      do iii = 1,n_size

        row_index = mod(iii-1,4) + 1
        col_index = (iii-1)/4 + 2

        Q_n(iii) = aux_small_1(row_index,col_index)

      end do

    sol_el(1:8,:) = auxiliary_sol(1:8,:)

    dual_iteration = 1

  ELSE

      aux_small_1 = sol_el(1:4,:)

      do iii = 1,n_size

        row_index = mod(iii-1,4) + 1
        col_index = (iii-1)/4 + 2

        Q_n(iii) = aux_small_1(row_index,col_index)

      end do

  end if


  end subroutine set_q_ns

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE IMPLICIT_ELECTRONS(sol)
  ! The subroutine updates the solution "sol" using implicit scheme for electrons
  ! and usual explicit scheme for Ions. Note that this method is available for
  ! multifluid simulations only.

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ, N_cells), INTENT(INOUT) :: sol
  REAL(KIND=8), DIMENSION(n_cells) :: E_SC
  INTEGER :: C_ID, I_INT
  REAL(KIND=8), DIMENSION(N_EQ) :: PDE_SRC_EM, PDE_SRC_COLL
  REAL(KIND=8), DIMENSION(N_EQ, N_cells-2) :: delta_cons
  REAL(KIND=8), DIMENSION(n_cells) :: Ey_field,b_field
  INTEGER :: INT_LEFT, INT_RIGHT
  REAL(KIND=8), DIMENSION(N_EQ) :: Electrons_RHS
  REAL(KIND=8), DIMENSION(N_EQ, N_int) :: F_n ! Fluxes at timestep n
  REAL(KIND=8), DIMENSION(3) :: coll_vect
  REAL(KIND=8), DIMENSION(n_cells) :: En_plus_1
  real(kind=8) :: p

  PDE_SRC_EM = 0.0d0
  PDE_SRC_COLL = 0.0d0

  E_SC = 0.0d0
  En_plus_1 = 0.0d0
  !En_plus_1(2:N_cells-1) = Predicted_E(2:N_cells-1)

  call set_Iterative_Method_parameters

  ! Compute fluxes at interfaces and save them into F_int
  CALL COMPUTE_NUMERICAL_FLUXES(sol, F_n)

  INTER_FLUX = 0.0D0

  ! Compute Jacobian
  call COMPUTE_JACOBIAN(SOL,F_N)

  ! Update self-consistent Electromagnetic field
  call compute_SC_EM_field(sol,E_SC(2:N_cells-1),Ey_field,B_field)


  ! Update solution (internal cells only, no ghost cells)
  DO C_ID = 2, N_cells-1

    CALL COMPUTE_PDE_SOURCES(sol(:,C_ID), E_SC(C_ID), Ey_field(C_ID),& !En_plus_1(C_ID)
                                    B_field(C_ID), PDE_SRC_EM, PDE_SRC_COLL,C_ID) ! Compute sources in the cell

    ! Recomputing collisonal frequencies at every subiterations for ions population only
    coll_elastic_vect_el(C_ID) = nu_implicit_el
    nu_ions_vect(2,C_ID) = nu_implicit_el_iz

    INT_LEFT  = C_ID - 1
    INT_RIGHT = C_ID

    ! Retrieving RHS of Final Linear System

    Electrons_RHS= ((F_n(1:4, INT_LEFT) - F_n(1:4, INT_RIGHT))+ (PDE_SRC_EM(1:4) + PDE_SRC_COLL(1:4))*L_cells(C_ID))

    ! Source Jacobian is assembled
    call UPDATE_IMPLICIT_TERMS(Electrons_RHS ,L_cells(C_ID),&
              E_SC(C_ID), Ey_field(C_ID), B_field(C_ID),C_ID,sol(:,C_ID),PDE_SRC_COLL,PDE_SRC_EM) !En_plus_1(C_ID)


    call update_residual(F_n(:, INT_LEFT),F_n(:, INT_RIGHT),(PDE_SRC_EM(1:4) + PDE_SRC_COLL(1:4)) + Electrons_RHS ,L_cells(C_ID))


 END DO

 call compute_residual()

 call implicit_step(delta_cons)


 sol(:,2:N_cells-1) = sol(:,2:N_cells-1) + delta_cons * relaxation_factor

 call check_non_physical_points(sol)

  !res_Rho_el = maxval(abs(delta_cons/maxval(d_tau_vect)))

  ! Then, impose Boundary conditions
  call IMPOSE_BC(SOL)

  END SUBROUTINE IMPLICIT_ELECTRONS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE MIDPOINT_EULER(sol)

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ,N_cells), INTENT(INOUT) :: sol
  REAL(KIND=8), DIMENSION(n_cells) :: E_SC
  INTEGER :: C_ID, I_INT,int_right,int_left
  REAL(KIND=8), DIMENSION(N_EQ) :: PDE_SRC,PDE_SRC_EM,PDE_SRC_COLL
  REAL(KIND=8), DIMENSION(N_EQ,N_cells) :: sol_half    ! Solution at half timestep
  REAL(KIND=8) :: dt_half
  REAL(KIND=8), DIMENSION(N_EQ) :: dt_vect
  REAL(KIND=8), DIMENSION(n_cells) :: Ey_field,b_field
  REAL(KIND=8), DIMENSION(N_EQ, N_int)  :: F_n, F_half ! Fluxes at timestep n and n+1/2
  real(kind=8), dimension(3) :: coll_vect
  REAL(KIND=8), DIMENSION(N_cells) :: sol_aux

  sol_aux = sol(1,:)


  ! Init
  sol_half(:,1)       = sol(:,1)       ! Set ghost cells!
  sol_half(:,N_cells) = sol(:,N_cells) ! Set ghost cells!

  E_SC = 0.0d0
  PDE_SRC = 0.0d0
  PDE_SRC_EM = 0.0d0
  PDE_SRC_COLL = 0.0d0

  ! Compute fluxes at interfaces and save them into F_int
  CALL COMPUTE_NUMERICAL_FLUXES(sol, F_n)

  ! Update self-consistent Electromagnetic field
  call compute_SC_EM_field(sol,E_SC(2:N_cells-1),Ey_field,B_field)


  ! Update solution (internal cells only, no ghost cells)
  dt_half = dt/2.0d0


  DO C_ID = 2, N_cells-1


    CALL COMPUTE_PDE_SOURCES(sol(:,C_ID), E_SC(C_ID), Ey_field(C_ID), B_field(C_ID), PDE_SRC_EM, PDE_SRC_COLL,C_ID) ! Compute sources in the cell

    PDE_SRC = PDE_SRC_EM + PDE_SRC_COLL

    INT_LEFT  = C_ID - 1
    INT_RIGHT = C_ID
    sol_half(:,C_ID) = sol(:, C_ID) + dt_half/(L_cells(C_ID))*(F_n(:, INT_LEFT) - F_n(:, INT_RIGHT)) &
                                    + dt_half*PDE_SRC

    call check_Temperature(sol_half(:,C_ID))

    call update_residual(F_n(:, INT_LEFT),F_n(:, INT_RIGHT),PDE_SRC,L_cells(C_ID))


  END DO

  call compute_residual()

  call check_non_physical_points(sol_half)


  ! Impose BCs for sol_half (before computing fluxes!)
  call IMPOSE_BC(SOL_HALF)

  ! And compute fluxes from updated solution
  CALL COMPUTE_NUMERICAL_FLUXES(sol_half, F_half)

  ! Update self-consistent Electromagnetic field
  call compute_SC_EM_field(sol_half,E_SC(2:N_cells-1),Ey_field,B_field)

  PDE_SRC = 0.0d0
  PDE_SRC_EM = 0.0d0
  PDE_SRC_COLL = 0.0d0


  ! ++++++++++++++++++++++ Compute final solution using stuff at half-step ++++++++++++++++

  DO C_ID = 2, N_cells-1

    CALL COMPUTE_PDE_SOURCES(sol_half(:,C_ID), E_SC(C_ID), Ey_field(C_ID), B_field(C_ID), PDE_SRC_EM, PDE_SRC_COLL,C_ID)  ! Compute sources in the cell from solution at half-step

    coll_elastic_vect_el(C_ID) = nu_implicit_el
    coll_elastic_vect_io(C_ID) = nu_implicit_io
    nu_ions_vect(1,C_ID) = nu_implicit_io
    nu_ions_vect(2,C_ID) = nu_implicit_el_iz
    nu_ions_vect(3,C_ID) = nu_implicit_el_ex

    PDE_SRC = PDE_SRC_EM + PDE_SRC_COLL

    INT_LEFT  = C_ID - 1
    INT_RIGHT = C_ID
    sol(:,C_ID) = sol(:,C_ID) + dt/(L_cells(C_ID))*(F_half(:, INT_LEFT) - F_half(:, INT_RIGHT)) &
                              + dt*PDE_SRC

    call check_Temperature(sol(:,C_ID))

! ONLY MIDPOINT SOURCES - TEST    sol(:,C_ID) = sol(:,C_ID) + dt/(L_cells(C_ID))*(F_n(:, INT_LEFT) - F_n(:, INT_RIGHT)) &
! ONLY MIDPOINT SOURCES - TEST                              + dt*PDE_SRC ! Midpoint Euler only for sources

  END DO

  call check_non_physical_points(sol)

  ! Enforcing boundary conditions
  call IMPOSE_BC(SOL)




  END SUBROUTINE MIDPOINT_EULER


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine check_non_physical_points(sol)

implicit NONE

REAL(KIND=8), DIMENSION(N_EQ,n_cells), intent(inout) :: sol
real(kind=8) :: P_ele,P_io
integer :: C_ID

DO C_ID = 2, N_cells-1

  P_ele  = (GAS_GAMMA_1 - 1)*(sol(4,C_ID) - sol(2,C_ID)**2/sol(1,C_ID)/2.0d0 - sol(3,C_ID)**2/sol(1,C_ID)/2.0d0)

  if (P_ele .lt. 0)then

    P_ele = 1d-3*P_init_electrons

  end if


  sol(4,C_ID) = sol(3,C_ID)**2/sol(1,C_ID)/2.0d0 + sol(2,C_ID)**2/2.0d0/sol(1,C_ID) + P_ele/(GAS_GAMMA_1 - 1)

  if(SYS_EULER_MF_BOOL)then

    P_io = (GAS_GAMMA_2 - 1)*(sol(8,C_ID) - sol(6,C_ID)**2/sol(5,C_ID)/2.0d0 - sol(7,C_ID)**2/sol(5,C_ID)/2.0d0)

    if (P_io .lt. 0)then

      P_io = 1d-3*P_init_ions

    end if

    sol(8,C_ID) = sol(7,C_ID)**2/sol(5,C_ID)/2.0d0 + sol(6,C_ID)**2/2.0d0/sol(5,C_ID) + P_io/(GAS_GAMMA_2 - 1)

  end if


end do

end subroutine check_non_physical_points

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine update_residual(F_L,F_R,SRC,dx)

implicit none

REAL(KIND=8), DIMENSION(N_EQ),intent(in) :: SRC
REAL(KIND=8), intent(in) ::dx
REAL(KIND=8), DIMENSION(N_EQ), intent(in) :: F_L, F_R

if ((SYS_EULER_BOOL .and. FLUID_ELECTRONS) .or. SYS_Euler_MF_BOOL)THEN

res_rho_el = res_rho_el + ((F_L(1) - F_R(1)) + SRC(1)*dx)**2
res_ene_el = res_ene_el + ((F_L(4) - F_R(4)) + SRC(4)*dx)**2

end if

if (SYS_EULER_BOOL .and. FLUID_IONS) THEN

  res_rho_i = res_rho_i + ((F_L(1) - F_R(1)) + SRC(1)*dx)**2
  res_ene_i = res_ene_i + ((F_L(4) - F_R(4)) + SRC(4)*dx)**2

end if

if (SYS_Euler_MF_BOOL)THEN

 res_rho_i = res_rho_i + ((F_L(5) - F_R(5)) + SRC(5)*dx)**2
 res_ene_i = res_ene_i + ((F_L(8) - F_R(8)) + SRC(8)*dx)**2

end if

end SUBROUTINE update_residual

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_residual()

implicit none

if ((SYS_EULER_BOOL .and. FLUID_ELECTRONS) .or. SYS_Euler_MF_BOOL)THEN

res_rho_el = sqrt(res_rho_el / (N_cells-2))
res_ene_el = sqrt(res_ene_el / (N_cells-2))


  if(isnan(res_rho_el))THEN
  print*, "                                                                 "
  print*, "                                                                 "
  print*, " ------------------------ WARNING -------------------------------"
  print*, " Simulation has diverged: NaN found in electrons density residual "
  print*, " ----------------------------------------------------------------"
  print*, "                                                                 "
  print*, "                                                               "
   !STOP
  end if

end if

if ((SYS_EULER_BOOL .and. FLUID_IONS) .or. SYS_Euler_MF_BOOL)THEN

  res_rho_i = sqrt(res_rho_i / (N_cells-2))
  res_ene_i = sqrt(res_ene_i / (N_cells-2))

  if(isnan(res_rho_i))THEN
    print*, "                                                                 "
    print*, "                                                                 "
    print*, " ------------------------ WARNING -------------------------------"
    print*, " Simulation has diverged: NaN found in ions density residual "
    print*, " ----------------------------------------------------------------"
    print*, "                                                                 "
    print*, "                                                                 "
   !STOP
  end if

end if


end subroutine compute_residual

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine IMPOSE_BC(SOL)

implicit NONE

REAL(KIND=8), DIMENSION(N_EQ,n_cells), intent(inout) :: sol

if(SYS_EULER_BOOL .and. FLUID_ELECTRONS)then

  call impose_electrons_BC(SOL)

else if(SYS_EULER_BOOL .and. FLUID_IONS )then

  call impose_ions_BC(SOL)

else if(SYS_EULER_MF_BOOL)THEN

  call impose_mf_BC(sol)

end if

end subroutine impose_BC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine check_Temperature(sol)

implicit none

REAL(KIND=8), DIMENSION(N_EQ), intent(inout) :: sol
real(kind=8) :: P2

if(SYS_Euler_MF_BOOL .and. IONS_Temperature_Control)THEN


  P2 = sol(5)*P_converter_ions

  sol(8) = sol(7)**2/sol(5)/2.0d0 + sol(6)**2/2.0d0/sol(5) + &
      P2/(GAS_GAMMA_2 - 1)

else if(SYS_Euler_MF_BOOL .and. ELECTRONS_Temperature_Control)THEN


  P2 = sol(1)*P_converter_electrons

  sol(4) = sol(3)**2/sol(1)/2.0d0 + sol(2)**2/2.0d0/sol(1) + &
      P2/(GAS_GAMMA_1 - 1)

else if (SYS_EULER_BOOL .and. FLUID_IONS .and. IONS_Temperature_Control)THEN

  P2 = sol(1)*P_converter_ions

  sol(4) = sol(3)**2/sol(1)/2.0d0 + sol(2)**2/2.0d0/sol(1) + &
      P2/(GAS_GAMMA_2 - 1)

end if

end subroutine check_Temperature

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE TIME_INTEGRATORS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
