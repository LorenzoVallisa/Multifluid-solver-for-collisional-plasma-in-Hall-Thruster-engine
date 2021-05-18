PROGRAM FORMICA

  USE GLOBAL
  USE INITIALIZATION
  USE GRID
  USE EM_FIELDS
  USE VARIOUS
  USE NUMERICAL_FLUXES
  USE TIME_INTEGRATORS
  !use omp_lib
  use poisson_solver

  IMPLICIT NONE

  INTEGER :: C_ID, t_ID, iteration,PRINT_EVERY
  real(kind=8) :: t1,start,intermediate
  ! DBDBDB
  REAL(KIND=8), DIMENSION(10) :: sol
  REAL(KIND=8) :: rho, rhoux, ux, Pxx, Prr, qx, Riijj, SIGMA_64, SIGMA_65, SIGMA_66, SIGMA_67, SIGMA_68, dummy
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: U ! Solution, vector of conserved variables
  logical :: passed

  print *, " ----------------------------------------------------------------------------- "
  print *, "                              FORMICA                  "
  print *, " ------------------------------------------------------------------------------ "



  CALL READ_INPUT_FILE
  call init_phys_quantities
  call init_single_fluid_variables
  call init_initial_conditions
  call init_BC
  CALL INIT_VARIOUS
  CALL CREATE_GRID
  CALL INIT_EM_FIELDS
  call INIT_nondim_quantities
  CALL INIT_MATRIX
  call init_numerical_CFL_Settings

  ! Allocate solution and put to zero
  ALLOCATE(U(N_EQ, N_cells))         ! Size: N_EQ x N_cells

  U = 0

  CALL INIT_SOL(U) ! Initialize solution (done at the PDE level)

  ! ============= Solve in time =================
  ! First dt is the one set through input file
  PRINT_EVERY = WRITE_EVERY
  t_ID = 0
  global_iter = 0
  current_time = t_old_ID
  dual_iteration = 0
  passed = .False.
  if(DEBUG_BOOL)THEN
    PRINT_EVERY = 1
  END IF


  write(*,*)"--------------------------------------------------------------------"
  write(*,*)"-----------STARTING SIMULATION", current_time, "[s]"
  write(*,*)"--------------------------------------------------------------------"


  call cpu_time(start)
  DO while (current_time < t_end ) ! 1d-7 1.417d-6

    res_ene_el = 0.0d0
    res_rho_i = 0.0d0
    res_ene_i = 0.0d0
    res_rho_el = 0.0d0

    !---------------------------------- ITERATE --------------------------------------

    if (ELECTRONS_IMPLCIT_BOOL) then

      call COMPUTE_DUAL_TIMESTEP(U)

    else

      CALL COMPUTE_TIMESTEP(U)

    END IF

    t_ID = t_ID + 1

    current_time = current_time + dt

    global_iter = t_ID

    !---------------------------- POSTPROCESSING --------------------------------------------------


    IF(MOD(T_id,PRINT_EVERY) .EQ. 0)THEN

      PRINT*, "Max ions CFL (estimated):     "    , CFL_max_ions
      PRINT*, "Max electrons CFL (estimated): ", CFL_max_electrons
      PRINT*, "Current time ",current_time, ";  with dt = ",dt, " [s] , and maximum d_tau = ",maxval(d_tau_vect)

    END IF

    IF(MOD(t_ID, WRITE_EVERY) .EQ. 0) then

       CALL WRITE_SOLUTION_TO_FILE(t_ID,U)

    end if




   !------------------------------ CONVERGENCE CHECK ----------------------------------

    IF(MOD(t_ID, 5*WRITE_EVERY) .EQ. 0)then

       if(SYS_EULER_MF_BOOL)THEN

       PRINT*, "-------- GLOBAL CONVERGENCE UPDATE --------------"

           WRITE(*,*) "Density residuals :"
           WRITE(*,*) "Electrons :",log10(res_rho_el)
             if(isnan(res_rho_el))THEN
            print*, " Simulation has diverged: NaN found in electrons density residual "
            STOP
            end if
           WRITE(*,*) "Ions :",log10(res_rho_i)
           !WRITE(*,*) "Energy residuals :"
           !WRITE(*,*) "Electrons :",log10(res_ene_el)
           !WRITE(*,*) "Ions :",log10(res_ene_i)
           WRITE(*,*) "Current time :",current_time

      PRINT*, "--------------------------------------------"

      else

      PRINT*, "-------- GLOBAL CONVERGENCE UPDATE --------------"

           WRITE(*,*) "Density residuals :"
           WRITE(*,*) "Electrons :",log10(res_rho_el)
             if(isnan(res_rho_el))THEN
            print*, " Simulation has diverged: NaN found in electrons density residual "
            STOP
            end if
           WRITE(*,*) "Ions :",log10(res_rho_i)
           !WRITE(*,*) "Energy residuals :"
           !WRITE(*,*) "Electrons :",log10(res_ene_el)
           !WRITE(*,*) "Ions :",log10(res_ene_i)
           WRITE(*,*) "Current time :",current_time

      PRINT*, "--------------------------------------------"

      end if

    end if

    !if(t_ID .eq. 2)THEN
    !  call cpu_time(intermediate)
    !  print*," Simulation lasted ", intermediate - start, " [s] "
    !end if

    !call cpu_time(intermediate)
    !print*," We are at ", intermediate - start, " [s] "

    !if((mod(t_ID,1) .eq. 0) .and. (CFL_max_electrons .lt. 15) .and. (PASSED)) then
    !  CFL_INP_IONS = CFL_INP_IONS * (1 + 1d-1)
    !else if ((mod(t_ID,1) .eq. 0) .and. (CFL_max_electrons .gt. 20)) then
    !  CFL_INP_IONS = CFL_INP_IONS * (1 - 5d-1)
    !  PASSED = .TRue.
    !end if


  end do
  call cpu_time(intermediate)
  print*," Simulation lasted ", intermediate - start, " [s] "



  !print *," Time elapsed, for dual time stepping x100 , to reach 1e-7 [s]",omp_get_wtime()-t1
  !print*, " Iterations needed: ", t_ID
  !stop


  ! ========= Write solution =========
  !DO C_ID = 2, N_cells-1
  !  WRITE(*,*) "SOL: ", x_cc(C_ID), U(:, C_ID)
  !END DO

END PROGRAM FORMICA
