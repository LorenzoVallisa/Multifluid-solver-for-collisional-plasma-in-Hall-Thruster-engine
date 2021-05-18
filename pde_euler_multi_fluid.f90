MODULE pde_euler_multi_fluid

  USE GLOBAL
  use rates_box

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!! Euler subroutines !!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE

  CONTAINS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE EULER_mf_RESTART_SOL(SOL)


  ! Notice that the mesh has to be EXACTLY THE SAME: same numebr of cells (otherwise it won't work)
  REAL(KIND=8), DIMENSION(N_EQ, N_cells), INTENT(INOUT) :: sol
  real(kind=8) :: dummy_time,dummy_x,dummy_E,dummy_rhoQ,dummy_phi
  integer(kind=8) :: IC,dummy_t


  OPEN(unit=123,file='restart_sol.dat', status='old', action='read')


      if (SC_EM_field) then

        DO IC = 2, N_cells-1 ! Write solution in physical cells (no ghost cells)

          ! Reading final time of previous simulation
          read(123, *) t_old_ID, dummy_x, sol(:, IC),dummy_E,dummy_rhoQ, dummy_phi

        END DO

      else

        DO IC = 2, N_cells-1 ! Write solution in physical cells (no ghost cells)
          read(123, *) t_old_ID, dummy_x, sol(:, IC)
        END DO

      end if

  close(123)

  write(*,*)" Starting from time   ",t_old_ID

  END SUBROUTINE EULER_mf_RESTART_SOL

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE EULER_MF_INIT_SOL(sol)
  ! Initializes the solution vector for the Euler PDEs.

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ, N_cells), INTENT(INOUT) :: sol
  REAL(KIND=8) :: rho1, ux1,uy1, P1, rho2, ux2,uy2, P2
  REAL(KIND=8) :: rhoL1, rhoR1, uxL1,uyL1, uxR1, uyR1, PL1, PR1
  REAL(KIND=8) :: rhoL2, rhoR2, uxL2,uyL2, uxR2, uyR2, PL2, PR2
  integer :: position
  real(kind=8) :: inlet_electrons_velocity
  position = 1

  ! Electrons
   rho1 = rho_init_electrons
   P1   = P_init_electrons
   ux1  = u_init_electrons
   uy1 = v_init_electrons

   ! Ions
   rho2 = rho_init_ions
   P2   = P_init_ions
   ux2  = u_init_ions
   uy2 = v_init_ions

   if(IONS_Temperature_Control)THEN

    P2 = rho2*P_converter_ions

   end if

   if(ELECTRONS_Temperature_Control)THEN

    P1 = rho1*P_converter_electrons

   end if


  if (cont_bool) then

   !!!! UNIFORM CONDITIONS !!!!

   ! Conserved variables equal on both sides
   sol(1,:) = rho1
   sol(2,:) = rho1*ux1
   sol(3,:) = rho1*uy1
   sol(4,:) = rho1*ux1**2/2.0d0 + rho1*uy1**2/2.0d0 + P1/(GAS_GAMMA_1 - 1)

   sol(5,:) = rho2
   sol(6,:) = rho2*ux2
   sol(7,:) = rho2*uy2
   sol(8,:) = rho2*ux2**2/2.0d0 + rho2*uy2**2/2.0d0 + P2/(GAS_GAMMA_2 - 1)

  else if(test_case_1) then

    !!! PLASMA EXPANSION !!!
    write(*,*) " Plasma expansion Test Case "


    do while (x_cc(position) .lt. 0)

      position = position + 1

    end do


    ! Conserved variables, left state
    sol(1,:) = rho1
    sol(2,:) = rho1*ux1
    sol(3,:) = rho1*uy1
    sol(4,:) = rho1*ux1**2/2.0d0 + rho1*uy1**2/2.0d0 + P1/(GAS_GAMMA_1 - 1)

    ! Electrons_vacuum case add-on
    sol(1,position:) = rho1*1d-3
    sol(2,position:) = 0.000d000
    sol(3,position:) = 0.0d0
    sol(4,position:) = P1*1d-3/(GAS_GAMMA_1 - 1)

    sol(5,:) = rho2
    sol(6,:) = rho2*ux2
    sol(7,:) = rho2*uy2
    sol(8,:) = rho2*ux2**2/2.0d0 + rho2*uy2**2/2.0d0 + P2/(GAS_GAMMA_2 - 1)

    ! Ions_vacuum case add-on
    sol(5,position:) = rho2*1d-3
    sol(6,position:) = 0.000d000
    sol(7,position:) = 0.0d0
    sol(8,position:) = P2*1d-3/(GAS_GAMMA_2 - 1)

  else

      !!! SOD SHOCK PROBLEM !!!

      rhoL1 = rho1
      PL1 =  P1
      uxL1  = ux1
      uyL1 =  0.0d0

      rhoR1 = rho1/4.0d0
      PR1   = P1/4.0d0
      uxR1  = ux1/4.0d0
      uyR1 = 0.0d0

      rhoL2 = rho2
      PL2   = P2
      uxL2  = ux2
      uyL2  = 0.0d0

      rhoR2 = rho2/4.0d0
      PR2   = P2/4.0d0
      uxR2  = ux2/4.0d0
      uyR2 = 0.0d0

      ! Conserved variables, left state
      sol(1,:) = rhoL1
      sol(2,:) = rhoL1*uxL1
      sol(3,:) = rhoL1*uyL1
      sol(4,:) = rho1*uyL1**2/2.0d0 + rhoL1*uxL1**2/2.0d0 + PL1/(GAS_GAMMA_1 - 1)

      sol(5,:) = rhoL2
      sol(6,:) = rhoL2*uxL2
      sol(7,:) = rhoL2*uyL2
      sol(8,:) = rhoL2*uyL2**2/2.0d0 + rhoL2*uxL2**2/2.0d0 + PL2/(GAS_GAMMA_2 - 1)

      !!! SOD SHOCK PROBLEM !!!   ! Conserved variables, right state
      sol(1,FLOOR(N_cells/2.0):) = rhoR1
      sol(2,FLOOR(N_cells/2.0):) = rhoR1*uxR1
      sol(3,FLOOR(N_cells/2.0):) = rhoR1*uyR1
      sol(4,FLOOR(N_cells/2.0):) = rhoR1*uyR1**2/2.0d0 + rhoR1*uxR1**2/2.0d0 + PR1/(GAS_GAMMA_1 - 1)

      sol(5,FLOOR(N_cells/2.0):) = rhoR2
      sol(6,FLOOR(N_cells/2.0):) = rhoR2*uxR2
      sol(7,FLOOR(N_cells/2.0):) = rhoR2*uyR2
      sol(8,FLOOR(N_cells/2.0):) = rhoR2*uyR2**2/2.0d0 + rhoR2*uxR2**2/2.0d0 + PR2/(GAS_GAMMA_2 - 1)

end if




  ! Note that in this way it works open BC and homogeneous Neumann only

  if (BOOL_RESTART) THEN

    call EULER_mf_RESTART_SOL(SOL)

  end if

  !------------------------ SETTING WEAK BC --------------------------------!

  call impose_mf_BC(sol)

  END SUBROUTINE EULER_MF_INIT_SOL

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine impose_mf_BC(sol)

  implicit NONE

  REAL(KIND=8), DIMENSION(N_EQ,N_Cells), intent(inout) :: sol
  real(kind=8) :: inlet_electrons_velocity


  IF (BCs_PERIODIC_BOOL_electrons) THEN

    sol(1:4, 1)       = sol(1:4, N_cells - 1) ! First ghost cell equal to last physical cell
    sol(1:4, N_cells) = sol(1:4, 2)           ! Last ghost cell equal to first physical cell

 end if

  if(BCs_OPEN_BOOL_Left_electrons)THEN

    sol(1:4,1) = BC_electrons_left

  end if


  if(BCs_OPEN_BOOL_Right_electrons)THEN

    sol(1:4,N_cells) = BC_electrons_right

  end if

  IF (BCs_NEU_BOOL_Left_electrons) then

    sol(:, 1)       = sol(:, 2)

  end if

  If (BCs_NEU_BOOL_Right_electrons) then

    sol(:, N_cells)       = sol(:,N_cells - 1)

  END IF

  If (BCs_ROBIN_BOOL_Right_electrons) then

    ! Density is equal across the Cathode
    sol(1, N_cells)       = sol(1,N_cells - 1)

    ! According to density value set an incoming velocity
    if(nondim_bool)THEN

    inlet_electrons_velocity = Electrons_discharge_current_cathode/&
                    (sol(1,N_cells)*rho_el_carac/GAS_M_1 * GAS_Q_1 * Engine_CS_Area)

      inlet_electrons_velocity = inlet_electrons_velocity / u_el_carac

    else

      inlet_electrons_velocity = 0.6/(sol(1,N_cells)/GAS_M_1 * GAS_Q_1 * 4d-3)

    end if

    sol(2,N_cells) = inlet_electrons_velocity * sol(1, N_cells)

    sol(3, N_cells)       = sol(3,N_cells - 1)
    sol(4, N_cells)       = sol(4,N_cells - 1)

  END IF

  IF (BCs_PERIODIC_BOOL_ions) THEN

    sol(5:8, 1)       = sol(5:8, N_cells - 1) ! First ghost cell equal to last physical cell
    sol(5:8, N_cells) = sol(5:8, 2)           ! Last ghost cell equal to first physical cell

  end if

  IF (BCs_OPEN_BOOL_Left_ions) then

    sol(5:8, 1)       = BC_ions_left

  end if

  if (BCs_OPEN_BOOL_Right_ions) then

    sol(5:8, N_cells) = BC_ions_right

  END IF

  IF (BCs_NEU_BOOL_Left_ions) then

    sol(5:8, 1)       = sol(5:8, 2)

  end if

  If (BCs_NEU_BOOL_Right_ions) then

    sol(5:8, N_cells)       = sol(5:8,N_cells - 1)

  END IF



  END SUBROUTINE impose_mf_BC

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine compute_Ustar_left_right_MF(sol,S1,S_star1,S2,S_star2,U_star)

  implicit none

  REAL(KIND=8), DIMENSION(N_EQ), intent(in) :: sol
  REAL(KIND=8), DIMENSION(N_EQ), intent(inout) :: u_star
  REAL(KIND=8), intent(in) :: S1,S_star1,S2,S_star2
  real(kind=8) :: rho1,ux1,uy1,E1,P1
  real(kind=8) :: rhoE1,rhoux1,rhouy1
  real(kind=8) :: rho2,ux2,uy2,E2,P2
  real(kind=8) :: rhoE2,rhoux2,rhouy2

  rho1    = sol(1)
  rhoux1  = sol(2)
  rhouy1  = sol(3)
  ux1     = rhoux1/(rho1 + 1.0e-35)
  uy1     = rhouy1/(rho1 + 1.0e-35)
  rhoE1   = sol(4)
  E1 = rhoE1/rho1
  P1      = (GAS_GAMMA_1 - 1)*(rhoE1 - rho1*ux1*ux1/2.0 - rho1*uy1*uy1/2.0)

  u_star(1) = 1

  u_star(2) = S_star1

  u_star(3) = uy1

  u_star(4) = E1 + (S_star1 - ux1) * (S_star1 + p1/(rho1 * (S1-ux1)))

  u_star(1:4) = u_star(1:4) * (rho1 * (S1 - ux1)/(S1 - S_star1))

  rho2    = sol(5)
  rhoux2  = sol(6)
  rhouy2  = sol(7)
  ux2     = rhoux2/(rho2 + 1.0e-35)
  uy2     = rhouy2/(rho2 + 1.0e-35)
  rhoE2   = sol(8)
  E2 = rhoE2/rho2
  P2      = (GAS_GAMMA_2 - 1)*(rhoE2 - rho2*ux2*ux2/2.0 - rho2*uy2*uy2/2.0)

  u_star(5) = 1

  u_star(6) = S_star2

  u_star(7) = uy2

  u_star(8) = E2 + (S_star2 - ux2) * (S_star2 + p2/(rho2 * (S2-ux2)))

  u_star(5:8) = u_star(5:8) * (rho2 * (S2 - ux2)/(S2 - S_star2))

  end subroutine compute_Ustar_left_right_MF


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE EULER_MF_fluxes(sol, Flux)

  ! Computes flux from a solution vector, of size N_EQ

  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(OUT) :: Flux

  REAL(KIND=8) :: rho1, rhoux1, rhouy1, ux1, uy1, rhoE1, P1, rho2, rhoux2, rhouy2, ux2, uy2, rhoE2, P2

  ! Unpack the solution and compute primitive variables
  !rho1    = sol(1)
  !rhoux1  = sol(2)
  !rhouy1  = sol(3)
  !ux1     = rhoux1/(rho1 + 1.0e-35)
  !uy1     = rhouy1/(rho1 + 1.0e-35)
  !rhoE1   = sol(4)
  !P1      = (GAS_GAMMA_1 - 1.0d0)*(rhoE1 - rhoux1*ux1/2.0d0 - rhouy1*uy1/2.0d0)

  !rho2    = sol(5)
  !rhoux2  = sol(6)
  !rhouy2  = sol(7)
  !ux2     = rhoux2/(rho2 + 1.0e-35)
  !uy2     = rhouy2/(rho2 + 1.0e-35)
  !rhoE2   = sol(8)
  !P2      = (GAS_GAMMA_2 - 1.0d0)*(rhoE2 - rhoux2*ux2/2.0d0 - rhouy2*uy2/2.0d0)



  ! Compute fluxes
  !Flux(1) = rhoux1
  !Flux(2) = rhoux1*ux1 + P1
  !Flux(3) = rhoux1*uy1
  !Flux(4) = rhoE1*ux1 + P1*ux1

  !Flux(5) = rhoux2
  !Flux(6) = rhoux2*ux2 + P2
  !Flux(7) = rhoux2*uy2
  !Flux(8) = rhoE2*ux2 + P2*ux2

  ! Unpack the solution and compute primitive variables (electrons)

    rho1    = sol(1)
    ux1     = sol(2)/(rho1 + 1.0e-35)
    uy1     = sol(3)/(rho1 + 1.0e-35)
    P1     = (GAS_GAMMA_1 - 1)*(sol(4) - rho1*ux1*ux1/2.0 - rho1*uy1*uy1/2.0)

    ! Compute fluxes
    !Flux(1) = rho1*ux1
    !Flux(2) = rho1*ux1*ux1 + P1
    !Flux(3) = rho1*ux1*uy1
    !Flux(4) = sol(4)*ux1 + P1*ux1

    ! Compute fluxes (electrons )

    Flux(1) = rho1*ux1
    Flux(2) = rho1*ux1*ux1 + P1
    Flux(3) = rho1*ux1*uy1
    Flux(4) = rho1*(ux1)**(3.0d0)/2.0d0 + P1*ux1/(GAS_GAMMA_1 - 1) + P1*ux1

    ! Unpack the solution and compute primitive variables (ions)

    rho2    = sol(5)
    ux2     = sol(6)/(rho2 + 1.0e-35)
    uy2     = sol(7)/(rho2 + 1.0e-35)
    P2     = (GAS_GAMMA_2 - 1)*(sol(8) - rho2*ux2*ux2/2.0 - rho2*uy2*uy2/2.0)

    ! Compute fluxes (ions )

    Flux(5) = rho2*ux2
    Flux(6) = rho2*ux2*ux2 + P2
    Flux(7) = rho2*ux2*uy2
    Flux(8) = rho2*(ux2)**(3.0d0)/2.0d0 + P2*ux2/(GAS_GAMMA_2 - 1) + P2*ux2



  END SUBROUTINE EULER_MF_fluxes

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE EULER_MF_sources(sol, Ex, Ey, B, Src_EM, Src_COLL)


  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8),                  INTENT(IN)  :: Ex, Ey, B
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(OUT) :: Src_EM
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(OUT) :: Src_COLL

  Src_COLL = 0.0d0
  Src_EM = 0.0d0


  if(nondim_bool)THEN

    CALL EULER_EM_SOURCES_MF_ad(sol, Ex, Ey, B,  Src_EM)

    CALL euler_collisional_sources_MF_ad(sol,  Src_COLL)

  else

    CALL EULER_EM_SOURCES_MF(sol, Ex, Ey, B,  Src_EM)

    CALL euler_collisional_sources_MF(sol, Src_COLL)

  end if


  END SUBROUTINE EULER_MF_sources

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE euler_collisional_sources_MF(sol, Src)


  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(inOUT) :: Src
  real(kind=8) :: nu_coll_elastic_electrons,nu_coll_ionizing,nu_coll_excitation
  REAL(KIND=8) :: rho1, ux1, uy1, rhoE1, P1, T1
  REAL(KIND=8) :: rho2, ux2, uy2, rhoE2, P2, T2

  real(kind=8) :: nu_coll_elastic_ions
  real(kind=8), parameter :: pi = 3.141593
  real(kind=8) :: ionized_k, excited_k , ions_elastic_k ,electrons_elastic_k
  real(kind=8) ::electron_number_density, T_neutrals, Pressure_PC_electrons

  nu_coll_elastic_electrons = 0.0d0
  nu_coll_excitation = 0.0d0
  nu_coll_ionizing = 0.d0
  src  = 0.0d0

  ! Reconstruct solution
  rho1  = sol(1)
  ux1   = sol(2)/rho1
  uy1   = sol(3)/rho1
  P1    = (GAS_GAMMA_1 - 1)*(sol(4) - rho1*ux1*ux1/2.0d0 - rho1*uy1*uy1/2.0d0)

  rho2  = sol(5)
  ux2    = sol(6)/rho2
  uy2    = sol(7)/rho2
  P2    = (GAS_GAMMA_2 - 1)*(sol(8) - rho2*ux2*ux2/2.0d0 - rho2*uy2*uy2/2.0d0)

  ! Compute temperatures (ideal gas)
  T1 = P1/(rho1*1.38064852d-23/GAS_M_1)
  T2 = P2/(rho2*1.38064852d-23/GAS_M_2)


  if (coll_bool_full) then


    !!!!!!!!!!!!!!!!!!!!!!!------------ ELECTRONS -----------------!!!!!!!!!!!!!!!!!!!!!!

    electron_number_density = rho1/GAS_M_1;

    !write(*,*)" Elastic electrons "
    call get_elastic_k_rate(electrons_elastic_k,T1)

    nu_coll_elastic_electrons = electrons_elastic_k*neutrals_density

    !write(*,*)" Ionization "
    call get_ionized_k_rate(ionized_k,T1)

    ! For the moment using same nu_coll_ionizing for both electrons and ions

    nu_coll_ionizing = ionized_k*neutrals_density

    ! For the moment no contribution to excitation mentioned within ions, so only electrons part

    !write(*,*)" Excitation "
    call get_excited_k_rate(excited_k,T1)

    nu_coll_excitation = excited_k*neutrals_density

    ! In case of elastic collision electron-neutrals total energy is assumed to be conserved
    Pressure_PC_electrons = 2.0d0/3.0d0*sol(4)



    src(1) = + rho1*nu_coll_ionizing

    ! Only momentum contribution for electrons

    src(2) =  - rho1*ux1*nu_coll_elastic_electrons
    src(3) =  - rho1*uy1*nu_coll_elastic_electrons

    src(4) =  - threshold_energy_xe*electron_number_density*nu_coll_ionizing  &                      ! Ionization energetic contribution                                                                           ! Elastic scattering contribution (void) - nu_coll_elastic_electrons*( sol(4) - 3.0d0/2.0d0*Pressure_PC_electrons)
                          - excitation_threshold*electron_number_density*nu_coll_excitation                ! Exctation energetic contribution


    !!!!!!!!!!!!!!!!!!!!!!!------------ IONS -----------------!!!!!!!!!!!!!!!!!!!!!!


    !write(*,*)" Elastic - Ions "
    call get_ions_elastic_k_rate(ions_elastic_k,T2)


    nu_coll_elastic_ions = ions_elastic_k*neutrals_density

    src(5) =  + GAS_M_2*electron_number_density*nu_coll_ionizing


    ! Only momentum contribution for ions

    src(6) =  - rho2*ux2*nu_coll_elastic_ions + GAS_M_2*electron_number_density*u_neutrals*nu_coll_ionizing
    src(7) =  - rho2*uy2*nu_coll_elastic_ions

    T_neutrals = pi/8.0d0*u_neutrals**2*GAS_M_2/k_boltzmann

    src(8) =  + nu_coll_ionizing*GAS_M_2*electron_number_density*(u_neutrals**2/2.0d0 + 3.0/2.0*k_boltzmann* &
              T_neutrals/GAS_M_2)


  else if(coll_bool_elastic) then

      !!!!!!!!!!!!!!!!!!!!!!!------------ ELECTRONS - ONLY ELASTIC -----------------!!!!!!!!!!!!!!!!!!!!!!

      electron_number_density = rho1/GAS_M_1;

      !write(*,*)" Elastic electrons "
      call get_elastic_k_rate(electrons_elastic_k,T1)

      nu_coll_elastic_electrons = electrons_elastic_k*neutrals_density



      ! Only momentum contribution for electrons

      src(2) =  - rho1*ux1*nu_coll_elastic_electrons
      src(3) =  - rho1*uy1*nu_coll_elastic_electrons



      !!!!!!!!!!!!!!!!!!!!!!!------------ IONS - ONLY ELASTIC -----------------!!!!!!!!!!!!!!!!!!!!!!


      !write(*,*)" Elastic - Ions "
      call get_ions_elastic_k_rate(ions_elastic_k,T2)


      nu_coll_elastic_ions = ions_elastic_k*neutrals_density



      ! Only momentum contribution for ions

      src(6) = src(6) - rho2*ux2*nu_coll_elastic_ions
      src(7) = src(7) - rho2*uy2*nu_coll_elastic_ions



  end if

  nu_implicit_el = nu_coll_elastic_electrons
  nu_implicit_io = nu_coll_elastic_ions
  nu_implicit_el_iz = nu_coll_ionizing
  nu_implicit_el_ex = nu_coll_excitation


  END SUBROUTINE euler_collisional_sources_MF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! NON DIMENSIONAL COLLISIONS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE euler_collisional_sources_Mf_ad(sol, Src)


  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(inOUT) :: Src
  real(kind=8) :: nu_coll_elastic_electrons,nu_coll_ionizing,nu_coll_excitation

  REAL(KIND=8) :: rho1, ux1, uy1, rhoE1, P1, T1
  REAL(KIND=8) :: rho2, ux2, uy2, rhoE2, P2, T2

  real(kind=8) :: nu_coll_elastic_ions
  real(kind=8), parameter :: pi = 3.141593
  real(kind=8) :: ionized_k, excited_k , ions_elastic_k ,electrons_elastic_k
  real(kind=8) ::electron_number_density, T_neutrals

  nu_coll_elastic_electrons = 0.0d0
  nu_coll_excitation = 0.0d0
  nu_coll_ionizing = 0.d0
  src  = 0.0d0


  ! Reconstruct of dimensional quantities
  rho1  = sol(1)*rho_el_carac
  ux1   = sol(2)/sol(1) * u_el_carac
  uy1   = sol(3)/sol(1) * uy_el_carac
  P1    = (GAS_GAMMA_1 - 1)*((sol(4)*Energy_el_carac*rho_el_carac) - rho1*ux1*ux1/2.0d0 - rho1*uy1*uy1/2.0d0)

  rho2  = sol(5)*rho_i_carac
  ux2    = sol(6)/sol(5)*u_i_carac
  uy2    = sol(7)/sol(5)*uy_i_carac
  P2    = (GAS_GAMMA_2 - 1)*((sol(8)*Energy_i_carac*rho_i_carac) - rho2*ux2*ux2/2.0d0 - rho2*uy2*uy2/2.0d0)

  ! Compute temperatures (ideal gas)
  T1 = P1/(rho1*1.38064852d-23/GAS_M_1)
  T2 = P2/(rho2*1.38064852d-23/GAS_M_2)


    if (coll_bool_full) then



    !!!!!!!!!!!!!!!!!!!!!!!------------ ELECTRONS -----------------!!!!!!!!!!!!!!!!!!!!!!

    !write(*,*)" Elastic electrons "
    call get_elastic_k_rate(electrons_elastic_k,T1)

    nu_coll_elastic_electrons = electrons_elastic_k*neutrals_density

    nu_coll_elastic_electrons = nu_coll_elastic_electrons/nu_coll_carac

    !write(*,*)" Ionization "
    call get_ionized_k_rate(ionized_k,T1)

    ! For the moment using same nu_coll_ionizing for both electrons and ions

    nu_coll_ionizing = ionized_k*neutrals_density

    nu_coll_ionizing = nu_coll_ionizing/nu_coll_carac

    ! For the moment no contribution to excitation mentioned within ions, so only electrons part

    !write(*,*)" Excitation "
    call get_excited_k_rate(excited_k,T1)

    nu_coll_excitation = excited_k*neutrals_density

    nu_coll_excitation = nu_coll_excitation/nu_coll_carac

    src(1) =  + sol(1)*nu_coll_ionizing*beta

    ! Only momentum contribution for electrons

    src(2) =  - sol(2)*nu_coll_elastic_electrons*beta
    src(3) =  - sol(3)*nu_coll_elastic_electrons*beta

    src(4) =    - threshold_energy_xe*sol(1)*nu_coll_ionizing*cp3  &                                         ! Ionization energetic contribution
              - excitation_threshold*sol(1)*nu_coll_excitation*cp3                           ! Exctation energetic contribution



    !!!!!!!!!!!!!!!!!!!!!!!------------ IONS -----------------!!!!!!!!!!!!!!!!!!!!!!


    !write(*,*)" Electrons Ions "
    call get_ions_elastic_k_rate(ions_elastic_k,T2)

    nu_coll_elastic_ions = ions_elastic_k*neutrals_density

    nu_coll_elastic_ions = nu_coll_elastic_ions/nu_coll_carac

    src(5) =  + sol(1)*nu_coll_ionizing*beta*tmc


    ! Only momentum contribution for ions

    src(6) =  - sol(6)*nu_coll_elastic_ions*beta + tmc*beta*sol(1)*u_neutrals_to_ions_ratio*nu_coll_ionizing
    src(7) =  - sol(7)*nu_coll_elastic_ions*beta

    src(8) =  + tmc*beta*sol(1)*Energy_neutrals_to_ions_ratio*nu_coll_ionizing


  else if (coll_bool_elastic) then


    !!!!!!!!!!!!!!!!!!!!!!!------------ ELECTRONS -----------------!!!!!!!!!!!!!!!!!!!!!!

    !write(*,*)" Elastic electrons "
    call get_elastic_k_rate(electrons_elastic_k,T1)

    nu_coll_elastic_electrons = electrons_elastic_k*neutrals_density

    nu_coll_elastic_electrons = nu_coll_elastic_electrons/nu_coll_carac


    ! Only momentum contribution for electrons

    src(2) =  - sol(2)*nu_coll_elastic_electrons*beta
    src(3) =  - sol(3)*nu_coll_elastic_electrons*beta




    !!!!!!!!!!!!!!!!!!!!!!!------------ IONS -----------------!!!!!!!!!!!!!!!!!!!!!!


    !write(*,*)" Electrons Ions "
    call get_ions_elastic_k_rate(ions_elastic_k,T2)

    nu_coll_elastic_ions = ions_elastic_k*neutrals_density

    nu_coll_elastic_ions = nu_coll_elastic_ions/nu_coll_carac

    ! Only momentum contribution for ions

    src(6) =  - sol(6)*nu_coll_elastic_ions*beta
    src(7) =  - sol(7)*nu_coll_elastic_ions*beta


  end if

  nu_implicit_el = nu_coll_elastic_electrons
  nu_implicit_io = nu_coll_elastic_ions
  nu_implicit_el_iz = nu_coll_ionizing
  nu_implicit_el_ex = nu_coll_excitation


  END SUBROUTINE euler_collisional_sources_MF_ad

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE EULER_EM_SOURCES_MF(sol, Ex, Ey, B, EM_SRC)
  ! Computes EM sources

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8),                  INTENT(IN)  :: Ex, Ey, B
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(inOUT) :: EM_SRC

  REAL(KIND=8) :: ux1, uy1, ux2, uy2

  EM_SRC = 0.0d0

  ux1  = sol(2)/(sol(1)+1.0e-35)
  uy1  = sol(3)/(sol(1)+1.0e-35)

  EM_SRC(1) = 0.0d0              ! Mass equation
  EM_SRC(2) = GAS_Q_1/GAS_M_1*(sol(1)*( Ex + uy1*B ) )  ! rho * q / m * ( Ex + uy B)
  EM_SRC(3) = GAS_Q_1/GAS_M_1*(sol(1)*( Ey - ux1*B ) )
  EM_SRC(4) = GAS_Q_1/GAS_M_1*(sol(2)*Ex + sol(3)*Ey)  ! rho * q / m * E * ux = (rho ux) * q / m * E


  ux2  = sol(6)/(sol(5)+1.0e-35)
  uy2  = sol(7)/(sol(5)+1.0e-35)

  EM_SRC(5) = 0.0d0              ! Mass equation
  EM_SRC(6) = GAS_Q_2/GAS_M_2*(sol(5)*( Ex + uy2*B ) )  ! rho * q / m * ( Ex + uy B)
  EM_SRC(7) = GAS_Q_2/GAS_M_2*(sol(5)*( Ey - ux2*B ) )
  EM_SRC(8) = GAS_Q_2/GAS_M_2*(sol(6)*Ex + sol(7)*Ey)  ! rho * q / m * E * ux = (rho ux) * q / m * E

  END SUBROUTINE EULER_EM_SOURCES_MF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! NON DIMENSIONAL EM

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE EULER_EM_SOURCES_MF_ad(sol, Ex, Ey, B, EM_SRC)
  ! Computes EM sources

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8),                  INTENT(IN)  :: Ex, Ey, B
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(inOUT) :: EM_SRC

  REAL(KIND=8) :: ux1, uy1, ux2, uy2

  EM_SRC = 0.0d0

  ux1  = sol(2)/(sol(1)+1.0e-35)
  uy1  = sol(3)/(sol(1)+1.0e-35)

  EM_SRC(1) = 0.0d0              ! Mass equation
  EM_SRC(2) = sol(1)*( E_x_e_AD * Ex + B_x_e_AD * uy1 * B )   ! rho * q / m * ( Ex + uy B)
  EM_SRC(3) = sol(1)*( jrswish * Ey - B_y_e_AD * ux1 * B )
  EM_SRC(4) = E_x_e_ene_AD * sol(2) * Ex + kcp_y * sol(3) * Ey  ! rho * q / m * E * ux = (rho ux) * q / m * E


  ux2  = sol(6)/(sol(5)+1.0e-35)
  uy2  = sol(7)/(sol(5)+1.0e-35)

  EM_SRC(5) = 0.0d0              ! Mass equation
  EM_SRC(6) = sol(5)*( E_x_i_AD * Ex + B_x_i_AD * uy2 * B )  ! rho * q / m * ( Ex + uy B)
  EM_SRC(7) = sol(5)*( rj * Ey - B_y_i_AD* ux2 *B )
  EM_SRC(8) = E_x_i_ene_AD * sol(6) * Ex + kg_y * sol(7) * Ey  ! rho * q / m * E * ux = (rho ux) * q / m * E

  END SUBROUTINE EULER_EM_SOURCES_MF_ad

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE EULER_2FLUIDS_MIN_MAX_WAVESPEEDS(sol, sMin1, sMax1, sMin2, sMax2)

  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8),                  INTENT(OUT) :: sMin1, sMax1, sMin2, sMax2

  REAL(KIND=8) :: rho1, rhoux1, rhouy1, ux1, uy1, rhoE1, P1
  REAL(KIND=8) :: rho2, rhoux2, rhouy2, ux2, uy2, rhoE2, P2


  ! Unpack the solution and compute primitive variables
  rho1    = sol(1)
  rhoux1  = sol(2)
  rhouy1  = sol(3)
  ux1     = rhoux1/(rho1 + 1.0e-35)
  uy1     = rhouy1/(rho1 + 1.0e-35)
  rhoE1   = sol(4)
  P1      = (GAS_GAMMA_1 - 1.0d0)*(rhoE1 - rhoux1*ux1/2.0d0 - rhouy1*uy1/2.0d0)

  rho2    = sol(5)
  rhoux2  = sol(6)
  rhouy2  = sol(7)
  ux2     = rhoux2/(rho2 + 1.0e-35)
  uy2     = rhouy2/(rho2 + 1.0e-35)
  rhoE2   = sol(8)
  P2      = (GAS_GAMMA_2 - 1.0d0)*(rhoE2 - rhoux2*ux2/2.0d0 - rhouy2*uy2/2.0d0)


  sMin1 = ux1 - SQRT(GAS_GAMMA_1*P1/(rho1 + 1.0e-35))
  sMax1 = ux1 + SQRT(GAS_GAMMA_1*P1/(rho1 + 1.0e-35))

  sMin2 = ux2 - SQRT(GAS_GAMMA_2*P2/(rho2 + 1.0e-35))
  sMax2 = ux2 + SQRT(GAS_GAMMA_2*P2/(rho2 + 1.0e-35))



  END SUBROUTINE EULER_2FLUIDS_MIN_MAX_WAVESPEEDS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE COMPUTE_MAX_MIN_IONS_WAVESPEEDS(sol, sMin2, sMax2)

  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8),                  INTENT(OUT) :: sMin2, sMax2


  REAL(KIND=8) :: rho2, rhoux2, rhouy2, ux2, uy2, rhoE2, P2

  rho2    = sol(5)
  rhoux2  = sol(6)
  rhouy2  = sol(7)
  ux2     = rhoux2/(rho2 + 1.0e-35)
  uy2     = rhouy2/(rho2 + 1.0e-35)
  rhoE2   = sol(8)
  P2      = (GAS_GAMMA_2 - 1.0d0)*(rhoE2 - rhoux2*ux2/2.0d0 - rhouy2*uy2/2.0d0)

  sMin2 = ux2 - SQRT(GAS_GAMMA_2*P2/(rho2 + 1.0e-35))
  sMax2 = ux2 + SQRT(GAS_GAMMA_2*P2/(rho2 + 1.0e-35))

END SUBROUTINE COMPUTE_MAX_MIN_IONS_WAVESPEEDS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE pde_euler_multi_fluid
