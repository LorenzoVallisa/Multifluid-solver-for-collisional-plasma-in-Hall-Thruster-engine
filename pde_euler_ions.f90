MODULE PDE_EULER_IONS

  USE GLOBAL
  use rates_box

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!! Euler subroutines !!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE IONS_INIT_SOL(sol)


    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(N_EQ, N_cells), INTENT(INOUT) :: sol
    REAL(KIND=8) :: rhoL, rhoR, uxL, uxR, uyL, uyR, PL, PR

    rhoL = rho_init_ions
    uxl = u_init_ions
    uyL  = v_init_ions
    PL = P_init_ions

    if(IONS_Temperature_Control)THEN

     PL = rhoL*P_converter_ions

    end if


  if (cont_bool) then

    sol(1,:) = rhoL
    sol(2,:) = rhoL*uxL
    sol(3,:) = rhoL*uyL
    sol(4,:) = rhoL*(uxL*uxL + uyL*uyL)/2 + PL/(GAS_GAMMA - 1)

  else if (disp_bool) then

    sol(1,:) = rhoL*(1 + 1e-4*sin(8*3.14159265359*x_cc/(x_max-x_min)))
    sol(2,:) = rhoL*uxL
    sol(3,:) = rhoL*uyL
    sol(4,:) = rhoL*(uxL*uxL + uyL*uyL)/2 + PL/(GAS_GAMMA - 1)

    !sol(1,:) = rhoL*(1 + 1e-4*sin(16*3.14159265359*x_cc/(x_max-x_min))) ! 1e-4*
    !sol(2,:) = sol(1,:) *uxL
    !sol(3,:) = sol(1,:) *uyL
    !sol(4,:) = sol(1,:) *(uxL*uxL + uyL*uyL)/2 + PL/(GAS_GAMMA - 1)

  else

     rhoR = rho_init_ions/4.0d0
     uxR  = u_init_ions/4.0d0
     uyR  = v_init_ions/4.0d0
     PR   = P_init_ions/4.0d0

    ! Conserved variables, left state
    sol(1,:) = rhoL
    sol(2,:) = rhoL*uxL
    sol(3,:) = rhoL*uyL
    sol(4,:) = rhoL*(uxL*uxL + uyL*uyL)/2 + PL/(GAS_GAMMA - 1)

    sol(1,FLOOR(N_cells/2.0):) = rhoR
    sol(2,FLOOR(N_cells/2.0):) = rhoR*uxR
    sol(3,FLOOR(N_cells/2.0):) = rhoR*uyR
    sol(4,FLOOR(N_cells/2.0):) = rhoR*(uxR*uxR + uyR*uyR)/2 + PR/(GAS_GAMMA - 1)

  end if


  ! Note that in this way it works for open BC and homogeneous Neumann only

  if (BOOL_RESTART) THEN

    call EULER_IONS_RESTART_SOL(SOL)

  end if

  !--------------------- SETTING WEAK BC ---------------------------!

  call impose_ions_BC(sol)

  END SUBROUTINE IONS_INIT_SOL


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine impose_ions_BC(sol)

  implicit none

  REAL(KIND=8), DIMENSION(N_EQ,N_cells), INTENT(INout)  :: sol
  real(kind=8) :: inlet_electrons_velocity


  IF (BCs_PERIODIC_BOOL_ions) THEN

    sol(:, 1)       = sol(:, N_cells - 1) ! First ghost cell equal to last physical cell
    sol(:, N_cells) = sol(:, 2)           ! Last ghost cell equal to first physical cell

  end if

  if(BCs_OPEN_BOOL_Left_ions)THEN

    sol(1:4,1) = BC_ions_left

  end if


  if(BCs_OPEN_BOOL_Right_ions)THEN

    sol(1:4,N_cells) = BC_ions_right

  end if

  IF (BCs_NEU_BOOL_Left_ions) then

    sol(:, 1)       = sol(:, 2)

  end if

  If (BCs_NEU_BOOL_Right_ions) then

    sol(:, N_cells)       = sol(:,N_cells - 1)

  END IF



  END SUBROUTINE impose_ions_BC


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE EULER_IONS_RESTART_SOL(SOL)


  ! Notice that the mesh has to be EXACTLY THE SAME: same numebr of cells (otherwise it won't work)
  REAL(KIND=8), DIMENSION(N_EQ, N_cells), INTENT(INOUT) :: sol
  real(kind=8) :: dummy_time,dummy_x,dummy_E,dummy_rhoQ,dummy_phi
  integer(kind=8) :: IC,dummy_t


  OPEN(unit=123,file='restart_sol.dat', status='old', action='read')


      ! Possibilità di riscriverlo come restart evitando le velocità lungo y
      if (SC_EM_field) then

        DO IC = 2, N_cells-1 ! Write solution in physical cells (no ghost cells)
          read(123, *) t_old_ID, dummy_x, sol(:, IC),dummy_E,dummy_phi,dummy_rhoQ
        END DO

      else

        DO IC = 2, N_cells-1 ! Write solution in physical cells (no ghost cells)
          read(123, *) t_old_ID, dummy_x, sol(:, IC)
        END DO

      end if

  close(123)

  write(*,*)" Starting from time   ",t_old_ID

  END SUBROUTINE EULER_IONS_RESTART_SOL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE IONS_EM_SOURCES(sol, Ex, Ey, B, Src_EM)

  ! Computes EM sources for ions

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8),                  INTENT(IN)  :: Ex, Ey, B
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(OUT) :: Src_EM
  real(kind=8) :: t1
  REAL(KIND=8) :: ux1, uy1

  Src_EM = 0.0d0

  ux1  = sol(2)/(sol(1)+1.0e-35)
  uy1  = sol(3)/(sol(1)+1.0e-35)


  if(nondim_bool)THEN !!!!!!!!!!!!!!!!!!!! NON DIM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Src_EM(1) = 0.0d0              ! Mass equation
    Src_EM(2) = sol(1)*( E_x_i_AD * Ex + B_x_i_AD * uy1 * B )  ! rho * q / m * ( Ex + uy B)
    Src_EM(3) = sol(1)*( rj * Ey - B_y_i_AD* ux1 *B )
    Src_EM(4) = E_x_i_ene_AD * sol(2) * Ex + kg_y * sol(3) * Ey  ! rho * q / m * E * ux = (rho ux) * q / m * E


  else


    Src_EM(1) = 0.0d0              ! Mass equation
    Src_EM(2) = GAS_Q/GAS_M*(sol(1)*( Ex + uy1*B ) )  ! rho * q / m * ( Ex + uy B)
    Src_EM(3) = GAS_Q/GAS_M*(sol(1)*( Ey - ux1*B ) )
    Src_EM(4) = GAS_Q/GAS_M*(sol(2)*Ex + sol(2)*Ey)  ! rho * q / m * E * ux = (rho ux) * q / m * E

  end if


  END SUBROUTINE IONS_EM_SOURCES

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE IONS_COLLISIONAL_SOURCES(sol, Src)


  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(inOUT) :: Src
  real(kind=8) :: nu_coll_elastic_ions,nu_coll_ionizing

  ! Ions local primitive variables
  REAL(KIND=8) :: Temperature,ux,uy,Pressure,rho

  real(kind=8), parameter :: pi = 3.141593
  real(kind=8) :: ionized_k, ions_elastic_k
  real(kind=8) :: electron_number_density, T_neutrals

  nu_coll_elastic_ions = 0.0d0
  nu_coll_ionizing = 0.d0
  src  = 0.0d0

  if (nondim_bool)THEN

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!! NON DIM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Reconstruct of dimensional quantities
    rho  = sol(1)*rho_el_carac

    Pressure    = (GAS_GAMMA - 1)*((sol(4)*Energy_el_carac*rho_el_carac) - sol(2)**2/sol(1)/2.0d0*&
          rho_el_carac*u_el_carac*u_el_carac   - sol(3)**2/sol(1)/2.0d0*&
                rho_el_carac*uy_el_carac*uy_el_carac)

    ! Compute temperatures (ideal gas)
    Temperature = Pressure/(rho*1.38064852d-23/GAS_M)

    if (coll_bool_full) then


      !write(*,*)" Ionization "
      call get_ionized_k_rate(ionized_k,Temperature_electrons)

      ! For the moment using same nu_coll_ionizing for both electrons and ions

      nu_coll_ionizing = ionized_k*neutrals_density

      nu_coll_ionizing = nu_coll_ionizing/nu_coll_carac


      !write(*,*)" Electrons Ions "
      call get_ions_elastic_k_rate(ions_elastic_k,Temperature)

      nu_coll_elastic_ions = ions_elastic_k*neutrals_density

      nu_coll_elastic_ions = nu_coll_elastic_ions/nu_coll_carac

      src(1) =  + Density_electrons*nu_coll_ionizing*beta*tmc


      ! Only momentum contribution for ions

      src(2) =  - sol(2)*nu_coll_elastic_ions*beta + tmc*beta*Density_electrons*u_neutrals_to_ions_ratio*nu_coll_ionizing
      src(3) =  - sol(3)*nu_coll_elastic_ions*beta

      src(4) =  + tmc*beta*Density_electrons*Energy_neutrals_to_ions_ratio*nu_coll_ionizing


    else if (coll_bool_elastic) then

      !write(*,*)" Electrons Ions "
      call get_ions_elastic_k_rate(ions_elastic_k,Temperature)

      nu_coll_elastic_ions = ions_elastic_k*neutrals_density

      nu_coll_elastic_ions = nu_coll_elastic_ions/nu_coll_carac

      ! Only momentum contribution for ions

      src(2) =  - sol(2)*nu_coll_elastic_ions*beta
      src(3) =  - sol(3)*nu_coll_elastic_ions*beta

    end if

  else

      rho = sol(1)
      ux = sol(2)/rho
      uy = sol(3)/rho
      Pressure = (GAS_GAMMA - 1)*(sol(4) - ux*ux*rho/2.0 - uy*uy*rho/2.0)
      Temperature = Pressure/(rho*1.38064852d-23/GAS_M)
      electron_number_density = Density_electrons/GAS_M;

    if (coll_bool_full) then


        !write(*,*)" Ionization "
        call get_ionized_k_rate(ionized_k,Temperature_electrons)

        nu_coll_ionizing = ionized_k*neutrals_density

        call get_ions_elastic_k_rate(ions_elastic_k,Temperature)

        nu_coll_elastic_ions = ions_elastic_k*neutrals_density

        src(1) =  + GAS_M*electron_number_density*nu_coll_ionizing


        ! Only momentum contribution for ions

        src(2) =  - rho*ux*nu_coll_elastic_ions + GAS_M*electron_number_density*u_neutrals*nu_coll_ionizing
        src(3) =  - rho*uy*nu_coll_elastic_ions

        T_neutrals = pi/8.0d0*u_neutrals**2*GAS_M/k_boltzmann

        src(4) =  + nu_coll_ionizing*GAS_M*electron_number_density*(u_neutrals**2/2.0d0 + 3.0/2.0*k_boltzmann* &
                  T_neutrals/GAS_M)



      else if (coll_bool_elastic) then


            call get_ions_elastic_k_rate(ions_elastic_k,Temperature)

            !WRITE(*,*)" --------------------- Elastic rate ------------------ "
            !WRITE(*,*)elastic_k

            nu_coll_elastic_ions = ions_elastic_k*neutrals_density

            ! In case of elastic collision electron-neutrals total energy is assumed to be conserved
            src(2) =  - sol(1)*ux*nu_coll_elastic_ions
            src(3) =  - sol(1)*uy*nu_coll_elastic_ions

      end if

  end if


  nu_implicit_el_iz = nu_coll_ionizing


  end SUBROUTINE IONS_COLLISIONAL_SOURCES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_electrons_primitive(sol)

! sol would be the solution of electrons conservatives only

implicit none

REAL(KIND=8), DIMENSION(4,N_cells), INTENT(IN)  :: sol
real(kind=8), DIMENSION(N_cells) :: rho, Pressure, Temperature, ux, uy

Density_electrons_vector = sol(1,:)

if(nondim_bool)THEN

  rho  = sol(1,:)*rho_el_carac
  Pressure   = (GAS_GAMMA - 1)*((sol(4,:)*Energy_el_carac*rho_el_carac) - sol(2,:)**2/sol(1,:)/2.0d0*&
        rho_el_carac*u_el_carac*u_el_carac   - sol(3,:)**2/sol(1,:)/2.0d0*&
              rho_el_carac*uy_el_carac*uy_el_carac)


  Temperature = Pressure/(rho*1.38064852d-23/GAS_M_1)

else

  ux = sol(2,:)/sol(1,:)
  uy = sol(3,:)/sol(1,:)
  Pressure = (GAS_GAMMA - 1)*(sol(4,:) - ux*ux*sol(1,:)/2.0 - uy*uy*sol(1,:)/2.0)

  Temperature = Pressure/(rho*1.38064852d-23/GAS_M_1)

end if

Temperature_electrons_vector = Temperature

end subroutine set_electrons_primitive

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



END MODULE PDE_EULER_IONS
