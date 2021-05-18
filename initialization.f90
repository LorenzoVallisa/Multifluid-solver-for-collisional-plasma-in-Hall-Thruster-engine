MODULE INITIALIZATION

  USE GLOBAL
  use implicit_electrons_utilities
  use rates_box
  use EM_fields


  CONTAINS


  SUBROUTINE READ_INPUT_FILE

    OPEN(unit=1,file='test_case.inp')

    READ(1,*)N_cells
    READ(1,*)x_min
    READ(1,*)x_max
    READ(1,*)GRID_TYPE
    READ(1,*)SYS_NAME
    READ(1,*)t_start
    READ(1,*)t_end
    READ(1,*)N_time
    READ(1,*)TIME_INTEGRATOR_NAME
    READ(1,*)RECONSTRUCTION_ORDER
    READ(1,*)WRITE_EVERY
    READ(1,*)type_e
    read(1,*)type_init
    read(1,*)type_coll
    READ(1,*)type_dim
    READ(1,*)type_restart
    read(1,*)flux_scheme
    read(1,*)type_debug


    CLOSE(1)

    write(*,*)"------------------------------"
    write(*,*)"       NUMERICAL SETTINGS"
    write(*,*)"------------------------------"
    PRINT *,  "Space"
    print*,   "Grid ",GRID_TYPE
    PRINT *, X_MIN," -- ",N_CELLS," --",X_MAX
    print *, " "
    PRINT *,  "Time"
    print *, " "
    PRINT *, T_START,"[s] --------> ",T_END, "[s]"


  END SUBROUTINE READ_INPUT_FILE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_numerical_CFL_Settings

  OPEN(unit=1990,file='cfl_input.inp')

  read(1990,*)CFL_INP_IONS
  read(1990,*)CFL_PHYS_MAX
  read(1990,*)CFL_INP_ELECTRONS
  read(1990,*)CFL_NUM_EL
  read(1990,*)toll
  read(1990,*)relaxation_factor
  read(1990,*)type_splitting
  read(1990,*)ions_subiterations
  read(1990,*)alfa_1
  read(1990,*)alfa_2


  CLOSE(1990)

  write(*,*)"------------------------------"
  write(*,*)"       CFL SETTINGS"
  write(*,*)"------------------------------"
  print*," Imposed CFL for ions population : ",CFL_INP_IONS
  print*," Imposed CFL for electrons population : ",CFL_INP_ELECTRONS
  print*," Imposed starting CFL for implicit electrons subiteration : ",CFL_NUM_EL
  print *," Relaxation factor set to ", relaxation_factor
  print *," alfa 1 ", alfa_1
  print *," alfa 2 ", alfa_2

  end subroutine init_numerical_CFL_Settings
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE init_single_fluid_variables

    OPEN(unit=69,file='single_fluid_vars.inp')

    read(69,*)n_ions_background
    read(69,*)Density_electrons
    read(69,*)Temperature_electrons
    read(69,*)type_control_ions
    read(69,*)T_fixed_ions
    read(69,*)type_control_electrons
    read(69,*)T_fixed_electrons

    CLOSE(69)

  END SUBROUTINE init_single_fluid_variables

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE init_phys_quantities

    OPEN(unit=2,file='physical_constants.inp')

    READ(2,*)GAS_GAMMA_1
    READ(2,*)Gas_M_1
    READ(2,*)Gas_q_1
    READ(2,*)k_boltzmann
    READ(2,*)eps_0
    READ(2,*)GAS_GAMMA_2
    READ(2,*)Gas_M_2
    READ(2,*)Gas_q_2
    READ(2,*)neutrals_density
    READ(2,*)u_neutrals
    READ(2,*)threshold_energy_xe
    read(2,*)excitation_threshold

    CLOSE(2)


  END SUBROUTINE init_phys_quantities

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE init_initial_conditions

    OPEN(unit=3,file='initial_conditions.inp')

    READ(3,*)n_init_electrons
    READ(3,*)T_init_electrons
    READ(3,*)u_init_electrons
    READ(3,*)v_init_electrons
    READ(3,*)n_init_ions
    READ(3,*)T_init_ions
    READ(3,*)u_init_ions
    READ(3,*)v_init_ions


    CLOSE(3)

    rho_init_electrons = n_init_electrons * GAS_M_1
    rho_init_ions = n_init_ions * GAS_M_2

    P_init_electrons = rho_init_electrons/GAS_M_1*k_boltzmann*T_init_electrons
    P_init_ions = rho_init_ions/GAS_M_2*k_boltzmann*T_init_ions

  END SUBROUTINE init_initial_conditions

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

  SUBROUTINE INIT_BC

  implicit NONE


    OPEN(unit=95,file='boundary_conditions.inp')

    read(95,*)BCs_TYPE_Right_electrons
    read(95,*)n_el_r
    READ(95,*)T_el_r
    READ(95,*)u_el_r
    READ(95,*)uy_el_r
    read(95,*)BCs_TYPE_Left_electrons
    read(95,*)n_el_l
    READ(95,*)T_el_l
    READ(95,*)u_el_l
    READ(95,*)uy_el_l
    read(95,*)BCs_TYPE_Right_ions
    read(95,*)n_i_r
    READ(95,*)T_i_r
    READ(95,*)u_i_r
    READ(95,*)uy_i_r
    read(95,*)BCs_TYPE_Left_ions
    read(95,*)n_i_l
    READ(95,*)T_i_l
    READ(95,*)u_i_l
    READ(95,*)uy_i_l

    close(95)

    rho_i_r = n_i_r * GAS_M_2
    rho_i_l = n_i_l * GAS_M_2
    rho_el_r = n_el_r * GAS_M_1
    rho_el_l = n_el_l * GAS_M_1

    allocate(BC_electrons_right(4),BC_electrons_left(4),BC_ions_right(4),BC_ions_left(4))

    BC_electrons_right(1)=rho_el_r
    BC_electrons_right(2)=rho_el_r*u_el_r
    BC_electrons_right(3)=rho_el_r*uy_el_r
    BC_electrons_right(4)=rho_el_r*u_el_r**2/2.0d0 + rho_el_r*uy_el_r**2/2.0d0 + &
      (rho_el_r/GAS_M_1*k_boltzmann*T_el_r)/(Gas_gamma_1 - 1)

    BC_electrons_left(1)=rho_el_l
    BC_electrons_left(2)=rho_el_l*u_el_l
    BC_electrons_left(3)=rho_el_l*uy_el_l
    BC_electrons_left(4)=rho_el_l*u_el_l**2/2.0d0 + rho_el_l*uy_el_l**2/2.0d0 + &
      (rho_el_l/GAS_M_1*k_boltzmann*T_el_l)/(Gas_gamma_1 - 1)

    BC_ions_right(1)=rho_i_r
    BC_ions_right(2)=rho_i_r*u_i_r
    BC_ions_right(3)=rho_i_r*uy_i_r
    BC_ions_right(4)=rho_i_r*u_i_r**2/2.0d0 + rho_i_r*uy_i_r**2/2.0d0 + &
     (rho_i_r/GAS_M_2*k_boltzmann*T_i_r)/(Gas_gamma_2 - 1)

    BC_ions_left(1)=rho_i_l
    BC_ions_left(2)=rho_i_l*u_i_l
    BC_ions_left(3)=rho_i_l*uy_i_l
    BC_ions_left(4)=rho_i_l*u_i_l**2/2.0d0 + rho_i_l*uy_i_l**2/2.0d0 + &
     (rho_i_l/GAS_M_2*k_boltzmann*T_i_l)/(Gas_gamma_2 - 1)

  end SUBROUTINE INIT_BC

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

  SUBROUTINE INIT_nondim_quantities

  implicit NONE

  real(kind=8), parameter :: pi = 3.141593
  real(kind=8) :: T_el, T_i,P_ad_el_l,P_ad_el_r,P_ad_i_r,P_ad_i_l

  if(Nondim_bool) then

    OPEN(unit=90,file='characteristics_quantities.inp')

    read(90,*)rho_el_carac
    read(90,*)rho_i_carac
    read(90,*)x_carac
    read(90,*)nu_coll_carac
    read(90,*)phi_carac
    read(90,*)B_mag_carac
    read(90,*)Ey_carac
    read(90,*)uy_el_carac
    read(90,*)uy_i_carac
    read(90,*)u_el_carac


    close(90)


    !!!!!!!!!!!!! POISSON SOLVER !!!!!!!!!!!!!!!!!!!!!

    ! Hypothesis for electric field to be checked: indeed when characteristic for phi is given
    E_carac = phi_carac/x_carac
    kobe = rho_el_carac/rho_i_carac*GAS_M_2/GAS_M_1
    shaq = phi_carac*GAS_M_2/rho_i_carac*eps_0/GAS_Q_1/x_carac/x_carac  ! Not fully non-dimensional

    ! Hypothesis for phi to be checked: indeed when characteristic for E_el is given
    ! phi_carac = E_carac*x_carac

    !!!!!!!!!!! ELECTRONS !!!!!!!!!!!!!
    P_el_carac = rho_el_carac * u_el_carac * u_el_carac
    Energy_el_carac = P_el_carac/rho_el_carac
    t_carac = x_carac/u_el_carac


    ! First 2 non-dimensional quantites within MASS-ELECTRONS
    beta = nu_coll_carac*t_carac

    ! Togheter with alfa and beta those are the non-dim quantities for X-MOMENTUM-ELECTRONS
    E_x_e_AD = GAS_Q_1 * E_carac * x_carac / GAS_M_1 / u_el_carac /u_el_carac
    B_x_e_AD = GAS_Q_1/GAS_M_1*b_mag_carac*uy_el_carac/u_el_carac*t_carac
    kb8 = P_el_carac/rho_el_carac/u_el_carac/u_el_carac

    ! Non-dim quantities for Y-MOMENTUM-ELECTRONS
    jrswish = GAS_Q_1/GAS_M_1*Ey_carac/uy_el_carac*t_carac
    B_y_e_AD = GAS_Q_1/GAS_M_1*b_mag_carac*t_carac*u_el_carac/uy_el_carac

    ! Togheter with alfa and beta those are the non-dim quantities for ENERGY-ELECTRONS
    cp3 = nu_coll_carac*t_carac/GAS_M_1/Energy_el_carac         ! Not fully non-dimensional
    E_x_e_ene_AD = GAS_Q_1 * E_carac * x_carac / GAS_M_1 / Energy_el_carac
    kcp_y =  GAS_Q_1/GAS_M_1*Ey_carac*t_carac*uy_el_carac/Energy_el_carac
    kp = P_el_carac/rho_el_carac/Energy_el_carac


    !!!!!!!!!!! IONS !!!!!!!!!!!!!

    u_i_carac = u_el_carac
    Energy_i_carac = u_i_carac*u_i_carac
    P_i_carac = rho_i_carac * u_i_carac * u_i_carac


    u_neutrals_to_ions_ratio = u_neutrals/u_i_carac
    T_neutrals = pi/8.0d0*u_neutrals**2*GAS_M_2/k_boltzmann
    Energy_neutrals = 0.5*u_neutrals*u_neutrals + 3.0d0/2.0d0*k_boltzmann*T_neutrals/GAS_M_2    ! No y component of neutrals velocity
    Energy_neutrals_to_ions_ratio = Energy_neutrals/Energy_i_carac

    ! Parameter introduced to account for fact that during ionization the source term is
    ! mass_neutrals*electron_number_density*un_coll_iz ; mass_neutrals = mass_ion
    tmc = rho_el_carac/rho_i_carac*GAS_M_2/GAS_M_1

    ! Togheter with alfa and beta those are the non-dim quantities for X-MOMENTUM-IONS
    E_x_i_AD = GAS_Q_2 * E_carac * x_carac / GAS_M_2 / u_i_carac /u_i_carac
    B_x_i_AD = GAS_Q_2/GAS_M_2*b_mag_carac*t_carac*uy_i_carac/u_i_carac
    kb24 = P_i_carac/rho_i_carac/u_i_carac/u_i_carac


    ! Togheter with alfa and beta those are the non-dim quantities for Y-MOMENTUM-IONS
    rj = GAS_Q_2/GAS_M_2*Ey_carac*t_carac/uy_i_carac
    B_y_i_AD = GAS_Q_2/GAS_M_2*b_mag_carac*t_carac*u_i_carac/uy_i_carac

    ! Non-dim quantities for ENERGY-IONS
    E_x_i_ene_AD = GAS_Q_2 * E_carac * x_carac / GAS_M_2 / Energy_i_carac
    kg_y = GAS_Q_2/GAS_M_2*Ey_carac*t_carac*uy_i_carac*t_carac/Energy_i_carac
    thj = P_i_carac/rho_i_carac/Energy_i_carac

    !!!!! Some DEBUG !!!

    write(*,*)" --------------Reference variables derived-------------"
    write(*,*)" Velocity-x ",u_el_carac,"  -   ",u_i_carac
    write(*,*)" Velocity-y ",uy_el_carac,"  -   ",uy_i_carac
    write(*,*)" Pressure ",P_el_carac,"  -   ",P_i_carac
    write(*,*)" Energy ",Energy_el_carac,"  -   ",Energy_i_carac
    write(*,*)" Electric_x_field", E_carac
    write(*,*)" Potential ", phi_carac
    write(*,*)" Time ",t_carac
    write(*,*)" Space ",x_carac

    ! Nondimensionalization of intitial conditions

    rho_init_electrons = rho_init_electrons/rho_el_carac
    P_init_electrons = P_init_electrons/P_el_carac
    u_init_electrons = u_init_electrons/u_el_carac
    v_init_electrons = v_init_electrons/uy_el_carac

    rho_init_ions = rho_init_ions/rho_i_carac
    P_init_ions = P_init_ions/P_i_carac
    u_init_ions = u_init_ions/u_i_carac
    v_init_ions = v_init_ions/uy_i_carac

    ! Nondimensionalization of boundary conditions
    P_ad_el_r = (rho_el_r/GAS_M_1*k_boltzmann*T_el_r)/P_el_carac
    rho_el_r = rho_el_r/rho_el_carac
    u_el_r = u_el_r/u_el_carac
    uy_el_r = uy_el_r/uy_el_carac

    BC_electrons_right(1)=rho_el_r
    BC_electrons_right(2)=rho_el_r*u_el_r
    BC_electrons_right(3)=rho_el_r*uy_el_r
    BC_electrons_right(4)=rho_el_r*u_el_r**2/2.0d0 + rho_el_r*uy_el_r**2/2.0d0 + P_ad_el_r/(Gas_gamma_1 - 1)

    P_ad_el_l = (rho_el_l/GAS_M_1*k_boltzmann*T_el_l)/P_el_carac
    rho_el_l = rho_el_l/rho_el_carac
    u_el_l = u_el_l/u_el_carac
    uy_el_l = uy_el_l/uy_el_carac

    BC_electrons_left(1)=rho_el_l
    BC_electrons_left(2)=rho_el_l*u_el_l
    BC_electrons_left(3)=rho_el_l*uy_el_l
    BC_electrons_left(4)=rho_el_l*u_el_l**2/2.0d0+ rho_el_l*uy_el_l**2/2.0d0 + P_ad_el_l/(Gas_gamma_1 - 1)

    P_ad_i_l = (rho_i_l/GAS_M_2*k_boltzmann*T_i_l)/P_i_carac
    rho_i_l = rho_i_l/rho_i_carac
    u_i_l = u_i_l/u_i_carac
    uy_i_l = uy_i_l/uy_i_carac

    BC_ions_left(1)=rho_i_l
    BC_ions_left(2)=rho_i_l*u_i_l
    BC_ions_left(3)=rho_i_l*uy_i_l
    BC_ions_left(4)=rho_i_l*u_i_l**2/2.0d0 + rho_i_l*uy_i_l**2/2.0d0 + P_ad_i_l/(Gas_gamma_2 - 1)

    P_ad_i_r = (rho_i_r/GAS_M_2*k_boltzmann*T_i_r)/P_i_carac
    rho_i_r = rho_i_r/rho_i_carac
    u_i_r = u_i_r/u_i_carac
    uy_i_r = uy_i_r/uy_i_carac

    BC_ions_right(1)=rho_i_r
    BC_ions_right(2)=rho_i_r*u_i_r
    BC_ions_right(3)=rho_i_r*uy_i_r
    BC_ions_right(4)=rho_i_r*u_i_r**2/2.0d0  +rho_i_r*uy_i_r**2/2.0d0 + P_ad_i_r/(Gas_gamma_2 - 1)

    ! Nondimensionalization of numerical settings

    dt = dt / t_carac
    d_tau = d_tau / t_carac
    t_end = t_end / t_carac
    t_old_id = t_old_id / t_carac
    x_min = x_min  /x_carac
    x_max = x_max /x_carac
    x_exhaust = x_exhaust/ x_carac
    l_cells = l_cells / x_carac

    ! Nondimensionalization of Em variables

    phi_0 = phi_0/phi_carac
    phi_end = phi_end/phi_carac
    Neu_left = Neu_left / E_carac
    Neu_right = Neu_right /E_carac
    B_field = B_field/B_mag_carac

    ! Nondimensional parameters for Poisson predictor
    Predictor_1 = t_carac /eps_0 * rho_el_carac * (Gas_q_1/GAS_M_1)**2 /nu_coll_carac
    Predictor_2 = t_carac /eps_0 * rho_i_carac * (Gas_q_2/GAS_M_2)**2 /nu_coll_carac
    Predictor_3 = GAS_Q_1 * rho_el_carac * x_carac * x_carac /GAS_M_1 / phi_carac
    Predictor_3 = GAS_Q_2 * rho_i_carac * x_carac * x_carac /GAS_M_2 / phi_carac

    P_converter_ions = rho_i_carac/GAS_M_2*T_fixed_ions*k_boltzmann/P_i_carac
    P_converter_electrons = rho_el_carac/GAS_M_1*T_fixed_electrons*k_boltzmann/P_el_carac


    write(*,*)"-------------------INITIAL NON-DIMENSIONALIZED QUANTITIES----------"
    write(*,*)" Density:  ",rho_init_electrons,"   -   ",rho_init_ions
    write(*,*)" Velocity - x:  ",u_init_electrons,"   -   ",u_init_ions
    write(*,*)" Velocity - y:  ",v_init_electrons,"   -   ",v_init_ions
    write(*,*)" Pressure:  ",P_init_electrons,"   -   ",P_init_ions
    write(*,*)" dt: ",dt
    write(*,*)" Potential phi BC:  ",phi_0,"  -   ",phi_end
    write(*,*)" Neumann: ",Neu_left,"   -   ",Neu_right
    write(*,*)" B (at exhaust): ",maxval(B_field)


    write(*,*)" --------------Non-dimensional quantities: team electrons -------------"
    write(*,*)" ------------ Convective parameters ----------------"
    write(*,*)kb8,"     -> Pe/(rhoe * ue^2)"
    write(*,*)kp,"      -> Pe/(rhoe * Ee)"
    write(*,*)" ------------ EM parameters ----------------"
    write(*,*)B_y_e_AD,"         ->qe / me * B * t * ue / ve"
    write(*,*)E_x_e_AD,"       -> qe * Ex * x / (me * ue * ue)"
    write(*,*)E_x_e_ene_AD,"         -> qe * Ex * x /( me * Ee)"
    write(*,*)B_x_e_AD,"         ->qe / me * B * ve / ue * t"


    write(*,*)" --------------Non-dimensional quantities: team ions -------------"
    write(*,*)" ------------ Convective parameters ----------------"
    write(*,*)kb24,"     -> Pi/(rhoi * ui^2)"
    write(*,*)thj,"      -> Pi/(rhoi * Ei)"
    write(*,*)" ------------ EM parameters ----------------"
    write(*,*)E_x_i_AD,"       -> qi * Ex * x / (mi * ui * ui)"
    write(*,*)E_x_i_ene_AD,"         -> qi * Ex *x /( mi * Ei)"

    write(*,*)" ------------ Collisional parameters (In common)----------------"
    write(*,*)beta,"         ->nu_coll * t "
    write(*,*)cp3,"         ->nu_coll * t /(me * Ee)"
    write(*,*)tmc,"         ->rho_el /rho_i * mi / me"

    !!!!! End of Some DEBUG !!!

  end if

  end subroutine INIT_nondim_quantities


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE INIT_VARIOUS

  implicit NONE

#include "lisf.h"

  integer :: i,triple_index_iii
  LIS_INTEGER :: ierr,comm
  CHARACTER(LEN=512) :: options

  if(type_splitting .eq. "TS")THEN
    TS_BOOL = .True.
  end if


    ! Default problem with dimensional quantities
    if (type_dim .EQ. "Nondim") then
        Nondim_bool = .True.
    end if

    if (flux_scheme .EQ. "HLL") then
        HLL_BOOL = .True.
    else if(flux_scheme .EQ. "HLLC")THEN
        HLLC_BOOL = .True.
    end if

    write(*,*)"------------------------------"
    write(*,*)"       SELF CONSISTENT EM "
    write(*,*)"------------------------------"
    ! SC EM default not included
    if (type_e .EQ. "esc") then
        SC_EM_field = .True.
        PRINT *, " Simulation with Poisson solver "
    end if



    ! Set Restart bool
    if (type_restart .EQ. "Restart") then
        BOOL_RESTART = .True.
    end if

    write(*,*)"------------------------------"
    write(*,*)"       COLLISIONAL PART "
    write(*,*)"------------------------------"
    ! Default no collision
    if (type_coll .EQ. "collision_full") then
        coll_bool_full = .True.
        PRINT *," Elastic/Inelastic: Boltzmann operator derivation"
    else if (type_coll .EQ. "collisionA") then
        coll_boolA = .True.
        PRINT *," Elastic: thermal velocity approximation "
    else if (type_coll .EQ. "collision_elastic") then
        coll_bool_elastic = .True.
        PRINT *," Elastic: Boltzmann operator derivation"
    end if

    if(type_control_electrons .eq. "ELECTRONS_CHECK")THEN
      ELECTRONS_Temperature_Control = .True.
    end if

    write(*,*)"------------------------------"
    write(*,*)"       TEST CASE "
    write(*,*)"------------------------------"
    ! Default Riemann problem
    if (type_init .EQ. "continous") then
        cont_bool = .True.
        PRINT *," Uniform initial conditions: Complete Hall-Thruster test case"
    else if (type_init .EQ. "dispersion") THEN
        disp_bool = .True.
        PRINT *," Testing dispersion relations "
    else if (type_init .EQ. "plasma_expansion") Then
        test_case_1 = .True.
        if ( ELECTRONS_Temperature_Control ) THEN
          PRINT *," Plasma expansion into vacuum, with electrons Temperature fixed at ",T_fixed_electrons
        else
          PRINT *," Plasma expansion into vacuum"
        end if
    else
        PRINT *," Sod-Shock problem"
    end if



    write(*,*)"------------------------------"
    write(*,*)"       MODEL TYPE"
    write(*,*)"------------------------------"
    ! It's very procedural.... but set here the properties of the system
    IF (SYS_NAME .EQ. "Euler_Electrons") THEN
      SYS_EULER_BOOL = .True.
      FLUID_ELECTRONS = .True.
      N_EQ = 4
      GAS_M=GAS_M_1
      GAS_Q=GAS_Q_1
      GAS_GAMMA=GAS_GAMMA_1
      PRINT *," Electrons as single fluid: Euler set of equations "
    else iF (SYS_NAME .EQ. "Euler_Ions") THEN
      SYS_EULER_BOOL = .True.
      FLUID_IONS = .True.
      N_EQ = 4
      GAS_M=GAS_M_2
      GAS_Q=GAS_Q_2
      GAS_GAMMA=GAS_GAMMA_2
      PRINT *," Ions as single fluid: Euler set of equations "
    ELSE IF (SYS_NAME .EQ. "Euler_multi_fluid") THEN
      SYS_Euler_MF_BOOL = .True.
      N_EQ = 8
      PRINT *," Multi fluid Euler "
    ELSE IF (SYS_NAME .EQ. "5mom") THEN
      SYS_MOM5_BOOL = .True.
      N_EQ = 5
    ELSE IF (SYS_NAME .EQ. "14mom") THEN
      SYS_MOM14_BOOL = .True.
      N_EQ = 14
    ELSE IF (SYS_NAME .EQ. "14mom_axi") THEN
      SYS_MOM14_AXI_BOOL = .True.
      N_EQ = 6
    ELSE
      WRITE(*,*) "Attention! System type not recognized! Check input file. ABORTING!"
      STOP
    END IF

    write(*,*)"------------------------------"
    write(*,*)"       NONDIMENSIONAL OPTION"
    write(*,*)"------------------------------"
    ! Check coeherence of input
    !if (SYS_EULER_BOOL .and. Nondim_bool)THEN

    !  write(*,*)" ERROR: Nondimensional analysis for single fluid solver not yet implemented, please try later "
    !  STOP

    !end if

    IF(NONDIM_BOOL)Then
      PRINT *," Nondimensional simulation "
    end if

    write(*,*)"------------------------------"
    write(*,*)"       RESTART OPTION"
    write(*,*)"------------------------------"
    if(BOOL_RESTART)THEN
      if(SYS_Euler_MF_BOOL)then
        write(*,*)" Restarting multi-fluid solution "
      else
        write(*,*)" Restarting single-fluid solution "
      end if
    end if

    write(*,*)"------------------------------"
    write(*,*)"       CFL OPTION "
    write(*,*)"------------------------------"
    ! CFL defualt changes
    if (type_dt .EQ. "CFL") then
        CFL_BOOL = .True.
        PRINT *, " Simulation with dynamic dt, CFL to be respected is ",CFL_INPUT
    end if

    write(*,*)"------------------------------"
    write(*,*)"       TIME INTEGRATOR"
    write(*,*)"------------------------------"
    ! Time integrator
    IF (TIME_INTEGRATOR_NAME .EQ. "ForwardEuler") THEN
      FORWARD_EULER_BOOL = .True.
      PRINT *,"ForwardEuler"
    ELSE IF (TIME_INTEGRATOR_NAME .EQ. "MidpointEuler") THEN
      MIDPOINT_EULER_BOOL = .True.
      PRINT *,"MidpointEuler"
    ELSE IF (TIME_INTEGRATOR_NAME .EQ. "ElectronsImplicit") THEN
      ELECTRONS_IMPLCIT_BOOL = .True.
      PRINT *,"Simulation with only electrons population being treated implicitly"
    ELSE
      WRITE(*,*) "Attention! Time integrator not recognized! Check input file. ABORTING!"
      STOP
    END IF

    write(*,*)"------------------------------"
    write(*,*)"       DUAL TIME STEPPING OPTION"
    write(*,*)"------------------------------"
    IF (ELECTRONS_IMPLCIT_BOOL) THEN
      PRINT *," Implicit simulation for electrons population only with dual time stepping strategy"
      !print *," Pseudo time step set to ", d_tau
    end if




    !if((ELECTRONS_IMPLCIT_BOOL .and. SYS_EULER_BOOL) .or. (ELECTRONS_IMPLCIT_BOOL .and. SPLIT_TS)) THEN

    !  write(*,*) "ERROR: Implicit simulation for electrons implemented for multifluid only! "
    !  STOP

    !end if


    write(*,*)"------------------------------"
    write(*,*)"       CONVERGENCE CRITERIA "
    write(*,*)"------------------------------"
    PRINT *, " Tolerance set on density in log10 scale ",toll

    if (type_debug .eq. "Debug")THEN
      DEBUG_BOOL = .True.
    end if

    if(type_control_ions .eq. "IONS_CHECK")THEN
      IONS_Temperature_Control = .True.
    end if




    dt = (t_end - t_start)/N_time
    t_old_ID = 0.0d0




    ! Allocate grid quantities and put them to zero
    N_int = N_cells - 1 ! First and last cells are ghost cells
    ALLOCATE(x_cc(N_cells))
    ALLOCATE(L_cells(N_cells))
    ALLOCATE(x_int(N_int))

    x_cc  = 0
    x_int = 0
    x_exhaust = 0.025d0

    if ((SYS_Euler_BOOL .and. FLUID_ELECTRONS ) .or. SYS_Euler_MF_BOOL) THEN

    write(*,*)"------------------------------"
    write(*,*)"       INITIAL SPECIES VALUE"
    write(*,*)"------------------------------"
    PRINT *, " ELECTRONS "
    PRINT *, " DENSITY  ",rho_init_electrons
    PRINT *, " TEMPERATURE  ",T_init_electrons
    PRINT *, " PRESSURE  ",p_init_electrons
    PRINT *, " VELOCITY X ",u_init_electrons
    PRINT *, " VELOCITY Y ",v_init_electrons

    end if

    if (SYS_Euler_MF_BOOL .or. (SYS_Euler_BOOL .and. IONS_Temperature_Control .and. FLUID_IONS))then

      PRINT *, " --------------------------------- "
      PRINT *, " IONS "
      PRINT *, " DENSITY  ",rho_init_ions
      PRINT *, " TEMPERATURE (FIXED) ",T_fixed_ions
      PRINT *, " PRESSURE  ",p_init_ions
      PRINT *, " VELOCITY X ",u_init_ions
      PRINT *, " VELOCITY Y ",v_init_ions

    else if (SYS_Euler_MF_BOOL .or. (SYS_Euler_BOOL .and. FLUID_IONS)) THEN

      PRINT *, " --------------------------------- "
      PRINT *, " IONS "
      PRINT *, " DENSITY  ",rho_init_ions
      PRINT *, " TEMPERATURE  ",T_init_ions
      PRINT *, " PRESSURE  ",p_init_ions
      PRINT *, " VELOCITY X ",u_init_ions
      PRINT *, " VELOCITY Y ",v_init_ions
    END IF

    PRINT *, " --------------------------------- "
    PRINT *, " NEUTRALS "
    print*,neutrals_density

    if ((SYS_Euler_BOOL .and. FLUID_ELECTRONS ) .or. SYS_Euler_MF_BOOL) THEN
    write(*,*)"------------------------------"
    write(*,*)"       ELECTRONS FLUID BC"
    write(*,*)"------------------------------"
    ! Boundary conditions type Right (fluids variables)
    IF (BCs_TYPE_Left_electrons .EQ. "open") THEN
      BCs_OPEN_BOOL_Left_electrons = .True.
      write(*,*) " Weak Dirichlet on the left "
      write(*,*)rho_el_l
      write(*,*)T_el_l
      write(*,*)u_el_l
      write(*,*)uy_el_l
    ELSE IF (BCs_TYPE_Left_electrons .EQ. "Neumann") THEN
      BCs_NEU_BOOL_Left_electrons = .True.
      write(*,*) " Weak Neumann on the left "
    else if(BCs_type_Left_electrons .eq. "periodic") THEN
      BCs_PERIODIC_BOOL_electrons = .True.
    ELSE
      WRITE(*,*) "Attention! BCs type not recognized! Check input file. ABORTING!"
      STOP
    END IF

    ! Boundary conditions type Right for electrons (fluids variables)
    IF (BCs_TYPE_Right_electrons .EQ. "open") THEN
      BCs_OPEN_BOOL_Right_electrons = .True.
      write(*,*) " Weak Dirichlet on the right "
      write(*,*)rho_el_r
      write(*,*)T_el_r
      write(*,*)u_el_r
      write(*,*)uy_el_r
    ELSE IF (BCs_TYPE_Right_electrons .EQ. "Neumann") THEN
      BCs_NEU_BOOL_Right_electrons = .True.
      write(*,*) " Weak Neumann on the right "
    else if(BCs_TYPE_Right_electrons .eq. "Robin")THEN
      BCs_ROBIN_BOOL_Right_electrons = .True.
      write(*,*) " Specific Robin condition for Hall Thruster engine at cathode: "
      write(*,*) " Current 0.6 Ampere enforced at cathode, adapting speed to number density and temperature values"
    else if(BCs_TYPE_Right_electrons .eq. "periodic") THEN
      if(BCs_PERIODIC_BOOL_electrons)then
        write(*,*) " Periodic boundary conditions set successfully "
      ELSE
        WRITE(*,*) "Attention! Periodic BCs can't be choosen on one side only "
        STOP
      end if
    ELSE
      WRITE(*,*) "Attention! BCs type not recognized! Check input file. ABORTING!"
      STOP
    END IF

end if


    if ((SYS_Euler_BOOL .and. FLUID_IONS ) .or. SYS_Euler_MF_BOOL) THEN



        write(*,*)"------------------------------"
        write(*,*)"       IONS FLUID BC"
        write(*,*)"------------------------------"

        ! Boundary conditions type Right (fluids variables)
        IF (BCs_TYPE_Left_ions .EQ. "open") THEN
          BCs_OPEN_BOOL_Left_ions = .True.
          write(*,*) " Weak Dirichlet on the left "
          if(IONS_Temperature_Control)THEN
            write(*,*)rho_i_l
            write(*,*)"(Fixed Temperature)",T_fixed_ions
            write(*,*)u_i_l
            write(*,*)uy_i_l
          else
            write(*,*)rho_i_l
            write(*,*)T_i_l
            write(*,*)u_i_l
            write(*,*)uy_i_l
          end if
        ELSE IF (BCs_TYPE_Left_ions .EQ. "Neumann") THEN
          BCs_NEU_BOOL_Left_ions = .True.
          write(*,*) " Weak Neumann on the left "
        else if(BCs_type_Left_ions .eq. "periodic") THEN
          BCs_PERIODIC_BOOL_ions = .True.
        ELSE
          WRITE(*,*) "Attention! BCs type not recognized! Check input file. ABORTING!"
          STOP
        END IF



    ! Boundary conditions type Right for ions (fluids variables)
    IF (BCs_TYPE_Right_ions .EQ. "open") THEN
      BCs_OPEN_BOOL_Right_ions = .True.
      write(*,*) " Weak Dirichlet on the right "
        if(IONS_Temperature_Control)THEN
          write(*,*)rho_i_r
          write(*,*)"(Fixed Temperature)",T_fixed_ions
          write(*,*)u_i_r
          write(*,*)uy_i_r
        else
          write(*,*)rho_i_r
          write(*,*)T_i_r
          write(*,*)u_i_r
          write(*,*)uy_i_r
        end if
    ELSE IF (BCs_TYPE_Right_ions .EQ. "Neumann") THEN
      BCs_NEU_BOOL_Right_ions = .True.
      write(*,*) " Weak Neumann on the right "
    else if(BCs_TYPE_Right_ions .eq. "periodic") THEN
      if(BCs_PERIODIC_BOOL_ions)then
        write(*,*) " Periodic boundary conditions set successfully "
      ELSE
        WRITE(*,*) "Attention! Periodic BCs can't be choosen on one side only "
        STOP
      end if
    ELSE
      WRITE(*,*) "Attention! BCs type not recognized! Check input file. ABORTING!"
      STOP
    END IF

  end if


  ! Initializing parameters for assembling block-diagonal double tensor
  ! using features of LIS library

  if(ELECTRONS_IMPLCIT_BOOL)THEN

    bnnz = 3*(N_cells-2) - 2

    allocate(bptr(0:N_cells-2))
    allocate(bindex(0:bnnz-1))
    allocate(value_vector(0:bnnz*16-1))


    ! Building bptr structure
    bptr(0) = 0
    bptr(1) = bptr(0) + 2

    do i=2,N_cells-3

      bptr(i) = bptr(i-1) + 3

    end do

    bptr(N_cells - 2)= bptr(N_cells - 3) + 2


    ! Building bindex structure

    bindex(0) = 0
    bindex(1) = 1

    triple_index_iii = 0

    do i=2,bnnz-3

      bindex(i) = triple_index_iii

      triple_index_iii = triple_index_iii + 1

      if (mod(i-1,3) .eq. 0) THEN

        triple_index_iii = triple_index_iii - 2

      end if


    end do

    bindex(bnnz-2) = N_cells-4
    bindex(bnnz-1) = N_cells-3

    allocate(LS_rhs((N_cells-2)*4))
    allocate(Q_n((N_cells-2)*4))
    allocate(Q_n_minus_1((N_cells-2)*4))

    bnnz = 3*(N_cells-2) - 2
    n_size = (N_cells-2) * 4
    bnr = 4
    bnc = 4
    nr =(n_size-1)/bnr+1

    ! Initializaing matrix
    call lis_matrix_create(comm,A_MATRIX,ierr)
    call lis_matrix_set_size(A_MATRIX,0,n_size,ierr)

    ! Initializing solution
    call lis_vector_duplicate(A_MATRIX,x_solution,ierr)
    call lis_vector_duplicate(A_MATRIX,IMPLICIT_RHS,ierr)



    ! Initializing solver
    call lis_solver_create(solver,ierr)
    call lis_precon_create(solver,my_prec,ierr)
    call lis_solver_set_option('-i 9 -p 2 -maxiter 100',solver,ierr)

    ALLOCATE(INTER_FLUX(N_EQ,4))



  end if


  allocate(d_tau_vect(N_cells))
  allocate(sol_n_minus1(N_eq,N_Cells))
  sol_n_minus1 = 0.0d0
  allocate(coll_elastic_vect_el(N_cells))
  coll_elastic_vect_el = 0.0d0
  allocate(coll_elastic_vect_io(N_cells))
  coll_elastic_vect_io = 0.0d0
  allocate(nu_ions_vect(3,N_cells))
  nu_ions_vect = 0.0d0

  allocate(Density_electrons_vector(N_cells))
  allocate(Temperature_electrons_vector(N_cells))

  Density_electrons_vector = Density_electrons
  Temperature_electrons_vector = Temperature_electrons

  P_converter_ions = T_fixed_ions/GAS_M_2*k_boltzmann
  P_converter_electrons = T_fixed_electrons/GAS_M_1*k_boltzmann


  ! Building matrix of elastic-scattering rates vs Temperature
  if (coll_bool_elastic) THEN

    T_samples = 1000;

    allocate(T_K_elastic(T_samples,2),T_K_elastic_ions(T_samples,2))

    OPEN(unit=7,file='elastic_rates.dat', status='old', action='read')

    do i=1,T_samples

      read(7,*) T_K_elastic(i,:)
      !write(*,*) T_K_elastic(i,:)

    end do

    close(7)

    OPEN(unit=10,file='elastic_ions_rates.dat', status='old', action='read')

    do i=1,T_samples

      read(10,*) T_K_elastic_ions(i,:)
      !write(*,*) T_K_elastic_ions(i,:)

    end do

    close(10)



  else if (coll_bool_full) then

    T_samples = 1000;

    allocate(T_K_elastic(T_samples,2),T_K_ionized(T_samples,2),T_K_excited(T_samples,2),T_K_elastic_ions(T_samples,2))

    OPEN(unit=7,file='elastic_rates.dat', status='old', action='read')

    do i=1,T_samples

      read(7,*) T_K_elastic(i,:)
      !write(*,*) T_K_elastic(i,:)

    end do

    close(7)

    OPEN(unit=8,file='ionization_rates.dat', status='old', action='read')

    do i=1,T_samples

      read(8,*) T_K_ionized(i,:)
      !write(*,*) T_K_ionized(i,:)

    end do

    close(8)

    OPEN(unit=9,file='excitation_rates.dat', status='old', action='read')

    do i=1,T_samples

      read(9,*) T_K_excited(i,:)
      !write(*,*) T_K_excited(i,:)

    end do

    close(9)

    OPEN(unit=10,file='elastic_ions_rates.dat', status='old', action='read')

    do i=1,T_samples

      read(10,*) T_K_elastic_ions(i,:)
      !write(*,*) T_K_elastic_ions(i,:)

    end do

    close(10)



   end if

  END SUBROUTINE INIT_VARIOUS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE INITIALIZATION
