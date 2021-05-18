MODULE GLOBAL
  ! This module contains the most important variables of the program.
  ! All other module will load this one.

  IMPLICIT NONE

  CHARACTER(LEN=512) :: SYS_NAME
  LOGICAL :: SYS_EULER_BOOL = .False. ! Init to false!
  LOGICAL :: SYS_MOM5_BOOL  = .False. ! Init to false!
  LOGICAL :: SYS_MOM14_BOOL = .False. ! Init to false!
  LOGICAL :: SYS_MOM14_AXI_BOOL = .False. ! Init to false!

  ! Thermodynamic variables
  REAL(KIND=8) :: GAS_GAMMA, GAS_M, GAS_Q


  ! Time integration
  INTEGER :: N_time
  REAL(KIND=8) :: CURRENT_TIME
  REAL(KIND=8) :: t_start, t_end, dt
  CHARACTER(LEN=512) :: TIME_INTEGRATOR_NAME
  LOGICAL :: FORWARD_EULER_BOOL  = .False. ! Init to false!
  LOGICAL :: MIDPOINT_EULER_BOOL = .False. ! Init to false!
  real(kind=8),dimension(:),allocatable :: coll_elastic_vect_el,coll_elastic_vect_io
  real(kind=8) :: nu_implicit_el,nu_implicit_io,alfa_1,alfa_2


  ! MANGOTURB

  LOGICAL :: SC_EM_field = .False.
  REAL(kind=8) :: phi_0,phi_end,neutrals_density,k_boltzmann,eps_0,threshold_energy_xe,excitation_threshold
  REAL(kind=8) :: Neu_left, Neu_right
  CHARACTER(LEN=512) :: type_e,type_init,type_coll,type_dim,type_dt,type_timestep,flux_scheme,type_splitting
  CHARACTER(LEN=512) :: type_debug,type_Control_ions,type_control_electrons
  logical :: cont_bool = .False.
  logical :: disp_bool = .False.
  logical :: coll_bool_full = .False.
  logical :: coll_bool_elastic = .False.
  logical :: coll_boolA = .False.
  logical :: single_bool = .False.
  logical :: Nondim_bool = .False.
  logical :: NL_D_bool = .False.
  logical :: test_case_1 = .False.
  LOGICAL :: SYS_EULER_MF_BOOL = .False.
  LOGICAL :: ELECTRONS_IMPLCIT_BOOL = .False.
  logical :: CFL_BOOL = .False.
  INTEGER :: KL,KU,LDAB,INFO
  REAL(KIND=8) :: GAS_GAMMA_1, GAS_M_1, GAS_Q_1
  REAL(KIND=8) :: GAS_GAMMA_2, GAS_M_2, GAS_Q_2
  REAL(KIND=8) :: u_neutrals,Energy_neutrals, T_neutrals
  INTEGER :: T_samples
  real(kind=8) :: t_old_ID,CFL_LIMIT_EXPLICIT
  ! Quantities needed for nondimensional analysis: characteristics variables
  real(kind=8) :: rho_el_carac, rho_i_carac, u_el_carac, u_i_carac, uy_el_carac, uy_i_carac, P_el_carac, P_i_carac
  real(kind=8) :: x_carac, t_carac
  real(kind=8) :: phi_carac, E_carac, B_mag_carac,Ey_carac
  real(kind=8) :: nu_coll_carac
  real(kind=8) :: Energy_el_carac, Energy_i_carac
  ! Quantities needed for nondimensional analysis: nondimensional quantities
  real(kind=8) :: relaxation_factor,NUM_CFL_EVOLUTION_FACTOR
  integer :: dual_iteration,global_iter
  logical :: HLLC_BOOL = .False.
  logical :: HLL_BOOL = .False.
  integer :: ions_subiterations
  logical :: TS_BOOL = .False.
  logical :: DEBUG_BOOL = .False.
  real(kind=8) :: CFL_PHYS_MAX
  real(kind=8) :: n_el_l,n_el_r,n_i_l,n_i_r,n_init_ions,n_init_electrons
  logical :: NaN_BOOL = .False.
  logical :: BCs_ROBIN_BOOL_Right_electrons = .False.
  real(kind=8) :: Electrons_discharge_current_cathode,Engine_CS_Area
  real(kind=8), dimension(:),allocatable :: Density_electrons_vector, Temperature_electrons_vector
  real(kind=8) :: Density_electrons, Temperature_electrons
  real(kind=8) :: T_fixed_ions = 1100.0d0 ![K]
  real(kind=8) :: T_fixed_electrons
  real(kind=8) :: P_converter_ions, P_converter_electrons
  logical :: IONS_Temperature_Control = .False.
  logical :: ELECTRONS_Temperature_Control = .False.
  logical :: FLUID_ELECTRONS = .False.
  logical :: FLUID_IONS = .False.



  ! Team electrons

  real(kind=8) :: beta = 1.0d0 ! nu_coll_carac*t_carac
  real(kind=8) :: E_x_e_AD = 1.0d0 ! GAS_Q_1 * E_carac * x_carac / GAS_M_1 / u_el_carac /u_el_carac
  real(kind=8) :: B_x_e_AD = 1.0d0 ! GAS_Q_1/GAS_M_1*b_mag_carac*uy_el_carac/u_el_carac*t_carac
  real(kind=8) :: kb8 = 1.0d0 ! P_el_carac/rho_el_carac/u_el_carac/u_el_carac
  real(kind=8) :: jrswish = 1.0d0 ! GAS_Q_1/GAS_M_1*Ey_carac/uy_el_carac*t_carac
  real(kind=8) :: B_y_e_AD = 1.0d0 ! GAS_Q_1/GAS_M_1*b_mag_carac*t_carac*u_el_carac/uy_el_carac
  real(kind=8) :: cp3 = 1.0d0 ! nu_coll_carac*t_carac/GAS_M_1/Energy_el_carac
  real(kind=8) :: E_x_e_ene_AD = 1.0d0 ! GAS_Q_1 * E_carac * x_carac / GAS_M_1 / Energy_el_carac
  real(kind=8) :: kcp_y = 1.0d0 !  GAS_Q_1/GAS_M_1*Ey_carac*t_carac*uy_el_carac/Energy_el_carac
  real(kind=8) :: kp = 1.0d0 ! P_el_carac/rho_el_carac/Energy_el_carac

  ! Team ions

  real(kind=8) :: u_neutrals_to_ions_ratio = 1.0d0 ! 0.5*u_neutrals*u_neutrals + 3.0d0/2.0d0*k_boltzmann*T_neutrals/GAS_M_2
  real(kind=8) :: Energy_neutrals_to_ions_ratio =1.0d0 ! Energy_neutrals/Energy_i_carac
  real(kind=8) :: tmc = 1.0d0 ! rho_el_carac/rho_i_carac*GAS_M_2/GAS_M_1
  real(kind=8) :: E_x_i_AD = 1.0d0 ! GAS_Q_2 * E_carac * x_carac / GAS_M_2 / u_i_carac /u_i_carac
  real(kind=8) :: B_x_i_AD = 1.0d0 ! GAS_Q_2/GAS_M_2*b_mag_carac*t_carac*uy_i_carac/u_i_carac
  real(kind=8) :: kb24 = 1.0d0 ! P_i_carac/rho_i_carac/u_i_carac/u_i_carac
  real(kind=8) :: rj = 1.0d0 ! GAS_Q_2/GAS_M_2*Ey_carac*t_carac/uy_i_carac
  real(kind=8) :: B_y_i_AD = 1.0d0 ! GAS_Q_2/GAS_M_2*b_mag_carac*t_carac*u_i_carac/uy_i_carac
  real(kind=8) :: E_x_i_ene_AD = 1.0d0 ! GAS_Q_2 * E_carac * x_carac / GAS_M_2 / Energy_i_carac
  real(kind=8) :: kg_y = 1.0d0 ! GAS_Q_2/GAS_M_2*Ey_carac*t_carac*uy_i_carac*t_carac/Energy_i_carac
  real(kind=8) :: thj = 1.0d0 ! P_i_carac/rho_i_carac/Energy_i_carac

  ! Poisson refs
  real(kind=8) :: kobe = 1.0d0 ! rho_el_carac/rho_i_carac*GAS_M_2/GAS_M_1
  real(kind=8) :: shaq = 1.0d0 ! phi_carac*GAS_M_2/rho_i_carac*eps_0/GAS_Q/x_carac/x_carac  ! Not fully non-dimensional
  REAL(KIND=8) :: Predictor_1,Predictor_2,Predictor_3,Predictor_4


  ! Initial conditions
  REAL(KIND=8) :: rho_init_electrons, P_init_electrons, u_init_electrons, T_init_electrons,v_init_electrons
  REAL(KIND=8) :: rho_init_ions, P_init_ions, u_init_ions, T_init_ions,v_init_ions
  real(kind=8) :: n_ions_background
  real(kind=8) , allocatable, DIMENSION(:) :: BC_electrons_right,BC_electrons_left,BC_ions_right,BC_ions_left
  REAL(KIND=8) :: rho_el_r,rho_i_r,T_el_r,T_i_r,u_el_r,u_i_r,rho_el_l,rho_i_l,T_el_l,T_i_l,u_el_l,u_i_l
  REAL(KIND=8) :: uy_el_r,uy_i_r,uy_el_l,uy_i_l
  real(Kind=8) :: CFL_INPUT,res_rho_i,res_ene_i,res_ene_el,res_rho_el,toll
  real(kind=8) :: d_tau,nu_implicit_el_iz,nu_implicit_el_ex
  real(kind=8) :: CFL_max_ions,CFL_MAX_electrons,CFL_INP_IONS,CFL_INP_ELECTRONS,CFL_NUM_EL
  real(kind=8), dimension(:,:),allocatable :: nu_ions_vect



  ! Variables
  INTEGER :: N_EQ ! Number of equations (system size)

  ! Numerical
  REAL(KIND=8) :: CFL_max
  CHARACTER :: sigma_small_flag = ' '

  REAL(KIND=8) :: sig_min, sig_max

  ! Boundary conditions
  LOGICAL :: BCs_OPEN_BOOL_Left_electrons  = .False. ! Init to false!
  LOGICAL :: BCs_OPEN_BOOL_Right_electrons = .False.
  LOGICAL :: BCs_NEU_BOOL_Left_electrons   = .False. ! Init to false!
  LOGICAL :: BCs_NEU_BOOL_Right_electrons = .False.
  LOGICAL :: BCs_OPEN_BOOL_Left_ions  = .False. ! Init to false!
  LOGICAL :: BCs_OPEN_BOOL_Right_ions = .False.
  LOGICAL :: BCs_NEU_BOOL_Left_ions   = .False. ! Init to false!
  LOGICAL :: BCs_NEU_BOOL_Right_ions = .False.
  LOGICAL :: BCs_NEU_BOOL_Left_EM   = .False. ! Init to false!
  LOGICAL :: BCs_NEU_BOOL_Right_EM = .False.
  LOGICAL :: BCs_DIR_BOOL_Left_EM   = .False. ! Init to false!
  LOGICAL :: BCs_DIR_BOOL_Right_EM = .False.
  LOGICAL :: BCs_PERIODIC_BOOL_EM = .False.
  LOGICAL :: BCs_PERIODIC_BOOL_electrons = .False.
  LOGICAL :: BCs_PERIODIC_BOOL_ions = .False.
  LOGICAL :: BOOL_RESTART = .False.
  CHARACTER(LEN=512) :: BCs_TYPE_Left_electrons,BCs_TYPE_Right_electrons,BCs_TYPE_Left_EM,BCs_TYPE_Right_EM,type_restart
  CHARACTER(LEN=512) :: BCs_TYPE_Left_ions,BCs_TYPE_Right_ions

  ! Spatial scheme
  INTEGER :: RECONSTRUCTION_ORDER ! 0: first order scheme, 1: linear reconstruction -> second order

  ! Grid
  INTEGER :: N_cells, N_int
  CHARACTER(LEN=512) :: GRID_TYPE
  REAL(KIND=8) :: x_min, x_max, x_exhaust
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: x_cc, x_int, L_cells,d_tau_vect ! Center cells and interfaces


  ! Dump
  INTEGER :: WRITE_EVERY

END MODULE GLOBAL
