MODULE EM_FIELDS

  USE GLOBAL

  IMPLICIT NONE

  ! Electrical and magnetic fields
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Ex_field, Ey_field, B_field, rho_Q_nodes ,phi,E_SC
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: En_plus_1_plot,E_SC_plot

  CONTAINS

  SUBROUTINE INIT_EM_FIELDS

  implicit none

  real(kind=8) :: ex_in,ey_in,B0

  ! Allocate electric and magnatic fields on the grid
  ALLOCATE(Ex_field(N_cells))
  ALLOCATE(Ey_field(N_cells))
  ALLOCATE(B_field(N_cells))
  allocate(rho_Q_nodes(N_cells))
  allocate(phi(N_int))
  ALLOCATE(En_plus_1_plot(N_cells))
  allocate(E_SC_plot(N_cells))
  allocate(E_SC(N_cells-2))
  E_SC_plot = 0.0d0


   ! Read EM settings file
  OPEN(unit=94,file='EM_fields.dat',status='old',action='read')

     READ(94,*)Ex_in
     READ(94,*)Ey_in
     READ(94,*)B0
     READ(94,*)phi_0
     READ(94,*)phi_end
     READ(94,*)BCs_TYPE_Left_EM
     READ(94,*)BCs_TYPE_Right_EM
     READ(94,*)Neu_left
     READ(94,*)Neu_right
     READ(94,*)Electrons_discharge_current_cathode
     READ(94,*)Engine_CS_Area


 CLOSE(94)

  Ex_field = Ex_in
  Ey_field = ey_in
  !B_field = B0/1.2d0*(exp(-16*(x_cc/x_exhaust - 1.0d0)**2) + 0.2d0*x_cc/x_exhaust)
  B_field = B0*(exp(-16*(x_cc/x_exhaust - 1.0d0)**2))

  rho_Q_nodes = 0
  phi = 0

  IF(SC_EM_field)THEN

    write(*,*)"------------------------------"
    write(*,*)"       POISSON SOLVER (Engine Features) "
    write(*,*)"------------------------------"
    ! Boundary conditions type Right (EM variables)
    IF (BCs_TYPE_Left_EM .EQ. "Dirichlet") THEN
      BCs_DIR_BOOL_Left_EM = .True.
      write(*,*) " LHS Dirichlet     ",phi_0," V"
    ELSE IF (BCs_TYPE_Left_EM .EQ. "Neumann") THEN
      BCs_NEU_BOOL_Left_EM = .True.
      write(*,*) " LHS Neumann"
    else if(BCs_TYPE_Left_EM .eq. "periodic") THEN
      BCs_PERIODIC_BOOL_EM = .True.
    ELSE
      WRITE(*,*) "Attention! BCs type not recognized! Check input file. ABORTING!"
      STOP
    END IF

    ! Boundary conditions type Right (EM variables)
    IF (BCs_TYPE_Right_EM .EQ. "Dirichlet") THEN
      BCs_DIR_BOOL_Right_EM = .True.
      write(*,*) "RHS Dirichlet     ",phi_end," V"
    ELSE IF (BCs_TYPE_Right_EM .EQ. "Neumann") THEN
      if(BCs_NEU_BOOL_Left_EM) THEN
        WRITE(*,*) "Problem not pointwise defined, numerics do NOT like it!"
        WRITE(*,*) "Go and change the code yourself if you want to take this risk!"
        STOP
      else
        BCs_NEU_BOOL_Right_EM = .True.
        write(*,*) " RHS Neumann"
      end if
    else if(BCs_TYPE_Right_EM .eq. "periodic") THEN
      if(BCs_PERIODIC_BOOL_EM)then
        write(*,*) " Periodic boundary conditions set successfully "
      ELSE
        WRITE(*,*) "Attention! Periodic BCs can't be choosen on one side only "
        STOP
      end if
    ELSE
      WRITE(*,*) "Attention! BCs type not recognized! Check input file. ABORTING!"
      STOP
    END IF
    print*," Magnetic field at exhaust: ",B0, " [Tesla]"
    print*," Discharge current for electrons at Cathode: ",Electrons_discharge_current_cathode," [Ampere]"
    print*," Cross sectional area for engine: ",Engine_CS_Area," [m^2]"
    print*," -------------------------------------------------------------------------------------------"


  END IF

  END SUBROUTINE INIT_EM_FIELDS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE COMPUTE_E_FIELD

  IMPLICIT NONE

  ! Hard-coded parameters
  REAL(KIND=8) :: E_SIN_AMPL      = 1000.0d0
  REAL(KIND=8) :: E_SIN_k         = 2513.2741228718d0
  REAL(KIND=8) :: E_SIN_omega     = 0
  REAL(KIND=8) :: E_SIN_PHASE_rad = 0

  REAL(KIND=8) :: B_field_uniform = 0.02 ! [Tesla]

  ! INTEGER :: IC

  ! Compute field at current time
  ! Ex_field = 10000.0d0*EXP(- 5.0d0*100000.0d0*(x_cc - 0.005)**2) ! Gaussian
  ! Ex_field = E_SIN_AMPL * SIN( E_SIN_k*x_cc - E_SIN_omega*CURRENT_TIME + E_SIN_PHASE_rad ) ! sin
  ! Ex_field = E_SIN_AMPL

  ! Set hard-coded value of Electric Field:
  ! Indeed within Hall-Thruster simulation, only Electric Field along y axis
  ! would be used as uniform, the one one the x always come from potential gap

  Ex_field = 0.0d00
  !Ex_field = 12000.0d00
  Ey_field = 0.0d0
  B_field = 0.0d0

!   B_field = 0.01d0
!   B_field = B_field_uniform

   B_field = 0.00d0

   if(Nondim_bool)THEN

     B_field = B_field/B_mag_carac

     Ex_field = Ex_field/E_carac

     Ey_field = Ey_field/Ey_carac

   end if


  END SUBROUTINE

END MODULE EM_FIELDS
