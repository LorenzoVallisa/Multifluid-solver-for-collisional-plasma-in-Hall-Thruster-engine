MODULE PDE_EULER

  USE GLOBAL
  use PDE_EULER_ELECTRONS
  use PDE_EULER_IONS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!! Euler subroutines !!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE EULER_INIT_SOL(sol)

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ, N_cells), INTENT(INOUT) :: sol

  if (FLUID_ELECTRONS)THEN

    CALL ELECTRONS_INIT_SOL(sol)

  ELSE IF (FLUID_IONS)Then

    CALL IONS_INIT_SOL(sol)

  end if


  END SUBROUTINE EULER_INIT_SOL


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine compute_Ustar_left_right(sol,S,S_star,U_star)

  implicit none

  REAL(KIND=8), DIMENSION(N_EQ), intent(in) :: sol
  REAL(KIND=8), DIMENSION(N_EQ), intent(inout) :: u_star
  REAL(KIND=8), intent(in) :: S,S_star
  real(kind=8) :: rho,ux,uy,E,P
  real(kind=8) :: rhoE,rhoux,rhouy

  rho    = sol(1)
  rhoux  = sol(2)
  rhouy  = sol(3)
  ux     = rhoux/(rho + 1.0e-35)
  uy     = rhouy/(rho + 1.0e-35)
  rhoE   = sol(4)
  E = rhoE/rho
  P      = (GAS_GAMMA - 1)*(rhoE - rho*ux*ux/2.0 - rho*uy*uy/2.0)

  u_star(1) = 1

  u_star(2) = S_star

  u_star(3) = uy

  u_star(4) = E + (S_star - ux) * (S_star + p/(rho * (S-ux)))

  u_star = u_star * (rho * (S - ux)/(S - S_star))

  end subroutine compute_Ustar_left_right

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE EULER_FLUXES(sol, Flux)

    ! Computes flux from a solution vector, of size N_EQ

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
    REAL(KIND=8), DIMENSION(N_EQ), INTENT(OUT) :: Flux

    REAL(KIND=8) :: rho, rhoux, rhouy, ux, uy, rhoE, P

    ! Unpack the solution and compute primitive variables
    rho    = sol(1)
    rhoux  = sol(2)
    rhouy  = sol(3)
    ux     = rhoux/(rho + 1.0e-35)
    uy     = rhouy/(rho + 1.0e-35)
    rhoE   = sol(4)
    P      = (GAS_GAMMA - 1)*(rhoE - rhoux**2/rho/2.0d0 - rhouy**2/rho/2.0d0)


    Flux(1) = rhoux
    Flux(2) = rhoux**2/rho + P
    Flux(3) = rhoux*rhouy/rho
    Flux(4) = rhoE*ux + P*ux



  END SUBROUTINE EULER_FLUXES

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE EULER_SOURCES(sol, Ex, Ey, B, Src_EM, Src_COLL)
  ! Computes all sources for Euler

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8),                  INTENT(IN)  :: Ex, Ey, B
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(OUT) :: Src_COLL
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(OUT) :: Src_EM

  Src_COLL = 0.0d0
  Src_EM = 0.0d0

  if (FLUID_ELECTRONS)THEN

    CALL ELECTRONS_EM_SOURCES(sol, Ex, Ey, B, Src_EM)

    CALL ELECTRONS_COLLISIONAL_SOURCES(sol,Src_COLL)

  ELSE IF (FLUID_IONS)Then

    CALL iONS_EM_SOURCES(sol, Ex, Ey, B, Src_EM)

    CALL IONS_COLLISIONAL_SOURCES(sol,Src_COLL)

  end if

  END SUBROUTINE EULER_SOURCES

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE EULER_MIN_MAX_WAVESPEEDS(sol, sMin, sMax)

  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8),                  INTENT(OUT) :: sMin, sMax

  REAL(KIND=8) :: rho, rhoux, rhouy, ux, uy, rhoE, P

  ! Unpack the solution and compute primitive variables
  rho    = sol(1)
  rhoux  = sol(2)
  rhouy  = sol(3)
  ux     = rhoux/(rho + 1.0e-35)
  uy     = rhouy/(rho + 1.0e-35)
  rhoE   = sol(4)
  P      = (GAS_GAMMA - 1)*(rhoE - rho*ux*ux/2.0 - rho*uy*uy/2.0)

  sMin = ux - SQRT(GAS_GAMMA*P/(rho + 1.0e-35))
  sMax = ux + SQRT(GAS_GAMMA*P/(rho + 1.0e-35))

  END SUBROUTINE EULER_MIN_MAX_WAVESPEEDS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE PDE_EULER
