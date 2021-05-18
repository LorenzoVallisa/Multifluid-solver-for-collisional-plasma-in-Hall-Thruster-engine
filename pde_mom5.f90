MODULE PDE_MOM5

  !################################################################################
  !
  ! This module implements the axisymmetric version of the 14 moment equations.
  ! Exploiting symmetries, the solution is made by 6 equations.
  ! This works for the 1D shockwave for example, where the VDF along x is 
  ! Mott-Smith-like, while along y and z is Maxwellian.
  ! 
  ! ------
  !
  ! The theory for the 14-moments problem and the interpolative closure can be 
  ! found in:
  !
  ! McDonald & Torrilhon, Affordable robust moment closures for CFD based on the 
  ! maximum-entropy hierarchy, Journal of Computational Physics (2013).
  ! 
  ! More theory in:
  ! 
  ! Levermore, Moment closure hierarchies for kinetic theories, Journal of
  ! Statistical Physics (1996).
  !
  !################################################################################

  USE GLOBAL
  USE VARIOUS

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE MOM5_INIT_SOL(sol)
  ! Initializes the solution vector for the 14-moments PDEs. 
  ! HARD-CODED FOR NOW!

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ, N_cells), INTENT(INOUT) :: sol

  REAL(KIND=8) :: rho5, u5, P5, q5, R5

  ! ============= Initialize solution (for 14 moments) ===============

  ! Uniform initial solution
  rho5 = 9.1094E-11
  u5   = 0.0d0
  ! P5   = 138.0d0
  P5   = 0.1380d0
  q5   = 0.0d0
  r5   = 3.0d0*(P5**2)/rho5 

  ! Conserved variables 
  sol(1,:) = rho5
  sol(2,:) = rho5*u5   
  sol(3,:) = rho5*u5**2 + P5
  sol(4,:) = rho5*u5**3 + 3*u5*P5 + q5
  sol(5,:) = rho5*u5**4 + 6*u5**2*P5 + 4*u5*q5 + r5

!!! SOD SHOCK STATE !!!  ! Left state
!!! SOD SHOCK STATE !!!  rho5 = 4.696d0
!!! SOD SHOCK STATE !!!  u5   = 500.0d0
!!! SOD SHOCK STATE !!!  P5   = 404400.0d0
!!! SOD SHOCK STATE !!!  q5   = 0.0d0
!!! SOD SHOCK STATE !!!  r5   = 3.0d0*P5**2/rho5
!!! SOD SHOCK STATE !!!
!!! SOD SHOCK STATE !!!  sol(1,:) = rho5
!!! SOD SHOCK STATE !!!  sol(2,:) = rho5*u5   
!!! SOD SHOCK STATE !!!  sol(3,:) = rho5*u5**2 + P5
!!! SOD SHOCK STATE !!!  sol(4,:) = rho5*u5**3 + 3*u5*P5 + q5
!!! SOD SHOCK STATE !!!  sol(5,:) = rho5*u5**4 + 6*u5**2*P5 + 4*u5*q5 + r5
!!! SOD SHOCK STATE !!!
!!! SOD SHOCK STATE !!!  ! Right state
!!! SOD SHOCK STATE !!!  ! rho5 = 1.408d0
!!! SOD SHOCK STATE !!!  u5   = -u5
!!! SOD SHOCK STATE !!!  ! P5   = 101100.0d0
!!! SOD SHOCK STATE !!!  ! q5   = 0.0d0
!!! SOD SHOCK STATE !!!  ! r5   = 3.0d0*P5**2/rho5
!!! SOD SHOCK STATE !!!
!!! SOD SHOCK STATE !!!  sol(1,FLOOR(N_cells/2.0d0):) = rho5
!!! SOD SHOCK STATE !!!  sol(2,FLOOR(N_cells/2.0d0):) = rho5*u5   
!!! SOD SHOCK STATE !!!  sol(3,FLOOR(N_cells/2.0d0):) = rho5*u5**2 + P5
!!! SOD SHOCK STATE !!!  sol(4,FLOOR(N_cells/2.0d0):) = rho5*u5**3 + 3*u5*P5 + q5
!!! SOD SHOCK STATE !!!  sol(5,FLOOR(N_cells/2.0d0):) = rho5*u5**4 + 6*u5**2*P5 + 4*u5*q5 + r5

  END SUBROUTINE MOM5_INIT_SOL

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  SUBROUTINE MOM5_FLUXES(sol, Flux)

  ! Computes flux from a solution vector, of size N_EQ

  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(OUT) :: Flux

  REAL(KIND=8) :: rho5, rhou5, u5, P5, q5, r5, SIG, s5
!  REAL(KIND=8) :: sig, Qxxx, Qxrr, Rxxjj, Sxiijj

  ! Reconstruct state
  rho5  = sol(1) + 1.0e-35
  rhou5 = sol(2)
  u5    = rhou5/rho5
  P5    = sol(3) - rhou5*u5
  q5    = sol(4) - (rhou5*u5*u5 + 3*u5*P5)
  r5    = sol(5) - (rho5*u5**4 + 6*u5**2*P5 + 4*u5*q5)

  ! Compute closing moments
  CALL MOM5_CLOSING_INTERPOLATIVE(rho5, P5, q5, r5, SIG, s5)

  Flux(1) = sol(2)
  Flux(2) = sol(3)
  Flux(3) = sol(4)
  Flux(4) = sol(5)
  Flux(5) = rhou5*u5*u5*u5*u5 + 10*u5**3*P5 + 10*u5**2*q5 + 5*u5*r5 + s5

  END SUBROUTINE MOM5_FLUXES

  ! ****************************************************************************************

  SUBROUTINE MOM5_SOURCES(sol, Ex, Ey, B, Src)

  ! Computes the source terms for this PDE. 
  ! Can be of different nature.

  IMPLICIT NONE
  
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8),                  INTENT(IN)  :: Ex, Ey, B
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(OUT) :: Src

  REAL(KIND=8), DIMENSION(N_EQ) :: Src_BGK, Src_EM

  CALL MOM5_BGK_SOURCES(sol, Src_BGK) 
  CALL MOM5_EM_SOURCES(sol, Ex, Ey, B, Src_EM)

  Src = Src_BGK + Src_EM

  END SUBROUTINE MOM5_SOURCES

  ! ****************************************************************************************

  SUBROUTINE MOM5_BGK_SOURCES(sol, BGK_Src)
  ! Computes BGK-style collision sources

  IMPLICIT NONE
  
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(OUT) :: BGK_Src
 
  ! REAL(KIND=8) :: Tau = 1.0e-11 ! HARD-CODED, mean time between collisions, 1/collis_freq
  REAL(KIND=8) :: Tau = 100.0d0 ! HUGE! ! HARD-CODED, mean time between collisions, 1/collis_freq 

  BGK_Src = 0.0d0
  
  END SUBROUTINE MOM5_BGK_SOURCES

  ! ****************************************************************************************

  SUBROUTINE MOM5_EM_SOURCES(sol, Ex, Ey, B, EM_SRC)
  ! Computes EM sources
 
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8),                  INTENT(IN)  :: Ex, Ey, B
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(OUT) :: EM_SRC

  EM_SRC(1) = 0.0d0
  EM_SRC(2) =   GAS_Q/GAS_M*Ex*sol(1)
  EM_SRC(3) = 2*GAS_Q/GAS_M*Ex*sol(2)
  EM_SRC(4) = 3*GAS_Q/GAS_M*Ex*sol(3)
  EM_SRC(5) = 4*GAS_Q/GAS_M*Ex*sol(4)

  END SUBROUTINE MOM5_EM_SOURCES

  ! ****************************************************************************************

  SUBROUTINE MOM5_CLOSING_INTERPOLATIVE(rho, P, q, r, SIG, s)
  ! Computes interpolative closing moments for the 1D1V 5-mom system

  IMPLICIT NONE
  
  REAL(KIND=8), INTENT(IN)  :: rho, P, q, R
  REAL(KIND=8), INTENT(OUT) :: SIG, s

  REAL(KIND=8) :: qNON, rNON, sNON, SIG_LOW_LIMIT, POW_SIG

  ! Non-dimensionalize variables
  qNON = q/rho*(rho/P)**1.5d0 
  rNON = r/rho*(rho/P)**2.0d0 

  CALL COMPUTE_SIGMA(qNON, rNON,  SIG)
 
  ! POW_SIG = 3.0d0/5.0d0 ! Old paper
  POW_SIG = 1.0d0/2.0d0

  sNON = qNON**3/SIG**2 + (10.0d0 - 8*SIG**POW_SIG)*qNON

  s = sNON*rho*(P/rho)**2.5d0
  
  END SUBROUTINE MOM5_CLOSING_INTERPOLATIVE
 
  ! *****************************************************************************************

  SUBROUTINE COMPUTE_SIGMA(qNON, rNON,  SIG)
  ! Computes sigma from nondimensional variables qNON and rNON

  IMPLICIT NONE

  REAL(KIND=8), INTENT(IN)  :: qNON, rNON
  REAL(KIND=8), INTENT(OUT) :: SIG

  REAL(KIND=8) :: SIG_LOW_LIMIT

  ! Compute sigma and limit it artificially
  SIG = (3 - rNON + SQRT( (3 - rNON)**2 + 8*qNON**2 ) )/4.0d0

  ! Artificial limit to sigma
  SIG_LOW_LIMIT = 1.0d-5
  sigma_small_flag = ' '

  IF (SIG .LT. SIG_LOW_LIMIT) THEN
    sigma_small_flag = '*'
    SIG = MAX(SIG, SIG_LOW_LIMIT) ! Artificial limitation on sigma
  END IF

  END SUBROUTINE COMPUTE_SIGMA

  ! *****************************************************************************************

  SUBROUTINE MOM5_MIN_MAX_WAVESPEEDS(sol, sMin, sMax)

  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8),                  INTENT(OUT) :: sMin, sMax

  REAL(KIND=8) :: rho5, rhou5, u5, P5, q5, r5, a5, MULT_CONST
  REAL(KIND=8) :: GAMMA_5_MOM

  REAL(KIND=8) :: A, B, C, D, E, X, Y ! working variables
  REAL(KIND=8) :: SIG, r, q

  ! Reconstruct state
  rho5  = sol(1) + 1.0e-35
  rhou5 = sol(2)
  u5    = rhou5/rho5
  P5    = sol(3) - rhou5*u5
  q5    = sol(4) - (rhou5*u5*u5 + 3*u5*P5)
  r5    = sol(5) - (rho5*u5**4 + 6*u5**2*P5 + 4*u5*q5)

  !!!!  ! We arbitrarily scale the speed of sound, by a constant.
  !!!!  ! Using something in the order of 1 - 10 should do, but notice that it can really 
  !!!!  ! be much different, depending on the problem! 
  !!!!  ! This will introduce unphysical dissipation anyway, so you should at least try different 
  !!!!  ! values.
  !!!!  MULT_CONST = 2.0d0 
  !!!!  
  !!!!  GAMMA_5_MOM = 3.0d0 ! Gamma is 3 for 1V gases!
  !!!!  a5     = SQRT(GAMMA_5_MOM*P5/rho5) ! Speed of sound
  !!!! 
  !!!!  sMin = u5 - a5*MULT_CONST
  !!!!  sMax = u5 + a5*MULT_CONST

  ! Non-dimensionalize variables
  q = q5/rho5*(rho5/P5)**1.5d0 
  r = r5/rho5*(rho5/P5)**2.0d0 

  CALL COMPUTE_SIGMA(q, r,  sig) ! q,r nondimensional

  ! Use fits for the wavespeeds, from Baradaran (see Baradaran's master thesis or paper with J. McDonald)
  A = 5 - 4*sqrt(sig) - sqrt(10 - 16*sqrt(sig) + 6*sig)
  B = 5 - 4*sqrt(sig) + sqrt(10 - 16*sqrt(sig) + 6*sig)
  C = sqrt(3 - 3*sig)
  D = 3.0d0/10.0d0*sqrt(3 - 3*sig)
  E = 8.0d0/10.0d0*sqrt(3 - 3*sig)
  X = A + D**2 + 2*sqrt(A)*D
  Y = B + E**2 - 2*sqrt(B)*E

  ! Lambda 0 and Lambda 4 (in Baradaran's THESIS notation, not the paper) are the max and min eigenvalues
  ! respectively of the NONDIMENSIONAL Jacobian. They need to be scaled
  sMax = (q + sqrt(q**2 - 4.0d0/5.0d0*q*sig*C + 4*sig**2*Y))/(2.0d0*sig) + E
  sMin = (q - sqrt(q**2 + 4.0d0/5.0d0*q*sig*C + 4*sig**2*Y))/(2.0d0*sig) - E

  ! Rescale wavespeeds
  sMin = sMin*SQRT(P5/rho5) + u5
  sMax = sMax*SQRT(P5/rho5) + u5

  END SUBROUTINE MOM5_MIN_MAX_WAVESPEEDS

END MODULE PDE_MOM5
