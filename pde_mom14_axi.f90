MODULE PDE_MOM14_AXI

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

  SUBROUTINE MOM14_AXI_INIT_SOL(sol)
  ! Initializes the solution vector for the 14-moments PDEs. 
  ! HARD-CODED FOR NOW!

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ, N_cells), INTENT(INOUT) :: sol

  REAL(KIND=8) :: P14, rho14, ux14, Pxx14, Prr14, qx14, Riijj14

  ! ============= Initialize solution (for 14 moments) ===============

  ! Uniform initial solution
  rho14   = 9.1094E-11
  ux14    = 0.0d0
  P14     = 138.0d0
  Pxx14   = P14
  Prr14   = 2.0d0*P14 ! Pyy + Pzz
  qx14    = 0.0d0
  Riijj14 = 15.0d0*(P14**2)/rho14

  ! Conserved variables, left state
  sol(1,:) = rho14   
  sol(2,:) = rho14*ux14   
  sol(3,:) = rho14*ux14**2 + Pxx14
  sol(4,:) = Prr14
  sol(5,:) = rho14*ux14**3 + ux14*(3*Pxx14 + Prr14) + qx14
  sol(6,:) = rho14*ux14**4 + ux14**2*(6*Pxx14 + 2*Prr14) + 4*qx14 + Riijj14

  !!! SOD SHOCK !!! ! Initial solution for Riemann problem with two beams, one going left one right.
  !!! SOD SHOCK !!! rho14   = 4.696d0
  !!! SOD SHOCK !!! P14     = 404400.0d0
  !!! SOD SHOCK !!! ux14    = 0.0d0
  !!! SOD SHOCK !!! Pxx14   = P14
  !!! SOD SHOCK !!! Prr14   = 2.0d0*P14 ! Pyy + Pzz
  !!! SOD SHOCK !!! qx14    = 0.0d0
  !!! SOD SHOCK !!! Riijj14 = 15.0d0*(P14**2)/rho14

  !!! SOD SHOCK !!! ! Conserved variables, left state
  !!! SOD SHOCK !!! sol(1,:) = rho14   
  !!! SOD SHOCK !!! sol(2,:) = rho14*ux14   
  !!! SOD SHOCK !!! sol(3,:) = rho14*ux14**2 + Pxx14
  !!! SOD SHOCK !!! sol(4,:) = Prr14
  !!! SOD SHOCK !!! sol(5,:) = rho14*ux14**3 + ux14*(3*Pxx14 + Prr14) + qx14
  !!! SOD SHOCK !!! sol(6,:) = rho14*ux14**4 + ux14**2*(6*Pxx14 + 2*Prr14) + 4*qx14 + Riijj14

  !!! SOD SHOCK !!! ! Conserved variables, right state
  !!! SOD SHOCK !!! ! ux14 = - ux14

  !!! SOD SHOCK !!! rho14 = 1.408d0
  !!! SOD SHOCK !!! P14   = 101100.0d0

  !!! SOD SHOCK !!! ux14    = 0.0d0
  !!! SOD SHOCK !!! Pxx14   = P14
  !!! SOD SHOCK !!! Prr14   = 2.0d0*P14 ! Pyy + Pzz
  !!! SOD SHOCK !!! qx14    = 0.0d0
  !!! SOD SHOCK !!! Riijj14 = 15.0d0*(P14**2)/rho14

  !!! SOD SHOCK !!! sol(1,FLOOR(N_CELLS/2.0d0):) = rho14   
  !!! SOD SHOCK !!! sol(2,FLOOR(N_CELLS/2.0d0):) = rho14*ux14   
  !!! SOD SHOCK !!! sol(3,FLOOR(N_CELLS/2.0d0):) = rho14*ux14**2 + Pxx14
  !!! SOD SHOCK !!! sol(4,FLOOR(N_CELLS/2.0d0):) = Prr14
  !!! SOD SHOCK !!! sol(5,FLOOR(N_CELLS/2.0d0):) = rho14*ux14**3 + ux14*(3*Pxx14 + Prr14) + qx14
  !!! SOD SHOCK !!! sol(6,FLOOR(N_CELLS/2.0d0):) = rho14*ux14**4 + ux14**2*(6*Pxx14 + 2*Prr14) + 4*ux14*qx14 + Riijj14

  END SUBROUTINE MOM14_AXI_INIT_SOL

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  SUBROUTINE MOM14_AXI_FLUXES(sol, Flux)

  ! Computes flux from a solution vector, of size N_EQ

  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(OUT) :: Flux

  REAL(KIND=8) :: rho, rhoux, ux, Pxx, Prr, qx, Riijj
  REAL(KIND=8) :: sig, Qxxx, Qxrr, Rxxjj, Sxiijj

  ! Reconstruct state
  rho   = sol(1) + 1.0e-35
  rhoux = sol(2)
  ux    = rhoux/rho
  Pxx   = sol(3) - rhoux*ux
  Prr   = sol(4)
  qx    = sol(5) - (rhoux*ux*ux + ux*(3*Pxx + Prr))
  Riijj = sol(6) - (rhoux*ux*ux*ux + ux*ux*(6*Pxx + 2*Prr) + 4*ux*qx)

  ! Compute closing moments
  CALL MOM14_AXI_CLOSING_INTERPOLATIVE(rho, Pxx, Prr, qx, Riijj, sig, Qxxx, Qxrr, Rxxjj, Sxiijj)

  Flux(1) = rhoux
  Flux(2) = rhoux*ux + Pxx
  Flux(3) = rhoux*ux*ux + 3*ux*Pxx + Qxxx
  Flux(4) = ux*Prr + Qxrr
  Flux(5) = rhoux*ux*ux*ux + ux*ux*(6*Pxx + Prr) + 2*ux*(qx + Qxxx) + Rxxjj
  Flux(6) = rhoux*ux*ux*ux*ux + ux*ux*ux*(10*Pxx + 2*Prr) + ux*ux*(6*qx + 4*Qxxx) &
            + ux*(Riijj + 4*Rxxjj) + Sxiijj

  END SUBROUTINE MOM14_AXI_FLUXES

  ! ****************************************************************************************

  SUBROUTINE MOM14_AXI_SOURCES(sol, Ex, Ey, B, Src)

  ! Computes the source terms for this PDE. 
  ! Can be of different nature.

  IMPLICIT NONE
  
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8),                  INTENT(IN)  :: Ex, Ey, B
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(OUT) :: Src

  REAL(KIND=8), DIMENSION(N_EQ) :: Src_BGK, Src_EM

  CALL MOM14_AXI_BGK_SOURCES(sol, Src_BGK) 
  CALL MOM14_AXI_EM_SOURCES(sol, Ex, Ey, B, Src_EM)

  Src = Src_BGK + Src_EM

  END SUBROUTINE MOM14_AXI_SOURCES

  ! ****************************************************************************************

  SUBROUTINE MOM14_AXI_BGK_SOURCES(sol, BGK_Src)
  ! Computes BGK-style collision sources

  IMPLICIT NONE
  
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(OUT) :: BGK_Src
 
  ! REAL(KIND=8) :: Tau = 1.0e-7 ! HARD-CODED, mean time between collisions, 1/collis_freq
  REAL(KIND=8) :: Tau = 10000.0d0 ! HUGE! ! HARD-CODED, mean time between collisions, 1/collis_freq 
  REAL(KIND=8) :: rho, rhoux, ux, Pxx, Prr, qx, Riijj
  
  ! Reconstruct state
  rho   = sol(1) + 1.0e-35
  rhoux = sol(2)
  ux    = rhoux/rho
  Pxx   = sol(3) - rhoux*ux
  Prr   = sol(4)
  qx    = sol(5) - (rhoux*ux*ux + ux*(3*Pxx + Prr))
  Riijj = sol(6) - (rhoux*ux*ux*ux + ux*ux*(6*Pxx + 2*Prr) + 4*ux*qx)

  ! Build BGK sources vector
  BGK_Src(1) = 0.0d0
  BGK_Src(2) = 0.0d0
  BGK_Src(3) = (Prr - 2*Pxx)/(3.0d0*Tau + 1.0e-35)
  BGK_Src(4) = (2*Pxx - Prr)/(3.0d0*Tau + 1.0e-35)
  BGK_Src(5) = (2*ux*(Prr - 2*Pxx) - 3*qx)/(3.0d0*Tau)
  BGK_Src(6) = (4*ux**2*(Prr - 2*Pxx) - 12*ux*qx + 5*(Pxx+Prr)**2/rho - 3*Riijj)/(3.0d0*Tau)

  END SUBROUTINE MOM14_AXI_BGK_SOURCES

  ! ****************************************************************************************

  SUBROUTINE MOM14_AXI_EM_SOURCES(sol, Ex, Ey, B, EM_SRC)
  ! Computes EM sources
 
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8),                  INTENT(IN)  :: Ex, Ey, B
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(OUT) :: EM_SRC

  EM_SRC(1) = 0.0d0
  EM_SRC(2) =       GAS_Q/GAS_M*Ex*sol(1)
  EM_SRC(3) = 2.0d0*GAS_Q/GAS_M*Ex*sol(2)
  EM_SRC(4) = 0.0d0
  EM_SRC(5) =       GAS_Q/GAS_M*Ex*(3*sol(3) + sol(4))
  EM_SRC(6) = 4.0d0*GAS_Q/GAS_M*Ex*sol(5)

  END SUBROUTINE MOM14_AXI_EM_SOURCES

  ! ****************************************************************************************

  SUBROUTINE MOM14_AXI_CLOSING_INTERPOLATIVE(rho, Pxx, Prr, qx, Riijj, sig, Qxxx, Qxrr, Rxxii, Sxiijj)
  ! Computes interpolative closing moments for the axisymmetric 14-mom closure

  IMPLICIT NONE
  
  REAL(KIND=8), INTENT(IN)  :: rho, Pxx, Prr, qx, Riijj
  REAL(KIND=8), INTENT(OUT) :: sig, Qxxx, Qxrr, Rxxii, Sxiijj
  
  ! Working variables
  REAL(KIND=8) :: A, B_tmp, POW_SIG, SIG_LOW_LIMIT

  ! Compute sigma
  B_TMP = 2*Pxx**2 + Prr**2 + (Pxx + Prr)**2 - rho*Riijj 
  
  SIG = (B_TMP + SQRT(B_TMP**2 + 4*rho*(2*Pxx**2 + Prr**2)*qx**2/Pxx))/(2.0d0*(2.0d0*Pxx**2 + Prr**2))

  ! Artificial limit to sigma
  SIG_LOW_LIMIT = 2.0d-4

  IF (SIG .LT. SIG_LOW_LIMIT) THEN
    sigma_small_flag = '*'
    SIG = MAX(SIG, SIG_LOW_LIMIT) ! Artificial limitation on sigma
  END IF

  ! Compute closing moments
  POW_SIG = 1.0d0/2.0d0
  A      = (6*Pxx**2)/(Prr**2 + 6*Pxx**2)

  Qxxx   = A*qx
  Qxrr   = (1-A)*qx

  Rxxii  = (A/SIG)*qx**2/Pxx + (2*(1 - SIG)*Pxx**2 + (Pxx + Prr)*Pxx)/rho
  
  Sxiijj = A/(SIG**2)*qx**3/(Pxx**2)    & 
           + 2.0d0/rho*(Pxx + Prr + (1-SIG**POW_SIG)*(Prr**3 + 2*(Prr**2)*Pxx + 24*Pxx**3) &
                                    /(Prr**2 + 6*Pxx**2))*qx
 
  END SUBROUTINE MOM14_AXI_CLOSING_INTERPOLATIVE
 
  ! *****************************************************************************************

  SUBROUTINE MOM14_AXI_MIN_MAX_WAVESPEEDS(sol, sMin, sMax)

  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8),                  INTENT(OUT) :: sMin, sMax

  REAL(KIND=8) :: rho, rhoux, ux, Pxx, ax, MULT_CONST

  REAL(KIND=8), DIMENSION(N_EQ,N_EQ) :: JMAT
  REAL(KIND=8), DIMENSION(N_EQ) :: sol_nondim

  INTEGER :: INFO
  REAL(KIND=8), DIMENSION(N_EQ) :: WR, WI, W
  
  INTEGER :: II
 
  ! Unpack the solution and compute primitive variables
  rho    = sol(1)
  rhoux  = sol(2)
  ux     = rhoux/(rho + 1.0d-35)
  Pxx    = sol(3)  - rhoux*ux
  ax     = SQRT(GAS_GAMMA*Pxx/(rho+1.0d-35)) ! Speed of sound along x - IS IT CORRECT TO USE Pxx??

  !!!  ! We arbitrarily scale the speed of sound, by a constant.
  !!!  ! Using something in the order of 1 - 10 should do, but notice that it can really 
  !!!  ! be much different, depending on the problem! 
  !!!  ! This will introduce unphysical dissipation anyway, so you should at least try different 
  !!!  ! values.
  !!!  MULT_CONST = 5.0d0 
  !!!  
  !!!  sMin = ux - ax*MULT_CONST
  !!!  sMax = ux + ax*MULT_CONST
 
  ! Compute Jacobian and its maximum and minimum wavespeeds
  ! NONDIMENSIONALIZE SOLUTION (using Pxx)
  sol_nondim = sol/rho

  sol_nondim(1) = sol_nondim(1) ! ok
  sol_nondim(2) = sol_nondim(2)/(SQRT(Pxx/rho)) 
  sol_nondim(3) = sol_nondim(3)/(SQRT(Pxx/rho)**2)  ! Pxx
  sol_nondim(4) = sol_nondim(4)/(SQRT(Pxx/rho)**2)  ! Prr
  sol_nondim(5) = sol_nondim(5)/(SQRT(Pxx/rho)**3)  ! qx
  sol_nondim(6) = sol_nondim(6)/(SQRT(Pxx/rho)**4)  ! Riijj

  CALL MOM14_AXI_COMPUTE_JACOBIAN_FD(sol_nondim, JMAT)

  CALL COMPUTE_EIGS_LAPACK_DGEEV(JMAT, N_EQ, WR, WI, INFO)

  ! The 3D moment closure is not strictly hyperbolic.  
  ! Take the modulus of real and imag parts, and use the sign of the real part
  ! W = SIGN(Re)*SQRT( Re**2 + Im**2 )
  DO II = 1,N_EQ
    W(II) = DSIGN( SQRT(WR(II)**2 + WI(II)**2)  , WR(II))
  END DO

  sMin = MINVAL(W)
  sMax = MAXVAL(W)

  sMin = sMin*SQRT(Pxx/rho) + ux
  sMax = sMax*SQRT(Pxx/rho) + ux

  IF (INFO .GT. 0) THEN ! Algorithm failed
    PRINT*, "Eigenvalues of Jacobian could not be computed. Using wavespeed estimates instead!"

    MULT_CONST = 5.0d0
    sMin = ux - ax*MULT_CONST
    sMax = ux + ax*MULT_CONST

  END IF

  END SUBROUTINE MOM14_AXI_MIN_MAX_WAVESPEEDS

  ! ******************************************************************************************

  SUBROUTINE MOM14_AXI_COMPUTE_JACOBIAN(sol, JMAT)

  ! Computes the Jacobian from the solution, using a bit of analytical derivation and finite
  ! differences for the rest.
  ! This is done in a frame of reference where the velocities are zero.
  !
  ! First, primitive variables are reconstructed from the solution. Then, assuming that ux is zero
  ! the jacobian is constructed.

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ),       INTENT(IN)  :: sol
  REAL(KIND=8), DIMENSION(N_EQ, N_EQ), INTENT(OUT) :: JMAT

  REAL(KIND=8) :: rho, rhoux, ux, Pxx, Prr, qx, Riijj

  REAL(KIND=8) :: A, dummy, PERT, TOL
  REAL(KIND=8) :: dQxxxdPxx, dQxxxdPrr, dQxxxdqx
  REAL(KIND=8) :: dQxrrdPxx, dQxrrdPrr, dQxrrdqx
  REAL(KIND=8) :: dRxxjjdrho,  dRxxjjdPxx,  dRxxjjdPrr,  dRxxjjdqx,  dRxxjjdRiijj
  REAL(KIND=8) :: dSxiijjdrho,  dSxiijjdPxx,  dSxiijjdPrr,  dSxiijjdqx,  dSxiijjdRiijj

  REAL(KIND=8) :: Rxxjj, Rxxjj_pert, Sxiijj, Sxiijj_pert
  REAL(KIND=8) :: rho_pert, Pxx_pert, Prr_pert, qx_pert, Riijj_pert

  ! Reconstruct state
  rho   = sol(1) + 1.0e-35
  rhoux = sol(2)
  ux    = rhoux/rho
  Pxx   = sol(3) - rhoux*ux
  Prr   = sol(4)
  qx    = sol(5) - (rhoux*ux*ux + ux*(3*Pxx + Prr))
  Riijj = sol(6) - (rhoux*ux*ux*ux + ux*ux*(6*Pxx + 2*Prr) + 4*ux*qx)

  ! Compute derivatives of Q

  A = 6*Pxx**2/(Prr**2 + 6*Pxx**2)

  dQxxxdPxx = 12*qx*Pxx*Prr**2/((Prr**2 + 6*Pxx**2)**2)
  dQxxxdPrr = -12*qx*Pxx**2*Prr/((Prr**2 + 6*Pxx**2)**2)
  dQxxxdqx  = A

  dQxrrdPxx = -12*qx*Pxx*Prr**2/((Prr**2 + 6*Pxx**2)**2)
  dQxrrdPrr = 12*qx*Pxx**2*Prr/((Prr**2 + 6*Pxx**2)**2)
  dQxrrdqx  = 1-A

  ! --- Compute other derivatives by perturbing the state

  ! Unperturbed moments
  CALL MOM14_AXI_CLOSING_INTERPOLATIVE(rho, Pxx, Prr, qx, Riijj, dummy, dummy, dummy, Rxxjj, Sxiijj)

  PERT = 1.0d-6
  TOL  = 1.0d-15

  ! Perturb density <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  rho_pert = rho + PERT*MAX(ABS(rho), TOL)
  CALL MOM14_AXI_CLOSING_INTERPOLATIVE(rho_pert, Pxx, Prr, qx, Riijj, dummy, dummy, dummy, Rxxjj_pert, Sxiijj_pert)

  dRxxjjdrho  = (Rxxjj_pert  - Rxxjj)/(rho_pert - rho)
  dSxiijjdrho = (Sxiijj_pert - Sxiijj)/(rho_pert - rho)
  
  ! Perturb Pxx <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  Pxx_pert = Pxx + PERT*MAX(ABS(Pxx), TOL)
  CALL MOM14_AXI_CLOSING_INTERPOLATIVE(rho, Pxx_pert, Prr, qx, Riijj, dummy, dummy, dummy, Rxxjj_pert, Sxiijj_pert)

  dRxxjjdPxx  = (Rxxjj_pert  - Rxxjj)/(Pxx_pert - Pxx)
  dSxiijjdPxx = (Sxiijj_pert - Sxiijj)/(Pxx_pert - Pxx)
 
  ! Perturb Prr <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  Prr_pert = Prr + PERT*MAX(ABS(Prr), TOL)
  CALL MOM14_AXI_CLOSING_INTERPOLATIVE(rho, Pxx, Prr_pert, qx, Riijj, dummy, dummy, dummy, Rxxjj_pert, Sxiijj_pert)

  dRxxjjdPrr  = (Rxxjj_pert  - Rxxjj)/(Prr_pert - Prr)
  dSxiijjdPrr = (Sxiijj_pert - Sxiijj)/(Prr_pert - Prr)
 
  ! Perturb qx <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  qx_pert = qx + PERT*MAX(ABS(qx), TOL)
  CALL MOM14_AXI_CLOSING_INTERPOLATIVE(rho, Pxx, Prr, qx_pert, Riijj, dummy, dummy, dummy, Rxxjj_pert, Sxiijj_pert)

  dRxxjjdqx  = (Rxxjj_pert  - Rxxjj)/(qx_pert - qx)
  dSxiijjdqx = (Sxiijj_pert - Sxiijj)/(qx_pert - qx)
 
  ! Perturb Riijj <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  Riijj_pert = Riijj + PERT*MAX(ABS(Riijj), TOL)
  CALL MOM14_AXI_CLOSING_INTERPOLATIVE(rho, Pxx, Prr, qx, Riijj_pert, dummy, dummy, dummy, Rxxjj_pert, Sxiijj_pert)

  dRxxjjdRiijj  = (Rxxjj_pert  - Rxxjj)/(Riijj_pert - Riijj)
  dSxiijjdRiijj = (Sxiijj_pert - Sxiijj)/(Riijj_pert - Riijj)

  ! ++++++++ Compose the Jacobian ++++++++
  JMAT = 0.0d0

  JMAT(2, 3) = 1.0d0

  JMAT(3,3) = dQxxxdPxx
  JMAT(3,4) = dQxxxdPrr
  JMAT(3,5) = dQxxxdqx

  JMAT(4,3) = dQxrrdPxx
  JMAT(4,4) = dQxrrdPrr
  JMAT(4,5) = dQxrrdqx

  JMAT(5,1) = dRxxjjdrho
  JMAT(5,3) = dRxxjjdPxx
  JMAT(5,4) = dRxxjjdPrr
  JMAT(5,5) = dRxxjjdqx
  JMAT(5,6) = dRxxjjdRiijj
  
  JMAT(6,1) = dSxiijjdrho
  JMAT(6,3) = dSxiijjdPxx
  JMAT(6,4) = dSxiijjdPrr
  JMAT(6,5) = dSxiijjdqx
  JMAT(6,6) = dSxiijjdRiijj

!   PRINT*, "DBDBDBD >>>>>>>>>>>>>>>>>>>>>>>>>>> JACOBIAN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
!   PRINT*, JMAT
!   PRINT*, "DBDBDBD >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

  END SUBROUTINE MOM14_AXI_COMPUTE_JACOBIAN


  ! ******************************************************************************************

  SUBROUTINE MOM14_AXI_COMPUTE_JACOBIAN_FD(sol, JMAT)

  ! Computes the Jacobian from the solution, using finite differences.

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ),       INTENT(IN)  :: sol
  REAL(KIND=8), DIMENSION(N_EQ, N_EQ), INTENT(OUT) :: JMAT

  REAL(KIND=8), DIMENSION(N_EQ) :: sol_p ! Perturbed solution
  REAL(KIND=8), DIMENSION(N_EQ) :: flux, flux_p ! Perturbed solution

  INTEGER :: I_EQ
  REAL(KIND=8) :: PERT, ETA, TOL

  CALL MOM14_AXI_FLUXES(sol, flux) ! Compute flux from unperturbed sol

  ! Perturb solution
  ETA  = 1.0d-6
  TOL  = 1.0d-15
 
  DO I_EQ = 1, N_EQ ! For every entry in the solution vector

    sol_p = sol ! Reset

    PERT  =  ETA*DSIGN(MAX(ABS(sol(I_EQ)), TOL), sol(I_EQ)) ! Perturbation
    sol_p(I_EQ) = sol(I_EQ) + PERT        ! Perturbed state

    CALL MOM14_AXI_FLUXES(sol_p, flux_p) ! Perturbed flux

    JMAT(:, I_EQ) = (flux_p - flux)/PERT

  END DO

  !!!  PRINT*, "DBDBDBD >>>>>>>>>>>>>>>>>>>>>>>>>>> JACOBIAN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
  !!!  PRINT*, JMAT
  !!!  PRINT*, "DBDBDBD >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

  END SUBROUTINE MOM14_AXI_COMPUTE_JACOBIAN_FD

END MODULE PDE_MOM14_AXI
