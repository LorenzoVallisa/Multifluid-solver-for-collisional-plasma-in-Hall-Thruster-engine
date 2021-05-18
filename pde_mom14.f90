MODULE PDE_MOM14

  !################################################################################
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

  SUBROUTINE MOM14_INIT_SOL(sol)
  ! Initializes the solution vector for the 14-moments PDEs. 
  ! HARD-CODED FOR NOW!

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ, N_cells), INTENT(INOUT) :: sol

  REAL(KIND=8) :: P14, rho14, ux14, uy14, uz14, u2_14
  REAL(KIND=8) :: Pxx14, Pxy14, Pxz14, Pyy14, Pyz14, Pzz14
  REAL(KIND=8) :: qx14, qy14, qz14, Riijj14

  ! ============= Initialize solution (for 14 moments) ===============
!!! UNIFORM STATE !!!  ! Uniform state
!!! UNIFORM STATE !!!
!!! UNIFORM STATE !!!  rho14   = 9.1094E-11
!!! UNIFORM STATE !!!  P14     = 138.0d0
!!! UNIFORM STATE !!!
!!! UNIFORM STATE !!!  ux14    = 0.0d0
!!! UNIFORM STATE !!!  uy14    = 0.0d0
!!! UNIFORM STATE !!!  uz14    = 0.0d0
!!! UNIFORM STATE !!!  Pxx14   = P14
!!! UNIFORM STATE !!!  Pxy14   = 0.0d0
!!! UNIFORM STATE !!!  Pxz14   = 0.0d0
!!! UNIFORM STATE !!!  Pyy14   = P14
!!! UNIFORM STATE !!!  Pyz14   = 0.0d0
!!! UNIFORM STATE !!!  Pzz14   = P14
!!! UNIFORM STATE !!!  qx14    = 0.0d0    
!!! UNIFORM STATE !!!  qy14    = 0.0d0    
!!! UNIFORM STATE !!!  qz14    = 0.0d0    
!!! UNIFORM STATE !!!  Riijj14 = 15.0d0*(P14**2)/rho14
!!! UNIFORM STATE !!!
!!! UNIFORM STATE !!!  u2_14 = SQRT(ux14**2 + uy14**2 + uz14**2)
!!! UNIFORM STATE !!!
!!! UNIFORM STATE !!!  ! Conserved variables, left state
!!! UNIFORM STATE !!!  sol(1,:) = rho14   ! density
!!! UNIFORM STATE !!!
!!! UNIFORM STATE !!!  sol(2,:) = rho14 * ux14   ! momentum x
!!! UNIFORM STATE !!!  sol(3,:) = rho14 * uy14   ! momentum y
!!! UNIFORM STATE !!!  sol(4,:) = rho14 * uz14   ! momentum z
!!! UNIFORM STATE !!!
!!! UNIFORM STATE !!!  sol(5,:)  = rho14*ux14*ux14 + Pxx14
!!! UNIFORM STATE !!!  sol(6,:)  = rho14*ux14*uy14 + Pxy14
!!! UNIFORM STATE !!!  sol(7,:)  = rho14*ux14*uz14 + Pxz14
!!! UNIFORM STATE !!!  sol(8,:)  = rho14*uy14*uy14 + Pyy14
!!! UNIFORM STATE !!!  sol(9,:)  = rho14*uy14*uz14 + Pyz14
!!! UNIFORM STATE !!!  sol(10,:) = rho14*uz14*uz14 + Pzz14 
!!! UNIFORM STATE !!!
!!! UNIFORM STATE !!!  sol(11,:) = rho14*ux14*u2_14 + 3*ux14*P14 + 2*(ux14*Pxx14 + uy14*Pxy14 + uz14*Pxz14) + qx14
!!! UNIFORM STATE !!!  sol(12,:) = rho14*uy14*u2_14 + 3*uy14*P14 + 2*(ux14*Pxy14 + uy14*Pyy14 + uz14*Pyz14) + qy14
!!! UNIFORM STATE !!!  sol(13,:) = rho14*uz14*u2_14 + 3*uz14*P14 + 2*(ux14*Pxz14 + uy14*Pyz14 + uz14*Pzz14) + qz14
!!! UNIFORM STATE !!!
!!! UNIFORM STATE !!!  sol(14,:) = rho14*u2_14*u2_14 + 2*u2_14*3*P14 &
!!! UNIFORM STATE !!!          + 4*( ux14*ux14*Pxx14 + uy14*uy14*Pyy14 + uz14*uz14*Pzz14 &
!!! UNIFORM STATE !!!                + 2*( ux14*uy14*Pxy14 + uy14*uz14*Pyz14 + ux14*uz14*Pxz14 ) ) &
!!! UNIFORM STATE !!!          + 4*( ux14*qx14 + uy14*qy14 + uz14*qz14 ) + Riijj14
!!! UNIFORM STATE !!!

  ! Left and right states
  P14 = 404400.0d0

  rho14   = 4.696d0
  ux14    = 5000.0d0
  uy14    = 0.0d0
  uz14    = 0.0d0
  Pxx14   = P14
  Pxy14   = 0.0d0
  Pxz14   = 0.0d0
  Pyy14   = P14
  Pyz14   = 0.0d0
  Pzz14   = P14
  qx14    = 0.0d0    
  qy14    = 0.0d0    
  qz14    = 0.0d0    
  Riijj14 = 15.0d0*(P14**2)/rho14

  u2_14 = ux14**2 + uy14**2 + uz14**2

  ! Conserved variables, left state
  sol(1,:) = rho14   ! density

  sol(2,:) = rho14 * ux14   ! momentum x
  sol(3,:) = rho14 * uy14   ! momentum y
  sol(4,:) = rho14 * uz14   ! momentum z

  sol(5,:)  = rho14*ux14*ux14 + Pxx14
  sol(6,:)  = rho14*ux14*uy14 + Pxy14
  sol(7,:)  = rho14*ux14*uz14 + Pxz14
  sol(8,:)  = rho14*uy14*uy14 + Pyy14
  sol(9,:)  = rho14*uy14*uz14 + Pyz14
  sol(10,:) = rho14*uz14*uz14 + Pzz14 

  sol(11,:) = rho14*ux14*u2_14 + 3*ux14*P14 + 2*(ux14*Pxx14 + uy14*Pxy14 + uz14*Pxz14) + qx14
  sol(12,:) = rho14*uy14*u2_14 + 3*uy14*P14 + 2*(ux14*Pxy14 + uy14*Pyy14 + uz14*Pyz14) + qy14
  sol(13,:) = rho14*uz14*u2_14 + 3*uz14*P14 + 2*(ux14*Pxz14 + uy14*Pyz14 + uz14*Pzz14) + qz14

  sol(14,:) = rho14*u2_14*u2_14 + 2*u2_14*3*P14 &
            + 4*( ux14*ux14*Pxx14 + uy14*uy14*Pyy14 + uz14*uz14*Pzz14 &
                  + 2*( ux14*uy14*Pxy14 + uy14*uz14*Pyz14 + ux14*uz14*Pxz14 ) ) &
            + 4*( ux14*qx14 + uy14*qy14 + uz14*qz14 ) + Riijj14

  ! Right state
  ! P14 = 101100.0d0

  ! rho14   = 1.408d0
  ux14    = -ux14
  ! uy14    = 0.0d0
  ! uz14    = 0.0d0
  ! Pxx14   = P14
  ! Pxy14   = 0.0d0
  ! Pxz14   = 0.0d0
  ! Pyy14   = P14
  ! Pyz14   = 0.0d0
  ! Pzz14   = P14
  ! qx14    = 0.0d0    
  ! qy14    = 0.0d0    
  ! qz14    = 0.0d0    
  ! Riijj14 = 15.0d0*(P14**2)/rho14

  u2_14 = ux14**2 + uy14**2 + uz14**2

  ! Conserved variables, left state
  sol(1, FLOOR(N_cells/2.0d0):) = rho14   ! density

  sol(2, FLOOR(N_cells/2.0d0):) = rho14 * ux14   ! momentum x
  sol(3, FLOOR(N_cells/2.0d0):) = rho14 * uy14   ! momentum y
  sol(4, FLOOR(N_cells/2.0d0):) = rho14 * uz14   ! momentum z

  sol(5, FLOOR(N_cells/2.0d0):)  = rho14*ux14*ux14 + Pxx14
  sol(6, FLOOR(N_cells/2.0d0):)  = rho14*ux14*uy14 + Pxy14
  sol(7, FLOOR(N_cells/2.0d0):)  = rho14*ux14*uz14 + Pxz14
  sol(8, FLOOR(N_cells/2.0d0):)  = rho14*uy14*uy14 + Pyy14
  sol(9, FLOOR(N_cells/2.0d0):)  = rho14*uy14*uz14 + Pyz14
  sol(10,FLOOR(N_cells/2.0d0):) = rho14*uz14*uz14 + Pzz14 

  sol(11,FLOOR(N_cells/2.0d0):) = rho14*ux14*u2_14 + 3*ux14*P14 + 2*(ux14*Pxx14 + uy14*Pxy14 + uz14*Pxz14) + qx14
  sol(12,FLOOR(N_cells/2.0d0):) = rho14*uy14*u2_14 + 3*uy14*P14 + 2*(ux14*Pxy14 + uy14*Pyy14 + uz14*Pyz14) + qy14
  sol(13,FLOOR(N_cells/2.0d0):) = rho14*uz14*u2_14 + 3*uz14*P14 + 2*(ux14*Pxz14 + uy14*Pyz14 + uz14*Pzz14) + qz14

  sol(14,FLOOR(N_cells/2.0d0):) = rho14*u2_14*u2_14 + 2*u2_14*3*P14 &
            + 4*( ux14*ux14*Pxx14 + uy14*uy14*Pyy14 + uz14*uz14*Pzz14 &
                  + 2*( ux14*uy14*Pxy14 + uy14*uz14*Pyz14 + ux14*uz14*Pxz14 ) ) &
            + 4*( ux14*qx14 + uy14*qy14 + uz14*qz14 ) + Riijj14


  END SUBROUTINE MOM14_INIT_SOL


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  SUBROUTINE MOM14_FLUXES(sol, Flux)

  ! Computes flux from a solution vector, of size N_EQ

  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(OUT) :: Flux

  REAL(KIND=8) :: rho, rhoux, rhouy, rhouz, ux, uy, uz, u2
  REAL(KIND=8) :: Pxx, Pxy, Pxz, Pyy, Pyz, Pzz, qx, qy, qz, Riijj
  REAL(KIND=8), DIMENSION(3,3) :: Pij_MAT
  REAL(KIND=8), DIMENSION(3)   :: uVEC, qVEC 

  ! Closing moments
  REAL(KIND=8), DIMENSION(10) :: Qijk_VECT  
  REAL(KIND=8), DIMENSION(6)  :: Rijkk_VECT 
  REAL(KIND=8), DIMENSION(3)  :: Sijjkk_VECT
  REAL(KIND=8) :: Qxxx, Qxxy, Qxyy, Qyyy, Qyyz, Qyzz, Qzzz, Qxxz, Qxzz, Qxyz
  REAL(KIND=8) :: Rxxkk, Rxykk, Ryykk, Ryzkk, Rzzkk, Rxzkk, Sxjjkk, Syjjkk, Szjjkk

  REAL(KIND=8) :: sig
  REAL(KIND=8), DIMENSION(3,3) :: Qij_x

  ! Reconstruct state
  rho   = sol(1)

  rhoux = sol(2)
  rhouy = sol(3)
  rhouz = sol(4)
  ux = rhoux/(rho + 1.0e-35)
  uy = rhouy/(rho + 1.0e-35)
  uz = rhouz/(rho + 1.0e-35)

  Pxx = sol(5)  - rhoux*ux
  Pxy = sol(6)  - rhoux*uy
  Pxz = sol(7)  - rhoux*uz
  Pyy = sol(8)  - rhouy*uy
  Pyz = sol(9)  - rhouy*uz
  Pzz = sol(10) - rhouz*uz

  u2 = ux**2 + uy**2 + uz**2

  uVEC(1) = ux
  uVEC(2) = uy
  uVEC(3) = uz

  Pij_MAT(1,1) = Pxx
  Pij_MAT(1,2) = Pxy
  Pij_MAT(1,3) = Pxz

  Pij_MAT(2,1) = Pxy
  Pij_MAT(2,2) = Pyy
  Pij_MAT(2,3) = Pyz

  Pij_MAT(3,1) = Pxz
  Pij_MAT(3,2) = Pyz
  Pij_MAT(3,3) = Pzz

  qx = sol(11) - ( rhoux*u2 + (Pxx+Pyy+Pzz)*ux + 2*(ux*Pxx + uy*Pxy + uz*Pxz) )
  qy = sol(12) - ( rhouy*u2 + (Pxx+Pyy+Pzz)*uy + 2*(ux*Pxy + uy*Pyy + uz*Pyz) )
  qz = sol(13) - ( rhouz*u2 + (Pxx+Pyy+Pzz)*uz + 2*(ux*Pxz + uy*Pyz + uz*Pzz) )

  qVEC(1) = qx
  qVEC(2) = qy
  qVEC(3) = qz

  Riijj = sol(14) - (rho*u2*u2                                   & 
                     + 2*u2*(Pxx+Pyy+Pzz)                        &
                     + 4*DOT_PRODUCT(uVEC, MATMUL(Pij_MAT,uVEC)) &
                     + 4*DOT_PRODUCT(uVEC, qVEC))

  ! Compute closing moments
  CALL MOM14_CLOSING_INTERPOLATIVE(rho, Pij_MAT, qVEC, Riijj, sig, Qijk_VECT, Rijkk_VECT, Sijjkk_VECT)

  sig_min = MIN(sig_min, sig)
  sig_max = MAX(sig_max, sig)

  Qxxx = Qijk_VECT(1) 
  Qxxy = Qijk_VECT(2) 
  Qxyy = Qijk_VECT(3) 
  Qyyy = Qijk_VECT(4) 
  Qyyz = Qijk_VECT(5) 
  Qyzz = Qijk_VECT(6) 
  Qzzz = Qijk_VECT(7) 
  Qxxz = Qijk_VECT(8) 
  Qxzz = Qijk_VECT(9) 
  Qxyz = Qijk_VECT(10)
    
  Rxxkk = Rijkk_VECT(1) 
  Rxykk = Rijkk_VECT(2) 
  Rxzkk = Rijkk_VECT(3) 
  Ryykk = Rijkk_VECT(4) 
  Ryzkk = Rijkk_VECT(5) 
  Rzzkk = Rijkk_VECT(6)
  
  Sxjjkk = Sijjkk_VECT(1) 
  Syjjkk = Sijjkk_VECT(2) 
  Szjjkk = Sijjkk_VECT(3)


!   PRINT*,   Qxxx 
!   PRINT*,   Qxxy 
!   PRINT*,   Qxyy 
!   PRINT*,   Qyyy 
!   PRINT*,   Qyyz 
!   PRINT*,   Qyzz 
!   PRINT*,   Qzzz 
!   PRINT*,   Qxxz 
!   PRINT*,   Qxzz 
!   PRINT*,   Qxyz 
!   PRINT*
!   PRINT*,   Rxxkk 
!   PRINT*,   Rxykk 
!   PRINT*,   Rxzkk 
!   PRINT*,   Ryykk 
!   PRINT*,   Ryzkk 
!   PRINT*,   Rzzkk 
!   PRINT*
!   PRINT*,   "Riijj: ", Riijj
!   PRINT*
!   PRINT*,   Sxjjkk
!   PRINT*,   Syjjkk
!   PRINT*,   Szjjkk
!   PRINT*
! 
!   WRITE(*,*) "DBDBDBDBDB SIGMA SHOULD BE AROUND 0: ", sig

  ! Compute fluxes
  Flux(1) = sol(2) ! rho ux
 
  Flux(2) = sol(5) ! rho ux ux + Pxx
  Flux(3) = sol(6) ! rho uy ux + Pyx
  Flux(4) = sol(7) ! rho uz ux + Pzx

  Flux(5)  = rhoux*ux*ux + ux*Pxx + ux*Pxx + ux*Pxx + Qxxx  
  Flux(6)  = rhoux*uy*ux + ux*Pxy + uy*Pxx + ux*Pxy + Qxxy  
  Flux(7)  = rhoux*uz*ux + ux*Pxz + uz*Pxx + ux*Pxz + Qxxz  
  Flux(8)  = rhouy*uy*ux + uy*Pxy + uy*Pxy + ux*Pyy + Qxyy  
  Flux(9)  = rhouy*uz*ux + uy*Pxz + uz*Pxy + ux*Pyz + Qxyz  
  Flux(10) = rhouz*uz*ux + uz*Pxz + uz*Pxz + ux*Pzz + Qxzz  

  Flux(11) = rhoux*ux*u2 + ux*ux*(Pxx+Pyy+Pzz) + 2*ux*(ux*Pxx + uy*Pxy + uz*Pxz) &
           + 2*ux*(Pxx*ux + Pxy*uy + Pxz*uz) + u2*Pxx + ux*qx + ux*qx            &
           + 2*(ux*Qxxx + uy*Qxxy + uz*Qxxz) + Rxxkk

  Flux(12) = rhouy*ux*u2 + uy*ux*(Pxx+Pyy+Pzz) + 2*uy*(ux*Pxx + uy*Pxy + uz*Pxz) &
           + 2*ux*(Pxy*ux + Pyy*uy + Pyz*uz) + u2*Pxy + uy*qx + ux*qy            & 
           + 2*(ux*Qxxy + uy*Qxyy + uz*Qxyz) + Rxykk

  Flux(13) = rhouz*ux*u2 + uz*ux*(Pxx+Pyy+Pzz) + 2*uz*(ux*Pxx + uy*Pxy + uz*Pxz) &
           + 2*ux*(Pxz*ux + Pyz*uy + Pzz*uz) + u2*Pxz + uz*qx + ux*qz            &
           + 2*(ux*Qxxz + uy*Qxyz + uz*Qxzz) + Rxzkk

  Qij_x(1,:) = (/ Qxxx, Qxxy, Qxxz /) 
  Qij_x(2,:) = (/ Qxxy, Qxyy, Qxyz /) 
  Qij_x(3,:) = (/ Qxxz, Qxyz, Qxzz /) 

  Flux(14) = rhoux*u2*u2 + 2*ux*u2*(Pxx+Pyy+Pzz) + 4*u2*(ux*Pxx + uy*Pxy + uz*Pxz)    &
           + 4*DOT_PRODUCT(uVEC, MATMUL(Pij_MAT,uVEC))*ux + 2*u2*qx                   &
           + 4*ux*(ux*qx+uy*qy+uz*qz) + 4*DOT_PRODUCT(uVEC, MATMUL(Qij_x, uVEC))      &
           + 4*(ux*Rxxkk + uy*Rxykk + uz*Rxzkk) + ux*(Rxxkk + Ryykk + Rzzkk) + Sxjjkk

  END SUBROUTINE MOM14_FLUXES

  ! ****************************************************************************************

  SUBROUTINE MOM14_SOURCES(sol, Ex, Ey, B, Src)

  ! Computes the source terms for this PDE. 
  ! Can be of different nature.

  IMPLICIT NONE
  
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8),                  INTENT(IN)  :: Ex, Ey, B
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(OUT) :: Src

  REAL(KIND=8), DIMENSION(N_EQ) :: Src_BGK, Src_EM

  CALL MOM14_BGK_SOURCES(sol, Src_BGK) 
  CALL MOM14_EM_SOURCES(sol, Ex, Ey, B, Src_EM)

  Src = Src_BGK + Src_EM

  END SUBROUTINE MOM14_SOURCES

  ! ****************************************************************************************

  SUBROUTINE MOM14_BGK_SOURCES(sol, BGK_Src)
  ! Computes BGK-style collision sources

  IMPLICIT NONE
  
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(OUT) :: BGK_Src
 
  ! REAL(KIND=8) :: Tau = 1.0e-11 ! HARD-CODED, mean time between collisions, 1/collis_freq
  REAL(KIND=8) :: Tau = 10000.0d0 ! HUGE! ! HARD-CODED, mean time between collisions, 1/collis_freq 

  REAL(KIND=8) :: rho,rhoux,rhouy,rhouz,ux,uy,uz,u2,Pxx,Pyy,Pzz, P_BGK

  ! Reconstruct state
  rho   = sol(1)

  rhoux = sol(2)
  rhouy = sol(3)
  rhouz = sol(4)
  ux = rhoux/(rho + 1.0e-35)
  uy = rhouy/(rho + 1.0e-35)
  uz = rhouz/(rho + 1.0e-35)

  u2 = ux**2 + uy**2 + uz**2

  Pxx = sol(5)  - rhoux*ux
  Pyy = sol(8)  - rhouy*uy
  Pzz = sol(10) - rhouz*uz

  P_BGK = (Pxx+Pyy+Pzz)/3.0d0

  ! Assemble BGK source
  BGK_Src(1) = 0.0d0 ! Density
  
  BGK_Src(2) = 0.0d0 ! rho ux
  BGK_Src(3) = 0.0d0 ! rho uy
  BGK_Src(4) = 0.0d0 ! rho uz

  BGK_Src(5)  = -(sol(5)  - rhoux*ux - P_BGK)/Tau
  BGK_Src(6)  = -(sol(6)  - rhoux*uy        )/Tau
  BGK_Src(7)  = -(sol(7)  - rhoux*uz        )/Tau
  BGK_Src(8)  = -(sol(8)  - rhouy*uy - P_BGK)/Tau
  BGK_Src(9)  = -(sol(9)  - rhouy*uz        )/Tau
  BGK_Src(10) = -(sol(10) - rhouz*uz - P_BGK)/Tau

  BGK_Src(11) = -(sol(11) - rhoux*u2 - 5.0d0*ux*P_BGK)/Tau
  BGK_Src(12) = -(sol(12) - rhouy*u2 - 5.0d0*uy*P_BGK)/Tau
  BGK_Src(13) = -(sol(13) - rhouz*u2 - 5.0d0*uz*P_BGK)/Tau

  BGK_Src(14) = -(sol(14) - rho*u2**2 - u2*P_BGK - 15.0d0*P_BGK**2/(rho + 1.0e-35))/Tau

  END SUBROUTINE MOM14_BGK_SOURCES

  ! ****************************************************************************************

  SUBROUTINE MOM14_EM_SOURCES(sol, Ex, Ey, B, EM_SRC)
  ! Computes EM sources
 
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8),                  INTENT(IN)  :: Ex, Ey, B
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(OUT) :: EM_SRC

  EM_SRC = 0.0d0

! PRINT*, "ATTENTION!!! CHECK THE EQUATIONS - DERIVE THEM ONCE MORE TO BE SURE!!!"
! PRINT*, "ATTENTION!!! CHECK THE EQUATIONS - DERIVE THEM ONCE MORE TO BE SURE!!!"
! PRINT*, "ATTENTION!!! CHECK THE EQUATIONS - DERIVE THEM ONCE MORE TO BE SURE!!!"
! PRINT*, "ATTENTION!!! CHECK THE EQUATIONS - DERIVE THEM ONCE MORE TO BE SURE!!!"
! PRINT*, "ATTENTION!!! CHECK THE EQUATIONS - DERIVE THEM ONCE MORE TO BE SURE!!!"
! PRINT*, "ATTENTION!!! CHECK THE EQUATIONS - DERIVE THEM ONCE MORE TO BE SURE!!!"
! PRINT*, "ATTENTION!!! CHECK THE EQUATIONS - DERIVE THEM ONCE MORE TO BE SURE!!!"
! PRINT*, "ATTENTION!!! CHECK THE EQUATIONS - DERIVE THEM ONCE MORE TO BE SURE!!!"

  EM_SRC(1)  = GAS_Q/GAS_M*( 0.0d0                                                      )    ! Mass equation

  EM_SRC(2)  = GAS_Q/GAS_M*( Ex*sol(1) + B*sol(3)                                       )    ! vx
  EM_SRC(3)  = GAS_Q/GAS_M*( Ey*sol(1) - B*sol(2)                                       )    ! vy
  EM_SRC(4)  = GAS_Q/GAS_M*( 0.0d0                                                      )    ! vz

  EM_SRC(5)  = GAS_Q/GAS_M*( 2*Ex*sol(2) + 2*B*sol(6)                                   )    ! vx vx
  EM_SRC(6)  = GAS_Q/GAS_M*( Ex*sol(3) + Ey*sol(2) + B*(sol(8) - sol(5))                )    ! vx vy
  EM_SRC(7)  = GAS_Q/GAS_M*( Ex*sol(4) + B*sol(9)                                       )    ! vx vz
  EM_SRC(8)  = GAS_Q/GAS_M*( 2*Ey*sol(3) - 2*B*sol(6)                                   )    ! vy vy
  EM_SRC(9)  = GAS_Q/GAS_M*( Ey*sol(4) - B*sol(7)                                       )    ! vy vz
  EM_SRC(10) = GAS_Q/GAS_M*( 0.0d0                                                      )    ! vz vz
  
  EM_SRC(11) = GAS_Q/GAS_M*( Ex*(3*sol(5) + sol(8) + sol(10)) + 2*Ey*sol(6) + B*sol(12) )    ! vx v^2
  EM_SRC(12) = GAS_Q/GAS_M*( 2*Ex*sol(6) + Ey*(sol(5) + 3*sol(8) + sol(10)) - B*sol(11) )    ! vy v^2
  EM_SRC(13) = GAS_Q/GAS_M*( 2*Ex*sol(7) + 2*Ey*sol(9)                                  )    ! vz v^2

  EM_SRC(14) = GAS_Q/GAS_M*( 4*Ex*sol(11) + 4*Ey*sol(12)                                )    ! v^4

  END SUBROUTINE MOM14_EM_SOURCES

  ! ****************************************************************************************

  SUBROUTINE MOM14_CLOSING_INTERPOLATIVE(rho, PMAT, qVEC, Riijj, sig, Qijk_VECT, Rijkk_VECT, Sijjkk_VECT)
  ! Computes interpolative closing moments for the 14-mom closure

  IMPLICIT NONE
  
  REAL(KIND=8),                 INTENT(IN)  :: rho, Riijj
  REAL(KIND=8), DIMENSION(3,3), INTENT(IN)  :: PMAT
  REAL(KIND=8), DIMENSION(3),   INTENT(IN)  :: qVEC
  REAL(KIND=8),                 INTENT(OUT) :: sig

  REAL(KIND=8), DIMENSION(10), INTENT(OUT) :: Qijk_VECT  
  REAL(KIND=8), DIMENSION(6),  INTENT(OUT) :: Rijkk_VECT 
  REAL(KIND=8), DIMENSION(3),  INTENT(OUT) :: Sijjkk_VECT

  ! Working variables
  REAL(KIND=8) :: A, TR_PMAT, TR_PMAT2, TR_PMAT3, q_invP_q, POW_SIG
  REAL(KIND=8), DIMENSION(3,3) :: invPMAT, BMAT, invBMAT, WMAT, PMAT2, PMAT3, PMAT4
  REAL(KIND=8), DIMENSION(3)   :: invP_q, Wq
  REAL(KIND=8), DIMENSION(3,3) :: Q1_kl, Q2_kl, Q3_kl

  LOGICAL :: DUMMY_BOOL

  ! Closing moments
  REAL(KIND=8) :: Qxxx, Qxxy, Qxyy, Qyyy, Qyyz, Qyzz, Qzzz, Qxxz, Qxzz, Qxyz
  REAL(KIND=8) :: Rxxkk, Rxykk, Ryykk, Ryzkk, Rzzkk, Rxzkk, Sxjjkk, Syjjkk, Szjjkk

  ! Compute sigma, parameter of the parabolic mapping

  ! Define some stuff
  PMAT2    = MATMUL(PMAT, PMAT)
  PMAT3    = MATMUL(PMAT, PMAT2)
  PMAT4    = MATMUL(PMAT, PMAT3)

  TR_PMAT  = PMAT(1,1) + PMAT(2,2) + PMAT(3,3)
  TR_PMAT2 = PMAT2(1,1) + PMAT2(2,2) + PMAT2(3,3)
  TR_PMAT3 = PMAT3(1,1) + PMAT3(2,2) + PMAT3(3,3)

  CALL M33INV(PMAT, invPMAT, DUMMY_BOOL) ! Compute inverse of 3x3 matrix

  IF (DUMMY_BOOL .EQV. .False.) THEN
    PRINT*, "ATTENTION! COULD NOT INVERT PRESSURE TENSOR!!!!!!!"
  END IF

  invP_q   = MATMUL(invPMAT, qVEC)
  q_invP_q = DOT_PRODUCT(qVEC, invP_q)

  A   = 2*TR_PMAT2 + TR_PMAT**2 - rho*Riijj
  SIG = (A + SQRT(A**2 + 8*rho*TR_PMAT2*q_invP_q))/(4.0d0*TR_PMAT2)
  SIG = MAX(SIG, 1.0E-5) ! Artificial limitation on sigma
  SIG = MIN(SIG, 1.0d0)  ! Limit sigma to 1.0 (physical realizability boundary)

  ! Check physical realizability
  IF (Riijj .LT. (q_invP_q + TR_PMAT**2/rho) ) THEN
    WRITE(*,*) "ATTENTION!! Riijj APPEARS TO BE OUT OF PHYSICAL REALIZABILITY BOUNDARY!!!!"
    WRITE(*,*) "PROBABLY THE SOLUTION IS JUST CRASHING...."
  END IF

  ! Compute working stuff
  BMAT = 2*PMAT*TR_PMAT2 + 4*PMAT3
  CALL M33INV(BMAT, invBMAT, DUMMY_BOOL)

  IF (DUMMY_BOOL .EQV. .False.) THEN
    PRINT*, "ATTENTION! COULD NOT INVERT BMAT!!!"
  END IF

  ! Compute Qijk
  Qxxx = Q_ijk(PMAT, PMAT2, invBMAT, qVEC, 1, 1, 1);
  Qxxy = Q_ijk(PMAT, PMAT2, invBMAT, qVEC, 1, 1, 2);
  Qxyy = Q_ijk(PMAT, PMAT2, invBMAT, qVEC, 1, 2, 2);
  Qyyy = Q_ijk(PMAT, PMAT2, invBMAT, qVEC, 2, 2, 2);
  Qyyz = Q_ijk(PMAT, PMAT2, invBMAT, qVEC, 2, 2, 3);
  Qyzz = Q_ijk(PMAT, PMAT2, invBMAT, qVEC, 2, 3, 3);
  Qzzz = Q_ijk(PMAT, PMAT2, invBMAT, qVEC, 3, 3, 3);
  Qxxz = Q_ijk(PMAT, PMAT2, invBMAT, qVEC, 1, 1, 3);
  Qxzz = Q_ijk(PMAT, PMAT2, invBMAT, qVEC, 1, 3, 3);
  Qxyz = Q_ijk(PMAT, PMAT2, invBMAT, qVEC, 1, 2, 3);

  ! Compute Rijkk
  Rxxkk = R_ijkk(rho, SIG, Qxxx, Qxxy, Qxxz, invP_q, PMAT, PMAT2, 1, 1);
  Rxykk = R_ijkk(rho, SIG, Qxxy, Qxyy, Qxyz, invP_q, PMAT, PMAT2, 1, 2);
  Ryykk = R_ijkk(rho, SIG, Qxyy, Qyyy, Qyyz, invP_q, PMAT, PMAT2, 2, 2);
  Ryzkk = R_ijkk(rho, SIG, Qxyz, Qyyz, Qyzz, invP_q, PMAT, PMAT2, 2, 3);
  Rzzkk = R_ijkk(rho, SIG, Qxzz, Qyzz, Qzzz, invP_q, PMAT, PMAT2, 3, 3);
  Rxzkk = R_ijkk(rho, SIG, Qxxz, Qxyz, Qxzz, invP_q, PMAT, PMAT2, 1, 3);

  ! Compute Sijjkk
!!!!!!!!!!  !!!! WMAT = 2*(TR_PMAT**3)*MATMUL(PMAT, invBMAT) &
!!!!!!!!!!  !!!!      + 12*(TR_PMAT3)*MATMUL(PMAT, invBMAT)  &
!!!!!!!!!!  !!!!      + 14*(TR_PMAT2)*MATMUL(PMAT2, invBMAT) & 
!!!!!!!!!!  !!!!      + 20*TR_PMAT*MATMUL(PMAT3, invBMAT)    &
!!!!!!!!!!  !!!!      + .......
!!!!!!!!!!
!!!!!!!!!!  !!!! WMAT = WMAT/(rho+1.0e-35) ! And.. divide!
!!!!!!!!!!
!!!!!!!!!!! AAAAAAAAAAAA OLD OLD OL
!!!!!!!!!!  ! WRONG!!!  WMAT = 1.0d0/(rho+1.0e-35)*(2*(TR_PMAT**3)* PMAT + 12*TR_PMAT3*PMAT + 14*TR_PMAT2*PMAT2  &
!!!!!!!!!!  ! WRONG!!!              + 20*TR_PMAT*PMAT3 + 20*PMAT4 - 2*TR_PMAT2*TR_PMAT*PMAT &
!!!!!!!!!!  ! WRONG!!!              - 6*(TR_PMAT**2)*PMAT2)*invBMAT;
  WMAT = 1.0d0/(rho+1.0e-35) *   &
          MATMUL(2*(TR_PMAT**3)* PMAT + 12*TR_PMAT3*PMAT + 14*TR_PMAT2*PMAT2  &
                 + 20*TR_PMAT*PMAT3 + 20*PMAT4 - 2*TR_PMAT2*TR_PMAT*PMAT         &
                 - 6*(TR_PMAT**2)*PMAT2, invBMAT)

  Wq = MATMUL(WMAT, qVEC);

  Q1_kl(1,:) = (/ Qxxx, Qxxy, Qxxz /) 
  Q1_kl(2,:) = (/ Qxxy, Qxyy, Qxyz /)
  Q1_kl(3,:) = (/ Qxxz, Qxyz, Qxzz /)
 
  Q2_kl(1,:) = (/ Qxxy, Qxyy, Qxyz /)
  Q2_kl(2,:) = (/ Qxyy, Qyyy, Qyyz /)
  Q2_kl(3,:) = (/ Qxyz, Qyyz, Qyzz /)

  Q3_kl(1,:) = (/ Qxxz, Qxyz, Qxzz /)
  Q3_kl(2,:) = (/ Qxyz, Qyyz, Qyzz /)
  Q3_kl(3,:) = (/ Qxzz, Qyzz, Qzzz /)
 
  ! pow_SIG = 3.0d0/5.0d0; ! old value
  POW_SIG = 1.0d0/2.0d0; ! Watch out for divisions!
 
  Sxjjkk = 1.0d0/(SIG**2.0d0)*DOT_PRODUCT(invP_q, MATMUL(Q1_kl, invP_q)) &
         + 2.0d0*(SIG**pow_SIG)*TR_PMAT*qVEC(1)/(rho+1.0e-35)            &
         + (1.0d0 - SIG**pow_SIG)*Wq(1)

  Syjjkk = 1.0d0/(SIG**2.0d0)*DOT_PRODUCT(invP_q, MATMUL(Q2_kl, invP_q)) &
         + 2.0d0*(SIG**pow_SIG)*TR_PMAT*qVEC(2)/(rho+1.0e-35)            &
         + (1.0d0 - SIG**pow_SIG)*Wq(2)

  Szjjkk = 1.0d0/(SIG**2.0d0)*DOT_PRODUCT(invP_q, MATMUL(Q3_kl, invP_q)) &
         + 2.0d0*(SIG**pow_SIG)*TR_PMAT*qVEC(3)/(rho+1.0e-35)            &
         + (1.0d0 - SIG**pow_SIG)*Wq(3)

  ! Compose arrays to output solution
  Qijk_VECT   = (/ Qxxx, Qxxy, Qxyy, Qyyy, Qyyz, Qyzz, Qzzz, Qxxz, Qxzz, Qxyz /)
  Rijkk_VECT  = (/ Rxxkk, Rxykk, Rxzkk, Ryykk, Ryzkk, Rzzkk /)
  Sijjkk_VECT = (/ Sxjjkk, Syjjkk, Szjjkk /)

  END SUBROUTINE MOM14_CLOSING_INTERPOLATIVE
 
  ! *****************************************************************************************

  FUNCTION R_ijkk(rho, SIG, Q_ij1, Q_ij2, Q_ij3, invPq, PMAT, P2MAT, i, j)
    ! Computes the element of the vector Rijkk, with indices (i,j) 
    ! (kk is a sum in the Einstein notation), for the interpolative closure 
    ! of McDonald & Torrilhon, 2013

    REAL(KIND=8) :: rho, SIG, Q_ij1, Q_ij2, Q_ij3
    REAL(KIND=8), DIMENSION(3), INTENT(IN)   :: invPq
    REAL(KIND=8), DIMENSION(3,3), INTENT(IN) :: PMAT
    REAL(KIND=8), DIMENSION(3,3), INTENT(IN) :: P2MAT
    INTEGER, INTENT(IN) :: i, j

    REAL(KIND=8) :: R_ijkk

    R_ijkk = 1.0d0/(SIG+1.0e-35)*(Q_ij1*invPq(1) + Q_ij2*invPq(2) + Q_ij3*invPq(3)) &
             + (2*(1.0d0-SIG)*P2MAT(i,j) + PMAT(i,j)*(PMAT(1,1)+PMAT(2,2)+PMAT(3,3)))/(rho + 1.0e-35);

    RETURN

  END FUNCTION R_ijkk

  ! *****************************************************************************************

  FUNCTION Q_ijk(PMAT, P2MAT, invBMAT, qVEC, i, j, k)
    ! Computes the matrix element Q(i,j,k) for the interpolative closure - McDonald & Torrilhon, 2013

    REAL(KIND=8), DIMENSION(3,3), INTENT(IN) :: PMAT
    REAL(KIND=8), DIMENSION(3,3), INTENT(IN) :: P2MAT
    REAL(KIND=8), DIMENSION(3,3), INTENT(IN) :: invBMAT
    REAL(KIND=8), DIMENSION(3), INTENT(IN)   :: qVEC
    INTEGER, INTENT(IN) :: i, j, k

    REAL(KIND=8) :: Q_ijk

    Q_ijk =   K_ijkm(PMAT, P2MAT, invBMAT, i,j,k,1)*qVEC(1) &
            + K_ijkm(PMAT, P2MAT, invBMAT, i,j,k,2)*qVEC(2) &
            + K_ijkm(PMAT, P2MAT, invBMAT, i,j,k,3)*qVEC(3)

    RETURN

  END FUNCTION Q_ijk

  ! *****************************************************************************************

  FUNCTION K_ijkm(PMAT, P2MAT, invBMAT, i, j, k, m)
    ! Computes the matrix element K(i,j,k,m) for the interpolative closure - McDonald & Torrilhon, 2013

    REAL(KIND=8), DIMENSION(3,3), INTENT(IN) :: PMAT
    REAL(KIND=8), DIMENSION(3,3), INTENT(IN) :: P2MAT
    REAL(KIND=8), DIMENSION(3,3), INTENT(IN) :: invBMAT
    INTEGER, INTENT(IN) :: i, j, k, m

    REAL(KIND=8) :: K_ijkm

    K_ijkm =  (2*PMAT(i,1)*P2MAT(j,k) + 2*PMAT(k,1)*P2MAT(i,j) + 2*PMAT(j,1)*P2MAT(i,k))*invBMAT(1,m) & 
             +(2*PMAT(i,2)*P2MAT(j,k) + 2*PMAT(k,2)*P2MAT(i,j) + 2*PMAT(j,2)*P2MAT(i,k))*invBMAT(2,m) & 
             +(2*PMAT(i,3)*P2MAT(j,k) + 2*PMAT(k,3)*P2MAT(i,j) + 2*PMAT(j,3)*P2MAT(i,k))*invBMAT(3,m)

    RETURN

  END FUNCTION K_ijkm

  ! *****************************************************************************************

  SUBROUTINE MOM14_MIN_MAX_WAVESPEEDS(sol, sMin, sMax)

  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8),                  INTENT(OUT) :: sMin, sMax

  REAL(KIND=8) :: rho, rhoux, ux, Pxx, ax, MULT_CONST
 
  ! Unpack the solution and compute primitive variables
  rho    = sol(1)
  rhoux  = sol(2)
  ux     = rhoux/(rho + 1.0e-35)
  Pxx    = sol(5)  - rhoux*ux
  ax     = SQRT(GAS_GAMMA*Pxx/(rho+1.0e-35)) ! Speed of sound along x - IS IT CORRECT TO USE Pxx??

  ! We arbitrarily scale the speed of sound, by a constant.
  ! Using something in the order of 1 - 10 should do, but notice that it can really 
  ! be much different, depending on the problem! 
  ! This will introduce unphysical dissipation anyway, so you should at least try different 
  ! values.
  MULT_CONST = 5.0d0 

  sMin = ux - ax*MULT_CONST
  sMax = ux + ax*MULT_CONST

  END SUBROUTINE MOM14_MIN_MAX_WAVESPEEDS

END MODULE PDE_MOM14
