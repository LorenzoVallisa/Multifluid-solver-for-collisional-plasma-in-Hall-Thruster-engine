MODULE PDE

  ! This module provides subroutines that can be called elsewhere, and redirects to the proper
  ! PDEs. With this, one doesn't have to modify everywhere when a new PDE is added, but only
  ! create a module with the proper name.

  USE GLOBAL
  USE PDE_EULER
  USE PDE_MOM5
  USE PDE_MOM14
  USE PDE_MOM14_AXI
  USE VARIOUS
  use pde_euler_multi_fluid

  ! This module contains stuff needed to solve the PDE

  CONTAINS

  ! ==============================================================
  ! === General subroutines (called by other parts of program) ===
  ! ==============================================================

  SUBROUTINE INIT_SOL(sol)
  ! Initializes the solution calling the proper PDE

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(N_EQ, N_cells), INTENT(INOUT) :: sol

  ! Do only one of the following
  IF (SYS_EULER_BOOL) THEN

    CALL EULER_INIT_SOL(sol)

  ELSE IF (SYS_EULER_MF_BOOL) THEN

    CALL EULER_MF_INIT_SOL(sol)

  ELSE IF (SYS_MOM5_BOOL) THEN

    CALL MOM5_INIT_SOL(sol)

  ELSE IF (SYS_MOM14_BOOL) THEN

    CALL MOM14_INIT_SOL(sol)

  ELSE IF (SYS_MOM14_AXI_BOOL) THEN

    CALL MOM14_AXI_INIT_SOL(sol)

  ELSE

    WRITE(*,*) "ATTENTION! NO BOOL SET FOR THE PDE, I DON'T KNOW WHICH FLUX TO COMPUTE! ABORTING"
    STOP

  END IF

  END SUBROUTINE INIT_SOL

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE COMPUTE_PDE_FLUXES(sol, Flux)

  ! This subroutine computes the fluxes for a given solution, for the requested PDE.

  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(OUT) :: Flux

  ! Do only one of the following
  IF (SYS_EULER_BOOL) THEN

    CALL EULER_FLUXES(sol, Flux)

  ELSE IF (SYS_EULER_MF_BOOL) THEN

    CALL EULER_MF_fluxes(sol,Flux)

  ELSE IF (SYS_MOM5_BOOL) THEN

    CALL MOM5_FLUXES(sol, Flux)

  ELSE IF (SYS_MOM14_BOOL) THEN

    CALL MOM14_FLUXES(sol, Flux)

  ELSE IF (SYS_MOM14_AXI_BOOL) THEN

    CALL MOM14_AXI_FLUXES(sol, Flux)

  ELSE

    WRITE(*,*) "ATTENTION! NO BOOL SET FOR THE PDE, I DON'T KNOW WHICH FLUX TO COMPUTE! ABORTING"
    STOP

  END IF

  END SUBROUTINE COMPUTE_PDE_FLUXES

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE COMPUTE_PDE_SOURCES(sol, Ex, Ey, B, Src_EM, Src_COLL,cell_index)

  ! This subroutine computes the sources for a given solution, for the requested PDE.
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
  REAL(KIND=8),                  INTENT(IN)  :: Ex, Ey,B
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(OUT) :: Src_EM
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(OUT) :: Src_COLL
  integer, intent(in) :: cell_index

  Density_electrons = Density_electrons_vector(cell_index)
  Temperature_electrons = Temperature_electrons_vector(cell_index)

  ! Do only one of the following
  IF (SYS_EULER_BOOL) THEN

    CALL EULER_SOURCES(sol, Ex, Ey, B, Src_EM, Src_COLL)

  ELSE IF (SYS_EULER_MF_BOOL) THEN

    call EULER_MF_sources(sol, Ex, Ey, B, Src_EM, Src_COLL)

  ELSE

    WRITE(*,*) "ATTENTION! NO BOOL SET FOR THE PDE, I DON'T KNOW WHICH SOURCES TO COMPUTE! ABORTING"
    STOP

  END IF

  END SUBROUTINE COMPUTE_PDE_SOURCES

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE COMPUTE_MAX_MIN_WAVESPEEDS(sol, sMin, sMax)

      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
      REAL(KIND=8),                  INTENT(inOUT) :: sMin, sMax

      ! Do only one of the following
      IF (SYS_EULER_BOOL) THEN

        CALL EULER_MIN_MAX_WAVESPEEDS(sol, sMin, sMax)

      ELSE IF (SYS_MOM5_BOOL) THEN

        CALL MOM5_MIN_MAX_WAVESPEEDS(sol, sMin, sMax)

      ELSE IF (SYS_MOM14_BOOL) THEN

        CALL MOM14_MIN_MAX_WAVESPEEDS(sol, sMin, sMax)

      ELSE IF (SYS_MOM14_AXI_BOOL) THEN

        CALL MOM14_AXI_MIN_MAX_WAVESPEEDS(sol, sMin, sMax)

      ELSE

        WRITE(*,*) "ATTENTION! NO BOOL SET FOR THE PDE, I DON'T KNOW WHICH SOURCES TO COMPUTE! ABORTING"
        STOP

      END IF

  END SUBROUTINE COMPUTE_MAX_MIN_WAVESPEEDS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE COMPUTE_MAX_MIN_WAVESPEEDS_MF(sol, sMin1, sMax1, sMin2, sMax2)

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: sol
    REAL(KIND=8),                  INTENT(OUT) :: sMin1, sMax1, sMin2, sMax2

    call EULER_2FLUIDS_MIN_MAX_WAVESPEEDS(sol, sMin1, sMax1, sMin2, sMax2)

  end subroutine COMPUTE_MAX_MIN_WAVESPEEDS_MF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine COMPUTE_Perturbed_FLUXES(index,Cons,Pert_cons, Flux)

  implicit NONE

  REAL(KIND=8), DIMENSION(4), INTENT(IN)  :: Cons,Pert_cons
  REAL(KIND=8), DIMENSION(4), INTENT(INout)  :: Flux
  REAL(KIND=8), DIMENSION(4) :: inter_cons
  integer :: index

  inter_cons = cons
  inter_cons(index) = inter_cons(index) + Pert_cons(index)

  CALL EULER_FLUXES(inter_cons, Flux)

  end subroutine COMPUTE_Perturbed_FLUXES

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE PDE
