MODULE VARIOUS

  USE GLOBAL
  USE EM_FIELDS

  IMPLICIT NONE
  CONTAINS

  SUBROUTINE PRINT_LOGO

    ! CIAO

  END SUBROUTINE PRINT_LOGO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE WRITE_SOLUTION_TO_FILE(t_ID,U)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: t_ID
  REAL(KIND=8), DIMENSION(N_EQ,N_cells), intent(inout) :: U
  CHARACTER(LEN=512) :: file_name
  INTEGER :: IC
  real(kind=8) :: eps_0 = 8.854187d-12
  real(kind=8) :: max_E,ion_front_vel

  ! Be careful with this implementation that the restarted simulation overrides the old one
  WRITE(file_name,'(A, I8.8, A)') './dumps/solution_', t_ID, '.dat'

  OPEN(12399, FILE=file_name, STATUS="replace", FORM="formatted")




    ! In this way at WRITE_EVERY iterations the solution at time CURRENT_TIME is actually saved

    IF(MOD(t_ID, WRITE_EVERY) .EQ. 0) then

    WRITE(12399, '(A, ES14.7,A)') '# Solution at time ', current_time, ' s.'
    WRITE(12399, '(A)')         '# Simulation time [s]  |  Cell center [m]  |  Solution in conserved variables '

        if (SC_EM_field) then

          DO IC = 2, N_cells-1 ! Write solution in physical cells (no ghost cells)
            WRITE(12399, *) current_time, x_cc(IC), U(:, IC),E_SC_plot(IC), rho_Q_nodes(IC), phi(IC),nu_ions_vect(2,IC)   ! B_field(IC) !, nu_ions_vect(1,IC),nu_ions_vect(2,IC),nu_ions_vect(3,IC) !
          END DO

        else

          DO IC = 2, N_cells-1 ! Write solution in physical cells (no ghost cells)
            WRITE(12399, *) current_time, x_cc(IC), U(:, IC)
          END DO

        end if

   end if

  !  else

  !      ion_front_vel = maxval(U(6,:)/U(5,:))
  !      max_E = maxval(E_total_x)
  !      write(12399,*)  (t_ID+t_old_ID)*dt,max_E,ion_front_vel

  !  end if




  CLOSE(12399)

  END SUBROUTINE WRITE_SOLUTION_TO_FILE



  !***********************************************************************************************************************************
  !  Programmer:   David G. Simpson
  !                NASA Goddard Space Flight Center
  !                Greenbelt, Maryland  20771
  !
  !  M33INV  -  Compute the inverse of a 3x3 matrix.
  !
  !  A       = input 3x3 matrix to be inverted
  !  AINV    = output 3x3 inverse of matrix A
  !  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
  !***********************************************************************************************************************************

  SUBROUTINE M33INV (A, AINV, OK_FLAG)

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(3,3), INTENT(IN)  :: A
  REAL(KIND=8), DIMENSION(3,3), INTENT(OUT) :: AINV
  LOGICAL, INTENT(OUT) :: OK_FLAG

  REAL(KIND=8), PARAMETER :: EPS = 1.0D-15
  REAL(KIND=8) :: DET
  REAL(KIND=8), DIMENSION(3,3) :: COFACTOR


  DET =   A(1,1)*A(2,2)*A(3,3)  &
        - A(1,1)*A(2,3)*A(3,2)  &
        - A(1,2)*A(2,1)*A(3,3)  &
        + A(1,2)*A(2,3)*A(3,1)  &
        + A(1,3)*A(2,1)*A(3,2)  &
        - A(1,3)*A(2,2)*A(3,1)

  IF (ABS(DET) .LE. EPS) THEN
     AINV = 0.0D0
     OK_FLAG = .FALSE.
     RETURN
  END IF

  COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
  COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
  COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
  COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
  COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
  COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
  COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
  COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
  COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

  AINV = TRANSPOSE(COFACTOR) / DET

  OK_FLAG = .TRUE.

  RETURN

  END SUBROUTINE M33INV

  ! ************************************************************************************

  SUBROUTINE MAX_EIG_POWER_ITER(MAT, N, mu)
  ! Computes maximum eigenvalue using power iterations.
  ! Note that it doesn't work too well! JUST A QUICK TEST!!!
  ! JUST A QUICK TEST!!!!
  ! JUST A QUICK TEST!!!!
  ! JUST A QUICK TEST!!!!
  ! JUST A QUICK TEST!!!!
  ! JUST A QUICK TEST!!!!
  ! JUST A QUICK TEST!!!!

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: N
  REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: MAT
  REAL(KIND=8),                 INTENT(OUT) :: mu

  REAL(KIND=8), DIMENSION(N) :: bVEC
  INTEGER :: I
  REAL(KIND=8) :: TOL, ERRR, mu_new

  ! -----------------------------------------
  ! Initialize bVEC randomly, so that it's not aligned to any
  ! eigenvector probably

  bVEC = 0.0d0 ! Init

  DO I = 1, N
    bVEC(I) = bVEC(I) + RAND()
  END DO

  bVEC = bVEC/NORM2(bVEC,1) ! Normalize
  ! ------------------------------------------

  TOL  = 1.0d-5  ! Tolerance on the eigenvalue
  ERRR = 100.0d0 ! Initialize error
  mu   = 0.0d0   ! Initialize

  DO WHILE (ERRR .GT. TOL)

    PRINT*, "Iterating..."

    mu_new = DOT_PRODUCT(bVEC, MATMUL(MAT,bVEC))/NORM2(bVEC) ! Compute mu

    bVEC   = MATMUL(MAT,bVEC) ! Update bvec
    bVEC   = bVEC/NORM2(bVEC) ! Normalize it

    ERRR = ABS(mu_new - mu)   ! Compute error

    mu = mu_new ! Update mu

  END DO

  END SUBROUTINE MAX_EIG_POWER_ITER

  ! ****************************************************************

  SUBROUTINE COMPUTE_MAX_MIN_EIG(JMAT, N, mu_min, mu_max)

  ! mu_min: minimum eigenvalue (most negative one)
  ! mu_max: maximum eigenvalue

  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: JMAT
  INTEGER,      INTENT(IN)  :: N
  REAL(KIND=8), INTENT(OUT) :: mu_min, mu_max

  INTEGER      :: ii
  REAL(KIND=8) :: TMP, EIG_MULT_FACTOR, SHIFT_MAT
  REAL(KIND=8), DIMENSION(N,N) :: JMAT_shift

  ! Shift matrix JMAT by SHIFT_MAT
  SHIFT_MAT = 10.0d0

  JMAT_shift = jMAT
  DO ii = 1, N
    JMAT_shift(ii,ii) = JMAT_shift(ii,ii) + SHIFT_MAT
  END DO

  ! Find maximum eigenvalue
  CALL MAX_EIG_POWER_ITER(JMAT_shift, N, mu_max)

  mu_max = mu_max - SHIFT_MAT

PRINT*, "DONE COMPUTING MAX EIG!", mu_max

  ! Now remove it from the diagonal, to shift all eigenvalues.
  ! Actually, remove more: remove EIG_MULT_FACTOR times it.
  EIG_MULT_FACTOR = 10.0d0

  JMAT_shift = JMAT
  DO ii = 1, N
    JMAT_shift(ii,ii) = JMAT_shift(ii,ii) - EIG_MULT_FACTOR*mu_max - SHIFT_MAT
  END DO

  CALL MAX_EIG_POWER_ITER(JMAT_shift, N, mu_min)

  mu_min = mu_min + EIG_MULT_FACTOR*mu_max + SHIFT_MAT

  ! Order them
  IF (mu_min .GT. mu_max) THEN
    TMP    = mu_max
    mu_max = mu_min
    mu_min = TMP
  END IF

  END SUBROUTINE COMPUTE_MAX_MIN_EIG

  ! **************************************************************

  SUBROUTINE COMPUTE_EIGS_LAPACK_DGEEV(MAT, N, WR, WI, INFO)

  ! INPUT - MAT: ..
  ! INPUT - N: dimension of square matrix MAT
  !
  ! OUTPUT - WR: real part of eigenvalues
  ! OUTPUT - WI: imaginary part of eigenvalues
  ! OUTPUT - INFO: if larger than zero, the algorithm has failed

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: N ! Dimension of matrix "MAT"
  REAL(KIND=8), DIMENSION(N, N), INTENT(IN)  :: MAT
  REAL(KIND=8), DIMENSION(N),    INTENT(OUT) :: WR, WI
  INTEGER, INTENT(OUT) :: INFO

  ! Local working variables
  INTEGER      ::  LDA, LDVL, LDVR
  INTEGER      ::  LWMAX
  PARAMETER        ( LWMAX = 1000 )

  INTEGER      ::  LWORK

  REAL(KIND=8) :: VL(N,N), VR(N,N), WORK( LWMAX )

  LDA  = N
  LDVL = N
  LDVR = N

  ! Query the optimal workspace
  LWORK = -1
  CALL DGEEV( 'N', 'N', N, MAT, LDA, WR, WI, VL, LDVL, &
                VR, LDVR, WORK, LWORK, INFO )
  LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

  ! Solve eigenproblem
  CALL DGEEV( 'N', 'N', N, MAT, LDA, WR, WI, VL, LDVL, &
             VR, LDVR, WORK, LWORK, INFO )

  ! Check for convergence
  IF( INFO.GT.0 ) THEN
     WRITE(*,*)'The algorithm failed to compute eigenvalues.'
  END IF

  END SUBROUTINE COMPUTE_EIGS_LAPACK_DGEEV

END MODULE VARIOUS
