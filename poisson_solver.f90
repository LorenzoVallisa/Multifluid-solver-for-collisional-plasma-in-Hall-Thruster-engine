MODULE poisson_solver

  USE GLOBAL
  use EM_FIELDS

  implicit none

  REAL(KIND=8), allocatable, DIMENSION(:) :: DL,DU,IPIV
  REAL(KIND=8), allocatable, DIMENSION(:) :: D
  REAL(KIND=8), allocatable, DIMENSION(:,:) :: AB
  REAL(KIND=8), allocatable, DIMENSION(:) :: n_background ! Background ions density
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Predicted_E



  CONTAINS

  ! Auxialiary function used into Poisson Matrix assembly

  function alfa(index) result(out)

    implicit none

    real(kind=8) :: out
    integer, intent(in) :: index

    out = 2/(l_cells(index+1)+l_cells(index+2))

  end function alfa

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_matrix

    implicit none

    integer :: N,i
    real(kind=8) :: dx

    if (BCs_PERIODIC_BOOL_EM) then

        ! Common Formulation tridiagonal matrix with values on interfaces only of phi-Formulation with tridiagonality broken

        N = N_int
        KL = N-1
        KU = N-1
        LDAB = 2*KL + KU + 1
        allocate(AB(LDAB,N),IPIV(N),D(N),DL(N-1),DU(N-1))

        DO i = 1, N-1

          D(i) = -alfa(i-1)*(1/l_cells(i) + 1/l_cells(i+1))
          DU(i) = alfa(i-1)/l_cells(i+1)
          DL(i) = alfa(i)/l_cells(i+1)

        END DO

        D(N) = -alfa(N-1)*(1/l_cells(N) + 1/l_cells(N+1))

        ! Formulation with tridiagonality broken

        AB = 0

        DO i = 1, N-1

          AB(KL+KU,i+1) = DU(i)
          AB(KL+KU+1,i) = D(i)
          AB(KL+KU+2,i) = DL(i)

        end do

        AB(KL+KU+1,N) = -1
        AB(KL+KU+N,1) = 1
        AB(KL+KU+2-N,N) = D(N)
        AB(KL+KU+2-(N-1),N-1) = DL(N-1)
        AB(KL+KU+2,N-1) = 0

    else if ((BCs_DIR_BOOL_Right_EM .and. BCs_NEU_BOOL_Left_EM) .or. (BCs_DIR_BOOL_Left_EM .and. BCs_NEU_BOOL_Right_EM)) then


      N = N_int
      allocate(D(N),DL(N-1),DU(N-1))


      DO i = 1, N-1

        D(i) = -alfa(i-1)*(1/l_cells(i) + 1/l_cells(i+1))
        DU(i) = alfa(i-1)/l_cells(i+1)
        DL(i) = alfa(i)/l_cells(i+1)

      END DO

      D(N) = -alfa(N-1)*(1/l_cells(N) + 1/l_cells(N+1))


    else ! Both Dirichlet conditions


      N = N_int - 2
      allocate(D(N),DL(N-1),DU(N-1))


      DO i = 1, N-1


      ! Note that l_cells(1) is the ghost_cell size, so indexes are shifted by 1

        D(i) = -alfa(i)*(1/l_cells(i+1) + 1/l_cells(i+2))
        DU(i) = alfa(i)/l_cells(i+2)
        DL(i) = alfa(i+1)/l_cells(i+2)

        END DO

      D(N) = -alfa(N)*(1/l_cells(N+1) + 1/l_cells(N+2))

    end if

      !else if (BCs_OPEN_BOOL .AND. (GRID_TYPE .eq. "uniform")) THEN

      !N = N_int - 2
      !allocate(D(N),DL(N-1),DU(N-1))

      ! Taking any random interval
      !dx = l_cells(3)-l_cells(2);

      !DO i = 1, N-1

      !  D(i) = -2.0/(dx*dx)
      !  DU(i) = 1.0/(dx*dx)
      !  DL(i) = 1.0/(dx*dx)

      !END DO

      !D(N) = -2.0/(dx*dx)

      !end if

      ! Allocating n_background as well
      allocate(n_background(N_cells))
      n_background = n_ions_background
      allocate(Predicted_E(N_cells))
      Predicted_E = 0.0d0

      if (nondim_bool) then

          ! Nondimensionalization of background ions density

          n_background = n_background/rho_i_carac

      end if

    end subroutine init_matrix


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remove_numerical_noise(force_term_nodes)

  implicit none

  REAL(KIND=8), DIMENSION(N_int), INTENT(INOUT) :: force_term_nodes
  integer :: i
  real(kind=8) :: sum = 0.00d0

  DO i = 1, N_cells

    sum = sum + force_term_nodes(i)

  END DO

  sum = sum/(real(N_cells,8))

  force_term_nodes = force_term_nodes - sum

  end subroutine remove_numerical_noise

  ! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  subroutine Poisson_prediction(sol,sol_n_minus1)

  implicit none

  REAL(KIND=8), DIMENSION(8,N_cells), INTENT(IN) :: sol,sol_n_minus1
  REAL(KIND=8), DIMENSION(N_cells) :: E_n,force_term_nodes
  REAL(KIND=8), DIMENSION(N_int) :: E_n_inter,nu_coll_minus1_el,force_term_inter,nu_coll_minus1_io
  REAL(KIND=8), DIMENSION(N_int) ::  n_n_ions,n_n_electrons,n_nminus1_ions,n_nminus1_electrons,Grad_rho_electrons
  REAL(KIND=8), DIMENSION(N_int) ::  rho_n_ions,rho_n_electrons,rho_nminus1_ions,rho_nminus1_electrons,Grad_rho_ions
  REAL(KIND=8), DIMENSION(N_int) :: mu_mobility_inter_el,mu_mobility_inter_io,Grad_n_electrons,Grad_n_ions
  REAL(KIND=8), DIMENSION(N_int) :: CHI_E,phi_nplus1,Grad_CHI_E
  REAL(KIND=8), DIMENSION(N_int) :: CHI_EL_AD,CHI_IO_AD,Grad_CHI_EL_AD,Grad_CHI_IO_AD
  REAL(KIND=8), allocatable, DIMENSION(:) :: BIG_DL,BIG_DU
  REAL(KIND=8), allocatable, DIMENSION(:) :: BIG_D,B
  INTEGER :: NRHS = 1
  INTEGER :: i,LDB,INFO
  integer :: N
  REAL(KIND=8), DIMENSION(N_int) :: BIG_RHS

  ! Intializing Electric field at timestep n
  E_n = 0.0d0
  E_n(2:N_cells-1) = E_SC

  ! Electirc field at interface ready for LS rhs
  call interp_interface(E_n,E_n_inter)

  ! nu_coll_el_minus1 at interface ready for LS rhs
  call interp_interface(coll_elastic_vect_el,nu_coll_minus1_el)
  call interp_interface(coll_elastic_vect_io,nu_coll_minus1_io)

  if (NONDIM_BOOL)THEN

      ! Retrieve nondimensional densities ad gradients for both pops on the domain interfaces
      call interp_interface(sol(1,:),rho_n_electrons)
      call interp_interface(sol(5,:),rho_n_ions)
      call interp_interface(sol_n_minus1(1,:),rho_nminus1_electrons)
      call interp_interface(sol_n_minus1(5,:),rho_nminus1_ions)
      call compute_Grad_n(sol(1,:),Grad_rho_electrons)
      call compute_Grad_n(sol(5,:),Grad_rho_ions)

      CHI_EL_AD = dt/nu_coll_minus1_el * rho_n_electrons
      Grad_CHI_EL_AD = dt/nu_coll_minus1_el * Grad_rho_electrons

      CHI_IO_AD = dt/nu_coll_minus1_io * rho_n_ions
      Grad_CHI_IO_AD = dt/nu_coll_minus1_io * Grad_rho_ions

      ! Assembling Big_A

      N = N_int - 2
      allocate(BIG_D(N),BIG_DL(N-1),BIG_DU(N-1),B(N))


      DO i = 1, N-1

        BIG_D(i) = (-alfa(i)*(1/l_cells(i+1) + 1/l_cells(i+2)))*(eps_0 + Predictor_1 * CHI_EL_AD(i+1) + &
                    Predictor_2 * CHI_IO_AD(i+1))

        BIG_DU(i) = alfa(i)/l_cells(i+2)*(eps_0 + Predictor_1 * CHI_EL_AD(i+1) + Predictor_2 * CHI_IO_AD(i+1) ) - &
                    (Grad_CHI_EL_AD(i+1)* Predictor_1 + Grad_CHI_IO_AD(i+1)* Predictor_2)/(2*l_cells(i+2))

        BIG_DL(i) = alfa(i+1)/l_cells(i+2)*(eps_0 + Predictor_1 * CHI_EL_AD(i+2) + Predictor_2 * CHI_IO_AD(i+2)) + &
                    (Grad_CHI_EL_AD(i+2)* Predictor_1 + Grad_CHI_IO_AD(i+2)* Predictor_2)/(2*l_cells(i+2))

      END DO

      BIG_D(N) = (-alfa(N)*(1/l_cells(N+1) + 1/l_cells(N+2)))*(eps_0 + Predictor_1 * CHI_EL_AD(N_int-1) + &
                  Predictor_2 * CHI_IO_AD(N_int-1))

      ! Assembling Big_RHS

      BIG_RHS = - (Predictor_3 * (2 * rho_n_electrons - rho_nminus1_electrons) + Predictor_4 * ( 2 * rho_n_ions - &
            rho_nminus1_ions) + &
            ( Predictor_1 * Grad_CHI_EL_AD + Predictor_2 * Grad_CHI_IO_AD) * 1 * E_n_inter + &
            ( Predictor_1 * CHI_EL_AD + Predictor_2 * CHI_IO_AD) * (kobe * rho_n_electrons - rho_n_ions)/shaq)


      B = BIG_RHS(2:N_int-1)

      B(1) = B(1) - (alfa(1)/l_cells(2) * (eps_0 + Predictor_1 * CHI_EL_AD(2) + Predictor_2 * CHI_IO_AD(2) ) + &
                    (Grad_CHI_EL_AD(2)* Predictor_1 + Grad_CHI_IO_AD(2)* Predictor_2)/(2*l_cells(2)))* phi_0

      B(N) = B(N) - (alfa(N)/l_cells(N+2) * (eps_0 + Predictor_1 * CHI_EL_AD(N+1) + Predictor_2 * CHI_IO_AD(N+1) ) + &
                    (Grad_CHI_EL_AD(N+1)* Predictor_1 + Grad_CHI_IO_AD(N+1)* Predictor_2)/(2*l_cells(N+1)))*phi_end


  else



      mu_mobility_inter_el = Gas_q_1 /( GAS_M_1 * nu_coll_minus1_el)
      mu_mobility_inter_io = GAS_Q_2 /(Gas_M_2 * nu_coll_minus1_io)

      ! Retrieve number densities for both pops
      call interp_interface(sol(1,:)/GAS_M_1,n_n_electrons)
      call interp_interface(sol(5,:)/GAS_M_2,n_n_ions)
      call interp_interface(sol_n_minus1(1,:)/GAS_M_1,n_nminus1_electrons)
      call interp_interface(sol_n_minus1(5,:)/GAS_M_2,n_nminus1_ions)

      ! Rho_Q at timestep n
      rho_Q_nodes = (GAS_Q_1*sol(1,:) /GAS_M_1 + GAS_Q_2*sol(5,:) /GAS_M_2)
      force_term_nodes = rho_Q_nodes/eps_0;
      call interp_interface(force_term_nodes,force_term_inter)

      call compute_Grad_n(sol(1,:)/GAS_M_1,Grad_n_electrons)
      call compute_Grad_n(sol(5,:)/GAS_M_2,Grad_n_ions)

      CHI_E = dt/eps_0 * (GAS_Q_1 * mu_mobility_inter_el * n_n_electrons + GAS_Q_2 * mu_mobility_inter_io * n_n_ions)
      Grad_CHI_E = dt/eps_0 *( GAS_Q_1 * mu_mobility_inter_el * Grad_n_electrons + GAS_Q_2 * mu_mobility_inter_io * Grad_n_ions)

      ! Assembling Big_A

      N = N_int - 2
      allocate(BIG_D(N),BIG_DL(N-1),BIG_DU(N-1),B(N))


      DO i = 1, N-1

        BIG_D(i) = (-alfa(i)*(1/l_cells(i+1) + 1/l_cells(i+2)))*(1+CHI_E(i+1))
        BIG_DU(i) = alfa(i)/l_cells(i+2)*(1+CHI_E(i+1)) - Grad_CHI_E(i+1)/(2*l_cells(i+2))
        BIG_DL(i) = alfa(i+1)/l_cells(i+2)*(1 + CHI_E(i+2)) + Grad_CHI_E(i+2)/(2*l_cells(i+2))
        !BIG_D(i) = (-alfa(i)*(1/l_cells(i+1) + 1/l_cells(i+2)))*(1+CHI_E(i+1))
        !BIG_DU(i) = alfa(i)/l_cells(i+2)*(1 + CHI_E(i+2))               ! - Grad_CHI_E(i+1)/(2*l_cells(i+2)))
        !BIG_DL(i) = alfa(i+1)/l_cells(i+2)*(1 + CHI_E(i+1))           !+ Grad_CHI_E(i+2)/(2*l_cells(i+2))

      END DO

      !BIG_D(N) = (-alfa(N)*(1/l_cells(N+1) + 1/l_cells(N+2)))*(1+CHI_E(N_int-1))
      BIG_D(N) = (-alfa(N)*(1/l_cells(N+1) + 1/l_cells(N+2)))*(1+CHI_E(N_int-1))


      ! Assembling Big_RHS

      BIG_RHS = -1/eps_0 * (GAS_Q_1 * (2 * n_n_electrons - n_nminus1_electrons) + GAS_Q_2 * ( 2 * n_n_ions - n_nminus1_ions) + &
            eps_0 * ( CHI_E * force_term_inter + Grad_CHI_E * E_n_inter))




      B = BIG_RHS(2:N_int-1)

      !B(1) = B(1) - (alfa(1)/l_cells(2) * (1 + CHI_E(1)))* phi_0
      !B(N) = B(N) - (alfa(N)/l_cells(N+2) * (1 + CHI_E(N+2)))*phi_end

      B(1) = B(1) - (alfa(1)/l_cells(2) * (1 + CHI_E(2)) + Grad_CHI_E(2)/(2*l_cells(2)))* phi_0
      B(N) = B(N) - (alfa(N)/l_cells(N+2) * (1 + CHI_E(N+1)) - Grad_CHI_E(N+1)/(2*l_cells(N+1)))*phi_end



  end if



  LDB = N

  ! N is order of Matrix, in our case has to be number of internal interfaces (N_cells -3 )
  ! NRHS is the number of column of RHS, indeed 1
  ! DL contains elements of lower-diagonal of A
  ! DU contains elements of upper-diagonal of A
  ! D contains elements of diagonal of A
  ! B is the RHS, which will come out as solution
  ! LDB dimension of array B


  CALL dgtsv( N, NRHS, BIG_DL, BIG_D, BIG_DU, B, LDB, INFO )

  phi_nplus1(1)=phi_0
  phi_nplus1(N_int)=phi_end
  phi_nplus1(2:N_int-1) = B

  call compute_E(phi_nplus1,Predicted_E(2:N_cells-1))

  En_plus_1_plot = Predicted_E

  end subroutine Poisson_prediction


  ! |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


    ! Subroutine used to compute electric field

    subroutine compute_Grad_n(n,Grad_n)

    implicit none

    real(kind=8), dimension(n_int), intent(inout) :: Grad_n
    real(kind=8), dimension(n_cells), intent(in) :: n
    integer :: j

    DO j = 1, n_int

      Grad_n(j) = (n(j+1) - n(j))/l_cells(j)

    END DO

    end subroutine compute_Grad_n

    ! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ! Subroutine to generically interpolate potential on interfaces

  subroutine interp_interface(term_nodes , term_inter)

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(N_cells), INTENT(IN) :: term_nodes
    REAL(KIND=8), DIMENSION(N_int), INTENT(INOUT) :: term_inter
    INTEGER :: i
    REAL(KIND=8) :: dx_r,dx_l,dx

    DO i = 1, N_int

    !    GC----|------°------|-----°-----|----GC
    !                 2 dx_l   dx_r 3


      dx_l = l_cells(i)/2
      dx_r = l_cells(i+1)/2
      dx = dx_l+dx_r
      term_inter(i) = ( term_nodes(i)*dx_l/dx + term_nodes(i+1)*dx_r/dx )


    END DO

  end subroutine interp_interface


! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ! Subroutine to compute potential using lapack linear system solver for poisson problem

  subroutine compute_phi(force_term_inter ,force_term_nodes ,phi)

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(N_int), INTENT(INOUT) :: phi
    REAL(KIND=8), DIMENSION(N_int), INTENT(IN) :: force_term_inter
    real(KIND=8), dimension(n_cells), intent(in) :: force_term_nodes
    INTEGER :: NRHS = 1
    REAL(KIND=8), allocatable, DIMENSION(:) :: DL_,DU_
    REAL(KIND=8), allocatable, DIMENSION(:) :: D_,B
    INTEGER :: i,LDB,INFO
    REAL(KIND=8), allocatable, DIMENSION(:,:) :: ABC
    integer :: n


    if (BCs_PERIODIC_BOOL_EM) then

        ! Common Formulation tridiagonal matrix with values on interfaces only of phi-Formulation with tridiagonality broken

        N = N_int
        allocate(B(N))
        LDB = N
        B = force_term_inter
        B(1) = B(1) + B(N)
        B(N) = 0
        allocate(ABC(LDAB,N))
        ABC = AB

        ! N is order of Matrix, in our case has to be number of interfaces (N_int)
        ! NRHS is the number of column of RHS, indeed 1
        ! D contains elements of diagonal of A
        ! B is the RHS, which will come out as solution
        ! LDB dimension of array B
        ! KL subdiag
        ! KU upperdiag
        ! ABC has matrix A stored somehow inside
        ! IPIV pivoting array


        call dgbsv (N, KL, KU, NRHS, ABC, LDAB, IPIV, B, LDB, INFO)

    else if (BCs_DIR_BOOL_Left_EM .and. BCs_NEU_BOOL_Right_EM) then

        !!!! Neumann on the right - Dirichlet on the left

        N = N_int
        allocate(D_(N),DL_(N-1),DU_(N-1),B(N))
        B = force_term_inter(2:N_int-1)

        LDB = N
        DL_ = DL
        D_ = D
        DU_ = DU


        ! Neumann Right - Dirichlet Left
        B(1) = phi_0
        B(N) = Neu_right * l_cells(N)


        D(1) = 1.00d0
        DU(1) = 0.0d0
        D(N) = -1.0d0
        DL(N-1) = 1.0d0

        ! N is order of Matrix, in our case has to be number of internal interfaces (N_cells -3 )
        ! NRHS is the number of column of RHS, indeed 1
        ! DL contains elements of lower-diagonal of A
        ! DU contains elements of upper-diagonal of A
        ! D contains elements of diagonal of A
        ! B is the RHS, which will come out as solution
        ! LDB dimension of array B


        CALL dgtsv( N, NRHS, DL_, D_, DU_, B, LDB, INFO )

    else if (BCs_DIR_BOOL_Right_EM .and. BCs_NEU_BOOL_Left_EM) then

            !!!! Neumann on the left - Dirichlet on the right

            N = N_int
            allocate(D_(N),DL_(N-1),DU_(N-1),B(N))
            B = force_term_inter(2:N_int-1)

            LDB = N
            DL_ = DL
            D_ = D
            DU_ = DU

            ! Neumann Left - Dirichlet Right
            B(1) = Neu_left * l_cells(2)
            B(N) = phi_end

            D(1) = 1.00d0
            DU(1) = -1.0d0
            D(N) = 1.00d0
            DL(N-1) = 0.000d000


            ! N is order of Matrix, in our case has to be number of internal interfaces (N_cells -3 )
            ! NRHS is the number of column of RHS, indeed 1
            ! DL contains elements of lower-diagonal of A
            ! DU contains elements of upper-diagonal of A
            ! D contains elements of diagonal of A
            ! B is the RHS, which will come out as solution
            ! LDB dimension of array B


            CALL dgtsv( N, NRHS, DL_, D_, DU_, B, LDB, INFO )

    else


        N = N_int-2
        allocate(D_(N),DL_(N-1),DU_(N-1),B(N))
        B = force_term_inter(2:N_int-1)
        B(1) = B(1) - alfa(1)/l_cells(2)*phi_0
        B(N) = B(N) - alfa(N)/l_cells(N+2)*phi_end

        LDB = N
        DL_ = DL
        D_ = D
        DU_ = DU

        ! N is order of Matrix, in our case has to be number of internal interfaces (N_cells -3 )
        ! NRHS is the number of column of RHS, indeed 1
        ! DL contains elements of lower-diagonal of A
        ! DU contains elements of upper-diagonal of A
        ! D contains elements of diagonal of A
        ! B is the RHS, which will come out as solution
        ! LDB dimension of array B


        CALL dgtsv( N, NRHS, DL_, D_, DU_, B, LDB, INFO )

    end if


    if ( (info.LT.0) .or. (info.GT.0)  ) then
        WRITE(*,*) "ERROR TYPE: IMMINENT EXPLOSION "
        stop
    end if

    if (BCs_PERIODIC_BOOL_EM) then

      phi = B

    else if ((BCs_DIR_BOOL_Right_EM .and. BCs_NEU_BOOL_Left_EM) .or. (BCs_DIR_BOOL_Left_EM .and. BCs_NEU_BOOL_Right_EM)) then

      phi = B

    ELSE

      phi(1)=phi_0
      phi(N_int)=phi_end
      phi(2:N_int-1) = B

    end if

    !write(*,*)" ------------------------------------------------------------------"
    !write(*,*)" BC left ------------>  ",phi_0
    !write(*,*)" BC right ----------->  ",phi_end
    !write(*,*)" ------------------------------------------------------------------"


    end subroutine compute_phi


! |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||



  ! Subroutine used to compute electric field

  subroutine compute_E(phi,E_SelfConsistent)

  implicit none

  real(kind=8), dimension(n_int), intent(in) :: phi
  real(kind=8), dimension(n_cells-2), intent(inout) :: E_SelfConsistent
  integer :: j

  DO j = 1, n_cells-2

    E_SelfConsistent(j) = -(phi(j+1) - phi(j))/l_cells(j+1)

  END DO

  end subroutine compute_E

  ! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


  ! Subroutine used into time_integrators.f90 to access Electric field computation


  subroutine compute_SC_EM_field(sol,E_SelfConsistent,Ey,B)

    implicit none

  REAL(KIND=8), DIMENSION(N_EQ,N_cells), INTENT(IN) :: sol
  REAL(KIND=8), DIMENSION(N_cells) :: sol_e,sol_i
  REAL(KIND=8), DIMENSION(N_cells-2), INTENT(INOUT) :: E_SelfConsistent
  REAL(KIND=8), DIMENSION(N_cells), INTENT(INOUT) :: b,Ey

  E_SelfConsistent = 0

  ! Fixed parameters updated from EM_FIELDS module
  Ey = Ey_field
  B = B_field

  if(SYS_Euler_MF_BOOL)THEN

    sol_i = sol(5,:)
    sol_e = sol(1,:)

  else if(SYS_EULER_BOOL .and. FLUID_IONS)THEN

    sol_i = sol(1,:)
    sol_e = Density_electrons_vector

  else if(SYS_EULER_BOOL .and. FLUID_ELECTRONS)THEN

    sol_i = n_background
    sol_e = sol(1,:)

  end if

  if( SC_EM_field ) then

    call update_SC_EM(sol_e,E_SelfConsistent,sol_i)

  end if

  E_SelfConsistent = E_SelfConsistent + Ex_field(2:n_cells-1)

  E_SC_plot(2:N_cells-1) = E_SelfConsistent
  E_SC = E_SelfConsistent

  end subroutine compute_SC_EM_field

  ! |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


  ! Subroutine to compute self_consistent_electric field through poisson solver.

  ! Firstly note that output is the electric field over the entire domain (vector long n_cells)
  ! but self-consistent is computed only on internal one: I only need values of variables on
  ! ghost cells for those which are used to compute fluxes

  ! Secondly note that the returned Electric field is the force term PLUS the self-consistent
  ! electric field: this will be the only self-consistent one in case the force term is void.
  ! Therefore at the end of this beautiful method, the self-consistent electric field, freshly computed, will be added to
  ! the input one (that is why it is Initialized as INOUT)


  subroutine update_SC_EM(sol_e,E_SelfConsistent,sol_i)

  IMPLICIT NONE

  ! Input parameter needed is just first conservative variable (Density)


  REAL(KIND=8), DIMENSION(N_cells), INTENT(IN) :: sol_e,sol_i
  REAL(KIND=8), DIMENSION(N_cells-2), INTENT(INOUT) :: E_SelfConsistent
  real(kind=8), dimension(n_cells) :: force_term_nodes
  REAL(KIND=8), DIMENSION(N_int) :: force_term_inter
  integer :: i


  ! Build force term for Poisson Problem


  if(nondim_bool) THEN

    rho_Q_nodes = (kobe*sol_e - sol_i)/shaq

    force_term_nodes = -rho_Q_nodes;

  ELSE

    rho_Q_nodes = (GAS_Q_1*sol_e /GAS_M_1 + GAS_Q_2*sol_i /GAS_M_2)

    force_term_nodes = -rho_Q_nodes/eps_0;

 end if

  ! Remove background noise due to numerical errors due to missed non-dimensionalisation process

  ! call remove_numerical_noise(force_term_nodes)

  ! Interpolate charge density on interface according to grid structure

  call interp_interface(force_term_nodes, force_term_inter)


  ! Solve linear system using lapack in order to obtain potential on all interfaces
  ! DEBUG Verifica1 force_term_inter = 0
  ! DEBUG Verifica2 force_term_inter = -4*3.14159265359*3.14159265359*sin(2*3.14159265359*x_int)
  ! DEBUG Periodic force_term_inter  = -3.75
  ! DEBUG force_term_inter  =  - (-1.602179512121615e-04 + sin(3.14159265359*(x_int+1)))/eps_0  ! - (-1.602179512121615e-04 + 1.602146506484098e-08*sin(3.14159265359*(x_int+1)))/eps_0
  ! DEBUG force_term_inter = + 1.6d-4/eps_0
  ! force_term_inter = -4*3.14159265359*3.14159265359*sin(2*3.14159265359*x_int)


  call compute_phi(force_term_inter,force_term_nodes,phi)

  ! Method used to obtain self-consistent electric field in all internal nodes

  call compute_E(phi,E_SelfConsistent)

  end subroutine update_SC_EM


END MODULE poisson_solver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
