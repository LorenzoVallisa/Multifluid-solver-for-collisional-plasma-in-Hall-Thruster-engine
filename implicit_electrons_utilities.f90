MODULE IMPLICIT_ELECTRONS_UTILITIES

  USE GLOBAL
  !use omp_lib
  use pde_euler_electrons
  use pde_euler

  implicit none

#include "lisf.h"

  LIS_INTEGER,allocatable :: bptr(:),bindex(:)
  LIS_SCALAR,allocatable :: value_vector(:)
  real(kind=8),dimension(:),allocatable :: LS_rhs,Q_n,Q_n_minus_1
  LIS_MATRIX :: A_MATRIX
  LIS_VECTOR :: x_solution,IMPLICIT_RHS
  LIS_SOLVER :: solver
  LIS_INTEGER :: n_size,bnr,bnc,bnnz,nr
  LIS_PRECON :: my_prec
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: INTER_FLUX
  real(kind=8),dimension(:,:),allocatable :: sol_n_minus1


  CONTAINS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_Iterative_Method_parameters

  implicit none

  ! Setting to 0 all the values of the Jacobian Tensor
  value_vector = 0.0d0

  ! RHS will have dimension of 4 (number of equations) * block rows of implicit matrix
  LS_rhs = 0.0d0

  end subroutine set_Iterative_Method_parameters

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine implicit_step(delta_cons)

  implicit NONE

  REAL(KIND=8), DIMENSION(4, N_cells-2), intent(inout) :: delta_cons
  LIS_INTEGER :: ierr,comm
  integer :: local_index,row_DC_index,col_DC_index,iii
  CHARACTER(LEN=512) :: matname,vecname,vecname1


  matname = 'MATRIX'
  vecname = 'SOLUTION'
  vecname1 = 'RHS'



  ! Assembling main matrix
  call lis_matrix_set_bsr(bnr,bnc,bnnz,bptr,bindex,value_vector,A_MATRIX,ierr)
  call lis_matrix_assemble(A_MATRIX,ierr)


  !call lis_output_matrix(A_MATRIX,2,matname,ierr);



  do local_index=1,n_size

    call lis_vector_set_value(LIS_INS_VALUE, local_index, LS_rhs(local_index), IMPLICIT_RHS, ierr)

  end do

  !call lis_output_vector(IMPLICIT_RHS,2,vecname1,ierr);




  !Creating solver
  !call lis_solve(A_MATRIX,IMPLICIT_RHS,x_solution,solver,ierr)
  call lis_solve_kernel(A_MATRIX,IMPLICIT_RHS,x_solution,solver,my_prec,ierr)



  do iii=1,n_size

    row_DC_index = mod(iii-1,4) + 1
    col_DC_index = (iii-1)/4 + 1

    call lis_vector_get_value(x_solution, iii, delta_cons(row_DC_index, col_DC_index), ierr)

  end do

  !call lis_output_vector(x_solution,2,vecname,ierr);


  end subroutine implicit_step

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine compute_collisional_jacobian(coll_tensor,sol,Res_coll)

  implicit none

  REAL(KIND=8), DIMENSION(4),intent(in) :: sol,Res_coll
  REAL(KIND=8), DIMENSION(4,4),intent(inout) :: coll_tensor
  REAL(KIND=8), DIMENSION(4) :: perturbed_cons, pert
  REAL(KIND=8), DIMENSION(N_EQ) :: Res_coll_pert
  real(kind=8) :: dx
  integer :: cons_var,ii

  Res_coll_pert = 0.0d0
  coll_tensor = 0.0d0

  call COMPUTE_PERTURBED_CONS(sol,PERTURBED_CONS)

    do CONS_VAR=1,4

      call apply_perturbation(sol,perturbed_cons,pert,cons_var)

      call ELECTRONS_COLLISIONAL_SOURCES(pert,Res_coll_pert)

      do ii = 1,4

        coll_tensor(ii,cons_var) = (Res_coll_pert(ii) - Res_Coll(ii)) / &
                    ( perturbed_cons(cons_var) )

      end do

  end do

  end subroutine compute_collisional_jacobian

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine compute_EM_jacobian(EM_tensor,sol,Res_EM,Ex,Ey,B)

  implicit none

  REAL(KIND=8), DIMENSION(4),intent(in) :: sol,Res_EM
  REAL(KIND=8), DIMENSION(4,4),intent(inout) :: EM_tensor
  REAL(KIND=8), DIMENSION(4) :: perturbed_cons, pert
  REAL(KIND=8), DIMENSION(N_EQ) :: Res_EM_pert
  real(kind=8),intent(in) :: ex,ey,B
  real(kind=8) :: dx
  integer :: cons_var,ii

  Res_EM_pert = 0.0d0
  EM_tensor = 0.0d0

  call COMPUTE_PERTURBED_CONS(sol,PERTURBED_CONS)

    do CONS_VAR=1,4

      call apply_perturbation(sol,perturbed_cons,pert,cons_var)

      call ELECTRONS_EM_SOURCES(pert, Ex, Ey, B, Res_EM_pert)

      do ii = 1,4

        EM_tensor(ii,cons_var) = (Res_EM_pert(ii) - Res_EM(ii)) / &
                    ( perturbed_cons(cons_var) )

      end do

  end do

  end subroutine compute_EM_jacobian

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine UPDATE_IMPLICIT_TERMS(Electrons_RHS, dx, Ex, Ey, B,cell_index,Q_el_star,Res_coll,Res_EM)

  implicit NONE

  REAL(KIND=8), DIMENSION(4),intent(inout) :: Electrons_RHS,Q_el_star
  real(kind=8), dimension(4, 4) :: Diagonal_Source_Tensor
  real(kind=8) :: nu_coll_elastic_electrons,nu_coll_ionizing,nu_coll_excitation
  real(kind=8), intent(in) :: dx,Ex,Ey,B
  integer, intent(in) :: cell_index
  integer :: diag_index,rhs_index,iii
  REAL(KIND=8), DIMENSION(N_EQ), INTENT(IN)  :: Res_Coll,Res_EM
  LIS_INTEGER :: ierr
  REAL(KIND=8), DIMENSION(4) :: perturbed_cons, pert
  real(kind=8) :: density_perturbed,Coeff_EM_Momentum,Coeff_EM_Ene
  real(kind=8), dimension(4,4) :: Collisional_perturbed_tensor,EM_perturbed_tensor


  Collisional_perturbed_tensor = 0.0d0

  ! Compute Jacobian relative to collisional part
  call compute_collisional_jacobian(Collisional_perturbed_tensor,Q_el_star,Res_coll)
  call compute_EM_jacobian(EM_perturbed_tensor,Q_el_star,Res_EM,Ex,Ey,B)


  ! Internal cells start with index 2, but matrix starts with index 1
  diag_index = cell_index - 1

  Diagonal_Source_Tensor = 0.0d0

  Diagonal_Source_Tensor = - Collisional_perturbed_tensor - EM_perturbed_tensor

  Diagonal_Source_Tensor(1,1) = Diagonal_Source_Tensor(1,1) + (2.0d0/2.0d0/dt + 1/d_tau_vect(cell_index))
  Diagonal_Source_Tensor(2,2) = Diagonal_Source_Tensor(2,2) + (2.0d0/2.0d0/dt + 1/d_tau_vect(cell_index))
  Diagonal_Source_Tensor(3,3) = Diagonal_Source_Tensor(3,3) + (2.0d0/2.0d0/dt + 1/d_tau_vect(cell_index))
  Diagonal_Source_Tensor(4,4) = Diagonal_Source_Tensor(4,4) + (2.0d0/2.0d0/dt + 1/d_tau_vect(cell_index))



  ! Updating the final tensor with the source Jacobian contribution of the cell_index cell

  call Insert_in_Final_Tensor(dx*Diagonal_Source_Tensor,diag_index,diag_index)

  ! Updating the RHS with the contribution of the cell_index cell

  RHS_INDEX = (diag_index-1)*4


  do iii=1,4

    LS_rhs(rhs_index + iii) = Electrons_RHS(iii) - dx/dt * (Q_el_star(iii)- Q_n(rhs_index + iii))


  end do

  do iii=1,4

    Electrons_RHS(iii) = - 1.0d0/dt * (Q_el_star(iii)- Q_n(rhs_index + iii))

  end do


  end subroutine UPDATE_IMPLICIT_TERMS


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE COMPUTE_PERTURBED_CONS(CONS_VECT,PERTURBED_CONS)

  implicit NONE
  real(kind=8), dimension(n_eq), intent(in) :: cons_vect
  real(kind=8), dimension(n_eq), intent(inout) :: perturbed_cons
  real(Kind=8) :: tolerance
  integer :: i
  tolerance = 1d-15

  perturbed_cons = 0.00d0


  ! Taken the value of the conservative variables on the left/right cell of the interface
  ! goal is to compute a perturbation of these conservative variables, but it shall be done for
  ! electrons only.

  do i=1,n_eq

    perturbed_cons(i) =  + 1d-4*max(abs(cons_vect(i)),abs(BC_electrons_left(i)), &
                        abs(BC_electrons_right(i)),tolerance) * sign(1.0d0,cons_vect(i))

  end do


  end SUBROUTINE COMPUTE_PERTURBED_CONS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine apply_perturbation(U,PERTURBATION,U_PERTURBED,INDEX)

  implicit NONE
  real(kind=8), dimension(N_EQ), intent(in) :: PERTURBATION,U
  real(kind=8), dimension(N_EQ), intent(inout) :: U_PERTURBED
  integer, INTENT(IN) :: INDEX

  U_PERTURBED = U

  U_PERTURBED(INDEX) = U_PERTURBED(INDEX) + PERTURBATION(INDEX)

  end SUBROUTINE apply_perturbation

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine Insert_in_Final_Tensor(Tensor_Flux,row_index,col_index)

  implicit NONE

  REAL(KIND=8), DIMENSION(N_eq,4),intent(in) ::Tensor_Flux
  integer, intent(in) :: row_index,col_index
  integer :: position,starting_index,i_col,i_row


  ! Finding the position in order to properly access the
  ! value_vector field

  if(row_index .eq. 1)THEN

    position = col_index - 1

  ELSE

    position = (2 + 3*(row_index - 2)) + (1 + (col_index - row_index))

  end if

  ! value_vector is initialized starting from 0-index, each block has 16 elements
  starting_index = position * 16

  do i_col = 1,4

    do i_row = 1,N_eq

      value_vector(starting_index) = value_vector(starting_index) + Tensor_Flux(i_row,i_col)

      starting_index = starting_index + 1

    end do

  end do

  end SUBROUTINE Insert_in_Final_Tensor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module IMPLICIT_ELECTRONS_UTILITIES
