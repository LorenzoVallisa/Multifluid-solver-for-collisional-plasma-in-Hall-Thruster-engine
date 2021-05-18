!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE GRID

  ! The grid is structured like this: the first cell is a ghost cell (GC),
  ! so the first interface is at its right.
  !
  !      cell1     int1   cell2    int2   cell3   int3          intN     cellN
  ! ----(GC_1)------|------(C)------|------(C)-----|--- (...) ---|------(GC_N)----
  !               x_min                                        x_max
  !

  USE GLOBAL

  IMPLICIT NONE

  INTEGER :: ii
  REAL(KIND=8) :: L_start, L_end, L_GHOST

  CONTAINS

  SUBROUTINE CREATE_GRID

    implicit none

    REAL(kind=8), PARAMETER :: Pi = 3.1415927
    real(kind=8), allocatable, dimension(:) :: chebs
    real(kind=8), allocatable ,dimension(:) :: i_vect
    integer :: I,ii
    allocate(chebs(n_int),i_vect(n_int-2))


    chebs(1) = 1
    chebs(n_int) = - 1
    i_vect = (/ (I, I = 1, n_int-2) /)
    chebs(2:n_int-1) = cos((2*i_vect-1)/(2*(N_int-2))*Pi)

    ! Create interfaces (LINEAR SPACING)
    IF (GRID_TYPE .EQ. "uniform") THEN

      DO ii = 1, N_int
        x_int(ii) = x_min + (ii-1)*(x_max - x_min)/(N_cells-2) ! N_cells - 2 = # of internal cells
      END DO

    ELSE if (GRID_TYPE .EQ. "Chebyshev1") then

    x_int(1) = x_min
    x_int(n_int) = x_max

    DO ii = 2, N_int-1
      x_int(ii) = x_int(ii-1) + abs(chebs(ii)-chebs(ii-1))*(x_int(n_int)-x_int(1))/2
    END DO


    END IF


    ! Now create cell centers array and length of cells
    L_start = x_int(2) - x_int(1)           ! Length of first cell
    L_end   = x_int(N_int) - x_int(N_int-1) ! Length of last cell
    L_GHOST = (L_start + L_end)/2 ! We use this for both ghost cells, so they are equal.

    x_cc(1)          = x_min - L_GHOST/2 ! Left ghost cell
    x_cc(N_cells)    = x_max + L_GHOST/2 ! Right ghost cell

    L_cells(1)       = L_GHOST ! Length of fist ghost cell
    L_cells(N_cells) = L_GHOST ! Length of last ghost cell

    DO ii = 2, N_cells-1 ! Internal cells
      x_cc(ii)    = (x_int(ii) + x_int(ii-1))/2
      L_cells(ii) = (x_int(ii) - x_int(ii-1))
    END DO


    !write(*,*)" dx from l_cells: ",(l_cells(3)),"   -   ",(l_cells(11)),"   -   ",(l_cells(57))
    !write(*,*)" dx from x_cc ND: ",(x_cc(24)-x_cc(23)), "dx from x_cc D: ",(x_cc(24)-x_cc(23))*x_carac
    !write(*,*)" dx from x_int: ",(x_int(45)-x_int(44))
    !write(*,*) x_int(1),"   -   ",x_int(1)*x_carac
    !write(*,*) x_int(N_int),"   -   ",x_int(N_int)*x_carac

  END SUBROUTINE CREATE_GRID

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE LEFT_CELL_FROM_INT(INT_ID, C_ID_LEFT)
    ! From the interface with id "INT_ID", the subroutine computes the ID of the
    ! cell at its LEFT.

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: INT_ID    ! ID of interface
    INTEGER, INTENT(OUT) :: C_ID_LEFT ! ID of left cell

    C_ID_LEFT = INT_ID
  END SUBROUTINE

  SUBROUTINE RIGHT_CELL_FROM_INT(INT_ID, C_ID_RIGHT)
    ! From the interface with id "INT_ID", the subroutine computes the ID of the
    ! cell at its RIGHT.

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: INT_ID     ! ID of interface
    INTEGER, INTENT(OUT) :: C_ID_RIGHT ! ID of right cell

    C_ID_RIGHT = INT_ID + 1
  END SUBROUTINE

END MODULE GRID
