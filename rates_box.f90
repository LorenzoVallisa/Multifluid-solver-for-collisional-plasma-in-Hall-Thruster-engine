MODULE rates_box

USE GLOBAL

implicit none

REAL(KIND=8), allocatable, DIMENSION(:,:) :: T_K_elastic,T_K_elastic_ions
real(kind=8), allocatable, DIMENSION(:,:) :: T_K_ionized,T_K_excited


CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_elastic_k_rate(elastic_k,Temperature)

implicit NONE

REAL(KIND=8), INTENT(INOUT) :: elastic_k
REAL(KIND=8), INTENT(IN) :: temperature

call get_da_rate_1(elastic_k,T_K_elastic,Temperature)

end subroutine get_elastic_k_rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_ionized_k_rate(ionized_k,Temperature)

implicit NONE

REAL(KIND=8), INTENT(INOUT) :: ionized_k
REAL(KIND=8), INTENT(IN) :: temperature


call get_da_rate_1(ionized_k,T_K_ionized,Temperature)

end subroutine get_ionized_k_rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_excited_k_rate(excited_k,Temperature)

implicit NONE

REAL(KIND=8), INTENT(INOUT) :: excited_k
REAL(KIND=8), INTENT(IN) :: temperature


call get_da_rate_1(excited_k,T_K_excited,Temperature)

end subroutine get_excited_k_rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_ions_elastic_k_rate(ions_elastic_k,Temperature)

implicit NONE

REAL(KIND=8), INTENT(INOUT) :: ions_elastic_k
REAL(KIND=8), INTENT(IN) :: temperature

call get_da_rate_1(ions_elastic_k,T_K_elastic_ions,Temperature)

end subroutine get_ions_elastic_k_rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_da_rate_1(local_k,local_T_K,Temperature)


implicit NONE

REAL(KIND=8), INTENT(INOUT) :: local_k
real(kind=8), dimension(T_samples,2), intent(in) :: local_T_K
REAL(KIND=8), INTENT(IN)  :: Temperature
integer :: ID_Left,ID_Right


! If fluid temperature higher than the one sampled than assessing the maximum reaction rate

if (Temperature .ge. local_T_K(T_samples,1)) THEN

  local_k = maxval(local_T_K(:,2))

else if (local_T_K(1,1) .gt. Temperature ) THEN

  local_k = minval(local_T_K(:,2))

else

  ID_Left  = floor((Temperature - local_T_K(1,1))/(local_T_K(T_samples,1)-local_T_K(1,1))*(T_samples-1)) + 1
  ID_Right = ID_Left + 1

  local_k = local_T_K(ID_Left,2) + (Temperature - local_T_K(ID_Left,1))/(local_T_k(ID_Right,1) - &
                    local_T_K(ID_Left,1))*(local_T_k(ID_Right,2)-local_T_k(ID_Left,2))
end if


if (local_k .gt. maxval(local_T_K(:,2))) then

  write(*,*) "Reaction rate out of bounds: ", local_k, "greater than ", maxval(local_T_K(:,2))
  write(*,*) " Temperature: ", Temperature ," K "
  write(*,*) local_T_k(ID_Left,1),"   -   ",local_T_k(ID_Right,1)
  write(*,*)(Temperature - local_T_K(1,1))
  write(*,*)(local_T_K(T_samples,1)-local_T_K(1,1))*T_samples
  STOP

end if

end subroutine get_da_rate_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



END MODULE rates_box
