! This source file is part of GAME-DA, which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME-DA

! This module contains linear algebra functions.

module lin_algebra

  use iso_c_binding

  implicit none

  ! precision of doubles
  integer, parameter        :: wp = c_double
  
  contains
    
  subroutine gauss(to_be_inverted, inv, matrix_size) &
  bind(c, name = "gauss")

    ! This subroutine computes the inverse inv of the matrix to_be_inverted, using the Gauss scheme.
    ! CAUTION: in the process, to_be_inverted will be modified.
    
    real(wp), intent(inout) :: to_be_inverted(matrix_size,matrix_size)
    real(wp), intent(inout) :: inv           (matrix_size,matrix_size)
    integer,  intent(in)    :: matrix_size

    ! local variables
    integer                 :: permute_index_found, permute_index_counter
    real(wp)                :: factor
    ! loop indices
    integer                 :: i, j, k

    ! firstly, the inverse is initialized with the unity matrix
    !$omp parallel
    !$omp do private(i)
    do i = 1,matrix_size
      inv(i,i) = 1._wp
    enddo
    !$omp end do
    !$omp end parallel

    ! Gaussian downwards
    ! ------------------
    ! we will start to modify to_be_inverted now (misuse of name)
    do i = 1,matrix_size-1
      ! checking if a permutation is necessary
      ! Firstly, the permutation index has to be found.
      permute_index_found = 0
      permute_index_counter = i
      do while (permute_index_found == 0)
        if (to_be_inverted(permute_index_counter,i) /= 0._wp) then
          permute_index_found = 1
        else
          permute_index_counter = permute_index_counter + 1
        endif
      enddo
      ! actually performing the permutation
      if (permute_index_counter > i) then
        call permute_lines(to_be_inverted, i, permute_index_counter, matrix_size)
        call permute_lines(inv, i, permute_index_counter, matrix_size)
      endif

      ! permutation is done, now comes the actual calculation
      ! dividing the line by to_be_inverted(i,i)
      factor = 1._wp/to_be_inverted(i,i)
      !$omp parallel
      !$omp do private(j)
      do j=i,matrix_size
        to_be_inverted(i,j) = factor*to_be_inverted(i,j)
      enddo
      !$omp end do
      !$omp end parallel
      !$omp parallel
      !$omp do private(j)
      do j = 1,matrix_size
        inv(i,j) = factor*inv(i,j)
      enddo
      !$omp end do
      !$omp end parallel
      ! loop over all the lines that are below the current line
      !$omp parallel
      !$omp do private(factor,j,k)
      do j=i+1,matrix_size
        factor = -to_be_inverted(j,i)
        do k=i,matrix_size
          to_be_inverted(j,k) = to_be_inverted(j,k) + factor*to_be_inverted(i,k)
        enddo
        inv(j,:) = inv(j,:) + factor*inv(i,:)
      enddo
      !$omp end do
      !$omp end parallel
    enddo

    !$omp parallel
    !$omp do private(j)
    do j = 1,matrix_size
      inv(matrix_size,j) = inv(matrix_size,j)/to_be_inverted(matrix_size,matrix_size)
    enddo
    !$omp end do
    !$omp end parallel
    to_be_inverted(matrix_size,matrix_size) = 1._wp

    ! Gaussian upwards
    ! ----------------
    do i = matrix_size,2,-1
      !$omp parallel
      !$omp do private(j)
      do j = i-1,1,-1
        inv(j,:) = inv(j,:) - to_be_inverted(j,i)*inv(i,:)
      enddo
      !$omp end do
      !$omp end parallel
    enddo
  
  end subroutine gauss

  subroutine permute_lines(matrix, line_a, line_b, matrix_size)

    ! This subroutine permutes line_a with line_b of matrix.

    ! arguments
    real(wp), intent(inout) :: matrix(matrix_size,matrix_size)
    integer,  intent(in)    :: line_a
    integer,  intent(in)    :: line_b
    integer,  intent(in)    :: matrix_size

    ! local variables
    real(wp)                :: line_a_pre(matrix_size)

    line_a_pre(:) = matrix(line_a,:)
    matrix(line_a,:) = matrix(line_b,:)
    matrix(line_b,:) = line_a_pre(:)

  end subroutine permute_lines

end module lin_algebra








