! This source file is part of real2GAME, which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/real2GAME

module this_does_something

  ! This module contains some spatial operators.
  
  use definitions, only: wp
  
  implicit none
  
  contains

  subroutine useful(in_field,out_field)
  
    ! This subroutine does something useful.
    
    real(wp), intent(in)  :: in_field(10)
    real(wp), intent(out) :: out_field(10)
    
    ! local variables
    integer :: ji
    
    !$omp parallel do private(ji)
    do ji=1,10
      out_field(ji) = in_field(ji) + 1._wp
    enddo
    !$omp end parallel do
    
  end subroutine useful

end module this_does_something

















