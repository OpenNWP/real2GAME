! This source file is part of real2GAME, which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/real2GAME

module index_helpers
  
  ! This module contains helper functions concerned with simple algebraic operations on vectors.

  use iso_c_binding
  
  implicit none
  
  public :: calculate_distance_h
  public :: find_min_index
  
  contains
  
  function calculate_distance_h(latitude_a,longitude_a,latitude_b,longitude_b,radius) &
  bind(c,name = "calculate_distance_h")
  
    ! This function returns the geodetic distance of two points given their geographical coordinates.
    
    real(c_double), intent(in) :: latitude_a,longitude_a,latitude_b,longitude_b,radius
    real(c_double)             :: calculate_distance_h
    
    calculate_distance_h = 2.0*radius*asin(sqrt(0.5-0.5*(cos(latitude_a)*cos(latitude_b) &
    *cos(longitude_b-longitude_a)+sin(latitude_a)*sin(latitude_b))))
    
  end function calculate_distance_h

  function find_min_index(vector,vector_length) &
  bind(c,name = "find_min_index")
  
    ! This function returns the index where a vector has its minimum.
    
    integer(c_int), intent(in) :: vector_length
    real(c_double), intent(in) :: vector(vector_length)
    integer(c_int)             :: find_min_index
    
    ! local variables
    integer        :: ji
    real(c_double) :: current_min
    
    find_min_index = 1
    current_min = vector(1)
    
    do ji=2,vector_length
      if (vector(ji)<current_min) then
        current_min = vector(ji) 
        find_min_index = ji
      endif
    enddo
    
    find_min_index = find_min_index - 1
    
  end function find_min_index

end module index_helpers










