! This source file is part of real2GAME,which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/real2GAME

module mo_shared

  ! This file contains some definitions and helper functions.

  use netcdf
  
  implicit none
  
  ! setting the floating pointeger :: precision
  ! single precision
  integer,  parameter :: ps = 6
  integer,  parameter :: rs = 37
  
  ! double :: precision
  integer,  parameter :: pd = 12
  integer,  parameter :: rd = 37
  
  integer,  parameter :: sp = selected_real_kind(ps,rs)                 ! single precission
  integer,  parameter :: dp = selected_real_kind(pd,rd)                 ! double precission
  
  integer,  parameter :: wp = dp                                        ! working precission
  
  integer,  parameter :: n_layers_input = 12                            ! number of layers from the input model that are used
  integer,  parameter :: n_avg_points = 20                              ! number of points used for averaging from the source model (such as ICON) to the target model (such as GAME)
  integer,  parameter :: n_sst_points = 259200                          ! the number of points of the SST grid
  integer,  parameter :: n_points_per_layer_input_icon_global = 2949120 ! the number of points per layer of the global ICON model
  integer,  parameter :: n_points_per_layer_input_icon_d2 = 1084249     ! the number of points per layer of the regional ICON-D2 model
  real(wp), parameter :: M_PI = 4._wp*atan(1._wp)                       ! pi
  
  contains
  
  function calculate_distance_h(latitude_a,longitude_a,latitude_b,longitude_b,radius)
  
    ! This function returns the geodetic distance of two points given their geographical coordinates.
    
    real(wp), intent(in) :: latitude_a           ! latitude of the first point (rad)
    real(wp), intent(in) :: longitude_a          ! longitude of the first point (rad)
    real(wp), intent(in) :: latitude_b           ! latitude of the second point (rad)
    real(wp), intent(in) :: longitude_b          ! longitude of the second point (rad)
    real(wp), intent(in) :: radius               ! radius of the sphere to work with (m)
    real(wp)             :: calculate_distance_h ! result (m)
    
    ! computing the result
    calculate_distance_h = 2._wp*radius*asin(sqrt(0.5_wp-0.5_wp*(cos(latitude_a)*cos(latitude_b) &
    *cos(longitude_b-longitude_a)+sin(latitude_a)*sin(latitude_b))))
    
  end function calculate_distance_h

  function find_min_index(vector,vector_length)
  
    ! This function returns the index where a vector has its minimum.
    
    integer,  intent(in) :: vector_length         ! length of the vector of which the minimum has to be found
    real(wp), intent(in) :: vector(vector_length) ! vector of which to find the minimum
    integer              :: find_min_index        ! result
    
    ! local variables
    integer  :: ji          ! loop index
    real(wp) :: current_min ! the lowest value encountered so far
    
    ! initialization of the result and the current minimum
    find_min_index = 1
    current_min = vector(1)
    
    ! loop over all remaining elements of the vector
    do ji=2,vector_length
      if (vector(ji)<current_min) then
        current_min = vector(ji) 
        find_min_index = ji
      endif
    enddo
    
  end function find_min_index
  
  subroutine nc_check(i_status)
  
    ! This subroutine checks wether a NetCDF function threw an error.
  
    integer, intent(in) :: i_status ! status ID of a netCDF function

    if(i_status/=nf90_noerr) then 
      print *, trim(nf90_strerror(i_status))
      stop "Netcdf threw an error."
    end if
    
  end subroutine nc_check
  
  character(len=64) function int2string(input)
  
    ! This is a helper function which converts an integer to a string.
  
    integer, intent(in) :: input ! integer to convert to a string
    
    write(int2string,*) input
    int2string = adjustl(int2string)
    
  end function int2string

end module mo_shared





