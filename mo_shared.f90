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
  
  ! real(wp) :: precision
  integer,  parameter :: pd = 12
  integer,  parameter :: rd = 37
  
  integer,  parameter :: sp = selected_real_kind(ps,rs)     ! single precission
  integer,  parameter :: dp = selected_real_kind(pd,rd)     ! real(wp) :: precission
  
  integer,  parameter :: wp = dp                            ! working precission
  
  integer,  parameter :: n_layers_input = 12                ! number of levels from the input model that are used
  integer,  parameter :: n_avg_points = 20                  ! number of points used for averaging
  integer,  parameter :: n_sst_points = 259200              ! the number of points of the SST grid
  integer,  parameter :: n_points_per_layer_input = 2949120 ! the number of points per layer of the input model
  real(wp), parameter :: M_PI = 4._wp*atan(1._wp)           ! pi
  
  contains
  
  function calculate_distance_h(latitude_a,longitude_a,latitude_b,longitude_b,radius)
  
    ! This function returns the geodetic distance of two points given their geographical coordinates.
    
    real(wp), intent(in) :: latitude_a,longitude_a,latitude_b,longitude_b,radius
    real(wp)             :: calculate_distance_h
    
    calculate_distance_h = 2.0*radius*asin(sqrt(0.5-0.5*(cos(latitude_a)*cos(latitude_b) &
    *cos(longitude_b-longitude_a)+sin(latitude_a)*sin(latitude_b))))
    
  end function calculate_distance_h

  function find_min_index(vector,vector_length)
  
    ! This function returns the index where a vector has its minimum.
    
    integer,  intent(in) :: vector_length
    real(wp), intent(in) :: vector(vector_length)
    integer              :: find_min_index
    
    ! local variables
    integer   :: ji
    real(wp) :: current_min
    
    find_min_index = 1
    current_min = vector(1)
    
    do ji=2,vector_length
      if (vector(ji)<current_min) then
        current_min = vector(ji) 
        find_min_index = ji
      endif
    enddo
    
  end function find_min_index
  
  subroutine nc_check(i_status)
  
    ! This checks wether a NetCDF function threw an error.
  
    integer, intent(in) :: i_status

    if(i_status/=nf90_noerr) then 
      print *, trim(nf90_strerror(i_status))
      stop "Netcdf threw an error."
    end if
    
  end subroutine nc_check
  
  character(len=64) function int2string(input)
  
    ! This is a helper function which converts an integer to a string.
  
    integer, intent(in) :: input
    
    write(int2string,*) input
    int2string = adjustl(int2string)
    
  end function int2string

end module mo_shared





