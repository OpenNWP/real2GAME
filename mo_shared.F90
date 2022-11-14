! This source file is part of real2GAME,which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/real2GAME

module mo_shared

  ! This file contains some definitions and helper functions.

  use netcdf
  
  implicit none
  
  ! setting the floating point precision
  ! single precision
  integer,  parameter :: ps = 6                                         ! single decimal precision
  integer,  parameter :: rs = 37                                        ! single exponent precision
  ! double precision
  integer,  parameter :: pd = 12                                        ! double decimal precision
  integer,  parameter :: rd = 37                                        ! double exponent precision
  
  integer,  parameter :: sp = selected_real_kind(ps,rs)                 ! single precision
  integer,  parameter :: dp = selected_real_kind(pd,rd)                 ! double precision
  
  integer,  parameter :: wp = dp                                        ! working precision
  
  integer,  parameter :: n_layers_input = 12                            ! number of layers from the input model that are used
  integer,  parameter :: n_avg_points = 20                              ! number of points used for averaging from the source model (such as ICON) to the target model (such as GAME)
  integer,  parameter :: n_sst_points = 259200                          ! the number of points of the SST grid
  integer,  parameter :: n_points_per_layer_input_icon_global = 2949120 ! the number of points per layer of the global ICON model
  integer,  parameter :: n_points_per_layer_input_icon_d2 = 1084249     ! the number of points per layer of the regional ICON-D2 model
  ! constants
  real(wp), parameter :: M_PI = 4._wp*atan(1._wp)                       ! pi
  real(wp), parameter :: n_a =  6.02214076e23_wp                        ! Avogadro's number
  real(wp), parameter :: p_0 = 100000._wp                               ! reference pressure
  real(wp), parameter :: m_d = n_a*0.004810e-23_wp                      ! molar mass of dry air
  real(wp), parameter :: m_v = n_a*0.002991e-23_wp                      ! molar mass of water
  real(wp), parameter :: r_d = 287.057811_wp                            ! specific gas constant of dry air
  real(wp), parameter :: c_d_p = 1005._wp                               ! isobaric specific heat capacity of dry air
  real(wp), parameter :: scale_height = 8000._wp                        ! scale_height
  real(wp), parameter :: t_0 = 273.15_wp                                ! 273.15 K
  real(wp), parameter :: r_v = 461.524879_wp                            ! specific gas constant of water vapour
  
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

  function find_min_index(vector)
  
    ! This function returns the index where a vector has its minimum.
    
    real(wp), intent(in) :: vector(:)      ! vector of which to find the minimum
    integer              :: find_min_index ! result
    
    ! local variables
    integer  :: ji          ! loop index
    real(wp) :: current_min ! the lowest value encountered so far
    
    ! initialization of the result and the current minimum
    find_min_index = 1
    current_min = vector(1)
    
    ! loop over all remaining elements of the vector
    do ji=2,size(vector)
      if (vector(ji)<current_min) then
        current_min = vector(ji) 
        find_min_index = ji
      endif
    enddo
    
  end function find_min_index
  
  subroutine nc_check(i_status)
  
    ! This subroutine checks wether a netCDF function threw an error.
  
    integer, intent(in) :: i_status ! status ID of a netCDF function

    if(i_status/=nf90_noerr) then 
      print *, trim(nf90_strerror(i_status))
      stop "netCDF threw an error."
    end if
    
  end subroutine nc_check
  
  character(len=64) function int2string(input)
  
    ! This is a helper function which converts an integer to a string.
  
    integer, intent(in) :: input ! integer to convert to a string
    
    write(int2string,*) input
    int2string = adjustl(int2string)
    
  end function int2string
  
  function rel_humidity(abs_humidity,temperature)
    
    ! This function returns the relative humidity as a function of the absolute humidity in kg/m^3 and the temperature in K.
    
    real(wp), intent(in) :: abs_humidity ! absolute humidity (water vapour mass density (kg/m**3))
    real(wp), intent(in) :: temperature  ! temperature (K)
    real(wp)             :: rel_humidity ! result
    
    ! local variables
    real(wp)             :: vapour_pressure     ! actual water vapour pressure
    real(wp)             :: saturation_pressure ! saturation water vapour pressure
    
    ! calculation of the water vapour pressure according to the equation of state
    vapour_pressure = abs_humidity*r_v*temperature
    
    if (temperature>t_0) then
      saturation_pressure = saturation_pressure_over_water(temperature)
    endif
    if (temperature<=t_0) then
      saturation_pressure = saturation_pressure_over_ice(temperature)
    endif
    
    rel_humidity = vapour_pressure/saturation_pressure
    
  end function rel_humidity

  function saturation_pressure_over_water(temperature)

    ! This function returns the saturation pressure in Pa over liquid water as a function of the temperature in K.
    ! It uses the formula by Huang: A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice, 2018, DOI: 10.1175/JAMC-D-17-0334.1.
    
    real(wp), intent(in) :: temperature                    ! temperature
    real(wp)             :: saturation_pressure_over_water ! result
    
    ! local variables
    real(wp)  :: temp_c ! temperature in degrees Celsius

    temp_c = temperature - t_0
    ! clipping too extreme values for this approximation
    if (temp_c>100._wp) then
      temp_c = 100._wp
    endif
    
    if (temp_c>0._wp) then
      saturation_pressure_over_water = exp(34.494_wp-4924.99_wp/(temp_c+237.1_wp))/(temp_c+105._wp)**1.57_wp
    ! For super-cooled water we use the formula cited in Pruppacher and Klett (2010), p. 854, Eq. (A.4-1).
    else
      ! Clipping values that are too extreme for this approximation.
      if (temp_c<-50._wp) then
        temp_c = -50._wp
      endif
      saturation_pressure_over_water &
      = 6.107799961_wp &
      + 4.436518521e-1_wp*temp_c &
      + 1.428945805e-2_wp*temp_c**2 &
      + 2.650648471e-4_wp*temp_c**3 &
      + 3.031240396e-6_wp*temp_c**4 &
      + 2.034080948e-8_wp*temp_c**5 &
      + 6.136820929e-11_wp*temp_c**6
    endif

  end function saturation_pressure_over_water
  
  function saturation_pressure_over_ice(temperature)
    
    ! This function returns the saturation pressure in Pa over ice as a function of the temperature in K.
    ! It blends the two formulas of Huang and Murphy.
    
    real(wp), intent(in) :: temperature                  ! temperature
    real(wp)             :: saturation_pressure_over_ice ! result
    
    ! local variables
    real(wp) :: t_local      ! local copy of temperature
    real(wp) :: temp_c       ! temperature in degrees Celsius
    real(wp) :: huang_weight ! weight of the Huang formula

    t_local = temperature

    temp_c = t_local - t_0
    
    if (temp_c>=-80._wp) then
      ! at temperatures > 0 degrees Celsius ice cannot exist in equilibrium which is why this is clipped
      if (t_local>t_0) then
        t_local = t_0
      endif
      saturation_pressure_over_ice = saturation_pressure_ice_huang(t_local)
    elseif (temp_c>=-100._wp) then
      huang_weight = (temp_c+100._wp)/20._wp
      saturation_pressure_over_ice = huang_weight*saturation_pressure_ice_huang(t_local) &
                                     + (1._wp-huang_weight)+saturation_pressure_ice_murphy(t_local)
    else
      ! clipping too extreme values for this approximation
      if (t_local<110._wp) then
        t_local = 110._wp
      endif
      saturation_pressure_over_ice = saturation_pressure_ice_murphy(t_local)
    endif
    
  end function saturation_pressure_over_ice
  
  function saturation_pressure_ice_huang(temperature)
  
    ! This function computes the saturation pressure over ice.
    ! It follows the formula by Huang: A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice, 2018, DOI: 10.1175/JAMC-D-17-0334.1.
    
    real(wp), intent(in) :: temperature                   ! temperature
    real(wp)             :: saturation_pressure_ice_huang ! result
    
    ! local variables
    real(wp) :: temp_c ! temperature in degrees Celsius
    
    temp_c = temperature - t_0
    
    saturation_pressure_ice_huang = exp(43.494_wp-6545.8_wp/(temp_c+278._wp))/(temp_c+868._wp)**2
  
  end function saturation_pressure_ice_huang
  
  function saturation_pressure_ice_murphy(temperature)
  
    ! This function computes the saturation pressure over ice.
    ! It follows Eq. (7) in Murphy DM, Koop T. Review of the vapour pressures of ice and supercooled water for atmospheric applications.
    ! QUARTERLY JOURNAL OF THE ROYAL METEOROLOGICAL SOCIETY. 2005;131(608):1539-1565.

    real(wp), intent(in) :: temperature                    ! temperature
    real(wp)             :: saturation_pressure_ice_murphy ! result
    
    ! computing the result
    saturation_pressure_ice_murphy = exp(9.550426_wp-5723.265_wp/temperature+3.53068_wp*log(temperature)-0.00728332_wp*temperature)
  
  end function saturation_pressure_ice_murphy

end module mo_shared





