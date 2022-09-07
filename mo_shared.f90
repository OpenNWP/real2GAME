! This source file is part of real2GAME,which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/real2GAME

module mo_shared

  ! This file contains some definitions and helper functions.
  
  implicit none
  
  ! setting the floating pointeger :: precision
  ! single precision
  integer, parameter :: ps = 6
  integer, parameter :: rs = 37
  
  ! real(wp) :: precision
  integer, parameter :: pd = 12
  integer, parameter :: rd = 37
  
  integer, parameter :: sp = selected_real_kind(ps,rs) ! single precission
  integer, parameter :: dp = selected_real_kind(pd,rd) ! real(wp) :: precission
  
  integer, parameter :: wp = dp                        ! working precission
  
  integer, parameter :: n_levels_input = 12                ! number of levels from the input model that are used
  integer, parameter :: n_sst_points = 259200              ! the number of points of the SST grid
  integer, parameter :: N_POINTS_PER_LAYER_INPUT = 2949120 ! the number of points per layer of the input model

end module mo_shared
