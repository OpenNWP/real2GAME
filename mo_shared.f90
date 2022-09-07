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

end module mo_shared
