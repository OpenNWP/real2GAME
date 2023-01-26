! This source file is part of real2GAME, which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/real2GAME

program control
  
  ! This tool reads the output from other models / data assimilation systems and brings it into a standardized format.
  
  use netcdf
  use eccodes
  use mo_shared, only: wp,n_layers_input,n_sst_points,n_points_per_layer_input_icon_global,nc_check,int2string,M_PI, &
                       n_points_per_layer_input_icon_d2
  
  implicit none
  
  integer               :: jl                            ! layer index
  integer               :: ncid                          ! ID of the netCDF file
  integer               :: h_dimid                       ! netCDF ID of the horizontal dimension
  integer               :: v_dimid                       ! netCDF ID of the vertical dimension
  integer               :: z_surf_id                     ! netCDF ID of the surface height
  integer               :: sp_id                         ! netCDF ID of the surface pressure
  integer               :: sst_dimid                     ! netCDF ID of the dimension of the sea surface temperature grid
  integer               :: t_id                          ! netCDF ID of the temperature
  integer               :: spec_hum_id                   ! netCDF ID of the specific humidity
  integer               :: z_id                          ! netCDF ID of the z-coordinates
  integer               :: u_id                          ! netCDF ID of the zonal wind
  integer               :: v_id                          ! netCDF ID of the meridional wind
  integer               :: sst_id                        ! netCDF ID of the sea surface temperature
  integer               :: dim_vector(2)                 ! vector of two netCDF dimensions
  integer               :: layers_vector(n_layers_input) ! vector containing the layers of the input system used for the initialization
  integer               :: jfile                         ! ID of a grib file
  integer               :: jgrib                         ! ID of a grib message
  integer               :: n_points_per_layer_input      ! number of points per layer of the input system
  real(wp), allocatable :: z_height_amsl_one_layer(:)    ! vertical height above mean sea level on one layer
  real(wp), allocatable :: temperature_one_layer(:)      ! temperature on one layer
  real(wp), allocatable :: spec_hum_one_layer(:)         ! specific humidity on one layer
  real(wp), allocatable :: u_one_layer(:)                ! zonal wind on one layer
  real(wp), allocatable :: v_one_layer(:)                ! meridional wind on one layer
  real(wp), allocatable :: z_height_amsl(:,:)            ! vertical height above mean sea level
  real(wp), allocatable :: temperature(:,:)              ! temperature
  real(wp), allocatable :: spec_hum(:,:)                 ! specific humidity
  real(wp), allocatable :: u_wind(:,:)                   ! zonal wind
  real(wp), allocatable :: v_wind(:,:)                   ! meridional wind
  real(wp), allocatable :: sst(:)                        ! sea surface temperature
  real(wp), allocatable :: surface_height(:)             ! surface height
  real(wp), allocatable :: pressure_surface(:)           ! surface pressure
  character(len=4)      :: year_string                   ! year of the initialization time as a string (command line argument)
  character(len=2)      :: month_string                  ! month of the initialization time as a string (command line argument)
  character(len=2)      :: day_string                    ! day of the initialization time as a string (command line argument)
  character(len=2)      :: hour_string                   ! hour of the initialization time as a string (command line argument)
  character(len=128)    :: real2game_root_dir            ! root directory of real2GAME (command line argument)
  character(len=256)    :: z_input_model_file            ! name of file containing the z-coordinates of the input system
  character(len=256)    :: temperature_file              ! name of file containing the temperature of the input system
  character(len=256)    :: u_wind_file                   ! name of file containing the zonal wind of the input system
  character(len=256)    :: v_wind_file                   ! name of file containing the meridional wind of the input system
  character(len=256)    :: sfc_height_file               ! name of file containing the surface height of the input system
  character(len=256)    :: sfc_pres_file                 ! name of file containing the surface pressure of the input system
  character(len=256)    :: sst_file                      ! name of file containing the sea surface temperature of the input system
  character(len=256)    :: spec_hum_file                 ! name of file containing the specific humidity of the input system
  character(len=256)    :: output_file                   ! name of the output file
  
  ! defining the levels of the model we want to use
  layers_vector(1) = 20
  layers_vector(2) = 30
  layers_vector(3) = 40
  layers_vector(4) = 50
  layers_vector(5) = 60
  layers_vector(6) = 70
  layers_vector(7) = 80
  layers_vector(8) = 90
  layers_vector(9) = 100
  layers_vector(10) = 110
  layers_vector(11) = 120
  
  ! shell arguments
  call get_command_argument(1,year_string)
  call get_command_argument(2,month_string)
  call get_command_argument(3,day_string)
  call get_command_argument(4,hour_string)
  call get_command_argument(5,real2game_root_dir)
  
  n_points_per_layer_input = n_points_per_layer_input_icon_global
  
  ! single-layer arrays from grib
  allocate(z_height_amsl_one_layer(n_points_per_layer_input))
  allocate(temperature_one_layer(n_points_per_layer_input))
  allocate(spec_hum_one_layer(n_points_per_layer_input))
  allocate(u_one_layer(n_points_per_layer_input))
  allocate(v_one_layer(n_points_per_layer_input))
  
  !$omp parallel workshare
  z_height_amsl_one_layer = 0._wp
  temperature_one_layer = 0._wp
  spec_hum_one_layer = 0._wp
  u_one_layer = 0._wp
  v_one_layer = 0._wp
  !$omp end parallel workshare
  
  ! 2D-arrays for netCDF
  allocate(z_height_amsl(n_points_per_layer_input,n_layers_input))
  allocate(temperature(n_points_per_layer_input,n_layers_input))
  allocate(spec_hum(n_points_per_layer_input,n_layers_input))
  allocate(u_wind(n_points_per_layer_input,n_layers_input))
  allocate(v_wind(n_points_per_layer_input,n_layers_input))
  
  !$omp parallel workshare
  z_height_amsl = 0._wp
  temperature = 0._wp
  spec_hum = 0._wp
  u_wind = 0._wp
  v_wind = 0._wp
  !$omp end parallel workshare
  
  ! grib stuff
  
  ! reading the data from the free atmosphere
  ! loop over all relevant levels in the free atmosphere
  do jl=1,n_layers_input
    ! vertical position of the current layer
    z_input_model_file = trim(real2game_root_dir) // "/input/icon_global_icosahedral_time-invariant_" // year_string &
                         // month_string // day_string // hour_string // "_" &
                         // trim(int2string(layers_vector(jl))) // "_HHL.grib2"
    
    call codes_open_file(jfile,trim(z_input_model_file),"r")
    call codes_grib_new_from_file(jfile,jgrib)
    call codes_get(jgrib,"values",z_height_amsl_one_layer)
    call codes_release(jgrib)
    call codes_close_file(jfile)
    
    !$omp parallel workshare
    z_height_amsl(:,jl) = z_height_amsl_one_layer
    !$omp end parallel workshare
    
    ! reading the temperature
    temperature_file = trim(real2game_root_dir) // "/input/icon_global_icosahedral_model-level_" // year_string // month_string &
                       // day_string // hour_string // "_000_" // trim(int2string(layers_vector(jl))) // "_T.grib2"
    
    call codes_open_file(jfile,trim(temperature_file),"r")
    call codes_grib_new_from_file(jfile,jgrib)
    call codes_get(jgrib,"values",temperature_one_layer)
    call codes_release(jgrib)
    call codes_close_file(jfile)
    
    !$omp parallel workshare
    temperature(:,jl) = temperature_one_layer
    !$omp end parallel workshare
    
    ! reading the specific humidity
    spec_hum_file = trim(real2game_root_dir) // "/input/icon_global_icosahedral_model-level_" // year_string // month_string // &
                    day_string // hour_string // "_000_" // trim(int2string(layers_vector(jl))) // &
                    "_QV.grib2"
    
    call codes_open_file(jfile,trim(spec_hum_file),"r")
    call codes_grib_new_from_file(jfile,jgrib)
    call codes_get(jgrib,"values",spec_hum_one_layer)
    call codes_release(jgrib)
    call codes_close_file(jfile)
    
    !$omp parallel workshare
    spec_hum(:,jl) = spec_hum_one_layer
    !$omp end parallel workshare
    
    ! reading the u wind
    u_wind_file = trim(real2game_root_dir) // "/input/icon_global_icosahedral_model-level_" // year_string // month_string // &
                  day_string // hour_string // "_000_" // trim(int2string(layers_vector(jl))) // & 
                  "_U.grib2"
    
    call codes_open_file(jfile,trim(u_wind_file),"r")
    call codes_grib_new_from_file(jfile,jgrib)
    call codes_get(jgrib,"values",u_one_layer)
    call codes_release(jgrib)
    call codes_close_file(jfile)
    
    !$omp parallel workshare
    u_wind(:,jl) = u_one_layer
    !$omp end parallel workshare
    
    ! reading the v wind
    v_wind_file = trim(real2game_root_dir) // "/input/icon_global_icosahedral_model-level_" // year_string // &
                  month_string // day_string // hour_string // "_000_" // trim(int2string(layers_vector(jl))) &
                  // "_V.grib2"
    
    call codes_open_file(jfile,trim(v_wind_file),"r")
    call codes_grib_new_from_file(jfile,jgrib)
    call codes_get(jgrib,"values",v_one_layer)
    call codes_release(jgrib)
    call codes_close_file(jfile)
    
    !$omp parallel workshare
    v_wind(:,jl) = v_one_layer
    !$omp end parallel workshare
    
  enddo
  
  deallocate(z_height_amsl_one_layer)
  deallocate(temperature_one_layer)
  deallocate(spec_hum_one_layer)
  deallocate(u_one_layer)
  deallocate(v_one_layer)
  
  ! reading the surface height
  allocate(surface_height(n_points_per_layer_input))
  !$omp parallel workshare
  surface_height = 0._wp
  !$omp end parallel workshare
  sfc_height_file = trim(real2game_root_dir) // "/input/icon_global_icosahedral_time-invariant_" // year_string // month_string // &
                    day_string // hour_string // "_HSURF.grib2"
  
  call codes_open_file(jfile,trim(sfc_height_file),"r")
  call codes_grib_new_from_file(jfile,jgrib)
  call codes_get(jgrib,"values",surface_height)
  call codes_release(jgrib)
  call codes_close_file(jfile)
  
  ! reading the surface presure
  allocate(pressure_surface(n_points_per_layer_input))
  !$omp parallel workshare
  pressure_surface = 0._wp
  !$omp end parallel workshare
  
  sfc_pres_file = trim(real2game_root_dir) // "/input/icon_global_icosahedral_single-level_" // year_string // month_string // &
                  day_string // hour_string // "_000_PS.grib2"
  
  call codes_open_file(jfile,trim(sfc_pres_file),"r")
  call codes_grib_new_from_file(jfile,jgrib)
  call codes_get(jgrib,"values",pressure_surface)
  call codes_release(jgrib)
  call codes_close_file(jfile)
  
  ! reading the SST
  allocate(sst(n_sst_points))
  !$omp parallel workshare
  sst = 0._wp
  !$omp end parallel workshare
  
  sst_file = trim(real2game_root_dir) // "/input/rtgssthr_grb_0.5.grib2"
  
  call codes_open_file(jfile,trim(sst_file),"r")
  call codes_grib_new_from_file(jfile,jgrib)
  call codes_get(jgrib,"values",sst)
  call codes_release(jgrib)
  call codes_close_file(jfile)
  
  ! writing the observations to a netCDF file
  output_file = trim(real2game_root_dir) // "/input/obs_" // year_string // month_string // day_string // &
                hour_string // ".nc"
  
  call nc_check(nf90_create(trim(output_file),NF90_CLOBBER,ncid))
  
  ! defining the dimensions
  call nc_check(nf90_def_dim(ncid,"h_index",n_points_per_layer_input,h_dimid))
  call nc_check(nf90_def_dim(ncid,"v_index",n_layers_input,v_dimid))
  call nc_check(nf90_def_dim(ncid,"sst_index",n_sst_points,sst_dimid))
  dim_vector(1) = h_dimid
  dim_vector(2) = v_dimid
  call nc_check(nf90_def_var(ncid,"z_surface",NF90_REAL,h_dimid,z_surf_id))
  call nc_check(nf90_def_var(ncid,"pressure_surface",NF90_REAL,h_dimid,sp_id))
  call nc_check(nf90_def_var(ncid,"z_height",NF90_REAL,dim_vector,z_id))
  call nc_check(nf90_def_var(ncid,"temperature",NF90_REAL,dim_vector,t_id))
  call nc_check(nf90_def_var(ncid,"spec_humidity",NF90_REAL,dim_vector,spec_hum_id))
  call nc_check(nf90_def_var(ncid,"u_wind",NF90_REAL,dim_vector,u_id))
  call nc_check(nf90_def_var(ncid,"v_wind",NF90_REAL,dim_vector,v_id))
  call nc_check(nf90_def_var(ncid,"sst",NF90_REAL,sst_dimid,sst_id))
  call nc_check(nf90_enddef(ncid))
  
  ! setting the variables
  call nc_check(nf90_put_var(ncid,z_surf_id,surface_height))
  call nc_check(nf90_put_var(ncid,sp_id,pressure_surface))
  call nc_check(nf90_put_var(ncid,z_id,z_height_amsl))
  call nc_check(nf90_put_var(ncid,t_id,temperature))
  call nc_check(nf90_put_var(ncid,spec_hum_id,spec_hum))
  call nc_check(nf90_put_var(ncid,u_id,u_wind))
  call nc_check(nf90_put_var(ncid,v_id,v_wind))
  call nc_check(nf90_put_var(ncid,sst_id,sst))
  call nc_check(nf90_close(ncid))
  
  ! freeing the memory
  deallocate(z_height_amsl)
  deallocate(surface_height)
  deallocate(pressure_surface)
  deallocate(u_wind)
  deallocate(v_wind)
  deallocate(temperature)
  deallocate(spec_hum)
  deallocate(sst)
  
end program control










