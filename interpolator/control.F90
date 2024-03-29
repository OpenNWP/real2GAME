! This source file is part of real2GAME, which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/real2GAME

program control
  
  ! This program coordinates the data interpolation process.
  
  use netcdf
  use mo_shared, only: wp,n_layers_input,n_sst_points,n_points_per_layer_input_icon_global,n_avg_points, &
                       nc_check,int2string,find_min_index,rel_humidity,M_PI,n_a,p_0,m_d,m_v,r_d,c_d_p, &
                       scale_height,t_0
  
  implicit none
  
  logical               :: ltke_avail                         ! switch indicating if TKE (turbulent kinetic energy) is present in the background state
  logical               :: lt_soil_avail                      ! switch indicating if soil temperature is present in the background state
  logical               :: lno_hydrometeors_found             ! switch indicating if hydrometeors are missing from the background state
  integer               :: ji                                 ! cell index
  integer               :: jl                                 ! layer index
  integer               :: jm                                 ! interpolation index
  integer               :: jn                                 ! another layer index
  integer               :: latitudes_game_id                  ! netCDF ID of the latitudes of the cell centers of GAME
  integer               :: longitudes_game_id                 ! netCDF ID of the longitudes of the cell centers of GAME
  integer               :: z_game_id                          ! netCDF ID of the z-coordinates of the cell centers of GAME
  integer               :: z_game_wind_id                     ! netCDF ID of the z-coordinates of the edges of GAME
  integer               :: gravity_potential_game_id          ! netCDF ID of the gravity potential of GAME
  integer               :: constituent_dimid                  ! netCDF ID of the constituent dimension
  integer               :: cell_dimid                         ! netCDF ID of the cell dimension
  integer               :: directions_id                      ! netCDF ID of the directions of the vectors of GAME
  integer               :: ncid                               ! netCDF file ID
  integer               :: interpolation_indices_scalar_id    ! interpolation indices for scalar quantities
  integer               :: interpolation_weights_scalar_id    ! interpolation weights for scalar quantities
  integer               :: interpolation_indices_vector_id    ! interpolation indices for vector quantities
  integer               :: interpolation_weights_vector_id    ! interpolation weights for vector quantities
  integer               :: dimids_vector_2(2)                 ! vector of two netCDF dimensions
  integer               :: dimids_vector_3(3)                 ! vector of three netCDF dimensions
  integer               :: n_condensed_constituents           ! number of condensed constituents of GAME or L-GAME
  integer               :: n_constituents                     ! number of constituents of GAME or L-GAME
  integer               :: closest_index                      ! index needed for vertical interpolation/extrapolation
  integer               :: other_index                        ! index needed for vertical interpolation/extrapolation
  integer               :: sp_id                              ! netCDF ID of the surface pressure of the input system
  integer               :: z_surf_id                          ! netCDF ID of the z-coordinates of the surface of the input system
  integer               :: z_coords_id                        ! netCDF ID of the z-coordinates of the input system
  integer               :: t_in_id                            ! netCDF ID of the temperature of the input system
  integer               :: spec_hum_id                        ! netCDF ID of the specific humidity of the input system
  integer               :: u_id                               ! netCDF ID of the zonal wind of the input system
  integer               :: v_id                               ! netCDF ID of the meridional wind of the input system
  integer               :: sst_id                             ! netCDF ID of the sea surface temperature of the input system
  integer               :: densities_background_id            ! netCDF ID of the mass densities of the background state
  integer               :: tke_id                             ! netCDF ID of the specific turbulent kinetic energy (input and output)
  integer               :: t_soil_id                          ! netCDF ID of the soil temperature (input and output)
  integer               :: edge_dimid                         ! netCDF ID of the edge dimension
  integer               :: single_double_dimid                ! netCDF ID of a single double
  integer               :: densities_id                       ! netCDF ID of the mass densities
  integer               :: temperature_id                     ! netCDF ID of the temperature
  integer               :: wind_h_id                          ! netCDF ID of the horizontal wind
  integer               :: wind_v_id                          ! netCDF ID of the vertical wind
  integer               :: soil_layer_dimid                   ! netCDF ID of the soil layer dimension
  integer               :: level_dimid                        ! netCDF of the level dimension
  integer               :: layer_dimid                        ! netCDF of the layer dimension
  integer               :: nsoillays                          ! number of soil layers of the GAME/L-GAME grid
  integer               :: n_points_per_layer_input           ! number of points per layer of the input system (excluding SST)
  integer               :: n_pentagons                        ! number of pentagons of the GAME grid
  integer               :: n_hexagons                         ! number of hexagons of the GAME grid
  integer               :: n_cells                            ! number of cells of the GAME grid
  integer               :: n_edges                            ! number of edges of the GAME grid
  integer               :: n_layers                           ! number of layers of the GAME grid
  integer               :: n_levels                           ! number of levels of the GAME grid
  integer               :: res_id                             ! resolution ID of GAME
  integer               :: oro_id                             ! orography ID of GAME or L-GAME
  integer               :: interpolation_indices_sst_id       ! netCDF ID of the interpolation indices for the SST interpolation
  integer               :: interpolation_weights_sst_id       ! netCDF ID of the interpolation weights for the SST interpolation
  real(wp)              :: rh                                 ! relative humidity value
  real(wp)              :: maximum_cloud_water_content        ! maximum cloud water content in (kg cloud)/(kg dry air)
  real(wp)              :: closest_value                      ! vertical interpolation weight
  real(wp)              :: other_value                        ! vertical interpolation weight
  real(wp)              :: df                                 ! used for computing the vertical gradient of a quantity
  real(wp)              :: dz                                 ! used for computing the vertical gradient of a quantity
  real(wp)              :: gradient                           ! vertical gradient of a quantity
  real(wp)              :: delta_z                            ! used for vertical interpolation/extrapolation
  real(wp)              :: b                                  ! used for vertically integrating the hydrostatic equation
  real(wp)              :: c                                  ! used for vertically integrating the hydrostatic equation
  real(wp)              :: u_local                            ! zonal wind at one gridpoint
  real(wp)              :: v_local                            ! meridional wind at one gridpoint
  real(wp)              :: vector_to_minimize(n_layers_input) ! vector needed for vertical interpolations
  real(wp), allocatable :: latitudes_game(:)                  ! latitudes of the cell centers of GAME
  real(wp), allocatable :: longitudes_game(:)                 ! longitudes of the cell centers of GAME
  real(wp), allocatable :: z_game(:,:)                        ! z-coordinates of the scalar points of GAME
  real(wp), allocatable :: directions(:)                      ! directions of the horizontal vectors of GAME
  real(wp), allocatable :: z_game_wind(:,:)                   ! z-coordinates of the horizontal vector points of GAME
  real(wp), allocatable :: gravity_potential_game(:,:)        ! gravity potential of GAME
  real(wp), allocatable :: densities_background(:,:,:)        ! densities of the background state
  real(wp), allocatable :: tke(:,:)                           ! specific turbulent kinetic energy (input = output)
  real(wp), allocatable :: t_soil(:,:)                        ! soil temperature (input = output)
  real(wp), allocatable :: z_coords_input_model(:,:)          ! vertical coordinates of the input model's gridpoints
  real(wp), allocatable :: temperature_in(:,:)                ! input temperature
  real(wp), allocatable :: spec_hum_in(:,:)                   ! input specific humidity
  real(wp), allocatable :: u_wind_in(:,:)                     ! input zonal wind
  real(wp), allocatable :: v_wind_in(:,:)                     ! input meridional wind
  real(wp), allocatable :: z_surf_in(:)                       ! input surface elevation
  real(wp), allocatable :: p_surf_in(:)                       ! input surface pressure
  real(wp), allocatable :: sst_in(:)                          ! input sea surface temperature
  real(wp), allocatable :: temperature_out(:,:)               ! resulting sea surface temperature
  real(wp), allocatable :: spec_hum_out(:,:)                  ! resulting specific humidity
  real(wp), allocatable :: pressure_lowest_layer_out(:)       ! resulting pressure in the lowest layer
  real(wp), allocatable :: density_moist_out(:,:)             ! resulting mass density of the moist air
  real(wp), allocatable :: temperature_v(:,:)                 ! resulting virtual temperature (only needed as a helper variable)
  real(wp), allocatable :: densities_out(:,:,:)               ! resulting mass densities
  real(wp), allocatable :: exner(:,:)                         ! resulting Exner pressure, only needed as a helper variables
  real(wp), allocatable :: wind_out_h(:,:)                    ! resulting horizontal wind
  real(wp), allocatable :: wind_out_v(:,:)                    ! resulting vertical wind
  real(wp), allocatable :: sst_out(:)                         ! resulting sea surface temperature
  real(wp), allocatable :: interpolation_weights_scalar(:,:)  ! interpolation weights for scalar quantities
  real(wp), allocatable :: interpolation_weights_vector(:,:)  ! interpolation weights for vector quantities
  real(wp), allocatable :: interpolation_weights_sst(:,:)     ! interpolation weights for the SST
  integer,  allocatable :: interpolation_indices_scalar(:,:)  ! interpolation indices for scalar quantities
  integer,  allocatable :: interpolation_indices_vector(:,:)  ! interpolation indices for vector quantities
  integer,  allocatable :: interpolation_indices_sst(:,:)     ! interpolation indices for the SST
  character(len=4)      :: year_string                        ! year of the initialization time as a string (command line argument)
  character(len=4)      :: n_layers_string                    ! number of layers of GAME or L-GAME as a string (command line argument)
  character(len=2)      :: month_string                       ! month of the initialization time as a string (command line argument)
  character(len=2)      :: day_string                         ! day of the initialization time as a string (command line argument)
  character(len=2)      :: hour_string                        ! hour of the initialization time as a string (command line argument)
  character(len=2)      :: res_id_string                      ! resolution ID of GAME as a string (command line argument)
  character(len=2)      :: nsoillays_string                   ! number of soil layer of GAME or L-GAME as a string (command line argument)
  character(len=128)    :: oro_id_string                      ! orography ID as a string (command line argument)
  character(len=128)    :: real2game_root_dir                 ! root directory of real2GAME (command line argument)
  character(len=128)    :: model_home_dir                     ! root directory of GAME or L-GAME (command line argument)
  character(len=256)    :: background_state_file              ! netCDF file containing the background state (command line argument)
  character(len=256)    :: geo_prop_file                      ! grid file of GAME or L-GAME
  character(len=256)    :: input_file                         ! file containing the data to interpolate to the model grid
  character(len=256)    :: interpol_file                      ! nnetCDF file containing interpolation indices and weights
  character(len=256)    :: output_file                        ! output filename
  
  call get_command_argument(1,res_id_string)
  read(res_id_string,*) res_id
  call get_command_argument(2,n_layers_string)
  read(n_layers_string,*) n_layers
  ! grid properties
  n_pentagons = 12
  n_hexagons = 10*(2**(2*res_id)-1)
  n_cells = n_pentagons+n_hexagons
  n_edges = (5*n_pentagons/2 + 6/2*n_hexagons)
  n_levels = n_layers+1
  n_condensed_constituents = 4
  n_constituents = 6
  call get_command_argument(3,nsoillays_string)
  read(nsoillays_string,*) nsoillays
  call get_command_argument(4,year_string)
  call get_command_argument(5,month_string)
  call get_command_argument(6,day_string)
  call get_command_argument(7,hour_string)
  call get_command_argument(8,model_home_dir)
  call get_command_argument(9,oro_id_string)
  read(oro_id_string,*) oro_id
  call get_command_argument(10,background_state_file)
  call get_command_argument(11,real2game_root_dir)
  
  n_points_per_layer_input = n_points_per_layer_input_icon_global
  
  write(*,*) "Background state file: ",trim(background_state_file)
  
  ! allocating memory for the grid properties
  allocate(latitudes_game(n_cells))
  allocate(longitudes_game(n_cells))
  allocate(z_game(n_cells,n_layers))
  allocate(directions(n_edges))
  allocate(z_game_wind(n_edges,n_layers))
  allocate(gravity_potential_game(n_cells,n_layers))
  
  ! initializing all array values with zeros
  !$omp parallel workshare
  latitudes_game = 0._wp
  longitudes_game = 0._wp
  z_game = 0._wp
  directions = 0._wp
  z_game_wind = 0._wp
  gravity_potential_game = 0._wp
  !$omp end parallel workshare
  
  ! reading the grid properties
  geo_prop_file = trim(model_home_dir) // "/grid_generator/grids/RES" // trim(int2string(res_id)) // "_L" // &
                  trim(int2string(n_layers)) // "_ORO" // trim(int2string(oro_id)) // ".nc"
  write(*,*) "Grid file: ",trim(geo_prop_file)
  write(*,*) "Reading grid file ..."
  call nc_check(nf90_open(trim(geo_prop_file),NF90_NOWRITE,ncid))
  call nc_check(nf90_inq_varid(ncid,"lat_c",latitudes_game_id))
  call nc_check(nf90_inq_varid(ncid,"lon_c",longitudes_game_id))
  call nc_check(nf90_inq_varid(ncid,"z_scalar",z_game_id))
  call nc_check(nf90_inq_varid(ncid,"direction",directions_id))
  call nc_check(nf90_inq_varid(ncid,"z_vector_h",z_game_wind_id))
  call nc_check(nf90_inq_varid(ncid,"gravity_potential",gravity_potential_game_id))
  call nc_check(nf90_get_var(ncid,latitudes_game_id,latitudes_game))
  call nc_check(nf90_get_var(ncid,longitudes_game_id,longitudes_game))
  call nc_check(nf90_get_var(ncid,z_game_id,z_game))
  call nc_check(nf90_get_var(ncid,directions_id,directions))
  call nc_check(nf90_get_var(ncid,z_game_wind_id,z_game_wind))
  call nc_check(nf90_get_var(ncid,gravity_potential_game_id,gravity_potential_game))
  call nc_check(nf90_close(ncid))
  write(*,*) "Grid file read."
  
  ! constructing the filename of the input file for GAME
  output_file = trim(model_home_dir) // "/nwp_init/" // year_string // month_string // day_string // hour_string // ".nc"
  
  ! These are the arrays of the background state.
  allocate(tke(n_cells,n_layers))
  allocate(t_soil(n_cells,nsoillays))
  allocate(densities_background(n_cells,n_layers,n_constituents))
  !$omp parallel workshare
  tke = 0._wp
  t_soil = 0._wp
  densities_background = 0._wp
  !$omp end parallel workshare
  
  ! Reading the background state.
  write(*,*) "Reading background state ..."
  call nc_check(nf90_open(background_state_file,NF90_NOWRITE,ncid))
  call nc_check(nf90_inq_varid(ncid,"densities",densities_background_id))
  ltke_avail = .false.
  if (nf90_inq_varid(ncid,"tke",tke_id)==0) then
    ltke_avail = .true.
    call nc_check(nf90_inq_varid(ncid,"tke",tke_id))
    write(*,*) "TKE found in background state file."
  else
    write(*,*) "TKE not found in background state file."
  endif
  lt_soil_avail = .false.
  if (nf90_inq_varid(ncid,"t_soil",t_soil_id)==0) then
    lt_soil_avail = .true.
    call nc_check(nf90_inq_varid(ncid,"t_soil",t_soil_id))
    write(*,*) "Soil temperature found in background state file."
  else
    write(*,*) "Soil temperature not found in background state file."
  endif
  call nc_check(nf90_get_var(ncid,densities_background_id,densities_background))
  if (ltke_avail) then
    call nc_check(nf90_get_var(ncid,tke_id,tke))
  endif
  if (lt_soil_avail) then
    call nc_check(nf90_get_var(ncid,t_soil_id,t_soil))
  endif
  call nc_check(nf90_close(ncid))
  write(*,*) "Background state read."
  
  ! allocating the memory for the analysis of the other model
  allocate(z_coords_input_model(n_points_per_layer_input,n_layers_input))
  allocate(temperature_in(n_points_per_layer_input,n_layers_input))
  allocate(spec_hum_in(n_points_per_layer_input,n_layers_input))
  allocate(u_wind_in(n_points_per_layer_input,n_layers_input))
  allocate(v_wind_in(n_points_per_layer_input,n_layers_input))
  allocate(z_surf_in(n_points_per_layer_input))
  allocate(p_surf_in(n_points_per_layer_input))
  allocate(sst_in(n_sst_points))
  !$omp parallel workshare
  z_coords_input_model = 0._wp
  temperature_in = 0._wp
  spec_hum_in = 0._wp
  u_wind_in = 0._wp
  v_wind_in = 0._wp
  z_surf_in = 0._wp
  p_surf_in = 0._wp
  sst_in = 0._wp
  !$omp end parallel workshare
  
  ! determining the name of the input file
  input_file = trim(real2game_root_dir) // "/input/obs_" // year_string // month_string // day_string // hour_string // ".nc"
  
  ! reading the analysis of the other model
  write(*,*) "Reading input ..."
  call nc_check(nf90_open(trim(input_file),NF90_NOWRITE,ncid))
  ! Defining the variables.
  call nc_check(nf90_inq_varid(ncid,"z_height",z_coords_id))
  call nc_check(nf90_inq_varid(ncid,"temperature",t_in_id))
  call nc_check(nf90_inq_varid(ncid,"spec_humidity",spec_hum_id))
  call nc_check(nf90_inq_varid(ncid,"u_wind",u_id))
  call nc_check(nf90_inq_varid(ncid,"v_wind",v_id))
  call nc_check(nf90_inq_varid(ncid,"z_surface",z_surf_id))
  call nc_check(nf90_inq_varid(ncid,"pressure_surface",sp_id))
  call nc_check(nf90_inq_varid(ncid,"sst",sst_id))
  call nc_check(nf90_get_var(ncid,z_coords_id,z_coords_input_model))
  call nc_check(nf90_get_var(ncid,t_in_id,temperature_in))
  call nc_check(nf90_get_var(ncid,spec_hum_id,spec_hum_in))
  call nc_check(nf90_get_var(ncid,u_id,u_wind_in))
  call nc_check(nf90_get_var(ncid,v_id,v_wind_in))
  call nc_check(nf90_get_var(ncid,z_surf_id,z_surf_in))
  call nc_check(nf90_get_var(ncid,sp_id,p_surf_in))
  call nc_check(nf90_get_var(ncid,sst_id,sst_in))
  call nc_check(nf90_close(ncid))
  write(*,*) "Input read."
  
  ! memory alloction for the interpolation indices and weights
  allocate(interpolation_indices_scalar(n_avg_points,n_cells))
  allocate(interpolation_weights_scalar(n_avg_points,n_cells))
  allocate(interpolation_indices_vector(n_avg_points,n_edges))
  allocate(interpolation_weights_vector(n_avg_points,n_edges))
  allocate(interpolation_indices_sst(n_avg_points,n_cells))
  allocate(interpolation_weights_sst(n_avg_points,n_cells))
  !$omp parallel workshare
  interpolation_indices_scalar = 0
  interpolation_weights_scalar = 0._wp
  interpolation_indices_vector = 0
  interpolation_weights_vector = 0._wp
  interpolation_indices_sst = 0
  interpolation_weights_sst = 0._wp
  !$omp end parallel workshare
  
  write(*,*) "Reading the interpolation indices and weights."
  
  ! constructing the name of the interpolation indices and weights file
  interpol_file  = trim(real2game_root_dir) // "/interpolation_files/icon-global2game" // trim(int2string(res_id)) // ".nc"
  
  ! reading the interpolation file
  call nc_check(nf90_open(trim(interpol_file),NF90_NOWRITE,ncid))
  call nc_check(nf90_inq_varid(ncid,"interpolation_indices_scalar",interpolation_indices_scalar_id))
  call nc_check(nf90_inq_varid(ncid,"interpolation_weights_scalar",interpolation_weights_scalar_id))
  call nc_check(nf90_inq_varid(ncid,"interpolation_indices_vector",interpolation_indices_vector_id))
  call nc_check(nf90_inq_varid(ncid,"interpolation_weights_vector",interpolation_weights_vector_id))
  call nc_check(nf90_inq_varid(ncid,"interpolation_indices_sst",interpolation_indices_sst_id))
  call nc_check(nf90_inq_varid(ncid,"interpolation_weights_sst",interpolation_weights_sst_id))
  call nc_check(nf90_get_var(ncid,interpolation_indices_scalar_id,interpolation_indices_scalar))
  call nc_check(nf90_get_var(ncid,interpolation_weights_scalar_id,interpolation_weights_scalar))
  call nc_check(nf90_get_var(ncid,interpolation_indices_vector_id,interpolation_indices_vector))
  call nc_check(nf90_get_var(ncid,interpolation_weights_vector_id,interpolation_weights_vector))
  call nc_check(nf90_get_var(ncid,interpolation_indices_sst_id,interpolation_indices_sst))
  call nc_check(nf90_get_var(ncid,interpolation_weights_sst_id,interpolation_weights_sst))
  call nc_check(nf90_close(ncid))
  write(*,*) "Interpolation indices and weights read."
  
  ! Begin of the actual interpolation.
  
  ! INTERPOLATION OF SCALAR QUANTITIES
  ! ----------------------------------
  
  write(*,*) "Starting the interpolation of scalar quantities ..."
  
  ! These are the arrays for the result of the interpolation process.
  allocate(temperature_out(n_cells,n_layers))
  allocate(spec_hum_out(n_cells,n_layers))
  !$omp parallel workshare
  temperature_out = 0._wp
  spec_hum_out = 0._wp
  !$omp end parallel workshare
  
  !$omp parallel workshare
  temperature_out = 0._wp
  spec_hum_out = 0._wp
  !$omp end parallel workshare
  
  !$omp parallel do private(ji,jl,jm,jn,vector_to_minimize,closest_index,closest_value,other_index,other_value,df,dz, &
  !$omp gradient,delta_z)
  do jl=1,n_layers
    do ji=1,n_cells
      
      ! loop over all points over which the averaging is executed
      do jm=1,n_avg_points
        ! computing linear vertical interpolation
        ! vertical distance vector
        do jn=1,n_layers_input
          vector_to_minimize(jn) = abs(z_game(ji,jl) - z_coords_input_model(interpolation_indices_scalar(jm,ji),jn))
        enddo
        
        ! closest vertical index
        closest_index = find_min_index(vector_to_minimize)
        
        ! value at the closest vertical index
        closest_value = temperature_in(interpolation_indices_scalar(jm,ji),closest_index)
        
        other_index = closest_index-1
        if (z_game(ji,jl)<z_coords_input_model(interpolation_indices_scalar(jm,ji),closest_index)) then
          other_index = closest_index+1
        endif
        
        ! avoiding array excess
        if (other_index==n_layers_input+1) then
          other_index = closest_index-1
        endif
        
        ! the value at the second point used for vertical interpolation
        other_value = temperature_in(interpolation_indices_scalar(jm,ji),other_index)
        
        ! computing the vertical gradient of u in the input model
        df = closest_value - other_value
        dz = z_coords_input_model(interpolation_indices_scalar(jm,ji),closest_index) &
             - z_coords_input_model(interpolation_indices_scalar(jm,ji),other_index)
        gradient = df/dz
        
        delta_z = z_game(ji,jl) - z_coords_input_model(interpolation_indices_scalar(jm,ji),closest_index)
        
        ! vertical interpolation of the temperature
        temperature_out(ji,jl) = temperature_out(ji,jl) + interpolation_weights_scalar(jm,ji)*(closest_value + delta_z*gradient)
        
        ! vertical interpolation of the specific humidity
        closest_value = spec_hum_in(interpolation_indices_scalar(jm,ji),closest_index)
        other_value = spec_hum_in(interpolation_indices_scalar(jm,ji),other_index)
        ! computing the vertical gradient of the specific humidity in the input model
        df = closest_value - other_value
        gradient = df/dz
        
        ! specific humidity
        spec_hum_out(ji,jl) = spec_hum_out(ji,jl) + interpolation_weights_scalar(jm,ji)*(closest_value + delta_z*gradient)
      enddo
      ! avoiding negative humidities, which can occur through vertical extrapolation
      if (spec_hum_out(ji,jl)<0._wp) then
        spec_hum_out(ji,jl) = 0._wp
      endif
    enddo
  enddo
  !$omp end parallel do
  ! these array are now interpolated and not needed any further
  
  deallocate(spec_hum_in)
  deallocate(temperature_in)
  ! surface pressure interpolation
  allocate(pressure_lowest_layer_out(n_cells))
  !$omp parallel workshare
  pressure_lowest_layer_out = 0._wp
  !$omp end parallel workshare
  
  !$omp parallel do private(ji,jm)
  do ji=1,n_cells
    do jm=1,n_avg_points
      pressure_lowest_layer_out(ji) = pressure_lowest_layer_out(ji) &
      ! horizontal component of the interpolation
      + interpolation_weights_scalar(jm,ji)*p_surf_in(interpolation_indices_scalar(jm,ji)) &
      ! vertical component of the interpolation according to the barometric height formula
      *exp(-(z_game(ji,n_layers) - z_surf_in(interpolation_indices_scalar(jm,ji)))/scale_height)
    enddo
  enddo
  !$omp end parallel do
  deallocate(z_game)
  deallocate(z_surf_in)
  deallocate(p_surf_in)
  ! no more interpolation of scalar quantities will be executed,this is why these interpolation indices and weights can be freed
  deallocate(interpolation_indices_scalar)
  deallocate(interpolation_weights_scalar)
  
  ! density is determined out of the hydrostatic equation
  allocate(density_moist_out(n_cells,n_layers))
  ! firstly setting the virtual temperature
  allocate(temperature_v(n_cells,n_layers))
  !$omp parallel workshare
  density_moist_out = 0._wp
  temperature_v = temperature_out*(1._wp + spec_hum_out*(m_d/m_v-1._wp))
  !$omp end parallel workshare
  
  ! the Exner pressure is just a temporarily needed helper variable here to integrate the hydrostatic equation
  allocate(exner(n_cells,n_layers))
  !$omp parallel workshare
  exner = 0._wp
  !$omp end parallel workshare
  !$omp parallel do private(ji,jl,b,c)
  do ji=1,n_cells
    do jl=n_layers,1,-1
      if (jl==n_layers) then
        density_moist_out(ji,jl) = pressure_lowest_layer_out(ji)/(r_d*temperature_v(ji,jl))
        exner(ji,jl) = (density_moist_out(ji,jl)*r_d*temperature_v(ji,jl)/p_0)**(r_d/c_d_p)
      else
        ! solving a quadratic equation for the Exner pressure
        b = -0.5_wp*exner(ji,jl+1)/temperature_v(ji,jl+1)*(temperature_v(ji,jl) - temperature_v(ji,jl+1) &
        + 2._wp/c_d_p*(gravity_potential_game(ji,jl)-gravity_potential_game(ji,jl+1)))
        c = exner(ji,jl+1)**2*temperature_v(ji,jl)/temperature_v(ji,jl+1)
        exner(ji,jl) = b + (b**2 + c)**0.5_wp
        density_moist_out(ji,jl) = p_0*exner(ji,jl)**(c_d_p/r_d)/(r_d*temperature_v(ji,jl))
      endif
    enddo
  enddo
  !$omp end parallel do
  
  deallocate(temperature_v)
  deallocate(pressure_lowest_layer_out)
  deallocate(gravity_potential_game)
  ! the Exner pressure was only needed to integrate the hydrostatic equation
  deallocate(exner)
  ! end of the interpolation of the dry thermodynamic state
  write(*,*) "Interpolation of scalar quantities completed."
  
  ! WIND INTERPOLATION
  ! ------------------
  
  write(*,*) "Starting the wind interpolation ..."
  allocate(wind_out_h(n_edges,n_layers))
  allocate(wind_out_v(n_cells,n_levels))
  !$omp parallel workshare
  wind_out_h = 0._wp
  wind_out_v = 0._wp
  !$omp end parallel workshare
  
  ! loop over all horizontal vector points
  !$omp parallel do private(ji,jl,jm,jn,vector_to_minimize,closest_index,other_index,closest_value, &
  !$omp other_value,df,dz,gradient,delta_z,u_local,v_local)
  do jl=1,n_layers
    do ji=1,n_edges
      
      ! the u- and v-components of the wind at the gridpoint of GAME
      u_local = 0._wp
      v_local = 0._wp
      ! loop over all horizontal points that are used for averaging
      do jm=1,n_avg_points
        ! computing linear vertical interpolation
        ! vertical distance vector
        do jn=1,n_layers_input
          vector_to_minimize(jn) = abs(z_game_wind(ji,jl) - z_coords_input_model(interpolation_indices_vector(jm,ji),jn))
        enddo
        ! closest vertical index
        closest_index = find_min_index(vector_to_minimize)
        ! value at the closest vertical index
        closest_value = u_wind_in(interpolation_indices_vector(jm,ji),closest_index)
        
        other_index = closest_index-1
        if (z_game_wind(ji,jl)<z_coords_input_model(interpolation_indices_vector(jm,ji),closest_index)) then
          other_index = closest_index+1
        endif
        ! avoiding array excess
        if (other_index==n_layers_input+1) then
          other_index = closest_index-1
        endif
        
        ! the value at the second point used for vertical interpolation
        other_value = u_wind_in(interpolation_indices_vector(jm,ji),other_index)
        
        ! computing the vertical gradient of u in the input model
        df = closest_value - other_value
        dz = z_coords_input_model(interpolation_indices_vector(jm,ji),closest_index) &
             - z_coords_input_model(interpolation_indices_vector(jm,ji),other_index)
        gradient = df/dz
        
        delta_z = z_game_wind(ji,jl) - z_coords_input_model(interpolation_indices_vector(jm,ji),closest_index)
        
        u_local = u_local + interpolation_weights_vector(jm,ji)*(closest_value + gradient*delta_z)
        
        ! vertical interpolation of v
        closest_value = v_wind_in(interpolation_indices_vector(jm,ji),closest_index)
        other_value = v_wind_in(interpolation_indices_vector(jm,ji),other_index)
        ! computing the vertical gradient of v in the input model
        df = closest_value - other_value
        gradient = df/dz
        
        v_local = v_local + interpolation_weights_vector(jm,ji)*(closest_value + gradient*delta_z)
        
      enddo
      
      ! projection onto the direction of the vector in GAME
      wind_out_h(ji,jl) = u_local*cos(directions(ji)) + v_local*sin(directions(ji))
      
    enddo
  enddo
  !$omp end parallel do
  
  deallocate(u_wind_in)
  deallocate(v_wind_in)
  deallocate(z_coords_input_model)
  deallocate(directions)
  deallocate(z_game_wind)
  deallocate(interpolation_indices_vector)
  deallocate(interpolation_weights_vector)
  
  write(*,*) "Wind interpolation completed."
  
  ! INTERPOLATION OF THE SST
  ! ------------------------
  
  write(*,*) "Interpolating the SST to the model grid ..."
  allocate(sst_out(n_cells))
  !$omp parallel workshare
  sst_out = 0._wp
  !$omp end parallel workshare
  !$omp parallel do private(ji,jm)
  do ji=1,n_cells
    do jm=1,n_avg_points
      sst_out(ji) = sst_out(ji) + interpolation_weights_sst(jm,ji)*sst_in(interpolation_indices_sst(jm,ji))
    enddo
  enddo
  !$omp end parallel do
  deallocate(interpolation_indices_sst)
  deallocate(interpolation_weights_sst)
  deallocate(sst_in)
  deallocate(latitudes_game)
  deallocate(longitudes_game)
  
  write(*,*) "Interpolation of the SST completed."
  
  ! PREPARING THE OUTPUT
  ! --------------------
  
  ! setting the mass densities of the result
  allocate(densities_out(n_cells,n_layers,n_constituents))
  !$omp parallel workshare
  ! clouds and precipitation are set equal to the background state
  densities_out(:,:,1:n_condensed_constituents) = densities_background(:,:,1:n_condensed_constituents)
  densities_out(:,:,n_condensed_constituents+1) = density_moist_out
  densities_out(:,:,n_condensed_constituents+2) = spec_hum_out*density_moist_out
  !$omp end parallel workshare
  
  !$omp parallel workshare
  lno_hydrometeors_found = minval(densities_out(:,:,1:n_condensed_constituents)) == 0._wp &
                           .and. maxval(densities_out(:,:,1:n_condensed_constituents)) == 0._wp
  !$omp end parallel workshare
  ! If no hydrometeors are present in the background state, we add hydrometeors according to a very simple ansatz.
  if (lno_hydrometeors_found) then
    write(*,*) "No hydrometeors found in the background state."
    
    maximum_cloud_water_content = 0.2e-3_wp
    
    !$omp parallel do private(ji,jl,rh)
    do jl=1,n_layers
      do ji=1,n_cells
        ! computing the relative humidity at the gridpoint
        rh = rel_humidity(densities_out(ji,jl,n_condensed_constituents+2),temperature_out(ji,jl))
        ! in the case of a relative humidity above 95 %, we add clouds and precipitation
        if (rh>=0.95_wp) then
          if (temperature_out(ji,jl)>=t_0) then
            ! water clouds
            densities_out(ji,jl,4) = maximum_cloud_water_content &
            *(densities_out(ji,jl,n_condensed_constituents+1) - densities_out(ji,jl,n_condensed_constituents+2))
            ! rain
            densities_out(ji,jl,2) = 0.5_wp*densities_out(ji,jl,4)
          else
            ! ice clouds
            densities_out(ji,jl,3) = maximum_cloud_water_content &
            *(densities_out(ji,jl,n_condensed_constituents+1) - densities_out(ji,jl,n_condensed_constituents+2))
            ! snow
            densities_out(ji,jl,1) = 0.5_wp*densities_out(ji,jl,3)
          endif
        ! in the case of a relative humidity above 90 %, we add clouds
        elseif (rh>=0.9_wp) then
          if (temperature_out(ji,jl)>=t_0) then
            ! water clouds
            densities_out(ji,jl,4) = 0.7_wp*maximum_cloud_water_content &
            *(densities_out(ji,jl,n_condensed_constituents+1) - densities_out(ji,jl,n_condensed_constituents+2))
          else
            ! ice clouds
            densities_out(ji,jl,3) = 0.7_wp*maximum_cloud_water_content &
            *(densities_out(ji,jl,n_condensed_constituents+1) - densities_out(ji,jl,n_condensed_constituents+2))
          endif
        endif
      enddo
    enddo
    !$omp end parallel do
    write(*,*) "Hydrometeors set."
  endif
  
  deallocate(density_moist_out)
  deallocate(spec_hum_out)
  deallocate(densities_background)
  
  ! writing the result to a netCDF file
  ! -----------------------------------
  
  write(*,*) "Output file: ",trim(output_file)
  write(*,*) "Writing result to output file ..."
  call nc_check(nf90_create(trim(output_file),NF90_CLOBBER,ncid))
  
  ! declaring dimensions
  call nc_check(nf90_def_dim(ncid,"cell_index",n_cells,cell_dimid))
  call nc_check(nf90_def_dim(ncid,"edge_index",n_edges,edge_dimid))
  call nc_check(nf90_def_dim(ncid,"layer_index",n_layers,layer_dimid))
  call nc_check(nf90_def_dim(ncid,"level_index",n_levels,level_dimid))
  call nc_check(nf90_def_dim(ncid,"constituent_index",6,constituent_dimid))
  call nc_check(nf90_def_dim(ncid,"soil_layer_index",nsoillays,soil_layer_dimid))
  call nc_check(nf90_def_dim(ncid,"single_double_dimid_index",1,single_double_dimid))
  
  ! defining variables
  ! mass densities
  dimids_vector_3(1) = cell_dimid
  dimids_vector_3(2) = layer_dimid
  dimids_vector_3(3) = constituent_dimid
  call nc_check(nf90_def_var(ncid,"densities",NF90_REAL,dimids_vector_3,densities_id))
  call nc_check(nf90_put_att(ncid,densities_id,"units","kg/m^3"))
  
  ! temperature
  dimids_vector_2(1) = cell_dimid
  dimids_vector_2(2) = layer_dimid
  call nc_check(nf90_def_var(ncid,"temperature",NF90_REAL,dimids_vector_2,temperature_id))
  call nc_check(nf90_put_att(ncid,temperature_id,"units","K"))
  
  ! horizontal wind
  dimids_vector_2(1) = edge_dimid
  dimids_vector_2(2) = layer_dimid
  call nc_check(nf90_def_var(ncid,"wind_h",NF90_REAL,dimids_vector_2,wind_h_id))
  call nc_check(nf90_put_att(ncid,wind_h_id,"units","m/s"))
  
  ! vertical wind
  dimids_vector_2(1) = cell_dimid
  dimids_vector_2(2) = level_dimid
  call nc_check(nf90_def_var(ncid,"wind_v",NF90_REAL,dimids_vector_2,wind_v_id))
  call nc_check(nf90_put_att(ncid,wind_v_id,"units","m/s"))
  
  ! SST
  call nc_check(nf90_def_var(ncid,"sst",NF90_REAL,cell_dimid,sst_id))
  call nc_check(nf90_put_att(ncid,sst_id,"units","K"))
  
  ! TKE
  if (ltke_avail) then
    dimids_vector_2(1) = cell_dimid
    dimids_vector_2(2) = layer_dimid
    call nc_check(nf90_def_var(ncid,"tke",NF90_REAL,dimids_vector_2,tke_id))
    call nc_check(nf90_put_att(ncid,tke_id,"units","J/kg"))
  endif
  
  ! soil temperature
  if (lt_soil_avail) then
    dimids_vector_2(1) = cell_dimid
    dimids_vector_2(2) = soil_layer_dimid
    call nc_check(nf90_def_var(ncid,"t_soil",NF90_REAL,dimids_vector_2,t_soil_id))
    call nc_check(nf90_put_att(ncid,t_soil_id,"units","K"))
  endif
  
  call nc_check(nf90_enddef(ncid))
  call nc_check(nf90_put_var(ncid,densities_id,densities_out))
  call nc_check(nf90_put_var(ncid,temperature_id,temperature_out))
  call nc_check(nf90_put_var(ncid,wind_h_id,wind_out_h))
  call nc_check(nf90_put_var(ncid,wind_v_id,wind_out_v))
  call nc_check(nf90_put_var(ncid,sst_id,sst_out))
  if (ltke_avail) then
    call nc_check(nf90_put_var(ncid,tke_id,tke))
  endif
  if (lt_soil_avail) then
    call nc_check(nf90_put_var(ncid,t_soil_id,t_soil))
  endif
  call nc_check(nf90_close(ncid))
  write(*,*) "Result successfully written."
  
  ! freeing the stil occupied memory
  deallocate(densities_out)
  deallocate(temperature_out)
  deallocate(wind_out_h)
  deallocate(wind_out_v)
  deallocate(sst_out)
  deallocate(tke)
  deallocate(t_soil)
  
end program control

















