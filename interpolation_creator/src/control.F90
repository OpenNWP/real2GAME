! This source file is part of real2GAME, which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/real2GAME

program control
  
  ! This program computes the horizontal interpolation indices and weights for interpolating from a foreign model to GAME or L-GAME.
  
  use netcdf
  use eccodes
  use mo_shared, only: wp,n_avg_points,M_PI,nc_check,int2string,find_min_index,calculate_distance_h, &
                       n_points_per_layer_input_icon_global,n_points_per_layer_input_icon_d2,n_sst_points
  
  implicit none
  
  integer               :: ji                                        ! horizontal index
  integer               :: jk                                        ! horizontal index
  integer               :: jm                                        ! interpolation index
  integer               :: oro_id                                    ! orography ID of GAME (only needed for reading the correct grid file)
  integer               :: model_target_id                           ! the ID of the model to interpolate to (1: GAME, 2: L-GAME)
  integer               :: ny                                        ! number of points in the y-direction of the L-GAME grid
  integer               :: nx                                        ! number of points in the x-direction of the L-GAME grid
  integer               :: ncid                                      ! netCDF file ID
  integer               :: lat_game_id                               ! netCDF ID of the latitudes of the cell centers of GAME
  integer               :: lat_game_wind_id                          ! netCDF ID of the latitudes of the edges of GAME
  integer               :: lon_game_id                               ! netCDF ID of the longitudes of the cell centers of GAME
  integer               :: lon_game_wind_id                          ! netCDF ID of the longitudes of the edges of GAME
  integer               :: lat_lgame_id                              ! netCDF ID of the latitudes of the cell centers of L-GAME
  integer               :: lon_lgame_id                              ! netCDF ID of the longitudes of the cell centers of L-GAME
  integer               :: lat_lgame_wind_u_id                       ! netCDF ID of the latitudes of the u-vector points of L-GAME
  integer               :: lon_lgame_wind_u_id                       ! netCDF ID of the longitudes of the u-vector points of L-GAME
  integer               :: lat_lgame_wind_v_id                       ! netCDF ID of the latitudes of the v-vector points of L-GAME
  integer               :: lon_lgame_wind_v_id                       ! netCDF ID of the longitudes of the v-vector points of L-GAME
  integer               :: dim_vector_2(2)                           ! vector of two netCDF dimensions
  integer               :: interpolation_indices_scalar_id           ! interpolation indices of the scalar gridpoints of L-GAME
  integer               :: interpolation_weights_scalar_id           ! interpolation weights of the scalar gridpoints of L-GAME
  integer               :: interpolation_indices_vector_u_id         ! interpolation indices of the u-vector gridpoints of L-GAME
  integer               :: interpolation_weights_vector_u_id         ! interpolation weights of the u-vector gridpoints of L-GAME
  integer               :: interpolation_indices_vector_v_id         ! interpolation indices of the v-vector gridpoints of L-GAME
  integer               :: interpolation_weights_vector_v_id         ! interpolation weights of the v-vector gridpoints of L-GAME
  integer               :: cell_dimid                                ! netCDF cell dimension ID
  integer               :: avg_dimid                                 ! netCDF averaging dimension ID
  integer               :: edge_dimid                                ! netCDF edge dimension ID
  integer               :: jfile                                     ! eccodes grib file ID
  integer               :: jgrib                                     ! eccodes grib message ID
  integer               :: interpolation_indices_vector_id           ! interpolation inidces of the vector gridpoints of GAME
  integer               :: interpolation_weights_vector_id           ! interpolation weights of the vector gridpoints of GAME
  integer               :: res_id                                    ! resolution ID of the GAME grid
  integer               :: n_layers                                  ! number of layers of the GAME grid
  integer               :: n_pentagons                               ! number of pentagons of the GAME grid
  integer               :: n_hexagons                                ! number of hexagons of the GAME grid
  integer               :: n_cells                                   ! number of cells of the GAME grid
  integer               :: n_edges                                   ! number of edges of the GAME grid
  integer               :: model_source_id                           ! ID of the input system (1: ICON-global, 2: GAME, 3: ICON-D2)
  integer               :: dim_vector_3(3)                           ! vector of three netCDF dimensions
  integer               :: y_dimid                                   ! netCDF dimid of the y-direction of L-GAME
  integer               :: x_dimid                                   ! netCDF dimid of the x-direction of L-GAME
  integer               :: yp1_dimid                                 ! netCDF dimid of the y-direction of L-GAME + 1 gridpoint
  integer               :: xp1_dimid                                 ! netCDF dimid of the x-direction of L-GAME + 1 gridpoint
  integer               :: n_points_per_layer_input                  ! number of points per layer of the input system
  integer               :: interpolation_indices_sst_id              ! netCDF ID of the SST interpolation indices
  integer               :: interpolation_weights_sst_id              ! netCDF ID of the SST interpolation weights
  real(wp)              :: sum_of_weights                            ! sum of interpolation weights used for normalization
  real(wp)              :: interpol_exp                              ! interpolation exponent (1/distance**interpol_exp-interpolation)
  integer,  allocatable :: interpolation_indices_scalar_game(:,:)    ! interpolation indices of scalar quantities to the GAME grid
  integer,  allocatable :: interpolation_indices_vector(:,:)         ! interpolation indices of vector quantities to the GAME grid
  integer,  allocatable :: interpolation_indices_vector_u(:,:,:)     ! interpolation indices of vector quantities to the u-vector points of the L-GAME grid
  integer,  allocatable :: interpolation_indices_vector_v(:,:,:)     ! interpolation indices of vector quantities to the v-vector points of the L-GAME grid
  integer,  allocatable :: interpolation_indices_scalar_lgame(:,:,:) ! interpolation indices of scalar quantities to the scalar data points of the L-GAME grid
  integer,  allocatable :: interpolation_indices_sst_game(:,:)       ! interpolation indices of the sea surface temperature to the GAME grid
  integer,  allocatable :: interpolation_indices_sst_lgame(:,:,:)    ! interpolation indices of the sea surface temperature to the L-GAME grid
  real(wp), allocatable :: interpolation_weights_sst_game(:,:)       ! interpolation weights of the sea surface temperature to the GAME grid
  real(wp), allocatable :: interpolation_weights_sst_lgame(:,:,:)    ! interpolation weights of the sea surface temperature to the L-GAME grid
  real(wp), allocatable :: lat_sst(:)                                ! latitudes of the sea surface temperature grid
  real(wp), allocatable :: lon_sst(:)                                ! longitudes of the sea surface temperature grid
  real(wp), allocatable :: lat_input_model(:)                        ! latitudes of the input system's gridpoints
  real(wp), allocatable :: lon_input_model(:)                        ! longitudes of the input system's gridpoints
  real(wp), allocatable :: lat_game(:)                               ! latitudes of the cell centers of GAME
  real(wp), allocatable :: lon_game(:)                               ! longitudes of the cell centers of GAME
  real(wp), allocatable :: lat_game_wind(:)                          ! latitudes of the edges of GAME
  real(wp), allocatable :: lon_game_wind(:)                          ! longitudes of the edges of GAME
  real(wp), allocatable :: interpolation_weights_scalar_game(:,:)    ! interpolation weights for scalar quantities of the GAME grid
  real(wp), allocatable :: interpolation_weights_vector(:,:)         ! interpolation weights for vector quantities of the GAME grid
  real(wp), allocatable :: distance_vector(:)                        ! vector of geodetic distances needed for computing the interpolation
  real(wp), allocatable :: lat_lgame(:,:)                            ! latitudes of the scalar gridpoints of L-GAME
  real(wp), allocatable :: lon_lgame(:,:)                            ! longitudes of the scalar gridpoints of L-GAME
  real(wp), allocatable :: lat_lgame_wind_u(:,:)                     ! latitudes of the u-vector gridpoints of L-GAME
  real(wp), allocatable :: lon_lgame_wind_u(:,:)                     ! longitudes of the u-vector gridpoints of L-GAME
  real(wp), allocatable :: lat_lgame_wind_v(:,:)                     ! latitudes of the v-vector gridpoints of L-GAME
  real(wp), allocatable :: lon_lgame_wind_v(:,:)                     ! longitudes of the v-vector gridpoints of L-GAME
  real(wp), allocatable :: interpolation_weights_scalar_lgame(:,:,:) ! interpolation weights for the scalar gridpoints of the L-GAME grid
  real(wp), allocatable :: interpolation_weights_vector_u(:,:,:)     ! interpolation weights for the u-vector gridpoints of the L-GAME grid
  real(wp), allocatable :: interpolation_weights_vector_v(:,:,:)     ! interpolation weights for the v-vector gridpoints of the L-GAME grid
  character(len=4)      :: ny_string                                 ! number of gridpoints in y-direction of the L-GAME grid (command line argument)
  character(len=4)      :: nx_string                                 ! number of gridpoints in x-direction of the L-GAME grid (command line argument)
  character(len=4)      :: n_layers_string                           ! number of layers of GAME or L-GAME as a string (command line argument)
  character(len=8)      :: interpol_exp_string                       ! the interpolation exponent as a string (command line argument)
  character(len=4)      :: year_string                               ! year of the grid file of the input system as a string (command line argument)
  character(len=2)      :: month_string                              ! month of the grid file of the input system as a string (command line argument)
  character(len=2)      :: day_string                                ! day of the grid file of the input system as a string (command line argument)
  character(len=2)      :: hour_string                               ! hour of the grid file of the input system as a string (command line argument)
  character(len=2)      :: oro_id_string                             ! orography ID of GAME or L-GAME as a string (command line argument)
  character(len=2)      :: model_target_id_string                    ! ID of the model to interpolate to as a string (1: GAME, 2: L-GAME; command line argument)
  character(len=2)      :: res_id_string                             ! resolution ID of the GAME grid as a string (command line argument)
  character(len=2)      :: model_source_id_string                    ! ID of the input system as a string (1: ICON-global, 2: GAME, 3: ICON-D2; command line argument)
  character(len=128)    :: real2game_root_dir                        ! root directory of real2GAME
  character(len=128)    :: model_home_dir                            ! root directory of GAME or L-GAME
  character(len=128)    :: lgame_grid                                ! name of the grid file of L-GAME (excluding directory)
  character(len=256)    :: lat_obs_file                              ! name of the file containing the latitudes of the input system
  character(len=256)    :: lon_obs_file                              ! name of the file containing the longitudes of the input system
  character(len=256)    :: geo_pro_file                              ! name of the grid file of GAME or L-GAME (including directory)
  character(len=256)    :: output_file                               ! name of the output file
  character(len=256)    :: sst_file                                  ! file containing the sea surface temperature
  
  ! shell arguments
  call get_command_argument(1,year_string)
  call get_command_argument(2,month_string)
  call get_command_argument(3,day_string)
  call get_command_argument(4,hour_string)
  call get_command_argument(5,real2game_root_dir)
  call get_command_argument(6,model_home_dir)
  call get_command_argument(7,oro_id_string)
  read(oro_id_string,*) oro_id
  call get_command_argument(8,model_target_id_string)
  read(model_target_id_string,*) model_target_id
  call get_command_argument(9,ny_string)
  read(ny_string,*) ny
  call get_command_argument(10,nx_string)
  read(nx_string,*) nx
  call get_command_argument(11,interpol_exp_string)
  read(interpol_exp_string,*) interpol_exp
  call get_command_argument(12,lgame_grid)
  call get_command_argument(13,res_id_string)
  read(res_id_string,*) res_id
  call get_command_argument(14,n_layers_string)
  read(n_layers_string,*) n_layers
  call get_command_argument(15,model_source_id_string)
  read(model_source_id_string,*) model_source_id
  
  ! grid properties
  n_pentagons = 12
  n_hexagons = 10*(2**(2*res_id)-1)
  n_cells = n_pentagons+n_hexagons
  n_edges = (5*n_pentagons/2 + 6/2*n_hexagons)
  
  ! Properties of the input model's grid.
  ! latitudes of the grid
  ! ICON-global
  if (model_source_id==1) then
    lat_obs_file = trim(real2game_root_dir) // "/interpolation_creator/icon_global_icosahedral_time-invariant_" // year_string &
                   // month_string // day_string // hour_string // "_CLAT.grib2"
    n_points_per_layer_input = n_points_per_layer_input_icon_global
  endif
  ! ICON-D2
  if (model_source_id==3) then
    lat_obs_file = trim(real2game_root_dir) // "/interpolation_creator/icon-d2_germany_icosahedral_time-invariant_" // year_string &
                   // month_string // day_string // hour_string // "_000_0_clat.grib2"
    n_points_per_layer_input = n_points_per_layer_input_icon_d2
  endif
  
  allocate(lat_input_model(n_points_per_layer_input))
  allocate(lon_input_model(n_points_per_layer_input))
  !$omp parallel workshare
  lat_input_model = 0._wp
  lon_input_model = 0._wp
  !$omp end parallel workshare
  
  call codes_open_file(jfile,trim(lat_obs_file),"r")
  call codes_grib_new_from_file(jfile,jgrib)
  call codes_get(jgrib,"values",lat_input_model)
  call codes_release(jgrib)
  call codes_close_file(jfile)
  
  ! transforming the latitude coordinates of the input model from degrees to radians
  !$omp parallel workshare
  lat_input_model = 2._wp*M_PI*lat_input_model/360._wp
  !$omp end parallel workshare
  
  ! longitudes of the grid
  ! ICON-global
  if (model_source_id==1) then
    lon_obs_file = trim(real2game_root_dir) // "/interpolation_creator/icon_global_icosahedral_time-invariant_" // year_string &
                   // month_string // day_string // hour_string // "_CLON.grib2"
  endif
  ! ICON-D2
  if (model_source_id==3) then
    lon_obs_file = trim(real2game_root_dir) // "/interpolation_creator/icon-d2_germany_icosahedral_time-invariant_" // year_string &
                   // month_string // day_string // hour_string // "_000_0_clon.grib2"
  endif
  
  call codes_open_file(jfile,trim(lon_obs_file),"r")
  call codes_grib_new_from_file(jfile,jgrib)
  call codes_get(jgrib,"values",lon_input_model)
  call codes_release(jgrib)
  call codes_close_file(jfile)
  
  ! transforming the longitude coordinates of the input model from degrees to radians
  !$omp parallel workshare
  lon_input_model = 2._wp*M_PI*lon_input_model/360._wp
  !$omp end parallel workshare
  
  ! SST points
  allocate(lat_sst(n_sst_points))
  allocate(lon_sst(n_sst_points))
  !$omp parallel workshare
  lat_sst = 0._wp
  lon_sst = 0._wp
  !$omp end parallel workshare
  
  ! SST filename
  sst_file = trim(real2game_root_dir) // "/interpolation_creator/rtgssthr_grb_0.5.grib2"
  
  call codes_open_file(jfile,trim(sst_file),"r")
  call codes_grib_new_from_file(jfile,jgrib)
  call codes_get(jgrib,"values",lat_sst)
  call codes_get(jgrib,"latitudes",lat_sst)
  call codes_get(jgrib,"longitudes",lon_sst)
  call codes_release(jgrib)
  call codes_close_file(jfile)
  
  ! transforming the coordinates of the SST file from degrees to radians
  !$omp parallel workshare
  lat_sst = 2._wp*M_PI*lat_sst/360._wp
  lon_sst = 2._wp*M_PI*lon_sst/360._wp
  !$omp end parallel workshare
  
  ! GAME
  if (model_target_id==1) then
    
    ! reading the horizontal coordinates of the grid of GAME
    allocate(lat_game(n_cells))
    allocate(lon_game(n_cells))
    allocate(lat_game_wind(n_edges))
    allocate(lon_game_wind(n_edges))
    !$omp parallel workshare
    lat_game = 0._wp
    lon_game = 0._wp
    lat_game_wind = 0._wp
    lon_game_wind = 0._wp
    !$omp end parallel workshare
    geo_pro_file = trim(model_home_dir) // "/grid_generator/grids/RES" // trim(int2string(res_id)) //"_L" &
                   // trim(int2string(n_layers)) // "_ORO" // trim(int2string(oro_id)) // ".nc"
    write(*,*) "Grid file: ",trim(geo_pro_file)
    write(*,*) "Reading grid file of GAME ..."
    call nc_check(nf90_open(trim(geo_pro_file),NF90_NOWRITE,ncid))
    
    call nc_check(nf90_inq_varid(ncid,"lat_c",lat_game_id))
    call nc_check(nf90_inq_varid(ncid,"lon_c",lon_game_id))
    call nc_check(nf90_inq_varid(ncid,"lat_e",lat_game_wind_id))
    call nc_check(nf90_inq_varid(ncid,"lon_e",lon_game_wind_id))
    call nc_check(nf90_get_var(ncid,lat_game_id,lat_game))
    call nc_check(nf90_get_var(ncid,lon_game_id,lon_game))
    call nc_check(nf90_get_var(ncid,lat_game_wind_id,lat_game_wind))
    call nc_check(nf90_get_var(ncid,lon_game_wind_id,lon_game_wind))
    call nc_check(nf90_close(ncid))
    write(*,*) "Grid file of GAME read."
    
    ! allocating memory for the result arrays
    allocate(interpolation_indices_scalar_game(n_avg_points,n_cells))
    allocate(interpolation_weights_scalar_game(n_avg_points,n_cells))
    allocate(interpolation_indices_vector(n_avg_points,n_edges))
    allocate(interpolation_weights_vector(n_avg_points,n_edges))
    allocate(interpolation_indices_sst_game(n_avg_points,n_cells))
    allocate(interpolation_weights_sst_game(n_avg_points,n_cells))
    !$omp parallel workshare
    interpolation_indices_scalar_game = 0
    interpolation_weights_scalar_game = 0._wp
    interpolation_indices_vector = 0
    interpolation_weights_vector = 0._wp
    interpolation_indices_sst_game = 0
    interpolation_weights_sst_game = 0._wp
    !$omp end parallel workshare
    
    ! executing the actual interpolation
    write(*,*) "Calculating interpolation indices and weights ..."
    
    allocate(distance_vector(n_points_per_layer_input))
    !$omp parallel workshare
    distance_vector = 0._wp
    !$omp end parallel workshare
    
    ! cell centers
    !$omp parallel do private(ji,jk,distance_vector,sum_of_weights)
    do ji=1,n_cells
      do jk=1,n_points_per_layer_input
        distance_vector(jk) = calculate_distance_h(lat_game(ji),lon_game(ji), &
                                                   lat_input_model(jk),lon_input_model(jk),1._wp)
      enddo
      sum_of_weights = 0._wp
      do jk=1,n_avg_points
        interpolation_indices_scalar_game(jk,ji) = find_min_index(distance_vector)
        interpolation_weights_scalar_game(jk,ji) = 1._wp/distance_vector(interpolation_indices_scalar_game(jk,ji))**interpol_exp
        distance_vector(interpolation_indices_scalar_game(jk,ji)) = 2._wp*M_PI
        sum_of_weights = sum_of_weights + interpolation_weights_scalar_game(jk,ji)
      enddo
      interpolation_weights_scalar_game(:,ji) = interpolation_weights_scalar_game(:,ji)/sum_of_weights
    enddo
    !$omp end parallel do
    
    ! edges
    !$omp parallel do private(ji,jk,sum_of_weights,distance_vector)
    do ji=1,n_edges
      do jk=1,n_points_per_layer_input
        distance_vector(jk) = calculate_distance_h(lat_game_wind(ji),lon_game_wind(ji), &
                                                   lat_input_model(jk),lon_input_model(jk),1._wp)
      enddo
      sum_of_weights = 0._wp
      do jk=1,n_avg_points
        interpolation_indices_vector(jk,ji) = find_min_index(distance_vector)
        interpolation_weights_vector(jk,ji) = 1._wp/distance_vector(interpolation_indices_vector(jk,ji))**interpol_exp
        distance_vector(interpolation_indices_vector(jk,ji)) = 2._wp*M_PI
        sum_of_weights = sum_of_weights + interpolation_weights_vector(jk,ji)
      enddo
      interpolation_weights_vector(:,ji) = interpolation_weights_vector(:,ji)/sum_of_weights
    enddo
    !$omp end parallel do
    
    deallocate(distance_vector)
    allocate(distance_vector(n_sst_points))
    
    ! SST points
    !$omp parallel do private(ji,jk,distance_vector,sum_of_weights)
    do ji=1,n_cells
      do jk=1,n_sst_points
        distance_vector(jk) = calculate_distance_h(lat_game(ji),lon_game(ji),lat_sst(jk),lon_sst(jk),1._wp)
      enddo
      sum_of_weights = 0._wp
      do jk=1,n_avg_points
        interpolation_indices_sst_game(jk,ji) = find_min_index(distance_vector)
        interpolation_weights_sst_game(jk,ji) = 1._wp/distance_vector(interpolation_indices_sst_game(jk,ji))**interpol_exp
        distance_vector(interpolation_indices_sst_game(jk,ji)) = 2._wp*M_PI
        sum_of_weights = sum_of_weights + interpolation_weights_sst_game(jk,ji)
      enddo
      interpolation_weights_sst_game(:,ji) = interpolation_weights_sst_game(:,ji)/sum_of_weights
    enddo
    !$omp end parallel do
    
    deallocate(distance_vector)
    deallocate(lat_sst)
    deallocate(lon_sst)
    
    write(*,*) "Calculating interpolation indices and weights finished."
    
    ! freeing memory we do not need further
    deallocate(lat_input_model)
    deallocate(lon_input_model)
    deallocate(lat_game)
    deallocate(lon_game)
    deallocate(lat_game_wind)
    deallocate(lon_game_wind)
    
    ! writing the result to a netCDF file
    output_file = trim(real2game_root_dir) // "/interpolation_files/icon-global2game" // trim(int2string(res_id)) // ".nc"
    write(*,*) "Starting to write to output file ..."
    call nc_check(nf90_create(trim(output_file),NF90_CLOBBER,ncid))
    call nc_check(nf90_def_dim(ncid,"cell_index",n_cells,cell_dimid))
    call nc_check(nf90_def_dim(ncid,"vector_index",n_edges,edge_dimid))
    call nc_check(nf90_def_dim(ncid,"interpol_index",n_avg_points,avg_dimid))
    dim_vector_2(1) = avg_dimid
    dim_vector_2(2) = cell_dimid
    call nc_check(nf90_def_var(ncid,"interpolation_indices_scalar",NF90_INT,dim_vector_2,interpolation_indices_scalar_id))
    call nc_check(nf90_def_var(ncid,"interpolation_weights_scalar",NF90_REAL,dim_vector_2,interpolation_weights_scalar_id))
    call nc_check(nf90_def_var(ncid,"interpolation_indices_sst",NF90_INT,dim_vector_2,interpolation_indices_sst_id))
    call nc_check(nf90_def_var(ncid,"interpolation_weights_sst",NF90_REAL,dim_vector_2,interpolation_weights_sst_id))
    dim_vector_2(1) = avg_dimid
    dim_vector_2(2) = edge_dimid
    call nc_check(nf90_def_var(ncid,"interpolation_indices_vector",NF90_INT,dim_vector_2,interpolation_indices_vector_id))
    call nc_check(nf90_def_var(ncid,"interpolation_weights_vector",NF90_REAL,dim_vector_2,interpolation_weights_vector_id))
    call nc_check(nf90_enddef(ncid))
    call nc_check(nf90_put_var(ncid,interpolation_indices_scalar_id,interpolation_indices_scalar_game))
    call nc_check(nf90_put_var(ncid,interpolation_weights_scalar_id,interpolation_weights_scalar_game))
    call nc_check(nf90_put_var(ncid,interpolation_indices_vector_id,interpolation_indices_vector))
    call nc_check(nf90_put_var(ncid,interpolation_weights_vector_id,interpolation_weights_vector))
    call nc_check(nf90_put_var(ncid,interpolation_indices_sst_id,interpolation_indices_sst_game))
    call nc_check(nf90_put_var(ncid,interpolation_weights_sst_id,interpolation_weights_sst_game))
    call nc_check(nf90_close(ncid))
    
    ! freeing the memory
    deallocate(interpolation_indices_scalar_game)
    deallocate(interpolation_weights_scalar_game)
    deallocate(interpolation_indices_vector)
    deallocate(interpolation_weights_vector)
    deallocate(interpolation_indices_sst_game)
    deallocate(interpolation_weights_sst_game)
    
  endif
  
  ! L-GAME
  if (model_target_id==2) then
    
    ! reading the horizontal coordinates of the grid of L-GAME
    allocate(lat_lgame(ny,nx))
    allocate(lon_lgame(ny,nx))
    allocate(lat_lgame_wind_u(ny,nx+1))
    allocate(lon_lgame_wind_u(ny,nx+1))
    allocate(lat_lgame_wind_v(ny+1,nx))
    allocate(lon_lgame_wind_v(ny+1,nx))
    !$omp parallel workshare
    lat_lgame = 0._wp
    lon_lgame = 0._wp
    lat_lgame_wind_u = 0._wp
    lon_lgame_wind_u = 0._wp
    lat_lgame_wind_v = 0._wp
    lon_lgame_wind_v = 0._wp
    !$omp end parallel workshare
    geo_pro_file = trim(model_home_dir) // "/grids/" // trim(lgame_grid)
    write(*,*) "Grid file:",trim(geo_pro_file)
    write(*,*) "Reading grid file of L-GAME ..."
    call nc_check(nf90_open(geo_pro_file,NF90_NOWRITE,ncid))
    
    call nc_check(nf90_inq_varid(ncid,"lat_geo",lat_lgame_id))
    call nc_check(nf90_inq_varid(ncid,"lon_geo",lon_lgame_id))
    call nc_check(nf90_inq_varid(ncid,"lat_geo_u",lat_lgame_wind_u_id))
    call nc_check(nf90_inq_varid(ncid,"lon_geo_u",lon_lgame_wind_u_id))
    call nc_check(nf90_inq_varid(ncid,"lat_geo_v",lat_lgame_wind_v_id))
    call nc_check(nf90_inq_varid(ncid,"lon_geo_v",lon_lgame_wind_v_id))
    call nc_check(nf90_get_var(ncid,lat_lgame_id,lat_lgame))
    call nc_check(nf90_get_var(ncid,lon_lgame_id,lon_lgame))
    call nc_check(nf90_get_var(ncid,lat_lgame_wind_u_id,lat_lgame_wind_u))
    call nc_check(nf90_get_var(ncid,lon_lgame_wind_u_id,lon_lgame_wind_u))
    call nc_check(nf90_get_var(ncid,lat_lgame_wind_v_id,lat_lgame_wind_v))
    call nc_check(nf90_get_var(ncid,lon_lgame_wind_v_id,lon_lgame_wind_v))
    call nc_check(nf90_close(ncid))
    write(*,*) "Grid file of L-GAME read."
    
    ! allocating memory for the result arrays
    allocate(interpolation_indices_scalar_lgame(n_avg_points,ny,nx))
    allocate(interpolation_weights_scalar_lgame(n_avg_points,ny,nx))
    allocate(interpolation_indices_vector_u(n_avg_points,ny,nx+1))
    allocate(interpolation_weights_vector_u(n_avg_points,ny,nx+1))
    allocate(interpolation_indices_vector_v(n_avg_points,ny+1,nx))
    allocate(interpolation_weights_vector_v(n_avg_points,ny+1,nx))
    allocate(interpolation_indices_sst_lgame(n_avg_points,ny,nx))
    allocate(interpolation_weights_sst_lgame(n_avg_points,ny,nx))
    !$omp parallel workshare
    interpolation_indices_scalar_lgame = 0
    interpolation_weights_scalar_lgame = 0._wp
    interpolation_indices_vector_u = 0
    interpolation_weights_vector_u = 0._wp
    interpolation_indices_vector_v = 0
    interpolation_weights_vector_v = 0._wp
    interpolation_indices_sst_lgame = 0
    interpolation_weights_sst_lgame = 0._wp
    !$omp end parallel workshare
    
    ! executing the actual interpolation
    write(*,*) "Calculating interpolation indices and weights ..."
    
    allocate(distance_vector(n_points_per_layer_input))
    !$omp parallel workshare
    distance_vector = 0._wp
    !$omp end parallel workshare
    
    ! scalar points
    !$omp parallel do private(ji,jk,jm,sum_of_weights,distance_vector)
    do jk=1,nx
      do ji=1,ny
        do jm=1,n_points_per_layer_input
          distance_vector(jm) = calculate_distance_h(lat_lgame(ji,jk),lon_lgame(ji,jk), &
                                                     lat_input_model(jm),lon_input_model(jm),1._wp)
        enddo
        sum_of_weights = 0._wp
        do jm=1,n_avg_points
          interpolation_indices_scalar_lgame(jm,ji,jk) = find_min_index(distance_vector)
          interpolation_weights_scalar_lgame(jm,ji,jk) = 1._wp/ &
                                      distance_vector(interpolation_indices_scalar_lgame(jm,ji,jk))**interpol_exp
          distance_vector(interpolation_indices_scalar_lgame(jm,ji,jk)) = 2._wp*M_PI
          sum_of_weights = sum_of_weights + interpolation_weights_scalar_lgame(jm,ji,jk)
        enddo
        interpolation_weights_scalar_lgame(:,ji,jk) = interpolation_weights_scalar_lgame(:,ji,jk)/sum_of_weights
      enddo
    enddo
    !$omp end parallel do
    
    ! u vector points
    !$omp parallel do private(ji,jk,jm,sum_of_weights,distance_vector)
    do jk=1,nx+1
      do ji=1,ny
        do jm=1,n_points_per_layer_input
          distance_vector(jm) = calculate_distance_h(lat_lgame_wind_u(ji,jk), &
                                lon_lgame_wind_u(ji,jk),lat_input_model(jm),lon_input_model(jm),1._wp)
        enddo
        sum_of_weights = 0._wp
        do jm=1,n_avg_points
          interpolation_indices_vector_u(jm,ji,jk) = find_min_index(distance_vector)
          interpolation_weights_vector_u(jm,ji,jk) = 1._wp/distance_vector(interpolation_indices_vector_u(jm,ji,jk))**interpol_exp
          distance_vector(interpolation_indices_vector_u(jm,ji,jk)) = 2._wp*M_PI
          sum_of_weights = sum_of_weights + interpolation_weights_vector_u(jm,ji,jk)
        enddo
        interpolation_weights_vector_u(:,ji,jk) = interpolation_weights_vector_u(:,ji,jk)/sum_of_weights
      enddo
    enddo
    !$omp end parallel do
    
    ! v vector points
    !$omp parallel do private(ji,jk,jm,sum_of_weights,distance_vector)
    do jk=1,nx
      do ji=1,ny+1
        do jm=1,n_points_per_layer_input
          distance_vector(jm) = calculate_distance_h(lat_lgame_wind_v(ji,jk),lon_lgame_wind_v(ji,jk), &
                                                     lat_input_model(jm),lon_input_model(jm),1._wp)
        enddo
        sum_of_weights = 0._wp
        do jm=1,n_avg_points
          interpolation_indices_vector_v(jm,ji,jk) = find_min_index(distance_vector)
          interpolation_weights_vector_v(jm,ji,jk) = 1._wp/distance_vector(interpolation_indices_vector_v(jm,ji,jk))**interpol_exp
          distance_vector(interpolation_indices_vector_v(jm,ji,jk)) = 2._wp*M_PI
          sum_of_weights = sum_of_weights + interpolation_weights_vector_v(jm,ji,jk)
        enddo
        interpolation_weights_vector_v(:,ji,jk) = interpolation_weights_vector_v(:,ji,jk)/sum_of_weights
      enddo
    enddo
    !$omp end parallel do
    
    deallocate(distance_vector)
    allocate(distance_vector(n_sst_points))
    
    ! SST
    !$omp parallel do private(ji,jk,jm,sum_of_weights,distance_vector)
    do jk=1,nx
      do ji=1,ny
        do jm=1,n_sst_points
          distance_vector(jm) = calculate_distance_h(lat_lgame(ji,jk),lon_lgame(ji,jk),lat_sst(jm),lon_sst(jm),1._wp)
        enddo
        sum_of_weights = 0._wp
        do jm=1,n_avg_points
          interpolation_indices_sst_lgame(jm,ji,jk) = find_min_index(distance_vector)
          interpolation_weights_sst_lgame(jm,ji,jk) = 1._wp/ &
                                      distance_vector(interpolation_indices_sst_lgame(jm,ji,jk))**interpol_exp
          distance_vector(interpolation_indices_sst_lgame(jm,ji,jk)) = 2._wp*M_PI
          sum_of_weights = sum_of_weights + interpolation_weights_sst_lgame(jm,ji,jk)
        enddo
        interpolation_weights_sst_lgame(:,ji,jk) = interpolation_weights_sst_lgame(:,ji,jk)/sum_of_weights
      enddo
    enddo
    !$omp end parallel do
    
    deallocate(distance_vector)
    deallocate(lat_sst)
    deallocate(lon_sst)
    
    write(*,*) "Calculating interpolation indices and weights finished."
    
    ! freeing memory we do not need further
    deallocate(lat_input_model)
    deallocate(lon_input_model)
    deallocate(lat_lgame)
    deallocate(lon_lgame)
    deallocate(lat_lgame_wind_u)
    deallocate(lon_lgame_wind_u)
    deallocate(lat_lgame_wind_v)
    deallocate(lon_lgame_wind_v)
    
    ! writing the result to a netCDF file
    output_file = trim(real2game_root_dir) // "/interpolation_files/icon-d22lgame_" // trim(lgame_grid)
    write(*,*) "Starting to write to output file ..."
    call nc_check(nf90_create(trim(output_file),NF90_CLOBBER,ncid))
    call nc_check(nf90_def_dim(ncid,"y_index",ny,y_dimid))
    call nc_check(nf90_def_dim(ncid,"x_index",nx,x_dimid))
    call nc_check(nf90_def_dim(ncid,"yp1_index",ny+1,yp1_dimid))
    call nc_check(nf90_def_dim(ncid,"xp1_index",nx+1,xp1_dimid))
    call nc_check(nf90_def_dim(ncid,"interpol_index",n_avg_points,avg_dimid))
    dim_vector_3(1) = avg_dimid
    dim_vector_3(2) = y_dimid
    dim_vector_3(3) = x_dimid
    call nc_check(nf90_def_var(ncid,"interpolation_indices_scalar",NF90_INT,dim_vector_3,interpolation_indices_scalar_id))
    call nc_check(nf90_def_var(ncid,"interpolation_weights_scalar",NF90_REAL,dim_vector_3,interpolation_weights_scalar_id))
    call nc_check(nf90_def_var(ncid,"interpolation_indices_sst",NF90_INT,dim_vector_3,interpolation_indices_sst_id))
    call nc_check(nf90_def_var(ncid,"interpolation_weights_sst",NF90_REAL,dim_vector_3,interpolation_weights_sst_id))
    dim_vector_3(1) = avg_dimid
    dim_vector_3(2) = y_dimid
    dim_vector_3(3) = xp1_dimid
    call nc_check(nf90_def_var(ncid,"interpolation_indices_vector_u",NF90_INT,dim_vector_3,interpolation_indices_vector_u_id))
    call nc_check(nf90_def_var(ncid,"interpolation_weights_vector_u",NF90_REAL,dim_vector_3,interpolation_weights_vector_u_id))
    dim_vector_3(1) = avg_dimid
    dim_vector_3(2) = yp1_dimid
    dim_vector_3(3) = x_dimid
    call nc_check(nf90_def_var(ncid,"interpolation_indices_vector_v",NF90_INT,dim_vector_3,interpolation_indices_vector_v_id))
    call nc_check(nf90_def_var(ncid,"interpolation_weights_vector_v",NF90_REAL,dim_vector_3,interpolation_weights_vector_v_id))
    call nc_check(nf90_enddef(ncid))
    call nc_check(nf90_put_var(ncid,interpolation_indices_scalar_id,interpolation_indices_scalar_lgame))
    call nc_check(nf90_put_var(ncid,interpolation_weights_scalar_id,interpolation_weights_scalar_lgame))
    call nc_check(nf90_put_var(ncid,interpolation_indices_vector_u_id,interpolation_indices_vector_u))
    call nc_check(nf90_put_var(ncid,interpolation_weights_vector_u_id,interpolation_weights_vector_u))
    call nc_check(nf90_put_var(ncid,interpolation_indices_vector_v_id,interpolation_indices_vector_v))
    call nc_check(nf90_put_var(ncid,interpolation_weights_vector_v_id,interpolation_weights_vector_v))
    call nc_check(nf90_put_var(ncid,interpolation_indices_sst_id,interpolation_indices_sst_lgame))
    call nc_check(nf90_put_var(ncid,interpolation_weights_sst_id,interpolation_weights_sst_lgame))
    call nc_check(nf90_close(ncid))
    
    ! freeing the memory
    deallocate(interpolation_indices_scalar_lgame)
    deallocate(interpolation_weights_scalar_lgame)
    deallocate(interpolation_indices_vector_u)
    deallocate(interpolation_weights_vector_u)
    deallocate(interpolation_indices_vector_v)
    deallocate(interpolation_weights_vector_v)
    deallocate(interpolation_indices_sst_lgame)
    deallocate(interpolation_weights_sst_lgame)
    
  endif
  
  write(*,*) "Finished."
  
end program control







