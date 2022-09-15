! This source file is part of real2GAME, which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/real2GAME

program control

  ! This file prepares the horizontal interpolation from the foreign model to GAME.
  
  use netcdf
  use eccodes
  use mo_shared, only: wp,n_points_per_layer_input,n_avg_points,M_PI, &
                       nc_check,int2string,find_min_index,calculate_distance_h
  
  implicit none
  
  integer               :: ji,jk,jm,oro_id,model_target_id,nlat,nlon,ncid,lat_game_id,lat_game_wind_id, &
                           lon_game_id,lon_game_wind_id,lat_lgame_id,lon_lgame_id, &
                           lon_lgame_wind_u_id,lat_lgame_wind_u_id,lon_lgame_wind_v_id, &
                           lat_lgame_wind_v_id,dim_vector(2),interpolation_indices_scalar_id, &
                           interpolation_weights_scalar_id,interpolation_indices_vector_u_id, &
                           interpolation_weights_vector_u_id,interpolation_indices_vector_v_id, &
                           interpolation_weights_vector_v_id,scalar_dimid,u_dimid,v_dimid,avg_dimid, &
                           vector_dimid,jfile,jgrib,interpolation_indices_vector_id,interpolation_weights_vector_id, &
                           res_id,n_layers,n_pentagons,n_hexagons,n_cells,n_vectors_h
  real(wp)              :: sum_of_weights,interpol_exp
  integer,  allocatable :: interpolation_indices_scalar(:,:),interpolation_indices_vector(:,:), &
                           interpolation_indices_vector_u(:,:),interpolation_indices_vector_v(:,:)
  real(wp), allocatable :: lat_input_model(:),lon_input_model(:),lat_game(:),lon_game(:), &
                           lat_game_wind(:),lon_game_wind(:),interpolation_weights_scalar(:,:), &
                           interpolation_weights_vector(:,:),distance_vector(:),lat_lgame(:,:), &
                           lon_lgame(:,:),lat_lgame_wind_u(:,:),lon_lgame_wind_u(:,:), &
                           lat_lgame_wind_v(:,:),lon_lgame_wind_v(:,:), &
                           interpolation_weights_vector_u(:,:),interpolation_weights_vector_v(:,:)
  character(len=4)      :: year_string,nlat_string,nlon_string,n_layers_string
  character(len=8)      :: interpol_exp_string
  character(len=2)      :: month_string,day_string,hour_string,oro_id_string,model_target_id_string, &
                           res_id_string
  character(len=128)    :: real2game_root_dir,model_home_dir,lgame_grid
  character(len=256)    :: lat_obs_file,lon_obs_file,geo_pro_file,output_file
  
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
  call get_command_argument(9,nlat_string)
  read(nlat_string,*) nlat
  call get_command_argument(10,nlon_string)
  read(nlon_string,*) nlon
  call get_command_argument(11,interpol_exp_string)
  read(interpol_exp_string,*) interpol_exp
  call get_command_argument(12,lgame_grid)
  call get_command_argument(13,res_id_string)
  read(res_id_string,*) res_id
  call get_command_argument(14,n_layers_string)
  read(n_layers_string,*) n_layers
  ! grid properties
  n_pentagons = 12
  n_hexagons = 10*(2**(2*res_id)-1)
  n_cells = n_pentagons+n_hexagons
  n_vectors_h = (5*n_pentagons/2 + 6/2*n_hexagons)

  allocate(lat_input_model(n_points_per_layer_input))
  allocate(lon_input_model(n_points_per_layer_input))
  
  ! Properties of the input model's grid.
  ! latitudes of the grid
  lat_obs_file = trim(real2game_root_dir) // "/interpolation_creator/icon_global_icosahedral_time-invariant_" // year_string &
                 // month_string // day_string // hour_string // "_CLAT.grib2"
    
  call codes_open_file(jfile,trim(lat_obs_file),"r")
  call codes_grib_new_from_file(jfile,jgrib)
  call codes_get(jgrib,"values",lat_input_model)
  call codes_release(jgrib)
  call codes_close_file(jfile)
  
  ! transforming the latitude coordinates of the input model from degrees to radians
  !$omp parallel do private(ji)
  do ji=1,n_points_per_layer_input
    lat_input_model(ji) = 2._wp*M_PI*lat_input_model(ji)/360._wp
  enddo
  !$omp end parallel do
  
  ! longitudes of the grid
  lon_obs_file = trim(real2game_root_dir) // "/interpolation_creator/icon_global_icosahedral_time-invariant_" // year_string &
                 // month_string // day_string // hour_string // "_CLON.grib2"
  
  call codes_open_file(jfile,trim(lon_obs_file),"r")
  call codes_grib_new_from_file(jfile,jgrib)
  call codes_get(jgrib,"values",lon_input_model)
  call codes_release(jgrib)
  call codes_close_file(jfile)
  
  ! transforming the longitude coordinates of the input model from degrees to radians
  !$omp parallel do private(ji)
  do ji=1,n_points_per_layer_input
    lon_input_model(ji) = 2._wp*M_PI*lon_input_model(ji)/360._wp
  enddo
  !$omp end parallel do
  
  ! GAME
  if (model_target_id==1) then
  
    ! reading the horizontal coordinates of the grid of GAME
    allocate(lat_game(n_cells))
    allocate(lon_game(n_cells))
    allocate(lat_game_wind(n_vectors_h))
    allocate(lon_game_wind(n_vectors_h))
    geo_pro_file = trim(model_home_dir) // "/grid_generator/grids/RES" // trim(int2string(res_id)) //"_L" &
                   // trim(int2string(n_layers)) // "_ORO" // trim(int2string(oro_id)) // ".nc"
    write(*,*) "Grid file: ",trim(geo_pro_file)
    write(*,*) "Reading grid file of GAME ..."
    call nc_check(nf90_open(geo_pro_file,NF90_NOWRITE,ncid))
     
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
    allocate(interpolation_indices_scalar(n_cells,n_avg_points))
    allocate(interpolation_weights_scalar(n_cells,n_avg_points))
    allocate(interpolation_indices_vector(n_vectors_h,n_avg_points))
    allocate(interpolation_weights_vector(n_vectors_h,n_avg_points))
    
    ! executing the actual interpolation
    write(*,*) "Calculating interpolation indices and weights ..."
    
    allocate(distance_vector(n_points_per_layer_input))
    !$omp parallel do private(ji,jk,distance_vector,sum_of_weights)
    do ji=1,n_cells
      do jk=1,n_points_per_layer_input
        distance_vector(jk) = calculate_distance_h(lat_game(ji),lon_game(ji), &
                                                   lat_input_model(jk),lon_input_model(jk),1._wp)
      enddo
      sum_of_weights = 0._wp
      do jk=1,n_avg_points
        interpolation_indices_scalar(ji,jk) = find_min_index(distance_vector,n_points_per_layer_input)
        interpolation_weights_scalar(ji,jk) = 1._wp/distance_vector(interpolation_indices_scalar(ji,jk))**interpol_exp
        distance_vector(interpolation_indices_scalar(ji,jk)) = 2._wp*M_PI
        sum_of_weights = sum_of_weights + interpolation_weights_scalar(ji,jk)
      enddo
      do jk=1,n_avg_points
        interpolation_weights_scalar(ji,jk) = interpolation_weights_scalar(ji,jk)/sum_of_weights
      enddo
    enddo
    !$omp end parallel do
  
    !$omp parallel do private(ji,jk,sum_of_weights,distance_vector)
    do ji=1,n_vectors_h
      do jk=1,n_points_per_layer_input
        distance_vector(jk) = calculate_distance_h(lat_game_wind(ji),lon_game_wind(ji), &
                                                   lat_input_model(jk),lon_input_model(jk),1._wp)
      enddo
      sum_of_weights = 0._wp
      do jk=1,n_avg_points
        interpolation_indices_vector(ji,jk) = find_min_index(distance_vector,n_points_per_layer_input)
        interpolation_weights_vector(ji,jk) = 1._wp/distance_vector(interpolation_indices_vector(ji,jk))**interpol_exp
        distance_vector(interpolation_indices_vector(ji,jk)) = 2._wp*M_PI
        sum_of_weights = sum_of_weights + interpolation_weights_vector(ji,jk)
      enddo
      do jk=1,n_avg_points
        interpolation_weights_vector(ji,jk) = interpolation_weights_vector(ji,jk)/sum_of_weights
      enddo
    enddo
    !$omp end parallel do
    
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
    call nc_check(nf90_def_dim(ncid,"cell_index",n_cells,scalar_dimid))
    call nc_check(nf90_def_dim(ncid,"vector_index",n_vectors_h,vector_dimid))
    call nc_check(nf90_def_dim(ncid,"interpol_index",n_avg_points,avg_dimid))
    dim_vector(1) = scalar_dimid
    dim_vector(2) = avg_dimid
    call nc_check(nf90_def_var(ncid,"interpolation_indices_scalar",NF90_INT,dim_vector,interpolation_indices_scalar_id))
    call nc_check(nf90_def_var(ncid,"interpolation_weights_scalar",NF90_REAL,dim_vector,interpolation_weights_scalar_id))
    dim_vector(1) = vector_dimid
    call nc_check(nf90_def_var(ncid,"interpolation_indices_vector",NF90_INT,dim_vector,interpolation_indices_vector_id))
    call nc_check(nf90_def_var(ncid,"interpolation_weights_vector",NF90_REAL,dim_vector,interpolation_weights_vector_id))
    call nc_check(nf90_enddef(ncid))
    call nc_check(nf90_put_var(ncid,interpolation_indices_scalar_id,interpolation_indices_scalar))
    call nc_check(nf90_put_var(ncid,interpolation_weights_scalar_id,interpolation_weights_scalar))
    call nc_check(nf90_put_var(ncid,interpolation_indices_vector_id,interpolation_indices_vector))
    call nc_check(nf90_put_var(ncid,interpolation_weights_vector_id,interpolation_weights_vector))
    call nc_check(nf90_close(ncid))
      
    ! freeing the memory
    deallocate(interpolation_indices_scalar)
    deallocate(interpolation_weights_scalar)
    deallocate(interpolation_indices_vector)
    deallocate(interpolation_weights_vector)
    
  endif
  
  ! L-GAME
  if (model_target_id==2) then
  
    ! reading the horizontal coordinates of the grid of L-GAME
    allocate(lat_lgame(nlon,nlat))
    allocate(lon_lgame(nlon,nlat))
    allocate(lat_lgame_wind_u(nlon,nlat+1))
    allocate(lon_lgame_wind_u(nlon,nlat+1))
    allocate(lat_lgame_wind_v(nlon+1,nlat))
    allocate(lon_lgame_wind_v(nlon+1,nlat))
    geo_pro_file = trim(model_home_dir) // "/grids/" // trim(lgame_grid)
    write(*,*) "Grid file:",geo_pro_file
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
    allocate(interpolation_indices_scalar(nlat*nlon,n_avg_points))
    allocate(interpolation_weights_scalar(nlat*nlon,n_avg_points))
    allocate(interpolation_indices_vector_u(nlat*(nlon+1),n_avg_points))
    allocate(interpolation_weights_vector_u(nlat*(nlon+1),n_avg_points))
    allocate(interpolation_indices_vector_v((nlat+1)*nlon,n_avg_points))
    allocate(interpolation_weights_vector_v((nlat+1)*nlon,n_avg_points))
    
    ! executing the actual interpolation
    write(*,*) "Calculating interpolation indices and weights ..."
    
    !$omp parallel do private(ji,jk,jm,sum_of_weights,distance_vector)
    do ji=1,nlat
      do jk=1,nlon
        do jm=1,n_points_per_layer_input
          distance_vector(jm) = calculate_distance_h(lat_lgame(ji,jk),lon_lgame(ji,jk), &
                                                     lat_input_model(jm),lon_input_model(jm),1._wp)
        enddo
        sum_of_weights = 0._wp
        do jm=1,n_avg_points
          interpolation_indices_scalar((jk-1)*nlat+ji,jm) = find_min_index(distance_vector,n_points_per_layer_input)
          interpolation_weights_scalar((jk-1)*nlat+ji,jm) = 1._wp/ &
                                      distance_vector(interpolation_indices_scalar((jk-1)*nlat+ji,jm))**interpol_exp
          distance_vector(interpolation_indices_scalar((jk-1)*nlat+ji,jm)) = 2._wp*M_PI
          sum_of_weights = sum_of_weights + interpolation_weights_scalar((jk-1)*nlat+ji,jm)
        enddo
        do jm=1,n_avg_points
          interpolation_weights_scalar((jk-1)*nlat+ji,jm) = interpolation_weights_scalar((jk-1)*nlat+ji,jm)/sum_of_weights
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    !$omp parallel do private(ji,jk,jm,sum_of_weights,distance_vector)
    do ji=1,nlat
      do jk=1,nlon+1
        do jm=1,n_points_per_layer_input
          distance_vector(jm) = calculate_distance_h(lat_lgame_wind_u(ji,jk), &
                                lon_lgame_wind_u(ji,jk),lat_input_model(jm),lon_input_model(jm),1._wp)
        enddo
        sum_of_weights = 0._wp
        do jm=1,n_avg_points
          interpolation_indices_vector_u((jk-1)*nlat+ji,jm) = find_min_index(distance_vector,n_points_per_layer_input)
          interpolation_weights_vector_u((jk-1)*nlat+ji,jm) = 1._wp/ &
                                        distance_vector(interpolation_indices_vector_u((jk-1)*nlat+ji,jm))**interpol_exp
          distance_vector(interpolation_indices_vector_u((jk-1)*nlat+ji,jm)) = 2._wp*M_PI
          sum_of_weights = sum_of_weights + interpolation_weights_vector_u((jk-1)*nlat+ji,jm)
        enddo
        do jm=1,n_avg_points
          interpolation_weights_vector_u((jk-1)*nlat+ji,jm) = interpolation_weights_vector_u((jk-1)*nlat+ji,jm)/sum_of_weights
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    !$omp parallel do private(ji,jk,jm,sum_of_weights,distance_vector)
    do ji=1,nlat+1
      do jk=1,nlon
        do jm=1,n_points_per_layer_input
          distance_vector(jm) = calculate_distance_h(lat_lgame_wind_v(ji,jk),lon_lgame_wind_v(ji,jk), &
                                                     lat_input_model(jm),lon_input_model(jm),1._wp)
        enddo
        sum_of_weights = 0._wp
        do jm=1,n_avg_points
          interpolation_indices_vector_v((jk-1)*(nlat+1)+ji,jm) = find_min_index(distance_vector,n_points_per_layer_input)
          interpolation_weights_vector_v((jk-1)*(nlat+1)+ji,jm) = 1._wp/distance_vector( &
                                                       interpolation_indices_vector_v((jk-1)*(nlat+1)+ji,jm))**interpol_exp
          distance_vector(interpolation_indices_vector_v((jk-1)*(nlat+1)+ji,jm)) = 2._wp*M_PI
          sum_of_weights = sum_of_weights + interpolation_weights_vector_v((jk-1)*(nlat+1)+ji,jm)
        enddo
        do jm=1,n_avg_points
          interpolation_weights_vector_v((jk-1)*(nlat+1)+ji,jm) &
          = interpolation_weights_vector_v((jk-1)*(nlat+1)+ji,jm)/sum_of_weights
        enddo
      enddo
    enddo
    !$omp end parallel do
    
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
    
    ! writing the result to a NetCDF file
    output_file = trim(real2game_root_dir) // "/interpolation_files/icon-d22lgame_" // trim(lgame_grid)
    write(*,*) "Starting to write to output file ..."
    call nc_check(nf90_create(trim(output_file),NF90_CLOBBER,ncid))
    call nc_check(nf90_def_dim(ncid,"cell_index",nlat*nlon,scalar_dimid))
    call nc_check(nf90_def_dim(ncid,"u_index",nlat*(nlon+1),u_dimid))
    call nc_check(nf90_def_dim(ncid,"v_index",(nlat+1)*nlon,v_dimid))
    call nc_check(nf90_def_dim(ncid,"interpol_index",n_avg_points,avg_dimid))
    dim_vector(1) = scalar_dimid
    dim_vector(2) = avg_dimid
    call nc_check(nf90_def_var(ncid,"interpolation_indices_scalar",NF90_INT,dim_vector,interpolation_indices_scalar_id))
    call nc_check(nf90_def_var(ncid,"interpolation_weights_scalar",NF90_REAL,dim_vector,interpolation_weights_scalar_id))
    dim_vector(1) = u_dimid
    call nc_check(nf90_def_var(ncid,"interpolation_indices_vector_u",NF90_INT,dim_vector,interpolation_indices_vector_u_id))
    call nc_check(nf90_def_var(ncid,"interpolation_weights_vector_u",NF90_REAL,dim_vector,interpolation_weights_vector_u_id))
    dim_vector(1) = v_dimid
    call nc_check(nf90_def_var(ncid,"interpolation_indices_vector_v",NF90_INT,dim_vector,interpolation_indices_vector_v_id))
    call nc_check(nf90_def_var(ncid,"interpolation_weights_vector_v",NF90_REAL,dim_vector,interpolation_weights_vector_v_id))
    call nc_check(nf90_enddef(ncid))
    call nc_check(nf90_put_var(ncid,interpolation_indices_scalar_id,interpolation_indices_scalar))
    call nc_check(nf90_put_var(ncid,interpolation_weights_scalar_id,interpolation_weights_scalar))
    call nc_check(nf90_put_var(ncid,interpolation_indices_vector_u_id,interpolation_indices_vector_u))
    call nc_check(nf90_put_var(ncid,interpolation_weights_vector_u_id,interpolation_weights_vector_u))
    call nc_check(nf90_put_var(ncid,interpolation_indices_vector_v_id,interpolation_indices_vector_v))
    call nc_check(nf90_put_var(ncid,interpolation_weights_vector_v_id,interpolation_weights_vector_v))
    call nc_check(nf90_close(ncid))
    
    ! freeing the memory
    deallocate(interpolation_indices_scalar)
    deallocate(interpolation_weights_scalar)
    deallocate(interpolation_indices_vector_u)
    deallocate(interpolation_weights_vector_u)
    deallocate(interpolation_indices_vector_v)
    deallocate(interpolation_weights_vector_v)
  
  endif
  
  write(*,*) "Finished."
  
end program control







