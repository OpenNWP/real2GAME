! This source file is part of real2GAME,which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/real2GAME

program control

  ! This file prepares the horizontal interpolation from the foreign model to GAME.
  
  use netcdf
  
  implicit none
  
  ! shell arguments
  char year_string(strlen(argv(1)) + 1)
  strcpy(year_string,argv(1))
  char month_string(strlen(argv(2)) + 1)
  strcpy(month_string,argv(2))
  char day_string(strlen(argv(3)) + 1)
  strcpy(day_string,argv(3))
  char hour_string(strlen(argv(4)) + 1)
  strcpy(hour_string,argv(4))
  char real2game_root_dir(strlen(argv(5)) + 1)
  strcpy(real2game_root_dir,argv(5))
  char model_home_dir(strlen(argv(6)) + 1)
  strcpy(model_home_dir,argv(6))
  int oro_id = strtod(argv(7),NULL)
  int model_target_id = strtod(argv(8),NULL)
  int nlat = strtod(argv(9),NULL)
  int nlon = strtod(argv(10),NULL)
  int interpol_exp = strtod(argv(11),NULL)
  char lgame_grid(strlen(argv(12)) + 1)
  strcpy(lgame_grid,argv(12))
  
  int no_of_points_per_layer_input_model
  no_of_points_per_layer_input_model = 2949120

  double *latitudes_input_model = malloc(no_of_points_per_layer_input_model*sizeof(double))
  double *longitudes_input_model = malloc(no_of_points_per_layer_input_model*sizeof(double))
  
  int err
  codes_handle *handle = NULL
  
  ! Properties of the input model's grid.
  ! latitudes of the grid
  char lat_obs_file_pre(200)
  sprintf(lat_obs_file_pre ,"%s/interpolation_creator/icon_global_icosahedral_time-invariant_%s%s%s%s_CLAT.grib2",
  real2game_root_dir,year_string,month_string,day_string,hour_string)
  char lat_obs_file(strlen(lat_obs_file_pre) + 1)
  strcpy(lat_obs_file,lat_obs_file_pre)
  FILE *ECC_FILE
  ECC_FILE = fopen(lat_obs_file,"r")
  handle = codes_handle_new_from_file(NULL,ECC_FILE,PRODUCT_GRIB,err)
  if (err != 0) ECCERR(err)
  size_t no_of_points_per_layer_input_model_SIZE_T
  no_of_points_per_layer_input_model_SIZE_T = (size_t) no_of_points_per_layer_input_model
  ECCCHECK(codes_get_double_array(handle,"values",latitudes_input_model,no_of_points_per_layer_input_model_SIZE_T))
  codes_handle_delete(handle)
  fclose(ECC_FILE)
  
  ! transforming the latitude coordinates of the input model from degrees to radians
  #pragma omp parallel for
  do (int i = 0 i < no_of_points_per_layer_input_model ++i)
    latitudes_input_model(i) = 2._wp*M_PI*latitudes_input_model(i)/360._wp
  enddo
  
  ! longitudes of the grid
  char lon_obs_file_pre(200)
  sprintf(lon_obs_file_pre ,"%s/interpolation_creator/icon_global_icosahedral_time-invariant_%s%s%s%s_CLON.grib2",
  real2game_root_dir,year_string,month_string,day_string,hour_string)
  char lon_obs_file(strlen(lon_obs_file_pre) + 1)
  strcpy(lon_obs_file,lon_obs_file_pre)
  
  ECC_FILE = fopen(lon_obs_file,"r")
  handle = codes_handle_new_from_file(NULL,ECC_FILE,PRODUCT_GRIB,err)
  if (err != 0) ECCERR(err)
  no_of_points_per_layer_input_model_SIZE_T = (size_t) no_of_points_per_layer_input_model
  ECCCHECK(codes_get_double_array(handle,"values",longitudes_input_model,no_of_points_per_layer_input_model_SIZE_T))
  codes_handle_delete(handle)
  fclose(ECC_FILE)
  
  ! transforming the longitude coordinates of the input model from degrees to radians
  !$omp parallel do private(ji)
  do (int i = 0 i < no_of_points_per_layer_input_model ++i)
    longitudes_input_model(i) = 2._wp*M_PI*longitudes_input_model(i)/360._p
  enddo
  !$omp end parallel do
  
  ! GAME
  if (model_target_id==1) then
  
    ! reading the horizontal coordinates of the grid of GAME
    double *latitudes_game = malloc(n_scalars_h*sizeof(double))
    double *longitudes_game = malloc(n_scalars_h*sizeof(double))
    double *latitudes_game_wind = malloc(n_vectors_h*sizeof(double))
    double *longitudes_game_wind = malloc(n_vectors_h*sizeof(double))
    char geo_pro_file_pre(200)
    sprintf(geo_pro_file_pre,"%s/grid_generator/grids/RES%d_L%d_ORO%d.nc",model_home_dir,RES_ID,N_LAYERS,oro_id)
    char geo_pro_file(strlen(geo_pro_file_pre) + 1)
    strcpy(geo_pro_file,geo_pro_file_pre)
    printf("Grid file: %s\n",geo_pro_file)
    write(*,*) "Reading grid file of GAME ..."
    int ncid
    call nc_check(nc_open(geo_pro_file,NC_NOWRITE,ncid))
    
    int latitudes_game_id,latitudes_game_wind_id,longitudes_game_id,longitudes_game_wind_id  
    call nc_check(nc_inq_varid(ncid,"latitude_scalar",latitudes_game_id))
    call nc_check(nc_inq_varid(ncid,"longitude_scalar",longitudes_game_id))
    call nc_check(nc_inq_varid(ncid,"latitude_vector",latitudes_game_wind_id))
    call nc_check(nc_inq_varid(ncid,"longitude_vector",longitudes_game_wind_id))
    call nc_check(nc_get_var_double(ncid,latitudes_game_id,latitudes_game))
    call nc_check(nc_get_var_double(ncid,longitudes_game_id,longitudes_game))
    call nc_check(nc_get_var_double(ncid,latitudes_game_wind_id,latitudes_game_wind))
    call nc_check(nc_get_var_double(ncid,longitudes_game_wind_id,longitudes_game_wind))
    call nc_check(nc_close(ncid))
    printf("Grid file of GAME read.\n")
  
    ! allocating memory for the result arrays
    int (*interpolation_indices_scalar,n_avg_points) = malloc(sizeof(int(n_scalars_h,n_avg_points)))
    double (*interpolation_weights_scalar,n_avg_points) = malloc(sizeof(double(n_scalars_h,n_avg_points)))
    int (*interpolation_indices_vector,n_avg_points) = malloc(sizeof(int(n_vectors_h,n_avg_points)))
    double (*interpolation_weights_vector,n_avg_points) = malloc(sizeof(double(n_vectors_h,n_avg_points)))
    
    ! executing the actual interpolation
    write(*,*) "Calculating interpolation indices and weights ..."
    
    !$omp parallel do private(ji)
    do (int i = 0 i < n_scalars_h ++i)
      double *distance_vector = malloc(no_of_points_per_layer_input_model*sizeof(double))
      do (int j = 0 j < no_of_points_per_layer_input_model ++j)
        distance_vector(j) = calculate_distance_h(&latitudes_game(i),longitudes_game(i),latitudes_input_model(j),longitudes_input_model(j),1._wp)
      enddo
      sum_of_weights = 0._wp
      do (int j = 0 j < n_avg_points ++j)
        interpolation_indices_scalar(i,j) = find_min_index(distance_vector,no_of_points_per_layer_input_model)
        interpolation_weights_scalar(i,j) = 1._wp/(pow(distance_vector(interpolation_indices_scalar(i,j)),interpol_exp))
        distance_vector(interpolation_indices_scalar(i,j)) = 2._wp*M_PI
        sum_of_weights += interpolation_weights_scalar(i,j)
      enddo
      deallocate(distance_vector)
      do (int j = 0 j < n_avg_points ++j)
        interpolation_weights_scalar(i,j) = interpolation_weights_scalar(i,j)/sum_of_weights
      enddo
    enddo
    !$omp end parallel do
  
    !$omp parallel do private(ji)
    do (int i = 0 i < n_vectors_h ++i)
      double *distance_vector = malloc(no_of_points_per_layer_input_model*sizeof(double))
      do (int j = 0 j < no_of_points_per_layer_input_model ++j)
        distance_vector(j) =  calculate_distance_h(&latitudes_game_wind(i),longitudes_game_wind(i),latitudes_input_model(j),longitudes_input_model(j),1._wp)
      enddo
      sum_of_weights = 0._wp
      do (int j = 0 j < n_avg_points ++j)
        interpolation_indices_vector(i,j) = find_min_index(distance_vector,no_of_points_per_layer_input_model)
        interpolation_weights_vector(i,j) = 1._wp/(pow(distance_vector(interpolation_indices_vector(i,j)),interpol_exp))
        distance_vector(interpolation_indices_vector(i,j)) = 2._wp*M_PI
        sum_of_weights += interpolation_weights_vector(i,j)
      enddo
      deallocate(distance_vector)
      do (int j = 0 j < n_avg_points ++j)
        interpolation_weights_vector(i,j) = interpolation_weights_vector(i,j)/sum_of_weights
      enddo
    enddo
    !$omp end parallel do
    
    write(*,*) "Calculating interpolation indices and weights finished."
    
    ! freeing memory we do not need further
    deallocate(latitudes_input_model)
    deallocate(longitudes_input_model)
    deallocate(latitudes_game)
    deallocate(longitudes_game)
    deallocate(latitudes_game_wind)
    deallocate(longitudes_game_wind)
    
    ! writing the result to a NetCDF file
    char output_file_pre(200)
    sprintf(output_file_pre,"%s/interpolation_files/icon-global2game%d.nc",real2game_root_dir,RES_ID)
    char output_file(strlen(output_file_pre) + 1)
    strcpy(output_file,output_file_pre)
    printf("Starting to write to output file ...\n")
    int interpolation_indices_scalar_id,interpolation_weights_scalar_id,interpolation_indices_vector_id,interpolation_weights_vector_id,
    scalar_dimid,vector_dimid,avg_dimid
    int dim_vector(2)
    call nc_check(nc_create(output_file,NC_CLOBBER,ncid))
    call nc_check(nc_def_dim(ncid,"scalar_index",n_scalars_h,scalar_dimid))
    call nc_check(nc_def_dim(ncid,"vector_index",n_vectors_h,vector_dimid))
    call nc_check(nc_def_dim(ncid,"interpol_index",n_avg_points,avg_dimid))
    dim_vector(1) = scalar_dimid
    dim_vector(2) = avg_dimid
    call nc_check(nc_def_var(ncid,"interpolation_indices_scalar",NC_INT,2,dim_vector,interpolation_indices_scalar_id))
    call nc_check(nc_def_var(ncid,"interpolation_weights_scalar",NC_DOUBLE,2,dim_vector,interpolation_weights_scalar_id))
    dim_vector(1) = vector_dimid
    call nc_check(nc_def_var(ncid,"interpolation_indices_vector",NC_INT,2,dim_vector,interpolation_indices_vector_id))
    call nc_check(nc_def_var(ncid,"interpolation_weights_vector",NC_DOUBLE,2,dim_vector,interpolation_weights_vector_id))
    call nc_check(nc_enddef(ncid))
    call nc_check(nc_put_var_int(ncid,interpolation_indices_scalar_id,interpolation_indices_scalar))
    call nc_check(nc_put_var_double(ncid,interpolation_weights_scalar_id,interpolation_weights_scalar))
    call nc_check(nc_put_var_int(ncid,interpolation_indices_vector_id,interpolation_indices_vector))
    call nc_check(nc_put_var_double(ncid,interpolation_weights_vector_id,interpolation_weights_vector))
    call nc_check(nc_close(ncid))
      
    ! freeing the memory
    deallocate(interpolation_indices_scalar)
    deallocate(interpolation_weights_scalar)
    deallocate(interpolation_indices_vector)
    deallocate(interpolation_weights_vector)
    
  endif
  
  ! L-GAME
  if (model_target_id==2) then
  
    ! reading the horizontal coordinates of the grid of L-GAME
    double (*latitudes_lgame,nlat) = malloc(sizeof(double(nlon,nlat)))
    double (*longitudes_lgame,nlat) = malloc(sizeof(double(nlon,nlat)))
    double (*latitudes_lgame_wind_u,nlat) = malloc(sizeof(double(nlon,nlat+1)))
    double (*longitudes_lgame_wind_u,nlat) = malloc(sizeof(double(nlon,nlat+1)))
    double (*latitudes_lgame_wind_v,nlat) = malloc(sizeof(double(nlon+1,nlat)))
    double (*longitudes_lgame_wind_v,nlat) = malloc(sizeof(double(nlon+1,nlat)))
    char geo_pro_file_pre(200)
    sprintf(geo_pro_file_pre,"%s/grids/%s",model_home_dir,lgame_grid)
    char geo_pro_file(strlen(geo_pro_file_pre) + 1)
    strcpy(geo_pro_file,geo_pro_file_pre)
    printf("Grid file: %s\n",geo_pro_file)
    printf("Reading grid file of L-GAME ...\n")
    int ncid
    call nc_check(nc_open(geo_pro_file,NC_NOWRITE,ncid))
    
    int latitudes_lgame_id,longitudes_lgame_id,longitudes_lgame_wind_u_id,latitudes_lgame_wind_u_id,
    longitudes_lgame_wind_v_id,latitudes_lgame_wind_v_id
    call nc_check(nc_inq_varid(ncid,"lat_geo",latitudes_lgame_id))
    call nc_check(nc_inq_varid(ncid,"lon_geo",longitudes_lgame_id))
    call nc_check(nc_inq_varid(ncid,"lat_geo_u",latitudes_lgame_wind_u_id))
    call nc_check(nc_inq_varid(ncid,"lon_geo_u",longitudes_lgame_wind_u_id))
    call nc_check(nc_inq_varid(ncid,"lat_geo_v",latitudes_lgame_wind_v_id))
    call nc_check(nc_inq_varid(ncid,"lon_geo_v",longitudes_lgame_wind_v_id))
    call nc_check(nc_get_var_double(ncid,latitudes_lgame_id,latitudes_lgame))
    call nc_check(nc_get_var_double(ncid,longitudes_lgame_id,longitudes_lgame))
    call nc_check(nc_get_var_double(ncid,latitudes_lgame_wind_u_id,latitudes_lgame_wind_u))
    call nc_check(nc_get_var_double(ncid,longitudes_lgame_wind_u_id,longitudes_lgame_wind_u))
    call nc_check(nc_get_var_double(ncid,latitudes_lgame_wind_v_id,latitudes_lgame_wind_v))
    call nc_check(nc_get_var_double(ncid,longitudes_lgame_wind_v_id,longitudes_lgame_wind_v))
    call nc_check(nc_close(ncid))
    printf("Grid file of L-GAME read.\n")
    
    ! allocating memory for the result arrays
    int (*interpolation_indices_scalar,n_avg_points) = malloc(sizeof(int(nlat*nlon,n_avg_points)))
    double (*interpolation_weights_scalar,n_avg_points) = malloc(sizeof(double(nlat*nlon,n_avg_points)))
    int (*interpolation_indices_vector_u,n_avg_points) = malloc(sizeof(int(nlat*(nlon+1),n_avg_points)))
    double (*interpolation_weights_vector_u,n_avg_points) = malloc(sizeof(double(nlat*(nlon+1),n_avg_points)))
    int (*interpolation_indices_vector_v,n_avg_points) = malloc(sizeof(int((nlat+1)*nlon,n_avg_points)))
    double (*interpolation_weights_vector_v,n_avg_points) = malloc(sizeof(double((nlat+1)*nlon,n_avg_points)))
    
    ! executing the actual interpolation
    printf("Calculating interpolation indices and weights ...\n")
    
    !$omp parallel do private(ji)
    do (int i = 0 i < nlat ++i)
      do (int j = 0 j < nlon ++j)
        double *distance_vector = malloc(no_of_points_per_layer_input_model*sizeof(double))
        for (int k = 0 k < no_of_points_per_layer_input_model ++k)
          distance_vector(k) = calculate_distance_h(&latitudes_lgame(i,j),longitudes_lgame(i,j),latitudes_input_model(k),longitudes_input_model(k),one)
        enddo
        sum_of_weights = 0._wp
        do (int k = 0 k < n_avg_points ++k)
          interpolation_indices_scalar(j*nlat + i,k) = find_min_index(distance_vector,no_of_points_per_layer_input_model)
          interpolation_weights_scalar(j*nlat + i,k) = 1._wp/(pow(distance_vector(interpolation_indices_scalar(j*nlat + i,k)),interpol_exp))
          distance_vector(interpolation_indices_scalar(j*nlat + i,k)) = 2._wp*M_PI
          sum_of_weights = sum_of_weights + interpolation_weights_scalar(j*nlat + i,k)
        enddo
        deallocate(distance_vector)
        do (int k = 0 k < n_avg_points ++k)
          interpolation_weights_scalar(j*nlat + i,k) = interpolation_weights_scalar(j*nlat + i,k)/sum_of_weights
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    !$omp parallel do private(ji)
    do (int i = 0 i < nlat ++i)
      do (int j = 0 j < nlon+1 ++j)
        double *distance_vector = malloc(no_of_points_per_layer_input_model*sizeof(double))
        do (int k = 0 k < no_of_points_per_layer_input_model ++k)
          distance_vector(k) = calculate_distance_h(&latitudes_lgame_wind_u(i,j),longitudes_lgame_wind_u(i,j),latitudes_input_model(k),longitudes_input_model(k),1._wp)
        enddo
        sum_of_weights = 0._wp
        do (int k = 0 k < n_avg_points ++k)
          interpolation_indices_vector_u(j*nlat + i,k) = find_min_index(distance_vector,no_of_points_per_layer_input_model)
          interpolation_weights_vector_u(j*nlat + i,k) = 1._wp/(pow(distance_vector(interpolation_indices_vector_u(j*nlat + i,k)),interpol_exp))
          distance_vector(interpolation_indices_vector_u(j*nlat + i,k)) = 2._wp*M_PI
          sum_of_weights = sum_of_weights + interpolation_weights_vector_u(j*nlat + i,k)
        enddo
        deallocate(distance_vector)
        do (int k = 0 k < n_avg_points ++k)
          interpolation_weights_vector_u(j*nlat + i,k) = interpolation_weights_vector_u(j*nlat + i,k)/sum_of_weights
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    !$omp parallel do private(ji)
    do (int i = 0 i < nlat+1 ++i)
      do (int j = 0 j < nlon ++j)
      double *distance_vector = malloc(no_of_points_per_layer_input_model*sizeof(double))
        do (int k = 0 k < no_of_points_per_layer_input_model ++k)
          distance_vector(k) = calculate_distance_h(&latitudes_lgame_wind_v(i,j),longitudes_lgame_wind_v(i,j),latitudes_input_model(k),longitudes_input_model(k),1._wp)
        enddo
        double sum_of_weights = 0._wp
        do (int k = 0 k < n_avg_points ++k)
          interpolation_indices_vector_v(j*(nlat+1) + i,k) = find_min_index(distance_vector,no_of_points_per_layer_input_model)
          interpolation_weights_vector_v(j*(nlat+1) + i,k) = 1._wp/(pow(distance_vector(interpolation_indices_vector_v(j*(nlat+1) + i,k)),interpol_exp))
          distance_vector(interpolation_indices_vector_v(j*(nlat+1) + i,k)) = 2._wp*M_PI
          sum_of_weights += interpolation_weights_vector_v(j*(nlat+1) + i,k)
        enddo
        deallocate(distance_vector)
        do (int k = 0 k < n_avg_points ++k)
          interpolation_weights_vector_v(j*(nlat+1) + i,k) = interpolation_weights_vector_v(j*(nlat+1) + i,k)/sum_of_weights
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    write(*,*) "Calculating interpolation indices and weights finished."
    
    ! freeing memory we do not need further
    deallocate(latitudes_input_model)
    deallocate(longitudes_input_model)
    deallocate(latitudes_lgame)
    deallocate(longitudes_lgame)
    deallocate(latitudes_lgame_wind_u)
    deallocate(longitudes_lgame_wind_u)
    deallocate(latitudes_lgame_wind_v)
    deallocate(longitudes_lgame_wind_v)
    
    ! writing the result to a NetCDF file
    char output_file_pre(200)
    sprintf(output_file_pre,"%s/interpolation_files/icon-d22lgame_%s",real2game_root_dir,lgame_grid)
    char output_file(strlen(output_file_pre) + 1)
    strcpy(output_file,output_file_pre)
    printf("Starting to write to output file ...\n")
    int interpolation_indices_scalar_id,interpolation_weights_scalar_id,interpolation_indices_vector_u_id,interpolation_weights_vector_u_id,
    interpolation_indices_vector_v_id,interpolation_weights_vector_v_id,scalar_dimid,u_dimid,v_dimid,avg_dimid
    int dim_vector(2)
    call nc_check(nc_create(output_file,NC_CLOBBER,ncid))
    call nc_check(nc_def_dim(ncid,"scalar_index",nlat*nlon,scalar_dimid))
    call nc_check(nc_def_dim(ncid,"u_index",nlat*(nlon+1),u_dimid))
    call nc_check(nc_def_dim(ncid,"v_index",(nlat+1)*nlon,v_dimid))
    call nc_check(nc_def_dim(ncid,"interpol_index",n_avg_points,avg_dimid))
    dim_vector(1) = scalar_dimid
    dim_vector(2) = avg_dimid
    call nc_check(nc_def_var(ncid,"interpolation_indices_scalar",NC_INT,2,dim_vector,interpolation_indices_scalar_id))
    call nc_check(nc_def_var(ncid,"interpolation_weights_scalar",NC_DOUBLE,2,dim_vector,interpolation_weights_scalar_id))
    dim_vector(1) = u_dimid
    call nc_check(nc_def_var(ncid,"interpolation_indices_vector_u",NC_INT,2,dim_vector,interpolation_indices_vector_u_id))
    call nc_check(nc_def_var(ncid,"interpolation_weights_vector_u",NC_DOUBLE,2,dim_vector,interpolation_weights_vector_u_id))
    dim_vector(1) = v_dimid
    call nc_check(nc_def_var(ncid,"interpolation_indices_vector_v",NC_INT,2,dim_vector,interpolation_indices_vector_v_id))
    call nc_check(nc_def_var(ncid,"interpolation_weights_vector_v",NC_DOUBLE,2,dim_vector,interpolation_weights_vector_v_id))
    call nc_check(nc_enddef(ncid))
    call nc_check(nc_put_var_int(ncid,interpolation_indices_scalar_id,interpolation_indices_scalar))
    call nc_check(nc_put_var_double(ncid,interpolation_weights_scalar_id,interpolation_weights_scalar))
    call nc_check(nc_put_var_int(ncid,interpolation_indices_vector_u_id,interpolation_indices_vector_u))
    call nc_check(nc_put_var_double(ncid,interpolation_weights_vector_u_id,interpolation_weights_vector_u))
    call nc_check(nc_put_var_int(ncid,interpolation_indices_vector_v_id,interpolation_indices_vector_v))
    call nc_check(nc_put_var_double(ncid,interpolation_weights_vector_v_id,interpolation_weights_vector_v))
    call nc_check(nc_close(ncid))
    
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







