! This source file is part of real2GAME,which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/real2GAME

program formatter

  ! This tool reads the output from other models / data assimilation systems and brings it into a standardized format.

  use netcdf
  use eccodes
  use mo_shared, only: wp

  implicit none

  integer, parameter    :: n_levels_input = 12
  integer               :: ji,ncid,h_dimid,v_dimid,z_surf_id,sp_id,sst_dimid,t_id,spec_hum_id,z_id,u_id,v_id, &
                           lat_sst_id,lon_sst_id,sst_id,dim_vector(2),levels_vector(n_levels_input)
  real(wp), allocatable :: z_height_amsl_one_layer(:),temperature_one_layer(:),spec_hum_one_layer(:), &
                           u_one_layer(:),v_one_layer(:),z_height_amsl(:,:),temperature(:,:),spec_hum(:,:), &
                           u_wind(:,:),v_wind(:,:)
  character(len=4)      :: year_string
  character(len=2)      :: month_string,day_string,hour_string
  character(len=128)    :: real2game_root_dir

  ! defining the levels of the model we want to use
  levels_vector(1) = 1
  levels_vector(2) = 10
  levels_vector(3) = 19
  levels_vector(4) = 27
  levels_vector(5) = 35
  levels_vector(6) = 43
  levels_vector(7) = 51
  levels_vector(8) = 59
  levels_vector(9) = 67
  levels_vector(10) = 75
  levels_vector(11) = 83
  levels_vector(12) = 90
  
  ! shell arguments
  call get_command_argument(1,year_string)
  call get_command_argument(2,month_string)
  call get_command_argument(3,day_string)
  call get_command_argument(4,hour_string)
  call get_command_argument(5,real2game_root_dir)
  
  ! single-layer arrays from grib
  allocate(z_height_amsl_one_layer(n_points_per_layer_input))
  allocate(temperature_one_layer(n_points_per_layer_input))
  allocate(spec_hum_one_layer(n_points_per_layer_input))
  allocate(u_one_layer(n_points_per_layer_input))
  allocate(v_one_layer(n_points_per_layer_input))
  
  ! 2D-arrays for NetCDF
  allocate(z_height_amsl(n_points_per_layer_input,n_levels_input))
  allocate(temperature(n_points_per_layer_input,n_levels_input))
  allocate(spec_hum(n_points_per_layer_input,n_levels_input))
  allocate(u_wind(n_points_per_layer_input,n_levels_input))
  allocate(v_wind(n_points_per_layer_input,n_levels_input))
  
  ! grib stuff
  FILE *ecc_file
  size_t no_of_points_per_layer_input_size_t
  int err
  codes_handle *handle = NULL
  
  ! reading the data from the free atmosphere
  ! loop over all relevant levels in the free atmosphere
  do level_index=1,n_levels_input
    ! vertical position of the current layer
    z_inpput_model_file = real2game_root_dir // "/input/icon_global_icosahedral_time-invariant_" // year_string &
                          // month_string // day_string // hour_string // "_" &
                          // trim(int2string(levels_vector(level_index))) // "_HHL.grib2"
    
    ecc_file = fopen(z_input_model_file,"r")
    handle = codes_handle_new_from_file(NULL,ecc_file,PRODUCT_GRIB,err)
    if (err != 0) ECCERR(err)
    no_of_points_per_layer_input_size_t = (size_t) n_points_per_layer_input
    ECCCHECK(codes_get_double_array(handle,"values",z_height_amsl_one_layer,no_of_points_per_layer_input_size_t))
    codes_handle_delete(handle)
    fclose(ecc_file)
    
    !$omp parallel do private(ji)
    do ji=1,n_points_per_layer_input
      z_height_amsl(ji,level_index) = z_height_amsl_one_layer(ji)
    enddo
    !$omp end parallel do
    
    ! reading the temperature
    temperature_file = real2game_root_dir // "/input/icon_global_icosahedral_model-level_" year_string // month_string // day_string &
                       // hour_string "_000_" // int2string(trim(levels_vector(level_index))) // "_T.grib2"
    
    ecc_file = fopen(temperature_file,"r")
    handle = codes_handle_new_from_file(NULL,ecc_file,PRODUCT_GRIB,err)
    if (err != 0) ECCERR(err)
    ECCCHECK(codes_get_double_array(handle,"values",temperature_one_layer(0),no_of_points_per_layer_input_size_t))
    codes_handle_delete(handle)
    fclose(ecc_file)
    
    !$omp parallel do private(ji)
    do ji=1,n_points_per_layer_input
      temperature(ji,level_index) = temperature_one_layer(ji)
    enddo
    !$omp end parallel do
    
    ! reading the specific humidity
    spec_hum_file = real2game_root_dir // "/input/icon_global_icosahedral_model-level_" // year_string // month_string // &
                    month_string // day_string // hour_string // "_000_" // trim(int2string(levels_vector(level_index))) // &
                    "_QV.grib2"
    
    ecc_file = fopen(spec_hum_file,"r")
    handle = codes_handle_new_from_file(NULL,ecc_file,PRODUCT_GRIB,err)
    if (err != 0) ECCERR(err)
    ECCCHECK(codes_get_double_array(handle,"values",spec_hum_one_layer,no_of_points_per_layer_input_size_t))
    codes_handle_delete(handle)
    fclose(ecc_file)
    
    !$omp parallel do private(ji)
    do ji=1,n_points_per_layer_input
      spec_hum(ji,level_index) = spec_hum_one_layer(ji)
    enddo
    !$omp end parallel do
    
    ! reading the u wind
    u_wind_file = real2game_root_dir // "/input/icon_global_icosahedral_model-level_" // year_string // month_string // &
                  day_string // hour_string // "_000_" // trim(int2string(levels_vector(level_index))) & 
                  // "_U.grib2",
    
    ecc_file = fopen(u_wind_file,"r")
    handle = codes_handle_new_from_file(NULL,ecc_file,PRODUCT_GRIB,err)
    if (err != 0) ECCERR(err)
    ECCCHECK(codes_get_double_array(handle,"values",u_one_layer,no_of_points_per_layer_input_size_t))
    codes_handle_delete(handle)
    fclose(ecc_file)
    
    !$omp parallel do private(ji)
    do ji=1,n_points_per_layer_input
      u_wind(ji,level_index) = u_one_layer(ji)
    enddo
    !$omp end parallel do
    
    ! reading the v wind
    v_wind_file = real2game_root_dir // "/input/icon_global_icosahedral_model-level_" // year_string // &
                  month_string // day_string // hour_string // "%s_000_" // trim(int2string(levels_vector(level_index))) &
                  // "_V.grib2",
    
    ecc_file = fopen(v_wind_file,"r")
    handle = codes_handle_new_from_file(NULL,ecc_file,PRODUCT_GRIB,err)
    if (err != 0) ECCERR(err)
    ECCCHECK(codes_get_double_array(handle,"values",v_one_layer,no_of_points_per_layer_input_size_t))
    codes_handle_delete(handle)
    fclose(ecc_file)
    
    !$omp parallel do private(ji)
    do (int i = 0 i < n_points_per_layer_input ++i)
      v_wind(ji,level_index) = v_one_layer(ji)
    enddo
    !$omp end parallel do
  
  enddo
  
  deallocate(z_height_amsl_one_layer)
  deallocate(temperature_one_layer)
  deallocate(spec_hum_one_layer)
  deallocate(u_one_layer)
  deallocate(v_one_layer)
  
  ! reading the surface height
  double *surface_height = malloc(n_points_per_layer_input)
  sfc_obs_file = real2game_root_dir // "/input/icon_global_icosahedral_time-invariant_" // year_string // month_string // &
                 day_string // hour_string // "%s_HSURF.grib2"
  
  ecc_file = fopen(sfc_obs_file,"r")
  handle = codes_handle_new_from_file(NULL,ecc_file,PRODUCT_GRIB,err)
  if (err != 0) ECCERR(err)
  ECCCHECK(codes_get_double_array(handle,"values",surface_height,no_of_points_per_layer_input_size_t))
  codes_handle_delete(handle)
  fclose(ecc_file)
  
  ! reading the surface presure
  double *pressure_surface = malloc(n_points_per_layer_input)
  
  sfc_pres_file = real2game_root_dir // "/input/icon_global_icosahedral_single-level_" // year_string // month_string // &
                  day_string // hour_string // "_000_PS.grib2"
  
  ecc_file = fopen(sfc_pres_file,"r")
  handle = codes_handle_new_from_file(NULL,ecc_file,PRODUCT_GRIB,err)
  if (err != 0) ECCERR(err)
  ECCCHECK(codes_get_double_array(handle,"values",pressure_surface,no_of_points_per_layer_input_size_t))
  codes_handle_delete(handle)
  fclose(ecc_file)
  
  ! reading the SST
  double *latitudes_sst = malloc(N_SST_POINTS)
  double *longitudes_sst = malloc(N_SST_POINTS)
  double *sst = malloc(N_SST_POINTS)
  
  sst_file = real2game_root_dir // "/input/rtgssthr_grb_0.5.grib2"
  
  ecc_file = fopen(sst_file,"r")
  handle = codes_handle_new_from_file(NULL,ecc_file,PRODUCT_GRIB,err)
  if (err /= 0) ECCERR(err)
  size_t no_of_sst_points_size_t = (size_t) N_SST_POINTS
  ECCCHECK(codes_get_double_array(handle,"values",sst,no_of_sst_points_size_t))
  ECCCHECK(codes_get_double_array(handle,"latitudes",latitudes_sst,no_of_sst_points_size_t))
  ECCCHECK(codes_get_double_array(handle,"longitudes",longitudes_sst,no_of_sst_points_size_t))
  codes_handle_delete(handle)
  fclose(ecc_file)
  
  ! transforming the coordinates of the SST grid from degrees to radians
  !$omp parallel do private(ji)
  do ji=1,n_sst_oints
    latitudes_sst(ji) = 2._wp*M_PI*latitudes_sst(ji)/360._wp
    longitudes_sst(ji) = 2._wp*M_PI*longitudes_sst(ji)/360._wp
  enddo
  !$omp end parallel do
    
  ! Writing the observations to a netcdf file.
  output_file = real2game_root_dir // "/input/obs_" // year_string // month_string // day_string // &
                hour_string // "%s%s%s.nc"
    
  nc_check(nc_create(output_file,NC_CLOBBER,ncid))
  ! Defining the dimensions.
  nc_check(nc_def_dim(ncid,"h_index",n_points_per_layer_input,h_dimid))
  nc_check(nc_def_dim(ncid,"v_index",n_levels_input,v_dimid))
  nc_check(nc_def_dim(ncid,"sst_index",N_SST_POINTS,sst_dimid))
  dim_vector(1) = h_dimid
  dim_vector(2) = v_dimid
  nc_check(nf90_def_var(ncid,"z_surface",NC_DOUBLE,1,h_dimid,z_surf_id))
  nc_check(nf90_def_var(ncid,"pressure_surface",NC_DOUBLE,1,h_dimid,sp_id))
  nc_check(nf90_def_var(ncid,"z_height",NC_DOUBLE,2,dim_vector,z_id))
  nc_check(nf90_def_var(ncid,"temperature",NC_DOUBLE,2,dim_vector,t_id))
  nc_check(nf90_def_var(ncid,"spec_humidity",NC_DOUBLE,2,dim_vector,spec_hum_id))
  nc_check(nf90_def_var(ncid,"u_wind",NC_DOUBLE,2,dim_vector,u_id))
  nc_check(nf90_def_var(ncid,"v_wind",NC_DOUBLE,2,dim_vector,v_id))
  nc_check(nf90_def_var(ncid,"lat_sst",NC_DOUBLE,1,sst_dimid,lat_sst_id))
  nc_check(nf90_def_var(ncid,"lon_sst",NC_DOUBLE,1,sst_dimid,lon_sst_id))
  nc_check(nf90_def_var(ncid,"sst",NC_DOUBLE,1,sst_dimid,sst_id))
  nc_check(nc_enddef(ncid))
  ! Setting the variables.
  nc_check(nf90_put_var(ncid,z_surf_id,surface_height))
  nc_check(nf90_put_var(ncid,sp_id,pressure_surface))
  nc_check(nf90_put_var(ncid,z_id,z_height_amsl))
  nc_check(nf90_put_var(ncid,t_id,temperature))
  nc_check(nf90_put_var(ncid,spec_hum_id,spec_hum))
  nc_check(nf90_put_var(ncid,u_id,u_wind))
  nc_check(nf90_put_var(ncid,v_id,v_wind))
  nc_check(nf90_put_var(ncid,lat_sst_id,latitudes_sst))
  nc_check(nf90_put_var(ncid,lon_sst_id,longitudes_sst))
  nc_check(nf90_put_var(ncid,sst_id,sst))
  nc_check(nc_close(ncid))
    
  ! freeing the memory
  deallocate(z_height_amsl)
  deallocate(surface_height)
  deallocate(pressure_surface)
  deallocate(u_wind)
  deallocate(v_wind)
  deallocate(temperature)
  deallocate(spec_hum)
  deallocate(sst)
  deallocate(latitudes_sst)
  deallocate(longitudes_sst)

end program formatter










