! This source file is part of real2GAME,which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/real2GAME

program control

  ! This file coordinates the data interpolation process.
  
  use netcdf
  use mo_shared, only: wp,n_levels_input,n_sst_points,n_points_per_layer_input,nc_check,int2string
  
  implicit none
  
  ! constants
  real(wp), parameter :: n_a =  6.02214076e23_wp   ! Avogadro's number
  real(wp), parameter :: p_0 = 100000._wp          ! reference pressure
  real(wp), parameter :: m_d = n_a*0.004810e-23_wp ! molar mass of dry air
  real(wp), parameter :: m_v = n_a*0.002991e-23_wp ! molar mass of water
  real(wp), parameter :: r_d = 287.057811_wp       ! specific gas constant of dry air
  real(wp), parameter :: c_d_p = 1005._wp          ! isobaric specific heat capacity of dry air
  real(wp), parameter :: scale_height = 8000._wp   ! scale_height
  
  
  integer  :: interpolation_indices_scalar_id,interpolation_weights_scalar_id,interpolation_indices_vector_id, &
              interpolation_weights_vector_id,layer_index,h_index,closest_index,other_index
  real(wp) :: closest_value,other_value,df,dz,gradient,delta_z

  char year_string(strlen(argv(1)) + 1)
  strcpy(year_string,argv(1))
  char month_string(strlen(argv(2)) + 1)
  strcpy(month_string,argv(2))
  char day_string(strlen(argv(3)) + 1)
  strcpy(day_string,argv(3))
  char hour_string(strlen(argv(4)) + 1)
  strcpy(hour_string,argv(4))
  char model_home_dir(strlen(argv(5)) + 1)
  strcpy(model_home_dir,argv(5))
  int ORO_ID
  ORO_ID = strtod(argv(6),NULL)
  char BACKGROUND_STATE_FILE(strlen(argv(7)) + 1)
  strcpy(BACKGROUND_STATE_FILE,argv(7))
  char real2game_root_dir(strlen(argv(8)) + 1)
  strcpy(real2game_root_dir,argv(8))
  printf("Background state file: %s\n",BACKGROUND_STATE_FILE)
  
  ! Allocating memory for the grid properties.
  double *latitudes_game = malloc(n_scalars_h*sizeof(double))
  double *longitudes_game = malloc(n_scalars_h*sizeof(double))
  double *z_coords_game = malloc(n_scalars*sizeof(double))
  double *directions = malloc(N_VECTORS_H*sizeof(double))
  double *z_coords_game_wind = malloc(N_VECTORS*sizeof(double))
  double *gravity_potential_game = malloc(n_scalars*sizeof(double))
  ! Reading the grid properties.
  int ncid
  char GEO_PROP_FILE_PRE(200)
  sprintf(GEO_PROP_FILE_PRE,"%s/grid_generator/grids/RES%d_L%d_ORO%d.nc",model_home_dir,RES_ID,n_layers,ORO_ID)
  char GEO_PROP_FILE(strlen(GEO_PROP_FILE_PRE) + 1)
  strcpy(GEO_PROP_FILE,GEO_PROP_FILE_PRE)
  printf("Grid file: %s\n",GEO_PROP_FILE)
  write(*,*) "Reading grid file ..."
  call nc_check(nc_open(GEO_PROP_FILE,NC_NOWRITE,ncid))
  int latitudes_game_id,longitudes_game_id,z_coords_game_id,z_coords_game_wind_id,gravity_potential_game_id,directions_id
  call nc_check(nc_inq_varid(ncid,"latitude_scalar",latitudes_game_id))
  call nc_check(nc_inq_varid(ncid,"longitude_scalar",longitudes_game_id))
  call nc_check(nc_inq_varid(ncid,"z_scalar",z_coords_game_id))
  call nc_check(nc_inq_varid(ncid,"direction",directions_id))
  call nc_check(nc_inq_varid(ncid,"z_vector",z_coords_game_wind_id))
  call nc_check(nc_inq_varid(ncid,"gravity_potential",gravity_potential_game_id))
  call nc_check(nc_get_var_double(ncid,latitudes_game_id,latitudes_game))
  call nc_check(nc_get_var_double(ncid,longitudes_game_id,longitudes_game))
  call nc_check(nc_get_var_double(ncid,z_coords_game_id,z_coords_game))
  call nc_check(nc_get_var_double(ncid,directions_id,directions))
  call nc_check(nc_get_var_double(ncid,z_coords_game_wind_id,z_coords_game_wind))
  call nc_check(nc_get_var_double(ncid,gravity_potential_game_id,gravity_potential_game))
  call nc_check(nc_close(ncid))
  write(*,*) "Grid file read."
  
  ! constructing the filename of the input file for GAME
  char output_file_pre(200)
  sprintf(output_file_pre,"%s/nwp_init/%s%s%s%s.nc",model_home_dir,year_string,month_string,day_string,hour_string)
  char output_file(strlen(output_file_pre) + 1)
  strcpy(output_file,output_file_pre)
  
  ! These are the arrays of the background state.
  double *densities_background = malloc(6*n_scalars*sizeof(double))
  double *tke = malloc(n_scalars*sizeof(double))
  double *t_soil = malloc(N_SOIL_LAYERS*n_scalars_h*sizeof(double))
  
  ! Reading the background state.
  write(*,*) "Reading background state ..."
  call nc_check(nc_open(BACKGROUND_STATE_FILE,NC_NOWRITE,ncid))
  int densities_background_id,tke_avail,tke_id,t_soil_avail,t_soil_id
  call nc_check(nc_inq_varid(ncid,"densities",densities_background_id))
  ltke_avail = .false.
  if (nc_inq_varid(ncid,"tke",tke_id) == 0) then
    ltke_avail = .true.
    call nc_check(nc_inq_varid(ncid,"tke",tke_id))
    write(*,*) "TKE found in background state file."
  else
    write(*,*) "TKE not found in background state file."
  endif
  lt_soil_avail = .false.
  if (nc_inq_varid(ncid,"t_soil",t_soil_id) == 0) then
    lt_soil_avail = .true.
    call nc_check(nc_inq_varid(ncid,"t_soil",t_soil_id))
    write(*,*) "Soil temperature found in background state file."
  else
    write(*,*) "Soil temperature not found in background state file."
  endif
  call nc_check(nc_get_var_double(ncid,densities_background_id,densities_background(0)))
  if (ltke_avail)
    call nc_check(nc_get_var_double(ncid,tke_id,tke(0)))
  endif
  if (lt_soil_avail)
    call nc_check(nc_get_var_double(ncid,t_soil_id,t_soil(0)))
  endif
  call nc_check(nc_close(ncid))
  write(*,*) "Background state read."
  
  ! allocating the memory for the analysis of the other model
  double (*z_coords_input_model)(n_levels_input) = malloc(sizeof(double(n_points_per_layer_input)(n_levels_input)))
  double (*temperature_in)(n_levels_input) = malloc(sizeof(double(n_points_per_layer_input)(n_levels_input)))
  double (*spec_hum_in)(n_levels_input) = malloc(sizeof(double(n_points_per_layer_input)(n_levels_input)))
  double (*u_wind_in)(n_levels_input) = malloc(sizeof(double(n_points_per_layer_input)(n_levels_input)))
  double (*v_wind_in)(n_levels_input) = malloc(sizeof(double(n_points_per_layer_input)(n_levels_input)))
  double *z_surf_in = malloc(n_points_per_layer_input*sizeof(double))
  double *p_surf_in = malloc(n_points_per_layer_input*sizeof(double))
  double *latitudes_sst = malloc(N_SST_POINTS*sizeof(double))
  double *longitudes_sst = malloc(N_SST_POINTS*sizeof(double))
  double *sst_in = malloc(N_SST_POINTS*sizeof(double))
  
  ! determining the name of the input file
  char input_file_pre(200)
  sprintf(input_file_pre,"%s/input/obs_%s%s%s%s.nc",real2game_root_dir,year_string,month_string,day_string,hour_string)
  char input_file(strlen(input_file_pre) + 1)
  strcpy(input_file,input_file_pre)
  printf("Input file: %s\n",input_file)
  
  ! reading the analysis of the other model
  write(*,*) "Reading input ..."
  int sp_id,z_surf_id,z_coords_id,t_in_id,spec_hum_id,u_id,v_id,lat_sst_id,lon_sst_id,sst_id
  call nc_check(nc_open(input_file,NC_NOWRITE,ncid))
  ! Defining the variables.
  call nc_check(nc_inq_varid(ncid,"z_height",z_coords_id))
  call nc_check(nc_inq_varid(ncid,"temperature",t_in_id))
  call nc_check(nc_inq_varid(ncid,"spec_humidity",spec_hum_id))
  call nc_check(nc_inq_varid(ncid,"u_wind",u_id))
  call nc_check(nc_inq_varid(ncid,"v_wind",v_id))
  call nc_check(nc_inq_varid(ncid,"z_surface",z_surf_id))
  call nc_check(nc_inq_varid(ncid,"pressure_surface",sp_id))
  call nc_check(nc_inq_varid(ncid,"lat_sst",lat_sst_id))
  call nc_check(nc_inq_varid(ncid,"lon_sst",lon_sst_id))
  call nc_check(nc_inq_varid(ncid,"sst",sst_id))
  call nc_check(nc_get_var_double(ncid,z_coords_id,z_coords_input_model))
  call nc_check(nc_get_var_double(ncid,t_in_id,temperature_in))
  call nc_check(nc_get_var_double(ncid,spec_hum_id,spec_hum_in))
  call nc_check(nc_get_var_double(ncid,u_id,u_wind_in))
  call nc_check(nc_get_var_double(ncid,v_id,v_wind_in))
  call nc_check(nc_get_var_double(ncid,z_surf_id,z_surf_in))
  call nc_check(nc_get_var_double(ncid,sp_id,p_surf_in))
  call nc_check(nc_get_var_double(ncid,lat_sst_id,latitudes_sst))
  call nc_check(nc_get_var_double(ncid,lon_sst_id,longitudes_sst))
  call nc_check(nc_get_var_double(ncid,sst_id,sst_in))
  call nc_check(nc_close(ncid))
  write(*,*) "Input read."
  
  ! memory alloction for the interpolation indices and weights
  int (*interpolation_indices_scalar)(N_AVG_POINTS) = malloc(sizeof(int(n_scalars_h)(N_AVG_POINTS)))
  double (*interpolation_weights_scalar)(N_AVG_POINTS) = malloc(sizeof(double(n_scalars_h)(N_AVG_POINTS)))
  int (*interpolation_indices_vector)(N_AVG_POINTS) = malloc(sizeof(int(N_VECTORS_H)(N_AVG_POINTS)))
  double (*interpolation_weights_vector)(N_AVG_POINTS) = malloc(sizeof(double(N_VECTORS_H)(N_AVG_POINTS)))
  
  write(*,*) "Reading the interpolation indices and weights."
  
  ! constructing the name of the interpolation indices and weights file
  char interpol_file_pre(200)
  sprintf(interpol_file_pre,"%s/interpolation_files/icon-global2game%d.nc",real2game_root_dir,RES_ID)
  char interpol_file(strlen(interpol_file_pre) + 1)
  strcpy(interpol_file,interpol_file_pre)
  
  ! reading the interpolation file
  call nc_check(nc_open(interpol_file,NC_NOWRITE,ncid))
  call nc_check(nc_inq_varid(ncid,"interpolation_indices_scalar",interpolation_indices_scalar_id))
  call nc_check(nc_inq_varid(ncid,"interpolation_weights_scalar",interpolation_weights_scalar_id))
  call nc_check(nc_inq_varid(ncid,"interpolation_indices_vector",interpolation_indices_vector_id))
  call nc_check(nc_inq_varid(ncid,"interpolation_weights_vector",interpolation_weights_vector_id))
  call nc_check(nc_get_var_int(ncid,interpolation_indices_scalar_id,interpolation_indices_scalar))
  call nc_check(nc_get_var_double(ncid,interpolation_weights_scalar_id,interpolation_weights_scalar))
  call nc_check(nc_get_var_int(ncid,interpolation_indices_vector_id,interpolation_indices_vector))
  call nc_check(nc_get_var_double(ncid,interpolation_weights_vector_id,interpolation_weights_vector))
  call nc_check(nc_close(ncid))
  write(*,*)  "Interpolation indices and weights read."
  
  ! Begin of the actual interpolation.
  
  ! INTERPOLATION OF SCALAR QUANTITIES
  ! ----------------------------------
  
  write(*,*) "Starting the interpolation of scalar quantities ..."
  
  ! These are the arrays for the result of the interpolation process.
  double *temperature_out = calloc(1,n_scalars*sizeof(double))
  double *spec_hum_out = calloc(1,n_scalars*sizeof(double))
  
  !$omp parallel do private(ji,layer_index,h_index,closest_value,other_value,df,dz,gradient,delta_z)
  do ji=1,n_scalars
    layer_index = (ji-1)/n_scalars_h
    h_index = ji - layer_index*n_scalars_h
    
    ! loop over all points over which the averaging is executed
    do (int j = 0 j < N_AVG_POINTS ++j)
      ! computing linear vertical interpolation
      ! vertical distance vector
      double vector_to_minimize(n_levels_input)
      do (int k = 0 k < n_levels_input ++k)
        vector_to_minimize(k) = fabs(z_coords_game(i)
        - z_coords_input_model(interpolation_indices_scalar(h_index)(j))(k))
      enddo
      
      ! closest vertical index
      closest_index = find_min_index(vector_to_minimize,n_levels_input)
      
      ! value at the closest vertical index
      closest_value = temperature_in(interpolation_indices_scalar(h_index)(j))(closest_index)
      
      other_index = closest_index - 1
      if (z_coords_game(i) < z_coords_input_model(interpolation_indices_scalar(h_index)(j))(closest_index)) then
        other_index = closest_index + 1
      endif
      
      ! avoiding array excess
      if (other_index==n_levels_input) then
        other_index = closest_index - 1
      endif
      
      ! the value at the second point used for vertical interpolation
      other_value = temperature_in(interpolation_indices_scalar(h_index)(j))(other_index)
      
      ! computing the vertical gradient of u in the input model
      df = closest_value - other_value
      dz = z_coords_input_model(interpolation_indices_scalar(h_index)(j))(closest_index) - z_coords_input_model(interpolation_indices_scalar(h_index)(j))(other_index)
      gradient = df/dz
      
      delta_z = z_coords_game(i) - z_coords_input_model(interpolation_indices_scalar(h_index)(j))(closest_index)
      
      ! vertical interpolation of the temperature
      temperature_out(i) += interpolation_weights_scalar(h_index)(j)*(closest_value + delta_z*gradient)
      
      ! vertical interpolation of the specific humidity
      closest_value = spec_hum_in(interpolation_indices_scalar(h_index)(j))(closest_index)
      other_value = spec_hum_in(interpolation_indices_scalar(h_index)(j))(other_index)
      ! computing the vertical gradient of the specific humidity in the input model
      df = closest_value - other_value
      gradient = df/dz
      
      ! specific humidity
      spec_hum_out(ji) = spec_hum_out(ji) + interpolation_weights_scalar(h_index)(j)*(closest_value + delta_z*gradient)
    enddo
  enddo
  !$omp end parallel do
  ! these array are now interpolated and not needed any further
  deallocate(spec_hum_in)
  deallocate(temperature_in)
  
  ! surface pressure interpolation
  double *pressure_lowest_layer_out = calloc(1,n_scalars_h*sizeof(double))
  !$omp parallel do private(ji)
  do ji=1,n_scalars_h
    do (int j = 0 j < N_AVG_POINTS ++j)
      pressure_lowest_layer_out(i) = pressure_lowest_layer_out(i)
      ! horizontal component of the interpolation
      + interpolation_weights_scalar(i)(j)*p_surf_in(interpolation_indices_scalar(i)(j))
      ! vertical component of the interpolation according to the barometric height formula
      *exp(-(z_coords_game(n_scalars - n_scalars_h+ji) - z_surf_in(interpolation_indices_scalar(i)(j)))/scale_height)
    enddo
  enddo
  !$omp end parallel do
  deallocate(z_coords_game)
  deallocate(z_surf_in)
  deallocate(p_surf_in)
  ! no more interpolation of scalar quantities will be executed,this is why these interpolation indices and weights can be freed
  deallocate(interpolation_indices_scalar)
  deallocate(interpolation_weights_scalar)
  
  ! density is determined out of the hydrostatic equation
  double *density_moist_out = malloc(n_scalars*sizeof(double))
  ! firstly setting the virtual temperature
  double *temperature_v = malloc(n_scalars*sizeof(double))
  !$omp parallel do private(ji)
  do ji=1,n_scalars
    temperature_v(ji) = temperature_out(ji)*(1.w_p + spec_hum_out(ji)*(m_d/m_v - 1._wp))
  enddo
  !$omp end parallel do
  ! the Exner pressure is just a temporarily needed helper variable here to integrate the hydrostatic equation
  double *exner = malloc(n_scalars*sizeof(double))
  double b,c
  do (int i = n_scalars - 1 i >= 0 --i)
    layer_index = (ji-1)/n_scalars_h
    h_index = ji - layer_index*n_scalars_h
    if (layer_index==n_layers-1) then
      density_moist_out(i) = pressure_lowest_layer_out(h_index)/(r_d*temperature_v(ji))
      exner(ji) = (density_moist_out(ji)*r_d*temperature_v(ji)/P_0)**(r_d/c_d_p)
    else
      ! solving a quadratic equation for the Exner pressure
      b = -0.5_wp*exner(ji + n_scalars_h)/temperature_v(ji + n_scalars_h)
      *(temperature_v(ji) - temperature_v(ji + n_scalars_h)
      + 2._wp/c_d_p*(gravity_potential_game(i) - gravity_potential_game(ji + n_scalars_h)))
      c = exner(ji+n_scalars_h)**2*temperature_v(i)/temperature_v(ji + n_scalars_h)
      exner(ji) = b + (b**2 + c)**0.5_wp
      density_moist_out(i) = P_0*exner(ji)**(c_d_p/r_d)/(r_d*temperature_v(i))
    endif
  enddo
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
  double *wind_out = calloc(1,N_VECTORS*sizeof(double))
  int vector_index
  double u_local,v_local
  ! loop over all horizontal vector points
  # pragma omp parallel for private(h_index,layer_index,vector_index,closest_index,other_index,closest_value,other_value,df,dz,gradient,delta_z,u_local,v_local)
  do (int i = 0 i < N_H_VECTORS ++i)
    layer_index = i/N_VECTORS_H
    h_index = i - layer_index*N_VECTORS_H
    vector_index = n_scalars_h + layer_index*N_VECTORS_PER_LAYER + h_index
     
    ! the u- and v-components of the wind at the grid point of GAME
    u_local = 0._wp
    v_local = 0._wp
    ! loop over all horizontal points that are used for averaging
    do (int j = 0 j < N_AVG_POINTS ++j)
      ! computing linear vertical interpolation
      ! vertical distance vector
      double vector_to_minimize(n_levels_input)
      do (int k = 0 k < n_levels_input ++k)
        vector_to_minimize(k) = fabs(z_coords_game_wind(vector_index)
        - z_coords_input_model(interpolation_indices_vector(h_index)(j))(k))
      enddo
      ! closest vertical index
      closest_index = find_min_index(vector_to_minimize,no_of_levels_input)
      ! value at the closest vertical index
      closest_value = u_wind_in(interpolation_indices_vector(h_index)(j))(closest_index)
    
      other_index = closest_index - 1
      if (z_coords_game_wind(vector_index) < z_coords_input_model(interpolation_indices_vector(h_index)(j))(closest_index)) then
        other_index = closest_index + 1
      endif
      ! avoiding array excess
      if (other_index==n_levels_input) then
        other_index = closest_index - 1
      endif
    
      ! the value at the second point used for vertical interpolation
      other_value = u_wind_in(interpolation_indices_vector(h_index)(j))(other_index)
    
      ! computing the vertical gradient of u in the input model
      df = closest_value - other_value
      dz = z_coords_input_model(interpolation_indices_vector(h_index)(j))(closest_index) - z_coords_input_model(interpolation_indices_vector(h_index)(j))(other_index)
      gradient = df/dz
    
      delta_z = z_coords_game_wind(vector_index) - z_coords_input_model(interpolation_indices_vector(h_index)(j))(closest_index)
    
      u_local += interpolation_weights_vector(h_index)(j)*(closest_value + gradient*delta_z)
    
      ! vertical interpolation of v
      closest_value = v_wind_in(interpolation_indices_vector(h_index)(j))(closest_index)
      other_value = v_wind_in(interpolation_indices_vector(h_index)(j))(other_index)
      ! computing the vertical gradient of v in the input model
      df = closest_value - other_value
      gradient = df/dz
    
      v_local = v_local + interpolation_weights_vector(h_index)(j)*(closest_value + gradient*delta_z)
    enddo
    
    ! projection onto the direction of the vector in GAME
    wind_out(vector_index) = u_local*cos(directions(h_index)) + v_local*sin(directions(h_index))
  enddo
  deallocate(u_wind_in)
  deallocate(v_wind_in)
  deallocate(z_coords_input_model)
  deallocate(directions)
  deallocate(z_coords_game_wind)
  deallocate(interpolation_indices_vector)
  deallocate(interpolation_weights_vector)
  
  write(*,*) "Wind interpolation completed."
  
  ! INTERPOLATION OF THE SST
  ! ------------------------
  
  int no_of_sst_points = N_SST_POINTS
  write(*,*) "Interpolating the SST to the model grid ..."
  double *sst_out = malloc(n_scalars_h*sizeof(double))
  int min_index
  #pragma omp parallel for private(min_index)
  do (int i = 0 i < n_scalars_h ++i)
    double *distance_vector = malloc(N_SST_POINTS*sizeof(double))
    do (int j = 0 j < N_SST_POINTS ++j)
      distance_vector(j) = calculate_distance_h(latitudes_sst(j),longitudes_sst(j),latitudes_game(i),longitudes_game(i),1._wp)
    enddo
  min_index = find_min_index(distance_vector,no_of_sst_points)
  sst_out(i) = sst_in(min_index)
  deallocate(distance_vector)
  enddo
  deallocate(latitudes_sst)
  deallocate(longitudes_sst)
  deallocate(sst_in)
  deallocate(latitudes_game)
  deallocate(longitudes_game)
  
  write(*,*) "Interpolation of the SST completed."

  ! PREPARING THE OUTPUT
  ! --------------------
  
  ! clouds and precipitation are set equal to the background state
  double *densities = malloc(6*n_scalars*sizeof(double))
  !$omp parallel do private(ji)
  do ji=1,n_scalars
    ! setting the mass densities of the result
    ! condensate densities are not assimilated
    densities(i) = densities_background(i)
    densities(n_scalars+ji) = densities_background(n_scalars+ji)
    densities(2*n_scalars+ji) = densities_background(2*n_scalars+ji)
    densities(3*n_scalars+ji) = densities_background(3*n_scalars+ji)
    densities(4*n_scalars+ji) = density_moist_out(i)
    densities(5*n_scalars+ji) = spec_hum_out(ji)*density_moist_out(ji)
    if (densities(5*n_scalars+ji)<0._wp) then
      densities(5*n_scalars+ji) = 0._wp
    endif
  enddo
  !$omp end parallel do
  deallocate(density_moist_out)
  deallocate(spec_hum_out)
  deallocate(densities_background)
  
  ! writing the result to a netCDF file
  ! -----------------------------------
  
  printf("Output file: %s\n",output_file)
  write(*,*) "Writing result to output file ..."
  int densities_dimid,scalar_dimid,vector_dimid,scalar_h_dimid,single_double_dimid,
  densities_id,temperature_id,wind_id,soil_dimid
  call nc_check(nc_create(output_file,NC_CLOBBER,ncid))
  call nc_check(nc_def_dim(ncid,"densities_index",6*n_scalars,densities_dimid))
  call nc_check(nc_def_dim(ncid,"vector_index",N_VECTORS,vector_dimid))
  call nc_check(nc_def_dim(ncid,"scalar_index",n_scalars,scalar_dimid))
  call nc_check(nc_def_dim(ncid,"soil_index",N_SOIL_LAYERS*n_scalars_h,soil_dimid))
  call nc_check(nc_def_dim(ncid,"scalar_h_index",n_scalars_h,scalar_h_dimid))
  call nc_check(nc_def_dim(ncid,"single_double_dimid_index",1,single_double_dimid))
  call nc_check(nc_def_var(ncid,"densities",NC_DOUBLE,1,densities_dimid,densities_id))
  call nc_check(nc_put_att_text(ncid,densities_id,"units",strlen("kg/m^3"),"kg/m^3"))
  call nc_check(nc_def_var(ncid,"temperature",NC_DOUBLE,1,scalar_dimid,temperature_id))
  call nc_check(nc_put_att_text(ncid,temperature_id,"units",strlen("K"),"K"))
  call nc_check(nc_def_var(ncid,"wind",NC_DOUBLE,1,vector_dimid,wind_id))
  call nc_check(nc_put_att_text(ncid,wind_id,"units",strlen("m/s"),"m/s"))
  call nc_check(nc_def_var(ncid,"sst",NC_DOUBLE,1,scalar_h_dimid,sst_id))
  call nc_check(nc_put_att_text(ncid,sst_id,"units",strlen("K"),"K"))
  if (ltke_avail) then
    call nc_check(nc_def_var(ncid,"tke",NC_DOUBLE,1,scalar_dimid,tke_id))
    call nc_check(nc_put_att_text(ncid,tke_id,"units",strlen("J/kg"),"J/kg"))
  endif
  if (lt_soil_avail) then
    call nc_check(nc_def_var(ncid,"t_soil",NC_DOUBLE,1,soil_dimid,t_soil_id))
    call nc_check(nc_put_att_text(ncid,t_soil_id,"units",strlen("K"),"K"))
  endif
  call nc_check(nc_enddef(ncid))
  call nc_check(nc_put_var_double(ncid,densities_id,densities))
  call nc_check(nc_put_var_double(ncid,temperature_id,temperature_out))
  call nc_check(nc_put_var_double(ncid,wind_id,wind_out))
  call nc_check(nc_put_var_double(ncid,sst_id,sst_out))
  if (ltke_avail) then
    call nc_check(nc_put_var_double(ncid,tke_id,tke))
  endif
  if (lt_soil_avail) then
    call nc_check(nc_put_var_double(ncid,t_soil_id,t_soil))
  endif
  call nc_check(nc_close(ncid))
  write(*,*) "Result successfully written."
  
  ! freeing the stil occupied memory
  deallocate(densities)
  deallocate(temperature_out)
  deallocate(wind_out)
  deallocate(sst_out)
  deallocate(tke)
  deallocate(t_soil)
  
end program control

















