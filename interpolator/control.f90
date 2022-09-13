! This source file is part of real2GAME, which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/real2GAME

program control

  ! This file coordinates the data interpolation process.
  
  use netcdf
  use mo_shared, only: wp,n_layers_input,n_sst_points,n_points_per_layer_input,n_avg_points, &
                       nc_check,int2string,find_min_index,calculate_distance_h
  
  implicit none
  
  ! constants
  real(wp), parameter :: n_a =  6.02214076e23_wp   ! Avogadro's number
  real(wp), parameter :: p_0 = 100000._wp          ! reference pressure
  real(wp), parameter :: m_d = n_a*0.004810e-23_wp ! molar mass of dry air
  real(wp), parameter :: m_v = n_a*0.002991e-23_wp ! molar mass of water
  real(wp), parameter :: r_d = 287.057811_wp       ! specific gas constant of dry air
  real(wp), parameter :: c_d_p = 1005._wp          ! isobaric specific heat capacity of dry air
  real(wp), parameter :: scale_height = 8000._wp   ! scale_height
  
  logical               :: ltke_avail,lt_soil_avail
  integer               :: ji,jk,jl,latitudes_game_id,longitudes_game_id,z_coords_game_id,z_coords_game_wind_id, &
                           gravity_potential_game_id, &
                           directions_id,ncid,interpolation_indices_scalar_id,interpolation_weights_scalar_id,& 
                           interpolation_indices_vector_id, &
                           interpolation_weights_vector_id,layer_index,h_index,closest_index,other_index, sp_id,z_surf_id, &
                           z_coords_id,t_in_id,spec_hum_id,u_id,v_id,lat_sst_id,lon_sst_id,sst_id,densities_background_id, &
                           tke_id,t_soil_id,vector_index,min_index,densities_dimid,scalar_dimid, &
                           vector_dimid,scalar_h_dimid,single_double_dimid,densities_id,temperature_id,wind_id,soil_dimid, &
                           res_id,oro_id,n_pentagons,n_hexagons,n_scalars_h,n_scalars,n_vectors_h,n_layers,n_h_vectors, &
                           n_levels,n_v_vectors,n_vectors_per_layer,n_vectors,nsoillays
  real(wp)              :: closest_value,other_value,df,dz,gradient,delta_z,b,c,u_local,v_local,vector_to_minimize(n_layers_input)
  real(wp), allocatable :: latitudes_game(:),longitudes_game(:),z_coords_game(:),directions(:),z_coords_game_wind(:), &
                           gravity_potential_game(:),densities_background(:),tke(:),t_soil(:),z_coords_input_model(:,:), &
                           temperature_in(:,:),spec_hum_in(:,:),u_wind_in(:,:),v_wind_in(:,:),z_surf_in(:), &
                           p_surf_in(:),latitudes_sst(:),longitudes_sst(:),sst_in(:),interpolation_weights_scalar(:,:), &
                           interpolation_weights_vector(:,:),temperature_out(:),spec_hum_out(:),pressure_lowest_layer_out(:), &
                           density_moist_out(:),temperature_v(:),distance_vector(:),densities(:),exner(:),wind_out(:),sst_out(:)
  integer,  allocatable :: interpolation_indices_scalar(:,:),interpolation_indices_vector(:,:)
  character(len=4)      :: year_string,n_layers_string
  character(len=2)      :: month_string,day_string,hour_string,res_id_string,nsoillays_string,oro_id_string
  character(len=128)    :: real2game_root_dir,model_home_dir
  character(len=256)    :: background_state_file,geo_prop_file,input_file,interpol_file,output_file

  call get_command_argument(1,res_id_string)
  read(res_id_string,*) res_id
  call get_command_argument(2,n_layers_string)
  read(n_layers_string,*) n_layers
  ! grid properties
  n_pentagons = 12
  n_hexagons = 10*(2**(2*res_id)-1)
  n_scalars_h = n_pentagons+n_hexagons
  n_scalars = n_layers*n_scalars_h
  n_vectors_h = (5*n_pentagons/2 + 6/2*n_hexagons)
  n_h_vectors = n_layers*n_vectors_h
  n_levels = n_layers+1
  n_v_vectors = n_levels*n_scalars_h
  n_vectors_per_layer = n_vectors_h+n_scalars_h
  n_vectors = n_h_vectors+n_v_vectors
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
  
  write(*,*) "Background state file:",trim(background_state_file)
  
  ! Allocating memory for the grid properties.
  allocate(latitudes_game(n_scalars_h))
  allocate(longitudes_game(n_scalars_h))
  allocate(z_coords_game(n_scalars))
  allocate(directions(n_vectors_h))
  allocate(z_coords_game_wind(n_vectors))
  allocate(gravity_potential_game(n_scalars))
  ! Reading the grid properties.
  geo_prop_file = trim(model_home_dir) // "/grid_generator/grids/RES" // trim(int2string(res_id)) // "_L" // &
                  trim(int2string(n_layers)) // "_ORO" // trim(int2string(oro_id)) // ".nc"
  write(*,*) "Grid file:",trim(geo_prop_file)
  write(*,*) "Reading grid file ..."
  call nc_check(nf90_open(trim(geo_prop_file),NF90_NOWRITE,ncid))
  call nc_check(nf90_inq_varid(ncid,"latitude_scalar",latitudes_game_id))
  call nc_check(nf90_inq_varid(ncid,"longitude_scalar",longitudes_game_id))
  call nc_check(nf90_inq_varid(ncid,"z_scalar",z_coords_game_id))
  call nc_check(nf90_inq_varid(ncid,"direction",directions_id))
  call nc_check(nf90_inq_varid(ncid,"z_vector",z_coords_game_wind_id))
  call nc_check(nf90_inq_varid(ncid,"gravity_potential",gravity_potential_game_id))
  call nc_check(nf90_get_var(ncid,latitudes_game_id,latitudes_game))
  call nc_check(nf90_get_var(ncid,longitudes_game_id,longitudes_game))
  call nc_check(nf90_get_var(ncid,z_coords_game_id,z_coords_game))
  call nc_check(nf90_get_var(ncid,directions_id,directions))
  call nc_check(nf90_get_var(ncid,z_coords_game_wind_id,z_coords_game_wind))
  call nc_check(nf90_get_var(ncid,gravity_potential_game_id,gravity_potential_game))
  call nc_check(nf90_close(ncid))
  write(*,*) "Grid file read."
  
  ! constructing the filename of the input file for GAME
  output_file = trim(model_home_dir) // "/nwp_init/" // year_string // month_string // day_string // hour_string // ".nc"
  
  ! These are the arrays of the background state.
  allocate(densities_background(6*n_scalars))
  allocate(tke(n_scalars))
  allocate(t_soil(nsoillays*n_scalars_h))
  
  ! Reading the background state.
  write(*,*) "Reading background state ..."
  call nc_check(nf90_open(background_state_file,NF90_NOWRITE,ncid))
  call nc_check(nf90_inq_varid(ncid,"densities",densities_background_id))
  ltke_avail = .false.
  if (nf90_inq_varid(ncid,"tke",tke_id) == 0) then
    ltke_avail = .true.
    call nc_check(nf90_inq_varid(ncid,"tke",tke_id))
    write(*,*) "TKE found in background state file."
  else
    write(*,*) "TKE not found in background state file."
  endif
  lt_soil_avail = .false.
  if (nf90_inq_varid(ncid,"t_soil",t_soil_id) == 0) then
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
  allocate(latitudes_sst(n_sst_points))
  allocate(longitudes_sst(n_sst_points))
  allocate(sst_in(n_sst_points))
  
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
  call nc_check(nf90_inq_varid(ncid,"lat_sst",lat_sst_id))
  call nc_check(nf90_inq_varid(ncid,"lon_sst",lon_sst_id))
  call nc_check(nf90_inq_varid(ncid,"sst",sst_id))
  call nc_check(nf90_get_var(ncid,z_coords_id,z_coords_input_model))
  call nc_check(nf90_get_var(ncid,t_in_id,temperature_in))
  call nc_check(nf90_get_var(ncid,spec_hum_id,spec_hum_in))
  call nc_check(nf90_get_var(ncid,u_id,u_wind_in))
  call nc_check(nf90_get_var(ncid,v_id,v_wind_in))
  call nc_check(nf90_get_var(ncid,z_surf_id,z_surf_in))
  call nc_check(nf90_get_var(ncid,sp_id,p_surf_in))
  call nc_check(nf90_get_var(ncid,lat_sst_id,latitudes_sst))
  call nc_check(nf90_get_var(ncid,lon_sst_id,longitudes_sst))
  call nc_check(nf90_get_var(ncid,sst_id,sst_in))
  call nc_check(nf90_close(ncid))
  write(*,*) "Input read."
  
  ! memory alloction for the interpolation indices and weights
  allocate(interpolation_indices_scalar(n_scalars_h,n_avg_points))
  allocate(interpolation_weights_scalar(n_scalars_h,n_avg_points))
  allocate(interpolation_indices_vector(n_vectors_h,n_avg_points))
  allocate(interpolation_weights_vector(n_vectors_h,n_avg_points))
  
  write(*,*) "Reading the interpolation indices and weights."
  
  ! constructing the name of the interpolation indices and weights file
  interpol_file  = trim(real2game_root_dir) // "/interpolation_files/icon-global2game" // trim(int2string(res_id)) // ".nc"
  
  ! reading the interpolation file
  call nc_check(nf90_open(trim(interpol_file),NF90_NOWRITE,ncid))
  call nc_check(nf90_inq_varid(ncid,"interpolation_indices_scalar",interpolation_indices_scalar_id))
  call nc_check(nf90_inq_varid(ncid,"interpolation_weights_scalar",interpolation_weights_scalar_id))
  call nc_check(nf90_inq_varid(ncid,"interpolation_indices_vector",interpolation_indices_vector_id))
  call nc_check(nf90_inq_varid(ncid,"interpolation_weights_vector",interpolation_weights_vector_id))
  call nc_check(nf90_get_var(ncid,interpolation_indices_scalar_id,interpolation_indices_scalar))
  call nc_check(nf90_get_var(ncid,interpolation_weights_scalar_id,interpolation_weights_scalar))
  call nc_check(nf90_get_var(ncid,interpolation_indices_vector_id,interpolation_indices_vector))
  call nc_check(nf90_get_var(ncid,interpolation_weights_vector_id,interpolation_weights_vector))
  call nc_check(nf90_close(ncid))
  write(*,*) "Interpolation indices and weights read."
  
  ! Begin of the actual interpolation.
  
  ! INTERPOLATION OF SCALAR QUANTITIES
  ! ----------------------------------
  
  write(*,*) "Starting the interpolation of scalar quantities ..."
  
  ! These are the arrays for the result of the interpolation process.
  allocate(temperature_out(n_scalars))
  allocate(spec_hum_out(n_scalars))
  
  !$omp parallel workshare
  temperature_out = 0._wp
  spec_hum_out = 0._wp
  !$omp end parallel workshare
  
  !$omp parallel do private(ji,jk,jl,vector_to_minimize,layer_index,h_index,closest_value,other_value,df,dz,gradient,delta_z)
  do ji=1,n_scalars
    layer_index = (ji-1)/n_scalars_h
    h_index = ji - layer_index*n_scalars_h
    
    ! loop over all points over which the averaging is executed
    do jk=1,n_avg_points
      ! computing linear vertical interpolation
      ! vertical distance vector
      do jl=1,n_layers_input
        vector_to_minimize(jl) = abs(z_coords_game(ji) - z_coords_input_model(interpolation_indices_scalar(h_index,jk),jl))
      enddo
      
      ! closest vertical index
      closest_index = find_min_index(vector_to_minimize,n_layers_input)
      
      ! value at the closest vertical index
      closest_value = temperature_in(interpolation_indices_scalar(h_index,jk),closest_index)
      
      other_index = closest_index-1
      if (z_coords_game(ji)<z_coords_input_model(interpolation_indices_scalar(h_index,jk),closest_index)) then
        other_index = closest_index+1
      endif
      
      ! avoiding array excess
      if (other_index==n_layers_input+1) then
        other_index = closest_index-1
      endif
      
      ! the value at the second point used for vertical interpolation
      other_value = temperature_in(interpolation_indices_scalar(h_index,jk),other_index)
      
      ! computing the vertical gradient of u in the input model
      df = closest_value - other_value
      dz = z_coords_input_model(interpolation_indices_scalar(h_index,jk),closest_index) &
           - z_coords_input_model(interpolation_indices_scalar(h_index,jk),other_index)
      gradient = df/dz
      
      delta_z = z_coords_game(ji) - z_coords_input_model(interpolation_indices_scalar(h_index,jk),closest_index)
      
      ! vertical interpolation of the temperature
      temperature_out(ji) = temperature_out(ji) + interpolation_weights_scalar(h_index,jk)*(closest_value + delta_z*gradient)
      
      ! vertical interpolation of the specific humidity
      closest_value = spec_hum_in(interpolation_indices_scalar(h_index,jk),closest_index)
      other_value = spec_hum_in(interpolation_indices_scalar(h_index,jk),other_index)
      ! computing the vertical gradient of the specific humidity in the input model
      df = closest_value - other_value
      gradient = df/dz
      
      ! specific humidity
      spec_hum_out(ji) = spec_hum_out(ji) + interpolation_weights_scalar(h_index,jk)*(closest_value + delta_z*gradient)
    enddo
  enddo
  !$omp end parallel do
  ! these array are now interpolated and not needed any further
  
  deallocate(spec_hum_in)
  deallocate(temperature_in)
  ! surface pressure interpolation
  allocate(pressure_lowest_layer_out(n_scalars_h))
  !$omp parallel workshare
  pressure_lowest_layer_out = 0._wp
  !$omp end parallel workshare
  
  !$omp parallel do private(ji,jk)
  do ji=1,n_scalars_h
    do jk=1,n_avg_points
      pressure_lowest_layer_out(ji) = pressure_lowest_layer_out(ji) &
      ! horizontal component of the interpolation
      + interpolation_weights_scalar(ji,jk)*p_surf_in(interpolation_indices_scalar(ji,jk)) &
      ! vertical component of the interpolation according to the barometric height formula
      *exp(-(z_coords_game(n_scalars-n_scalars_h+ji) - z_surf_in(interpolation_indices_scalar(ji,jk)))/scale_height)
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
  allocate(density_moist_out(n_scalars))
  ! firstly setting the virtual temperature
  allocate(temperature_v(n_scalars))
  !$omp parallel do private(ji)
  do ji=1,n_scalars
    temperature_v(ji) = temperature_out(ji)*(1._wp + spec_hum_out(ji)*(m_d/m_v-1._wp))
  enddo
  !$omp end parallel do
  ! the Exner pressure is just a temporarily needed helper variable here to integrate the hydrostatic equation
  allocate(exner(n_scalars))
  do ji=n_scalars,1,-1
    layer_index = (ji-1)/n_scalars_h
    h_index = ji - layer_index*n_scalars_h
    if (layer_index==n_layers-1) then
      density_moist_out(ji) = pressure_lowest_layer_out(h_index)/(r_d*temperature_v(ji))
      exner(ji) = (density_moist_out(ji)*r_d*temperature_v(ji)/p_0)**(r_d/c_d_p)
    else
      ! solving a quadratic equation for the Exner pressure
      b = -0.5_wp*exner(ji+n_scalars_h)/temperature_v(ji+n_scalars_h) &
      *(temperature_v(ji) - temperature_v(ji+n_scalars_h) &
      + 2._wp/c_d_p*(gravity_potential_game(ji) - gravity_potential_game(ji+n_scalars_h)))
      c = exner(ji+n_scalars_h)**2*temperature_v(ji)/temperature_v(ji+n_scalars_h)
      exner(ji) = b + (b**2 + c)**0.5_wp
      density_moist_out(ji) = p_0*exner(ji)**(c_d_p/r_d)/(r_d*temperature_v(ji))
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
  allocate(wind_out(n_vectors))
  !$omp parallel workshare
  wind_out = 0._wp
  !$omp end parallel workshare
  ! loop over all horizontal vector points
  !$omp parallel do private(ji,jk,jl,vector_to_minimize,h_index,layer_index,vector_index,closest_index,other_index,closest_value, &
  !$omp other_value,df,dz,gradient,delta_z,u_local,v_local)
  do ji=1,n_h_vectors
    layer_index = (ji-1)/n_vectors_h
    h_index = ji - layer_index*n_vectors_h
    vector_index = n_scalars_h + layer_index*n_vectors_per_layer + h_index
     
    ! the u- and v-components of the wind at the grid point of GAME
    u_local = 0._wp
    v_local = 0._wp
    ! loop over all horizontal points that are used for averaging
    do jk=1,n_avg_points
      ! computing linear vertical interpolation
      ! vertical distance vector
      do jl=1,n_layers_input
        vector_to_minimize(jl) = abs(z_coords_game_wind(vector_index) &
                                    - z_coords_input_model(interpolation_indices_vector(h_index,jk),jl))
      enddo
      ! closest vertical index
      closest_index = find_min_index(vector_to_minimize,n_layers_input)
      ! value at the closest vertical index
      closest_value = u_wind_in(interpolation_indices_vector(h_index,jk),closest_index)
    
      other_index = closest_index-1
      if (z_coords_game_wind(vector_index) < z_coords_input_model(interpolation_indices_vector(h_index,jk),closest_index)) then
        other_index = closest_index+1
      endif
      ! avoiding array excess
      if (other_index==n_layers_input+1) then
        other_index = closest_index-1
      endif
    
      ! the value at the second point used for vertical interpolation
      other_value = u_wind_in(interpolation_indices_vector(h_index,jk),other_index)
    
      ! computing the vertical gradient of u in the input model
      df = closest_value - other_value
      dz = z_coords_input_model(interpolation_indices_vector(h_index,jk),closest_index) &
           - z_coords_input_model(interpolation_indices_vector(h_index,jk),other_index)
      gradient = df/dz
    
      delta_z = z_coords_game_wind(vector_index) - z_coords_input_model(interpolation_indices_vector(h_index,jk),closest_index)
    
      u_local = u_local + interpolation_weights_vector(h_index,jk)*(closest_value + gradient*delta_z)
    
      ! vertical interpolation of v
      closest_value = v_wind_in(interpolation_indices_vector(h_index,jk),closest_index)
      other_value = v_wind_in(interpolation_indices_vector(h_index,jk),other_index)
      ! computing the vertical gradient of v in the input model
      df = closest_value - other_value
      gradient = df/dz
    
      v_local = v_local + interpolation_weights_vector(h_index,jk)*(closest_value + gradient*delta_z)
    enddo
    
    ! projection onto the direction of the vector in GAME
    wind_out(vector_index) = u_local*cos(directions(h_index)) + v_local*sin(directions(h_index))
  enddo
  !$omp end parallel do
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
  
  write(*,*) "Interpolating the SST to the model grid ..."
  allocate(sst_out(n_scalars_h))
  allocate(distance_vector(n_sst_points))
  !$omp parallel do private(ji,jk,distance_vector,min_index)
  do ji=1,n_scalars_h
    do jk=1,n_sst_points
      distance_vector(jk) = calculate_distance_h(latitudes_sst(jk),longitudes_sst(jk),latitudes_game(ji),longitudes_game(ji),1._wp)
    enddo
    min_index = find_min_index(distance_vector,n_sst_points)
    sst_out(ji) = sst_in(min_index)
  enddo
  !$omp end parallel do
  deallocate(distance_vector)
  deallocate(latitudes_sst)
  deallocate(longitudes_sst)
  deallocate(sst_in)
  deallocate(latitudes_game)
  deallocate(longitudes_game)
  
  write(*,*) "Interpolation of the SST completed."

  ! PREPARING THE OUTPUT
  ! --------------------
  
  ! clouds and precipitation are set equal to the background state
  allocate(densities(6*n_scalars))
  !$omp parallel do private(ji)
  do ji=1,n_scalars
    ! setting the mass densities of the result
    ! condensate densities are not assimilated
    densities(ji) = densities_background(ji)
    densities(n_scalars+ji) = densities_background(n_scalars+ji)
    densities(2*n_scalars+ji) = densities_background(2*n_scalars+ji)
    densities(3*n_scalars+ji) = densities_background(3*n_scalars+ji)
    densities(4*n_scalars+ji) = density_moist_out(ji)
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
  
  write(*,*) "Output file:",trim(output_file)
  write(*,*) "Writing result to output file ..."
  call nc_check(nf90_create(trim(output_file),NF90_CLOBBER,ncid))
  call nc_check(nf90_def_dim(ncid,"densities_index",6*n_scalars,densities_dimid))
  call nc_check(nf90_def_dim(ncid,"vector_index",n_vectors,vector_dimid))
  call nc_check(nf90_def_dim(ncid,"scalar_index",n_scalars,scalar_dimid))
  call nc_check(nf90_def_dim(ncid,"soil_index",nsoillays*n_scalars_h,soil_dimid))
  call nc_check(nf90_def_dim(ncid,"scalar_h_index",n_scalars_h,scalar_h_dimid))
  call nc_check(nf90_def_dim(ncid,"single_double_dimid_index",1,single_double_dimid))
  call nc_check(nf90_def_var(ncid,"densities",NF90_REAL,densities_dimid,densities_id))
  call nc_check(nf90_put_att(ncid,densities_id,"units","kg/m^3"))
  call nc_check(nf90_def_var(ncid,"temperature",NF90_REAL,scalar_dimid,temperature_id))
  call nc_check(nf90_put_att(ncid,temperature_id,"units","K"))
  call nc_check(nf90_def_var(ncid,"wind",NF90_REAL,vector_dimid,wind_id))
  call nc_check(nf90_put_att(ncid,wind_id,"units","m/s"))
  call nc_check(nf90_def_var(ncid,"sst",NF90_REAL,scalar_h_dimid,sst_id))
  call nc_check(nf90_put_att(ncid,sst_id,"units","K"))
  if (ltke_avail) then
    call nc_check(nf90_def_var(ncid,"tke",NF90_REAL,scalar_dimid,tke_id))
    call nc_check(nf90_put_att(ncid,tke_id,"units","J/kg"))
  endif
  if (lt_soil_avail) then
    call nc_check(nf90_def_var(ncid,"t_soil",NF90_REAL,soil_dimid,t_soil_id))
    call nc_check(nf90_put_att(ncid,t_soil_id,"units","K"))
  endif
  call nc_check(nf90_enddef(ncid))
  call nc_check(nf90_put_var(ncid,densities_id,densities))
  call nc_check(nf90_put_var(ncid,temperature_id,temperature_out))
  call nc_check(nf90_put_var(ncid,wind_id,wind_out))
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
  deallocate(densities)
  deallocate(temperature_out)
  deallocate(wind_out)
  deallocate(sst_out)
  deallocate(tke)
  deallocate(t_soil)
  
end program control

















