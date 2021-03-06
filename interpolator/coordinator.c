/*
This source file is part of real2GAME, which is released under the MIT license.
Github repository: https://github.com/OpenNWP/real2GAME
*/

/*
This file coordinates the data interpolation process.
*/

#include <stdlib.h>
#include "../game_properties.h"
#include "../header.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <netcdf.h>
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(1);}
#define NCCHECK(e) {if(e != 0) NCERR(e)}
#define SCALE_HEIGHT 8000.0
#define P_0 100000.0
#define R_D 287.057811
#define C_D_P 1005.0
#define N_A (6.02214076e23)
#define M_D (N_A*0.004810e-23)
#define M_V (N_A*0.002991e-23)

int main(int argc, char *argv[])
{
    char year_string[strlen(argv[1]) + 1];
    strcpy(year_string, argv[1]);
    char month_string[strlen(argv[2]) + 1];
    strcpy(month_string, argv[2]);
    char day_string[strlen(argv[3]) + 1];
    strcpy(day_string, argv[3]);
    char hour_string[strlen(argv[4]) + 1];
    strcpy(hour_string, argv[4]);
    char model_home_dir[strlen(argv[5]) + 1];
    strcpy(model_home_dir, argv[5]);
	int ORO_ID;
    ORO_ID = strtod(argv[6], NULL);
    char BACKGROUND_STATE_FILE[strlen(argv[7]) + 1];
    strcpy(BACKGROUND_STATE_FILE, argv[7]);
    char real2game_root_dir[strlen(argv[8]) + 1];
    strcpy(real2game_root_dir, argv[8]);
	printf("Background state file: %s\n", BACKGROUND_STATE_FILE);
    
    // Allocating memory for the grid properties.
    double *latitudes_game = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *longitudes_game = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *z_coords_game = malloc(NO_OF_SCALARS*sizeof(double));
    double *directions = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *z_coords_game_wind = malloc(NO_OF_VECTORS*sizeof(double));
    double *gravity_potential_game = malloc(NO_OF_SCALARS*sizeof(double));
    // Reading the grid properties.
    int ncid;
    char GEO_PROP_FILE_PRE[200];
    sprintf(GEO_PROP_FILE_PRE, "%s/grid_generator/grids/RES%d_L%d_ORO%d.nc", model_home_dir, RES_ID, NO_OF_LAYERS, ORO_ID);
    char GEO_PROP_FILE[strlen(GEO_PROP_FILE_PRE) + 1];
    strcpy(GEO_PROP_FILE, GEO_PROP_FILE_PRE);
	printf("Grid file: %s\n", GEO_PROP_FILE);
	printf("Reading grid file ...\n");
    NCCHECK(nc_open(GEO_PROP_FILE, NC_NOWRITE, &ncid));
    int latitudes_game_id, longitudes_game_id, z_coords_game_id, z_coords_game_wind_id, gravity_potential_game_id, directions_id;
    NCCHECK(nc_inq_varid(ncid, "latitude_scalar", &latitudes_game_id));
    NCCHECK(nc_inq_varid(ncid, "longitude_scalar", &longitudes_game_id));
    NCCHECK(nc_inq_varid(ncid, "z_scalar", &z_coords_game_id));
    NCCHECK(nc_inq_varid(ncid, "direction", &directions_id));
    NCCHECK(nc_inq_varid(ncid, "z_vector", &z_coords_game_wind_id));
    NCCHECK(nc_inq_varid(ncid, "gravity_potential", &gravity_potential_game_id));
    NCCHECK(nc_get_var_double(ncid, latitudes_game_id, &latitudes_game[0]));
    NCCHECK(nc_get_var_double(ncid, longitudes_game_id, &longitudes_game[0]));
    NCCHECK(nc_get_var_double(ncid, z_coords_game_id, &z_coords_game[0]));
    NCCHECK(nc_get_var_double(ncid, directions_id, &directions[0]));
    NCCHECK(nc_get_var_double(ncid, z_coords_game_wind_id, &z_coords_game_wind[0]));
    NCCHECK(nc_get_var_double(ncid, gravity_potential_game_id, &gravity_potential_game[0]));
    NCCHECK(nc_close(ncid));
	printf("Grid file read.\n");
	
	// constructing the filename of the input file for GAME
    char output_file_pre[200];
    sprintf(output_file_pre, "%s/nwp_init/%s%s%s%s.nc", model_home_dir, year_string, month_string, day_string, hour_string);
    char output_file[strlen(output_file_pre) + 1];
    strcpy(output_file, output_file_pre);
    
    // These are the arrays of the background state.
    double *densities_background = malloc(6*NO_OF_SCALARS*sizeof(double));
    double *tke = malloc(NO_OF_SCALARS*sizeof(double));
    double *t_soil = malloc(NO_OF_SOIL_LAYERS*NO_OF_SCALARS_H*sizeof(double));
    
    // Reading the background state.
	printf("Reading background state ...\n");
    NCCHECK(nc_open(BACKGROUND_STATE_FILE, NC_NOWRITE, &ncid));
    int densities_background_id, tke_avail, tke_id, t_soil_avail, t_soil_id;
    NCCHECK(nc_inq_varid(ncid, "densities", &densities_background_id));
    tke_avail = 0;
    if (nc_inq_varid(ncid, "tke", &tke_id) == 0)
    {
    	tke_avail = 1;
		NCCHECK(nc_inq_varid(ncid, "tke", &tke_id));
		printf("TKE found in background state file.\n");
    }
    else
    {
    	printf("TKE not found in background state file.\n");
    }
    t_soil_avail = 0;
    if (nc_inq_varid(ncid, "t_soil", &t_soil_id) == 0)
    {
    	t_soil_avail = 1;
		NCCHECK(nc_inq_varid(ncid, "t_soil", &t_soil_id));
		printf("Soil temperature found in background state file.\n");
    }
    else
    {
    	printf("Soil temperature not found in background state file.\n");
    }
    NCCHECK(nc_get_var_double(ncid, densities_background_id, &densities_background[0]));
    if (tke_avail == 1)
    {
		NCCHECK(nc_get_var_double(ncid, tke_id, &tke[0]));
    }
    if (t_soil_avail == 1)
    {
		NCCHECK(nc_get_var_double(ncid, t_soil_id, &t_soil[0]));
    }
    NCCHECK(nc_close(ncid));
	printf("Background state read.\n");	
	
	// allocating the memory for the analysis of the other model
	double (*z_coords_input_model)[NO_OF_LEVELS_INPUT] = malloc(sizeof(double[NO_OF_POINTS_PER_LAYER_INPUT][NO_OF_LEVELS_INPUT]));
	double (*temperature_in)[NO_OF_LEVELS_INPUT] = malloc(sizeof(double[NO_OF_POINTS_PER_LAYER_INPUT][NO_OF_LEVELS_INPUT]));
	double (*spec_hum_in)[NO_OF_LEVELS_INPUT] = malloc(sizeof(double[NO_OF_POINTS_PER_LAYER_INPUT][NO_OF_LEVELS_INPUT]));
	double (*u_wind_in)[NO_OF_LEVELS_INPUT] = malloc(sizeof(double[NO_OF_POINTS_PER_LAYER_INPUT][NO_OF_LEVELS_INPUT]));
	double (*v_wind_in)[NO_OF_LEVELS_INPUT] = malloc(sizeof(double[NO_OF_POINTS_PER_LAYER_INPUT][NO_OF_LEVELS_INPUT]));
	double *z_surf_in = malloc(NO_OF_POINTS_PER_LAYER_INPUT*sizeof(double));
	double *p_surf_in = malloc(NO_OF_POINTS_PER_LAYER_INPUT*sizeof(double));
	double *latitudes_sst = malloc(NO_OF_SST_POINTS*sizeof(double));
	double *longitudes_sst = malloc(NO_OF_SST_POINTS*sizeof(double));
	double *sst_in = malloc(NO_OF_SST_POINTS*sizeof(double));
    
    // determining the name of the input file
    char input_file_pre[200];
    sprintf(input_file_pre, "%s/input/obs_%s%s%s%s.nc", real2game_root_dir, year_string, month_string, day_string, hour_string);
    char input_file[strlen(input_file_pre) + 1];
    strcpy(input_file, input_file_pre);
	printf("Input file: %s\n", input_file);
    
    // reading the analysis of the other model
	printf("Reading input ...\n");
    int sp_id, z_surf_id, z_coords_id, t_in_id, spec_hum_id, u_id, v_id, lat_sst_id, lon_sst_id, sst_id;
    NCCHECK(nc_open(input_file, NC_NOWRITE, &ncid));
    // Defining the variables.
    NCCHECK(nc_inq_varid(ncid, "z_height", &z_coords_id));
    NCCHECK(nc_inq_varid(ncid, "temperature", &t_in_id));
    NCCHECK(nc_inq_varid(ncid, "spec_humidity", &spec_hum_id));
    NCCHECK(nc_inq_varid(ncid, "u_wind", &u_id));
    NCCHECK(nc_inq_varid(ncid, "v_wind", &v_id));
    NCCHECK(nc_inq_varid(ncid, "z_surface", &z_surf_id));
    NCCHECK(nc_inq_varid(ncid, "pressure_surface", &sp_id));
    NCCHECK(nc_inq_varid(ncid, "lat_sst", &lat_sst_id));
    NCCHECK(nc_inq_varid(ncid, "lon_sst", &lon_sst_id));
    NCCHECK(nc_inq_varid(ncid, "sst", &sst_id));
    NCCHECK(nc_get_var_double(ncid, z_coords_id, &z_coords_input_model[0][0]));
    NCCHECK(nc_get_var_double(ncid, t_in_id, &temperature_in[0][0]));
    NCCHECK(nc_get_var_double(ncid, spec_hum_id, &spec_hum_in[0][0]));
    NCCHECK(nc_get_var_double(ncid, u_id, &u_wind_in[0][0]));
    NCCHECK(nc_get_var_double(ncid, v_id, &v_wind_in[0][0]));
    NCCHECK(nc_get_var_double(ncid, z_surf_id, &z_surf_in[0]));
    NCCHECK(nc_get_var_double(ncid, sp_id, &p_surf_in[0]));
    NCCHECK(nc_get_var_double(ncid, lat_sst_id, &latitudes_sst[0]));
    NCCHECK(nc_get_var_double(ncid, lon_sst_id, &longitudes_sst[0]));
    NCCHECK(nc_get_var_double(ncid, sst_id, &sst_in[0]));
    NCCHECK(nc_close(ncid));
	printf("Input read.\n");
	
	// memory alloction for the interpolation indices and weights
	int (*interpolation_indices_scalar)[NO_OF_AVG_POINTS] = malloc(sizeof(int[NO_OF_SCALARS_H][NO_OF_AVG_POINTS]));
	double (*interpolation_weights_scalar)[NO_OF_AVG_POINTS] = malloc(sizeof(double[NO_OF_SCALARS_H][NO_OF_AVG_POINTS]));
	int (*interpolation_indices_vector)[NO_OF_AVG_POINTS] = malloc(sizeof(int[NO_OF_VECTORS_H][NO_OF_AVG_POINTS]));
	double (*interpolation_weights_vector)[NO_OF_AVG_POINTS] = malloc(sizeof(double[NO_OF_VECTORS_H][NO_OF_AVG_POINTS]));
	
	printf("Reading the interpolation indices and weights.\n");
	
	// constructing the name of the interpolation indices and weights file
	char interpol_file_pre[200];
	sprintf(interpol_file_pre, "%s/interpolation_files/icon-global2game%d.nc", real2game_root_dir, RES_ID);
    char interpol_file[strlen(interpol_file_pre) + 1];
    strcpy(interpol_file, interpol_file_pre);
	
	// reading the interpolation file
	int interpolation_indices_scalar_id, interpolation_weights_scalar_id, interpolation_indices_vector_id, interpolation_weights_vector_id;
    NCCHECK(nc_open(interpol_file, NC_NOWRITE, &ncid));
    NCCHECK(nc_inq_varid(ncid, "interpolation_indices_scalar", &interpolation_indices_scalar_id));
    NCCHECK(nc_inq_varid(ncid, "interpolation_weights_scalar", &interpolation_weights_scalar_id));
	NCCHECK(nc_inq_varid(ncid, "interpolation_indices_vector", &interpolation_indices_vector_id));
	NCCHECK(nc_inq_varid(ncid, "interpolation_weights_vector", &interpolation_weights_vector_id));
    NCCHECK(nc_get_var_int(ncid, interpolation_indices_scalar_id, &interpolation_indices_scalar[0][0]));
    NCCHECK(nc_get_var_double(ncid, interpolation_weights_scalar_id, &interpolation_weights_scalar[0][0]));
	NCCHECK(nc_get_var_int(ncid, interpolation_indices_vector_id, &interpolation_indices_vector[0][0]));
	NCCHECK(nc_get_var_double(ncid, interpolation_weights_vector_id, &interpolation_weights_vector[0][0]));
    NCCHECK(nc_close(ncid));
    printf("Interpolation indices and weights read.\n");
	
	// Begin of the actual interpolation.
    
    /*
    INTERPOLATION OF SCALAR QUANTITIES
    ----------------------------------
    */
	
	printf("Starting the interpolation of scalar quantities ...\n");
    
    // These are the arrays for the result of the interpolation process.
    double *temperature_out = calloc(1, NO_OF_SCALARS*sizeof(double));
	double *spec_hum_out = calloc(1, NO_OF_SCALARS*sizeof(double));
	
	int no_of_levels_input = NO_OF_LEVELS_INPUT;
	
    int layer_index, h_index, closest_index, other_index;
    double closest_value, other_value, df, dz, gradient, delta_z;
    #pragma omp parallel for private(layer_index, h_index, closest_value, other_value, df, dz, gradient, delta_z)
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	layer_index = i/NO_OF_SCALARS_H;
    	h_index = i - layer_index*NO_OF_SCALARS_H;
    	
    	// loop over all points over which the averaging is executed
    	for (int j = 0; j < NO_OF_AVG_POINTS; ++j)
    	{
    		// computing linear vertical interpolation
    		// vertical distance vector
    		double vector_to_minimize[NO_OF_LEVELS_INPUT];
    		for (int k = 0; k < NO_OF_LEVELS_INPUT; ++k)
    		{
    			vector_to_minimize[k] = fabs(z_coords_game[i]
    			- z_coords_input_model[interpolation_indices_scalar[h_index][j]][k]);
    		}
    		
    		// closest vertical index
    		closest_index = find_min_index(vector_to_minimize, &no_of_levels_input);
    		
    		// value at the closest vertical index
    		closest_value = temperature_in[interpolation_indices_scalar[h_index][j]][closest_index];
    		
    		other_index = closest_index - 1;
    		if (z_coords_game[i] < z_coords_input_model[interpolation_indices_scalar[h_index][j]][closest_index])
    		{
    			other_index = closest_index + 1;
    		}
    		
    		// avoiding array excess
    		if (other_index == NO_OF_LEVELS_INPUT)
    		{
    			other_index = closest_index - 1;
    		}
    		
    		// the value at the second point used for vertical interpolation
    		other_value = temperature_in[interpolation_indices_scalar[h_index][j]][other_index];
    		
    		// computing the vertical gradient of u in the input model
    		df = closest_value - other_value;
    		dz = z_coords_input_model[interpolation_indices_scalar[h_index][j]][closest_index] - z_coords_input_model[interpolation_indices_scalar[h_index][j]][other_index];
    		gradient = df/dz;
    		
    		delta_z = z_coords_game[i] - z_coords_input_model[interpolation_indices_scalar[h_index][j]][closest_index];
    		
    		// vertical interpolation of the temperature
    		temperature_out[i] += interpolation_weights_scalar[h_index][j]*(closest_value + delta_z*gradient);
    		
    		// vertical interpolation of the specific humidity
    		closest_value = spec_hum_in[interpolation_indices_scalar[h_index][j]][closest_index];
    		other_value = spec_hum_in[interpolation_indices_scalar[h_index][j]][other_index];
    		// computing the vertical gradient of the specific humidity in the input model
    		df = closest_value - other_value;
    		gradient = df/dz;
    		
    		// specific humidity
    		spec_hum_out[i] += interpolation_weights_scalar[h_index][j]*(closest_value + delta_z*gradient);
    	}
    }
    // these array are now interpolated and not needed any further
    free(spec_hum_in);
    free(temperature_in);
    
    // surface pressure interpolation
    double *pressure_lowest_layer_out = calloc(1, NO_OF_SCALARS_H*sizeof(double));
    # pragma omp parallel for
    for (int i = 0; i < NO_OF_SCALARS_H; ++i)
    {
    	for (int j = 0; j < NO_OF_AVG_POINTS; ++j)
    	{
    		pressure_lowest_layer_out[i]
    		// horizontal component of the interpolation
    		+= interpolation_weights_scalar[i][j]*p_surf_in[interpolation_indices_scalar[i][j]]
    		// vertical component of the interpolation according to the barometric height formula
    		*exp(-(z_coords_game[NO_OF_SCALARS - NO_OF_SCALARS_H + i] - z_surf_in[interpolation_indices_scalar[i][j]])/SCALE_HEIGHT);
    	}
    }
    free(z_coords_game);
	free(z_surf_in);
	free(p_surf_in);
	// no more interpolation of scalar quantities will be executed, this is why these interpolation indices and weights can be freed
	free(interpolation_indices_scalar);
	free(interpolation_weights_scalar);
    
    // density is determined out of the hydrostatic equation
    double *density_moist_out = malloc(NO_OF_SCALARS*sizeof(double));
    // firstly setting the virtual temperature
    double *temperature_v = malloc(NO_OF_SCALARS*sizeof(double));
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	temperature_v[i] = temperature_out[i]*(1.0 + spec_hum_out[i]*(M_D/M_V - 1.0));
    }
    // the Exner pressure is just a temporarily needed helper variable here to integrate the hydrostatic equation
    double *exner = malloc(NO_OF_SCALARS*sizeof(double));
    double b, c;
    for (int i = NO_OF_SCALARS - 1; i >= 0; --i)
    {
    	layer_index = i/NO_OF_SCALARS_H;
    	h_index = i - layer_index*NO_OF_SCALARS_H;
    	if (layer_index == NO_OF_LAYERS - 1)
    	{
        	density_moist_out[i] = pressure_lowest_layer_out[h_index]/(R_D*temperature_v[i]);
        	exner[i] = pow((density_moist_out[i]*R_D*temperature_v[i])/P_0, R_D/C_D_P);
        }
        else
        {
			// solving a quadratic equation for the Exner pressure
			b = -0.5*exner[i + NO_OF_SCALARS_H]/temperature_v[i + NO_OF_SCALARS_H]
			*(temperature_v[i] - temperature_v[i + NO_OF_SCALARS_H]
			+ 2.0/C_D_P*(gravity_potential_game[i] - gravity_potential_game[i + NO_OF_SCALARS_H]));
			c = pow(exner[i + NO_OF_SCALARS_H], 2.0)*temperature_v[i]/temperature_v[i + NO_OF_SCALARS_H];
			exner[i] = b + pow((pow(b, 2.0) + c), 0.5);
        	density_moist_out[i] = P_0*pow(exner[i], C_D_P/R_D)/(R_D*temperature_v[i]);
        }
    }
    free(temperature_v);
    free(pressure_lowest_layer_out);
	free(gravity_potential_game);
	// the Exner pressure was only needed to integrate the hydrostatic equation
    free(exner);
    // end of the interpolation of the dry thermodynamic state
	printf("Interpolation of scalar quantities completed.\n");
	
	/*
	WIND INTERPOLATION
	------------------
	*/
	
	printf("Starting the wind interpolation ...\n");
    double *wind_out = calloc(1, NO_OF_VECTORS*sizeof(double));
    int vector_index;
    double u_local, v_local;
    // loop over all horizontal vector points
    # pragma omp parallel for private(h_index, layer_index, vector_index, closest_index, other_index, closest_value, other_value, df, dz, gradient, delta_z, u_local, v_local)
    for (int i = 0; i < NO_OF_H_VECTORS; ++i)
    {
    	layer_index = i/NO_OF_VECTORS_H;
    	h_index = i - layer_index*NO_OF_VECTORS_H;
   		vector_index = NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index;
   		
   		// the u- and v-components of the wind at the grid point of GAME
   		u_local = 0.0;
   		v_local = 0.0;
   		// loop over all horizontal points that are used for averaging
    	for (int j = 0; j < NO_OF_AVG_POINTS; ++j)
    	{
    		// computing linear vertical interpolation
    		// vertical distance vector
    		double vector_to_minimize[NO_OF_LEVELS_INPUT];
    		for (int k = 0; k < NO_OF_LEVELS_INPUT; ++k)
    		{
    			vector_to_minimize[k] = fabs(z_coords_game_wind[vector_index]
    			- z_coords_input_model[interpolation_indices_vector[h_index][j]][k]);
    		}
    		// closest vertical index
    		closest_index = find_min_index(vector_to_minimize, &no_of_levels_input);
    		// value at the closest vertical index
    		closest_value = u_wind_in[interpolation_indices_vector[h_index][j]][closest_index];
    		
    		other_index = closest_index - 1;
    		if (z_coords_game_wind[vector_index] < z_coords_input_model[interpolation_indices_vector[h_index][j]][closest_index])
    		{
    			other_index = closest_index + 1;
    		}
    		// avoiding array excess
    		if (other_index == NO_OF_LEVELS_INPUT)
    		{
    			other_index = closest_index - 1;
    		}
    		
    		// the value at the second point used for vertical interpolation
    		other_value = u_wind_in[interpolation_indices_vector[h_index][j]][other_index];
    		
    		// computing the vertical gradient of u in the input model
    		df = closest_value - other_value;
    		dz = z_coords_input_model[interpolation_indices_vector[h_index][j]][closest_index] - z_coords_input_model[interpolation_indices_vector[h_index][j]][other_index];
    		gradient = df/dz;
    		
    		delta_z = z_coords_game_wind[vector_index] - z_coords_input_model[interpolation_indices_vector[h_index][j]][closest_index];
    		
    		u_local += interpolation_weights_vector[h_index][j]*(closest_value + gradient*delta_z);
    		
    		// vertical interpolation of v
    		closest_value = v_wind_in[interpolation_indices_vector[h_index][j]][closest_index];
    		other_value = v_wind_in[interpolation_indices_vector[h_index][j]][other_index];
    		// computing the vertical gradient of v in the input model
    		df = closest_value - other_value;
    		gradient = df/dz;
    		
    		v_local += interpolation_weights_vector[h_index][j]*(closest_value + gradient*delta_z);
    	}
    	
    	// projection onto the direction of the vector in GAME
		wind_out[vector_index] = u_local*cos(directions[h_index]) + v_local*sin(directions[h_index]);
    }
    free(u_wind_in);
    free(v_wind_in);
    free(z_coords_input_model);
	free(directions);
	free(z_coords_game_wind);
	free(interpolation_indices_vector);
	free(interpolation_weights_vector);
	
	printf("Wind interpolation completed.\n");
	
	/*
	INTERPOLATION OF THE SST
	------------------------
	*/
	int no_of_sst_points = NO_OF_SST_POINTS;
	printf("Interpolating the SST to the model grid ...\n");
	double *sst_out = malloc(NO_OF_SCALARS_H*sizeof(double));
	int min_index;
	#pragma omp parallel for private(min_index)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
    	double *distance_vector = malloc(NO_OF_SST_POINTS*sizeof(double));
    	for (int j = 0; j < NO_OF_SST_POINTS; ++j)
    	{
    		double one = 1.0;
    		distance_vector[j] = calculate_distance_h(&latitudes_sst[j], &longitudes_sst[j], &latitudes_game[i], &longitudes_game[i], &one);
    	}
		min_index = find_min_index(distance_vector, &no_of_sst_points);
		sst_out[i] = sst_in[min_index];
		free(distance_vector);
	}
	free(latitudes_sst);
	free(longitudes_sst);
	free(sst_in);
	free(latitudes_game);
	free(longitudes_game);
	
	printf("Interpolation of the SST completed.\n");

	/*
	PREPARING THE OUTPUT
	--------------------
	*/
	
	// clouds and precipitation are set equal to the background state
	double *densities = malloc(6*NO_OF_SCALARS*sizeof(double));
    #pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		// setting the mass densities of the result
		// condensate densities are not assimilated
		densities[i] = densities_background[i];
		densities[NO_OF_SCALARS + i] = densities_background[NO_OF_SCALARS + i];
		densities[2*NO_OF_SCALARS + i] = densities_background[2*NO_OF_SCALARS + i];
		densities[3*NO_OF_SCALARS + i] = densities_background[3*NO_OF_SCALARS + i];
		densities[4*NO_OF_SCALARS + i] = density_moist_out[i];
		densities[5*NO_OF_SCALARS + i] = spec_hum_out[i]*density_moist_out[i];
		if (densities[5*NO_OF_SCALARS + i] < 0.0)
		{
			densities[5*NO_OF_SCALARS + i] = 0.0;
		}
    }
    free(density_moist_out);
	free(spec_hum_out);
	free(densities_background);
    
    /*
    writing the result to a NetCDF file
    -----------------------------------
    */
    
    printf("Output file: %s\n", output_file);
    printf("Writing result to output file ...\n");
    int densities_dimid, scalar_dimid, vector_dimid, scalar_h_dimid, single_double_dimid,
    densities_id, temperature_id, wind_id, soil_dimid;
    NCCHECK(nc_create(output_file, NC_CLOBBER, &ncid));
    NCCHECK(nc_def_dim(ncid, "densities_index", 6*NO_OF_SCALARS, &densities_dimid));
    NCCHECK(nc_def_dim(ncid, "vector_index", NO_OF_VECTORS, &vector_dimid));
    NCCHECK(nc_def_dim(ncid, "scalar_index", NO_OF_SCALARS, &scalar_dimid));
    NCCHECK(nc_def_dim(ncid, "soil_index", NO_OF_SOIL_LAYERS*NO_OF_SCALARS_H, &soil_dimid));
    NCCHECK(nc_def_dim(ncid, "scalar_h_index", NO_OF_SCALARS_H, &scalar_h_dimid));
    NCCHECK(nc_def_dim(ncid, "single_double_dimid_index", 1, &single_double_dimid));
    NCCHECK(nc_def_var(ncid, "densities", NC_DOUBLE, 1, &densities_dimid, &densities_id));
    NCCHECK(nc_put_att_text(ncid, densities_id, "units", strlen("kg/m^3"), "kg/m^3"));
    NCCHECK(nc_def_var(ncid, "temperature", NC_DOUBLE, 1, &scalar_dimid, &temperature_id));
    NCCHECK(nc_put_att_text(ncid, temperature_id, "units", strlen("K"), "K"));
    NCCHECK(nc_def_var(ncid, "wind", NC_DOUBLE, 1, &vector_dimid, &wind_id));
    NCCHECK(nc_put_att_text(ncid, wind_id, "units", strlen("m/s"), "m/s"));
    NCCHECK(nc_def_var(ncid, "sst", NC_DOUBLE, 1, &scalar_h_dimid, &sst_id));
    NCCHECK(nc_put_att_text(ncid, sst_id, "units", strlen("K"), "K"));
    if (tke_avail == 1)
    {
		NCCHECK(nc_def_var(ncid, "tke", NC_DOUBLE, 1, &scalar_dimid, &tke_id));
		NCCHECK(nc_put_att_text(ncid, tke_id, "units", strlen("J/kg"), "J/kg"));
    }
    if (t_soil_avail == 1)
    {
		NCCHECK(nc_def_var(ncid, "t_soil", NC_DOUBLE, 1, &soil_dimid, &t_soil_id));
		NCCHECK(nc_put_att_text(ncid, t_soil_id, "units", strlen("K"), "K"));
    }
    NCCHECK(nc_enddef(ncid));
    NCCHECK(nc_put_var_double(ncid, densities_id, &densities[0]));
    NCCHECK(nc_put_var_double(ncid, temperature_id, &temperature_out[0]));
    NCCHECK(nc_put_var_double(ncid, wind_id, &wind_out[0]));
    NCCHECK(nc_put_var_double(ncid, sst_id, &sst_out[0]));
    if (tke_avail == 1)
    {
		NCCHECK(nc_put_var_double(ncid, tke_id, &tke[0]));
    }
    if (t_soil_avail == 1)
    {
		NCCHECK(nc_put_var_double(ncid, t_soil_id, &t_soil[0]));
    }
    NCCHECK(nc_close(ncid));
    printf("Result successfully written.\n");
    
    // freeing the stil occupied memory
	free(densities);
	free(temperature_out);
	free(wind_out);
	free(sst_out);
	free(tke);
	free(t_soil);
	
	// that's it
    return 0;
}

















