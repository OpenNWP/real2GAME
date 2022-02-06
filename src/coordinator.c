/*
This source file is part of GAME-DA, which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME-DA
*/

/*
This file coordinates the data assimilation process.
*/

#include <stdlib.h>
#include "enum.h"
#include "game-da.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <netcdf.h>
#include <geos95.h>
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}
#define EPSILON 1e-4
#define SCALE_HEIGHT 8000.0
#define P_0 100000
#define R_D 287.057811

int obs_op_setup(double [], double [][NO_OF_REL_MODEL_DOFS_PER_OBS], int [][NO_OF_REL_MODEL_DOFS_PER_OBS], double [], double [], double [], double [], double [], double [], double []);
int obs_op_setup_wind(double [], double [][NO_OF_REL_MODEL_DOFS_PER_OBS], int [][NO_OF_REL_MODEL_DOFS_PER_OBS], double [], double [], double [], double [],
double [], double [], double [], double []);

int main(int argc, char *argv[])
{
	if (fmod(NO_OF_REL_MODEL_DOFS_PER_OBS, 2) == 1)
	{
		printf("NO_OF_REL_MODEL_DOFS_PER_OBS must be even.\n");
		printf("Aborting.\n");
		exit(1);
	}
	double C_D_P = 1005.0;
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
    char game_da_root_dir[strlen(argv[8]) + 1];
    strcpy(game_da_root_dir, argv[8]);
	printf("Background state file: %s\n", BACKGROUND_STATE_FILE);
    
    // Allocating memory for the grid properties.
    double *latitudes_model = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *longitudes_model = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *z_coords_model = malloc(NO_OF_SCALARS*sizeof(double));
    double *latitudes_model_wind = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *longitudes_model_wind = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *directions = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *z_coords_model_wind = malloc(NO_OF_VECTORS*sizeof(double));
    double *gravity_potential_model = malloc(NO_OF_SCALARS*sizeof(double));
    double *normal_distance = malloc(NO_OF_VECTORS*sizeof(double));
    int *from_index = malloc(NO_OF_VECTORS_H*sizeof(int));
    int *to_index = malloc(NO_OF_VECTORS_H*sizeof(int));
    int *adjacent_vector_indices_h = malloc(6*NO_OF_SCALARS_H*sizeof(int));
    // Reading the grid properties.
    int ncid_grid, retval;
    char GEO_PROP_FILE_PRE[200];
    sprintf(GEO_PROP_FILE_PRE, "%s/grid_generator/grids/RES%d_L%d_ORO%d.nc", model_home_dir, RES_ID, NO_OF_LAYERS, ORO_ID);
    char GEO_PROP_FILE[strlen(GEO_PROP_FILE_PRE) + 1];
    strcpy(GEO_PROP_FILE, GEO_PROP_FILE_PRE);
	printf("Grid file: %s\n", GEO_PROP_FILE);
	printf("Reading grid file ...\n");
    if ((retval = nc_open(GEO_PROP_FILE, NC_NOWRITE, &ncid_grid)))
        NCERR(retval);
    int latitudes_model_id, latitudes_model_wind_id, longitudes_model_id, longitudes_model_wind_id, z_coords_model_id,
    z_coords_model_wind_id, gravity_potential_model_id, normal_distance_id, from_index_id, to_index_id, adjacent_vector_indices_h_id, directions_id;
    if ((retval = nc_inq_varid(ncid_grid, "latitude_scalar", &latitudes_model_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "longitude_scalar", &longitudes_model_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "z_scalar", &z_coords_model_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "latitude_vector", &latitudes_model_wind_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "longitude_vector", &longitudes_model_wind_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "direction", &directions_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "z_vector", &z_coords_model_wind_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "gravity_potential", &gravity_potential_model_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "normal_distance", &normal_distance_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "from_index", &from_index_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "to_index", &to_index_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "adjacent_vector_indices_h", &adjacent_vector_indices_h_id)))
        NCERR(retval);
    if ((retval = nc_get_var_int(ncid_grid, from_index_id, &from_index[0])))
        NCERR(retval);
    if ((retval = nc_get_var_int(ncid_grid, to_index_id, &to_index[0])))
        NCERR(retval);
    if ((retval = nc_get_var_int(ncid_grid, adjacent_vector_indices_h_id, &adjacent_vector_indices_h[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, latitudes_model_id, &latitudes_model[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, longitudes_model_id, &longitudes_model[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, z_coords_model_id, &z_coords_model[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, latitudes_model_wind_id, &latitudes_model_wind[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, longitudes_model_wind_id, &longitudes_model_wind[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, directions_id, &directions[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, z_coords_model_wind_id, &z_coords_model_wind[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, gravity_potential_model_id, &gravity_potential_model[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, normal_distance_id, &normal_distance[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid_grid)))
        NCERR(retval);
	printf("Grid file read.\n");
	
    char output_file_pre[200];
    sprintf(output_file_pre, "%s/nwp_init/%s%s%s%s.nc", model_home_dir, year_string, month_string, day_string, hour_string);
    char output_file[strlen(output_file_pre) + 1];
    strcpy(output_file, output_file_pre);
    
    // These are the arrays of the background state.
    double *densities_background = malloc(6*NO_OF_SCALARS*sizeof(double));
    double *temperatures_background = malloc(5*NO_OF_SCALARS*sizeof(double));
    double *wind_background = malloc(NO_OF_VECTORS*sizeof(double));
    
    // Reading the background state.
	printf("Reading background state ...\n");
    int ncid;
    if ((retval = nc_open(BACKGROUND_STATE_FILE, NC_NOWRITE, &ncid)))
        NCERR(retval);
    int densities_background_id, temperatures_background_id, wind_background_id;
    if ((retval = nc_inq_varid(ncid, "densities", &densities_background_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "temperatures", &temperatures_background_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "wind", &wind_background_id)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, densities_background_id, &densities_background[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, temperatures_background_id, &temperatures_background[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, wind_background_id, &wind_background[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid)))
        NCERR(retval);
	printf("Background state read.\n");

	// saving the relevant part of the background state in one array
	double *background_dry = malloc(NO_OF_MODEL_DOFS_DRY*sizeof(double));
    #pragma omp parallel for
	for (int i = 0; i < NO_OF_MODEL_DOFS_DRY; ++i)
	{
		// temperature of the gas phase of the background state
		if (i < NO_OF_SCALARS)
		{
			background_dry[i] = temperatures_background[4*NO_OF_SCALARS + i];
		}
		// dry air density in the lowest layer of the background state
		else
		{
			background_dry[i] = densities_background[4*NO_OF_SCALARS + i - NO_OF_SCALARS_H];
		}
	}
    free(temperatures_background);	
	
	// Allocating the memory for the observations.
	double *latitude_vector_obs = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(double));
	double *longitude_vector_obs = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(double));
	double *z_coords_obs = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(double));
	double *observations_vector = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(double));
    
    char OBSERVATIONS_FILE_PRE[200];
    sprintf(OBSERVATIONS_FILE_PRE, "%s/input/obs_%s%s%s%s.nc", game_da_root_dir, year_string, month_string, day_string, hour_string);
    char OBSERVATIONS_FILE[strlen(OBSERVATIONS_FILE_PRE) + 1];
    strcpy(OBSERVATIONS_FILE, OBSERVATIONS_FILE_PRE);
	printf("Observations file: %s\n", OBSERVATIONS_FILE);
    
    // Reading the observations.
	printf("Reading observations ...\n");
    if ((retval = nc_open(OBSERVATIONS_FILE, NC_NOWRITE, &ncid)))
        NCERR(retval);
    int latitude_obs_id, longitude_obs_id, z_coords_id, obervations_id;
    // Defining the variables.
    if ((retval = nc_inq_varid(ncid, "latitude_vector", &latitude_obs_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "longitude_vector", &longitude_obs_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "z_coords_obs", &z_coords_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "observations_vector", &obervations_id)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, latitude_obs_id, &latitude_vector_obs[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, longitude_obs_id, &longitude_vector_obs[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, z_coords_id, &z_coords_obs[0])))
        NCERR(retval);  
    if ((retval = nc_get_var_double(ncid, obervations_id, &observations_vector[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid)))
    	NCERR(retval);
	printf("Observations read.\n");
	
	// Begin of the actual assimilation.
    
    /*
    DRY THERMODYNAMIC STATE ASSIMILATION
    ------------------------------------
    */
	
	printf("Starting the dry assimilation ...\n");
    
	// setting up the observations operator
	double *interpolated_model_dry = malloc(NO_OF_CHOSEN_OBSERVATIONS_DRY*sizeof(double));
	double (*obs_op_jacobian_reduced_matrix_dry)[NO_OF_REL_MODEL_DOFS_PER_OBS] = malloc(sizeof(double[NO_OF_CHOSEN_OBSERVATIONS_DRY][NO_OF_REL_MODEL_DOFS_PER_OBS]));
	int (*relevant_model_dofs_matrix_dry)[NO_OF_REL_MODEL_DOFS_PER_OBS] = malloc(sizeof(int[NO_OF_CHOSEN_OBSERVATIONS_DRY][NO_OF_REL_MODEL_DOFS_PER_OBS]));
	
	// setting up the obervations operator
	obs_op_setup(interpolated_model_dry, obs_op_jacobian_reduced_matrix_dry, relevant_model_dofs_matrix_dry,
	latitude_vector_obs, longitude_vector_obs, z_coords_obs, latitudes_model, longitudes_model, z_coords_model, background_dry);
    
    free(z_coords_model);
	// now, all the constituents of the gain matrix are known
	
	double *observations_vector_dry = malloc(NO_OF_CHOSEN_OBSERVATIONS_DRY*sizeof(double));
	// setting up the dry observations vector
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS_DRY; ++i)
	{
    	if (i < NO_OF_CHOSEN_OBSERVATIONS_DRY - NO_OF_CHOSEN_POINTS_PER_LAYER_OBS)
    	{
			observations_vector_dry[i] = observations_vector[i];
		}
		else
		{
			observations_vector_dry[i] = observations_vector[NO_OF_CHOSEN_OBSERVATIONS_MOIST + i];
		}
	}
	
    // setting up the measurement error covariance matrix
    double *obs_error_cov_dry = malloc(sizeof(double[NO_OF_CHOSEN_OBSERVATIONS_DRY]));
    double temperature_error_obs = 0.25;
    double pressure_error_obs = 100;
	#pragma omp parallel for
    for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS_DRY; ++i)
    {
    	if (i < NO_OF_CHOSEN_OBSERVATIONS_DRY - NO_OF_CHOSEN_POINTS_PER_LAYER_OBS)
    	{
			obs_error_cov_dry[i] = pow(temperature_error_obs, 2);
		}
		else
		{
			obs_error_cov_dry[i] = pow(pressure_error_obs, 2);
		}
    }
    
    // setting up the background error covariance matrix (only the diagonal)
    double (*bg_error_cov_dry)[7] = calloc(1, sizeof(double[NO_OF_MODEL_DOFS_DRY][7]));
    double temperature_error_model = 1;
    double pressure_error_model = 400;
    double e_folding_length, distance;
    e_folding_length = 750e3;
    int no_of_edges, layer_index, h_index;;
    #pragma omp parallel for private(distance, layer_index, h_index, no_of_edges)
    for (int i = 0; i < NO_OF_MODEL_DOFS_DRY; ++i)
    {
    	layer_index = i/NO_OF_SCALARS_H;
    	h_index = i - layer_index*NO_OF_SCALARS_H;
    	if (layer_index == NO_OF_LAYERS)
    	{
    		layer_index = NO_OF_LAYERS - 1;
    	}
    	no_of_edges = 6;
    	if (h_index < NO_OF_PENTAGONS)
    	{
    		no_of_edges = 5;
    	}
    	// diagonal terms
    	if (i < NO_OF_SCALARS)
    	{
			bg_error_cov_dry[i][0] = pow(temperature_error_model, 2);
    	}
    	else
    	{
    		// density = p/(R_D*T) (Gauss'ian error propagation)
    		bg_error_cov_dry[i][0] = pow(pressure_error_model/(R_D*background_dry[i - NO_OF_SCALARS_H]), 2) + pow(background_dry[i]/background_dry[i - NO_OF_SCALARS_H]*temperature_error_model, 2);
    	}
    	// non-diagonal terms
    	for (int j = 1; j < no_of_edges + 1; ++j)
    	{
    		distance = normal_distance[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*h_index + j - 1]];
    		bg_error_cov_dry[i][j] = bg_error_cov_dry[i][0]*exp(-distance/e_folding_length);
    	}
    }
	
	/*
	Calling the optimum interpolation
	---------------------------------
	*/
	double *model_vector_dry = malloc((NO_OF_SCALARS + NO_OF_SCALARS_H)*sizeof(double));
	oi(obs_error_cov_dry, obs_op_jacobian_reduced_matrix_dry, relevant_model_dofs_matrix_dry, bg_error_cov_dry, interpolated_model_dry,
	background_dry, observations_vector_dry, model_vector_dry, NO_OF_CHOSEN_OBSERVATIONS_DRY, NO_OF_MODEL_DOFS_DRY);
	// freeing the memory
	free(obs_error_cov_dry);
	free(bg_error_cov_dry);
	free(interpolated_model_dry);
	free(background_dry);
	free(observations_vector_dry);
	
	// data assimilation is finished at this point
    
    // These are the arrays for the result of the assimilation process.
    double *density_dry = malloc(NO_OF_SCALARS*sizeof(double));
    double *wind = malloc(NO_OF_VECTORS*sizeof(double));
    double *exner = malloc(NO_OF_SCALARS*sizeof(double));
    
    // density is determined out of the hydrostatic equation
    double b, c;
    for (int i = NO_OF_SCALARS - 1; i >= 0; --i)
    {
    	layer_index = i/NO_OF_SCALARS_H;
    	h_index = i - layer_index*NO_OF_SCALARS_H;
    	// at the lowest layer the density is part of the model_vector_dry
    	if (layer_index == NO_OF_LAYERS - 1)
    	{
        	density_dry[i] = model_vector_dry[NO_OF_SCALARS + h_index];
        	exner[i] = pow((density_dry[i]*R_D*model_vector_dry[i])/P_0, R_D/C_D_P);
        }
        else
        {
			// solving a quadratic equation for the Exner pressure
			b = -0.5*exner[i + NO_OF_SCALARS_H]/model_vector_dry[i + NO_OF_SCALARS_H]
			*(model_vector_dry[i] - model_vector_dry[i + NO_OF_SCALARS_H]
			+ 2/C_D_P*(gravity_potential_model[i] - gravity_potential_model[i + NO_OF_SCALARS_H]));
			c = pow(exner[i + NO_OF_SCALARS_H], 2)*model_vector_dry[i]/model_vector_dry[i + NO_OF_SCALARS_H];
			exner[i] = b + pow((pow(b, 2) + c), 0.5);
        	density_dry[i] = P_0*pow(exner[i], C_D_P/R_D)/(R_D*model_vector_dry[i]);
        }
    }
	free(gravity_potential_model);
    free(exner);
    // end of the assimilation of the dry thermodynamic state
	printf("Dry assimilation completed.\n");
    
    /*
    MOISTURE ASSIMILATION
    ---------------------
    */
    
	printf("Starting the moist assimilation ...\n");
	
	// setting up the background error covariance matrix (only the diagonal)
	double (*bg_error_cov_moist)[7] = calloc(1, sizeof(double[NO_OF_MODEL_DOFS_MOIST][7]));
	double spec_hum_error_model = 0.01;
	#pragma omp parallel for private(distance, layer_index, h_index, no_of_edges)
	for (int i = 0; i < NO_OF_MODEL_DOFS_MOIST; ++i)
	{
    	// diagonal terms
		bg_error_cov_moist[i][0] = pow(spec_hum_error_model, 2);
		
    	layer_index = i/NO_OF_SCALARS_H;
    	h_index = i - layer_index*NO_OF_SCALARS_H;
    	if (layer_index == NO_OF_LAYERS)
    	{
    		layer_index = NO_OF_LAYERS - 1;
    	}
    	no_of_edges = 6;
    	if (h_index < NO_OF_PENTAGONS)
    	{
    		no_of_edges = 5;
    	}
    	// non-diagonal terms
    	for (int j = 1; j < no_of_edges + 1; ++j)
    	{
    		distance = normal_distance[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + adjacent_vector_indices_h[6*h_index + j - 1]];
    		bg_error_cov_moist[i][j] = bg_error_cov_moist[i][0]*exp(-distance/e_folding_length);
    	}
	}
    free(normal_distance);
    free(from_index);
    free(to_index);
    free(adjacent_vector_indices_h);
	
	// writing the background state into a single vector
    double *background_moist = malloc(NO_OF_MODEL_DOFS_MOIST*sizeof(double));
    // the data assimilation is being calculated with the specific humidity for pragmatic reasons
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_MODEL_DOFS_MOIST; ++i)
	{
		background_moist[i] = densities_background[5*NO_OF_SCALARS + i]/(densities_background[4*NO_OF_SCALARS + i] + densities_background[5*NO_OF_SCALARS + i]);
	}
	
	// setting up the observations operator
	double *interpolated_model_moist = malloc(NO_OF_CHOSEN_OBSERVATIONS_MOIST*sizeof(double));
	double (*obs_op_jacobian_reduced_matrix_moist)[NO_OF_REL_MODEL_DOFS_PER_OBS] = malloc(sizeof(double[NO_OF_CHOSEN_OBSERVATIONS_MOIST][NO_OF_REL_MODEL_DOFS_PER_OBS]));
	int (*relevant_model_dofs_matrix_moist)[NO_OF_REL_MODEL_DOFS_PER_OBS] = malloc(sizeof(int[NO_OF_CHOSEN_OBSERVATIONS_MOIST][NO_OF_REL_MODEL_DOFS_PER_OBS]));
	// setting up the moist observations operator using the dry observations operator
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS_MOIST; ++i)
	{
		interpolated_model_moist[i] = 0;
		for (int j = 0; j < NO_OF_REL_MODEL_DOFS_PER_OBS; ++j)
		{
			obs_op_jacobian_reduced_matrix_moist[i][j] = obs_op_jacobian_reduced_matrix_dry[i][j];
			relevant_model_dofs_matrix_moist[i][j] = relevant_model_dofs_matrix_dry[i][j];
			interpolated_model_moist[i] += obs_op_jacobian_reduced_matrix_moist[i][j]*background_moist[relevant_model_dofs_matrix_moist[i][j]];
		}
	}
	free(obs_op_jacobian_reduced_matrix_dry);
	free(relevant_model_dofs_matrix_dry);
	
    // writing the moist observations into a single vector
    double *observations_vector_moist = malloc(NO_OF_CHOSEN_OBSERVATIONS_MOIST*sizeof(double));
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS_MOIST; ++i)
	{
		observations_vector_moist[i] = observations_vector[NO_OF_CHOSEN_OBSERVATIONS_MOIST + i];
	}
	
	// setting up the measurement error covariance matrix
	double *obs_error_cov_moist = malloc(sizeof(double[NO_OF_CHOSEN_OBSERVATIONS_MOIST]));
	double spec_hum_error_obs = 0.0025;
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS_MOIST; ++i)
	{
		obs_error_cov_moist[i] = pow(spec_hum_error_obs, 2);
	}
	
	/*
	Calling the optimum interpolation
	---------------------------------
	*/
	double *model_vector_moist = malloc(NO_OF_MODEL_DOFS_MOIST*sizeof(double));
	oi(obs_error_cov_moist, obs_op_jacobian_reduced_matrix_moist, relevant_model_dofs_matrix_moist, bg_error_cov_moist, interpolated_model_moist,
	background_moist, observations_vector_moist, model_vector_moist, NO_OF_CHOSEN_OBSERVATIONS_MOIST, NO_OF_MODEL_DOFS_MOIST);
	
	// freeing arrays we do not need anymore
	free(obs_op_jacobian_reduced_matrix_moist);
	free(background_moist);
	free(relevant_model_dofs_matrix_moist);
	free(obs_error_cov_moist);
	free(bg_error_cov_moist);
	free(interpolated_model_moist);
	free(observations_vector_moist);
	printf("Moist assimilation completed.\n");
	
	/*
	WIND ASSIMILATION
	-----------------
	*/
	
	printf("Starting the wind assimilation ...\n");
	// writing the background state into a single vector
	double *background_wind = malloc(NO_OF_H_VECTORS*sizeof(double));
    #pragma omp parallel for private(layer_index, h_index)
	for (int i = 0; i < NO_OF_H_VECTORS; ++i)
	{
		layer_index = i/NO_OF_VECTORS_H;
		h_index = i - layer_index*NO_OF_VECTORS_H;
		background_wind[i] = wind_background[NO_OF_SCALARS_H + layer_index*NO_OF_VECTORS_PER_LAYER + h_index];
	}

	double *interpolated_model_wind = malloc(NO_OF_CHOSEN_OBSERVATIONS_WIND*sizeof(double));
	double (*obs_op_jacobian_reduced_matrix_wind)[NO_OF_REL_MODEL_DOFS_PER_OBS] = malloc(sizeof(double[NO_OF_CHOSEN_OBSERVATIONS_WIND][NO_OF_REL_MODEL_DOFS_PER_OBS]));
	int (*relevant_model_dofs_matrix_wind)[NO_OF_REL_MODEL_DOFS_PER_OBS] = malloc(sizeof(int[NO_OF_CHOSEN_OBSERVATIONS_WIND][NO_OF_REL_MODEL_DOFS_PER_OBS]));
	
	// setting up the obervations operator
	obs_op_setup_wind(interpolated_model_wind, obs_op_jacobian_reduced_matrix_wind, relevant_model_dofs_matrix_wind,
	latitude_vector_obs, longitude_vector_obs, z_coords_obs, latitudes_model_wind, longitudes_model_wind, z_coords_model_wind, directions, background_wind);
    
    free(z_coords_model_wind);
    free(directions);
	free(z_coords_obs);
    free(latitudes_model_wind);
    free(longitudes_model_wind);
	
    // writing the observations into one vector
	double *observations_vector_wind = malloc(NO_OF_CHOSEN_OBSERVATIONS_WIND*sizeof(double));
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS_WIND; ++i)
	{
		observations_vector_wind[i] = observations_vector[NO_OF_CHOSEN_OBSERVATIONS_DRY + NO_OF_CHOSEN_OBSERVATIONS_MOIST + i];
	}
	
	// setting the wind observations error covariance matrix
    double *obs_error_cov_wind = malloc(sizeof(double[NO_OF_CHOSEN_OBSERVATIONS_WIND]));
    double wind_error_obs = 0.5;
	#pragma omp parallel for
    for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS_WIND; ++i)
    {
		obs_error_cov_wind[i] = pow(wind_error_obs, 2);
    }
    
	// setting the background error covariance matrix
	double (*bg_error_cov_wind)[7] = calloc(1, sizeof(double[NO_OF_H_VECTORS][7]));
    double wind_error_model = 0.5;
	#pragma omp parallel for
    for (int i = 0; i < NO_OF_H_VECTORS; ++i)
    {
    	// diagonal terms
		bg_error_cov_wind[i][0] = pow(wind_error_model, 2);
    	
    	// non-diagonal terms
    	no_of_edges = 6;
    	for (int j = 1; j < no_of_edges + 1; ++j)
    	{
    		distance = 240e3;
    		bg_error_cov_wind[i][j] = bg_error_cov_wind[i][0]*exp(-distance/e_folding_length);
    	}
    }
    
	/*
	Calling the optimum interpolation
	---------------------------------
	*/
    
	double *model_vector_wind = malloc(NO_OF_H_VECTORS*sizeof(double));
	oi(obs_error_cov_wind, obs_op_jacobian_reduced_matrix_wind, relevant_model_dofs_matrix_wind, bg_error_cov_wind, interpolated_model_wind,
	background_wind, observations_vector_wind, model_vector_wind, NO_OF_CHOSEN_OBSERVATIONS_WIND, NO_OF_H_VECTORS);
	// freeing the memory
	free(obs_op_jacobian_reduced_matrix_wind);
	free(obs_error_cov_wind);
	free(bg_error_cov_wind);
	free(interpolated_model_wind);
	free(background_wind);
	free(observations_vector_wind);
	
	printf("Wind assimilation completed.\n");
	
	/*
	INTERPOLATION OF THE SST
	------------------------
	*/
	printf("Interpolating the SST to the model grid ...\n");
	double *sst = malloc(NO_OF_SCALARS_H*sizeof(double));
	int min_index;
	#pragma omp parallel for private(min_index)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
    	double *distance_vector = malloc(NO_OF_SST_POINTS*sizeof(double));
    	for (int j = 0; j < NO_OF_SST_POINTS; ++j)
    	{
    		distance_vector[j] = calculate_distance_h(latitude_vector_obs[NO_OF_CHOSEN_OBSERVATIONS - NO_OF_SST_POINTS + j],
    		longitude_vector_obs[NO_OF_CHOSEN_OBSERVATIONS - NO_OF_SST_POINTS + j], latitudes_model[i], longitudes_model[i], 1);
    	}
		min_index = find_min_index(distance_vector, NO_OF_SST_POINTS);
		sst[i] = observations_vector[NO_OF_CHOSEN_OBSERVATIONS - NO_OF_SST_POINTS + min_index];
		free(distance_vector);
	}
    free(latitudes_model);
    free(longitudes_model);
	free(latitude_vector_obs);
	free(longitude_vector_obs);
	free(observations_vector);
	
	printf("Interpolation of the SST completed.\n");

	/*
	PREPARING THE OUTPUT
	--------------------
	*/
	
	// individual condensate temperatures are for higher resolutions, not yet implemented
	// clouds and precipitation are set equal to the background state
	double *densities = malloc(6*NO_OF_SCALARS*sizeof(double));
	double *temperatures = malloc(5*NO_OF_SCALARS*sizeof(double));
    #pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		// setting the mass densities of the result
		// condensate densities are not assimilated
		densities[i] = densities_background[i];
		densities[NO_OF_SCALARS + i] = densities_background[NO_OF_SCALARS + i];
		densities[2*NO_OF_SCALARS + i] = densities_background[2*NO_OF_SCALARS + i];
		densities[3*NO_OF_SCALARS + i] = densities_background[3*NO_OF_SCALARS + i];
		densities[4*NO_OF_SCALARS + i] = density_dry[i];
		densities[5*NO_OF_SCALARS + i] = model_vector_moist[i]/(1 - model_vector_moist[i])*density_dry[i];
		if (densities[5*NO_OF_SCALARS + i] < 0)
		{
			densities[5*NO_OF_SCALARS + i] = 0;
		}
		// setting the temperatures of the result
		// assuming an LTE (local thermodynamic equilibrium)
		temperatures[i] = model_vector_dry[i];
		temperatures[NO_OF_SCALARS + i] = model_vector_dry[i];
		temperatures[2*NO_OF_SCALARS + i] = model_vector_dry[i];
		temperatures[3*NO_OF_SCALARS + i] = model_vector_dry[i];
		temperatures[4*NO_OF_SCALARS + i] = model_vector_dry[i];
    }
    free(model_vector_dry);
    free(density_dry);
	free(model_vector_moist);
	free(densities_background);
    
    // writing the result of the wind data assimilation to the resulting wind field
    #pragma omp parallel for private(layer_index, h_index)
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
    	layer_index = i/NO_OF_VECTORS_PER_LAYER;
    	h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
    	if (h_index < NO_OF_SCALARS_H)
    	{
    		wind[i] = wind_background[i];
    	}
    	else
    	{
    		wind[i] = model_vector_wind[layer_index*NO_OF_VECTORS_H + h_index - NO_OF_SCALARS_H];
    	}
    }
    free(model_vector_wind);
    free(wind_background);
    
    /*
    writing the result to a netcdf file
    -----------------------------------
    */
    
    printf("Output file: %s\n", output_file);
    printf("Writing result to output file ...\n");
    int densities_dimid, temperatures_dimid, vector_dimid, scalar_h_dimid, single_double_dimid, densities_id, temperatures_id, wind_id, sst_id;
    if ((retval = nc_create(output_file, NC_CLOBBER, &ncid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "densities_index", 6*NO_OF_SCALARS, &densities_dimid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "temperatures_index", 5*NO_OF_SCALARS, &temperatures_dimid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "vector_index", NO_OF_VECTORS, &vector_dimid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "scalar_h_index", NO_OF_SCALARS_H, &scalar_h_dimid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "single_double_dimid_index", 1, &single_double_dimid)))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "densities", NC_DOUBLE, 1, &densities_dimid, &densities_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, densities_id, "units", strlen("kg/m^3"), "kg/m^3")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "temperatures", NC_DOUBLE, 1, &temperatures_dimid, &temperatures_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, temperatures_id, "units", strlen("K"), "K")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "wind", NC_DOUBLE, 1, &vector_dimid, &wind_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, wind_id, "units", strlen("m/s"), "m/s")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "sst", NC_DOUBLE, 1, &scalar_h_dimid, &sst_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, sst_id, "units", strlen("K"), "K")))
        NCERR(retval);
    if ((retval = nc_enddef(ncid)))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, densities_id, &densities[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, temperatures_id, &temperatures[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, wind_id, &wind[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, sst_id, &sst[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid)))
    	NCERR(retval);
    printf("Result successfully written.\n");
    
    // freeing the stil occupied memory
	free(densities);
	free(temperatures);
	free(wind);
	free(sst);
	
	// that's it
    return 0;
}

int obs_op_setup(double interpolated_model[], double obs_op_jacobian_reduced_matrix[][NO_OF_REL_MODEL_DOFS_PER_OBS], int relevant_model_dofs_matrix[][NO_OF_REL_MODEL_DOFS_PER_OBS], double lat_used_obs[], double lon_used_obs[], double z_used_obs[], double lat_model[], double lon_model[], double z_model[], double background[])
{
	/*
	This functions calculates the observations operator.
	It is the background state, interpolated to the observations + the derivative of this function, which will be used to calculate
	the perturbation induced by the observations.
	*/
	
	// finding the NO_OF_REL_MODEL_DOFS_PER_OBS/2 closest grid points (horizontally) for each observation
	// the vector containing the relevant horizontal model indices for each observation
	int (*rel_h_index_vector)[NO_OF_REL_MODEL_DOFS_PER_OBS/2] = malloc(sizeof(int[NO_OF_CHOSEN_POINTS_PER_LAYER_OBS][NO_OF_REL_MODEL_DOFS_PER_OBS/2]));
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_CHOSEN_POINTS_PER_LAYER_OBS; ++i)
	{
		// the vector containing the horizontal distances between the observation at hand and each horizontal model gridpoint
		double *dist_vector = malloc(NO_OF_SCALARS_H*sizeof(double));
		// filling up the dist_vector
		for (int j = 0; j < NO_OF_SCALARS_H; ++j)
		{
			dist_vector[j] = calculate_distance_h(lat_used_obs[i], lon_used_obs[i], lat_model[j], lon_model[j], 1);
		}
		// finding the NO_OF_REL_MODEL_DOFS_PER_OBS/2 closest points
		for (int j = 0; j < NO_OF_REL_MODEL_DOFS_PER_OBS/2; ++j)
		{
			rel_h_index_vector[i][j] = find_min_index(dist_vector, NO_OF_SCALARS_H);
			dist_vector[rel_h_index_vector[i][j]] = M_PI + EPSILON;
		}
		free(dist_vector);
	}
	
	int layer_index, obs_index_h;
	// finally setting up the reduced observations operator
	#pragma omp parallel for private(layer_index, obs_index_h)
	for (int obs_index = 0; obs_index < NO_OF_CHOSEN_OBSERVATIONS_DRY; ++obs_index)
	{
		layer_index = obs_index/NO_OF_CHOSEN_POINTS_PER_LAYER_OBS;
		obs_index_h = obs_index - layer_index*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS;
		// the vector containing the vertical distance between the observation at hand and the model gridpoints
		double vert_distance_vector[NO_OF_LAYERS];
		// the vector containing preliminary interpolation weights
		double weights_vector[NO_OF_REL_MODEL_DOFS_PER_OBS];
		// the closest vertical indices
		int closest_vert_index, other_vert_index;
		double sum_of_interpol_weights, distance, closest_vert_weight, other_vert_weight;
		// free atmosphere quantities (temperature)
		if (obs_index < NO_OF_CHOSEN_OBSERVATIONS_DRY - NO_OF_CHOSEN_POINTS_PER_LAYER_OBS)
		{
			sum_of_interpol_weights = 0;
			interpolated_model[obs_index] = 0;
			// loop over all relevant horizontal model gridpoints
			for (int j = 0; j < NO_OF_REL_MODEL_DOFS_PER_OBS/2; ++j)
			{
				// finding out which layer is the closest to the observation
				for (int k = 0; k < NO_OF_LAYERS; ++k)
				{
					vert_distance_vector[k] = fabs(z_model[k*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j]] - z_used_obs[obs_index]);
				}
				closest_vert_index = find_min_index(vert_distance_vector, NO_OF_LAYERS);
				// vertical interpolation
				// firstly setting for the other vertical index
				other_vert_index = closest_vert_index + 1;
				// if the the closest model point is below the observation, the next higher point is taken into account for the interpolation
				if (z_model[closest_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j]] < z_used_obs[obs_index])
				{
					other_vert_index = closest_vert_index - 1;
				}
				// if the observation is below the lowest layer of the model
				if (other_vert_index == NO_OF_LAYERS)
				{
					other_vert_index = NO_OF_LAYERS - 2;
					closest_vert_weight = 1 - (z_used_obs[obs_index] - z_model[closest_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j]])
					/(z_model[other_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j]] - z_model[closest_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j]]);
				}
				else
				{
					closest_vert_weight = fabs(z_model[other_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j]] - z_used_obs[obs_index])
					/fabs(z_model[closest_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j]] - z_model[other_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j]]);
				}
				other_vert_weight = 1 - closest_vert_weight;
				// now we know which gridpoint is relevant to this observation
				// the closest vertical point
				relevant_model_dofs_matrix[obs_index][j] = closest_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j];
				// the second closest vertical point
				relevant_model_dofs_matrix[obs_index][j + NO_OF_REL_MODEL_DOFS_PER_OBS/2] = other_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j];
				// radius does not matter here
				distance = calculate_distance_h(lat_used_obs[obs_index], lon_used_obs[obs_index], lat_model[rel_h_index_vector[obs_index_h][j]], lon_model[rel_h_index_vector[obs_index_h][j]], 1);
				// 1/r-interpolation
				weights_vector[j] = closest_vert_weight/pow(distance + EPSILON, INTERPOL_EXP);
				weights_vector[j + NO_OF_REL_MODEL_DOFS_PER_OBS/2] = other_vert_weight/pow(distance + EPSILON, INTERPOL_EXP);
				interpolated_model[obs_index] += weights_vector[j]*background[relevant_model_dofs_matrix[obs_index][j]];
				interpolated_model[obs_index] += weights_vector[j + NO_OF_REL_MODEL_DOFS_PER_OBS/2]*background[relevant_model_dofs_matrix[obs_index][j + NO_OF_REL_MODEL_DOFS_PER_OBS/2]];
				sum_of_interpol_weights += weights_vector[j];
				sum_of_interpol_weights += weights_vector[j + NO_OF_REL_MODEL_DOFS_PER_OBS/2];
			}
			for (int k = 0; k < NO_OF_REL_MODEL_DOFS_PER_OBS; ++k)
			{
				// we have to divide by the sum of weights here
				obs_op_jacobian_reduced_matrix[obs_index][k] = weights_vector[k]/sum_of_interpol_weights;
			}
			interpolated_model[obs_index] = interpolated_model[obs_index]/sum_of_interpol_weights;
		}
		// surface quantities (only surface pressure for now)
		else
		{
			sum_of_interpol_weights = 0;
			interpolated_model[obs_index] = 0;
			// loop over all relevant model degrees of freedom
			for (int j = 0; j < NO_OF_REL_MODEL_DOFS_PER_OBS; ++j)
			{
				// we pick the lowest layer here (independant on wether we look at the temperature or the density)
				closest_vert_index = NO_OF_LAYERS - 1;
				// How is the suface pressure affected by the temperature in the lowest layer?
				if (j < NO_OF_REL_MODEL_DOFS_PER_OBS/2)
				{
					// radius does not matter here
					distance = calculate_distance_h(lat_used_obs[obs_index], lon_used_obs[obs_index],
					lat_model[rel_h_index_vector[obs_index_h][j]], lon_model[rel_h_index_vector[obs_index_h][j]], 1);
					// now we know which gridpoint is relevant to this observation
					relevant_model_dofs_matrix[obs_index][j] = closest_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j];
					// 1/r-interpolation
					weights_vector[j] = 1/pow(distance + EPSILON, INTERPOL_EXP)
					*R_D*background[relevant_model_dofs_matrix[obs_index][j] + NO_OF_SCALARS_H]
					*exp(-(z_used_obs[NO_OF_CHOSEN_OBSERVATIONS_MOIST + obs_index] - z_model[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j]])/SCALE_HEIGHT);
					sum_of_interpol_weights += 1/pow(distance + EPSILON, INTERPOL_EXP);
					// the result
					if (j == NO_OF_REL_MODEL_DOFS_PER_OBS/2 - 1)
					{
						// loop over all relevant gridpoints
						for (int k = 0; k < NO_OF_REL_MODEL_DOFS_PER_OBS/2; ++k)
						{
							// we have to divide by the sum of weights here
							obs_op_jacobian_reduced_matrix[obs_index][k] = weights_vector[k]/sum_of_interpol_weights;
						}
					}
				}
				// How is the surface pressure affected by the density in the lowest layer?
				else
				{
					// as a new interpolation will be conducted now, the sum_of_interpol_weights variable has to be reset to zero
					if (j == NO_OF_REL_MODEL_DOFS_PER_OBS/2)
					{
						sum_of_interpol_weights = 0;
					}
					// radius does not matter here
					distance = calculate_distance_h(lat_used_obs[obs_index], lon_used_obs[obs_index],
					lat_model[rel_h_index_vector[obs_index_h][j - NO_OF_REL_MODEL_DOFS_PER_OBS/2]], lon_model[rel_h_index_vector[obs_index_h][j - NO_OF_REL_MODEL_DOFS_PER_OBS/2]], 1);
					// now we know which gridpoint is relevant to this observation
					relevant_model_dofs_matrix[obs_index][j] = (closest_vert_index + 1)*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j - NO_OF_REL_MODEL_DOFS_PER_OBS/2];
					// interpolation weights
					weights_vector[j] = 1/pow(distance + EPSILON, INTERPOL_EXP)
					*R_D*background[relevant_model_dofs_matrix[obs_index][j] - NO_OF_SCALARS_H]
					*exp(-(z_used_obs[NO_OF_CHOSEN_OBSERVATIONS_MOIST + obs_index]
					- z_model[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j - NO_OF_REL_MODEL_DOFS_PER_OBS/2]])/SCALE_HEIGHT);
					// interpolation to the surface pressure
					interpolated_model[obs_index] += weights_vector[j]*background[relevant_model_dofs_matrix[obs_index][j]];
					sum_of_interpol_weights += 1/pow(distance + EPSILON, INTERPOL_EXP);
				}
			}
			// the result
			// the interpolation to the surface pressure
			interpolated_model[obs_index] = interpolated_model[obs_index]/sum_of_interpol_weights;
			// loop over all relevant gridpoints
			for (int k = NO_OF_REL_MODEL_DOFS_PER_OBS/2; k < NO_OF_REL_MODEL_DOFS_PER_OBS; ++k)
			{
				// we have to divide by the sum of weights here
				obs_op_jacobian_reduced_matrix[obs_index][k] = weights_vector[k]/sum_of_interpol_weights;
			}
		}
	}
	free(rel_h_index_vector);
	
	// returning 0 indicating success
	return 0;
}

int obs_op_setup_wind(double interpolated_model[], double obs_op_jacobian_reduced_matrix[][NO_OF_REL_MODEL_DOFS_PER_OBS], int relevant_model_dofs_matrix[][NO_OF_REL_MODEL_DOFS_PER_OBS], double lat_used_obs[], double lon_used_obs[], double z_used_obs[], double lat_model[], double lon_model[], double z_model[], double directions[], double background[])
{
	/*
	same as obs_op_setup, only for the wind
	*/
	
	// finding the NO_OF_REL_MODEL_DOFS_PER_OBS/2 closest grid points (horizontally) for each observation
	// the vector containing the relevant horizontal model indices for each observation
	int (*rel_h_index_vector)[NO_OF_REL_MODEL_DOFS_PER_OBS/2] = malloc(sizeof(int[NO_OF_CHOSEN_WIND_POINTS_PER_LAYER_OBS][NO_OF_REL_MODEL_DOFS_PER_OBS/2]));
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_CHOSEN_WIND_POINTS_PER_LAYER_OBS; ++i)
	{
		// the vector containing the horizontal distances between the observation at hand and each horizontal model gridpoint
		double *dist_vector = malloc(NO_OF_VECTORS_H*sizeof(double));
		// filling up the dist_vector
		for (int j = 0; j < NO_OF_VECTORS_H; ++j)
		{
			dist_vector[j] = calculate_distance_h(lat_used_obs[NO_OF_CHOSEN_OBSERVATIONS_DRY + NO_OF_CHOSEN_OBSERVATIONS_MOIST + i],
			lon_used_obs[NO_OF_CHOSEN_OBSERVATIONS_DRY + NO_OF_CHOSEN_OBSERVATIONS_MOIST + i], lat_model[j], lon_model[j], 1);
		}
		// finding the NO_OF_REL_MODEL_DOFS_PER_OBS/2 closest points
		for (int j = 0; j < NO_OF_REL_MODEL_DOFS_PER_OBS/2; ++j)
		{
			rel_h_index_vector[i][j] = find_min_index(dist_vector, NO_OF_VECTORS_H);
			dist_vector[rel_h_index_vector[i][j]] = M_PI + EPSILON;
		}
		free(dist_vector);
	}

	int layer_index, obs_index_h;
	// finally setting up the reduced observations operator
	#pragma omp parallel for private(layer_index, obs_index_h)
	for (int obs_index = 0; obs_index < NO_OF_CHOSEN_OBSERVATIONS_WIND; ++obs_index)
	{
		// just a helper variable for calculating obs_index_h; gets >= NO_OF_LAYERS but won't be used
		layer_index = obs_index/NO_OF_CHOSEN_WIND_POINTS_PER_LAYER_OBS;
		obs_index_h = obs_index - layer_index*NO_OF_CHOSEN_WIND_POINTS_PER_LAYER_OBS;
		// the vector containing the vertical distance between the observation at hand and the model gridpoints
		double vert_distance_vector[NO_OF_LAYERS];
		// the vector containing preliminary interpolation weights
		double weights_vector[NO_OF_REL_MODEL_DOFS_PER_OBS];
		// the vector containing preliminary interpolation weights with sin or cos prefactors accounting for the direction of the
		// vectors on the C-grid
		double weights_vector_with_dir[NO_OF_REL_MODEL_DOFS_PER_OBS];
		// the closest vertical indices
		int closest_vert_index, other_vert_index;
		double sum_of_interpol_weights, distance, closest_vert_weight, other_vert_weight, distance_factor, direction;
		sum_of_interpol_weights = 0;
		interpolated_model[obs_index] = 0;
		// loop over all relevant horizontal model gridpoints
		for (int j = 0; j < NO_OF_REL_MODEL_DOFS_PER_OBS/2; ++j)
		{
			// finding out which layer is the closest to the observation
			for (int k = 0; k < NO_OF_LAYERS; ++k)
			{
				vert_distance_vector[k] = fabs(z_model[NO_OF_SCALARS_H + k*NO_OF_VECTORS_PER_LAYER + rel_h_index_vector[obs_index_h][j]]
				- z_used_obs[NO_OF_CHOSEN_OBSERVATIONS_DRY + NO_OF_CHOSEN_OBSERVATIONS_MOIST + obs_index]);
			}
			closest_vert_index = find_min_index(vert_distance_vector, NO_OF_LAYERS);
			// vertical interpolation
			// firstly setting for the other vertical index
			other_vert_index = closest_vert_index + 1;
			// if the the closest model point is below the observation, the next higher point is taken into account for the interpolation
			if (z_model[NO_OF_SCALARS_H + closest_vert_index*NO_OF_VECTORS_PER_LAYER + rel_h_index_vector[obs_index_h][j]]
			< z_used_obs[NO_OF_CHOSEN_OBSERVATIONS_DRY + NO_OF_CHOSEN_OBSERVATIONS_MOIST + obs_index])
			{
				other_vert_index = closest_vert_index - 1;
			}
			// if the observation is below the lowest layer of the model
			if (other_vert_index == NO_OF_LAYERS)
			{
				other_vert_index = NO_OF_LAYERS - 2;
				closest_vert_weight = 1 - (z_used_obs[NO_OF_CHOSEN_OBSERVATIONS_DRY + NO_OF_CHOSEN_OBSERVATIONS_MOIST + obs_index]
				- z_model[NO_OF_SCALARS_H + closest_vert_index*NO_OF_VECTORS_PER_LAYER + rel_h_index_vector[obs_index_h][j]])
				/(z_model[NO_OF_SCALARS_H + other_vert_index*NO_OF_VECTORS_PER_LAYER + rel_h_index_vector[obs_index_h][j]]
				- z_model[NO_OF_SCALARS_H + closest_vert_index*NO_OF_VECTORS_PER_LAYER + rel_h_index_vector[obs_index_h][j]]);
			}
			else
			{
				closest_vert_weight = fabs(z_model[NO_OF_SCALARS_H + other_vert_index*NO_OF_VECTORS_PER_LAYER + rel_h_index_vector[obs_index_h][j]]
				- z_used_obs[NO_OF_CHOSEN_OBSERVATIONS_DRY + NO_OF_CHOSEN_OBSERVATIONS_MOIST + obs_index])
				/fabs(z_model[NO_OF_SCALARS_H + closest_vert_index*NO_OF_VECTORS_PER_LAYER + rel_h_index_vector[obs_index_h][j]]
				- z_model[NO_OF_SCALARS_H + other_vert_index*NO_OF_VECTORS_PER_LAYER + rel_h_index_vector[obs_index_h][j]]);
			}
			other_vert_weight = 1 - closest_vert_weight;
			// now we know which gridpoint is relevant to this observation
			// the closest vertical point
			relevant_model_dofs_matrix[obs_index][j] = closest_vert_index*NO_OF_VECTORS_H + rel_h_index_vector[obs_index_h][j];
			// the second closest vertical point
			relevant_model_dofs_matrix[obs_index][j + NO_OF_REL_MODEL_DOFS_PER_OBS/2] = other_vert_index*NO_OF_VECTORS_H + rel_h_index_vector[obs_index_h][j];
			// radius does not matter here
			distance = calculate_distance_h(lat_used_obs[NO_OF_CHOSEN_OBSERVATIONS_DRY + NO_OF_CHOSEN_OBSERVATIONS_MOIST + obs_index],
			lon_used_obs[NO_OF_CHOSEN_OBSERVATIONS_DRY + NO_OF_CHOSEN_OBSERVATIONS_MOIST + obs_index],
			lat_model[rel_h_index_vector[obs_index_h][j]], lon_model[rel_h_index_vector[obs_index_h][j]], 1);
			// interpolation weights
			distance_factor = 1.0/pow(distance + EPSILON, INTERPOL_EXP);
			direction = directions[rel_h_index_vector[obs_index_h][j]];
			// u
			if (obs_index < NO_OF_CHOSEN_OBSERVATIONS_WIND/2)
			{
				weights_vector[j] = closest_vert_weight*distance_factor;
				weights_vector[j + NO_OF_REL_MODEL_DOFS_PER_OBS/2] = other_vert_weight*distance_factor;
				weights_vector_with_dir[j] = weights_vector[j]*cos(direction);
				weights_vector_with_dir[j + NO_OF_REL_MODEL_DOFS_PER_OBS/2] = weights_vector[j + NO_OF_REL_MODEL_DOFS_PER_OBS/2]*cos(direction);
			}
			// v
			else
			{
				weights_vector[j] = closest_vert_weight*distance_factor;
				weights_vector[j + NO_OF_REL_MODEL_DOFS_PER_OBS/2] = other_vert_weight*distance_factor;
				weights_vector_with_dir[j] = weights_vector[j]*sin(direction);
				weights_vector_with_dir[j + NO_OF_REL_MODEL_DOFS_PER_OBS/2] = weights_vector[j + NO_OF_REL_MODEL_DOFS_PER_OBS/2]*sin(direction);
			}
			interpolated_model[obs_index] += weights_vector_with_dir[j]*background[relevant_model_dofs_matrix[obs_index][j]];
			interpolated_model[obs_index] += weights_vector_with_dir[j + NO_OF_REL_MODEL_DOFS_PER_OBS/2]*background[relevant_model_dofs_matrix[obs_index][j + NO_OF_REL_MODEL_DOFS_PER_OBS/2]];
			sum_of_interpol_weights += weights_vector[j];
			sum_of_interpol_weights += weights_vector[j + NO_OF_REL_MODEL_DOFS_PER_OBS/2];
		}
		for (int k = 0; k < NO_OF_REL_MODEL_DOFS_PER_OBS; ++k)
		{
			// we have to divide by the sum of weights here
			obs_op_jacobian_reduced_matrix[obs_index][k] = weights_vector_with_dir[k]/(0.5*sum_of_interpol_weights);
		}
		interpolated_model[obs_index] = interpolated_model[obs_index]/(0.5*sum_of_interpol_weights);
	}
	free(rel_h_index_vector);
	
	// returning 0 indicating success
	return 0;
}

















