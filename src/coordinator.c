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
#include "geos95.h"
#include "atmostracers.h"
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}
#define EPSILON 1e-4
#define SCALE_HEIGHT 8000.0
#define P_0 100000

int obs_op_setup(double [], double [][NO_OF_REL_MODEL_DOFS_PER_OBS], int [][NO_OF_REL_MODEL_DOFS_PER_OBS], double [], double [], double [], double [], double [], double [], double []);

int main(int argc, char *argv[])
{
	if (fmod(NO_OF_REL_MODEL_DOFS_PER_OBS, 2) == 1)
	{
		printf("NO_OF_REL_MODEL_DOFS_PER_OBS must be even.\n");
		printf("Aborting.\n");
		exit(1);
	}
	double R_D = specific_gas_constants_lookup(0);
	double C_D_P = spec_heat_capacities_p_gas_lookup(0);
    size_t len = strlen(argv[1]);
    char *year_string = malloc((len + 1)*sizeof(char));
    strcpy(year_string, argv[1]);
    len = strlen(argv[2]);
    char *month_string = malloc((len + 1)*sizeof(char));
    strcpy(month_string, argv[2]);
    len = strlen(argv[3]);
    char *day_string = malloc((len + 1)*sizeof(char));
    strcpy(day_string, argv[3]);
    len = strlen(argv[4]);
    char *hour_string = malloc((len + 1)*sizeof(char));
    strcpy(hour_string, argv[4]);
    len = strlen(argv[5]);
    char *model_home_dir = malloc((len + 1)*sizeof(char));
    strcpy(model_home_dir, argv[5]);
	int ORO_ID;
    ORO_ID = strtod(argv[6], NULL);
    len = strlen(argv[7]);
    char *BACKGROUND_STATE_FILE = malloc((len + 1)*sizeof(char));
    strcpy(BACKGROUND_STATE_FILE, argv[7]);
    len = strlen(argv[8]);
    char *game_da_root_dir = malloc((len + 1)*sizeof(char));
    strcpy(game_da_root_dir, argv[8]);
	int NO_OF_ORO_LAYERS;
    NO_OF_ORO_LAYERS = strtod(argv[9], NULL);
	int TOA = strtod(argv[10], NULL);
	printf("background state file: %s\n", BACKGROUND_STATE_FILE);
    
    // Allocating memory for the grid properties.
    double *latitudes_model = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *longitudes_model = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *z_coords_model = malloc(NO_OF_SCALARS*sizeof(double));
    double *gravity_potential_model = malloc(NO_OF_SCALARS*sizeof(double));
    
    // Reading the grid properties.
    int ncid_grid, retval;
    char GEO_PROP_FILE_PRE[200];
    sprintf(GEO_PROP_FILE_PRE, "%s/grid_generator/grids/B%dL%dT%d_O%d_OL%d_SCVT.nc", model_home_dir, RES_ID, NO_OF_LAYERS, TOA, ORO_ID, NO_OF_ORO_LAYERS);
    char GEO_PROP_FILE[strlen(GEO_PROP_FILE_PRE) + 1];
    strcpy(GEO_PROP_FILE, GEO_PROP_FILE_PRE);
	printf("grid file: %s\n", GEO_PROP_FILE);
	printf("reading grid file ...\n");
    if ((retval = nc_open(GEO_PROP_FILE, NC_NOWRITE, &ncid_grid)))
        NCERR(retval);
    int latitudes_model_id, longitudes_model_id, z_coords_model_id, gravity_potential_model_id, stretching_parameter_grid_id;
    double stretching_parameter_grid;
    if ((retval = nc_inq_varid(ncid_grid, "latitude_scalar", &latitudes_model_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "longitude_scalar", &longitudes_model_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "z_scalar", &z_coords_model_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "gravity_potential", &gravity_potential_model_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "stretching_parameter", &stretching_parameter_grid_id)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, stretching_parameter_grid_id, &stretching_parameter_grid)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, latitudes_model_id, &latitudes_model[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, longitudes_model_id, &longitudes_model[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, z_coords_model_id, &z_coords_model[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, gravity_potential_model_id, &gravity_potential_model[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid_grid)))
        NCERR(retval);
	printf("Grid file read.\n");
	
    char OUTPUT_FILE_PRE[200];
    sprintf(OUTPUT_FILE_PRE, "%s/nwp_init/%s%s%s%s_B%dL%dT%d_O%d_OL%d_SCVT.nc", model_home_dir, year_string, month_string, day_string, hour_string, RES_ID, NO_OF_LAYERS, TOA, ORO_ID, NO_OF_ORO_LAYERS);
    free(model_home_dir);
    char OUTPUT_FILE[strlen(OUTPUT_FILE_PRE) + 1];
    strcpy(OUTPUT_FILE, OUTPUT_FILE_PRE);
    
    // These are the arrays of the background state.
    double *temperature_gas_background = malloc(NO_OF_SCALARS*sizeof(double));
    double *density_dry_background = malloc(NO_OF_SCALARS*sizeof(double));
    double *wind_background = malloc(NO_OF_VECTORS*sizeof(double));
    double *water_vapour_density_background = malloc(NO_OF_SCALARS*sizeof(double));
    double *liquid_water_density_background = malloc(NO_OF_SCALARS*sizeof(double));
    double *solid_water_density_background = malloc(NO_OF_SCALARS*sizeof(double));
    double *liquid_water_temperature_background = malloc(NO_OF_SCALARS*sizeof(double));
    double *solid_water_temperature_background = malloc(NO_OF_SCALARS*sizeof(double));
    
    // Reading the background state.
	printf("reading background state ...\n");
    int ncid;
    if ((retval = nc_open(BACKGROUND_STATE_FILE, NC_NOWRITE, &ncid)))
        NCERR(retval);
    free(BACKGROUND_STATE_FILE);
    int temperature_gas_background_id, density_dry_background_id, wind_background_id, density_vapour_background_id, density_liquid_background_id, density_solid_background_id, temperature_liquid_background_id, temperature_solid_background_id, stretching_parameter_state_id;
    double stretching_parameter_state;
    if ((retval = nc_inq_varid(ncid, "density_dry", &density_dry_background_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "temperature_gas",&temperature_gas_background_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "wind", &wind_background_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "density_vapour", &density_vapour_background_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "density_liquid", &density_liquid_background_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "density_solid", &density_solid_background_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "temperature_liquid", &temperature_liquid_background_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "temperature_solid", &temperature_solid_background_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "stretching_parameter", &stretching_parameter_state_id)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, stretching_parameter_state_id, &stretching_parameter_state)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, temperature_gas_background_id, &temperature_gas_background[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, density_dry_background_id, &density_dry_background[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, wind_background_id, &wind_background[0])))
        NCERR(retval);    
    if ((retval = nc_get_var_double(ncid, density_vapour_background_id, &water_vapour_density_background[0])))
        NCERR(retval);    
    if ((retval = nc_get_var_double(ncid, density_liquid_background_id, &liquid_water_density_background[0])))
        NCERR(retval);    
    if ((retval = nc_get_var_double(ncid, density_solid_background_id, &solid_water_density_background[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, temperature_liquid_background_id, &liquid_water_temperature_background[0])))
        NCERR(retval);    
    if ((retval = nc_get_var_double(ncid, temperature_solid_background_id, &solid_water_temperature_background[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid)))
        NCERR(retval);
	printf("Background state read.\n");

	// saving the relevant part of the background state in one array
	double *background_dry = malloc(NO_OF_MODEL_DOFS_DRY*sizeof(double));
    #pragma omp parallel for
	for (int i = 0; i < NO_OF_MODEL_DOFS_DRY; ++i)
	{
		if (i < NO_OF_SCALARS)
		{
			background_dry[i] = temperature_gas_background[i];
		}
		else
		{
			background_dry[i] = density_dry_background[i - NO_OF_SCALARS_H];
		}
	}
    free(temperature_gas_background);

	// Comparing the stretching parameters of the grid and the background state.
	if (stretching_parameter_grid != stretching_parameter_state)
	{
		printf("stretching_parameters of grid and background state do not conform.\n");
		printf("Aborting.\n");
		exit(1);
	}
	double stretching_parameter = stretching_parameter_grid;	
	
	// Allocating the memory for the observations.
	double *latitude_vector_obs = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(double));
	double *longitude_vector_obs = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(double));
	double *z_coords_obs = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(double));
	double *observations_vector = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(double));
    
    char OBSERVATIONS_FILE_PRE[200];
    sprintf(OBSERVATIONS_FILE_PRE, "%s/input/obs_%s%s%s%s.nc", game_da_root_dir, year_string, month_string, day_string, hour_string);
	free(game_da_root_dir);
	free(year_string);
	free(month_string);
	free(day_string);
	free(hour_string);
    char OBSERVATIONS_FILE[strlen(OBSERVATIONS_FILE_PRE) + 1];
    strcpy(OBSERVATIONS_FILE, OBSERVATIONS_FILE_PRE);
	printf("observations file: %s\n", OBSERVATIONS_FILE);
    
    // Reading the observations.
	printf("reading observations ...\n");
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
    
    // setting up the measurement error covariance matrix
    double *obs_error_cov_dry = malloc(sizeof(double[NO_OF_CHOSEN_OBSERVATIONS_DRY]));
    double temperature_error_obs = 1;
    double pressure_error_obs = 100;
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
    double *bg_error_cov_dry = malloc(NO_OF_MODEL_DOFS_DRY*sizeof(double));
    double temperature_error_model = 1;
    double pressure_error_model = 100;
    #pragma omp parallel for
    for (int i = 0; i < NO_OF_MODEL_DOFS_DRY; ++i)
    {
    	if (i < NO_OF_SCALARS)
    	{
			bg_error_cov_dry[i] = pow(temperature_error_model, 2);
    	}
    	else
    	{
    		// density = p/(R_D*T) (Gauss'ian error propagation)
    		bg_error_cov_dry[i] = pow(pressure_error_model/(R_D*background_dry[i - NO_OF_SCALARS_H]), 2) + pow(background_dry[i]/background_dry[i - NO_OF_SCALARS_H]*temperature_error_model, 2);
    	}
    }
    
	// setting up the observations operator
	double *interpolated_model_dry = malloc(NO_OF_CHOSEN_OBSERVATIONS_DRY*sizeof(double));
	double (*obs_op_jacobian_reduced_matrix_dry)[NO_OF_REL_MODEL_DOFS_PER_OBS] = malloc(sizeof(double[NO_OF_CHOSEN_OBSERVATIONS_DRY][NO_OF_REL_MODEL_DOFS_PER_OBS]));
	int (*relevant_model_dofs_matrix_dry)[NO_OF_REL_MODEL_DOFS_PER_OBS] = malloc(sizeof(int[NO_OF_CHOSEN_OBSERVATIONS_DRY][NO_OF_REL_MODEL_DOFS_PER_OBS]));
	obs_op_setup(interpolated_model_dry, obs_op_jacobian_reduced_matrix_dry, relevant_model_dofs_matrix_dry,
	latitude_vector_obs, longitude_vector_obs, z_coords_obs, latitudes_model, longitudes_model, z_coords_model, background_dry);
    free(z_coords_model);
    free(latitudes_model);
    free(longitudes_model);
	free(latitude_vector_obs);
	free(longitude_vector_obs);
	free(z_coords_obs);
	
	double *observations_vector_dry = malloc(NO_OF_CHOSEN_OBSERVATIONS_DRY*sizeof(double));
	// setting up the dry observations vector
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
	
	// now, all the constituents of the gain matrix are known
	double *model_vector_dry = malloc((NO_OF_SCALARS + NO_OF_SCALARS_H)*sizeof(double));
	oi(obs_error_cov_dry, obs_op_jacobian_reduced_matrix_dry, relevant_model_dofs_matrix_dry, bg_error_cov_dry, interpolated_model_dry,
	background_dry, observations_vector_dry, model_vector_dry, NO_OF_CHOSEN_OBSERVATIONS_DRY, NO_OF_MODEL_DOFS_DRY);
	
	// data assimilation is finished at this point
	// freeing the memory
	free(obs_error_cov_dry);
	free(bg_error_cov_dry);
	free(interpolated_model_dry);
	free(background_dry);
	free(observations_vector_dry);
    
    // These are the arrays for the result of the assimilation process.
    double *temperature = malloc(NO_OF_SCALARS*sizeof(double));
    double *density_dry = malloc(NO_OF_SCALARS*sizeof(double));
    double *wind = malloc(NO_OF_VECTORS*sizeof(double));
    double *water_vapour_density = malloc(NO_OF_SCALARS*sizeof(double));
    double *liquid_water_density = malloc(NO_OF_SCALARS*sizeof(double));
    double *solid_water_density = malloc(NO_OF_SCALARS*sizeof(double));
    double *liquid_water_temp = malloc(NO_OF_SCALARS*sizeof(double));
    double *solid_water_temp = malloc(NO_OF_SCALARS*sizeof(double));
    double *exner = malloc(NO_OF_SCALARS*sizeof(double));
    
    // the temperature comes first in the model_vector_dry
    #pragma omp parallel for
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	temperature[i] = model_vector_dry[i];
    }
    
    // density is determined out of the hydrostatic equation
    int layer_index, h_index;
    double b, c;
    for (int i = NO_OF_SCALARS - 1; i >= 0; --i)
    {
    	layer_index = i/NO_OF_SCALARS_H;
    	h_index = i - layer_index*NO_OF_SCALARS_H;
    	// at the lowest layer the density is part of the model_vector_dry
    	if (layer_index == NO_OF_LAYERS - 1)
    	{
        	density_dry[i] = model_vector_dry[NO_OF_SCALARS + h_index];
        	exner[i] = pow((density_dry[i]*R_D*temperature[i])/P_0, R_D/C_D_P);
        }
        else
        {
			// solving a quadratic equation for the Exner pressure
			b = -0.5*exner[i + NO_OF_SCALARS_H]/temperature[i + NO_OF_SCALARS_H]
			*(temperature[i] - temperature[i + NO_OF_SCALARS_H]
			+ 2/C_D_P*(gravity_potential_model[i] - gravity_potential_model[i + NO_OF_SCALARS_H]));
			c = pow(exner[i + NO_OF_SCALARS_H], 2)*temperature[i]/temperature[i + NO_OF_SCALARS_H];
			exner[i] = b + pow((pow(b, 2) + c), 0.5);
        	density_dry[i] = P_0*pow(exner[i], C_D_P/R_D)/(R_D*temperature[i]);
        }
    }
	free(gravity_potential_model);
    free(exner);
    free(model_vector_dry);
    
    // Wind is set equal to the background wind for now. Later it will be derived from the balance equation.
    #pragma omp parallel for
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
    	wind[i] = wind_background[i];
    }
    free(wind_background);
    // end of the assimilation of the dry state
    
    // separate moisture assimilation
    double *observations_vector_moist = malloc(NO_OF_CHOSEN_OBSERVATIONS_MOIST*sizeof(double));
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS_MOIST; ++i)
	{
		observations_vector_moist[i] = observations_vector[NO_OF_CHOSEN_OBSERVATIONS_MOIST + i];
	}
	free(observations_vector);
	
    double *background_moist = malloc(NO_OF_MODEL_DOFS_MOIST*sizeof(double));
    // the data assimilation is being calculated with the specific humidity for pragmatic reasons
	for (int i = 0; i < NO_OF_MODEL_DOFS_MOIST; ++i)
	{
		background_moist[i] = water_vapour_density_background[i]/(density_dry_background[i] + water_vapour_density_background[i]);
	}
    free(density_dry_background);
    free(water_vapour_density_background);
	
	// setting up the measurement error covariance matrix
	double *obs_error_cov_moist = malloc(sizeof(double[NO_OF_CHOSEN_OBSERVATIONS_MOIST]));
	double abs_moisture_error_obs = 0.005;
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS_MOIST; ++i)
	{
		obs_error_cov_moist[i] = pow(abs_moisture_error_obs, 2);
	}
	
	// setting up the background error covariance matrix (only the diagonal)
	double *bg_error_cov_moist = malloc(NO_OF_MODEL_DOFS_MOIST*sizeof(double));
	double abs_moisture_error_model = 0.005;
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_MODEL_DOFS_MOIST; ++i)
	{
		bg_error_cov_moist[i] = pow(abs_moisture_error_model, 2);
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
	
	// now, all the constituents of the gain matrix are known
	double *model_vector_moist = malloc(NO_OF_MODEL_DOFS_MOIST*sizeof(double));
	oi(obs_error_cov_moist, obs_op_jacobian_reduced_matrix_moist, relevant_model_dofs_matrix_moist, bg_error_cov_moist, interpolated_model_moist,
	background_moist, observations_vector_moist, model_vector_moist, NO_OF_CHOSEN_OBSERVATIONS_MOIST, NO_OF_MODEL_DOFS_MOIST);
	
	free(obs_op_jacobian_reduced_matrix_moist);
	free(background_moist);
	free(relevant_model_dofs_matrix_moist);
	free(obs_error_cov_moist);
	free(bg_error_cov_moist);
	free(interpolated_model_moist);
	free(observations_vector_moist);
	
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		// cloud water density
		water_vapour_density[i] = model_vector_moist[i]/(1 - model_vector_moist[i])*density_dry[i];
		if (water_vapour_density[i] < 0)
		{
			water_vapour_density[i] = 0;
		}
	}
	free(model_vector_moist);
	
	// individual condensate temperatures are for higher resolutions, not yet implemented
	// clouds and precipitation are set equal to the background state
    #pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		solid_water_density[i] = solid_water_density_background[i];
		liquid_water_density[i] = liquid_water_density_background[i];
		solid_water_temp[i] = solid_water_temperature_background[i];
		liquid_water_temp[i] = liquid_water_temperature_background[i];
    }
    free(solid_water_density_background);
    free(liquid_water_density_background);
    free(solid_water_temperature_background);
    free(liquid_water_temperature_background);
    
    // Writing the result to a netcdf file.
    printf("output file: %s\n", OUTPUT_FILE);
    printf("writing result to output file ...\n");
    int scalar_dimid, vector_dimid, temp_id, density_dry_id, wind_id, density_vapour_id, density_liquid_id, density_solid_id, temperature_liquid_id, temperature_solid_id, single_double_dimid, stretching_parameter_id;
    if ((retval = nc_create(OUTPUT_FILE, NC_CLOBBER, &ncid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "scalar_index", NO_OF_SCALARS, &scalar_dimid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "vector_index", NO_OF_VECTORS, &vector_dimid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "single_double_dimid_index", 1, &single_double_dimid)))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "stretching_parameter", NC_DOUBLE, 1, &single_double_dimid, &stretching_parameter_id)))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "temperature_gas", NC_DOUBLE, 1, &scalar_dimid, &temp_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, temp_id, "units", strlen("K"), "K")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "density_dry", NC_DOUBLE, 1, &scalar_dimid, &density_dry_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, density_dry_id, "units", strlen("kg/m^3"), "kg/m^3")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "wind", NC_DOUBLE, 1, &vector_dimid, &wind_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, wind_id, "units", strlen("m/s"), "m/s")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "density_vapour", NC_DOUBLE, 1, &scalar_dimid, &density_vapour_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, density_vapour_id, "units", strlen("kg/m^3"), "kg/m^3")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "density_liquid", NC_DOUBLE, 1, &scalar_dimid, &density_liquid_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, density_liquid_id, "units", strlen("kg/m^3"), "kg/m^3")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "density_solid", NC_DOUBLE, 1, &scalar_dimid, &density_solid_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, density_solid_id, "units", strlen("kg/m^3"), "kg/m^3")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "temperature_liquid", NC_DOUBLE, 1, &scalar_dimid, &temperature_liquid_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, temperature_liquid_id, "units", strlen("T"), "T")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "temperature_solid", NC_DOUBLE, 1, &scalar_dimid, &temperature_solid_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, temperature_solid_id, "units", strlen("T"), "T")))
        NCERR(retval);
    if ((retval = nc_enddef(ncid)))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, stretching_parameter_id, &stretching_parameter)))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, temp_id, &temperature[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, density_dry_id, &density_dry[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, wind_id, &wind[0])))
        NCERR(retval);    
    if ((retval = nc_put_var_double(ncid, density_vapour_id, &water_vapour_density[0])))
        NCERR(retval);    
    if ((retval = nc_put_var_double(ncid, density_liquid_id, &liquid_water_density[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, density_solid_id, &solid_water_density[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, temperature_liquid_id, &liquid_water_temp[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, temperature_solid_id, &solid_water_temp[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid)))
    	NCERR(retval);
    printf("Result successfully written.\n");
	free(density_dry);
	free(temperature);
	free(wind);
	free(water_vapour_density);
	free(solid_water_density);
	free(liquid_water_density);
	free(solid_water_temp);
	free(liquid_water_temp);
    return 0;
}

int obs_op_setup(double interpolated_model_dry[], double obs_op_jacobian_reduced_matrix[][NO_OF_REL_MODEL_DOFS_PER_OBS], int relevant_model_dofs_matrix_dry[][NO_OF_REL_MODEL_DOFS_PER_OBS], double lat_used_obs[], double lon_used_obs[], double z_used_obs[], double lat_model[], double lon_model[], double z_model[], double background[])
{
	/*
	This functions calculates the observations operator.
	It is the background state, interpolated to the observations
	+ the derivative of this function, which will be used to calculate
	the perturbation induced by the observations.
	*/
	double R_D = specific_gas_constants_lookup(0);
	// finding the NO_OF_REL_MODEL_DOFS_PER_OBS closest grid points (horizontally) for each observation
	int (*rel_h_index_vector)[NO_OF_REL_MODEL_DOFS_PER_OBS/2] = malloc(sizeof(int[NO_OF_CHOSEN_POINTS_PER_LAYER_OBS][NO_OF_REL_MODEL_DOFS_PER_OBS/2])); // the vector containing the relevant horizontal model indices for each observation
	double *dist_vector = malloc(NO_OF_SCALARS_H*sizeof(double)); // the vector containing the horizontal distances between the observation at hand and each horizontal model gridpoint
	for (int i = 0; i < NO_OF_CHOSEN_POINTS_PER_LAYER_OBS; ++i)
	{
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
	}
	// dist_vector is no longer needed
	free(dist_vector);
	
	int layer_index, obs_index_h;
	// finally setting up the reduced observations operator
	for (int obs_index = 0; obs_index < NO_OF_CHOSEN_OBSERVATIONS_DRY; ++obs_index)
	{
		layer_index = obs_index/NO_OF_CHOSEN_POINTS_PER_LAYER_OBS;
		obs_index_h = obs_index - layer_index*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS;
		// the vector containing the vertical distance between the observation at hand and the model gridpoints
		double vert_distance_vector[NO_OF_LAYERS];
		// the vector containing preliminary interpolation weights
		double weights_vector[NO_OF_REL_MODEL_DOFS_PER_OBS];
		// the closest vertical index
		int closest_vert_index, other_vert_index;
		double sum_of_interpol_weights, distance, closest_vert_weight, other_vert_weight;
		// free atmosphere quantities (temperature)
		if (obs_index < NO_OF_CHOSEN_OBSERVATIONS_DRY - NO_OF_CHOSEN_POINTS_PER_LAYER_OBS)
		{
			sum_of_interpol_weights = 0;
			interpolated_model_dry[obs_index] = 0;
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
				// first setting for the other vertical index
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
					closest_vert_weight = fabs(z_model[other_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j]] - z_used_obs[obs_index])/fabs(z_model[closest_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j]] - z_model[other_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j]]);
				}
				other_vert_weight = 1 - closest_vert_weight;
				// now we know which gridpoint is relevant to this observation
				// the closest vertical point
				relevant_model_dofs_matrix_dry[obs_index][j] = closest_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j];
				// the second closest vertical point
				relevant_model_dofs_matrix_dry[obs_index][j + NO_OF_REL_MODEL_DOFS_PER_OBS/2] = other_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j];
				// radius does not matter here
				distance = calculate_distance_h(lat_used_obs[obs_index], lon_used_obs[obs_index], lat_model[rel_h_index_vector[obs_index_h][j]], lon_model[rel_h_index_vector[obs_index_h][j]], 1);
				// 1/r-interpolation
				weights_vector[j] = closest_vert_weight/pow(distance + EPSILON, T_INTERPOL_EXP);
				weights_vector[j + NO_OF_REL_MODEL_DOFS_PER_OBS/2] = other_vert_weight/pow(distance + EPSILON, T_INTERPOL_EXP);
				interpolated_model_dry[obs_index] += weights_vector[j]*background[relevant_model_dofs_matrix_dry[obs_index][j]];
				interpolated_model_dry[obs_index] += weights_vector[j + NO_OF_REL_MODEL_DOFS_PER_OBS/2]*background[relevant_model_dofs_matrix_dry[obs_index][j + NO_OF_REL_MODEL_DOFS_PER_OBS/2]];
				sum_of_interpol_weights += weights_vector[j];
				sum_of_interpol_weights += weights_vector[j + NO_OF_REL_MODEL_DOFS_PER_OBS/2];
			}
			for (int k = 0; k < NO_OF_REL_MODEL_DOFS_PER_OBS; ++k)
			{
				// we have to divide by the sum of weights here
				obs_op_jacobian_reduced_matrix[obs_index][k] = weights_vector[k]/sum_of_interpol_weights;
			}
			interpolated_model_dry[obs_index] = interpolated_model_dry[obs_index]/sum_of_interpol_weights;
		}
		// surface quantities (only surface pressure for now)
		else
		{
			sum_of_interpol_weights = 0;
			interpolated_model_dry[obs_index] = 0;
			// loop over all relevant model degrees of freedom
			for (int j = 0; j < NO_OF_REL_MODEL_DOFS_PER_OBS; ++j)
			{
				// we pick the lowest layer here (independant on wether we look at the temperature or the density)
				closest_vert_index = NO_OF_LAYERS - 1;
				// How is the suface pressure affected by the temperature in the lowest layer?
				if (j < NO_OF_REL_MODEL_DOFS_PER_OBS/2)
				{
					// radius does not matter here
					distance = calculate_distance_h(lat_used_obs[obs_index], lon_used_obs[obs_index], lat_model[rel_h_index_vector[obs_index_h][j]], lon_model[rel_h_index_vector[obs_index_h][j]], 1);
					// now we know which gridpoint is relevant to this observation
					relevant_model_dofs_matrix_dry[obs_index][j] = closest_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j];
					// 1/r-interpolation
					weights_vector[j] = 1/pow(distance + EPSILON, SP_INTERPOL_EXP)
					*R_D*background[relevant_model_dofs_matrix_dry[obs_index][j] + NO_OF_SCALARS_H]
					*exp(-(z_used_obs[obs_index] - z_model[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j]])/SCALE_HEIGHT);
					sum_of_interpol_weights += 1/pow(distance + EPSILON, SP_INTERPOL_EXP);
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
				// How is the suface pressure affected by the density in the lowest layer?
				else
				{
					// as a new interpolation will be conducted now, the sum_of_interpol_weights variable has to be reset to zero
					if (j == NO_OF_REL_MODEL_DOFS_PER_OBS/2)
					{
						sum_of_interpol_weights = 0;
					}
					// radius does not matter here
					distance = calculate_distance_h(lat_used_obs[obs_index], lon_used_obs[obs_index], lat_model[rel_h_index_vector[obs_index_h][j - NO_OF_REL_MODEL_DOFS_PER_OBS/2]], lon_model[rel_h_index_vector[obs_index_h][j - NO_OF_REL_MODEL_DOFS_PER_OBS/2]], 1);
					// now we know which gridpoint is relevant to this observation
					relevant_model_dofs_matrix_dry[obs_index][j] = (closest_vert_index + 1)*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j - NO_OF_REL_MODEL_DOFS_PER_OBS/2];
					// 1/r-interpolation
					weights_vector[j] = 1/pow(distance + EPSILON, SP_INTERPOL_EXP)
					*R_D*background[relevant_model_dofs_matrix_dry[obs_index][j] - NO_OF_SCALARS_H]
					*exp(-(z_used_obs[NO_OF_CHOSEN_OBSERVATIONS_MOIST + obs_index]
					- z_model[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + rel_h_index_vector[obs_index_h][j - NO_OF_REL_MODEL_DOFS_PER_OBS/2]])/SCALE_HEIGHT);
					// interpolation to the surface pressure
					interpolated_model_dry[obs_index] += weights_vector[j]*background[relevant_model_dofs_matrix_dry[obs_index][j]];
					sum_of_interpol_weights += 1/pow(distance + EPSILON, SP_INTERPOL_EXP);
				}
			}
			// the result
			// the interpolation to the surface pressure
			interpolated_model_dry[obs_index] = interpolated_model_dry[obs_index]/sum_of_interpol_weights;
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




