/*
This source file is part of ndvar, which is released under the MIT license.
Github repository: https://github.com/OpenNWP/ndvar
*/
/*
This file coordinates the data assimilation process.
*/

#include <stdlib.h>
#include "enum.h"
#include "ndvar.h"
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
	if (NO_OF_SURFACE_FIELDS_OBS < 1)
	{
		printf("NO_OF_SURFACE_FIELDS_OBS must be larger than zero.\n");
		printf("Aborting.\n");
		exit(2);
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
    char *ndvar_root_dir = malloc((len + 1)*sizeof(char));
    strcpy(ndvar_root_dir, argv[8]);
	int NO_OF_ORO_LAYERS;
    NO_OF_ORO_LAYERS = strtod(argv[9], NULL);
	int TOA = strtod(argv[10], NULL);
	int OI_SOLUTION_METHOD = strtod(argv[11], NULL);
	printf("background state file: %s\n", BACKGROUND_STATE_FILE);
    
    // Allocating memory for the grid properties.
    double *direction = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *latitude_scalar = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *longitude_scalar = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *latitude_vector = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *longitude_vector = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *z_scalar = malloc(NO_OF_SCALARS*sizeof(double));
    double *z_vector = malloc(NO_OF_VECTORS*sizeof(double));
    double *gravity_potential = malloc(NO_OF_VECTORS*sizeof(double));
    
    // Reading the grid properties.
    int ncid_grid, retval;
    int GEO_PROP_FILE_LENGTH = 200;
    char *GEO_PROP_FILE_PRE = malloc((GEO_PROP_FILE_LENGTH + 1)*sizeof(char));
    sprintf(GEO_PROP_FILE_PRE, "%s/grid_generator/grids/B%dL%dT%d_O%d_OL%d_SCVT.nc", model_home_dir, RES_ID, NO_OF_LAYERS, TOA, ORO_ID, NO_OF_ORO_LAYERS);
    GEO_PROP_FILE_LENGTH = strlen(GEO_PROP_FILE_PRE);
    free(GEO_PROP_FILE_PRE);
    char *GEO_PROP_FILE = malloc((GEO_PROP_FILE_LENGTH + 1)*sizeof(char));
    sprintf(GEO_PROP_FILE, "%s/grid_generator/grids/B%dL%dT%d_O%d_OL%d_SCVT.nc", model_home_dir, RES_ID, NO_OF_LAYERS, TOA, ORO_ID, NO_OF_ORO_LAYERS);
	printf("grid file: %s\n", GEO_PROP_FILE);
	printf("reading grid file ...\n");
    if ((retval = nc_open(GEO_PROP_FILE, NC_NOWRITE, &ncid_grid)))
        NCERR(retval);
    int direction_id, latitude_scalar_id, longitude_scalar_id, latitude_vector_id, longitude_vector_id, z_scalar_id, z_vector_id, gravity_potential_id, stretching_parameter_grid_id;
    double stretching_parameter_grid;
    if ((retval = nc_inq_varid(ncid_grid, "direction", &direction_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "latitude_scalar", &latitude_scalar_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "longitude_scalar", &longitude_scalar_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "latitude_vector", &latitude_vector_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "longitude_vector", &longitude_vector_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "z_scalar", &z_scalar_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "z_vector", &z_vector_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "gravity_potential", &gravity_potential_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "stretching_parameter", &stretching_parameter_grid_id)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, stretching_parameter_grid_id, &stretching_parameter_grid)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, direction_id, &direction[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, latitude_scalar_id, &latitude_scalar[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, longitude_scalar_id, &longitude_scalar[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, latitude_vector_id, &latitude_vector[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, longitude_vector_id, &longitude_vector[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, z_scalar_id, &z_scalar[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, z_vector_id, &z_vector[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, gravity_potential_id, &gravity_potential[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid_grid)))
        NCERR(retval);
	printf("Grid file read.\n");
	
    int OUTPUT_FILE_LENGTH = 100;
    char *OUTPUT_FILE_PRE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE_PRE, "%s/input/%s%s%s%s_nwp_B%dL%dT%d_O%d_OL%d_SCVT.nc", model_home_dir, year_string, month_string, day_string, hour_string, RES_ID, NO_OF_LAYERS, TOA, ORO_ID, NO_OF_ORO_LAYERS);
    OUTPUT_FILE_LENGTH = strlen(OUTPUT_FILE_PRE);
    free(OUTPUT_FILE_PRE);
    char *OUTPUT_FILE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE, "%s/input/%s%s%s%s_nwp_B%dL%dT%d_O%d_OL%d_SCVT.nc", model_home_dir, year_string, month_string, day_string, hour_string, RES_ID, NO_OF_LAYERS, TOA, ORO_ID, NO_OF_ORO_LAYERS);
    
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
	double *background = malloc(NO_OF_MODEL_DOFS_DRY*sizeof(double));
    #pragma omp parallel for
	for (int i = 0; i < NO_OF_MODEL_DOFS_DRY; ++i)
	{
		if (i < NO_OF_SCALARS)
		{
			background[i] = temperature_gas_background[i];
		}
		else
		{
			background[i] = density_dry_background[i - NO_OF_SCALARS_H];
		}
	}

	// Comparing the stretching parameters of the grid and the background state.
	if (stretching_parameter_grid != stretching_parameter_state)
	{
		printf("stretching_parameters of grid and background state do not conform.\n");
		printf("Aborting.\n");
		exit(1);
	}
    free(GEO_PROP_FILE);
	double stretching_parameter = stretching_parameter_grid;	
	
	// Allocating the memory for the observations.
	double *latitude_vector_obs = malloc(NO_OF_CHOSEN_OBSERVATIONS_DRY*sizeof(double));
	double *longitude_vector_obs = malloc(NO_OF_CHOSEN_OBSERVATIONS_DRY*sizeof(double));
	double *vert_vector = malloc(NO_OF_CHOSEN_OBSERVATIONS_DRY*sizeof(double));
	double *observations_vector = malloc(NO_OF_CHOSEN_OBSERVATIONS_DRY*sizeof(double));
	int *type_vector = malloc(NO_OF_CHOSEN_OBSERVATIONS_DRY*sizeof(int));
    
    int OBSERVATIONS_FILE_LENGTH = 100;
    char *OBSERVATIONS_FILE_PRE = malloc((OBSERVATIONS_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OBSERVATIONS_FILE_PRE, "%s/input/obs_%s%s%s%s.nc", ndvar_root_dir, year_string, month_string, day_string, hour_string);
    OBSERVATIONS_FILE_LENGTH = strlen(OBSERVATIONS_FILE_PRE);
    free(OBSERVATIONS_FILE_PRE);
    char *OBSERVATIONS_FILE = malloc((OBSERVATIONS_FILE_LENGTH + 1)*sizeof(char));
	sprintf(OBSERVATIONS_FILE, "%s/input/obs_%s%s%s%s.nc", ndvar_root_dir, year_string, month_string, day_string, hour_string);
	printf("observations file: %s\n", OBSERVATIONS_FILE);
    
    // Reading the observations.
	printf("reading observations ...\n");
    if ((retval = nc_open(OBSERVATIONS_FILE, NC_NOWRITE, &ncid)))
        NCERR(retval);
    int latitude_obs_id, longitude_obs_id, vert_id, obervations_id, type_id;
    // Defining the variables.
    if ((retval = nc_inq_varid(ncid, "latitude_vector", &latitude_obs_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "longitude_vector", &longitude_obs_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "vert_vector", &vert_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "observations_vector", &obervations_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "type_vector", &type_id)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, latitude_obs_id, &latitude_vector_obs[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, longitude_obs_id, &longitude_vector_obs[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, vert_id, &vert_vector[0])))
        NCERR(retval);  
    if ((retval = nc_get_var_double(ncid, obervations_id, &observations_vector[0])))
        NCERR(retval);    
    if ((retval = nc_get_var_int(ncid, type_id, &type_vector[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid)))
    	NCERR(retval);
	printf("Observations read.\n");
	
	// Begin of the actual assimilation.
    
    // setting up the measurement error covariance matrix
    double *obs_error_cov = malloc(sizeof(double[NO_OF_CHOSEN_OBSERVATIONS_DRY]));
    double temperature_error_obs = 1;
    double pressure_error_obs = 100;
    for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS_DRY; ++i)
    {
    	if (i < NO_OF_CHOSEN_OBSERVATIONS_DRY - NO_OF_SURFACE_FIELDS_OBS*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS)
    	{
			obs_error_cov[i] = pow(temperature_error_obs, 2);
		}
		else
		{
			obs_error_cov[i] = pow(pressure_error_obs, 2);
		}
    }
    
    // setting up the background error covariance matrix (only the diagonal)
    double *bg_error_cov = malloc(NO_OF_MODEL_DOFS_DRY*sizeof(double));
    double temperature_error_model = 1;
    double pressure_error_model = 50;
    #pragma omp parallel for
    for (int i = 0; i < NO_OF_MODEL_DOFS_DRY; ++i)
    {
    	if (i < NO_OF_SCALARS)
    	{
			bg_error_cov[i] = pow(temperature_error_model, 2);
    	}
    	else
    	{
    		// density = p/(R_D*T) (Gauss'ian error propagation)
    		bg_error_cov[i] = pow(pressure_error_model/(R_D*background[i - NO_OF_SCALARS_H]), 2) + pow(background[i]/background[i - NO_OF_SCALARS_H]*temperature_error_model, 2);
    	}
    }
    
	// setting up the observations operator
	double *interpolated_model = malloc(NO_OF_CHOSEN_OBSERVATIONS_DRY*sizeof(double));
	double (*obs_op_jacobian_reduced_matrix_dry)[NO_OF_REL_MODEL_DOFS_PER_OBS] = malloc(sizeof(double[NO_OF_CHOSEN_OBSERVATIONS_DRY][NO_OF_REL_MODEL_DOFS_PER_OBS]));
	int (*relevant_model_dofs_matrix_dry)[NO_OF_REL_MODEL_DOFS_PER_OBS] = malloc(sizeof(int[NO_OF_CHOSEN_OBSERVATIONS_DRY][NO_OF_REL_MODEL_DOFS_PER_OBS]));
	obs_op_setup(interpolated_model, obs_op_jacobian_reduced_matrix_dry, relevant_model_dofs_matrix_dry, latitude_vector_obs, longitude_vector_obs, vert_vector, latitude_scalar, longitude_scalar, z_scalar, background);
	
	double *observations_vector_dry = malloc(NO_OF_CHOSEN_OBSERVATIONS_DRY*sizeof(double));
	// setting up the dry observations vector
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS_DRY; ++i)
	{
		observations_vector_dry[i] = observations_vector[i];
	}
	
	// now, all the constituents of the gain matrix are known
	double *model_vector = malloc((NO_OF_SCALARS + NO_OF_SCALARS_H)*sizeof(double));
	oi(obs_error_cov, obs_op_jacobian_reduced_matrix_dry, relevant_model_dofs_matrix_dry, bg_error_cov, interpolated_model, background, observations_vector_dry, model_vector, OI_SOLUTION_METHOD);
	
	// data assimilation is finished at this point
	// freeing the memory
	free(obs_error_cov);
	free(bg_error_cov);
	free(interpolated_model);
	free(background);
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
    
    // the temperature comes first in the model_vector
    #pragma omp parallel for
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	temperature[i] = model_vector[i];
    }
    
    // density is determined out of the hydrostatic equation
    int layer_index, h_index;
    double b, c;
    for (int i = NO_OF_SCALARS - 1; i >= 0; --i)
    {
    	layer_index = i/NO_OF_SCALARS_H;
    	h_index = i - layer_index*NO_OF_SCALARS_H;
    	// at the lowest layer the density is part of the model_vector
    	if (layer_index == NO_OF_LAYERS - 1)
    	{
        	density_dry[i] = model_vector[NO_OF_SCALARS + h_index];
        	exner[i] = pow((model_vector[NO_OF_SCALARS + h_index]*R_D*temperature[i])/P_0, R_D/C_D_P);
        }
        else
        {
			// solving a quadratic equation for the Exner pressure
			b = -0.5*exner[i + NO_OF_SCALARS_H]/temperature[i + NO_OF_SCALARS_H]
			*(temperature[i] - temperature[i + NO_OF_SCALARS_H]
			+ 2/C_D_P*(gravity_potential[i] - gravity_potential[i + NO_OF_SCALARS_H]));
			c = pow(exner[i + NO_OF_SCALARS_H], 2)*temperature[i]/temperature[i + NO_OF_SCALARS_H];
			exner[i] = b + pow((pow(b, 2) + c), 0.5);
        	density_dry[i] = P_0*pow(exner[i], C_D_P/R_D)/(R_D*temperature[i]);
        }
    }
    
    free(exner);
    free(model_vector);
    
    // Wind is set equal to the background wind for now. Later it will be derived from the balance equation.
    #pragma omp parallel for
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
    	wind[i] = wind_background[i];
    }
    
    // end of the assimilation of the dry state
    
    // separate moisture assimilation
    double *observations_vector_moist = malloc(NO_OF_CHOSEN_OBSERVATIONS_MOIST*sizeof(double));
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS_MOIST; ++i)
	{
		observations_vector_moist[i] = observations_vector[i];
	}
    
    // if moisture is not assimilated, it is set equally to the background state
    if (NO_OF_FIELDS_PER_LAYER_OBS == 1)
    {
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			water_vapour_density[i] = water_vapour_density_background[i];
			liquid_water_density[i] = liquid_water_density_background[i];
			solid_water_density[i] = solid_water_density_background[i];
		}
	}
	// this is the case where moisture is assimilated
    if (NO_OF_FIELDS_PER_LAYER_OBS == 2)
    {
	    double *background_moist = malloc(NO_OF_SCALARS*sizeof(double));
	    // the data assimilation is being calculated with the specific humidity for pragmatic reasons
    	for (int i = 0; i < NO_OF_SCALARS; ++i)
    	{
    		background_moist[i] = water_vapour_density_background[i]/(density_dry_background[i] + water_vapour_density_background[i]);
    	}
    	
		// setting up the measurement error covariance matrix
		double *obs_error_cov = malloc(sizeof(double[NO_OF_CHOSEN_OBSERVATIONS_MOIST]));
		double abs_moisture_error_obs = 0.01;
		for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS_MOIST; ++i)
		{
			obs_error_cov[i] = pow(abs_moisture_error_obs, 2);
		}
		
		// setting up the background error covariance matrix (only the diagonal)
		double *bg_error_cov = malloc(NO_OF_MODEL_DOFS_MOIST*sizeof(double));
		double abs_moisture_error_model = 0.005;
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_MODEL_DOFS_MOIST; ++i)
		{
			bg_error_cov[i] = pow(abs_moisture_error_model, 2);
		}
		
		// setting up the observations operator
		double *interpolated_model = malloc(NO_OF_CHOSEN_OBSERVATIONS_MOIST*sizeof(double));
		double (*obs_op_jacobian_reduced_matrix_moist)[NO_OF_REL_MODEL_DOFS_PER_OBS] = malloc(sizeof(double[NO_OF_CHOSEN_OBSERVATIONS_MOIST][NO_OF_REL_MODEL_DOFS_PER_OBS]));
		int (*relevant_model_dofs_matrix_moist)[NO_OF_REL_MODEL_DOFS_PER_OBS] = malloc(sizeof(int[NO_OF_CHOSEN_OBSERVATIONS_MOIST][NO_OF_REL_MODEL_DOFS_PER_OBS]));
		// setting up the moist obs operator using the dry obs operator
		for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS_MOIST; ++i)
		{
			interpolated_model[i] = 0;
			for (int j = 0; j < NO_OF_REL_MODEL_DOFS_PER_OBS; ++j)
			{
				obs_op_jacobian_reduced_matrix_moist[i][j] = obs_op_jacobian_reduced_matrix_dry[i][j];
				relevant_model_dofs_matrix_moist[i][j] = relevant_model_dofs_matrix_dry[i][j];
				interpolated_model[i] += obs_op_jacobian_reduced_matrix_moist[i][j] + background_moist[relevant_model_dofs_matrix_dry[i][j]];
			}
		}
		
		// now, all the constituents of the gain matrix are known
		double *model_vector = malloc(NO_OF_SCALARS*sizeof(double));
		oi(obs_error_cov, obs_op_jacobian_reduced_matrix_moist, relevant_model_dofs_matrix_moist, bg_error_cov, interpolated_model, background_moist, observations_vector_moist, model_vector, OI_SOLUTION_METHOD);
		
		free(obs_op_jacobian_reduced_matrix_moist);
		free(background_moist);
		free(relevant_model_dofs_matrix_moist);
		free(obs_error_cov);
		free(bg_error_cov);
		free(interpolated_model);
		
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_SCALARS; ++i)
		{
			water_vapour_density[i] = model_vector[i]/(1 - model_vector[i])*density_dry[i];
		}
		
		free(model_vector);
	}
	free(relevant_model_dofs_matrix_dry);
	free(obs_op_jacobian_reduced_matrix_dry);
	free(observations_vector_moist);
	
	// individual condensate temperatures are for higher resolutions, not yet implemented
    #pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		liquid_water_temp[i] = temperature[i];
		solid_water_temp[i] = temperature[i];
    }
    
    // precipitation is not assimilated for now
    
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

	free(gravity_potential);
	free(BACKGROUND_STATE_FILE);
	free(ndvar_root_dir);
	free(OBSERVATIONS_FILE);
    free(model_home_dir);
	free(temperature);
	free(density_dry);
	free(wind);
	free(water_vapour_density);
	free(liquid_water_density);
	free(solid_water_density);
	free(liquid_water_temp);
	free(solid_water_temp);
	free(year_string);
	free(month_string);
	free(day_string);
	free(hour_string);
    free(wind_background);
    free(temperature_gas_background);
    free(density_dry_background);
    free(water_vapour_density_background);
    free(liquid_water_density_background);
    free(solid_water_density_background);
    free(liquid_water_temperature_background);
    free(solid_water_temperature_background);
    free(z_scalar);
	free(latitude_vector_obs);
	free(longitude_vector_obs);
	free(vert_vector);
	free(type_vector);
    free(z_vector);
    free(latitude_scalar);
    free(longitude_scalar);
    free(direction);
    free(latitude_vector);
    free(longitude_vector);
    free(OUTPUT_FILE);
    return 0;
}

int obs_op_setup(double interpolated_model[], double obs_op_jacobian_reduced_matrix[][NO_OF_REL_MODEL_DOFS_PER_OBS], int relevant_model_dofs_matrix_dry[][NO_OF_REL_MODEL_DOFS_PER_OBS], double lat_used_obs[], double lon_used_obs[], double z_used_obs[], double lat_model[], double lon_model[], double z_model[], double background[])
{
	/*
	This functions calculates the observations operator.
	It is the background state, interpolated to the observations
	+ the derivative of this function, which will be used to calculate
	the perturbation induced by the observations.
	*/
	double R_D = specific_gas_constants_lookup(0);
	// finding the NO_OF_REL_MODEL_DOFS_PER_OBS closest grid points (horizontally) for each observation
	int (*rel_h_index_vector)[NO_OF_REL_MODEL_DOFS_PER_OBS/2] = malloc(sizeof(int[NO_OF_CHOSEN_OBSERVATIONS][NO_OF_REL_MODEL_DOFS_PER_OBS/2])); // the vector containing the relevant horizontal model indices for each observation
	double *dist_vector = malloc(NO_OF_SCALARS_H*sizeof(double)); // the vector containing the horizontal distances between the observation at hand and each horizontal model gridpoint
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS; ++i)
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
	
	// the vector containing the vertical distance between the observation at hand and the model gridpoints
	double vert_distance_vector[NO_OF_LAYERS];
	// the vector containing preliminary interpolation weights
	double weights_vector[NO_OF_REL_MODEL_DOFS_PER_OBS];
	// the closest vertical index
	int closest_vert_index, other_vert_index;
	double sum_of_interpol_weights, distance, closest_vert_weight, other_vert_weight;
	// finally setting up the reduced observations operator
	#pragma omp parallel for private(vert_distance_vector, weights_vector, closest_vert_index, other_vert_index, sum_of_interpol_weights, distance, closest_vert_weight, other_vert_weight)
	for (int obs_index = 0; obs_index < NO_OF_CHOSEN_OBSERVATIONS; ++obs_index)
	{
		// free atmosphere quantities (temperature, specific humidity)
		if (obs_index < NO_OF_CHOSEN_OBSERVATIONS - NO_OF_SURFACE_FIELDS_OBS*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS)
		{
			sum_of_interpol_weights = 0;
			interpolated_model[obs_index] = 0;
			// loop over all relevant horizontal model gridpoints
			for (int j = 0; j < NO_OF_REL_MODEL_DOFS_PER_OBS/2; ++j)
			{
				// finding out which layer is the closest to the observation
				for (int k = 0; k < NO_OF_LAYERS; ++k)
				{
					vert_distance_vector[k] = fabs(z_model[k*NO_OF_SCALARS_H + rel_h_index_vector[obs_index][j]] - z_used_obs[obs_index]);
				}
				closest_vert_index = find_min_index(vert_distance_vector, NO_OF_LAYERS);
				// vertical interpolation
				// first setting for the other vertical index
				other_vert_index = closest_vert_index + 1;
				// if the the closest model point is below the observation, the next higher point is taken into account for the interpolation
				if (z_model[closest_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index][j]] < z_used_obs[obs_index])
				{
					other_vert_index = closest_vert_index - 1;
				}
				// if the observation is below the lowest layer of the model
				if (other_vert_index == NO_OF_LAYERS)
				{
					other_vert_index = NO_OF_LAYERS - 2;
					closest_vert_weight = 1 - (z_used_obs[obs_index] - z_model[closest_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index][j]])
					/(z_model[other_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index][j]] - z_model[closest_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index][j]]);
					other_vert_weight = (z_used_obs[obs_index] - z_model[closest_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index][j]])
					/(z_model[other_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index][j]] - z_model[closest_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index][j]]);
				}
				else
				{
					closest_vert_weight = fabs(z_model[other_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index][j]] - z_used_obs[obs_index])/fabs(z_model[closest_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index][j]] - z_model[other_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index][j]]);
					other_vert_weight = 1 - closest_vert_weight;
				}
				// now we know which gridpoint is relevant to this observation
				// the closest vertical point
				relevant_model_dofs_matrix_dry[obs_index][j] = closest_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index][j];
				// the second closest vertical point
				relevant_model_dofs_matrix_dry[obs_index][j + NO_OF_REL_MODEL_DOFS_PER_OBS/2] = other_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index][j];
				// radius does not matter here
				distance = calculate_distance_h(lat_used_obs[obs_index], lon_used_obs[obs_index], lat_model[rel_h_index_vector[obs_index][j]], lon_model[rel_h_index_vector[obs_index][j]], 1);
				// 1/r-interpolation
				weights_vector[j] = closest_vert_weight/pow(distance + EPSILON, T_INTERPOL_EXP);
				weights_vector[j + NO_OF_REL_MODEL_DOFS_PER_OBS/2] = other_vert_weight/pow(distance + EPSILON, T_INTERPOL_EXP);
				interpolated_model[obs_index] += weights_vector[j]*background[relevant_model_dofs_matrix_dry[obs_index][j]];
				interpolated_model[obs_index] += weights_vector[j + NO_OF_REL_MODEL_DOFS_PER_OBS/2]*background[relevant_model_dofs_matrix_dry[obs_index][j + NO_OF_REL_MODEL_DOFS_PER_OBS/2]];
				sum_of_interpol_weights += weights_vector[j];
				sum_of_interpol_weights += weights_vector[j + NO_OF_REL_MODEL_DOFS_PER_OBS/2];
				if (j == NO_OF_REL_MODEL_DOFS_PER_OBS/2 - 1)
				{
					for (int k = 0; k < NO_OF_REL_MODEL_DOFS_PER_OBS; ++k)
					{
						// we have to divide by the sum of weights here
						obs_op_jacobian_reduced_matrix[obs_index][k] = weights_vector[k]/sum_of_interpol_weights;
					}
					interpolated_model[obs_index] = interpolated_model[obs_index]/sum_of_interpol_weights;
				}
			}
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
					distance = calculate_distance_h(lat_used_obs[obs_index], lon_used_obs[obs_index], lat_model[rel_h_index_vector[obs_index][j]], lon_model[rel_h_index_vector[obs_index][j]], 1);
					// now we know which gridpoint is relevant to this observation
					relevant_model_dofs_matrix_dry[obs_index][j] = closest_vert_index*NO_OF_SCALARS_H + rel_h_index_vector[obs_index][j];
					// 1/r-interpolation
					weights_vector[j] = 1/pow(distance + EPSILON, SP_INTERPOL_EXP)
					*R_D*background[relevant_model_dofs_matrix_dry[obs_index][j] + NO_OF_SCALARS_H]
					*exp(-(z_used_obs[obs_index] - z_model[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + rel_h_index_vector[obs_index][j]])/SCALE_HEIGHT);
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
					distance = calculate_distance_h(lat_used_obs[obs_index], lon_used_obs[obs_index], lat_model[rel_h_index_vector[obs_index][j - NO_OF_REL_MODEL_DOFS_PER_OBS/2]], lon_model[rel_h_index_vector[obs_index][j - NO_OF_REL_MODEL_DOFS_PER_OBS/2]], 1);
					// now we know which gridpoint is relevant to this observation
					relevant_model_dofs_matrix_dry[obs_index][j] = (closest_vert_index + 1)*NO_OF_SCALARS_H + rel_h_index_vector[obs_index][j - NO_OF_REL_MODEL_DOFS_PER_OBS/2];
					// 1/r-interpolation
					weights_vector[j] = 1/pow(distance + EPSILON, SP_INTERPOL_EXP)
					*R_D*background[relevant_model_dofs_matrix_dry[obs_index][j] - NO_OF_SCALARS_H]
					*exp(-(z_used_obs[obs_index] - z_model[(NO_OF_LAYERS - 1)*NO_OF_SCALARS_H + rel_h_index_vector[obs_index][j - NO_OF_REL_MODEL_DOFS_PER_OBS/2]])/SCALE_HEIGHT);
					// interpolation to the surface pressure
					interpolated_model[obs_index] += weights_vector[j]*background[relevant_model_dofs_matrix_dry[obs_index][j]];
					sum_of_interpol_weights += 1/pow(distance + EPSILON, SP_INTERPOL_EXP);
					// the result
					if (j == NO_OF_REL_MODEL_DOFS_PER_OBS - 1)
					{
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
			}
		}
	}
	free(rel_h_index_vector);
	// returning 0 indicating success
	return 0;
}




