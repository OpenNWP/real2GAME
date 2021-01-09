/*
This source file is part of ndvar, which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/ndvar
*/

#include <stdlib.h>
#include "enum.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <netcdf.h>
#include "geos95.h"
#include "atmostracers.h"
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}
#define N_A (6.0221409e23)
#define K_B (1.380649e-23)
#define M_D 0.028964420
#define R (N_A*K_B)
#define R_D (R/M_D)
#define P_0 100000.0
#define OMEGA (7.292115e-5)
#define C_D_P 1005.0
#define EPSILON 1e-4

const int NO_OF_LEVELS_OBS = 3;
const int NO_OF_FIELDS_PER_LAYER_OBS = 1;
const int NO_OF_SURFACE_FIELDS_OBS = 0;
const int NO_OF_POINTS_PER_LAYER_OBS = 10;
const int NO_OF_OBSERVATIONS = (NO_OF_LEVELS_OBS*NO_OF_FIELDS_PER_LAYER_OBS + NO_OF_SURFACE_FIELDS_OBS)*NO_OF_POINTS_PER_LAYER_OBS;
// the number of observations that will actually be used (must be <= NO_OF_OBSERVATIONS)

int interpolate_bg_to_obs(double [], double [], double [], double [], double [], double [], double [], double [], double [][NO_OF_SCALARS]);

int main(int argc, char *argv[])
{	
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
	int TOA;
	TOA = strtod(argv[10], NULL);
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
    int GEO_PROP_FILE_LENGTH = 100;
    char *GEO_PROP_FILE_PRE = malloc((GEO_PROP_FILE_LENGTH + 1)*sizeof(char));
    sprintf(GEO_PROP_FILE_PRE, "%s/grids/B%dL%dT%d_O%d_OL%d_SCVT.nc", model_home_dir, RES_ID, NO_OF_LAYERS, TOA, ORO_ID, NO_OF_ORO_LAYERS);
    GEO_PROP_FILE_LENGTH = strlen(GEO_PROP_FILE_PRE);
    free(GEO_PROP_FILE_PRE);
    char *GEO_PROP_FILE = malloc((GEO_PROP_FILE_LENGTH + 1)*sizeof(char));
    sprintf(GEO_PROP_FILE, "%s/grids/B%dL%dT%d_O%d_OL%d_SCVT.nc", model_home_dir, RES_ID, NO_OF_LAYERS, TOA, ORO_ID, NO_OF_ORO_LAYERS);
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
	double *latitude_vector_obs = malloc(NO_OF_OBSERVATIONS*sizeof(double));
	double *longitude_vector_obs = malloc(NO_OF_OBSERVATIONS*sizeof(double));
	double *vert_vector = malloc(NO_OF_OBSERVATIONS*sizeof(double));
	double *observations_vector = malloc(NO_OF_OBSERVATIONS*sizeof(double));
	int *type_vector = malloc(NO_OF_OBSERVATIONS*sizeof(int));
    
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
    double (*obs_error_cov)[NO_OF_OBSERVATIONS] = malloc(sizeof(double[NO_OF_OBSERVATIONS][NO_OF_OBSERVATIONS]));
    double temperature_error_obs = 0.2;
    for (int i = 0; i < NO_OF_OBSERVATIONS; ++i)
    {
    	for (int j = 0; j < NO_OF_OBSERVATIONS; ++j)
    	{
    		if (i == j)
    		{
    			obs_error_cov[i][j] = pow(temperature_error_obs, 2);
    		}
    		else
    		{
    			obs_error_cov[i][j] = 0;
    		}
    	}
    }
    
    // setting up the background error covariance matrix (only the diagonal)
    double *bg_error_cov = malloc(NO_OF_SCALARS*sizeof(double));
    double temperature_error_model = 0.2;
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
		bg_error_cov[i] = pow(temperature_error_model, 2);
    }
    
    // setting up the observations_operator
    double (*obs_op)[NO_OF_SCALARS] = malloc(sizeof(double[NO_OF_OBSERVATIONS][NO_OF_SCALARS]));
	
	// this vector will contain the values expected for the observations, assuming the background state
	double *interpolated_model = malloc(NO_OF_OBSERVATIONS*sizeof(double));
	// interpolating the background state to the observations
	interpolate_bg_to_obs(interpolated_model, latitude_vector_obs, longitude_vector_obs, vert_vector, latitude_scalar, longitude_scalar, z_scalar, temperature_gas_background, obs_op);
	
	// now, all the constituents of the gain matrix are known
	// short notation: b: background error covariance, h: observations operator; r: observations error covariance
	double (*b_ht)[NO_OF_OBSERVATIONS] = malloc(sizeof(double[NO_OF_SCALARS][NO_OF_OBSERVATIONS]));
	
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		for (int j = 0; j < NO_OF_OBSERVATIONS; ++j)
		{
			b_ht[i][j] = bg_error_cov[i]*obs_op[j][i];
		}
	}
	double (*h_b_ht_plus_r)[NO_OF_OBSERVATIONS] = malloc(sizeof(double[NO_OF_OBSERVATIONS][NO_OF_OBSERVATIONS]));
	
	for (int i = 0; i < NO_OF_OBSERVATIONS; ++i)
	{
		for (int j = 0; j < NO_OF_OBSERVATIONS; ++j)
		{
			h_b_ht_plus_r[i][j] = 0;
			for (int k = 0; k < NO_OF_SCALARS; ++k)
			{
				h_b_ht_plus_r[i][j] = h_b_ht_plus_r[i][j] + obs_op[i][k]*b_ht[k][j];
			}
			h_b_ht_plus_r[i][j] = h_b_ht_plus_r[i][j] + obs_error_cov[i][j];
		}
	}
	
	// h_b_ht_plus_r needs to be inversed in order to calculate the gain matrix
	// this is actually the main task of OI
	double (*h_b_ht_plus_r_inv)[NO_OF_OBSERVATIONS] = malloc(sizeof(double[NO_OF_OBSERVATIONS][NO_OF_OBSERVATIONS]));
	// firstly, the inverse is initialized with the unity matrix
	for (int i = 0; i < NO_OF_OBSERVATIONS; ++i)
	{
		for (int j = 0; j < NO_OF_OBSERVATIONS; ++j)
		{
			h_b_ht_plus_r_inv[i][j] = 0;
			if (i == j)
			{
				h_b_ht_plus_r_inv[i][j] = 1;
			}
		}
	}
	// we will start to modify h_b_ht_plus_r now (misuse of name)
	// Gaussian downwards
	double factor = 0;
	// starting with the first line, down to the second but last line
	for (int i = 0; i < NO_OF_OBSERVATIONS - 1; ++i)
	{
		// dividing the line by h_b_ht_plus_r[i][i]
		for (int j = 0; j < NO_OF_OBSERVATIONS; ++j)
		{
			h_b_ht_plus_r[i][j] = h_b_ht_plus_r[i][j]/h_b_ht_plus_r[i][i];
			h_b_ht_plus_r_inv[i][j] = h_b_ht_plus_r_inv[i][j]/h_b_ht_plus_r[i][i];
		}
		for (int j = i + 1; j < NO_OF_OBSERVATIONS; ++j)
		{
			factor = -h_b_ht_plus_r[j][i];
			for (int k = 0; k < NO_OF_OBSERVATIONS; ++k)
			{
				h_b_ht_plus_r[j][k] = h_b_ht_plus_r[j][k] + factor*h_b_ht_plus_r[i][k];
				h_b_ht_plus_r_inv[j][k] = h_b_ht_plus_r_inv[j][k] + factor*h_b_ht_plus_r_inv[i][k];
			}
		}
	}
	for (int j = 0; j < NO_OF_OBSERVATIONS; ++j)
	{
		h_b_ht_plus_r[NO_OF_OBSERVATIONS - 1][j] = h_b_ht_plus_r[NO_OF_OBSERVATIONS - 1][j]/h_b_ht_plus_r[NO_OF_OBSERVATIONS - 1][NO_OF_OBSERVATIONS - 1];
		h_b_ht_plus_r_inv[NO_OF_OBSERVATIONS - 1][j] = h_b_ht_plus_r_inv[NO_OF_OBSERVATIONS - 1][j]/h_b_ht_plus_r[NO_OF_OBSERVATIONS - 1][NO_OF_OBSERVATIONS - 1];
	}
	// Gaussian upwards
	// starting with the last line, then going up to the last but first
	for (int i = NO_OF_OBSERVATIONS - 1; i >= 1; --i)
	{
		for (int j = i - 1; j >= 0; --j)
		{
			factor = -h_b_ht_plus_r[j][i];
			for (int k = 0; k < NO_OF_OBSERVATIONS; ++k)
			{
				h_b_ht_plus_r[j][k] = h_b_ht_plus_r[j][k] + factor*h_b_ht_plus_r[i][k];
				h_b_ht_plus_r_inv[j][k] = h_b_ht_plus_r_inv[j][k] + factor*h_b_ht_plus_r_inv[i][k];
			}
		}
	}
	
	// now, the gain matrix can finally be put together
    double (*gain)[NO_OF_OBSERVATIONS] = malloc(sizeof(double[NO_OF_SCALARS][NO_OF_OBSERVATIONS]));
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{	
		for (int j = 0; j < NO_OF_OBSERVATIONS; ++j)
		{
			gain[i][j] = 0;
			for (int k = 0; k < NO_OF_OBSERVATIONS; ++k)
			{
				gain[i][j] = gain[i][j] + b_ht[i][k]*h_b_ht_plus_r_inv[k][j];
			}
		}
	}
	// now, the main job is already done
	
	free(h_b_ht_plus_r_inv);
	free(b_ht);
	free(h_b_ht_plus_r);
	free(obs_error_cov);
	free(bg_error_cov);
	free(obs_op);
	
	// this vector will contain the product of the model forecast error and the gain matrix
	double *prod_with_gain_matrix = malloc(NO_OF_SCALARS*sizeof(double));
	// multiplying (obs - (interpolated model)) by the gain matrix
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{	
		prod_with_gain_matrix[i] = 0;
		for (int j = 0; j < NO_OF_OBSERVATIONS; ++j)
		{
			prod_with_gain_matrix[i] += gain[i][j]*(observations_vector[j] - interpolated_model[j]);
		}
	}
	free(interpolated_model);
	
	double *model_vector = malloc(NO_OF_SCALARS*sizeof(double));
	
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		model_vector[i] = temperature_gas_background[i] + prod_with_gain_matrix[i];
	}
	
	free(gain);
	free(prod_with_gain_matrix);
	free(observations_vector);
	// End of the actual assimilation
    
    // These are the arrays for the result of the assimilation process.
    double *temperature = malloc(NO_OF_SCALARS*sizeof(double));
    double *density_dry = malloc(NO_OF_SCALARS*sizeof(double));
    double *wind = malloc(NO_OF_VECTORS*sizeof(double));
    double *water_vapour_density = malloc(NO_OF_SCALARS*sizeof(double));
    double *liquid_water_density = malloc(NO_OF_SCALARS*sizeof(double));
    double *solid_water_density = malloc(NO_OF_SCALARS*sizeof(double));
    double *liquid_water_temp = malloc(NO_OF_SCALARS*sizeof(double));
    double *solid_water_temp = malloc(NO_OF_SCALARS*sizeof(double));
    
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	temperature[i] = model_vector[i];
    	density_dry[i] = density_dry_background[i];
		water_vapour_density[i] = water_vapour_density_background[i];
		liquid_water_density[i] = liquid_water_density_background[i];
		solid_water_density[i] = solid_water_density_background[i];
		liquid_water_temp[i] = liquid_water_temperature_background[i];
		solid_water_temp[i] = solid_water_temperature_background[i];
    }
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
    	wind[i] = wind_background[i];
    }
    
    free(model_vector);
    
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
    if ((retval = nc_put_att_text(ncid, wind_background_id, "units", strlen("m/s"), "m/s")))
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


int interpolate_bg_to_obs(double interpolated_model[], double lat_used_obs[], double lon_used_obs[], double z_used_obs[], double lat_model[], double lon_model[], double z_model[], double background[], double obs_op [][NO_OF_SCALARS])
{
	// this functions calculates the expected values for the observations from the background state
	double weight, sum_of_weights, distance;
	double vert_distance_vector[NO_OF_LAYERS];
	int min_index;
	// loop over all observations to which we want to interpolate
	for (int i = 0; i < NO_OF_OBSERVATIONS; ++i)
	{
		// initializing the sum of interpolation weights as well as the interpolated value
		sum_of_weights = 0;
		interpolated_model[i] = 0;
		// loop over all horizontal model gridpoints
		for (int j = 0; j < NO_OF_SCALARS_H; ++j)
		{
			// finding out which layer is the closest to the observation
			for (int k = 0; k < NO_OF_LAYERS; ++k)
			{
				vert_distance_vector[k] = fabs(z_model[k*NO_OF_SCALARS_H + j] - z_used_obs[i]);
				// initializing the observations operator
				obs_op[i][k*NO_OF_SCALARS_H + j] = 0;
			}
			min_index = find_min_index(vert_distance_vector, NO_OF_LAYERS);
			// radius does not matter here
			distance = calculate_distance_h(lat_used_obs[i], lon_used_obs[i], lat_model[j], lon_model[j], 1);
			// 1/r-interpolation
			weight = 1/(distance + EPSILON);
			interpolated_model[i] += weight*background[min_index*NO_OF_SCALARS_H + j];
			sum_of_weights += weight;
		}
		interpolated_model[i] = interpolated_model[i]/sum_of_weights;
		
		// determining the derivative of the interpolation (the observations operator)
		for (int j = 0; j < NO_OF_SCALARS_H; ++j)
		{
			// finding out which layer is the closest to the observation
			for (int k = 0; k < NO_OF_LAYERS; ++k)
			{
				vert_distance_vector[k] = fabs(z_model[k*NO_OF_SCALARS_H + j] - z_used_obs[i]);
			}
			min_index = find_min_index(vert_distance_vector, NO_OF_LAYERS);
			// radius does not matter here
			distance = calculate_distance_h(lat_used_obs[i], lon_used_obs[i], lat_model[j], lon_model[j], 1);
			// 1/r-interpolation
			weight = 1/(distance + EPSILON);
			obs_op[i][min_index*NO_OF_SCALARS_H + j] = weight/sum_of_weights;
		}
	}
	return 0;
}










