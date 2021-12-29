/*
This source file is part of GAME-DA, which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME-DA
*/

#include <netcdf.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "eccodes.h"
#include "../../src/game-da.h"
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}
#define ERRCODE 3
#define ECCERR(e) {printf("Error: Eccodes failed with error code %d. See http://download.ecmwf.int/test-data/eccodes/html/group__errors.html for meaning of the error codes.\n", e); exit(ERRCODE);}

int main(int argc, char *argv[])
{
	// defining the levels of the model we want to use
	int levels_vector[NO_OF_LEVELS_OBS];
	levels_vector[0] = 24;
	levels_vector[1] = 37;
	levels_vector[2] = 50;
	levels_vector[3] = 63;
	levels_vector[4] = 77;
	levels_vector[5] = 90;
	
	// shell arguments
    char year_string[strlen(argv[1]) + 1];
    strcpy(year_string, argv[1]);
    char month_string[strlen(argv[2]) + 1];
    strcpy(month_string, argv[2]);
    char day_string[strlen(argv[3]) + 1];
    strcpy(day_string, argv[3]);
    char hour_string[strlen(argv[4]) + 1];
    strcpy(hour_string, argv[4]);
    char game_da_root_dir[strlen(argv[5]) + 1];
    strcpy(game_da_root_dir, argv[5]);
	
	// Properties of the input model's grid.
	double *latitudes_one_layer = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
	double *longitudes_one_layer = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
	double *latitudes_sst = malloc(NO_OF_SST_POINTS*sizeof(double));
	double *longitudes_sst = malloc(NO_OF_SST_POINTS*sizeof(double));
    
	int retval, err;
	codes_handle *handle = NULL;
	
	// latitudes of the grid
	char LAT_OBS_FILE_PRE[200];
    sprintf(LAT_OBS_FILE_PRE , "%s/input/icon_global_icosahedral_time-invariant_%s%s%s%s_CLAT.grib2", game_da_root_dir, year_string, month_string, day_string, hour_string);
	char LAT_OBS_FILE[strlen(LAT_OBS_FILE_PRE) + 1];
	strcpy(LAT_OBS_FILE, LAT_OBS_FILE_PRE);
    
	FILE *ECC_FILE;
	ECC_FILE = fopen(LAT_OBS_FILE, "r");
	handle = codes_handle_new_from_file(NULL, ECC_FILE, PRODUCT_GRIB, &err);
	if (err != 0)
		ECCERR(err);
	size_t NO_OF_POINTS_PER_LAYER_OBS_SIZE_T;
	NO_OF_POINTS_PER_LAYER_OBS_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_OBS;
	if ((retval = codes_get_double_array(handle, "values", &latitudes_one_layer[0], &NO_OF_POINTS_PER_LAYER_OBS_SIZE_T)))
		ECCERR(retval);
	codes_handle_delete(handle);
    fclose(ECC_FILE);
    
    for (int i = 0; i < NO_OF_POINTS_PER_LAYER_OBS; ++i)
    {
    	latitudes_one_layer[i] = 2*M_PI*latitudes_one_layer[i]/360;
    }
    
    // longitudes of the grid
	char LON_OBS_FILE_PRE[200];
    sprintf(LON_OBS_FILE_PRE , "%s/input/icon_global_icosahedral_time-invariant_%s%s%s%s_CLON.grib2", game_da_root_dir, year_string, month_string, day_string, hour_string);
	char LON_OBS_FILE[strlen(LON_OBS_FILE_PRE) + 1];
	strcpy(LON_OBS_FILE, LON_OBS_FILE_PRE);
    
	ECC_FILE = fopen(LON_OBS_FILE, "r");
	handle = codes_handle_new_from_file(NULL, ECC_FILE, PRODUCT_GRIB, &err);
	if (err != 0)
		ECCERR(err);
	NO_OF_POINTS_PER_LAYER_OBS_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_OBS;
	if ((retval = codes_get_double_array(handle, "values", &longitudes_one_layer[0], &NO_OF_POINTS_PER_LAYER_OBS_SIZE_T)))
		ECCERR(retval);
	codes_handle_delete(handle);
    fclose(ECC_FILE);
    
    for (int i = 0; i < NO_OF_POINTS_PER_LAYER_OBS; ++i)
    {
    	longitudes_one_layer[i] = 2*M_PI*longitudes_one_layer[i]/360;
    }
    
	// Allocating the memory for the final result.
	double *latitude_vector = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(double));
	double *longitude_vector = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(double));
	double *z_coords_amsl = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(double));
	double *observations_vector = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(double));
	
	double *temperature_one_layer = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
	double *u_wind_one_layer = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
	double *v_wind_one_layer = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
	double *spec_hum_one_layer = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
	
	// the indices of the chosen points for the scalar assimilation
	int *chosen_indices = malloc(NO_OF_CHOSEN_POINTS_PER_LAYER_OBS*sizeof(int));
	for (int i = 0; i < NO_OF_CHOSEN_POINTS_PER_LAYER_OBS; ++i)
	{
		chosen_indices[i] = NO_OF_POINTS_PER_LAYER_OBS/NO_OF_CHOSEN_POINTS_PER_LAYER_OBS*i;
	}
	// the indices of the chosen points for the wind assimilation
	int *chosen_indices_wind = malloc(NO_OF_CHOSEN_WIND_POINTS_PER_LAYER_OBS*sizeof(int));
	for (int i = 0; i < NO_OF_CHOSEN_WIND_POINTS_PER_LAYER_OBS; ++i)
	{
		chosen_indices_wind[i] = NO_OF_POINTS_PER_LAYER_OBS/NO_OF_CHOSEN_WIND_POINTS_PER_LAYER_OBS*i;
	}
	
	// reading the data from the free atmosphere
	double *z_height_amsl = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
	// loop over all relevant level in the free atmosphere
	for (int level_index = 0; level_index < NO_OF_LEVELS_OBS; ++level_index)
	{
		// vertical position of the current layer
		char Z_OBS_FILE_PRE[200];
		sprintf(Z_OBS_FILE_PRE , "%s/input/icon_global_icosahedral_time-invariant_%s%s%s%s_%d_HHL.grib2", game_da_root_dir, year_string, month_string, day_string, hour_string, levels_vector[level_index]);
		char Z_OBS_FILE[strlen(Z_OBS_FILE_PRE) + 1];
		strcpy(Z_OBS_FILE, Z_OBS_FILE_PRE);
		
		ECC_FILE = fopen(Z_OBS_FILE, "r");
		handle = codes_handle_new_from_file(NULL, ECC_FILE, PRODUCT_GRIB, &err);
		if (err != 0)
			ECCERR(err);
		NO_OF_POINTS_PER_LAYER_OBS_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_OBS;
		if ((retval = codes_get_double_array(handle, "values", &z_height_amsl[0], &NO_OF_POINTS_PER_LAYER_OBS_SIZE_T)))
			ECCERR(retval);
		codes_handle_delete(handle);
		fclose(ECC_FILE);
		
	   	// reading the temperature
		char TEMPERATURE_FILE_PRE[200];
		sprintf(TEMPERATURE_FILE_PRE , "%s/input/icon_global_icosahedral_model-level_%s%s%s%s_000_%d_T.grib2",
		game_da_root_dir, year_string, month_string, day_string, hour_string, levels_vector[level_index]);
		char TEMPERATURE_FILE[strlen(TEMPERATURE_FILE_PRE) + 1];
		strcpy(TEMPERATURE_FILE, TEMPERATURE_FILE_PRE);
		
		ECC_FILE = fopen(TEMPERATURE_FILE, "r");
		handle = codes_handle_new_from_file(NULL, ECC_FILE, PRODUCT_GRIB, &err);
		if (err != 0)
			ECCERR(err);
		NO_OF_POINTS_PER_LAYER_OBS_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_OBS;
		if ((retval = codes_get_double_array(handle, "values", &temperature_one_layer[0], &NO_OF_POINTS_PER_LAYER_OBS_SIZE_T)))
			ECCERR(retval);
		codes_handle_delete(handle);
		fclose(ECC_FILE);
		
	   	// reading the specific humidity
		char SPEC_HUM_FILE_PRE[200];
		sprintf(SPEC_HUM_FILE_PRE , "%s/input/icon_global_icosahedral_model-level_%s%s%s%s_000_%d_QV.grib2",
		game_da_root_dir, year_string, month_string, day_string, hour_string, levels_vector[level_index]);
		char SPEC_HUM_FILE[strlen(SPEC_HUM_FILE_PRE) + 1];
		strcpy(SPEC_HUM_FILE, SPEC_HUM_FILE_PRE);
		
		ECC_FILE = fopen(SPEC_HUM_FILE, "r");
		handle = codes_handle_new_from_file(NULL, ECC_FILE, PRODUCT_GRIB, &err);
		if (err != 0)
			ECCERR(err);
		NO_OF_POINTS_PER_LAYER_OBS_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_OBS;
		if ((retval = codes_get_double_array(handle, "values", &spec_hum_one_layer[0], &NO_OF_POINTS_PER_LAYER_OBS_SIZE_T)))
			ECCERR(retval);
		codes_handle_delete(handle);
		fclose(ECC_FILE);
		
	   	// reading the u wind
		char U_WIND_FILE_PRE[200];
		sprintf(U_WIND_FILE_PRE , "%s/input/icon_global_icosahedral_model-level_%s%s%s%s_000_%d_U.grib2",
		game_da_root_dir, year_string, month_string, day_string, hour_string, levels_vector[level_index]);
		char U_WIND_FILE[strlen(U_WIND_FILE_PRE) + 1];
		strcpy(U_WIND_FILE, U_WIND_FILE_PRE);
		
		ECC_FILE = fopen(U_WIND_FILE, "r");
		handle = codes_handle_new_from_file(NULL, ECC_FILE, PRODUCT_GRIB, &err);
		if (err != 0)
			ECCERR(err);
		NO_OF_POINTS_PER_LAYER_OBS_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_OBS/2;
		if ((retval = codes_get_double_array(handle, "values", &u_wind_one_layer[0], &NO_OF_POINTS_PER_LAYER_OBS_SIZE_T)))
			ECCERR(retval);
		codes_handle_delete(handle);
		fclose(ECC_FILE);
		
	   	// reading the v wind
		char V_WIND_FILE_PRE[200];
		sprintf(V_WIND_FILE_PRE , "%s/input/icon_global_icosahedral_model-level_%s%s%s%s_000_%d_V.grib2",
		game_da_root_dir, year_string, month_string, day_string, hour_string, levels_vector[level_index]);
		char V_WIND_FILE[strlen(V_WIND_FILE_PRE) + 1];
		strcpy(V_WIND_FILE, V_WIND_FILE_PRE);
		
		ECC_FILE = fopen(V_WIND_FILE, "r");
		handle = codes_handle_new_from_file(NULL, ECC_FILE, PRODUCT_GRIB, &err);
		if (err != 0)
			ECCERR(err);
		NO_OF_POINTS_PER_LAYER_OBS_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_OBS/2;
		if ((retval = codes_get_double_array(handle, "values", &v_wind_one_layer[0], &NO_OF_POINTS_PER_LAYER_OBS_SIZE_T)))
			ECCERR(retval);
		codes_handle_delete(handle);
		fclose(ECC_FILE);
		
		// formatting the observations
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_CHOSEN_POINTS_PER_LAYER_OBS; ++i)
		{
			// temperature
			latitude_vector[level_index*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS + i] = latitudes_one_layer[chosen_indices[i]];
			longitude_vector[level_index*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS + i] = longitudes_one_layer[chosen_indices[i]];
			z_coords_amsl[level_index*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS + i] = z_height_amsl[chosen_indices[i]];
			observations_vector[level_index*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS + i] = temperature_one_layer[chosen_indices[i]];
			
			// specific humidity
			latitude_vector[NO_OF_CHOSEN_OBSERVATIONS_MOIST + level_index*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS + i] = latitudes_one_layer[chosen_indices[i]];
			longitude_vector[NO_OF_CHOSEN_OBSERVATIONS_MOIST + level_index*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS + i] = longitudes_one_layer[chosen_indices[i]];
			z_coords_amsl[NO_OF_CHOSEN_OBSERVATIONS_MOIST + level_index*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS + i] = z_height_amsl[chosen_indices[i]];
			observations_vector[NO_OF_CHOSEN_OBSERVATIONS_MOIST + level_index*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS + i] = spec_hum_one_layer[chosen_indices[i]];
			
			// u wind
			latitude_vector[NO_OF_CHOSEN_OBSERVATIONS_DRY + NO_OF_CHOSEN_OBSERVATIONS_MOIST + level_index*NO_OF_CHOSEN_WIND_POINTS_PER_LAYER_OBS + i]
			= latitudes_one_layer[chosen_indices_wind[i]];
			longitude_vector[NO_OF_CHOSEN_OBSERVATIONS_DRY + NO_OF_CHOSEN_OBSERVATIONS_MOIST + level_index*NO_OF_CHOSEN_WIND_POINTS_PER_LAYER_OBS + i]
			= longitudes_one_layer[chosen_indices_wind[i]];
			z_coords_amsl[NO_OF_CHOSEN_OBSERVATIONS_DRY + NO_OF_CHOSEN_OBSERVATIONS_MOIST + level_index*NO_OF_CHOSEN_WIND_POINTS_PER_LAYER_OBS + i]
			= z_height_amsl[chosen_indices_wind[i]];
			observations_vector[NO_OF_CHOSEN_OBSERVATIONS_DRY + NO_OF_CHOSEN_OBSERVATIONS_MOIST + level_index*NO_OF_CHOSEN_WIND_POINTS_PER_LAYER_OBS + i]
			= u_wind_one_layer[chosen_indices_wind[i]];
			
			// v wind
			latitude_vector[NO_OF_CHOSEN_OBSERVATIONS_DRY + NO_OF_CHOSEN_OBSERVATIONS_MOIST + NO_OF_CHOSEN_OBSERVATIONS_WIND/2 + level_index*NO_OF_CHOSEN_WIND_POINTS_PER_LAYER_OBS + i]
			= latitudes_one_layer[chosen_indices_wind[i]];
			longitude_vector[NO_OF_CHOSEN_OBSERVATIONS_DRY + NO_OF_CHOSEN_OBSERVATIONS_MOIST + NO_OF_CHOSEN_OBSERVATIONS_WIND/2 + level_index*NO_OF_CHOSEN_WIND_POINTS_PER_LAYER_OBS + i]
			= longitudes_one_layer[chosen_indices_wind[i]];
			z_coords_amsl[NO_OF_CHOSEN_OBSERVATIONS_DRY + NO_OF_CHOSEN_OBSERVATIONS_MOIST + NO_OF_CHOSEN_OBSERVATIONS_WIND/2 + level_index*NO_OF_CHOSEN_WIND_POINTS_PER_LAYER_OBS + i]
			= z_height_amsl[chosen_indices_wind[i]];
			observations_vector[NO_OF_CHOSEN_OBSERVATIONS_DRY + NO_OF_CHOSEN_OBSERVATIONS_MOIST + NO_OF_CHOSEN_OBSERVATIONS_WIND/2 + level_index*NO_OF_CHOSEN_WIND_POINTS_PER_LAYER_OBS + i]
			= v_wind_one_layer[chosen_indices_wind[i]];
		}
	}
	free(z_height_amsl);
	free(temperature_one_layer);
	free(u_wind_one_layer);
	free(v_wind_one_layer);
	free(chosen_indices_wind);
	
	// reading the surface height
	double *surface_height = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
	char SFC_OBS_FILE_PRE[200];
	sprintf(SFC_OBS_FILE_PRE , "%s/input/icon_global_icosahedral_time-invariant_%s%s%s%s_HSURF.grib2", game_da_root_dir, year_string, month_string, day_string, hour_string);
	char SFC_OBS_FILE[strlen(SFC_OBS_FILE_PRE) + 1];
	strcpy(SFC_OBS_FILE, SFC_OBS_FILE_PRE);
	
	ECC_FILE = fopen(SFC_OBS_FILE, "r");
	handle = codes_handle_new_from_file(NULL, ECC_FILE, PRODUCT_GRIB, &err);
	if (err != 0)
		ECCERR(err);
	NO_OF_POINTS_PER_LAYER_OBS_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_OBS;
	if ((retval = codes_get_double_array(handle, "values", &surface_height[0], &NO_OF_POINTS_PER_LAYER_OBS_SIZE_T)))
		ECCERR(retval);
	codes_handle_delete(handle);
	fclose(ECC_FILE);
	
	// surface pressure
	// reading the surface presure
	double *pressure_one_layer = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
	
	char SFC_PRES_FILE_PRE[200];
	sprintf(SFC_PRES_FILE_PRE , "%s/input/icon_global_icosahedral_single-level_%s%s%s%s_000_PS.grib2", game_da_root_dir, year_string, month_string, day_string, hour_string);
	char SFC_PRES_FILE[strlen(SFC_PRES_FILE_PRE) + 1];
	strcpy(SFC_PRES_FILE, SFC_PRES_FILE_PRE);
	
	ECC_FILE = fopen(SFC_PRES_FILE, "r");
	handle = codes_handle_new_from_file(NULL, ECC_FILE, PRODUCT_GRIB, &err);
	if (err != 0)
		ECCERR(err);
	NO_OF_POINTS_PER_LAYER_OBS_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_OBS;
	if ((retval = codes_get_double_array(handle, "values", &pressure_one_layer[0], &NO_OF_POINTS_PER_LAYER_OBS_SIZE_T)))
		ECCERR(retval);
	codes_handle_delete(handle);
	fclose(ECC_FILE);
	
	// writing the surface pressure to the observations
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_CHOSEN_POINTS_PER_LAYER_OBS; ++i)
	{
		latitude_vector[NO_OF_LEVELS_OBS*2*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS + i] = latitudes_one_layer[chosen_indices[i]];
		longitude_vector[NO_OF_LEVELS_OBS*2*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS + i] = longitudes_one_layer[chosen_indices[i]];
		z_coords_amsl[NO_OF_LEVELS_OBS*2*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS + i] = surface_height[chosen_indices[i]];
		observations_vector[NO_OF_LEVELS_OBS*2*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS + i] = pressure_one_layer[chosen_indices[i]];
	}
	
	// SST
	// reading the SST
	double *sst = malloc(NO_OF_SST_POINTS*sizeof(double));
	
	char SST_FILE_PRE[200];
	sprintf(SST_FILE_PRE , "%s/input/rtgssthr_grb_0.5.grib2", game_da_root_dir);
	char SST_FILE[strlen(SST_FILE_PRE) + 1];
	strcpy(SST_FILE, SST_FILE_PRE);
	
	ECC_FILE = fopen(SST_FILE, "r");
	handle = codes_handle_new_from_file(NULL, ECC_FILE, PRODUCT_GRIB, &err);
	if (err != 0)
		ECCERR(err);
	size_t NO_OF_SST_POINTS_SIZE_T = (size_t) NO_OF_SST_POINTS;
	if ((retval = codes_get_double_array(handle, "values", &sst[0], &NO_OF_SST_POINTS_SIZE_T)))
		ECCERR(retval);
	if ((retval = codes_get_double_array(handle, "latitudes", &latitudes_sst[0], &NO_OF_SST_POINTS_SIZE_T)))
		ECCERR(retval);
	if ((retval = codes_get_double_array(handle, "longitudes", &longitudes_sst[0], &NO_OF_SST_POINTS_SIZE_T)))
		ECCERR(retval);
	codes_handle_delete(handle);
	fclose(ECC_FILE);
	
	// writing the SST to the observations
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SST_POINTS; ++i)
	{
		latitude_vector[NO_OF_CHOSEN_OBSERVATIONS - NO_OF_SST_POINTS + i] = 2*M_PI*latitudes_sst[i]/360;
		longitude_vector[NO_OF_CHOSEN_OBSERVATIONS - NO_OF_SST_POINTS + i] = 2*M_PI*longitudes_sst[i]/360;
		z_coords_amsl[NO_OF_CHOSEN_OBSERVATIONS - NO_OF_SST_POINTS + i] = 0;
		observations_vector[NO_OF_CHOSEN_OBSERVATIONS - NO_OF_SST_POINTS + i] = sst[i];
	}
	
	// freeing the memory
	free(sst);
	free(latitudes_sst);
	free(longitudes_sst);
	free(latitudes_one_layer);
	free(longitudes_one_layer);
	free(surface_height);
	free(pressure_one_layer);
    free(chosen_indices);
    
    // Writing the observations to a netcdf file.
    char OUTPUT_FILE_PRE[200];
    sprintf(OUTPUT_FILE_PRE, "%s/input/obs_%s%s%s%s.nc", game_da_root_dir, year_string, month_string, day_string, hour_string);
	char OUTPUT_FILE[strlen(OUTPUT_FILE_PRE) + 1];
	strcpy(OUTPUT_FILE, OUTPUT_FILE_PRE);
    
    int ncid, observation_dimid, latitude_id, longitude_id, vert_id, obervations_id;
    if ((retval = nc_create(OUTPUT_FILE, NC_CLOBBER, &ncid)))
        NCERR(retval);
    // Defining the dimensions.
    if ((retval = nc_def_dim(ncid, "observation_index", NO_OF_CHOSEN_OBSERVATIONS, &observation_dimid)))
        NCERR(retval);
    // Defining the variables.
    if ((retval = nc_def_var(ncid, "latitude_vector", NC_DOUBLE, 1, &observation_dimid, &latitude_id)))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "longitude_vector", NC_DOUBLE, 1, &observation_dimid, &longitude_id)))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "z_coords_obs", NC_DOUBLE, 1, &observation_dimid, &vert_id)))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "observations_vector", NC_DOUBLE, 1, &observation_dimid, &obervations_id)))
        NCERR(retval);
    
    if ((retval = nc_enddef(ncid)))
        NCERR(retval);
        
    // Setting the variables.
    if ((retval = nc_put_var_double(ncid, latitude_id, &latitude_vector[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, longitude_id, &longitude_vector[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, vert_id, &z_coords_amsl[0])))
        NCERR(retval);    
    if ((retval = nc_put_var_double(ncid, obervations_id, &observations_vector[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid)))
    	NCERR(retval);
    
    // Freeing the memory.
	free(spec_hum_one_layer);
	free(latitude_vector);
	free(longitude_vector);
	free(z_coords_amsl);
	free(observations_vector);

	return 0;
}
