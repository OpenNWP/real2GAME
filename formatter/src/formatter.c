/*
This source file is part of ndvar, which is released under the MIT license.
Github repository: https://github.com/OpenNWP/ndvar
*/

#include <netcdf.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "eccodes.h"
#include "../../src/ndvar.h"
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
    char *ndvar_root_dir = malloc((len + 1)*sizeof(char));
    strcpy(ndvar_root_dir, argv[5]);
	
	// Properties of the input model's grid.
	double *latitudes_one_layer = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
	double *longitudes_one_layer = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
    
	int retval, err;
	codes_handle *handle = NULL;
	
	// latitudes of the grid
	int LAT_OBS_FILE_LENGTH = 100;
    char *LAT_OBS_FILE_PRE = malloc((LAT_OBS_FILE_LENGTH + 1)*sizeof(char));
    sprintf(LAT_OBS_FILE_PRE , "%s/input/icon_global_icosahedral_time-invariant_%s%s%s%s_CLAT.grib2", ndvar_root_dir, year_string, month_string, day_string, hour_string);
    LAT_OBS_FILE_LENGTH = strlen(LAT_OBS_FILE_PRE);
	free(LAT_OBS_FILE_PRE);
    char *LAT_OBS_FILE = malloc((LAT_OBS_FILE_LENGTH + 1)*sizeof(char));
    sprintf(LAT_OBS_FILE, "%s/input/icon_global_icosahedral_time-invariant_%s%s%s%s_CLAT.grib2", ndvar_root_dir, year_string, month_string, day_string, hour_string);
    
	FILE *ECC_FILE;
	ECC_FILE = fopen(LAT_OBS_FILE, "r");
	free(LAT_OBS_FILE);
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
	int LON_OBS_FILE_LENGTH = 100;
    char *LON_OBS_FILE_PRE = malloc((LON_OBS_FILE_LENGTH + 1)*sizeof(char));
    sprintf(LON_OBS_FILE_PRE , "%s/input/icon_global_icosahedral_time-invariant_%s%s%s%s_CLON.grib2", ndvar_root_dir, year_string, month_string, day_string, hour_string);
    LON_OBS_FILE_LENGTH = strlen(LON_OBS_FILE_PRE);
	free(LON_OBS_FILE_PRE);
    char *LON_OBS_FILE = malloc((LON_OBS_FILE_LENGTH + 1)*sizeof(char));
    sprintf(LON_OBS_FILE, "%s/input/icon_global_icosahedral_time-invariant_%s%s%s%s_CLON.grib2", ndvar_root_dir, year_string, month_string, day_string, hour_string);
    
	ECC_FILE = fopen(LON_OBS_FILE, "r");
	free(LON_OBS_FILE);
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
	double *vert_vector = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(double));
	double *observations_vector = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(double));
	
	double *temperature_one_layer = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
	double *spec_hum_one_layer = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
	
	// the indices of the chosen points
	int *chosen_indices = malloc(NO_OF_CHOSEN_POINTS_PER_LAYER_OBS*sizeof(int));
	
	for (int i = 0; i < NO_OF_CHOSEN_POINTS_PER_LAYER_OBS; ++i)
	{
		chosen_indices[i] = NO_OF_POINTS_PER_LAYER_OBS/NO_OF_CHOSEN_POINTS_PER_LAYER_OBS*i;
	}
	
	// reading the data from the free atmosphere
	double *z_height_amsl = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
	// loop over all relevant level in the free atmosphere
	for (int level_index = 0; level_index < NO_OF_LEVELS_OBS; ++level_index)
	{
		// vertical position of the current layer
		int Z_OBS_FILE_LENGTH = 100;
		char *Z_OBS_FILE_PRE = malloc((Z_OBS_FILE_LENGTH + 1)*sizeof(char));
		sprintf(Z_OBS_FILE_PRE , "%s/input/icon_global_icosahedral_time-invariant_%s%s%s%s_%d_HHL.grib2", ndvar_root_dir, year_string, month_string, day_string, hour_string, levels_vector[level_index]);
		Z_OBS_FILE_LENGTH = strlen(Z_OBS_FILE_PRE);
		free(Z_OBS_FILE_PRE);
		char *Z_OBS_FILE = malloc((Z_OBS_FILE_LENGTH + 1)*sizeof(char));
		sprintf(Z_OBS_FILE, "%s/input/icon_global_icosahedral_time-invariant_%s%s%s%s_%d_HHL.grib2", ndvar_root_dir, year_string, month_string, day_string, hour_string, levels_vector[level_index]);
		
		ECC_FILE = fopen(Z_OBS_FILE, "r");
		free(Z_OBS_FILE);
		handle = codes_handle_new_from_file(NULL, ECC_FILE, PRODUCT_GRIB, &err);
		if (err != 0)
			ECCERR(err);
		NO_OF_POINTS_PER_LAYER_OBS_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_OBS;
		if ((retval = codes_get_double_array(handle, "values", &z_height_amsl[0], &NO_OF_POINTS_PER_LAYER_OBS_SIZE_T)))
			ECCERR(retval);
		codes_handle_delete(handle);
		fclose(ECC_FILE);
		
	   	// reading the temperature
		int TEMPERATURE_FILE_LENGTH = 100;
		char *TEMPERATURE_FILE_PRE = malloc((TEMPERATURE_FILE_LENGTH + 1)*sizeof(char));
		sprintf(TEMPERATURE_FILE_PRE , "%s/input/icon_global_icosahedral_model-level_%s%s%s%s_000_%d_T.grib2", ndvar_root_dir, year_string, month_string, day_string, hour_string, levels_vector[level_index]);
		TEMPERATURE_FILE_LENGTH = strlen(TEMPERATURE_FILE_PRE);
		free(TEMPERATURE_FILE_PRE);
		char *TEMPERATURE_FILE = malloc((TEMPERATURE_FILE_LENGTH + 1)*sizeof(char));
		sprintf(TEMPERATURE_FILE, "%s/input/icon_global_icosahedral_model-level_%s%s%s%s_000_%d_T.grib2", ndvar_root_dir, year_string, month_string, day_string, hour_string, levels_vector[level_index]);
		
		ECC_FILE = fopen(TEMPERATURE_FILE, "r");
		free(TEMPERATURE_FILE);
		handle = codes_handle_new_from_file(NULL, ECC_FILE, PRODUCT_GRIB, &err);
		if (err != 0)
			ECCERR(err);
		NO_OF_POINTS_PER_LAYER_OBS_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_OBS;
		if ((retval = codes_get_double_array(handle, "values", &temperature_one_layer[0], &NO_OF_POINTS_PER_LAYER_OBS_SIZE_T)))
			ECCERR(retval);
		codes_handle_delete(handle);
		fclose(ECC_FILE);
		
	   	// reading the specific humidity
		int SPEC_HUM_FILE_LENGTH = 100;
		char *SPEC_HUM_FILE_PRE = malloc((SPEC_HUM_FILE_LENGTH + 1)*sizeof(char));
		sprintf(SPEC_HUM_FILE_PRE , "%s/input/icon_global_icosahedral_model-level_%s%s%s%s_000_%d_QV.grib2", ndvar_root_dir, year_string, month_string, day_string, hour_string, levels_vector[level_index]);
		SPEC_HUM_FILE_LENGTH = strlen(SPEC_HUM_FILE_PRE);
		free(SPEC_HUM_FILE_PRE);
		char *SPEC_HUM_FILE = malloc((SPEC_HUM_FILE_LENGTH + 1)*sizeof(char));
		sprintf(SPEC_HUM_FILE, "%s/input/icon_global_icosahedral_model-level_%s%s%s%s_000_%d_QV.grib2", ndvar_root_dir, year_string, month_string, day_string, hour_string, levels_vector[level_index]);
		
		ECC_FILE = fopen(SPEC_HUM_FILE, "r");
		free(SPEC_HUM_FILE);
		handle = codes_handle_new_from_file(NULL, ECC_FILE, PRODUCT_GRIB, &err);
		if (err != 0)
			ECCERR(err);
		NO_OF_POINTS_PER_LAYER_OBS_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_OBS;
		if ((retval = codes_get_double_array(handle, "values", &spec_hum_one_layer[0], &NO_OF_POINTS_PER_LAYER_OBS_SIZE_T)))
			ECCERR(retval);
		codes_handle_delete(handle);
		fclose(ECC_FILE);
		
		// formatting the observations
		for (int i = 0; i < NO_OF_CHOSEN_POINTS_PER_LAYER_OBS; ++i)
		{
			for (int j = 0; j < 2; ++j)
			{
				latitude_vector[(level_index*2 + j)*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS + i] = latitudes_one_layer[chosen_indices[i]];
				longitude_vector[(level_index*2 + j)*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS + i] = longitudes_one_layer[chosen_indices[i]];
				vert_vector[(level_index*2 + j)*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS + i] = z_height_amsl[chosen_indices[i]];
				
				if (j == 0)
				{
					observations_vector[(level_index*2 + j)*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS + i] = temperature_one_layer[chosen_indices[i]];
				}
				if (j == 1)
				{
					observations_vector[(level_index*2 + j)*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS + i] = spec_hum_one_layer[chosen_indices[i]];
				}
			}
		}
	}
	
	// reading the surface height
	double *surface_height = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
	int SFC_OBS_FILE_LENGTH = 100;
	char *SFC_OBS_FILE_PRE = malloc((SFC_OBS_FILE_LENGTH + 1)*sizeof(char));
	sprintf(SFC_OBS_FILE_PRE , "%s/input/icon_global_icosahedral_time-invariant_%s%s%s%s_HSURF.grib2", ndvar_root_dir, year_string, month_string, day_string, hour_string);
	SFC_OBS_FILE_LENGTH = strlen(SFC_OBS_FILE_PRE);
	free(SFC_OBS_FILE_PRE);
	char *SFC_OBS_FILE = malloc((SFC_OBS_FILE_LENGTH + 1)*sizeof(char));
	sprintf(SFC_OBS_FILE, "%s/input/icon_global_icosahedral_time-invariant_%s%s%s%s_HSURF.grib2", ndvar_root_dir, year_string, month_string, day_string, hour_string);
	
	ECC_FILE = fopen(SFC_OBS_FILE, "r");
	free(SFC_OBS_FILE);
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
	
	int SFC_PRES_FILE_LENGTH = 100;
	char *SFC_PRES_FILE_PRE = malloc((SFC_PRES_FILE_LENGTH + 1)*sizeof(char));
	sprintf(SFC_PRES_FILE_PRE , "%s/input/icon_global_icosahedral_single-level_%s%s%s%s_000_PS.grib2", ndvar_root_dir, year_string, month_string, day_string, hour_string);
	SFC_PRES_FILE_LENGTH = strlen(SFC_PRES_FILE_PRE);
	free(SFC_PRES_FILE_PRE);
	char *SFC_PRES_FILE = malloc((SFC_PRES_FILE_LENGTH + 1)*sizeof(char));
	sprintf(SFC_PRES_FILE, "%s/input/icon_global_icosahedral_single-level_%s%s%s%s_000_PS.grib2", ndvar_root_dir, year_string, month_string, day_string, hour_string);
	
	ECC_FILE = fopen(SFC_PRES_FILE, "r");
	free(SFC_PRES_FILE);
	handle = codes_handle_new_from_file(NULL, ECC_FILE, PRODUCT_GRIB, &err);
	if (err != 0)
		ECCERR(err);
	NO_OF_POINTS_PER_LAYER_OBS_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_OBS;
	if ((retval = codes_get_double_array(handle, "values", &pressure_one_layer[0], &NO_OF_POINTS_PER_LAYER_OBS_SIZE_T)))
		ECCERR(retval);
	codes_handle_delete(handle);
	fclose(ECC_FILE);
	
	// writing the surface pressure to the observations
	for (int i = 0; i < NO_OF_CHOSEN_POINTS_PER_LAYER_OBS; ++i)
	{
		latitude_vector[NO_OF_LEVELS_OBS*2*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS + i] = latitudes_one_layer[chosen_indices[i]];
		longitude_vector[NO_OF_LEVELS_OBS*2*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS + i] = longitudes_one_layer[chosen_indices[i]];
		vert_vector[NO_OF_LEVELS_OBS*2*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS + i] = surface_height[chosen_indices[i]];
		observations_vector[NO_OF_LEVELS_OBS*2*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS + i] = pressure_one_layer[chosen_indices[i]];
	}
	
	free(pressure_one_layer);
    free(chosen_indices);
    
    // Writing the observations to a netcdf file.
    int OUTPUT_FILE_LENGTH = 100;
    char *OUTPUT_FILE_PRE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE_PRE, "%s/input/obs_%s%s%s%s.nc", ndvar_root_dir, year_string, month_string, day_string, hour_string);
    OUTPUT_FILE_LENGTH = strlen(OUTPUT_FILE_PRE);
    free(OUTPUT_FILE_PRE);
    char *OUTPUT_FILE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE, "%s/input/obs_%s%s%s%s.nc", ndvar_root_dir, year_string, month_string, day_string, hour_string);
    
    int ncid, observation_dimid, latitude_id, longitude_id, vert_id, obervations_id;
    if ((retval = nc_create(OUTPUT_FILE, NC_CLOBBER, &ncid)))
        NCERR(retval);
    free(OUTPUT_FILE);
    // Defining the dimensions.
    if ((retval = nc_def_dim(ncid, "observation_index", NO_OF_CHOSEN_OBSERVATIONS, &observation_dimid)))
        NCERR(retval);
    // Defining the variables.
    if ((retval = nc_def_var(ncid, "latitude_vector", NC_DOUBLE, 1, &observation_dimid, &latitude_id)))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "longitude_vector", NC_DOUBLE, 1, &observation_dimid, &longitude_id)))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "vert_vector", NC_DOUBLE, 1, &observation_dimid, &vert_id)))
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
    if ((retval = nc_put_var_double(ncid, vert_id, &vert_vector[0])))
        NCERR(retval);    
    if ((retval = nc_put_var_double(ncid, obervations_id, &observations_vector[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid)))
    	NCERR(retval);
    
    // Freeing the memory.
	free(temperature_one_layer);
	free(spec_hum_one_layer);
    free(ndvar_root_dir);
	free(latitude_vector);
	free(longitude_vector);
	free(vert_vector);
	free(observations_vector);
	free(year_string);
	free(month_string);
	free(day_string);
	free(hour_string);
	free(latitudes_one_layer);
	free(longitudes_one_layer);
	free(z_height_amsl);
	free(surface_height);

	return 0;
}
