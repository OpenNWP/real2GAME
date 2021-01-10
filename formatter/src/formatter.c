/*
This source file is part of ndvar, which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/ndvar
*/

// types of observations
// 0: temperature
// 1: pressure
// 2: specific humidity
// 3: total precipitation rate

#include <netcdf.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "eccodes.h"
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}
#define ERRCODE 3
#define ECCERR(e) {printf("Error: Eccodes failed with error code %d. See http://download.ecmwf.int/test-data/eccodes/html/group__errors.html for meaning of the error codes.\n", e); exit(ERRCODE);}

// number of levels in the free atmosphere from which we want to use data
const int NO_OF_OBS_LEVELS_OBS = 1;
// number of fields per layer which we want to use
const int NO_OF_FIELDS_PER_LAYER_OBS = 1;
// no of fields at the surface we want to use
const int NO_OF_SURFACE_FIELDS_OBS = 0;
// the number of points per layer of the input model
const int NO_OF_POINTS_PER_LAYER_OBS = 2949120;
// the number of points per layer that are actually picked for the assimilation process
const int NO_OF_CHOSEN_POINTS_PER_LAYER = 2;
// the total number of available observations
const int NO_OF_OBSERVATIONS = (NO_OF_OBS_LEVELS_OBS*NO_OF_FIELDS_PER_LAYER_OBS + NO_OF_SURFACE_FIELDS_OBS)*NO_OF_POINTS_PER_LAYER_OBS;
// the number of observations we want to actually use
const int NO_OF_CHOSEN_OBSERVATIONS = (NO_OF_OBS_LEVELS_OBS*NO_OF_FIELDS_PER_LAYER_OBS + NO_OF_SURFACE_FIELDS_OBS)*NO_OF_CHOSEN_POINTS_PER_LAYER;

int main(int argc, char *argv[])
{
	// defining the levels of the model we want to use
	int levels_vector[NO_OF_OBS_LEVELS_OBS];
	levels_vector[0] = 90;
	
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
	double *latitude_vector_one_layer = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
	double *longitude_vector_one_layer = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
    
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
	if ((retval = codes_get_double_array(handle, "values", &latitude_vector_one_layer[0], &NO_OF_POINTS_PER_LAYER_OBS_SIZE_T)))
		ECCERR(retval);
	codes_handle_delete(handle);
    fclose(ECC_FILE);
    
    for (int i = 0; i < NO_OF_POINTS_PER_LAYER_OBS; ++i)
    {
    	latitude_vector_one_layer[i] = 2*M_PI*latitude_vector_one_layer[i]/360;
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
	if ((retval = codes_get_double_array(handle, "values", &longitude_vector_one_layer[0], &NO_OF_POINTS_PER_LAYER_OBS_SIZE_T)))
		ECCERR(retval);
	codes_handle_delete(handle);
    fclose(ECC_FILE);
    
    for (int i = 0; i < NO_OF_POINTS_PER_LAYER_OBS; ++i)
    {
    	longitude_vector_one_layer[i] = 2*M_PI*longitude_vector_one_layer[i]/360;
    }
    
	// Allocating the memory for the final result.
	double *latitude_vector = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(double));
	double *longitude_vector = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(double));
	double *vert_vector = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(double));
	double *observations_vector = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(double));
	int *type_vector = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(int));
	
	double *temperature_one_layer = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
	double *spec_hum_one_layer = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
	
	// the indices of the chosen points
	int *chosen_indices = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(int));
	
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS; ++i)
	{
		chosen_indices[i] = NO_OF_OBSERVATIONS/NO_OF_CHOSEN_OBSERVATIONS*i;
	}
	
	// reading the data from the free atmosphere
	double *vert_vector_free_atmosphere = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
	// loop over all relevant level in the free atmosphere
	for (int level_index = 0; level_index < NO_OF_OBS_LEVELS_OBS; ++level_index)
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
		if ((retval = codes_get_double_array(handle, "values", &vert_vector_free_atmosphere[0], &NO_OF_POINTS_PER_LAYER_OBS_SIZE_T)))
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
		for (int i = 0; i < NO_OF_CHOSEN_POINTS_PER_LAYER; ++i)
		{
			for (int j = 0; j < NO_OF_FIELDS_PER_LAYER_OBS; ++j)
			{
				latitude_vector[(level_index*NO_OF_FIELDS_PER_LAYER_OBS + j)*NO_OF_CHOSEN_POINTS_PER_LAYER + i] = latitude_vector_one_layer[chosen_indices[i]];
				longitude_vector[(level_index*NO_OF_FIELDS_PER_LAYER_OBS + j)*NO_OF_CHOSEN_POINTS_PER_LAYER + i] = longitude_vector_one_layer[chosen_indices[i]];
				vert_vector[(level_index*NO_OF_FIELDS_PER_LAYER_OBS + j)*NO_OF_CHOSEN_POINTS_PER_LAYER + i] = vert_vector_free_atmosphere[chosen_indices[i]];
				
				if (j == 0)
				{
					observations_vector[(level_index*NO_OF_FIELDS_PER_LAYER_OBS + j)*NO_OF_CHOSEN_POINTS_PER_LAYER + i] = temperature_one_layer[chosen_indices[i]];
				}
				if (j == 1)
				{
					observations_vector[(level_index*NO_OF_FIELDS_PER_LAYER_OBS + j)*NO_OF_CHOSEN_POINTS_PER_LAYER + i] = spec_hum_one_layer[chosen_indices[i]];
				}
				
				if (j == 0)
				{
					type_vector[(level_index*NO_OF_FIELDS_PER_LAYER_OBS + j)*NO_OF_CHOSEN_POINTS_PER_LAYER + i] = 0;
				}
				if (j == 1)
				{
					type_vector[(level_index*NO_OF_FIELDS_PER_LAYER_OBS + j)*NO_OF_CHOSEN_POINTS_PER_LAYER + i] = 2;
				}
			}
		}
	}
	
	// reading the surface height
	double *surface_height = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
	// this only needs to be done if we have surface fields
	if (NO_OF_SURFACE_FIELDS_OBS > 0)
	{
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
    }
	
	
	if (NO_OF_SURFACE_FIELDS_OBS == 2)
	{
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
		for (int i = 0; i < NO_OF_CHOSEN_POINTS_PER_LAYER; ++i)
		{
			latitude_vector[NO_OF_OBS_LEVELS_OBS*NO_OF_FIELDS_PER_LAYER_OBS*NO_OF_CHOSEN_POINTS_PER_LAYER + i] = latitude_vector_one_layer[chosen_indices[i]];
			
			longitude_vector[NO_OF_OBS_LEVELS_OBS*NO_OF_FIELDS_PER_LAYER_OBS*NO_OF_CHOSEN_POINTS_PER_LAYER + i] = longitude_vector_one_layer[chosen_indices[i]];
			
			vert_vector[NO_OF_OBS_LEVELS_OBS*NO_OF_FIELDS_PER_LAYER_OBS*NO_OF_CHOSEN_POINTS_PER_LAYER + i] = surface_height[chosen_indices[i]];
			
			observations_vector[NO_OF_OBS_LEVELS_OBS*NO_OF_FIELDS_PER_LAYER_OBS*NO_OF_CHOSEN_POINTS_PER_LAYER + i] = pressure_one_layer[chosen_indices[i]];
			
			type_vector[NO_OF_OBS_LEVELS_OBS*NO_OF_FIELDS_PER_LAYER_OBS*NO_OF_CHOSEN_POINTS_PER_LAYER + i] = 1;
		}
		
		// reading the total precipitation rate
		double *tot_prec = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
		
		int TOT_PREC_FILE_LENGTH = 100;
		char *TOT_PREC_FILE_PRE = malloc((TOT_PREC_FILE_LENGTH + 1)*sizeof(char));
		sprintf(TOT_PREC_FILE_PRE , "%s/input/icon_global_icosahedral_single-level_%s%s%s%s_001_TOT_PREC.grib2", ndvar_root_dir, year_string, month_string, day_string, hour_string);
		TOT_PREC_FILE_LENGTH = strlen(TOT_PREC_FILE_PRE);
		free(TOT_PREC_FILE_PRE);
		char *TOT_PREC_FILE = malloc((TOT_PREC_FILE_LENGTH + 1)*sizeof(char));
		sprintf(TOT_PREC_FILE, "%s/input/icon_global_icosahedral_single-level_%s%s%s%s_001_TOT_PREC.grib2", ndvar_root_dir, year_string, month_string, day_string, hour_string);
		
		ECC_FILE = fopen(TOT_PREC_FILE, "r");
		free(TOT_PREC_FILE);
		handle = codes_handle_new_from_file(NULL, ECC_FILE, PRODUCT_GRIB, &err);
		if (err != 0)
			ECCERR(err);
		NO_OF_POINTS_PER_LAYER_OBS_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_OBS;
		if ((retval = codes_get_double_array(handle, "values", &tot_prec[0], &NO_OF_POINTS_PER_LAYER_OBS_SIZE_T)))
			ECCERR(retval);
		codes_handle_delete(handle);
		fclose(ECC_FILE);
		
		// writing the total precipitation rate to the observations
		for (int i = 0; i < NO_OF_CHOSEN_POINTS_PER_LAYER; ++i)
		{
			latitude_vector[(NO_OF_OBS_LEVELS_OBS*NO_OF_FIELDS_PER_LAYER_OBS + 1)*NO_OF_CHOSEN_POINTS_PER_LAYER + i] = latitude_vector_one_layer[chosen_indices[i]];
			
			longitude_vector[(NO_OF_OBS_LEVELS_OBS*NO_OF_FIELDS_PER_LAYER_OBS + 1)*NO_OF_CHOSEN_POINTS_PER_LAYER + i] = longitude_vector_one_layer[chosen_indices[i]];
			
			vert_vector[(NO_OF_OBS_LEVELS_OBS*NO_OF_FIELDS_PER_LAYER_OBS + 1)*NO_OF_CHOSEN_POINTS_PER_LAYER + i] = surface_height[chosen_indices[i]];
			
			observations_vector[(NO_OF_OBS_LEVELS_OBS*NO_OF_FIELDS_PER_LAYER_OBS + 1)*NO_OF_CHOSEN_POINTS_PER_LAYER + i] = tot_prec[chosen_indices[i]];
			
			type_vector[(NO_OF_OBS_LEVELS_OBS*NO_OF_FIELDS_PER_LAYER_OBS + 1)*NO_OF_CHOSEN_POINTS_PER_LAYER + i] = 1;
		}
    
    	free(tot_prec);
		free(pressure_one_layer);
    
    }
    
    // Writing the observations to a netcdf file.
    int OUTPUT_FILE_LENGTH = 100;
    char *OUTPUT_FILE_PRE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE_PRE, "%s/input/obs_%s%s%s%s.nc", ndvar_root_dir, year_string, month_string, day_string, hour_string);
    OUTPUT_FILE_LENGTH = strlen(OUTPUT_FILE_PRE);
    free(OUTPUT_FILE_PRE);
    char *OUTPUT_FILE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE, "%s/input/obs_%s%s%s%s.nc", ndvar_root_dir, year_string, month_string, day_string, hour_string);
    
    int ncid, observation_dimid, latitude_id, longitude_id, vert_id, obervations_id, type_id;
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
    if ((retval = nc_def_var(ncid, "type_vector", NC_INT, 1, &observation_dimid, &type_id)))
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
    if ((retval = nc_put_var_int(ncid, type_id, &type_vector[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid)))
    	NCERR(retval);
    
    // Freeing the memory.
    free(chosen_indices);
	free(temperature_one_layer);
	free(spec_hum_one_layer);
    free(ndvar_root_dir);
	free(latitude_vector);
	free(longitude_vector);
	free(vert_vector);
	free(observations_vector);
	free(type_vector);
	free(year_string);
	free(month_string);
	free(day_string);
	free(hour_string);
	free(latitude_vector_one_layer);
	free(longitude_vector_one_layer);
	free(vert_vector_free_atmosphere);
	free(surface_height);

	return 0;
}
