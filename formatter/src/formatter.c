/*
This source file is part of ndvar, which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/ndvar
*/

// types of observations
// 0: temperature
// 1: pressure
// 2: specific humidity

#include <netcdf.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}

const int NO_OF_OBS_LEVELS_OBS = 3;
const int NO_OF_FIELDS_PER_LAYER_OBS = 3;
const int NO_OF_SURFACE_FIELDS_OBS = 1;
const int NO_OF_POINTS_PER_LAYER_OBS = 2949120;
const int NO_OF_OBSERVATIONS = (NO_OF_OBS_LEVELS_OBS*NO_OF_FIELDS_PER_LAYER_OBS + NO_OF_SURFACE_FIELDS_OBS)*NO_OF_POINTS_PER_LAYER_OBS;

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
    char *ndvar_root_dir = malloc((len + 1)*sizeof(char));
    strcpy(ndvar_root_dir, argv[5]);
    
	// Allocating the memory.
	double *latitude_vector = malloc(NO_OF_OBSERVATIONS*sizeof(double));
	double *longitude_vector = malloc(NO_OF_OBSERVATIONS*sizeof(double));
	double *vert_vector = malloc(NO_OF_OBSERVATIONS*sizeof(double));
	double *observations_vector = malloc(NO_OF_OBSERVATIONS*sizeof(double));
	int *type_vector = malloc(NO_OF_OBSERVATIONS*sizeof(int));
	double *pressure_on_layer = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
	double *temperature_on_layer = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
	double *spec_hum_on_layer = malloc(NO_OF_POINTS_PER_LAYER_OBS*sizeof(double));
	
	// Freeing parts of the memory.
	free(pressure_on_layer);
	free(temperature_on_layer);
	free(spec_hum_on_layer);
	
    int OUTPUT_FILE_LENGTH = 100;
    char *OUTPUT_FILE_PRE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE_PRE, "%s/input/obs_%s%s%s%s.nc", ndvar_root_dir, year_string, month_string, day_string, hour_string);
    OUTPUT_FILE_LENGTH = strlen(OUTPUT_FILE_PRE);
    free(OUTPUT_FILE_PRE);
    char *OUTPUT_FILE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE, "%s/input/obs_%s%s%s%s.nc", ndvar_root_dir, year_string, month_string, day_string, hour_string);
    
    // Writing the observations to a netcdf file.
    int ncid, retval, observation_dimid, latitude_id, longitude_id, vert_id, obervations_id, type_id;
    if ((retval = nc_create(OUTPUT_FILE, NC_CLOBBER, &ncid)))
        NCERR(retval);
    // Defining the dimensions.
    if ((retval = nc_def_dim(ncid, "observation_index", NO_OF_OBSERVATIONS, &observation_dimid)))
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

	return 0;
}
