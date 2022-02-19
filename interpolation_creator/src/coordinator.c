/*
This source file is part of real2GAME, which is released under the MIT license.
Github repository: https://github.com/OpenNWP/real2GAME
*/

/*
This file prepares the horizontal interpolation from the foreign model to GAME.
*/

#include <netcdf.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "eccodes.h"
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}
#define ERRCODE 3
#define ECCERR(e) {printf("Error: Eccodes failed with error code %d. See http://download.ecmwf.int/test-data/eccodes/html/group__errors.html for meaning of the error codes.\n", e); exit(ERRCODE);}

const int NO_OF_POINTS_PER_LAYER_ICON = 2949120;

int main(int argc, char *argv[])
{
	
	// shell arguments
    char year_string[strlen(argv[1]) + 1];
    strcpy(year_string, argv[1]);
    char month_string[strlen(argv[2]) + 1];
    strcpy(month_string, argv[2]);
    char day_string[strlen(argv[3]) + 1];
    strcpy(day_string, argv[3]);
    char hour_string[strlen(argv[4]) + 1];
    strcpy(hour_string, argv[4]);
    char real2game_root_dir[strlen(argv[5]) + 1];
    strcpy(real2game_root_dir, argv[5]);
    
	double *latitudes_model = malloc(NO_OF_POINTS_PER_LAYER_ICON*sizeof(double));
	double *longitudes_model = malloc(NO_OF_POINTS_PER_LAYER_ICON*sizeof(double));
	
	int retval, err;
	codes_handle *handle = NULL;
	
	// Properties of the input model's grid.
	// latitudes of the grid
	char LAT_OBS_FILE_PRE[200];
    sprintf(LAT_OBS_FILE_PRE , "%s/input/icon_global_icosahedral_time-invariant_%s%s%s%s_CLAT.grib2", real2game_root_dir, year_string, month_string, day_string, hour_string);
	char LAT_OBS_FILE[strlen(LAT_OBS_FILE_PRE) + 1];
	strcpy(LAT_OBS_FILE, LAT_OBS_FILE_PRE);
    
	FILE *ECC_FILE;
	ECC_FILE = fopen(LAT_OBS_FILE, "r");
	handle = codes_handle_new_from_file(NULL, ECC_FILE, PRODUCT_GRIB, &err);
	if (err != 0)
		ECCERR(err);
	size_t NO_OF_POINTS_PER_LAYER_ICON_SIZE_T;
	NO_OF_POINTS_PER_LAYER_ICON_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_ICON;
	if ((retval = codes_get_double_array(handle, "values", &latitudes_model[0], &NO_OF_POINTS_PER_LAYER_ICON_SIZE_T)))
		ECCERR(retval);
	codes_handle_delete(handle);
    fclose(ECC_FILE);
    
    for (int i = 0; i < NO_OF_POINTS_PER_LAYER_ICON; ++i)
    {
    	latitudes_model[i] = 2*M_PI*latitudes_model[i]/360;
    }
    
    // longitudes of the grid
	char LON_OBS_FILE_PRE[200];
    sprintf(LON_OBS_FILE_PRE , "%s/input/icon_global_icosahedral_time-invariant_%s%s%s%s_CLON.grib2", real2game_root_dir, year_string, month_string, day_string, hour_string);
	char LON_OBS_FILE[strlen(LON_OBS_FILE_PRE) + 1];
	strcpy(LON_OBS_FILE, LON_OBS_FILE_PRE);
    
	ECC_FILE = fopen(LON_OBS_FILE, "r");
	handle = codes_handle_new_from_file(NULL, ECC_FILE, PRODUCT_GRIB, &err);
	if (err != 0)
		ECCERR(err);
	NO_OF_POINTS_PER_LAYER_ICON_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_ICON;
	if ((retval = codes_get_double_array(handle, "values", &longitudes_model[0], &NO_OF_POINTS_PER_LAYER_ICON_SIZE_T)))
		ECCERR(retval);
	codes_handle_delete(handle);
    fclose(ECC_FILE);
    
    for (int i = 0; i < NO_OF_POINTS_PER_LAYER_ICON; ++i)
    {
    	longitudes_model[i] = 2*M_PI*longitudes_model[i]/360;
    }
    
    
    
	free(latitudes_model);
	free(longitudes_model);
	
	return 0;
}


