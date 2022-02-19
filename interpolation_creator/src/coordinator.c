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
#include "../../interpolator/enum.h"
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}
#define ERRCODE 3
#define ECCERR(e) {printf("Error: Eccodes failed with error code %d. See http://download.ecmwf.int/test-data/eccodes/html/group__errors.html for meaning of the error codes.\n", e); exit(ERRCODE);}

const int NO_OF_POINTS_PER_LAYER_INPUT_MODEL = 2949120;

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
    char model_home_dir[strlen(argv[6]) + 1];
    strcpy(model_home_dir, argv[6]);
	int ORO_ID;
    ORO_ID = strtod(argv[7], NULL);
    
	double *latitudes_model = malloc(NO_OF_POINTS_PER_LAYER_INPUT_MODEL*sizeof(double));
	double *longitudes_model = malloc(NO_OF_POINTS_PER_LAYER_INPUT_MODEL*sizeof(double));
	
	int retval, err;
	codes_handle *handle = NULL;
	
	// Properties of the input model's grid.
	// latitudes of the grid
	char LAT_OBS_FILE_PRE[200];
    sprintf(LAT_OBS_FILE_PRE , "%s/interpolation_creator/icon_global_icosahedral_time-invariant_%s%s%s%s_CLAT.grib2", real2game_root_dir, year_string, month_string, day_string, hour_string);
	char LAT_OBS_FILE[strlen(LAT_OBS_FILE_PRE) + 1];
	strcpy(LAT_OBS_FILE, LAT_OBS_FILE_PRE);
	FILE *ECC_FILE;
	ECC_FILE = fopen(LAT_OBS_FILE, "r");
	handle = codes_handle_new_from_file(NULL, ECC_FILE, PRODUCT_GRIB, &err);
	if (err != 0)
		ECCERR(err);
	size_t NO_OF_POINTS_PER_LAYER_INPUT_MODEL_SIZE_T;
	NO_OF_POINTS_PER_LAYER_INPUT_MODEL_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_INPUT_MODEL;
	if ((retval = codes_get_double_array(handle, "values", &latitudes_model[0], &NO_OF_POINTS_PER_LAYER_INPUT_MODEL_SIZE_T)))
		ECCERR(retval);
	codes_handle_delete(handle);
    fclose(ECC_FILE);
    
    for (int i = 0; i < NO_OF_POINTS_PER_LAYER_INPUT_MODEL; ++i)
    {
    	latitudes_model[i] = 2*M_PI*latitudes_model[i]/360;
    }
    
    // longitudes of the grid
	char LON_OBS_FILE_PRE[200];
    sprintf(LON_OBS_FILE_PRE , "%s/interpolation_creator/icon_global_icosahedral_time-invariant_%s%s%s%s_CLON.grib2", real2game_root_dir, year_string, month_string, day_string, hour_string);
	char LON_OBS_FILE[strlen(LON_OBS_FILE_PRE) + 1];
	strcpy(LON_OBS_FILE, LON_OBS_FILE_PRE);
    
	ECC_FILE = fopen(LON_OBS_FILE, "r");
	handle = codes_handle_new_from_file(NULL, ECC_FILE, PRODUCT_GRIB, &err);
	if (err != 0)
		ECCERR(err);
	NO_OF_POINTS_PER_LAYER_INPUT_MODEL_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_INPUT_MODEL;
	if ((retval = codes_get_double_array(handle, "values", &longitudes_model[0], &NO_OF_POINTS_PER_LAYER_INPUT_MODEL_SIZE_T)))
		ECCERR(retval);
	codes_handle_delete(handle);
    fclose(ECC_FILE);
    
    for (int i = 0; i < NO_OF_POINTS_PER_LAYER_INPUT_MODEL; ++i)
    {
    	longitudes_model[i] = 2*M_PI*longitudes_model[i]/360;
    }
    
    
    // reading the horizontal coordinates of the grid of GAME
    double *latitudes_game = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *longitudes_game = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *latitudes_game_wind = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *longitudes_game_wind = malloc(NO_OF_VECTORS_H*sizeof(double));
    char GEO_PROP_FILE_PRE[200];
    sprintf(GEO_PROP_FILE_PRE, "%s/grid_generator/grids/RES%d_L%d_ORO%d.nc", model_home_dir, RES_ID, NO_OF_LAYERS, ORO_ID);
    char GEO_PROP_FILE[strlen(GEO_PROP_FILE_PRE) + 1];
    strcpy(GEO_PROP_FILE, GEO_PROP_FILE_PRE);
	printf("Grid file: %s\n", GEO_PROP_FILE);
	printf("Reading grid file of GAME ...\n");
	int ncid_grid;
    if ((retval = nc_open(GEO_PROP_FILE, NC_NOWRITE, &ncid_grid)))
        NCERR(retval);
    
    int latitudes_game_id, latitudes_game_wind_id, longitudes_game_id, longitudes_game_wind_id;
    if ((retval = nc_inq_varid(ncid_grid, "latitude_scalar", &latitudes_game_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "longitude_scalar", &longitudes_game_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "latitude_vector", &latitudes_game_wind_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "longitude_vector", &longitudes_game_wind_id)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, latitudes_game_id, &latitudes_game[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, longitudes_game_id, &longitudes_game[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, latitudes_game_wind_id, &latitudes_game_wind[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, longitudes_game_wind_id, &longitudes_game_wind[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid_grid)))
        NCERR(retval);
	printf("Grid file of GAME read.\n");
	
	// writing the result to a NetCDF file
	
    
	free(latitudes_model);
	free(longitudes_model);
    free(latitudes_game);
    free(longitudes_game);
    free(latitudes_game_wind);
    free(longitudes_game_wind);
	
	return 0;
}


