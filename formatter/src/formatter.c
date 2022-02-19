/*
This source file is part of real2GAME, which is released under the MIT license.
Github repository: https://github.com/OpenNWP/real2GAME
*/

#include <netcdf.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "eccodes.h"
#include "../../header.h"
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}
#define ERRCODE 3
#define ECCERR(e) {printf("Error: Eccodes failed with error code %d. See http://download.ecmwf.int/test-data/eccodes/html/group__errors.html for meaning of the error codes.\n", e); exit(ERRCODE);}

int main(int argc, char *argv[])
{
	// defining the levels of the model we want to use
	int levels_vector[NO_OF_LEVELS_INPUT];
	levels_vector[0] = 1;
	levels_vector[1] = 10;
	levels_vector[2] = 19;
	levels_vector[3] = 27;
	levels_vector[4] = 35;
	levels_vector[5] = 43;
	levels_vector[6] = 51;
	levels_vector[7] = 59;
	levels_vector[8] = 67;
	levels_vector[9] = 75;
	levels_vector[10] = 83;
	levels_vector[11] = 90;
	
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
    
	// Allocating the memory for the final result.
	double *z_height_amsl_single_layer = malloc(NO_OF_POINTS_PER_LAYER_INPUT*sizeof(double));
	
	// single-layer arrays from grib
	double *temperature_one_layer = malloc(NO_OF_POINTS_PER_LAYER_INPUT*sizeof(double));
	double *spec_hum_one_layer = malloc(NO_OF_POINTS_PER_LAYER_INPUT*sizeof(double));
	double *u_one_layer = malloc(NO_OF_POINTS_PER_LAYER_INPUT*sizeof(double));
	double *v_one_layer = malloc(NO_OF_POINTS_PER_LAYER_INPUT*sizeof(double));
	
	// 2D-arrays for NetCDF
	double (*temperature)[NO_OF_LEVELS_INPUT] = malloc(sizeof(double[NO_OF_POINTS_PER_LAYER_INPUT][NO_OF_LEVELS_INPUT]));
	double (*spec_hum)[NO_OF_LEVELS_INPUT] = malloc(sizeof(double[NO_OF_POINTS_PER_LAYER_INPUT][NO_OF_LEVELS_INPUT]));
	double (*u_wind)[NO_OF_LEVELS_INPUT] = malloc(sizeof(double[NO_OF_POINTS_PER_LAYER_INPUT][NO_OF_LEVELS_INPUT]));
	double (*v_wind)[NO_OF_LEVELS_INPUT] = malloc(sizeof(double[NO_OF_POINTS_PER_LAYER_INPUT][NO_OF_LEVELS_INPUT]));
	double (*z_height_amsl)[NO_OF_LEVELS_INPUT] = malloc(sizeof(double[NO_OF_POINTS_PER_LAYER_INPUT][NO_OF_LEVELS_INPUT]));
	
	// grib stuff
	FILE *ecc_file;
	size_t NO_OF_POINTS_PER_LAYER_INPUT_SIZE_T;
	int retval, err;
	codes_handle *handle = NULL;
	
	// reading the data from the free atmosphere
	// loop over all relevant levels in the free atmosphere
	for (int level_index = 0; level_index < NO_OF_LEVELS_INPUT; ++level_index)
	{
		// vertical position of the current layer
		char Z_OBS_FILE_PRE[200];
		sprintf(Z_OBS_FILE_PRE , "%s/input/icon_global_icosahedral_time-invariant_%s%s%s%s_%d_HHL.grib2",
		real2game_root_dir, year_string, month_string, day_string, hour_string, levels_vector[level_index]);
		char Z_OBS_FILE[strlen(Z_OBS_FILE_PRE) + 1];
		strcpy(Z_OBS_FILE, Z_OBS_FILE_PRE);
		
		ecc_file = fopen(Z_OBS_FILE, "r");
		handle = codes_handle_new_from_file(NULL, ecc_file, PRODUCT_GRIB, &err);
		if (err != 0)
			ECCERR(err);
		NO_OF_POINTS_PER_LAYER_INPUT_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_INPUT;
		if ((retval = codes_get_double_array(handle, "values", &z_height_amsl_single_layer[0], &NO_OF_POINTS_PER_LAYER_INPUT_SIZE_T)))
			ECCERR(retval);
		codes_handle_delete(handle);
		fclose(ecc_file);
		
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_POINTS_PER_LAYER_INPUT; ++i)
		{
			z_height_amsl[i][level_index] = z_height_amsl_single_layer[i];
		}
		
	   	// reading the temperature
		char TEMPERATURE_FILE_PRE[200];
		sprintf(TEMPERATURE_FILE_PRE , "%s/input/icon_global_icosahedral_model-level_%s%s%s%s_000_%d_T.grib2",
		real2game_root_dir, year_string, month_string, day_string, hour_string, levels_vector[level_index]);
		char TEMPERATURE_FILE[strlen(TEMPERATURE_FILE_PRE) + 1];
		strcpy(TEMPERATURE_FILE, TEMPERATURE_FILE_PRE);
		
		ecc_file = fopen(TEMPERATURE_FILE, "r");
		handle = codes_handle_new_from_file(NULL, ecc_file, PRODUCT_GRIB, &err);
		if (err != 0)
			ECCERR(err);
		NO_OF_POINTS_PER_LAYER_INPUT_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_INPUT;
		if ((retval = codes_get_double_array(handle, "values", &temperature_one_layer[0], &NO_OF_POINTS_PER_LAYER_INPUT_SIZE_T)))
			ECCERR(retval);
		codes_handle_delete(handle);
		fclose(ecc_file);
		
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_POINTS_PER_LAYER_INPUT; ++i)
		{
			temperature[i][level_index] = temperature_one_layer[i];
		}
		
	   	// reading the specific humidity
		char SPEC_HUM_FILE_PRE[200];
		sprintf(SPEC_HUM_FILE_PRE , "%s/input/icon_global_icosahedral_model-level_%s%s%s%s_000_%d_QV.grib2",
		real2game_root_dir, year_string, month_string, day_string, hour_string, levels_vector[level_index]);
		char SPEC_HUM_FILE[strlen(SPEC_HUM_FILE_PRE) + 1];
		strcpy(SPEC_HUM_FILE, SPEC_HUM_FILE_PRE);
		
		ecc_file = fopen(SPEC_HUM_FILE, "r");
		handle = codes_handle_new_from_file(NULL, ecc_file, PRODUCT_GRIB, &err);
		if (err != 0)
			ECCERR(err);
		NO_OF_POINTS_PER_LAYER_INPUT_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_INPUT;
		if ((retval = codes_get_double_array(handle, "values", &spec_hum_one_layer[0], &NO_OF_POINTS_PER_LAYER_INPUT_SIZE_T)))
			ECCERR(retval);
		codes_handle_delete(handle);
		fclose(ecc_file);
		
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_POINTS_PER_LAYER_INPUT; ++i)
		{
			spec_hum[i][level_index] = spec_hum_one_layer[i];
		}
		
	   	// reading the u wind
		char U_WIND_FILE_PRE[200];
		sprintf(U_WIND_FILE_PRE , "%s/input/icon_global_icosahedral_model-level_%s%s%s%s_000_%d_U.grib2",
		real2game_root_dir, year_string, month_string, day_string, hour_string, levels_vector[level_index]);
		char U_WIND_FILE[strlen(U_WIND_FILE_PRE) + 1];
		strcpy(U_WIND_FILE, U_WIND_FILE_PRE);
		
		ecc_file = fopen(U_WIND_FILE, "r");
		handle = codes_handle_new_from_file(NULL, ecc_file, PRODUCT_GRIB, &err);
		if (err != 0)
			ECCERR(err);
		NO_OF_POINTS_PER_LAYER_INPUT_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_INPUT;
		if ((retval = codes_get_double_array(handle, "values", &u_one_layer[0], &NO_OF_POINTS_PER_LAYER_INPUT_SIZE_T)))
			ECCERR(retval);
		codes_handle_delete(handle);
		fclose(ecc_file);
		
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_POINTS_PER_LAYER_INPUT; ++i)
		{
			u_wind[i][level_index] = u_one_layer[i];
		}
		
	   	// reading the v wind
		char V_WIND_FILE_PRE[200];
		sprintf(V_WIND_FILE_PRE , "%s/input/icon_global_icosahedral_model-level_%s%s%s%s_000_%d_V.grib2",
		real2game_root_dir, year_string, month_string, day_string, hour_string, levels_vector[level_index]);
		char V_WIND_FILE[strlen(V_WIND_FILE_PRE) + 1];
		strcpy(V_WIND_FILE, V_WIND_FILE_PRE);
		
		ecc_file = fopen(V_WIND_FILE, "r");
		handle = codes_handle_new_from_file(NULL, ecc_file, PRODUCT_GRIB, &err);
		if (err != 0)
			ECCERR(err);
		NO_OF_POINTS_PER_LAYER_INPUT_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_INPUT;
		if ((retval = codes_get_double_array(handle, "values", &v_one_layer[0], &NO_OF_POINTS_PER_LAYER_INPUT_SIZE_T)))
			ECCERR(retval);
		codes_handle_delete(handle);
		fclose(ecc_file);
		
		#pragma omp parallel for
		for (int i = 0; i < NO_OF_POINTS_PER_LAYER_INPUT; ++i)
		{
			v_wind[i][level_index] = v_one_layer[i];
		}
		
	}
	free(z_height_amsl_single_layer);
	free(temperature_one_layer);
	free(spec_hum_one_layer);
	free(u_one_layer);
	free(v_one_layer);
	
	// reading the surface height
	double *surface_height = malloc(NO_OF_POINTS_PER_LAYER_INPUT*sizeof(double));
	char SFC_OBS_FILE_PRE[200];
	sprintf(SFC_OBS_FILE_PRE , "%s/input/icon_global_icosahedral_time-invariant_%s%s%s%s_HSURF.grib2", real2game_root_dir, year_string, month_string, day_string, hour_string);
	char SFC_OBS_FILE[strlen(SFC_OBS_FILE_PRE) + 1];
	strcpy(SFC_OBS_FILE, SFC_OBS_FILE_PRE);
	
	ecc_file = fopen(SFC_OBS_FILE, "r");
	handle = codes_handle_new_from_file(NULL, ecc_file, PRODUCT_GRIB, &err);
	if (err != 0)
		ECCERR(err);
	NO_OF_POINTS_PER_LAYER_INPUT_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_INPUT;
	if ((retval = codes_get_double_array(handle, "values", &surface_height[0], &NO_OF_POINTS_PER_LAYER_INPUT_SIZE_T)))
		ECCERR(retval);
	codes_handle_delete(handle);
	fclose(ecc_file);
	
	// reading the surface presure
	double *pressure_surface = malloc(NO_OF_POINTS_PER_LAYER_INPUT*sizeof(double));
	
	char SFC_PRES_FILE_PRE[200];
	sprintf(SFC_PRES_FILE_PRE , "%s/input/icon_global_icosahedral_single-level_%s%s%s%s_000_PS.grib2", real2game_root_dir, year_string, month_string, day_string, hour_string);
	char SFC_PRES_FILE[strlen(SFC_PRES_FILE_PRE) + 1];
	strcpy(SFC_PRES_FILE, SFC_PRES_FILE_PRE);
	
	ecc_file = fopen(SFC_PRES_FILE, "r");
	handle = codes_handle_new_from_file(NULL, ecc_file, PRODUCT_GRIB, &err);
	if (err != 0)
		ECCERR(err);
	NO_OF_POINTS_PER_LAYER_INPUT_SIZE_T = (size_t) NO_OF_POINTS_PER_LAYER_INPUT;
	if ((retval = codes_get_double_array(handle, "values", &pressure_surface[0], &NO_OF_POINTS_PER_LAYER_INPUT_SIZE_T)))
		ECCERR(retval);
	codes_handle_delete(handle);
	fclose(ecc_file);
	
	// reading the SST
	double *latitudes_sst = malloc(NO_OF_SST_POINTS*sizeof(double));
	double *longitudes_sst = malloc(NO_OF_SST_POINTS*sizeof(double));
	double *sst = malloc(NO_OF_SST_POINTS*sizeof(double));
	
	char SST_FILE_PRE[200];
	sprintf(SST_FILE_PRE , "%s/input/rtgssthr_grb_0.5.grib2", real2game_root_dir);
	char SST_FILE[strlen(SST_FILE_PRE) + 1];
	strcpy(SST_FILE, SST_FILE_PRE);
	
	ecc_file = fopen(SST_FILE, "r");
	handle = codes_handle_new_from_file(NULL, ecc_file, PRODUCT_GRIB, &err);
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
	fclose(ecc_file);
	
	// transforming the coordinates of the SST grid from degrees to radians
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_SST_POINTS; ++i)
	{
		latitudes_sst[i] = 2.0*M_PI*latitudes_sst[i]/360.0;
		longitudes_sst[i] = 2.0*M_PI*longitudes_sst[i]/360.0;
	}
    
    // Writing the observations to a netcdf file.
    char OUTPUT_FILE_PRE[200];
    sprintf(OUTPUT_FILE_PRE, "%s/input/obs_%s%s%s%s.nc", real2game_root_dir, year_string, month_string, day_string, hour_string);
	char OUTPUT_FILE[strlen(OUTPUT_FILE_PRE) + 1];
	strcpy(OUTPUT_FILE, OUTPUT_FILE_PRE);
    
    int ncid, h_dimid, v_dimid, sst_dimid, t_id, spec_hum_id, sp_id, z_id, z_surf_id, u_id, v_id, lat_sst_id, lon_sst_id, sst_id;
    int dim_vector[2];
    if ((retval = nc_create(OUTPUT_FILE, NC_CLOBBER, &ncid)))
        NCERR(retval);
    // Defining the dimensions.
    if ((retval = nc_def_dim(ncid, "h_index", NO_OF_POINTS_PER_LAYER_INPUT, &h_dimid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "v_index", NO_OF_LEVELS_INPUT, &v_dimid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "sst_index", NO_OF_SST_POINTS, &sst_dimid)))
        NCERR(retval);
    // Defining the variables."interpol_index", NO_OF_AVG_POINTS, &avg_dimid)))
        NCERR(retval);
    dim_vector[0] = h_dimid;
    dim_vector[1] = v_dimid;
    if ((retval = nc_def_var(NC_DOUBLE, "pressure_surface", NC_DOUBLE, 1, &h_dimid, &sp_id)))
        NCERR(retval);
    if ((retval = nc_def_var(NC_DOUBLE, "z_surface", NC_DOUBLE, 1, &h_dimid, &z_surf_id)))
        NCERR(retval);
    if ((retval = nc_def_var(NC_DOUBLE, "z_height", NC_DOUBLE, 2, dim_vector, &z_id)))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "temperature", NC_DOUBLE, 2, dim_vector, &t_id)))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "spec_humidity", NC_DOUBLE, 2, dim_vector, &spec_hum_id)))
        NCERR(retval);
	if ((retval = nc_def_var(ncid, "u_wind", NC_DOUBLE, 2, dim_vector, &u_id)))
	  	NCERR(retval);
	if ((retval = nc_def_var(ncid, "v_wind", NC_DOUBLE, 2, dim_vector, &v_id)))
	  	NCERR(retval);
	if ((retval = nc_def_var(ncid, "lat_sst", NC_DOUBLE, 1, &sst_dimid, &lat_sst_id)))
	  	NCERR(retval);
	if ((retval = nc_def_var(ncid, "lon_sst", NC_DOUBLE, 1, &sst_dimid, &lon_sst_id)))
	  	NCERR(retval);
	if ((retval = nc_def_var(ncid, "sst", NC_DOUBLE, 1, &sst_dimid, &sst_id)))
	  	NCERR(retval);
    if ((retval = nc_enddef(ncid)))
        NCERR(retval);
    // Setting the variables.
	if ((retval = nc_put_var_double(ncid, sp_id, &pressure_surface[0])))
	  	NCERR(retval);
	if ((retval = nc_put_var_double(ncid, z_id, &z_height_amsl[0][0])))
	  	NCERR(retval);
	if ((retval = nc_put_var_double(ncid, z_surf_id, &surface_height[0])))
	  	NCERR(retval);
	if ((retval = nc_put_var_double(ncid, t_id, &temperature[0][0])))
	  	NCERR(retval);
	if ((retval = nc_put_var_double(ncid, spec_hum_id, &spec_hum[0][0])))
	  	NCERR(retval);
	if ((retval = nc_put_var_double(ncid, u_id, &u_wind[0][0])))
	  	NCERR(retval);
	if ((retval = nc_put_var_double(ncid, v_id, &v_wind[0][0])))
	  	NCERR(retval);
	if ((retval = nc_put_var_double(ncid, lat_sst_id, &latitudes_sst[0])))
	  	NCERR(retval);
	if ((retval = nc_put_var_double(ncid, lon_sst_id, &longitudes_sst[0])))
	  	NCERR(retval);
	if ((retval = nc_put_var_double(ncid, sst_id, &sst[0])))
	  	NCERR(retval);
    if ((retval = nc_close(ncid)))
    	NCERR(retval);
    
    // Freeing the memory.
    free(z_height_amsl);
    free(surface_height);
	free(pressure_surface);
    free(u_wind);
    free(v_wind);
    free(temperature);
    free(spec_hum);
	free(sst);
	free(latitudes_sst);
	free(longitudes_sst);

	return 0;
}










