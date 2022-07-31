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
#include "../../header.h"
#include "../../game_properties.h"
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(1);}
#define NCCHECK(e) {if(e != 0) NCERR(e)}
#define ECCERR(e) {printf("Error: Eccodes failed with error code %d. See http://download.ecmwf.int/test-data/eccodes/html/group__errors.html for meaning of the error codes.\n", e); exit(1);}
#define ECCCHECK(e) {if(e != 0) ECCERR(e)}

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
	int oro_id = strtod(argv[7], NULL);
	int model_target_id = strtod(argv[8], NULL);
	int nlat = strtod(argv[9], NULL);
	int nlon = strtod(argv[10], NULL);
	int interpol_exp = strtod(argv[11], NULL);
    char lgame_grid[strlen(argv[12]) + 1];
    strcpy(lgame_grid, argv[12]);
    
    int no_of_points_per_layer_input_model;
	no_of_points_per_layer_input_model = 2949120;

	double *latitudes_input_model = malloc(no_of_points_per_layer_input_model*sizeof(double));
	double *longitudes_input_model = malloc(no_of_points_per_layer_input_model*sizeof(double));
	
	int err;
	codes_handle *handle = NULL;
	
	double one = 1.0;
	
	// Properties of the input model's grid.
	// latitudes of the grid
	char lat_obs_file_pre[200];
	sprintf(lat_obs_file_pre , "%s/interpolation_creator/icon_global_icosahedral_time-invariant_%s%s%s%s_CLAT.grib2",
	real2game_root_dir, year_string, month_string, day_string, hour_string);
	char lat_obs_file[strlen(lat_obs_file_pre) + 1];
	strcpy(lat_obs_file, lat_obs_file_pre);
	FILE *ECC_FILE;
	ECC_FILE = fopen(lat_obs_file, "r");
	handle = codes_handle_new_from_file(NULL, ECC_FILE, PRODUCT_GRIB, &err);
	if (err != 0) ECCERR(err);
	size_t no_of_points_per_layer_input_model_SIZE_T;
	no_of_points_per_layer_input_model_SIZE_T = (size_t) no_of_points_per_layer_input_model;
	ECCCHECK(codes_get_double_array(handle, "values", &latitudes_input_model[0], &no_of_points_per_layer_input_model_SIZE_T));
	codes_handle_delete(handle);
    fclose(ECC_FILE);
    
    // transforming the latitude coordinates of the input model from degrees to radians
    #pragma omp parallel for
    for (int i = 0; i < no_of_points_per_layer_input_model; ++i)
    {
    	latitudes_input_model[i] = 2.0*M_PI*latitudes_input_model[i]/360.0;
    }
    
    // longitudes of the grid
	char lon_obs_file_pre[200];
	sprintf(lon_obs_file_pre , "%s/interpolation_creator/icon_global_icosahedral_time-invariant_%s%s%s%s_CLON.grib2",
	real2game_root_dir, year_string, month_string, day_string, hour_string);
	char lon_obs_file[strlen(lon_obs_file_pre) + 1];
	strcpy(lon_obs_file, lon_obs_file_pre);
    
	ECC_FILE = fopen(lon_obs_file, "r");
	handle = codes_handle_new_from_file(NULL, ECC_FILE, PRODUCT_GRIB, &err);
	if (err != 0) ECCERR(err);
	no_of_points_per_layer_input_model_SIZE_T = (size_t) no_of_points_per_layer_input_model;
	ECCCHECK(codes_get_double_array(handle, "values", &longitudes_input_model[0], &no_of_points_per_layer_input_model_SIZE_T));
	codes_handle_delete(handle);
    fclose(ECC_FILE);
    
    // transforming the longitude coordinates of the input model from degrees to radians
    #pragma omp parallel for
    for (int i = 0; i < no_of_points_per_layer_input_model; ++i)
    {
    	longitudes_input_model[i] = 2.0*M_PI*longitudes_input_model[i]/360.0;
    }
    
    // GAME
    if (model_target_id == 0)
    {
		// reading the horizontal coordinates of the grid of GAME
		double *latitudes_game = malloc(N_SCALARS_H*sizeof(double));
		double *longitudes_game = malloc(N_SCALARS_H*sizeof(double));
		double *latitudes_game_wind = malloc(N_VECTORS_H*sizeof(double));
		double *longitudes_game_wind = malloc(N_VECTORS_H*sizeof(double));
		char geo_pro_file_pre[200];
		sprintf(geo_pro_file_pre, "%s/grid_generator/grids/RES%d_L%d_ORO%d.nc", model_home_dir, RES_ID, N_LAYERS, oro_id);
		char geo_pro_file[strlen(geo_pro_file_pre) + 1];
		strcpy(geo_pro_file, geo_pro_file_pre);
		printf("Grid file: %s\n", geo_pro_file);
		printf("Reading grid file of GAME ...\n");
		int ncid;
		NCCHECK(nc_open(geo_pro_file, NC_NOWRITE, &ncid));
		
		int latitudes_game_id, latitudes_game_wind_id, longitudes_game_id, longitudes_game_wind_id;
		NCCHECK(nc_inq_varid(ncid, "latitude_scalar", &latitudes_game_id));
		NCCHECK(nc_inq_varid(ncid, "longitude_scalar", &longitudes_game_id));
		NCCHECK(nc_inq_varid(ncid, "latitude_vector", &latitudes_game_wind_id));
		NCCHECK(nc_inq_varid(ncid, "longitude_vector", &longitudes_game_wind_id));
		NCCHECK(nc_get_var_double(ncid, latitudes_game_id, &latitudes_game[0]));
		NCCHECK(nc_get_var_double(ncid, longitudes_game_id, &longitudes_game[0]));
		NCCHECK(nc_get_var_double(ncid, latitudes_game_wind_id, &latitudes_game_wind[0]));
		NCCHECK(nc_get_var_double(ncid, longitudes_game_wind_id, &longitudes_game_wind[0]));
		NCCHECK(nc_close(ncid));
		printf("Grid file of GAME read.\n");
		
		// allocating memory for the result arrays
		int (*interpolation_indices_scalar)[N_AVG_POINTS] = malloc(sizeof(int[N_SCALARS_H][N_AVG_POINTS]));
		double (*interpolation_weights_scalar)[N_AVG_POINTS] = malloc(sizeof(double[N_SCALARS_H][N_AVG_POINTS]));
		int (*interpolation_indices_vector)[N_AVG_POINTS] = malloc(sizeof(int[N_VECTORS_H][N_AVG_POINTS]));
		double (*interpolation_weights_vector)[N_AVG_POINTS] = malloc(sizeof(double[N_VECTORS_H][N_AVG_POINTS]));
		
		// executing the actual interpolation
		printf("Calculating interpolation indices and weights ...\n");
		
		#pragma omp parallel for
		for (int i = 0; i < N_SCALARS_H; ++i)
		{
			double *distance_vector = malloc(no_of_points_per_layer_input_model*sizeof(double));
			for (int j = 0; j < no_of_points_per_layer_input_model; ++j)
			{
				distance_vector[j] = calculate_distance_h(&latitudes_game[i], &longitudes_game[i], &latitudes_input_model[j], &longitudes_input_model[j], &one);
			}
			double sum_of_weights = 0.0;
			for (int j = 0; j < N_AVG_POINTS; ++j)
			{
				interpolation_indices_scalar[i][j] = find_min_index(distance_vector, &no_of_points_per_layer_input_model);
				interpolation_weights_scalar[i][j] = 1.0/(pow(distance_vector[interpolation_indices_scalar[i][j]], interpol_exp));
				distance_vector[interpolation_indices_scalar[i][j]] = 2.0*M_PI;
				sum_of_weights += interpolation_weights_scalar[i][j];
			}
			free(distance_vector);
			for (int j = 0; j < N_AVG_POINTS; ++j)
			{
				interpolation_weights_scalar[i][j] = interpolation_weights_scalar[i][j]/sum_of_weights;
			}
		}
		
		#pragma omp parallel for
		for (int i = 0; i < N_VECTORS_H; ++i)
		{
			double *distance_vector = malloc(no_of_points_per_layer_input_model*sizeof(double));
			for (int j = 0; j < no_of_points_per_layer_input_model; ++j)
			{
				distance_vector[j] =  calculate_distance_h(&latitudes_game_wind[i], &longitudes_game_wind[i], &latitudes_input_model[j], &longitudes_input_model[j], &one);
			}
			double sum_of_weights = 0.0;
			for (int j = 0; j < N_AVG_POINTS; ++j)
			{
				interpolation_indices_vector[i][j] = find_min_index(distance_vector, &no_of_points_per_layer_input_model);
				interpolation_weights_vector[i][j] = 1.0/(pow(distance_vector[interpolation_indices_vector[i][j]], interpol_exp));
				distance_vector[interpolation_indices_vector[i][j]] = 2.0*M_PI;
				sum_of_weights += interpolation_weights_vector[i][j];
			}
			free(distance_vector);
			for (int j = 0; j < N_AVG_POINTS; ++j)
			{
				interpolation_weights_vector[i][j] = interpolation_weights_vector[i][j]/sum_of_weights;
			}
		}
		printf("Calculating interpolation indices and weights finished.\n");
		
		// freeing memory we do not need further
		free(latitudes_input_model);
		free(longitudes_input_model);
		free(latitudes_game);
		free(longitudes_game);
		free(latitudes_game_wind);
		free(longitudes_game_wind);
		
		// writing the result to a NetCDF file
		char output_file_pre[200];
		sprintf(output_file_pre, "%s/interpolation_files/icon-global2game%d.nc", real2game_root_dir, RES_ID);
		char output_file[strlen(output_file_pre) + 1];
		strcpy(output_file, output_file_pre);
		printf("Starting to write to output file ...\n");
		int interpolation_indices_scalar_id, interpolation_weights_scalar_id, interpolation_indices_vector_id, interpolation_weights_vector_id,
		scalar_dimid, vector_dimid, avg_dimid;
		int dim_vector[2];
		NCCHECK(nc_create(output_file, NC_CLOBBER, &ncid));
		NCCHECK(nc_def_dim(ncid, "scalar_index", N_SCALARS_H, &scalar_dimid));
		NCCHECK(nc_def_dim(ncid, "vector_index", N_VECTORS_H, &vector_dimid));
		NCCHECK(nc_def_dim(ncid, "interpol_index", N_AVG_POINTS, &avg_dimid));
		dim_vector[0] = scalar_dimid;
		dim_vector[1] = avg_dimid;
		NCCHECK(nc_def_var(ncid, "interpolation_indices_scalar", NC_INT, 2, dim_vector, &interpolation_indices_scalar_id));
		NCCHECK(nc_def_var(ncid, "interpolation_weights_scalar", NC_DOUBLE, 2, dim_vector, &interpolation_weights_scalar_id));
		dim_vector[0] = vector_dimid;
		NCCHECK(nc_def_var(ncid, "interpolation_indices_vector", NC_INT, 2, dim_vector, &interpolation_indices_vector_id));
		NCCHECK(nc_def_var(ncid, "interpolation_weights_vector", NC_DOUBLE, 2, dim_vector, &interpolation_weights_vector_id));
		NCCHECK(nc_enddef(ncid));
		NCCHECK(nc_put_var_int(ncid, interpolation_indices_scalar_id, &interpolation_indices_scalar[0][0]));
		NCCHECK(nc_put_var_double(ncid, interpolation_weights_scalar_id, &interpolation_weights_scalar[0][0]));
		NCCHECK(nc_put_var_int(ncid, interpolation_indices_vector_id, &interpolation_indices_vector[0][0]));
		NCCHECK(nc_put_var_double(ncid, interpolation_weights_vector_id, &interpolation_weights_vector[0][0]));
		NCCHECK(nc_close(ncid));
		    
		// freeing the memory
		free(interpolation_indices_scalar);
		free(interpolation_weights_scalar);
		free(interpolation_indices_vector);
		free(interpolation_weights_vector);
		    
    }
    
    // L-GAME
    if (model_target_id == 1)
    {
		// reading the horizontal coordinates of the grid of L-GAME
		double (*latitudes_lgame)[nlat] = malloc(sizeof(double[nlon][nlat]));
		double (*longitudes_lgame)[nlat] = malloc(sizeof(double[nlon][nlat]));
		double (*latitudes_lgame_wind_u)[nlat] = malloc(sizeof(double[nlon][nlat+1]));
		double (*longitudes_lgame_wind_u)[nlat] = malloc(sizeof(double[nlon][nlat+1]));
		double (*latitudes_lgame_wind_v)[nlat] = malloc(sizeof(double[nlon+1][nlat]));
		double (*longitudes_lgame_wind_v)[nlat] = malloc(sizeof(double[nlon+1][nlat]));
		char geo_pro_file_pre[200];
		sprintf(geo_pro_file_pre, "%s/grids/%s", model_home_dir, lgame_grid);
		char geo_pro_file[strlen(geo_pro_file_pre) + 1];
		strcpy(geo_pro_file, geo_pro_file_pre);
		printf("Grid file: %s\n", geo_pro_file);
		printf("Reading grid file of L-GAME ...\n");
		int ncid;
		NCCHECK(nc_open(geo_pro_file, NC_NOWRITE, &ncid));
		
		int latitudes_lgame_id, longitudes_lgame_id, longitudes_lgame_wind_u_id, latitudes_lgame_wind_u_id,
		longitudes_lgame_wind_v_id, latitudes_lgame_wind_v_id;
		NCCHECK(nc_inq_varid(ncid, "lat_geo", &latitudes_lgame_id));
		NCCHECK(nc_inq_varid(ncid, "lon_geo", &longitudes_lgame_id));
		NCCHECK(nc_inq_varid(ncid, "lat_geo_u", &latitudes_lgame_wind_u_id));
		NCCHECK(nc_inq_varid(ncid, "lon_geo_u", &longitudes_lgame_wind_u_id));
		NCCHECK(nc_inq_varid(ncid, "lat_geo_v", &latitudes_lgame_wind_v_id));
		NCCHECK(nc_inq_varid(ncid, "lon_geo_v", &longitudes_lgame_wind_v_id));
		NCCHECK(nc_get_var_double(ncid, latitudes_lgame_id, &latitudes_lgame[0][0]));
		NCCHECK(nc_get_var_double(ncid, longitudes_lgame_id, &longitudes_lgame[0][0]));
		NCCHECK(nc_get_var_double(ncid, latitudes_lgame_wind_u_id, &latitudes_lgame_wind_u[0][0]));
		NCCHECK(nc_get_var_double(ncid, longitudes_lgame_wind_u_id, &longitudes_lgame_wind_u[0][0]));
		NCCHECK(nc_get_var_double(ncid, latitudes_lgame_wind_v_id, &latitudes_lgame_wind_v[0][0]));
		NCCHECK(nc_get_var_double(ncid, longitudes_lgame_wind_v_id, &longitudes_lgame_wind_v[0][0]));
		NCCHECK(nc_close(ncid));
		printf("Grid file of L-GAME read.\n");
		
		// allocating memory for the result arrays
		int (*interpolation_indices_scalar)[N_AVG_POINTS] = malloc(sizeof(int[nlat*nlon][N_AVG_POINTS]));
		double (*interpolation_weights_scalar)[N_AVG_POINTS] = malloc(sizeof(double[nlat*nlon][N_AVG_POINTS]));
		int (*interpolation_indices_vector_u)[N_AVG_POINTS] = malloc(sizeof(int[nlat*(nlon+1)][N_AVG_POINTS]));
		double (*interpolation_weights_vector_u)[N_AVG_POINTS] = malloc(sizeof(double[nlat*(nlon+1)][N_AVG_POINTS]));
		int (*interpolation_indices_vector_v)[N_AVG_POINTS] = malloc(sizeof(int[(nlat+1)*nlon][N_AVG_POINTS]));
		double (*interpolation_weights_vector_v)[N_AVG_POINTS] = malloc(sizeof(double[(nlat+1)*nlon][N_AVG_POINTS]));
		
		// executing the actual interpolation
		printf("Calculating interpolation indices and weights ...\n");
		
		#pragma omp parallel for
		for (int i = 0; i < nlat; ++i)
		{
			for (int j = 0; j < nlon; ++j)
			{
				double *distance_vector = malloc(no_of_points_per_layer_input_model*sizeof(double));
				for (int k = 0; k < no_of_points_per_layer_input_model; ++k)
				{
					distance_vector[k] = calculate_distance_h(&latitudes_lgame[i][j], &longitudes_lgame[i][j], &latitudes_input_model[k], &longitudes_input_model[k], &one);
				}
				double sum_of_weights = 0.0;
				for (int k = 0; k < N_AVG_POINTS; ++k)
				{
					interpolation_indices_scalar[j*nlat + i][k] = find_min_index(distance_vector, &no_of_points_per_layer_input_model);
					interpolation_weights_scalar[j*nlat + i][k] = 1.0/(pow(distance_vector[interpolation_indices_scalar[j*nlat + i][k]], interpol_exp));
					distance_vector[interpolation_indices_scalar[j*nlat + i][k]] = 2.0*M_PI;
					sum_of_weights += interpolation_weights_scalar[j*nlat + i][k];
				}
				free(distance_vector);
				for (int k = 0; k < N_AVG_POINTS; ++k)
				{
					interpolation_weights_scalar[j*nlat + i][k] = interpolation_weights_scalar[j*nlat + i][k]/sum_of_weights;
				}
			}
		}
		
		#pragma omp parallel for
		for (int i = 0; i < nlat; ++i)
		{
			for (int j = 0; j < nlon+1; ++j)
			{
				double *distance_vector = malloc(no_of_points_per_layer_input_model*sizeof(double));
				for (int k = 0; k < no_of_points_per_layer_input_model; ++k)
				{
					distance_vector[k] = calculate_distance_h(&latitudes_lgame_wind_u[i][j], &longitudes_lgame_wind_u[i][j], &latitudes_input_model[k], &longitudes_input_model[k], &one);
				}
				double sum_of_weights = 0.0;
				for (int k = 0; k < N_AVG_POINTS; ++k)
				{
					interpolation_indices_vector_u[j*nlat + i][k] = find_min_index(distance_vector, &no_of_points_per_layer_input_model);
					interpolation_weights_vector_u[j*nlat + i][k] = 1.0/(pow(distance_vector[interpolation_indices_vector_u[j*nlat + i][k]], interpol_exp));
					distance_vector[interpolation_indices_vector_u[j*nlat + i][k]] = 2.0*M_PI;
					sum_of_weights += interpolation_weights_vector_u[j*nlat + i][k];
				}
				free(distance_vector);
				for (int k = 0; k < N_AVG_POINTS; ++k)
				{
					interpolation_weights_vector_u[j*nlat + i][k] = interpolation_weights_vector_u[j*nlat + i][k]/sum_of_weights;
				}
			}
		}
		
		#pragma omp parallel for
		for (int i = 0; i < nlat+1; ++i)
		{
			for (int j = 0; j < nlon; ++j)
			{
				double *distance_vector = malloc(no_of_points_per_layer_input_model*sizeof(double));
				for (int k = 0; k < no_of_points_per_layer_input_model; ++k)
				{
					distance_vector[k] = calculate_distance_h(&latitudes_lgame_wind_v[i][j], &longitudes_lgame_wind_v[i][j], &latitudes_input_model[k], &longitudes_input_model[k], &one);
				}
				double sum_of_weights = 0.0;
				for (int k = 0; k < N_AVG_POINTS; ++k)
				{
					interpolation_indices_vector_v[j*(nlat+1) + i][k] = find_min_index(distance_vector, &no_of_points_per_layer_input_model);
					interpolation_weights_vector_v[j*(nlat+1) + i][k] = 1.0/(pow(distance_vector[interpolation_indices_vector_v[j*(nlat+1) + i][k]], interpol_exp));
					distance_vector[interpolation_indices_vector_v[j*(nlat+1) + i][k]] = 2.0*M_PI;
					sum_of_weights += interpolation_weights_vector_v[j*(nlat+1) + i][k];
				}
				free(distance_vector);
				for (int k = 0; k < N_AVG_POINTS; ++k)
				{
					interpolation_weights_vector_v[j*(nlat+1) + i][k] = interpolation_weights_vector_v[j*(nlat+1) + i][k]/sum_of_weights;
				}
			}
		}
		printf("Calculating interpolation indices and weights finished.\n");
		
		// freeing memory we do not need further
		free(latitudes_input_model);
		free(longitudes_input_model);
		free(latitudes_lgame);
		free(longitudes_lgame);
		free(latitudes_lgame_wind_u);
		free(longitudes_lgame_wind_u);
		free(latitudes_lgame_wind_v);
		free(longitudes_lgame_wind_v);
		
		// writing the result to a NetCDF file
		char output_file_pre[200];
		sprintf(output_file_pre, "%s/interpolation_files/icon-d22lgame_%s", real2game_root_dir, lgame_grid);
		char output_file[strlen(output_file_pre) + 1];
		strcpy(output_file, output_file_pre);
		printf("Starting to write to output file ...\n");
		int interpolation_indices_scalar_id, interpolation_weights_scalar_id, interpolation_indices_vector_u_id, interpolation_weights_vector_u_id,
        interpolation_indices_vector_v_id, interpolation_weights_vector_v_id, scalar_dimid, u_dimid, v_dimid, avg_dimid;
		int dim_vector[2];
		NCCHECK(nc_create(output_file, NC_CLOBBER, &ncid));
		NCCHECK(nc_def_dim(ncid, "scalar_index", nlat*nlon, &scalar_dimid));
		NCCHECK(nc_def_dim(ncid, "u_index", nlat*(nlon+1), &u_dimid));
		NCCHECK(nc_def_dim(ncid, "v_index", (nlat+1)*nlon, &v_dimid));
		NCCHECK(nc_def_dim(ncid, "interpol_index", N_AVG_POINTS, &avg_dimid));
		dim_vector[0] = scalar_dimid;
		dim_vector[1] = avg_dimid;
		NCCHECK(nc_def_var(ncid, "interpolation_indices_scalar", NC_INT, 2, dim_vector, &interpolation_indices_scalar_id));
		NCCHECK(nc_def_var(ncid, "interpolation_weights_scalar", NC_DOUBLE, 2, dim_vector, &interpolation_weights_scalar_id));
		dim_vector[0] = u_dimid;
		NCCHECK(nc_def_var(ncid, "interpolation_indices_vector_u", NC_INT, 2, dim_vector, &interpolation_indices_vector_u_id));
		NCCHECK(nc_def_var(ncid, "interpolation_weights_vector_u", NC_DOUBLE, 2, dim_vector, &interpolation_weights_vector_u_id));
		dim_vector[0] = v_dimid;
		NCCHECK(nc_def_var(ncid, "interpolation_indices_vector_v", NC_INT, 2, dim_vector, &interpolation_indices_vector_v_id));
		NCCHECK(nc_def_var(ncid, "interpolation_weights_vector_v", NC_DOUBLE, 2, dim_vector, &interpolation_weights_vector_v_id));
		NCCHECK(nc_enddef(ncid));
		NCCHECK(nc_put_var_int(ncid, interpolation_indices_scalar_id, &interpolation_indices_scalar[0][0]));
		NCCHECK(nc_put_var_double(ncid, interpolation_weights_scalar_id, &interpolation_weights_scalar[0][0]));
		NCCHECK(nc_put_var_int(ncid, interpolation_indices_vector_u_id, &interpolation_indices_vector_u[0][0]));
		NCCHECK(nc_put_var_double(ncid, interpolation_weights_vector_u_id, &interpolation_weights_vector_u[0][0]));
		NCCHECK(nc_put_var_int(ncid, interpolation_indices_vector_v_id, &interpolation_indices_vector_v[0][0]));
		NCCHECK(nc_put_var_double(ncid, interpolation_weights_vector_v_id, &interpolation_weights_vector_v[0][0]));
		NCCHECK(nc_close(ncid));
	
		// freeing the memory
		free(interpolation_indices_scalar);
		free(interpolation_weights_scalar);
		free(interpolation_indices_vector_u);
		free(interpolation_weights_vector_u);
		free(interpolation_indices_vector_v);
		free(interpolation_weights_vector_v);
	
    }
    
    printf("Finished.\n");
	
	return 0;
}







