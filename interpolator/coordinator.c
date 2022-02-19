/*
This source file is part of real2GAME, which is released under the MIT license.
Github repository: https://github.com/OpenNWP/real2GAME
*/

/*
This file coordinates the data interpolation process.
*/

#include <stdlib.h>
#include "../game_properties.h"
#include "../header.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <netcdf.h>
#include <geos95.h>
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}
#define EPSILON 1e-4
#define SCALE_HEIGHT 8000.0
#define P_0 100000
#define R_D 287.057811
#define C_D_P 1005.0

int main(int argc, char *argv[])
{
    char year_string[strlen(argv[1]) + 1];
    strcpy(year_string, argv[1]);
    char month_string[strlen(argv[2]) + 1];
    strcpy(month_string, argv[2]);
    char day_string[strlen(argv[3]) + 1];
    strcpy(day_string, argv[3]);
    char hour_string[strlen(argv[4]) + 1];
    strcpy(hour_string, argv[4]);
    char model_home_dir[strlen(argv[5]) + 1];
    strcpy(model_home_dir, argv[5]);
	int ORO_ID;
    ORO_ID = strtod(argv[6], NULL);
    char BACKGROUND_STATE_FILE[strlen(argv[7]) + 1];
    strcpy(BACKGROUND_STATE_FILE, argv[7]);
    char real2game_root_dir[strlen(argv[8]) + 1];
    strcpy(real2game_root_dir, argv[8]);
	printf("Background state file: %s\n", BACKGROUND_STATE_FILE);
    
    // Allocating memory for the grid properties.
    double *latitudes_game = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *longitudes_game = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *z_coords_game = malloc(NO_OF_SCALARS*sizeof(double));
    double *latitudes_game_wind = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *longitudes_game_wind = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *directions = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *z_coords_game_wind = malloc(NO_OF_VECTORS*sizeof(double));
    double *gravity_potential_game = malloc(NO_OF_SCALARS*sizeof(double));
    // Reading the grid properties.
    int ncid_grid, retval;
    char GEO_PROP_FILE_PRE[200];
    sprintf(GEO_PROP_FILE_PRE, "%s/grid_generator/grids/RES%d_L%d_ORO%d.nc", model_home_dir, RES_ID, NO_OF_LAYERS, ORO_ID);
    char GEO_PROP_FILE[strlen(GEO_PROP_FILE_PRE) + 1];
    strcpy(GEO_PROP_FILE, GEO_PROP_FILE_PRE);
	printf("Grid file: %s\n", GEO_PROP_FILE);
	printf("Reading grid file ...\n");
    if ((retval = nc_open(GEO_PROP_FILE, NC_NOWRITE, &ncid_grid)))
        NCERR(retval);
    int latitudes_game_id, longitudes_game_id, z_coords_game_id, z_coords_game_wind_id, gravity_potential_game_id, adjacent_vector_indices_h_id, directions_id;
    if ((retval = nc_inq_varid(ncid_grid, "latitude_scalar", &latitudes_game_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "longitude_scalar", &longitudes_game_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "z_scalar", &z_coords_game_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "direction", &directions_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "z_vector", &z_coords_game_wind_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "gravity_potential", &gravity_potential_game_id)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, latitudes_game_id, &latitudes_game[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, longitudes_game_id, &longitudes_game[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, z_coords_game_id, &z_coords_game[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, directions_id, &directions[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, z_coords_game_wind_id, &z_coords_game_wind[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, gravity_potential_game_id, &gravity_potential_game[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid_grid)))
        NCERR(retval);
	printf("Grid file read.\n");
	
    char output_file_pre[200];
    sprintf(output_file_pre, "%s/nwp_init/%s%s%s%s.nc", model_home_dir, year_string, month_string, day_string, hour_string);
    char output_file[strlen(output_file_pre) + 1];
    strcpy(output_file, output_file_pre);
    
    // These are the arrays of the background state.
    double *densities_background = malloc(6*NO_OF_SCALARS*sizeof(double));
    double *temperatures_background = malloc(5*NO_OF_SCALARS*sizeof(double));
    double *wind_background = malloc(NO_OF_VECTORS*sizeof(double));
    double *tke = malloc(NO_OF_SCALARS*sizeof(double));
    double *t_soil = malloc(NO_OF_SOIL_LAYERS*NO_OF_SCALARS_H*sizeof(double));
    
    // Reading the background state.
	printf("Reading background state ...\n");
    int ncid;
    if ((retval = nc_open(BACKGROUND_STATE_FILE, NC_NOWRITE, &ncid)))
        NCERR(retval);
    int densities_background_id, temperatures_background_id, wind_background_id, tke_avail, tke_id, t_soil_avail, t_soil_id;
    if ((retval = nc_inq_varid(ncid, "densities", &densities_background_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "temperatures", &temperatures_background_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "wind", &wind_background_id)))
        NCERR(retval);
    tke_avail = 0;
    if (nc_inq_varid(ncid, "tke", &tke_id) == 0)
    {
    	tke_avail = 1;
		if ((retval = nc_inq_varid(ncid, "tke", &tke_id)))
		    NCERR(retval);
		printf("TKE found in background state file.\n");
    }
    else
    {
    	printf("TKE not found in background state file.\n");
    }
    t_soil_avail = 0;
    if (nc_inq_varid(ncid, "t_soil", &t_soil_id) == 0)
    {
    	t_soil_avail = 1;
		if ((retval = nc_inq_varid(ncid, "t_soil", &t_soil_id)))
		    NCERR(retval);
		printf("Soil temperature found in background state file.\n");
    }
    else
    {
    	printf("Soil temperature not found in background state file.\n");
    }
    if ((retval = nc_get_var_double(ncid, densities_background_id, &densities_background[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, temperatures_background_id, &temperatures_background[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, wind_background_id, &wind_background[0])))
        NCERR(retval);
    if (tke_avail == 1)
    {
		if ((retval = nc_get_var_double(ncid, tke_id, &tke[0])))
		    NCERR(retval);
    }
    if (t_soil_avail == 1)
    {
		if ((retval = nc_get_var_double(ncid, t_soil_id, &t_soil[0])))
		    NCERR(retval);
    }
    if ((retval = nc_close(ncid)))
        NCERR(retval);
	printf("Background state read.\n");	
	
	// Allocating the memory for the observations.
	double (*z_coords_input_model)[NO_OF_LEVELS_INPUT] = malloc(sizeof(double[NO_OF_POINTS_PER_LAYER_INPUT][NO_OF_LEVELS_INPUT]));
	double (*temperature_in)[NO_OF_LEVELS_INPUT] = malloc(sizeof(double[NO_OF_POINTS_PER_LAYER_INPUT][NO_OF_LEVELS_INPUT]));
    
    char input_file_pre[200];
    sprintf(input_file_pre, "%s/input/obs_%s%s%s%s.nc", real2game_root_dir, year_string, month_string, day_string, hour_string);
    char input_file[strlen(input_file_pre) + 1];
    strcpy(input_file, input_file_pre);
	printf("Observations file: %s\n", input_file);
    
    // Reading the observations.
	printf("Reading observations ...\n");
    if ((retval = nc_open(input_file, NC_NOWRITE, &ncid)))
        NCERR(retval);
    int z_coords_id, t_in_id;
    // Defining the variables.
    if ((retval = nc_inq_varid(ncid, "z_height", &z_coords_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "temperature", &t_in_id)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, z_coords_id, &z_coords_input_model[0][0])))
        NCERR(retval);  
    if ((retval = nc_get_var_double(ncid, t_in_id, &temperature_in[0][0])))
        NCERR(retval);
    if ((retval = nc_close(ncid)))
    	NCERR(retval);
	printf("Observations read.\n");
	
	// Begin of the actual interpolation.
    
    /*
    DRY THERMODYNAMIC STATE INTERPOLATION
    -------------------------------------
    */
	
	printf("Starting the dry interpolation ...\n");
	
	// dry data interpolation is finished at this point
    
    // These are the arrays for the result of the interpolation process.
    double *pressure_surface_out = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *density_dry = malloc(NO_OF_SCALARS*sizeof(double));
    double *wind = malloc(NO_OF_VECTORS*sizeof(double));
    double *exner = malloc(NO_OF_SCALARS*sizeof(double));
    double *temperature = malloc(NO_OF_SCALARS*sizeof(double));
    
    // density is determined out of the hydrostatic equation
    int layer_index, h_index;
    double b, c;
    for (int i = NO_OF_SCALARS - 1; i >= 0; --i)
    {
    	layer_index = i/NO_OF_SCALARS_H;
    	h_index = i - layer_index*NO_OF_SCALARS_H;
    	if (layer_index == NO_OF_LAYERS - 1)
    	{
        	density_dry[i] = pressure_surface_out[h_index]/(R_D*temperature[i]);
        	exner[i] = pow((density_dry[i]*R_D*temperature[i])/P_0, R_D/C_D_P);
        }
        else
        {
			// solving a quadratic equation for the Exner pressure
			b = -0.5*exner[i + NO_OF_SCALARS_H]/temperature[i + NO_OF_SCALARS_H]
			*(temperature[i] - temperature[i + NO_OF_SCALARS_H]
			+ 2.0/C_D_P*(gravity_potential_game[i] - gravity_potential_game[i + NO_OF_SCALARS_H]));
			c = pow(exner[i + NO_OF_SCALARS_H], 2.0)*temperature[i]/temperature[i + NO_OF_SCALARS_H];
			exner[i] = b + pow((pow(b, 2.0) + c), 0.5);
        	density_dry[i] = P_0*pow(exner[i], C_D_P/R_D)/(R_D*temperature[i]);
        }
    }
	free(gravity_potential_game);
    free(exner);
    // end of the interpolation of the dry thermodynamic state
	printf("Dry interpolation completed.\n");
    
    /*
    MOISTURE INTERPOLATION
    ----------------------
    */
    
	printf("Starting the moist interpolation ...\n");
	
	
	
	printf("Moist interpolation completed.\n");
	
	/*
	WIND INTERPOLATION
	------------------
	*/
	
	printf("Starting the wind interpolation ...\n");
	
	
	
	printf("Wind interpolation completed.\n");
	
	/*
	INTERPOLATION OF THE SST
	------------------------
	*/
	printf("Interpolating the SST to the model grid ...\n");
	double *sst_out = malloc(NO_OF_SCALARS_H*sizeof(double));
	int min_index;
	#pragma omp parallel for private(min_index)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
    	double *distance_vector = malloc(NO_OF_SST_POINTS*sizeof(double));
    	for (int j = 0; j < NO_OF_SST_POINTS; ++j)
    	{
    		distance_vector[j] = calculate_distance_h(latitudes_sst[j], longitudes_sst[j], latitudes_game[i], longitudes_game[i], 1);
    	}
		min_index = find_min_index(distance_vector, NO_OF_SST_POINTS);
		sst_out[i] = sst_in[min_index];
		free(distance_vector);
	}
	free(latitudes_sst);
	free(longitudes_sst);
	free(latitudes_game);
	free(longitudes_game);
	free(sst_in);
	
	printf("Interpolation of the SST completed.\n");

	/*
	PREPARING THE OUTPUT
	--------------------
	*/
	
	// individual condensate temperatures are for higher resolutions, not yet implemented
	// clouds and precipitation are set equal to the background state
	double *densities = malloc(6*NO_OF_SCALARS*sizeof(double));
	double *temperatures = malloc(5*NO_OF_SCALARS*sizeof(double));
    #pragma omp parallel for
	for (int i = 0; i < NO_OF_SCALARS; ++i)
	{
		// setting the mass densities of the result
		// condensate densities are not assimilated
		densities[i] = densities_background[i];
		densities[NO_OF_SCALARS + i] = densities_background[NO_OF_SCALARS + i];
		densities[2*NO_OF_SCALARS + i] = densities_background[2*NO_OF_SCALARS + i];
		densities[3*NO_OF_SCALARS + i] = densities_background[3*NO_OF_SCALARS + i];
		densities[4*NO_OF_SCALARS + i] = density_dry[i];
		densities[5*NO_OF_SCALARS + i] = model_vector_moist[i]/(1 - model_vector_moist[i])*density_dry[i];
		if (densities[5*NO_OF_SCALARS + i] < 0)
		{
			densities[5*NO_OF_SCALARS + i] = 0;
		}
		// setting the temperatures of the result
		// assuming an LTE (local thermodynamic equilibrium)
		temperatures[i] = temperature[i];
		temperatures[NO_OF_SCALARS + i] = temperature[i];
		temperatures[2*NO_OF_SCALARS + i] = temperature[i];
		temperatures[3*NO_OF_SCALARS + i] = temperature[i];
		temperatures[4*NO_OF_SCALARS + i] = temperature[i];
    }
    free(temperature);
    free(density_dry);
	free(model_vector_moist);
	free(densities_background);
    
    // writing the result of the wind data interpolation to the resulting wind field
    #pragma omp parallel for private(layer_index, h_index)
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
    	layer_index = i/NO_OF_VECTORS_PER_LAYER;
    	h_index = i - layer_index*NO_OF_VECTORS_PER_LAYER;
    	if (h_index < NO_OF_SCALARS_H)
    	{
    		wind[i] = wind_background[i];
    	}
    	else
    	{
    		wind[i] = wind_h_out[layer_index*NO_OF_VECTORS_H + h_index - NO_OF_SCALARS_H];
    	}
    }
    free(abs_humidity_game);
    free(wind_background);
    
    /*
    writing the result to a netcdf file
    -----------------------------------
    */
    
    printf("Output file: %s\n", output_file);
    printf("Writing result to output file ...\n");
    int densities_dimid, temperatures_dimid, scalar_dimid, vector_dimid, scalar_h_dimid, single_double_dimid,
    densities_id, temperatures_id, wind_id, sst_id, soil_dimid;
    if ((retval = nc_create(output_file, NC_CLOBBER, &ncid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "densities_index", 6*NO_OF_SCALARS, &densities_dimid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "temperatures_index", 5*NO_OF_SCALARS, &temperatures_dimid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "vector_index", NO_OF_VECTORS, &vector_dimid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "scalar_index", NO_OF_SCALARS, &scalar_dimid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "soil_index", NO_OF_SOIL_LAYERS*NO_OF_SCALARS_H, &soil_dimid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "scalar_h_index", NO_OF_SCALARS_H, &scalar_h_dimid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "single_double_dimid_index", 1, &single_double_dimid)))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "densities", NC_DOUBLE, 1, &densities_dimid, &densities_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, densities_id, "units", strlen("kg/m^3"), "kg/m^3")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "temperatures", NC_DOUBLE, 1, &temperatures_dimid, &temperatures_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, temperatures_id, "units", strlen("K"), "K")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "wind", NC_DOUBLE, 1, &vector_dimid, &wind_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, wind_id, "units", strlen("m/s"), "m/s")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "sst", NC_DOUBLE, 1, &scalar_h_dimid, &sst_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, sst_id, "units", strlen("K"), "K")))
        NCERR(retval);
    if (tke_avail == 1)
    {
		if ((retval = nc_def_var(ncid, "tke", NC_DOUBLE, 1, &scalar_dimid, &tke_id)))
		    NCERR(retval);
		if ((retval = nc_put_att_text(ncid, tke_id, "units", strlen("J/kg"), "J/kg")))
		    NCERR(retval);
    }
    if (t_soil_avail == 1)
    {
		if ((retval = nc_def_var(ncid, "t_soil", NC_DOUBLE, 1, &soil_dimid, &t_soil_id)))
		    NCERR(retval);
		if ((retval = nc_put_att_text(ncid, t_soil_id, "units", strlen("K"), "K")))
		    NCERR(retval);
    }
    if ((retval = nc_enddef(ncid)))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, densities_id, &densities[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, temperatures_id, &temperatures[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, wind_id, &wind[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, sst_id, &sst_out[0])))
        NCERR(retval);
    if (tke_avail == 1)
    {
		if ((retval = nc_put_var_double(ncid, tke_id, &tke[0])))
		    NCERR(retval);
    }
    if (t_soil_avail == 1)
    {
		if ((retval = nc_put_var_double(ncid, t_soil_id, &t_soil[0])))
		    NCERR(retval);
    }
    if ((retval = nc_close(ncid)))
    	NCERR(retval);
    printf("Result successfully written.\n");
    
    // freeing the stil occupied memory
	free(densities);
	free(temperatures);
	free(wind);
	free(sst_out);
	free(tke);
	free(t_soil);
	
	// that's it
    return 0;
}

















