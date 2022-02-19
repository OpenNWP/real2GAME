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

int main(int argc, char *argv[])
{
	double C_D_P = 1005.0;
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
    double *latitudes_model = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *longitudes_model = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *z_coords_model = malloc(NO_OF_SCALARS*sizeof(double));
    double *latitudes_model_wind = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *longitudes_model_wind = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *directions = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *z_coords_model_wind = malloc(NO_OF_VECTORS*sizeof(double));
    double *gravity_potential_model = malloc(NO_OF_SCALARS*sizeof(double));
    double *normal_distance = malloc(NO_OF_VECTORS*sizeof(double));
    int *from_index = malloc(NO_OF_VECTORS_H*sizeof(int));
    int *to_index = malloc(NO_OF_VECTORS_H*sizeof(int));
    int *adjacent_vector_indices_h = malloc(6*NO_OF_SCALARS_H*sizeof(int));
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
    int latitudes_model_id, latitudes_model_wind_id, longitudes_model_id, longitudes_model_wind_id, z_coords_model_id,
    z_coords_model_wind_id, gravity_potential_model_id, normal_distance_id, from_index_id, to_index_id, adjacent_vector_indices_h_id, directions_id;
    if ((retval = nc_inq_varid(ncid_grid, "latitude_scalar", &latitudes_model_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "longitude_scalar", &longitudes_model_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "z_scalar", &z_coords_model_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "latitude_vector", &latitudes_model_wind_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "longitude_vector", &longitudes_model_wind_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "direction", &directions_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "z_vector", &z_coords_model_wind_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "gravity_potential", &gravity_potential_model_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "normal_distance", &normal_distance_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "from_index", &from_index_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "to_index", &to_index_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "adjacent_vector_indices_h", &adjacent_vector_indices_h_id)))
        NCERR(retval);
    if ((retval = nc_get_var_int(ncid_grid, from_index_id, &from_index[0])))
        NCERR(retval);
    if ((retval = nc_get_var_int(ncid_grid, to_index_id, &to_index[0])))
        NCERR(retval);
    if ((retval = nc_get_var_int(ncid_grid, adjacent_vector_indices_h_id, &adjacent_vector_indices_h[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, latitudes_model_id, &latitudes_model[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, longitudes_model_id, &longitudes_model[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, z_coords_model_id, &z_coords_model[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, latitudes_model_wind_id, &latitudes_model_wind[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, longitudes_model_wind_id, &longitudes_model_wind[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, directions_id, &directions[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, z_coords_model_wind_id, &z_coords_model_wind[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, gravity_potential_model_id, &gravity_potential_model[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, normal_distance_id, &normal_distance[0])))
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

	// saving the relevant part of the background state in one array
	double *background_dry = malloc(NO_OF_MODEL_DOFS_DRY*sizeof(double));
    #pragma omp parallel for
	for (int i = 0; i < NO_OF_MODEL_DOFS_DRY; ++i)
	{
		// temperature of the gas phase of the background state
		if (i < NO_OF_SCALARS)
		{
			background_dry[i] = temperatures_background[4*NO_OF_SCALARS + i];
		}
		// dry air density in the lowest layer of the background state
		else
		{
			background_dry[i] = densities_background[4*NO_OF_SCALARS + i - NO_OF_SCALARS_H];
		}
	}
    free(temperatures_background);	
	
	// Allocating the memory for the observations.
	double *latitude_vector_obs = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(double));
	double *longitude_vector_obs = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(double));
	double *z_coords_obs = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(double));
	double *observations_vector = malloc(NO_OF_CHOSEN_OBSERVATIONS*sizeof(double));
    
    char OBSERVATIONS_FILE_PRE[200];
    sprintf(OBSERVATIONS_FILE_PRE, "%s/input/obs_%s%s%s%s.nc", real2game_root_dir, year_string, month_string, day_string, hour_string);
    char OBSERVATIONS_FILE[strlen(OBSERVATIONS_FILE_PRE) + 1];
    strcpy(OBSERVATIONS_FILE, OBSERVATIONS_FILE_PRE);
	printf("Observations file: %s\n", OBSERVATIONS_FILE);
    
    // Reading the observations.
	printf("Reading observations ...\n");
    if ((retval = nc_open(OBSERVATIONS_FILE, NC_NOWRITE, &ncid)))
        NCERR(retval);
    int latitude_obs_id, longitude_obs_id, z_coords_id, obervations_id;
    // Defining the variables.
    if ((retval = nc_inq_varid(ncid, "latitude_vector", &latitude_obs_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "longitude_vector", &longitude_obs_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "z_coords_obs", &z_coords_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "observations_vector", &obervations_id)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, latitude_obs_id, &latitude_vector_obs[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, longitude_obs_id, &longitude_vector_obs[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, z_coords_id, &z_coords_obs[0])))
        NCERR(retval);  
    if ((retval = nc_get_var_double(ncid, obervations_id, &observations_vector[0])))
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
    
	// setting up the observations operator
	double *interpolated_model_dry = malloc(NO_OF_CHOSEN_OBSERVATIONS_DRY*sizeof(double));
    
    free(z_coords_model);
	
	
	// data interpolation is finished at this point
    
    // These are the arrays for the result of the interpolation process.
    double *density_dry = malloc(NO_OF_SCALARS*sizeof(double));
    double *wind = malloc(NO_OF_VECTORS*sizeof(double));
    double *exner = malloc(NO_OF_SCALARS*sizeof(double));
    
    // density is determined out of the hydrostatic equation
    double b, c;
    for (int i = NO_OF_SCALARS - 1; i >= 0; --i)
    {
    	layer_index = i/NO_OF_SCALARS_H;
    	h_index = i - layer_index*NO_OF_SCALARS_H;
    	// at the lowest layer the density is part of the model_vector_dry
    	if (layer_index == NO_OF_LAYERS - 1)
    	{
        	density_dry[i] = model_vector_dry[NO_OF_SCALARS + h_index];
        	exner[i] = pow((density_dry[i]*R_D*model_vector_dry[i])/P_0, R_D/C_D_P);
        }
        else
        {
			// solving a quadratic equation for the Exner pressure
			b = -0.5*exner[i + NO_OF_SCALARS_H]/model_vector_dry[i + NO_OF_SCALARS_H]
			*(model_vector_dry[i] - model_vector_dry[i + NO_OF_SCALARS_H]
			+ 2.0/C_D_P*(gravity_potential_model[i] - gravity_potential_model[i + NO_OF_SCALARS_H]));
			c = pow(exner[i + NO_OF_SCALARS_H], 2.0)*model_vector_dry[i]/model_vector_dry[i + NO_OF_SCALARS_H];
			exner[i] = b + pow((pow(b, 2.0) + c), 0.5);
        	density_dry[i] = P_0*pow(exner[i], C_D_P/R_D)/(R_D*model_vector_dry[i]);
        }
    }
	free(gravity_potential_model);
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
	double *sst = malloc(NO_OF_SCALARS_H*sizeof(double));
	int min_index;
	#pragma omp parallel for private(min_index)
	for (int i = 0; i < NO_OF_SCALARS_H; ++i)
	{
    	double *distance_vector = malloc(NO_OF_SST_POINTS*sizeof(double));
    	for (int j = 0; j < NO_OF_SST_POINTS; ++j)
    	{
    		distance_vector[j] = calculate_distance_h(latitude_vector_obs[NO_OF_CHOSEN_OBSERVATIONS - NO_OF_SST_POINTS + j],
    		longitude_vector_obs[NO_OF_CHOSEN_OBSERVATIONS - NO_OF_SST_POINTS + j], latitudes_model[i], longitudes_model[i], 1);
    	}
		min_index = find_min_index(distance_vector, NO_OF_SST_POINTS);
		sst[i] = observations_vector[NO_OF_CHOSEN_OBSERVATIONS - NO_OF_SST_POINTS + min_index];
		free(distance_vector);
	}
	free(latitude_vector_obs);
	free(longitude_vector_obs);
	free(observations_vector);
	
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
		temperatures[i] = model_vector_dry[i];
		temperatures[NO_OF_SCALARS + i] = model_vector_dry[i];
		temperatures[2*NO_OF_SCALARS + i] = model_vector_dry[i];
		temperatures[3*NO_OF_SCALARS + i] = model_vector_dry[i];
		temperatures[4*NO_OF_SCALARS + i] = model_vector_dry[i];
    }
    free(model_vector_dry);
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
    		wind[i] = model_vector_wind[layer_index*NO_OF_VECTORS_H + h_index - NO_OF_SCALARS_H];
    	}
    }
    free(model_vector_wind);
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
    if ((retval = nc_put_var_double(ncid, sst_id, &sst[0])))
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
	free(sst);
	free(tke);
	free(t_soil);
	
	// that's it
    return 0;
}

















