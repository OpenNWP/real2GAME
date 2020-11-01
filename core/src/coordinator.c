/*
This source file is part of ndvar, which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/ndvar
*/

#include <stdlib.h>
#include "enum.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <netcdf.h>
#include "geos95.h"
#include "atmostracers.h"
#define NCERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}
#define N_A (6.0221409e23)
#define K_B (1.380649e-23)
#define M_D 0.028964420
#define R (N_A*K_B)
#define R_D (R/M_D)
#define P_0 100000.0
#define OMEGA (7.292115e-5)
#define C_D_P 1005.0

// constants specifying the grid
const double TOA = 30000;

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
    char *hour = malloc((len + 1)*sizeof(char));
    strcpy(hour, argv[4]);
    len = strlen(argv[5]);
    char *model_home_dir = malloc((len + 1)*sizeof(char));
    strcpy(model_home_dir, argv[5]);
	int ORO_ID;
    ORO_ID = strtod(argv[6], NULL);
    len = strlen(argv[7]);
    char *BACKGROUND_STATE_FILE = malloc((len + 1)*sizeof(char));
    strcpy(BACKGROUND_STATE_FILE, argv[7]);
    
    // The grid properties.
    double *direction = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *latitude_scalar = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *longitude_scalar = malloc(NO_OF_SCALARS_H*sizeof(double));
    double *latitude_vector = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *longitude_vector = malloc(NO_OF_VECTORS_H*sizeof(double));
    double *z_scalar = malloc(NO_OF_SCALARS*sizeof(double));
    double *z_vector = malloc(NO_OF_VECTORS*sizeof(double));
    double *gravity_potential = malloc(NO_OF_VECTORS*sizeof(double));
    int ncid_grid, retval;
    int GEO_PROP_FILE_LENGTH = 100;
    char *GEO_PROP_FILE_PRE = malloc((GEO_PROP_FILE_LENGTH + 1)*sizeof(char));
    sprintf(GEO_PROP_FILE_PRE, "%s/grids/B%dL%dT%d_O%d_OL%d_SCVT.nc", model_home_dir, RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
    GEO_PROP_FILE_LENGTH = strlen(GEO_PROP_FILE_PRE);
    free(GEO_PROP_FILE_PRE);
    char *GEO_PROP_FILE = malloc((GEO_PROP_FILE_LENGTH + 1)*sizeof(char));
    sprintf(GEO_PROP_FILE, "%s/grids/B%dL%dT%d_O%d_OL%d_SCVT.nc", model_home_dir, RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
    if ((retval = nc_open(GEO_PROP_FILE, NC_NOWRITE, &ncid_grid)))
        NCERR(retval);
    free(GEO_PROP_FILE);
    int direction_id, latitude_scalar_id, longitude_scalar_id, latitude_vector_id, longitude_vector_id, z_scalar_id, z_vector_id, gravity_potential_id;
    if ((retval = nc_inq_varid(ncid_grid, "direction", &direction_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "latitude_scalar", &latitude_scalar_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "longitude_scalar", &longitude_scalar_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "latitude_vector", &latitude_vector_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "longitude_vector", &longitude_vector_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "z_scalar", &z_scalar_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "z_vector", &z_vector_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid_grid, "gravity_potential", &gravity_potential_id)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, direction_id, &direction[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, latitude_scalar_id, &latitude_scalar[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, longitude_scalar_id, &longitude_scalar[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, latitude_vector_id, &latitude_vector[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, longitude_vector_id, &longitude_vector[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, z_scalar_id, &z_scalar[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, z_vector_id, &z_vector[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid_grid, gravity_potential_id, &gravity_potential[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid_grid)))
        NCERR(retval);
    int OUTPUT_FILE_LENGTH = 100;
    char *OUTPUT_FILE_PRE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE_PRE, "%s/input/%s%s%s%s_nwp_B%dL%dT%d_O%d_OL%d_SCVT.nc", model_home_dir, year_string, month_string, day_string, hour, RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
    OUTPUT_FILE_LENGTH = strlen(OUTPUT_FILE_PRE);
    free(OUTPUT_FILE_PRE);
    char *OUTPUT_FILE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE, "%s/input/%s%s%s%s_nwp_B%dL%dT%d_O%d_OL%d_SCVT.nc", model_home_dir, year_string, month_string, day_string, hour, RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
    
    // These are the arrays of the background state.
    double *temperature_gas_background = malloc(NO_OF_SCALARS*sizeof(double));
    double *density_dry_background = malloc(NO_OF_SCALARS*sizeof(double));
    double *wind_background = malloc(NO_OF_VECTORS*sizeof(double));
    double *water_vapour_density_background = malloc(NO_OF_SCALARS*sizeof(double));
    double *liquid_water_density_background = malloc(NO_OF_SCALARS*sizeof(double));
    double *solid_water_density_background = malloc(NO_OF_SCALARS*sizeof(double));
    double *liquid_water_temperature_background = malloc(NO_OF_SCALARS*sizeof(double));
    double *solid_water_temperature_background = malloc(NO_OF_SCALARS*sizeof(double));
    int ncid;
    if ((retval = nc_open(BACKGROUND_STATE_FILE, NC_NOWRITE, &ncid)))
        NCERR(retval);
    int temperature_gas_background_id, density_dry_background_id, wind_background_id, density_vapour_background_id, density_liquid_background_id, density_solid_background_id, temperature_liquid_background_id, temperature_solid_background_id;
    if ((retval = nc_inq_varid(ncid, "density_dry", &density_dry_background_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "temperature_gas",&temperature_gas_background_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "wind", &wind_background_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "density_vapour", &density_vapour_background_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "density_liquid", &density_liquid_background_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "density_solid", &density_solid_background_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "temperature_liquid", &temperature_liquid_background_id)))
        NCERR(retval);
    if ((retval = nc_inq_varid(ncid, "temperature_solid", &temperature_solid_background_id)))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, temperature_gas_background_id, &temperature_gas_background[0])))
        NCERR(retval);    
    if ((retval = nc_get_var_double(ncid, density_dry_background_id, &density_dry_background[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, wind_background_id, &wind_background[0])))
        NCERR(retval);    
    if ((retval = nc_get_var_double(ncid, density_vapour_background_id, &water_vapour_density_background[0])))
        NCERR(retval);    
    if ((retval = nc_get_var_double(ncid, density_liquid_background_id, &liquid_water_density_background[0])))
        NCERR(retval);    
    if ((retval = nc_get_var_double(ncid, density_solid_background_id, &solid_water_density_background[0])))
        NCERR(retval);
    if ((retval = nc_get_var_double(ncid, temperature_liquid_background_id, &liquid_water_temperature_background[0])))
        NCERR(retval);    
    if ((retval = nc_get_var_double(ncid, temperature_solid_background_id, &solid_water_temperature_background[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid)))
        NCERR(retval);
    
    // These are the arrays for the result of the assimilation process.
    double *temperature = malloc(NO_OF_SCALARS*sizeof(double));
    double *density_dry = malloc(NO_OF_SCALARS*sizeof(double));
    double *wind = malloc(NO_OF_VECTORS*sizeof(double));
    double *water_vapour_density = malloc(NO_OF_SCALARS*sizeof(double));
    double *liquid_water_density = malloc(NO_OF_SCALARS*sizeof(double));
    double *solid_water_density = malloc(NO_OF_SCALARS*sizeof(double));
    double *liquid_water_temp = malloc(NO_OF_SCALARS*sizeof(double));
    double *solid_water_temp = malloc(NO_OF_SCALARS*sizeof(double));
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	temperature[i] = temperature_gas_background[i];
    	density_dry[i] = density_dry_background[i];
		water_vapour_density[i] = water_vapour_density_background[i];
		liquid_water_density[i] = liquid_water_density_background[i];
		solid_water_density[i] = solid_water_density_background[i];
		liquid_water_temp[i] = liquid_water_temperature_background[i];
		solid_water_temp[i] = solid_water_temperature_background[i];
    }
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
    	wind[i] = wind_background[i];
    }
    
    // Freeing the grid properties.
    free(z_vector);
    free(latitude_scalar);
    free(longitude_scalar);
    free(direction);
    free(latitude_vector);
    free(longitude_vector);
    
    // Writing the result to a netcdf file.
    int scalar_dimid, vector_dimid, temp_id, density_dry_id, wind_id, density_vapour_id, density_liquid_id, density_solid_id, temperature_liquid_id, temperature_solid_id;
    if ((retval = nc_create(OUTPUT_FILE, NC_CLOBBER, &ncid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "scalar_index", NO_OF_SCALARS, &scalar_dimid)))
        NCERR(retval);
    if ((retval = nc_def_dim(ncid, "vector_index", NO_OF_VECTORS, &vector_dimid)))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "temperature_gas", NC_DOUBLE, 1, &scalar_dimid, &temp_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, temp_id, "units", strlen("K"), "K")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "density_dry", NC_DOUBLE, 1, &scalar_dimid, &density_dry_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, density_dry_id, "units", strlen("kg/m^3"), "kg/m^3")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "wind", NC_DOUBLE, 1, &vector_dimid, &wind_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, wind_background_id, "units", strlen("m/s"), "m/s")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "density_vapour", NC_DOUBLE, 1, &scalar_dimid, &density_vapour_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, density_vapour_id, "units", strlen("kg/m^3"), "kg/m^3")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "density_liquid", NC_DOUBLE, 1, &scalar_dimid, &density_liquid_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, density_liquid_id, "units", strlen("kg/m^3"), "kg/m^3")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "density_solid", NC_DOUBLE, 1, &scalar_dimid, &density_solid_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, density_solid_id, "units", strlen("kg/m^3"), "kg/m^3")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "temperature_liquid", NC_DOUBLE, 1, &scalar_dimid, &temperature_liquid_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, temperature_liquid_id, "units", strlen("T"), "T")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "temperature_solid", NC_DOUBLE, 1, &scalar_dimid, &temperature_solid_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, temperature_solid_id, "units", strlen("T"), "T")))
        NCERR(retval);
    if ((retval = nc_enddef(ncid)))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, temp_id, &temperature[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, density_dry_id, &density_dry[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, wind_id, &wind[0])))
        NCERR(retval);    
    if ((retval = nc_put_var_double(ncid, density_vapour_id, &water_vapour_density[0])))
        NCERR(retval);    
    if ((retval = nc_put_var_double(ncid, density_liquid_id, &liquid_water_density[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, density_solid_id, &solid_water_density[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, temperature_liquid_id, &liquid_water_temp[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, temperature_solid_id, &solid_water_temp[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid)))
    	NCERR(retval);

    free(model_home_dir);
	free(temperature);
	free(density_dry);
	free(wind);
	free(water_vapour_density);
	free(liquid_water_density);
	free(solid_water_density);
	free(liquid_water_temp);
	free(solid_water_temp);
	free(year_string);
	free(month_string);
	free(day_string);
	free(hour);
    free(wind_background);
    free(temperature_gas_background);
    free(density_dry_background);
    free(water_vapour_density_background);
    free(liquid_water_density_background);
    free(solid_water_density_background);
    free(liquid_water_temperature_background);
    free(solid_water_temperature_background);
    free(z_scalar);
    free(OUTPUT_FILE);
    return 0;
}












