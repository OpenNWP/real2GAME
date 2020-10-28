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

// constants needed for the standard atmosphere
const double G = 9.80616;
const double TROPO_HEIGHT = 12e3;
const double T_SFC = 273.15 + 15;
double T_0 = 288;
const double TEMP_GRADIENT = -0.65/100;
double GAMMA = 0.005;

int main(int argc, char *argv[])
{	
    size_t len = strlen(argv[1]);
    char *year = malloc((len + 1)*sizeof(char));
    strcpy(year, argv[1]);
    len = strlen(argv[2]);
    char *month = malloc((len + 1)*sizeof(char));
    strcpy(month, argv[2]);
    len = strlen(argv[3]);
    char *day = malloc((len + 1)*sizeof(char));
    strcpy(day, argv[3]);
    len = strlen(argv[4]);
    char *hour = malloc((len + 1)*sizeof(char));
    strcpy(hour, argv[4]);
	int ORO_ID = 3;
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
    sprintf(GEO_PROP_FILE_PRE, "/home/max/compiled/game_dev/grids/B%dL%dT%d_O%d_OL%d_SCVT.nc", RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
    GEO_PROP_FILE_LENGTH = strlen(GEO_PROP_FILE_PRE);
    free(GEO_PROP_FILE_PRE);
    char *GEO_PROP_FILE = malloc((GEO_PROP_FILE_LENGTH + 1)*sizeof(char));
    sprintf(GEO_PROP_FILE, "/home/max/compiled/game_dev/grids/B%dL%dT%d_O%d_OL%d_SCVT.nc", RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
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
    sprintf(OUTPUT_FILE_PRE, "/home/max/compiled/game_dev/input/%s%s%s%s_nwp_B%dL%dT%d_O%d_OL%d_SCVT.nc", year, month, day, hour, RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
    OUTPUT_FILE_LENGTH = strlen(OUTPUT_FILE_PRE);
    free(OUTPUT_FILE_PRE);
    char *OUTPUT_FILE = malloc((OUTPUT_FILE_LENGTH + 1)*sizeof(char));
    sprintf(OUTPUT_FILE, "/home/max/compiled/game_dev/input/%s%s%s%s_nwp_B%dL%dT%d_O%d_OL%d_SCVT.nc", year, month, day, hour, RES_ID, NO_OF_LAYERS, (int) TOA, ORO_ID, NO_OF_ORO_LAYERS);
    double *pressure_background = malloc(NO_OF_SCALARS*sizeof(double));
    double *temperature_background = malloc(NO_OF_SCALARS*sizeof(double));
    double *density_dry_background = malloc(NO_OF_SCALARS*sizeof(double));
    double *rel_humidity_background = malloc(NO_OF_SCALARS*sizeof(double));
    double *wind_background = malloc(NO_OF_VECTORS*sizeof(double));
    double *water_vapour_density_background = malloc(NO_OF_SCALARS*sizeof(double));
    double *liquid_water_density_background = malloc(NO_OF_SCALARS*sizeof(double));
    double *solid_water_density_background = malloc(NO_OF_SCALARS*sizeof(double));
    double *liquid_water_temp_background = malloc(NO_OF_SCALARS*sizeof(double));
    double *solid_water_temp_background = malloc(NO_OF_SCALARS*sizeof(double));
    const double TROPO_TEMP = T_SFC + TROPO_HEIGHT*TEMP_GRADIENT;
    double z_height;
    // 3D scalar fields determined here, apart from density
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
        z_height = z_scalar[i];
        rel_humidity_background[i] = 0;
        // standard atmosphere
        if (z_height < TROPO_HEIGHT)
        {
            temperature_background[i] = T_SFC + z_height*TEMP_GRADIENT;
            pressure_background[i] = 101325*pow(1 + TEMP_GRADIENT*z_height/T_SFC, -G/(R_D*TEMP_GRADIENT));
        }
        else
        {
            temperature_background[i] = TROPO_TEMP;
            pressure_background[i] = 101325*pow(1 + TEMP_GRADIENT*TROPO_HEIGHT/T_SFC, -G/(R_D*TEMP_GRADIENT))*exp(-G*(z_height - TROPO_HEIGHT)/(R_D*TROPO_TEMP));
        }
        liquid_water_density_background[i] = 0;
        solid_water_density_background[i] = 0;
        liquid_water_temp_background[i] = temperature_background[i];
        solid_water_temp_background[i] = temperature_background[i];
    }
    // density is determined out of the hydrostatic equation
    int layer_index;
    double entropy_value, temperature_background_mean, delta_temperature_background, delta_gravity_potential, pot_temp_value, lower_entropy_value, pressure_background_value;
    for (int i = NO_OF_SCALARS - 1; i >= 0; --i)
    {
    	layer_index = i/NO_OF_SCALARS_H;
    	// at the lowest layer the density is set using the equation of state, can be considered a boundary condition
    	if (layer_index == NO_OF_LAYERS - 1)
    	{
        	density_dry_background[i] = pressure_background[i]/(R_D*temperature_background[i]);
        }
        else
        {
        	lower_entropy_value = C_D_P*log(temperature_background[i + NO_OF_SCALARS_H]*pow(P_0/(density_dry_background[i + NO_OF_SCALARS_H]*R_D*temperature_background[i + NO_OF_SCALARS_H]), R_D/C_D_P));
        	temperature_background_mean = 0.5*(temperature_background[i] + temperature_background[i + NO_OF_SCALARS_H]);
        	delta_temperature_background = temperature_background[i] - temperature_background[i + NO_OF_SCALARS_H];
        	delta_gravity_potential = gravity_potential[i] - gravity_potential[i + NO_OF_SCALARS_H];
        	entropy_value = lower_entropy_value + (delta_gravity_potential + C_D_P*delta_temperature_background)/temperature_background_mean;
        	pot_temp_value = exp(entropy_value/C_D_P);
        	pressure_background_value = P_0*pow(temperature_background[i]/pot_temp_value, C_D_P/R_D);
        	density_dry_background[i] = pressure_background_value/(R_D*temperature_background[i]);
        }
        water_vapour_density_background[i] = water_vapour_density_from_rel_humidity(rel_humidity_background[i], temperature_background[i], density_dry_background[i]);
        if (water_vapour_density_background[i] < 0)
        	printf("water_vapour_density_background negative.\n.");
    }
    free(gravity_potential);
    for (int i = 0; i < NO_OF_LAYERS; ++i)
    {
        // horizontal wind_background fields are determind here
        for (int j = 0; j < NO_OF_VECTORS_H; ++j)
        {
            z_height = z_vector[NO_OF_SCALARS_H + j + i*NO_OF_VECTORS_PER_LAYER];
            // standard atmosphere: no wind_background
            wind_background[NO_OF_SCALARS_H + i*NO_OF_VECTORS_PER_LAYER + j] = 0;
        }
    }
    for (int i = 0; i < NO_OF_LEVELS; ++i)
    {
        for (int j = 0; j < NO_OF_SCALARS_H; ++j)
        {
            z_height = z_vector[j + i*NO_OF_VECTORS_PER_LAYER];
            wind_background[i*NO_OF_VECTORS_PER_LAYER + j] = 0;
        }
    }
    double *temperature = malloc(NO_OF_SCALARS*sizeof(double));
    double *density_dry = malloc(NO_OF_SCALARS*sizeof(double));
    double *rel_humidity = malloc(NO_OF_SCALARS*sizeof(double));
    double *wind = malloc(NO_OF_VECTORS*sizeof(double));
    double *water_vapour_density = malloc(NO_OF_SCALARS*sizeof(double));
    double *liquid_water_density = malloc(NO_OF_SCALARS*sizeof(double));
    double *solid_water_density = malloc(NO_OF_SCALARS*sizeof(double));
    double *liquid_water_temp = malloc(NO_OF_SCALARS*sizeof(double));
    double *solid_water_temp = malloc(NO_OF_SCALARS*sizeof(double));
    for (int i = 0; i < NO_OF_SCALARS; ++i)
    {
    	temperature[i] = temperature_background[i];
    	density_dry[i] = density_dry_background[i];
    	rel_humidity[i] = rel_humidity_background[i];
		water_vapour_density[i] = water_vapour_density_background[i];
		liquid_water_density[i] = liquid_water_density_background[i];
		solid_water_density[i] = solid_water_density_background[i];
		liquid_water_temp[i] = liquid_water_temp_background[i];
		solid_water_temp[i] = solid_water_temp_background[i];
    }
    for (int i = 0; i < NO_OF_VECTORS; ++i)
    {
    	wind[i] = wind_background[i];
    }
    free(z_vector);
    free(latitude_scalar);
    free(longitude_scalar);
    free(direction);
    free(latitude_vector);
    free(longitude_vector);
    int scalar_dimid, vector_dimid, temp_id, density_dry_id, wind_background_id, density_vapour_id, density_liquid_id, density_solid_id, temperature_background_liquid_id, temperature_background_solid_id, ncid;
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
    if ((retval = nc_def_var(ncid, "wind", NC_DOUBLE, 1, &vector_dimid, &wind_background_id)))
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
    if ((retval = nc_def_var(ncid, "temperature_liquid", NC_DOUBLE, 1, &scalar_dimid, &temperature_background_liquid_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, temperature_background_liquid_id, "units", strlen("T"), "T")))
        NCERR(retval);
    if ((retval = nc_def_var(ncid, "temperature_solid", NC_DOUBLE, 1, &scalar_dimid, &temperature_background_solid_id)))
        NCERR(retval);
    if ((retval = nc_put_att_text(ncid, temperature_background_solid_id, "units", strlen("T"), "T")))
        NCERR(retval);
    if ((retval = nc_enddef(ncid)))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, temp_id, &temperature[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, density_dry_id, &density_dry[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, wind_background_id, &wind[0])))
        NCERR(retval);    
    if ((retval = nc_put_var_double(ncid, density_vapour_id, &water_vapour_density[0])))
        NCERR(retval);    
    if ((retval = nc_put_var_double(ncid, density_liquid_id, &liquid_water_density[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, density_solid_id, &solid_water_density[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, temperature_background_liquid_id, &liquid_water_temp[0])))
        NCERR(retval);
    if ((retval = nc_put_var_double(ncid, temperature_background_solid_id, &solid_water_temp[0])))
        NCERR(retval);
    if ((retval = nc_close(ncid)))
    	NCERR(retval);
	free(temperature);
	free(density_dry);
	free(rel_humidity);
	free(wind);
	free(water_vapour_density);
	free(liquid_water_density);
	free(solid_water_density);
	free(liquid_water_temp);
	free(solid_water_temp);
	free(year);
	free(month);
	free(day);
	free(hour);
    free(wind_background);
    free(pressure_background);
    free(temperature_background);
    free(density_dry_background);
    free(water_vapour_density_background);
    free(liquid_water_density_background);
    free(solid_water_density_background);
    free(liquid_water_temp_background);
    free(solid_water_temp_background);
    free(z_scalar);
    free(rel_humidity_background);
    free(OUTPUT_FILE);
    return 0;
}












