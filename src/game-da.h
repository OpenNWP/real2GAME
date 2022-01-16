/*
This source file is part of GAME-DA, which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME-DA
*/

/*
These things can be considered the config of the data assimilation:
------------------------------------------------------------------
*/
// the number of levels from which we use observations
#define NO_OF_LEVELS_OBS 6
// the number of points on each layer
#define NO_OF_CHOSEN_POINTS_PER_LAYER_OBS 1200
// the number of wind points on each layer
#define NO_OF_CHOSEN_WIND_POINTS_PER_LAYER_OBS 600
// the number of model degrees of freedom that are used to interpolate to an observation (must be even)
#define NO_OF_REL_MODEL_DOFS_PER_OBS 14
// values smaller than two smooth the fields
#define INTERPOL_EXP 2.001

/*
These things should rarely need to be modified by the user:
----------------------------------------------------------
*/
// the number of points per layer of the input model
#define NO_OF_POINTS_PER_LAYER_OBS 2949120
// the number of points of the SST grid
#define NO_OF_SST_POINTS 259200

/*
These things should never need to be modified by the user:
---------------------------------------------------------
*/
// the total number of observations we take into account for the dry thermodynamic assimilation
#define NO_OF_CHOSEN_OBSERVATIONS_DRY ((NO_OF_LEVELS_OBS + 1)*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS)
// the total number of observations we take into account for the moist assimilation
#define NO_OF_CHOSEN_OBSERVATIONS_MOIST (NO_OF_LEVELS_OBS*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS)
// the total number of observations we take into account for the wind assimilation (factor 2 for two directions in the horizontal)
#define NO_OF_CHOSEN_OBSERVATIONS_WIND (NO_OF_LEVELS_OBS*2*NO_OF_CHOSEN_WIND_POINTS_PER_LAYER_OBS)
// the total number of observations including SST
#define NO_OF_CHOSEN_OBSERVATIONS (NO_OF_CHOSEN_OBSERVATIONS_DRY + NO_OF_CHOSEN_OBSERVATIONS_MOIST + NO_OF_CHOSEN_OBSERVATIONS_WIND + NO_OF_SST_POINTS)
// the total number of model degree of freedoms for the dry assimilation
#define NO_OF_MODEL_DOFS_DRY (NO_OF_SCALARS + NO_OF_SCALARS_H)
// the total number of model degree of freedoms for the moisture assimilation
#define NO_OF_MODEL_DOFS_MOIST NO_OF_SCALARS

/*
function declarations
---------------------
*/

int oi(double [], double [][NO_OF_REL_MODEL_DOFS_PER_OBS], int [][NO_OF_REL_MODEL_DOFS_PER_OBS], double [][7], double [], double [], double [], double [], int, int);
int inv_gauss_dry(double [][NO_OF_CHOSEN_OBSERVATIONS_DRY], double [][NO_OF_CHOSEN_OBSERVATIONS_DRY]);
int inv_gauss_moist(double [][NO_OF_CHOSEN_OBSERVATIONS_MOIST], double [][NO_OF_CHOSEN_OBSERVATIONS_MOIST]);
int inv_gauss_wind(double [][NO_OF_CHOSEN_OBSERVATIONS_WIND], double [][NO_OF_CHOSEN_OBSERVATIONS_WIND]);















