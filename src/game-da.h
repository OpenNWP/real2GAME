/*
This source file is part of GAME-DA, which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME-DA
*/

// the number of levels from which we use observations
#define NO_OF_LEVELS_OBS 6
// the number of points per layer of the input model
#define NO_OF_POINTS_PER_LAYER_OBS 2949120
// the number of points on each layer
#define NO_OF_CHOSEN_POINTS_PER_LAYER_OBS 1000
// the total number of observations we take into account for the dry assimilation
#define NO_OF_CHOSEN_OBSERVATIONS_DRY ((NO_OF_LEVELS_OBS + 1)*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS)
// the total number of observations we take into account for the moist assimilation
#define NO_OF_CHOSEN_OBSERVATIONS_MOIST (NO_OF_LEVELS_OBS*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS)
#define NO_OF_CHOSEN_OBSERVATIONS (NO_OF_CHOSEN_OBSERVATIONS_DRY + NO_OF_CHOSEN_OBSERVATIONS_MOIST)
// the total number of model degree of freedoms for the dry assimilation
#define NO_OF_MODEL_DOFS_DRY (NO_OF_SCALARS + NO_OF_SCALARS_H)
// the total number of model degree of freedoms for the moisture assimilation
#define NO_OF_MODEL_DOFS_MOIST NO_OF_SCALARS
// the number of model degrees of freedom that are used to interpolate to an observation
#define NO_OF_REL_MODEL_DOFS_PER_OBS 14
// values smaller than one smooth the MSLP field
#define SP_INTERPOL_EXP 0.5
// values smaller than one smooth the temperature field
#define T_INTERPOL_EXP 0.75

int oi(double [], double [][NO_OF_REL_MODEL_DOFS_PER_OBS], int [][NO_OF_REL_MODEL_DOFS_PER_OBS], double [], double [], double [], double [], double [], int, int);
int inv_gauss_dry(double [][NO_OF_CHOSEN_OBSERVATIONS_DRY], double [][NO_OF_CHOSEN_OBSERVATIONS_DRY]);
int inv_gauss_moist(double [][NO_OF_CHOSEN_OBSERVATIONS_MOIST], double [][NO_OF_CHOSEN_OBSERVATIONS_MOIST]);
