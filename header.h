/*
This source file is part of real2GAME, which is released under the MIT license.
Github repository: https://github.com/OpenNWP/real2GAME
*/

/*
These things can be considered the config of the data assimilation:
------------------------------------------------------------------
*/
// the number of levels from which we use observations
#define NO_OF_LEVELS_OBS 12
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















