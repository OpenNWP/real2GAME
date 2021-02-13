// the number of levels from which we use observations
#define NO_OF_LEVELS_OBS 6
// the number of fields we use on each layer
#define NO_OF_FIELDS_PER_LAYER_OBS 1
// the number of surface variables (order: surface pressure, precipitation rate)
#define NO_OF_SURFACE_FIELDS_OBS 0
// the number of points on each layer
#define NO_OF_POINTS_PER_LAYER_OBS 1020
// the total number of observations
#define NO_OF_OBSERVATIONS ((NO_OF_LEVELS_OBS*NO_OF_FIELDS_PER_LAYER_OBS + NO_OF_SURFACE_FIELDS_OBS)*NO_OF_POINTS_PER_LAYER_OBS)
// the number of model grid points that are used to interpolate to an observation
#define NO_OF_REL_MODEL_DOFS 10
// the total number of model degree of freedoms
#define NO_OF_MODEL_DOFS (NO_OF_SCALARS + NO_OF_SURFACE_FIELDS_OBS*NO_OF_SCALARS_H)

int oi(double [], double [][NO_OF_REL_MODEL_DOFS], int [][NO_OF_REL_MODEL_DOFS], double [], double [], double [], double [], double []);
