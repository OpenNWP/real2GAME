// the number of levels from which we use observations
#define NO_OF_LEVELS_OBS 6
// the number of fields we use on each layer
#define NO_OF_FIELDS_PER_LAYER_OBS 1
// the number of surface variables (order: surface pressure, precipitation rate)
#define NO_OF_SURFACE_FIELDS_OBS 1
// the number of points per layer of the input model
#define NO_OF_POINTS_PER_LAYER_OBS 2949120
// the number of points on each layer
#define NO_OF_CHOSEN_POINTS_PER_LAYER_OBS 200
// the total number of observations we take into account
#define NO_OF_CHOSEN_OBSERVATIONS ((NO_OF_LEVELS_OBS*NO_OF_FIELDS_PER_LAYER_OBS + NO_OF_SURFACE_FIELDS_OBS)*NO_OF_CHOSEN_POINTS_PER_LAYER_OBS)
// the total number of model degree of freedoms
#define NO_OF_MODEL_DOFS (NO_OF_SCALARS + NO_OF_SURFACE_FIELDS_OBS*NO_OF_SCALARS_H)
// the number of model degrees of freedom that are used to interpolate to an observation
#define NO_OF_REL_MODEL_DOFS_PER_OBS 400
// values smaller than one smooth the MSLP
#define SP_INTERPOL_EXP 0.5
// values smaller than one smooth the temperature field
#define T_INTERPOL_EXP 0.75

int oi(double [], double [][NO_OF_REL_MODEL_DOFS_PER_OBS], int [][NO_OF_REL_MODEL_DOFS_PER_OBS], double [], double [], double [], double [], double [], int);
int inv_gauss(double [][NO_OF_CHOSEN_OBSERVATIONS], double [][NO_OF_CHOSEN_OBSERVATIONS]);
int inv_lu(double [][NO_OF_CHOSEN_OBSERVATIONS], double [][NO_OF_CHOSEN_OBSERVATIONS]);
