/*
This source file is part of ndvar, which is released under the MIT license.
Github repository: https://github.com/OpenNWP/ndvar
*/

/*
Optimum interpolation.
*/

#include "ndvar.h"
#include "enum.h"
#include <stdlib.h>
#include "geos95.h"

int oi(double obs_error_cov[], double obs_op_jacobian_reduced_matrix[][NO_OF_REL_MODEL_DOFS_PER_OBS], int relevant_model_dofs_matrix[][NO_OF_REL_MODEL_DOFS_PER_OBS], double bg_error_cov[], double interpolated_model[], double background[], double observations_vector[], double model_vector[], int no_of_obs, int no_of_model_dofs)
{
	// short notation: b: background error covariance, h: derivative of the observations operator; r: observations error covariance
	double (*h_b_ht_plus_r)[no_of_obs] = malloc(sizeof(double[no_of_obs][no_of_obs]));
	int index_found;
	int check_vector[NO_OF_REL_MODEL_DOFS_PER_OBS];
    //#pragma omp parallel for private(index_found, check_vector)
	for (int i = 0; i < no_of_obs; ++i)
	{
		for (int l = 0; l < NO_OF_REL_MODEL_DOFS_PER_OBS; ++l)
		{
			check_vector[l] = relevant_model_dofs_matrix[i][l];
		}
		for (int j = 0; j < no_of_obs; ++j)
		{
			h_b_ht_plus_r[i][j] = 0;
			for (int k = 0; k < NO_OF_REL_MODEL_DOFS_PER_OBS; ++k)
			{
				if (in_bool_calculator(relevant_model_dofs_matrix[j][k], check_vector, NO_OF_REL_MODEL_DOFS_PER_OBS) == 1)
				{
					index_found = 0;
					for (int l = 0; l < NO_OF_REL_MODEL_DOFS_PER_OBS; ++l)
					{
						if (relevant_model_dofs_matrix[j][l] == relevant_model_dofs_matrix[i][k])
						{
							index_found = l;
						}
					}
					h_b_ht_plus_r[i][j]
					+= bg_error_cov[relevant_model_dofs_matrix[i][k]]
					*obs_op_jacobian_reduced_matrix[i][k]
					*obs_op_jacobian_reduced_matrix[j][index_found];
				}
			}
		}
		h_b_ht_plus_r[i][i] = h_b_ht_plus_r[i][i] + obs_error_cov[i];
	}
	
	// h_b_ht_plus_r needs to be inversed in order to calculate the gain matrix
	// this is actually the main task of OI
	double (*h_b_ht_plus_r_inv)[no_of_obs] = calloc(1, sizeof(double[no_of_obs][no_of_obs]));
	if (no_of_obs == NO_OF_CHOSEN_OBSERVATIONS_DRY)
	{
		inv_gauss_dry(h_b_ht_plus_r, h_b_ht_plus_r_inv);
	}
	if (no_of_obs == NO_OF_CHOSEN_OBSERVATIONS_MOIST)
	{
		inv_gauss_moist(h_b_ht_plus_r, h_b_ht_plus_r_inv);
	}
	
	// now, the main job is already done
	free(h_b_ht_plus_r);
	
	// this vector will contain the product of the model forecast error and the gain matrix
	double *prod_with_gain_matrix = calloc(no_of_model_dofs, sizeof(double));
	// multiplying (obs - (interpolated model)) by the gain matrix
    //#pragma omp parallel for
	for (int i = 0; i < NO_OF_REL_MODEL_DOFS_PER_OBS; ++i)
	{
		for (int j = 0; j < no_of_obs; ++j)
		{
			for (int k = 0; k < no_of_obs; ++k)
			{
				prod_with_gain_matrix[relevant_model_dofs_matrix[k][i]]
				+= (h_b_ht_plus_r_inv[k][j]
				*bg_error_cov[relevant_model_dofs_matrix[k][i]]
				*obs_op_jacobian_reduced_matrix[k][i])
				*(observations_vector[j] - interpolated_model[j]);
			}
		}
	}
	
	free(h_b_ht_plus_r_inv);
	
    #pragma omp parallel for
	for (int i = 0; i < no_of_model_dofs; ++i)
	{
		model_vector[i] = background[i] + prod_with_gain_matrix[i];
	}
	free(prod_with_gain_matrix);
	
	return 0;
}















