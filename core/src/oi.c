/*
This source file is part of ndvar, which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/ndvar
*/
/*
Optimum interpolation.
*/

#include "ndvar.h"
#include "enum.h"
#include <stdlib.h>

int oi(double obs_error_cov[], double obs_op_reduced_matrix[][NO_OF_REL_MODEL_DOFS], int relevant_model_dofs_matrix[][NO_OF_REL_MODEL_DOFS], double bg_error_cov[], double interpolated_model[], double background[], double observations_vector[], double model_vector[])
{
	// short notation: b: background error covariance, h: observations operator; r: observations error covariance
	double (*h_b_ht_plus_r)[NO_OF_OBSERVATIONS] = malloc(sizeof(double[NO_OF_OBSERVATIONS][NO_OF_OBSERVATIONS]));
	for (int i = 0; i < NO_OF_OBSERVATIONS; ++i)
	{
		for (int j = 0; j < NO_OF_OBSERVATIONS; ++j)
		{
			h_b_ht_plus_r[i][j] = 0;
			for (int k = 0; k < NO_OF_REL_MODEL_DOFS; ++k)
			{
				h_b_ht_plus_r[i][j]
				= h_b_ht_plus_r[i][j]
				+ obs_op_reduced_matrix[i][k]
				*bg_error_cov[relevant_model_dofs_matrix[i][k]]
				*obs_op_reduced_matrix[j][k];
			}
			if (i == j)
			{
				h_b_ht_plus_r[i][j] = h_b_ht_plus_r[i][j] + obs_error_cov[i];
			}
		}
	}
	
	// h_b_ht_plus_r needs to be inversed in order to calculate the gain matrix
	// this is actually the main task of OI
	double (*h_b_ht_plus_r_inv)[NO_OF_OBSERVATIONS] = malloc(sizeof(double[NO_OF_OBSERVATIONS][NO_OF_OBSERVATIONS]));
	// firstly, the inverse is initialized with the unity matrix
	for (int i = 0; i < NO_OF_OBSERVATIONS; ++i)
	{
		for (int j = 0; j < NO_OF_OBSERVATIONS; ++j)
		{
			h_b_ht_plus_r_inv[i][j] = 0;
			if (i == j)
			{
				h_b_ht_plus_r_inv[i][j] = 1;
			}
		}
	}
	// we will start to modify h_b_ht_plus_r now (misuse of name)
	// Gaussian downwards
	double factor = 0;
	// starting with the first line, down to the second but last line
	for (int i = 0; i < NO_OF_OBSERVATIONS - 1; ++i)
	{
		// dividing the line by h_b_ht_plus_r[i][i]
		factor = 1/h_b_ht_plus_r[i][i];
		for (int j = 0; j < NO_OF_OBSERVATIONS; ++j)
		{
			h_b_ht_plus_r[i][j] = factor*h_b_ht_plus_r[i][j];
			h_b_ht_plus_r_inv[i][j] = factor*h_b_ht_plus_r_inv[i][j];
		}
		for (int j = i + 1; j < NO_OF_OBSERVATIONS; ++j)
		{
			factor = -h_b_ht_plus_r[j][i];
			for (int k = 0; k < NO_OF_OBSERVATIONS; ++k)
			{
				h_b_ht_plus_r[j][k] = h_b_ht_plus_r[j][k] + factor*h_b_ht_plus_r[i][k];
				h_b_ht_plus_r_inv[j][k] = h_b_ht_plus_r_inv[j][k] + factor*h_b_ht_plus_r_inv[i][k];
			}
		}
	}
	factor = 1/h_b_ht_plus_r[NO_OF_OBSERVATIONS - 1][NO_OF_OBSERVATIONS - 1];
	for (int j = 0; j < NO_OF_OBSERVATIONS; ++j)
	{
		h_b_ht_plus_r[NO_OF_OBSERVATIONS - 1][j] = factor*h_b_ht_plus_r[NO_OF_OBSERVATIONS - 1][j];
		h_b_ht_plus_r_inv[NO_OF_OBSERVATIONS - 1][j] = factor*h_b_ht_plus_r_inv[NO_OF_OBSERVATIONS - 1][j];
	}
	// Gaussian upwards
	// starting with the last line, then going up to the last but first
	for (int i = NO_OF_OBSERVATIONS - 1; i >= 1; --i)
	{
		for (int j = i - 1; j >= 0; --j)
		{
			factor = -h_b_ht_plus_r[j][i];
			for (int k = 0; k < NO_OF_OBSERVATIONS; ++k)
			{
				h_b_ht_plus_r[j][k] = h_b_ht_plus_r[j][k] + factor*h_b_ht_plus_r[i][k];
				h_b_ht_plus_r_inv[j][k] = h_b_ht_plus_r_inv[j][k] + factor*h_b_ht_plus_r_inv[i][k];
			}
		}
	}
	
	// now, the main job is already done
	
	free(h_b_ht_plus_r);
	
	// this vector will contain the product of the model forecast error and the gain matrix
	double *prod_with_gain_matrix = calloc(NO_OF_MODEL_DOFS, sizeof(double));
	// multiplying (obs - (interpolated model)) by the gain matrix
	for (int i = 0; i < NO_OF_REL_MODEL_DOFS; ++i)
	{	
		for (int j = 0; j < NO_OF_OBSERVATIONS; ++j)
		{
			for (int k = 0; k < NO_OF_OBSERVATIONS; ++k)
			{
				prod_with_gain_matrix[relevant_model_dofs_matrix[k][i]]
				+= (bg_error_cov[relevant_model_dofs_matrix[k][i]]
				*obs_op_reduced_matrix[k][i]
				*h_b_ht_plus_r_inv[k][j])
				*(observations_vector[j] - interpolated_model[j]);
			}
		}
	}	
	
	free(h_b_ht_plus_r_inv);
	
	for (int i = 0; i < NO_OF_MODEL_DOFS; ++i)
	{
		model_vector[i] = background[i] + prod_with_gain_matrix[i];
	}
	free(prod_with_gain_matrix);
	
	return 0;
}















