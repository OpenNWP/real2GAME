/*
This source file is part of ndvar, which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/ndvar
*/

/*
Linear algebra functions.
*/

#include "ndvar.h"

int inv_gauss(double to_be_inverted[][NO_OF_CHOSEN_OBSERVATIONS], double inv[][NO_OF_CHOSEN_OBSERVATIONS])
{
	/*
	This function computes the inverse inv of the matrix to_be_inverted, using the Gauss scheme.
	*/
	// firstly, the inverse is initialized with the unity matrix
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS; ++i)
	{
		for (int j = 0; j < NO_OF_CHOSEN_OBSERVATIONS; ++j)
		{
			inv[i][j] = 0;
			if (i == j)
			{
				inv[i][j] = 1;
			}
		}
	}
	// we will start to modify to_be_inverted now (misuse of name)
	// Gaussian downwards
	double factor = 0;
	// starting with the first line, down to the second but last line
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS - 1; ++i)
	{
		// dividing the line by to_be_inverted[i][i]
		factor = 1/to_be_inverted[i][i];
		for (int j = 0; j < NO_OF_CHOSEN_OBSERVATIONS; ++j)
		{
			to_be_inverted[i][j] = factor*to_be_inverted[i][j];
			inv[i][j] = factor*inv[i][j];
		}
		for (int j = i + 1; j < NO_OF_CHOSEN_OBSERVATIONS; ++j)
		{
			factor = -to_be_inverted[j][i];
			for (int k = 0; k < NO_OF_CHOSEN_OBSERVATIONS; ++k)
			{
				to_be_inverted[j][k] = to_be_inverted[j][k] + factor*to_be_inverted[i][k];
				inv[j][k] = inv[j][k] + factor*inv[i][k];
			}
		}
	}
	factor = 1/to_be_inverted[NO_OF_CHOSEN_OBSERVATIONS - 1][NO_OF_CHOSEN_OBSERVATIONS - 1];
	for (int j = 0; j < NO_OF_CHOSEN_OBSERVATIONS; ++j)
	{
		to_be_inverted[NO_OF_CHOSEN_OBSERVATIONS - 1][j] = factor*to_be_inverted[NO_OF_CHOSEN_OBSERVATIONS - 1][j];
		inv[NO_OF_CHOSEN_OBSERVATIONS - 1][j] = factor*inv[NO_OF_CHOSEN_OBSERVATIONS - 1][j];
	}
	// Gaussian upwards
	// starting with the last line, then going up to the last but first
	for (int i = NO_OF_CHOSEN_OBSERVATIONS - 1; i >= 1; --i)
	{
		for (int j = i - 1; j >= 0; --j)
		{
			factor = -to_be_inverted[j][i];
			for (int k = 0; k < NO_OF_CHOSEN_OBSERVATIONS; ++k)
			{
				inv[j][k] = inv[j][k] + factor*inv[i][k];
			}
		}
	}
	return 0;
}

int inv_lu(double to_be_inverted[][NO_OF_CHOSEN_OBSERVATIONS], double inv[][NO_OF_CHOSEN_OBSERVATIONS])
{
	/*
	This function computes the inverse inv of the matrix to_be_inverted, using the LU decomposition.
	*/
	// firstly, the inverse is initialized with the unity matrix
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS; ++i)
	{
		for (int j = 0; j < NO_OF_CHOSEN_OBSERVATIONS; ++j)
		{
			inv[i][j] = 0;
			if (i == j)
			{
				inv[i][j] = 1;
			}
		}
	}
	// we will start to modify to_be_inverted now (misuse of name)
	// Gaussian downwards
	double factor = 0;
	// starting with the first line, down to the second but last line
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS - 1; ++i)
	{
		// dividing the line by to_be_inverted[i][i]
		factor = 1/to_be_inverted[i][i];
		for (int j = 0; j < NO_OF_CHOSEN_OBSERVATIONS; ++j)
		{
			to_be_inverted[i][j] = factor*to_be_inverted[i][j];
			inv[i][j] = factor*inv[i][j];
		}
		for (int j = i + 1; j < NO_OF_CHOSEN_OBSERVATIONS; ++j)
		{
			factor = -to_be_inverted[j][i];
			for (int k = 0; k < NO_OF_CHOSEN_OBSERVATIONS; ++k)
			{
				to_be_inverted[j][k] = to_be_inverted[j][k] + factor*to_be_inverted[i][k];
				inv[j][k] = inv[j][k] + factor*inv[i][k];
			}
		}
	}
	factor = 1/to_be_inverted[NO_OF_CHOSEN_OBSERVATIONS - 1][NO_OF_CHOSEN_OBSERVATIONS - 1];
	for (int j = 0; j < NO_OF_CHOSEN_OBSERVATIONS; ++j)
	{
		to_be_inverted[NO_OF_CHOSEN_OBSERVATIONS - 1][j] = factor*to_be_inverted[NO_OF_CHOSEN_OBSERVATIONS - 1][j];
		inv[NO_OF_CHOSEN_OBSERVATIONS - 1][j] = factor*inv[NO_OF_CHOSEN_OBSERVATIONS - 1][j];
	}
	// Gaussian upwards
	// starting with the last line, then going up to the last but first
	for (int i = NO_OF_CHOSEN_OBSERVATIONS - 1; i >= 1; --i)
	{
		for (int j = i - 1; j >= 0; --j)
		{
			factor = -to_be_inverted[j][i];
			for (int k = 0; k < NO_OF_CHOSEN_OBSERVATIONS; ++k)
			{
				inv[j][k] = inv[j][k] + factor*inv[i][k];
			}
		}
	}
	return 0;
}














