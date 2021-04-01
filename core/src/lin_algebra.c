/*
This source file is part of ndvar, which is released under the MIT license.
Github repository: https://github.com/AUN4GFD/ndvar
*/

/*
Linear algebra functions.
*/

#include <stdlib.h>
#include "ndvar.h"

int inv_gauss(double to_be_inverted[][NO_OF_CHOSEN_OBSERVATIONS], double inv[][NO_OF_CHOSEN_OBSERVATIONS])
{
	/*
	This function computes the inverse inv of the matrix to_be_inverted, using the Gauss scheme.
	CAUTION: in the process, to_be_inverted will be modified.
	*/
	
	// firstly, the inverse is initialized with the unity matrix
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS; ++i)
	{
		inv[i][i] = 1;
	}
	
	/*
	Gaussian downwards
	------------------
	we will start to modify to_be_inverted now (misuse of name)
	*/
	double factor;
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS - 1; ++i)
	{
		// dividing the line by to_be_inverted[i][i]
		factor = 1/to_be_inverted[i][i];
		for (int j = i; j < NO_OF_CHOSEN_OBSERVATIONS; ++j)
		{
			to_be_inverted[i][j] = factor*to_be_inverted[i][j];
		}
		for (int j = 0; j <= i; ++j)
		{
			inv[i][j] = factor*inv[i][j];
		}
		// loop over all the lines that are below the current line
		for (int j = i + 1; j < NO_OF_CHOSEN_OBSERVATIONS; ++j)
		{
			factor = -to_be_inverted[j][i];
			for (int k = i; k < NO_OF_CHOSEN_OBSERVATIONS; ++k)
			{
				to_be_inverted[j][k] = to_be_inverted[j][k] + factor*to_be_inverted[i][k];
			}
			for (int k = 0; k <= i; ++k)
			{
				inv[j][k] = inv[j][k] + factor*inv[i][k];
			}
		}
	}
	
	for (int j = 0; j < NO_OF_CHOSEN_OBSERVATIONS; ++j)
	{
		inv[NO_OF_CHOSEN_OBSERVATIONS - 1][j] = inv[NO_OF_CHOSEN_OBSERVATIONS - 1][j]/to_be_inverted[NO_OF_CHOSEN_OBSERVATIONS - 1][NO_OF_CHOSEN_OBSERVATIONS - 1];
	}
	to_be_inverted[NO_OF_CHOSEN_OBSERVATIONS - 1][NO_OF_CHOSEN_OBSERVATIONS - 1] = 1;
	
	/*
	Gaussian upwards
	----------------
	*/
	for (int i = NO_OF_CHOSEN_OBSERVATIONS - 1; i >= 1; --i)
	{
		for (int j = i - 1; j >= 0; --j)
		{
			for (int k = 0; k < NO_OF_CHOSEN_OBSERVATIONS; ++k)
			{
				inv[j][k] = inv[j][k] - to_be_inverted[j][i]*inv[i][k];
			}
		}
	}
	return 0;
}

int inv_lu(double to_be_inverted[][NO_OF_CHOSEN_OBSERVATIONS], double inv[][NO_OF_CHOSEN_OBSERVATIONS])
{
	/*
	This function computes the inverse inv of the matrix to_be_inverted, using the LU decomposition.
	CAUTION: in the process, to_be_inverted will be modified.
	*/
	double (*l_matrix)[NO_OF_CHOSEN_OBSERVATIONS] = calloc(1, sizeof(double[NO_OF_CHOSEN_OBSERVATIONS][NO_OF_CHOSEN_OBSERVATIONS]));
	/*
	downward sweep
	--------------
	to_be_inverted will become the r_matrix now (misuse of name)
	*/
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS - 1; ++i)
	{
		l_matrix[i][i] = 1;
		for (int j = i + 1; j < NO_OF_CHOSEN_OBSERVATIONS; ++j)
		{
			l_matrix[j][i] = to_be_inverted[j][i]/to_be_inverted[i][i];
			for (int k = i; k < NO_OF_CHOSEN_OBSERVATIONS; ++k)
			{
				to_be_inverted[j][k] = to_be_inverted[j][k] - l_matrix[j][i]*to_be_inverted[i][k];
			}
		}
	}
	
	/*
	The LU decomposition is already done at this point.
	We now use the LU decomposition for the inversion.
	We know LU = A. We want to solve AA^-1 = 1, a.k.a. LUA^-1 = 1.
	Therefore, we firstly solve LB = 1 with a downward sweep.
	*/
	double (*b_matrix)[NO_OF_CHOSEN_OBSERVATIONS] = calloc(1, sizeof(double[NO_OF_CHOSEN_OBSERVATIONS][NO_OF_CHOSEN_OBSERVATIONS]));
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS - 1; ++i)
	{
		b_matrix[i][i] = 1/l_matrix[i][i];
		for (int j = 0; j < i; ++j)
		{
			for (int k = 0; k < i; ++k)
			{
				b_matrix[i][j] -= l_matrix[i][k]*b_matrix[k][j]/l_matrix[i][i];
			}
		}
	}
	// l_matrix is not needed anymore
	free(l_matrix);
	
	/*
	Now we have to solve UA^-1 = B with an upward sweep.
	U is to_be_inverted (see above).
	*/
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS; ++i)
	{
		inv[NO_OF_CHOSEN_OBSERVATIONS - 1][i] = b_matrix[NO_OF_CHOSEN_OBSERVATIONS - 1][i]/to_be_inverted[NO_OF_CHOSEN_OBSERVATIONS - 1][NO_OF_CHOSEN_OBSERVATIONS - 1];
		for (int j = NO_OF_CHOSEN_OBSERVATIONS - 2; j >= 0; --j)
		{
			for (int k = j + 1; k < NO_OF_CHOSEN_OBSERVATIONS ; ++j)
			{
				inv[j][i] = b_matrix[j][i] - to_be_inverted[j][k]*inv[k][i]/to_be_inverted[j][j];
			}
		}
	}
	
	// that's it, b_matrix is not needed anymore
	free(b_matrix);
	return 0;
}










