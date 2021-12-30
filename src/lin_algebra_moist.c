/*
This source file is part of GAME-DA, which is released under the MIT license.
Github repository: https://github.com/OpenNWP/GAME-DA
*/

/*
linear algebra functions for the moist assimilation process
*/

#include <stdlib.h>
#include <stdio.h>
#include "game-da.h"

int permute_lines_moist(double [][NO_OF_CHOSEN_OBSERVATIONS_MOIST], int, int);

int inv_gauss_moist(double to_be_inverted[][NO_OF_CHOSEN_OBSERVATIONS_MOIST], double inv[][NO_OF_CHOSEN_OBSERVATIONS_MOIST])
{
	/*
	This function computes the inverse inv of the matrix to_be_inverted, using the Gauss scheme.
	CAUTION: in the process, to_be_inverted will be modified.
	*/
	
	// firstly, the inverse is initialized with the unity matrix
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS_MOIST; ++i)
	{
		inv[i][i] = 1;
	}
	
	/*
	Gaussian downwards
	------------------
	we will start to modify to_be_inverted now (misuse of name)
	*/
	int permute_index_found, permute_index_counter;
	double factor;
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS_MOIST - 1; ++i)
	{
		/*
		checking if a permutation is necessary
		*/
		// Firstly, the permutation index has to be found.
		permute_index_found = 0;
		permute_index_counter = i;
		while (permute_index_found == 0)
		{
			if (to_be_inverted[permute_index_counter][i] != 0)
			{
				permute_index_found = 1;
			}
			else
			{
				permute_index_counter += 1;
			}
		}
		// actually performing the permutation
		if (permute_index_counter > i)
		{
			permute_lines_moist(to_be_inverted, i, permute_index_counter);
			permute_lines_moist(inv, i, permute_index_counter);
		}
		
		// permutation is done, now comes the actual calculation
		// dividing the line by to_be_inverted[i][i]
		factor = 1/to_be_inverted[i][i];
		#pragma omp parallel for
		for (int j = i; j < NO_OF_CHOSEN_OBSERVATIONS_MOIST; ++j)
		{
			to_be_inverted[i][j] = factor*to_be_inverted[i][j];
		}
		#pragma omp parallel for
		for (int j = 0; j < NO_OF_CHOSEN_OBSERVATIONS_MOIST; ++j)
		{
			inv[i][j] = factor*inv[i][j];
		}
		// loop over all the lines that are below the current line
		#pragma omp parallel for private(factor)
		for (int j = i + 1; j < NO_OF_CHOSEN_OBSERVATIONS_MOIST; ++j)
		{
			factor = -to_be_inverted[j][i];
			for (int k = i; k < NO_OF_CHOSEN_OBSERVATIONS_MOIST; ++k)
			{
				to_be_inverted[j][k] = to_be_inverted[j][k] + factor*to_be_inverted[i][k];
			}
			for (int k = 0; k < NO_OF_CHOSEN_OBSERVATIONS_MOIST; ++k)
			{
				inv[j][k] = inv[j][k] + factor*inv[i][k];
			}
		}
	}
	
	#pragma omp parallel for
	for (int j = 0; j < NO_OF_CHOSEN_OBSERVATIONS_MOIST; ++j)
	{
		inv[NO_OF_CHOSEN_OBSERVATIONS_MOIST - 1][j] = inv[NO_OF_CHOSEN_OBSERVATIONS_MOIST - 1][j]/to_be_inverted[NO_OF_CHOSEN_OBSERVATIONS_MOIST - 1][NO_OF_CHOSEN_OBSERVATIONS_MOIST - 1];
	}
	to_be_inverted[NO_OF_CHOSEN_OBSERVATIONS_MOIST - 1][NO_OF_CHOSEN_OBSERVATIONS_MOIST - 1] = 1;
	
	/*
	Gaussian upwards
	----------------
	*/
	for (int i = NO_OF_CHOSEN_OBSERVATIONS_MOIST - 1; i >= 1; --i)
	{
		#pragma omp parallel for
		for (int j = i - 1; j >= 0; --j)
		{
			for (int k = 0; k < NO_OF_CHOSEN_OBSERVATIONS_MOIST; ++k)
			{
				inv[j][k] = inv[j][k] - to_be_inverted[j][i]*inv[i][k];
			}
		}
	}
	return 0;
}

int permute_lines_moist(double matrix[][NO_OF_CHOSEN_OBSERVATIONS_MOIST], int line_a, int line_b)
{
	/*
	Permutes line_a with line_b of matrix.
	*/
	double line_a_pre[NO_OF_CHOSEN_OBSERVATIONS_MOIST];
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS_MOIST; ++i)
	{
		line_a_pre[i] = matrix[line_a][i];
	}
	#pragma omp parallel for
	for (int i = 0; i < NO_OF_CHOSEN_OBSERVATIONS_MOIST; ++i)
	{
		matrix[line_a][i] = matrix[line_b][i];
		matrix[line_b][i] = line_a_pre[i];
	}
	return 0;
}









