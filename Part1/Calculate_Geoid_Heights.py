from time import time
import math as mt
from progressbar import progressbar
import numpy as np
from Part1 import Sub_Functions_Part1 as SFP1


# Main function used for computation and saving of gravimetric geoidal heights.
def calculate_geoidal_heights_world(n_max, egm=True):
    start_time = time()

    # Check to determine what model to calculate and where to save the results.
    print('Creating a dictionary containing necessary coefficients for each pair of n and m...')
    if egm:
        data_dict = SFP1.read_file(True)
        file_name = 'egm_results_test.txt'
        file_path = '../../geoprog/Part1/Results/' + file_name
        model_name = 'EGM2008'

    else:
        data_dict = SFP1.read_file(False)
        file_name = 'ggm_results_test.txt'
        file_path = '../../geoprog/Part1/Results/' + file_name
        model_name = 'GGM03S'

    # Creation of matricces containing values of r and q for each pair om n and m up to n_max.
    print('Creating r_ and q_matrices...')
    r_matrix, q_matrix = SFP1.get_rq_bar_matrix(data_dict, n_max)

    # Creates an empty dictionary to contain height values for each pair of phi and lmd
    n_dict = {}

    # Computing the height values
    print('Calculating geoidal heights all over the world using EGM2008 values and writing these to file...')
    # Looping through phi on top as the p values are the same for each phi. In this way, the construction of the full
    # p_matrix is limited to the least amount of constructions, and hence being more time effective.
    for phi in progressbar(np.arange(-90, 90.5, 0.5)):

        # Creating the full p_matrix consisting of p_values for each pair of n and m for the specific phi.
        p_matrix = np.zeros((n_max + 1) ** 2).reshape(n_max + 1, n_max + 1)
        p_value_dict = {}
        for n in range(0, n_max + 1):
            for m in range(0, n_max + 1):
                p_matrix[n][m] = SFP1.get_p_bar(n, m, mt.radians(phi), p_value_dict)

        # Looping through each lambda in order to compute the height. A copy of the r and q matrices are being made in
        # order to secure the original matrices.
        for lmd in np.arange(-180, 180.5, 0.5):
            r = np.copy(r_matrix)
            q = np.copy(q_matrix)
            # Performing trigonometric operations on the whole matrices in order to save computational time.
            for m in range(0, n_max + 1):
                r[:, m] *= mt.cos(m * mt.radians(lmd))
                q[:, m] *= mt.sin(m * mt.radians(lmd))
            # Compute the height value for the specific pair of phi and lambda.
            n_dict[(phi, lmd)] = SFP1.get_n_grv(lmd, n_max, p_matrix, r, q)

    # Saving the computed dictionary of height values to file. This helps making the geoidal maps faster,
    # and not having to run the entire computations every time.
    file = open(file_path, 'w')
    file.write('LAT\tLON\tGeoidal Height\n')
    for key in n_dict:
        file.write(str(key[0]) + '\t' + str(key[1]) + '\t' + str(n_dict[key]) + '\n')
    file.close()

    print('Calculated geoidal heights using ' + str(model_name) + ' values of degree n=' + str(n_max) + ' in ' +
          str(time() - start_time) + ' seconds. The results are saved as ' + file_name + ' found in the results folder')
