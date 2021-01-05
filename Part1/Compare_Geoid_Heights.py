from Part1 import Sub_Functions_Part1 as SFP1
from progressbar import progressbar
import numpy as np


# Function used to calculate the geometric heights from provided GNSS leveling data.
# Returns a dictionary of geometric height values for each available pair of phi and lambda.
def calculate_n_geometric(gnss_plot=False):
    if gnss_plot:
        file = open("../../Datafiles/Part1/GNSS_data.txt")
    else:
        file = open("../../geoprog/Datafiles/Part1/GNSS_data.txt")

    n_dict = {}

    # Looping through the file, computing the height and saving them to a dictionary
    for line in file:
        line_values = line.split()
        lat, long = float(line_values[1]), float(line_values[2])
        h_ort, h_ell = float(line_values[3]), float(line_values[4])
        n_dict[(lat, long)] = h_ell - h_ort

    file.close()

    return n_dict


# Function used to calculate the gravimetric heights for the locations provided in the GNSS data file returned as a dict
def calculate_n_gravimetric(gnss_plot=False, egm=True):
    if gnss_plot:
        file = open("../../Datafiles/Part1/GNSS_data.txt")
    else:
        file = open("../../geoprog/Datafiles/Part1/GNSS_data.txt")

    n_dict = {}
    gnss_dict = {}

    # Check to determine what kind of model to ues in the coputations, and the degree n of these.
    print('Creating value dictionaries based on EGM2008 or GGM03S data...')
    if egm:
        value_dict = SFP1.read_file(True)
        n_max = 2050
    else:
        value_dict = SFP1.read_file(False)
        n_max = 180

    # Creates a dictionary to obtain every pair of available pairings of phi and lambda.
    print('Reading GNSS file, creating a dictionary on format (phi, lmd) = geoid height...')
    for line in file:
        line_values = line.split()
        phi, lmd = float(line_values[1]), float(line_values[2])
        n = float(line_values[5])
        gnss_dict[(phi, lmd)] = n
    file.close()

    # Loop through each available key, and computes the gravimetric height at the given location, saving these in a dict
    # The p_value_dictionary is emptied for every computation in order for the code not to crash because of memory error
    print('Calculating geoid height based on EGM or GGM model, saving them to a dictionary on format (phi, lmd) = N...')
    for key in progressbar(gnss_dict):
        p_value_dict = {}
        n_dict[(key[0], key[1])] = SFP1.get_n_grv_modified(key[0], key[1], value_dict, n_max, p_value_dict)

    return n_dict


# Function used to save a dictionary to a text file
def save_n_gravimetric_gnss(egm=True):

    # Computes the dictionary of gravimetric height values.
    n_dict = calculate_n_gravimetric(egm)

    # Determine where to save the resultant file.
    if egm:
        file = open('../../geoprog/Part1/Results/Calculated_GNSS_heights_EGM2008_1', 'w')
    else:
        file = open('../../geoprog/Part1/Results/Calculated_GNSS_heights_GGM03S_1', 'w')

    # Writes the values to file on format Latitude - Longitude - Height
    file.write('LAT\tLON\tGeoidal Height\n')

    for key in n_dict:
        file.write(str(key[0]) + '\t' + str(key[1]) + '\t' + str(n_dict[key]) + '\n')

    file.close()


# Function used to read a pre processed file of gravimetric heights computed from the GNSS file and return the values
# as a dictionary
def get_n_gravimetric(gnss_plot=False, egm=True):

    n_dict = {}

    # Determines where to collect the pre processed values
    if gnss_plot:
        if egm:
            path = '../../Part1/Results/Calculated_GNSS_heights_EGM2008'
        else:
            path = '../../Part1/Results/Calculated_GNSS_heights_GGM03S'
    else:
        if egm:
            path = '../../geoprog/Part1/Results/Calculated_GNSS_heights_EGM2008'
        else:
            path = '../../geoprog/Part1/Results/Calculated_GNSS_heights_GGM03S'

    file = open(path, 'r')
    counter = 0

    # Loops through the text file, saving the values to a dictionary
    for line in file:
        if counter > 0:
            line_values = line.split()
            lat, lon = float(line_values[0]), float(line_values[1])
            n = float(line_values[2])
            n_dict[(lat, lon)] = n
        counter += 1
    file.close()

    return n_dict


# Function used to create and return a dictionary containing differences in values between geometric and gravimetric
# height values
def get_diff_dict(gnss_plot=False, egm=True):

    # Creates the two dictionaries
    geo_dict = calculate_n_geometric(gnss_plot=gnss_plot)
    grav_dict = get_n_gravimetric(gnss_plot=gnss_plot, egm=egm)

    diff_dict = {}

    # Creates a dictionary containing differences
    for key in geo_dict:
        diff_dict[key] = geo_dict[key] - grav_dict[key]

    return diff_dict


# Function used to create and return relevant statistics of the differences
def get_statistics(dictionary):

    # Two lists to contain the differences, with one being the difference in absolute value
    diffs = []
    abs_diffs = []

    # Appending the values to the empty lists
    for key in dictionary:
        diffs.append(dictionary[key])
        abs_diffs.append(abs(dictionary[key]))

    # Saving the extremum values
    max_diff = max(abs_diffs)
    min_diff = min(abs_diffs)
    mean_diff = np.mean(abs_diffs)
    std = np.std(diffs)

    return max_diff, min_diff, mean_diff, std


# Main function used to present the statistics for both computational methods.
def main():

    print('Creating statistics related to the comparison between gravimetric and geometric heights based on EGM2008...')
    diff_dict_egm = get_diff_dict(egm=True)

    max_diff_egm, min_diff_egm, mean_diff_egm, std_egm = get_statistics(diff_dict_egm)

    print('Maximum difference: ', max_diff_egm)
    print('Minimum difference: ', min_diff_egm)
    print('Mean of differences: ', mean_diff_egm)
    print('Standard deviation: ', std_egm)

    print('-----------------------------------------------------------------------------------------------------------')

    print('Creating statistics related to the comparison between gravimetric and geometric heights based on GGM03S...')
    diff_dict_ggm = get_diff_dict(egm=False)

    max_diff_ggm, min_diff_ggm, mean_diff_ggm, std_ggm = get_statistics(diff_dict_ggm)

    print('Maximum difference: ', max_diff_ggm)
    print('Minimum difference: ', min_diff_ggm)
    print('Mean of differences: ', mean_diff_ggm)
    print('Standard deviation: ', std_ggm)


if __name__ == '__main__':
    main()
