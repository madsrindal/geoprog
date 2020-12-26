from Part1 import Gravimetric_geoid_model as Ggm
from Part1 import p_bar_testing as pbt
from progressbar import progressbar
import numpy as np


def calculate_n_geometric():
    file = open("../../Datafiles/Part1/GNSS_data.txt")

    n_dict = {}

    for line in file:
        line_values = line.split()
        lat, long = float(line_values[1]), float(line_values[2])
        h_ort, h_ell = float(line_values[3]), float(line_values[4])
        n_dict[(lat, long)] = h_ell - h_ort

    file.close()

    return n_dict


def calculate_n_gravimetric(egm=True):
    file = open("../../geoprog/Datafiles/Part1/GNSS_data.txt")

    n_dict = {}
    gnss_dict = {}

    print('Creating value dictionaries based on EGM2008 or GGM03S data...')
    if egm:
        value_dict = Ggm.read_file(True)
        n_max = 2050
    else:
        value_dict = Ggm.read_file(False)
        n_max = 180

    print('Reading GNSS file, creating a dictionary on format (phi, lmd) = geoid height...')
    for line in file:
        line_values = line.split()
        phi, lmd = float(line_values[1]), float(line_values[2])
        n = float(line_values[5])
        gnss_dict[(phi, lmd)] = n
    file.close()

    print('Calculating geoid height based on EGM or GGM model, saving them to a dictionary on format (phi, lmd) = N...')
    for key in progressbar(gnss_dict):
        p_value_dict = {}
        n_dict[(key[0], key[1])] = pbt.get_n_grv(key[0], key[1], value_dict, n_max, p_value_dict)

    return n_dict


def save_n_gravimetric_gnss(egm=True):

    n_dict = calculate_n_gravimetric(egm)

    if egm:
        file = open('../../geoprog/Part1/Results/Calculated_GNSS_heights_EGM2008_1', 'w')
    else:
        file = open('../../geoprog/Part1/Results/Calculated_GNSS_heights_GGM03S_1', 'w')

    file.write('LAT\tLON\tGeoidal Height\n')

    for key in n_dict:
        file.write(str(key[0]) + '\t' + str(key[1]) + '\t' + str(n_dict[key]) + '\n')

    file.close()


def get_n_gravimetric(egm=True):

    n_dict = {}

    if egm:
        path = '../../Part1/Results/Calculated_GNSS_heights_EGM2008'
    else:
        path = '../../Part1/Results/Calculated_GNSS_heights_GGM03S'

    file = open(path, 'r')
    counter = 0

    for line in file:
        if counter > 0:
            line_values = line.split()
            lat, lon = float(line_values[0]), float(line_values[1])
            n = float(line_values[2])
            n_dict[(lat, lon)] = n
        counter += 1
    file.close()

    return n_dict


def get_diff_dict(egm=True):

    geo_dict = calculate_n_geometric()
    grav_dict = get_n_gravimetric(egm)

    diff_dict = {}

    for key in geo_dict:
        diff_dict[key] = geo_dict[key] - grav_dict[key]

    return diff_dict


def get_statistics(dictionary):

    diffs = []
    abs_diffs = []

    for key in dictionary:
        diffs.append(dictionary[key])
        abs_diffs.append(abs(dictionary[key]))

    max_diff = max(abs_diffs)
    min_diff = min(abs_diffs)
    mean_diff = np.mean(abs_diffs)
    std = np.std(diffs)

    return max_diff, min_diff, mean_diff, std


def main():

    diff_dict_egm = get_diff_dict(egm=True)

    max_diff_egm, min_diff_egm, mean_diff_egm, std_egm = get_statistics(diff_dict_egm)

    print('Maximum difference: ', max_diff_egm)
    print('Minimum difference: ', min_diff_egm)
    print('Mean of differences: ', mean_diff_egm)
    print('Standard deviation: ', std_egm)

    print('-------------------------------------------------')

    diff_dict_ggm = get_diff_dict(egm=False)

    max_diff_ggm, min_diff_ggm, mean_diff_ggm, std_ggm = get_statistics(diff_dict_ggm)

    print('Maximum difference: ', max_diff_ggm)
    print('Minimum difference: ', min_diff_ggm)
    print('Mean of differences: ', mean_diff_ggm)
    print('Standard deviation: ', std_ggm)


if __name__ == '__main__':
    main()
