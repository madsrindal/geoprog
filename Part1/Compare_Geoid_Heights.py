from Part1 import Gravimetric_geoid_model as Ggm
from Part1 import p_bar_testing as pbt
from progressbar import progressbar
import numpy as np


def calculate_n_geometric():
    file = open("../../geoprog/Datafiles/Part1/GNSS_data.txt")

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


def main():
    print("Starting...")

    print('Calculating gravimetric heights based on GNSS locations...')
    grav_dict = calculate_n_gravimetric(False)

    print('Calculating geometric heights based on provided GNSS data...')
    geo_dict = calculate_n_geometric()

    diffs = []

    print('Calculating differences between gravimetric and geometric heights...')
    for key in grav_dict:
        diffs.append(grav_dict[key] - geo_dict[key])

    abs_diffs = []
    for i in diffs:
        abs_diffs.append(abs(i))

    print('---------------------------------------------------------------------')

    print('Maximum difference: ' + str(max(abs_diffs)))
    print('Minimum difference: ' + str(min(abs_diffs)))
    print('Mean of differences: ' + str(np.mean(diffs)))
    print('Standard deviation: ' + str(np.std(diffs)))


if __name__ == '__main__':
    save_n_gravimetric_gnss(True)
