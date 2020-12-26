import numpy as np
from time import time
from progressbar import progressbar
from Part2 import Stokes_coeff as Sc


def pre_processing():
    print("Program starting...")

    print('Reading file containing Love numbers')
    Sc.read_love_file()

    print('Creating two lists of GLDAS and ECCO dictionaries, one dictionary for each GLDAS and ECCO file...')
    for year_num in range(5, 15):
        gldas_dict = Sc.read_file('../../geoprog/Datafiles/Part2/GLDAS/GLDAS_Oct_' + str(year_num) + '.txt')
        ecco_dict = Sc.read_file('../../geoprog/Datafiles/Part2/ECCO/ECCO_Oct_' + str(year_num) + '.txt')
        Sc.gldas_list.append(gldas_dict)
        Sc.ecco_list.append(ecco_dict)

    print('---------------------------------------------------------------------------------------------------------')


def calculate_all_coefficients():
    start_time = time()
    print('Calculating coefficients for all GLDAS files from 2005 to 2014...')
    print('The resulting files containing these coefficients can later be found in Part2/Results/GLDAS')
    for year_num in progressbar(range(5, 15)):
        Sc.get_all_coefficients(100, year_num, True)
    print('All GLDAS coefficients calculated in ' + str(time()-start_time) + ' seconds.')

    print('---------------------------------------------------------------------------------------------------------')

    second_time = time()
    print('Calculating coefficients for all ECCO files from 2005 to 2014...')
    print('The resulting files containing these coefficients can later be found in Part2/Results/ECCO')
    for year_num in progressbar(range(5, 15)):
        Sc.get_all_coefficients(100, year_num, False)
    print('All ECCO coefficients calculated in ' + str(time() - second_time) + ' seconds.')

    print('---------------------------------------------------------------------------------------------------------')

    print('Program finished in a total of ' + str(time()-start_time) + ' seconds.')


def calculate_sigmas_for_one_file(year_num, gldas=True):
    if gldas:
        true_file = '../../geoprog/Datafiles/Part2/GLDAS/GLDAS_Oct_' + str(year_num) + '.txt'
        coef_file = '../../geoprog/Part2/Results/GLDAS/GLDAS_Coefficients_Oct_' + str(year_num) + '.txt'
        res_file = '../../geoprog/Part2/Results/GLDAS/GLDAS_totH20_Oct_' + str(year_num) + '.txt'
    else:
        true_file = '../../geoprog/Datafiles/Part2/ECCO/ECCO_Oct_' + str(year_num) + '.txt'
        coef_file = '../../geoprog/Part2/Results/ECCO/ECCO_Coefficients_Oct_' + str(year_num) + '.txt'
        res_file = '../../geoprog/Part2/Results/ECCO/ECCO_OBP_Oct_' + str(year_num) + '.txt'

    f = open(coef_file, 'r')
    res_dict = {}
    counter = 0
    for line in f:
        if counter > 0:
            line_values = line.split()
            n, m = float(line_values[0]), float(line_values[1])
            r, q = float(line_values[2]), float(line_values[3])
            res_dict[(n, m)] = r, q
        counter += 1
    f.close()

    key_dict = Sc.read_file(true_file)

    f = open(res_file, 'w')
    f.write('LON\tLAT\tVALUE\n')
    for key in key_dict:
        sigma = Sc.get_sigma(key[0], key[1], 100, year_num, res_dict, gldas=True, rq_dict=True)
        f.write(str(key[1]) + '\t' + str(key[0]) + '\t' + str(sigma/1000) + '\n')
    f.close()


def calculate_sigmas_for_every_file():
    time_one = time()
    print('Calculating Sigma over Rho_water for every phi and lamda having a known value in original GLDAS files...')
    print('The resulting files containing calculated total H20 values can later be found in Part2/Results/GLDAS')

    for year_num in progressbar(range(5, 15)):
        calculate_sigmas_for_one_file(year_num, True)

    print('Total h20 values for every year from 2005 to 2014 calculated in ' + str(time() - time_one) + ' seconds')

    print('---------------------------------------------------------------------------------------------------------')

    time_two = time()
    print('Calculating Sigma over Rho_water for every phi and lamda having a known value in original ECCO files...')
    print('The resulting files containing calculated Ocean Bottom Pressure values can later be found in '
          'Part2/Results/ECCO')

    for year_num in progressbar(range(5, 15)):
        calculate_sigmas_for_one_file(year_num, False)

    print('OBP values for every year from 2005 to 2014 calculated in ' + str(time() - time_two) + ' seconds')

    print('---------------------------------------------------------------------------------------------------------')

    print('Program finished in a total of ' + str(time() - time_one) + ' seconds.')


def compare_sigmas(year_num, gldas=True):
    if gldas:
        file = open('../../geoprog/Part2/Results/GLDAS/GLDAS_totH20_Oct_' + str(year_num) + '.txt')
        true_sigma_dict = Sc.gldas_list[year_num-5]
        name = 'GLDAS'
    else:
        file = open('../../geoprog/Part2/Results/ECCO/ECCO_OBP_Oct_' + str(year_num) + '.txt')
        true_sigma_dict = Sc.ecco_list[year_num - 5]
        name = 'ECCO'

    calculated_sigma_dict = {}

    counter = 0
    for line in file:
        if counter > 0:
            line_values = line.split()
            lat, lon = float(line_values[1]), float(line_values[0])
            value = float(line_values[2])
            calculated_sigma_dict[(lat, lon)] = value
        counter += 1

    diffs = []

    for key in calculated_sigma_dict:
        diffs.append(calculated_sigma_dict[key] - true_sigma_dict[key])

    abs_diffs = []
    for i in diffs:
        abs_diffs.append(abs(i))

    print('Statistics regarding the comparison between given sigma and calculated sigma based on file ' + name + '_Oct_'
          + str(year_num))
    print('---------------------------------------------------------------------------------------------------------')
    print('Standard deviation: ' + str(np.std(diffs)))
    print('Mean of differences: ' + str(np.mean(abs_diffs)))
    print('Maximum difference: ' + str(max(abs_diffs)))
    print('Minimum difference: ' + str(min(abs_diffs)))


if __name__ == '__main__':
    pre_processing()
    # calculate_all_coefficients()
    # calculate_sigmas_for_every_file()
    compare_sigmas(5, True)
