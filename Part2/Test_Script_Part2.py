from Part2 import Calculate_Coefficients_and_Sigmas as CCS
from Part2 import Sub_Functions_Part2 as SFP2


def main(n_max, year_num, gldas=True):
    print('The following values has been chosen for this test: ')
    print('n_max: ', n_max)
    print('Year number: ', year_num)

    if gldas:
        file_type = 'GLDAS'
        phi = -54.5
        lmd = 289.5
    else:
        file_type = 'ECCO'
        phi = -72.5
        lmd = 170.5

    print('Type of data to be used: ', file_type)
    print('---------------------------------------')

    CCS.pre_processing()

    if gldas:
        true_sigma_dict = SFP2.gldas_list[year_num-5]
    else:
        true_sigma_dict = SFP2.ecco_list[year_num-5]

    print('Computing coefficients up to degree n_max=' + str(n_max) + ' for the ' + file_type + ' file from October_'
          + str(year_num) + ' and saving them to the results folder as ' + file_type + '_Coefficients_Oct_'
          + str(year_num) + '_test.txt...\n')

    print('These values will be identical to the ones already computed located in the same folder, but with "test" '
          'removed from the file name...\n')

    SFP2.get_all_coefficients(n_max, year_num, gldas)

    print('Calculating a sigma value based on the computed coefficients at latitude=' + str(phi) + ' and longitude=' +
          str(lmd) + '...\n')
    computed_sigma = SFP2.get_sigma(phi, lmd, n_max, year_num, {}, gldas, False)
    print('The computed sigma divided by rho_water for comparison purposes is: ', computed_sigma/1000)

    measured_sigma = true_sigma_dict[phi, lmd]

    print('The measured value for sigma divided by rho_water at the same location is: ', measured_sigma)

    print('Most likely, this sigma is way off the measured value, and this is because of the low degree, n, chosen at '
          'the start. However, you have now tested the functions and hopefully with great success!\n')


if __name__ == '__main__':
    main(10, 14, gldas=True)

