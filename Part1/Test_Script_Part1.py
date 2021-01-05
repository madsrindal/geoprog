from Part1 import Calculate_Geoid_Heights as CGH
from Part1.Plots import Plot_geoid_model as Plot


def main(n_max, egm=True):
    if egm:
        model = 'EGM2008'
    else:
        model = 'GGM03S'

    CGH.calculate_geoidal_heights_world(n_max=n_max, egm=egm)

    print('Creating a geoid model based based on ' + model + '...')

    Plot.plot_pre_processed_height_file_world(egm=egm, test_script=True)

    print('The presented geoid model is made of coefficients up to a degree, n=' + str(n_max) +
          ' resulting in a less accurate model, but the most common shapes are already starting to appear.\n')

    print('This script was written to present the whole process in short, '
          'and may be altered by the user if the degree n_max or model is changed. '
          'Please note that higher degrees results in significant longer run time.')


if __name__ == '__main__':
    main(n_max=20, egm=True)
