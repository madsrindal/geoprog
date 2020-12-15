import numpy as np
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.pyplot as plt


def read_processed_file(filename):
    file = open(filename, 'r')

    resulting_dict = {}
    counter = 0

    for line in file:
        if counter > 0:
            line_values = line.split()
            lat, long = float(line_values[0]), float(line_values[1])
            n_value = float(line_values[2])
            resulting_dict[(lat, long)] = n_value
        counter += 1

    return resulting_dict


@np.vectorize
def get_n(phi, lmd, dictionary):
    return dictionary[(phi, lmd)]


def plot_pre_processed_height_file_world(egm=True):
    lon = np.arange(-180, 180.5, 0.5)
    lat = np.arange(-90, 90.5, 0.5)

    if egm:
        dictionary = read_processed_file('../../../geoprog/Part1/Results/egm_results_n100.txt')
        model_name = 'EGM2008'
    else:
        dictionary = read_processed_file('../../../geoprog/Part1/Results/ggm_results_nmax.txt')
        model_name = 'GGM03S'

    lon2d, lat2d = np.meshgrid(lon, lat)

    mycmap = plt.get_cmap('jet')

    data = get_n(lat2d, lon2d, dictionary)
    heights = plt.contourf(lon2d, lat2d, data, cmap=mycmap)

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines()
    ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter()
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    ax.contourf(lon, lat, data, cmap=mycmap)

    bar = plt.colorbar(heights, ax=ax, shrink=0.55)
    bar.ax.set_title('[m]')
    plt.title('Geoidal Heights for Global Gravity Model ' + model_name, fontdict={'fontsize': 13})

    plt.show()


def display_all_plots():
    print('Plotting global geoid model based on EGM2008 values...')
    plot_pre_processed_height_file_world(True)

    print('Plotting global geoid model based on GGM03S values...')
    plot_pre_processed_height_file_world(False)


if __name__ == '__main__':
    display_all_plots()
