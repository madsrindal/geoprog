import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from Part1 import Compare_Geoid_Heights as Cgh


@np.vectorize
def get_n(phi, lmd, dictionary):
    return dictionary[phi, lmd]


def get_degree_interval(key_values=True):
    n_dict = Cgh.calculate_n_geometric()
    lat_list = []
    lon_list = []

    for key in n_dict:
        lat_list.append(key[0])
        lon_list.append(key[1])

    lat_min = min(lat_list)
    lat_max = max(lat_list)
    lon_min = min(lon_list)
    lon_max = max(lon_list)

    if key_values:
        return lat_min, lat_max, lon_min, lon_max
    else:
        return lat_list, lon_list


def plot_height_differences_fail(egm=True):
    lat_min, lat_max, lon_min, lon_max = get_degree_interval()

    lon = np.arange(lon_min, lon_max)
    lat = np.arange(lat_min, lat_max)

    if egm:
        model_name = 'EGM2008'
        diff_dict = Cgh.get_diff_dict(True)
    else:
        model_name = 'GGM03S'
        diff_dict = Cgh.get_diff_dict(False)

    lon2d, lat2d = np.meshgrid(lon, lat)
    mycmap = plt.get_cmap('jet')

    data = get_n(lat2d, lon2d, diff_dict)

    heights = plt.contourf(lon2d, lat2d, data, cmap=mycmap)

    ax = plt.axes(projection=ccrs.EuroPP())
    extent = [lon_min, lon_max, lat_min, lat_max]
    ax.set_extent(extent)
    ax.coastlines()
    # ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.EuroPP())
    # ax.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.EuroPP())
    lon_formatter = LongitudeFormatter()
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    ax.contourf(lon, lat, data, cmap=mycmap)

    bar = plt.colorbar(heights, ax=ax, shrink=0.55)
    bar.ax.set_title('[m]')
    plt.title('Geometric and gravimetric height differences based on ' + model_name, fontdict={'fontsize': 10})

    plt.show()


def plot_height_differences(egm=True):
    lat_min, lat_max, lon_min, lon_max = get_degree_interval(True)

    lat_list, lon_list = get_degree_interval(False)

    if egm:
        model_name = 'EGM2008'
        diff_dict = Cgh.get_diff_dict(True)
    else:
        model_name = 'GGM03S'
        diff_dict = Cgh.get_diff_dict(False)

    values = []

    for i in range(0, len(lat_list)):
        values.append([diff_dict[lat_list[i], lon_list[i]], lat_list[i], lon_list[i]])

    mycmap = plt.get_cmap('jet')

    # data = get_n(lat_list, lon_list, diff_dict)

    heights = plt.contourf(lat_list, lon_list, values, cmap=mycmap)

    ax = plt.axes(projection=ccrs.EuroPP())
    extent = [lon_min, lon_max, lat_min, lat_max]
    ax.set_extent(extent)
    ax.coastlines()
    # ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.EuroPP())
    # ax.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.EuroPP())
    lon_formatter = LongitudeFormatter()
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    ax.contourf(lon_list, lat_list, values, cmap=mycmap)

    bar = plt.colorbar(heights, ax=ax, shrink=0.55)
    bar.ax.set_title('[m]')
    plt.title('Geometric and gravimetric height differences based on ' + model_name, fontdict={'fontsize': 10})

    plt.show()


def display_plot():
    # print('Plotting differences between true sigma values and calculated sigma values...')
    plot_height_differences(True)


if __name__ == '__main__':
    display_plot()
