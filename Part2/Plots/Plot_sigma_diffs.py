import numpy as np
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.pyplot as plt


def get_sigma_dict(year_num, gldas=True, true_sigma=True):
    if true_sigma:
        counter = 0
        if gldas:
            path = '../../../geoprog/Datafiles/Part2/GLDAS/GLDAS_Oct_' + str(year_num) + '.txt'
        else:
            path = '../../../geoprog/Datafiles/Part2/ECCO/ECCO_Oct_' + str(year_num) + '.txt'
    else:
        counter = -1
        if gldas:
            path = '../../../geoprog/Part2/Results/GLDAS/GLDAS_totH20_Oct_' + str(year_num) + '.txt'
        else:
            path = '../../../geoprog/Part2/Results/ECCO/ECCO_OBP_Oct_' + str(year_num) + '.txt'

    sigma_dict = {}

    file = open(path)

    for line in file:
        if counter >= 0:
            line_values = line.split()
            lon, lat = float(line_values[0]), float(line_values[1])
            value = float(line_values[2])
            if value != 32767.000:
                sigma_dict[(lat, lon)] = value
        counter += 1
    file.close()
    return sigma_dict


@np.vectorize
def get_sigma(phi, lmd, dictionary):
    try:
        return dictionary[phi, lmd]
    except:
        return 0


def plot_pre_processed_sigma_file_world(year_num, gldas=True):
    lon = np.arange(0.5, 360.5, 1)
    lat = np.arange(-90.5, 90.5, 1)

    if gldas:
        true_dict = get_sigma_dict(year_num, gldas, True)
        calculated_dict = get_sigma_dict(year_num, gldas, False)
        model_name = 'GLDAS'
    else:
        true_dict = get_sigma_dict(year_num, gldas, True)
        calculated_dict = get_sigma_dict(year_num, gldas, False)
        model_name = 'ECCOS'

    diff_dict = {}

    good_count = 0
    total_count = 0

    for key in calculated_dict:
        if abs(true_dict[key] - calculated_dict[key]) < 100:
            good_count += 1
        total_count += 1

    print('Good count: ', good_count)
    print('Total count: ', total_count)
    print('Percentage good values: ', (good_count/total_count)*100)

    for key in calculated_dict:
        if abs(true_dict[key] - calculated_dict[key]) < 1000:
            diff_dict[key] = true_dict[key] - calculated_dict[key]

    lon2d, lat2d = np.meshgrid(lon, lat)
    mycmap = plt.get_cmap('jet')

    data = get_sigma(lat2d, lon2d, diff_dict)

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
    bar.ax.set_title('[cmH20]')
    plt.title('Sigma differences based on ' + model_name + ' October_' + str(year_num), fontdict={'fontsize': 10})

    plt.show()


def display_plot():
    print('Plotting differences between true sigma values and calculated sigma values...')
    plot_pre_processed_sigma_file_world(8, True)


if __name__ == '__main__':
    display_plot()
