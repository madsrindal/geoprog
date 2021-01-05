import numpy as np
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.pyplot as plt


# Function used to return a dictionary of the available sigma values related to a pair of latitude and longitude
def get_sigma_dict(year_num, gldas=True, true_sigma=True):
    # Check used to determine to return a dictionary of measured sigma values or the ones computed with coefficients
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

    # Looping through the file adding the non-fill (32767.000) values to a dictionary on format (lat, lon) = sigma
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


# A simple function returning an already calculated sigma. The function is vectorised in order to return arrays used
# to making grids in the plot.
@np.vectorize
def get_sigma(phi, lmd, dictionary):
    try:
        return dictionary[phi, lmd]
    except:
        return 0


# Function plotting the differnecees on a global map
def plot_pre_processed_sigma_file_world(year_num, gldas=True):
    lon = np.arange(0.5, 360.5, 1)
    lat = np.arange(-90.5, 90.5, 1)

    # Function used to save two dictionaries, one for true measured values and one for computed values.
    if gldas:
        true_dict = get_sigma_dict(year_num, gldas, True)
        calculated_dict = get_sigma_dict(year_num, gldas, False)
        model_name = 'GLDAS'
    else:
        true_dict = get_sigma_dict(year_num, gldas, True)
        calculated_dict = get_sigma_dict(year_num, gldas, False)
        model_name = 'ECCO'

    # Numbers used to determine the amount of error being inside a specific interval.
    # In this case, the maximum error value is set to +- 13 cmH20
    good_count = 0
    total_count = 0
    max_diff = 0
    limit_value = 13

    # Creates an empty dictionary to contain the error values for each available coordinate pair.
    diff_dict = {}
    for key in calculated_dict:
        if abs(true_dict[key] - calculated_dict[key]) < limit_value:
            diff_dict[key] = true_dict[key] - calculated_dict[key]
            good_count += 1

            # Check to find the max error within the maximum error value set above
            if abs(true_dict[key] - calculated_dict[key]) > max_diff:
                max_diff = abs(true_dict[key] - calculated_dict[key])

        total_count += 1

    # Prints for readability when running the code.
    print('Errors in range +-' + str(limit_value) + ' cmH20:')
    print('Error count included: ', good_count)
    print('Total error count: ', total_count)
    print('Percentage good values: ', (good_count/total_count)*100)

    # Creates the grid
    lon2d, lat2d = np.meshgrid(lon, lat)

    # Color palette
    mycmap = plt.get_cmap('jet')

    # Determine what kind of data to be visualised
    data = get_sigma(lat2d, lon2d, diff_dict)

    # Save the plotted heights
    heights = plt.contourf(lon2d, lat2d, data, cmap=mycmap)

    # Set the map axis, projection, axis labels and degree format.
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines()
    ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter()
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    # Creates the plot across the map
    ax.contourf(lon, lat, data, cmap=mycmap)

    # Set colorbar labels, map title and value unit
    bar = plt.colorbar(heights, ax=ax, shrink=0.55)
    bar.ax.set_title('[cmH20]')
    plt.title('Sigma differences based on ' + model_name + ' October_' + str(year_num), fontdict={'fontsize': 10})

    plt.show()


# Main function used to plot the differences of either GLDAS or ECCO, in a given year between 5 and 14.
def display_plots(year_num):
    print('Plotting differences between true sigma values and calculated sigma values for both ECCO and GLDAS file...')
    plot_pre_processed_sigma_file_world(year_num, True)
    print('----------------------------------------------')
    plot_pre_processed_sigma_file_world(year_num, False)


if __name__ == '__main__':
    display_plots(year_num=14)
