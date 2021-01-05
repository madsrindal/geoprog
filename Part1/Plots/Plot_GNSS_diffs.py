import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from Part1 import Compare_Geoid_Heights as Cgh


# Function used to return the max and min value of latitude and longitude related to provided GNSS data.
# Can also be used to return two lists containing every value of latitude and longitude.
def get_degree_interval(key_values=True):

    # Creates the dictionary containing the available values for lat and long, and two empty lists
    n_dict = Cgh.calculate_n_geometric(gnss_plot=True)
    lat_list = []
    lon_list = []

    # Appends every value of lat and long to the empty lists
    for key in n_dict:
        lat_list.append(key[0])
        lon_list.append(key[1])

    # Saves the extremum values
    lat_min = min(lat_list)
    lat_max = max(lat_list)
    lon_min = min(lon_list)
    lon_max = max(lon_list)

    # Returns the desired values being either two lists or the extremum values
    if key_values:
        return lat_min, lat_max, lon_min, lon_max
    else:
        return lat_list, lon_list


# Function used to plot the height difference between GNSS leveling data and computed gravimetric geoidal heights
def plot_height_differences(egm=True):
    # Saves the extremum values with the use of the sub function above
    lat_min, lat_max, lon_min, lon_max = get_degree_interval(True)

    # Creates a dictionary depending on the model used and saves the name (EGM or GGM)
    if egm:
        model_name = 'EGM2008'
        diff_dict = Cgh.get_diff_dict(gnss_plot=True, egm=True)
    else:
        model_name = 'GGM03S'
        diff_dict = Cgh.get_diff_dict(gnss_plot=True, egm=False)

    # Chosen color palette
    mycmap = plt.get_cmap('jet')

    # Creates the map projection, axis labels and step size, and set the degree to lat and long.
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    step_size_lon = (lon_max-lon_min)/3
    step_size_lat = (lat_max-lat_min)/3
    ax.set_xticks([lon_min, lon_min + step_size_lon, lon_min + step_size_lon*2, lon_max], crs=ccrs.PlateCarree())
    ax.set_yticks([lat_min, lat_min + step_size_lat, lat_min + step_size_lat*2, lat_max], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter()
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    # Creates two emtpy lists to contain x and y coordinates of the provided GNSS data locations
    # Creates two empty lists to contain height value and size of location point.
    x = []
    y = []
    heights = []
    point_sizes = []

    # Loop to add values to the created lists above
    for key in diff_dict.keys():
        x.append(key[1])
        y.append(key[0])
        heights.append(diff_dict[key])
        point_sizes.append(abs(diff_dict[key] * 5))  # Multiplied by five to make the dots bigger for visual reasons

    # Set the colorbar lables, map title and creates the plot
    plt.scatter(x, y, alpha=1, c=heights, s=point_sizes, cmap=mycmap)
    bar = plt.colorbar(shrink=0.50)
    bar.ax.set_title('[m]')
    plt.title('Geometric and gravimetric height differences based on ' + model_name, fontdict={'fontsize': 10})
    plt.show()


# Main function used to display a plot for each model (EGM2008 and GGM03S)
def display_plot():
    plot_height_differences(True)
    plot_height_differences(False)


if __name__ == '__main__':
    display_plot()
