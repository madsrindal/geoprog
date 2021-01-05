from progressbar import progressbar
import math as mt
from Part1 import Sub_Functions_Part1 as SFP1

# Constants provided in the assignment
love_numbers = []
rho_w = 1000  # density of water [kg/m^3]
rho_ave = 5517  # average density of the Earth [kg/m^3]
a = 637100000  # radius of the Earth [m]


# Function used to modify a string in order to making it more suitable for rading.
def clean_string(string):
    return string.strip()[:-4]


# Function made to read and fix the provided love file
def read_love_file():
    counter = 0
    file = open('../Datafiles/Part2/love_number.txt', 'r')
    for line in file:
        line_values = line.split(' ')
        if counter >= 16:
            love_numbers.append(float(clean_string(line_values[1]))*10**-2)
        else:
            love_numbers.append(float(line_values[1]))
        counter += 1
    file.close()


# Function used to easily access the provided Love numbers.
def get_love_n(n):
    return love_numbers[n]


# Reads file on format (lon lat value) and return a dict on format {(lat, long): value}
def read_file(filename):
    file = open(filename, 'r')
    res_dict = {}

    # Does not append fill values being 32767.000
    for line in file:
        line_values = line.split()
        lon, lat = float(line_values[0]), float(line_values[1])
        value = float(line_values[2])
        if value != float('32767.000'):
            res_dict[(lat, lon)] = value
    file.close()

    return res_dict


# Function used to compute a single coefficient (either r or q) for a specific pairing of n and m.
def get_coefficient(n, m, dictionary, r_bar=True):
    constant = 1/(4*mt.pi) * ((1 + get_love_n(n))/(2*n + 1)) * ((3*rho_w)/(a * rho_ave))

    # Implementation of the summing formula, ranging from available phi and lamda values.
    tot_sum = 0
    for key in dictionary:
        if 0 <= mt.radians(key[1]) <= 2*mt.pi and -mt.pi/2 <= mt.radians(key[0]) <= mt.pi/2:

            # check of what kind of coefficient is calculated (R_bar or q_bar)
            if r_bar:
                res = dictionary[key] * mt.cos(m * mt.radians(key[1]))
            else:
                res = dictionary[key] * mt.sin(m * mt.radians(key[1]))

            tot_sum += res * SFP1.get_p_bar(n, m, mt.radians(key[0]), {}) * ((mt.pi / 180) ** 2) * mt.cos(mt.radians(key[0]))

    return constant*tot_sum


# Function usedd to calculate every coefficient for a specific file up to the max degree decided by n_max.
def get_all_coefficients(n_max, year_num, gldas=True):
    # Check to determine what kind of data and file is being computed.
    if gldas:
        dictionary = gldas_list[year_num-5]
        file = open('../Part2/Results/GLDAS/GLDAS_Coefficients_Oct_' + str(year_num) + '_test.txt', 'w')
    else:
        dictionary = ecco_list[year_num-5]
        file = open('../Part2/Results/ECCO/ECCO_Coefficients_Oct_' + str(year_num) + '_test.txt', 'w')

    file.write('N\tM\tR_coefficient\tQ_coefficient\n')

    # Saves the computed coefficients to a text file on format N M R Q.
    for n in range(0, n_max+1):
        for m in range(0, n+1):
            r_bar = get_coefficient(n, m, dictionary, True)
            q_bar = get_coefficient(n, m, dictionary, False)
            file.write(str(n) + '\t' + str(m) + '\t' + str(r_bar) + '\t' + str(q_bar) + '\n')

    file.close()


# Creates two list that will contain 10 dictionaries each, one for each year from 2005-2014
ecco_list = []
gldas_list = []


#  Function used to compute a value for sigma based on a specific location and data type.
def get_sigma(phi, lmd, n_max, year_number, dictionary, gldas=True, rq_dict=True):
    phi = mt.radians(phi)
    lmd = mt.radians(lmd)
    constant = ((a * rho_ave)/3)
    tot_sum = 0

    # The implemented formula adding every combination of n and m. The function is structured in a way making it
    # possible to compute the value in two different ways. Either from scratch or with an already computed dictionary
    # of relevant coefficients.
    for n in range(0, n_max + 1):
        part_sum = 0
        for m in range(0, n + 1):
            res = (2*n + 1)/(1 + get_love_n(n))

            if rq_dict:
                r_bar = dictionary[(n, m)][0]
                q_bar = dictionary[(n, m)][1]
            else:
                if gldas:
                    r_bar = get_coefficient(n, m, gldas_list[year_number-5], True)
                    q_bar = get_coefficient(n, m, gldas_list[year_number-5], False)
                else:
                    r_bar = get_coefficient(n, m, ecco_list[year_number-5], True)
                    q_bar = get_coefficient(n, m, ecco_list[year_number-5], False)

            part_sum += res * (r_bar * mt.cos(m * lmd) + q_bar * mt.sin(m * lmd)) * SFP1.get_p_bar(n, m, mt.sin(phi), {})
        tot_sum += part_sum
    return constant*tot_sum
