import math as mt
import numpy as np
import sys

# Constants
a = 6378137.0000
b = 6356752.3141
f_inv = 298.257222101
e_sqd = 0.006694380023
R = 6371000.7900
GM = 3986005 * 10**8
gamma = 9.81
j_2 = 0.108263*(10**(-2))

sys.setrecursionlimit(10**6)


# Takes a file written in format [Key, N, M, C, S, Sigma_C, Sigma_S] separated by tab, starting on N=2, M=0
def create_dict(filename):
    res = {}
    file = open(filename, "r")
    for line in file:
        line_values = line.split()
        res[(int(line_values[1]), int(line_values[2]))] = [float(line_values[3]), float(line_values[4])]
    file.close()
    return res


# Reads the provided file and returns a dictionary on format {(N, M) = [C, S]}
def read_file(egm2008=True):
    if egm2008:
        path = 'Datafiles/Part1/EGM2008.gfc'
    else:
        path = 'Datafiles/Part1/GGM03S.gfc'
    return create_dict(path)


# Boolean check to determine whether or not a number is even
def is_even(number):
    if number % 2 == 0:
        return True
    else:
        return False


def get_j_bar(number):
    ans = 0
    if is_even(number):
        ans = (-1)**(number/2) * ((3*mt.sqrt(e_sqd)**number)*(1-(number/2)+((5/2)*(j_2/e_sqd)*number))) / \
              ((number+1)*(number+3)*(2*number + 1)**0.5)
        return ans
    else:
        return ans


p_value_dict = {}


def get_initial_p_values(n, m, phi):
    t = mt.sin(phi)

    if n == 0 and m == 0:
        return 1

    if n == 1 and m == 0:
        return t*(3**0.5)

    if n == 1 and m == 1:
        return (3**0.5)*((1-t**2)**0.5)

    if n == 2 and m == 0:
        return (5**0.5)*(((3/2)*t**2)-0.5)

    if n == 2 and m == 1:
        return t*(15**0.5)*((1-t**2)**0.5)


def get_p_bar(n, m, phi):
    res = 0

    if (n, m, phi) in p_value_dict:
        return p_value_dict.get((n, m, phi))

    else:
        t = mt.sin(phi)

        if n <= 2 and m < 2:
            res = get_initial_p_values(n, m, phi)

        if n >= 2 and m == 0:
            res = -((((2*n+1)**0.5)/n)*((n-1)/((2*n - 3)**0.5)))*get_p_bar(n-2, m, phi) + \
                   t*((((2*n+1)*(2*n-1))/((n+m)*(n-m)))**0.5)*get_p_bar(n-1, m, phi)

        if n >= 3 and 1 <= m <= n-2:
            res = -((((2*n+1)*(n+m-1)*(n-m-1))/((2*n-3)*(n+m)*(n-m)))**0.5)*get_p_bar(n-2, m, phi) + \
                   t*((((2*n+1)*(2*n-1))/((n+m)*(n-m)))**0.5)*get_p_bar(n-1, m, phi)

        if n >= 1 and m == n-1:
            res = t*((2*n+1)**0.5)*get_p_bar(n-1, n-1, phi)

        if n >= 2 and n == m:
            res = ((2*n+1)**0.5)/((2*n)**0.5)*(1-t**2)**0.5*get_p_bar(n-1, n-1, phi)

        p_value_dict[(n, m, phi)] = res
        return res


def get_r_bar(n, m, dictionary):
    if m == 0:
        return dictionary.get((n, m))[0] - get_j_bar(n)
    else:
        return dictionary.get((n, m))[0]


def get_q_bar(n, m, dictionary):
    return dictionary.get((n, m))[1]


# phi = lat, lmd = long
# @np.vectorize
def get_n_grv(phi_in, lmd_in, dictionary, n_max):

    phi = mt.radians(phi_in)
    lmd = mt.radians(lmd_in)
    total_sum = 0
    for n in range(2, n_max + 1):
        part_sum = 0
        for m in range(0, n + 1):
            part_sum += get_p_bar(n, m, phi) * (
                        get_r_bar(n, m, dictionary) * mt.cos(m * lmd) + get_q_bar(n, m, dictionary) * mt.sin(m * lmd))
        total_sum += part_sum * (a / R) ** n

    return total_sum * (GM / (R * gamma))


# Lagre N verdien for en spesifikk n og m i en dict.
# I stede for å kalle get_p for hver loop, sjekk om du har N-verdi for (n, m-1) og pluss på en til