import math as mt
import numpy as np

# Constants provided in the assignment text
a = 6378137.0000
b = 6356752.3141
f_inv = 298.257222101
e_sqd = 0.006694380023
R = 6371000.7900
GM = 3986005 * 10**8
gamma = 9.81
j_2 = 0.108263*(10**(-2))


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
        path = '../Datafiles/Part1/EGM2008.gfc'
    else:
        path = '../Datafiles/Part1/GGM03S.gfc'
    return create_dict(path)


# Boolean check to determine whether or not a number is even
def is_even(number):
    if number % 2 == 0:
        return True
    else:
        return False


# Function used to return a j_value of degree n
def get_j_bar(number):
    ans = 0
    if is_even(number):
        ans = (-1)**(number/2) * ((3*mt.sqrt(e_sqd)**number)*(1-(number/2)+((5/2)*(j_2/e_sqd)*number))) / \
              ((number+1)*(number+3)*(2*number + 1)**0.5)
        return ans
    else:
        return ans


# Function used to return initial values of p, given a specific pairing of n, m and phi.
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


# The implemented function used to return the Normalized Legendre's Polynomials for a specific pair of n, m and phi.
# The function is recursive, making it essential to implement with time and memory complexity in mind.
def get_p_bar(n, m, phi, p_value_dict):
    res = 0

    # If the value is calculated before (and hence saved to the dictionary), the function returns the value directly.
    if (n, m) in p_value_dict:
        return p_value_dict.get((n, m))

    # If not, the function computes it based on a number of requirements, and then saves the value to the dictionary.
    else:
        t = mt.sin(phi)

        if n <= 2 and m < 2:
            res = get_initial_p_values(n, m, phi)

        if n >= 2 and m == 0:
            res = -((((2*n+1)**0.5)/n)*((n-1)/((2*n - 3)**0.5)))*get_p_bar(n-2, m, phi, p_value_dict) + \
                   t*((((2*n+1)*(2*n-1))/((n+m)*(n-m)))**0.5)*get_p_bar(n-1, m, phi, p_value_dict)

        if n >= 3 and 1 <= m <= n-2:
            res = -((((2*n+1)*(n+m-1)*(n-m-1))/((2*n-3)*(n+m)*(n-m)))**0.5)*get_p_bar(n-2, m, phi, p_value_dict) + \
                   t*((((2*n+1)*(2*n-1))/((n+m)*(n-m)))**0.5)*get_p_bar(n-1, m, phi, p_value_dict)

        if n >= 1 and m == n-1:
            res = t*((2*n+1)**0.5)*get_p_bar(n-1, n-1, phi, p_value_dict)

        if n >= 2 and n == m:
            res = ((2*n+1)**0.5)/((2*n)**0.5)*(1-t**2)**0.5*get_p_bar(n-1, n-1, phi, p_value_dict)

        p_value_dict[(n, m)] = res
        return res


# Function used to return a value of r given a specific pair of n and m.
def get_r_bar(n, m, dictionary):
    if m == 0:
        return dictionary.get((n, m))[0] - get_j_bar(n)
    else:
        return dictionary.get((n, m))[0]


# Function used to return a value of q given a specific number of n and m
def get_q_bar(n, m, dictionary):
    return dictionary.get((n, m))[1]


# A function used to create two matrices of every possible value of r and q up to the given n_max as input.
def get_rq_bar_matrix(dictionary, n_max):
    # Creates two matrices containing only zeroes with of desired size and shape as provided by n_max
    r_matrix = np.zeros((n_max+1)**2).reshape(n_max+1, n_max+1)
    q_matrix = np.zeros((n_max + 1) ** 2).reshape(n_max + 1, n_max + 1)

    # Fills the matrices with correct values
    for n in range(2, n_max+1):
        for m in range(0, n+1):
            r_matrix[n][m] = get_r_bar(n, m, dictionary)
            q_matrix[n][m] = get_q_bar(n, m, dictionary)

    return r_matrix, q_matrix


# Function used to compute the gravimetric height for a given pair of phi and lambda.
def get_n_grv(lmd_in, n_max, p_matrix, r_matrix, q_matrix):

    lmd = mt.radians(lmd_in)
    total_sum = 0
    # The implemented formula summing over every combination of n and m
    for n in range(2, n_max + 1):
        part_sum = 0
        for m in range(0, n + 1):
            part_sum += p_matrix[n][m] * (
                        r_matrix[n][m] + q_matrix[n][m])
        total_sum += part_sum * (a / R) ** n

    return total_sum * (GM / (R * gamma))


# A modified function for computations of gravimetric heights.
# This function does not utilize matrices, and is somewhat slower in operation,
# but may initially look more alike the provided formula
def get_n_grv_modified(phi_in, lmd_in, dictionary, n_max, p_dict):

    phi = mt.radians(phi_in)
    lmd = mt.radians(lmd_in)
    total_sum = 0
    for n in range(2, n_max + 1):
        part_sum = 0
        for m in range(0, n + 1):
            part_sum += get_p_bar(n, m, phi, p_dict) * (
                        get_r_bar(n, m, dictionary) * mt.cos(m * lmd) + get_q_bar(n, m, dictionary) * mt.sin(m * lmd))
        total_sum += part_sum * (a / R) ** n

    return total_sum * (GM / (R * gamma))
