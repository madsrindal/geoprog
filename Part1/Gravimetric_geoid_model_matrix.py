import numpy as np
import time as t
import math as mt
from Part1 import Gravimetric_geoid_model as Ggm
from Part1 import p_bar_testing as pbt

# Constants
a = 6378137.0000
b = 6356752.3141
f_inv = 298.257222101
e_sqd = 0.006694380023
R = 6371000.7900
GM = 3986005 * 10**8
gamma = 9.81
j_2 = 0.108263*(10**(-2))


# Need matrices for P_nm, R_nm and Q_nm
def get_prq_matrices(n_max, phi, egm2008=True):
    print('leser fil...')
    cs_dict = Ggm.read_file(egm2008)
    print('ferdig...')

    phi = mt.radians(phi)

    p_nm = np.zeros((n_max+1) ** 2).reshape(n_max+1, n_max+1)
    r_nm = np.zeros((n_max+1) ** 2).reshape(n_max+1, n_max+1)
    q_nm = np.zeros((n_max+1) ** 2).reshape(n_max+1, n_max+1)

    p_dict = {}

    for n in range(2, n_max+1):
        for m in range(0, n+1):
            p_nm[n][m] = pbt.get_p_bar(n, m, phi, p_dict)
            # p_nm[n][m] = Ggm.get_p_bar(n, m, phi)
            r_nm[n][m] = Ggm.get_r_bar(n, m, cs_dict)
            q_nm[n][m] = Ggm.get_q_bar(n, m, cs_dict)

    return p_nm, r_nm, q_nm


def get_cos_sin_matrices(n_max, lmd):
    cos_m = np.zeros(n_max + 1)
    sin_m = np.zeros(n_max + 1)

    lmd = mt.radians(lmd)

    for m in range(0, n_max + 1):
        cos_m[m] = mt.cos(m*lmd)
        sin_m[m] = mt.sin(m*lmd)

    return cos_m, sin_m


def get_first_sum_matrix(phi, lmd, n_max, egm2008=True):

    p_matrix, r_matrix, q_matrix = get_prq_matrices(n_max, phi, egm2008)
    cos_matrix, sin_matrix = get_cos_sin_matrices(n_max, lmd)

    a_matrix = np.multiply(r_matrix, cos_matrix) + np.multiply(q_matrix, sin_matrix)
    b_matrix = np.multiply(a_matrix, p_matrix)
    return b_matrix


def sum_of_matrix(matrix):
    dim = np.shape(matrix)[0]
    matrix_vector = np.zeros(dim-2)

    for n in range(2, dim):
        matrix_vector[n-2] = sum(matrix[n])

    return matrix_vector


def get_ar_vector(n_max):
    res = np.zeros(n_max - 1)
    for n in range(2, n_max + 1):
        res[n-2] = (a / R) ** n

    return res


start_time = t.time()
get_prq_matrices(2050, -90, True)
# sum1 = get_first_sum_matrix(45, 45, 2150, True)
print('Calculated first sum matrix in ' + str(t.time() - start_time) + ' seconds')

'''second_time = t.time()
vector_sum = sum_of_matrix(sum1)
print('Created the vector sum in ' + str(t.time() - second_time) + ' seconds')

third_time = t.time()
const_vector = get_ar_vector(2150)
print('Created the constant vector in ' + str(t.time() - third_time) + ' seconds')

fourth_time = t.time()
sum2 = np.multiply(vector_sum, const_vector)
print('Multiplied sum and constant vectors in ' + str(t.time() - fourth_time) + ' seconds')

fifth_time = t.time()
tot_sum = sum(sum2) * (GM/(R*gamma))
print(sum(sum2))
print('Multiplied two values in ' + str(t.time() - fifth_time) + ' seconds')

print('Matrix calculated N: ', tot_sum)
print('Matrix calculations used: ' + str(t.time()-start_time) + ' seconds.')

print('--------------------------------------')'''

'''second_time = t.time()
egm_dict = Ggm.read_file(True)
print('Originally calculated N: ', Ggm.get_n_grv(45, 45, egm_dict, 2150))
print('Originally calculations used: ' + str(t.time()-second_time) + ' seconds.')'''
