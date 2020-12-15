from time import time
import math as mt
from progressbar import progressbar
import numpy as np
from Part1 import Gravimetric_geoid_model as Ggm
from Part1 import p_bar_testing as pbt
import multiprocessing as mp


def calculate_ggm_model_world():
    start_time = time()
    print("Program starting...")

    print('Creating dictionaries from the GGM03S file')
    ggm03s_dict = Ggm.read_file(False)
    file = open('../../geoprog/Part1/Results/ggm_results.txt', 'w')
    file.write('LAT\tLON\tGeoidal Height\n')

    print('Calculating geoidal heights all over the world using GGM03S values and writing these to file...')
    for phi in progressbar(np.arange(-90, 90.5, 0.5)):
        for lmd in np.arange(-180, 180.5, 0.5):
            file.write(str(phi) + '\t' + str(lmd) + '\t' + str(Ggm.get_n_grv(phi, lmd, ggm03s_dict, n_max=180)) + '\n')

    file.close()
    print('---------------------------------------------------------------------------------------------------------')
    print('Calculated geoidal heights using GGM03S values in ' + str(time() - start_time) + ' seconds')


def calculate_egm_model_world():
    start_time = time()
    n_max = 2050

    print('Creating dictionaries from the EGM2008 file')
    egm2008_dict = Ggm.read_file()

    print('Creating r_ and q_matrices, multiplying each matrix with cos and sin of radians(lambda)*m')
    r_matrix, q_matrix = pbt.get_rq_bar_matrix(egm2008_dict, n_max)

    file = open('Results/egm_results_nmax.txt', 'w')
    file.write('LAT\tLON\tGeoidal Height\n')

    #  pool = mp.Pool(mp.cpu_count())

    r_matrix_list = []
    q_matrix_list = []

    print('Creating two lists containing of different matrices based on different lambdas...')
    for lmd in np.arange(-180, 180.5, 0.5):
        r = np.copy(r_matrix)
        q = np.copy(q_matrix)
        for m in range(0, n_max+1):
            r[:, m] *= mt.cos(m * mt.radians(lmd))
            q[:, m] *= mt.sin(m * mt.radians(lmd))
        r_matrix_list.append(r)
        q_matrix_list.append(q)

    print('Calculating geoidal heights all over the world using EGM2008 values and writing these to file...')
    for phi in progressbar(np.arange(-90, 90.5, 0.5)):

        p_matrix = np.zeros((n_max+1) ** 2).reshape(n_max+1, n_max+1)
        p_value_dict = {}
        for n in range(0, n_max+1):
            for m in range(0, n_max+1):
                p_matrix[n][m] = pbt.get_p_bar(n, m, mt.radians(phi), p_value_dict)

        lmd_counter = 0
        for lmd in progressbar(np.arange(-180, 180.5, 0.5)):
            # file.write(str(phi)+'\t'+str(lmd)+'\t'+str(pbt.get_n_grv_new(lmd, egm2008_dict, n_max, p_matrix))+'\n')
            file.write(str(phi)+'\t'+str(lmd)+'\t'+str(pbt.get_n_grv_new1(lmd, n_max, p_matrix,
                                                                          r_matrix_list[lmd_counter],
                                                                          q_matrix_list[lmd_counter])) + '\n')

            lmd_counter += 1
        '''lambdas = np.arange(-180, 180.5, 0.5)
        results = pool.starmap(pbt.get_n_grv_new, [(lmd, egm2008_dict, n_max, p_matrix) for lmd in lambdas])

        for i in range(len(results)):
            file.write(str(phi)+'\t'+str(lambdas[i])+'\t'+str(results[i])+'\n')'''

    file.close()
    print('Calculated geoidal heights using EGM2008 values in ' + str(time()-start_time) + ' seconds')


'''def calculate_egm_model_world1():
    start_time = time()
    n_max = 50

    print('Creating dictionaries from the EGM2008 file')
    egm2008_dict = Ggm.read_file()

    print('Creating r_ and q_matrices, multiplying each matrix with cos and sin of radians(lambda)*m')
    r_matrix, q_matrix = pbt.get_rq_bar_matrix(egm2008_dict, n_max)

    file = open('Results/egm_results_nmax.txt', 'w')
    file.write('LAT\tLON\tGeoidal Height\n')

    print('Creating two lists containing of different matrices based on different lambdas...')

    print('Calculating geoidal heights all over the world using EGM2008 values and writing these to file...')
    for phi in progressbar(np.arange(-90, 90.5, 0.5)):

        p_matrix = np.zeros((n_max+1) ** 2).reshape(n_max+1, n_max+1)
        p_value_dict = {}
        for n in range(0, n_max+1):
            for m in range(0, n_max+1):
                p_matrix[n][m] = pbt.get_p_bar(n, m, mt.radians(phi), p_value_dict)

        for lmd in np.arange(-180, 180.5, 0.5):
            r = np.copy(r_matrix)
            q = np.copy(q_matrix)

            for m in range(0, n_max + 1):
                r[:, m] *= mt.cos(m * mt.radians(lmd))
                q[:, m] *= mt.sin(m * mt.radians(lmd))

            # file.write(str(phi)+'\t'+str(lmd)+'\t'+str(pbt.get_n_grv_new(lmd, egm2008_dict, n_max, p_matrix))+'\n')
            file.write(str(phi)+'\t'+str(lmd)+'\t'+str(pbt.get_n_grv_new1(lmd, n_max, p_matrix, r, q)) + '\n')

        """lambdas = np.arange(-180, 180.5, 0.5)
        results = pool.starmap(pbt.get_n_grv_new, [(lmd, egm2008_dict, n_max, p_matrix) for lmd in lambdas])

        for i in range(len(results)):
            file.write(str(phi)+'\t'+str(lambdas[i])+'\t'+str(results[i])+'\n')"""

    file.close()
    print('Calculated geoidal heights using EGM2008 values in ' + str(time()-start_time) + ' seconds')'''


def test():
    tid = time()
    # egm2008_dict = Ggm.read_file()
    ggm03s_dict = Ggm.read_file(False)
    # print(egm2008_dict)
    print(Ggm.get_n_grv(61.9308563192723, 5.12764703841812, ggm03s_dict, 180))
    # print(Ggm.get_p_bar(2, 0, mt.asin(0.5)))
    # print(Ggm.get_r_bar(100, 89, ggm03s_dict))
    print('Calculated one N in ' + str(time()-tid) + ' seconds.')

    # lengde: 10.592871601647
    # bredde: 59.4218454378833
    # full n


'''def tull(n_max):
    egm_dict = Ggm.read_file(True)

    r_matrix, q_matrix = pbt.get_rq_bar_matrix(egm_dict, n_max)

    r_matrix_list = []
    q_matrix_list = []

    for lmd in progressbar(np.arange(-180, 180, 0.5)):
        r = np.copy(r_matrix)
        q = np.copy(q_matrix)
        for m in range(0, n_max + 1):
            r[:, m] *= mt.cos(m * mt.radians(lmd))
            q[:, m] *= mt.sin(m * mt.radians(lmd))
        r_matrix_list.append(r)
        q_matrix_list.append(q)

    return r_matrix_list, q_matrix_list'''


def test_script(n_max, phi, lmd):
    s1 = time()
    egm_dict = Ggm.read_file(True)
    print('Sekunder brukt på å lage egm_dict: ', time() - s1)

    s2 = time()
    # p_mtrix, p_dict = tull(phi)
    print('Sekunder brukt på å lage p_matrix: ', time() - s2)

    s3 = time()
    r_mtrix, q_mtrix = pbt.get_rq_bar_matrix(egm_dict, n_max)
    r_mtrix1, q_mtrix1 = pbt.get_rq_bar_matrix(egm_dict, n_max)
    print('Sekunder brukt på å lage 2x r_ og q_matrix: ', time() - s3)

    s8 = time()
    for m in range(0, n_max+1):
        r_mtrix1[:, m] *= mt.cos(m * mt.radians(lmd))
        q_mtrix1[:, m] *= mt.sin(m * mt.radians(lmd))
    print('Sekunder brukt på å gange matrisene med sin og cos: ', time() - s8)

    start = time()
    # print(pbt.get_n_grv_new(lmd, n_max, p_mtrix, r_mtrix, q_mtrix))
    print('Sekunder brukt på å regne ut n: ', time() - start)

    start = time()
    # print(pbt.get_n_grv_new1(lmd, n_max, p_mtrix, r_mtrix1, q_mtrix1))
    print('Sekunder brukt på å regne ut n_modified: ', time() - start)


if __name__ == '__main__':
    # start = time()
    calculate_egm_model_world()
    # print('Sekunder brukt på n=50, og full liste med matriser: ', time()-start)

    # start = time()
    # test_script(2050, -90, -90)
    # tull(2050)
    # print('Sekunder brukt på å regne ut matrise-listene: ', time()-start)

    # England: (phi: 20[45. 65]), (lmd: 30[-20, 10])