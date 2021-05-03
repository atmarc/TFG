import matplotlib.pyplot as plt
from matplotlib import rc

import numpy as np
from scipy.optimize import curve_fit
from collections import defaultdict
import os
import matplotlib

DIMENSION = np.log(8) / np.log(3)

def energy_to_inf(n_rec, print_graph=False):
    data = []
    with open(f'data/good_data_per_iteration/{n_rec}_data.txt', 'r') as f:
        for line in f.readlines():
            if line[-1] == '\n': line = line[:-1] 
            n, rec, e1, error_flag = line.split()
            data.append((int(n), int(rec), float(e1), int(error_flag)))

    if len(data) == 1:
        # If we have only one sample of this recursive level, we return that value
        print(f'Value at 0 ({n_rec}):', data[0][2])
        return data[0][2], 0

    MAX_SAMPLE_NUM = 5

    X = [1/x for (x, _, _, _) in data[-MAX_SAMPLE_NUM:]]        
    Y = [y for (_, _, y, _) in data[-MAX_SAMPLE_NUM:]]        

    recta = np.poly1d(np.polyfit(X, Y, 1))
    value_at_zero = recta(0)
    print(f'Value at 0 ({n_rec}):', value_at_zero)

    if print_graph:
        X_pred = np.arange(0, max(X) + 0.001, 0.0001)
        Y_pred = recta(X_pred)

        plt.title(f'Getting energy when N --> inf with rec = {n_rec}')
        # plt.xscale("log")
        
        plt.plot(X, Y, 'x', color="black")
        plt.plot(X_pred, Y_pred, '--', label=n_rec)
        plt.show()

    error = value_at_zero - data[-1][2]

    return value_at_zero, error 

    # --------------------------------------------------------------------------------------    


def rec_vs_energy():
    energies = []
    rec_iter = []
    errors = []
    for i in range(1,9):
        en, err = energy_to_inf(i, print_graph=True)
        energies.append(en)
        errors.append(err)
        rec_iter.append(i)


    def func(x, a, b):
        return b * np.sinh(x*a) 

    popt, pcov = curve_fit(func, rec_iter[0:5], energies[0:5], maxfev=10000)
    X_pred = np.arange(1, rec_iter[-1], 0.0001)
    Y_pred = func(X_pred, *popt)
    
    params = list(map(lambda x: round(x,4), popt))
    print('')
    print(f'y = {params[0]} * exp({params[1]}*x)')
    
    plt.plot(X_pred, Y_pred, '--', label="Aproximation")
    # plt.plot(rec_iter[0:5], energies[0:5], 'x')
    plt.errorbar(rec_iter[0:5], energies[0:5], yerr=errors[0:5], fmt='.',capsize=2)
    # plt.plot(rec_iter[5:], energies[5:], 'x')
    plt.errorbar(rec_iter[5:], energies[5:], yerr=errors[5:], fmt='.',capsize=2)

    X_upper_bound = list(range(9))
    Y_upper_bound = [np.pi**2 * 9**i for i in range(9)]
    plt.plot(X_upper_bound, Y_upper_bound, label="upper bound")

    # plt.xscale("log")
    plt.yscale("log")
    plt.title("Ground energy vs recursion iteration of the fractal")
    plt.xlabel("Iteration of the fractal")
    plt.ylabel("Energy")
    plt.legend()
    plt.show()


def execution_time():
    def ex_time_n(n_rec):
        data = []
        with open(f'data/execution_time/ex_time_{n_rec}.txt', 'r') as f:
            for line in f.readlines():
                if line[-1] == '\n': line = line[:-1] 
                n, rec, e1, ex_time = line.split()
                data.append((int(n), int(rec), float(e1), float(ex_time)))

        X = [n**4 for (n, _, _, _) in data]
        Y = [t/60 for (_, _, _, t) in data]

        recta = np.poly1d(np.polyfit(X, Y, 2))
        X_pred = np.arange(0, max(X), max(X)/100)
        Y_pred = recta(X_pred)

        plt.title(f'Matrix size vs execution time')

        plt.plot(X, Y, 'x', color="black")
        plt.plot(X_pred, Y_pred, '--', label=f"ajuste rec = {n_rec}")
    
    
    ex_time_n(4)
    ex_time_n(5)
    plt.xlabel("Number of elements of matrix (N^4)")
    plt.ylabel("Execution time (min)")
    plt.legend()
    plt.show()


def IPR_states():
    def read_data(filename):
        states = []
        energies = []
        IPRs = []
        with open(filename, 'r') as f:
            for line in f.readlines():
                state, energy, IPR = line[:-1].split(',')
                states.append(int(state)) 
                energies.append(float(energy)) 
                IPRs.append(float(IPR))

        return np.array(states), np.array(energies), np.array(IPRs) 

    def normalize(l):
        l_max = max(l)
        l_min = min(l)
        return list(map(lambda x: (x - l_min)/(l_max - l_min), l))


    X0, Y0, Z0 = read_data('data/IPR_data/IPR_data_rec0_243')
    X1, Y1, Z1 = read_data('data/IPR_data/IPR_data_rec1_243')
    X2, Y2, Z2 = read_data('data/IPR_data/IPR_data_rec2_243')
    X3, Y3, Z3 = read_data('data/IPR_data/IPR_data_rec3_243')
    X4, Y4, Z4 = read_data('data/IPR_data/IPR_data_rec4_243')
    X5, Y5, Z5 = read_data('data/IPR_data/IPR_data_rec5_243')

    # plt.plot(Y0, Z0, 'x', label='iteration 0')
    # plt.plot(Y1, Z1, 'x', label='iteration 1')
    # plt.plot(Y2, Z2, 'x', label='iteration 2')
    # plt.plot(Y3, Z3, 'x', label='iteration 3')
    # plt.plot(Y4, Z4, 'x', label='iteration 4')
    # plt.plot(Y5, Z5, 'x', label='iteration 5')

    # plt.plot(X0, Z0, '--', label='iteration 0')
    # plt.plot(X1, Z1, '--', label='iteration 1')
    # plt.plot(X2, Z2, '--', label='iteration 2')
    # plt.plot(X3, Z3, '--', label='iteration 3')
    # plt.plot(X4, Z4, '--', label='iteration 4')
    # plt.plot(X5, Z5, '--', label='iteration 5')

    # plt.plot(X0, Y0 / (DIMENSION**0), 'x', label='iteration 0')
    # plt.plot(X1, Y1 / (DIMENSION**1), 'x', label='iteration 1')
    # plt.plot(X2, Y2 / (DIMENSION**2), 'x', label='iteration 2')
    # plt.plot(X3, Y3 / (DIMENSION**3), 'x', label='iteration 3')
    # plt.plot(X4, Y4 / (DIMENSION**4), 'x', label='iteration 4')
    # plt.plot(X5, Y5 / (DIMENSION**5), 'x', label='iteration 5')

    plt.xlabel("Energy")
    plt.ylabel("IPR")
    # plt.yscale("log")
    plt.title("Energy of the different eigenstates")
    plt.legend()
    plt.show()


def min_size_energy_vs_rec():
    data = []
    folder = 'data/good_data_per_iteration/'
    for i in range(1, 10):
        with open(f'{folder}{i}_data.txt') as f:
            for line in f.readlines():
                if line[-1] == '\n': line = line[:-1]                
                n, rec, energy, error_flag = line.split()
                data.append((int(n), int(rec), float(energy), int(error_flag)))

    for fact_min in range(1, 7):
        Y = [e for (n, rec, e, _) in data if n == (3**rec * fact_min)]
        X = list(range(1, len(Y) + 1))
        print(X)
        print(Y)
        def func(x, b, c):
            return b * np.exp(x*c) 

        popt, pcov = curve_fit(func, X, Y)
        X_pred = np.arange(0, max(X) + 0.001, 0.0001)
        Y_pred = func(X_pred, *popt)
        
        params = list(map(lambda x: round(x,4), popt))
        print('')
        print(f'y = {params[0]} * exp({params[1]}*x)')
        
        plt.plot(X_pred, Y_pred, '--', label="Aproximation")
        plt.plot(X, Y, 'x')
        plt.yscale('log')
        plt.title(f'Energy vs iteration (with min factor = {fact_min})')
        plt.xlabel('Iteration of the fractal')
        plt.xlabel('Energy')
        plt.show()


def random_walks():
    data = []
    files = os.listdir('data/random_walks')

    diff_coefficients = []

    for f in [x for x in files if x.startswith('distances')]:
        with open('data/random_walks/' + f) as f:
            params = f.readline()
            params = list(map(float, (params[:-1].split())))
            iterations, n_walkers, rec_lvl, L, N, tam_min, dt, data_size = params

            values = f.read().split(',')
            values = list(map(float, values))

        
        X = np.array(list(range(len(values)))) * dt * iterations/data_size
        recta = np.poly1d(np.polyfit(X, values, 1))
        X_pred = np.arange(0, max(X), max(X)/100)
        Y_pred = recta(X_pred)
        print(f'Ajuste para rec_lvl {int(rec_lvl)}:')
        print(recta)
        diff_coefficients.append(recta[1])
        # plt.plot(X_pred, Y_pred, '--',label=f"Ajuste lineal rec_lvl={rec_lvl}")
        plt.plot(X, values, label=f"Real values rec_lvl={rec_lvl}")
        # plt.xscale("log")
        # plt.yscale("log")
    plt.legend()
    plt.title('Random walks')
    plt.xlabel('Time')
    plt.ylabel('Distance')
    plt.show()

    plt.plot(list(range(9)), [2] + diff_coefficients, 'x')
    plt.title('Diffusion coeficients')
    plt.xlabel('Iteration of fractal')
    plt.ylabel('Diffusion coefficient')
    plt.show()

if __name__ == "__main__":
    # execution_time()
    IPR_states()
    # rec_vs_energy()
    # energy_to_inf(6, print_graph=True)
    # min_size_energy_vs_rec()
    # random_walks()