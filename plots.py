from types import LambdaType
import matplotlib.pyplot as plt
from matplotlib import rc

import numpy as np
from scipy.optimize import curve_fit
from scipy import constants as const
from collections import defaultdict
import os

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

    recta = np.poly1d(np.polyfit(X, Y, 2))
    value_at_zero = recta(0)
    print(f'Value at 0 ({n_rec}):', value_at_zero)

    if print_graph:
        X_pred = np.arange(0, max(X) + 0.001, 0.0001)
        Y_pred = recta(X_pred)

        plt.title(f'Ground state energy with external potential \n with shape of Sierpinski carpet of iteration {n_rec}')
        # plt.xscale("log")
        
        plt.plot(X, Y, 'x', color="black")
        plt.plot(X_pred, Y_pred, '--', label=n_rec)
        plt.xlabel('1/N')
        plt.ylabel('Ground state energy $\left ( \\frac{\\hbar ^2}{2m} \\right )$')
        plt.show()

    error = value_at_zero - data[-1][2]

    return value_at_zero, error 

    # --------------------------------------------------------------------------------------    


def rec_vs_energy():
    energies = []
    rec_iter = []
    errors = []
    max_i = 6
    for i in range(1, max_i + 1):
        en, err = energy_to_inf(i, print_graph=False)
        # Lmin = 1 / (3**i)
        # en_inc = en / (1 / (Lmin**2))
        energies.append(en)
        errors.append(err)
        rec_iter.append(i)


    def func(x, a, b):
        return b * np.exp(x*a) 

    popt, pcov = curve_fit(func, rec_iter[0:5], energies[0:5], maxfev=10000)
    X_pred = np.arange(1, rec_iter[-1], 0.0001)
    Y_pred = func(X_pred, 1.93, 1.79)
    
    params = list(map(lambda x: x, popt))
    print('')
    # print(params)
    print(f'y = {params[1]} * exp({params[0]}*x)')
    
    plt.plot(X_pred, Y_pred, '--')
    plt.plot(rec_iter[0:5], energies[0:5], 'x')
    # plt.errorbar(rec_iter[0:5], energies[0:5], yerr=errors[0:5], fmt='.',capsize=2)
    plt.plot(rec_iter[5:], energies[5:], 'x')
    # plt.errorbar(rec_iter[5:], energies[5:], yerr=errors[5:], fmt='.',capsize=2)

    X_upper_bound = list(range(9))
    Y_upper_bound = [np.pi**2 * 9**i for i in range(9)]
    # plt.plot(X_upper_bound, Y_upper_bound, label="upper bound")

    plt.yscale("log")
    plt.title("Ground state energy vs iteration of the fractal")
    plt.xlabel("Iteration of the fractal")
    plt.ylabel("Ground state energy $\left ( \\frac{\\hbar ^2}{2m} \\right )$")
    plt.xticks(list(range(1, max_i + 1)))
    # plt.legend()
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
                state, energy, IPR = line[:-1].split()
                states.append(int(state)) 
                energies.append(float(energy)) 
                IPRs.append(float(IPR))

        return np.array(states), np.array(energies), np.array(IPRs) 

    def normalize(l):
        l_max = max(l)
        l_min = min(l)
        return list(map(lambda x: (x - l_min)/(l_max - l_min), l))


    X0, Y0, Z0 = read_data('data/IPR_data/IPR_data_rec0')
    X1, Y1, Z1 = read_data('data/IPR_data/IPR_data_rec1')
    X2, Y2, Z2 = read_data('data/IPR_data/IPR_data_rec2')
    X3, Y3, Z3 = read_data('data/IPR_data/IPR_data_rec3')
    X4, Y4, Z4 = read_data('data/IPR_data/IPR_data_rec4')
    # X4_pbc, Y4_pbc, Z4_pbc = read_data('data/IPR_data/IPR_data_rec4_pbc')
    # X5, Y5, Z5 = read_data('data/IPR_data/IPR_data_rec5_243')

    plt.plot(X0, Y0, 'x', label='No potential')
    plt.plot(X1, Y1, 'x', label='iteration 1')
    plt.plot(X2, Y2, 'x', label='iteration 2')
    plt.plot(X3, Y3, 'x', label='iteration 3')
    plt.plot(X4, Y4, 'x', label='iteration 4')

    # plt.plot(X4[:30], Y4[:30], 'x', label='zbc iteration 4')
    # plt.plot(X4_pbc, Y4_pbc, 'x', label='pbc iteration 4')
    # plt.plot(Y5, Z5, 'x', label='iteration 5')

    # plt.plot(X0, Z0, '--', label='iteration 0')

    # Y0_diff = [Y0[i] - Y0[i - 1] for i in range(1, len(Y0))]
    # Y1_diff = [Y1[i] - Y1[i - 1] for i in range(1, len(Y1))]
    # Y4_diff = [Y4[i] - Y4[i - 1] for i in range(1, len(Y4))]
    
    # plt.plot(list(range(len(Y0_diff))), Y0_diff, 'x', label='iteration 0')
    # plt.plot(list(range(len(Y0_diff))), Y4_diff, 'x', label='iteration 4')
    # plt.ylabel("Energy")
    # plt.show()
    
    # plt.hist(Y0_diff, bins=40)
    # plt.hist(Y4_diff, bins=40)
    

    plt.ylabel("Energy $\left ( \\frac{\\hbar ^2}{2m} \\right )$")
    plt.xlabel("Eigenstate")
    # plt.ylabel("IPR")
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
    # plt.legend()
    plt.title('Random walks')
    plt.xlabel('Time')
    plt.ylabel('Distance')
    plt.show()

    plt.plot(list(range(9)), [2] + diff_coefficients, 'x')
    plt.title('Diffusion coeficients')
    plt.xlabel('Iteration of fractal')
    plt.ylabel('Diffusion coefficient')
    plt.show()


def pbc_zbc():
    with open('3_zbc.txt') as f:
        X = []
        Y = []
        for line in f.readlines():
            N, _, energy, _ = line.split()
            X.append(1/int(N))
            Y.append(float(energy))
    
    recta = np.poly1d(np.polyfit(X, Y, 2))
    value_at_zero = recta(0)
    print(f'Value at 0 (zbc):', value_at_zero)

    X_pred = np.arange(0, max(X) + 0.001, 0.0001)
    Y_pred = recta(X_pred)
    
    plt.plot(X, Y, 'x')
    plt.plot(X_pred, Y_pred, '--', label="Zero bc")


    with open('3_pbc.txt') as f:
        X = []
        Y = []
        for line in f.readlines():
            N, _, energy, _ = line.split()
            X.append(1/int(N))
            Y.append(float(energy))
    
    recta = np.poly1d(np.polyfit(X, Y, 2))
    value_at_zero = recta(0)
    print(f'Value at 0 (pbc):', value_at_zero)

    X_pred = np.arange(0, max(X) + 0.001, 0.0001)
    Y_pred = recta(X_pred)
    
    plt.plot(X, Y, 'x')
    plt.plot(X_pred, Y_pred, '--', label="Periodic bc")
    plt.xlabel('1/N')
    plt.ylabel('Energy')
    plt.legend()
    plt.show()


def time_execution_eigs():
    T1 = []
    T2 = []
    N = []
    with open('data/execution_time/time_eigs.txt') as f:
        for line in f.readlines():
            n, t1, t2 = line.split()
            N.append(int(n)**4)
            T1.append(float(t1))
            T2.append(float(t2))

    plt.plot(N, T2, '--x', label="eigs method")
    plt.plot(N, T1, '--x', label="LOBPCG algorithm")
    plt.title('Comparing methods for obtaining first 20 eigenpairs')
    plt.xlabel('Size of the matrix')
    plt.ylabel('Execution time (s)')
    plt.legend()
    plt.show()


def time_evolution_displacement():
    def load_data(filename):
        Y = []
        with open(filename) as f:
            for line in f.readlines():
                x_2, y_2 = line.split()
                disp = float(x_2)
                # disp = (float(x_2) + float(y_2))/2
                Y.append(disp)

        return Y[:40] # There are 50 samples on the files, we can take the first ones
    
    # disp_0 = load_data('data/time_evolution_disp/disp_0_smallt.txt') 
    # disp_0 = load_data('data/time_evolution_disp/disp_0.txt') 
    disp_1 = load_data('data/time_evolution_disp/disp_1.txt') 
    disp_2 = load_data('data/time_evolution_disp/disp_2.txt') 
    disp_3 = load_data('data/time_evolution_disp/disp_3.txt') 
    disp_4 = load_data('data/time_evolution_disp/disp_4.txt') 

    X = np.array(list(range(len(disp_1)))) * 0.005
    # X = np.array(list(range(len(disp_0)))) * 0.0002

    # plt.plot(X, disp_0, label="no potential")
    plt.plot(X, disp_1, label="iteration 1")
    plt.plot(X, disp_2, label="iteration 2")
    plt.plot(X, disp_3, label="iteration 3")
    plt.plot(X, disp_4, label="iteration 4")

    # plt.plot([0, X[-1]], np.array([0.5/3, 0.5/3])**2, '--')
    # plt.plot([0, X[-1]], np.array([0.5, 0.5])**2, '--', color="tab:blue", label="box border")
    # plt.plot([0, X[-1]], np.array([1/9, 1/9])**2, '--', color="tab:orange")
    # plt.plot([0, X[-1]], np.array([1/(9**2), 1/(9**2)])**2, '--', color="tab:green")
    # plt.plot([0, X[-1]], np.array([0.3876, 0.3876])**2, '--')

    func = lambda x, a, c: np.power(x, a) + c 

    # for i, y_data in enumerate([disp_0, disp_1, disp_2, disp_3, disp_4]):
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
    ranges = [(3, 35), (3, 40), (3, 40), (3, 40)]
    for i, y_data in enumerate([disp_1, disp_2, disp_3, disp_4]):
        X_pred = X[ranges[i][0]:ranges[i][1]]
        Y_data = y_data[ranges[i][0]:ranges[i][1]]
        # X_pred = X[3:]
        # Y_data = y_data[3:]

        bounds = ([0, -0.5], [1, 0.5])
        # bounds = ([0, 0, -0.5], [1, 1, 0.5])
        popt, pcov = curve_fit(func, X_pred, Y_data, maxfev=10000, bounds=bounds)
        # if i == 0: popt = [0.50, 0.45, -0.006]
        Y_pred = func(X_pred, *popt)
        print(f'Iteration {i}: y = {round(popt[1], 4)} +  t^{round(popt[0], 4)}')
        # print(f'Iteration {i}: y = {round(popt[2], 4)} + {round(popt[1], 4)} * t^{round(popt[0], 4)}')
        plt.plot(X_pred, Y_pred, '--', color=colors[i])

    sigma = 0.1
    X_pred = np.arange(0, 0.05, 0.001)
    # ballistic_exp = (sigma**2/2) + ((X_pred**2) / (sigma**2))
    # ballistic_exp = (sigma**2/2) * (1 + (((const.h * X) / (sigma**2))**2) )

    # plt.plot(X_pred, ballistic_exp, '--', label="ballistic")
    # plt.plot(X[:11], ballistic_exp[:11], '--', label="ballistic")

    # plt.plot([0, X[-1]], np.array([sigma**2/2, sigma**2/2]), '--', label="$\sigma^2/2$")


    plt.legend(loc='upper left')
    plt.xlabel('Time')
    plt.ylabel('MSD')
    plt.title('Expected value of MSD for different external potential')
    # plt.yscale('log')
    plt.show()





if __name__ == "__main__":
    # execution_time()
    # IPR_states()
    # rec_vs_energy()
    # energy_to_inf(3, print_graph=True)
    # min_size_energy_vs_rec()
    # random_walks()
    # pbc_zbc()
    time_execution_eigs()
    # time_evolution_displacement()