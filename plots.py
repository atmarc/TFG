import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from collections import defaultdict


# ---------------------------- 1/N vs energy ----------------------------
def rec_vs_energy():
    def energy_to_inf(n_rec, print_graph=False):
        data = []
        with open(f'data/{n_rec}_data.txt', 'r') as f:
            for line in f.readlines():
                if line[-1] == '\n': line = line[:-1] 
                n, rec, e1, error_flag = line.split()
                data.append((int(n), int(rec), float(e1), int(error_flag)))

        if len(data) == 1:
            # If we have only one sample of this recursive level, we return that value
            print(f'Value at 0 ({n_rec}):', data[0][2])
            return data[0][2]

        MAX_SAMPLE_NUM = 5

        X = [1/x for (x, _, _, _) in data[-MAX_SAMPLE_NUM:]]        
        Y = [y for (_, _, y, _) in data[-MAX_SAMPLE_NUM:]]        

        recta = np.poly1d(np.polyfit(X, Y, 2))
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

        return value_at_zero

    # --------------------------------------------------------------------------------------    

    energies = []
    rec_iter = []
    for i in range(1,9):
        e = energy_to_inf(i, print_graph=False)
        energies.append(e)
        rec_iter.append(i)


    def func(x, b, c):
        return b * np.exp(x*c) 

    popt, pcov = curve_fit(func, rec_iter[0:5], energies[0:5])
    X_pred = np.arange(0, 8.001, 0.0001)
    Y_pred = func(X_pred, *popt)
    
    params = list(map(lambda x: round(x,4), popt))
    print('')
    print(f'y = {params[0]} * exp({params[1]}*x)')
    
    plt.plot(X_pred, Y_pred, '--', label="Aproximation")
    plt.plot(rec_iter[0:5], energies[0:5], 'x')
    plt.plot(rec_iter[5:], energies[5:], 'x')

    
    
    X_upper_bound = list(range(9))
    Y_upper_bound = [np.pi**2 * 9**i for i in range(9)]
    plt.plot(X_upper_bound, Y_upper_bound, label="upper bound")

    # plt.xscale("log")
    plt.yscale("log")
    plt.title("Ground energy vs recursion iteration of the fractal")
    plt.legend()
    plt.show()


def execution_time():
    def ex_time_n(n_rec):
        data = []
        with open(f'data/ex_time_{n_rec}.txt', 'r') as f:
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
        data = []
        with open(filename, 'r') as f:
            for line in f.readlines():
                state, energy, IPR = line[:-1].split()
                data.append((int(state), float(energy), float(IPR)))
        return data

    data_rec2 = read_data('IPR_data_rec2')
    data_rec3 = read_data('IPR_data_rec3')
    data_rec4 = read_data('IPR_data_rec4')

    X_r2 = [x for (x, _, _) in data_rec2]
    Y_r2 = [y for (_, y, _) in data_rec2]
    Z_r2 = [z for (_, _, z) in data_rec2]

    X_r3 = [x for (x, _, _) in data_rec3]
    Y_r3 = [y for (_, y, _) in data_rec3]
    Z_r3 = [z for (_, _, z) in data_rec3]
   
    X_r4 = [x for (x, _, _) in data_rec4]
    Y_r4 = [y for (_, y, _) in data_rec4]
    Z_r4 = [z for (_, _, z) in data_rec4]
   
    # Y_norm = list(map(lambda x: (x - min(Y))/(max(Y) - min(Y)), Y))
    # Z_norm = list(map(lambda x: (x - min(Z))/(max(Z) - min(Z)), Z))
   
    plt.plot(X_r2, Z_r2, '--', label='rec 2')

    plt.plot(X_r3, Z_r3, '--', label='rec 3')
    # plt.plot(X_r3, Z_r3, '.')
    plt.plot(X_r4, Z_r4, '--', label='rec 4')
    # plt.plot(X_r4, Z_r4, '.')
    
    # plt.plot(X, Z, 'x')
    plt.legend()
    plt.show()



if __name__ == "__main__":
    # execution_time()
    # IPR_states()
    rec_vs_energy()