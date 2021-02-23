import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from collections import defaultdict


# ---------------------------- 1/N vs energy ----------------------------
def n_vs_energy():
    def add_file(filename, data_label, pattern='-'):
        data = []
        with open(filename,'r') as f:
            for line in f.readlines():
                n, pot, rec, e1, e2 = line[:-1].split()
                if int(n) > 117: data.append((int(n), int(pot), int(rec), float(e1), float(e2)))
        
        X = np.array([1/x for (x, _, _, _, _) in data])
        Y = [y for (_, _, _, y, _) in data]
        Y2 = [y for (_, _, _, _, y) in data]

        recta = np.poly1d(np.polyfit(X, Y, 2))
        X_pred = np.arange(0, max(X) + 0.001, 0.001)
        Y_pred = recta(X_pred)

        plt.plot(X, Y, 'x', color="black")
        plt.plot(X_pred, Y_pred, pattern, label=data_label)
        print('Value at 0:', recta(0))
    
        recta2 = np.poly1d(np.polyfit(X, Y2, 2))
        X_pred2 = np.arange(0, max(X) + 0.001, 0.001)
        Y_pred2 = recta2(X_pred2)

        plt.plot(X, Y2, 'x', color="black")
        plt.plot(X_pred2, Y_pred2, pattern, label=data_label)
        print('Value at 0:', recta2(0))
    

    # add_file('data_pbc_pot_10.txt', "Ground energy with pbc and pot = 10", '--')
    # add_file('data_zbc_pot_10.txt', "Ground energy with zbc and pot = 10")
    
    # add_file('data_pbc_pot_100.txt', "Ground energy with pbc and pot = 100", '--')
    # add_file('data_zbc_pot_100.txt', "Ground energy with zbc and pot = 100")
    
    # add_file('data_pbc_pot_1000.txt', "Ground energy with pbc and pot = 1000", "--")
    # add_file('data_zbc_pot_1000.txt', "Ground energy with zbc and pot = 1000")
    
    add_file('data_pbc_pot_1000000.txt', "Ground energy with pbc and pot = 1000000")
    # add_file('data_zbc_pot_1000000.txt', "Ground energy with zbc and pot = 1000000")

    
    
    # ----- Show the correct version of 0 boundary condition --------
    # add_file('zbc_v1_no_pot.txt', "correct version")
    # add_file('zbc_v2_no_pot.txt', "version with 0s")

    # --------------------------------------------------------------------------------------    
    # plt.xscale("log")
    # plt.yscale("log")
    
    plt.legend()
    plt.xlabel("1/N")
    plt.ylabel("Energy")
    plt.title("1/N vs energy")

    plt.show()

# ---------------------------- rec vs energy ----------------------------
def rec_vs_energy():
    X = []
    E1 = []
    E2 = []
    data = []
    with open('data6.txt','r') as f:
        for line in f.readlines():
            rec, pot, n, e1, e2 = line[:-1].split()
            data.append((int(pot), int(rec), int(n), float(e1), float(e2)))
            X.append(int(rec))
            E1.append(float(e1))
            E2.append(float(e2))

    # Semilogaritmic
    # plt.xscale("log")
    plt.yscale("log")

    # recta2 = np.poly1d(np.polyfit(X, E1, 1))
    
    X = [1,2,3,4,5]
    for i in [1, 10, 100, 1000, 10000, 100000]:
        E = [e for (pot, _, _, e, _) in data if pot == i]
        plt.plot(X, E, label=i)
    # plt.plot(X, recta2(X), label="Regression")
    # plt.plot(X, E2,'.')

    # (pi^2)/(1/3^i)^2
    datos = [(np.pi**2)/(1/3**i)**2 for i in X]
    plt.plot(X, datos, label="1/3")
    
    # datos2 = [10/(1/3**i) for i in X]
    # plt.plot(X, datos2)

    plt.xlabel("Nivel de recursión")
    plt.ylabel("Energía")
    plt.legend()
    plt.title("Recursion level vs energy with sierpinski external potential with N = 3^5")
    plt.show()

# ---------------------------- Potential vs energy ----------------------------
def pot_vs_energy():
    X = []
    E1 = []
    E2 = []
    with open('data5.txt','r') as f:
        for line in f.readlines():
            pot, n, e1, e2 = line[:-1].split()
            if int(pot) > 10:
                X.append(1/int(pot))
                E1.append(float(e1))
                E2.append(float(e2))

    # ajuste = np.poly1d(np.polyfit(X, E1, 3))
    # X_pred = np.arange(0, 0.0125, 0.0001)
    # plt.plot(X_pred, ajuste(X_pred))

    def func(x, a, b, c, d):
        return a - b*np.arctan(x*c + d) 

    popt, pcov = curve_fit(func, X, E1)

    X_pred = np.array(X + [0])
    plt.plot(X_pred, func(X_pred, *popt), '--', label="Aproximation")
    print('Value at 0:', func(0, *popt))

    plt.plot(X, E1, '.', label="Energy 0")
    # plt.plot(X, E2, 'o', label="Energy 1")
    
    # plt.xscale("log")
   
    plt.title("Sierpinski external potential vs energy with N = 3^5 and rec = 2")
    plt.xlabel("1/Potencial")
    plt.ylabel("Energía")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    n_vs_energy()