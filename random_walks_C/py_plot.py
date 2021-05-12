import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


with open('track_mat.txt') as f:
    data = []
    for i, line in enumerate(f.readlines()):
        line_values = line[:-1].split(',')
        values = np.array(list(map(float, line_values)))
        data.append(values)

plt.imshow(data)
plt.show()


with open('distances.txt') as f:
    data = []
    params = f.readline()
    params = list(map(float, (params[:-1].split())))

    iterations, n_walkers, rec_lvl, L, N, tam_min, dt, data_size = params

    values = f.read().split(',')
    values = list(map(float, values))

X = np.array(list(range(len(values)))) * dt * iterations/data_size

recta = np.poly1d(np.polyfit(X, values, 1))
X_pred2 = np.arange(0, max(X), max(X)/100)
Y_pred2 = recta(X_pred2)

def func(x, a, b, c):
    return b + c*np.power(x, a) 

popt, pcov = curve_fit(func, X, values)
X_pred = np.arange(0, max(X), max(X)/100)
Y_pred = func(X_pred, *popt)

params_names = ['iterations', 'n_walkers', 'rec_lvl', 'L', 'N', 'tam_min', 'dt', 'data_size']
print('PARAMETERS:')
print('--------------------')
for val, name in zip(params, params_names):
    print(name + ':', val)

# print('--------------------')
# print(f'Ajuste: dt^({popt[0]})')


# plt.plot(X_pred, Y_pred, '--', label="Ajuste diff")
plt.plot(X_pred2, Y_pred2, '--', label="Ajuste recta")
plt.plot(X, values, label="Desplazamiento")
# plt.plot(X, X * 2, label="Dispersión libre")
# plt.xscale("log")
# plt.yscale("log")
plt.legend()
plt.show()