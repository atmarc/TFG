import matplotlib.pyplot as plt
import numpy as np

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
X_pred = np.arange(0, max(X), max(X)/100)
Y_pred = recta(X_pred)

params_names = ['iterations', 'n_walkers', 'rec_lvl', 'L', 'N', 'tam_min', 'dt', 'data_size']
print('PARAMETERS:')
print('--------------------')
for val, name in zip(params, params_names):
    print(name + ':', val)

print('--------------------')
print('Ajuste:', recta)


plt.plot(X, values, label="Desplazamiento")
# plt.plot(X, X * 2, label="Dispersión libre")
# plt.xscale("log")
# plt.yscale("log")
plt.plot(X_pred, Y_pred, label="Ajuste desplz")
plt.legend()
plt.show()