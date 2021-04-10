import matplotlib.pyplot as plt
import numpy as np

with open('track_mat.txt') as f:
    data = []
    for i, line in enumerate(f.readlines()):
        line_values = line[:-2].split(',')
        values = np.array(list(map(float, line_values)))
        data.append(values)
        if i % 2000 == 0: print(f'line {i} computed')

plt.imshow(data)
plt.show()

