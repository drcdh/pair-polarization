import sys
import numpy as np
import csv
OMEGA = float(sys.argv[2])

Phi, Psi, sigma = [], [], []
with open(sys.argv[1], 'r') as f:
    header = True
    for row in csv.reader(f):
        if header:
            header = False
            continue
        if float(row[0]) != OMEGA:
            continue
        Phi.append(float(row[1]))
        Psi.append(float(row[2]))
        sigma.append(float(row[3]))

import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(Phi, Psi, sigma)
ax.set_xlabel('Phi')
ax.set_ylabel('Psi')
#ax.set_xlim(3, np.pi)#; ax.set_zlim(4, 16)
plt.show()
