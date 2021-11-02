import sys
import numpy as np
import csv

from collections import defaultdict

CSV_FILE = sys.argv[1]

PSI_0_BY_OMEGA = defaultdict(list)
PSI_1p57_BY_OMEGA = defaultdict(list)

with open(CSV_FILE, 'r') as f:
    header = True
    for row in csv.reader(f):
        if header:
            header = False
            continue
        _omega = float(row[0])
        _phi = float(row[1])
        _psi = float(row[2])
        _sigma = float(row[3])
        _rerr = float(row[4])
        _D = PSI_0_BY_OMEGA if _psi == 0. else PSI_1p57_BY_OMEGA
        _D[_omega] += [[_phi, _sigma, _rerr]]

import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111)

COLORS = ['b', 'g', 'r', 'm', 'k', 'y', 'c']
for _omega, _color in zip(PSI_0_BY_OMEGA.keys(), COLORS):
    psi0_data = np.array(PSI_0_BY_OMEGA[_omega])
    psi1p57_data = np.array(PSI_1p57_BY_OMEGA[_omega])
#    ax.plot(psi0_data[:, 0], psi0_data[:, 1], label=f"{_omega} MeV, Psi = 0", c=_color)
#    ax.plot(psi1p57_data[:, 0], psi1p57_data[:, 1], label=f"{_omega} MeV, ψ = π/2", c=_color, linestyle='--')
    ax.errorbar(psi0_data[:, 0], psi0_data[:, 1], yerr=psi0_data[:, 2], label=f"ω = {_omega:.0f} MeV, ψ = 0", c=_color)
    ax.errorbar(psi1p57_data[:, 0], psi1p57_data[:, 1], yerr=psi1p57_data[:, 2], label=f"ω = {_omega:.0f} MeV, ψ = π/2", c=_color, linestyle='--')

ax.set_xlabel('φ [rad]')
ax.set_ylabel('dσ/dφ/dψ/α(Zr)^2')
ax.set_xlim(np.pi*3./4., np.pi)
#; ax.set_zlim(4, 16)
ax.yaxis.tick_right()
#ax.set_yscale('log')
ax.yaxis.grid(True, which='minor', linestyle=':')
ax.legend(loc='best')
ax.grid()
plt.show()
