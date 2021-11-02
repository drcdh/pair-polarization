
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

npz = np.load('lambda_omega_phi.npz')

lambdas = npz['lambdas']
omegas = npz['omegas']
phis = npz['phis']

levels = np.arange(-.5, .6, .1)

max_lambda = abs(lambdas).max()
norm = cm.colors.Normalize(vmax=max_lambda, vmin=-max_lambda)
cmap = cm.PRGn

fig, _axs = plt.subplots(nrows=3, ncols=1)
#_axs[-1,-1].remove()
axs = _axs.flatten()
for i, ax in enumerate(axs):
    if i == 3:
        break
    CS = ax.contourf(phis, omegas, lambdas[i].T, levels, norm=norm, cmap=cm.get_cmap(cmap, len(levels)-1))
#ax.clabel(CS, inline=True, fontsize=14)
fig.colorbar(CS)
ax.set_yscale('log')
ax.set_ylabel('ω/MeV')
ax.set_xlabel('ψ')
plt.show()
