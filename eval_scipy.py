import numpy as np
from scipy import integrate

from depaola import get_integrand, _I


def get_sp_integrand(w, phi, psi):
    _I = get_integrand(w, phi, psi)

    def blah(x, th_p, th_m):
        return _I([th_m, th_p, x])

    return blah


#integrand = get_sp_integrand(100., np.pi, np.pi / 2)
#arglebargle = integrate.tplquad(integrand, 0, 1, 0, np.pi, 0, np.pi)
#print(arglebargle)

fq = integrate.nquad(lambda x, tp, tm, w, phi, psi: get_integrand(w, phi, psi)([x, tp, tm]),
                     [[0., 1],
                      [0., np.pi],
                      [0., np.pi]],
                     args=(100., np.pi/2, np.pi/2),
                     )
print(fq)
