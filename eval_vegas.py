
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import vegas

from depaola import *

def eval_at(psi: float, phi: float, w: float = 100., verbose: bool = False,
            adapt_nitn: int = 10, adapt_neval: int = 10000,
            nitn: int = 10, neval: int = 2000,
            alpha: float = .5,
            sdev: bool = False,
            parts: str = 'ABC') -> float:
    if verbose:
        print('\nEvaluating parts {:s} at ϕ, ψ, ω = ({:.2f}, {:.2f}, {:.1f})'.format(parts, phi, psi, w))

    _sigma = get_integrand(w, phi, psi, parts)

    integ = vegas.Integrator([[0, 1], [0, pi], [0, pi]])
    # step 1 -- adapt to sigma; discard results
    adapt = integ(_sigma, nitn=adapt_nitn, neval=adapt_neval, alpha=alpha)
    if verbose:
        print('ADAPTATION:')
        print(adapt.summary())
    # step 2 -- integ has adapted to f; keep results
    result = integ(_sigma, nitn=nitn, neval=neval, alpha=alpha)
    if verbose:
        print('ESTIMATION:')
        print(result.summary())
        print('    Result = %f +- %f    Q = %.2f' % (result.mean, result.sdev, result.Q))
    if sdev:
        return result.mean, result.sdev
    return result.mean


def main():
    _X0 = dict(psi=0.4, phi=3.)
    _X1 = dict(psi=3.12, phi=3.12)
    parts = 'ABC'
    for _p in parts:
        break
        eval_at(**_X0, verbose=True, parts=_p)
        eval_at(**_X1, verbose=True, parts=_p)

    #return

    eval_at(**_X0, verbose=True)
    eval_at(**_X1, verbose=True)

    return

    N = 20
    tr = list(range(N))
    I0 = numpy.array([eval_at(**_X0, sdev=True, parts=parts) for _ in tr])
    I1 = numpy.array([eval_at(**_X1, sdev=True, parts=parts) for _ in tr])
    plt.errorbar(tr, I0[:, 0], yerr=I0[:, 1], label='ϕ, ψ = ({phi:.2f}, {psi:.2f})'.format(**_X0))
    plt.errorbar(tr, I1[:, 0], yerr=I1[:, 1], label='ϕ, ψ = ({phi:.2f}, {psi:.2f})'.format(**_X1))
    plt.grid(True)
    plt.legend()
    plt.show()

    return

    for i, _psi in enumerate(PSI):
        print('{:3d}/{:3d}: {:8f}'.format(i, len(PSI), _psi))
        for j, _phi in enumerate(PHI):
            PSI_PHI[i, j] = eval_at(psi=_psi, phi=_phi)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(*numpy.meshgrid(PSI, PHI), PSI_PHI, c=PSI_PHI)
    plt.show()


if __name__ == '__main__':
    main()
