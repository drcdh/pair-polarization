from functools import partial
from sys import argv

import matplotlib.pyplot as plt
import vegas
from numpy import pi as Pi, sin, cos, sqrt, power
from pycuba import Vegas
from scipy import integrate

MASS_MEV = 0.511  # mass, electron/positron


# int Integral(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
def _sig(X, omega, phi, psi):
    en = X[0]*omega  # energy, electron
    t  = X[1]*Pi/2.  # theta, electron
    tp = X[2]*Pi/2.  # theta, positron
#def _sig(en, t, tp, omega, phi, psi):
    enp = omega - en  # energy, positron

    b  = sqrt(1. - power(MASS_MEV/en,  2)) if en > MASS_MEV else 0
    bp = sqrt(1. - power(MASS_MEV/enp, 2)) if en > MASS_MEV else 0

    qa = en*enp*(1. - b*bp*(sin(t)*sin(tp)*cos(phi)
                   + cos(t)*cos(tp)))
    qb = (omega*en *(b *cos(t)  - 1)
        + omega*enp*(bp*cos(tp) - 1) + MASS_MEV**2)
    q2 = -2.*(qa + qb)

    a  = b  * sin(t)  / (1. - b  * cos(t))
    ap = bp * sin(tp) / (1. - bp * cos(tp))

    a1  = a  * cos(phi + psi)
    a1p = ap * cos(psi)

    dd = sin(t) * sin(tp) * b * bp

    s1 = (
            (en * enp * dd / pow(q2, 2) / pow(omega, 3)) * (4.*pow((en*a1p + enp*a1), 2)
                                                          - q2 * pow((a1p - a1), 2))
          )

    s2 = (
            a * ap / omega / pow(q2, 2) * pow(en * b * sin(t), 2)
        + pow(enp * bp * sin(tp), 2)
        + 2 * en * enp * b * bp * cos(phi) * sin(t) * sin(tp)
    )

    f = 0.  # f and Z are unused

    sig = (2./pow((2.*Pi), 2))*(s2 - s1)*pow((1.-f), 2)
    sig = sig*omega*pow(Pi, 2)/4.  # normalización de la integración...

    if sig < 0:
        sig = 0.

    return sig


def vegas_sigma(omega, phi, psi, verbose: bool = False,
        adapt_nitn: int = 10, adapt_neval: int = 10000,
        nitn: int = 10, neval: int = 2000,
        alpha: float = .5,
        sdev: bool = False):
    vegas_integrand = partial(_sig, omega=omega, phi=phi, psi=psi)
    integ = vegas.Integrator([[0, 1], [0, Pi], [0, Pi]])
    # step 1 -- adapt to sigma; discard results
    adapt = integ(vegas_integrand, nitn=adapt_nitn, neval=adapt_neval, alpha=alpha,
                  rtol=1e-3, atol=1e-12,)
    if verbose:
        print('ADAPTATION:')
        print(adapt.summary())
    # step 2 -- integ has adapted to f; keep results
    result = integ(vegas_integrand, nitn=nitn, neval=neval, alpha=alpha)
    if verbose:
        print('ESTIMATION:')
        print(result.summary())
        print('    Result = %f +- %f    Q = %.2f' % (result.mean, result.sdev, result.Q))
    if sdev:
        return result.mean, result.sdev
    return result.mean


def PyCube_print_results(name, results):
    keys = ['nregions', 'neval', 'fail']
    keys = list(filter(results.has_key, keys))
    text = ["%s %d" % (k, results[k]) for k in keys]
    print("%s RESULT:\t" % name.upper() + "\t".join(text))
    for comp in results['results']:
      print("%s RESULT:\t" % name.upper() + \
    "%(integral).8f +- %(error).8f\tp = %(prob).3f\n" % comp)


def PyCubaIntegrand(ndim, xx, ncomp, ff, userdata):
    omega, phi, psi = userdata
    ff[0] = _sig(xx, *userdata)
    return 0


def PyCuba_sigma(omega, phi, psi):
    PyCube_print_results("Vegas", Vegas(PyCubaIntegrand, 3, verbose=2))


def main():
    omega_v = [60., 70., 80., 90., 100.,
                                   200., 300., 400., 600., 700., 800., 900., 1000.]

    phi_break = 2.9  # Change the phi sampling rate at this value
    n1, n = 10, 20  # n is "paso"
    phi_v = [phi_break * i / n1 for i in range(n1)]
    phi_v += [(Pi - phi_break) * i / (n - n1) + phi_break for i in range(n - n1)]
    phi_v = phi_v[::-1]

    psi_v = [j*Pi/n for j in range(n)]

    res = []
    for omega in [100.]:  #omega_v:
        _res = []
        for phi in phi_v:
            __res = []
            for psi in psi_v:
                fq = integrate.nquad(_sig,
                                     [[0., 1],  # en
                                      [0., Pi],  # t
                                      [0., Pi]],  # tp
                                     args=(omega, phi, psi),
                                     )
                print("%.1f  %.2f  %.2f  %.6e +- %.6e" % (omega, phi, psi, fq[0], fq[1]))
                __res.append(fq)
            _res.append(__res)
        res.append(_res)
    # res is now a 3D grid indexed like [omega][phi][psi]
    res = res[0]
    plt.plot(res, phi_v, psi_v)


if __name__ == '__main__':
    #vegas_sigma(100, 3.12, 0., verbose=True)
    PyCuba_sigma(100., 3.12, 0.)
    #main()
