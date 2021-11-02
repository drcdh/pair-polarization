import csv
from functools import partial
from sys import argv

import matplotlib.pyplot as plt
#import vegas
from numpy import (pi as Pi, sin, cos, sqrt,
 #power as pow,
  linspace)
from cycuba import Vegas, Suave, Divonne, Cuhre
from scipy import integrate

MASS_MEV = 0.511  # mass, electron/positron


# int Integral(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
#def _sig(X, omega, phi, psi):
#    en = X[0]*omega  # energy, electron
#    t  = X[1]*Pi/2.  # theta, electron
#    tp = X[2]*Pi/2.  # theta, positron
def _sig(en, t, tp, omega, phi, psi):
    # Scale energy from 0,1 to 0,omega
    en *= omega
    # Scale thetas from 0,1 to 0,Pi
    t *= Pi/2
    tp *= Pi/2

    enp = omega - en  # energy, positron

    b  = 0. if en  < MASS_MEV else sqrt(1. - pow(MASS_MEV/en,  2))
    bp = 0. if enp < MASS_MEV else sqrt(1. - pow(MASS_MEV/enp, 2))

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
            a * ap / omega / pow(q2, 2) * (pow(en * b * sin(t), 2)
                                         + pow(enp * bp * sin(tp), 2)
                + 2 * en * enp * b * bp * cos(phi) * sin(t) * sin(tp))
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

METHODS = {
    'c': Cuhre,
    'd': Divonne,
    's': Suave,
    'v': Vegas,
}
PARAM = dict(
    ranges=None,  # None for unit hypercube
#    nstart=1000,
#    nincrease=500,
#    nbatch=1000,
#    gridno=0,
    verbosity=0,
#    last_samples_only=False,
#    do_not_smooth=False,
    retain_state_file=False,
#    file_grid_only=False,
#    level=0,
    epsrel=1e-3,
    epsabs=1e-12,
#    seed=1959,
    mineval=0,
    maxeval=1e6,
)
def cycuba_sigma(omega, phi, psi, method: str ='v', verbose=False):
    method = METHODS[method.lower()[0]]
    def cuba_integrand(en, t, tp):
        return [_sig(en, t, tp, omega=omega, phi=phi, psi=psi)]
    cycuba_result = method(cuba_integrand, **PARAM)
    if verbose:
        print('>> sigma = %.6g\n'
              '        +- %.6g\n'
              '   prob = %.2f' % (cycuba_result[0][0], cycuba_result[1][0], cycuba_result[2][0]))
    return cycuba_result

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


OMEGA_V = [
#    5., 10.,
    20.,  50.,
    #60., 70., 80., 90.,
    100.,
    200.,# 300., 400., 500.,
    #600., 700., 800., 900., 1000.
]

def main(csvfile, method, verbose=False):
    writer = csv.writer(csvfile)
    writer.writerow('omega Phi Psi sigma error'.split(' '))

    #phi_break = 2.9  # Change the phi sampling rate at this value
    #n1, n = 10, 20  # n is "paso"
    phi_break, n1, n = Pi*3/4., 0, 50
    phi_v = [phi_break * i / n1 for i in range(n1)]
    phi_v += [(Pi - phi_break) * i / (n - n1) + phi_break for i in range(n - n1)]
    phi_v += [Pi]

    #psi_v = [j*Pi/n for j in range(n)]
    #psi_v += [Pi]
    psi_v = [0., Pi/2.]

    for omega in OMEGA_V:
        for phi in phi_v:
            for psi in psi_v:
                cyc = cycuba_sigma(omega, phi=phi, psi=psi, method=method, verbose=False)
                if verbose:
                    print("%.1f  %.2f  %.2f  %.6e +- %.6e" % (omega, phi, psi, cyc[0][0], cyc[1][0]))
                writer.writerow(['%.6g'%v for v in [omega, phi, psi, cyc[0][0], cyc[1][0]]])

def demo_methods():
    for m in 'vsdc':
        print('\n@@ METHOD = %s' % METHODS[m].__name__)
        for psi_, phi_ in ((0., 3.), (Pi/2, 3.), (Pi/2, 3.12), (0., 3.12)):
            print(f'\n== psi = {psi_}\n== phi = {phi_}')
            cycuba_sigma(100., phi=phi_, psi=psi_, method=m, verbose=True)

if __name__ == '__main__':
    #vegas_sigma(100, 3.12, 0., verbose=True)
#    for psi_ in linspace(0.01, 3.12, 10, endpoint=True):
#        print('\n###  PSI = %4.2f  ###\n' % psi_)
#        cycuba_vegas_sigma(100, psi_, 0.)
#    for phi_ in linspace(0.01, 3.14, 10, endpoint=True):
#        print('\n###  PHI = %4.2f  ###\n' % phi_)
#        cycuba_vegas_sigma(100, 3.12, phi_)
    #PyCuba_sigma(100., 3.12, 0.)
    with open(argv[1], 'w', newline='') as csvfile:
        main(csvfile, method=argv[2], verbose=True)
#    demo_methods()
