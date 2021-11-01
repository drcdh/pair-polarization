# Depoala Figure 2

import math
import sys

from typing import Iterable, List, Tuple

import numpy
from numpy import pi, sin, cos

SHAPE = 10, 50
PSI_PHI = numpy.zeros(SHAPE)
PSI = numpy.linspace(0, pi, SHAPE[0])
PHI = numpy.linspace(3, pi, SHAPE[1])

m_e = 0.511  # electron rest mass [MeV]
m2 = m_e ** 2
# a = 1 / 137.


def _q2(x, w, ct_p, ct_m, st_p, st_m, cph):
    q2 = -2 * (x * (1 - x) * (1 - st_p * st_m * cph - ct_p * ct_m)
               + x * (ct_p - 1)
               + (1 - x) * (ct_m - 1)
               + (m_e / w) ** 2)
    return q2


def _stuff(s, x, w, ct_p, ct_m, st_p, st_m, cph):
    q2 = _q2(x, w, ct_p, ct_m, st_p, st_m, cph)
    s *= (x * (1 - x)) * st_p * st_m / (q2 ** 2)
    s *= -2 * (1 / 2 / pi) ** 2
    s *= (m_e / w) ** 2
    return s


def _I_A(x, ct_p, ct_m, st_p, st_m, cps, cpsph):
    return 4 * (x * cps * (st_m / (1 - ct_m))
                + (1 - x) * cpsph * (st_p / (1 - ct_p))
                ) ** 2


def _I_B(x, w, ct_p, ct_m, st_p, st_m, cph, cps, cpsph):
    q2 = _q2(x, w, ct_p, ct_m, st_p, st_m, cph)
    return -q2 * (cps * (st_m / (1 - ct_m))
                  - cpsph * (st_p / (1 - ct_p))
                  ) ** 2


def _I_C(x, ct_p, ct_m, st_p, st_m, cph):
    return (-(st_m / (1 - ct_m)) * (st_p / (1 - ct_p))
            * ((x / (1 - x)) * st_p / st_m
               + ((1 - x) / x) * st_m / st_p
               + 2 * cph
               )
            )


def _I(X, w, phi, psi, part: str):
    x, th_p, th_m = X  #[:, 0], X[:, 1], X[:, 2]
    st_p = sin(th_p)
    ct_p = cos(th_p)
    st_m = sin(th_m)
    ct_m = cos(th_m)
    cph = cos(phi)
    cps = cos(psi)
    cpsph = cos(psi + phi)
    if part == 'A':
        blah = _I_A(x, ct_p, ct_m, st_p, st_m, cps, cpsph)
    elif part == 'B':
        blah = _I_B(x, w, ct_p, ct_m, st_p, st_m, cph, cps, cpsph)
    elif part == 'C':
        blah = _I_C(x, ct_p, ct_m, st_p, st_m, cph)
    else:
        return 0.
    return _stuff(blah, x, w, ct_p, ct_m, st_p, st_m, cph)


def get_integrand(phi, psi, w=100., parts='ABC'):
    def _sigma(X):
        return sum((_I(X, w=w, phi=phi, psi=psi, part=_p) for _p in parts))
    # _sigma = vegas.batchintegrand(_sigma)
    return _sigma

