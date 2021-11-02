import numpy as np

import csv
from functools import partial
from sys import argv

from depaola_cport import PARAM

import matplotlib.pyplot as plt

from numpy import (    pi, sin, cos, sqrt,
  linspace)
from cycuba import Vegas, Suave, Divonne, Cuhre

OMEGA_V = [100.]

MASS_MEV = 0.511  # mass, electron/positron

PARAM['epsrel'] = 1e-2

def PhiXu(E, Ep, p, pp, t, tp, phi, omega):
    st = sin(t)
    stp = sin(tp)
    ct = cos(t)
    ctp = cos(tp)

    b  = 0. if E  < MASS_MEV else sqrt(1. - pow(MASS_MEV/E,  2))
    bp = 0. if Ep < MASS_MEV else sqrt(1. - pow(MASS_MEV/Ep, 2))
    qa = E*Ep*(1. - b*bp*(sin(t)*sin(tp)*cos(phi)
                   + cos(t)*cos(tp)))
    qb = (omega*E *(b *cos(t)  - 1)
        + omega*Ep*(bp*cos(tp) - 1) + MASS_MEV**2)
    q2 = -2.*(qa + qb)

    Xu = (
        (pp*stp/(Ep - pp*ctp))**2 * (4*E **2 - q2)
      + ( p*st /(E  -  p*ct ))**2 * (4*Ep**2 - q2)
      + 2*(4*E*Ep + q2 - 2*omega**2)/(E  - p *ct )/(Ep - pp*ctp)*( pp*stp     *  p*st * cos(phi))
      -   (              2*omega**2)/(Ep - pp*ctp)/(E  - p *ct )*((pp*stp)**2 + (p*st)**2)
    )

    Phi = -MASS_MEV**2 * p * pp / (2*pi)**2 / omega**3 / q2**2
    
    return Phi*Xu*st*stp
    
def PhiXp(E, Ep, p, pp, t, tp, phi, psi, omega):
    st = sin(t)
    stp = sin(tp)
    ct = cos(t)
    ctp = cos(tp)
    
    b  = 0. if E  < MASS_MEV else sqrt(1. - pow(MASS_MEV/E,  2))
    bp = 0. if Ep < MASS_MEV else sqrt(1. - pow(MASS_MEV/Ep, 2))
    qa = E*Ep*(1. - b*bp*(sin(t)*sin(tp)*cos(phi)
                   + cos(t)*cos(tp)))
    qb = (omega*E *(b *cos(t)  - 1)
        + omega*Ep*(bp*cos(tp) - 1) + MASS_MEV**2)
    q2 = -2.*(qa + qb)

    Xp = (
        (pp*stp/(Ep - pp*ctp))**2 * (4*E **2 - q2) * cos(2*psi)
      + ( p*st /(E  - p *ct ))**2 * (4*Ep**2 - q2) * cos(2*(phi + psi))
      + 2*(4*E*Ep + q2             )/(E  - p *ct )/(Ep - pp*ctp)*( pp*stp     *  p*st * cos(phi + 2*psi))
    )

    Phi = -MASS_MEV**2 * p * pp / (2*pi)**2 / omega**3 / q2**2
    
    return Phi*Xp*st*stp
    
def MOM(E, m=MASS_MEV):
    if not m:
        return E
    if E <= m:
        return 0.
    return sqrt(E**2 - m**2)

def main(method, verbose=False):
    phi_break, n1, n = pi*3/4., 0, 19
    phi_v = [phi_break * i / n1 for i in range(n1)]
    phi_v += [(pi - phi_break) * i / (n - n1) + phi_break for i in range(n - n1)]
    phi_v += [pi]

    #n = 10; psi_v = [j*pi/n for j in range(n//2)]
    #psi_v += [pi]
    psi_v = [0., pi/2.]
    #psi_v = []

    PhiXu_points = np.zeros((len(phi_v), 2))
    PhiXp_points = np.zeros((len(phi_v), len(psi_v), 2))
    for omega in OMEGA_V:
        print(f'omega = {omega:.1f}')
        for i, phi in enumerate(phi_v):
            print(f'{phi=:.2f}')
            def integrand(E, t, tp):
                # Scale energy from 0,1 to 0,omega
                E *= omega
                # Scale thetas from 0,1 to 0,Pi
                t *= pi/2
                tp *= pi/2
                return [PhiXu(E, omega-E, MOM(E), MOM(omega-E), t, tp, phi, omega)]
            cycuba_result = method(integrand, **PARAM)
            PhiXu_points[i, 0] = cycuba_result[0][0]
            PhiXu_points[i, 1] = cycuba_result[1][0]
            for j, psi in enumerate(psi_v):
                print(f'{psi=:.2f}')
                def integrand(E, t, tp):
                    # Scale energy from 0,1 to 0,omega
                    E *= omega
                    # Scale thetas from 0,1 to 0,Pi
                    t *= pi/2
                    tp *= pi/2
                    return [PhiXp(E, omega-E, MOM(E), MOM(omega-E), t, tp, phi, psi, omega)]
                cycuba_result = method(integrand, **PARAM)
                PhiXp_points[i, j, 0] = cycuba_result[0][0]
                PhiXp_points[i, j, 1] = cycuba_result[1][0]
    fig, ax = plt.subplots()
    ax.errorbar(phi_v, PhiXu_points[:, 0], yerr=PhiXu_points[:, 1], label=f"ω = {omega:.0f} MeV, Unpol.")
    for j, _psi in enumerate(psi_v):
        ax.errorbar(phi_v, PhiXp_points[:, j, 0], yerr=PhiXp_points[:, j, 1], label=f'ω = {omega:.0f} MeV, ψ = {_psi:.2f}')
    ax.grid()
    ax.legend(loc='best')
    plt.show()


if __name__ == '__main__':
    main(Cuhre, True)
