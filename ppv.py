# On the stopping of fast particles and on the creation of positive electrons
# H. Bethe and W. Heitler
# February 27, 1934
# 1934.0140

# Methods and parameters are scaled as follows:
# Incoming photon energy k and electron/positron energies are in units
# of the electron mass.
# dΦ is in units of α*Z**2*e**4.
# β := p/E = sqrt(E**2 - m**2)/E

import argparse
import csv
from functools import partial
from sys import argv

from numpy import (pi, sin, cos, sqrt)

import vegas


MASS_MEV = 0.511  # mass, electron/positron

def _q2(E, cθ, cθp, k, φ) -> float:
    # sθ, sθp = sin(θ), sin(θp)
    # cθ, cθp = cos(θ), cos(θp)
    sθ, sθp = sqrt(1 - cθ**2), sqrt(1 - cθp**2)
    cφ = cos(φ)
    Ep = k - E
    β, βp = sqrt(E**2  - 1)/E, sqrt(Ep**2 - 1)/Ep
    ###
    _q2 = (
        -2*(
            E*Ep*(1 - β*βp*(sθ*sθp*cφ
                            + cθ*cθp))  # qa in depaola_cport
            - k*E*(1 - β*cθ)
            - k*Ep*(1 - βp*cθp)
            + 1                         # qb in depaola_cport
        )
    )
    return _q2

def _dΦ(E, θ, θp, k, φ, ψ) -> float:
    """BH Equation 20"""
    sθ, sθp = sin(θ), sin(θp)
    cθ, cθp = cos(θ), cos(θp)
    q2 = _q2(E, cθ, cθp, k, φ)
    q4 = q2**2
    # sθ, sθp = sqrt(1 - cθ**2), sqrt(1 - cθp**2)
    cφ = cos(φ)
    cψ = cos(ψ)
    cφψ = cos(φ + ψ)
    Ep = k - E
    β, βp = sqrt(E**2  - 1)/E, sqrt(Ep**2 - 1)/Ep
    β2, βp2 = β**2, βp**2
    α, αp = β*sθ/(1 - β*cθ), βp*sθp/(1 - βp*cθp)
    α2, αp2 = α**2, αp**2
    #r, rp = E/Ep, Ep/E
    ### A, B, C, and D are the terms in the curly braces in BH Eq. 20
    A =    αp2 * cψ **2 * (4*E **2 - q2)  # βp2 * sθp**2 * (4*E **2 - q2) / (1 - βp*cθp)**2
    B =    α2  * cφψ**2 * (4*Ep**2 - q2)  # β2  * sθ **2 * (4*Ep**2 - q2) / (1 - β *cθ )**2
    C = 2*α*αp*cψ*cφψ   * (4*E*Ep  + q2)  # 2*β*βp*sθ*sθp*cφ*(4*E*Ep + q2)/(1 - β*cθ)/(1 - βp*cθp)
    #fixme D = -2*k**2 * (rp*βp2*sθp**2 + r*β2*sθ**2 + 2*β*βp*sθ*sθp*cφ)/(1 - β*cθ)/(1 - βp*cθp)
    #fixme _dΦ = -(E*Ep*β*βp/k**3)/(2*pi)/q2**2 * (A+B+C+D)
    ###
    # s1 and s2 are from depaola_cport
    s1 = E*Ep*β*βp*sθ*sθp/k**3/q4*(A+B+C)
    s2 = α*αp/k/q4*(E**2*β2*sθ**2 + Ep**2*βp2*sθp**2 + 2*E*Ep*β*βp*sθ*sθp*cφ)
    #print(f"{s1=:.3g}\n{s2=:.3g}")
    _dΦ = (s2 - s1)/4/pi**2
    return _dΦ

def _Φ(k, φ, ψ,
        adapt_nitn: int = 10, adapt_neval: int = 10000,
        nitn: int = 10, neval: int = 2000,
        alpha: float = .5,
        verbose: bool = False,
) -> (float, float):
    def integrand(X):
        E, θ, θp = X
        return [_dΦ(E, θ, θp, k=k, φ=φ, ψ=ψ)]
    integ = vegas.Integrator([[1, k-1], [0, pi], [0, pi]])
    # step 1 -- adapt to sigma; discard results
    adapt = integ(integrand, nitn=adapt_nitn, neval=adapt_neval, alpha=alpha)
    if verbose:
        print('ADAPTATION:')
        print(adapt.summary())
    # step 2 -- integ has adapted to f; keep results
    result = integ(integrand, nitn=nitn, neval=neval, alpha=alpha)
    if verbose:
        print('ESTIMATION:')
        print(result.summary())
        print('    Result = %s    Q = %.2f' % (result, result.Q))
        print('    Result = %f +- %f' % (result[0].mean, result[0].sdev))

    mu, sig, Q = result[0].mean, result[0].sdev, result.Q
    if verbose:
        print('>> sigma = %.6g\n'
              '        +- %.6g\n'
              '   Q = %.2g' % (mu, sig, Q))
    return mu, sig

def _main(omega=None, nphi=60, phi_min=0, phi_break=3, verbose=False):
    PHI = [phi_break * i / (nphi//2) for i in range(nphi//2)]
    PHI += [(pi - phi_break) * i / (nphi//2) + phi_break for i in range(nphi//2)]
    PHI += [pi]
    PSI = [0, pi/2]
    for ω in omega or [100]:
        for φ in PHI:
            for ψ in PSI:
                result = _Φ(k=ω/MASS_MEV, φ=φ, ψ=ψ, verbose=verbose)
                σ, ε = result[0]/MASS_MEV**3, result[1]/MASS_MEV**3
                if verbose:
                    print(f"{ω=:>3d}   {φ=:<4.2f}   {ψ=:<4.2f}   {σ=:<9.6g} +- {ε:<8.6e}")
                yield ω, φ, ψ, σ, ε

def main(omega=None, nphi=60, phi_min=0, phi_break=3, verbose=False):
    with open("ppv.csv", 'w', newline='', encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow('ω φ ψ σ ε'.split(' '))
        for ω, φ, ψ, σ, ε in _main(omega, nphi, phi_min, phi_break, verbose):
            writer.writerow(['%g' % _v for _v in [ω, φ, ψ, σ, ε]])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate the cross-section for electron/positron pair-production by polarized photons.")
    parser.add_argument("-o", "--omega", type=float, action="append", help="Incoming photon energy(s) in units of MeV. Can be repeated. If not given, 100 MeV will be used")
    parser.add_argument("-n", "--nphi", type=int, default=60, help="Number of values of phi (φ) in the range from phi-min to pi")
    parser.add_argument("-p", "--phi-min", type=float, default=0, help="Minimum value of phi (φ) in the range")
    parser.add_argument("-b", "--phi-break", type=float, default=3, help="Break in the range of values of phi (φ). Will calculate N/2 points in [phi-min, phi-break] and N/2 in [phi-break, pi]")
    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    main(**vars(parser.parse_args()))
