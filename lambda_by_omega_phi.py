import csv

import numpy as np

from numpy import pi

from cycuba import Cuhre as CyCubaMethod

from bernard import MOM, PhiXp, PhiXu
from depaola_cport import PARAM

def __gen_integrand_PhiXp(omega, phi, psi):
    def integrand(E, t, tp):
        # Scale energy from 0,1 to 0,omega
        E *= omega
        # Scale thetas from 0,1 to 0,Pi
        t *= pi/2
        tp *= pi/2
        return [PhiXp(E, omega-E, MOM(E), MOM(omega-E), t, tp, phi, psi, omega)]
    return integrand

def __gen_integrand_PhiXu(omega, phi):
    def integrand(E, t, tp):
        # Scale energy from 0,1 to 0,omega
        E *= omega
        # Scale thetas from 0,1 to 0,Pi
        t *= pi/2
        tp *= pi/2
        return [PhiXu(E, omega-E, MOM(E), MOM(omega-E), t, tp, phi, omega)]
    return integrand

def int_PhiXp(omega, phi, psi):
    cycuba_result = CyCubaMethod(__gen_integrand_PhiXp(omega, phi, psi), **PARAM)
    return cycuba_result[0][0], cycuba_result[1][0]

def int_PhiXu(omega, phi):
    cycuba_result = CyCubaMethod(__gen_integrand_PhiXu(omega, phi), **PARAM)
    return cycuba_result[0][0], cycuba_result[1][0]

@np.vectorize
def lambda_factor(omega, phi):
    print(f'{omega=:.2f} MeV, {phi=:.2f}')
    u, _ = int_PhiXu(omega, phi)
    pp, _ = int_PhiXp(omega, phi, pi/2)
    pm, _ = int_PhiXp(omega, phi, 0)
    lambda_plus = pp/u
    lambda_minus = pm/u
    return lambda_minus, lambda_plus, lambda_plus-lambda_minus

@np.vectorize
def int_PhiXu_PhiXp(omega, phi):
    print(f'{omega=:.2f} MeV, {phi=:.2f}')
    u, u_err = int_PhiXu(omega, phi)
    pp, pp_err = int_PhiXp(omega, phi, pi/2)
    pm, pm_err = int_PhiXp(omega, phi, 0)
    return u, u_err, pp, pp_err, pm, pm_err

omegas = np.logspace(1, 3, num=100)
phis = np.linspace(0, pi, num=100)

mesh_omegas, mesh_phis = np.meshgrid(omegas, phis)
#lambda_by_omega_phi = lambda_factor(mesh_omegas, mesh_phis)
#lambda_by_omega_phi = np.stack(lambda_by_omega_phi)
# indexed like [:, phi, omega]

#np.savez('lambda_omega_phi.npz', lambdas=lambda_by_omega_phi, omegas=mesh_omegas, phis=mesh_phis)

#f = open('lambda_by_omega_phi.csv', 'w')
f = open('PhiXu_PhiXp_by_omega_phi.csv', 'w')
csvf = csv.writer(f)
#csvf.writerow('omega[MeV] phi lambda'.split())
csvf.writerow('omega[MeV] phi PhiXu PhiXuErr PhiXpp PhiXppErr PhiXpm PhiXpmErr'.split())
for i, omega in enumerate(omegas):
    for j, phi in enumerate(phis):
        row = [omega, phi]
        row.extend(list(int_PhiXu_PhiXp(omega, phi)))
#        csvf.writerow(f'{mesh_omegas[i,j]} {mesh_phis[i,j]} {lamb}'.split())
        csvf.writerow([f'{_r}' for _r in row])
f.close()
