# On the stopping of fast particles and on the creation of positive electrons
# H. Bethe and W. Heitler
# February 27, 1934
# 1934.0140

# Dunder methods and parameters are scaled as follows:
# Incoming photon energy k is in units of the electron mass m.
# Electron/positron energies are in units of k (in units of m).
# Φ is in units of α*Z**2*e**4.

def q2(k, E, θ, θp, φa) -> float:
    q2 = (
        -2*(
            E*(1 - E)*(1 - sin(θ)*sin(θp)*cos(φ) - cos(θ)*cos(θp))
            - E*(1 - cos(θp))
            - (1 - E)*(1 - cos(θ))
            + 1
        )
    )
    return q2

def _dΦ(k, E, θ, θp, φ, ψ) -> float:
    """BH Equation 20"""
    
