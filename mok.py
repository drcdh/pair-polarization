# Motz, Olsen & Koch Formula 3D-3132
# Born approximation formula with dependance on (P_L, e): Screened point nucleus
# Conditions of validity, Table 6.01: A, E
# References May (1951), Gluckstern et al. (1951), Gluckstern & Hull (1953),
#  Bobel (1957), Claesson (1957), McVoy (1957), Fronsdale & Uberali (1958),
#  Banergee (1958)

import vegas
import math


def ds3(x):
    dx2 = 0
    for d in range(4):
        dx2 += (x[d] - 0.5) ** 2
    return math.exp(-dx2 * 100.) * 1013.2118364296088
