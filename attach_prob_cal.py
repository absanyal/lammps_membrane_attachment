import numpy as np
import matplotlib.pyplot as plt
import math
from numpy import exp
from math import gamma

T = 110
E = 500
n = 100 # number of linkers
continuous_val = 0

#--------------------
T0 = 310
kB = 1
kBT = kB*T0

T_entered = T

T = T/T0
beta = 1.0 / (kB * T)
E = -E / kBT
#--------------------

def attach_probability(n, m, continuous=False):
    """
    Calculate the attachment probability.

    Parameters:
    n (int): Total number of particles.
    m (int): Number of attached particles.
    continuous (bool, optional): Flag indicating whether continuous or discrete attachment is considered. 
                                 Defaults to False.

    Returns:
    float: The attachment probability.

    """
    if (continuous == True):
        nCm = gamma(n + 1) / ((gamma(m + 1)) * (gamma(n - m + 1)))
    else:
        nCm = math.comb(n, m)
    return (exp(-E * m * beta) / (-1 + pow(1 + exp(-E * beta), n))) * nCm


print("E =", E, "kB T")
print("T =", T_entered, "kelvin")
print("n =", n)

if (continuous_val == True):
    m = np.linspace(0, n, 5000, endpoint=True)
else:
    m = list(range(n+1))


# m = list(range(n+1))
m = np.array(m)

P = np.zeros(len(m))

for i in range(len(m)):
    P[i] = attach_probability(n, m[i], continuous_val)

plt.figure(tight_layout=True)

if (continuous_val == True):
    plt.plot(m, P)
else:
    plt.bar(m, P)

plt.xlabel("Number of attached linkers")
plt.ylabel("Probability of attachment")

# plt.yscale("log")

plt.ylim(bottom=0)

# plt.show()
plt.savefig("plots/attach_prob.pdf", bbox_inches='tight')
