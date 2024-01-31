import numpy as np
import matplotlib.pyplot as plt


def LJcut126(epsilon, sigma, Rc, r):
    """
    Calculate the Lennard-Jones 12-6 potential for a given distance.

    Parameters:
    epsilon (float): The depth of the potential well.
    sigma (float): The distance at which the potential energy is zero.
    Rc (float): The cutoff distance for the potential calculation.
    r (float): The distance at which to calculate the potential.

    Returns:
    float: The Lennard-Jones 12-6 potential energy at distance r, or 0 if r >= Rc.
    """
    if (r >= Rc):
        return 0
    else:
        return 4*epsilon*(((sigma/r)**12) - ((sigma/r)**6))


def LJcut93(epsilon, sigma, Rc, r):
    """
    Compute the Lennard-Jones potential 9-3 energy at a given distance.

    Args:
        epsilon (float): The depth of the potential well.
        sigma (float): The distance at which the potential energy is zero.
        Rc (float): The cutoff distance for the potential calculation.
        r (float): The distance at which to calculate the potential energy.

    Returns:
        float: The Lennard-Jones 9-3 potential energy at distance r, or 0 if r >= Rc.
    """
    if (r >= Rc):
        return 0
    else:
        return epsilon*((2/15)*((sigma/r)**9) - ((sigma/r)**3))


T = 310
kB = 1

kBT = kB*T

epsilon = 1500
sigma = 2.5
Rc = 10

threshold = 2.0

r = np.linspace(0.01, 15, 1000)

LJ126 = np.zeros((len(r)))
LJ93 = np.zeros((len(r)))

for i in range(len(r)):
    LJ126[i] = (1/kBT)*LJcut126(epsilon, sigma, Rc, r[i])
    LJ93[i] = (1/kBT)*LJcut93(epsilon, sigma, Rc, r[i])

# plt.figure(tight_layout=True, figsize=(6, 4))

# plt.xlabel("r (nm)")
# plt.ylabel("V(r) (kB T)")

# plt.xlim(left=0, right=max(r))
# plt.ylim(top=1, bottom=-10)

# plt.axhline(0, color="black", linewidth=1.0)

# plt.plot(r, LJ93, label="LJ 9-3")
# plt.plot(r, LJ126, '--', label="LJ 12-6")

# plt.legend(loc="lower right")
# plt.savefig("wall_potential.pdf", bbox_inches='tight')

##########################################################

plt.clf()

plt.figure(tight_layout=True, figsize=(6, 4))

plt.xlabel("r (nm)")
plt.ylabel("V(r) (kB T)")

plt.xlim(left=0, right=max(r))
plt.ylim(top=1, bottom=-6)

plt.axhline(0, color="black", linewidth=1.0)

threshold_line = []
y_vals = np.linspace(-10, 1, 1000)
for y in y_vals:
    threshold_line.append(threshold)

plt.plot(r, LJ93, label="LJ 9-3")
plt.plot(r, LJ126, 'g--', label="LJ 12-6")
plt.plot(threshold_line, y_vals, 'r--', label="Threshold", linewidth=0.7)

plt.minorticks_on()

plt.legend(loc="lower right")
plt.savefig("wall_potential.pdf", bbox_inches='tight')
