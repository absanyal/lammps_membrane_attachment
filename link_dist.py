import numpy as np
import matplotlib.pyplot as plt


class cylinder():
    def __init__(self, radius, rA, rB):
        self.radius = radius
        self.rA = rA
        self.rB = rB


def norm(A):
    sum = 0
    for i in range(len(A)):
        sum += A[i] * A[i]
    return np.sqrt(sum)


def distance_from_axis(cyl, rP):
    e = cyl.rA - cyl.rB
    d = norm((np.cross(e, rP - cyl.rA)))/norm(e)
    return d


def distance_from_surface(cyl, rP):
    return cyl.radius - distance_from_axis(cyl, rP)


rA = np.array([0.0, 175.0, 175.0])
rB = np.array([1000.0, 175.0, 175.0])
c = cylinder(175.0, rA, rB)

data = np.loadtxt("link_pos.txt", unpack=True)

# print(len(data))

num_linkers = int((len(data) - 1) / 3)

# print("Number of linkers =", num_linkers)

t = data[0]

avg_dist = np.zeros((len(t),))

plt.figure(tight_layout=True)

for i in range(1, num_linkers+1):
    
    dist_list = []
    
    x = data[3*(i-1)+1]
    y = data[3*(i-1)+2]
    z = data[3*(i-1)+3]

    for t_step in range(len(t)):
        rP = np.array([x[t_step], y[t_step], z[t_step]])
        dist = distance_from_surface(c, rP)
        
        dist_list.append(dist)
    
    dist_list = np.array(dist_list)
    plt.plot(t, dist_list, label="Linker {}".format(i), linewidth=0.5)
    avg_dist += dist_list
    

avg_dist *= (1/num_linkers)

plt.xlabel(r'$t/\tau$', fontsize=18)
plt.ylabel("Distance from cell membrane (nm)", fontsize=16)
plt.ylim(bottom=0)

plt.plot(t, avg_dist, 'k--', label="Average distance", linewidth=1.0)
plt.legend()

plt.savefig("link_dist.pdf")

# plt.show()
