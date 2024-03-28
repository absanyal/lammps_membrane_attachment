import numpy as np
import matplotlib.pyplot as plt
import cylindermath as cm

rA = np.array([0.0, 175.0, 175.0])
rB = np.array([1000.0, 175.0, 175.0])
c = cm.cylinder(175.0, rA, rB)

t_list, comx, comy, comz = np.loadtxt('data/com_pos.txt', unpack=True)

com_dist_list_str = "data/com_dist_list.txt"

window = 11
print("Smoothing window = {}".format(window))

detection_window = 21
print("Detection window = {}".format(detection_window))

threshold = 0.25
print("Threshold = {:.2f}".format(threshold))

hitting_distance = 2.5
print("Hitting distance = {:.2f}".format(hitting_distance))

info = np.loadtxt('data/info.txt')

num_monomers = int(info[0])
num_linkers = int(info[1])
num_skip = int(info[2])

dist_list = np.zeros_like(t_list)

for t_i in range(len(t_list)):
    rP = np.array([comx[t_i], comy[t_i], comz[t_i]])
    dist_list[t_i] = cm.distance_from_surface(c, rP)

plt.figure(tight_layout=True)
plt.plot(t_list, dist_list, label='CoM', color='k', linewidth=0.6)
plt.xlabel(r'$t/\tau$', fontsize=18)
plt.ylabel(r'Distance from cell membrane (nm)', fontsize=16)
plt.legend()
plt.ylim(bottom=0.0)
plt.savefig('plots/com_dist.pdf')
    
