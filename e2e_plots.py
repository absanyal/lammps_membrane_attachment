import numpy as np
import matplotlib.pyplot as plt

info = np.loadtxt('data/info.txt')

num_monomers = int(info[0])
num_linkers = int(info[1])
num_skip = int(info[2])

t, e1x, e1y, e1z, e2x, e2y, e2z = np.loadtxt('data/e2e_pos.txt', unpack=True)

hv0 = np.array([e1x[0] - e2x[0], e1y[0] - e2y[0], e1z[0] - e2z[0]])

ang_deviation_list = np.zeros(len(t))

for i in range(len(t)):
    hv = np.array([e1x[i] - e2x[i], e1y[i] - e2y[i], e1z[i] - e2z[i]])
    angle = np.arccos(np.dot(hv0, hv) / (np.linalg.norm(hv0) * np.linalg.norm(hv)))
    ang_deviation_list[i] = angle

with open('data/ang_deviation.txt', 'w') as f:
    f.write('#t ang_deviation(rad) ang_deviation(deg)\n')
    for i in range(len(t)):
        f.write('{}\t{:.4f}\t{:.4f}\n'.format(t[i], ang_deviation_list[i], np.rad2deg(ang_deviation_list[i])))

t, ang_deviation_list_rad, ang_deviation_list_deg = np.loadtxt('data/ang_deviation.txt', unpack=True)

plt.figure(tight_layout=True)

plt.plot(t, ang_deviation_list_deg, 'k-', linewidth=1, label='Angle deviation')
plt.xlabel(r'$t/\tau$', fontsize=18)
plt.ylabel("Angular deviation", fontsize=18)

plt.savefig('plots/ang_deviation.pdf')

with open('data/e2e_dist.txt', 'w') as f:
    f.write('#t e2e_dist\n')
    for i in range(len(t)):
        hv = np.array([e1x[i] - e2x[i], e1y[i] - e2y[i], e1z[i] - e2z[i]])
        f.write('{}\t{:.4f}\n'.format(t[i], np.linalg.norm(hv)))

t, e2e_dist_list = np.loadtxt('data/e2e_dist.txt', unpack=True)
print("End to end distance: ")
print("Mean: {:.4f}".format(np.mean(e2e_dist_list)))
print("Std: {:.4f}".format(np.std(e2e_dist_list)))

e2e_mean_line = np.mean(e2e_dist_list) * np.ones(len(t))

plt.clf()
plt.cla()

plt.figure(tight_layout=True)

plt.plot(t, e2e_dist_list, 'k-', linewidth=1, label='End to end distance')
plt.plot(t, e2e_mean_line, 'r--', linewidth=2, label='Mean distance')
plt.xlabel(r'$t/\tau$', fontsize=18)
plt.ylabel("End to end distance", fontsize=18)
plt.legend()

plt.savefig('plots/e2e_dist.pdf')

    