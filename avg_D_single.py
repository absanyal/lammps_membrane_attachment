import numpy as np
import matplotlib.pyplot as plt

D_list = np.loadtxt('runs_measurements/D_single_avg_list.txt')

# print(D_list)

D_avg = np.mean(D_list)

with open('runs_measurements/D_single_avg.txt', 'w') as fname:
    fname.write('#D_avg\n')
    fname.write('{}\n'.format(D_avg))

print("D_avg = {}".format(D_avg))

plt.figure(tight_layout=True)

plt.hist(D_list, bins='auto', color='b', alpha=0.7, edgecolor='k', density=True, rwidth=0.85)
plt.xlabel(r'$D/\sigma^2\tau^{-1}$', fontsize=18)
plt.ylabel('Frequency', fontsize=18)
plt.savefig('runs_measurements/D_single_avg_hist.pdf')