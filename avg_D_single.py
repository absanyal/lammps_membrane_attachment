import numpy as np
import matplotlib.pyplot as plt

D_list = np.loadtxt('runs_measurements/D_single_avg_list.txt')

# print(D_list)

D_avg = np.mean(D_list)
D_std = np.std(D_list)
n_trials = len(D_list)

with open('runs_measurements/D_single_avg.txt', 'w') as fname:
    fname.write('#D_avg\n')
    fname.write('{}\n'.format(D_avg))

print("Number of trials = {}".format(n_trials))
print("D_avg = {}".format(D_avg))
print("D_std = {}".format(D_std))

plt.figure(tight_layout=True)

D_min = np.min(D_list)
D_max = np.max(D_list)
D_step = 10

bins_list = np.arange(D_min, D_max, D_step)

plt.hist(D_list, bins=bins_list, color='b', alpha=0.7, edgecolor='k', density=True, rwidth=0.85)
plt.xlabel(r'$D/\sigma^2\tau^{-1}$', fontsize=18)
plt.ylabel('Frequency', fontsize=18)
plt.title(r'Trials = {}, $\langle D \rangle = {:.2f}$, $\sigma = {:.2f}$'.format(n_trials, D_avg, D_std), fontsize=18)
plt.savefig('runs_measurements/D_single_avg_hist.pdf')