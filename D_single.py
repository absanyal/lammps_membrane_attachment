import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def f(x, m, c):
    return m * x + c

t, x, y, z = np.loadtxt('data/single_pos.txt', unpack=True)

ds2 = np.zeros(len(t))

for i in range(len(t)):
    ds2[i] = (x[i]-x[0])**2 + (y[i]-y[0])**2 + (z[i]-z[0])**2

with open('data/single_displacement.txt', 'w') as file:
    file.write('#t ds2\n')
    for i in range(len(t)):
        file.write('{}\t{}\n'.format(t[i], ds2[i]))

t, ds2 = np.loadtxt('data/single_displacement.txt', unpack=True)

plt.figure(tight_layout=True)
plt.plot(t, ds2, 'k-', linewidth=0.5)
plt.xlabel(r'$t/\tau$', fontsize=18)
plt.ylabel(r'$\Delta s^2$', fontsize=18)
plt.savefig('plots/single_displacement.pdf')

t, x, y, z = np.loadtxt('data/single_pos.txt', unpack=True)

t_max = len(t)
sample_window_fraction = 0.02
sample_window = int(t_max * sample_window_fraction)
print("Total iterations: {}".format(t_max))
print("Sample window: {} iterations".format(sample_window))

ds2_avg_sampled = np.zeros(sample_window)
t0_iter_list = np.arange(0, t_max-sample_window)
print("Samples taken: {}".format(len(t0_iter_list)))

t_shortened = t[:sample_window]

dx = x - x[0]
dy = y - y[0]
dz = z - z[0]

for t0 in t0_iter_list:
    init_dx = dx[t0]
    init_dy = dy[t0]
    init_dz = dz[t0]
    ds2_avg_sampled += (dx[t0:t0+sample_window]-init_dx)**2 + (dy[t0:t0+sample_window]-init_dy)**2 + (dz[t0:t0+sample_window]-init_dz)**2

ds2_avg_sampled /= len(t0_iter_list)

with open('data/single_disp_sampled.txt', 'w') as file:
    file.write("# t0 \t <ds^2>\n")
    for i in range(sample_window):
        file.write("{} \t {}\n".format(t_shortened[i], ds2_avg_sampled[i]))

fit_params, fit_cov = curve_fit(f, t_shortened, ds2_avg_sampled)

print('Fitted parameters: m = {:.4e}, c = {:.4e}'.format(fit_params[0], fit_params[1]))
D = fit_params[0] / 6
print('Diffusion coefficient: D = {:.4e}'.format(D))

with open('data/D_single_value.txt', 'w') as file:
    file.write('# D\n')
    file.write('{}\n'.format(D))

fitline = f(t_shortened, *fit_params)

plt.figure(tight_layout=True)
plt.plot(t_shortened, ds2_avg_sampled, 'k-', label='Data')
plt.plot(t_shortened, fitline, 'r--', label='Fit')
plt.xlabel(r'$t/\tau$', fontsize=18)
plt.ylabel(r'$\langle \Delta s^2 \rangle$', fontsize=18)
plt.legend(fontsize=14)
plt.title(r'$D = {:.2f}$'.format(D))
plt.savefig('plots/single_disp_sampled.pdf')