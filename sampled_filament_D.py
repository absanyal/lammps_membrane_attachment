import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def f(x, m, c):
    return m * x + c

t, x, y, z = np.loadtxt('data/com_pos.txt', unpack=True)

dx = x - x[0]
dy = y - y[0]
dz = z - z[0]

t_max = len(t)
sample_window_fraction = 0.02
sample_window = int(t_max * sample_window_fraction)
print("Total iterations: {}".format(t_max))
print("Sample window: {} iterations".format(sample_window))

ds_sq_avg_sampled = np.zeros(sample_window)
t0_iter_list = np.arange(0, t_max-sample_window)
print("Samples taken: {}".format(len(t0_iter_list)))

t_shortened = t[:sample_window]

for t0_i in t0_iter_list:
    init_dx = dx[t0_i]
    init_dy = dy[t0_i]
    init_dz = dz[t0_i]
    
    ds_sq_avg_sampled += (dx[t0_i:t0_i+sample_window]-init_dx)**2 + (dy[t0_i:t0_i+sample_window]-init_dy)**2 + (dz[t0_i:t0_i+sample_window]-init_dz)**2
    
ds_sq_avg_sampled /= len(t0_iter_list)

with open('data/D_com_sampled.txt', 'w') as file:
    file.write("# t0 \t <ds^2>\n")
    for i in range(sample_window):
        file.write("{} \t {}\n".format(t_shortened[i], ds_sq_avg_sampled[i]))

fit_params, fit_cov = curve_fit(f, t_shortened, ds_sq_avg_sampled)
print('Fitted parameters: m = {:.4e}, c = {:.4e}'.format(fit_params[0], fit_params[1]))
D = fit_params[0]/(6)

print('Diffusion coefficient: D = {:.4e}'.format(D))
print('or, approximately D = {:.4f}'.format(D))

fitline = f(t_shortened, *fit_params)

plt.figure(tight_layout=True)
plt.plot(t_shortened, ds_sq_avg_sampled, 'k-', label='Data')
plt.plot(t_shortened, fitline, 'r--', label='Fit')
plt.xlabel(r'$t/\tau$', fontsize=18)
plt.ylabel(r'$\langle \Delta s^2 \rangle$', fontsize=18)
plt.title(r"$D_{{com}} = {:.4f}$".format(D))
plt.legend(fontsize=14)
plt.savefig('plots/D_com_sampled.pdf')


