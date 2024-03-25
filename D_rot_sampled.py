import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def f(x, m, c):
    return m * x + c

t, ang_dev_rad, ang_dev_deg = np.loadtxt(
    'data/ang_deviation.txt', unpack=True)

t_max = len(t)
sample_window_fraction = 0.02
sample_window = int(t_max * sample_window_fraction)
print("Total iterations: {}".format(t_max))
print("Sample window: {} iterations".format(sample_window))

ang_dev_sq_avg_sampled = np.zeros(sample_window)
t0_iter_list = np.arange(0, t_max-sample_window)
print("Samples taken: {}".format(len(t0_iter_list)))

t_shortened = t[:sample_window]

# print("t0_iter_list: {}".format(t0_iter_list))

for t0_i in t0_iter_list:
    init_angle = ang_dev_rad[t0_i]
    ang_dev_sq_avg_sampled += (ang_dev_rad[t0_i:t0_i+sample_window]-init_angle)**2

ang_dev_sq_avg_sampled /= len(t0_iter_list)

# plt.plot(t_shortened, ang_dev_sq_avg_sampled, 'k-', label='Data')
# plt.show()

with open('data/D_rot_sampled.txt', 'w') as file:
    file.write("# t0 \t <theta^2>\n")
    for i in range(sample_window):
        file.write("{} \t {}\n".format(t_shortened[i], ang_dev_sq_avg_sampled[i]))

fit_params, fit_cov = curve_fit(f, t_shortened, ang_dev_sq_avg_sampled)
print('Fitted parameters: m = {:.4e}, c = {:.4e}'.format(fit_params[0], fit_params[1]))
print('Diffusion coefficient: D = {:.4e}'.format(fit_params[0]/(2)))
print('or, approximately D = {:.4f}'.format(fit_params[0]/(2)))

fitline = f(t_shortened, *fit_params)

plt.figure(tight_layout=True)
plt.plot(t_shortened, ang_dev_sq_avg_sampled, 'k-', label='Data')
plt.plot(t_shortened, fitline, 'r--', label='Fit')
plt.xlabel(r'$t/\tau$', fontsize=18)
plt.ylabel(r'$\langle \theta^2 \rangle$', fontsize=18)
plt.title(r"$D_{{rot}} = {:.4f}$".format(fit_params[0]/(2)))
plt.legend(fontsize=14)
plt.savefig('plots/D_rot_sampled.pdf')

with open('data/D_rot_value.txt', 'w') as file:
    file.write('# D_rot (directly measured)\n')
    file.write('{}\n'.format(fit_params[0]/(2)))
