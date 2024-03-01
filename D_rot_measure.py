import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def f(x, m, c):
    return m*x + c


n_trials = 10

# Read trial 1 to know what the time array is
t, ang_dev_rad, ang_dev_deg = np.loadtxt(
    'runs/run_{}/data/ang_deviation.txt'.format(1), unpack=True)

avg_ang_dev_rad = np.zeros(len(t))

for n in range(1, n_trials+1):
    # Load data
    print('Loading angular deviation data from run_{}'.format(n))
    t, ang_dev_rad, ang_dev_deg = np.loadtxt(
        'runs/run_{}/data/ang_deviation.txt'.format(n), unpack=True)

    for t_iter in range(len(t)):
        avg_ang_dev_rad[t_iter] += (ang_dev_rad[t_iter])**2


avg_ang_dev_rad = avg_ang_dev_rad/n_trials

fit_params, fit_cov = curve_fit(f, t, avg_ang_dev_rad, method='lm')

fitted_curve = np.zeros(len(t))
for i in range(len(t)):
    fitted_curve[i] = f(t[i], fit_params[0], fit_params[1])

print('Fitted parameters: m = {:.4e}, c = {:.4e}'.format(fit_params[0], fit_params[1]))
print('Diffusion coefficient: D = {:.4e}'.format(fit_params[0]/(2)))

plt.figure(tight_layout=True)

plt.plot(t, avg_ang_dev_rad, 'k-', label='Data')
plt.plot(t, fitted_curve, 'r--', label='Fitted curve')
plt.xlabel(r'$t/\tau$', fontsize=18)
plt.ylabel(r'$\langle \theta^2 \rangle$', fontsize=18)
plt.legend()

plt.savefig('runs_measurements/avg_ang_dev_sq.pdf')

with open('runs_measurements/D_rot_measurements.txt', 'w') as f:
    f.write("# Slope \t D_rot(measured, unscaled)\n")
    f.write("{:.4e} \t {:.4e}\n".format(fit_params[0], fit_params[0]/(2)))

avg_e2e = np.zeros(len(t))
mean_e2e = 0

for n in range(1, n_trials+1):
    # Load data
    print('Loading e2e data from run_{}'.format(n))
    t, e2e = np.loadtxt(
        'runs/run_{}/data/e2e_dist.txt'.format(n), unpack=True)

    for t_iter in range(len(t)):
        avg_e2e[t_iter] += (e2e[t_iter])
    
    mean_e2e += np.mean(e2e)

avg_e2e = avg_e2e/n_trials
mean_e2e = mean_e2e/n_trials

mean_line = mean_e2e*np.ones(len(t))

plt.figure(tight_layout=True)

plt.plot(t, avg_e2e, 'k-', label='Averaged length', linewidth=1.0)
plt.plot(t, mean_line, 'r--', label='Mean length ({:.1f})'.format(mean_e2e), linewidth=1.0)

plt.xlabel(r'$t/\tau$', fontsize=18)
plt.ylabel(r'Averaged end to end distance', fontsize=18)

plt.legend()

plt.savefig('runs_measurements/avg_e2e_dist.pdf')

