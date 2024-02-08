"""
This module provides functionality for analyzing the distance of linkers from a cell membrane.

The module uses the cylindermath module to calculate distances from a cylinder representing the cell membrane. It reads linker positions from a text file, calculates their distances from the cell membrane, and applies a moving average smoothing to the distance data.

The module also includes functionality for detecting when linkers hit the cell membrane, calculating the degree of attachment of the linkers to the cell membrane, and counting the number of linkers attached to the cell membrane.

The module provides options for plotting the distance data, hit detection data, degree of attachment data, and number of attached linkers data.

The module writes the smoothed distance data, hit detection data, degree of attachment data, and number of attached linkers data to text files.

This module requires the numpy and matplotlib libraries.

Example:
    To run the module, ensure that the required data files are in the 'data' directory and the required output directories ('plots' and 'data') exist. Then, simply run the module with Python:

    ```python
    python3 link_dist.py
    ```

    The module will print information about the linkers and their attachment to the cell membrane to the console, and will create plots and data files in the 'plots' and 'data' directories, respectively.
"""

import numpy as np
import matplotlib.pyplot as plt
import cylindermath as cm


rA = np.array([0.0, 175.0, 175.0])
rB = np.array([1000.0, 175.0, 175.0])
c = cm.cylinder(175.0, rA, rB)

data = np.loadtxt("data/link_pos.txt", unpack=True)

dist_list_str = "data/dist_list.txt"
hit_detection_str = "data/hit_detection.txt"
degree_str = "data/degree.txt"
two_attached_str = "data/two_attached.txt"
num_attached_str = "data/num_attached.txt"

two_first_attached_time_str = "data/two_first_attached_time.txt"

# When this value of attachment is reached, filament will be considered attached
attachment_threshold = 0.2

# Plotting toggles
toggle_plot_traces = 1
toggle_plot_hit_detection = 1
toggle_plot_degree = 1
toggle_plot_two_attached = 1
toggle_plot_num_attached = 1

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

print("Number of linkers = {}".format(num_linkers))

t = data[0]

avg_dist = np.zeros((len(t),))

plt.figure(tight_layout=True)

smooth_all_linkers_dist = np.zeros((num_linkers, len(t)))
hit_detection_all_linkers = np.zeros((num_linkers, len(t)))

full_dist_list = np.zeros((len(t), num_linkers))

for i in range(1, num_linkers+1):
    dist_list = []

    x = data[3*(i-1)+1]
    y = data[3*(i-1)+2]
    z = data[3*(i-1)+3]

    for t_step in range(len(t)):
        rP = np.array([x[t_step], y[t_step], z[t_step]])
        dist = cm.distance_from_surface(c, rP)
        dist_list.append(dist)

    dist_list = np.array(dist_list)
    dist_list_s = cm.moving_average(dist_list, window, padding="edge")

    for d in range(len(dist_list)):
        full_dist_list[d][i-1] = dist_list_s[d]

    if (toggle_plot_traces == 1):
        plt.plot(t, dist_list_s, label="Linker {}".format(i), linewidth=0.5)
    avg_dist += dist_list

    for t_step in range(len(t)):
        smooth_all_linkers_dist[i-1][t_step] = dist_list_s[t_step]
        if (abs(hitting_distance - smooth_all_linkers_dist[i-1][t_step])/hitting_distance < threshold):
            hit_detection_all_linkers[i-1][t_step] = 1
        else:
            hit_detection_all_linkers[i-1][t_step] = 0

# --- Writing output files---

with open(dist_list_str, 'w') as f:
    for i in range(len(full_dist_list)):
        f.write('{}\t'.format(t[i]))
        for j in range(len(full_dist_list[i])):
            f.write('{:.4f}\t'.format(full_dist_list[i][j]))
        f.write('\n')

hit_detection_all_linkers_t = np.array(hit_detection_all_linkers).transpose()
with open(hit_detection_str, 'w') as f:
    for i in range(len(hit_detection_all_linkers_t)):
        f.write('{}\t'.format(t[i]))
        for j in range(len(hit_detection_all_linkers_t[i])):
            f.write('{}\t'.format(hit_detection_all_linkers_t[i][j]))
        f.write('\n')

# ---Degree of attachment---

is_attached = 0
attached_at_t = 0
detached_at_t = max(t)
one_linker_first_attached = 0
one_linker_first_attached_at_t = max(t)
two_linkers_first_attached = 0
two_linkers_first_attached_at_t = max(t)
num_linkers_attached = np.zeros((len(t)))

two_attached = np.zeros((len(t)))

attachment_line = attachment_threshold * np.ones((len(t)))

comparison_one = ""
comparison_two = ""
degree_of_attachment = np.zeros((len(t)))

for t_step in range(len(t)):
    monomers_currently_attached = []
    for i in range(num_linkers):
        if (hit_detection_all_linkers_t[t_step][i] == 1):
            j = i * num_skip
            j = num_monomers - j
            monomers_currently_attached.append(j)

    num_linkers_attached[t_step] = len(monomers_currently_attached)

    current_degree = 0
    if (len(monomers_currently_attached) > 0):
        start = min(monomers_currently_attached)
        stop = max(monomers_currently_attached)
        current_degree = (stop - start + 1) / num_monomers
        degree_of_attachment[t_step] = current_degree

    if (is_attached == 0):
        if (current_degree >= attachment_threshold):
            attached_at_t = t[t_step]
            is_attached = 1
    else:
        if (current_degree < attachment_threshold):
            detached_at_t = t[t_step]

    if (len(monomers_currently_attached) >= 2):
        two_attached[t_step] = 1
        if (two_linkers_first_attached == 0):
            two_linkers_first_attached = 1
            two_linkers_first_attached_at_t = t[t_step]
        if (two_linkers_first_attached_at_t == max(t)):
            comparison_two = ">"
        else:
            comparison_two = "="

    if (len(monomers_currently_attached) >= 1):
        if (one_linker_first_attached == 0):
            one_linker_first_attached = 1
            one_linker_first_attached_at_t = t[t_step]
        if (one_linker_first_attached_at_t == max(t)):
            comparison_one = ">"
        else:
            comparison_one = "="

if (comparison_one == ">" or one_linker_first_attached == 0):
    print("One linker never attached")
else:
    print("One linker first attached at t = {}".format(
        one_linker_first_attached_at_t))

if (comparison_two == ">" or two_linkers_first_attached == 0):
    print("Two linkers never attached")
else:
    print("Two linkers first attached at t = {}".format(
        two_linkers_first_attached_at_t))
with open(two_first_attached_time_str, 'w') as f:
    f.write("# First line: one linker first attached at t, second line: two linkers first attached at t\n")
    f.write('{}\n{}'.format(one_linker_first_attached_at_t, two_linkers_first_attached_at_t))

print("-"*20)

print("t = {}\nattached monomers = {}\ndegree of attachment = {:.4f}".format(
    max(t), monomers_currently_attached, degree_of_attachment[t_step]))

print("Number of linkers attached = {}".format(
    len(monomers_currently_attached)))

if (detached_at_t == max(t)):
    comparison = ">"
else:
    comparison = "="

print("First attached at t = {}".format(attached_at_t, comparison))

if (comparison == "="):
    print("Last detached at t = {}".format(detached_at_t))
else:
    print("Filament never detached")

with open(degree_str, 'w') as f:
    for i in range(len(degree_of_attachment)):
        f.write('{}\t{}\n'.format(t[i], degree_of_attachment[i]))

with open(two_attached_str, 'w') as f:
    for i in range(len(t)):
        f.write('{}\t{}\n'.format(t[i], two_attached[i]))

with open(num_attached_str, 'w') as f:
    for i in range(len(t)):
        f.write('{}\t{}\n'.format(t[i], num_linkers_attached[i]))

# ----Average Distance from Cell Membrane---
avg_dist *= (1/num_linkers)
avg_dist_s = cm.moving_average(avg_dist, window, padding="edge")

# print("Minimum avg point = {:.4f}".format(min(avg_dist)))
# print("Minimum avg(s) point = {:.4f}".format(min(avg_dist_s)))

if (toggle_plot_traces == 1):
    plt.plot(t, avg_dist_s, 'k--', label="Average", linewidth=1.0)
    # plt.plot(t, avg_dist, 'k--', label="Average", linewidth=1.0)

plt.xlabel(r'$t/\tau$', fontsize=18)
plt.ylabel("Distance from cell membrane (nm)", fontsize=16)
plt.ylim(bottom=0)
plt.tight_layout()

if (toggle_plot_traces == 1):
    plt.legend(fancybox=True)
    plt.savefig("plots/link_dist.pdf", bbox_inches='tight')

# plt.show()

# ---Hit Detection plot---

plt.clf()

offset_scale = 0.2 / num_linkers
offset = offset_scale * np.ones(len(t))

plt.xlabel(r'$t/\tau$', fontsize=18)
plt.ylabel("Hit detection", fontsize=18)
plt.ylim(bottom=-0.1, top=offset_scale * num_linkers + 2)
plt.tight_layout()

for i in range(1, num_linkers+1):
    hit_detection_all_linkers[i -
                              1] = hit_detection_all_linkers[i-1] + (i-1) * offset
    if (toggle_plot_hit_detection == 1):
        plt.scatter(t, hit_detection_all_linkers[i-1],
                    label="Linker {}".format(i), marker=".", s=0.25)

hit_detection_avg = []
smooth_detection = cm.moving_average(
    avg_dist, detection_window, padding="edge")

for i in range(len(smooth_detection)):
    if abs(hitting_distance - smooth_detection[i])/hitting_distance < threshold:
        hit_detection_avg.append(1)
    else:
        hit_detection_avg.append(0)

hit_detection_avg = hit_detection_avg + ((num_linkers)) * offset
if (toggle_plot_hit_detection == 1):
    plt.scatter(t, hit_detection_avg, label="Average",
                marker=".", color="black", s=0.3)

if (toggle_plot_hit_detection == 1):
    plt.legend(loc="best", markerscale=10)
    plt.savefig("plots/hit_detection.pdf", bbox_inches='tight')

# --- Degree of attachment plot ---

if (toggle_plot_degree == 1):
    plt.clf()
    plt.cla()
    plt.xlabel(r'$t/\tau$', fontsize=18)
    plt.ylabel("Degree of attachment", fontsize=18)
    plt.ylim(bottom=-0.1, top=1.1)
    plt.tight_layout()
    plt.plot(t, degree_of_attachment, 'k', linewidth=0.6)
    # plt.grid(axis='y', linestyle='--', linewidth=0.5, which='both')
    # plt.plot(t, attachment_line, 'r--', linewidth=0.5)
    plt.savefig("plots/degree_of_attachment.pdf", bbox_inches='tight')
    # plt.show()

# --- Two linkers attached plot ---

if (toggle_plot_two_attached == 1):
    plt.clf()
    plt.cla()
    plt.xlabel(r'$t/\tau$', fontsize=18)
    plt.ylabel("Two monomers attached", fontsize=18)
    plt.ylim(bottom=-0.1, top=1.1)
    plt.tight_layout()
    plt.plot(t, two_attached, 'k', linewidth=0.6)
    plt.grid(axis='y', linestyle='--', linewidth=0.5, which='major')
    plt.savefig("plots/two_attached.pdf", bbox_inches='tight')

# --- Number of linkers attached plot ---

line_at_two = 2 * np.ones(len(t))
if (toggle_plot_num_attached == 1):
    plt.clf()
    plt.cla()
    plt.xlabel(r'$t/\tau$', fontsize=18)
    plt.ylabel("Number of linkers attached", fontsize=18)
    plt.ylim(bottom=-0.1, top=num_linkers + 1)
    plt.tight_layout()
    plt.grid(axis='y', linestyle='--', linewidth=0.5)
    plt.plot(t, line_at_two, 'r--', linewidth=0.7)
    plt.plot(t, num_linkers_attached, 'k', linewidth=0.6)
    plt.savefig("plots/num_attached.pdf", bbox_inches='tight')
