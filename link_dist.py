import numpy as np
import matplotlib.pyplot as plt
import cylindermath as cm


rA = np.array([0.0, 175.0, 175.0])
rB = np.array([1000.0, 175.0, 175.0])
c = cm.cylinder(175.0, rA, rB)

data = np.loadtxt("link_pos.txt", unpack=True)

dist_list_str = "dist_list.txt"
hit_detection_str = "hit_detection.txt"
degree_str = "degree.txt"

# When this value of attachment is reached, filament will be considered attached
attachment_threshold = 0.2

# Plotting toggles
plot_traces = 1
plot_hit_detection = 1
plot_degree = 1

window = 11
print("Smoothing window ={}".format(window))

detection_window = 21
print("Detection window = {}".format(detection_window))

threshold = 0.25
print("Threshold = {:.2f}".format(threshold))

hitting_distance = 2.5
print("Hitting distance = {:.2f}".format(hitting_distance))

# print(len(data))

info = np.loadtxt('info.txt')

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
    
    if (plot_traces == 1):
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

attachment_line = attachment_threshold * np.ones((len(t)))

degree_of_attachment = np.zeros((len(t)))
for t_step in range(len(t)):
    monomers_currently_attached = []
    for i in range(num_linkers):
        if (hit_detection_all_linkers_t[t_step][i] == 1):
            j = i * num_skip
            j = num_monomers - j
            monomers_currently_attached.append(j)
    
    
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
        
    
print("t = {}\nattached monomers = {}\ndegree of attachment = {:.4f}".format(max(t), monomers_currently_attached, degree_of_attachment[t_step]))

if (detached_at_t == max(t)):
    comparison = ">"
else:
    comparison = "="

print("First attached at t = {}\nLast detached at t {} {}".format(attached_at_t, comparison, detached_at_t))

with open(degree_str, 'w') as f:
    for i in range(len(degree_of_attachment)):
        f.write('{}\t{}\n'.format(t[i], degree_of_attachment[i]))
    

# ----Average Distance from Cell Membrane---
avg_dist *= (1/num_linkers)
avg_dist_s = cm.moving_average(avg_dist, window, padding="edge")

# print("Minimum avg point = {:.4f}".format(min(avg_dist)))
# print("Minimum avg(s) point = {:.4f}".format(min(avg_dist_s)))

if (plot_traces == 1):
    plt.plot(t, avg_dist_s, 'k--', label="Average", linewidth=1.0)
    # plt.plot(t, avg_dist, 'k--', label="Average", linewidth=1.0)

plt.xlabel(r'$t/\tau$', fontsize=18)
plt.ylabel("Distance from cell membrane (nm)", fontsize=16)
plt.ylim(bottom=0)
plt.tight_layout()

if (plot_traces == 1):
    plt.legend()
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
    if (plot_hit_detection == 1):
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
if (plot_hit_detection == 1):
    plt.scatter(t, hit_detection_avg, label="Average",
                marker=".", color="black", s=0.3)

if (plot_hit_detection == 1):
    plt.legend(loc="best", markerscale=10)
    plt.savefig("plots/hit_detection.pdf", bbox_inches='tight')
    
# --- Degree of attachment plot ---

plt.clf()

plt.xlabel(r'$t/\tau$', fontsize=18)
plt.ylabel("Degree of attachment", fontsize=18)
plt.ylim(bottom=0, top=1.1)
plt.tight_layout()

if (plot_degree == 1):
    plt.plot(t, degree_of_attachment, 'k', linewidth=0.6)
    plt.plot(t, attachment_line, 'r--', linewidth=0.5)
    plt.savefig("plots/degree_of_attachment.pdf", bbox_inches='tight')
    # plt.show()



