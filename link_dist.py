import numpy as np
import matplotlib.pyplot as plt
import cylindermath as cm


rA = np.array([0.0, 175.0, 175.0])
rB = np.array([1000.0, 175.0, 175.0])
c = cm.cylinder(175.0, rA, rB)

data = np.loadtxt("link_pos.txt", unpack=True)

window = 21
print("Smoothing window =", window)

detection_window = 21
print("Detection window =", detection_window)

threshold = 0.25
print("Threshold =", threshold)

hitting_distance = 2.5

# print(len(data))

info = np.loadtxt('info.txt')

num_monomers = int(info[0])
num_linkers = int(info[1])

# print("Number of linkers =", num_linkers)

t = data[0]

avg_dist = np.zeros((len(t),))

plt.figure(tight_layout=True)

smooth_all_linkers_dist = np.zeros((num_linkers, len(t)))
hit_detection_all_linkers = np.zeros((num_linkers, len(t)))

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
    plt.plot(t, dist_list_s, label="Linker {}".format(i), linewidth=0.5)
    # plt.plot(t, dist_list, label="Linker {}".format(i), linewidth=0.5)
    avg_dist += dist_list

    for t_step in range(len(t)):
        smooth_all_linkers_dist[i-1][t_step] = dist_list_s[t_step]
        if (abs(hitting_distance - smooth_all_linkers_dist[i-1][t_step])/hitting_distance < threshold):
            hit_detection_all_linkers[i-1][t_step] = 1
        else:
            hit_detection_all_linkers[i-1][t_step] = 0


avg_dist *= (1/num_linkers)

print("Minimum avg point =", min(avg_dist))

avg_dist_s = cm.moving_average(avg_dist, window, padding="edge")

print("Minimum avg(s) point =", min(avg_dist_s))

plt.plot(t, avg_dist_s, 'k--', label="Average", linewidth=1.0)
# plt.plot(t, avg_dist, 'k--', label="Average", linewidth=1.0)

plt.xlabel(r'$t/\tau$', fontsize=18)
plt.ylabel("Distance from cell membrane (nm)", fontsize=16)
plt.ylim(bottom=0)
plt.tight_layout()

plt.legend()
plt.savefig("link_dist.pdf", bbox_inches='tight')

# plt.show()

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
plt.scatter(t, hit_detection_avg, label="Average", marker=".", color="black", s=0.3)

plt.legend(loc="best", markerscale=10)
plt.savefig("hit_detection.pdf", bbox_inches='tight')
