import numpy as np
import matplotlib.pyplot as plt

data = open('test.interactions.dump', 'r')
lines = data.readlines()

attachment_d = 2.7
print("Attachment distance: {}".format(attachment_d))

t_list_original = np.loadtxt('data/com_pos.txt', usecols=0)

t_list = []
d_list = []

t_attach_list = []
num_attach_list = []

timestep = 0
attachments_found = 0

for line in lines:
    if 'ITEM: TIMESTEP' in line:
        timestep += 1

        t_attach_list.append(t_list_original[timestep - 1])
        num_attach_list.append(attachments_found)

        # print(f'Processing timestep {timestep}')
        attachments_found = 0
    else:
        line = line.split()
        # t = t_list_original[timestep - 1]
        # d = line[3]
        d = round(float(line[3]), 2)

        # t_list.append(t)
        # d_list.append(d)

        if (d < 3.0):
            attachments_found += 1

# for i in range(len(t_list)):
#     print("{} \t {}".format(t_list[i], d_list[i]))

# plt.scatter(t_list, d_list)
# plt.xlabel(r'$t/\tau$', fontsize=18)
# plt.ylabel(r'Distance of interacting pairs (nm)', fontsize=14)
# plt.show()

# plt.clf()
# plt.cla()

ax = plt.figure(tight_layout=True).gca()
ax.yaxis.get_major_locator().set_params(integer=True)
# plt.figure(tight_layout=True)

plt.plot(t_attach_list, num_attach_list, color='k', linewidth=1.0)
plt.xlabel(r'$t/\tau$', fontsize=18)
plt.ylabel(r'Number of attached monomers', fontsize=14)
plt.grid(axis='y', linestyle='--', linewidth=0.8)
plt.axhline(y=2, color='r', linestyle='--', linewidth=1.0)
plt.savefig('plots/mono_attachment.pdf')
