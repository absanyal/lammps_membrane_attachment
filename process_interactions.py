import numpy as np
import matplotlib.pyplot as plt

data = open('test.interactions.dump', 'r')
lines = data.readlines()

attachment_d_list = [2.5 * pow(2, 1/6), 3.10, 3.57, 4.06]

color_list = ['red', 'blue', 'green', 'magenta', 'purple', 'brown', 'pink', 'gray']

print("Attachment distances:")
for attachment_d in attachment_d_list:
    print("{:.2f}   ".format(attachment_d), end=' ')
print()

attachment_minimum = 2
attachment_minimum = int(attachment_minimum)
print("Minimum monomers to be attatched: {}".format(attachment_minimum))

t_list_original = np.loadtxt('data/com_pos.txt', usecols=0)

plt.figure(tight_layout=True)

for i, attachment_d in enumerate(attachment_d_list):

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

            if (d < attachment_d):
                attachments_found += 1


    max_attachments = max(num_attach_list)
    max_attachments = max_attachments + 1

    yticklist = np.arange(0, max_attachments, 1)

    plt.plot(t_attach_list, num_attach_list, linewidth=1.0, label='{:.2f}'.format(attachment_d), color=color_list[i], alpha=0.5)

plt.xlabel(r'$t/\tau$', fontsize=18)
plt.ylabel(r'Number of attached monomers', fontsize=14)
plt.yticks(yticklist)
plt.grid(axis='y', linestyle='--', linewidth=0.8)
plt.axhline(y=attachment_minimum, color='r', linestyle='--', linewidth=1.0)
# plt.title("Minimum distance for attachment: {:.2f} nm".format(attachment_d))
plt.legend(title=r'$d_{{attach}}$ (nm)', fontsize=12)
plt.savefig('plots/mono_attachment.pdf')
