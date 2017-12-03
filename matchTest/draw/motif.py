import numpy as np
import matplotlib.pyplot as plt

f = open('../result')
d = dict()

for line in f:
    ele = line.split(' ')
    '''family = ele[0]
    if ele[0][-2:] == '10':
        family = ele[0][:-2]
    else:
        family = ele[0][:-1]
    if family not in d:
        d[family] = []
    match = float(ele[3].split('/')[0])
    RMSD = float(ele[6])
    d[family].append([match, RMSD])'''
    for i in range(9, len(ele) - 1, 2):
        motif = ele[i][-4:-2]
        if motif == '(8':
            motif = '8'
        if motif not in d:
            d[motif] = 0
        d[motif] += 1

f.close()

data = []
for i in d:
    data.append([i, d[i]])

data = sorted(data, key = lambda x: -x[1])
print data


'''families = []
match_means = []
RMSD_means = []

data = []

for i in d:
    families.append(i)
    sum1 = 0
    sum2 = 0
    for j in d[i]:
        sum1 += j[0]
        sum2 += j[1]
    match_means.append(sum1 / len(d[i]))
    RMSD_means.append(sum2 / float(len(d[i])))
    data.append([i, sum1 / len(d[i]), sum2 / len(d[i])])

data = sorted(data, key = lambda x: (-x[1], x[2]))
print data
families = []
match_means = []
RMSD_means = []

for i in data:
    families.append(i[0])
    match_means.append(i[1])
    RMSD_means.append(i[2])'''

N = len(data)

motifs = []
means = []
for i in data:
    motifs.append(i[0])
    means.append(i[1])

ind = np.arange(N)  # the x locations for the groups
width = 0.35       # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind, means, width, color='r')

# add some text for labels, title and axes ticks
ax.set_title('Matching Results for different Superfamilies')
ax.set_xticks(ind)
ax.set_xticklabels(motifs)

def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%d' % int(height),
                ha='center', va='bottom')
        
#autolabel(rects1)
#autolabel(rects2)

plt.show()
