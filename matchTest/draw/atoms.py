import pandas as pd
import numpy as np
import matplotlib
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

motifs = []
means = []
for i in data:
    motifs.append(i[0])
    means.append(i[1])

speed_count = [i*10 for i in range(len(motifs))]

fig = plt.figure(figsize = (13, 7))
ax = fig.add_subplot(111)
ax.set_xlabel("Atoms Number", fontsize = 20)
ax.set_ylabel("Matches Count", fontsize = 20)
matplotlib.rc('xtick', labelsize=20)
matplotlib.rc('ytick', labelsize=20)

bar_width = 4

ax.grid()
#ax.set_xlim(0, 130)
#ax.set_ylim(0, 125)

matches = ax.bar(speed_count, means, bar_width, alpha=0.8,color= 'r',label = 'Matches',hatch = '//')

#ax.legend(handles=[matches],loc = 'upper right',fontsize = 17)

ax.set_xticks(np.array(speed_count))
ax.set_xticklabels(motifs, {'fontsize':20})
    
plt.tight_layout()
plt.savefig('handoverresult4.eps', format = 'eps', dpi = 10000)
plt.show()

