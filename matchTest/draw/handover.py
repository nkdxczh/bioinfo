import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

f = open('../result')
d = dict()

for line in f:
    ele = line.split(' ')
    family = ele[0]
    '''if ele[0][-2:] == '10':
        family = ele[0][:-2]
    else:
        family = ele[0][:-1]'''
    if family not in d:
        d[family] = []
    match = float(ele[3].split('/')[0])
    RMSD = float(ele[6])
    d[family].append([match, RMSD])

f.close()

families = []
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
    RMSD_means.append(i[2])

print families

speed_count = [i*10 for i in range(len(families))]

fig = plt.figure(figsize = (13, 7))
ax = fig.add_subplot(111)
ax.set_xlabel("Families", fontsize = 20)
ax.set_ylabel("Average Matched Atoms Count", fontsize = 20)
matplotlib.rc('xtick', labelsize=20)
matplotlib.rc('ytick', labelsize=20)

ax2 = ax.twinx()
ax2.set_ylabel('Average LRMSD (A)', fontsize = 20)

bar_width = 2

ax.grid()
#ax.set_xlim(0, 130)
#ax.set_ylim(0, 125)

matches = ax.bar(speed_count, match_means, bar_width, alpha=0.8,color= 'r',label = 'Matches',hatch = '//')
LRMSDs = ax.bar(np.array(speed_count) + bar_width +0.5, RMSD_means, bar_width, alpha=0.8,color= 'b',label = 'LRMSD')

ax.legend(handles=[matches,LRMSDs],loc = 'upper right',fontsize = 17)

ax.set_xticks(np.array(speed_count) + bar_width / 2)
ax.set_xticklabels(families, {'fontsize':20})
    
plt.tight_layout()
plt.savefig('handoverresult4.eps', format = 'eps', dpi = 10000)
plt.show()

