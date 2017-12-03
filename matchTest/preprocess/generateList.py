from os import listdir
from os.path import isfile, join

mypath = '../pdb'
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
print onlyfiles

f = open("../blist", "w")
f.write('SIZE ' + str(len(onlyfiles)) + '\n')

for i in onlyfiles:
    f.write('pdb/' + i + '\n')

f.close()
