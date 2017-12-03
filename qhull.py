
f = open('anti-CRISPOR.pdb')
s = f.readlines()

res = []
for i in s:
    contents = i.split()
    if len(contents) == 10 and contents[0] == 'ATOM':
        res.append([float(contents[5]), float(contents[6]), float(contents[7])])

print '3'
print len(res)
for i in res:
    print i[0], i[1], i[2]
