f = open('../motif')

motif = dict()
temd = dict()

for line in f:
    tem = line.split(' ')
    motif[tem[0] + ' ' + tem[1] + " " + tem[2]] = int(tem[-1])
    temd[tem[1]] = int(tem[-1])

f.close()

print motif


f = open("../anti-CRISPOR.pdb")

output = open("../anti-CRISPOR.FAT", 'w')
count = 0
index = []

output.write('NUMSOURCES= 1\nNUMARRAYS= 20\nAMINOACIDS= true\n\nSOURCEPATTERN 0 anti-CRISPOR.pdb 710 87\n');

for line in f:
    ele = line.split()
    if ele[0] != 'ATOM':
        continue
    label = ele[3] + ' ' + ele[5] + ' ' + ele[2]
    if label in motif:
        print line
        index.append(count)
    count += 1
    temp = ele[9].split('.')
    #s = 'SOURCEATOM ' + ele[6] + ' ' + ele[7] + ' ' + ele[8] + ' ' + ele[5] + ' ' + temp[1][2:] + ' ' + ele[2] + ' ' + ele[3] + ' ' + ele[4] + ele[10] + ' ' + str(count) + '\n'
    if ele[5] in temd:
        s = 'SOURCEATOM ' + ele[6] + ' ' + ele[7] + ' ' + ele[8] + ' ' + ele[5] + ' ' + str(temd[ele[5]]) + ' ' + ele[2] + ' ' + ele[3] + ' ' + 'CA' + ' ' + str(count) + '\n'
    else:
        s = 'SOURCEATOM ' + ele[6] + ' ' + ele[7] + ' ' + ele[8] + ' ' + ele[5] + ' ' + '1' + ' ' + ele[2] + ' ' + ele[3] + ' ' + 'CA' + ' ' + str(count) + '\n'
    output.write(s)
    output.write('TEMPORARY FILE LINE\n')

output.write('MOTIF: ')
for i in index:
    output.write(str(i) + ' ')

output.write('\n')

output.write('ENDSOURCEPATTERN 0\nTARGETPATTERN foo.pdb 1\nTARGETATOM -0.569000 23.355000 6.050000 1 -1 CA ALA -1\nATOM      1  CA  ALA     1      -0.569  23.355   6.050  1.00113.78      1AGI 106\nMOTIF: 0 \nENDTARGETPATTERN -1');

output.close()
f.close()

