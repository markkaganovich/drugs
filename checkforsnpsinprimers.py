
import simplejson
'''
file = open('../1000GenomesData/partc')
lines = file.readlines()
file.close()

snpsites = {}
for i in range(1,22): snpsites[i] = []

for line in lines: 
    if line.startswith('#'):
        continue;
    pos = int(line.split('\t')[1])
    chr = int(line.split('\t')[0])
    if chr in snpsites.keys():
        snpsites[chr].append(pos)

file = open('snpsitesc','w')
simplejson.dump(snpsites, file)
file.close()

'''
file = open('snpsitesaa')
snpsites = simplejson.load(file)
file.close()
file = open('snpsitesb')
snpsites2 = simplejson.load(file)
file.close()
file = open('snpsitesc')
snpsites3 = simplejson.load(file)
file.close()
file = open('snpsitesd')
snpsites4 = simplejson.load(file)
file.close()
file = open('snpsitese')
snpsites5 = simplejson.load(file)
file.close()
file = open('snpsitesf')
snpsites6 = simplejson.load(file)
file.close()

i = 0
for snpdic in [snpsites2, snpsites3, snpsites4, snpsites5, snpsites6]:
    print i
    i = i+1
    for key in snpdic.keys():
        snpsites[key].extend(snpdic[key])








