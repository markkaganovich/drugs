barcodes = ['TCGATG', 'GGTCTA', 'TGGATC', 'TGCTAG', 'GTCAGT'] 

file = open('./pool31')
lines = file.readlines()
lcl = []
for l in lines:
    lcl.append('NA' + l.strip('\n'))

notinpools = ['NA19225','NA19238', 'NA19239', 'NA19240', 'NA18856', 'NA18858', 'NA18562', 'NA18563', 'NA18853', 'NA18861']

file = open('./orderprobes1115')
lines = file.readlines()
newseqs = {}
names = []
newseqslist = []
for l in lines:
    cell = l.split('\t')[0]
    if cell in lcl and not cell in notinpools:
        lcl.remove(cell)
        seq = l.split('\t')[1]
        ext = seq[57:]
        for i in range(0,len(barcodes)):
            newseqs[cell+'bar' + str(i)] = seq[0:56] + barcodes[i] + ext
            names.append(cell +'bar' + str(i))
            newseqslist.append(seq[0:56] + barcodes[i] + ext)

file = open('pool31barcodedprobes1208','w')
for k in range(0,len(names)):
    file.write(names[k] + '\t' + newseqslist[k])
file.close()







