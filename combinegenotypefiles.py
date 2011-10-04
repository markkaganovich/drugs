

file1 = open('../CEU.low_coverage.2010_09.genotypes.vcf')
lines1 = file1.readlines()
file2 = open('../CHBJPT.low_coverage.2010_09.genotypes.vcf')
lines2 = file2.readlines()
outputfile = open('newgenotypefile','w')


chroms1 = map(lambda x: x.split('\t')[0], lines1[18:])
chroms2 = map(lambda x: x.split('\t')[0], lines2[18:])

pos1 = map(lambda x: x.split('\t')[1], lines1[18:])
pos2 = map(lambda x: x.split('\t')[1], lines2[18:])

for i in range(0, len(pos1)):
    if (pos1[i] in pos2 
	and int(chroms1[i]) == int(chroms2[pos2.index(pos1[i])])):
	print "same SNP"
	newline = lines1[i+18]+lines2[pos2.index(pos1[i])+18]
	#for j in newline:
        outputfile.write(newline)
	    #outputfile.write('\t')
	outputfile.write('\n')




