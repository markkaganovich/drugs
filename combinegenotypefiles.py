

file1 = open('../CEU.low_coverage.2010_09.genotypes.vcf')
lines1 = file1.readlines()
file2 = open('../CHBJPT.low_coverage.2010_09.genotypes.vcf')
lines2 = file2.readlines()
outputfile = open('newgenotypefile','w')




for i in range(18,len(lines1)):
    a = lines1[i].split('\t')
    b = lines2[i].split('\t')
    if a[4] == b[4] and len(a[4]) == 1:
	newline = a+b[9:]
	for j in newline:
	    outputfile.write(j)
	    outputfile.write('\t')
	outputfile.write('\n')




