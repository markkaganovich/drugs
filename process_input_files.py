#----------------------------------------------------------
#
# split SNP genotype input files by chrom
#
#
#
#-----------------------------------------------------------

file = open('./CEU.genotypes.vcf')

chrom = 1
outputfile = open('./inputCEU/SNPs_chrom1', 'w')
while(1):
    line = file.readline()
    if not line:
        break
    tokens = line.split('\t')
    print int(tokens[0])
    if chrom == int(tokens[0]):
        outputfile.write(line)
    else:
        chrom = int(tokens[0])
        outputfile.close()
        outputfile = open('./inputCEU/SNPs_chrom'+str(chrom),'w')



