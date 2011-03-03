# ----------------------------------------------------------------------------------
#
# Mark Kaganovich,  Jan, 2011
#
# find SNPs that will identify 1000genome lines uniquely based on
# the recently released SNP calls
#
# -----------------------------------------------------------------------------------

import SNPhelpers
import simplejson 
import logging
import info
import os

def findbarcode(population, file = '../1000GenomesData/YRI.low_coverage.2010_09.genotypes.vcf'):   
	pool = info.pool().YRItest
	POS, HOMO, REF, ALT, FILTER, INFO = 1,2,3,4,6, 7
	file = open(file)

	barcode_snps1 = []
	barcode_snps2 = [] 
        t = info.lineinfo()
	individuals = t.individuals[population]
	celllines = t.YRIfilecells

	lines = file.readlines()
	i=0
	for line in lines:
	    i=i+1
	    print i
	    if line.startswith('#'):
		continue
	    l = line.split('\t')
	    a=l[INFO].split(';')
	    AC = filter(lambda x: 'AC' in x, a)
	    AN = filter(lambda x: 'AN' in x, a)
	    if ('PASS' not in l[FILTER] or len(l[ALT]) != 1 or 'MP'  in l[INFO]): 
		continue
	    genotype = SNPhelpers.vcf_to_geno_struct([line])
	    if genotype:
		[num_genos, cell_line_ID] = SNPhelpers.num_genotypes(genotype)
		hg18pos = int(l[POS])
		hg18seq = str(l[REF])
		varseq = str(l[ALT])
		chr = int(l[0])
		snp = [cell_line_ID, chr, hg18pos, hg18seq, varseq]
		if (num_genos == 1 and snp not in barcode_snps1 and len([j for j in cell_line_ID if individuals[j] in pool]) == 1):
		    barcode_snps1.append([cell_line_ID, chr, hg18pos, hg18seq, varseq])
		if (num_genos == 2 and snp not in barcode_snps2 and len([j for j in cell_line_ID if individuals[j] in pool]) == 1):
		    barcode_snps2.append([cell_line_ID, chr, hg18pos, hg18seq, varseq])
	file.close()

	return [barcode_snps1, barcode_snps2]

if __name__ == "__main__":
	
	dir = '../1000GenomesData/YRI'
	barcode1 = []
	barcode2 = []
	for file in os.listdir(dir):
	    b = findbarcode('YRI', dir+'/'+file)
	    barcode1.extend(b[0])
	    barcode2.extend(b[1])
	
	outputfile = open('./YRIoutputUniqueSNPs', 'w')
	simplejson.dump([barcode1, barcode2], outputfile)
	outputfile.close()





