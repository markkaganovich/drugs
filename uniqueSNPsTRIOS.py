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

def findbarcode(file = '../1000GenomesData/YRI.low_coverage.2010_09.genotypes.vcf'):   
#pool = info.pool().YRItest
    POS, HOMO, REF, ALT, FILTER, INFO = 1,2,3,4,6, 7
    file = open(file)
    barcodesnps = []
    lines = file.readlines(10000000)
    while lines != []:
        print file.tell()
        for i,line in enumerate(lines):
            if line.startswith('#') or (i > 0 and lines[i-1].startswith('#')):
                continue
            l = line.split('\t')
            a=l[INFO].split(';')
            AC = filter(lambda x: 'AC' in x, a)
            AN = filter(lambda x: 'AN' in x, a)
            if ('PASS' not in l[FILTER] or len(l[ALT]) != 1 or 'MP'  in l[INFO]): 
                continue
            l0 = lines[i-1].split('\t')
            if not(int(l[POS]) - int(l0[POS]) > 0 and int(l[POS]) -int(l0[POS]) <= 50):
                pass;
            genotype = SNPhelpers.vcf_to_geno_struct([line , lines[i-1]])
            if genotype:
                cell_line_ID = SNPhelpers.num_genotypes_trios(genotype)
                hg18pos = int(l[POS])
                hg18seq = str(l[REF])
                varseq = str(l[ALT])
                chr = int(l[0])
                snp = [cell_line_ID, chr, hg18pos, hg18seq, varseq]
            if snp not in barcodesnps:
                barcodesnps.append(snp)
        lines = file.readlines(1000000000)
    file.close()

    return barcodesnps

if __name__ == "__main__":
	
    file = '../1000GenomesData/YRI.low_coverage.2010_09.genotypes.vcf'
    b = findbarcode(file)
    outputfile = open(file + '.outputUniqueSNPs', 'w')
    simplejson.dump(b, outputfile)
    outputfile.close()





