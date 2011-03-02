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

#def findbarcode(population, file, dnastretch = 300):
population = 'YRI'
pool = info.pool().YRItest
POS, HOMO, REF, ALT, FILTER, INFO = 1,2,3,4,6, 7
file = '../1000GenomesData/YRI.low_coverage.2010_09.genotypes.vcf'   
file = open(file)
outputfile = open('./'+population+'output', 'w')

barcode_snps = [] 
individuals = info.lineinfo().individuals[population]

lines = file.readlines(1000000)
for line in lines:
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
        if (num_genos == 2 and snp not in barcode_snps and not (l[cell_line_ID[0]] in info.lineinfo().trios[population]
                and l[cell_line_ID[1]] in info.lineinfo().trios[population]) and len([j for j in cell_line_ID if individuals[j] in pool]) == 1):
            barcode_snps.append([cell_line_ID, chr, hg18pos, hg18seq, varseq])

#[genomic_positions, regions] = get_rid_of_redundancy(regions, genomic_positions)

simplejson.dump(barcode_snps, outputfile)                
file.close()
outputfile.close()

#return barcode_snps



