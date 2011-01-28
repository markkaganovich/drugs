# ----------------------------------------------------------------------------------
#
# Mark Kaganovich,  Jan, 2011
#
# find SNPs that will identify 1000genome lines uniquely based on
# the recently released SNP calls
#
# -----------------------------------------------------------------------------------

from helperfuns_findbarcode import *
import simplejson
import logging

#def findbarcode(population, file, dnastretch = 300):
population = 'CEU'
DATABYTES = 10000000
POS = 1
REF = 3
ALT = 4
FILTER = 6
INFO = 7
HOMO = 2
file = '../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf'   
file = open(file)
outputfile = open('./'+population+'output', 'w')

# read through header
lines = file.readlines(DATABYTES)
for z in range(0, len(lines)):
    print z
    if lines[z].startswith('#CHROM'):
        print lines[z]
        NUM_GENOMES = len(list(lines[z].split('\t')[9:]))
        print NUM_GENOMES         
        lines = lines[z+1:]
        break
    else:
        pass

barcode_snps=[]
regions=[]
genomic_positions=[]
j = 0

while(1):
    if not lines:
        break
#print j
    
    if len(lines) < j+1:
        del lines
        lines = file.readlines(DATABYTES)
        j=0
        if len(lines) <= j+1:                                # if reading more lines does not extend LINES, done
            break
    
    l = list(lines[j].split('\t'))
    a=l[INFO].split(';')
    AC = filter(lambda x: 'AC' in x, a)
    AN = filter(lambda x: 'AN' in x, a)
    if ('PASS' not in l[FILTER] or len(l[ALT]) != 1 or 'MP'  in l[INFO] or
        int(AC[0][AC[0].index('=')+1:]) != HOMO or int(AN[0][AN[0].index('=')+1:]) != NUM_GENOMES*2):
        j=j+1
        continue

    genotype = vcf_to_geno_struct([lines[j]])
    if genotype:
        [num_genos, cell_line_ID] =num_genotypes(genotype)
        
        hg18pos = int(l[POS])
        hg18seq = str(l[REF])
        varseq = str(l[ALT])
        chr = int(l[0])
            
        if num_genos == 1:
            barcode_snps.append([j, cell_line_ID, chr, hg18pos, hg18seq, varseq])

    j = j+1
#[genomic_positions, regions] = get_rid_of_redundancy(regions, genomic_positions)

simplejson.dump([genomic_positions, regions], outputfile)                
file.close()
outputfile.close()

#return barcode_snps



