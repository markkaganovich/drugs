# ----------------------------------------------------------------------------------
#
# Mark Kaganovich,  May 3, 2010
#
# find stretches of the genome that will identify 1000genome lines uniquely based on
# the recently released SNP calls
#
# -----------------------------------------------------------------------------------

from helperfuns_findbarcode import *
import simplejson

def findbarcode(population, dnastretch = 300):
#FILENAME = './'+population+'.genotypes.vcf'
#FILENAME = './test'
#FILENAME = './partial_input'
#FILENAME = './input' + population +'/SNPs_chrom1'
#   FILENAME = '../Desktop/release_20100608.CEU-merge.genotypes.phased_haplotypes.no_trio.vcf'
#    FILENAME = './temp'
    FILENAME = '../Desktop/CEU.low_coverage.indel.vcf'
    DATABYTES = 10000000
    DNASTRETCH = dnastretch
    POS = 1
    REF = 3
    ALT = 4
    INDEL_SIZE = 20

    file = open(FILENAME)
    outputfile_all_genos = open('./'+population+'allgenos_outputORIG', 'w')
    outputfile = open('./'+population+'output', 'w')

# read through header
    lines = file.readlines(DATABYTES)
    for z in range(0, len(lines)):
        print z
        if lines[z].startswith('#CHROM'):
            print lines[z]
            NUM_GENOMES = len(list(lines[z].split('\t')[9:-1]))
            print NUM_GENOMES         
            break
        else:
            pass
    
    barcode_indels=[]
    j=z+1 
    regions=[]
    genomic_positions=[]
# read blocks of data at a time: DATABYTES at a time
    while(1):
        if not lines:
            break
        print j
        print 'lines:  '+str(lines[j].split('\t')[0])
        genotype = vcf_to_geno_struct([lines[j]])
        if genotype:
            [num_genos, cell_line_ID] =num_genotypes(genotype)
            print num_genos
            hg18pos_start = int(lines[j].split('\t')[POS])
            hg18seq = str(lines[j].split('\t')[REF])
            varseq = str(lines[j].split('\t')[ALT])
            chr = int(lines[j].split('\t')[0])
            hg18pos_end = hg18pos_start + len(hg18seq)
#print hg18seq
#           print len(hg18seq)
            if num_genos == 1 and len(hg18seq) <= INDEL_SIZE and len(varseq) <= INDEL_SIZE:
                barcode_indels.append([j, cell_line_ID, chr, hg18pos_start, hg18pos_end, hg18seq, varseq])
#simplejson.dump([genomic_positions, regions], outputfile) 

        if len(lines) <= j+1:
            lines.extend(file.readlines(DATABYTES))
            if len(lines) <= j+1:                                # if reading more lines does not extend LINES, done
                break
        j = j+1
#[genomic_positions, regions] = get_rid_of_redundancy(regions, genomic_positions)
#   simplejson.dump([genomic_positions, regions], outputfile)                
    file.close()
    outputfile.close()
    outputfile_all_genos.close()
    
    return barcode_indels



