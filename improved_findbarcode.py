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
    FILENAME = './release_20100608.CEU-merge.genotypes.phased_haplotypes.no_trio.vcf'
    DATABYTES = 10000000
    DNASTRETCH = dnastretch
    NUM_GENOMES = 59
    
    file = open(FILENAME)
    outputfile_all_genos = open('./'+population+'allgenos_outputORIG', 'w')
    outputfile = open('./'+population+'output', 'w')

    lines = file.readlines(DATABYTES)
    for z in range(0, len(lines)):
        print z
        if lines[z].startswith('#CHROM'):
            print lines[z]
            print len(list(lines[z].split('\t')[9:-1]))
            lines = lines[z+1:-1]
            break
        else:
            pass
            
    startchrom = int(lines[0].split('\t')[0])
    startpos = int(lines[0].split('\t')[1])
    j=1
    regions=[]
    genomic_positions=[]
# read blocks of data at a time: DATABYTES at a time
    while(1):
        if not lines:
            break
        print j
        print 'lines:  '+str(lines[j].split('\t')[0])
        if int(lines[j].split('\t')[0]) == startchrom and int(lines[j].split('\t')[1]) - startpos < DNASTRETCH:
                   j=j+1
        else:
            genotype = vcf_to_geno_struct(lines[0:j])
            if genotype:
                [diff_genos, lcls_unphased] = gethaps(genotype)
                if len(diff_genos) > 58:
                    lcl_groups = group_lcls(diff_genos, lines[0:j], lcls_unphased)
                    regions.append(lcl_groups)
                    genomic_positions.append([startchrom, startpos, int(lines[j-1].split('\t')[1])])
            del lines[0]                                        # remove previous start line from lines and
            j=j-1
            startpos = int(lines[0].split('\t')[1])             # update startpos on same chromosome to next line,
            
            if not int(lines[j].split('\t')[0]) == startchrom:
                startchrom = int(lines[j].split('\t')[0])       # not on same chrom, so start new startpos and new genotype
                print "CHROM"
                print startchrom
                startpos = int(lines[j].split('\t')[1])
                outputfile.close()
                outputfile = open(population+'output'+str(startchrom), 'w')
                simplejson.dump([genomic_positions, regions], outputfile) 

        if len(lines) < j+1:
            lines.extend(file.readlines(DATABYTES))
            if len(lines) < j+1:                                # if reading more lines does not extend LINES, done
                break
    [genomic_positions, regions] = get_rid_of_redundancy(regions, genomic_positions)
    simplejson.dump([genomic_positions, regions], outputfile)                
    file.close()
    outputfile.close()
    outputfile_all_genos.close()
    
    return [genomic_positions, regions]



