# ----------------------------------------------------------------------------------
#
# Mark Kaganovich,  May 3, 2010
#
# find stretches of the genome that will identify 1000genome lines uniquely based on
# the recently released SNP calls
#
# -----------------------------------------------------------------------------------

from helperfuns_findbarcode import *
population = 'CEU'
dnastretch = 100
#FILENAME = './'+population+'.genotypes.vcf'
#FILENAME = './test'
#FILENAME = './partial_input'
FILENAME = './input' + population +'/SNPs_chrom1'
DATABYTES = 10000000
DNASTRETCH = dnastretch
NUM_GENOMES = 59

file = open(FILENAME)
outputfile_all_genos = open('./'+population+'allgenos_output', 'w')
outputfile = open('./'+population+'output', 'w')

lines = file.readlines(DATABYTES)
startchrom = int(lines[0].split('\t')[0])
startpos = int(lines[0].split('\t')[1])
j=1

# read blocks of data at a time: DATABYTES at a time
while(1):
#print j
    if not lines:
        break
    if int(lines[j].split('\t')[0]) == startchrom and int(lines[j].split('\t')[1]) - startpos < DNASTRETCH:
               j=j+1
    else:
        genotype = vcf_to_geno_struct(lines[0:j])
        if genotype:
            [diff_genos, lcls_unphased] = gethaps(genotype)
            if len(diff_genos) > 20:
                print diff_genos 
                print "\n"
                print len(diff_genos[0])
                print len(diff_genos[1])
                print j
                temp_lines = lines[0:j]
                temp_genotype = genotype
                temp_lcl = lcls_unphased
                temp_diff_Genos = diff_genos
                lcl_groups = group_lcls(diff_genos, lines[0:j], lcls_unphased)
                outputfile_all_genos.write(str(startchrom) + "\t" + str(startpos) + "\t" + lines[j-1].split('\t')[1] + "\t" + str(len(diff_genos)))
                for cells in range(0,len(lcl_groups)):
                    outputfile_all_genos.write("\t")
                    for item in range(0,len(lcl_groups[cells])):
                        if not item == 0:
                            outputfile_all_genos.write(",")
                        outputfile_all_genos.write(str(lcl_groups[cells][item]))
                outputfile_all_genos.write("\n")
            if len(diff_genos) == NUM_GENOMES:
                outputfile.write(str(startchrom) + "\t" + str(startpos) + "\t" + lines[j-1].split('\t')[1] + "\t" + str(len(diff_genos)) + "\n")

        del lines[0]                                        # remove previous start line from lines and
        j=j-1
        startpos = int(lines[0].split('\t')[1])             # update startpos on same chromosome to next line,
        
        if not int(lines[j].split('\t')[0]) == startchrom:
            startchrom = int(lines[j].split('\t')[0])       # not on same chrom, so start new startpos and new genotype
            print startchrom
            startpos = int(lines[j].split('\t')[1])
            outputfile.close()
            outputfile = open(population+'output'+str(startchrom), 'w')

    if len(lines) < j+1:
        lines.extend(file.readlines(DATABYTES))
        if len(lines) < j+1:                                # if reading more lines does not extend LINES, done
            break
            
file.close()
outputfile.close()
outputfile_all_genos.close()






