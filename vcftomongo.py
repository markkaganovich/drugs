# convert VCF genotype files to MongoDB

from pymongo import Connection 

connection = Connection()
db = connection['1000Genomes']
SNPs = db.SNPs
file = open('../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf')
lines = []

while(not lines):
    lines = file.readlines(1000000000)
    if not lines: break;
    print lines[0]
    for i in range(0, len(lines)):
        line = lines[i].split('\t')
        if str(line[0]) == '#CHROM':
            cell_lines = line[9:]
            parameters = line[0:8]
        if '#' not in str(line[0]):
            record = {}
            for par in parameters:
                record[str(par)] = line[parameters.index(par)]
            individuals_record = {}
            for individual in cell_lines:
                individuals_record[str(individual)] = line[cell_lines.index(individual)+9]
        
            record['individuals'] = individuals_record    
            SNPs.insert(record)
    lines = []
    print i 
            
file.close() 

# add info from other files

file = open('../1000GenomesData/CHBJPT.low_coverage.2010_09.genotypes.vcf')
lines = file.readlines()

CHROM = 0
POS = 1

newSNPs = db.newSNPs

for line in lines:
	line = line.split('\t')
	if str(line[0]) == '#CHROM':
        cell_lines = line[9:]
    if '#' not in str(line[0]):
        tmp = SNPs.find({'#CHROM': line[CHROM]})
        if tmp:
            record = {}
            for par in parameters:
                record[str(par)] = line[parameters.index(par)]
            individuals_record = {}
            for individual in cell_lines:
                individuals_record[str(individual)] = line[cell_lines.index(individual)+9]
            record['individuals'] = individuals_record    
            newSNPs.insert(record)
        else:
            indtemp = tmp['individuals']
            for individual in cell_lines:
                indtemp[str(individual)] = line[cell_lines.index(individual)+9]
            tmp['individuals'] = indtemp
            



