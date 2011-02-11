# convert VCF genotype files to MongoDB

from pymongo import Connection 

connection = Connection()
db = connection['1000Genomes']
SNPs = db.SNPs2

def createdb():
    file = open('../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf')
    lines = file.readlines()
    for line in lines:
        line = lines.strip('\n').split('\t')
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

SNPs = db.SNPs

for line in lines:
    line = line.split('\t')
    if str(line[0]) == '#CHROM':
        cell_lines = line[9:]
        parameters = line[0:8]
    if '#' not in str(line[0]):
        orig = SNPs.find({'#CHROM': line[CHROM]})
        if orig.count() == 0: 
            record = {}
            for par in parameters:
                record[str(par)] = line[parameters.index(par)]
            individuals_record = {}
            for individual in cell_lines:
                individuals_record[str(individual)] = line[cell_lines.index(individual)+9]
            record['individuals'] = individuals_record    
            SNPs.insert(record)
        else:
            indtemp = orig['individuals']
            newrecord = orig.copy()
            for individual in cell_lines:
                if str(individual) in indtemp.keys():
                    print "error: why the shit are there the same cell lines here"
                indtemp[str(individual)] = line[cell_lines.index(individual)+9]
            newrecord['individuals'] = indtemp
            SNPs.update(orig, newrecord, upsert = True)





