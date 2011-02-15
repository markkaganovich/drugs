# convert VCF genotype files to MongoDB

from pymongo import Connection 

connection = Connection()
db = connection['1000Genomes']

def createdb( file = '../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf'):
    f = open(file)
    lines = f.readlines()
    for line in lines:
        line = line.strip('\n')
        line = line.split('\t')
        if str(line[0]) == '#CHROM':
            cell_lines = line[9:]
            parameters = line[0:8]
        if '#' not in str(line[0]):
            record={}
            individuals={}
            for par in parameters:
                record[str(par)] = line[parameters.index(par)]
            for individual in cell_lines:
                individuals[str(individual)] = line[cell_lines.index(individual)+9]
            record.update(individuals)
            print eval('db.SNPs'+line[0]).count()
            print line[0]
            poslist= eval('db.SNPs'+line[0]).distinct('POS')
            if line[1] in poslist:
                eval('db.SNPs'+line[0]).update({'POS':line[1]}, {'$set':individuals})
            else:
                eval('db.SNPs'+line[0]).insert(record)

    f.close() 




