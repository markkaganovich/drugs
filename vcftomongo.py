# convert VCF genotype files to MongoDB

from pymongo import Connection 

def createdb(group = 'CEU', filename = '../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf'):
    global db
    f = open(filename)
    keys = []
    lines = f.readlines(100000)
#for line in lines:
    while lines!=[]:
        for j in range(0, len(lines)):
            line = lines[j]
            line = line.strip('\n').split('\t')
            if line[0] == "#CHROM":
                keys = line
            if line[0].startswith('#'):
                continue
            record = {}
            individuals ={}
            for i, name in enumerate(keys):
                if name.startswith('NA'):
                    individuals[name] = line[i]
                else:
                    record[name.strip('#')] = line[i]
            record['individuals'] = individuals
            db[group].insert(record)
            print record['POS']
            if j+1 == len(lines):
                j = 0
                lines = f.readlines(100000)
    f.close()

if __name__ == "__main__":
	con = Connection()
	db = con['1000Genomes']
	
	print 'CEU'
	createdb()
	print 'CHBJPT'
	createdb('CHBJPT', '../1000GenomesData/CHBJPT.low_coverage.2010_09.genotypes.vcf')
	print 'YRI'
	createdb('YRI', '../1000GenomesData/YRI.low_coverage.2010_09.genotypes.vcf')
 
'''
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
'''



