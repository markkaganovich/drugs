# convert VCF genotype files to MongoDB

from pymongo import Connection 

def createdb(filename = '../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf'):
    global db
    f = open(filename)
    lines = f.readlines()
    keys = []
    for line in lines:
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
        db[filename.split('.')[1]].insert(record)
        print filename.split('.')[1]
    f.close()

if __name__ == "__main__":
	con = Connection()
	db = con['1000Genomes']
	
	print 'CEU'
	createdb()
	print 'CHBJPT'
	createdb('../1000GenomesData/CHBJPT.low_coverage.2010_09.genotypes.vcf')
	print 'YRI'
	createdb('../1000GenomesData/YRI.low_coverage.2010_09.genotypes.vcf')
 
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



