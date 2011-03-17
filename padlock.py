from pygr import *
from Bio.SeqUtils import MeltingTemp
from types import *
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline


hg18 = worldbase.Bio.Seq.Genome.HUMAN.hg18()
OPT_TEMP_L = 58
OPT_TEMP_E = 53
BACKBONE = 32
LENGTH = 24

'''
11 bp gap between ligation and extension strands
look for primers that are within 1 degree of optimal Tms
'''

class Primer:
    def __init__(self, chr, pos):
    	self.chr = chr
        self.pos = pos    
        self.lpos = self.ligpos()
        self.epos = self.extpos()
        self.seqlig = hg18['chr'+str(self.chr)][self.lpos[0]:self.lpos[1]] 
        self.seqext = hg18['chr'+str(self.chr)][self.epos[0]:self.epos[1]]

    def ligpos(self):
        return [self.pos-30, self.pos-5]
    
    def extpos(self):
        return [self.pos+5,self.pos+28]

    def checktemp(self):
        if self.pos-28 <= 1:
            return False
        tmlig = MeltingTemp.Tm_staluc(str(self.seqlig))
        tmext = MeltingTemp.Tm_staluc(str(self.seqext))

        print tmlig, tmext

        if abs(tmlig - OPT_TEMP_L) <=1 and abs(tmext - OPT_TEMP_E) <= 1:
            return True
        else:
            return False

    '''
    run local blast
    '''
    def blastcheck(self):
        seq = (str(hg18['chr'+str(self.chr)][self.lpos[0]:self.lpos[1]]) +
         'NNNNNNNNNNN' + str(hg18['chr'+str(self.chr)][self.epos[0]:self.epos[1]]))
        print seq
        blastn_cline = NcbiblastxCommandline(query = seq, db="human_genomic", evalue=.01, outfmt=5, out="barcodeblast.xml")
        stdout, stderr = blastn_cline()
        handle = open("barcodeblast.xml")
        record = NCBIXML.read(handle)
        num_results  = len(filter(lambda x: 'GRCh37.p2' in x.title, record.alignments)) 
        print num_results
        if num_results == 1:
            return True
        else:
            return False
    '''       
    def blastprint(record):
        for alignment in record.alignments:
            for hsp in alignments.hsps:
    '''
'''
check if barcode set includes everything in the pool
'''
def checkset(barcodes):
    import info
    
    pool = info.pool().YRItest
    cells = info.lineinfo().YRIfilecells
    set = []
    for snp in barcodes:
        if cells[snp[0][0]] not in set:
            set.append(cells[snp[0][0]])

    #covered = filter(lambda x:  x in set, pool)
  
    return set 
	

def mapextend(ls):
    result=[]
    for e in ls:
        result.extend(e[0])
    return result

if __name__ == "__main__":
    import simplejson
    import info
    import os
    
    file = open('./trueprobes2')
    probes = simplejson.load(file)
    file.close()
    
    probesetfile = './probes'
    alreadyblastedfile = './blastresults/blasted'
    blastsetfile = './blastresults/blastset.fasta'
    blastoutput = './blastresults/blastoutput.xml'

    if probesetfile in os.listdir('./'):
        file = open(probsetfile)
        alreadyfound = simplejson.load(file)
        file.close()
    else:
        alreadyfound = []

    if alreadyblastedfile in os.listdir('./blastresults'):
        file = open(alreadyblastedfile)
        alreadyblasted = simplejson.load(file)
        file.close()
    else:
        alreadyblasted=[] 
    
    blastset = []

    for x in probes:
        if (list(set(x[0]) & set(mapextend(alreadyfound))) == [] and 
                x not in alreadyblasted and 
                list(set(x[0]) & set(mapextend(blastset))) == []):
            blastset.append(x)

    alreadyblasted.append(blastset)
    file = open(alreadyblastedfile,'w')
    simplejson.dump(alreadyblasted, file)
    file.close()

    pool = info.pool()

    file = open(blastsetfile,'w')
    for snp in blastset:
        file.write('>' + str(snp)+'\n')
        p = Primer(snp[1], snp[2])
        seq = (str(hg18['chr'+str(p.chr)][p.lpos[0]:p.lpos[1]]) 
        + str(hg18['chr'+str(p.chr)][p.epos[0]:p.epos[1]])) 
        print snp[1], snp[2], seq
        file.write(seq)
        file.write('\n')
    
    file.close()

    blastn_cline = NcbiblastnCommandline(query = blastsetfile, db="human_genomic", evalue=.1, outfmt=5, out=blastoutput, word_size=20)
    stdout, stderr = blastn_cline()






''' 
run all the tests
'''
def runalltests():
    import simplejson
    file = open('./YRIoutputUniqueSNPs')
    b = simplejson.load(file)
    file.close()

    for j in range(0,len(b)):
        trueprobes = []
        i = 0
        for snp in b[j]:
            i = i+1
            print i
            if Primer(snp[1], snp[2]).checktemp():
                trueprobes.append(snp)

        file = open('./checkedpadlocks'+str(j),'w')
        simplejson.dump(trueprobes, file)
        file.close()

    	
	

