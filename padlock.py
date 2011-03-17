from pygr import *
from Bio.SeqUtils import MeltingTemp
from types import *
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
 check the results of the local blast run; pick probes that have only two 
 significant matches
'''
def blastcheck(blastoutputfile):
    EXPECT_THRESH = .02
    handle = open(blastoutputfile)
    records = list(NCBIXML.parse(handle))
    goodprobes = []

    for r in range(0,len(records)):
        sa = 0
        for al in records[r].alignments:
            if 'GRCh37.p2' in al.title:
                for hsp in al.hsps:
                    if hsp.expect < EXPECT_THRESH:
                        sa = sa+1
        if sa == 2:
            goodprobes.append(r)

    return goodprobes
            
    num_results  = len(filter(lambda x: 'GRCh37.p2' in x.title, record.alignments)) 
    print num_results
    if num_results == 1:
        return True
    else:
            return False

'''
check if barcode set includes everything in the pool

one of the cell lines covered by probes for >1 cell line
is assumed to be in the desired pool

if we are to pool cell lines covered by the same barcode
then there will need to an additional list of cell lines that have
been covered by barcodes as these checks are performed
'''
def checkset(barcodes):
    import info
    
    lines = info.lineinfo()
    YRIpool = info.pool().YRItest
    pool = [p for p in YRIpool if p not in lines.trios['YRI']]
    cells = lines.YRIfilecells
    padlockset = []
    for snp in barcodes:
        c = filter(lambda x: cells[x] in pool, snp[0])
        if c[0] not in padlockset:
            padlockset.append(cells[c[0]])
    
    if len(list(set(padlockset) & set(pool))) == len(pool):
        return True
    else:
        return False

def mapextend(ls):
    result=[]
    for e in ls:
        result.extend(e[0])
    return result

'''
probeset : the probes that pass the blast test
alreadyblasted : the probes that have already been blasted
blastset : the current set of sequences to be blasted this round
blastoutput : the output of the current round of blast

probecandidates (trueprobes) are the inputs that pass the Tm screens
'''
if __name__ == "__main__":
    import simplejson
    import info
    import os
    
    file = open('./trueprobes1')
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
    
    pool = info.pool()
    while not checkset(alreadyfound):
        print str(len(alreadyblasted)) + '  /  ' + str(len(probes))
        blastset = []
        for x in probes:
            if (list(set(x[0]) & set(mapextend(alreadyfound))) == [] and 
                    x not in alreadyblasted and 
                    list(set(x[0]) & set(mapextend(blastset))) == []):
                blastset.append(x)
        alreadyblasted.extend(blastset)
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

        blastn_cline = (NcbiblastnCommandline(query = blastsetfile, 
                    db="human_genomic", evalue=.1, outfmt=5, out=blastoutput, word_size=20))
        stdout, stderr = blastn_cline()
        alreadyfound.extend(map(lambda x: blastset[x], blastcheck(blastoutput)))


    file = open(alreadyblastedfile,'w')
    simplejson.dump(alreadyblasted, file)
    file.close()
    file = open(probesetfile,'w')
    simplejson.dump(alreadyfound, file)
    file.close()




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

    	
	

