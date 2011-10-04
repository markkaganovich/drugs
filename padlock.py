from pygr import *
from Bio.SeqUtils import MeltingTemp
from types import *
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import Seq


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
    def __init__(self, chr, pos, snpline):
    	import info
        lines = info.lineinfo()
        cells = lines.individuals['CEU'] + lines.individuals['CHBJPT'] + lines.individuals['YRI']
        self.chr = chr
        self.pos = pos    
        self.lpos = self.ligpos()
        self.epos = self.extpos()
        self.seqlig = str(hg18['chr'+str(self.chr)][self.lpos[0]:self.lpos[1]]) 
        self.seqext = str(hg18['chr'+str(self.chr)][self.epos[0]:self.epos[1]])
        self.backbone = 'AGATCGGAAGAGCGTCGCATGTTATCGAGGTC'
        self.line = cells[snpline]  

    def ligpos(self):
        return [self.pos-30, self.pos-5]
    
    def extpos(self):
        return [self.pos+5,self.pos+28]

    def checktemp(self):
        if self.pos-28 <= 1:
            return False
        tmlig = MeltingTemp.Tm_staluc(self.seqlig)
        tmext = MeltingTemp.Tm_staluc(self.seqext)

        print tmlig, tmext

        if abs(tmlig - OPT_TEMP_L) <=1 and abs(tmext - OPT_TEMP_E) <= 1:
            return True
        else:
            return False

    def checkGCcontent(self):
        lenlig = float(len(self.seqlig))
        lenext = float(len(self.seqext))
        GClig = (len(filter(lambda x: x == 'G' or x == 'C' or x == 'c' or x == 'g', 
                    self.seqlig)))
        GCext = (len(filter(lambda x: x == 'G' or x == 'C' or x == 'c' or x == 'g', 
                        self.seqext)))

        print lenlig, lenext, GClig, GCext

        if (GClig/lenlig > .3 and GClig/lenlig <.7 and GCext/lenext > .3 and
                GCext/lenext < .7):
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

def printprobesandorderthem(probesfile, probesoutputfile):
    import os
    import simplejson
    import Bio.Seq
    import info

    if probesoutputfile in os.listdir('./'):
        print "WARNING: file exists"

    file = open(probesfile)
    probes = simplejson.load(file)
    file.close()
    lines = info.lineinfo()
    YRIpool = info.pool().YRItest
    pool = [p for p in YRIpool if p not in lines.trios['YRI']]

    primers = map(lambda x: Primer(x[1], x[2], x[0]), probes)
    file = open(probesoutputfile,'w')
    file.write('NAME'+','+ 'SEQUENCE' + '\n')
    for p in primers:
        print p.seqlig
        c = filter(lambda x: x in pool, p.lines)
        poolcell = str(c[0])
        (file.write(poolcell[2:] + str('Pl') + str(len(p.lines)) + ',' + 
                    Bio.Seq.reverse_complement(p.seqlig) + p.backbone
                    + Bio.Seq.reverse_complement(p.seqext) + '\n'))

    file.close()
    

def getsnpregionseq(probe):
    return str(hg18['chr'+str(probe[1])][int(probe[2])-30 : int(probe[2])+28])

''' 
run all the tests
'''
def alleles(snps, numalleles):
    newsnps = filter(lambda x: x[5] == numalleles, snps)
    return newsnps

def runTests(uniquesnps):
    trueprobes = []
    i = 0
    for snp in uniquesnps:
        i = i+1
        print i
        p = Primer(snp[1], snp[2], snp[0][0])
        if p.checktemp() and p.checkGCcontent():
            trueprobes.append(p)
    return trueprobes

def addtosetCheck(primer, padlockset):
    if primer.line not in padlockset.keys():
        return True
    elif padlockset[primer.line].__len__() >= 5:
        return False
    else:
        return True

def blast(primer):
    blastsetfile = './blastresults/blastset.fasta'
    blastoutput = './blastresults/blastoutput.xml'
    file = open(blastsetfile,'w')
    file.write('>' + 'Chr' + str(primer.chr) + 'Pos' + str(primer.pos) + str(primer.line) +'\n')
    seq = (str(hg18['chr'+str(primer.chr)][primer.lpos[0]:primer.lpos[1]]) 
    + str(hg18['chr'+str(primer.chr)][primer.epos[0]:primer.epos[1]])) 
    print primer.chr, primer.pos, seq
    file.write(seq)
    file.close()
    blastn_cline = (NcbiblastnCommandline(query = blastsetfile, 
                db="human_genomic", evalue=.1, outfmt=5, out=blastoutput, word_size=20))
    stdout, stderr = blastn_cline()
    if blastcheck(blastoutput):
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
            
#num_results  = len(filter(lambda x: 'GRCh37.p2' in x.title, record.alignments)) 
#   print num_results
#   if num_results == 1:
#       return True
#   else:
#       return False


if __name__ == "__main__":
    import simplejson
    import info
    import os
    
#file = open('../1000GenomesData/low_coverage.merged.vcf.uniquesnps')
    file = open('./probes1')
    candidateloci = simplejson.load(file)
    file.close()

    # take only num alleles == 2
#a = alleles(candidateloci, 2)
    a = candidateloci

    # check Tm
    primersbeforeblast = runTests(a)
   
    lines = info.lineinfo()
    cells = lines.individuals['CEU'] + lines.individuals['CHBJPT'] + lines.individuals['YRI']
    padlockset = {}
    for p in primersbeforeblast:
        if addtosetCheck(p, padlockset) and blast(p):
            if p.line in padlockset.keys():
                padlockset[p.line] = padlockset[p.line].append(p)
            else:
                padlockset[p.line] = [p]


    file = open('./padlockset','w')
    simplejson.dump(padlockset, file)
    file.close()
