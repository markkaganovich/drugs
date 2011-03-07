from pygr import *
from Bio.SeqUtils import MeltingTemp
from types import *

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

    def ligpos(self):
        return [self.pos-30, self.pos-5]
    
    def extpos(self):
        return [self.pos+5,self.pos+28]

    def checktemp(self):
        if self.pos-28 <= 1:
            return False
        lpos = self.ligpos()
        epos = self.extpos()
        tmlig = MeltingTemp.Tm_staluc(str(hg18['chr'+str(self.chr)][lpos[0]:lpos[1]]))
        tmext = MeltingTemp.Tm_staluc(str(hg18['chr'+str(self.chr)][epos[0]:epos[1]]))
        print tmlig, tmext

        if abs(tmlig - OPT_TEMP_L) <=1 and abs(tmext - OPT_TEMP_E) <= 1:
            return True
        else:
            return False
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
	

if __name__ == "__main__":
    import simplejson
    import padlock
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

    	
	

