from pygr import *
from Bio.SeqUtils import MeltingTemp
from types import *

hg18 = worldbase.Bio.Seq.Genome.HUMAN.hg18()
OPT_TEMP_L = 58
OPT_TEMP_E = 53
BACKBONE = 32
LENGTH = 24

class Primer:
    def __init__(self, chr, pos):
    	self.chr = chr
	self.pos = pos    


    def checktemp(self):
	if self.pos-28 <= 1:
	    return False
	tmlig = MeltingTemp.Tm_staluc(str(hg18['chr'+str(self.chr)][self.pos-29:self.pos-4]))
	tmext = MeltingTemp.Tm_staluc(str(hg18['chr'+str(self.chr)][self.pos+4:self.pos+27]))
	print tmlig, tmext

	if abs(tmlig - OPT_TEMP_L) <=1 and abs(tmext - OPT_TEMP_E) <= 1:
	    return True
	else:
	    return False

if __name__ == "__main__":
    import simplejson
    import padlock
    file = open('./YRIoutputUniqueSNPs')
    b = simplejson.load(file)
    file.close()

    trueprobes1 = []
    trueprobes2 = []
    i = 0
    for snp in b[0]:
	i=i+1
	print i
	if Primer(snp[1],snp[2]).checktemp():
            trueprobes1.append(snp)
	
    file = open('./trueprobes1')
    simpejson.dump(trueprobes1, file)
    file.close()

    	
	

