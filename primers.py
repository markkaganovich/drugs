#################################################################
#
# 1000Genomes uses hg18 (like hapmap); pygr genomic info is 0 indexed
#
#
#
################################################################
from pygr import *
from Bio.SeqUtils import MeltingTemp
from types import *

hg18 = worldbase.Bio.Seq.Genome.HUMAN.hg18()
OPT_TEMP_L = 58
OPT_TEMP_E = 53
REGION_LENGTH = 100

'''
input position (of start of barcode), return positions of primers
in the region, with Tms... extending up or down depending on whether
its a ligation (l) or extension (e) arm
'''
def helper1(chrom, pos, type = 'ligation'):
    
    if type == 'extension':
        opt_temp = OPT_TEMP_E
    else:
        opt_temp = OPT_TEMP_L

    pos = pos - 1 # indexing...
    i = 0
    temps=[]
    while(abs(i) <= 24):
        if type == 'extension':
            i = i+1
        else:
            i = i-1
        if abs(i) >= 17:
            Tm = MeltingTemp.Tm_staluc(str(hg18['chr'+str(chrom)][pos:pos+i]))
            temps.append([i, Tm])
        
    temps_sorted = sorted( temps , key=lambda temp: abs(temp[1] - opt_temp))
   
    return temps_sorted
'''
take genomic positions, call the helper1 function, find the oligo >16 <25 that has the closest 
Tm to the OPTIMAL_TMs
'''
def calc_primers(chr, genomic_position, var_length = 100):
    primers_L = []
    primers_E = []
    L = helper2(chr, genomic_position,'ligation')
    E = helper2(chr, genomic_position+var_length, 'extension')
    if type(E) != NoneType and type(L) != NoneType:
        primers_L.extend(L)
        primers_E.extend(E)
    
    return [primers_L, primers_E]

def SNPbarcode_primers(genomic_positions, var_length = 100):
    for region in genomic_positions:
        [primers_L, primers_E] = calc_primers(region[0], region[1], var_length)
    return [primers_L, primers_E]

def helper2(chr, pos, type):
    p = 0
    while( p<4):
        if type == 'ligation':
            newpos = pos-2*p
            temp = helper1(chr, newpos, type)[0]
            opt_temp = OPT_TEMP_L
        else:
            newpos = pos+2*p
            temp = helper1(chr, newpos, type)[0]
            opt_temp = OPT_TEMP_E
        i= temp[0]
        Tm = temp[1]
        if abs(Tm) > opt_temp - 3 and abs(Tm) < opt_temp + 3:
            return [chr, newpos,  newpos+i, Tm]
        else:
            p = p+1
        if p ==4:
            return None


sequence = hg18['chr1'][13104990:13105015]
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
blast_handle = NCBIWWW.qblast('blastn','nt', sequence)
blastrecord = NCBIXML.read(blast_handle)

for alignment in blastrecord.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < .1:
            print alignment.title
            print hsp.expect
            print hsp.query[0:75]
            print hsp.match[0:75]
            print hsp.sbjct[:75]

