from pygr import *
hg18 = worldbase.Bio.Seq.Genome.HUMAN.hg18()
hg19 = worldbase.Bio.Seq.Genome.HUMAN.hg19()
from Bio import pairwise2
import random

def checkilluminasnps(probesfile, probetype = 1, readsfile =
'110617_MONK_0184_B81N3EABXX_L4_eland_extended_pf.txt'):
    import simplejson

    file = open(probesfile)
    probes = file.readlines()
    file.close()
    alts = [0] * len(probes)
    refs = [0] * len(probes)
    others = [0] * len(probes)
    varinds = [[]] * len(probes)

    reads = open('./'+readsfile)
    lines = reads.readlines(100000)

    while(lines != []):
        for l in lines:
            r = random.random()
            if r < .5:
                continue;
            read = l.strip('\n').split('\t')
            readseq = read[1]
            readmappedinfo = read[3]
            if len(readmappedinfo) > 1:
                readchr = readmappedinfo.split('.')[0]
                if 'F' in readmappedinfo.split(':')[1] and ',' not in readmappedinfo:
                    readpos = int(readmappedinfo.split(':')[1].split('F')[0]) - 1
                else:
                    continue;
                for i in range(0,len(probes)):
                    prob = probes[i].strip('\n').split('\t')
                    probchr = prob[0].split(':')[0]
                    probpos = int(prob[0].split(':')[1])-1
                    REF = prob[1]
                    ALT = prob[2]
                    if probchr == readchr and probpos > readpos and probpos < readpos + 101:
                        print readseq[probpos - readpos]
                        print prob[1]
                        print prob[2]
                        print i
                        if readseq[probpos - readpos] == REF:
                            refs[i] = refs[i] + 1
                        elif readseq[probpos - readpos] == ALT:
                            alts[i] = alts[i] + 1
                        else:
                            others[i] = others[i] + 1
        lines = reads.readlines(100000)

    return [refs, alts, others]

def getsnpregionseq(probe):
        return str(hg18['chr'+str(probe[0])][int(probe[1])-30 : int(probe[1])+28])
        
def findnumalleles(probe, genotypefile):
    probechr = int(probe[1])
    probepos = int(probe[2])

    print probechr
    print probepos

    file = open(genotypefile)
    lines = file.readlines()
    file.close()

    for l in lines:
        line = l.split('\t')
        if l.startswith('#') != True and int(line[0]) == probechr and int(line[1]) == probepos:
            return int(line[7].split(';')[1].split('=')[1])

'''
file = open('./primersbeforeblast1007')
primerlist = simplejson.load(file)
for p in primerlist:
    primerseq = getsnpregionseq(p)
    for read in readlist:
        alignments = pairwise2.align.globalxx(primerseq, read)
'''     







