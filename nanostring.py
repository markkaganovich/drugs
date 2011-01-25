from helperfuns_findbarcode import *
from pygr import *
import simplejson as json
hg18 = worldbase.Bio.Seq.Genome.HUMAN.hg18()

def deletedregion(NUMGENOS = 1):
    file = open('./1000GenomesData/CEU.low_coverage.deletions.vcf')
    CEUlines = file.readlines()
    file.close()
    file= open('./1000GenomesData/YRI.low_coverage.deletions.vcf')
    YRIlines = file.readlines()
    file.close()
    file= open('./1000GenomesData/CHBJPT.low_coverage.deletions.vcf')
    CHBlines = file.readlines()
    file.close()
    
    LCLregions=[]
    file = open('./ALLcelllines')
    cell_lines = json.load(file)
    file.close()

    CHR = 0
    POS = 1
    REF = 3

    largeregions=[]
    largeandvalidated=[]

    for i in range(30,len(CEUlines)):
        l = []
        l.extend(CEUlines[i].split('\t'))
        l.extend(CHBlines[i].split('\t')[9:len(CHBlines[i].split('\t'))])
        l.extend(YRIlines[i].split('\t')[9:len(YRIlines[i].split('\t'))])
        if len(str(l[REF])) >=90 and 'VALIDATED' in l[7] and 'IMPRECISE' in l[7]:
            print l
        if len(str(l[REF])) >=300 and 'LowQual' not in CEUlines[i]:
            genotype = vcf_to_geno_justtokens(l[9:len(l)], 0 , 1)
            if genotype:
                celllines = onlynumgenos(genotype, NUMGENOS)
                if celllines != [-1]:
                    print i
                    LCLregions.append([celllines, int(l[CHR]), int(l[POS]), len(str(l[REF]))])
    
    LCLregions = sorted(LCLregions, key = lambda x: x[3], reverse=True)  
    return LCLregions
    
def checkgroups(proberegions, groups=[[],[],[]], regionsused =[]):

    for pair in proberegions:
        for x in pair[0]:
            for group in groups:
                if cangoingroup(x, group, proberegions):
                    if x not in group and pair not in regionsused and len(group) < 70:
                        group.append(x)
                        regionsused.append(pair)
                    break;
                
    return [regionsused, groups]

def cangoingroup(x, group, twoproberegions):
    
    if group == []:
        return True
    for y in group:
        if x == y:
            return True
        else:
            for pair in twoproberegions:
                if x in pair[0] and y in pair[0]:
                    return False

    return True

def checkgroups2(proberegions, groups):
    celllines = groups[0]+groups[1]+groups[2]
   
    print celllines
    for cell in celllines:
        print cell
        regions = filter(lambda x: cell in x[0], proberegions)
        regions = map(lambda x: x[0], regions)
        if [cell] in regions:
            return True
        t = filter(lambda x: len(x) == 2, regions)
        for i in t:
            g1 = [y for y in range(0,2) if i[0] in groups[y]]
            g2 = [y for y in range(0,2) if i[1] in groups[y]]
            if g1 != g2:
                return True
        t = filter(lambda x: len(x) == 3, regions)
        g = [y for y in range(0,2) if cell in groups[y]]
        for i in t:
            ind = [indeces for indeces in range(0,2) if i[indeces] != cell]
            for j in ind:
                g1 = [y for y in range(0,2) if j in groups[y]]
                if g != g1:
                    return True

    return False


# cover as many cell lines as possible with two probe and three probe groups, 
# then add the unique probe cell lines to the groups; this way some will have backups
# -- this can also be done in the reverse order
def adduniqueregions(uniqueregions, groups, regionsused):

    for pair in uniqueregions:
        for group in groups:
            if len(group) > 70 and pair[0][0]!=42:    # need to miake an exception for the trio.... CEU trio will be in one pool
                pass;
            else:
                if filter(lambda x: pair[0][0] in x, groups) == []:
                    group.append(pair[0][0])
                    regionsused.append(pair)

    return [regionsused, groups]


a = deletedregion()
b = deletedregion(2)
c = deletedregion(3)

d = []
d.extend(b)
d.extend(c)

#[regionsused, groups] = checkgroups(d)
#[regionsused, groups] = adduniqueregions(a, groups, regionsused)

### uniques first
[regionsused, groups] = adduniqueregions(a, [[],[],[]], [])
[regionsused, groups] = checkgroups(d, groups, regionsused)

## manually modify some groups, and make them into ~20 each
groups[2].extend(groups[0][63:])
groups[0] = groups[0][0:63]

newgroups = []
for j in range(0, 3):
    for i in range(0,3):
        newgroups.append(groups[j][i*21:(i+1)*21])
        if 42 not in newgroups[-1]:
            newgroups[-1].append(42)

newgroups[5].extend(newgroups[7])
newgroups[5] = newgroups[5][0:-1]
del newgroups[7]
del newgroups[7]

#################################################################
#
# retrieve sequences
#     input: regions dictionary, cell lines 
#
#################################################################
def retrieveseqs(regions, additionalregions):
    file = open('./ALLcelllines')
    cell_lines = json.load(file)
    file.close()
    file = open('nanostring_deletions_160probes','w')
    file.write('cell line\tdel length\tchr\tbeggining-50\tend+50\n')
    for region in regions:
        lcls = region[0]
        rg = region[1:]
        print lcls
#t = sorted(regions[lcl], key=lambda x: x[2], reverse=True)
        for lcl in lcls: file.write(cell_lines[lcl]+'\n')
        file.write( '\t' + str(rg[2]) + '\t' + str(rg[0]) + '\t' + str(hg18['chr'+str(rg[0])][rg[1]-50:rg[1]]) +'\t' + str(hg18['chr'+str(rg[0])][rg[1]+rg[2]:rg[1]+rg[2]+50]) +'\n')

    file.write('\n\n\n')
    file.write('backup probes\n')
    for region in additionalregions:
        lcls = region[0]
        print lcls
        rg = region[1:]
        for lcl in lcls: file.write(cell_lines[lcl]+'\n')
        file.write( '\t' + str(rg[2]) + '\t' + str(rg[0]) + '\t' + str(hg18['chr'+str(rg[0])][rg[1]-50:rg[1]]) +'\t' + str(hg18['chr'+str(rg[0])][rg[1]+rg[2]:rg[1]+rg[2]+50]) +'\n')
        
    file.close()






