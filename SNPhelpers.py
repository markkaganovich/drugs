def combinefiles(file1, file2):
    f1 = open(file1)
    f2 = open(file2)
    lines1 = f1.readlines()
    lines2 = f2.readlines()
   
    i = 0
    for line1 in lines1:
        print i
        l1 = line1.split('\t')
        i = i+1
        for line2 in lines2:
            l2 = line2.split('\t')
            if int(l2[1]) > int(l1[1]) or (l2[1] == l1[1] and l2[2]>l1[2]):
                break
            

#################################################################
#
# for indel barcode: return number of genomes that have the variant genotype
# Oct 11, 2010
########################################################################
def num_genotypes(genotypes):
    num_genos = 0
    lineID = []
    num_alleles = 0
    genotypes = genotypes[0]
    print "length of genotypes:" + str(len(genotypes))
    for lcl in range(0, len(genotypes)):
        if 1 in genotypes[lcl] or -1 in genotypes[lcl]:
            num_genos = num_genos+1
            lineID.append(lcl)
    num_alleles = sum(map(lambda x: sum(genotypes[x]), lineID))
        
    return [num_genos, num_alleles, lineID]


def num_genotypes_trios(genotypes):
    import info

    num_genos   = []
    num_alleles = []
    lineID      = []
    for geno in genotypes:
        print geno
        a = num_genotypes([geno])
        num_genos.append(a[0])
        num_alleles.append(a[1])
        lineID.append(a[2])

    if num_genos[0] == 2 and num_genos[1] == 2:
        i = info.lineinfo()
        snp1cells = map(lambda x: i.all[x], lineID[0])
        snp2cells = map(lambda x: i.all[x], lineID[1])
        if i.child['YRI'] in snp1cells and i.child['YRI'] in snp2cells:
            if i.parents['YRI'][0] in snp1cells+snp2cells and i.parents['YRI'][1] in snp1cells+snp2cells:
                return i.child['YRI']
    
#---------------------------------------------------------------------
# take vcf format LINES and output genotype struc:
#   NA123  [0] -> [0,1], [0,0], [1,0]
#   NA143  [1] ->
#
#---------------------------------------------------------------------
def vcf_to_geno_struct(lines, ignoreunphased = 0, countunphased=0):
    genotype = []
    print lines
    for line in lines:
        tokens = line.split('\t')
        tokens = tokens[9:len(tokens)]
        print "len of tokens:" + str(len(tokens))
        genotype.append([])
        for token in tokens:
            if token == '.':
                genotype[-1].append([0,0])
            elif token[1] == '|' or countunphased == 1:
                if token[0] == '.' or token[2] == '.':
                    genotype[-1].append([-2,-2])
                else:
                    genotype[-1].append([int(token[0]), int(token[2])])
            else:
               if ignoreunphased:
                   pass
               else:
                   if token[1] == '/':
                       genotype[-1].append([-1*int(token[0]), -1*int(token[2])])
               if token[1] == '\\':
                    print "\\"
                    pass
    return genotype

