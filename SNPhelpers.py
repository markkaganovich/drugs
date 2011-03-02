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
    cell_line_ID = []
    for lcl in range(0, len(genotypes[0])):
        if 1 in genotypes[0][lcl] or -1 in genotypes[0][lcl]:
            num_genos = num_genos+1
            cell_line_ID.append(lcl)

    return [num_genos, cell_line_ID]


#---------------------------------------------------------------------
# take vcf format LINES and output genotype struc:
#   NA123  [0] -> [0,1], [0,0], [1,0]
#   NA143  [1] ->
#
#---------------------------------------------------------------------
def vcf_to_geno_struct(lines, ignoreunphased = 0, countunphased=0):
    genotype = []
    for line in lines:
        tokens = line.split('\t')
        tokens = tokens[9:len(tokens)]
        genotype.append([])
        for token in tokens:
            if token[1] == '|' or countunphased == 1:
                if token[0] == '.' or token[2] == '.':
                    genotype[-1].append([-2,-2])
                else:
                    genotype[-1].append([int(token[0]), int(token[2])])
            else:
               if ignoreunphased:
                   pass
               else:
                   if token[1] == '/':
                       genotype[-1].append([-1,-1])
               if token[1] == '\\':
                    print "\\"
                    pass
    return genotype

