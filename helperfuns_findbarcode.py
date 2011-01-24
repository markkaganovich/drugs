def comphaps(genotype_list):
    diff_genotypes = []
    for geno in genotype_list:
        if not geno in diff_genotypes:
            diff_genotypes.append(geno)

    return diff_genotypes
#---------------------------------------------------------------------
# take chrom and pos and return the vcf format lines (the SNP 
#  genotype)
#
#---------------------------------------------------------------------
def lookup_geno(chrom, start, stop, surr = 0, outputfile="lookup_output"):
    inputdirs = ('./inputCEU/', './inputCHB+JPT/', './inputYRI/')
    
    SURROUNDING_REGION = 5000

    results =[]

    for dir in inputdirs:
        file = open(dir+'SNPs_chrom'+str(chrom))
#lines.append([])
#output = open(outputfile + dir[7:-1], 'w')
        result=[]
        while(1):
            line = file.readline()
            if not line:
                break
            tokens = line.split('\t')
            if int(tokens[1]) >= start - surr*SURROUNDING_REGION:
#               output.write(line)
                result.append(line)
            if int(tokens[1]) >= stop + surr * SURROUNDING_REGION:
                break
        results.append(result)
    return results
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

def vcf_to_geno_justtokens(tokens, ignoreunphased = 0, countunphased = 0):

    genotype = []
    for token in tokens:
        if token[1] == '|' or countunphased == 1:
            if token[0] == '.' or token[2] == '.':
                genotype.append([-2,-2])
            else:
                genotype.append([int(token[0]), int(token[2])])
        else:
           if ignoreunphased:
               pass
           else:
               if token[1] == '/':
                   genotype.append([-1,-1])
           if token[1] == '\\':
                print "\\"
                pass
    
    return genotype


######################################################################
#
# input: genotypes
# output: (1,0,1) (0,1,0)
#               
######################################################################

def gethaps(genotypes):
    # look for -1s
    lcls_unphased = []
    

    for lcl in range(0,len(genotypes[0])):
        for snp in range(0,len(genotypes)):
            if -1 in genotypes[snp][lcl]:
                lcls_unphased.append(lcl)

    haps = []

    for lcl in range(0,len(genotypes[0])):
        for allele in (0,1):
            temp = []
            for snp in range(0,len(genotypes)):
                temp.append(genotypes[snp][lcl][allele])
            if not temp in haps and not lcl in lcls_unphased:
                haps.append(temp)
    return [haps, lcls_unphased]
                
#################################################################
#
# for indel barcode: return number of genomes that have the variant genotype
# Oct 11, 2010
########################################################################
def num_genotypes(genotypes):
    num_genos = 0
    cell_line_ID = 0
    for lcl in range(0, len(genotypes[0])):
        if 1 in genotypes[0][lcl] or -1 in genotypes[0][lcl]:
            num_genos = num_genos+1
            cell_line_ID = lcl

    return [num_genos, cell_line_ID]
########################################################################
#
# for nanostring, check based on the input 'genotype', as with above function,
# that only one lcl has the deletion (i.e. has genotype [1,0],[0,1], or [1,1]) 
#- what does genotype .,. mean, and should LowQual be discarded?
#
####################################################################
def onlynumgenos(genotype, NUM = 1):
    onlyNUM = 0
    cellline = []
    for lcl in range(0,len(genotype)):
        if genotype[lcl] == [1,1] or genotype[lcl] == [0,1] or genotype[lcl]==[1,0]:
            onlyNUM = onlyNUM +1
            cellline.append(lcl)
            if onlyNUM > NUM:
                return [-1]
    if len(cellline) < NUM:
        return [-1]
    else:
        return cellline

def onlytwo(genotypes):
    genotype = genotypes[0]
    onlytwo = 0
    celllines = []
    for lcl in genotype:
        if lcl == [1,1] or lcl == [0,1] or lcl==[1,0]:
            if onlytwo == 2:
                return -1
            else:
                onlytwo = onlytwo+1
                celllines.append(genotype.index(lcl))
    return celllines
##################################################################
# Oct 11, 2010
# check that indels are <=50 bp 
#
#
#################################################################

#
#########################################################
#
# June 1 : group cell lines by the genotype that they have
#   
# input: genos: ([1,1,0,1],[1,0,1,1]...)    
##########################################################
def group_lcls(genos, lines, lcls_unphased):
    truelines = []
    lcl_groups =[]
    for l in range(0,len(lines)):
        truelines.append(lines[l].split('\t')[9:])
    for j in range(0,len(genos)):        
        lcl_groups.append([])
        for cells in range(0, len(truelines[0])):
            for allele in (0,2):
                hap = 1
                for snp in range(0, len(truelines)):
                    if int(truelines[snp][cells][allele]) == genos[j][snp]:
                        pass
                    else:
                        hap = 0
                if hap == 1 and cells not in lcls_unphased:
                    lcl_groups[-1].append(cells)

    return lcl_groups
#for i in range(0,len(lines)):
#       for j in 

def max_overlap(lcl_group1, lcl_group2):
    overlaps = []
    for item in map(lambda x: map(lambda y:list(set(x) & set(y)), lcl_group1),lcl_group2):
        overlaps.extend(item)
    return max(map(lambda x: len(x), overlaps))


# get rid of redundancy (regions that are subsets of other ones)            
def get_rid_of_redundancy(regions, genomic_positions):
    temp = []
    temp.extend(regions)
    temp_gp = []
    temp_gp.extend(genomic_positions)
    for i in range(1, len(regions)):
        if genomic_positions[i][1] == genomic_positions[i-1][1] or genomic_positions[i][2] == genomic_positions[i-1][2] or (genomic_positions[i][1] > genomic_positions[i-1][1] and genomic_positions[i][1] < genomic_positions[i-1][2]) or (genomic_positions[i-1][1] > genomic_positions[i][1] and genomic_positions[i-1][1] < genomic_positions[i][2]) :
            if abs(genomic_positions[i][1] - genomic_positions[i][2]) > abs(genomic_positions[i-1][1] - genomic_positions[i-1][2]):
                if regions[i-1] in temp and genomic_positions[i-1] in temp_gp:
                    temp.remove(regions[i-1])
                    temp_gp.remove(genomic_positions[i-1])
            else:
                if regions[i] in temp and genomic_positions[i] in temp_gp:
                    temp.remove(regions[i])
                    temp_gp.remove(genomic_positions[i])

    regions = temp
    genomic_positions = temp_gp
    return [genomic_positions, regions]




























