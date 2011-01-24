import simplejson

def check_combined(regions, indeces, num_of_lcls = 59):
    combo = []
    map(lambda x: map(lambda y: combo.extend(y),regions[x]) ,indeces)

    combo = list(set(combo))

    if len(combo)>= num_of_lcls:
        return True 
    else:
        return False 

######################################################################
# returns lcl_id, len: # of cell lines, each one is a list of [x,y] where 
# x is region group  index, y is which haplotype in the region group
#
# the specificity check asks if each of the cell lines is determined by
# a unique combination of haplotypes (of a region)
#
# the input is an array of regions, and indeces which is the group of
# regions upon which the test is being performed
#
# the output is the number of groups that this set of regions splits the
# cell lines into.... this needs to be the same as the number of cell lines
#
#######################################################################

#######################################################################
# for each cell line, return the regions i and the specific group defined
# within the region that includes that cell line (there will likely be multiple different
# ones) -- there will be at least two because there are two alleles
######################################################################       
def check_specificity(regions, indeces, num_of_lcls = 59):

    lcl_id=[]
    for lcl in range(0,num_of_lcls):
        lcl_id.append([])
        for i in indeces:
            for j in range(0,len(regions[i])):
                if lcl in regions[i][j]:
                    lcl_id[-1].append([i,j])
   
    a = len(check_specificity_helper1(lcl_id))
    if a == num_of_lcls:
        return True
    else: 
        if a > num_of_lcls:
            print "ERROR: wrong num of lcls"
        return False


########################################################################
# the helper function asks which lcls have a unique region+haplotype signature
##########################################################################
def check_specificity_helper1(lcl_id):

    unique_lcls = []
    for i in range(0,len(lcl_id)):
        for hap in lcl_id[i]:
            if is_unique(lcl_id,i, hap):
                unique_lcls.append(i)
                break
    return unique_lcls

def is_unique(lcl_id, index, hap):

    for j in range(0,len(lcl_id)):
        if not j==index:
            if hap in lcl_id[j]:
                return False

    return True

#########################################################################################
# read in regions, and get rid of duplicates, i.e. subsets of the regions (take the longest
# of the subsets)
# keep track of changes genomic_positions
# return regions - list of regions found, and each element is a list of cell line groups
    # that the region groups into
##########################################################################################
def read_regions(filename = './CEUallgenos_output'):
    file = open(filename)
    [genomic_positions, regions] = simplejson.load(file)

    return [genomic_positions, regions]


'''
find which region combinations satistfy the specificity and complete coverage of all cell lines criteria

return the indeces that indicate the best region combinations, the indeces correspond to the genomic_positions 
and regions lists returned by above function read_regions
'''
def group_barcodes(regions, num_lcls = 59):
    barcode_groups =[]

    #initial positions    
    start1 = 0
    start2 = 0
    end = 0

    while(start1 < len(regions)):
        start2 = start1
        while(start2 < len(regions)):
            start2 = start2 +1
            end = start2
            while(end < len(regions)):
                end = end+1
                indeces = [start1]+range(start2, end)
                if check_combined(regions, indeces):
                    if check_specificity(regions, indeces):
                        barcode_groups.append(indeces)
                        break
        start1=start1+1

# get rid of subsets
    temp = []
    temp.extend(barcode_groups)
    for i in range(1, len(barcode_groups)):
        if abs(len(barcode_groups[i]) - len(barcode_groups[i-1])) == 1:
            if len(list(set(barcode_groups[i]+barcode_groups[i-1]))) == max(len(barcode_groups[i]), len(barcode_groups[i-1])):
                if len(barcode_groups[i]) > len(barcode_groups[i-1]):
                    temp.remove(barcode_groups[i])
                else:
                    temp.remove(barcode_groups[i-1])

    barcode_groups = temp                

    return barcode_groups

