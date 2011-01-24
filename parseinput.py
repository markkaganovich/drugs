

def parseinput(inputfile, outputname):
    file = open(inputfile)
    
    line =1 
    chr = 1
    output = open(str(outputname)+str(chr), 'w')
    while(line):
        line = file.readline()
        if line =='':
            break;
        print line.split('\t')[0]
        print str(outputname)+str(chr)
        if not str(line.split('\t')[0][0]) == '#':
            chrom = int(line.split('\t')[0])
            if chr == chrom:
                output.write(line)
            else:
                chr = chrom
                output.close()
                output = open(str(outputname)+str(chr), 'w')
                output.write(line)



