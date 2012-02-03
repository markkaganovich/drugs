def getvcfmatrix(filename):
    file = open(filename)
    outputfile = open(filename+'matrixoutput', 'w')
    total = 1
    lines = file.readlines(1000000)
    while(lines!=[]):
        for line in lines:
            if not line.startswith('#'):
                total=+1
                l = line.split('\t')
                for g in l[9:]:
                    matrixentry = int(g[0]) + int(g[2])
                    outputfile.write(str(matrixentry) +'\t')
                outputfile.write('\n')
    print total
    file.close()
    outputfile.close()


