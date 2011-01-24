
file = open('CGEPgenes.txt')
betterfile = open('./CGEPgenesFormatted','w')
while(1):
    line = file.readline()
    if not line:
        break
    toks = line.split('\t')
    i=0
    newline=[]
    for tok in toks:
        if i < 6:
            newline.append(tok)
            print newline
            i=i+1
        else:
            for item in newline:
                betterfile.write(str(item) +'\t')
            betterfile.write('\n')
            newline=[]
            i=0



