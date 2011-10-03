import simplejson

def chrlinessites(sitesfile):
    file =       open(sitesfile)
    outputfile = open(str(sitesfile + 'chrlines'), 'write')
    lines =      file.readlines()
    file.close()
    chrlineshash = {}
    for i,l in enumerate(lines):
        p =   l.split('\t')
        chr = p[0]
        if chr not in chrlineshash.keys() and not chr.startswith('#'):
            chrlineshash[chr] = (i,)
            if str(int(chr)-1) in chrlineshash.keys():
                chrlineshash[str(int(chr)-1)] = (chrlineshash[chr][0],i)
            else:
                continue;
        else: 
            continue;
    lastchr = chrlineshash.keys()[-1]
    chrlineshash[lastchr] = (chrlineshash[lastchr][0],len(lines))
    simplejson.dump(chrlineshash, outputfile)
    outputfile.close()

