def checkilluminasnps(probesfile, probetype = 1, readsfile =
'110617_MONK_0184_B81N3EABXX_L4_eland_extended_pf.txt'):
    import simplejson

    file = open(probesfile)
    probes = simplejson.load(file)
    file.close()
    vars = [0] * len(probes)
    wts = [0] * len(probes)
    others = [0] * len(probes)
    varinds = [[]] * len(probes)

    reads = open('./'+readsfile)
    lines = reads.readlines(1000000)

    while(lines != []):
        for i in range(0,len(probes)):
            probeseq = getsnpregionseq(probes[i])
            probeseq = probeseq.upper()
            for l in lines:
                read = l.split('\t')[1]
                mismatches = len(probeseq) - len(filter(lambda x: probeseq[x] ==
read[x], range(0,len(probeseq))))
                if probeseq[29] == read[29] and mismatches <= 2:
                    wts[i] = wts[i]+1
                if probeseq[29] != read[29] and mismatches <= 3:
                    vars[i] = vars[i] +1
                    varinds[i].append(l)
                if mismatches <= 2:
                 others[i] = others[i] + 1
        lines = reads.readlines(1000000)

    return [wts, vars, varinds, others]

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

