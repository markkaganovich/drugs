# Sept 10, 2010: Download Phase 2 aligned reads

import os
import simplejson as json
file = open('./pilot_data.alignment.index')
lines = file.readlines()
alreadydownloaded =[]
harddrivename = 'MarkHD'
jpositions = json.load(open('jpositions'))


for line in lines:
    name =line.split('\t')[0]
    if str(name.split('.')[3]) == 'SRP000031' and str(name.split('.')[5]) != 'unmapped':  #and str(name.split('.')[1]) == 'chrom14': 
        print name
#os.popen("curl ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/" + name +" -o ~/../../Volumes/" +harddrivename +"/" + name.split('/')[4])
        alreadydownloaded.append(name)
        for j in range(0,len(jpositions)):
            os.popen("~/libraries/samtools-0.1.8/samtools view -b ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/" + name + ' ' + str(14)+':' + str(jpositions[j][0]) + '-' + str(jpositions[j][1]) +"  > "+ 'J'+str(j)+'_'+name.split('/')[4] )
            print ("~/libraries/samtools-0.1.8/samtools view -b ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/" + name + ' ' +str(14)+':' + str(jpositions[j][0]) + '-' + str(jpositions[j][1]) +"  > "+ 'J'+str(j)+'_'+name )


                






