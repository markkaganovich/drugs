file_reorder = open('./finalreorderlistDec2011.txt')
file_firstorder = open('./FirstLCLorderTotal.txt')

file_allpanel = open('../Dropbox/phase1_1KG_all_panel.txt')

alreadyhave = []

lines = file_reorder.readlines()
reorder_lines = str(lines[0]).split('\r')

lines = file_firstorder.readlines()
firstorder_lines = str(lines[0]).split('\r')

for l in firstorder_lines:
    if l not in reorder_lines:
        alreadyhave.append(l)

alreadyhavenum = [x[2:] for x in alreadyhave]
firstordernum = [x[2:] for x in firstorder_lines]
reordernum = [x[2:] for x in reorder_lines]

lines = file_allpanel.readlines()
allpanelnum = []
for l in lines:
    allpanelnum.append(l.split('\t')[0][2:])

order = []

for a in allpanelnum:
    if a not in alreadyhavenum:
        order.append(a)

cellorder = []
for i in range(0,len(order)):
    if i <411:
        cellorder.append('HG' + order[i])
    else:
        cellorder.append('GM' + order[i])

file = open('Dec2011cellorder.txt','w')
for o in cellorder:
    file.write(o + '\n') 

file.close()
