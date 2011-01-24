
from xlrd import *

File = open('celllines.txt', 'w')

book = open_workbook('./1000Samples_6Oct2009.xls')
sheet = book.sheet_by_index(0)
cellnum = 0
for i in range(0, sheet.nrows):
    if sheet.row_values(i)[2] != '' and sheet.row_values(i)[10] != '' and sheet.row_values(i)[2] != 'Sample\nID':
        cellnum = cellnum+1
        File.write(str(cellnum) + '   ')
        File.write('GM' + str(sheet.row_values(i)[2])[2:len(str(sheet.row_values(i)[2]))])
        File.write('\n')

File.write('\n'+'\n')

for i in range(29, sheet.nrows):
    if sheet.row_values(i)[2] != '' and sheet.row_values(i)[15] != '':
        cellnum = cellnum+1
        File.write(str(cellnum) + '   ')
        File.write('GM' + str(sheet.row_values(i)[2])[2:len(str(sheet.row_values(i)[2]))] + '     ' + str(sheet.row_values(i)[15]))
        File.write('\n')


