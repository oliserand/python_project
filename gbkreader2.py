import sys, os, re
from urllib import request


'''
The menu/ output
GBK Reader ­ (file.gbk)DNA  
====================================================================== 
Sequence: U49845 (SCU49845, 5028bp) 
Description: Saccharomyces cerevisiae (baker's yeast) 
Source: Saccharomyces cerevisiae (baker's yeast) 
Number of References: 3 
Number of Features: 6 
====================================================================== 
R:References S:Sequence M:Motif T:Translate F:Features E:Export Q:Quit
>
'''

#Sections
#R:References S:Sequence M:Motif T:Translate F:Features E:Export Q:Quit

#Codon table
codonTable = {}
url = 'https://popoolation.googlecode.com/svn-history/\
        r179/trunk/syn-nonsyn/codon-table.txt'
if not os.path.exists('./codon-table.txt'):
    webpage = request.urlopen(url)
    webpage = webpage.read().decode('utf8')
    for line in webpage.split('\n'):
        if not line.startswith('#'):
            codon, aa = line.split(':')
            codonTable[codon] = aa.lstrip()
else:
    page_handle = open('codon-table.txt','r')
    page = page_handle.readlines()
    for line in page:
        if not line.startswith('#'):
            line = line.rstrip()
            codon, aa = line.split(':')
            codonTable[codon] = aa.lstrip()

def multiLines(pattern, lineIdx, line, level1, level2):
    ref = {}
    space = ' '
    #Takes a string (REFERENCE or FEATURE) and returns a dict
    #Checking the first level. If it is level1, add it as a key
    if line.startswith(pattern):
        counter1 = 1
        while (lines[lineIdx+counter1].startswith(level1[0]*space) or lines[lineIdx+counter1].startswith(level1[1]*space)) and (not lines[lineIdx+counter1].startswith(level2*space)):
            #Extract the key
            nextLineLvl1 = lines[lineIdx+counter1]
            #The key is added to the REFERENCE or FEATURE in the genbank object
            tmpKey = nextLineLvl1[:12].strip()
            tmpVal = nextLineLvl1[12:].rstrip()
            counter1 += 1
            counter2 = 0
            ref[tmpKey] = tmpVal
            while lines[lineIdx+counter1+counter2].startswith(level2*space):
                newTmpVal = lines[lineIdx+counter1+counter2].strip()
                tmpVal += space + newTmpVal
                ref[tmpKey] = tmpVal
                counter2 += 1
    #The ref dictionary contains all pattern (i.e pattern, or feature) information 
    #Can be added to the list in the reference object
 #   print(ref)
    record.references.append(ref)
            


def seqHandler():
    pass

def motifHandler():
    pass

class GBRecord:
    #The genbank record container
    def __init__(self):
        self.locus = None
        self.source= None
        self.description = None
        self.references = []
        self.features = {}
        self.sequence = ''
        
    


records = []
filename = 'sequence.gb'
with open(filename, 'r') as handle:
    lines = list(map(lambda x: x.rstrip(), handle.readlines()))
    record = GBRecord()
    for lineIdx, line in enumerate(lines):
        if line.startswith('LOCUS'):
            tmpLine = line.split()
            record.locus = tmpLine[1]
            record.length = ''.join(tmpLine[2:4])

        if line.startswith('ACCESSION'):
            record.accession = line.rstrip().split()[1]
        
        if line.startswith('SOURCE'):
            record.source = line[12:].rstrip()
            record.description = record.source

        if line.startswith('REFERENCE'):
#            print(line)
            multiLines('REFERENCE', lineIdx, line, [2,3], 12)

        if line.startswith(''):
            pass
        if line.startswith(''):
            pass
        #Read ahead to get multilines
    
    
print(record.references)
