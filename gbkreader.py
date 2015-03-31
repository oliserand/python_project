import sys, os
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

def translateSeq(sequence,frame=0):
    #Translates sequence
    sequence = sequence[frame:]
    seqLen = len(sequence)
    protein = ''
    codons = [sequence[i:i+3] for i in range(0, seqLen-3+1, 3)]
    for codon in codons:
        protein += codonTable[codon]
    return protein


filename = 'sequence.gb'
#Parsing the genbank file
with open(filename) as handle:
    gbData = {}
    #Try to consider multiple records, in a later dict implementation
    locusInfo = handle.readline().split()
    locus = locusInfo[1]
    length = locusInfo[2] + locusInfo[3]
    lines = handle.readlines()
    readSeq = False
    readFeatures = False
    readRefs = False
    sequence = ''
    for line in lines:
        #Extracting sequence details
        if line.startswith('ACCESSION'):
            accession = line.split()[1]
        elif line.startswith('SOURCE'):
            description = line.split()[1:]
            description = ' '.join(description)
            source = description
        elif line.startswith('REFERENCE'):
            readRefs = True
            ref = {} #Create a dict for the references
        elif line.startswith('FEATURES'):
            #The features
            features = {} #Create a dict for the features
            readFeatures = True
            pass
        if (readRefs == True) and ('REFERENCE' not in line):
                #Get Reference details and put them in a ref dict
                #With key = number
                pass
        if line.startswith('ORIGIN'):
            readSeq = True
        if (readSeq == True) and ('ORIGIN' not in line) and ('//' not in line):
            tmpSeq = line.strip().split()
            tmpSeq = tmpSeq[1:]
            tmpSeq = ''.join(tmpSeq)
            sequence += tmpSeq
        #Accession
        #Description (Source)
        #References
        #Features
    readFeatures = False #Stop adding to string, in case of multiple records
    readSeq = False #Stop adding to string, in case of multiple records
    readRefs = False
    sequence = '' #Re-initialise sequence, in case of multiple records
  
    line1 = 'GBK Reader - (' + filename +') DNA'
    line2 = '='*70
    line3 = 'Sequence: ' + accession + ' (' + locus + ', ' + length + ')'
    line4 = 'Description: ' + description
    line5 = 'Source: ' + source

    print(line1)
    print(line2)
    print(line3)
    print(line4)
    print(line5)


#References
#Sequence
#Motif
#Y=Translate
#Features
#Export
#Quit
