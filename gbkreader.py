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

def translateSeq(sequence,frame=0):
    #Translates sequence
    sequence = sequence[frame:]
    seqLen = len(sequence)
    protein = ''
    codons = [sequence[i:i+3] for i in range(0, seqLen-3+1, 3)]
    for codon in codons:
        protein += codonTable[codon]
    return protein

#def motifMatch(motif, seq):
#    #Takes a motif and a string seq. Returns a match
#    #Allowed characters: *, ?, ATCG, [, ], (, )
#    pattern = re.compile(motif)
#    match = pattern.search(seq)
#    if match != None:
#        start = match.start()
#        match.stop()
#        result = seq[start:stop]
#        #More work to add lower/ uppercase functionality
#        return result

globList = []
def multiline(curPos, pattern='REFERENCE', level1=2, level2=3, nspace=12):
    currDict = {pattern:[]}
    #Give current line (level0), finds the next line
    if line.startswith(pattern):
        #Get the first key (for level0)
        data = line.split()
        key = data[0]
        #Adding the reference number if entry if reference
        if 'REFERENCE' in line:
            refnum = data[1]
            currDict[pattern].append({'REFNUM': refnum})
        #While the next line starts with 2, 3 or 5 spaces (levels 1 and 2 are treated as one level)
        counter1 = 1
        while (lines[curPos+counter1].startswith(' '*level1)) or (lines[curPos+counter1].startswith(' '*level2)):
            currLine = lines[curPos+counter1]
            #Getting the second level
            tmpKey = currLine[:nspace].strip()
            tmpVal = currLine[nspace:].strip()
            #Looking ahead for entry continuity across lines
            counter2 = 1
            while lines[curPos+counter1+counter2].startswith(' '*nspace):
                tmpLine = ' '.join((lines[curPos+counter1+counter2].split()[0:]))
                tmpVal += ' ' + tmpLine #Otherwise space is truncated
                counter2 += 1
            #It seems that an empty key is added.. Removed with if statement
            if tmpKey != '':
                currDict[pattern].append({tmpKey: tmpVal}) 
            counter1 += 1
    #Record position of dict
    globList.append(currDict)
    return currDict
    

filename = 'sequence.gb'
#Parsing the genbank file
with open(filename) as handle:
    #Try to consider multiple records, in a later dict implementation
    lines = handle.readlines()
    readSeq = False
    readFeatures = False
    sequence = ''
    for lineIdx, line in enumerate(lines):
        #Extracting sequence details
        if line.startswith('LOCUS'):
            locusInfo = line.split()
            locus = locusInfo[1]
            length = locusInfo[2] + locusInfo[3]
        
        elif line.startswith('ACCESSION'):
            accession = line.split()[1]
        
        elif line.startswith('SOURCE'):
            description = line.split()[1:]
            description = ' '.join(description)
            source = description
        
        elif line.startswith('REFERENCE'):
            #multiline(lineIdx, pattern='REFERENCE')
            print(multiline(lineIdx, pattern='REFERENCE'))

        elif line.startswith('FEATURES'):
            readFeatures = True
            print(multiline(lineIdx, pattern='FEATURES', level1=5, level2=5, nspace=21))
            
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
