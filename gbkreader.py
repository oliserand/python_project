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


globList = []
def multiline(curPos, pattern='REFERENCE', level1=2, level2=3, nspace=12):
    #Does the multiline stuffs.. Returns a dict for the current pattern
    currDict = {pattern:{}}
    #Given current line (level0) -> find the next line
    if line.startswith(pattern):
        #Get the first key (for level0)
        data = line.split()
        key = data[0]
        #Adding the reference number if entry if reference
        if 'REFERENCE' in line:
            refnum = data[1]
            currDict[pattern]['REFNUM'] = refnum
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
                if '/' in tmpVal:
                    tmpVal = tmpVal.split('/')
                    #Add a map strip operation here
                    tmpVal = list(map(lambda x: x.strip(), tmpVal))
                    #Removing extraneous spaces from translation (temp fix?)
                    for j, k in enumerate(tmpVal):
                        if 'translation' in k:
                            k = k.replace(' ', '')
                            tmpVal[j] = k
                
                currDict[pattern][tmpKey] = tmpVal
            counter1 += 1
    #Record position of dict. Maybe add it to the object somewhere around here?
    globList.append(currDict)
    
def refHandler():
    '''Chek possible redundant additions'''
    numRef = 0
    authors = []
    titles = []
    journal = ''
    for dictionary in globList:
        if list(dictionary.keys())[0] == 'REFERENCE':
            numRef += 1
            authors.append(dictionary['REFERENCE']['AUTHORS'])
            titles.append(dictionary['REFERENCE']['TITLE'])
            journal += dictionary['REFERENCE']['JOURNAL']

        #print('test',dictionary)
    print('There are', numRef, 'articles reported for the sequence', accession)
    for idx, author in enumerate(authors):
        print('['+str(idx+1)+']', author)
    refDetails = input('Input the number of a reference for details (M for the Menu) :')
    if refDetails == 'M':
        return 'continue'
    else:
        refDetails = int(refDetails)
        print('Title:\n\t'+titles[refDetails-1])
        print('Authors')
        currAuthors = authors[refDetails-1].replace('.,', '.').replace(' and ', ' ').split(' ')
        for auth in currAuthors:
            print('\t'+auth)
        print('Journal:\n\t'+journal)
        mprompt = input('Input the number of a reference for details (M for the Menu) :')
        if mprompt == 'M':
            return 'continue'
    return None

class GenBank:
    #Container for genbank record
    def __init__(self, sequence, references, features, accession, length):
        self.sequence = sequence
        self.reference = reference
        self.features = features
        self.accession = accession
        self.length = length


    def translateSeq(self, frame=0):
        #nucleotide self.sequence -> protein
        self.sequence = self.sequence[frame:]
        seqLen = len(self.sequence)
        protein = ''
        codons = [self.sequence[i:i+3] for i in range(0, seqLen-3+1, 3)]
        for codon in codons:
            aa = codonTable[codon]
            #Don't print the stop codon
            if aa == '-':
                break 
            else:
                protein += aa 

        return protein

    

def seqHandler():
    rangeDisp = input('Range: ')
    if rangeDisp == 'M':
        return 'continue'
    else:
        leftParen = rangeDisp[0]
        rightParen = rangeDisp[-1]
        rangeDisp = rangeDisp[1:-1].split(',')
        lowerLim = int(rangeDisp[0])
        upperLim = int(rangeDisp[1]) 
        '''To correct limits'''
        if leftParen == '(':
            #Exclude lim
            lowerLim += 1
        elif leftParen == '[':
            #Include lim
            lowerLim -= 1
        if rightParen == ')':
            #Exclude lim
            upperLim -= 1
        elif rightParen == ']':
            #Include lim
            pass
        #Build the query
        for dictionary in globList:
            if list(dictionary.keys())[0] == 'SEQUENCE':
                #String formatting
                currSeq = dictionary['SEQUENCE'][lowerLim:upperLim].upper()
                for i in range(0,len(currSeq), 60):
                    print(currSeq[i:i+60])
        return None
 

def motifHandler():
    pass



filename = 'sequence.gb'
#Parsing the genbank file
with open(filename) as handle:
    #Try to consider multiple records, in a later dict implementation
    lines = handle.readlines()
    recordLen = len(lines)
    readSeq = False
    sequence = ''
    running = True
    while running:
        #The main loop
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
                multiline(lineIdx, pattern='REFERENCE')

            elif line.startswith('FEATURES'):
                multiline(lineIdx, pattern='FEATURES', level1=5, level2=5, nspace=21)
                
            elif line.startswith('ORIGIN'):
                readSeq = True
            if (readSeq == True) and ('ORIGIN' not in line) and lineIdx <= recordLen: #Careful of the record separator
                tmpSeq = line.strip().split()
                tmpSeq = tmpSeq[1:]
                tmpSeq = ''.join(tmpSeq)
                sequence += tmpSeq
        globList.append({'SEQUENCE':sequence})
        readSeq = False #Stop adding to string, in case of multiple records
        sequence = '' #Re-initialise sequence, in case of multiple records

      
        line1 = 'GBK Reader - (' + filename +') DNA'
        line2 = '='*70
        line3 = 'Sequence: ' + accession + ' (' + locus + ', ' + length + ')'
        line4 = 'Description: ' + description
        line5 = 'Source: ' + source
        line6 = 'Number of References', None
        line7 = 'Number of Features', None
        line8 = '='*70

        print(line1)
        print(line2)
        print(line3)
        print(line4)
        print(line5)
        print(line6)
        print(line7)
        print(line8)
        print()
        prompt = input('R:References S:Sequence M:Motif T:Translate F:Features E:Export Q: Quit\n>')
        if prompt == 'Q':
            running = False
        elif prompt == 'R':
            #Do reference stuff
            '''Chek possible redundant additions'''
            if refHandler() == 'continue':
                continue

        elif prompt == 'S':
            #Do sequence stuff
            if seqHandler() == 'continue':
                continue
        
        elif prompt == 'M':
            #Do motif stuff
            motifPrompt = input('Motif:')
            while len(motifPrompt) < 5:
                motifPrompt = input('Motif:')
            
            print('Searching for the motif', motifPrompt.lower()+'...')
            '''Not finished'''
            pass
        
        elif prompt == 'T':
            #Do translation routine
            pass
        elif prompt == 'F':
            #Do features stuff
            pass
        elif prompt == 'E':
            #Export to fasta
            pass
        elif prompt == '':
            print(globList)
#
#
#References
#Sequence
#Motif
#Y=Translate
#Features
#Export
#Quit
