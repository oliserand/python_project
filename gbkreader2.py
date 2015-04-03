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

def getRefs(pattern, lineIdx, line, level1, level2):
    ref = {}
    space = ' '
    #Takes a string (REFERENCE) and returns a dict containing all pattern information
    #Checking the first level. If it is level1, add it as a key
    counter1 = 1
    if line.startswith(pattern):
        while lines[lineIdx+counter1].startswith(level1*space):
            #Extract the key
            nextLineLvl1 = lines[lineIdx+counter1]
            #The key is added to the REFERENCE or FEATURE in the genbank object
            tmpKey = nextLineLvl1[:12].strip()
            tmpVal = nextLineLvl1[level2:].rstrip()
            if tmpKey != '':
                ref[tmpKey] = tmpVal
            counter2 = 1
            while lines[lineIdx+counter1+counter2].startswith(level2*space):
                newTmpVal = lines[lineIdx+counter1+counter2].strip()
                tmpVal += (space + newTmpVal)
                ref[tmpKey] = tmpVal
                counter2 += 1
            counter1 += 1
        #Appending to the list in the reference object
        record.references.append(ref)
            

def getFeats(pattern, lineIdx, line, level1, level2):
    space = ' '
    #Takes a string (FEATURE) and returns a dict containing all pattern information
    #Checking the first level. If it is level1, add it as a key
    counter1 = 1
    if line.startswith(pattern):
        while lines[lineIdx+counter1].startswith(level1*space):
            feature = {}
            #Extract the key
            nextLineLvl1 = lines[lineIdx+counter1]
            #The key is added to the FEATURE in the genbank object
            tmpKey = nextLineLvl1[:12].strip()
            tmpVal = nextLineLvl1[level2:].rstrip()
            if tmpKey != '':
                feature[tmpKey] = tmpVal
            counter2 = 1
            while lines[lineIdx+counter1+counter2].startswith(level2*space) and tmpKey != '':
                newTmpVal = lines[lineIdx+counter1+counter2].strip()
                tmpVal += (space + newTmpVal)
                feature[tmpKey] = tmpVal
                counter2 += 1
            counter1 += 1
            #Somehow empty dicts were also added. Braching avoids that
            if feature != {}:
                #Cleaning up some fields
                if 'translation' in feature[tmpKey]: feature[tmpKey] = feature[tmpKey].replace(' ', '')
                if '/' in feature[tmpKey]: feature[tmpKey] = feature[tmpKey].split('/')
                record.features.append(feature)
            

def getSeq(lineIdx):
    #Assigns the sequence to the record
    counter1 = 1
    sequence = ''
    while lines[lineIdx+counter1].split()[0].strip().isnumeric():
        tmpSeq = lines[lineIdx+counter1].split()[1:]
        tmpSeq = ''.join(tmpSeq)
        sequence += tmpSeq
        counter1 += 1
    #Add sequence to record
    record.sequence = sequence

def gb2fasta(record):
    #Converts genbank to fasta format
    sequence = record.sequence
    formattedSeq = ''
    for i in range(0, len(sequence), 80):
        formattedSeq += sequence[i:i+80]+'\n'
    '''Do I have to get the right description? for the header'''
    header = '>'+record.description+'\n' 
    
    destfile = input('Filename :')
    while os.path.exists(destfile):
        print('File exists. Provide a different name')
        destfile = input('Filename :')
    with open(destfile, 'w') as fasta_handle:
        input('File ('+destfile+') created.')
        fasta_handle.writelines(header)
        fasta_handle.writelines(formattedSeq)

def dispMotif():
    '''To add Camel case for regex matches to *...'''
    tmpSeq = record.sequence
    motif = input('Motif:')
    while motif != '':
        while len(motif) < 5 and motif != '':
            motif = input('Motif:')
        if motif != '':
            print('Searching for the motif', motif+'...' )
            motifC = re.compile(motif)
            match = motifC.search(tmpSeq)
            tmpSeq2 = tmpSeq[match.start()-5:match.end()+5]
            print(tmpSeq2[:5] + '['+str(match.start())+']' + tmpSeq[match.start():match.end()]+ '['+str(match.end())+']' + tmpSeq2[-5:] )
            motif = input('Motif:')

def dispRef():
    #Displays reference details
    print('There are', len(record.references), 'articles reported for the sequence', record.accession)
    for refIdx, ref in enumerate(record.references):
        print('['+str(refIdx+1)+']', ref['AUTHORS'])
    refNum = input('Input the number of a reference for details (M for the Menu):')
    while refNum != 'M':
        refNum = int(refNum) - 1
        print('Title:')
        print('\t'+record.references[refNum]['TITLE'])
        print('Authors:')
        print('\t'+record.references[refNum]['AUTHORS'])
        print('Journal:')
        print('\t'+record.references[refNum]['JOURNAL'])
        refNum = input('Input the number of a reference for details (M for the Menu):')

def dispSeq():
    #Displays sequence ranges
    seqRange = input('Range: ')
    while seqRange != 'M' and seqRange != '':
        leftParen = seqRange[0]
        rightParen = seqRange[-1]
        seqRange = seqRange[1:-1].split(',')
        lowerLim = int(seqRange[0])
        upperLim = int(seqRange[1]) 
        
        if leftParen == '(':
            #Exclude lim
            pass
        elif leftParen == '[':
            #Include lim
            lowerLim -= 1
        if rightParen == ')':
            #Exclude lim
            pass
        elif rightParen == ']':
            #Include lim
            upperLim += 1
        #Build query
        currSeq = record.sequence[lowerLim:upperLim].upper()
        for i in range(0,len(currSeq), 60):
            print(currSeq[i:i+60])
        seqRange = input('Range: ')

def dispFeat():
    #Displays features
    def countFeat(featName):
        #Returns feature count
        totCount = 0
        for featDict in record.features:
            if featName in featDict.keys():
                totCount += 1
        return totCount
   '''To continue here... Task isolate ranges from tmpKey, by cleaning numbers from strings''' 
    def rangeHandler():
        for featDict in record.features:
            for feat in featDict.keys():
                tmpKey = featDict[feat][0].split('..')
                print(feat, tmpKey)

    featPrompt = input('Choose the type of feature query P(Position) or N(Name) :')
    while featPrompt != '':
        if featPrompt == 'P':
            #Check the first elements of the lists, mapped to the features in every dict
            rangeHandler()

        elif featPrompt == 'N':
            featName = input('Name:')
            countFeat(featName)
            print('Searching for features with name', featName+'...')
            featCount = 0
            for featDict in record.features:
                if featName in featDict.keys():
                    featCount += 1
                    tmpFeat = list(map(lambda x: '\t/'+x, featDict[featName]))
                    print('Feature', featCount, 'of', countFeat(featName))
                    print(featName, '\n'.join(tmpFeat))
        featPrompt = input('Choose the type of feature query P(Position) or N(Name) :')


def dispMenu():
    #Displays record summary
    print('GBK Reader - ('+filename+')'+'DNA')  
    print('======================================================================')
    print('Sequence:', record.accession, '('+record.locus+', '+record.length+')') 
    print('Description:', record.description) 
    print('Source:', record.source)
    print('Number of References:', len(record.references)) 
    print('Number of Features:', len(record.features)) 
    print('======================================================================')
    
    #Displays menu
    prompt = input('R:References S:Sequence M:Motif T:Translate F:Features E:Export Q:Quit\n>')
    if prompt == 'Q':
        confirm = input('Do you want to exit(E) or do you want load another file(F)')
        if confirm == 'E':
            global finished
            finished = True
    else:
        if prompt == 'R':
            #Reference handling
            dispRef()

        elif prompt == 'S':
            #Sequence handling
            dispSeq()
        
        elif prompt == 'M':
            #Motif handling
            dispMotif()
        
        elif prompt == 'T':
            #Translation
            pass
        
        elif prompt == 'F':
            #Feature handling
            dispFeat()
        
        elif prompt == 'E':
            #Export to fasta
            gb2fasta(record)


class GBRecord:
    #The genbank record container
    def __init__(self):
        self.accession = None
        self.locus = None
        self.length = None
        self.source= None
        self.description = None
        self.references = []
        self.features = []
        self.sequence = None

    def translate(self, frame=0):
        #Translates nucleotide -> protein, in different frames
        self.peptide = ''
        self.frame = frame
        self.tmpSeq = self.sequence[frame:]
        for i in range(0, len(self.tmpSeq)-3 + 1, 3):
            codon = self.sequence[i:i+3]
            print(codon)


finished = False
while not finished:
    #The container for all records
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
                getRefs('REFERENCE', lineIdx, line, 2, 12)

            if line.startswith('FEATURES'):
                getFeats('FEATURES', lineIdx, line, 5, 21)
                
            if line.startswith('ORIGIN'):
                getSeq(lineIdx)

    #Queries
    dispMenu() 
    
#print(record)
