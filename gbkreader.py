import sys, os, re
from urllib import request


#Codon table
codonTable = {}
url = 'https://popoolation.googlecode.com/svn-history/r179/trunk/syn-nonsyn/codon-table.txt'
if not os.path.exists('./codon-table.txt'):
    webpage = request.urlopen(url)
    webpage = webpage.read().decode('utf8')
    for line in webpage.split('\n'):
        if not line.startswith('#'):
            codon, aa = line.split(':')
            codonTable[codon] = aa.lstrip()
    #Also getting a local copy of the file
    request.urlretrieve(url, filename='./codon-table.txt')
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

def getDef(lineIdx, line):
    #Get the definition for the fasta defline
    counter = 1
    spaces = ' '*12
    output = line[12:]
    while lines[lineIdx+counter].startswith(spaces):
        output += lines[lineIdx+counter][12:]
        counter += 1
    return output
        
def gb2fasta(record):
    #Converts genbank to fasta format
    sequence = record.sequence
    formattedSeq = ''
    for i in range(0, len(sequence), 80):
        formattedSeq += sequence[i:i+80]+'\n'
    header = '>'+record.definition+'\n' 
    
    destfile = input('Filename :')
    while os.path.exists(destfile):
        print('File exists. Provide a different name')
        destfile = input('Filename :')
    with open(destfile, 'w') as fasta_handle:
        input('File ('+destfile+') created.')
        fasta_handle.writelines(header)
        fasta_handle.writelines(formattedSeq)

def dispMotif():

    def formatOutput(motif, match):
        #takes an motif (re object) and its equivalent string
        #split the string
        motif  = motif.replace('?',' ').replace('*', ' ').split(' ')
        groups = match.groups()
        output = ''
        try:
            for i in range(len(motif)):
                output += (motif[i] + groups[i].upper())
        except IndexError:
            #If there are fewer non-wildcard elements
            output += motif[i]
        return output
        
    tmpSeq = record.sequence
    #Request for motif. Make search case sensitive
    motif = input('Motif:')
    #Do not accept an empty string or a string containing only non-alphabet chars
    while motif != '' and (any(map(lambda x:x.isalpha(), motif))):
        #If motif is too short or empty, keep on asking for motif
        while len(motif) < 5 and motif != '':
            motif = input('Motif:')
        motif = motif.lower()
        if motif != '':
            print('Searching for the motif', motif+'...' )
            #Replacing wildcards by re wildcards
            motifReformatted = motif.replace('?','(.)').replace('*','(.*)')
            motifObj = re.compile(motifReformatted)
            match = motifObj.search(tmpSeq)
            tmpSeq2 = tmpSeq[match.start()-5:match.end()+5]
            print(tmpSeq2[:5] + '['+str(match.start())+']' + formatOutput(motif, match) + '['+str(match.end())+']' + tmpSeq2[-5:])
            motif = input('Motif:')

def dispRef():
    #Displays reference details
    print('There are', len(record.references), 'articles reported for the sequence', record.accession)
    for refIdx, ref in enumerate(record.references):
        print('['+str(refIdx+1)+']', ref['AUTHORS'])
    refNum = input('Input the number of a reference for details (M for the Menu):')
    while refNum != 'M' and refNum != '':
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
        
        if leftParen == '[':
            #Include lim
            lowerLim -= 1
        if rightParen == ']':
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
    
    def getLimitsAndTotal(givenRange, featureCounter, total, output):
        #Extracts the positions from the features and prepares total number of matching features
        for featDict in record.features:
            for feat in featDict.keys():
                #Every feature (feat) in featDict has its span as first element of a list
                tmpRange = featDict[feat][0]
                pattern = re.compile('(\d+)\.\.[<>]?(\d+)')
                match = pattern.search(tmpRange)
                #Limits are extracted
                lowerLim = int(match.groups()[0])
                upperLim = int(match.groups()[1])
                if (lowerLim >= givenRange[0]) and (upperLim <= givenRange[1]):
                    total.append(1)
                    output.append((feat,featDict[feat]))

    def rangeHandler():
        rangePrompt = input('Range:')
        givenRange = list(map(int, rangePrompt[1:-1].split(',')))
        leftBracket = rangePrompt[0]
        rightBracket = rangePrompt[-1]
        lowerLim = rangePrompt[1:-1].split(',')[0]
        upperLim = rangePrompt[1:-1].split(',')[1]
        print('Searching for features in the range', leftBracket+lowerLim+'-'+upperLim+rightBracket+'...')
        
        featureCounter = 0
        total = []
        output = []
        getLimitsAndTotal(givenRange, featureCounter, total, output)
        total = sum(total)

        #Display the feature
        for j in output:
            featureCounter += 1
            print('Feature', featureCounter, 'of', total)
            rangeToDisplay = '('+j[1][0].replace('..', ',').replace('<','').replace('>','').rstrip()+'):'
            print('\t'+j[0] + rangeToDisplay)
            for k in j[1][1:]: print('\t\t/'+k)

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
                    #Make a list of the output instead of printing
                    print('Feature', featCount, 'of', countFeat(featName))
                    print(featName +'('+tmpFeat[0].replace('..',',').replace('/','').replace('<','').replace('>','').strip()+'):'+'\n'+ '\n'.join(tmpFeat[1:]))
        featPrompt = input('Choose the type of feature query P(Position) or N(Name) :')


def dispMenu():
    global filename
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
        elif confirm == 'F':
            filename = input('Enter filename: ')
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
            dispPeptide()
        
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
        self.definition = None
        self.references = []
        self.features = []
        self.sequence = None

def translate(sequence, lowerLim, upperLim, frame):
    '''Re-test translation functions in several cases.. Something seems odd'''
    #Translates nucleotide -> protein, in different frames
    sequence = sequence[lowerLim:upperLim].upper()    
    
    def formatSeq(seq):
        #Formats sequence to 60chars per line
        output = ''
        for i in range(0, len(seq), 60):
            output += (seq[i:i+60]+'\n')
        return output

    def outputPep(sequence, lowerLim, upperLim, frame):
        peptide = ''
        tmpSeq = sequence[frame:]
        if frame == 0: fdisp = '1st'
        elif frame == 1: fdisp = '2nd'
        elif frame == 2: fdisp = '3rd'
        
        #If not None, return peptide
        if pattern.search(tmpSeq):
            match = pattern.search(tmpSeq)
            ntStart = match.start()
            ntEnd = match.end()
            hit = tmpSeq[ntStart:ntEnd]
        
            for i in range(0, len(hit), 3):
                codon = hit[i:i+3]
                peptide += codonTable[codon].rstrip()

            print('Amino acids sequence from nucleotide', ntStart+lowerLim+1, 'to', ntEnd+lowerLim, 'in the', fdisp,'ORF')
            return formatSeq(peptide).replace('-','')

    #Taking ATG -> start; ATG/TAA/TGA -> stop
    pattern = re.compile('ATG([ATGC]{3})*?(TAG|TAA|TGA)')
    if type(frame) == int:
        print(outputPep(sequence, lowerLim, upperLim, frame))
            
    elif type(frame) == list:
        for i in frame:
            print(outputPep(sequence, lowerLim, upperLim, i))
            
def dispPeptide():
    #Displays ORF translations
    seqRange = input('Range: ')
    frame = input('ORF:')
    while seqRange != '':
        if seqRange == 'FULL':
            lowerLim = 0
            upperLim = len(record.sequence)
        else:
            leftParen = seqRange[0]
            rightParen = seqRange[-1]
            seqRange = seqRange[1:-1].split(',')
            lowerLim = int(seqRange[0])
            upperLim = int(seqRange[1]) 
            
            if leftParen == '[':
                #Include lim
                lowerLim -= 1
            
            if rightParen == ']':
                #Include lim
                upperLim += 1

        if frame == '1':
            peptide = translate(record.sequence, lowerLim, upperLim, frame=0)
        elif frame == '2':
            peptide = translate(record.sequence, lowerLim, upperLim, frame=1)
        elif frame == '3':
            peptide = translate(record.sequence, lowerLim, upperLim, frame=2)
        elif frame == '':
            frame = [0,1,2]
            peptide = translate(record.sequence, lowerLim, upperLim, frame=frame)

        seqRange = input('Range: ')
        frame = input('ORF:')

finished = False
filename = sys.argv[1]
while not finished:
    with open(filename, 'r') as handle:
        lines = list(map(lambda x: x.rstrip(), handle.readlines()))
        record = GBRecord()
        for lineIdx, line in enumerate(lines):
            if line.startswith('LOCUS'):
                tmpLine = line.split()
                record.locus = tmpLine[1]
                record.length = ''.join(tmpLine[2:4])
            if line.startswith('DEFINITION'):
                record.definition = getDef(lineIdx, line)
            
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

    #Menu with prompts
    dispMenu() 
