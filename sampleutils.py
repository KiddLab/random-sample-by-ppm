import sys
import random


#####################################################################
def init_blank_list(listLen,element):
    myList = []
    for i in range(listLen):
        myList.append(element)
    return myList
#####################################################################
##############################################################################
# Returns complement of a bp.  If not ACGT then return same char
def complement(c):
    if c == 'A':
        return 'T'
    if c == 'T':
        return 'A'
    if c == 'C':
        return 'G'
    if c == 'G':
        return 'C'
    if c == 'a':
        return 't'
    if c == 't':
        return 'a'
    if c == 'c':
        return 'g'
    if c == 'g':
        return 'c'
    # If not ACTGactg simply return same character
    return c   
##############################################################################
##############################################################################
# Returns the reverse compliment of sequence 
def revcomp(seq):
    c = ''    
    seq = seq[::-1] #reverse
    # Note, this good be greatly sped up using list operations
    seq = [complement(i) for i in seq]
    c = ''.join(seq)
    return c
##############################################################################
###############################################################################
def read_fasta_file_to_list(fastaFile):
    myDict = {}
    inFile = open(fastaFile,'r')
    line = inFile.readline()
    line = line.rstrip()
    if line[0] != '>':
        print 'ERROR, FILE DOESNNOT START WITH >'
        sys.exit()
    myName = line[1:]
    myDict[myName] = {}
    myDict[myName]['seq'] = ''
    myDict[myName]['seqLen'] = 0    
    mySeq = ''
    while True:
        line = inFile.readline()
        if line == '':
            myDict[myName]['seq'] = mySeq
            myDict[myName]['seqLen'] = len(myDict[myName]['seq'])         
            break
        line = line.rstrip()
        if line[0] == '>':
            myDict[myName]['seq'] = mySeq
            myDict[myName]['seqLen'] = len(myDict[myName]['seq'])         
            myName = line[1:]
            myDict[myName] = {}
            myDict[myName]['seq'] = ''
            myDict[myName]['seqLen'] = 0    
            mySeq = ''
            continue
        mySeq += line
    inFile.close()
    return myDict
###############################################################################
#setup the genome sequence information for doing random samples
# gets sequence in memory, and also list of cname, start,end,len
def initialize_genome_sequences(data):
    print 'Reading in genome fasta file from',data['genomeFasta']
    data['genomeSeqDict'] = read_fasta_file_to_list(data['genomeFasta'])
    tmpNameList = data['genomeSeqDict'].keys()
    tmpNameList.sort()
    # will also be now random from 0 to seqlen, with start
    startPos = 0
    data['chromInfo'] = []
    for i in range(len(tmpNameList)):
        seqLen = data['genomeSeqDict'][tmpNameList[i]]['seqLen']
        n = []
        n.append(tmpNameList[i])
        n.append(seqLen)
        n.append(startPos) # from previos
        endPos = startPos + seqLen -1
        n.append(endPos)
        startPos = endPos + 1
        data['chromInfo'].append(n)
    lastPos = endPos
    data['lastPos'] = lastPos
    print 'Have initialized sequences for sampling'
    for i in data['chromInfo']:
        print i
    print data['lastPos']
    
    if data['useExclusionRegions'] is True:
        initialize_exclusion_regions(data)
        
###############################################################################
# set up set of 0,1 vectors of exclusion regions for comparison
def initialize_exclusion_regions(data):
    print 'setting up exclusions!'

    data['exclusionLists'] = {}
    for i in data['chromInfo']:
        chromName = i[0]
        chromLen = i[1]
        data['exclusionLists'][chromName] = init_blank_list(chromLen,0)
    print 'Initialzed blank exclusion lists for %i chromosomes' % len(data['exclusionLists'])

    print 'Reading in exclusion data from',data['exclusionBedFile']
    inFile = open(data['exclusionBedFile'],'r')
    for line in inFile:
       line = line.rstrip()
       line = line.split()
       cName = line[0]
       b = int(line[1])
       e = int(line[2])
       #if not one of the chroms we are doing, just skip on over'
       if cName in data['exclusionLists']:
           for i in range(b,e):
               data['exclusionLists'][cName][i] = 1
    inFile.close()
    print 'Exclusion list setup!'


###############################################################################
#random.randint(a, b)
#Return a random integer N such that a <= N <= b.
# select random into from 0 to lastPos
# then conver it to proper chrom,pos (zero based)
def select_random_position(data):
    myInt = random.randint(0,data['lastPos'])
    (chrom,pos)= rand_int_to_chrom(myInt,data)
    randSel = {}
    randSel['chrom'] = chrom
    randSel['pos'] = pos # zero based
    randSel['strand'] = random.choice(['+','-'])
    return randSel
###############################################################################
# is zero based
def rand_int_to_chrom(myInt,data):
    for chrom in data['chromInfo']:
        if myInt >= chrom[2] and myInt <= chrom[3]:
            cName = chrom[0]
            cPos = myInt - chrom[2]
            return (cName,cPos)
###############################################################################
# setup ppm file
def initialize_ppm(data):
    data['ppmMat'] = {} #go by row, col
    inFile = open(data['ppmFile'],'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        nuc = line[0]
        line = line[1:]
        data['ppmWidth'] = len(line)
        for i in range(len(line)):
            data['ppmMat'][nuc,i] = float(line[i])
    inFile.close()
###############################################################################
def score_seq_ppm(mySeq,myPPM):
    score = 1.0
    for i in range(len(mySeq)):
        c= mySeq[i]
        c = c.upper()
        score = score * myPPM[c,i]
    return score
###############################################################################
# select random position, proportional to ppm weights
def select_random_position_with_weights(data):
    # keep going until get one that passess
    numTry = 0
    while True:
        randSel = select_random_position(data)
        score_selected(data,randSel)
        numTry += 1
        if randSel['isChosen'] is True:
            return randSel
#    print 'rand sel was',randSel
#    print 'Took %i tries to get one!' % numTry
###############################################################################
# extract the sequence, and figure what the score is
def score_selected(data,randSel):
    randSel['isChosen'] = False
    
    # check to see if is in exclusion region
    if data['useExclusionRegions'] is True:
        m = data['exclusionLists'][randSel['chrom']][randSel['pos']]
        if m != 0:
            return

    if randSel['strand'] == '+':  #insert is in genome orientaiton, cut is after nuc, motif on opposite strand
        startBp = randSel['pos']
        endBp = startBp + data['ppmWidth'] # does not include last bp
        if endBp >= data['genomeSeqDict'][randSel['chrom']]['seqLen']:
            return  #off end of chrom, is false
        extractedSeq = data['genomeSeqDict'][randSel['chrom']]['seq'][startBp:endBp]
        extractedSeq = extractedSeq.upper()
        extractedSeq = revcomp(extractedSeq)
        randSel['extractedSeq'] = extractedSeq
        if 'N' in extractedSeq:
            return # has an N in it -- cannot score
        seqScore = score_seq_ppm(extractedSeq,data['ppmMat'])
        r = random.random()
        if r < seqScore:
            randSel['isChosen'] = True
            return
        else:
            return    
    else: # is minus
        endBp = randSel['pos']
        startBp = endBp - data['ppmWidth'] # does not include last bp    
        if startBp < 0:
            return  #off end of chrom, is false
        # so that we include the endBp    
        extractedSeq = data['genomeSeqDict'][randSel['chrom']]['seq'][startBp+1:endBp+1]
        extractedSeq = extractedSeq.upper()
        randSel['extractedSeq'] = extractedSeq
        if 'N' in extractedSeq:
            return # has an N in it -- cannot score
        seqScore = score_seq_ppm(extractedSeq,data['ppmMat'])
        r = random.random()
        if r < seqScore:
            randSel['isChosen'] = True
            return
        else:
            return    
###############################################################################


