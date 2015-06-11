import sys
from optparse import  OptionParser

###############################################################################
# returns as string
def matrix_to_string(matrix,alphabet):
    matW = len(matrix)
    matH = len(matrix[0])
    myString = []
    if matH != len(alphabet):
        print 'matrix height and alphabet length do not match!'
        sys.exit()
    for row_i in range(matH):
        s = []
        s.append(alphabet[row_i])
        for col_i in range(matW):
            s.append(matrix[col_i][row_i])
        s = [str(j) for j in s]
        s = ' '.join(s)
#        print s
        myString.append(s)
    myString = '\n'.join(myString)
    return myString
###############################################################################
###############################################################################
# returns as string
def convert_to_freqs(matrix):
    matW = len(matrix)
    matH = len(matrix[0])
    freqMat = []
    for col_i in range(matW):
        s = sum(matrix[col_i])
        s = float(s)
        n = []
        for row_i in range(matH):
            f = matrix[col_i][row_i]
            f = float(f)/s
            n.append(f)
        freqMat.append(n)
    return freqMat
###############################################################################
###############################################################################
USAGE = """
python create-ppm.py --in <input list of sequences>  --out <output for ppm file>

--pseduo <pseduo count for each column, default = 1>

    --in <input file of observed sequence motifs, 1 per line.  All must be same length>
    --out <output name for position probability matrix file>

"""
parser = OptionParser(USAGE)
parser.add_option('--in',dest='inputFile', help = 'input file of sequences')
parser.add_option('--out',dest='outFile', help = 'output file name')
parser.add_option('--pseduo',dest='pseduo',type='int',default='1', help = 'pseudo count total for each column [default 1]')

(options, args) = parser.parse_args()

if options.inputFile is None:
    parser.error('input observed motifs not given')
if options.outFile is None:
    parser.error('output ppm file name not given')

###############################################################################



seqCounts = {}
allSeqs = []
inFile = open(options.inputFile,'r')
for line in inFile:
    seq = line.rstrip()
    if seq not in seqCounts:
        seqCounts[seq] = 0
    seqCounts[seq] += 1
    allSeqs.append(seq)
inFile.close()

seqs = seqCounts.keys()
seqs.sort()

print 'Read in %i total sequences, of which %i are unique\n' % (len(allSeqs),len(seqCounts))

print 'motif\toccurences'
for s in seqs:
    print '%s\t%i' % (s,seqCounts[s])
    

seqLen = len(seqs[0])
# check that are all the same sequence
for seq in seqs:
    if len(seq) != seqLen:
        print 'ERROR, not all the same length %i %s %i' % (seqLen,seq,len(seq))
        sys.exit()


#default to DNA -- can change to other later if needed
# make dictioanry to go from 'A' --> 0, 'C' --> 1, etc
alphabet = ['A','C','G','T']
nucToPos = {}
for i,c in enumerate(alphabet):
    nucToPos[c] = i

psudeoPerChar = options.pseduo / float(len(alphabet))
print '\nAdding pseduo count of %f for each character in each position column' % psudeoPerChar

counts = []
for i in range(seqLen):
    empty = []
    for j in alphabet:
        empty.append(psudeoPerChar)  #start with pseduo counts of 0.5
    counts.append(empty)

print '\nInitial matrix'    
myS = matrix_to_string(counts,alphabet)    
print myS



for seq in allSeqs:
    for i in range(len(seq)):
        c = seq[i]
        counts[i][nucToPos[c]] += 1

print '\nRaw Counts'    
myS = matrix_to_string(counts,alphabet)    
print myS


freqMat = convert_to_freqs(counts)

print '\nRaw Freqs'    
myS = matrix_to_string(freqMat,alphabet)    
print myS

outFile = open(options.outFile,'w')
outFile.write(myS)
outFile.close()
