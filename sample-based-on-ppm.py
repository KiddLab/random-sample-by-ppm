#!/usr/bin/env python
# Jeff Kidd
# 1 June 2015
# make random samples from genome, weighted by some PSSM
# reads genome into memory, so be sure to have enough

import sampleutils
import sys
from optparse import  OptionParser

###############################################################################
USAGE = """
sample-based-on-ppm.py --ppm <position probability matrix>  --genome <genome file to extract from> 
                       --num <number of random samples to take>
                       --out <outfile for positions of random samples>

"""

parser = OptionParser(USAGE)
parser.add_option('--ppm',dest='ppm', help = 'position probability matrix')
parser.add_option('--genome', dest='genome', help = 'genome fasta file for selection')
parser.add_option('--num', dest='num',type='int', help = 'number of random samples to return')
parser.add_option('--out', dest='outFileName', help = 'name for output file')

(options, args) = parser.parse_args()

parser = OptionParser(USAGE)
if options.ppm is None:
    parser.error('input ppm file name not given')
if options.genome is None:
	parser.error('genome fasta file not given')
if options.outFileName is None:
	parser.error('output file not given')
if options.num is None:
	parser.error('number of random sample to return not given')
###############################################################################

myData = {}
# read in genome sequences to list of fasta files
myData['genomeFasta'] = options.genome
myData['ppmFile'] = options.ppm

sampleutils.initialize_ppm(myData)
sampleutils.initialize_genome_sequences(myData)


print 'Initialization complete, ready to sample'
print 'Sampled locations being written to',options.outFileName
outFile = open(options.outFileName,'w')
for sample_i in xrange(options.num):
    if sample_i % 100 == 0:
        print 'Did sample %i of %i...' % (sample_i,options.num)
    randSel = sampleutils.select_random_position_with_weights(myData)
    outFile.write('%s\t%i\t%s\t%s\n' % (randSel['chrom'],randSel['pos'],randSel['strand'],randSel['extractedSeq']))        
outFile.close()


