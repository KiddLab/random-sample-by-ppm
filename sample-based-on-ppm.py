#!/usr/bin/env python
# make random samples from genome, weighted by some PPM
# reads genome into memory, so be sure to have enough RAM!
# assumes that are looking at L1 endonuclease cleavage pattern

import sampleutils
import sys
from optparse import  OptionParser

###############################################################################
USAGE = """
sample-based-on-ppm.py --ppm <position probability matrix>  --genome <genome file to extract from> 
                       --outpre <prefix for random sample set files>
                       --n_per_set <number of random samples per set>
                       --n_sets <number of random sets to produce>
                       --exclusionbed <bed file of regions to exclude from sampling>

For example, to make 10 sets of 50 random samples, use --n_per_set 50 --n_sets 10.

"""

parser = OptionParser(USAGE)
parser.add_option('--ppm',dest='ppm', help = 'position probability matrix')
parser.add_option('--genome', dest='genome', help = 'genome fasta file for selection')
parser.add_option('--outpre', dest='outFilePrefix', help = 'prefix for output file')

parser.add_option('--n_sets', dest='numSet',type='int', help = 'number of random set to make')
parser.add_option('--n_per_set', dest='numPerSet',type='int', help = 'number of random samples per set to make')
parser.add_option('--exclusionbed', dest='exclusionBed', help = 'bed file of regions to exclude from sampling',default='NONE')


(options, args) = parser.parse_args()

parser = OptionParser(USAGE)
if options.ppm is None:
    parser.error('input ppm file name not given')
if options.genome is None:
	parser.error('genome fasta file not given')
if options.outFilePrefix is None:
	parser.error('output files prefix not given')
if options.numSet is None:
	parser.error('total number of sets not given')
if options.numPerSet is None:
	parser.error('number of random samples per set not given')
###############################################################################

myData = {}
# read in genome sequences to list of fasta files
myData['genomeFasta'] = options.genome
myData['ppmFile'] = options.ppm

if options.exclusionBed == 'NONE':
    myData['useExclusionRegions'] = False
else:
    myData['useExclusionRegions'] = True
    myData['exclusionBedFile'] = options.exclusionBed
    print 'Will exclude random sites that occur in file',options.exclusionBed


sampleutils.initialize_ppm(myData)
sampleutils.initialize_genome_sequences(myData)

print 'Initialization complete, ready to sample'

print 'Will make %i sets, each of which consists of %i random samples' % (options.numSet,options.numPerSet)

for setNum in range(options.numSet):
    outFileName = options.outFilePrefix + '.%i' %setNum
    print 'Sampled locations being written to',outFileName
    outFile = open(outFileName,'w')
    for sample_i in xrange(options.numPerSet):
        if sample_i % 100 == 0:
            print 'Did sample %i of %i...' % (sample_i,options.numPerSet)
        randSel = sampleutils.select_random_position_with_weights(myData)
        # write output as bed file for first 3 columns
        # pos+1, because pos is 0 zero based already....
        outFile.write('%s\t%i\t%i\t%s\t%s\n' % (randSel['chrom'],randSel['pos'],randSel['pos']+1,randSel['strand'],randSel['extractedSeq']))        
    outFile.close()


