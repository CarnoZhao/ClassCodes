'''
Todo:

    1. How many barcodes were analyzed

    2. How many reads per barcode

    3. How many total reads

    4. The distribution of read length per barcode

    5. The total number of bases sequenced per barcodes

    6. Taxonomy of all barcodes

    *7. Coverage to the genome

    *8. Assembly the genome and the assembly statistics
'''

from Bio import SeqIO, Seq
import gzip
from collections import defaultdict

class Genome():
    def __init__(self, filename):
        self.filename = filename

    def load(self):
        '''
        Load the file to a generator to save memory
        '''
        try:
            self.f = SeqIO.parse(self.filename, 'fastq')
        except:
            self.f = SeqIO.parse(gzip.open(self.filename, 'rt', encoding = 'utf-8'), 'fastq')

    def barcode_statistics(self):
        '''
        Count the number of barcodes and the number of reads per barcode.
        And do further analysis.
        (Todo 1-5: DONE !)
        '''
        self.load()
        barcodes = defaultdict(int)
        for rec in self.f:
            des = rec.description
            barcode = des.split(' ')[-1].split('=')[-1]
            barcodes[barcode] += 1
        print('In file %s, there are %d barcodes.' % (self.filename.split('/')[-1], len(barcodes)))
        for barcode, num in barcodes.items():
            print('For %s, there are %d reads' % (barcode, num))

    def 
        
