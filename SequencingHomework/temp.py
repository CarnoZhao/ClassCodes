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
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import numpy as np

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
        barcodes = defaultdict(list)
        for rec in self.f:
            des = rec.description
            barcode = des.split(' ')[-1].split('=')[-1]
            barcodes[barcode].append(len(rec.seq))
        print('In file %s, there are %d barcodes.' % (self.filename.split('/')[-1], len(barcodes))) # Todo 1
        print('The total number of reads is %d' % sum(len(reads) for reads in barcodes.values())) # Todo 3
        print('The total number of bases is %d' % sum(sum(reads) for reads in barcodes.values())) # Not mentioned in Todo
        pdf = PdfPages('./read_length_distribution.pdf')
        for barcode, reads in barcodes.items():
            print('For %s, there are %d reads and %d bases' % (barcode, len(reads), sum(reads))) # Todo 2, 5
            fig, ax = plt.subplots(figsize = (16, 9))
            plt.hist(reads, bins = 50)
            plt.yscale('log')
            plt.xlabel('Reads Length')
            plt.ylabel('Number of Reads')
            plt.title(barcode)
            pdf.savefig(fig) # Todo 4
        pdf.close()

    def taxonomy_analysis(self):
        pass     

if __name__ =='__main__':
    genome = Genome('/mnt/d/Grocery/DataSet/DNAlab/barcode1.fastq')
    genome.barcode_statistics()
        
