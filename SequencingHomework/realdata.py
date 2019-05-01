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

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from bs4 import BeautifulSoup
from collections import Counter, defaultdict
import gzip
import random
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
%matplotlib.inline

class Genome():
    def __init__(self, dirname):
        self.dirname = dirname

    def load(self):
        filelist = os.system('ls %s' % self.dirname)
        for filename in filelist:
            for read in SeqIO.parse(filename, 'fasta'):
                yield read

    def barcode_statistics(self, bins = 50):
        length = [len(read.seq) for read in self.load()]
        print('Ther are %d reads and %d bases in %s.' % (len(length), sum(length), self.dirname))
        plt.subplots(figsize = (16, 9))
        plt.hist(length, bins = bins)
        plt.yscale('log')
        plt.xlabel('Reads Length')
        plt.ylabel('Count of Reads')
        plt.title(self.dirname)
        plt.show()
        
    def __extract_fasta(self, length = 1000, squeeze = 10):
        os.system('mkdir -p xz2827_Personal')
        fasta = open('/%s.fa')
