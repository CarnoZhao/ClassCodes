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
import os
import gzip
import random
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
%matplotlib inline

class Genome():
    def __init__(self, dirname):
        self.dirname = dirname
        self.outpath = '../dirname'

    def load(self):
        filelist = os.listdir(self.outpath)
        for filename in filelist:
            for read in SeqIO.parse(self.outpath + '/' + filename, 'fastq'):
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
        os.system('mkdir -p fasta')
        fasta = open('fasta/%s.fa' % self.dirname)
        for read in self.load():
            seqlen = len(read.seq)
            if seqlen < length:
                continue
            else:
                for _ in range(1 + seqlen // (length * squeeze)):
                    idx = random.randint(0, seqlen - length)
                    seq = str(rec.seq[idx: idx + length])
                    fasta.write('>' + read.description + '\n' + seq + '\n')
        fasta.close() 
    
    def local_blast(self, method = 'blastn', db = 'nt'):
           pass
           
        
