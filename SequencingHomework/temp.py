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
from Bio.Blast import NCBIWWW
import random
from bs4 import BeautifulSoup

class Genome():
    def __init__(self, filename):
        self.filename = filename
        self.basename = filename.split('/')[-1].split('.')[0]

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
        print('In file %s, there are %d barcodes.' % (self.basename, len(barcodes))) # Todo 1
        print('The total number of reads is %d' % sum(len(reads) for reads in barcodes.values())) # Todo 3
        print('The total number of bases is %d' % sum(sum(reads) for reads in barcodes.values())) # Not mentioned in Todo
        pdf = PdfPages('./read_length_%s.pdf' % self.basename)
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

    def taxonomy_blast(self, lenth = 1000, squeeze = 10, method = 'blastn', db = 'nt'):
        '''
        qblast is too slow. A better solution is to blast manually in web browser
        Here is the code that can extract some 1kb fragments from reads. 
        The longer the read is, the more fragments will be extracted. The ratio is about 1 1k fragment per 10k read.
        '''
        self.load()
        query_seq = ''
        for rec in self.f:
            seqlen = len(rec.seq)
            if seqlen < lenth:
                continue
            for _ in range(1 + seqlen // (lenth * squeeze)): # roughly speaking, get 1k bases per 10k bases
                idx = random.randint(0, seqlen - lenth)
                seq = rec.seq[idx:idx + lenth]
                query_seq += '>' + rec.description + '\n' + str(seq) + '\n' # convert to fasta format
        res = NCBIWWW.qblast(method, db, query_seq)
        with open('%s.xml' % self.basename, 'w') as fw:
            fw.write(res.read())

    def taxonomy_read_xml(self):
        '''
        analyze the xml result from blast, using bs4 package (BeautifulSoup)
        Not sure that Bio.Blast.NCBIXML works for multi-reads results, so I use bs4 to deal with tags.
        available information: 
            hit_def, 
            hit_len, 
            hit_num, 
            hit_hsps:{
                hsp_bit-score, 
                hsp_evalue, 
                hsp_identity, 
                hsp_align-len}
        '''
        with open('%s.xml' % self.basename) as f:
            soup = BeautifulSoup(f.read(), 'html.parser')
        df = pd.DataFrame({})
        tag = 'hit_def' # species tag
        tlist = [t.text for t in soup.find_all(tag)]
        tlist = [' '.join(t.split(' ')[:2]) if 'PREDICTED:' not in t else ' '.join(t.split(' ')[1:3]) for t in tlist]
        df['species'] = tlist
        df.to_csv('%s.csv' % self.basename)


if __name__ =='__main__':
    genome = Genome('/mnt/d/Grocery/DataSet/DNAlab/barcode1.fastq')
    # genome.barcode_statistics()
    # genome.taxonomy_blast()
        
