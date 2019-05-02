from Bio import Entrez, SeqIO, Seq
import pandas as pd

Entrez.email = 'xz2827@columbia.edu'

def func():
    data = pd.read_csv('../../../DataLists/barcode01.out', sep = '\t')

    data.columns = ['qaccver', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
# print(data.head())

def get_org(accs):
    hdl = Entrez.efetch(db = 'nucleotide', id = accs, rettype = 'gb')
    recs = SeqIO.parse(hdl, 'gb')
    rec = next(recs)
    print(rec.__dir__())

get_org(['NC_012759.1'])

