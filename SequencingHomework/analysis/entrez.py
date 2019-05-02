from Bio import Entrez, SeqIO, Seq, GenBank
import pandas as pd
import collections

Entrez.email = 'xz2827@columbia.edu'

def func():
    data = pd.read_csv('../../../DataLists/barcode01.out', sep = '\t')

    data.columns = ['qaccver', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    ret = data['saccver'].tolist()
    print(ret)
    return ret

def get_org(accs):
    # phylo = 'kingdom-phylum-class-order-family-genus-species'.split('-')
    #df = dict((x, []) for x in phylo)
    orgs = []
    hdl = Entrez.efetch(db = 'nucleotide', id = accs, rettype = 'gb')
    recs = GenBank.parse(hdl)
    for i, rec in enumerate(recs):
        # tax = rec.taxonomy
        org = rec.organism
        orgs.append(org)
        # print(tax)
        # for j, t in enumerate(tax):
        #   df[phylo[j]].append(t)
        # df[phylo[j + 1]].append(org)
    # print(pd.DataFrame(df))
    print(collections.Counter(orgs))

get_org(func())

