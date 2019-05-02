from Bio import Entrez, SeqIO, Seq, GenBank
import pandas as pd
import collections

Entrez.email = 'xz2827@columbia.edu'

def func():
    data = pd.read_csv('../../../DataLists/barcode01.out', sep = '\t')
    data.columns = ['qaccver', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    data = data.groupby(data.qaccver).apply(lambda x: x[['saccver', 'bitscore']][x.bitscore == max(x.bitscore)])
    # ret = data['saccver'].tolist()
    return data

def get_org(accs):
    # phylo = 'kingdom-phylum-class-order-family-genus-species'.split('-')
    # df = dict((x, []) for x in phylo)
    orgs = []
    fw = open('orgs.txt', 'w')
    hdl = Entrez.efetch(db = 'nucleotide', id = accs.saccver.tolist(), rettype = 'gb')
    recs = GenBank.parse(hdl)
    for i, rec in enumerate(recs):
        # tax = rec.taxonomy
        org = rec.organism
        fw.write(org + '\n')
        # print(tax)
        print(org)
        orgs.append(org)
        # for j, t in enumerate(tax):
        #   df[phylo[j]].append(t)
        # df[phylo[j + 1]].append(org)
    # print(pd.DataFrame(df))
    print(collections.Counter(orgs))
    fw.close()

res = func()
# print(res.head(10))
get_org(res)

