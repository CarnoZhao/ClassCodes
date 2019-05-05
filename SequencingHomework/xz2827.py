#!/usr/bin/env python
# coding: utf-8

# # Field Trip Report
# 
# **PRINCIPLES & APPS MOD DNA SEQUENCING TEC**
# 
# Xun Zhao, xz2827
# 
# ## Background
# 
# In the field trip, some DNA samples from different species were sequenced, using NanoPore sequencing technology. The seqeunce data were given in FASTQ format, including 9 barcodes, which also means 9 species. The purpose of this report is to identify the taxonomy of these samples, and do some statistical analysis about the sequencing result, based on Python code and Bash script.
# 
# ## Import Packages
# 
# The code below is to import necessary packages and define a function that can read several FASTQ file that in one `barcode0*` directory, namely, all the FASTQ that have the same barcode. The function uses `SeqIO` package to parse the FASTQ format and returns a generator of each read. When iterating these reads, the description, sequence, quality and other information is available. 

# In[220]:


from Bio import SeqIO, GenBank, Entrez
from bs4 import BeautifulSoup
from collections import defaultdict
import os, json, time, random, matplotlib
import matplotlib.pyplot as plt
import pandas as pd
get_ipython().run_line_magic('matplotlib', 'inline')
global path
path = '/home/jovyan/ro-data/nanopore_raw_data/'
Entrez.email = 'xz2827@columbia.edu'

def load(barcode):
    filelist = os.listdir(path + barcode)
    for filename in filelist:
        if not filename.startswith('fastq_'):
            continue
        for read in SeqIO.parse(path + barcode + '/' + filename, 'fastq'):
            yield read


# ## Barcode Statistics
# 
# Here I iterate all the reads of one same barcode and store the length of each reads to a list. So the sum of the list is the total number of bases sequenced in this barcode directory, and the length of this list is the number of reads of this barcode.
# 
# Using the list that contains the length of each read, I drew the histogram of the distribution of reads length and calculate the total number of reads. 

# In[212]:


def barcode_statistics(barcode, bins = 50, length = None, text = None, if_return = True):
    '''
    1. read from files
    2. get a list of length of reads
    3. plot the distribution of reads length, sum of reads length (number of bases), and number of reads
    4. save all barcodes' list to a `total` list, plot the distribution of total sequencing result
    
    args:
        bins = 50:         number of bins in histogram
        length = None:     if `length` is given, use it instead of getting it from files
        text = None:       some notes that will show in the box in figures
        if_return = True:  define whether the function will return a list or None
    '''
    length = [len(read.seq) for read in load(barcode)] if length == None else length
    fig, ax = plt.subplots(figsize = (10, 7))
    ax.hist(length, bins = bins, color = 'skyblue')
    plt.yscale('log')
    plt.xlabel('Reads Length')
    plt.ylabel('Count of Reads')
    plt.title(barcode)
    label = '%s:\nReads: %dK (%d)\nBases: %dM (%d)' % (barcode,
            len(length) / 1000, len(length), 
            sum(length) / 1000000, sum(length))
    label += '\n' + text if text != None else ''
    plt.text(0.5, 
             0.75, 
             label,
            fontsize = 16,
            transform = ax.transAxes,
            bbox = dict(facecolor = 'grey', alpha = 0.5))
    plt.show()
    return length if if_return else None

total = []
cnt = 0
for barcode in sorted(x for x in os.listdir(path) if 'barcode0' in x):
    total.extend(barcode_statistics(barcode))
    cnt += 1
barcode_statistics('Total', length = total, text = 'There are %d barcodes' % cnt, if_return = False)


# ## Taxonomy Identification
# 
# The major tool that I use to identify the taxonomy of sequeces is BLAST. To prepare, the FASTQ file needs to be converted to FASTA file (`extract_fasta`). Then, I pass the FASTA file to local BLAST, using refseq database (`bash script`). The BLAST returns results that contain the quality information of alignment and the accession ID of each result sequence. So I define a filter function (`extract_high_score`) that remove low quality (especially short alignment length), because no matter which barcode I analysis, the majority results are E. coli. After filtering, I can get more information about the species and taxonomy from accession ID through Entrez API (`search_organisms`, `get_taxonomy`). Finally, I use weighted bit-score to calculated the probability of each possible result and visualize it (`organisms_histogram`, `plot_taxonomy`).

# ### Extract FASTA File
# 
# Both the local and remote BLAST needs FASTA file, and we need to extract small fraction of the total reads to BLAST. 
# 
# First, I exprect get 1% of the total reads, which equals to the `squeeze` parameter. Besides, I restrict all the smaller fragments to be exact 1kb long, which equals to the `length` parameter. So I ignore reads that are shorter, randomly choose reads using 0.01 as the probability. In addition, for a long read, I randomly choose 1kb part of it.
# 
# Then, I covert it to the FASTA format and store it in a new file.

# In[ ]:


def extract_fasta(barcode, length = 1000, squeeze = 0.001, outdir = 'fasta', minimum = 200):
    '''
    1. create an output directory
    2. read from files
    3. skip reads shorter than 1kb
    4. random choose a random 1kb fragment by `random.random() < squeeze`
    5. if number of reads < 200, repeat
    6. write them into FASTA
    
    args:
        length = 1000:     length of FASTA reads
        squeeze = 0.001:   get only 0.1% of reads from FASTQ
        outdir = 'fasta':  output dir
        minimum = 200:     at least get 200 reads
    '''
    faname = barcode + '.fasta'
    os.system('mkdir -p %s' % outdir)
    fasta = open(outdir + '/' + faname, 'w')
    cnt = 0
    while cnt < minimum:
        for read in load(barcode):
            seqlen = len(read.seq)
            if seqlen < length:
                continue
            elif random.random() < squeeze:
                idx = random.randint(0, seqlen - length)
                seq = str(read.seq[idx: idx + length])
                fasta.write('>' + read.description + '\n' + seq + '\n')
                cnt += 1
    fasta.close() 

for barcode in sorted(x for x in os.listdir(path) if 'barcode0' in x):
    extract_fasta(barcode)


# ### Local BLAST
# 
# The BLAST is using `refseq_genomic` database and tabular output format. It takes a lot of time to run, about 1 hour for a 1Mb FASTA file.

# In[ ]:


get_ipython().run_cell_magic('bash', '', '\n# args:\n#     -db:           refseq_database\n#     -query:        input FASTA file\n#     -out:          output file\n#     -outfmt:       output format, default 6 (tabular without header)\n#     -num_threads:  number of threads\noutdir="./blast_out/"\nmkdir -p $outdir\nfor name in ./fasta/barcode*\ndo\n    out=$outdir$(basename $name ".fa").out\n    echo $out\n    nohup blastn -db /home/jovyan/ro-data/blast_databases/refseq_genomic -query $name -out $out -outfmt 6 -num_threads 20\ndone')


# ### Filter BLAST Result
# 
# The analysis focuses on `BitScore`, which in some way is positively correlated with both alignment quality and alignment length. Both identity and E-value does not provide the alignment length information, which will consider some shorter but more accurate alignments to be better, instead of longer alignments.
# 
# This filter function is very customizable. The `query_choose` parameter means that for each query sequence, I only use the best ten results. And in total, I only keep 50 results with highest score, which is the `total_choose` parameter. Moreover, as my query is 1kb, the alignment is not informative if `BitScore` is lower than 100 (`threshold`), becuase there are many E. coli data that can map to query sequence, whose alignment length are shorter than 100.
# 
# For example:
# 
# `
# Query Result Score
#  1     1.1    1000
#  1     1.2    900      ==> Best 10 of Query1 
#  1     1.3    800                            
#  ...                                                 ==(sort)==>   Best 50
#  2     2.1    900                            
#  2     2.2    800      ==> Best 10 of Query2 
#  ...                      
#  3     3.1    90 (delete)
# `
# 
# Finally, I store the query name, result accession ID and bit score in a `csv` file.

# In[8]:


def extract_high_score(barcode, query_choose = 10, total_choose = 50, threshold = 100):
    '''
    1. read from BLAST output file (tabular format ==> pandas.DataFrame)
    2. choose columns (0: query name, 1: result ID, 11: bitscore)
    3. clean result ID (e.g. NC_12345.1 ==> NC_12345)
    4. filter data
        4.1. remove bitscore lower than threshold (default 100)
        4.2. make groups by query name
        4.3. for each group, choose the best 10 results (query_choose default 10)
        4.4. choose the best 50 from all chosen results (total_choose default 50)
    5. save to csv file
    
    args:
        query_choose = 10:  number of best results chosen for each query group
        total_choose = 50:  total number of results
        threshold = 100:    remove the results with bitscore lower than 100
    '''
    out = './blast_out/%s.out' % barcode
    csv = './blast_out/%s.csv' % barcode
    data = pd.read_csv(out, sep = '\t', header = None).iloc[:,[0, 1, 11]]
    data.columns = ['query', 'ret', 'score']
    data.ret = [ret.split('.')[0] for ret in data.ret]
    data = data         .loc[data.score > threshold,:]         .groupby('query')         .apply(lambda x: x             .sort_values(by = 'score')             .iloc[-query_choose:,])         .sort_values(by = 'score')         .iloc[-total_choose:,:]         .reset_index(drop = True)
    data.to_csv(csv, header = False)
    
for barcode in sorted(x for x in os.listdir(path) if 'barcode0' in x):
    extract_high_score(barcode)


# ### Search Species
# 
# I use the `Bio.Entrez.efetch` to search the accession ID obtained above, and use `Bio.GenBank.parser` to get species name from the `efetch` result.

# In[22]:


def search_organisms(barcode):
    '''
    1. read from csv 
    2. extract the column with accession ID of results
    3. get the species information from nucleotide database through Entrez API
    4. save to file 
    '''
    df = pd.read_csv('./blast_out/%s.csv' % barcode, sep = ',', header = None)
    ids = df.iloc[:,2].tolist()
    hdl = Entrez.efetch(db = 'nucleotide', id = ids, rettype = 'gb')
    recs = GenBank.parse(hdl)
    with open('./blast_out/%s.org' % barcode, 'w') as fw:
        fw.write('\n'.join(rec.organism for rec in recs))

for barcode in sorted(x for x in os.listdir(path) if 'barcode0' in x):
    search_organisms(barcode)


# ### Plot Species Histogram
# 
# The alignment length varies a lot, from hundreds to tens. So I dicide not to simply use count of reads to plot the graph. Instead, the bit scores are considered as the weight of each reads. That is, the longer alignments, which uasually have higher bit score, have higher weight. Conversely, those short alignments, which are more, have lower weight. 
# 
# For example:
# 
# `Species   BitScore                                 
#  E.coli    10           Species    WeightedScore    
#  E.coli    11     ==>   E.coli     33 = 10 + 11 + 12
#  E.coli    12           H.sapiens  250 = 100 + 150  
#  H.sapiens 100                                      
#  H.sapiens 150                                     .
# `
# 
# In the case above, even though there are 3 E.coli alignments, the score is too low. So we should consider H.sapiens with higher score.

# In[214]:


def organisms_histogram(barcode, num_species = 10):
    '''
    1. read from output files from last two cells
    2. extract bit-score column (`score`) and species column (`org`)
    3. make groups by different species names
    4. calculate the sum of scores as the final weighted score
    5. keep 10 species with highest weighed score
    6. plot the bar plot
    7. save all species names and scores to a new csv (for taxonomy)
    
    args:
        num_species = 10: only plot the 10 highest score species
    '''
    df = pd.read_csv('./blast_out/%s.csv' % barcode, sep = ',', names = ['query', 'ret', 'score'])
    orgs = pd.read_csv('./blast_out/%s.org' % barcode, sep = ',', names = ['orgs'])
    df['orgs'] = [' '.join(org.split(' ')[:2]) for org in orgs.orgs]
    species = df         .groupby('orgs')         .apply(lambda x: sum(x.score))         .sort_values(ascending = True)
    highest_10 = species.iloc[:num_species][::-1]
    plt.figure(figsize = (16, 9))
    graph = plt.barh(
        highest_10.index,
        highest_10.values,
        color = 'skyblue')
    for i, rect in enumerate(graph):
        label = '%s: %d' % (highest_10.index[i], highest_10.values[i])
        plt.text(0.1 * rect.get_width(), 
                 rect.get_y() + rect.get_height() * 0.5, 
                 label,
                 fontsize = 15,
                 verticalalignment = 'center')
    plt.yticks([])
    plt.ylabel('Species')
    plt.xlabel('Sum of scores')
    plt.title(barcode)
    plt.show()
    species.to_csv('./blast_out/%s.species' % barcode, header = False)

for barcode in sorted(x for x in os.listdir(path) if 'barcode0' in x):
    organisms_histogram(barcode)


# ### Search Taxonomy
# 
# In last cell, I saved the species into a new `csv` file. In NCBI's taxonomy database, we can get the taxonomy information. Here, I only use 7 normal ranks shown below, without sub- or super- ranks.
# 
# `kingdom, phylum, class, order, family, genus, species`
# 
# Like last cell, I still use weighted score to analysis the probability. 
# 
# I use `json` format to save the data instead of normal `pandas.DataFrame`.
# 
# The example structure is:
# 
# `{'kingdom':
#     {
#         'Metazoa': 100,
#         '...': ...
#     }
#   'phylum':
#     {
#         ...
#     }
#   ...
#  }`

# In[139]:


def get_taxonomy(barcode, which_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']):
    '''
    1. read from files
    2. for each line, get the species name and its weightes score
    3. using species name to get the species ID through `Entrez.esearch`
    4. using species ID to get the taxonomy information through `Entrez.efetch`
        4.1. already define the rank list with 7 common ranks
        4.2. from the returned HTML format, get the `taxon` tags
        4.3. skip the rank name that do not exist in rank list (`ranks`)
        4.4. get the scientfic nmae from each taxon
        4.5 add score to this name
    5. save the score dictionary to JSON format output file
    
    args:
        which_ranks: the rank names that will be used in this code
    '''
    ranks = dict((rank, defaultdict(int)) for rank in which_ranks)
    f = open('./blast_out/%s.species' % barcode)
    for line in f:
        org, score = line.strip().split(',')
        handle = Entrez.esearch(term = org, db = 'Taxonomy')
        id_list = Entrez.read(handle)['IdList']
        handle = Entrez.efetch(db = 'Taxonomy', id = id_list, rettype = 'gb')
        soup = BeautifulSoup(handle.read(), 'html.parser')
        for taxon in soup.find_all('taxon'):
            rank = taxon.find('rank').text
            name = taxon.find('scientificname').text
            name = name if rank != 'species' else name.replace(' ', '\n')
            if rank in ranks:
                ranks[rank][name] += eval(score)
    with open('./blast_out/%s.species.json' % barcode, 'w') as fw:
        fw.write(json.dumps(ranks))

for barcode in sorted(x for x in os.listdir(path) if 'barcode0' in x):
    get_taxonomy(barcode)
    time.sleep(120)


# ### Plot Taxonomy
# 
# In last cell, the returned data contains 7 ranks. In each rank, there is a `dict` whose keys are the names of taxa and the values are the score. So I normalize the scores to 1 to represent probability. The colors are randomly choosen by `matplotlib`, which do not indicate the taxonmy relationship.

# In[217]:


def plot_taxonomy(barcode):
    '''
    1. read from JSON file
    2. normalize the scores to probabilities 
    3. plot the probabilities as the figure shows
    '''
    with open('./blast_out/%s.species.json' % barcode) as f:
        ranks = json.load(fp = f)
    lenrank = len(ranks)
    maxdepth = max(len(names) for names in ranks.values())
    for rank in ranks:
        names = ranks[rank]
        names = list(names.items()) + [("", 0)] * (maxdepth - len(names))
        ranks[rank] = sorted(names, key = lambda x: x[1])
    datalist = [[name[1] for name in names] for names in ranks.values()]
    datalist = [[y / sum(x) for y in x] for x in datalist]
    namelist = [[name[0] for name in names] for names in ranks.values()]
    rankslist = [rank.capitalize() for rank in ranks]
    bottom = [0 for i in range(lenrank)]
    plt.figure(figsize = (18, 10))
    for i in range(maxdepth):
        rects = plt.barh(range(lenrank), 
                        [x[i] for x in datalist], 
                        left = bottom, 
                        tick_label = rankslist)
        for j, rect in enumerate(rects):
            width = rect.get_width()
            if namelist[j][i] == '' or width < 0.08:
                label = ''
            else:
                label = '%d%%\n%s' % (width * 100, namelist[j][i])
            plt.text(0.5 * width + bottom[j],
                     rect.get_y() + 0.5 * rect.get_height(),
                     label,
                    fontsize = 11 if width < 0.16 else 16,
                    horizontalalignment = 'center',
                    verticalalignment = 'center')
        bottom = [x + y for x, y in zip(bottom, [x[i] for x in datalist])]
    plt.xlabel('Sum of probability')
    plt.title(barcode)
    plt.show()
    
for barcode in sorted(x for x in os.listdir(path) if 'barcode0' in x):
    plot_taxonomy(barcode)

