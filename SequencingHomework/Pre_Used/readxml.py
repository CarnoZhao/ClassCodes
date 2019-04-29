from bs4 import BeautifulSoup
import pandas as pd
import numpy as np

for fileidx in (1, 2):
    f = open('./result%d.xml' % fileidx)
    soup = BeautifulSoup(f.read(), 'html.parser')
    df = pd.DataFrame({})
    for tag in ('hit_def', 'hit_len', 'hit_num', 'hit_hsps'):
        if tag != 'hit_hsps':
            tlist = [_.text for _ in soup.find_all(tag)]
            if tag == 'hit_def':
                tlist = [' '.join(_.split(' ')[:2]) if 'PREDICTED:' not in _ else ' '.join(_.split(' ')[1:3]) for _ in tlist]
            df[tag] = tlist
        else:
            tlist = soup.find_all(tag)
            for hsptag in ('hsp_bit-score', 'hsp_evalue', 'hsp_identity', 'hsp_align-len'):
                newtlist = [sum(eval(subtag.text) for subtag in hsp.find_all(hsptag)) / len(hsp.find_all(hsptag)) for hsp in tlist]
                df[hsptag] = newtlist
    df.to_csv('xml%d.csv' % fileidx)
    f.close()
