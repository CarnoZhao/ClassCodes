import pandas as pd
import numpy as np

filename = '/mnt/d/Codes/DataLists/barcode01.out'
data = pd.read_csv(filename, sep = '\t').iloc[:,[0, 1, 11]]
data.columns = ['query', 'ret', 'score']
data.ret = [ret.split('.')[0] for ret in data.ret]
data = data \
        .groupby(['query']) \
        .apply(lambda x: x.loc[x.score > np.quantile(x.score, 0.9),]) \
        .reset_index(drop = True)
# print(data.sort_values(by = 'score'))
# data = data.loc[data.score == data.groupby(data.query).apply(max, broadcast = True)]
data = data.loc[data.score > 100,]
orgs = pd.read_csv('./analysis/orgs.txt', sep = '\t', header = None)[0].tolist()
accs = pd.read_csv('./analysis/2.txt', sep = '\t', header = None)[0].tolist()


ret = data.groupby('ret').apply(len).sort_values()

# ret['org'] = ret.apply(lambda x: orgs[accs.index(x.ret)])
# ret.index = [orgs[accs.index(x)] if x in accs else None for x in ret.index]
print(ret)
