import pandas as pd
import numpy as np

filename = '/mnt/d/Codes/DataLists/barcode01.out'
data = pd.read_csv(filename, sep = '\t').iloc[:,[0, 1, 11]]
data.columns = ['query', 'ret', 'score']
data.ret = [ret.split('.')[0] for ret in data.ret]
data = data \
        .groupby(['query']) \
        .apply(lambda x: x.loc[x.score > np.quantile(x.score, 0.75),]) \
        .reset_index(drop = True)
print(data.sort_values(by = 'score'))
# data = data.loc[data.score == data.groupby(data.query).apply(max, broadcast = True)]
data = data.loc[data.score > 100,]
print(data.groupby(data.ret).apply(len).sort_values())
