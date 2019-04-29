import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import collections

for fileidx in (1, 2):
	d = pd.read_csv('xml%d.csv' % fileidx, index_col = 0)
	d.hit_def = d.hit_def.apply(lambda x: x.lower())
	counts = collections.Counter(d.hit_def)
	counts = dict(sorted(counts.items(), key = lambda x: x[1], reverse = True)[:10])

	plt.barh(width = list(counts.values()), y = [x.replace(' ', '\n') for x in counts.keys()])
	plt.xlabel('Counts of Reads')
	plt.title('Barcode%d' % fileidx)
	plt.show()
