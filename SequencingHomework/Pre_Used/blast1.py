from numpy import random
from Bio.Blast import NCBIWWW
from Bio import SeqIO

cnt = 100

for i in (1, 2):
    fileReads = SeqIO.parse('barcode%d.fastq' % i, 'fastq')
    reads = random.randint(0, 4000, size = cnt) if i == 1 else random.randint(0, 20000, size = cnt)
    reads = sorted(reads, reverse = True)
    
    start = 0
    while reads:
        onereads = next(fileReads)
        if start == reads[-1]:
            reads.pop()
            lenth = len(onereads.seq)
            if lenth < 1000:
                seq = onereads.seq
            else:
                idx = random.randint(0, lenth - 1000)
                seq = onereads.seq[idx: idx + 1000]
            print(seq[:10], '...', seq[-10:])
            blastResult = NCBIWWW.qblast('blastn', 'nt', seq)
            with open('./result/blastResult%d_%d.txt' % (i, start), 'w') as f:
                f.write(blastResult.read())
            print('Where:', i, start)
        else:
            start += 1
