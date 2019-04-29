import random
for i in (1, 2):
    f = open('barcode%d.fa' % i)
    ranline = [random.randint(0, 4000 if i == 1 else 20000) for _ in range(200)]
    ranline.sort()
    fw = open('clean%d.fa' % i, 'w')
    for j, line in enumerate(f):
        if j // 2 in ranline:
            if j % 2 == 0:
                pre = line
            elif len(line) > 1000:
                fw.write(pre)
                idx = random.randint(0, len(line) - 1000)
                fw.write(line.strip()[idx:idx + 1000] + '\n')
    f.close()
    fw.close()
