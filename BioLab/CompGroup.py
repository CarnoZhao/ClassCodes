import random

class CompGroup():
    def __init__(self, max_of_genes = 10):
        if max_of_genes < 4:
            raise Exception('Max number of genes too small')
        self.num_genes = random.randint(4, max_of_genes)
        self.num_colors = random.randint(4, num_genes)
        self.gene_name = [chr(ord('Z') - i) for i in range(num_genes - 1, -1, -1)]
        self.color_name = [chr(ord('A') + i) for i in range(num_colors)]

