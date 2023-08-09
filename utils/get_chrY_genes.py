import anndata as ad
from datetime import datetime
import pandas as pd
import utils.common as cm
import re
from params import Params
import numpy as np


p = Params()
chrY_genes = []
with open(p.folder + "genes.gtf", 'r') as file:
    for line in file:
        fields = line.strip().split('\t')
        if fields[0] == 'chrY' and fields[2] == 'gene':
            gene_name_pattern = r'gene_name "([^"]+)"'
            match = re.search(gene_name_pattern, fields[8])
            if match:
                gene_name = match.group(1)
                chrY_genes.append(match.group(1))
with open(p.folder + "chrY_genes.csv", 'w') as file:
    file.write("\n".join(chrY_genes))
