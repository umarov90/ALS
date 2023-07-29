import os
import numpy as np
import scipy.sparse as sp
from scipy.io import mmread
from params import Params
import pandas as pd
import utils.common as cm
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('macosx')

p = Params()

df = pd.read_csv(f'{p.folder}meta.tsv', sep='\t', index_col=False)
adatas = []
for index, row in df.iterrows():
    mtx_path = p.folder + f'cellranger_output/{row["sample_id"]}/outs/raw_feature_bc_matrix/matrix.mtx.gz'
    if not os.path.isfile(mtx_path):
        continue
    print(mtx_path)
    matrix = mmread(mtx_path)
    umi_counts = np.asarray(matrix.sum(axis=0)).flatten()
    sorted_indices = np.argsort(umi_counts)[::-1]
    sorted_umi_counts = umi_counts[sorted_indices]

    plt.figure(figsize=(8, 6))
    plt.plot(range(len(sorted_umi_counts)), sorted_umi_counts, marker='o')
    plt.xlabel('Droplet ID (Ranked by Count)')
    plt.ylabel('UMI Count per Droplet')
    plt.title('UMI Count per Droplet - Ranked by Count')
    plt.show()
    plt.savefig(cm.ensure(f'{p.folder}UMI_filtered/{row["sample_id"]}.png'))