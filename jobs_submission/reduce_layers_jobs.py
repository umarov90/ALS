import os
from pathlib import Path
import sys

with open("samples.txt", 'r') as file:
    sample_ids = file.readlines()
sample_ids = [line.strip() for line in sample_ids]
print(sample_ids)

wd = os.path.abspath(sys.argv[1])
sd = os.path.dirname(os.path.realpath(__file__))
Path(sd + "/temp").mkdir(parents=True, exist_ok=True)
i = 0
for sid in sample_ids:
    with open(f"{sd}/temp/job{i}", 'w+') as f:
        f.write("#!/bin/bash\n")
        f.write(f"#SBATCH --job-name=reduce{i}\n")
        f.write(f"#SBATCH --output={sd}/temp/job{i}.out\n")
        f.write(f"#SBATCH --error={sd}/temp/job{i}.err\n")
        f.write(f"#SBATCH --mem=64G\n")
        f.write(f"#SBATCH --partition=batch\n")
        f.write(f"python 10_batch_merge_velocyto_layers.py {wd} {sid}")
    os.system(f"sbatch {sd}/temp/job{i}")
    i += 1
