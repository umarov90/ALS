import os
from pathlib import Path
import sys

wd = os.path.abspath(sys.argv[1])
sd = os.path.dirname(os.path.realpath(__file__))
Path(sd + "/temp").mkdir(parents=True, exist_ok=True)
i = 0
for entry in os.scandir(wd):
    if entry.is_dir():
        directory_name = entry.name
        directory_path = entry.path
        h5_path = f"{wd}/../parts/{directory_name}.h5ad"
        if os.path.isfile(h5_path):
            print("Skipping")
            continue
        loom_path = f"{wd}/{directory_name}/velocyto/{directory_name}.loom"
        if os.path.isfile(loom_path):
            if os.path.getsize(loom_path) / (1024 * 1024) > 10:
                with open(f"{sd}/temp/job{i}", 'w+') as f:
                    f.write("#!/bin/bash\n")
                    f.write(f"#SBATCH --job-name=merge{i}\n")
                    f.write(f"#SBATCH --output={sd}/temp/job{i}.out\n")
                    f.write(f"#SBATCH --error={sd}/temp/job{i}.err\n")
                    f.write(f"#SBATCH --mem=64G\n")
                    # f.write(f"#SBATCH --partition=batch\n")
                    f.write(f"python -u 10_batch_merge_velocyto_layers.py {wd} {directory_name}\n")
                os.system(f"sbatch {sd}/temp/job{i}")
                i += 1
