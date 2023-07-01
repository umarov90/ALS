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
for entry in os.scandir(wd):
    if entry.is_dir():
        directory_name = entry.name
        directory_path = entry.path
        if directory_name not in sample_ids:
            print(f"Skipping {directory_name}, not in samples.txt")
            continue
        if os.path.isfile(f"{wd}/{directory_name}/velocyto/{directory_name}.loom"):
            print(f"Skipping {directory_name}, loom already exists")
            continue
        if not os.path.isfile(f"{wd}/{directory_name}/outs/possorted_genome_bam.bam"):
            print(f"Skipping {directory_name}, bam is not there")
            continue
        dir = f"{wd}/{directory_name}/outs"
        for file in os.listdir(dir):
            if file.startswith('cellsorted'):
                os.remove(os.path.join(dir, file))
        with open(f"{sd}/temp/job{i}", 'w+') as f:
            f.write("#!/bin/bash\n")
            f.write(f"#SBATCH --job-name=velo{i}\n")
            f.write(f"#SBATCH --output={sd}/temp/job{i}.out\n")
            f.write(f"#SBATCH --error={sd}/temp/job{i}.err\n")
            f.write(f"#SBATCH --mem=128G\n")
            # f.write(f"#SBATCH --cpus-per-task=4\n")
            # f.write(f"#SBATCH --partition=batch\n")
            f.write(f"samtools sort -t CB -O BAM -m 100G -@ 1 -o {wd}/{directory_name}/outs/cellsorted_possorted_genome_bam.bam {wd}/{directory_name}/outs/possorted_genome_bam.bam\n")
            f.write(f"velocyto run10x -m {sd}/hg38_rmsk.gtf {wd}/{directory_name} {sd}/refdata-gex-GRCh38-2020-A/genes/genes.gtf")
        os.system(f"sbatch {sd}/temp/job{i}")
        i += 1
