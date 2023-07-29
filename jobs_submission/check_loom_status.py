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
            # print(f"Skipping {directory_name}, not in samples.txt")
            continue
        bam_path = f"{wd}/{directory_name}/outs/possorted_genome_bam.bam"
        if os.path.isfile(bam_path):
            print(f"{directory_name}\t{os.path.getsize(bam_path)/(1024 * 1024)}", end="\t")
            c_bam_path = f"{wd}/{directory_name}/outs/cellsorted_possorted_genome_bam.bam"
            if os.path.isfile(c_bam_path):
                print(f"{os.path.getsize(c_bam_path)/(1024 * 1024)}", end="\t")
            loom_path = f"{wd}/{directory_name}/velocyto/{directory_name}.loom"
            if os.path.isfile(loom_path):
                print(f"{os.path.getsize(loom_path) / (1024 * 1024)}", end="\t")
            print("")
