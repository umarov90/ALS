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
        bam_path = f"{wd}/{directory_name}/outs/possorted_genome_bam.bam"
        if os.path.isfile(bam_path):
            print(f"{directory_name}\t{bam_path}")
