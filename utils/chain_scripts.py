import subprocess

python_files = ["7_obsm.py", "8_clustering.py", "9.1_deg_genes_umap.py"]

for file in python_files:
    subprocess.run(["python", file])