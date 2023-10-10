import subprocess

python_files = ["2.0_obsm.py", "3_clustering.py", "2.1_de_genes_umap.py"]

for file in python_files:
    subprocess.run(["python", file])