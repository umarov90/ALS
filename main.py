import subprocess

python_files = ["5_sub_obsm.py", "6_treatment_sub_obsm.py", "7_clustering.py"]

for file in python_files:
    subprocess.run(["python", file])