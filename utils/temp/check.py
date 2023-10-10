import joblib

a = joblib.load("data_split.gz")
c = joblib.load("hg38_regions.gz")
b = joblib.load("train_info.gz")
a = a["hg38"]["train"]
for i in range(len(a)):
    if a[i][0] != b[i][0] or a[i][1] != b[i][1]:
        print("Error")
        break
print(len(a))