import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler

ohsu = pd.read_csv(
    r"C:\Users\mba22ew\test\data_mrna_seq_rpkm.txt",  # Windows path
    sep="\t",
    index_col=0
)

target = pd.read_csv(
    r"C:\Users\mba22ew\test\data_mrna_seq_tpm.txt",  # Windows path
    sep="\t",
    index_col=0
)
mapping = pd.read_csv(
    r"C:\Users\mba22ew\test\entrez_to_hugo.tsv",
    sep="\t"
)


print("Original shapes:")
print("OHSU:", ohsu.shape)
print("TARGET:", target.shape)


# =========================
# 2. CHECK ORIENTATION
# =========================
# We expect: samples × genes
# If genes are rows, transpose

if ohsu.shape[0] > ohsu.shape[1]:
    print("Transposing OHSU")
    ohsu = ohsu.T

if target.shape[0] > target.shape[1]:
    print("Transposing TARGET")
    target = target.T


print("\nAfter orientation fix:")
print("OHSU:", ohsu.shape)
print("TARGET:", target.shape)
print("OHSU genes:", ohsu.columns[:5])
print("TARGET genes:", target.columns[:5])

# -----------------------------
# Map Entrez → Hugo
# -----------------------------
mapping = mapping.dropna()
mapping = mapping.drop_duplicates("Entrez_Gene_Id")

target["Hugo_Symbol"] = target.index.map(
    mapping.set_index("Entrez_Gene_Id")["Hugo_Symbol"]
)

target = target.dropna(subset=["Hugo_Symbol"])
target = target.set_index("Hugo_Symbol")

# -----------------------------
# Remove duplicate genes
# -----------------------------
target = target.groupby(target.index).mean()

# =========================
# 3. MATCH GENES (COLUMNS)
# =========================

common_genes = ohsu.columns.intersection(target.columns)

print("\nCommon genes:", len(common_genes))

ohsu = ohsu[common_genes]
target = target[common_genes]


# =========================
# 4. LOG2 TRANSFORM TARGET ONLY
# =========================

target = np.log2(target + 1)


# =========================
# 5. REMOVE LOW-VARIANCE GENES (OPTIONAL BUT RECOMMENDED)
# =========================

gene_variance = ohsu.var(axis=0)
high_var_genes = gene_variance[gene_variance > 1].index

ohsu = ohsu[high_var_genes]
target = target[high_var_genes]

print("\nAfter variance filter:")
print("OHSU:", ohsu.shape)
print("TARGET:", target.shape)


# =========================
# 6. SCALE USING OHSU ONLY
# =========================

scaler = StandardScaler()

X_ohsu = scaler.fit_transform(ohsu)
X_target = scaler.transform(target)


# =========================
# 7. SAVE CLEAN MATRICES
# =========================

np.save("X_ohsu.npy", X_ohsu)
np.save("X_target.npy", X_target)

ohsu.to_csv("ohsu_cleaned_expression.csv")
target.to_csv("target_cleaned_expression.csv")

print("\nPreprocessing complete ✅")