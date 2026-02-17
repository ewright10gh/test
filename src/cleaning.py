import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler

# =========================
# 1. LOAD DATA
# =========================

ohsu = pd.read_csv(
    r"C:\Users\mba22ew\test\data_mrna_seq_rpkm.txt",
    sep="\t",
    index_col=0
)

target = pd.read_csv(
    r"C:\Users\mba22ew\test\data_mrna_seq_tpm.txt",
    sep="\t",
    index_col=0
)

mapping = pd.read_csv(
    r"C:\Users\mba22ew\test\hgnc_complete_set.txt",
    sep="\t",
    low_memory=False
)

print("Original shapes:")
print("OHSU:", ohsu.shape)
print("TARGET:", target.shape)

# =========================
# 2. FIX ORIENTATION
# =========================

if ohsu.shape[0] > ohsu.shape[1]:
    print("Transposing OHSU")
    ohsu = ohsu.T

if target.shape[0] > target.shape[1]:
    print("Transposing TARGET")
    target = target.T

print("\nAfter transpose:")
print("OHSU:", ohsu.shape)
print("TARGET:", target.shape)

# =========================
# 3. MAP ENTREZ → HUGO (robust)
# =========================

mapping = mapping[["symbol", "entrez_id"]].dropna()

# Convert mapping Entrez IDs to int → str
mapping["entrez_id"] = mapping["entrez_id"].astype(float).astype(int).astype(str)

# Convert TARGET columns to str
target.columns = target.columns.astype(str)

# Create mapping dictionary
entrez_to_symbol = dict(zip(mapping["entrez_id"], mapping["symbol"]))

# Map TARGET columns
target.columns = target.columns.map(entrez_to_symbol)

# Drop unmapped genes
target = target.loc[:, target.columns.notna()]

# Collapse duplicate gene symbols by averaging
target = target.T.groupby(level=0).mean().T

print("\nTARGET after gene mapping:", target.shape)

# =========================
# 4. MATCH GENES
# =========================

common_genes = ohsu.columns.intersection(target.columns)

print("\nCommon genes:", len(common_genes))

if len(common_genes) == 0:
    raise ValueError("❌ No common genes found — check gene naming")

ohsu = ohsu[common_genes]
target = target[common_genes]

# =========================
# 5. FORCE NUMERIC
# =========================

ohsu = ohsu.apply(pd.to_numeric, errors="coerce")
target = target.apply(pd.to_numeric, errors="coerce")

# Drop genes with all NaN
ohsu = ohsu.dropna(axis=1, how="all")
target = target.dropna(axis=1, how="all")

# =========================
# 6. LOG2 TRANSFORM TARGET
# =========================

target = np.log2(target + 1)

# =========================
# 7. VARIANCE FILTER
# =========================

gene_variance = ohsu.var(axis=0)

high_var_genes = gene_variance[gene_variance > 1].index
high_var_genes = high_var_genes.intersection(target.columns)

print("Genes passing variance filter:", len(high_var_genes))

if len(high_var_genes) < 50:
    raise ValueError("❌ Too few genes after variance filter")

ohsu = ohsu[high_var_genes]
target = target[high_var_genes]

print("\nAfter variance filter:")
print("OHSU:", ohsu.shape)
print("TARGET:", target.shape)

# =========================
# 8. LOCK COLUMN ORDER
# =========================

target = target[ohsu.columns]
print("Column order identical:", (ohsu.columns == target.columns).all())

# =========================
# 9. SCALE
# =========================

scaler = StandardScaler()

X_ohsu = scaler.fit_transform(ohsu)
X_target = scaler.transform(target)

# =========================
# 10. SAVE
# =========================

np.save("X_ohsu.npy", X_ohsu)
np.save("X_target.npy", X_target)

ohsu.to_csv("ohsu_cleaned_expression.csv")
target.to_csv("target_cleaned_expression.csv")

# =========================
# 11. SUCCESS MESSAGE
# =========================

print("\n✅ PREPROCESSING COMPLETE")
print("Final OHSU matrix:", X_ohsu.shape)
print("Final TARGET matrix:", X_target.shape)
print("Gene count:", ohsu.shape[1])