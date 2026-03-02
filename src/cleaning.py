import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler

# -----------------------------
# 1. LOAD DATA
# -----------------------------
ohsu = pd.read_csv(r"C:\Users\mba22ew\test\data_mrna_seq_rpkm.txt", sep="\t", index_col=0)
target = pd.read_csv(r"C:\Users\mba22ew\test\data_mrna_seq_tpm.txt", sep="\t", index_col=0)
tcga = pd.read_csv(r"C:\Users\mba22ew\test\data_mrna_seq_v2_rsem.txt", sep="\t", index_col=0)
mapping = pd.read_csv(r"C:\Users\mba22ew\test\hgnc_complete_set.txt", sep="\t", low_memory=False)

print("Original shapes:")
print("OHSU:", ohsu.shape)
print("TARGET:", target.shape)
print("TCGA:", tcga.shape)

# -----------------------------
# 2. FIX ORIENTATION (samples in rows, genes in columns)
# -----------------------------
def transpose_if_needed(df):
    if df.shape[0] > df.shape[1]:
        return df.T
    return df

ohsu = transpose_if_needed(ohsu)
target = transpose_if_needed(target)
tcga = transpose_if_needed(tcga)

print("\nAfter transpose:")
print("OHSU:", ohsu.shape)
print("TARGET:", target.shape)
print("TCGA:", tcga.shape)

# -----------------------------
# 3. MAP ENTREZ → HUGO SYMBOL (TARGET dataset)
# -----------------------------
mapping = mapping[["symbol", "entrez_id"]].dropna()
mapping["entrez_id"] = mapping["entrez_id"].astype(float).astype(int).astype(str)

# Map TARGET columns
target.columns = target.columns.astype(str)
entrez_to_symbol = dict(zip(mapping["entrez_id"], mapping["symbol"]))
target.columns = target.columns.map(entrez_to_symbol)

# Drop unmapped genes
target = target.loc[:, target.columns.notna()]

# Collapse duplicates by averaging
target = target.T.groupby(level=0).mean().T

print("\nTARGET after gene mapping:", target.shape)

# -----------------------------
# 4. MATCH GENES ACROSS ALL DATASETS
# -----------------------------
# Make sure column names are comparable (strings, no surrounding whitespace)
ohsu.columns = ohsu.columns.astype(str).str.strip()
target.columns = target.columns.astype(str).str.strip()
tcga.columns = tcga.columns.astype(str).str.strip()

# collapse duplicate column names by averaging expression
# (TARGET was handled earlier but this is safe to call on all three)
def collapse_duplicates(df):
    if df.columns.duplicated().any():
        df = df.T.groupby(level=0).mean().T
    return df

ohsu = collapse_duplicates(ohsu)
target = collapse_duplicates(target)
tcga = collapse_duplicates(tcga)

# recompute intersection in case collapsing changed counts
common_genes = ohsu.columns.intersection(target.columns).intersection(tcga.columns)
print("\nCommon genes:", len(common_genes))

if len(common_genes) == 0:
    raise ValueError("❌ No common genes found — check gene naming")

# Subset datasets to common genes
ohsu = ohsu[common_genes]
target = target[common_genes]
tcga = tcga[common_genes]

# -----------------------------
# 5. FORCE NUMERIC AND DROP ALL-NaN GENES
# -----------------------------
def numeric_clean(df):
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(axis=1, how="all")
    return df

ohsu = numeric_clean(ohsu)
target = numeric_clean(target)
tcga = numeric_clean(tcga)

# report duplicates that might cause mismatched column numbers later
for name, df in [('OHSU', ohsu), ('TARGET', target), ('TCGA', tcga)]:
    total = len(df.columns)
    unique = df.columns.nunique()
    dup = df.columns[df.columns.duplicated()].unique()
    print(f"{name} columns: total={total}, unique={unique}, dup_count={len(dup)}")
    if len(dup) > 0:
        print(f"  duplicate names sample: {dup[:5]}")

# -----------------------------
# 6. LOG2 TRANSFORM TARGET & TCGA
# -----------------------------
target = np.log2(target + 1)
tcga = np.log2(tcga + 1)

# -----------------------------
# 7. VARIANCE FILTER (based on OHSU)
# -----------------------------
gene_variance = ohsu.var(axis=0)
high_var_genes = gene_variance[gene_variance > 1].index

# Intersect with other datasets to ensure all contain the same high-variance genes
high_var_genes = high_var_genes.intersection(target.columns).intersection(tcga.columns)
print("Genes passing variance filter in all datasets:", len(high_var_genes))

if len(high_var_genes) < 50:
    raise ValueError("❌ Too few genes after variance filter")

# Subset all datasets to high-variance genes
ohsu = ohsu[high_var_genes]
target = target[high_var_genes]
tcga = tcga[high_var_genes]

# --- sanity check: sometimes indexing quirks or dtypes can leave mismatches ---
common_after = ohsu.columns.intersection(target.columns).intersection(tcga.columns)
if len(common_after) != len(high_var_genes):
    print(f"WARNING: {len(high_var_genes) - len(common_after)} genes dropped when re-checking intersections")
    high_var_genes = common_after
    ohsu = ohsu[high_var_genes]
    target = target[high_var_genes]
    tcga = tcga[high_var_genes]

print("Genes retained after final intersection:", len(high_var_genes))

# -----------------------------
# 8. LOCK COLUMN ORDER (sort alphabetically)
# -----------------------------
ohsu = ohsu.sort_index(axis=1)
target = target[ohsu.columns]
tcga = tcga[ohsu.columns]

# Verify
print("\nAfter variance filter:")
print("OHSU:", ohsu.shape)
print("TARGET:", target.shape)
print("TCGA:", tcga.shape)
assert (ohsu.columns == target.columns).all() and (ohsu.columns == tcga.columns).all(), "Column mismatch!"

print("✅ Column order identical across all datasets")

# -----------------------------
# 9. SCALE
# -----------------------------
scaler = StandardScaler()
X_ohsu = scaler.fit_transform(ohsu)
X_target = scaler.transform(target)
X_tcga = scaler.transform(tcga)

# -----------------------------
# 10. SAVE
# -----------------------------
# ensure output directory exists
out_dir = "cleaned"
os.makedirs(out_dir, exist_ok=True)

np.save(os.path.join(out_dir, "X_ohsu.npy"), X_ohsu)
np.save(os.path.join(out_dir, "X_target.npy"), X_target)
np.save(os.path.join(out_dir, "X_tcga.npy"), X_tcga)

ohsu.to_csv(os.path.join(out_dir, "ohsu_cleaned_expression.csv"))
target.to_csv(os.path.join(out_dir, "target_cleaned_expression.csv"))
tcga.to_csv(os.path.join(out_dir, "tcga_cleaned_expression.csv"))

# -----------------------------
# 11. SUCCESS MESSAGE
# -----------------------------
print("\n✅ PREPROCESSING COMPLETE")
print("Final OHSU matrix:", X_ohsu.shape)
print("Final TARGET matrix:", X_target.shape)
print("Final TCGA matrix:", X_tcga.shape)
print("Gene count:", ohsu.shape[1])