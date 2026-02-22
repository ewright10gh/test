import pandas as pd

# Load clinical file
clinical = pd.read_csv(
    "aml_ohsu_2022_clinical_data.tsv",   # ← change if needed
    sep="\t",
    low_memory=False
)

# -----------------------------
# 1️⃣ Find the Cancer Type column safely
# -----------------------------
for col in clinical.columns:
    if "cancer" in col.lower() and "detailed" in col.lower():
        cancer_col = col

print("Using column:", cancer_col)

# -----------------------------
# 2️⃣ Show ALL unique values
# -----------------------------
unique_values = clinical[cancer_col].dropna().unique()

print("\nTotal unique values:", len(unique_values))
print("\n=== ALL VALUES ===")
for v in sorted(unique_values):
    print(v)

# -----------------------------
# 3️⃣ Show counts for each
# -----------------------------
print("\n=== VALUE COUNTS ===")
print(clinical[cancer_col].value_counts())

# -----------------------------
# 4️⃣ Show only fusion-related ones
# -----------------------------
fusion_mask = clinical[cancer_col].str.contains(
    "fusion|t\\(|inv\\(|cbfb|runx1|kmt2a|pml", 
    case=False,
    na=False
)

print("\n=== FUSION-LIKE LABELS ===")
print(clinical.loc[fusion_mask, cancer_col].value_counts())