import pandas as pd

# Load clinical file
clinical = pd.read_csv(
    "aml_ohsu_2022_clinical_data.tsv",   # ← change if needed
    sep="\t",
    low_memory=False
)

for col in clinical.columns:
    if "cancer" in col.lower() and "detailed" in col.lower():
        cancer_col = col

print("Using column:", cancer_col)


unique_values = clinical[cancer_col].dropna().unique()

print("\nTotal unique values:", len(unique_values))
print("\n=== ALL VALUES ===")
for v in sorted(unique_values):
    print(v)

print("\n=== VALUE COUNTS ===")
print(clinical[cancer_col].value_counts())


fusion_mask = clinical[cancer_col].str.contains(
    "fusion|t\\(|inv\\(|cbfb|runx1|kmt2a|pml", 
    case=False,
    na=False
)

print("\n=== FUSION-LIKE LABELS ===")
print(clinical.loc[fusion_mask, cancer_col].value_counts())