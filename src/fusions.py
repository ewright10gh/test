import pandas as pd
import os

BASE = r"C:\Users\mba22ew\test"

expr = pd.read_csv(os.path.join(BASE,"cleaned","ohsu_cleaned_expression.csv"), index_col=0)

clinical = pd.read_csv(
    os.path.join(BASE,"aml_ohsu_2022_clinical_data.tsv"),
    sep="\t",
    low_memory=False
).set_index("Sample ID")

clinical = clinical.loc[expr.index]

fusion_col = clinical["Cancer Type Detailed"].str.upper()

clinical["fusion_group"] = "OTHER"

clinical.loc[fusion_col.str.contains("PML-RARA", na=False), "fusion_group"] = "PML_RARA"
clinical.loc[fusion_col.str.contains("RUNX1-RUNX1T1", na=False), "fusion_group"] = "RUNX1_RUNX1T1"
clinical.loc[fusion_col.str.contains("CBFB-MYH11", na=False), "fusion_group"] = "CBFB_MYH11"

print(clinical["fusion_group"].value_counts())

groups = {}

for fusion in ["PML_RARA","RUNX1_RUNX1T1","CBFB_MYH11"]:
    groups[fusion] = expr[clinical["fusion_group"] == fusion]

other = expr[clinical["fusion_group"] == "OTHER"]

results = {}

for fusion in groups:

    diff = groups[fusion].mean() - other.mean()

    top = diff.sort_values(ascending=False).head(20)

    results[fusion] = top

    print("\nTop genes for", fusion)
    print(top.head(10))

    top.to_csv(os.path.join(BASE,"results",f"{fusion}_top_genes.csv"))