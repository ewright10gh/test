import os
import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

BASE = r"C:\Users\mba22ew\test"

# -------------------------------------------------
# FILE PATHS
# -------------------------------------------------

OHSU_EXPR = os.path.join(BASE,"cleaned","ohsu_cleaned_expression.csv")
OHSU_CLIN = os.path.join(BASE,"aml_ohsu_2022_clinical_data.tsv")

TCGA_EXPR = os.path.join(BASE,"cleaned","tcga_cleaned_expression.csv")
TCGA_PRED = os.path.join(BASE,"results","TCGA_favourable_fusion_predictions.csv")

TARGET_EXPR = os.path.join(BASE,"cleaned","target_cleaned_expression.csv")
TARGET_PRED = os.path.join(BASE,"results","TARGET_favourable_fusion_predictions.csv")

VALIDATA_EXPR = os.path.join(BASE,"cleaned","validata_cleaned_expression.csv")
VALIDATA_PRED = os.path.join(BASE,"results","validata_favourable_fusion_predictions.csv")

OUTDIR = os.path.join(BASE,"results","fusion_inference")
os.makedirs(OUTDIR,exist_ok=True)

# -------------------------------------------------
# LOAD DATA
# -------------------------------------------------

print("Loading OHSU expression...")
ohsu_expr = pd.read_csv(OHSU_EXPR,index_col=0)

print("Loading OHSU clinical...")
ohsu_clin = pd.read_csv(OHSU_CLIN,sep="\t",low_memory=False)
ohsu_clin = ohsu_clin.set_index("Sample ID")

print("Loading TCGA expression...")
tcga_expr = pd.read_csv(TCGA_EXPR,index_col=0)

print("Loading TCGA predictions...")
tcga_pred = pd.read_csv(TCGA_PRED,index_col=0)

print("Loading TARGET expression...")
target_expr = pd.read_csv(TARGET_EXPR,index_col=0)

print("Loading TARGET predictions...")
target_pred = pd.read_csv(TARGET_PRED,index_col=0)

print("Loading VALIDATA expression...")
validata_expr = pd.read_csv(VALIDATA_EXPR,index_col=0)

print("Loading VALIDATA predictions...")
validata_pred = pd.read_csv(VALIDATA_PRED,index_col=0)

# -------------------------------------------------
# ALIGN OHSU
# -------------------------------------------------

ohsu_clin = ohsu_clin.loc[ohsu_expr.index]

fusion_series = ohsu_clin["Cancer Type Detailed"].astype(str).str.upper()

fusion_group = pd.Series("OTHER",index=ohsu_expr.index)

fusion_group[fusion_series.str.contains("PML-RARA",na=False)] = "PML_RARA"
fusion_group[fusion_series.str.contains("RUNX1-RUNX1T1",na=False)] = "RUNX1_RUNX1T1"
fusion_group[fusion_series.str.contains("CBFB-MYH11",na=False)] = "CBFB_MYH11"

print("\nOHSU Fusion counts:")
print(fusion_group.value_counts())

# -------------------------------------------------
# COMPUTE FUSION CENTROIDS
# -------------------------------------------------

centroids = {}

for fusion in ["PML_RARA","RUNX1_RUNX1T1","CBFB_MYH11"]:

    subset = ohsu_expr[fusion_group==fusion]

    if len(subset)==0:
        continue

    centroids[fusion] = subset.mean()

centroids = pd.DataFrame(centroids)

print("\nCentroid matrix shape:",centroids.shape)

# -------------------------------------------------
# FUNCTION TO INFER FUSION PROGRAM
# -------------------------------------------------

def infer_fusion(expr, pred, dataset_name):

    print("\nProcessing", dataset_name)

    # -----------------------------
    # SET SAMPLE ID
    # -----------------------------

    if "Sample_ID" in pred.columns:
        pred = pred.set_index("Sample_ID")

    # -----------------------------
    # FIND PROBABILITY COLUMN
    # -----------------------------

    prob_col = None

    for c in pred.columns:
        if "prob" in c.lower():
            prob_col = c
            break

    if prob_col is None:
        raise ValueError("No probability column found")

    print("Using probability column:", prob_col)

    # -----------------------------
    # ALIGN SAMPLES
    # -----------------------------

    common_samples = expr.index.intersection(pred.index)

    expr = expr.loc[common_samples]
    pred = pred.loc[common_samples]

    print("Aligned samples:", len(common_samples))

    # -----------------------------
    # SELECT FAVOURABLE SAMPLES
    # -----------------------------

    fav_mask = pred[prob_col] > 0.675
    fav_expr = expr.loc[fav_mask]

    print(dataset_name, "favourable samples:", fav_expr.shape[0])

    if fav_expr.shape[0] == 0:
        print("No favourable samples found")
        return None

    # -----------------------------
    # ALIGN GENES
    # -----------------------------

    common_genes = fav_expr.columns.intersection(centroids.index)

    fav_expr = fav_expr[common_genes]
    cent = centroids.loc[common_genes]

    print("Genes used:", len(common_genes))

    # -----------------------------
    # COSINE SIMILARITY
    # -----------------------------

    sim = cosine_similarity(fav_expr, cent.T)

    sim_df = pd.DataFrame(
        sim,
        index=fav_expr.index,
        columns=["PML_RARA","RUNX1_RUNX1T1","CBFB_MYH11"]
    )

    # -----------------------------
    # ASSIGN FUSION PROGRAM
    # -----------------------------

    sim_df["Predicted_fusion_program"] = sim_df.idxmax(axis=1)

    # -----------------------------
    # SAVE TABLE
    # -----------------------------

    outfile = os.path.join(
        OUTDIR,
        f"{dataset_name}_fusion_inference.csv"
    )

    sim_df.to_csv(outfile)

    print("Saved:", outfile)

    print("\nFusion distribution:")
    print(sim_df["Predicted_fusion_program"].value_counts())

    return sim_df

# -------------------------------------------------
# RUN ANALYSIS
# -------------------------------------------------

tcga_results = infer_fusion(tcga_expr,tcga_pred,"TCGA")

target_results = infer_fusion(target_expr,target_pred,"TARGET")

validata_results = infer_fusion(validata_expr,validata_pred,"VALIDATA")

print("\nDone.")