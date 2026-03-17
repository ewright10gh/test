import streamlit as st
import pandas as pd
import numpy as np
import joblib
from sklearn.metrics.pairwise import cosine_similarity

st.set_page_config(page_title="AML Favourable Fusion Predictor", layout="wide")

st.title("AML Transcriptomic Favourable Fusion Classifier")

st.write(
"""
Upload a **TARGET/TCGA/OHSU RNA-seq expression matrix** to predict favourable fusion vs no favourable fusion.

Supported formats:
- TARGET TPM (Entrez IDs)
- TCGA RSEM (gene|entrez union IDs)
- OHSU RPKM (HUGO symbols)
"""
)

# -----------------------------
# LOAD TRAINED MODEL
# -----------------------------

@st.cache_resource
def load_model():

    model = joblib.load("models/ridge_eln_model.joblib")
    scaler = joblib.load("models/ridge_scaler.joblib")

    return model, scaler

model, scaler = load_model()

# -----------------------------
# LOAD OHSU CLEANED GENE LIST + FUSION CENTROIDS
# -----------------------------

ohsu_clean = pd.read_csv(
    "cleaned/ohsu_cleaned_expression.csv",
    index_col=0
)

# Build fusion centroids using OHSU clinical labels
ohsu_clin = pd.read_csv(
    "aml_ohsu_2022_clinical_data.tsv",
    sep="\t",
    low_memory=False
).set_index("Sample ID")

ohsu_clin = ohsu_clin.loc[ohsu_clean.index]

fusion_series = ohsu_clin["Cancer Type Detailed"].astype(str).str.upper()

fusion_group = pd.Series("OTHER", index=ohsu_clean.index)
fusion_group[fusion_series.str.contains("PML-RARA", na=False)] = "PML_RARA"
fusion_group[fusion_series.str.contains("RUNX1-RUNX1T1", na=False)] = "RUNX1_RUNX1T1"
fusion_group[fusion_series.str.contains("CBFB-MYH11", na=False)] = "CBFB_MYH11"

centroid_dict = {}
for fusion in ["PML_RARA", "RUNX1_RUNX1T1", "CBFB_MYH11"]:
    subset = ohsu_clean[fusion_group == fusion]
    if len(subset):
        centroid_dict[fusion] = subset.mean()

centroids = pd.DataFrame(centroid_dict)

model_genes = ohsu_clean.columns

# -----------------------------
# FILE UPLOAD
# -----------------------------

uploaded = st.file_uploader(
    "Upload TARGET RNA-seq file",
    type=["txt","tsv"]
)

threshold = st.slider(
    "Favourable probability threshold",
    0.5,
    0.9,
    0.66
)

if uploaded is not None:

    # -----------------------------
    # LOAD USER FILE
    # -----------------------------

    target = pd.read_csv(uploaded, sep="\t", index_col=0)

    st.write("Original shape:", target.shape)


    # -----------------------------
    # FIX ORIENTATION
    # -----------------------------

    # TARGET / TCGA files usually have genes as rows
    # We want genes as columns

    if target.shape[0] > target.shape[1]:
        target = target.T

    st.write("After transpose:", target.shape)


    # -----------------------------
    # AUTO DETECT DATASET TYPE
    # -----------------------------

    sample_genes = list(target.columns[:100])

    pipe_count = sum(["|" in str(g) for g in sample_genes])

    def is_number(x):
        try:
            float(x)
            return True
        except:
            return False

    numeric_count = sum([is_number(g) for g in sample_genes])

    if pipe_count > 20:
        dataset_type = "TCGA"

    elif numeric_count > 50:
        dataset_type = "TARGET"

    else:
        dataset_type = "OHSU"  # includes RPKM/OHSU-style gene symbols

    # Explicit fallback for Hugo_Symbol index (RPKM data_mrna_seq_rpkm style)
    if target.index.name and "hugo" in str(target.index.name).lower():
        dataset_type = "OHSU"

    st.write("Detected dataset type:", dataset_type)
    st.write("Example genes:", sample_genes[:10])

    # -----------------------------
    # DATASET-SPECIFIC PROCESSING
    # -----------------------------

    if dataset_type == "TCGA":

        # TCGA format: TP53|7157
        target.columns = target.columns.str.split("|").str[0]
        target = target.loc[:, target.columns != ""]

        # collapse duplicates
        target = target.T.groupby(level=0).mean().T


    elif dataset_type == "TARGET":

        # load HGNC mapping
        mapping = pd.read_csv(
            "hgnc_complete_set.txt",
            sep="\t",
            low_memory=False
        )

        mapping = mapping[["symbol","entrez_id"]].dropna()
        mapping["entrez_id"] = mapping["entrez_id"].astype(int).astype(str)

        entrez_to_symbol = dict(zip(mapping["entrez_id"], mapping["symbol"]))

        target.columns = target.columns.astype(str)
        target.columns = target.columns.map(entrez_to_symbol)

        # remove unmapped genes
        target = target.loc[:, target.columns.notna()]

        # collapse duplicates
        target = target.T.groupby(level=0).mean().T


    elif dataset_type == "OHSU":

        # already gene symbols
        target.columns = target.columns.astype(str).str.strip()


    st.write("After gene harmonisation:", target.shape)


    # -----------------------------
    # NUMERIC CLEAN
    # -----------------------------

    target = target.apply(pd.to_numeric, errors="coerce")
    target = target.dropna(axis=1, how="all")


    # -----------------------------
    # LOG TRANSFORM
    # -----------------------------

    if dataset_type in ["TARGET", "TCGA"]:
        target = np.log2(target + 1)
        st.write("Applied log2(x+1) transform (TARGET/TCGA)")
    else:
        st.write("No log transform applied (OHSU/RPKM)")


    # -----------------------------
    # MATCH MODEL GENES
    # -----------------------------

    common_genes = target.columns.intersection(model_genes)

    st.write("Genes overlapping with model:", len(common_genes))

    if len(common_genes) < 500:
        st.error("Too few genes overlap with the training dataset.")
        st.stop()

    target = target.loc[:, common_genes]

    # -----------------------------
    # REMOVE DUPLICATE GENES
    # -----------------------------

    if target.columns.duplicated().any():
        st.write("Duplicate genes detected — collapsing")
        target = target.T.groupby(level=0).mean().T

    # -----------------------------
    # ALIGN FEATURE SPACE
    # -----------------------------

    target = target.reindex(columns=model_genes, fill_value=0)

    # -----------------------------
    # HANDLE NaNs USING TRAINING MEANS
    # -----------------------------
    
    nan_count = target.isna().sum().sum()
    
    if nan_count > 0:
    
        st.write(f"NaNs detected: {nan_count}")
    
        # replace NaNs with the training mean expression for each gene
        gene_means = pd.Series(scaler.mean_, index=model_genes)
    
        target = target.fillna(gene_means)
    
        st.write("NaNs replaced with training gene means")
    # -----------------------------
    # SCALE
    # -----------------------------

    X = scaler.transform(target.values)

    # -----------------------------
    # PREDICT
    # -----------------------------

    prob = model.predict_proba(X)[:,1]
    pred = np.where(prob >= threshold, "Favourable fusion", "No favourable fusion")

    results = pd.DataFrame({

        "Sample_ID": target.index,
        "favourable_fusion_probability": prob,
        "favourable_fusion_predicted_class": pred

    })

    # -------------------------------------------------
    # FUSION INFERENCE FOR FAVOURABLE SAMPLES
    # -------------------------------------------------

    fusion_program = pd.Series("No favourable fusion", index=results.index)

    fav_mask = results["favourable_fusion_predicted_class"] == "Favourable fusion"
    fav_samples = results.loc[fav_mask, "Sample_ID"]

    if len(fav_samples) > 0 and not centroids.empty:
        fav_expr = target.loc[fav_samples]

        common_genes = fav_expr.columns.intersection(centroids.index)
        if len(common_genes) > 0:
            fav_expr = fav_expr[common_genes]
            cent = centroids.loc[common_genes]

            sim = cosine_similarity(fav_expr, cent.T)
            sim_df = pd.DataFrame(
                sim,
                index=fav_expr.index,
                columns=centroids.columns
            )

            pred_fusion = sim_df.idxmax(axis=1)
            fusion_program.loc[fav_mask] = pred_fusion.values
        else:
            st.warning("No overlapping genes between upload and fusion centroids; fusion program inference skipped.")
    else:
        if len(fav_samples) == 0:
            st.info("No favourable-fusion samples were predicted, so fusion-type inference is skipped.")

    results["favourable_fusion_program"] = fusion_program

    st.subheader("Prediction Results")

    st.dataframe(results)

    # -----------------------------
    # SUMMARY
    # -----------------------------

    st.subheader("Prediction counts")

    pred_counts = results["favourable_fusion_predicted_class"].value_counts()
    st.write(pred_counts)

    # Show total favourable fusion count and fusion-type breakdown
    fav_count = int(pred_counts.get("Favourable fusion", 0))
    st.write(f"Total favourable fusion samples: {fav_count}")

    fusion_counts = results.loc[results["favourable_fusion_program"] != "No favourable fusion", "favourable_fusion_program"].value_counts()
    if len(fusion_counts) > 0:
        st.write("Favourable fusion breakdown by inferred program:")
        st.write(fusion_counts)
    else:
        st.info("No fusion program inferred for favourable fusion samples.")

    # -----------------------------
    # DOWNLOAD
    # -----------------------------

    csv = results.to_csv(index=False).encode()

    st.download_button(
        "Download predictions",
        csv,
        "favourable_fusion_predictions.csv",
        "text/csv"
    )
    # ----------------------------------------
    # MODEL INTERPRETATION
    # ----------------------------------------

    st.subheader("Genes driving favourable fusion predictions (dataset-specific)")

    import matplotlib.pyplot as plt

    # load ridge coefficients (global model)
    coefs = pd.read_csv(
        "cleaned/ridge_feature_importance.csv",
        index_col=0
    )

    coefs.columns = ["Coefficient"]

    # compute dataset-specific mean expression stratified by prediction
    fav_samples = results[results["favourable_fusion_predicted_class"] == "Favourable fusion"]["Sample_ID"]
    nonfav_samples = results[results["favourable_fusion_predicted_class"] == "No favourable fusion"]["Sample_ID"]

    fav_expr = target.loc[fav_samples]
    nonfav_expr = target.loc[nonfav_samples]

    # handle empty classes safely
    if fav_expr.shape[0] == 0 or nonfav_expr.shape[0] == 0:
        st.warning("Not enough class diversity to compute dataset-specific top genes.")
        fav_expr = target
        nonfav_expr = target

    stats = pd.DataFrame({
        "Favourable_mean": fav_expr.mean(),
        "No_favourable_mean": nonfav_expr.mean(),
    })
    stats["mean_diff"] = stats["Favourable_mean"] - stats["No_favourable_mean"]

    # join with coefficients
    stats = stats.join(coefs, how="left")

    top_fav = stats.sort_values("mean_diff", ascending=False).head(5).reset_index()
    top_nonfav = stats.sort_values("mean_diff", ascending=True).head(5).reset_index()

    interpret_table = pd.concat([top_fav, top_nonfav], ignore_index=True)
    interpret_table["Associated risk"] = ["Favourable fusion"] * len(top_fav) + ["No favourable fusion"] * len(top_nonfav)

    interpret_table = interpret_table[["index", "Coefficient", "mean_diff", "Associated risk"]]
    interpret_table.columns = ["Gene", "Coefficient", "Mean expression difference (fav - nonfav)", "Associated risk"]

    # show table
    st.dataframe(interpret_table)

    # ----------------------------------------
    # BARPLOT
    # ----------------------------------------

    fig, ax = plt.subplots(figsize=(7,4))

    colors = ["#2ca02c" if c > 0 else "#d62728"
            for c in interpret_table["Coefficient"]]

    ax.barh(
        interpret_table["Gene"],
        interpret_table["Coefficient"],
        color=colors
    )

    ax.axvline(0, color="black", linewidth=1)

    ax.set_xlabel("Model coefficient")
    ax.set_title("Top genes driving favourable fusion predictions")

    plt.tight_layout()

    st.pyplot(fig)

else:
    st.info("Please upload a file to get started.")

# -----------------------------
# SIDEBAR INFO
# -----------------------------

st.sidebar.write(
"""
### Expected input

A **TARGET TPM expression matrix**

Format:

- tab separated
- rows = samples
- columns = genes
- values = TPM

Example file:
""")
