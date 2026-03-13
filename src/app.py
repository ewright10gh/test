import streamlit as st
import pandas as pd
import numpy as np
import joblib

st.set_page_config(page_title="AML ELN Predictor", layout="wide")

st.title("AML Transcriptomic Risk Classifier")

st.write(
"""
Upload a **TARGET-style RNA-seq TPM file** to predict ELN favourable AML risk.

The uploaded dataset replaces the TARGET dataset used in the analysis.
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
# LOAD OHSU CLEANED GENE LIST
# -----------------------------

ohsu_clean = pd.read_csv(
    "cleaned/ohsu_cleaned_expression.csv",
    index_col=0
)

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
    0.6
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

    target = np.log2(target + 1)


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
    pred = np.where(prob >= threshold, "Favourable", "Non-favourable")

    results = pd.DataFrame({

        "Sample_ID": target.index,
        "ELN_probability": prob,
        "ELN_predicted_class": pred

    })

    st.subheader("Prediction Results")

    st.dataframe(results)

    # -----------------------------
    # SUMMARY
    # -----------------------------

    st.subheader("Prediction counts")

    st.write(results["ELN_predicted_class"].value_counts())

    # -----------------------------
    # DOWNLOAD
    # -----------------------------

    csv = results.to_csv(index=False).encode()

    st.download_button(
        "Download predictions",
        csv,
        "eln_predictions.csv",
        "text/csv"
    )
    # ----------------------------------------
    # MODEL INTERPRETATION
    # ----------------------------------------

    st.subheader("Genes driving ELN predictions")

    import matplotlib.pyplot as plt

    # load ridge coefficients
    coefs = pd.read_csv(
        "cleaned/ridge_feature_importance.csv",
        index_col=0
    )

    coefs.columns = ["Coefficient"]

    # top genes
    top_fav = coefs.sort_values("Coefficient", ascending=False).head(5)
    top_adv = coefs.sort_values("Coefficient", ascending=True).head(5)

    interpret_table = pd.concat([top_fav, top_adv])

    interpret_table["Direction"] = [
        "Favourable"] * 5 + ["Non-favourable"] * 5

    interpret_table = interpret_table.reset_index()
    interpret_table.columns = ["Gene","Coefficient","Associated risk"]

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
    ax.set_title("Top genes driving AML ELN risk predictions")

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
