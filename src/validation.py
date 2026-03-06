import pandas as pd
from sklearn.metrics import classification_report, roc_auc_score

# =========================

# 1. LOAD MODEL PREDICTIONS

# =========================

results = pd.read_csv(r"C:\Users\mba22ew\test\results\tcga_eln_predictions.csv")

print("Predictions loaded:", results.shape)

# =========================

# 2. LOAD TCGA CLINICAL

# =========================

clinical = pd.read_csv(
r"C:\Users\mba22ew\test\laml_tcga_pub_clinical_data.tsv",
sep="\t",
low_memory=False
)

print("Clinical loaded:", clinical.shape)

# TCGA sample IDs

clinical = clinical.set_index("Sample ID")

# =========================

# 3. ALIGN SAMPLES

# =========================

common_samples = results["Sample_ID"].astype(str).isin(clinical.index)

results = results[common_samples]
clinical = clinical.loc[results["Sample_ID"]]

print("Samples overlapping:", len(results))

# =========================

# 4. MAP CYTOGENETICS → ELN

# =========================

def map_fusion(label):


    label = str(label).upper()

    if "15;17" in label:
        return 1

    if "8;21" in label:
        return 1

    if "INV(16)" in label:
        return 1

    return 0


y_true = clinical["Cytogenetics"].apply(map_fusion)
y_pred = results["ELN_predicted"]
y_prob = results["ELN_probability"]

print("\nTrue ELN distribution:")
print(y_true.value_counts())

# =========================

# 5. EVALUATION

# =========================

print("\n=== CLASSIFICATION REPORT ===")
print(classification_report(y_true, y_pred))

if y_true.nunique() > 1:


    roc = roc_auc_score(y_true, y_prob)
    print("\nROC AUC:", roc)


else:
    print("\nROC AUC cannot be computed (only one class present)")
