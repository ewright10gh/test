import pandas as pd
from sklearn.metrics import classification_report, roc_auc_score

import matplotlib.pyplot as plt

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

    import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
import seaborn as sns

cm = confusion_matrix(y_true, y_pred)

plt.figure(figsize=(5,4))
sns.heatmap(
    cm,
    annot=True,
    fmt="d",
    cmap="Blues",
    xticklabels=["Non-fav","Fav"],
    yticklabels=["Non-fav","Fav"]
)

plt.xlabel("Predicted")
plt.ylabel("True")
plt.title("TCGA ELN Classification Confusion Matrix")

plt.tight_layout()
plt.savefig("tcga_confusion_matrix.png", dpi=300)
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, average_precision_score

# load predictions
results = pd.read_csv(r"C:\Users\mba22ew\test\results\tcga_eln_predictions.csv")

# load clinical
clinical = pd.read_csv(
    r"C:\Users\mba22ew\test\laml_tcga_pub_clinical_data.tsv",
    sep="\t",
    low_memory=False
)

clinical = clinical.set_index("Sample ID")

# align
results = results[results["Sample_ID"].isin(clinical.index)]
clinical = clinical.loc[results["Sample_ID"]]

# map fusions
def map_fusion(label):
    label = str(label).upper()
    if "15;17" in label: return 1
    if "8;21" in label: return 1
    if "INV(16)" in label: return 1
    return 0

y_true = clinical["Cytogenetics"].apply(map_fusion)
y_prob = results["ELN_probability"]

precision, recall, thresholds = precision_recall_curve(y_true, y_prob)
ap = average_precision_score(y_true, y_prob)

plt.figure(figsize=(6,5))

plt.plot(thresholds, precision[:-1], linewidth=2, label="Precision")
plt.plot(thresholds, recall[:-1], linewidth=2, label="Recall")

plt.axvline(0.7, linestyle="--", color="black", label="Chosen threshold")

plt.xlabel("Prediction probability threshold", fontsize=12)
plt.ylabel("Score", fontsize=12)

plt.title("Precision–Recall Trade-off", fontsize=13)

plt.legend(frameon=False)

plt.tight_layout()
plt.savefig("precision_recall_threshold.png", dpi=400)
plt.show()

from sklearn.metrics import roc_curve, roc_auc_score

fpr, tpr, _ = roc_curve(y_true, y_prob)
auc = roc_auc_score(y_true, y_prob)

plt.figure(figsize=(6,5))
plt.plot(fpr, tpr, label=f"AUC = {auc:.3f}")
plt.plot([0,1],[0,1],'--')

plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve")
plt.legend()

plt.tight_layout()
plt.savefig("tcga_roc_curve.png", dpi=300)
plt.show()
import seaborn as sns

plot_df = pd.DataFrame({
    "Probability": y_prob,
    "True_Label": y_true
})

plt.figure(figsize=(6,5))
sns.kdeplot(data=plot_df, x="Probability", hue="True_Label", fill=True)

plt.title("ELN Probability Distribution")
plt.xlabel("Predicted ELN Favourable Probability")

plt.tight_layout()
plt.savefig("tcga_probability_distribution.png", dpi=300)
plt.show()
plt.figure(figsize=(6,5))

plt.plot(thresholds, precision[:-1], label="Precision")
plt.plot(thresholds, recall[:-1], label="Recall")

plt.axvline(0.7, linestyle="--")

plt.xlabel("Probability Threshold")
plt.ylabel("Score")
plt.title("Precision & Recall vs Threshold")

plt.legend()

plt.tight_layout()
plt.savefig("precision_recall_vs_threshold.png", dpi=300)
plt.show()
from sklearn.calibration import calibration_curve

prob_true, prob_pred = calibration_curve(y_true, y_prob, n_bins=5)

plt.figure(figsize=(5,5))
plt.plot(prob_pred, prob_true, marker='o')
plt.plot([0,1],[0,1],'--')

plt.xlabel("Predicted probability")
plt.ylabel("Observed frequency")
plt.title("Calibration Curve")

plt.tight_layout()
plt.savefig("calibration_curve.png", dpi=300)
plt.show()