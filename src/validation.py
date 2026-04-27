import pandas as pd
import numpy as np
from sklearn.metrics import classification_report, roc_auc_score
import os
import matplotlib.pyplot as plt
import seaborn as sns

# Get the directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)



# 1. LOAD MODEL PREDICTIONS



results = pd.read_csv(os.path.join(parent_dir, "results", "tcga_favourable_fusion_predictions.csv"))

print("Predictions loaded:", results.shape)



# 2. LOAD TCGA CLINICAL



clinical = pd.read_csv(
    os.path.join(parent_dir, "laml_tcga_pub_clinical_data.tsv"),
    sep="\t",
    low_memory=False
)

print("Clinical loaded:", clinical.shape)

# TCGA sample IDs

clinical = clinical.set_index("Sample ID")



# 3. ALIGN SAMPLES



common_samples = results["Sample_ID"].astype(str).isin(clinical.index)

results = results[common_samples]
clinical = clinical.loc[results["Sample_ID"]]

print("Samples overlapping:", len(results))



# 4. MAP CYTOGENETICS → ELN



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
y_pred = results["favourable_fusion_predicted"]
y_prob = results["favourable_fusion_probability"]

print("\nTrue favourable fusion distribution:")
print(y_true.value_counts())



# 5. EVALUATION



print("\n=== CLASSIFICATION REPORT ===")
print(classification_report(y_true, y_pred))

if y_true.nunique() > 1:


    roc = roc_auc_score(y_true, y_prob)
    print("\nROC AUC:", roc)


else:
    print("\nROC AUC cannot be computed (only one class present)")



# VALIDATA VALIDATION


validata_results = pd.read_csv(os.path.join(parent_dir, "results", "validata_favourable_fusion_predictions.csv"))
validata_clinical = pd.read_csv(
    os.path.join(parent_dir, "aml_ohsu_2018_clinical_data.tsv"),
    sep="\t",
    low_memory=False
)

if "sampleId" in validata_clinical.columns:
    validata_clinical = validata_clinical.set_index("sampleId")
else:
    validata_clinical = validata_clinical.set_index("Sample ID")

common_samples = validata_results["Sample_ID"].astype(str).isin(validata_clinical.index.astype(str))
validata_results = validata_results[common_samples]
validata_clinical = validata_clinical.loc[validata_results["Sample_ID"].astype(str)]

print("Validata samples overlapping:", len(validata_results))


def map_fusion_ohsu(label):
    label = str(label).upper()
    if "PML-RARA" in label or "15;17" in label:
        return 1
    if "RUNX1-RUNX1T1" in label or "8;21" in label:
        return 1
    if "CBFB-MYH11" in label or "INV(16)" in label:
        return 1
    return 0

if "CANCER_TYPE_DETAILED" in validata_clinical.columns:
    y_true_validata = validata_clinical["CANCER_TYPE_DETAILED"].apply(map_fusion_ohsu)
elif "CANCER_TYPE" in validata_clinical.columns:
    y_true_validata = validata_clinical["CANCER_TYPE"].apply(map_fusion_ohsu)
else:
    y_true_validata = pd.Series(0, index=validata_results.index)


y_pred_validata = validata_results["favourable_fusion_predicted"]
y_prob_validata = validata_results["favourable_fusion_probability"]

print("\nValidata true distribution:")
print(y_true_validata.value_counts())

print("\n=== VALIDATA CLASSIFICATION REPORT ===")
print(classification_report(y_true_validata, y_pred_validata))

if y_true_validata.nunique() > 1:
    roc_validata = roc_auc_score(y_true_validata, y_prob_validata)
    print("\nValidata ROC AUC:", roc_validata)
else:
    print("\nValidata ROC AUC cannot be computed (only one class present)")

# Load TARGET predictions for combined probability distributions
try:
    target_results = pd.read_csv(os.path.join(parent_dir, "results", "TARGET_favourable_fusion_predictions.csv"))
    target_results["favourable_fusion_probability"] = pd.to_numeric(
        target_results["favourable_fusion_probability"], errors="coerce"
    )
    target_results = target_results.dropna(subset=["favourable_fusion_probability"])
    target_probs = target_results["favourable_fusion_probability"].copy()
except FileNotFoundError:
    target_results = None
    target_probs = pd.Series(dtype=float)
    print("Warning: TARGET prediction file not found. Combined distribution figure will exclude TARGET.")

# Combine probabilities
plot_df = pd.DataFrame({
    "Probability": np.concatenate([
        y_prob.dropna().to_numpy(),
        y_prob_validata.dropna().to_numpy(),
        target_probs.to_numpy()
    ]),
    "Dataset": (
        ["TCGA"] * len(y_prob.dropna()) +
        ["OHSU 2018"] * len(y_prob_validata.dropna()) +
        ["TARGET"] * len(target_probs)
    )
}).dropna()

# Sanity check: proportions above threshold
threshold = 0.675

for ds in ["TCGA", "OHSU 2018", "TARGET"]:
    p = plot_df.loc[
        plot_df["Dataset"] == ds,
        "Probability"
    ]
    print(f"{ds}: {(p > threshold).mean():.3f} above threshold")

# Consistent palette
palette = {
    "TARGET": "#2C7FB8",
    "TCGA": "#D95F02",
    "OHSU 2018": "#FFC0CB"
}

# -----------------------------------------
# Two-panel figure
# -----------------------------------------

fig, axes = plt.subplots(
    1, 2,
    figsize=(12,5)
)


# Left panel: ECDF


sns.ecdfplot(
    data=plot_df,
    x="Probability",
    hue="Dataset",
    hue_order=["TARGET","TCGA","OHSU 2018"],
    palette=palette,
    linewidth=2,
    ax=axes[0]
)

axes[0].axvline(
    threshold,
    linestyle="--",
    color="black",
    linewidth=1.5
)

axes[0].annotate(
    "Threshold (0.675)",
    xy=(threshold,0.75),
    xytext=(0.76,0.60),
    arrowprops=dict(arrowstyle="->", lw=1.2),
    fontsize=10
)

axes[0].set_title("Cumulative probability distributions")
axes[0].set_xlabel("Predicted probability of favourable AML")
axes[0].set_ylabel("Cumulative fraction")
axes[0].set_xlim(0,1)


# Right panel: Histogram density


sns.histplot(
    data=plot_df,
    x="Probability",
    hue="Dataset",
    hue_order=["TARGET","TCGA","OHSU 2018"],
    palette=palette,
    bins=30,
    stat="density",
    common_norm=False,
    element="step",
    fill=True,
    alpha=0.25,
    ax=axes[1]
)

axes[1].axvline(
    threshold,
    linestyle="--",
    color="black",
    linewidth=1.5
)

axes[1].set_title("Probability density distributions")
axes[1].set_xlabel("Predicted probability of favourable AML")
axes[1].set_ylabel("Density")
axes[1].set_xlim(0,1)

plt.tight_layout()
plt.savefig(
    "final_probability_distribution_ecdf_hist.png",
    dpi=300
)
plt.show()

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
plt.title("TCGA favourable fusion Classification Confusion Matrix")

plt.tight_layout()
plt.savefig("tcga_confusion_matrix.png", dpi=300)
plt.show()

# VALIDATA confusion matrix
cm_validata = confusion_matrix(y_true_validata, y_pred_validata)

plt.figure(figsize=(5,4))
sns.heatmap(
    cm_validata,
    annot=True,
    fmt="d",
    cmap="Blues",
    xticklabels=["Non-fav","Fav"],
    yticklabels=["Non-fav","Fav"]
)

plt.xlabel("Predicted")
plt.ylabel("True")
plt.title("2018 OHSU Classification Confusion Matrix")

plt.tight_layout()
plt.savefig("validata_confusion_matrix.png", dpi=300)
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, average_precision_score

# load predictions
results = pd.read_csv(os.path.join(parent_dir, "results", "tcga_favourable_fusion_predictions.csv"))

# load clinical
clinical = pd.read_csv(
    os.path.join(parent_dir, "laml_tcga_pub_clinical_data.tsv"),
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
y_prob = results["favourable_fusion_probability"]

precision, recall, thresholds = precision_recall_curve(y_true, y_prob)
ap = average_precision_score(y_true, y_prob)

plt.figure(figsize=(6,5))

plt.plot(thresholds, precision[:-1], linewidth=2, label="Precision")
plt.plot(thresholds, recall[:-1], linewidth=2, label="Recall")

plt.axvline(0.675, linestyle="--", color="black")

plt.xlabel("Classification threshold (predicted probability of favourable AML)", fontsize=12)
plt.ylabel("Score", fontsize=12)

plt.title("Precision–Recall Trade-off", fontsize=13)

plt.legend(frameon=False)

# Annotation for chosen operating point
ymax = plt.gca().get_ylim()[1]
plt.annotate(
    "Chosen operating point",
    xy=(0.675, ymax * 0.8),
    xytext=(0.675 + 0.07, ymax * 0.9),
    arrowprops=dict(arrowstyle="->", lw=1.2),
    fontsize=10
)

# Note on sharp drop
plt.text(0.85, 0.05, "Sharp drop reflects few\nhigh-confidence samples", 
         fontsize=8, ha='center', va='bottom', transform=plt.gca().transAxes)

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

plt.title("Favourable Fusion Probability Distribution")
plt.xlabel("Predicted favourable fusion probability")

plt.tight_layout()
plt.savefig("tcga_probability_distribution.png", dpi=300)
plt.show()
plt.figure(figsize=(6,5))

plt.plot(thresholds, precision[:-1], label="Precision")
plt.plot(thresholds, recall[:-1], label="Recall")

plt.axvline(0.675, linestyle="--")

plt.xlabel("Probability Threshold")
plt.ylabel("Precision / Recall (proportion)")
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