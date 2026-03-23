import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, roc_auc_score
from joblib import dump

# Get the directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)

# =========================
# 1. LOAD EXPRESSION MATRIX
# =========================
expr = pd.read_csv(os.path.join(parent_dir, "cleaned", "ohsu_cleaned_expression.csv"), index_col=0)
X = expr.values
genes = expr.columns
samples = expr.index

print("Expression matrix:", X.shape)

# =========================
# 2. LOAD CLINICAL DATA
# =========================
clinical = pd.read_csv(
    os.path.join(parent_dir, "aml_ohsu_2022_clinical_data.tsv"),
    sep="\t",
    low_memory=False
).set_index("Sample ID")

# Align clinical data to expression samples
clinical = clinical.loc[samples]
print("Clinical data aligned:", clinical.shape)

# =========================
# 3. MAP TO FUSION LABELS
# =========================
def assign_fusion(label):
    if pd.isna(label):
        return "NONE"

    label = str(label).upper()

    # PML-RARA / APL
    if any(x in label for x in [
        "PML-RARA", "PML_RARA", "T(15;17)", "APL", "PROMYELOCYTIC"
    ]):
        return "PML_RARA"

    # RUNX1-RUNX1T1
    if any(x in label for x in [
        "RUNX1-RUNX1T1", "RUNX1_RUNX1T1", "AML1-ETO", "T(8;21)"
    ]):
        return "RUNX1_RUNX1T1"

    # CBFB-MYH11
    if any(x in label for x in [
        "CBFB-MYH11", "CBFB_MYH11", "INV(16)", "T(16;16)"
    ]):
        return "CBFB_MYH11"

    # KMT2A rearrangements
    if "KMT2A" in label or "MLL" in label:
        return "KMT2A"

    # MECOM / inv(3)
    if "INV(3)" in label or "EVI1" in label:
        return "MECOM_INV3"

    # DEK-NUP214
    if "DEK-NUP214" in label or "T(6;9)" in label:
        return "DEK_NUP214"

    return "NONE"

clinical["fusion_class"] = clinical["Cancer Type Detailed"].apply(assign_fusion)

# =========================
# 4. CREATE BINARY FAVOURABLE FUSION LABELS
# =========================
favorable_fusions = ["PML_RARA", "RUNX1_RUNX1T1", "CBFB_MYH11"]
y = clinical["fusion_class"].isin(favorable_fusions).astype(int)

print("\nBinary label distribution (1 = favourable fusion):")
print(y.value_counts())

# =========================
# 5. TRAIN/TEST SPLIT
# =========================
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, stratify=y, random_state=42
)

# =========================
# 6. SCALE DATA
# =========================
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)

# =========================
# 7. TRAIN RIDGE LOGISTIC REGRESSION WITH CV
# =========================
ridge = LogisticRegressionCV(
    Cs=np.logspace(-4, 4, 20),
    cv=5,
    penalty="l2",
    solver="lbfgs",
    max_iter=10000,
    class_weight="balanced",
    n_jobs=-1
)

ridge.fit(X_train, y_train)

# =========================
# 8. EVALUATION
# =========================
pred = ridge.predict(X_test)
prob = ridge.predict_proba(X_test)[:, 1]

print("\nBest C:", ridge.C_[0])
print("\n=== RIDGE CLASSIFICATION REPORT ===")
print(classification_report(y_test, pred))

roc_auc = roc_auc_score(y_test, prob)
print("\nROC-AUC:", roc_auc)

# =========================
# 9. FEATURE IMPORTANCE
# =========================
coefs = pd.Series(ridge.coef_[0], index=genes).sort_values(ascending=False)
coefs.to_csv(os.path.join(parent_dir, "cleaned", "ridge_feature_importance.csv"))

print("\nTop 20 genes driving favorable fusion classification:")
print(coefs.head(20))

# =========================
# 10. SAVE MODEL + SCALER + LABELS
# =========================
os.makedirs(os.path.join(parent_dir, "models"), exist_ok=True)
dump(ridge, os.path.join(parent_dir, "models", "ridge_eln_model.joblib"))
dump(scaler, os.path.join(parent_dir, "models", "ridge_scaler.joblib"))
np.save(os.path.join(parent_dir, "models", "ohsu_favorable_labels.npy"), y.values)
y.to_csv(os.path.join(parent_dir, "cleaned", "ohsu_favorable_labels.csv"))

print("\n✅ ELN Ridge model pipeline complete")
print("Samples:", X.shape[0])
print("Genes:", X.shape[1])
print("Favorable samples:", y.sum())
print("Non-favorable samples:", len(y) - y.sum())


import matplotlib.pyplot as plt

# 1. Get the range of C values (regularization strength) used in your CV
c_values = ridge.Cs_

# 2. Get the coefficients for each C (this requires a slightly different fit 
# or manual loop if using LogisticRegressionCV, but here is the logic:)
# For simplicity, we can plot the final coefficients sorted:
plt.figure(figsize=(10, 6))
plt.semilogx(c_values, ridge.coefs_paths_[1].mean(axis=0)) # Mean across CV folds

plt.title('Ridge Coefficient Paths')
plt.xlabel('C (Inverse Regularization Strength)')
plt.ylabel('Coefficients')
plt.axvline(ridge.C_, linestyle='--', color='k', label='Best C')
plt.legend()
plt.show()
