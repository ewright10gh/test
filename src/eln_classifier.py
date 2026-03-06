import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, roc_auc_score
from joblib import dump

# =========================
# 1. LOAD EXPRESSION MATRIX
# =========================
expr = pd.read_csv(r"C:\Users\mba22ew\test\cleaned\ohsu_cleaned_expression.csv", index_col=0)
X = expr.values
genes = expr.columns
samples = expr.index

print("Expression matrix:", X.shape)

# =========================
# 2. LOAD CLINICAL DATA
# =========================
clinical = pd.read_csv(
    r"C:\Users\mba22ew\test\aml_ohsu_2022_clinical_data.tsv",
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
    """Map free-text clinical fusion to concise label."""
    if pd.isna(label):
        return "NONE"
    label = label.upper()
    if "PML-RARA" in label:
        return "PML_RARA"
    elif "RUNX1-RUNX1T1" in label:
        return "RUNX1_RUNX1T1"
    elif "CBFB-MYH11" in label:
        return "CBFB_MYH11"
    elif "MLLT3-MLL" in label or "KMT2A" in label:
        return "KMT2A"
    elif "RPN1-EVI1" in label or "INV(3)" in label:
        return "MECOM_INV3"
    elif "DEK-NUP214" in label:
        return "DEK_NUP214"
    else:
        return "NONE"

clinical["fusion_class"] = clinical["Cancer Type Detailed"].apply(assign_fusion)

# =========================
# 4. CREATE BINARY ELN LABELS
# =========================
favorable_fusions = ["PML_RARA", "RUNX1_RUNX1T1", "CBFB_MYH11"]
y = clinical["fusion_class"].isin(favorable_fusions).astype(int)

print("\nBinary label distribution (1 = favorable):")
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
coefs.to_csv(r"C:\Users\mba22ew\test\cleaned\ridge_feature_importance.csv")

print("\nTop 20 genes driving favorable fusion classification:")
print(coefs.head(20))

# =========================
# 10. SAVE MODEL + SCALER + LABELS
# =========================
os.makedirs(r"C:\Users\mba22ew\test\models", exist_ok=True)
dump(ridge, r"C:\Users\mba22ew\test\models\ridge_eln_model.joblib")
dump(scaler, r"C:\Users\mba22ew\test\models\ridge_scaler.joblib")
np.save(r"C:\Users\mba22ew\test\models\ohsu_favorable_labels.npy", y.values)
y.to_csv(r"C:\Users\mba22ew\test\cleaned\ohsu_favorable_labels.csv")

print("\n✅ ELN Ridge model pipeline complete")
print("Samples:", X.shape[0])
print("Genes:", X.shape[1])
print("Favorable samples:", y.sum())
print("Non-favorable samples:", len(y) - y.sum())