import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score, cross_val_predict, train_test_split
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score
from sklearn.preprocessing import label_binarize
import joblib

# =========================
# 1. LOAD EXPRESSION MATRIX
# =========================
# Load preprocessed/scaled expression matrix (n_samples x n_genes) from cleaning.py.
X_ohsu = np.load(r"C:\Users\mba22ew\test\cleaned\X_ohsu.npy")

# Also load the dataframe to extract sample IDs and gene names for downstream reporting.
expr = pd.read_csv(r"C:\Users\mba22ew\test\cleaned\ohsu_cleaned_expression.csv", index_col=0)

ohsu_samples = expr.index  # sample identifiers
genes = expr.columns       # gene names (columns of X_ohsu)
print("Expression matrix:", X_ohsu.shape)

# =========================
# 2. LOAD CLINICAL DATA
# =========================
# Load clinical metadata with fusion annotations
clinical = pd.read_csv(r"C:\Users\mba22ew\test\aml_ohsu_2022_clinical_data.tsv",
                       sep="\t", low_memory=False).set_index("Sample ID")

# Ensure clinical samples align with expression matrix (same sample order and subset).
clinical = clinical.loc[ohsu_samples]
print("Clinical aligned:", clinical.shape)

# =========================
# 3. MAP TO FUSION LABELS (MULTICLASS)
# =========================
def assign_fusion(label):
    """Convert free-text clinical fusion annotation to standardized multiclass label.
    
    Known fusions (favorable or otherwise) are mapped to short identifiers.
    Missing or unrecognized labels are mapped to 'NONE'.
    """
    if pd.isna(label):
        return "NONE"
    
    label = label.upper()
    
    # Map common fusion patterns; order matters for overlapping patterns.
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
        # Anything not explicitly matched is classified as no fusion.
        return "NONE"

# Apply mapping to create multiclass labels.
clinical["fusion_class"] = clinical["Cancer Type Detailed"].apply(assign_fusion)
y_ohsu_multi = clinical["fusion_class"]  # used for CV and for downstream eln_classifier.py

print("\nFusion class distribution (multiclass):")
print(y_ohsu_multi.value_counts())

# =========================
# 4. CREATE BINARY LABELS (FAVORABLE vs. NON-FAVORABLE)
# =========================
# Define ELN-favorable fusions (core binding factor and PML-RARA).
favorable_fusions = ["CBFB_MYH11", "PML_RARA", "RUNX1_RUNX1T1"]

# Convert multiclass to binary: 1 = favorable, 0 = non-favorable.
y_ohsu = y_ohsu_multi.apply(lambda x: 1 if x in favorable_fusions else 0)

print("\nBinary label distribution (1=favorable, 0=non-favorable):")
print(y_ohsu.value_counts())

# =========================
# 5. SANITY CHECKS
# =========================
# Verify data integrity before training.
assert X_ohsu.shape[0] == len(y_ohsu), "Sample mismatch between expression and labels!"
assert not y_ohsu.isna().any(), "Missing labels detected!"

# =========================
# 6. DEFINE RANDOM FOREST MODEL
# =========================
# Create RandomForest with class weighting to handle imbalanced binary labels.
model = RandomForestClassifier(
    n_estimators=1000,
    class_weight="balanced",   # upweight minority class (favorable)
    random_state=42,            # reproducibility
    n_jobs=-1                   # use all CPU cores
)

# =========================
# 7. STRATIFIED CROSS-VALIDATION (BINARY)
# =========================
# Use stratified K-fold to preserve class proportions in each fold.
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# Compute F1 scores across folds to estimate binary classifier performance.
scores = cross_val_score(model, X_ohsu, y_ohsu, cv=cv, scoring="f1")
print("\nCV F1 scores (binary):", scores)
print("Mean CV F1:", scores.mean())

# Generate cross-validated predictions for a classification report.
y_pred_cv = cross_val_predict(model, X_ohsu, y_ohsu, cv=cv)
print("\nClassification report (Cross-Validation, binary):")
print(classification_report(y_ohsu, y_pred_cv))

# =========================
# 8. HOLD-OUT EVALUATION (BINARY)
# =========================
# Split data for a final hold-out test (80/20 split).
X_train, X_test, y_train, y_test = train_test_split(
    X_ohsu, y_ohsu, test_size=0.2, stratify=y_ohsu, random_state=42
)

# Train a separate model on the training split for evaluation.
model_split = RandomForestClassifier(
    n_estimators=1000,
    class_weight="balanced",
    random_state=42,
    n_jobs=-1
)

model_split.fit(X_train, y_train)
y_pred_split = model_split.predict(X_test)
y_prob_split = model_split.predict_proba(X_test)[:, 1]

# Report hold-out classification metrics.
print("\n=== HOLD-OUT EVALUATION (binary) ===")
print(classification_report(y_test, y_pred_split))

# Show confusion matrix.
cm = confusion_matrix(y_test, y_pred_split)
cm_df = pd.DataFrame(cm, index=["Non-favorable", "Favorable"], columns=["Pred Non-fav", "Pred Fav"])
print("\nConfusion Matrix:")
print(cm_df)

# Compute binary ROC-AUC.
roc_auc = roc_auc_score(y_test, y_prob_split)
print("\nROC-AUC:", roc_auc)

# =========================
# 9. FEATURE IMPORTANCE
# =========================
# Retrain the hold-out model on full dataset to get final feature importances.
model_split.fit(X_ohsu, y_ohsu)

# Extract and save feature importances (which genes drive the classification).
importance = pd.Series(model_split.feature_importances_, index=genes).sort_values(ascending=False)
importance.to_csv("fusion_feature_importance_binary.csv")

print("\nTop 20 genes driving binary classification:")
print(importance.head(20))

# =========================
# 10. SAVE MODEL AND LABELS
# =========================
# Save the trained binary model.
joblib.dump(model_split, "fusion_rf_binary_ohsu.joblib")

# Save binary labels (for reference/reuse).
np.save("ohsu_favorable_labels.npy", y_ohsu.values)
y_ohsu.to_csv(r"C:\Users\mba22ew\test\cleaned\ohsu_favorable_labels.csv")

# Save multiclass labels (IMPORTANT: needed by eln_classifier.py downstream).
# Construct a dataframe with multiclass labels and sample IDs as index.
multiclass_labels_df = pd.DataFrame({"fusion_class": y_ohsu_multi}, index=ohsu_samples)
multiclass_labels_df.to_csv(r"C:\Users\mba22ew\test\cleaned\ohsu_fusion_labels.csv")

print("\n✅ Binary classification training complete")
print("Samples:", X_ohsu.shape[0])
print("Genes:", X_ohsu.shape[1])
print("Favorable samples:", y_ohsu.sum())
print("Non-favorable samples:", len(y_ohsu) - y_ohsu.sum())
print("\nOutputs saved:")
print("  - fusion_rf_binary_ohsu.joblib")
print("  - fusion_feature_importance_binary.csv")
print("  - cleaned/ohsu_favorable_labels.csv (binary)")
print("  - cleaned/ohsu_fusion_labels.csv (multiclass, for eln_classifier.py)")