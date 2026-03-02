import numpy as np
import pandas as pd
import joblib

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    roc_auc_score
)
from sklearn.preprocessing import label_binarize


# LOAD DATA


X = np.load(r"C:\Users\mba22ew\test\cleaned\X_ohsu.npy")

y = pd.read_csv(
    r"C:\Users\mba22ew\test\cleaned\ohsu_fusion_labels.csv",
    index_col=0
)["fusion_class"]

print("Data loaded:")
print("X shape:", X.shape)
print("y shape:", y.shape)
print("\nClass distribution:\n", y.value_counts())


# TRAIN / TEST SPLIT


X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.2,
    stratify=y,
    random_state=42
)

print("\nTrain shape:", X_train.shape)
print("Test shape:", X_test.shape)


# MODEL


model = RandomForestClassifier(
    n_estimators=1000,
    class_weight="balanced",
    random_state=42,
    n_jobs=-1
)


# TRAIN


model.fit(X_train, y_train)


# PREDICT


y_pred = model.predict(X_test)
y_prob = model.predict_proba(X_test)


# METRICS


print("\n=== CLASSIFICATION REPORT ===")
print(classification_report(y_test, y_pred))

print("\n=== CONFUSION MATRIX ===")
cm = confusion_matrix(y_test, y_pred, labels=model.classes_)
cm_df = pd.DataFrame(cm, index=model.classes_, columns=model.classes_)
print(cm_df)


# MULTICLASS ROC-AUC


y_test_bin = label_binarize(y_test, classes=model.classes_)

roc_auc = roc_auc_score(
    y_test_bin,
    y_prob,
    multi_class="ovr",
    average="macro"
)

print("\nMacro ROC-AUC:", roc_auc)


# FEATURE IMPORTANCE


genes = pd.read_csv(
    r"C:\Users\mba22ew\test\cleaned\ohsu_cleaned_expression.csv",
    nrows=1,
    index_col=0
).columns

importance = pd.Series(
    model.feature_importances_,
    index=genes
).sort_values(ascending=False)

importance.head(30).to_csv("top_genes_final_model.csv")

print("\nTop 10 genes:")
print(importance.head(10))


# SAVE MODEL


joblib.dump(model, "fusion_rf_model.pkl")

print("\n✅ FINAL MODEL TRAINED & SAVED")