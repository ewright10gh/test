import numpy as np
import pandas as pd

X_ohsu = np.load(r"C:\Users\mba22ew\test\cleaned\X_ohsu.npy")

ohsu_samples = pd.read_csv(
    r"C:\Users\mba22ew\test\cleaned\ohsu_cleaned_expression.csv",
    index_col=0
).index

clinical = pd.read_csv("aml_ohsu_2022_clinical_data.tsv", sep="\t")
for v in clinical["Cancer Type Detailed"].unique():
    print(repr(v))


fusion_keywords = [
    "PML-RARA",
    "RUNX1-RUNX1T1",
    "CBFB-MYH11",
]

fusion_pattern = "|".join(fusion_keywords)

clinical["fusion_status"] = (
    clinical["Cancer Type Detailed"]
    .str.contains(fusion_pattern, case=False, na=False, regex=True)
    .astype(int)
)

print(clinical["fusion_status"].value_counts())

clinical = clinical.set_index("Sample ID")

y_ohsu = clinical.loc[ohsu_samples, "fusion_status"]

print(y_ohsu.value_counts())
print(y_ohsu.shape)
print(X_ohsu.shape)

# from sklearn.ensemble import RandomForestClassifier

# model = RandomForestClassifier(
#     n_estimators=500,
#     class_weight="balanced",
#     random_state=42,
#     n_jobs=-1
# )

# model.fit(X_ohsu, y_ohsu)

# from sklearn.model_selection import cross_val_score

# scores = cross_val_score(
#     model,
#     X_ohsu,
#     y_ohsu,
#     cv=5,
#     scoring="roc_auc"
# )

# print("CV AUC:", scores)
# print("Mean AUC:", scores.mean())


# from sklearn.model_selection import cross_val_predict

# y_prob = cross_val_predict(
#     model,
#     X_ohsu,
#     y_ohsu,
#     cv=5,
#     method="predict_proba"
# )[:, 1]

# genes = pd.read_csv("ohsu_cleaned_expression.csv", nrows=1).columns

# importance = pd.Series(
#     model.feature_importances_,
#     index=genes
# ).sort_values(ascending=False)

# print(importance.head(20))

# importance.to_csv("feature_importance.csv")

# np.save("ohsu_cv_probabilities.npy", y_prob)