# RNA-seq Data Cleaning & Integration

## 📁 Project Structure

```
test/
│── data_mrna_seq_rpkm.txt
│── data_mrna_seq_tpm.txt
│── requirements.txt
│
└── src/
    └── cleaning.py
```

---

## ⚙️ Environment Setup (Conda)

Open a terminal in VS Code and run:

```bash
conda create -n rna_seq_env python=3.9 pandas numpy -y
conda activate rna_seq_env
```

Verify installation:

```bash
python -c "import pandas, numpy; print('Environment ready')"
```

---

## 📦 Install Dependencies (optional if using requirements.txt)

```bash
pip install -r requirements.txt
```

---

## ▶️ Run the Cleaning Script

From the project root or inside `src`:

```bash
cd src
python cleaning.py
```

---

## 🧠 What the Script Does

* Loads OHSU and TARGET RNA-seq datasets
* Matches common genes between datasets
* Transposes matrices (genes → columns)
* Log₂ transforms TARGET expression values
* Outputs aligned data ready for downstream analysis

---

## 📊 Input Data

The script expects the following files in the project root:

* `data_mrna_seq_rpkm.txt`
* `data_mrna_seq_tpm.txt`

Both must be tab-separated with gene names in the first column.

---

## ❗ Notes

* File paths in `cleaning.py` are set for **Windows**.
* If you move the project, update the paths accordingly.
* Large files (>100MB) are not recommended for GitHub without Git LFS.

---

## 🚀 Future Improvements

* Save cleaned datasets to `/output`
* Add exploratory data analysis
* Add model training pipeline
* Make paths OS-independent using `pathlib`

---

## 👤 Author

MSc Bioinformatics Project
