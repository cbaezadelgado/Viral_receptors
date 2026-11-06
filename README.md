# Viral receptors

Scripts, pipelines, databases and models from the paper

**Cellular receptors for mammalian viruses**

_Ana Valero-Rello, Carlos Baeza-Delgado, Iván Andreu-Moreno, and Rafael Sanjuán_

_PLoS Pathog. 2024 Feb 20;20(2):e1012021. doi: [10.1371/journal.ppat.1012021](https://doi.org/10.1371/journal.ppat.1012021)_

---

## What is in this repository

The repository has two independent parts:

| Component | Language | What it does |
|---|---|---|
| [`automated_text_mining/`](automated_text_mining/) | Python + Bash | Pipeline (**virus_relations**) that reads PubMed abstracts and pulls out sentences where a **virus** is related to a **gene/protein** (candidate receptor), a drug or another entity type. It uses BERN2, NCBI Entrez, spaCy and Stanford OpenIE. |
| [`gradient_boosting_model/`](gradient_boosting_model/) | R | **XGBoost** classifier that predicts which human genes are likely viral receptors. It includes Bayesian hyper-parameter tuning, replicate-based model selection and feature importance. |

```
Viral_receptors/
├── README.md
├── automated_text_mining/            # virus–gene relation extraction from abstracts
│   ├── README.md                     # full documentation of the pipeline
│   ├── config.yaml                   # run configuration
│   ├── run_virus_relations.sh        # entry point
│   ├── extract_virus_relations.py    # main orchestration script
│   ├── requirements.txt
│   ├── scripts/                      # download, preprocessing, BERN2 control, mention parsing
│   ├── spacy/                        # spaCy relation extraction
│   ├── openie/                       # Stanford OpenIE relation extraction
│   ├── databases/alias_virus.json    # virus name → aliases dictionary
│   └── other_resources/virus_6034.txt
└── gradient_boosting_model/
    ├── XGBoost pipeline.R            # tuning, selection, training and feature importance
    └── Viral_Receptor_Prediction_Dataset.RData
```

---

## 1. Automated text mining (`virus_relations`)

You give the pipeline a list of PubMed IDs. It returns an Excel table of extracted relations such as

| term_1 | relation | term_2 |
|---|---|---|
| Angiotensin-converting enzyme 2 | is cellular receptor for | SARS coronavirus |
| Human signaling lymphocyte activation molecule | acts as | cellular receptor for Measles virus |

Each mention is normalised to NCBI Taxonomy (viruses) or NCBI Gene (genes). The pipeline runs in four steps:

1. Download the abstracts, expand acronyms and split the text into sentences.
2. Run named-entity recognition and normalisation with BERN2.
3. Curate the virus and gene mentions and keep only sentences where the requested entity types co-occur.
4. Extract relations with spaCy and OpenIE.

➡️ **Installation, configuration, usage, output format and preliminary validation are in [`automated_text_mining/README.md`](automated_text_mining/README.md).**

---

## 2. Gradient boosting receptor prediction model

`gradient_boosting_model/XGBoost pipeline.R` trains a binary XGBoost classifier (`binary:logistic`). The classifier separates known viral receptors (label `1`) from all other genes (label `0`).

### Input data

`Viral_Receptor_Prediction_Dataset.RData` loads two objects:

- **`Data`**: sparse matrix with one row per gene. Column 1 is the target label (known receptor = 1). The remaining columns are the predictive features. Missing values are `NA` and XGBoost handles them natively.
- **`Gene_info`**: gene names with their Ensembl and UniProt identifiers, row-aligned with `Data`.

### What the script does

1. **Bayesian hyper-parameter optimisation** with `ParBayesianOptimization::bayesOpt`, run in parallel on all cores.
   - Tuned parameters and their bounds: `eta` [0.01–0.5], `gamma` [0–10], `max_depth` [1–10], `subsample` [0.1–1], `colsample_bytree` [0.1–1], `alpha` [0–20], `lambda` [0–20].
   - Search size: 100 initial points plus 500 iterations, run 10 per epoch.
   - Each candidate is scored with stratified 10-fold `xgb.cv`, using early stopping (50 rounds) and up to 10 000 rounds.
   - Class imbalance is corrected with `scale_pos_weight = #neg / #pos`.
   - The objective is **Mean_PR = (precision + recall) / 2**, computed at a 0.5 threshold.
2. **Robust model selection**.
   - The 50 best parameter sets are each re-run in **100 cross-validation replicates** (5 000 CV runs).
   - For every set the script records TN/FP/FN/TP, precision, recall, F1, Mean_PR and AUC.
   - The best set is the one whose Mean_PR is significantly greater (one-sided t-test, p < 0.01) than the highest threshold.
   - `plot_meanPR` shows these per-model significance thresholds.
3. **Averaged predictions**. The 100 replicate CV predictions of the best model are averaged per gene into `Avg_Pred`, which contains `Gene_info` + `Avg_model_pred` + `label`. Final metrics, including the AUC from `pROC`, are stored in `Avg_best_model_metrics`.
4. **Final model and feature importance**.
   - An `xgboost` model is trained on the whole dataset with the best parameters and the average number of rounds.
   - Variable importance is computed with `xgb.importance` and plotted with `xgb.ggplot.importance`.

### Requirements

R (≥ 4.0) with the packages:

```r
install.packages(c("ParBayesianOptimization", "xgboost", "Matrix", "tidyverse",
                   "doParallel", "pROC", "Ckmeans.1d.dp"))
```

(`parallel` ships with base R.)

### Running

1. Open `XGBoost pipeline.R` in R/RStudio.
2. Make sure the `load("./Viral_Receptor_Prediction_Dataset.RData")` path points to the dataset. Either set the working directory to `gradient_boosting_model/` or edit the path.
3. Run the script section by section. The RStudio outline shows `#----` headers.

> ⚠️ The script is compute-intensive. It runs about 600 CV runs for the optimisation and 5 000 for the replicates, and it uses `detectCores()` workers by default. To limit resource use, lower `no_cores`, `iters.n`, `N` or `replicates`.

Main objects left in the R session:

| Object | Content |
|---|---|
| `opt_bayes$scoreSummary` | All Bayesian optimisation trials and their metrics |
| `XGBMetrics` / `XGBMeanMetrics` | Per-replicate / per-model metrics |
| `Best_model_ID` | Selected parameter set |
| `Avg_Pred` | Averaged receptor probability for every gene |
| `Avg_best_model_metrics` | Final performance of the selected model |
| `XGBOOST_Model`, `importance` | Final trained model and feature importance |

---

## How to cite

If you use any code or data from this repository, please cite:

> Valero-Rello A, Baeza-Delgado C, Andreu-Moreno I, Sanjuán R. Cellular receptors for mammalian viruses. *PLoS Pathog.* 2024;20(2):e1012021. https://doi.org/10.1371/journal.ppat.1012021

## Contact

Sanjuán Lab, Institute for Integrative Systems Biology (I2SysBio), Universitat de València–CSIC.
For questions about the code, please open an issue in this repository.
