# Product Requirements Document
## Histone Mark Co-occurrence Network Analysis Pipeline

**Project:** histone_mark_mining  
**Author:** Dr. Arli Aditya Parikesit  
**Version:** 1.2  
**Date:** 2026-06-16  
**Status:** Active development

---

## 1. Problem Statement

Histone modifications are key epigenetic regulators that govern chromatin accessibility, gene expression, and cell identity. Understanding how multiple histone marks co-occur across genomic regions, and whether their interaction topology exhibits scale-free properties, has direct implications for identifying regulatory hubs in cancer, development, and disease.

No lightweight, reproducible Python pipeline currently exists that integrates histone mark co-occurrence network construction, scale-free topology testing, statistical power-law fitting, and chromatin state clustering in a single script suitable for both exploratory analysis and HPC batch execution.

---

## 2. Goals

- Build a reproducible, end-to-end pipeline for histone mark co-occurrence network analysis.
- Test whether the resulting interaction network exhibits scale-free topology using statistical (not visual-only) methods.
- Cluster genomic regions into chromatin states and evaluate cluster quality with scientifically valid intrinsic metrics.
- Produce publication-ready output figures and quantitative summary statistics.
- Support both live database ingestion and simulated fallback data.

## 3. Non-Goals

- Real-time or streaming epigenomic data processing.
- Full ChIP-seq peak calling or read alignment (upstream of this pipeline).
- ATAC-seq, bisulfite-seq, or RNA-seq integration (separate pipelines).
- Supervised classification of chromatin states without an independent biological label source (flagged for future work).
- Production-grade deployment (authentication, multi-user concurrency, persistent storage).

---

## 4. Users

| User | Context |
|------|---------|
| Bioinformatician (primary) | Runs pipeline on processed ChIP-seq count tables; interprets network topology |
| Wet-lab collaborator | Reads output figures and summary statistics to prioritize experimental follow-up |
| HPC operator | Submits batch jobs; needs non-interactive execution and parameterized output paths |

---

## 5. Functional Requirements

### 5.1 Data Ingestion

| ID | Requirement |
|----|-------------|
| F1 | Pipeline attempts to fetch live data from Histome2 (ACTREC); falls back to reproducible simulated dataset on network failure or unexpected HTML structure. |
| F2 | Simulated dataset: 100 genomic regions, 3 histone marks (mark1--mark3), Poisson-distributed read counts; seeded for reproducibility (`np.random.seed(42)`). |
| F3 | If live data arrives without a `count` column, pipeline synthesizes one from Poisson(lam=25) rather than crashing. |
| F4 | All numeric histone mark columns (prefix `mark`) are auto-detected; pipeline is not hard-coded to exactly 3 marks. |

### 5.2 Normalization

| ID | Requirement |
|----|-------------|
| F5 | Histone mark columns are standardized (zero mean, unit variance) using `StandardScaler` before any downstream analysis. |

### 5.3 Graph Construction

| ID | Requirement |
|----|-------------|
| F6 | Directed graph (`nx.DiGraph`) with one node per genomic region. |
| F7 | Edges added between regions whose pairwise cosine similarity of standardized histone mark vectors exceeds a configurable threshold (default 0.95). |
| F8 | Pipeline reports node count, edge count, and top-5 regions by degree, betweenness, and closeness centrality. |

### 5.4 Correlation Analysis

| ID | Requirement |
|----|-------------|
| F9 | Pearson correlation matrix computed and printed for all histone mark pairs. |

### 5.5 Scale-Free Network Testing

| ID | Requirement |
|----|-------------|
| F10 | Average clustering coefficient of the constructed network compared against a size-matched Erdos-Renyi random graph. |
| F11 | Degree distribution plotted on log-log axes and saved to disk. |
| F12 | KS test performed comparing observed degree distribution to a Poisson distribution of matching mean. |
| F13 | Power-law exponent gamma estimated by log-log OLS (`scipy.optimize.curve_fit`); report gamma, standard error, and R² on log-log scale. Flag if gamma falls outside the scale-free range [2, 3]. Note: MLE via `powerlaw` package (Clauset et al. 2009) required for publication-grade claims. |

### 5.6 Clustering and Dimensionality Reduction

| ID | Requirement |
|----|-------------|
| F14 | PCA applied to histone mark features; components capped at `min(2, n_samples, n_features)`. |
| F15 | KMeans clustering (k=3, or fewer if n_samples < 3) applied to PCA-reduced data to define chromatin states. |

### 5.7 Cluster Quality Evaluation

| ID | Requirement |
|----|-------------|
| F16 | Silhouette score and Davies-Bouldin index computed on original histone mark features against KMeans labels. Neither metric uses KMeans labels as pseudo ground truth for a supervised model. |
| F17 | Pipeline prints interpretation guidance alongside scores (range, direction of improvement). |

### 5.8 Statistical Testing

| ID | Requirement |
|----|-------------|
| F18 | Per-cluster Poisson significance test using the `count` column mean as lambda; survival function p-value reported for each cluster. |

### 5.10 Streamlit Dashboard

| ID | Requirement |
|----|-------------|
| F22 | Dashboard launched via `streamlit run app.py`; opens on `localhost:8501`. |
| F23 | Sidebar controls: cosine similarity threshold (slider 0.50-1.00), number of clusters k (slider 2-5), Histome2 live data toggle (checkbox). All controls trigger live recomputation via `st.cache_data`. |
| F24 | Six tabs matching pipeline steps: Data Overview, Network Graph, Correlation, Scale-Free Analysis, Clustering, Statistics. Each tab displays figures (via `st.pyplot`), metrics (`st.metric`), and tables (`st.dataframe`). |
| F25 | Power-law OLS fit overlaid on log-log degree distribution plot; gamma value and [2,3] range check displayed as colored status indicator. |
| F26 | All matplotlib figures closed after rendering (`plt.close`) to prevent memory accumulation across Streamlit reruns. |

### 5.11 Output

| ID | Requirement |
|----|-------------|
| F19 | Two figures saved: `scale_free_network.png` (network graph, hubs highlighted in red) and `degree_distribution.png` (log-log degree distribution). DPI=300. |
| F20 | Output directory configurable via `--output-dir` CLI argument (default: current working directory). Directory created if absent. |
| F21 | All `plt.show()` calls replaced with `plt.close()`; execution never blocks in non-interactive environments. |

---

## 6. Non-Functional Requirements

| ID | Requirement |
|----|-------------|
| NF1 | Reproducible: fixed random seeds throughout; identical outputs given same data. |
| NF2 | Non-interactive: no blocking calls; safe to run under Slurm/SGE without a display. |
| NF3 | Self-contained: single Python script plus test file; all dependencies pip-installable. |
| NF4 | Testable: 9 smoke tests covering data generation, network fallback, graph centrality correctness; all pass before merge. |
| NF5 | Documented: script header with module, purpose, author, date, parameters, and references. |
| NF6 | Interactive: Streamlit dashboard exposes all six pipeline steps with live parameter controls; no page reload required on parameter change. |

---

## 7. Data Sources

| Source | Type | Status |
|--------|------|--------|
| Histome2 (ACTREC) | HTML scrape | Fragile; HTML structure unverified against live site. Functional only if `div.region` / `div.mark` classes exist. |
| Simulated data | `np.random` | Active fallback; seed=42 guarantees reproducibility. |
| ENCODE / GEO | ChIP-seq count tables | Intended real-data source for production use; not yet integrated. |
| EpiFactors (autosome.org) | REST API | Removed (no usable API endpoint identified); candidate for future reintegration via official download link. |

---

## 8. Technical Architecture

```
histone_mark_mining.py
│
├── data_mining()               # Simulated fallback dataset
├── fetch_histome2_data()       # Live Histome2 scraper (with fallback)
│
└── main(output_dir)
    ├── Step 1: Ingest + normalize
    ├── Step 2: Graph construction + centrality
    ├── Step 3: Pearson correlation
    ├── Step 4: Scale-free comparison + visualization
    ├── Step 5: Power-law fit (OLS) + KS test
    ├── Step 6: PCA + KMeans clustering
    ├── Step 7: Cluster quality (silhouette, Davies-Bouldin)
    └── Step 8: Poisson cluster significance

tests/
└── test_pipeline.py            # 9 smoke tests (pytest)
```

### 8.1 Knowledge Graph Analysis

Graph built via code-review-graph (2026-06-16): **16 nodes, 265 edges, 2 files, 1 flow, 2 communities.**

**Communities:**

| Community | Nodes | Cohesion | Contents |
|-----------|-------|----------|----------|
| `histone-mark-data` | 4 | 0.018 | `data_mining`, `fetch_histome2_data`, `_power_law_log`, `main` |
| `tests-data` | 10 | 0.000 | 9 test functions + `built_graph` fixture |

**Execution flow:**

| Flow | Criticality | Nodes | Notes |
|------|-------------|-------|-------|
| `main` | 0.36 | 4 | Only entry point; no sub-flows or CLI subcommands |

**Architecture diagram:**

```mermaid
graph TD
    subgraph HM["histone-mark-data (pipeline)"]
        DM["data_mining()"]
        FH["fetch_histome2_data()"]
        PL["_power_law_log()"]
        M["main(output_dir)"]
    end

    subgraph TD["tests-data (test suite)"]
        BG["built_graph [fixture]"]
        T1["test_data_mining_*  (x5)"]
        T2["test_fetch_histome2_*"]
        T3["test_graph_* / test_centrality_* (x3)"]
    end

    DM -->|simulated fallback| M
    FH -->|live data or None| M
    M -->|inline call| PL

    TD -. "sys.path import (runtime only; not graph-visible)" .-> HM
```

**Key structural findings:**

- `main()` concentrates all 8 pipeline steps inline. Cohesion of `histone-mark-data` is 0.018 (near-zero) because only `main` calls the other three functions; no inter-function calls exist outside `main`. Each step is a candidate for extraction into a discrete function to improve testability and graph cohesion.
- `tests-data` cohesion is 0 by design — flat pytest suite with no shared helper calls (the `built_graph` fixture is the only shared construct).
- Cross-community coupling is 0 in the static graph. The `sys.path.insert + import histone_mark_mining` pattern in `test_pipeline.py` is resolved at runtime but invisible to the static parser. Adding a `conftest.py` or proper packaging (`pyproject.toml`) would not resolve this; it is a fundamental limitation of static import tracing for path-manipulated imports.

**Key dependencies:**

| Package | Version constraint | Role |
|---------|--------------------|------|
| networkx | >=2.6 | Graph construction and centrality |
| scikit-learn | >=1.0 | StandardScaler, PCA, KMeans, silhouette, Davies-Bouldin, cosine similarity |
| scipy | >=1.7 | KS test, Poisson sf, curve_fit |
| matplotlib | >=3.4 | Figure generation |
| pandas | >=1.3 | DataFrame operations |
| numpy | >=1.21 | Numerical operations |
| requests | >=2.26 | Histome2 HTTP fetch |
| beautifulsoup4 | >=4.10 | Histome2 HTML parsing |

---

## 9. Deliverables

| Deliverable | Description |
|-------------|-------------|
| `histone_mark_mining.py` | Main analysis pipeline (CLI) |
| `app.py` | Streamlit interactive dashboard |
| `tests/test_pipeline.py` | Smoke test suite (pytest) |
| `scale_free_network.png` | Directed network graph; hubs in red |
| `degree_distribution.png` | Log-log degree distribution with OLS power-law fit |
| Console output | Centrality top-5, correlation matrix, KS test, gamma fit, silhouette/DB scores, Poisson p-values |

---

## 10. Validation Strategy

| Stage | Method | Pass Criterion |
|-------|--------|----------------|
| Unit | `pytest tests/test_pipeline.py` | 9/9 pass |
| Reproducibility | Run pipeline twice with seed=42 | Identical figures and statistics |
| Scale-free claim | Gamma in [2, 3] AND KS p < 0.05 vs. Poisson | Both conditions must hold before asserting scale-free topology in manuscript |
| Cluster validity | Silhouette > 0.5 AND Davies-Bouldin < 1.0 | Acceptable cluster separation for 3-state chromatin model |
| Live data path | Manual test against Histome2 URL | Returns non-empty DataFrame; mark columns auto-detected |

---

## 11. Known Limitations

1. **Simulated biology.** `data_mining()` hard-codes `mark2 = 0.6*mark1 + noise`. This artificial correlation inflates co-occurrence structure. All conclusions drawn from simulated data are illustrative only.
2. **HTML scraper fragility.** Histome2 scraper relies on undocumented CSS class names. A site redesign will silently trigger fallback.
3. **OLS power-law fit.** Log-log OLS underestimates gamma compared to MLE and cannot estimate x_min. Use `powerlaw` package for manuscript results.
4. **Unsupervised cluster labels.** Chromatin state labels come from KMeans, not from biological ground truth (e.g., ChIP-seq peak calls, ENCODE annotations). Supervised prediction requires an independent label source.
5. **O(n^2) similarity loop.** Graph edge construction iterates all region pairs. Scales poorly beyond ~5,000 regions; vectorized threshold masking needed for large datasets.
6. **EpiFactors not integrated.** No verified API endpoint for programmatic access; data must be manually downloaded for multi-source analysis.
7. **`main()` monolith.** All 8 pipeline steps are inline inside `main()`. Graph cohesion of the pipeline community is 0.018. Individual steps cannot be unit-tested, profiled, or swapped independently without refactoring. See roadmap.
8. **Static import tracing gap.** `sys.path.insert + import histone_mark_mining` in the test file is not resolved by static analysis tools. Cross-community coupling reports 0 edges despite runtime dependency. Proper packaging (`pyproject.toml` with `pip install -e .`) is the only full fix.

---

## 12. Open Questions

| # | Question | Owner |
|---|----------|-------|
| Q1 | Which real ChIP-seq dataset (ENCODE project, GEO accession) will replace simulated data? | Dr. Parikesit |
| Q2 | Should the similarity threshold (0.95) be tuned per dataset, or is a fixed value acceptable for the target publication? | Dr. Parikesit |
| Q3 | Is EpiFactors integration required for submission, or is Histome2 alone sufficient? | Dr. Parikesit |
| Q4 | What biological label source will be used for future supervised chromatin state classification? | Dr. Parikesit |
| Q5 | Target journal and its data/code availability requirements (e.g., GitHub release, GEO deposition)? | Dr. Parikesit |

---

## 13. Roadmap

| Priority | Item | Rationale |
|----------|------|-----------|
| HIGH | Replace simulated data with real ChIP-seq count table (ENCODE or GEO) | All biological conclusions require real data |
| HIGH | Verify or rewrite Histome2 scraper against actual live HTML structure | Current scraper untested against live site |
| HIGH | Replace OLS gamma estimate with MLE via `powerlaw` package for publication | OLS is insufficient for power-law claim in peer review |
| MEDIUM | Add independent biological label source for supervised chromatin state model | Resolves tautological evaluation from pre-merge review |
| MEDIUM | Vectorize edge construction with `np.where` masking | Required for datasets > 1,000 regions |
| MEDIUM | Add EpiFactors integration via CSV bulk download | Multi-source validation strengthens epigenetic factor coverage |
| LOW | Add `environments/environment.yml` and `environments/requirements.txt` | Reproducibility for collaborators |
| LOW | Add `--threshold` CLI argument for similarity cutoff | Sensitivity analysis without code edits |
| LOW | Refactor `main()` into discrete step functions (e.g., `build_graph`, `run_pca`, `evaluate_clusters`) | Raises graph cohesion from 0.018; enables per-step unit testing and profiling |
| LOW | Add `pyproject.toml` with `pip install -e .` entry point | Enables normal `import histone_mark_mining` without `sys.path` manipulation; makes test-to-module coupling graph-visible |
| DONE | Streamlit interactive dashboard (`app.py`) | Six-tab dashboard with sidebar parameter controls and live recomputation |
