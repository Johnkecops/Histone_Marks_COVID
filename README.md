# Histone Mark Co-occurrence Network Analysis Pipeline

A reproducible Python pipeline for histone mark co-occurrence network construction, scale-free topology testing, statistical power-law fitting, and chromatin state clustering. Designed for exploratory analysis and HPC batch execution.

**Author:** Dr. Arli Aditya Parikesit  
**Affiliation:** Department of Biotechnology, i3L University, Jakarta  
**ORCID:** https://orcid.org/0000-0001-8716-3926

---

## Background

Histone modifications regulate chromatin accessibility, gene expression, and cell identity. This pipeline asks whether histone mark co-occurrence across genomic regions forms a scale-free interaction network — a topology characterized by power-law degree distribution and hub regions with disproportionate connectivity. Identifying such hubs has direct implications for cancer epigenetics and transcriptional regulatory network analysis.

---

## Requirements

Python 3.8+. Install dependencies:

```bash
pip install networkx matplotlib pandas numpy requests beautifulsoup4 \
            scikit-learn scipy streamlit defusedxml
```

| Package | Version | Role |
|---------|---------|------|
| networkx | >=2.6 | Graph construction and centrality |
| scikit-learn | >=1.0 | Normalization, PCA, KMeans, cluster quality, cosine similarity |
| scipy | >=1.7 | KS test, Poisson sf, power-law curve fitting |
| matplotlib | >=3.4 | Figure output (non-interactive and Streamlit) |
| pandas | >=1.3 | DataFrame operations |
| numpy | >=1.21 | Numerical operations |
| requests | >=2.26 | HTTP fetch for all three live databases |
| beautifulsoup4 | >=4.10 | Histome2 HTML parsing |
| defusedxml | >=0.7 | Safe XML parsing of ProHistoneDB phylotree files (XXE protection) |
| streamlit | >=1.30 | Interactive web dashboard |

---

## Usage

### Command-line pipeline

```bash
# Run with default output directory (current working directory)
python histone_mark_mining.py

# Specify output directory (created if absent)
python histone_mark_mining.py --output-dir results/figures
```

Safe for non-interactive/HPC execution (no `plt.show()` blocking).

### Streamlit dashboard

```bash
streamlit run app.py
```

Opens at `http://localhost:8501`. Use the sidebar to adjust the similarity threshold, number of clusters, and data source. All six pipeline steps are displayed as interactive tabs with live parameter updates.

---

## Pipeline Steps

| Step | Description |
|------|-------------|
| 1 | Data ingestion: fetch Histome2 live data or fall back to reproducible simulated dataset (seed=42) |
| 2 | Graph construction: directed cosine-similarity graph of genomic regions; degree, betweenness, and closeness centrality |
| 3 | Pearson correlation matrix across histone mark pairs |
| 4 | Scale-free comparison: clustering coefficient vs. size-matched Erdos-Renyi random graph; network visualization |
| 5 | Power-law fit: log-log degree distribution plot, KS test vs. Poisson, OLS gamma estimate with stderr and R² |
| 6 | PCA + KMeans clustering (k=3) to define chromatin states |
| 7 | Cluster quality: silhouette score and Davies-Bouldin index (no tautological label reuse) |
| 8 | Poisson significance test per cluster using read count column |

---

## Output

| File | Description |
|------|-------------|
| `scale_free_network.png` | Directed network graph; hub nodes (degree > 80th percentile) highlighted in red |
| `degree_distribution.png` | Log-log degree distribution with OLS power-law fit |
| Console | Centrality top-5, correlation matrix, KS statistic, gamma ± stderr, R², silhouette/DB scores, per-cluster Poisson p-values |

---

## Data Sources

The pipeline attempts live data ingestion in this order:

1. **Histome2** (ACTREC, India) — HTML scrape of `https://www.actrec.gov.in/histome2/`. Falls back to simulated data if the site is unreachable or the HTML structure does not match expected CSS classes.
2. **Simulated data** — 100 genomic regions, 3 histone marks, Poisson read counts; `np.random.seed(42)` guarantees reproducibility.

For production use, replace `data_mining()` with a real ChIP-seq count table (ENCODE or GEO accession). See [Known Limitations](#known-limitations).

---

## Testing

```bash
pytest tests/test_pipeline.py -v
```

9 smoke tests covering: `data_mining()` output shape and value ranges, seed reproducibility, `fetch_histome2_data()` fallback safety, graph node count, centrality dict coverage and bounds.

---

## Known Limitations

1. **Simulated biology.** The default dataset hard-codes `mark2 = 0.6*mark1 + noise`. All scale-free and clustering results on simulated data are illustrative only.
2. **HTML scraper fragility.** Histome2 scraper relies on undocumented CSS class names; a site redesign silently triggers fallback.
3. **OLS power-law fit.** Log-log OLS underestimates gamma vs. MLE and cannot estimate x_min. Use the `powerlaw` package (Clauset et al. 2009) for publication-grade claims.
4. **Unsupervised cluster labels.** KMeans labels are not biological ground truth. Supervised chromatin state prediction requires an independent label source (e.g., ChIP-seq peak calls, ENCODE chromatin state annotations).
5. **O(n^2) edge construction.** Scales poorly beyond ~5,000 regions. Vectorized masking needed for large datasets.

---

## Project Structure

```
histone_mark_mining.py     # Main pipeline (CLI)
app.py                     # Streamlit interactive dashboard
tests/
└── test_pipeline.py       # Smoke tests (pytest)
docs/
└── PRD.md                 # Product Requirements Document
.gitignore
README.md
```

---

## References

- Clauset A, Shalizi CR, Newman MEJ (2009) Power-law distributions in empirical data. *SIAM Rev* 51:661-703. https://doi.org/10.1137/070710111
- Hagberg AA, Schult DA, Swart PJ (2008) Exploring network structure, dynamics, and function using NetworkX. *Proc SciPy 2008*, pp 11-15.

---

## License

For academic and research use. Contact the author for other uses.
