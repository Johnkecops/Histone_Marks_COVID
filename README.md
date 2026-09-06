# Histone Marks in Coronavirus Infection: Co-occurrence and Avoidance in the Human Nucleosome

A reproducible bioinformatics pipeline for histone mark co-occurrence and avoidance analysis, benchmarking three public histone databases (Histome2, HistoneDB 2.0, ProHistoneDB) alongside simulated nucleosome data.

**Author:** Dr. Arli Aditya Parikesit  
**Affiliation:** Department of Biotechnology, i3L University, Jakarta  
**ORCID:** https://orcid.org/0000-0001-8716-3926

---

## Background

Coronaviruses (SARS-CoV-2, MERS-CoV, SARS-CoV, HCoV-229E) manipulate host cell epigenetics through histone modifications to suppress antiviral responses. The SARS-CoV-2 ORF8 protein functions as a histone mimic of the ARKS motif in Histone H3, disrupting host post-translational modifications and inducing chromatin compaction. This pipeline analyzes histone mark co-occurrence and avoidance patterns across multiple data sources and contextualizes them within the molecular mechanisms of histone-mediated tissue damage.

---

## Requirements

Python 3.8+. Install dependencies:

```bash
pip install networkx matplotlib pandas numpy requests beautifulsoup4 \
            scikit-learn scipy streamlit defusedxml statsmodels
```

| Package | Version | Role |
|---------|---------|------|
| networkx | >=2.6 | Graph construction and centrality |
| scikit-learn | >=1.0 | Normalization, clustering |
| scipy | >=1.7 | Poisson significance test |
| matplotlib | >=3.4 | Figures (Agg backend for non-interactive) |
| pandas | >=1.3 | DataFrame operations |
| numpy | >=1.21 | Numerical operations |
| requests | >=2.26 | HTTP fetch for live databases |
| beautifulsoup4 | >=4.10 | Histome2 HTML parsing |
| defusedxml | >=0.7 | Safe XML parsing (XXE protection) |
| statsmodels | >=0.14 | FDR multiple testing correction |
| streamlit | >=1.30 | Interactive web dashboard |

---

## Usage

### Database benchmarking

```bash
python benchmark_databases.py
```

Fetches all three databases live, computes occurrence, co-occurrence, and avoidance (Poisson test + FDR correction), and saves three comparative figures to `results/figures/`.

### Figure generation

```bash
python gen_figures.py
```

Generates 6 figures (NETosis pathway, ORF8 mimicry network, immunothrombosis cycle, occurrence barchart, co-occurrence heatmap, network graph).

### Streamlit dashboard

```bash
streamlit run app-covid.py
```

Interactive dashboard at `http://localhost:8501` exploring histone-virus interaction mechanisms, database benchmarking results, and literature synthesis.

### Command-line pipeline

```bash
python histone_mark_mining.py [--output-dir results/figures]
```

Original pipeline for graph construction, scale-free testing, and chromatin state clustering (see [Legacy](#legacy) section).

---

## Data Sources

| Source | Type | Records | Description |
|--------|------|---------|-------------|
| **Histome2** | Live HTML scrape | 62 PTM sites | Human histone writer/eraser/mod-type counts |
| **HistoneDB 2.0** | Live CSV download | 30 variants | HMM scores, sequence counts, taxonomic breadth |
| **ProHistoneDB** | Live XML fetch | 24 categories | Proto-histone clade count, taxonomic codes, sequence length |
| **Simulated** | Synthetic (rules-based) | 200 nucleosomes | 5 binary marks with embedded co-occurrence/avoidance rules |

All continuous variables binarized by median split. Co-occurrence and avoidance assessed via Poisson test with Benjamini-Hochberg FDR correction (alpha = 0.05).

---

## Results Summary

| Data Source | Significant Pairs | Best Effect Size |
|-------------|-------------------|------------------|
| Histome2 | 0 | has_writer vs has_eraser (p_adj = 0.194, LFC = 0.67) |
| HistoneDB 2.0 | 0 | many_seq vs broad_tax (p_adj = 0.065, LFC = 0.90) |
| ProHistoneDB | 0 | many_seq vs broad_tax (p_adj = 0.507, LFC = 0.64) |
| Simulated | 6 (3 co-occur, 3 avoid) | H3K4me3 vs H3K27ac (p_adj < 0.001, LFC = 1.04) |

No real database showed significant co-occurrence or avoidance after FDR correction. The simulated positive control correctly recovered all embedded relationships, confirming analytical sensitivity.

---

## Output

| File/Directory | Description |
|----------------|-------------|
| `results/figures/fig1_db_occurrence_comparison.png` | Occurrence bar chart across all 4 data sources |
| `results/figures/fig2_db_cooccurrence_heatmap.png` | Co-occurrence heatmap panel (4 x 3 grid) |
| `results/figures/fig3_db_network_comparison.png` | Co-occurrence/avoidance networks |
| `results/figures/fig4_netosis.png` | NETosis pathway diagram |
| `results/figures/fig5_orf8_mimicry.png` | ORF8 histone mimicry network |
| `results/figures/fig6_immunothrombosis.png` | Immunothrombosis cycle |
| `results/figures/db_comparison_summary.csv` | Full pairwise statistics table |

## Project Structure

```
├── app-covid.py                    # Streamlit dashboard (main application)
├── benchmark_databases.py          # Database benchmarking script
├── gen_figures.py                  # Figure generation from pipeline
├── histone_mark_mining.py          # Original pipeline (CLI)
├── app.py                          # Original Streamlit dashboard
├── CLAUDE.md                       # Agent context
├── README.md
├── Coronavirus_Vault/
│   ├── coronavirus26.bib           # Reference bibliography
│   ├── Molecular Sabotage_.md      # Histone damage mechanisms synthesis
│   └── Viral Histone Mimicry_.md   # Epigenetic suppression synthesis
├── results/
│   └── figures/                    # All generated figures
├── tests/
│   ├── test_pipeline.py            # Smoke tests
│   └── ...
├── docs/
│   └── PRD.md
└── .gitignore
```

---

## Legacy

The original pipeline (`histone_mark_mining.py`, `app.py`) covers graph construction, scale-free topology testing via power-law fitting, and KMeans chromatin state clustering. These are maintained for reproducibility.

---

## Testing

```bash
pytest tests/test_pipeline.py -v
```

9 smoke tests covering data mining, live fetch fallback, graph construction, and centrality bounds.

---

## Known Limitations

1. **Median-split binarization discards information.** Continuous variables lose magnitude information; low-variance features (e.g., Histome2 multi_mod) become degenerate.
2. **Small sample sizes.** 24-62 records per database limit statistical power for moderate effect sizes.
3. **Database-level vs. PTM-level features.** Writer counts and HMM scores capture different biology than the PTM-level marks (H3K4me3, H3K27ac) that show established co-occurrence in ChIP-seq data.
4. **HTML scraper fragility.** Histome2 scraper relies on undocumented CSS class names; a site redesign triggers fallback to simulated data.

---

AI Assistance Disclaimer: This codebase was developed with the assistance of Claude Code. While the AI provided code generation, debugging, and structural support, the human developer maintains full responsibility for reviewing, testing, and maintaining all content and functionality.

## References

- Kee J, et al. (2022) SARS-CoV-2 disrupts host epigenetic regulation via histone mimicry. *Nature* 610:381-388. doi:10.1038/s41586-022-05282-z
- de Vries F, et al. (2022) The role of extracellular histones in COVID-19. *J Intern Med* 293:275. doi:10.1111/joim.13585
- Ives AN, et al. (2026) Human coronavirus 229E infection alters histone proteoforms. *J Proteome Res* 25:2698-2708. doi:10.1021/acs.jproteome.5c01015
- Kaur R, et al. (2021) Histome2: a database of human histone PTMs. *Nucleic Acids Res* 49:D1024-D1030. doi:10.1093/nar/gkaa641
- Draizen EJ, et al. (2023) HistoneDB 2.0: a histone database with variants. *Nucleic Acids Res* 51:D348-D355. doi:10.1093/nar/gkac1011
- Henneman B, et al. (2018) ProHistoneDB: a database of proto-histones. *Nucleic Acids Res* 46:D348-D354. doi:10.1093/nar/gkx1014

---

## License

For academic and research use. Contact the author for other uses.
