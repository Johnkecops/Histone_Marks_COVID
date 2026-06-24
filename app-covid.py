#!/usr/bin/env python3
"""
Module: app-covid
Purpose: Streamlit dashboard for histone-virus interaction research.
         Visualizes molecular mechanisms of histone-mediated tissue damage,
         histone mark alterations during coronavirus infection, NETosis
         pathways, therapeutic countermeasures, and full bibliography.
Author: Dr. Arli Aditya Parikesit
Date: 2026
References:
    - Kee et al. (2022) Nature 610:381-388. doi:10.1038/s41586-022-05282-z
    - de Vries et al. (2022) J Intern Med 293:275. doi:10.1111/joim.13585
    - Ives et al. (2026) J Proteome Res 25:2698-2708. doi:10.1021/acs.jproteome.5c01015
    - York (2022) Nat Rev Microbiol 20:703. doi:10.1038/s41579-022-00815-9
    - Yang et al. (2022) J Immunol 208(Suppl 1):111.19. doi:10.4049/jimmunol.208.Supp.111.19
    - Gupta et al. (2025) Front Immunol 16:1596135. doi:10.3389/fimmu.2025.1596135
    - Schafer & Baric (2017) Pathogens 6:8. doi:10.3390/pathogens6010008
"""
import os
import sys
from pathlib import Path
from itertools import combinations

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import streamlit as st
from scipy.stats import poisson
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

with st.spinner("Importing pipeline modules..."):
    from histone_mark_mining import (
        data_mining,
        fetch_histome2_data,
        fetch_histonedb_data,
        fetch_prohistonedb_data,
    )

st.set_page_config(
    page_title="Histone-Virus Interaction Research",
    page_icon=None,
    layout="wide",
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_CORONAVIRUS_VAULT = Path(__file__).resolve().parent / "Coronavirus_Vault"
_HISTONE_COVID_DIR = Path(__file__).resolve().parent / "Histone Mark of COVID-19"

VIRUSES = ["SARS-CoV-2", "MERS-CoV", "SARS-CoV", "HCoV-229E"]

HISTONE_MARKS_RECORDS = [
    {
        "type": "Methylation",
        "code": "H3K9me3, H3K27me3",
        "virus": "SARS-CoV-2, MERS-CoV",
        "annotation": "Transcriptional Repression / Gene Silencing: ORF8 acts as a histone mimic to increase these marks, promoting chromatin compaction and silencing antiviral genes.",
    },
    {
        "type": "Methylation",
        "code": "H3K4me3",
        "virus": "SARS-CoV-2, SARS-CoV",
        "annotation": "Transcriptional Activation: Found at promoters of highly active genes; depleted at ISG promoters in MERS-CoV to prevent an antiviral state.",
    },
    {
        "type": "Methylation",
        "code": "H3K9me2",
        "virus": "SARS-CoV-2",
        "annotation": "Repressive Mark: Inversely correlates with interferon expression; used by the virus to inhibit the innate immune response.",
    },
    {
        "type": "Acetylation",
        "code": "H3K9ac, H3K14ac",
        "virus": "SARS-CoV-2",
        "annotation": "Transcriptional Activation: Necessary to mount a potent antiviral state; decreased by ORF8 to weaken the host interferon (IFN) response.",
    },
    {
        "type": "Acetylation",
        "code": "H3K27ac",
        "virus": "SARS-CoV-2",
        "annotation": "Regulation of ACE2: Associated with the activation of ACE2-related genes, impacting viral entry.",
    },
    {
        "type": "Acetylation",
        "code": "H2AX-K6ac",
        "virus": "HCoV-229E",
        "annotation": "DNA Damage Response: Regulates H2AX release and chromatin reorganization after DNA double-strand breaks; decreased abundance during infection.",
    },
    {
        "type": "Citrullination",
        "code": "Cit-H3 (Citrullinated H3)",
        "virus": "SARS-CoV-2",
        "annotation": "NETosis Marker: Indicates neutrophil extracellular trap formation; elevated levels drive vascular damage, ARDS, and predict 28-day mortality.",
    },
    {
        "type": "Phosphorylation",
        "code": "H2AC11-T17ph",
        "virus": "HCoV-229E",
        "annotation": "Oxidative Stress: Identified as a specific marker for virus-induced oxidative stress in human lung fibroblasts.",
    },
    {
        "type": "Phosphorylation",
        "code": "H3S10ph",
        "virus": "SARS-CoV-2",
        "annotation": "Cell Cycle / Activation: Monitored during ORF8 expression studies.",
    },
    {
        "type": "Glutathionylation",
        "code": "H3C1 (at C97/C111)",
        "virus": "HCoV-229E",
        "annotation": "Structural Destabilization: Known to destabilize the nucleosome complex.",
    },
    {
        "type": "Truncation",
        "code": "C-terminally truncated H2A",
        "virus": "HCoV-229E",
        "annotation": "Structural Regulation: Truncation normally destabilizes nucleosomes and removes PTM sites; abundance is reduced during infection.",
    },
    {
        "type": "Truncation",
        "code": "N-terminally truncated H3",
        "virus": "HCoV-229E",
        "annotation": "Structural Regulation: Cleaved at R27 or K28; abundance reduced during infection to potentially modulate DNA accessibility.",
    },
]

ADDITIONAL_PTMS = [
    "H3K18ac", "H3K23ac", "H2AK5ac",
    "H3K4me1", "H3K36me3", "H3.3K27me3", "H4K20me1",
    "H4R20me", "H4K21me", "H4R24me", "H4R68me",
    "H3K28me", "H3K37me", "H3R50me", "H3K57me",
    "H3K80me", "H3R84me",
]

# ---------------------------------------------------------------------------
# Nucleosome co-occurrence simulation
# ---------------------------------------------------------------------------

NUCLEOSOME_MARKS = ["H3K4me3", "H3K27ac", "H3K9me3", "H3K36me3", "H2AK5ac"]

_NUCLEOSOME_CITATION = (
    "Co-occurrence and avoidance patterns modeled after known biology: "
    "H3K4me3 and H3K27ac co-occur at active promoters; "
    "H3K9me3 (repressive) avoids H3K4me3 and H3K27ac; "
    "H3K36me3 is independent (gene body); "
    "H2AK5ac weakly co-occurs with H3K27ac."
)

@st.cache_data
def generate_nucleosome_data(n_nuc: int = 200, seed: int = 42) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    regions = [f"nuc_{i}" for i in range(n_nuc)]

    h3k4me3 = rng.binomial(1, 0.35, n_nuc)
    h3k27ac = np.array([
        1 if (h3k4me3[i] == 1 and rng.random() < 0.75)
          or (h3k4me3[i] == 0 and rng.random() < 0.15)
        else 0
        for i in range(n_nuc)
    ])
    active = (h3k4me3 | h3k27ac).astype(bool)
    h3k9me3 = np.array([
        1 if (not active[i] and rng.random() < 0.40)
          or (active[i] and rng.random() < 0.06)
        else 0
        for i in range(n_nuc)
    ])
    h3k36me3 = rng.binomial(1, 0.25, n_nuc)
    h2ak5ac = np.array([
        1 if (h3k27ac[i] == 1 and rng.random() < 0.45)
          or (h3k27ac[i] == 0 and rng.random() < 0.10)
        else 0
        for i in range(n_nuc)
    ])

    return pd.DataFrame({
        "region": regions,
        "H3K4me3": h3k4me3.astype(int),
        "H3K27ac": h3k27ac.astype(int),
        "H3K9me3": h3k9me3.astype(int),
        "H3K36me3": h3k36me3.astype(int),
        "H2AK5ac": h2ak5ac.astype(int),
    })


@st.cache_data
def compute_cooccurrence_stats(
    df_binary: pd.DataFrame,
    marks: list[str],
    alpha: float = 0.05,
):
    n = len(df_binary)
    m = len(marks)
    occurrence = df_binary[marks].mean()

    obs_mat = np.zeros((m, m))
    exp_mat = np.zeros((m, m))
    lfc_mat = np.zeros((m, m))
    p_co = np.ones((m, m))
    p_av = np.ones((m, m))
    test_labels = []
    test_stats = []

    for i, mi in enumerate(marks):
        for j, mj in enumerate(marks):
            if i == j:
                obs_mat[i, j] = occurrence[mi]
                exp_mat[i, j] = occurrence[mi]
                continue
            if i > j:
                continue

            both = ((df_binary[mi] == 1) & (df_binary[mj] == 1)).sum()
            pi = occurrence[mi]
            pj = occurrence[mj]
            exp_cnt = n * pi * pj
            prop = both / n

            obs_mat[i, j] = prop
            obs_mat[j, i] = prop
            exp_mat[i, j] = exp_cnt / n
            exp_mat[j, i] = exp_cnt / n

            if exp_cnt > 0 and both > 0:
                lfc_mat[i, j] = np.log2(both / exp_cnt)
                lfc_mat[j, i] = lfc_mat[i, j]
            elif exp_cnt > 0 and both == 0:
                lfc_mat[i, j] = -10.0
                lfc_mat[j, i] = -10.0
            else:
                lfc_mat[i, j] = 0.0
                lfc_mat[j, i] = 0.0

            if exp_cnt > 0:
                p_co_val = poisson.sf(both - 1, exp_cnt)
                p_av_val = poisson.cdf(both, exp_cnt)
            else:
                p_co_val = 0.0 if both > 0 else 1.0
                p_av_val = 1.0
            p_co[i, j] = p_co_val
            p_co[j, i] = p_co_val
            p_av[i, j] = p_av_val
            p_av[j, i] = p_av_val

            test_labels.append(f"{mi} vs {mj}")
            test_stats.append({
                "pair": f"{mi} vs {mj}",
                "mark_i": mi,
                "mark_j": mj,
                "observed": int(both),
                "expected": round(exp_cnt, 2),
                "occurrence_i": round(float(pi), 3),
                "occurrence_j": round(float(pj), 3),
                "n_regions": n,
                "p_cooccurrence": round(p_co_val, 6),
                "p_avoidance": round(p_av_val, 6),
                "log2_fold_change": round(lfc_mat[i, j], 3),
            })

    df_test = pd.DataFrame(test_stats)
    pvals_co = df_test["p_cooccurrence"].values
    pvals_av = df_test["p_avoidance"].values
    reject_co, pcor_co, _, _ = multipletests(pvals_co, method="fdr_bh")
    reject_av, pcor_av, _, _ = multipletests(pvals_av, method="fdr_bh")
    df_test["p_cooccurrence_adj"] = pcor_co
    df_test["p_avoidance_adj"] = pcor_av
    df_test["cooccurrence_significant"] = (pcor_co < alpha) & (df_test["observed"] > df_test["expected"])
    df_test["avoidance_significant"] = (pcor_av < alpha) & (df_test["observed"] < df_test["expected"])

    return {
        "marks": marks,
        "occurrence": occurrence,
        "observed_matrix": obs_mat,
        "expected_matrix": exp_mat,
        "log2fc_matrix": lfc_mat,
        "p_cooccurrence": p_co,
        "p_avoidance": p_av,
        "test_table": df_test,
        "alpha": alpha,
    }


# ---------------------------------------------------------------------------
# Co-occurrence visualization helpers
# ---------------------------------------------------------------------------

def _draw_occurrence_barchart(occ: pd.Series):
    fig, ax = plt.subplots(figsize=(8, 4))
    colors = ["#3498db", "#e74c3c", "#2ecc71", "#f39c12", "#9b59b6"]
    bars = ax.bar(range(len(occ)), occ.values * 100, color=colors[:len(occ)],
                  alpha=0.8, edgecolor="white", linewidth=1.2)
    ax.set_xticks(range(len(occ)))
    ax.set_xticklabels(occ.index, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("Occurrence frequency (%)", fontsize=10)
    ax.set_title("Histone Mark Occurrence Across Nucleosomes", fontsize=12, fontweight="bold")
    for bar, val in zip(bars, occ.values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                f"{val * 100:.1f}%", ha="center", va="bottom", fontsize=8)
    ax.set_ylim(0, max(occ.values) * 100 * 1.2)
    st.pyplot(fig)
    plt.close(fig)


def _draw_cooccurrence_heatmap(stats: dict):
    marks = stats["marks"]
    obs = stats["observed_matrix"]
    exp = stats["expected_matrix"]
    p_co = stats["p_cooccurrence"]
    p_av = stats["p_avoidance"]
    alpha = stats["alpha"]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    titles = ["Observed co-occurrence", "Expected (independence)", "Significance direction"]
    matrices = [obs, exp, obs - exp]
    cmaps = ["YlOrRd", "YlOrRd", "RdYlBu_r"]
    vmins = [0, 0, None]
    vmaxs = [None, None, None]

    for ax, mat, title, cmap, vmin, vmax in zip(axes, matrices, titles, cmaps, vmins, vmaxs):
        im = ax.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax, aspect="equal")
        ax.set_xticks(range(len(marks)))
        ax.set_yticks(range(len(marks)))
        ax.set_xticklabels(marks, rotation=45, ha="right", fontsize=8)
        ax.set_yticklabels(marks, fontsize=8)
        ax.set_title(title, fontsize=10, fontweight="bold")
        plt.colorbar(im, ax=ax, shrink=0.8)

        for i in range(len(marks)):
            for j in range(len(marks)):
                val = mat[i, j]
                if i == j:
                    ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                            fontsize=7, color="gray")
                    continue
                stars = ""
                if p_co[i, j] < alpha:
                    stars += "*"
                if p_av[i, j] < alpha:
                    stars += "\u2020"
                txt = f"{val:.2f}{stars}" if stars else f"{val:.2f}"
                ax.text(j, i, txt, ha="center", va="center",
                        fontsize=7, color="black" if abs(val - 0.5) < 0.3 else "white")

    fig.suptitle("Histone Mark Co-occurrence Analysis (* p_co < alpha; \u2020 p_avoid < alpha)",
                 fontsize=11, fontweight="bold", y=1.02)
    fig.tight_layout()
    st.pyplot(fig)
    plt.close(fig)


def _draw_cooccurrence_network(stats: dict):
    marks = stats["marks"]
    obs = stats["observed_matrix"]
    p_co = stats["p_cooccurrence"]
    p_av = stats["p_avoidance"]
    alpha = stats["alpha"]

    G = nx.Graph()
    for m in marks:
        G.add_node(m, size=600)

    for i, mi in enumerate(marks):
        for j, mj in enumerate(marks):
            if i >= j:
                continue
            if p_co[i, j] < alpha and obs[i, j] > stats["expected_matrix"][i, j]:
                G.add_edge(mi, mj, color="#2ecc71", weight=2 + 5 * obs[i, j])
            elif p_av[i, j] < alpha and obs[i, j] < stats["expected_matrix"][i, j]:
                G.add_edge(mi, mj, color="#e74c3c", weight=2, style="dashed")

    if G.number_of_edges() == 0:
        st.info("No significant co-occurrence or avoidance relationships at the current threshold.")
        return

    fig, ax = plt.subplots(figsize=(8, 6))
    pos = nx.circular_layout(G)
    edge_colors = [G.edges[e]["color"] for e in G.edges()]
    edge_widths = [G.edges[e]["weight"] for e in G.edges()]
    edge_styles = [G.edges[e].get("style", "solid") for e in G.edges()]

    nx.draw_networkx_nodes(G, pos, ax=ax, node_color="#bdc3c7",
                           node_size=700, alpha=0.85)
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=9, font_weight="bold")
    for (u, v), color, width, style in zip(G.edges(), edge_colors, edge_widths, edge_styles):
        nx.draw_networkx_edges(
            G, pos, ax=ax, edgelist=[(u, v)],
            edge_color=color, width=width, style=style, alpha=0.7,
        )

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color="#2ecc71", lw=3, label="Co-occurrence (p_adj < alpha)"),
        Line2D([0], [0], color="#e74c3c", lw=2, linestyle="dashed", label="Avoidance (p_adj < alpha)"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=8)
    ax.set_title("Co-occurrence and Avoidance Network", fontsize=12, fontweight="bold")
    ax.axis("off")
    st.pyplot(fig)
    plt.close(fig)


def _draw_virus_mark_cooccurrence(df_marks: pd.DataFrame):
    mark_types = df_marks["type"].value_counts()
    pairs = list(combinations(mark_types.index, 2))

    n = len(df_marks)
    rows = []
    for t1, t2 in pairs:
        n1 = mark_types[t1]
        n2 = mark_types[t2]
        both = ((df_marks["type"] == t1) | (df_marks["type"] == t2)).sum()
        exp_cnt = n * (n1 / n) * (n2 / n)
        p_val = poisson.sf(both - 1, exp_cnt) if exp_cnt > 0 else 1.0
        rows.append({
            "Mark type 1": t1,
            "Mark type 2": t2,
            "Count (union)": both,
            "Expected (independence)": round(exp_cnt, 2),
            "Poisson p-value": round(p_val, 4),
            "Significant (p < 0.05)": p_val < 0.05,
        })

    df_virus_co = pd.DataFrame(rows).sort_values("Poisson p-value")
    st.dataframe(df_virus_co, use_container_width=True)


# ---------------------------------------------------------------------------
# Co-occurrence data source selection and preparation
# ---------------------------------------------------------------------------

COOC_DATA_SOURCES = ["Simulated", "Histome2", "HistoneDB 2.0", "ProHistoneDB"]

COOC_SOURCE_META = {
    "Simulated": {
        "marks": NUCLEOSOME_MARKS,
        "descriptions": {
            "H3K4me3": "Active promoter (trimethylation of H3 Lys4)",
            "H3K27ac": "Active enhancer/promoter (acetylation of H3 Lys27)",
            "H3K9me3": "Repressive heterochromatin (trimethylation of H3 Lys9)",
            "H3K36me3": "Gene body elongation (trimethylation of H3 Lys36)",
            "H2AK5ac": "Minor acetylation (acetylation of H2A Lys5)",
        },
        "binary_method": "Binomial simulation with known interaction parameters",
        "n_regions": 200,
        "citation": _NUCLEOSOME_CITATION,
    },
    "Histome2": {
        "marks": ["has_writer", "has_eraser", "multi_mod"],
        "descriptions": {
            "has_writer": "Above-median writer enzyme count at PTM site",
            "has_eraser": "Above-median eraser enzyme count at PTM site",
            "multi_mod": "Above-median number of modification types at site",
        },
        "binary_method": "Median split of discrete count variables (writer/eraser/mod_type counts)",
        "n_regions": "Variable (PTM sites from Histome2)",
        "citation": (
            "Human histone PTM sites fetched from Histome2. "
            "Continuous mark variables (n_writers, n_erasers, n_mod_types) "
            "binarized by median split."
        ),
    },
    "HistoneDB 2.0": {
        "marks": ["high_HMM_score", "many_sequences", "broad_taxonomy"],
        "descriptions": {
            "high_HMM_score": "Above-median HMM classification confidence score",
            "many_sequences": "Above-median number of sequences for this variant",
            "broad_taxonomy": "Above-median number of distinct taxonomic IDs",
        },
        "binary_method": "Median split of continuous/ordinal variables (HMM score, sequence count, taxon count)",
        "n_regions": "Variable (histone variants from HistoneDB)",
        "citation": (
            "Histone variant sequences with HMM classification scores from "
            "HistoneDB 2.0 (NCBI). Continuous marks binarized by median split."
        ),
    },
    "ProHistoneDB": {
        "marks": ["many_sequences", "broad_taxonomy", "long_sequence"],
        "descriptions": {
            "many_sequences": "Above-median number of sequences in proto-histone category",
            "broad_taxonomy": "Above-median number of distinct taxonomic groups",
            "long_sequence": "Above-median mean sequence length (amino acids)",
        },
        "binary_method": "Median split of continuous variables (sequence count, taxon count, mean length)",
        "n_regions": "Variable (proto-histone categories from ProHistoneDB)",
        "citation": (
            "Proto-histone categories from ProHistoneDB across Archaea, "
            "Bacteria, and viruses. Continuous marks binarized by median split."
        ),
    },
}

_COOC_LABEL_MAP = {
    "Histome2": {"mark1": "has_writer", "mark2": "has_eraser", "mark3": "multi_mod"},
    "HistoneDB 2.0": {"mark1": "high_HMM_score", "mark2": "many_sequences", "mark3": "broad_taxonomy"},
    "ProHistoneDB": {"mark1": "many_sequences", "mark2": "broad_taxonomy", "mark3": "long_sequence"},
}


@st.cache_data(show_spinner="Loading co-occurrence data...")
def prepare_cooccurrence_data(source: str):
    if source == "Simulated":
        df = generate_nucleosome_data()
        marks = [c for c in NUCLEOSOME_MARKS if c in df.columns]
        return df, marks, "Simulated (seed=42)"

    fetch_map = {
        "Histome2": fetch_histome2_data,
        "HistoneDB 2.0": fetch_histonedb_data,
        "ProHistoneDB": fetch_prohistonedb_data,
    }

    try:
        raw = fetch_map[source]()
        if raw is not None and not raw.empty:
            df = pd.DataFrame({"region": raw["region"].astype(str)})
            label_map = _COOC_LABEL_MAP[source]
            binary_marks = []
            for col_key, new_name in label_map.items():
                if col_key not in raw.columns:
                    continue
                vals = pd.to_numeric(raw[col_key], errors="coerce").dropna()
                if len(vals) < 2:
                    continue
                median = vals.median()
                df[new_name] = (raw[col_key].astype(float) > median).astype(int)
                binary_marks.append(new_name)
            if len(binary_marks) >= 2:
                return df, binary_marks, f"{source} (live, binarized by median split)"
    except Exception as e:
        st.warning(f"{source} fetch failed: {e}. Falling back to simulated data.")

    df = generate_nucleosome_data()
    marks = [c for c in NUCLEOSOME_MARKS if c in df.columns]
    return df, marks, f"Simulated ({source} unavailable)"


MECHANISMS = [
    {
        "id": 1,
        "title": "Direct Membrane Disruption (Physical Cytolysis)",
        "subtitle": "Histone H4 as a molecular detergent",
        "summary": (
            "Extracellular histones, specifically Histone H4, behave as molecular "
            "detergents. Their high positive charge density allows binding to negatively "
            "charged phospholipids of host membranes, forming lytic pores."
        ),
        "details": [
            "Plasma Membrane Permeabilization: Histone H4 intercalates into plasma "
            "membranes of vascular smooth muscle and endothelial cells, forming lytic "
            "pores that cause immediate cytolysis.",
            "Mitochondrial Outer Membrane Permeabilization (MOMP): Histones disrupt "
            "the mitochondrial boundary, causing leakage of cytochrome c into the cytosol.",
            "Apoptosome Activation: Cytochrome c translocation triggers apoptosome "
            "assembly, activating the caspase cascade (Caspase-9, then Caspase-3/7), "
            "mandating programmed cell death in uninfected bystander cells.",
        ],
        "color": "#e74c3c",
    },
    {
        "id": 2,
        "title": "Ionic Dyshomeostasis (Chemical Chaos)",
        "subtitle": "Histone H1 and calcium overload",
        "summary": (
            "Histones (particularly H1) disrupt cellular calcium (Ca2+) regulation, "
            "bypassing highly regulated signaling pathways and inducing chemical chaos."
        ),
        "details": [
            "Nonselective cation channel activation: Histones induce massive influx "
            "of extracellular Ca2+ through cell surface channels.",
            "ER calcium release: Histones trigger Ca2+ release from intracellular "
            "stores within the endoplasmic reticulum.",
            "Consequences: Signaling paralysis (bypassing endothelium-dependent "
            "vasodilation), ER stress (unfolded protein response), and contractile "
            "failure in cardiomyocytes contributing to arrhythmias.",
        ],
        "color": "#f39c12",
    },
    {
        "id": 3,
        "title": "Innate Immune Overdrive (TLR Activation)",
        "subtitle": "Histones as pseudo-pathogenic signals",
        "summary": (
            "The innate immune system recognizes extracellular histones via Toll-like "
            "Receptors (TLR2, TLR4, TLR9), mistakenly identifying them as viral threats."
        ),
        "details": [
            "Histone binding to TLRs recruits the adapter protein MyD88, which in turn "
            "recruits IRAK-4 and TRAF6.",
            "This signaling hub activates the IkappaB kinase (IKK) complex and mitogen-"
            "activated protein kinases (MAPK).",
            "Engagement of NF-kB and AP-1 transcription factors drives the cytokine storm "
            "(IL-6, TNF-alpha release).",
        ],
        "color": "#9b59b6",
    },
    {
        "id": 4,
        "title": "Inflammasome and Pyroptosis",
        "subtitle": "Explosive inflammatory cell death",
        "summary": (
            "Histones are sensed internally by NOD2 and NLRP3 receptors, leading to "
            "inflammasome assembly and pyroptosis, a highly inflammatory cell death."
        ),
        "details": [
            "NLRP3 inflammasome assembly facilitates Caspase-1 activation.",
            "Caspase-1 cleaves pro-IL-1beta and pro-IL-18 into their active forms.",
            "Pyroptosis results in cell swelling and violent cytolysis, releasing massive "
            "pro-inflammatory signals, unlike apoptosis which is immunologically quiet.",
        ],
        "color": "#e67e22",
    },
    {
        "id": 5,
        "title": "Complement System Sabotage",
        "subtitle": "Histones H3/H4 inhibit C3 and C5",
        "summary": (
            "Histones H3 and H4 bind complement protein C4, inhibiting functional "
            "activities of C3 and C5 while not preventing C4 cleavage."
        ),
        "details": [
            "The histone-C4 complex specifically inhibits C3 and C5 functional activities, "
            "impairing pathogen clearance.",
            "This creates a vicious circle: impaired viral clearance -> continued "
            "cell death -> more histone release -> immunothrombosis.",
            "Immunothrombosis: pathological clot formation in microvasculature driven "
            "by intersection of inflammation and coagulation factors.",
        ],
        "color": "#1abc9c",
    },
]

REFERENCES = [
    {
        "key": "Kee2022",
        "authors": "Kee J, Thudium S, Renner DM, Glastad K, Palozola K, Zhang Z, "
                   "Li Y, Lan Y, Cesare J, Poleshko A, Kiseleva AA, Truitt R, "
                   "Cardenas-Diaz FL, Zhang X, Xie X, Kotton DN, Alysandratos KD, "
                   "Epstein JA, Shi PY, Yang W, Morrisey E, Garcia BA, Berger SL, "
                   "Weiss SR, Korb E",
        "year": "2022",
        "title": "SARS-CoV-2 disrupts host epigenetic regulation via histone mimicry",
        "journal": "Nature",
        "volume": "610",
        "pages": "381-388",
        "doi": "10.1038/s41586-022-05282-z",
        "pmid": "36198800",
        "url": "https://www.nature.com/articles/s41586-022-05282-z",
    },
    {
        "key": "deVries2022",
        "authors": "de Vries F, Huckriede J, Wichapong K, Reutelingsperger C, Nicolaes GAF",
        "year": "2022",
        "title": "The role of extracellular histones in COVID-19",
        "journal": "Journal of Internal Medicine",
        "volume": "293",
        "pages": "275",
        "doi": "10.1111/joim.13585",
        "pmid": "36382685",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10108027/",
    },
    {
        "key": "Ives2026",
        "authors": "Ives AN, Thibert S, Berger MR, Gaffrey MJ, Mitchell HD, "
                   "Williams SM, Zhou M, Waters KM, Sims AC, Zhang T",
        "year": "2026",
        "title": "Human Coronavirus 229E Infection Alters Histone Proteoforms",
        "journal": "Journal of Proteome Research",
        "volume": "25",
        "pages": "2698-2708",
        "doi": "10.1021/acs.jproteome.5c01015",
        "pmid": "42047537",
        "url": "https://pubs.acs.org/doi/10.1021/acs.jproteome.5c01015",
    },
    {
        "key": "York2022",
        "authors": "York A",
        "year": "2022",
        "title": "Histone mimicry by SARS-CoV-2",
        "journal": "Nature Reviews Microbiology",
        "volume": "20",
        "pages": "703",
        "doi": "10.1038/s41579-022-00815-9",
        "pmid": "36207440",
        "url": "https://www.nature.com/articles/s41579-022-00815-9",
    },
    {
        "key": "Yang2022",
        "authors": "Yang X, Julian J, Schnell T, Albrecht H, Owens W, Yang XS",
        "year": "2022",
        "title": "Epigenetic regulation of inflammatory genes involving histone "
                 "methylation is associated with hyperimmune response in COVID-19 patients",
        "journal": "The Journal of Immunology",
        "volume": "208",
        "pages": "111.19",
        "doi": "10.4049/jimmunol.208.Supp.111.19",
        "url": "https://dx.doi.org/10.4049/jimmunol.208.Supp.111.19",
    },
    {
        "key": "Gupta2025",
        "authors": "Gupta S, Hemeg HA, Afrin F",
        "year": "2025",
        "title": "Immuno-epigenetic paradigms in coronavirus infection",
        "journal": "Frontiers in Immunology",
        "volume": "16",
        "pages": "1596135",
        "doi": "10.3389/fimmu.2025.1596135",
        "pmid": "40969755",
        "url": "https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2025.1596135",
    },
    {
        "key": "Schafer2017",
        "authors": "Schafer A, Baric RS",
        "year": "2017",
        "title": "Epigenetic landscape during coronavirus infection",
        "journal": "Pathogens",
        "volume": "6",
        "pages": "8",
        "doi": "10.3390/pathogens6010008",
        "pmid": "28212305",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5371896/",
    },
]

# ---------------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------------

df_marks = pd.DataFrame(HISTONE_MARKS_RECORDS)
df_refs = pd.DataFrame(REFERENCES)

# ---------------------------------------------------------------------------
# Page header
# ---------------------------------------------------------------------------

st.title("Histone-Virus Interaction Research")
st.caption(
    "Molecular mechanisms of histone-mediated tissue damage during coronavirus infection, "
    "epigenetic reprogramming by viral histone mimicry, and therapeutic countermeasures.  "
    "Integrates knowledge from Coronavirus_Vault literature synthesis."
)

# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------

st.sidebar.title("Navigation & Filters")

view_mode = st.sidebar.radio(
    "View mode",
    options=["Guided Tour", "Reference Browser", "Mechanism Explorer"],
    index=0,
    help=(
        "Guided Tour: stepwise through mechanisms. "
        "Reference Browser: full bibliography. "
        "Mechanism Explorer: interactive deep-dive."
    ),
)

st.sidebar.markdown("---")

filter_virus = st.sidebar.multiselect(
    "Filter by virus",
    options=VIRUSES,
    default=VIRUSES,
    help="Only show histone marks associated with selected viruses.",
)

st.sidebar.markdown("---")

cooc_source = st.sidebar.radio(
    "Co-occurrence data source",
    options=COOC_DATA_SOURCES,
    index=0,
    help=(
        "Simulated: synthetic nucleosomes with 5 binary marks and known biology. "
        "Histome2: human PTM site writer/eraser counts binarized by median split. "
        "HistoneDB 2.0: histone variant HMM scores and counts binarized. "
        "ProHistoneDB: proto-histone category metrics binarized."
    ),
)

st.sidebar.markdown("---")
st.sidebar.caption(
    "Dr. Arli Aditya Parikesit  \n"
    "i3L University Jakarta  \n"
    "[ORCID 0000-0001-8716-3926](https://orcid.org/0000-0001-8716-3926)"
)

# ---------------------------------------------------------------------------
# Filtered data
# ---------------------------------------------------------------------------

_mask = df_marks["virus"].apply(
    lambda v: any(fv in v for fv in filter_virus)
) if filter_virus else pd.Series(True, index=df_marks.index)
df_marks_filtered = df_marks[_mask].reset_index(drop=True)

# ---------------------------------------------------------------------------
# Co-occurrence analysis data (shared across tabs)
# ---------------------------------------------------------------------------

df_nuc, marks_nuc, cooc_source_label = prepare_cooccurrence_data(cooc_source)

cooc_alpha = st.sidebar.slider(
    "Co-occurrence alpha (FDR)", 0.001, 0.10, 0.05, 0.005,
    key="cooc_alpha",
    help="Adjusted p-value threshold for co-occurrence and avoidance significance.",
)

cooc_stats = compute_cooccurrence_stats(df_nuc, marks_nuc, alpha=cooc_alpha)
df_cooc_test = cooc_stats["test_table"]

# ---------------------------------------------------------------------------
# Helper: mechanism card
# ---------------------------------------------------------------------------

def _render_mechanism_card(m):
    with st.container(border=True):
        cols = st.columns([0.05, 0.95])
        with cols[0]:
            st.markdown(
                f"<div style='background:{m['color']};width:6px;height:100%;"
                f"border-radius:3px;'></div>",
                unsafe_allow_html=True,
            )
        with cols[1]:
            st.subheader(f"Mechanism {m['id']}: {m['title']}")
            st.caption(m["subtitle"])
            st.markdown(m["summary"])
            with st.expander("Detailed mechanism"):
                for d in m["details"]:
                    st.markdown(f"- {d}")


# ---------------------------------------------------------------------------
# Helper: NETosis pathway figure
# ---------------------------------------------------------------------------

def _draw_netosis_pathway():
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7)
    ax.axis("off")

    steps = [
        (1.5, 6.0, "1. Activation", "Neutrophil triggered\nby viral PAMPs or\npro-inflammatory cytokines", "#3498db"),
        (4.0, 6.0, "2. ROS Production", "NADPH oxidase generates\nReactive Oxygen Species", "#e74c3c"),
        (6.5, 6.0, "3. Citrullination", "PAD4 converts histone\narginine to citrulline\nchromatin decondenses", "#e67e22"),
        (1.5, 3.5, "4. Dissolution", "Nuclear envelope breaks;\nchromatin mixes with\nMPO and Neutrophil Elastase", "#9b59b6"),
        (4.0, 3.5, "5. Externalization", "Cell membrane ruptures;\nNET released into\nextracellular space", "#1abc9c"),
        (6.5, 3.5, "6. Unmasking", "Host DNase I digests\nNET DNA scaffold;\nhistones regain cytotoxicity", "#2c3e50"),
    ]

    for x, y, title, desc, color in steps:
        ax.add_patch(mpatches.FancyBboxPatch(
            (x - 0.7, y - 0.5), 1.4, 1.0,
            boxstyle="round,pad=0.15",
            facecolor=color, alpha=0.15, edgecolor=color, linewidth=2,
        ))
        ax.text(x, y + 0.15, title, ha="center", va="bottom", fontsize=9,
                fontweight="bold", color=color)
        ax.text(x, y - 0.15, desc, ha="center", va="top", fontsize=7,
                color="gray")

    arrows = [
        (2.2, 6.5, 3.3, 6.5),
        (4.7, 6.5, 5.8, 6.5),
        (2.2, 5.5, 2.2, 4.0),
        (4.7, 5.5, 4.7, 4.0),
        (6.5, 5.5, 6.5, 4.0),
        (2.2, 3.0, 3.3, 3.0),
        (4.7, 3.0, 5.8, 3.0),
    ]
    for x1, y1, x2, y2 in arrows:
        ax.annotate(
            "", xy=(x2, y2), xytext=(x1, y1),
            arrowprops=dict(arrowstyle="->", color="gray", lw=1.5, alpha=0.6),
        )

    ax.set_title("NETosis Pathway: The Double-Edged Sword", fontsize=13, fontweight="bold")
    st.pyplot(fig)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Helper: histone mimicry network graph
# ---------------------------------------------------------------------------

def _draw_mimicry_network():
    G = nx.DiGraph()

    G.add_node("SARS-CoV-2 ORF8", color="#e74c3c", size=800, type="viral")
    G.add_node("ARKS motif\n(H3 mimic)", color="#f39c12", size=600, type="motif")
    G.add_node("Host chromatin", color="#3498db", size=700, type="host")
    G.add_node("H3K9me3 / H3K27me3\n(repressive marks)", color="#2c3e50", size=500, type="mark")
    G.add_node("IFN response genes\n(silenced)", color="#95a5a6", size=500, type="target")
    G.add_node("H3K9ac / H3K14ac\n(active marks depleted)", color="#7f8c8d", size=500, type="mark")
    G.add_node("Chromatin\ncompaction", color="#e67e22", size=550, type="outcome")

    edges = [
        ("SARS-CoV-2 ORF8", "ARKS motif\n(H3 mimic)", "mimics"),
        ("ARKS motif\n(H3 mimic)", "Host chromatin", "associates with"),
        ("Host chromatin", "H3K9me3 / H3K27me3\n(repressive marks)", "increases"),
        ("H3K9me3 / H3K27me3\n(repressive marks)", "IFN response genes\n(silenced)", "silences"),
        ("Host chromatin", "H3K9ac / H3K14ac\n(active marks depleted)", "decreases"),
        ("H3K9ac / H3K14ac\n(active marks depleted)", "IFN response genes\n(silenced)", "contributes to"),
        ("Host chromatin", "Chromatin\ncompaction", "promotes"),
        ("Chromatin\ncompaction", "IFN response genes\n(silenced)", "enhances"),
    ]

    for u, v, label in edges:
        G.add_edge(u, v, label=label)

    fig, ax = plt.subplots(figsize=(10, 7))
    pos = nx.spring_layout(G, k=1.2, seed=42, iterations=50)

    node_colors = [G.nodes[n]["color"] for n in G.nodes()]
    node_sizes = [G.nodes[n]["size"] for n in G.nodes()]

    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=node_sizes, alpha=0.85)
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=8, font_weight="bold")
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="gray",
                           alpha=0.4, arrows=True, arrowstyle="-|>",
                           arrowsize=15, connectionstyle="arc3,rad=0.1")

    ax.set_title(
        "ORF8 Histone Mimicry: Molecular Mechanism of Epigenetic Subversion",
        fontsize=12, fontweight="bold",
    )
    ax.axis("off")
    st.pyplot(fig)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Helper: immunothrombosis cycle figure
# ---------------------------------------------------------------------------

def _draw_vicious_cycle():
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.axis("off")

    nodes = {
        "Viral Infection": (0.0, 1.5),
        "Cell Death\n& Histone Release": (1.5, 0.5),
        "Complement\nInhibition (C3/C5)": (1.5, -0.8),
        "Impaired\nViral Clearance": (0.0, -1.5),
        "More Cell Death": (-1.5, -0.8),
        "Histone Release": (-1.5, 0.5),
    }

    for label, (x, y) in nodes.items():
        ax.add_patch(mpatches.FancyBboxPatch(
            (x - 0.7, y - 0.3), 1.4, 0.6,
            boxstyle="round,pad=0.1",
            facecolor="#e74c3c" if "Histone" in label or "Cell Death" in label
                    else "#f39c12" if "Complement" in label
                    else "#3498db",
            alpha=0.15, edgecolor="gray", linewidth=1.5,
        ))
        ax.text(x, y, label, ha="center", va="center", fontsize=8, fontweight="bold")

    cycle_edges = [
        ("Viral Infection", "Cell Death\n& Histone Release"),
        ("Cell Death\n& Histone Release", "Complement\nInhibition (C3/C5)"),
        ("Complement\nInhibition (C3/C5)", "Impaired\nViral Clearance"),
        ("Impaired\nViral Clearance", "More Cell Death"),
        ("More Cell Death", "Histone Release"),
        ("Histone Release", "Viral Infection"),
    ]

    for u, v in cycle_edges:
        if u in nodes and v in nodes:
            ax.annotate(
                "", xy=nodes[v], xytext=nodes[u],
                arrowprops=dict(
                    arrowstyle="->", color="#e74c3c", lw=2,
                    alpha=0.6, connectionstyle="arc3,rad=0.3",
                ),
            )

    ax.set_title(
        "Immunothrombosis: The Vicious Circle",
        fontsize=13, fontweight="bold", color="#c0392b",
    )
    st.pyplot(fig)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Tabs
# ---------------------------------------------------------------------------

tabs = st.tabs([
    "Overview",
    "Histone Marks in Infection",
    "Occurrence",
    "Co-occurrence",
    "Avoidance",
    "Molecular Mechanisms",
    "NETosis & Pathways",
    "Clinical Countermeasures",
    "References",
])

# ===== Tab 1: Overview ===================================================
with tabs[0]:
    st.header("The Biological Jekyll and Hyde of Histones")

    st.markdown(
        "In molecular pathology, histones are recognized primarily as the structural "
        "'spools' of the nucleosome, essential for DNA packaging and epigenetic "
        "regulation. However, during severe viral infections such as COVID-19, "
        "histones undergo a transformation from essential nuclear proteins to "
        "cytotoxic agents."
    )

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Internal Role ('The Hero')")
        st.markdown(
            "- **Location:** Nucleus (bound to DNA)\n"
            "- **Function:** DNA packaging & epigenetic regulation\n"
            "- **Impact:** Essential for genetic stability and transcriptional control\n"
            "- **Regulation:** Post-translational modifications (PTMs) orchestrate "
            "chromatin state"
        )

    with col2:
        st.subheader("External Role ('The Villain')")
        st.markdown(
            "- **Location:** Extracellular space (interstitium / plasma)\n"
            "- **Function:** Cytotoxic DAMP & pro-inflammatory signal\n"
            "- **Impact:** Driver of cytolysis, ARDS, and immunothrombosis\n"
            "- **Mechanisms:** Membrane disruption, TLR activation, complement sabotage"
        )

    st.divider()

    st.subheader("Core Hypothesis")
    st.markdown(
        "SARS-CoV-2 ORF8 protein acts as a **histone mimic**, imitating the "
        "'ARKS' motif of Histone H3. This mimicry associates with chromatin, "
        "disrupts critical host PTMs, and induces chromatin compaction to suppress "
        "the host's antiviral transcriptional response. As infection progresses "
        "and cell death occurs, histones are externalized into the interstitium "
        "and intravascular compartment, functioning as cytotoxic Damage-Associated "
        "Molecular Patterns (DAMPs)."
    )

    st.info(
        "**Key Reference:** Kee et al. (2022) demonstrated that deletion of either "
        "the ORF8 gene or the histone mimic site attenuates SARS-CoV-2's ability to "
        "disrupt host cell chromatin and reduces viral genome copy number.",
        icon="\U0001f4d6",
    )

    st.subheader("Viruses Analyzed")
    vcols = st.columns(4)
    virus_info = [
        ("SARS-CoV-2", "COVID-19 pandemic; primary focus of ORF8 mimicry studies"),
        ("MERS-CoV", "Highly pathogenic; used to study ISG suppression"),
        ("SARS-CoV", "2002/2003 outbreak; comparative epigenetics"),
        ("HCoV-229E", "Model for histone proteoform analysis via top-down proteomics"),
    ]
    for i, (name, desc) in enumerate(virus_info):
        with vcols[i]:
            st.metric(name, desc[:50] + "...")

    st.divider()

    st.subheader("Data Sources")
    st.caption(
        "This dashboard integrates data from multiple live databases, literature "
        "syntheses, and simulated controls. Live fetches fall back to simulated "
        "data when a service is unreachable."
    )

    _data_sources = [
        {
            "name": "Histome 2.0",
            "url": "https://www.actrec.gov.in/histome2/",
            "institution": "ACTREC, India",
            "content": "Human histone PTM sites with writer enzyme, eraser enzyme, "
                       "and modification type counts across lysine methylation, "
                       "lysine acetylation, arginine methylation, and lysine "
                       "ubiquitination categories.",
            "used_in": "Simulated fallback in the original pipeline "
                       "(histone_mark_mining.py fetch_histome2_data).",
            "reference": "Histome2: Kaur et al. (2021) Nucleic Acids Res. 49(D1):D1024-D1030.",
        },
        {
            "name": "HistoneDB 2.0",
            "url": "https://www.ncbi.nlm.nih.gov/research/histonedb/",
            "institution": "NCBI / NIH, USA",
            "content": "Histone variant sequences with HMM classification scores. "
                       "Provides mean HMM scores, sequence counts, and taxonomic "
                       "diversity per variant across sequenced genomes.",
            "used_in": "Histone variant aggregation and taxonomy-based organism "
                       "classification via NCBI E-utilities.",
            "reference": "HistoneDB 2.0: Draizen et al. (2023) Nucleic Acids Res. 51(D1):D348-D355.",
        },
        {
            "name": "ProHistoneDB",
            "url": "https://prohistonedb.org/",
            "institution": "Max Planck Institute / University of Marburg",
            "content": "Proto-histone categories across Archaea, Bacteria, and "
                       "viruses. Contains 28 phylogenetic categories including "
                       "Nucleosomal, Viral singlet/doublet, Phage histone, "
                       "Bacterial H2A-H2B, and archaeal variants.",
            "used_in": "Proto-histone category fetch with per-category sequence "
                       "counts, taxonomic breadth, and mean sequence length.",
            "reference": "ProHistoneDB: Henneman et al. (2018) Nucleic Acids Res. 46(D1):D348-D354.",
        },
        {
            "name": "Coronavirus_Vault (Literature Synthesis)",
            "url": "Local markdown files in Coronavirus_Vault/",
            "institution": "Compiled from published literature",
            "content": "Two synthesized literature reviews: (1) 'Molecular Sabotage' "
                       "detailing five mechanisms of histone-mediated tissue damage "
                       "during viral infection; (2) 'Viral Histone Mimicry' "
                       "documenting specific histone marks altered by SARS-CoV-2, "
                       "MERS-CoV, SARS-CoV, and HCoV-229E.",
            "used_in": "Histone marks table (Tab 2), mechanism cards (Tab 6), "
                       "and NETosis pathway (Tab 7).",
            "reference": "References in coronavirus26.bib (7 entries) and 5 PDFs "
                       "in Histone Mark of COVID-19/.",
        },
        {
            "name": "Simulated Nucleosome Data",
            "url": "numpy.random (local, seed=42)",
            "institution": "Built-in",
            "content": "200 synthetic nucleosomes with 5 binary histone marks "
                       "(H3K4me3, H3K27ac, H3K9me3, H3K36me3, H2AK5ac) "
                       "generated with known co-occurrence and avoidance patterns "
                       "reflecting known biology.",
            "used_in": "Occurrence, co-occurrence, and avoidance analysis "
                       "(Tabs 3, 4, 5).",
            "reference": "Designed for this dashboard.",
        },
    ]

    for ds in _data_sources:
        with st.container(border=True):
            st.markdown(f"**{ds['name']}**")
            st.markdown(f"{ds['content']}")
            st.markdown(f"- **Website:** [{ds['url']}]({ds['url']})")
            st.markdown(f"- **Used in:** {ds['used_in']}")
            st.markdown(f"- **Reference:** {ds['reference']}")

    st.divider()

    st.subheader("External APIs Queried at Runtime")
    st.caption(
        "These endpoints are accessed live during data fetches. All connections "
        "use HTTPS with configurable timeouts (10-30 s) and fall back to "
        "simulated data on failure."
    )
    _external_apis = [
        ("Histome2 PTM pages",
         "https://www.actrec.gov.in/histome2/Human/ptm.php?mod_type={mod_type}",
         "HTML scrape of PTM tables for each modification type"),
        ("HistoneDB 2.0 CSV dumps",
         "https://www.ncbi.nlm.nih.gov/research/histonedb/static/browse/dumps/seqs.txt",
         "CSV download of histone variant sequences and HMM scores"),
        ("ProHistoneDB phyloXML trees",
         "https://prohistonedb.org/static/phylotrees/{codename}.xml",
         "PhyloXML files for each proto-histone category (28 files via ThreadPoolExecutor)"),
        ("NCBI Taxonomy E-utilities",
         "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={taxids}&retmode=xml",
         "Batch taxonomy lookup for organism classification"),
    ]
    for name, url_template, desc in _external_apis:
        with st.container(border=True):
            st.markdown(f"**{name}**")
            st.markdown(f"- Endpoint: `{url_template}`")
            st.markdown(f"- {desc}")

    _pdf_dir = _HISTONE_COVID_DIR
    if _pdf_dir.exists():
        pdf_files = sorted([p.name for p in _pdf_dir.iterdir() if p.suffix == ".pdf"])
        if pdf_files:
            st.divider()
            st.subheader("Related Literature (PDFs)")
            st.caption("Full-text PDFs in the Histone Mark of COVID-19/ corpus.")
            for f in pdf_files:
                st.markdown(f"- [{f}]({_pdf_dir / f})")

# ===== Tab 2: Histone Marks in Infection ==================================
with tabs[1]:
    st.header("Histone Marks Altered During Coronavirus Infection")

    st.caption(
        "Data synthesized from Viral Histone Mimicry literature review "
        "(Coronavirus_Vault). Marks are grouped by modification type and virus."
    )

    _type_order = df_marks_filtered["type"].value_counts().index.tolist()
    fig, ax = plt.subplots(figsize=(10, max(4, len(df_marks_filtered) * 0.35)))
    if not df_marks_filtered.empty:
        type_colors = {
            "Methylation": "#3498db",
            "Acetylation": "#e74c3c",
            "Citrullination": "#2ecc71",
            "Phosphorylation": "#f39c12",
            "Glutathionylation": "#9b59b6",
            "Truncation": "#1abc9c",
        }
        y_pos = range(len(df_marks_filtered))
        bar_colors = [type_colors.get(t, "#95a5a6") for t in df_marks_filtered["type"]]
        ax.barh(y_pos, [1] * len(df_marks_filtered), color=bar_colors, alpha=0.7, edgecolor="white")
        ax.set_yticks(y_pos)
        ax.set_yticklabels(df_marks_filtered["code"], fontsize=9)
        ax.set_xlim(0, 1.5)
        ax.tick_params(left=False, labelleft=True)
        for idx, row in df_marks_filtered.iterrows():
            ax.text(1.02, idx, f"  {row['virus']}  |  {row['annotation'][:80]}...",
                    va="center", fontsize=7, color="gray")
        ax.set_title("Histone Marks Affected by Coronavirus Infection", fontsize=11, fontweight="bold")
        ax.axis("off")
    else:
        ax.text(0.5, 0.5, "No marks match the current filter.", ha="center", va="center")
    st.pyplot(fig)
    plt.close(fig)

    st.subheader("Filtered Mark Table")
    st.dataframe(
        df_marks_filtered.rename(columns={
            "type": "Modification Type",
            "code": "Mark Code",
            "virus": "Virus",
            "annotation": "Functional Annotation",
        }),
        use_container_width=True,
        column_config={
            "Modification Type": st.column_config.TextColumn(width="small"),
            "Mark Code": st.column_config.TextColumn(width="medium"),
            "Virus": st.column_config.TextColumn(width="small"),
            "Functional Annotation": st.column_config.TextColumn(width="large"),
        },
    )

    st.subheader("Additional PTMs Identified via Proteomics")
    st.caption("Broader set of modifications from top-down proteomics studies (mined from literature).")
    cols_additional = st.columns(3)
    for i, ptm in enumerate(ADDITIONAL_PTMS):
        cols_additional[i % 3].markdown(f"- {ptm}")

# ===== Tab 3: Occurrence ==================================================
with tabs[2]:
    _cooc_meta = COOC_SOURCE_META[cooc_source]

    st.header("Histone Mark Occurrence")
    st.caption(
        "Occurrence measures the fraction of regions carrying each histone mark "
        "in the current dataset."
    )

    st.info(f"**Data source:** {cooc_source_label}", icon="\U0001f4ca")

    nuc_cols = st.columns(4)
    nuc_cols[0].metric("Regions", len(df_nuc))
    nuc_cols[1].metric("Histone marks", len(marks_nuc))
    nuc_cols[2].metric("Encoding", "Binary (present / absent)")
    nuc_cols[3].metric("Significance (FDR)", f"{cooc_alpha}")

    with st.expander("Data preview (first 10 rows)"):
        st.caption(_cooc_meta["citation"])
        st.dataframe(df_nuc.head(10), use_container_width=True)

    st.subheader("Occurrence Frequency")
    occ = df_nuc[marks_nuc].mean()
    _draw_occurrence_barchart(occ)

    st.subheader("Mark Descriptions")
    for mark, desc in _cooc_meta["descriptions"].items():
        freq = occ.get(mark, 0)
        st.markdown(f"- **{mark}**: {desc} (occurrence = {freq * 100:.1f}%)")
    st.caption(f"Binarization method: {_cooc_meta['binary_method']}")

# ===== Tab 4: Co-occurrence ==============================================
with tabs[3]:
    _cooc_meta = COOC_SOURCE_META[cooc_source]

    st.header("Histone Mark Co-occurrence")
    st.caption(
        "Co-occurrence measures whether two histone marks appear on the same "
        "region more often than expected under independence. "
        "The null expectation is modeled as Poisson(lambda = n * P(A) * P(B)). "
        "P-values are FDR-adjusted (Benjamini-Hochberg)."
    )

    st.info(f"**Data source:** {cooc_source_label}", icon="\U0001f4ca")
    st.caption(f"Binarization: {_cooc_meta['binary_method']}")

    st.subheader("Co-occurrence Matrices")
    _draw_cooccurrence_heatmap(cooc_stats)

    _n_sig_co = df_cooc_test["cooccurrence_significant"].sum()
    if _n_sig_co > 0:
        st.success(
            f"Found {_n_sig_co} significant co-occurrence relationships "
            f"(FDR-adjusted p < {cooc_alpha})."
        )
    else:
        st.info("No significant co-occurrence at the chosen threshold.")

    st.subheader("Significant Co-occurrence Pairs")
    _co_pairs = df_cooc_test[df_cooc_test["cooccurrence_significant"]]
    if not _co_pairs.empty:
        st.dataframe(
            _co_pairs[["pair", "observed", "expected", "log2_fold_change",
                       "p_cooccurrence_adj"]]
            .sort_values("p_cooccurrence_adj")
            .rename(columns={
                "pair": "Mark pair",
                "observed": "Observed co-occurrence",
                "expected": "Expected (independence)",
                "log2_fold_change": "Log2(obs/exp)",
                "p_cooccurrence_adj": "p-adj",
            }),
            use_container_width=True,
        )
    else:
        st.info("No pairs pass the significance threshold.")

    st.subheader("Co-occurrence Network")
    _draw_cooccurrence_network(cooc_stats)

    st.divider()

    st.subheader("Full Pairwise Results (all pairs)")
    st.dataframe(
        df_cooc_test.sort_values("p_cooccurrence_adj").rename(columns={
            "pair": "Mark pair",
            "observed": "Observed",
            "expected": "Expected",
            "log2_fold_change": "Log2(obs/exp)",
            "p_cooccurrence_adj": "p-adj (co-occurrence)",
            "p_avoidance_adj": "p-adj (avoidance)",
            "cooccurrence_significant": "Co-occurrence",
            "avoidance_significant": "Avoidance",
        }),
        use_container_width=True,
    )

    st.divider()

    if cooc_source == "Simulated":
        st.subheader("Mark Type Co-occurrence in Viral Literature")
        st.caption(
            "Tests whether pairs of histone mark modification types (methylation, "
            "acetylation, etc.) co-occur across viral literature entries more than "
            "expected by chance, using the same Poisson framework."
        )
        _draw_virus_mark_cooccurrence(df_marks)
        st.info(
            "**Interpretation:** A significant result suggests that certain mark "
            "types are preferentially co-regulated by coronaviruses, potentially "
            "reflecting coordinated epigenetic reprogramming strategies.",
            icon="\U0001f9ec",
        )

# ===== Tab 5: Avoidance ==================================================
with tabs[4]:
    _cooc_meta = COOC_SOURCE_META[cooc_source]

    st.header("Histone Mark Avoidance")
    st.caption(
        "Avoidance (mutual exclusion) occurs when two histone marks appear "
        "together less often than expected by chance. "
        "Significance is assessed via the left-tailed Poisson test "
        "P(observed <= k | lambda) with FDR correction."
    )

    st.info(f"**Data source:** {cooc_source_label}", icon="\U0001f4ca")
    st.caption(f"Binarization: {_cooc_meta['binary_method']}")

    st.subheader("Mark Descriptions")
    for mark, desc in _cooc_meta["descriptions"].items():
        st.markdown(f"- **{mark}**: {desc}")

    _n_sig_av = df_cooc_test["avoidance_significant"].sum()
    if _n_sig_av > 0:
        st.success(
            f"Found {_n_sig_av} significant avoidance relationships "
            f"(FDR-adjusted p < {cooc_alpha})."
        )
    else:
        st.info("No significant avoidance at the chosen threshold.")

    st.subheader("Significant Avoidance Pairs")
    _av_pairs = df_cooc_test[df_cooc_test["avoidance_significant"]]
    if not _av_pairs.empty:
        st.dataframe(
            _av_pairs[["pair", "observed", "expected", "log2_fold_change",
                       "p_avoidance_adj"]]
            .sort_values("p_avoidance_adj")
            .rename(columns={
                "pair": "Mark pair",
                "observed": "Observed co-occurrence",
                "expected": "Expected (independence)",
                "log2_fold_change": "Log2(obs/exp)",
                "p_avoidance_adj": "p-adj (avoidance)",
            }),
            use_container_width=True,
        )
    else:
        st.info("No pairs pass the significance threshold for avoidance.")

    st.subheader("Avoidance Network")
    _draw_cooccurrence_network(cooc_stats)
    st.caption(
        "Red dashed edges highlight significant avoidance (marks that rarely "
        "co-occur). Green edges show co-occurrence for comparison."
    )

    st.subheader("Summary Table")
    st.dataframe(
        df_cooc_test[["pair", "observed", "expected", "log2_fold_change",
                       "p_avoidance_adj", "avoidance_significant"]]
        .sort_values("p_avoidance_adj")
        .rename(columns={
            "pair": "Mark pair",
            "observed": "Observed co-occurrence",
            "expected": "Expected (independence)",
            "log2_fold_change": "Log2(obs/exp)",
            "p_avoidance_adj": "p-adj (avoidance)",
            "avoidance_significant": "Significant avoidance",
        }),
        use_container_width=True,
    )

# ===== Tab 6: Molecular Mechanisms =======================================
with tabs[5]:
    st.header("Five Mechanisms of Histone-Mediated Tissue Damage")

    if view_mode == "Mechanism Explorer":
        selected_ids = st.multiselect(
            "Select mechanisms to display",
            options=[m["id"] for m in MECHANISMS],
            format_func=lambda x: f"M{x}: {MECHANISMS[x-1]['title'][:50]}...",
            default=[m["id"] for m in MECHANISMS],
        )
        for m in MECHANISMS:
            if m["id"] in selected_ids:
                _render_mechanism_card(m)
    else:
        for m in MECHANISMS:
            _render_mechanism_card(m)

# ===== Tab 7: NETosis & Pathways =========================================
with tabs[6]:
    st.header("NETosis and the DNA-Shielding Paradox")

    st.markdown(
        "Neutrophil Extracellular Traps (NETs) are a double-edged sword in "
        "viral infection. They are released as a defense mechanism, but "
        "extracellular histones within NETs drive vascular damage, ARDS, "
        "and predict 28-day mortality in COVID-19 patients."
    )

    _draw_netosis_pathway()

    st.info(
        "**The Shielding Paradox:** DNA forms the scaffold of the trap, but also "
        "acts as a polyanionic shield that temporarily neutralizes histone toxicity. "
        "Digestion of this DNA by host nucleases (e.g., DNase I) 'unmasks' the "
        "histones, restoring their potent cytotoxicity toward the endothelium.",
        icon="\U0001f9ea",
    )

    st.divider()

    st.subheader("Apoptosis vs. Pyroptosis")
    apop_cols = st.columns([1, 2, 2])
    apop_cols[0].markdown("**Feature**")
    apop_cols[1].markdown("**Apoptosis**")
    apop_cols[2].markdown("**Pyroptosis**")
    apop_rows = [
        ("Primary Trigger", "Mitochondrial stress / Caspase 3, 7",
         "Inflammasome / Caspase 1"),
        ("Morphology", "Cell shrinkage & membrane blebbing",
         "Cell swelling & violent cytolysis"),
        ("Inflammatory Impact", "Minimal (ordered phagocytosis)",
         "Massive (secretion of IL-1beta and IL-18)"),
    ]
    for label, apop, pyro in apop_rows:
        c1, c2, c3 = st.columns([1, 2, 2])
        c1.markdown(f"**{label}**")
        c2.markdown(apop)
        c3.markdown(pyro)

    st.divider()

    st.subheader("ORF8 Histone Mimicry: Epigenetic Subversion Network")
    st.caption(
        "SARS-CoV-2 ORF8 mimics the ARKS motif of Histone H3, associating with "
        "chromatin to disrupt PTM regulation and promote gene silencing."
    )
    _draw_mimicry_network()

    st.divider()

    st.subheader("Immunothrombosis Vicious Circle")
    st.caption(
        "Histone-mediated complement inhibition creates a self-reinforcing cycle "
        "of impaired viral clearance, continued cell death, and escalating histone release."
    )
    _draw_vicious_cycle()

# ===== Tab 8: Clinical Countermeasures ===================================
with tabs[7]:
    st.header("Molecular Shields: Therapeutic Countermeasures")

    st.markdown(
        "Targeting extracellular histones offers a 'host-directed therapy' "
        "approach that could bridge the gap between antiviral treatments and "
        "management of hyperinflammatory sequelae."
    )

    shields = [
        {
            "name": "C-Reactive Protein (CRP)",
            "mechanism": "Binds directly to Histone H4, neutralizing "
                         "membrane-disruptive toxicity.",
            "evidence": "Natural pentraxin; acute-phase reactant.",
            "status": "Endogenous; therapeutic potential being explored.",
        },
        {
            "name": "Activated Protein C (APC)",
            "mechanism": "Proteolytic cleavage of histones, destroying "
                         "their cytotoxic activity.",
            "evidence": "Natural negative feedback loop in coagulation.",
            "status": "Recombinant APC studied in sepsis (drotrecogin alfa).",
        },
        {
            "name": "Low-Molecular-Weight Heparin (LMWH)",
            "mechanism": "Highly negative charge neutralizes histones via "
                         "electrostatic binding; protects endothelial "
                         "glycocalyx from heparan sulfate degradation.",
            "evidence": "Standard of care for COVID-19 thromboprophylaxis.",
            "status": "Approved; used clinically in COVID-19.",
        },
        {
            "name": "Designed Cyclic Peptides (e.g., HIPe)",
            "mechanism": "Specifically block N-terminal tails of H4 or H2A, "
                         "preventing interaction with lipid bilayer.",
            "evidence": "Rational design; in vitro validation.",
            "status": "Preclinical / experimental.",
        },
    ]

    cols = st.columns(2)
    for i, s in enumerate(shields):
        with cols[i % 2]:
            with st.container(border=True):
                st.subheader(s["name"])
                st.markdown(f"**Mechanism:** {s['mechanism']}")
                st.markdown(f"**Evidence:** {s['evidence']}")
                st.markdown(f"**Status:** {s['status']}")

# ===== Tab 9: References =================================================
with tabs[8]:
    st.header("Bibliography")

    st.caption(
        "References from the Coronavirus_Vault literature corpus. "
        "All DOIs link to the published articles."
    )

    ref_group = st.radio(
        "Group by",
        options=["Chronological", "By Author", "Flat list"],
        index=0,
        horizontal=True,
    )

    refs_sorted = list(REFERENCES)
    if ref_group == "Chronological":
        refs_sorted.sort(key=lambda x: x["year"], reverse=True)
    elif ref_group == "By Author":
        refs_sorted.sort(key=lambda x: x["authors"].split(",")[0].strip())

    for ref in refs_sorted:
        with st.container(border=True):
            cols = st.columns([0.8, 0.2])
            with cols[0]:
                st.markdown(
                    f"**{ref['authors']}** ({ref['year']}). "
                    f"*{ref['title']}*. "
                    f"{ref['journal']} **{ref['volume']}**: {ref['pages']}. "
                    f"DOI: [{ref['doi']}](https://doi.org/{ref['doi']})"
                )
            with cols[1]:
                st.markdown(f"PMID: {ref['pmid']}" if ref.get("pmid") else "", unsafe_allow_html=True)
                if ref.get("url"):
                    st.markdown(f"[Open]({ref['url']})")

    st.divider()
    st.caption(
        "Additional references embedded in the Coronavirus_Vault markdown files "
        "cite literature from the Histone Mark of COVID-19/ PDF corpus "
        "(full_text.pdf, full_text1.pdf, full_text2.pdf, full_textxx.pdf, "
        "human_coronavirus_229e_infection_alters_histone_proteoforms.pdf)."
    )
