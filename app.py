#!/usr/bin/env python3
"""
Module: app
Purpose: Streamlit interactive dashboard for histone mark network analysis
Author: Dr. Arli Aditya Parikesit
Date: 2026
Parameters:
    - threshold: cosine similarity threshold for graph edges (sidebar slider)
    - n_clusters: number of KMeans chromatin state clusters (sidebar slider)
    - db_choice: data source selection (sidebar radio: Histome2 / HistoneDB 2.0 / ProHistoneDB / Simulated)
References:
    - Clauset A, Shalizi CR, Newman MEJ (2009) Power-law distributions in
      empirical data. SIAM Rev 51:661-703. https://doi.org/10.1137/070710111
    - Hagberg AA, Schult DA, Swart PJ (2008) Exploring network structure,
      dynamics, and function using NetworkX. Proc SciPy 2008, pp 11-15.
"""
import os
import sys
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import streamlit as st
from scipy.optimize import curve_fit
from scipy.stats import ks_2samp, poisson
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import davies_bouldin_score, silhouette_score
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(__file__))
from histone_mark_mining import (
    data_mining,
    fetch_histome2_data,
    fetch_histonedb_data,
    fetch_prohistonedb_data,
)

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="Histone Mark Network Analysis",
    page_icon=None,
    layout="wide",
)

# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------

st.sidebar.title("Parameters")
db_choice = st.sidebar.radio(
    "Data source",
    options=["Histome2", "HistoneDB 2.0", "ProHistoneDB", "Simulated"],
    index=0,
    help=(
        "Histome2: human histone PTM sites (writers/erasers). "
        "HistoneDB 2.0: histone variant sequences with HMM scores (NCBI). "
        "ProHistoneDB: proto-histone categories across Archaea, Bacteria, viruses. "
        "Simulated: synthetic data for testing (seed=42)."
    ),
)
threshold = st.sidebar.slider("Cosine similarity threshold", 0.50, 1.00, 0.95, 0.01)
n_clusters = st.sidebar.slider("Chromatin clusters (k)", 2, 5, 3)

st.sidebar.markdown("---")
st.sidebar.caption(
    "Dr. Arli Aditya Parikesit  \n"
    "i3L University Jakarta  \n"
    "[ORCID 0000-0001-8716-3926](https://orcid.org/0000-0001-8716-3926)"
)

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

_MARK_DESCRIPTIONS = {
    "Histome2": {
        "mark1": "n_writer_enzymes — number of enzymes that add this PTM",
        "mark2": "n_eraser_enzymes — number of enzymes that remove this PTM",
        "mark3": "n_modification_types — how many PTM classes this site belongs to",
    },
    "HistoneDB 2.0": {
        "mark1": "mean_HMM_score — mean best-model HMM score across sequences of this variant",
        "mark2": "n_sequences — total sequences classified as this variant",
        "mark3": "n_taxa — number of distinct taxonomic IDs carrying this variant",
    },
    "ProHistoneDB": {
        "mark1": "n_sequences — number of proto-histone sequences in this category",
        "mark2": "n_taxa — number of distinct taxonomic groups in this category",
        "mark3": "mean_seq_len_aa — mean amino acid sequence length",
    },
    "Simulated": {
        "mark1": "H3K4me3 signal (active promoter mark, simulated)",
        "mark2": "H3K27ac signal (active enhancer mark, correlated with mark1, simulated)",
        "mark3": "H3K9me3 signal (repressive heterochromatin mark, simulated)",
    },
}

@st.cache_data(show_spinner="Loading data...")
def load_data(db_choice: str):
    if db_choice == "Histome2":
        df = fetch_histome2_data()
        if df is not None and not df.empty:
            if "count" not in df.columns:
                df["count"] = np.random.default_rng(42).poisson(lam=25, size=len(df))
            return df, "Histome2 (live)"
        return pd.DataFrame(data_mining()), "Simulated (Histome2 unavailable)"
    if db_choice == "HistoneDB 2.0":
        df = fetch_histonedb_data()
        if df is not None and not df.empty:
            return df, "HistoneDB 2.0 (live)"
        return pd.DataFrame(data_mining()), "Simulated (HistoneDB unavailable)"
    if db_choice == "ProHistoneDB":
        df = fetch_prohistonedb_data()
        if df is not None and not df.empty:
            return df, "ProHistoneDB (live)"
        return pd.DataFrame(data_mining()), "Simulated (ProHistoneDB unavailable)"
    return pd.DataFrame(data_mining()), "Simulated (seed=42)"

df_raw, source_label = load_data(db_choice)
st.sidebar.info(f"Source: **{source_label}**")
if "unavailable" in source_label:
    st.sidebar.warning(f"{db_choice} fetch failed — using simulated fallback.")

numeric_cols = df_raw.select_dtypes(include=[np.number]).columns
histone_marks = [c for c in numeric_cols if c.startswith("mark")]

if not histone_marks:
    st.error("No histone mark columns detected (expected columns prefixed with 'mark'). Check data source.")
    st.stop()

df = df_raw.copy()
scaler = StandardScaler()
df[histone_marks] = scaler.fit_transform(df[histone_marks])

# ---------------------------------------------------------------------------
# Graph construction (cached on data + threshold)
# ---------------------------------------------------------------------------

@st.cache_data(show_spinner="Building similarity graph...")
def build_graph(mark_values: np.ndarray, regions: list, sim_threshold: float):
    sim = cosine_similarity(mark_values)
    G = nx.DiGraph()
    for r in regions:
        G.add_node(r)
    for i in range(len(regions)):
        for j in range(len(regions)):
            if i != j and sim[i, j] > sim_threshold:
                G.add_edge(regions[i], regions[j])
    return G, sim

regions = list(df["region"])
G, sim_matrix = build_graph(df[histone_marks].values, regions, threshold)

deg_c = nx.degree_centrality(G)
between_c = nx.betweenness_centrality(G)
close_c = nx.closeness_centrality(G)

# ---------------------------------------------------------------------------
# PCA + KMeans (computed once; results used across tabs)
# ---------------------------------------------------------------------------

features = df[histone_marks].values
n_s, n_f = features.shape
n_comps = min(2, n_s, n_f)
pca_model = PCA(n_components=n_comps)
reduced = pca_model.fit_transform(features) if n_comps > 0 else features

k = min(n_clusters, n_s)
kmeans = KMeans(n_clusters=k, random_state=42, n_init="auto")
cluster_labels = kmeans.fit_predict(reduced)
df["cluster"] = cluster_labels

# ---------------------------------------------------------------------------
# Page header
# ---------------------------------------------------------------------------

st.title("Histone Mark Co-occurrence Network Analysis")
st.caption(
    "Interactive dashboard for epigenetic network construction, scale-free topology "
    "testing, chromatin state clustering, and statistical evaluation."
)

# ---------------------------------------------------------------------------
# Tabs
# ---------------------------------------------------------------------------

tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
    "Data Overview",
    "Network Graph",
    "Correlation",
    "Scale-Free Analysis",
    "Clustering",
    "Statistics",
])

# ---- Tab 1: Data Overview --------------------------------------------------
with tab1:
    c1, c2, c3 = st.columns(3)
    c1.metric("Regions", len(df))
    c2.metric("Histone marks", len(histone_marks))
    c3.metric("Data source", source_label)

    st.subheader("Raw data (first 10 rows)")
    display_cols = (
        ["region", "histone_modification"] + histone_marks + ["count"]
        if "histone_modification" in df_raw.columns
        else list(df_raw.columns)
    )
    st.dataframe(
        df_raw[display_cols].head(10).rename(
            columns={"histone_modification": "Histone modification (dominant mark)"}
        ),
        use_container_width=True,
    )

    if "histone_modification" in df_raw.columns:
        st.caption(
            "**Histone modification column**: dominant mark per region by raw signal intensity. "
            "H3K4me3 = trimethyl-Lys4 H3 (active promoter); "
            "H3K27ac = acetyl-Lys27 H3 (active enhancer/promoter); "
            "H3K9me3 = trimethyl-Lys9 H3 (repressive heterochromatin)."
        )
        mod_counts = df_raw["histone_modification"].value_counts().reset_index()
        mod_counts.columns = ["Modification", "Region count"]
        st.dataframe(mod_counts, use_container_width=True)

    _db_key = db_choice if db_choice in _MARK_DESCRIPTIONS else "Simulated"
    _mdesc = _MARK_DESCRIPTIONS[_db_key]
    st.subheader("Mark column meanings for this source")
    for col, desc in _mdesc.items():
        if col in histone_marks:
            st.markdown(f"- **{col}**: {desc}")

    st.subheader("Descriptive statistics (raw)")
    st.dataframe(df_raw[histone_marks].describe().round(4), use_container_width=True)

# ---- Tab 2: Network Graph --------------------------------------------------
with tab2:
    c1, c2, c3 = st.columns(3)
    c1.metric("Nodes", G.number_of_nodes())
    c2.metric("Edges", G.number_of_edges())
    c3.metric("Threshold", threshold)

    fig, ax = plt.subplots(figsize=(10, 8))
    pos = nx.spring_layout(G, k=0.15, seed=42)
    node_degrees = [G.degree(n) for n in G.nodes()]
    max_deg = max(node_degrees) if node_degrees else 1
    hub_threshold = np.percentile(node_degrees, 80) if node_degrees else 0
    colors = ["red" if G.degree(n) > hub_threshold else "lightblue" for n in G.nodes()]
    nx.draw(G, pos, ax=ax, with_labels=False, node_color=colors,
            edge_color="gray", node_size=100, alpha=0.8, width=0.5)
    ax.set_title("Epigenetic Network (red = hub nodes, degree > 80th percentile)")
    st.pyplot(fig)
    plt.close(fig)

    st.subheader("Centrality — top 10")
    _mod_map = (
        df_raw.set_index("region")["histone_modification"].to_dict()
        if "histone_modification" in df_raw.columns
        else {}
    )
    cent_df = pd.DataFrame({
        "Region": list(deg_c.keys()),
        "Histone modification": [_mod_map.get(r, "") for r in deg_c],
        "Degree centrality": [round(deg_c[r], 4) for r in deg_c],
        "Betweenness centrality": [round(between_c[r], 4) for r in deg_c],
        "Closeness centrality": [round(close_c[r], 4) for r in deg_c],
    }).sort_values("Degree centrality", ascending=False).head(10).reset_index(drop=True)
    st.dataframe(cent_df, use_container_width=True)

# ---- Tab 3: Correlation ----------------------------------------------------
with tab3:
    st.subheader("Pearson correlation matrix")
    corr = df[histone_marks].corr(method="pearson")

    fig, ax = plt.subplots(figsize=(max(5, len(histone_marks) + 1), max(4, len(histone_marks))))
    im = ax.imshow(corr.values, cmap="coolwarm", vmin=-1, vmax=1)
    plt.colorbar(im, ax=ax)
    ax.set_xticks(range(len(histone_marks)))
    ax.set_yticks(range(len(histone_marks)))
    ax.set_xticklabels(histone_marks, rotation=45, ha="right")
    ax.set_yticklabels(histone_marks)
    for i in range(len(histone_marks)):
        for j in range(len(histone_marks)):
            ax.text(j, i, f"{corr.values[i, j]:.2f}", ha="center", va="center", fontsize=11)
    ax.set_title("Pearson Correlation of Standardized Histone Marks")
    st.pyplot(fig)
    plt.close(fig)

    st.dataframe(corr.round(4), use_container_width=True)

# ---- Tab 4: Scale-Free Analysis --------------------------------------------
with tab4:
    avg_clust = nx.average_clustering(G)
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    p_rand = n_edges / (n_nodes * (n_nodes - 1)) if n_nodes > 1 else 0
    p_rand = min(max(p_rand, 0.0), 1.0)
    rand_G = nx.fast_gnp_random_graph(n_nodes, p_rand, directed=True, seed=42)
    avg_clust_rand = nx.average_clustering(rand_G)

    c1, c2 = st.columns(2)
    c1.metric("Clustering coeff (network)", f"{avg_clust:.4f}")
    c2.metric("Clustering coeff (Erdos-Renyi)", f"{avg_clust_rand:.4f}")

    degree_values = [d for _, d in G.degree()]
    degree_counts = Counter(degree_values)
    x_deg = sorted([d for d in degree_counts if d > 0])
    y_freq = [degree_counts[d] for d in x_deg]

    gamma_fit = gamma_se = r2 = None
    _fit_failed = False
    if len(x_deg) >= 3:
        log_x = np.log(np.array(x_deg, dtype=float))
        log_y = np.log(np.array(y_freq, dtype=float))
        def _pl(lk, lc, g): return lc - g * lk
        try:
            popt, pcov = curve_fit(_pl, log_x, log_y)
            gamma_fit = popt[1]
            _gamma_se_raw = np.sqrt(pcov[1, 1])
            gamma_se = _gamma_se_raw if np.isfinite(_gamma_se_raw) else None
            log_y_pred = _pl(log_x, *popt)
            ss_res = np.sum((log_y - log_y_pred) ** 2)
            ss_tot = np.sum((log_y - np.mean(log_y)) ** 2)
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        except Exception:
            _fit_failed = True

    fig, ax = plt.subplots(figsize=(8, 5))
    if x_deg:
        ax.loglog(x_deg, y_freq, "bo", label="Observed", markersize=6)
        if gamma_fit is not None:
            x_fit = np.linspace(min(x_deg), max(x_deg), 200)
            y_fit = np.exp(popt[0]) * x_fit ** (-popt[1])
            se_label = f" ± {gamma_se:.2f}" if gamma_se is not None else " (SE undefined)"
            ax.loglog(x_fit, y_fit, "r-",
                      label=f"Power-law fit  γ = {gamma_fit:.2f}{se_label}")
    ax.set_xlabel("Degree (log scale)")
    ax.set_ylabel("Frequency (log scale)")
    ax.set_title("Degree Distribution (Log-Log)")
    ax.legend()
    st.pyplot(fig)
    plt.close(fig)

    if _fit_failed:
        st.caption("Power-law fit did not converge.")

    st.subheader("KS test vs. Poisson")
    mean_deg = np.mean(degree_values) if degree_values else 0
    _rng = np.random.default_rng(42)
    pois_s = _rng.poisson(lam=mean_deg, size=len(degree_values)) if mean_deg > 0 else np.zeros(len(degree_values))
    ks_stat, ks_p = ks_2samp(degree_values, pois_s)
    c1, c2 = st.columns(2)
    c1.metric("KS statistic", f"{ks_stat:.4f}")
    c2.metric("p-value", f"{ks_p:.4f}")
    if ks_p < 0.05:
        st.success("Degree distribution differs significantly from Poisson — supports non-random topology.")
    else:
        st.warning("Degree distribution does not differ significantly from Poisson.")

    if gamma_fit is not None:
        st.subheader("Power-law fit (OLS)")
        c1, c2, c3 = st.columns(3)
        c1.metric("gamma", f"{gamma_fit:.3f} ± {gamma_se:.3f}")
        c2.metric("R² (log-log)", f"{r2:.4f}")
        if 2.0 <= gamma_fit <= 3.0:
            c3.success("gamma in [2, 3]: consistent with scale-free.")
        else:
            c3.warning(f"gamma = {gamma_fit:.3f} outside [2, 3].")
        st.caption(
            "OLS estimate only. Use `powerlaw` package (Clauset et al. 2009) for publication-grade claims."
        )

# ---- Tab 5: Clustering -----------------------------------------------------
with tab5:
    c1, c2, c3 = st.columns(3)
    c1.metric("PCA components", n_comps)
    c2.metric("Clusters (k)", k)
    if len(set(cluster_labels)) > 1:
        sil = silhouette_score(reduced, cluster_labels)
        db = davies_bouldin_score(reduced, cluster_labels)
        c3.metric("Silhouette score", f"{sil:.4f}")
        st.metric("Davies-Bouldin index", f"{db:.4f}",
                  help="Lower = more compact, better-separated clusters")
    else:
        st.warning("Only one cluster found — adjust k or threshold.")

    fig, ax = plt.subplots(figsize=(8, 6))
    for c in sorted(set(cluster_labels)):
        mask = cluster_labels == c
        x_vals = reduced[mask, 0]
        y_vals = reduced[mask, 1] if n_comps > 1 else np.zeros(mask.sum())
        ax.scatter(x_vals, y_vals, label=f"Cluster {c}", s=60, alpha=0.8)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2" if n_comps > 1 else "")
    ax.set_title(f"KMeans Clusters (k={k}) on PCA-Reduced Histone Marks")
    ax.legend()
    st.pyplot(fig)
    plt.close(fig)

    st.subheader("Cluster assignments (first 20 regions)")
    _cluster_cols = (
        ["region", "histone_modification"] + histone_marks + ["cluster"]
        if "histone_modification" in df.columns
        else ["region"] + histone_marks + ["cluster"]
    )
    st.dataframe(df[_cluster_cols].head(20), use_container_width=True)

# ---- Tab 6: Statistics -----------------------------------------------------
with tab6:
    st.subheader("Poisson cluster significance")
    st.caption(
        "Tests whether cluster size exceeds expectation under "
        "Poisson(lambda = mean read count for that cluster)."
    )
    rows = []
    cluster_counts = df["cluster"].value_counts()
    for c in sorted(cluster_counts.index):
        count = int(cluster_counts[c])
        mean_val = float(df.loc[df["cluster"] == c, "count"].mean()) if "count" in df.columns else 0.0
        p_val = float(poisson.sf(count, mean_val)) if mean_val > 0 else 1.0
        rows.append({
            "Cluster": c,
            "n regions": count,
            "Mean read count (lambda)": round(mean_val, 2),
            "Poisson p-value": round(p_val, 4),
            "Significant (p<0.05)": p_val < 0.05,
        })
    st.dataframe(pd.DataFrame(rows), use_container_width=True)
