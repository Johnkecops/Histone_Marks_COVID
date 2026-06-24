#!/usr/bin/env python3
"""
Benchmark three databases: Histome2, HistoneDB 2.0, ProHistoneDB
Compute occurrence, co-occurrence, avoidance for each.
Author: Dr. Arli Aditya Parikesit
Date: 2026
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import poisson
from statsmodels.stats.multitest import multipletests
from matplotlib.lines import Line2D
from histone_mark_mining import fetch_histome2_data, fetch_histonedb_data, fetch_prohistonedb_data

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "figures")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---- Fetch all three databases ----
print("Fetching Histome2 data...")
histome2_raw = fetch_histome2_data()

print("Fetching HistoneDB 2.0 data...")
histonedb_raw = fetch_histonedb_data()

print("Fetching ProHistoneDB data...")
prohistonedb_raw = fetch_prohistonedb_data()

# ---- Binarization and mark naming ----
DATABASES = {}

# Histome2
if histome2_raw is not None and not histome2_raw.empty:
    df = histome2_raw.copy()
    marks_map = {"mark1": "has_writer", "mark2": "has_eraser", "mark3": "multi_mod"}
    binary_marks = []
    for col_key, new_name in marks_map.items():
        vals = pd.to_numeric(df[col_key], errors="coerce").dropna()
        if len(vals) >= 2:
            median = vals.median()
            df[new_name] = (df[col_key].astype(float) > median).astype(int)
            binary_marks.append(new_name)
    DATABASES["Histome2"] = {
        "df": df,
        "marks": binary_marks,
        "n": len(df),
        "descriptions": {
            "has_writer": "Above-median writer enzyme count at PTM site",
            "has_eraser": "Above-median eraser enzyme count at PTM site",
            "multi_mod": "Above-median number of modification types at site",
        },
        "source_info": f"{len(df)} PTM sites from human histone database, writer/eraser/mod-type counts binarized by median split",
    }
    print(f"  Histome2: {len(df)} sites, marks={binary_marks}")

# HistoneDB 2.0
if histonedb_raw is not None and not histonedb_raw.empty:
    df = histonedb_raw.copy()
    marks_map = {"mark1": "high_HMM_score", "mark2": "many_sequences", "mark3": "broad_taxonomy"}
    binary_marks = []
    for col_key, new_name in marks_map.items():
        vals = pd.to_numeric(df[col_key], errors="coerce").dropna()
        if len(vals) >= 2:
            median = vals.median()
            df[new_name] = (df[col_key].astype(float) > median).astype(int)
            binary_marks.append(new_name)
    DATABASES["HistoneDB 2.0"] = {
        "df": df,
        "marks": binary_marks,
        "n": len(df),
        "descriptions": {
            "high_HMM_score": "Above-median HMM classification confidence",
            "many_sequences": "Above-median number of sequences per variant",
            "broad_taxonomy": "Above-median number of taxonomic IDs",
        },
        "source_info": f"{len(df)} histone variants from NCBI with HMM scores, sequence counts, and taxonomic breadth",
    }
    print(f"  HistoneDB: {len(df)} variants, marks={binary_marks}")

# ProHistoneDB
if prohistonedb_raw is not None and not prohistonedb_raw.empty:
    df = prohistonedb_raw.copy()
    marks_map = {"mark1": "many_sequences", "mark2": "broad_taxonomy", "mark3": "long_sequence"}
    binary_marks = []
    for col_key, new_name in marks_map.items():
        vals = pd.to_numeric(df[col_key], errors="coerce").dropna()
        if len(vals) >= 2:
            median = vals.median()
            df[new_name] = (df[col_key].astype(float) > median).astype(int)
            binary_marks.append(new_name)
    DATABASES["ProHistoneDB"] = {
        "df": df,
        "marks": binary_marks,
        "n": len(df),
        "descriptions": {
            "many_sequences": "Above-median sequence count per category",
            "broad_taxonomy": "Above-median taxonomic group count",
            "long_sequence": "Above-median mean sequence length (aa)",
        },
        "source_info": f"{len(df)} proto-histone categories across Archaea, Bacteria, and viruses",
    }
    print(f"  ProHistoneDB: {len(df)} categories, marks={binary_marks}")

# ---- Simulated nucleosome data (added alongside real databases) ----
NUCLEOSOME_MARKS = ["H3K4me3", "H3K27ac", "H3K9me3", "H3K36me3", "H2AK5ac"]

def generate_nucleosome_data(n_nuc=200, seed=42):
    rng = np.random.default_rng(seed)
    regions = [f"nuc_{i}" for i in range(n_nuc)]
    h3k4me3 = rng.binomial(1, 0.35, n_nuc)
    h3k27ac = np.array([
        1 if (h3k4me3[i] == 1 and rng.random() < 0.75)
          or (h3k4me3[i] == 0 and rng.random() < 0.15)
        else 0 for i in range(n_nuc)
    ])
    active = (h3k4me3 | h3k27ac).astype(bool)
    h3k9me3 = np.array([
        1 if (not active[i] and rng.random() < 0.40)
          or (active[i] and rng.random() < 0.06)
        else 0 for i in range(n_nuc)
    ])
    h3k36me3 = rng.binomial(1, 0.25, n_nuc)
    h2ak5ac = np.array([
        1 if (h3k27ac[i] == 1 and rng.random() < 0.45)
          or (h3k27ac[i] == 0 and rng.random() < 0.10)
        else 0 for i in range(n_nuc)
    ])
    return pd.DataFrame({
        "region": regions,
        "H3K4me3": h3k4me3.astype(int),
        "H3K27ac": h3k27ac.astype(int),
        "H3K9me3": h3k9me3.astype(int),
        "H3K36me3": h3k36me3.astype(int),
        "H2AK5ac": h2ak5ac.astype(int),
    })

df_sim = generate_nucleosome_data(n_nuc=200, seed=42)
sim_marks = NUCLEOSOME_MARKS
DATABASES["Simulated"] = {
    "df": df_sim,
    "marks": sim_marks,
    "n": len(df_sim),
    "descriptions": {
        "H3K4me3": "Active promoter (trimethylation of H3 Lys4)",
        "H3K27ac": "Active enhancer/promoter (acetylation of H3 Lys27)",
        "H3K9me3": "Repressive heterochromatin (trimethylation of H3 Lys9)",
        "H3K36me3": "Gene body elongation (trimethylation of H3 Lys36)",
        "H2AK5ac": "Minor acetylation (acetylation of H2A Lys5)",
    },
    "source_info": f"{len(df_sim)} synthetic nucleosomes with 5 binary marks encoding known co-occurrence and avoidance rules",
}
print(f"  Simulated: {len(df_sim)} nucleosomes, marks={sim_marks}")

# ---- Compute co-occurrence and avoidance for each ----
def compute_stats(df_binary, marks, alpha=0.05):
    n = len(df_binary)
    m = len(marks)
    occurrence = df_binary[marks].mean()
    obs_mat = np.zeros((m, m))
    exp_mat = np.zeros((m, m))
    p_co = np.ones((m, m))
    p_av = np.ones((m, m))
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
            p_co_val = poisson.sf(both - 1, exp_cnt) if exp_cnt > 0 else (0.0 if both > 0 else 1.0)
            p_av_val = poisson.cdf(both, exp_cnt) if exp_cnt > 0 else 1.0
            p_co[i, j] = p_co_val
            p_co[j, i] = p_co_val
            p_av[i, j] = p_av_val
            p_av[j, i] = p_av_val
            lfc = np.log2(both / exp_cnt) if (exp_cnt > 0 and both > 0) else (-10.0 if (exp_cnt > 0 and both == 0) else 0.0)
            test_stats.append({
                "pair": f"{mi} vs {mj}",
                "mark_i": mi, "mark_j": mj,
                "observed": int(both), "expected": round(exp_cnt, 2),
                "occurrence_i": round(float(pi), 3), "occurrence_j": round(float(pj), 3),
                "n_regions": n,
                "p_cooccurrence": round(p_co_val, 6), "p_avoidance": round(p_av_val, 6),
                "log2_fold_change": round(lfc, 3),
            })
    df_test = pd.DataFrame(test_stats)
    pvals_co = df_test["p_cooccurrence"].values
    pvals_av = df_test["p_avoidance"].values
    _, pcor_co, _, _ = multipletests(pvals_co, method="fdr_bh")
    _, pcor_av, _, _ = multipletests(pvals_av, method="fdr_bh")
    df_test["p_cooccurrence_adj"] = pcor_co
    df_test["p_avoidance_adj"] = pcor_av
    df_test["cooccurrence_significant"] = (pcor_co < alpha) & (df_test["observed"] > df_test["expected"])
    df_test["avoidance_significant"] = (pcor_av < alpha) & (df_test["observed"] < df_test["expected"])
    return {
        "marks": marks, "occurrence": occurrence,
        "observed_matrix": obs_mat, "expected_matrix": exp_mat,
        "p_cooccurrence": p_co, "p_avoidance": p_av,
        "test_table": df_test, "alpha": alpha,
    }

for db_name in DATABASES:
    db = DATABASES[db_name]
    db["stats"] = compute_stats(db["df"], db["marks"], alpha=0.05)
    tt = db["stats"]["test_table"]
    n_sig_co = tt["cooccurrence_significant"].sum()
    n_sig_av = tt["avoidance_significant"].sum()
    print(f"\n{db_name}:")
    print(f"  N={db['n']}, marks={db['marks']}")
    print(f"  Occurrence: {dict(round(db['stats']['occurrence']*100, 1))}")
    print(f"  Significant co-occurrence pairs: {n_sig_co}")
    print(f"  Significant avoidance pairs: {n_sig_av}")
    for _, row in tt.iterrows():
        if row["cooccurrence_significant"]:
            print(f"    CO-OCCUR: {row['pair']} (obs={row['observed']}, exp={row['expected']}, lfc={row['log2_fold_change']}, p_adj={row['p_cooccurrence_adj']:.4f})")
        if row["avoidance_significant"]:
            print(f"    AVOID:    {row['pair']} (obs={row['observed']}, exp={row['expected']}, lfc={row['log2_fold_change']}, p_adj={row['p_avoidance_adj']:.4f})")

# ---- FIGURE 1: Comparative occurrence bar chart ----
fig, axes = plt.subplots(1, 4, figsize=(20, 4))
colors = ["#3498db", "#e74c3c", "#2ecc71", "#f39c12", "#9b59b6"]

for idx, (db_name, db) in enumerate(DATABASES.items()):
    ax = axes[idx]
    occ = db["stats"]["occurrence"]
    marks = db["marks"]
    bar_colors = colors[:len(marks)]
    bars = ax.bar(range(len(occ)), occ.values * 100, color=bar_colors, alpha=0.8, edgecolor="white", linewidth=1.2)
    ax.set_xticks(range(len(occ)))
    ax.set_xticklabels(occ.index, rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("Occurrence (%)", fontsize=9)
    ax.set_title(f"{db_name} (n={db['n']})", fontsize=10, fontweight="bold")
    for bar, val in zip(bars, occ.values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                f"{val*100:.1f}%", ha="center", va="bottom", fontsize=7)
    ax.set_ylim(0, max(occ.values) * 100 * 1.25)

fig.suptitle("Histone Mark Occurrence: Database Comparison", fontsize=12, fontweight="bold")
plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "fig1_db_occurrence_comparison.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("\nSaved fig1_db_occurrence_comparison.png")

# ---- FIGURE 2: Co-occurrence heatmaps (3 databases x 2 columns = obs, sig-dir) ----
fig, axes = plt.subplots(4, 3, figsize=(15, 16))
row_labels = list(DATABASES.keys())
col_labels = ["Observed co-occurrence", "Expected (independence)", "Observed - Expected"]

for row, (db_name, db) in enumerate(DATABASES.items()):
    stats = db["stats"]
    marks = stats["marks"]
    obs = stats["observed_matrix"]
    exp_m = stats["expected_matrix"]
    p_co = stats["p_cooccurrence"]
    p_av = stats["p_avoidance"]
    alpha = stats["alpha"]
    matrices = [obs, exp_m, obs - exp_m]
    cmaps = ["YlOrRd", "YlOrRd", "RdYlBu_r"]

    for col, (mat, cmap, title) in enumerate(zip(matrices, cmaps, col_labels)):
        ax = axes[row, col]
        im = ax.imshow(mat, cmap=cmap, aspect="equal")
        ax.set_xticks(range(len(marks)))
        ax.set_yticks(range(len(marks)))
        ax.set_xticklabels(marks, rotation=45, ha="right", fontsize=7)
        ax.set_yticklabels(marks, fontsize=7)

        # Add significance markers
        for i in range(len(marks)):
            for j in range(len(marks)):
                val = mat[i, j]
                if i == j:
                    ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6, color="gray")
                    continue
                stars = ""
                if p_co[i, j] < alpha: stars += "*"
                if p_av[i, j] < alpha: stars += "\u2020"
                txt = f"{val:.2f}{stars}" if stars else f"{val:.2f}"
                ax.text(j, i, txt, ha="center", va="center", fontsize=6,
                        color="black" if abs(val - 0.5) < 0.3 else "white")

        if row == 0:
            ax.set_title(title, fontsize=9, fontweight="bold")
        ax.set_ylabel(db_name if col == 0 else "", fontsize=9, fontweight="bold")

fig.suptitle("Database Co-occurrence Comparison (* p_co<alpha; \u2020 p_avoid<alpha)", fontsize=11, fontweight="bold", y=1.01)
plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "fig2_db_cooccurrence_heatmap.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved fig2_db_cooccurrence_heatmap.png")

# ---- FIGURE 3: Network graphs for each database ----
fig, axes = plt.subplots(1, 4, figsize=(20, 5))
for idx, (db_name, db) in enumerate(DATABASES.items()):
    ax = axes[idx]
    stats = db["stats"]
    marks = stats["marks"]
    obs = stats["observed_matrix"]
    p_co = stats["p_cooccurrence"]
    p_av = stats["p_avoidance"]
    alpha = stats["alpha"]

    G = nx.Graph()
    for m in marks: G.add_node(m)
    for i, mi in enumerate(marks):
        for j, mj in enumerate(marks):
            if i >= j: continue
            if p_co[i, j] < alpha and obs[i, j] > stats["expected_matrix"][i, j]:
                G.add_edge(mi, mj, color="#2ecc71", weight=2 + 5 * obs[i, j])
            elif p_av[i, j] < alpha and obs[i, j] < stats["expected_matrix"][i, j]:
                G.add_edge(mi, mj, color="#e74c3c", weight=2, style="dashed")

    if G.number_of_edges() > 0:
        pos = nx.circular_layout(G)
        edge_colors = [G.edges[e]["color"] for e in G.edges()]
        edge_widths = [G.edges[e]["weight"] for e in G.edges()]
        edge_styles = [G.edges[e].get("style", "solid") for e in G.edges()]
        nx.draw_networkx_nodes(G, pos, ax=ax, node_color="#bdc3c7", node_size=600, alpha=0.85)
        nx.draw_networkx_labels(G, pos, ax=ax, font_size=8, font_weight="bold")
        for (u, v), color, width, style in zip(G.edges(), edge_colors, edge_widths, edge_styles):
            nx.draw_networkx_edges(G, pos, ax=ax, edgelist=[(u, v)],
                                   edge_color=color, width=width, style=style, alpha=0.7)
    else:
        ax.text(0.5, 0.5, "No significant\nco-occurrence\nor avoidance", ha="center", va="center",
                fontsize=9, color="gray", transform=ax.transAxes)

    ax.set_title(f"{db_name}\n(n={db['n']})", fontsize=10, fontweight="bold")
    ax.axis("off")

# Shared legend
legend_elements = [
    Line2D([0], [0], color="#2ecc71", lw=3, label="Co-occurrence"),
    Line2D([0], [0], color="#e74c3c", lw=2, linestyle="dashed", label="Avoidance"),
]
fig.legend(handles=legend_elements, loc="lower center", ncol=2, fontsize=9)
fig.suptitle("Co-occurrence and Avoidance Networks by Database", fontsize=12, fontweight="bold", y=1.02)
plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "fig3_db_network_comparison.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved fig3_db_network_comparison.png")

# ---- Save summary table for manuscript ----
print("\n\n=== SUMMARY TABLE ===")
rows = []
for db_name, db in DATABASES.items():
    stats = db["stats"]
    tt = stats["test_table"]
    occ_str = "; ".join([f"{k}={v*100:.1f}%" for k, v in stats["occurrence"].items()])
    for _, row in tt.iterrows():
        rows.append({
            "Database": db_name,
            "Pair": row["pair"],
            "Observed": row["observed"],
            "Expected": row["expected"],
            "Log2FC": row["log2_fold_change"],
            "Co-occurrence (p_adj)": row["p_cooccurrence_adj"],
            "Avoidance (p_adj)": row["p_avoidance_adj"],
            "Co-occurrence significant": "Yes" if row["cooccurrence_significant"] else "",
            "Avoidance significant": "Yes" if row["avoidance_significant"] else "",
        })

df_summary = pd.DataFrame(rows)
summary_path = os.path.join(OUTPUT_DIR, "db_comparison_summary.csv")
df_summary.to_csv(summary_path, index=False)
print(f"\nSummary saved to {summary_path}")
print(df_summary.to_string(index=False))
print("\nDone.")
