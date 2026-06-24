#!/usr/bin/env python3
"""
Generate figures from app-covid.py for manuscript use.
Author: Dr. Arli Aditya Parikesit
Date: 2026
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import poisson
from statsmodels.stats.multitest import multipletests
from matplotlib.lines import Line2D

# Constants from app-covid.py
NUCLEOSOME_MARKS = ["H3K4me3", "H3K27ac", "H3K9me3", "H3K36me3", "H2AK5ac"]

def generate_nucleosome_data(n_nuc=200, seed=42):
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

def compute_cooccurrence_stats(df_binary, marks, alpha=0.05):
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

output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "figures")
os.makedirs(output_dir, exist_ok=True)

df_nuc = generate_nucleosome_data(n_nuc=200, seed=42)
marks = [c for c in NUCLEOSOME_MARKS if c in df_nuc.columns]
cooc_stats = compute_cooccurrence_stats(df_nuc, marks, alpha=0.05)
occ = df_nuc[marks].mean()
obs = cooc_stats["observed_matrix"]
exp_m = cooc_stats["expected_matrix"]
p_co = cooc_stats["p_cooccurrence"]
p_av = cooc_stats["p_avoidance"]
alpha = cooc_stats["alpha"]

# === Figure 1: Occurrence barchart ===
fig, ax = plt.subplots(figsize=(8, 4))
colors = ["#3498db", "#e74c3c", "#2ecc71", "#f39c12", "#9b59b6"]
bars = ax.bar(range(len(occ)), occ.values * 100, color=colors[:len(occ)], alpha=0.8, edgecolor="white", linewidth=1.2)
ax.set_xticks(range(len(occ)))
ax.set_xticklabels(occ.index, rotation=30, ha="right", fontsize=9)
ax.set_ylabel("Occurrence frequency (%)", fontsize=10)
ax.set_title("Histone Mark Occurrence Across Nucleosomes", fontsize=12, fontweight="bold")
for bar, val in zip(bars, occ.values):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5, f"{val * 100:.1f}%", ha="center", va="bottom", fontsize=8)
ax.set_ylim(0, max(occ.values) * 100 * 1.2)
plt.tight_layout()
fig.savefig(os.path.join(output_dir, "fig1_occurrence.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved fig1_occurrence.png")

# === Figure 2: Co-occurrence heatmap ===
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
titles = ["Observed co-occurrence", "Expected (independence)", "Observed - Expected"]
matrices = [obs, exp_m, obs - exp_m]
cmaps = ["YlOrRd", "YlOrRd", "RdYlBu_r"]
for ax, mat, title, cmap in zip(axes, matrices, titles, cmaps):
    im = ax.imshow(mat, cmap=cmap, aspect="equal")
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
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7, color="gray")
                continue
            stars = ""
            if p_co[i, j] < alpha: stars += "*"
            if p_av[i, j] < alpha: stars += "\u2020"
            txt = f"{val:.2f}{stars}" if stars else f"{val:.2f}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=7, color="black" if abs(val - 0.5) < 0.3 else "white")
fig.suptitle("Histone Mark Co-occurrence Analysis", fontsize=11, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(output_dir, "fig2_cooccurrence_heatmap.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved fig2_cooccurrence_heatmap.png")

# === Figure 3: Co-occurrence network ===
G = nx.Graph()
for m in marks: G.add_node(m)
for i, mi in enumerate(marks):
    for j, mj in enumerate(marks):
        if i >= j: continue
        if p_co[i, j] < alpha and obs[i, j] > cooc_stats["expected_matrix"][i, j]:
            G.add_edge(mi, mj, color="#2ecc71", weight=2 + 5 * obs[i, j])
        elif p_av[i, j] < alpha and obs[i, j] < cooc_stats["expected_matrix"][i, j]:
            G.add_edge(mi, mj, color="#e74c3c", weight=2, style="dashed")

if G.number_of_edges() > 0:
    fig, ax = plt.subplots(figsize=(8, 6))
    pos = nx.circular_layout(G)
    edge_colors = [G.edges[e]["color"] for e in G.edges()]
    edge_widths = [G.edges[e]["weight"] for e in G.edges()]
    edge_styles = [G.edges[e].get("style", "solid") for e in G.edges()]
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color="#bdc3c7", node_size=700, alpha=0.85)
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=9, font_weight="bold")
    for (u, v), color, width, style in zip(G.edges(), edge_colors, edge_widths, edge_styles):
        nx.draw_networkx_edges(G, pos, ax=ax, edgelist=[(u, v)], edge_color=color, width=width, style=style, alpha=0.7)
    legend_elements = [
        Line2D([0], [0], color="#2ecc71", lw=3, label="Co-occurrence"),
        Line2D([0], [0], color="#e74c3c", lw=2, linestyle="dashed", label="Avoidance"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=8)
    ax.set_title("Co-occurrence and Avoidance Network", fontsize=12, fontweight="bold")
    ax.axis("off")
    fig.savefig(os.path.join(output_dir, "fig3_cooccurrence_network.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved fig3_cooccurrence_network.png")

# === Figure 4: NETosis pathway ===
fig, ax = plt.subplots(figsize=(10, 7))
ax.set_xlim(0, 10); ax.set_ylim(0, 7); ax.axis("off")
steps = [
    (1.5, 6.0, "1. Activation", "Neutrophil triggered by viral PAMPs", "#3498db"),
    (4.0, 6.0, "2. ROS Production", "NADPH oxidase generates ROS", "#e74c3c"),
    (6.5, 6.0, "3. Citrullination", "PAD4 citrullinates histones", "#e67e22"),
    (1.5, 3.5, "4. Dissolution", "Nuclear envelope breaks", "#9b59b6"),
    (4.0, 3.5, "5. Externalization", "NET released extracellularly", "#1abc9c"),
    (6.5, 3.5, "6. Unmasking", "DNase I digests DNA; histones cytotoxic", "#2c3e50"),
]
for x, y, title, desc, color in steps:
    ax.add_patch(mpatches.FancyBboxPatch((x - 0.7, y - 0.5), 1.4, 1.0, boxstyle="round,pad=0.15", facecolor=color, alpha=0.15, edgecolor=color, linewidth=2))
    ax.text(x, y + 0.15, title, ha="center", va="bottom", fontsize=9, fontweight="bold", color=color)
    ax.text(x, y - 0.15, desc, ha="center", va="top", fontsize=7, color="gray")
arrows = [(2.2, 6.5, 3.3, 6.5), (4.7, 6.5, 5.8, 6.5), (2.2, 5.5, 2.2, 4.0), (4.7, 5.5, 4.7, 4.0), (6.5, 5.5, 6.5, 4.0), (2.2, 3.0, 3.3, 3.0), (4.7, 3.0, 5.8, 3.0)]
for x1, y1, x2, y2 in arrows:
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="->", color="gray", lw=1.5, alpha=0.6))
ax.set_title("NETosis Pathway", fontsize=13, fontweight="bold")
fig.savefig(os.path.join(output_dir, "fig4_netosis.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved fig4_netosis.png")

# === Figure 5: ORF8 mimicry network ===
G2 = nx.DiGraph()
G2.add_node("SARS-CoV-2 ORF8", color="#e74c3c", size=800)
G2.add_node("ARKS motif (H3 mimic)", color="#f39c12", size=600)
G2.add_node("Host chromatin", color="#3498db", size=700)
G2.add_node("H3K9me3/H3K27me3", color="#2c3e50", size=500)
G2.add_node("IFN response genes", color="#95a5a6", size=500)
G2.add_node("H3K9ac/H3K14ac", color="#7f8c8d", size=500)
G2.add_node("Chromatin compaction", color="#e67e22", size=550)
edges2 = [("SARS-CoV-2 ORF8", "ARKS motif (H3 mimic)"), ("ARKS motif (H3 mimic)", "Host chromatin"), ("Host chromatin", "H3K9me3/H3K27me3"), ("H3K9me3/H3K27me3", "IFN response genes"), ("Host chromatin", "H3K9ac/H3K14ac"), ("H3K9ac/H3K14ac", "IFN response genes"), ("Host chromatin", "Chromatin compaction"), ("Chromatin compaction", "IFN response genes")]
for u, v in edges2: G2.add_edge(u, v)
fig, ax = plt.subplots(figsize=(10, 7))
pos2 = nx.spring_layout(G2, k=1.2, seed=42, iterations=50)
node_colors = [G2.nodes[n]["color"] for n in G2.nodes()]
node_sizes = [G2.nodes[n]["size"] for n in G2.nodes()]
nx.draw_networkx_nodes(G2, pos2, ax=ax, node_color=node_colors, node_size=node_sizes, alpha=0.85)
nx.draw_networkx_labels(G2, pos2, ax=ax, font_size=8, font_weight="bold")
nx.draw_networkx_edges(G2, pos2, ax=ax, edge_color="gray", alpha=0.4, arrows=True, arrowstyle="-|>", arrowsize=15, connectionstyle="arc3,rad=0.1")
ax.set_title("ORF8 Histone Mimicry: Epigenetic Subversion", fontsize=12, fontweight="bold")
ax.axis("off")
fig.savefig(os.path.join(output_dir, "fig5_orf8_mimicry.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved fig5_orf8_mimicry.png")

# === Figure 6: Immunothrombosis cycle ===
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_xlim(-2, 2); ax.set_ylim(-2, 2); ax.axis("off")
nodes_dict = {
    "Viral Infection": (0.0, 1.5), "Cell Death & Histone Release": (1.5, 0.5),
    "Complement Inhibition (C3/C5)": (1.5, -0.8), "Impaired Viral Clearance": (0.0, -1.5),
    "More Cell Death": (-1.5, -0.8), "Histone Release": (-1.5, 0.5),
}
for label, (x, y) in nodes_dict.items():
    fc = "#e74c3c" if "Histone" in label or "Cell Death" in label else "#f39c12" if "Complement" in label else "#3498db"
    ax.add_patch(mpatches.FancyBboxPatch((x - 0.7, y - 0.3), 1.4, 0.6, boxstyle="round,pad=0.1", facecolor=fc, alpha=0.15, edgecolor="gray", linewidth=1.5))
    ax.text(x, y, label, ha="center", va="center", fontsize=8, fontweight="bold")
cycle_edges = [("Viral Infection", "Cell Death & Histone Release"), ("Cell Death & Histone Release", "Complement Inhibition (C3/C5)"), ("Complement Inhibition (C3/C5)", "Impaired Viral Clearance"), ("Impaired Viral Clearance", "More Cell Death"), ("More Cell Death", "Histone Release"), ("Histone Release", "Viral Infection")]
for u, v in cycle_edges:
    if u in nodes_dict and v in nodes_dict:
        ax.annotate("", xy=nodes_dict[v], xytext=nodes_dict[u], arrowprops=dict(arrowstyle="->", color="#e74c3c", lw=2, alpha=0.6, connectionstyle="arc3,rad=0.3"))
ax.set_title("Immunothrombosis: The Vicious Circle", fontsize=13, fontweight="bold", color="#c0392b")
fig.savefig(os.path.join(output_dir, "fig6_immunothrombosis.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved fig6_immunothrombosis.png")

print(f"\nAll 6 figures saved to {output_dir}/")
