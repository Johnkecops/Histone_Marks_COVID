#!/usr/bin/env python3
"""
Generate manuscript figures from the app-covid.py pipeline
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
import numpy as np
import pandas as pd

OUTPUT = "Histone Mark/figures"
os.makedirs(OUTPUT, exist_ok=True)

# ===== Figure 1: Occurrence =====
marks = ["H3K4me3", "H3K27ac", "H3K9me3", "H3K36me3", "H2AK5ac"]
occ = [0.35, 0.36, 0.28, 0.25, 0.225]
colors = ["#3498db", "#e74c3c", "#2ecc71", "#f39c12", "#9b59b6"]

fig, ax = plt.subplots(figsize=(8, 4))
bars = ax.bar(range(len(marks)), [v*100 for v in occ], color=colors, alpha=0.8, edgecolor="white", linewidth=1.2)
ax.set_xticks(range(len(marks)))
ax.set_xticklabels(marks, rotation=30, ha="right", fontsize=10)
ax.set_ylabel("Occurrence frequency (%)", fontsize=11)
ax.set_title("Histone Mark Occurrence Across Nucleosomes", fontsize=13, fontweight="bold")
for bar, val in zip(bars, occ):
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.5, f"{val*100:.1f}%", ha="center", va="bottom", fontsize=9)
ax.set_ylim(0, max(occ)*100*1.25)
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT, "figure1_occurrence.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("Figure 1 saved")

# ===== Figure 2: Co-occurrence heatmaps =====
# Generate data like in app-covid.py
np.random.seed(42)
n_nuc = 200
h3k4me3 = np.random.binomial(1, 0.35, n_nuc)
h3k27ac = np.array([1 if (h3k4me3[i]==1 and np.random.random()<0.75) or (h3k4me3[i]==0 and np.random.random()<0.15) else 0 for i in range(n_nuc)])
active = (h3k4me3 | h3k27ac).astype(bool)
h3k9me3 = np.array([1 if (not active[i] and np.random.random()<0.40) or (active[i] and np.random.random()<0.06) else 0 for i in range(n_nuc)])
h3k36me3 = np.random.binomial(1, 0.25, n_nuc)
h2ak5ac = np.array([1 if (h3k27ac[i]==1 and np.random.random()<0.45) or (h3k27ac[i]==0 and np.random.random()<0.10) else 0 for i in range(n_nuc)])

df = pd.DataFrame({"H3K4me3": h3k4me3, "H3K27ac": h3k27ac, "H3K9me3": h3k9me3, "H3K36me3": h3k36me3, "H2AK5ac": h2ak5ac})
marks_list = list(df.columns)
m = len(marks_list)
occ_s = df.mean()

obs = np.zeros((m,m))
exp = np.zeros((m,m))
for i, mi in enumerate(marks_list):
    for j, mj in enumerate(marks_list):
        if i == j:
            obs[i,j] = occ_s[mi]
            exp[i,j] = occ_s[mi]
            continue
        both = ((df[mi]==1)&(df[mj]==1)).sum()
        pi, pj = occ_s[mi], occ_s[mj]
        obs[i,j] = both/n_nuc
        exp[i,j] = pi*pj

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
titles = ["Observed co-occurrence", "Expected (independence)", "Difference"]
matrices = [obs, exp, obs - exp]
cmaps = ["YlOrRd", "YlOrRd", "RdYlBu_r"]

for ax, mat, title, cmap in zip(axes, matrices, titles, cmaps):
    im = ax.imshow(mat, cmap=cmap, aspect="equal")
    ax.set_xticks(range(m))
    ax.set_yticks(range(m))
    ax.set_xticklabels(marks_list, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(marks_list, fontsize=8)
    ax.set_title(title, fontsize=10, fontweight="bold")
    plt.colorbar(im, ax=ax, shrink=0.8)
    for i in range(m):
        for j in range(m):
            ax.text(j, i, f"{mat[i,j]:.2f}", ha="center", va="center", fontsize=7,
                    color="black" if abs(mat[i,j]-0.5)<0.3 else "white")

fig.suptitle("Histone Mark Co-occurrence Analysis", fontsize=12, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT, "figure2_cooccurrence_heatmap.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("Figure 2 saved")

# ===== Figure 3: Co-occurrence network =====
G = nx.Graph()
for m in marks_list:
    G.add_node(m)

# Add edges based on co-occurrence/avoidance
for i, mi in enumerate(marks_list):
    for j, mj in enumerate(marks_list):
        if i >= j: continue
        both = ((df[mi]==1)&(df[mj]==1)).sum()
        pi, pj = occ_s[mi], occ_s[mj]
        exp_cnt = n_nuc * pi * pj
        if exp_cnt > 0:
            from scipy.stats import poisson
            p_co = poisson.sf(both-1, exp_cnt)
            p_av = poisson.cdf(both, exp_cnt)
            if p_co < 0.05 and both > exp_cnt:
                G.add_edge(mi, mj, color="#2ecc71", weight=2+5*both/n_nuc)
            elif p_av < 0.05 and both < exp_cnt:
                G.add_edge(mi, mj, color="#e74c3c", weight=2, style="dashed")

fig, ax = plt.subplots(figsize=(8, 6))
pos = nx.circular_layout(G)
edge_colors = [G.edges[e]["color"] for e in G.edges()]
edge_widths = [G.edges[e]["weight"] for e in G.edges()]
nx.draw_networkx_nodes(G, pos, ax=ax, node_color="#bdc3c7", node_size=700, alpha=0.85)
nx.draw_networkx_labels(G, pos, ax=ax, font_size=10, font_weight="bold")
for (u,v), color, width in zip(G.edges(), edge_colors, edge_widths):
    style = G.edges[(u,v)].get("style", "solid")
    nx.draw_networkx_edges(G, pos, ax=ax, edgelist=[(u,v)], edge_color=color, width=width, style=style, alpha=0.7)

from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0],[0], color="#2ecc71", lw=3, label="Co-occurrence"),
    Line2D([0],[0], color="#e74c3c", lw=2, linestyle="dashed", label="Avoidance"),
]
ax.legend(handles=legend_elements, loc="upper right", fontsize=9)
ax.set_title("Co-occurrence and Avoidance Network", fontsize=13, fontweight="bold")
ax.axis("off")
fig.savefig(os.path.join(OUTPUT, "figure3_network.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("Figure 3 saved")

# ===== Figure 4: ORF8 mimicry network =====
G2 = nx.DiGraph()
G2.add_node("SARS-CoV-2 ORF8", color="#e74c3c", size=800)
G2.add_node("ARKS motif\n(H3 mimic)", color="#f39c12", size=600)
G2.add_node("Host chromatin", color="#3498db", size=700)
G2.add_node("H3K9me3/H3K27me3\n(repressive)", color="#2c3e50", size=500)
G2.add_node("IFN response genes\n(silenced)", color="#95a5a6", size=500)
G2.add_node("H3K9ac/H3K14ac\n(active depleted)", color="#7f8c8d", size=500)
G2.add_node("Chromatin\ncompaction", color="#e67e22", size=550)

edges = [
    ("SARS-CoV-2 ORF8", "ARKS motif\n(H3 mimic)"),
    ("ARKS motif\n(H3 mimic)", "Host chromatin"),
    ("Host chromatin", "H3K9me3/H3K27me3\n(repressive)"),
    ("H3K9me3/H3K27me3\n(repressive)", "IFN response genes\n(silenced)"),
    ("Host chromatin", "H3K9ac/H3K14ac\n(active depleted)"),
    ("H3K9ac/H3K14ac\n(active depleted)", "IFN response genes\n(silenced)"),
    ("Host chromatin", "Chromatin\ncompaction"),
    ("Chromatin\ncompaction", "IFN response genes\n(silenced)"),
]
for u,v in edges:
    G2.add_edge(u,v)

fig, ax = plt.subplots(figsize=(10, 7))
pos = nx.spring_layout(G2, k=1.2, seed=42, iterations=50)
node_colors = [G2.nodes[n]["color"] for n in G2.nodes()]
node_sizes = [G2.nodes[n]["size"] for n in G2.nodes()]
nx.draw_networkx_nodes(G2, pos, ax=ax, node_color=node_colors, node_size=node_sizes, alpha=0.85)
nx.draw_networkx_labels(G2, pos, ax=ax, font_size=8, font_weight="bold")
nx.draw_networkx_edges(G2, pos, ax=ax, edge_color="gray", alpha=0.4, arrows=True, arrowstyle="-|>", arrowsize=15, connectionstyle="arc3,rad=0.1")
ax.set_title("ORF8 Histone Mimicry: Epigenetic Subversion", fontsize=13, fontweight="bold")
ax.axis("off")
fig.savefig(os.path.join(OUTPUT, "figure4_orf8_mimicry.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("Figure 4 saved")

# ===== Figure 5: NETosis pathway =====
fig, ax = plt.subplots(figsize=(10, 7))
ax.set_xlim(0, 10)
ax.set_ylim(0, 7)
ax.axis("off")

steps = [
    (1.5, 6.0, "1. Activation", "Neutrophil triggered\nby viral PAMPs", "#3498db"),
    (4.0, 6.0, "2. ROS Production", "NADPH oxidase\ngenerates ROS", "#e74c3c"),
    (6.5, 6.0, "3. Citrullination", "PAD4 citrullinates H3;\nchromatin decondenses", "#e67e22"),
    (1.5, 3.5, "4. Dissolution", "Nuclear envelope breaks;\nchromatin + MPO/NE", "#9b59b6"),
    (4.0, 3.5, "5. Externalization", "Membrane ruptures;\nNET released", "#1abc9c"),
    (6.5, 3.5, "6. Unmasking", "DNase I digests DNA;\nhistones regain toxicity", "#2c3e50"),
]

for x, y, title, desc, color in steps:
    ax.add_patch(mpatches.FancyBboxPatch((x-0.7, y-0.5), 1.4, 1.0, boxstyle="round,pad=0.15", facecolor=color, alpha=0.15, edgecolor=color, linewidth=2))
    ax.text(x, y+0.15, title, ha="center", va="bottom", fontsize=9, fontweight="bold", color=color)
    ax.text(x, y-0.15, desc, ha="center", va="top", fontsize=7, color="gray")

arrows = [(2.2,6.5,3.3,6.5), (4.7,6.5,5.8,6.5), (2.2,5.5,2.2,4.0), (4.7,5.5,4.7,4.0), (6.5,5.5,6.5,4.0), (2.2,3.0,3.3,3.0), (4.7,3.0,5.8,3.0)]
for x1,y1,x2,y2 in arrows:
    ax.annotate("", xy=(x2,y2), xytext=(x1,y1), arrowprops=dict(arrowstyle="->", color="gray", lw=1.5, alpha=0.6))

ax.set_title("NETosis Pathway: The Double-Edged Sword", fontsize=13, fontweight="bold")
fig.savefig(os.path.join(OUTPUT, "figure5_netosis.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("Figure 5 saved")

print(f"\nAll figures saved to {OUTPUT}/")
