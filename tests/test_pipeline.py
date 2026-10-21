#!/usr/bin/env python3
"""
Module: test_pipeline
Purpose: Smoke tests for histone_mark_mining.py
Author: Dr. Arli Aditya Parikesit
Date: 2026
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import histone_mark_mining as _mod

import networkx as nx
import numpy as np
import pandas as pd
import pytest
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.preprocessing import StandardScaler


# ---------------------------------------------------------------------------
# data_mining()
# ---------------------------------------------------------------------------

def test_data_mining_keys():
    result = _mod.data_mining()
    expected = {
        "region", "histone_modification", "organism_type", "species",
        "mark1", "mark2", "mark3", "count",
    }
    assert set(result.keys()) == expected


def test_data_mining_length():
    result = _mod.data_mining()
    for key in result:
        assert len(result[key]) == 100, f"Column '{key}' has wrong length"


def test_data_mining_mark_ranges():
    result = _mod.data_mining()
    for col in ("mark1", "mark3"):
        arr = np.array(result[col])
        assert arr.min() >= 0.1 and arr.max() <= 0.9, f"{col} out of [0.1, 0.9]"


def test_data_mining_counts_nonnegative():
    result = _mod.data_mining()
    assert all(c >= 0 for c in result["count"])


def test_data_mining_reproducible():
    r1 = _mod.data_mining()
    r2 = _mod.data_mining()
    assert r1["mark1"] == r2["mark1"], "data_mining() not reproducible with fixed seed"


# ---------------------------------------------------------------------------
# fetch_histome2_data()
# ---------------------------------------------------------------------------

@pytest.mark.network
def test_fetch_histome2_returns_none_or_dataframe():
    result = _mod.fetch_histome2_data()
    assert result is None or isinstance(result, pd.DataFrame)


# ---------------------------------------------------------------------------
# Graph construction and centrality
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def built_graph():
    df = pd.DataFrame(_mod.data_mining())
    marks = ["mark1", "mark2", "mark3"]
    df[marks] = StandardScaler().fit_transform(df[marks])
    sim = cosine_similarity(df[marks].values)
    G = nx.DiGraph()
    for r in df["region"]:
        G.add_node(r)
    for i in range(len(df)):
        for j in range(len(df)):
            if i != j and sim[i, j] > 0.95:
                G.add_edge(df.iloc[i]["region"], df.iloc[j]["region"])
    return G, df["region"].tolist()


def test_graph_has_all_nodes(built_graph):
    G, regions = built_graph
    assert G.number_of_nodes() == 100


def test_centrality_covers_all_nodes(built_graph):
    G, regions = built_graph
    deg_c = nx.degree_centrality(G)
    between_c = nx.betweenness_centrality(G)
    close_c = nx.closeness_centrality(G)
    for c in (deg_c, between_c, close_c):
        assert set(c.keys()) == set(regions)


def test_degree_centrality_bounded(built_graph):
    G, _ = built_graph
    deg_c = nx.degree_centrality(G)
    assert all(0.0 <= v <= 1.0 for v in deg_c.values())


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))


# ---------------------------------------------------------------------------
# compute_cooccurrence_stats() (extracted via ast: generate_figures.py runs on import)
# ---------------------------------------------------------------------------

def test_cooccurrence_drops_constant_and_keeps_small_p():
    import ast
    from scipy.stats import poisson
    from statsmodels.stats.multitest import multipletests
    path = os.path.join(os.path.dirname(__file__), "..", "generate_figures.py")
    tree = ast.parse(open(path).read())
    fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "compute_cooccurrence_stats")
    ns = {"np": np, "pd": pd, "poisson": poisson, "multipletests": multipletests}
    exec(compile(ast.Module([fn], []), path, "exec"), ns)

    a = np.array([1] * 100 + [0] * 900)
    df = pd.DataFrame({"a": a, "b": a, "const": 0})
    res = ns["compute_cooccurrence_stats"](df, ["a", "b", "const"])
    assert res["marks"] == ["a", "b"]
    p_adj = res["test_table"]["p_cooccurrence_adj"].iloc[0]
    assert 0 < p_adj < 1e-6  # not floored to 0 by pre-FDR rounding
