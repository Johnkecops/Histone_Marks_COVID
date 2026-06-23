#!/usr/bin/env python3
"""
Module: histone_mark_mining
Purpose: Network analysis of histone mark co-occurrence across genomic regions
Author: Dr. Arli Aditya Parikesit
Date: 2026
Parameters:
    - output_dir: Directory for saving output figures (default: current working directory)
References:
    - Clauset A, Shalizi CR, Newman MEJ (2009) Power-law distributions in empirical data.
      SIAM Rev 51:661-703. https://doi.org/10.1137/070710111
    - Hagberg AA, Schult DA, Swart PJ (2008) Exploring network structure, dynamics,
      and function using NetworkX. Proc SciPy 2008, pp 11-15.
"""
import os
import re
import argparse
import defusedxml.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import requests
from bs4 import BeautifulSoup
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score, davies_bouldin_score
from sklearn.metrics.pairwise import cosine_similarity
from scipy.stats import ks_2samp, poisson
from scipy.optimize import curve_fit
from collections import Counter

# Biological labels for each simulated mark column.
# mark1 (H3K4me3): trimethylation of H3 Lys4 — active promoter mark, independent.
# mark2 (H3K27ac): acetylation of H3 Lys27 — active enhancer/promoter, correlated
#                  with H3K4me3 (simulated via 0.6*mark1 + noise).
# mark3 (H3K9me3): trimethylation of H3 Lys9 — repressive heterochromatin, independent.
_MARK_LABELS = {
    'mark1': 'H3K4me3',
    'mark2': 'H3K27ac',
    'mark3': 'H3K9me3',
}

# Function to simulate data mining from a repository website
def data_mining():
    # Simulating 100 genomic regions
    np.random.seed(42)
    n_samples = 100
    regions = [f"region_{i}" for i in range(n_samples)]

    # Generate 3 histone marks with some correlations to make it interesting
    mark1 = np.random.uniform(0.1, 0.9, n_samples)
    mark2 = 0.6 * mark1 + np.random.uniform(0.0, 0.4, n_samples)
    mark3 = np.random.uniform(0.1, 0.9, n_samples)

    # Assign each region its dominant histone modification (highest raw signal)
    mark_array = np.column_stack([mark1, mark2, mark3])
    mark_keys = ['mark1', 'mark2', 'mark3']
    dominant_marks = [_MARK_LABELS[mark_keys[i]] for i in np.argmax(mark_array, axis=1)]

    # Generate read counts
    counts = np.random.poisson(lam=25, size=n_samples)

    return {
        'region': regions,
        'histone_modification': dominant_marks,
        'organism_type': ['Human'] * n_samples,
        'species': ['Homo sapiens'] * n_samples,
        'mark1': list(mark1),
        'mark2': list(mark2),
        'mark3': list(mark3),
        'count': list(counts)
    }

# ---------------------------------------------------------------------------
# Organism classification helpers (shared across all live fetch functions)
# ---------------------------------------------------------------------------

def _classify_organism(scientific_name, lineage, division):
    if scientific_name == 'Homo sapiens':
        return 'Human'
    if 'Viridiplantae' in lineage:
        return 'Plant'
    if any(k in lineage for k in ('Bacteria', 'Archaea', 'Viruses', 'Fungi', 'Apicomplexa')):
        return 'Microbe'
    if division in ('Primates', 'Rodents', 'Mammals', 'Vertebrates', 'Invertebrates'):
        return 'Animal'
    return 'Other'

def _lookup_ncbi_taxids(taxids):
    """Batch-lookup NCBI taxonomy for a list of taxids.
    Returns {taxid: (scientific_name, organism_type)}.
    """
    if not taxids:
        return {}
    ids = ','.join(str(int(t)) for t in taxids)
    url = (
        'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
        f'?db=taxonomy&id={ids}&retmode=xml'
    )
    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        root = ET.fromstring(r.content)
        result = {}
        for taxon in root.findall('Taxon'):
            tid = int(taxon.findtext('TaxId') or 0)
            name = taxon.findtext('ScientificName') or 'Unknown'
            lineage = taxon.findtext('Lineage') or ''
            div = taxon.findtext('Division') or ''
            result[tid] = (name, _classify_organism(name, lineage, div))
        return result
    except Exception as e:
        print(f"Warning: NCBI taxonomy lookup failed: {e}")
        return {}


# PTM modification types to fetch from Histome2 Human DB.
# mark1 = n_writers, mark2 = n_erasers, mark3 = n_mod_types (site cross-regulation breadth).
_HISTOME2_BASE = 'https://www.actrec.gov.in/histome2/Human/'
_HISTOME2_MOD_TYPES = [
    'lysine_methylation',
    'lysine_acetylation',
    'arginine_methylation',
    'lysine_ubiquitination',
]
_HISTONE_PROTEIN_RE = re.compile(r'^(H[^KRSTY]+)')

def fetch_histome2_data():
    site_data = {}
    for mod_type in _HISTOME2_MOD_TYPES:
        url = f'{_HISTOME2_BASE}ptm.php?mod_type={mod_type}'
        try:
            response = requests.get(url, timeout=10)
            response.raise_for_status()
            soup = BeautifulSoup(response.content, 'html.parser')
            table = soup.find('table')
            if table is None:
                continue
            tbody = table.find('tbody')
            if tbody is None:
                continue
            for row in tbody.find_all('tr'):
                cells = row.find_all('td')
                if len(cells) < 3:
                    continue
                site = cells[0].get_text(strip=True)
                if not site:
                    continue
                writers = [a.get_text(strip=True) for a in cells[1].find_all('a')]
                erasers = [a.get_text(strip=True) for a in cells[2].find_all('a')]
                if site not in site_data:
                    site_data[site] = {'writers': set(), 'erasers': set(), 'mod_types': set()}
                site_data[site]['writers'].update(writers)
                site_data[site]['erasers'].update(erasers)
                site_data[site]['mod_types'].add(mod_type)
        except Exception as e:
            print(f"Warning: could not fetch {mod_type}: {e}")
            continue

    if not site_data:
        print("Warning: Histome2 returned no PTM data. Using fallback.")
        return None

    rows = []
    for site, info in sorted(site_data.items()):
        m = _HISTONE_PROTEIN_RE.match(site)
        histone_protein = m.group(1) if m else 'unknown'
        rows.append({
            'region': site,
            'histone_modification': histone_protein,
            'organism_type': 'Human',
            'species': 'Homo sapiens',
            'mark1': len(info['writers']),
            'mark2': len(info['erasers']),
            'mark3': len(info['mod_types']),
            'count': len(info['writers']) + len(info['erasers']),
        })

    print(f"Fetched {len(rows)} histone modification sites from Histome2.")
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# HistoneDB 2.0 (NCBI) — histone variant sequences with HMM classification
# mark1 = mean HMM score (confidence of variant assignment)
# mark2 = n_sequences (number of sequences classified as this variant)
# mark3 = n_taxa (number of distinct taxonomic IDs for this variant)
# ---------------------------------------------------------------------------

_HISTONEDB_BASE = (
    'https://www.ncbi.nlm.nih.gov/research/histonedb/static/browse/dumps/'
)

def fetch_histonedb_data():
    try:
        seqs = pd.read_csv(_HISTONEDB_BASE + 'seqs.txt')
        scores = pd.read_csv(_HISTONEDB_BASE + 'scores.txt')
    except Exception as e:
        print(f"Error fetching HistoneDB data: {e}")
        return None

    best_scores = (
        scores[scores['used_for_classification'] == True]
        .groupby('accession')['score']
        .max()
        .reset_index()
        .rename(columns={'score': 'best_score'})
    )

    merged = seqs.merge(best_scores, on='accession', how='left')

    agg = (
        merged.groupby(['hist_var', 'hist_type'])
        .agg(
            mean_score=('best_score', 'mean'),
            n_sequences=('accession', 'count'),
            n_taxa=('taxid', 'nunique'),
            dominant_taxid=('taxid', lambda x: x.mode()[0]),
        )
        .reset_index()
    )

    taxid_map = _lookup_ncbi_taxids(agg['dominant_taxid'].dropna().unique().tolist())

    agg = agg.rename(columns={'hist_var': 'region', 'hist_type': 'histone_modification'})
    agg['mark1'] = agg['mean_score'].round(2)
    agg['mark2'] = agg['n_sequences'].astype(float)
    agg['mark3'] = agg['n_taxa'].astype(float)
    agg['count'] = agg['n_sequences']
    agg['species'] = agg['dominant_taxid'].apply(
        lambda t: taxid_map.get(int(t), ('Unknown', 'Other'))[0]
    )
    agg['organism_type'] = agg['dominant_taxid'].apply(
        lambda t: taxid_map.get(int(t), ('Unknown', 'Other'))[1]
    )

    result = (
        agg[['region', 'histone_modification', 'organism_type', 'species',
             'mark1', 'mark2', 'mark3', 'count']]
        .dropna(subset=['mark1'])
        .sort_values('region')
        .reset_index(drop=True)
    )
    print(f"Fetched {len(result)} histone variants from HistoneDB 2.0.")
    return result


# ---------------------------------------------------------------------------
# ProHistoneDB — proto-histone categories in Archaea, Bacteria, and viruses
# mark1 = n_sequences (number of sequences in category)
# mark2 = n_taxa (number of distinct taxonomic groups)
# mark3 = mean_seq_len (mean amino acid sequence length)
# ---------------------------------------------------------------------------

_PROHISTONEDB_NS = {'p': 'http://www.phyloxml.org'}
_PROHISTONEDB_CATEGORIES = {
    1: 'Pacman_doublet',
    2: 'Nucleosomal',
    3: 'Xer_histone',
    4: 'Halo_doublet',
    5: 'Viral_singlet',
    6: 'RdgC_histone',
    7: 'Theion_coiled-coil_histone',
    8: 'Transmembrane_histone',
    9: 'Face-to-face',
    10: 'Poseidoniia_doublet',
    11: 'Methanococcales_histone',
    12: 'Bacterial_dimer',
    13: 'Beta-propeller_histone',
    14: 'Nanohalo_coiled-coil_histone',
    15: 'Rab_GTPase_histone',
    16: 'ZZ_histone',
    17: 'IHF_histone',
    18: 'Viral_doublet',
    21: 'Phage_histone',
    24: 'Bacterial_H2A_H2B',
    25: 'DUF1931',
    26: 'Dimer',
    27: 'Coiled-coil_histone',
    28: 'Miscellaneous',
}

def _broad_category(codename):
    n = codename.lower()
    if 'nucleosomal' in n:
        return 'Nucleosomal'
    if 'viral' in n or 'phage' in n:
        return 'Viral/Phage'
    if 'bacterial' in n or 'ihf' in n:
        return 'Bacterial'
    if any(k in n for k in ['halo', 'methanococcales', 'poseidoniia', 'nanohalo', 'theion']):
        return 'Archaeal'
    return 'Other'

def _fetch_prohistonedb_category(codename):
    url = f'https://prohistonedb.org/static/phylotrees/{codename}.xml'
    try:
        r = requests.get(url, timeout=15)
        r.raise_for_status()
        root = ET.fromstring(r.content)
        clades = root.findall('.//p:clade[p:name]', _PROHISTONEDB_NS)
        tax_counts = Counter()
        seq_lens = []
        for c in clades:
            tax = c.find('p:taxonomy/p:code', _PROHISTONEDB_NS)
            if tax is not None and tax.text:
                tax_counts[tax.text] += 1
            seq = c.find('p:sequence/p:mol_seq', _PROHISTONEDB_NS)
            if seq is not None and seq.text:
                seq_lens.append(len(seq.text))
        dominant_taxon = tax_counts.most_common(1)[0][0] if tax_counts else 'Unknown'
        return {
            'region': codename.replace('_', ' '),
            'histone_modification': _broad_category(codename),
            'organism_type': 'Microbe',
            'species': dominant_taxon,
            'mark1': float(len(clades)),
            'mark2': float(len(tax_counts)),
            'mark3': float(np.mean(seq_lens)) if seq_lens else 0.0,
            'count': len(clades),
        }
    except Exception as e:
        print(f"Warning: could not fetch {codename}: {e}")
        return None

def fetch_prohistonedb_data():
    rows = []
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {
            executor.submit(_fetch_prohistonedb_category, name): name
            for name in _PROHISTONEDB_CATEGORIES.values()
        }
        for future in as_completed(futures):
            result = future.result()
            if result is not None:
                rows.append(result)

    if not rows:
        print("Warning: ProHistoneDB returned no data.")
        return None

    df = (
        pd.DataFrame(rows)
        .sort_values('region')
        .reset_index(drop=True)
    )
    result = df[['region', 'histone_modification', 'organism_type', 'species',
                 'mark1', 'mark2', 'mark3', 'count']].copy()
    print(f"Fetched {len(result)} proto-histone categories from ProHistoneDB.")
    return result


# Main function to integrate and analyze data
def main(output_dir="."):
    os.makedirs(output_dir, exist_ok=True)
    print("=== Step 1: Data Collection, Preparation, and Normalization ===")
    
    # Attempt to fetch Histome2 data; fall back to simulated dataset
    histome2_df = fetch_histome2_data()
    if histome2_df is not None and not histome2_df.empty:
        print("Using Histome2 live data.")
        df = histome2_df
        if 'count' not in df.columns:
            df['count'] = np.random.poisson(lam=25, size=len(df))
    else:
        print("Using simulated histone mark mining dataset.")
        df = pd.DataFrame(data_mining())
    
    # Define histone mark columns
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    histone_marks = [col for col in numeric_cols if col.startswith('mark')]
    
    # Normalize: Standardize the data to ensure all histone marks are on a comparable scale
    print(f"Standardizing histone marks: {histone_marks}")
    scaler = StandardScaler()
    df[histone_marks] = scaler.fit_transform(df[histone_marks])
    
    print("\nDataFrame Preview (First 5 rows):")
    print(df.head())

    # Step 2: Graph Construction
    print("\n=== Step 2: Graph Construction & Centrality Measures ===")
    G = nx.DiGraph()
    
    # Add nodes for each region
    for _, row in df.iterrows():
        G.add_node(row['region'])
        
    # Calculate pairwise similarity matrix between genomic regions using normalized histone marks
    similarity_matrix = cosine_similarity(df[histone_marks].values)
    
    # Add edges if similarity exceeds threshold (0.95 gives a robust network with clear hubs)
    similarity_threshold = 0.95
    for i in range(len(df)):
        for j in range(len(df)):
            if i != j and similarity_matrix[i, j] > similarity_threshold:
                G.add_edge(df.iloc[i]['region'], df.iloc[j]['region'])
                
    print(f"Graph constructed with {G.number_of_nodes()} nodes (genomic regions) and {G.number_of_edges()} edges (interactions).")

    # Compute network centrality measures to identify key hubs / hotspots
    deg_centrality = nx.degree_centrality(G)
    between_centrality = nx.betweenness_centrality(G)
    closeness_centrality = nx.closeness_centrality(G)
    
    print("\nTop 5 Regions by Degree Centrality:")
    print("  ", sorted(deg_centrality.items(), key=lambda x: x[1], reverse=True)[:5])
    print("Top 5 Regions by Betweenness Centrality:")
    print("  ", sorted(between_centrality.items(), key=lambda x: x[1], reverse=True)[:5])
    print("Top 5 Regions by Closeness Centrality:")
    print("  ", sorted(closeness_centrality.items(), key=lambda x: x[1], reverse=True)[:5])

    # Step 3: Correlation Analysis
    print("\n=== Step 3: Correlation Analysis ===")
    # Calculate Pearson's correlation coefficient to assess relationships between histone marks
    corr_matrix = df[histone_marks].corr(method='pearson')
    print("Pearson Correlation Matrix of Histone Marks:")
    print(corr_matrix)

    # Step 4: Scale-Free Network Identification & Visualization
    print("\n=== Step 4: Scale-Free Network Identification & Visualization ===")
    # Calculate clustering coefficient
    avg_clustering_G = nx.average_clustering(G)
    
    # Compare against a random Erdos-Renyi network with the same size
    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()
    p = num_edges / (num_nodes * (num_nodes - 1)) if num_nodes > 1 else 0
    p = min(max(p, 0.0), 1.0)
    random_G = nx.fast_gnp_random_graph(num_nodes, p, directed=True, seed=42)
    avg_clustering_rand = nx.average_clustering(random_G)
    
    print(f"Average Clustering Coefficient (our network): {avg_clustering_G:.4f}")
    print(f"Average Clustering Coefficient (random network): {avg_clustering_rand:.4f}")

    # Visualize the graph
    plt.figure(figsize=(10, 8))
    pos = nx.spring_layout(G, k=0.15, seed=42)
    # Highlight high-degree hubs in red, others in lightblue
    node_degrees = [G.degree(node) for node in G.nodes()]
    max_deg = max(node_degrees) if node_degrees else 1
    node_colors = ['red' if G.degree(node) > 0.8 * max_deg else 'lightblue' for node in G.nodes()]
    
    nx.draw(G, pos, with_labels=False, node_color=node_colors, edge_color='gray', node_size=100, alpha=0.8, width=0.5)
    plt.title("Epigenetic Network of Histone Marks Similarity")
    net_path = os.path.join(output_dir, "scale_free_network.png")
    plt.savefig(net_path, dpi=300, bbox_inches='tight')
    print(f"Saved network visualization to '{net_path}'.")
    plt.close()

    # Step 5: Power-Law Distribution Test
    print("\n=== Step 5: Power-Law & Poisson Distribution Test ===")
    degree_values = [d for node, d in G.degree()]
    
    plt.figure(figsize=(8, 6))
    degree_counts = Counter(degree_values)
    x_deg = [deg for deg in degree_counts.keys() if deg > 0]
    y_freq = [degree_counts[deg] for deg in x_deg]
    
    if x_deg:
        plt.loglog(x_deg, y_freq, 'bo', label='Observed Degree Distribution')
    else:
        plt.plot([], [], 'bo', label='Observed Degree Distribution (No positive degrees)')
        
    plt.xlabel('Degree (log)')
    plt.ylabel('Frequency (log)')
    plt.title('Degree Distribution of the Network (Log-Log)')
    plt.legend()
    deg_path = os.path.join(output_dir, "degree_distribution.png")
    plt.savefig(deg_path, dpi=300, bbox_inches='tight')
    print(f"Saved degree distribution plot to '{deg_path}'.")
    plt.close()

    # Perform KS test comparing observed degrees to a Poisson distribution of similar mean
    mean_degree = np.mean(degree_values) if degree_values else 0
    poisson_sample = np.random.poisson(lam=mean_degree, size=len(degree_values)) if mean_degree > 0 else np.zeros(len(degree_values))
    ks_statistic, ks_pvalue = ks_2samp(degree_values, poisson_sample)
    print(f"KS Goodness-of-Fit Test (comparing to Poisson/random graph distribution):")
    print(f"  KS Statistic: {ks_statistic:.4f}")
    print(f"  P-value:      {ks_pvalue:.4f}")
    if ks_pvalue < 0.05:
        print("  Conclusion: The degree distribution differs significantly from a Poisson (random) network.")
    else:
        print("  Conclusion: The degree distribution does not differ significantly from a Poisson (random) network.")

    # Power-law exponent estimation via log-log OLS (Clauset et al. 2009 MLE preferred for publication;
    # use the `powerlaw` package for that).
    print("\nPower-Law Fit (log-log OLS):")
    if len(x_deg) >= 3:
        log_x = np.log(np.array(x_deg, dtype=float))
        log_y = np.log(np.array(y_freq, dtype=float))

        def _power_law_log(log_k, log_c, gamma):
            return log_c - gamma * log_k

        try:
            popt, pcov = curve_fit(_power_law_log, log_x, log_y)
            gamma_fit = popt[1]
            gamma_stderr = np.sqrt(pcov[1, 1])
            log_y_pred = _power_law_log(log_x, *popt)
            ss_res = np.sum((log_y - log_y_pred) ** 2)
            ss_tot = np.sum((log_y - np.mean(log_y)) ** 2)
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
            print(f"  Exponent gamma: {gamma_fit:.3f} ± {gamma_stderr:.3f}")
            print(f"  R² (log-log):   {r2:.4f}")
            print(f"  Note: Scale-free networks typically have gamma in [2, 3].")
        except Exception as e:
            print(f"  Fit failed: {e}")
    else:
        print("  Insufficient unique positive degrees for fitting (need >= 3).")

    # Step 6: Clustering and PCA
    print("\n=== Step 6: Clustering and PCA ===")
    features_for_pca = df[histone_marks]
    n_samples, n_features = features_for_pca.shape
    n_comps = min(2, n_samples, n_features)
    
    if n_comps > 0:
        pca = PCA(n_components=n_comps)
        reduced_data = pca.fit_transform(features_for_pca)
        print(f"PCA reduced data to {n_comps} dimensions.")
    else:
        reduced_data = features_for_pca.values
        print("PCA skipped due to insufficient dimensions/samples.")

    # Perform KMeans clustering to classify chromatin states
    n_clusters = min(3, n_samples)
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init='auto')
    df['cluster'] = kmeans.fit_predict(reduced_data)
    print(f"K-Means grouped samples into {n_clusters} clusters (chromatin states).")

    # Step 7: Cluster Quality Evaluation
    # Silhouette and Davies-Bouldin scores measure intrinsic cluster quality
    # without treating KMeans labels as an independent ground truth.
    # For supervised prediction, supply biological labels (e.g., ChIP-seq
    # peak calls, gene expression quartile) that are independent of histone
    # mark values used as features.
    print("\n=== Step 7: Cluster Quality Evaluation ===")
    features_eval = df[histone_marks].values
    labels = df['cluster'].values

    if len(set(labels)) > 1:
        sil = silhouette_score(features_eval, labels)
        db = davies_bouldin_score(features_eval, labels)
        print(f"Silhouette Score:      {sil:.4f}  (range [-1, 1]; closer to 1 = better separation)")
        print(f"Davies-Bouldin Index:  {db:.4f}  (lower = more compact, better-separated clusters)")
    else:
        print("Only one cluster found; quality metrics require >= 2 clusters.")

    # Step 8: Statistical Testing
    print("\n=== Step 8: Statistical Testing (Poisson Cluster Significance) ===")
    cluster_counts = df['cluster'].value_counts()

    for cluster in sorted(cluster_counts.index):
        count = cluster_counts[cluster]
        # Calculate mean from the 'count' column for samples in this cluster
        cluster_samples_count = df.loc[df['cluster'] == cluster, 'count']
        mean_val = cluster_samples_count.mean() if not cluster_samples_count.empty else 0
        
        # Calculate Poisson significance (survival function)
        p_value = poisson.sf(count, mean_val) if mean_val > 0 else 1.0
        print(f"Cluster {cluster}: Count = {count}, Mean of 'count' = {mean_val:.2f}, Poisson P-value = {p_value:.4f}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Histone mark network analysis pipeline")
    parser.add_argument("--output-dir", default=".", help="Directory for output figures (default: .)")
    args = parser.parse_args()
    main(output_dir=args.output_dir)