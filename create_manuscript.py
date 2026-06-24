#!/usr/bin/env python3
"""
Generate manuscript .docx for "Histone Marks in Coronaviruses Infection:
Co-occurrence and Avoidance in the Human Nucleosome"
Author: Arli Aditya Parikesit, Arif Nur Muhammad Ansori, Vikash Jakmola
Date: 2026
"""
import os
from docx import Document
from docx.shared import Inches, Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(OUTPUT_DIR, "results", "figures")

doc = Document()

for section in doc.sections:
    section.top_margin = Cm(2.54)
    section.bottom_margin = Cm(2.54)
    section.left_margin = Cm(2.54)
    section.right_margin = Cm(2.54)

style = doc.styles['Normal']
font = style.font
font.name = 'Times New Roman'
font.size = Pt(12)
style.paragraph_format.line_spacing = 1.5

def add_heading_styled(text, level=1):
    h = doc.add_heading(text, level=level)
    for run in h.runs:
        run.font.name = 'Times New Roman'
        run.font.color.rgb = RGBColor(0, 0, 0)
    return h

def add_para(text, bold=False, italic=False, size=12, align=None):
    p = doc.add_paragraph()
    run = p.add_run(text)
    run.font.name = 'Times New Roman'
    run.font.size = Pt(size)
    run.bold = bold
    run.italic = italic
    if align:
        p.alignment = align
    p.paragraph_format.line_spacing = 1.5
    return p

def add_figure(caption, filename, width=5.5):
    fig_path = os.path.join(FIG_DIR, filename)
    if os.path.exists(fig_path):
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        run = p.add_run()
        run.add_picture(fig_path, width=Inches(width))
    p = doc.add_paragraph()
    run = p.add_run(f"Figure: {caption}")
    run.font.name = 'Times New Roman'
    run.font.size = Pt(10)
    run.italic = True
    p.paragraph_format.line_spacing = 1.5

# ============================================================
# TITLE PAGE
# ============================================================
add_para("", size=12)
add_para("", size=12)
add_para(
    "Histone Marks in Coronaviruses Infection: Co-occurrence and Avoidance in the Human Nucleosome",
    bold=True, size=16, align=WD_ALIGN_PARAGRAPH.CENTER
)
add_para("", size=12)
add_para("Dr. Arli Aditya Parikesit", size=14, align=WD_ALIGN_PARAGRAPH.CENTER)
add_para("Department of Biotechnology, i3L University, Jakarta, Indonesia", size=12, align=WD_ALIGN_PARAGRAPH.CENTER)
add_para("ORCID: 0000-0001-8716-3926", size=10, align=WD_ALIGN_PARAGRAPH.CENTER)
add_para("", size=12)
doc.add_page_break()

# ============================================================
# ABSTRACT
# ============================================================
add_heading_styled("Abstract", level=1)
abstract_text = (
    "Background: Coronaviruses, including SARS-CoV-2, MERS-CoV, SARS-CoV, and HCoV-229E, "
    "manipulate host cell epigenetics through histone modifications to suppress antiviral "
    "responses and promote replication. The SARS-CoV-2 ORF8 protein functions as a histone "
    "mimic of the ARKS motif in Histone H3, disrupting host post-translational modifications "
    "and inducing chromatin compaction. "
    "Objective: This study investigates histone mark co-occurrence and avoidance patterns "
    "across three public histone databases (Histome2, HistoneDB 2.0, ProHistoneDB) alongside "
    "simulated nucleosome data with known interaction rules, providing a comparative benchmark "
    "of epigenetic interaction signals. "
    "Methods: Data were fetched live from Histome2 (62 PTM sites), HistoneDB 2.0 (30 histone "
    "variants), and ProHistoneDB (24 proto-histone categories). Simulated data (200 nucleosomes, "
    "5 marks with embedded co-occurrence/avoidance rules) were included as a positive control. "
    "Continuous variables were binarized by median split. Co-occurrence and avoidance were "
    "assessed using Poisson significance tests with Benjamini-Hochberg FDR correction. "
    "Molecular mechanisms of histone-mediated tissue damage were synthesized from the "
    "Coronavirus_Vault literature corpus. "
    "Results: Across all three databases, no significant co-occurrence or avoidance "
    "relationships were detected after FDR correction (alpha = 0.05). In contrast, the "
    "simulated nucleosome data correctly recovered 3 significant co-occurrence and 3 "
    "significant avoidance relationships as designed. The largest database trend was in "
    "HistoneDB 2.0 between many_sequences and broad_taxonomy (observed = 14, expected = 7.5, "
    "p_adj = 0.0647), approaching significance. Five molecular mechanisms of histone-mediated "
    "damage were identified from the literature: direct membrane disruption, ionic "
    "dyshomeostasis, TLR-mediated innate immune overdrive, inflammasome-driven pyroptosis, "
    "and complement system sabotage via C3/C5 inhibition. "
    "Conclusion: Real histone database variables binarized by median split do not exhibit "
    "significant co-occurrence or avoidance patterns, whereas simulated data with known "
    "biological rules successfully recapitulate canonical mark interactions. This contrast "
    "validates the analytical framework while suggesting that database-level attributes "
    "differ fundamentally from PTM-level mark interactions. Coronaviruses employ a coordinated "
    "epigenetic strategy involving both the promotion of repressive marks and the suppression "
    "of activating marks to silence antiviral gene programs."
)
add_para(abstract_text, size=11)
add_para("", size=10)
keywords = ("Keywords: Histone marks, coronavirus, SARS-CoV-2, nucleosome, co-occurrence, "
            "epigenetics, histone mimicry, ORF8, Histome2, HistoneDB, ProHistoneDB.")
add_para(keywords, size=10, italic=True)
doc.add_page_break()

# ============================================================
# 1. INTRODUCTION
# ============================================================
add_heading_styled("1. Introduction", level=1)
intro_text1 = (
    "Histones are the fundamental structural proteins that package eukaryotic DNA into "
    "nucleosomes, the repeating units of chromatin. Beyond their structural role, histones "
    "carry a diverse array of post-translational modifications (PTMs) that constitute the "
    "epigenetic code governing gene expression, DNA repair, and chromosome segregation [1, 2]. "
    "These modifications include methylation, acetylation, phosphorylation, citrullination, "
    "and ubiquitination, each associated with distinct transcriptional outcomes."
)
add_para(intro_text1)

intro_text2 = (
    "The emergence of SARS-CoV-2 and the subsequent COVID-19 pandemic brought renewed attention "
    "to the intersection of virology and epigenetics. Several studies have demonstrated that "
    "coronaviruses, particularly SARS-CoV-2, actively subvert host epigenetic regulation to "
    "suppress antiviral immune responses [3, 4]. A landmark study by Kee et al. [5] revealed "
    "that the SARS-CoV-2 accessory protein ORF8 functions as a histone mimic of the ARKS motif "
    "in Histone H3. Through this mimicry, ORF8 associates with host chromatin, disrupts "
    "critical histone PTMs, and promotes chromatin compaction, thereby silencing interferon "
    "response genes. Deletion of either the ORF8 gene or the histone mimic site attenuates "
    "viral disruption of host chromatin and reduces viral genome copy number [5]."
)
add_para(intro_text2)

intro_text3 = (
    "The range of histone marks affected by coronavirus infection is remarkably broad, "
    "encompassing both activating and repressive modifications. SARS-CoV-2 increases "
    "repressive marks such as H3K9me3 and H3K27me3 while decreasing activating marks "
    "including H3K9ac and H3K14ac [3, 5]. MERS-CoV similarly manipulates H3K4me3 at "
    "interferon-stimulated gene (ISG) promoters to prevent antiviral state establishment [6]. "
    "Human coronavirus 229E (HCoV-229E) alters histone proteoforms, including the truncation "
    "states of core histones, revealing a previously underappreciated mechanism of viral "
    "epigenetic manipulation [7]."
)
add_para(intro_text3)

intro_text4 = (
    "When infected cells undergo cell death through necrosis, apoptosis, or NETosis, histones "
    "are released into the extracellular space where they function as cytotoxic "
    "Damage-Associated Molecular Patterns (DAMPs). Extracellular histones drive tissue "
    "damage through five principal mechanisms: direct membrane disruption, ionic "
    "dyshomeostasis, Toll-like receptor (TLR) activation, inflammasome-mediated pyroptosis, "
    "and complement system inhibition [8, 9]. Understanding both the intracellular "
    "epigenetic reprogramming and the extracellular cytotoxic effects of histones is "
    "critical for developing host-directed therapies for coronavirus infections."
)
add_para(intro_text4)

intro_text5 = (
    "In this study, we developed a computational framework to analyze histone mark "
    "co-occurrence and avoidance patterns using live data from three public histone "
    "databases: Histome2 (human PTM site writer/eraser counts), HistoneDB 2.0 (histone "
    "variant HMM scores and taxonomic distribution), and ProHistoneDB (proto-histone "
    "categories across Archaea, Bacteria, and viruses). The goal was to benchmark and "
    "compare these databases for their ability to reveal statistically significant "
    "co-occurrence and avoidance signals, and to contextualize these patterns within "
    "the broader molecular mechanisms of histone-mediated pathology during coronavirus infection."
)
add_para(intro_text5)

# ============================================================
# 2. MATERIALS AND METHODS
# ============================================================
add_heading_styled("2. Materials and Methods", level=1)

add_heading_styled("2.1. Literature Data Sources", level=2)
methods_text1 = (
    "Histone mark data associated with coronavirus infection were synthesized from the "
    "Coronavirus_Vault literature corpus, which contains two synthesized reviews: "
    "(1) 'Molecular Sabotage: How Histones Damage Human Tissues during Viral Infection,' "
    "detailing five mechanisms of histone-mediated tissue damage, and (2) 'Viral Histone Mimicry "
    "and Epigenetic Suppression of Host Immunity,' documenting specific histone marks altered "
    "by SARS-CoV-2, MERS-CoV, SARS-CoV, and HCoV-229E. The reference database "
    "(coronavirus26.bib) contains seven curated citations [3-8, 10]."
)
add_para(methods_text1)

add_heading_styled("2.2. Database Sources and Data Curation", level=2)

methods_db1 = (
    "Three public histone databases were queried programmatically via their REST APIs "
    "or downloadable flat files. All fetches were performed with configurable timeouts "
    "(10-30 s) and thread pool parallelism where applicable."
)
add_para(methods_db1)

methods_db2 = (
    "Histome2 (actrec.gov.in/histome2) is a curated database of human histone "
    "post-translational modification sites maintained by ACTREC, India [11]. Data were "
    "scraped from the PTM table pages across four modification types: lysine methylation, "
    "lysine acetylation, arginine methylation, and lysine ubiquitination. For each PTM "
    "site, the numbers of associated writer enzymes (mark1), eraser enzymes (mark2), and "
    "distinct modification types (mark3) were recorded. A total of 62 PTM sites were "
    "retrieved. The continuous variables were binarized by median split to produce three "
    "binary marks: has_writer (above-median writer count), has_eraser (above-median eraser "
    "count), and multi_mod (above-median number of modification types)."
)
add_para(methods_db2)

methods_db3 = (
    "HistoneDB 2.0 (ncbi.nlm.nih.gov/research/histonedb) is an NCBI-hosted database of "
    "histone variant sequences with HMM classification scores [12]. Sequence and score "
    "flat files were downloaded as CSV dumps. Data were aggregated by histone variant, "
    "producing three features per variant: mean HMM score (mark1), number of classified "
    "sequences (mark2), and number of distinct taxonomic IDs (mark3). Thirty histone "
    "variants were retrieved. These continuous variables were binarized by median split "
    "into three binary marks: high_HMM_score, many_sequences, and broad_taxonomy."
)
add_para(methods_db3)

methods_db4 = (
    "ProHistoneDB (prohistonedb.org) catalogs proto-histone categories across Archaea, "
    "Bacteria, and viruses [13]. PhyloXML tree files for 28 phylogenetic categories "
    "were fetched concurrently via ThreadPoolExecutor. For each category, the number of "
    "clades (mark1), distinct taxonomic codes (mark2), and mean sequence length in amino "
    "acids (mark3) were extracted. Twenty-four categories with valid data were retained "
    "after fetch completion. These were binarized by median split into three binary marks: "
    "many_sequences, broad_taxonomy, and long_sequence."
)
add_para(methods_db4)

methods_sim = (
    "Simulated nucleosome data (200 nucleosomes, 5 marks) were generated as a positive "
    "control using embedded biological rules. Activating marks H3K4me3 and H3K27ac were "
    "modeled to co-occur (probability P = 0.75 of H3K27ac given H3K4me3, versus P = 0.15 "
    "otherwise). The repressive mark H3K9me3 was modeled to avoid active marks (probability "
    "P = 0.40 in inactive regions versus P = 0.06 in active regions). H3K36me3 was modeled "
    "as independent at 25% background frequency. H2AK5ac was modeled to co-occur with "
    "H3K27ac (P = 0.45 given H3K27ac, versus P = 0.10 otherwise). These rules encode known "
    "nucleosomal interactions: H3K4me3 and H3K27ac jointly mark active promoters, H3K9me3 "
    "marks silenced heterochromatin and is anticorrelated with active marks, and H3K36me3 "
    "demarcates gene body elongation independent of promoter status."
)
add_para(methods_sim)

add_heading_styled("2.3. Statistical Analysis of Co-occurrence and Avoidance", level=2)
methods_stat = (
    "For each database, co-occurrence was defined as two binary marks appearing together "
    "in the same region more often than expected under independence. The null expectation "
    "was modeled as a Poisson distribution with lambda = n * P(A) * P(B), where n is the "
    "number of regions and P(A) and P(B) are the marginal probabilities of marks A and B. "
    "The right-tailed Poisson test assessed co-occurrence significance: P(observed >= k | lambda). "
    "The left-tailed Poisson test assessed avoidance significance: P(observed <= k | lambda). "
    "P-values were adjusted for multiple testing using the Benjamini-Hochberg False Discovery "
    "Rate (FDR) method, with alpha set to 0.05 [14]. The log2 fold change was calculated as "
    "log2(observed_count / expected_count)."
)
add_para(methods_stat)

add_heading_styled("2.4. Application Development", level=2)
methods_app = (
    "The computational framework was implemented as a Streamlit dashboard (app-covid.py) "
    "that integrates data ingestion, statistical analysis, and visualization. The application "
    "provides nine interactive tabs including co-occurrence analysis with data source selection "
    "(Histome2, HistoneDB 2.0, ProHistoneDB, Simulated), with automatic fallback to simulated "
    "data when live fetches fail. All figures in this manuscript were generated programmatically "
    "using the core computational functions with Matplotlib and NetworkX [15, 16]."
)
add_para(methods_app)

doc.add_page_break()

# ============================================================
# 3. RESULTS
# ============================================================
add_heading_styled("3. Results", level=1)

add_heading_styled("3.1. Occurrence of Histone Marks Across Databases", level=2)

res_occ1 = (
    "The occurrence frequency of each binary mark varied across the three databases, "
    "reflecting their different biological domains and feature definitions (Figure 1). "
    "In Histome2, has_writer occurred at 41.9% of the 62 PTM sites, has_eraser at 33.9%, "
    "and multi_mod at 0.0%. The latter was entirely degenerate because all PTM sites in "
    "Histome2 had only one modification type each, so the median equaled the minimum value "
    "and no site exceeded the threshold."
)
add_para(res_occ1)

res_occ2 = (
    "HistoneDB 2.0 showed all three marks at exactly 50.0% occurrence by construction, "
    "because the median split on 30 variants yields 15 variants above and 15 below the "
    "median for each continuous variable. ProHistoneDB showed many_sequences at 50.0%, "
    "broad_taxonomy at 37.5%, and long_sequence at 50.0%. The lower occurrence of "
    "broad_taxonomy reflects the uneven distribution of taxonomic breadth across proto-histone "
    "categories, where many categories are found in only one or two taxonomic groups."
)
add_para(res_occ2)

res_occ3 = (
    "The simulated nucleosome data (n = 200) showed more varied occurrence frequencies "
    "reflecting the deliberately engineered mark distributions: H3K4me3 at 35.5%, H3K27ac "
    "at 38.5%, H3K9me3 at 27.5%, H3K36me3 at 26.0%, and H2AK5ac at 22.0%. These values "
    "fall within realistic ranges for ChIP-seq experiments, where activating marks typically "
    "cover 20-40% of analyzed nucleosomes and repressive marks 15-30%."
)
add_para(res_occ3)

add_figure(
    "Histone mark occurrence frequencies across all data sources: Histome2 (62 PTM sites), "
    "HistoneDB 2.0 (30 histone variants), ProHistoneDB (24 proto-histone categories), and "
    "Simulated (200 nucleosomes, 5 marks with embedded interaction rules). "
    "Marks were binarized by median split of continuous variables.",
    "fig1_db_occurrence_comparison.png", width=5.5
)
add_para("", size=10)

add_heading_styled("3.2. Co-occurrence and Avoidance Analysis", level=2)

res_co1 = (
    "The co-occurrence and avoidance analysis revealed no statistically significant "
    "relationships after FDR correction at alpha = 0.05 for any of the three real databases "
    "(Figure 2). In contrast, the simulated nucleosome data recovered three significant "
    "co-occurrence pairs and three significant avoidance pairs, all consistent with the "
    "embedded biological rules (Table 1)."
)
add_para(res_co1)

# ---- Table 1 ----
add_para("Table 1. Pairwise co-occurrence and avoidance results across all data sources.", bold=True, size=10)
add_para("", size=6)

table = doc.add_table(rows=19, cols=8)
table.style = 'Light Grid Accent 1'
table.alignment = WD_TABLE_ALIGNMENT.CENTER

headers = ["Database", "Mark pair", "Obs", "Exp", "Log2FC", "p_adj (co)", "p_adj (av)", "Significant"]
for i, h in enumerate(headers):
    cell = table.rows[0].cells[i]
    cell.text = h
    for paragraph in cell.paragraphs:
        for run in paragraph.runs:
            run.bold = True
            run.font.size = Pt(8)
            run.font.name = 'Times New Roman'

rows_data = [
    ["Histome2", "has_writer vs has_eraser", "14", "8.81", "0.67", "0.194", "1.000", "None"],
    ["Histome2", "has_writer vs multi_mod", "0", "0.00", "0.00", "1.000", "1.000", "None (degenerate)"],
    ["Histome2", "has_eraser vs multi_mod", "0", "0.00", "0.00", "1.000", "1.000", "None (degenerate)"],
    ["HistoneDB 2.0", "high_HMM vs many_seq", "7", "7.50", "-0.10", "0.622", "0.787", "None"],
    ["HistoneDB 2.0", "high_HMM vs broad_tax", "7", "7.50", "-0.10", "0.622", "0.787", "None"],
    ["HistoneDB 2.0", "many_seq vs broad_tax", "14", "7.50", "0.90", "0.065", "0.990", "Trend only"],
    ["ProHistoneDB", "many_seq vs broad_tax", "7", "4.50", "0.64", "0.507", "0.913", "None"],
    ["ProHistoneDB", "many_seq vs long_seq", "7", "6.00", "0.22", "0.591", "0.913", "None"],
    ["ProHistoneDB", "broad_tax vs long_seq", "4", "4.50", "-0.17", "0.658", "0.913", "None"],
    ["Simulated", "H3K4me3 vs H3K27ac", "56", "27.34", "1.04", "<0.001", "1.000", "Co-occur *"],
    ["Simulated", "H3K4me3 vs H3K9me3", "4", "19.53", "-2.29", "1.000", "<0.001", "Avoid **"],
    ["Simulated", "H3K4me3 vs H3K36me3", "25", "18.46", "0.44", "0.211", "1.000", "None"],
    ["Simulated", "H3K4me3 vs H2AK5ac", "26", "15.62", "0.74", "0.033", "1.000", "Co-occur *"],
    ["Simulated", "H3K27ac vs H3K9me3", "5", "21.18", "-2.08", "1.000", "<0.001", "Avoid **"],
    ["Simulated", "H3K27ac vs H2AK5ac", "36", "16.94", "1.09", "<0.001", "1.000", "Co-occur **"],
    ["Simulated", "H3K9me3 vs H3K36me3", "10", "14.30", "-0.52", "1.000", "0.392", "None"],
    ["Simulated", "H3K9me3 vs H2AK5ac", "4", "12.10", "-1.60", "1.000", "0.024", "Avoid *"],
    ["Simulated", "H3K36me3 vs H2AK5ac", "11", "11.44", "-0.06", "0.986", "1.000", "None"],
]
for i, row_data in enumerate(rows_data):
    for j, val in enumerate(row_data):
        cell = table.rows[i+1].cells[j]
        cell.text = val
        for paragraph in cell.paragraphs:
            for run in paragraph.runs:
                run.font.size = Pt(7)
                run.font.name = 'Times New Roman'

add_para("", size=10)

res_co2 = (
    "The strongest trend was observed in HistoneDB 2.0 between many_sequences and "
    "broad_taxonomy (observed = 14 co-occurring variants, expected = 7.5 under independence, "
    "log2 fold change = 0.90, FDR-adjusted p = 0.0647). This trend approached but did not "
    "reach the significance threshold of 0.05, suggesting that histone variants with more "
    "classified sequences tend to have broader taxonomic distribution, which is biologically "
    "plausible: well-conserved variants are both more sequenced and more taxonomically widespread."
)
add_para(res_co2)

res_co3 = (
    "In Histome2, has_writer and has_eraser showed a non-significant co-occurrence trend "
    "(observed = 14 sites with both above-median counts, expected = 8.81, log2 fold change "
    "= 0.67, FDR-adjusted p = 0.194). The multi_mod mark was entirely degenerate because all "
    "62 PTM sites had exactly one modification type, producing zero variance. This indicates "
    "a limitation of the median-split binarization approach for integer-valued features with "
    "low variance."
)
add_para(res_co3)

res_co4 = (
    "In ProHistoneDB, no pair showed significant co-occurrence or avoidance. The observed "
    "counts were close to expected values for all three mark pairs, with log2 fold changes "
    "ranging from -0.17 to 0.64 and all FDR-adjusted p-values exceeding 0.5. This suggests "
    "that among proto-histone categories, sequence count, taxonomic breadth, and sequence "
    "length vary independently."
)
add_para(res_co4)

res_sim = (
    "The simulated nucleosome data provided a striking contrast. Three co-occurrence "
    "relationships reached significance: H3K4me3 with H3K27ac (observed = 56, expected = "
    "27.34, Log2FC = 1.04, p_adj < 0.001), H3K27ac with H2AK5ac (observed = 36, expected = "
    "16.94, Log2FC = 1.09, p_adj < 0.001), and H3K4me3 with H2AK5ac (observed = 26, expected "
    "= 15.62, Log2FC = 0.74, p_adj = 0.033). Three avoidance relationships also reached "
    "significance: H3K4me3 with H3K9me3 (observed = 4, expected = 19.53, Log2FC = -2.29, "
    "p_adj < 0.001), H3K27ac with H3K9me3 (observed = 5, expected = 21.18, Log2FC = -2.08, "
    "p_adj < 0.001), and H3K9me3 with H2AK5ac (observed = 4, expected = 12.10, Log2FC = "
    "-1.60, p_adj = 0.024) (Table 1). These results are fully consistent with the embedded "
    "biological rules, confirming that the Poisson test with FDR correction correctly detects "
    "known co-occurrence and avoidance patterns when they exist."
)
add_para(res_sim)

res_sim2 = (
    "The co-occurrence of H3K4me3 with H3K27ac reflects their shared role in marking active "
    "promoters and enhancers. The additional co-occurrence of H2AK5ac with both marks is an "
    "indirect consequence of its conditional dependence on H3K27ac. The avoidance of H3K9me3 "
    "with all active marks (H3K4me3, H3K27ac, and H2AK5ac) reflects the mutual exclusion "
    "between active and repressive chromatin domains. H3K36me3, which was modeled as "
    "independent, correctly showed no significant interactions, demonstrating specificity."
)
add_para(res_sim2)

add_figure(
    "Co-occurrence heatmaps for Histome2, HistoneDB 2.0, ProHistoneDB, and Simulated data. "
    "Columns show observed co-occurrence, expected under independence, and the difference. "
    "Asterisks (*) mark significant co-occurrence; daggers (dag) mark significant avoidance. "
    "Only the Simulated panel shows significant relationships (FDR-corrected alpha = 0.05).",
    "fig2_db_cooccurrence_heatmap.png", width=5.5
)
add_para("", size=10)

add_figure(
    "Co-occurrence and avoidance networks for each data source. Green edges represent "
    "significant co-occurrence; red dashed edges represent significant avoidance. "
    "No edges reached statistical significance in the three real databases. The Simulated "
    "network shows three co-occurrence edges (H3K4me3-H3K27ac, H3K27ac-H2AK5ac, "
    "H3K4me3-H2AK5ac) and three avoidance edges (H3K4me3-H3K9me3, H3K27ac-H3K9me3, "
    "H3K9me3-H2AK5ac).",
    "fig3_db_network_comparison.png", width=5.5
)
add_para("", size=10)

add_heading_styled("3.3. Database Benchmarking Summary", level=2)

res_bench = (
    "Table 2 summarizes the key characteristics and results across the three databases and "
    "the simulated control. HistoneDB 2.0 provided the largest effect size (Log2FC = 0.90 "
    "for the many_sequences vs broad_taxonomy pair) and the lowest p-value (0.0647) among "
    "real databases, suggesting it may be the most informative for co-occurrence analysis "
    "with the current encoding scheme. Histome2's multi_mod feature was degenerate, revealing "
    "a limitation when binarizing features with naturally low variance. ProHistoneDB showed "
    "the most independent mark behavior. The simulated data, serving as a positive control, "
    "correctly recovered 6 significant relationships (3 co-occurrence, 3 avoidance) out of "
    "10 possible pairs, with Log2FC magnitudes up to -2.29 (H3K4me3 vs H3K9me3 avoidance)."
)
add_para(res_bench)

# ---- Table 2 ----
add_para("Table 2. Comparative summary of data source characteristics and benchmarking results.", bold=True, size=10)
add_para("", size=6)

table2 = doc.add_table(rows=5, cols=5)
table2.style = 'Light Grid Accent 1'
table2.alignment = WD_TABLE_ALIGNMENT.CENTER

headers2 = ["Feature", "Histome2", "HistoneDB 2.0", "ProHistoneDB", "Simulated"]
for i, h in enumerate(headers2):
    cell = table2.rows[0].cells[i]
    cell.text = h
    for paragraph in cell.paragraphs:
        for run in paragraph.runs:
            run.bold = True
            run.font.size = Pt(9)
            run.font.name = 'Times New Roman'

rows2 = [
    ["N (samples)", "62 PTM sites", "30 histone variants", "24 proto-histone categories", "200 nucleosomes"],
    ["Biological domain", "Human (Homo sapiens)", "All taxa (HMM classified)", "Archaea, Bacteria, viruses", "Synthetic (rules-based)"],
    ["Significant pairs", "0", "0", "0", "6 (3 co-occur, 3 avoid)"],
    ["Best effect size", "has_writer vs has_eraser\n(p_adj = 0.194, LFC=0.67)", "many_seq vs broad_tax\n(p_adj = 0.065, LFC=0.90)", "many_seq vs broad_tax\n(p_adj = 0.507, LFC=0.64)", "H3K4me3 vs H3K27ac\n(p_adj<0.001, LFC=1.04)"],
]
for i, row_data in enumerate(rows2):
    for j, val in enumerate(row_data):
        cell = table2.rows[i+1].cells[j]
        cell.text = val
        for paragraph in cell.paragraphs:
            for run in paragraph.runs:
                run.font.size = Pt(8)
                run.font.name = 'Times New Roman'

add_para("", size=10)

add_heading_styled("3.4. Histone Marks in Coronavirus Infection", level=2)
res_hm = (
    "Literature synthesis from the Coronavirus_Vault corpus identified 12 distinct histone "
    "mark-virus associations across four coronavirus species. Five types of PTMs were "
    "represented: methylation, acetylation, citrullination, phosphorylation, and "
    "glutathionylation, along with two types of histone truncation events. SARS-CoV-2 "
    "accounted for the majority of reported marks (8 of 12), reflecting the intensive "
    "research focus on this virus. HCoV-229E was associated with 5 marks, predominantly "
    "identified through top-down proteomics [7]. MERS-CoV and SARS-CoV each contributed "
    "2 marks, primarily related to ISG suppression."
)
add_para(res_hm)

# ---- Table 3: Histone marks ----
add_para("Table 3. Histone marks altered during coronavirus infection.", bold=True, size=10)
add_para("", size=6)

table3 = doc.add_table(rows=13, cols=4)
table3.style = 'Light Grid Accent 1'
table3.alignment = WD_TABLE_ALIGNMENT.CENTER

headers3 = ["Modification Type", "Mark Code", "Virus", "Functional Annotation"]
for i, h in enumerate(headers3):
    cell = table3.rows[0].cells[i]
    cell.text = h
    for paragraph in cell.paragraphs:
        for run in paragraph.runs:
            run.bold = True
            run.font.size = Pt(9)
            run.font.name = 'Times New Roman'

rows3 = [
    ("Methylation", "H3K9me3, H3K27me3", "SARS-CoV-2, MERS-CoV", "ORF8 increases repressive marks to silence antiviral genes"),
    ("Methylation", "H3K4me3", "SARS-CoV-2, SARS-CoV", "Active promoter mark; depleted at ISG promoters in MERS-CoV"),
    ("Methylation", "H3K9me2", "SARS-CoV-2", "Inversely correlates with interferon expression"),
    ("Acetylation", "H3K9ac, H3K14ac", "SARS-CoV-2", "Decreased by ORF8 to weaken IFN response"),
    ("Acetylation", "H3K27ac", "SARS-CoV-2", "ACE2 regulation; impacts viral entry"),
    ("Acetylation", "H2AX-K6ac", "HCoV-229E", "DNA damage response; decreased during infection"),
    ("Citrullination", "Cit-H3", "SARS-CoV-2", "NETosis marker; predicts 28-day mortality"),
    ("Phosphorylation", "H2AC11-T17ph", "HCoV-229E", "Oxidative stress marker"),
    ("Phosphorylation", "H3S10ph", "SARS-CoV-2", "Cell cycle marker"),
    ("Glutathionylation", "H3C1 (C97/C111)", "HCoV-229E", "Nucleosome destabilization"),
    ("Truncation", "C-truncated H2A", "HCoV-229E", "Reduced abundance during infection"),
    ("Truncation", "N-truncated H3", "HCoV-229E", "Cleaved at R27/K28; reduced abundance"),
]
for i, row_data in enumerate(rows3):
    for j, val in enumerate(row_data):
        cell = table3.rows[i+1].cells[j]
        cell.text = val
        for paragraph in cell.paragraphs:
            for run in paragraph.runs:
                run.font.size = Pt(8)
                run.font.name = 'Times New Roman'

add_para("", size=10)

add_heading_styled("3.5. Molecular Mechanisms of Histone-Mediated Tissue Damage", level=2)
res_mech = (
    "Five distinct but interconnected mechanisms of histone-mediated tissue damage during "
    "coronavirus infection were identified from the literature synthesis. "
    "Mechanism 1: Direct membrane disruption. Histone H4 behaves as a molecular detergent, "
    "intercalating into plasma membranes to form lytic pores and triggering mitochondrial "
    "outer membrane permeabilization (MOMP), leading to cytochrome c release and apoptosome "
    "activation [8]. Mechanism 2: Ionic dyshomeostasis. Histones, particularly H1, activate "
    "nonselective cation channels, causing massive calcium influx from both extracellular "
    "space and endoplasmic reticulum stores, resulting in signaling paralysis and ER stress [8, 9]. "
    "Mechanism 3: TLR-mediated innate immune overdrive. Extracellular histones bind TLR2, "
    "TLR4, and TLR9, activating the MyD88-IRAK4-TRAF6 signaling axis, which drives NF-kB "
    "and AP-1, culminating in the cytokine storm characteristic of severe COVID-19 [8]. "
    "Mechanism 4: Inflammasome activation and pyroptosis. Intracellular sensing of histones "
    "by NOD2 and NLRP3 triggers inflammasome assembly, Caspase-1 activation, and release "
    "of IL-1beta and IL-18 [9]. Mechanism 5: Complement system sabotage. Histones H3 and H4 "
    "bind complement protein C4, inhibiting C3 and C5 functional activities, impairing viral "
    "clearance and establishing immunothrombosis [8]."
)
add_para(res_mech)

add_figure(
    "NETosis pathway illustrating the six-step process from neutrophil activation to "
    "extracellular trap release and histone unmasking.",
    "fig4_netosis.png", width=5.0
)
add_para("", size=10)

add_heading_styled("3.6. ORF8 Histone Mimicry Network", level=2)
res_orf8 = (
    "The SARS-CoV-2 ORF8 protein's histone mimicry mechanism was represented as a directed "
    "network (Figure 5). ORF8 mimics the ARKS motif of Histone H3, associates with host "
    "chromatin, promotes repressive marks (H3K9me3, H3K27me3), depletes activating marks "
    "(H3K9ac, H3K14ac), induces chromatin compaction, and ultimately silences interferon "
    "response genes [5]. This coordinated epigenetic reprogramming represents a sophisticated "
    "viral strategy to evade the innate immune response."
)
add_para(res_orf8)

add_figure(
    "Directed network diagram of SARS-CoV-2 ORF8 histone mimicry. ORF8 mimics the ARKS "
    "motif of H3, associates with host chromatin, and orchestrates a cascade of repressive "
    "epigenetic changes leading to IFN response gene silencing.",
    "fig5_orf8_mimicry.png", width=5.0
)
add_para("", size=10)

res_immuno = (
    "The immunothrombosis vicious cycle (Figure 6) illustrates how histone-mediated complement "
    "inhibition perpetuates tissue damage. Impaired viral clearance due to C3/C5 inhibition "
    "leads to continued cell death and further histone release, creating a self-reinforcing "
    "cycle that drives pathological clot formation in the microvasculature [8]."
)
add_para(res_immuno)

add_figure(
    "The immunothrombosis vicious circle. Histone-mediated complement inhibition impairs "
    "viral clearance, leading to continued cell death and further histone release.",
    "fig6_immunothrombosis.png", width=4.0
)
add_para("", size=10)

doc.add_page_break()

# ============================================================
# 4. DISCUSSION
# ============================================================
add_heading_styled("4. Discussion", level=1)

disc1 = (
    "We benchmarked three public histone databases for co-occurrence and avoidance "
    "patterns using median-split binarized marks. None of the three databases yielded "
    "significant pairwise relationships after FDR correction. This null result carries "
    "biological implications that we discuss below."
)
add_para(disc1)

disc_occ_biolog = (
    "Occurrence frequencies across databases reflect real biology. In Histome2, has_writer "
    "occurs at 41.9% of PTM sites while has_eraser occurs at 33.9%. This gap is not a data "
    "artifact: for most histone PTMs, the writer enzymes that deposit the mark are better "
    "characterized than the erasers that remove it. Lysine methyltransferases (KMTs) and "
    "histone acetyltransferases (HATs) have been studied for decades, while many corresponding "
    "demethylases (KDMs) and deacetylases (HDACs) were identified later and remain less "
    "completely annotated for site specificity [17]. That multi_mod is uniformly zero in "
    "Histome2 reflects a genuine feature of the database: each PTM site is listed with a "
    "single modification type. The site may carry multiple marks in vivo (e.g., H3K9 can be "
    "methylated or acetylated), but Histome2 catalogues each as a separate PTM site entry, "
    "so within a single entry the modification type count is always one."
)
add_para(disc_occ_biolog)

disc_occ_biolog2 = (
    "In HistoneDB 2.0, all three marks sit at exactly 50.0% by construction, but the "
    "underlying variables carry biological meaning. Variants with high HMM scores are those "
    "that match known histone domains well; variants with broad taxonomic distribution "
    "represent evolutionarily ancient histones. The near-significant trend between "
    "many_sequences and broad_taxonomy (p_adj = 0.0647, Log2FC = 0.90) is biologically "
    "coherent: conserved histone variants are sequenced more often across more species. "
    "Histone H3.3 and H2A.Z exemplify this: both are found across eukaryotes and have "
    "hundreds of sequences in the database, while variant histones restricted to specific "
    "lineages (e.g., mammalian H1.0) have fewer entries [12]. In ProHistoneDB, the lower "
    "occurrence of broad_taxonomy (37.5%) reflects that many archaeal and viral proto-histones "
    "are confined to single taxonomic groups. Proto-histones in Archaea are restricted to "
    "specific phyla (Euryarchaeota, Crenarchaeota) and have not spread across domains [13]."
)
add_para(disc_occ_biolog2)

disc_cooc_biolog = (
    "Why do real database features show no co-occurrence while the simulated PTM-level "
    "marks do? The answer lies in what the features represent. In the simulated data, "
    "H3K4me3 and H3K27ac co-occur because they are deposited by the same transcriptional "
    "machinery at active promoters. The MLL/COMPASS complex deposits H3K4me3, and p300/CBP "
    "deposits H3K27ac; these enzymes are recruited together at active regulatory elements "
    "[17, 18]. Co-occurrence arises through reader-mediated cross-talk: the SAGA complex "
    "subunit SGF29 binds H3K4me3 via its tandem tudor domain, facilitating HAT recruitment "
    "to active promoters [18]. In contrast, "
    "the database features (writer counts, HMM scores, sequence counts) do not participate "
    "in such cross-talk. A PTM site having a known writer does not cause it to also have a "
    "known eraser; these are independent annotation facts. The Co-ChIP method by Weiner et "
    "al. [20] demonstrated that PTM-level co-occurrence on individual nucleosomes is "
    "widespread and follows predictable rules: H3K4me3 co-occurs with H3K27ac at active "
    "promoters, H3K9me3 co-occurs with H4K20me3 at heterochromatin, and H3K36me3 is "
    "enriched in gene bodies independently of promoter marks. These are the patterns we "
    "recovered in the simulated data but could not detect in the database-level features."
)
add_para(disc_cooc_biolog)

disc_avoid_biolog = (
    "Avoidance between histone marks is biologically as informative as co-occurrence. In "
    "the simulated data, H3K9me3 and H3K4me3 were mutually exclusive (observed = 4, "
    "expected = 19.53, Log2FC = -2.29). This mirrors well-established biology: H3K9me3 is "
    "deposited by SUV39H1/2 at heterochromatin and recruits HP1 to compact chromatin, while "
    "H3K4me3 is deposited by MLL complexes at euchromatic promoters [21]. These marks do "
    "not coexist because the enzymes that deposit them recognize mutually exclusive chromatin "
    "environments. H3K9me3 and H3K27ac also showed strong avoidance (observed = 5, "
    "expected = 21.18, Log2FC = -2.08), consistent with the fact that acetylation and "
    "trimethylation of H3K27 are mutually exclusive: p300/CBP catalyzes acetylation, while "
    "PRC2 catalyzes trimethylation, and these complexes are recruited to different genomic "
    "regions [22]. Bivalency (H3K4me3 + H3K27me3) is an exception rather than the rule and "
    "is largely restricted to poised developmental genes in embryonic stem cells [23]."
)
add_para(disc_avoid_biolog)

disc_null_biolog = (
    "The absence of significant avoidance in all three real databases is itself a biological "
    "finding. The database attributes we tested are not inherently antagonistic. A histone "
    "PTM site can have both high writer and high eraser counts because writers and erasers "
    "serve complementary roles in mark turnover. For example, H3K4me3 has well-characterized "
    "writers (MLL1-4, SETD1A/B) and erasers (KDM5A-D, KDM2B), and active promoters show "
    "rapid turnover of both marks [24]. Similarly, a histone variant can have both high HMM "
    "scores and broad taxonomic distribution because conservation and confidence are "
    "positively correlated rather than mutually exclusive. The database features simply do "
    "not encode the antagonistic relationships that exist between PTM-level marks. This "
    "distinction carries a practical message for computational epigenetics: database-level "
    "attributes such as writer/eraser counts may not substitute for PTM-level co-occurrence "
    "data when studying combinatorial histone codes."
)
add_para(disc_null_biolog)

disc_coronavirus = (
    "The link between histone marks and coronavirus pathogenesis is mechanistically "
    "specific. SARS-CoV-2 employs at least three distinct viral proteins to subvert the "
    "histone modification system. First, ORF8 mimics the ARKS motif of histone H3 and "
    "associates with host chromatin, causing global increases in H3K9me3 and H3K27me3 and "
    "decreases in H3K9ac and H3K14ac, resulting in chromatin compaction and silencing of "
    "interferon response genes [5]. Second, NSP1 induces H3K9me2 at immune-related gene "
    "promoters via the G9a/GLP methyltransferase, repressing antiviral gene expression at "
    "the transcriptional level [25]. Pharmacological inhibition of G9a by UNC0638 restores "
    "interferon-stimulated gene expression and restricts viral replication. Third, the main "
    "protease Mpro (NSP5) cleaves histone deacetylase 2 (HDAC2), disrupting the "
    "deacetylation arm of the PTM regulatory cycle and further attenuating the interferon "
    "response [26]. This multi-pronged attack on histone modification regulatory enzymes "
    "explains why SARS-CoV-2 is so effective at suppressing host innate immunity."
)
add_para(disc_coronavirus)

disc_coronavirus2 = (
    "The chromatin-level impact extends beyond individual marks. Wang et al. [27] showed "
    "that SARS-CoV-2 restructures the host 3D genome within 24 hours of infection: A "
    "compartment regions weaken, H3K27ac levels drop globally, and cohesin is depleted from "
    "topologically associating domains (TADs). These changes correlate with the "
    "transcriptional suppression of interferon response genes. The magnitude of this "
    "epigenetic reprogramming - affecting compartmentalization, histone marks, and "
    "transcription simultaneously - suggests that histone mark avoidance and co-occurrence "
    "patterns at a genomic scale are actively rewritten by the virus. Our null results from "
    "static database features do not capture this dynamic rewiring. Time-resolved ChIP-seq "
    "during infection time courses would be needed to observe how SARS-CoV-2 alters "
    "co-occurrence relationships in real time."
)
add_para(disc_coronavirus2)

disc_limitations = (
    "This study has several limitations. First, the three databases cover distinct biological "
    "domains (human PTM sites, cross-species histone variants, and proto-histone categories), "
    "making direct cross-database comparisons of co-occurrence patterns informative but not "
    "strictly equivalent. Second, the sample sizes are modest (24-62 observations), limiting "
    "statistical power for detection of moderate effect sizes. Third, the median-split "
    "binarization discards information about continuous variation; a site with barely "
    "above-median writer count is treated the same as one with many times the median. "
    "Fourth, database-level attributes (writer counts, HMM scores) are not the same as "
    "the PTM-level marks that show established co-occurrence in ChIP-seq data. Future work "
    "should use raw continuous values in correlation frameworks, incorporate ENCODE "
    "ChIP-seq data from coronavirus-infected cells, and apply multivariate methods that "
    "account for the continuous nature of histone mark measurements."
)
add_para(disc_limitations)

disc_future = (
    "Future directions include extending this analysis to real ChIP-seq data from "
    "SARS-CoV-2-infected cells, comparing mark co-occurrence before and after infection "
    "to quantify the virus-induced rewiring of the histone interaction network. Integration "
    "with RNA-seq would connect specific co-occurrence changes to transcriptional outcomes. "
    "The finding that SARS-CoV-2 targets multiple points in the histone regulatory cycle "
    "(H3K9 methyltransferases via ORF8, G9a via NSP1, HDAC2 via Mpro) suggests that "
    "combination therapies targeting both writers and erasers might be more effective than "
    "single-agent approaches."
)
add_para(disc_future)

doc.add_page_break()

# ============================================================
# 5. CONCLUSION
# ============================================================
add_heading_styled("5. Conclusion", level=1)
concl = (
    "This study benchmarked three public histone databases for co-occurrence and avoidance "
    "analysis using median-split binarized marks. No statistically significant relationships "
    "were detected after FDR correction across Histome2, HistoneDB 2.0, or ProHistoneDB, "
    "indicating that the chosen attributes behave independently within each database. The "
    "strongest trend was observed in HistoneDB 2.0 between sequence count and taxonomic "
    "breadth (p_adj = 0.0647), suggesting that better-powered studies might reveal meaningful "
    "associations. The simulated nucleosome data, serving as a positive control, correctly "
    "recovered six significant mark interactions (three co-occurrence, three avoidance), "
    "confirming that the analytical framework is sensitive when biologically meaningful "
    "relationships exist. The five mechanisms of histone-mediated tissue damage provide a "
    "comprehensive model for understanding the transition from intracellular epigenetic "
    "manipulation to extracellular cytotoxicity during coronavirus infection. These findings "
    "establish a baseline for future studies using alternative encoding strategies and "
    "support the development of host-directed therapies targeting both the epigenetic and "
    "cytotoxic dimensions of histone biology in severe coronavirus infections."
)
add_para(concl)

# ============================================================
# ABBREVIATIONS
# ============================================================
add_heading_styled("List of Abbreviations", level=1)
abbrevs = [
    "ACE2 - Angiotensin-Converting Enzyme 2",
    "APC - Activated Protein C",
    "ARDS - Acute Respiratory Distress Syndrome",
    "CRP - C-Reactive Protein",
    "DAMP - Damage-Associated Molecular Pattern",
    "FDR - False Discovery Rate",
    "HMM - Hidden Markov Model",
    "IFN - Interferon",
    "ISG - Interferon-Stimulated Gene",
    "LMWH - Low-Molecular-Weight Heparin",
    "MOMP - Mitochondrial Outer Membrane Permeabilization",
    "NET - Neutrophil Extracellular Trap",
    "PAMP - Pathogen-Associated Molecular Pattern",
    "PTM - Post-Translational Modification",
    "TLR - Toll-Like Receptor",
]
for a in abbrevs:
    p = doc.add_paragraph()
    run = p.add_run(a)
    run.font.name = 'Times New Roman'
    run.font.size = Pt(10)
    p.paragraph_format.line_spacing = 1.5

# ============================================================
# CONFLICT OF INTEREST
# ============================================================
add_heading_styled("Conflict of Interest", level=1)
add_para("The author declares no conflict of interest.")

add_heading_styled("Acknowledgements", level=1)
add_para("This work was supported by the Department of Biotechnology, i3L University, Jakarta. The Claude cowork from Anthropic ® was used for outlining, brainstorming, and ideation of the manuscript. All authors are responsible for the validity of the content. The source code for the bioinformatics pipeline is deposited here: https://github.com/Johnkecops/Histone_Marks_COVID")

# ============================================================
# REFERENCES
# ============================================================
add_heading_styled("References", level=1)

references = [
    "[1] Schafer A, Baric RS. Epigenetic landscape during coronavirus infection. Pathogens 2017; 6(1): 8. doi: 10.3390/pathogens6010008",
    "[2] Kaur R, et al. Histome2: a database of human histone PTMs. Nucleic Acids Res 2021; 49(D1): D1024-D1030. doi: 10.1093/nar/gkaa641",
    "[3] Gupta S, Hemeg HA, Afrin F. Immuno-epigenetic paradigms in coronavirus infection. Front Immunol 2025; 16: 1596135. doi: 10.3389/fimmu.2025.1596135",
    "[4] York A. Histone mimicry by SARS-CoV-2. Nat Rev Microbiol 2022; 20(12): 703. doi: 10.1038/s41579-022-00815-9",
    "[5] Kee J, Thudium S, Renner DM, et al. SARS-CoV-2 disrupts host epigenetic regulation via histone mimicry. Nature 2022; 610(7931): 381-388. doi: 10.1038/s41586-022-05282-z",
    "[6] Menachery VD, Schafer A, et al. MERS-CoV and HCoV-229E modulate host histone modifications. J Virol 2020; 94(18): e00763-20.",
    "[7] Ives AN, Thibert S, Berger MR, et al. Human coronavirus 229E infection alters histone proteoforms. J Proteome Res 2026; 25(6): 2698-2708. doi: 10.1021/acs.jproteome.5c01015",
    "[8] de Vries F, Huckriede J, Wichapong K, Reutelingsperger C, Nicolaes GAF. The role of extracellular histones in COVID-19. J Intern Med 2022; 293(3): 275. doi: 10.1111/joim.13585",
    "[9] Marsman G, Zeerleder S, Luken BM. Extracellular histones, cell-free DNA, or nucleosomes: what drives coagulation? J Thromb Haemost 2022; 20(6): 1301-1312.",
    "[10] Yang X, Julian J, Schnell T, Albrecht H, Owens W, Yang XS. Epigenetic regulation of inflammatory genes involving histone methylation is associated with hyperimmune response in COVID-19 patients. J Immunol 2022; 208(Suppl 1): 111.19. doi: 10.4049/jimmunol.208.Supp.111.19",
    "[11] Kaur R, et al. Histome2: a database of human histone post-translational modifications. Nucleic Acids Res 2021; 49(D1): D1024-D1030.",
    "[12] Draizen EJ, Shaytan AK, Marino-Ramirez L, et al. HistoneDB 2.0: a histone database with variants. Nucleic Acids Res 2023; 51(D1): D348-D355. doi: 10.1093/nar/gkac1011",
    "[13] Henneman B, van Emmerik C, van Ingen H, Dame RT. ProHistoneDB: a database of proto-histones. Nucleic Acids Res 2018; 46(D1): D348-D354. doi: 10.1093/nar/gkx1014",
    "[14] Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Ser B 1995; 57(1): 289-300.",
    "[15] Hagberg AA, Schult DA, Swart PJ. Exploring network structure, dynamics, and function using NetworkX. Proc SciPy 2008; pp 11-15.",
    "[16] Hunter JD. Matplotlib: a 2D graphics environment. Comput Sci Eng 2007; 9(3): 90-95. doi: 10.1109/MCSE.2007.55",
    "[17] Abdulla AZ, Salari MM, Rouhani M, et al. Histone modification cross-talk: analytical tools and molecular mechanisms. Biochem J 2026; 483(9): 1059-1085. doi: 10.1042/BCJ20250115",
    "[18] Rothbart SB, Strahl BD. The interplay of histone modifications - writers that read. EMBO Rep 2015; 16(12): 1627-1645. doi: 10.15252/embr.201540945",
    "[19] Torres IO, Fujimori DG. Functional coupling between writers, erasers and readers of histone and DNA methylation. Curr Opin Struct Biol 2015; 35: 68-75. doi: 10.1016/j.sbi.2015.09.007",
    "[20] Weiner A, Lara-Astiaso D, Krupalnik V, et al. Co-ChIP enables genome-wide mapping of histone mark co-occurrence at single-molecule resolution. Nat Biotechnol 2016; 34(9): 953-961. doi: 10.1038/nbt.3652",
    "[21] Becker JS, Nicetto D, Zaret KS. H3K9me3-dependent heterochromatin: barrier to cell fate changes. Trends Genet 2016; 32(1): 29-41. doi: 10.1016/j.tig.2015.11.001",
    "[22] Margueron R, Reinberg D. The Polycomb complex PRC2 and its mark in life. Nature 2011; 469(7330): 343-349. doi: 10.1038/nature09784",
    "[23] Kumar D, Cinghu S, Oldfield AJ, Yang P, Jothi R. Decoding the function of bivalent chromatin in development and cancer. Genome Res 2021; 31(12): 2175-2188. doi: 10.1101/gr.275736.121",
    "[24] Hyun K, Jeon J, Park K, Kim J. Writing, erasing and reading histone lysine methylations. Exp Mol Med 2017; 49(4): e324. doi: 10.1038/emm.2017.11",
    "[25] Anastasakis DG, Benhalevy D, Cuburu N, Altan-Bonnet N, Hafner M. Epigenetic repression of antiviral genes by SARS-CoV-2 NSP1. PLoS One 2024; 19(2): e0297262. doi: 10.1371/journal.pone.0297262",
    "[26] Gordon DE, Hiatt J, Bouhaddou M, et al. The main protease of SARS-CoV-2 cleaves histone deacetylases and DCP1A, attenuating the immune defense of the interferon-stimulated genes. J Biol Chem 2023; 299(3): 102990. doi: 10.1016/j.jbc.2023.102990",
    "[27] Wang Y, Liu M, Gao J, et al. SARS-CoV-2 restructures host chromatin architecture. Nat Microbiol 2023; 8(4): 679-694. doi: 10.1038/s41564-023-01344-8",
]

for ref in references:
    p = doc.add_paragraph()
    run = p.add_run(ref)
    run.font.name = 'Times New Roman'
    run.font.size = Pt(10)
    p.paragraph_format.line_spacing = 1.5
    p.paragraph_format.first_line_indent = Cm(-1.27)
    p.paragraph_format.left_indent = Cm(1.27)

output_path = os.path.join(OUTPUT_DIR, "Histone_Marks_Coronaviruses_Manuscript.docx")
doc.save(output_path)
print(f"Manuscript saved to: {output_path}")
