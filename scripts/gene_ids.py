#!/usr/bin/env python3
"""
Gene/Transcript IDs from notes.md for mirna-kras project
Format: GENENAME-201 (symbol-version suffix for GenCode/Ensembl compatibility)

These are the genes mentioned in the project notes for microRNA-target coupling analysis.
"""

# KRAS pathway genes (from notes.md)
KRAS_PATHWAY_GENES = [
    ("KRAS", "ENSG00000133703"),
    ("EGFR", "ENSG00000146648"),
    ("MET", "ENSG00000105976"),
    ("MYC", "ENSG00000136997"),
    ("RAF1", "ENSG00000132155"),
    ("TBK1", "ENSG00000153540"),
    ("STK33", "ENSG00000104248"),
    ("CDK1", "ENSG00000170312"),
    ("CCND1", "ENSG00000110092"),
    ("PARP1", "ENSG00000138996"),
    ("KLF5", "ENSG00000103024"),
]

# Tumor suppressor genes
TUMOR_SUPPRESSOR_GENES = [
    ("TP53", "ENSG00000141510"),
    ("RB1", "ENSG00000139687"),
    ("PTEN", "ENSG00000197062"),
]

# MAPK pathway genes
MAPK_PATHWAY_GENES = [
    ("MAP2K1", "ENSG00000169032"),
    ("MAP2K2", "ENSG00000126934"),
    ("MAPK3", "ENSG00000102882"),
    ("MAPK1", "ENSG00000100030"),
    ("PIK3CA", "ENSG00000121858"),
    ("AKT1", "ENSG00000142208"),
    ("AKT2", "ENSG00000105952"),
    ("AKT3", "ENSG00000117020"),
    ("MTOR", "ENSG00000198746"),
    ("NFKB1", "ENSG00000109320"),
    ("RELA", "ENSG00000173039"),
]

# All categories combined
ALL_GENES = {
    "kras_pathway": KRAS_PATHWAY_GENES,
    "tumor_suppressors": TUMOR_SUPPRESSOR_GENES,
    "mapk_pathway": MAPK_PATHWAY_GENES,
}
