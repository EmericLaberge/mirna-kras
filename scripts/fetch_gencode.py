#!/usr/bin/env python3
"""Fetch transcript sequences from Ensembl REST API for genes in notes.md"""

import requests
import sys
import json
import os
from pathlib import Path

# Gene IDs from notes.md (Ensembl transcript IDs with -201 suffix)
# KRAS pathway
KRAS_PATHWAY = [
    "ENST00000256078",  # KRAS-201
    "ENST00000262020",  # EGFR-201
    "ENST00000389605",  # MET-201
    "ENST00000396095",  # MYC-201
    "ENST00000416153",  # RAF1-201
    "ENST00000342728",  # TBK1-201
    "ENST00000261954",  # STK33-201
    "ENST00000330567",  # CDK1-201
    "ENST00000303299",  # CCND1-201
    "ENST00000337655",  # PARP1-203
    "ENST00000237885",  # KLF5-201
]

# Tumor suppressors
TUMOR_SUPPRESSORS = [
    "ENST00000269305",  # TP53-201
    "ENST00000267157",  # RB1-201
    "ENST00000371953",  # PTEN-201
]

# MAPK pathway
MAPK_PATHWAY = [
    "ENST00000296788",  # MAP2K1-201
    "ENST00000398575",  # MAP2K2-201
    "ENST00000496030",  # MAPK3-201
    "ENST00000216448",  # MAPK1-201
    "ENST00000263967",  # PIK3CA-201
    "ENST00000281869",  # AKT1-201
    "ENST00000336047",  # AKT2-201
    "ENST00000345202",  # AKT3-201
    "ENST00000166971",  # MTOR-201
    "ENST00000628269",  # NFKB1-201
    "ENST00000164944",  # RELA-201
]

SERVER = "https://rest.ensembl.org"


def fetch_sequence(transcript_id, output_dir):
    """Fetch transcript sequence from Ensembl and save to FASTA file"""
    ext = f"/sequence/id/{transcript_id}"
    headers = {"Content-Type": "text/x-fasta"}
    params = {"species": "human"}
    r = requests.get(f"{SERVER}{ext}", headers=headers, params=params)

    if not r.ok:
        print(f"Error fetching {transcript_id}: {r.status_code} - {r.text[:200]}", file=sys.stderr)
        return False

    # Extract filename from transcript ID (remove version)
    safe_id = transcript_id.split('.')[0]
    output_file = output_dir / f"{safe_id}.fasta"

    with open(output_file, 'w') as f:
        f.write(r.text)

    print(f"Saved: {output_file}")
    return True


def fetch_gene_info(gene_id):
    """Fetch gene symbol and description from Ensembl"""
    ext = f"/lookup/id/{gene_id}"
    headers = {"Content-Type": "application/json"}
    r = requests.get(f"{SERVER}{ext}", headers=headers)

    if r.ok:
        return r.json()
    return None


def main():
    """CLI entry: fetch all KRas-pathway / tumor-suppressor / MAPK-pathway transcripts from Ensembl REST."""
    base_dir = Path(__file__).parent

    categories = {
        "kras_pathway": KRAS_PATHWAY,
        "tumor_suppressors": TUMOR_SUPPRESSORS,
        "mapk_pathway": MAPK_PATHWAY,
    }

    for category, transcripts in categories.items():
        output_dir = base_dir / category
        output_dir.mkdir(exist_ok=True)

        print(f"\n=== {category.upper()} ===")
        for transcript_id in transcripts:
            fetch_sequence(transcript_id, output_dir)


if __name__ == "__main__":
    main()
