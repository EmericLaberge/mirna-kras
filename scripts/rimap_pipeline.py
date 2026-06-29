#!/usr/bin/env python3
"""
RIMap-RISC Query Pipeline
Queries RIMap-RISC API for all transcript × microRNA combinations
Stores results in DuckDB for analysis

Usage:
    source venv/bin/activate
    python rimap_pipeline.py
"""

import requests
import duckdb
import json
import time
import re
import os
from pathlib import Path
from datetime import datetime
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

API_URL = "https://rimap-risc.api.major.iric.ca/api/data"
# DB_PATH can be overridden via env var to keep the production DB safe
# (e.g. `DB_PATH=/tmp/test/rimap_results.duckdb python scripts/rimap_pipeline.py`)
import os
DB_PATH = os.environ.get("DB_PATH", "rimap_results.duckdb")
CHECKPOINT_FILE = "pipeline_checkpoint.txt"
REQUEST_DELAY = 0  # no delay, let server response time be the rate limit
MAX_WORKERS = 4  # parallel API calls (conservative)


def get_mirna_name_from_fasta():
    """Parse mirnas.fa to get MIMAT ID -> name mapping"""
    mirna_map = {}
    with open("mirnas.fa", "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                match = re.match(r">(\S+)\s+miRNA\s+\(([^)]+)\)", line)
                if match:
                    mirna_map[match.group(1)] = match.group(2)
    return mirna_map


def query_rimap(transcript_name, mirna_name):
    """Query RIMap-RISC API for one transcript-miRNA pair"""
    payload = {
        "gencode": "v46",
        "transcriptName": transcript_name,
        "miRName": mirna_name,
        "minSeed": "6mer",
        "UTR5": False,
        "CDS": True,
        "UTR3": True,
        "supplementaryPaired": False,
        "seedAccessible": 0.0,
        "supplementaryAccessible": 0.0,
        "conserved": 0
    }
    try:
        resp = requests.post(API_URL, json=payload, timeout=60)
        if resp.status_code == 200:
            return resp.json()
        else:
            return {"error": f"HTTP {resp.status_code}", "text": resp.text[:200]}
    except Exception as e:
        return {"error": str(e)}


def parse_interactions(binding_landscape):
    """Parse all interactions from binding landscape text"""
    interactions = []
    blocks = re.split(r"(?=5'-)", binding_landscape.strip())

    for block in blocks:
        if not block.strip():
            continue

        lines = block.strip().split('\n')
        if len(lines) < 2:
            continue

        line0 = lines[0]
        pos_match = re.search(r'\[(\d+)\.\.(\d+)\]', line0)
        region_match = re.search(r'\((\w+)\)', line0)
        seed_match = re.search(r'target site:\s*(\S+)', line0)

        line1 = lines[1]
        dg_match = re.search(r'dG:\s*([-\d.]+)', line1)
        kd_match = re.search(r'Kd:\s*([\d,]+)', line1)
        phastcons_match = re.search(r'phastcons:\s*([-\d.]+)', line1)
        bridge_match = re.search(r'bridge:\s*(\d+)', line1)
        validated = 'TarBase' in block

        if pos_match and dg_match:
            interactions.append({
                'position_start': int(pos_match.group(1)),
                'position_end': int(pos_match.group(2)),
                'region': region_match.group(1) if region_match else None,
                'seed_type': seed_match.group(1) if seed_match else None,
                'binding_dG': float(dg_match.group(1)),
                'Kd': float(kd_match.group(1).replace(',', '')) if kd_match else None,
                'phastcons': float(phastcons_match.group(1)) if phastcons_match else None,
                'bridge': int(bridge_match.group(1)) if bridge_match else None,
                'validated': validated,
                'interaction_string': block.strip()
            })

    return interactions


def process_rimap_result(transcript_id, mirna_id, result):
    """Process RIMap-RISC result and return list of interactions"""
    if "error" in result:
        return []

    binding_landscape = result.get("bindingLandscape", "")
    if not binding_landscape or binding_landscape == "No interactions found":
        return []

    interactions = parse_interactions(binding_landscape)
    for inter in interactions:
        inter["transcript_id"] = transcript_id
        inter["mirna_id"] = mirna_id

    return interactions


def init_database():
    """Initialize DuckDB with schema"""
    conn = duckdb.connect(DB_PATH)

    conn.execute("""
        CREATE TABLE IF NOT EXISTS transcripts (
            transcript_id TEXT PRIMARY KEY,
            gene_name TEXT,
            category TEXT,
            enst_id TEXT
        )
    """)

    conn.execute("""
        CREATE TABLE IF NOT EXISTS mirnas (
            mirna_id TEXT PRIMARY KEY,
            name TEXT
        )
    """)

    conn.execute("""
        CREATE TABLE IF NOT EXISTS results (
            transcript_id TEXT,
            mirna_id TEXT,
            binding_dG REAL,
            Kd REAL,
            position_start BIGINT,
            position_end BIGINT,
            region TEXT,
            seed_type TEXT,
            phastcons REAL,
            bridge INTEGER,
            validated BOOLEAN,
            interaction_string TEXT,
            PRIMARY KEY (transcript_id, mirna_id, position_start)
        )
    """)

    return conn


def load_transcripts(conn):
    """Load transcript metadata"""
    transcript_data = [
        ("ENST00000256078", "KRAS", "kras_pathway"),
        ("ENST00000275493", "EGFR", "kras_pathway"),
        ("ENST00000318493", "MET", "kras_pathway"),
        ("ENST00000259523", "MYC", "kras_pathway"),
        ("ENST00000251849", "RAF1", "kras_pathway"),
        ("ENST00000331710", "TBK1", "kras_pathway"),
        ("ENST00000315204", "STK33", "kras_pathway"),
        ("ENST00000316629", "CDK1", "kras_pathway"),
        ("ENST00000227507", "CCND1", "kras_pathway"),
        ("ENST00000366790", "PARP1", "kras_pathway"),
        ("ENST00000377687", "KLF5", "kras_pathway"),
        ("ENST00000269305", "TP53", "tumor_suppressors"),
        ("ENST00000267163", "RB1", "tumor_suppressors"),
        ("ENST00000371953", "PTEN", "tumor_suppressors"),
        ("ENST00000307102", "MAP2K1", "mapk_pathway"),
        ("ENST00000262948", "MAP2K2", "mapk_pathway"),
        ("ENST00000263025", "MAPK3", "mapk_pathway"),
        ("ENST00000215832", "MAPK1", "mapk_pathway"),
        ("ENST00000263967", "PIK3CA", "mapk_pathway"),
        ("ENST00000349310", "AKT1", "mapk_pathway"),
        ("ENST00000311278", "AKT2", "mapk_pathway"),
        ("ENST00000263826", "AKT3", "mapk_pathway"),
        ("ENST00000361445", "MTOR", "mapk_pathway"),
        ("ENST00000226574", "NFKB1", "mapk_pathway"),
        ("ENST00000308639", "RELA", "mapk_pathway"),
    ]

    for tid, gene, cat in transcript_data:
        conn.execute(
            "INSERT OR IGNORE INTO transcripts VALUES (?, ?, ?, ?)",
            [tid, gene, cat, f"{gene}-201"]
        )


def load_mirnas(conn):
    """Load microRNA metadata"""
    mirna_map = get_mirna_name_from_fasta()
    for mid, name in mirna_map.items():
        conn.execute(
            "INSERT OR IGNORE INTO mirnas VALUES (?, ?)",
            [mid, name]
        )


def query_worker(args):
    """Worker function for parallel queries"""
    transcript_name, transcript_id, mirna_id, mirna_name = args
    time.sleep(REQUEST_DELAY)  # Rate limiting
    result = query_rimap(transcript_name, mirna_name)
    interactions = process_rimap_result(transcript_id, mirna_id, result)
    return interactions


def run_pipeline():
    """Main pipeline execution with tqdm progress and parallel workers"""
    print("=" * 60)
    print("RIMap-RISC Analysis Pipeline")
    print("=" * 60)
    print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Parallel workers: {MAX_WORKERS}")
    print(f"Estimated time: ~3 hours for 65,000 queries")
    print("=" * 60)

    # Clean start
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)
    if os.path.exists(CHECKPOINT_FILE):
        os.remove(CHECKPOINT_FILE)

    print(f"\nInitializing database...")
    conn = init_database()

    print("Loading transcripts...")
    load_transcripts(conn)

    print("Loading microRNAs...")
    load_mirnas(conn)

    # Get all transcripts and mirnas
    transcripts = conn.execute("SELECT transcript_id, gene_name FROM transcripts").fetchall()
    mirnas = conn.execute("SELECT mirna_id, name FROM mirnas").fetchall()

    # Create work items
    work_items = []
    for t_id, gene_name in transcripts:
        transcript_name = f"{gene_name}-201"
        for m_id, mirna_name in mirnas:
            work_items.append((transcript_name, t_id, m_id, mirna_name))

    total_combos = len(work_items)
    print(f"\nFound {len(transcripts)} transcripts and {len(mirnas)} microRNAs")
    print(f"Total combinations: {total_combos}")

    # Process with progress bar
    pbar = tqdm(total=total_combos, desc="Overall progress")
    total_results = 0

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(query_worker, item): item for item in work_items}

        for future in as_completed(futures):
            try:
                interactions = future.result()
                for inter in interactions:
                    try:
                        conn.execute("""
                            INSERT OR IGNORE INTO results VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                        """, [
                            inter["transcript_id"],
                            inter["mirna_id"],
                            inter.get("binding_dG"),
                            inter.get("Kd"),
                            inter["position_start"],
                            inter["position_end"],
                            inter["region"],
                            inter["seed_type"],
                            inter.get("phastcons"),
                            inter.get("bridge"),
                            inter["validated"],
                            inter.get("interaction_string")
                        ])
                        total_results += 1
                    except Exception:
                        pass  # Skip duplicates
            except Exception:
                pass  # Skip failed queries

            pbar.update(1)
            pbar.set_postfix({"results": total_results})

    pbar.close()

    print(f"\n\n{'=' * 60}")
    print("PIPELINE COMPLETE")
    print(f"{'=' * 60}")
    print(f"Total interactions stored: {total_results}")
    print(f"Completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Summary
    print("\n--- Summary by Pathway ---")
    summary = conn.execute("""
        SELECT
            t.category,
            COUNT(DISTINCT t.transcript_id) as transcripts,
            COUNT(*) as interactions,
            AVG(r.binding_dG) as avg_energy,
            MIN(r.binding_dG) as strongest
        FROM transcripts t
        LEFT JOIN results r ON t.transcript_id = r.transcript_id
        GROUP BY t.category
    """).fetchall()

    for row in summary:
        print(f"  {row[0]}: {row[1]} tx, {row[2]} hits, avg dG: {row[3]:.2f}, strongest: {row[4]:.2f}")

    conn.close()


if __name__ == "__main__":
    run_pipeline()
