#!/usr/bin/env python3
"""
RIMap-RISC Local Query Pipeline
Runs RIMap calculation locally using Java MiRScanS
Stores results in DuckDB for analysis

Usage:
    source venv/bin/activate
    python rimap_pipeline_local.py

Requirements:
    - MC-Flashfold installed at /home/emeric/MC-Flashfold-v37.0/
    - GENCODE v46 transcripts FASTA (gencode.v46.transcripts.fa.gz)
    - Java classes compiled (mvn compile)
"""

import subprocess
import duckdb
import json
import os
import re
import gzip
from pathlib import Path
from datetime import datetime
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

# Paths - use absolute paths
SCRIPT_DIR = Path(__file__).parent.resolve()
# Override via env var (used by the Docker image to point at the in-container
# GENCODE archive path); default to a sibling of the script for host runs.
GENCODE_FASTA = os.environ.get(
    "GENCODE_FASTA",
    str(SCRIPT_DIR / "gencode.v46.transcripts.fa.gz"),
)
# Paths can be overridden via env vars (used by the Docker image to point at
# the in-container locations); default to the host paths used during development.
MCFLASHFOLD_PATH = os.environ.get("MCFLASHFOLD_PATH", "/home/emeric/MC-Flashfold-v37.0")
JAVA_CLASSPATH = os.environ.get(
    "JAVA_CLASSPATH",
    str(SCRIPT_DIR / "rimap_java_code" / "target" / "classes"),
)
JACKSON_JARS = [
    f"{os.path.expanduser('~')}/.m2/repository/com/fasterxml/jackson/core/jackson-databind/2.15.3/jackson-databind-2.15.3.jar",
    f"{os.path.expanduser('~')}/.m2/repository/com/fasterxml/jackson/core/jackson-core/2.15.3/jackson-core-2.15.3.jar",
    f"{os.path.expanduser('~')}/.m2/repository/com/fasterxml/jackson/core/jackson-annotations/2.15.3/jackson-annotations-2.15.3.jar",
]
# Override via env var (used by the Docker image to point the new DB at a
# path outside the production DB); default to a sibling of the script.
DB_PATH = os.environ.get("DB_PATH", str(SCRIPT_DIR / "rimap_results.duckdb"))
JAVA_MAIN = "ca.iric.major.rinexus.rimaprisc.MiRScanS"
MAX_WORKERS = 4  # parallel Java processes


def get_transcript_sequence(enst_id, gencode_fasta):
    """Extract transcript sequence from compressed GENCODE FASTA"""
    cmd = [
        "zcat", gencode_fasta
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    lines = result.stdout.split('\n')

    seq_lines = []
    capture = False
    for line in lines:
        if line.startswith(f'>{enst_id}'):
            capture = True
            continue
        if capture:
            if line.startswith('>'):
                break
            seq_lines.append(line.strip())
    return ''.join(seq_lines) if seq_lines else None


def get_mirna_sequence(mirna_name, mirna_map):
    """Get mature miRNA sequence from our mirna_map (mimatic -> sequence)"""
    # mirna_map is already loaded: mirna_id -> (name, sequence)
    for mimat, (name, seq) in mirna_map.items():
        if name == mirna_name or mimat == mirna_name:
            return mimat, name, seq
    return None, None, None


def get_mirna_name_from_fasta():
    """Parse mirnas.fa to get MIMAT ID -> (name, sequence) mapping"""
    mirna_map = {}
    current_mimat = None
    current_seq = None

    # mirnas.fa can live in several places: SCRIPT_DIR (host default), or
    # repo_root/data (Docker default, where the Dockerfile copies it).
    # Allow override via MIRNAS_FA env var.
    _candidates = [
        os.environ.get("MIRNAS_FA"),
        str(SCRIPT_DIR / "mirnas.fa"),
        str(SCRIPT_DIR.parent / "data" / "mirnas.fa"),
    ]
    _mirnas_path = next((p for p in _candidates if p and os.path.exists(p)), None)
    if not _mirnas_path:
        raise FileNotFoundError(
            "mirnas.fa not found. Tried: " + ", ".join(str(p) for p in _candidates if p)
        )
    with open(_mirnas_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                match = re.match(r">(\S+)\s+miRNA\s+\(([^)]+)\)", line)
                if match:
                    current_mimat = match.group(1)
                    current_name = match.group(2)
                    current_seq = None
                else:
                    current_mimat = None
            else:
                if current_mimat and current_seq is None:
                    current_seq = line.strip().replace('T', 'U')
                    mirna_map[current_mimat] = (current_name, current_seq)

    return mirna_map


def run_mirscan(transcript_name, transcript_id, mirna_id, mirna_name, mirna_seq, target_seq):
    """Run MiRScanS locally and return parsed JSON results"""
    env = os.environ.copy()
    env['PATH'] = f"{MCFLASHFOLD_PATH}:{env.get('PATH', '')}"

    classpath = JAVA_CLASSPATH + ":" + ":".join(JACKSON_JARS)

    cmd = [
        "java", "-cp", classpath, JAVA_MAIN,
        mirna_name,           # guide name
        mirna_id,             # guide MIMAT
        "Homo sapiens",       # organism
        mirna_seq,            # guide sequence
        transcript_name,       # target name (gene)
        transcript_id,        # target ENST
        target_seq            # target sequence
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120,
            env=env,
            cwd=MCFLASHFOLD_PATH,  # mcff looks for its tables/ in CWD
        )

        if result.returncode != 0:
            return []

        interactions = []
        for line in result.stdout.strip().split('\n'):
            line = line.strip()
            if not line or not line.startswith('{'):
                continue
            try:
                data = json.loads(line)
                inter = parse_mirscan_result(data, transcript_id, mirna_id)
                if inter:
                    interactions.append(inter)
            except json.JSONDecodeError:
                continue

        return interactions
    except subprocess.TimeoutExpired:
        return []
    except Exception as e:
        return []


def parse_mirscan_result(data, transcript_id, mirna_id):
    """Parse MiRScanS JSON output into interaction dict"""
    try:
        # Parse position: "pst" is start, "pen" is end
        position_start = data.get('pst')
        position_end = data.get('pen')

        if position_start is None or position_end is None:
            return None

        # Seed type from "sty" field, e.g., "5mer-m2.6", "A1-7mer-b(5.6)"
        seed_type = data.get('sty', '')

        # Region from "rgn" field
        region = data.get('rgn', '')

        # dG (free energy)
        binding_dG = data.get('sen')  # seed energy

        # Kd
        kd = data.get('kdv')

        # Target sequence for interaction_string
        tseq = data.get('tseq', '')
        alt = data.get('alt', '')
        alp = data.get('alp', '')
        alg = data.get('alg', '')

        # Build interaction string similar to API format
        interaction_string = f"[{position_start}..{position_end}] ({region}) target site: {seed_type}\n"
        interaction_string += f"dG: {binding_dG} Kd: {kd} phastcons: {data.get('pcs', 0)} bridge: {data.get('brl', 0)}\n"
        interaction_string += f"{tseq}\n{alt}\n{alp}\n{alg}"

        return {
            'transcript_id': transcript_id,
            'mirna_id': mirna_id,
            'position_start': position_start,
            'position_end': position_end,
            'region': region,
            'seed_type': seed_type,
            'binding_dG': binding_dG,
            'Kd': kd,
            'phastcons': data.get('pcs', 0),
            'bridge': data.get('brl', 0),
            'validated': False,  # local calculation doesn't have TarBase info
            'interaction_string': interaction_string
        }
    except Exception as e:
        return None


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
    for mid, (name, seq) in mirna_map.items():
        conn.execute(
            "INSERT OR IGNORE INTO mirnas VALUES (?, ?)",
            [mid, name]
        )


def worker(args):
    """Worker function for parallel local queries"""
    transcript_id, gene_name, mirna_id, mirna_name, mirna_seq, target_seq = args

    # Build transcript name (gene-201 format used by GENCODE)
    transcript_name = f"{gene_name}-201"

    interactions = run_mirscan(
        transcript_name, transcript_id,
        mirna_id, mirna_name, mirna_seq,
        target_seq
    )
    return interactions


def run_pipeline():
    """Main pipeline execution with tqdm progress and parallel workers"""
    print("=" * 60)
    print("RIMap-RISC Local Analysis Pipeline")
    print("=" * 60)
    print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Parallel workers: {MAX_WORKERS}")
    print(f"GENCODE FASTA: {GENCODE_FASTA}")
    print(f"MC-Flashfold: {MCFLASHFOLD_PATH}")
    print("=" * 60)

    # Clean start
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)

    print(f"\nInitializing database...")
    conn = init_database()

    print("Loading transcripts...")
    load_transcripts(conn)

    print("Loading microRNAs...")
    mirna_map = get_mirna_name_from_fasta()
    load_mirnas(conn)

    # Get all transcripts and mirnas
    transcripts = conn.execute("SELECT transcript_id, gene_name FROM transcripts").fetchall()
    mirnas = conn.execute("SELECT mirna_id, name FROM mirnas").fetchall()

    # Create mirna_id -> sequence map
    mirna_seqs = {}
    for mimat, (name, seq) in mirna_map.items():
        mirna_seqs[mimat] = (name, seq)

    total_combos = len(transcripts) * len(mirnas)
    print(f"\nFound {len(transcripts)} transcripts and {len(mirnas)} microRNAs")
    print(f"Total combinations: {total_combos}")

    # Extract all transcript sequences first (one-time cost)
    print("\nExtracting transcript sequences from GENCODE FASTA...")
    transcript_seqs = {}
    for t_id, gene_name in tqdm(transcripts, desc="Extracting transcripts"):
        seq = get_transcript_sequence(t_id, GENCODE_FASTA)
        if seq:
            transcript_seqs[t_id] = seq

    print(f"Extracted {len(transcript_seqs)} transcript sequences")

    # Create work items
    work_items = []
    for t_id, gene_name in transcripts:
        if t_id not in transcript_seqs:
            continue
        for m_id, m_name in mirnas:
            if m_id in mirna_seqs:
                name, seq = mirna_seqs[m_id]
                work_items.append((t_id, gene_name, m_id, name, seq, transcript_seqs[t_id]))

    print(f"Work items with valid sequences: {len(work_items)}")

    # Process with progress bar
    pbar = tqdm(total=len(work_items), desc="Overall progress")
    total_results = 0

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(worker, item): item for item in work_items}

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
        if row[3]:
            print(f"  {row[0]}: {row[1]} tx, {row[2]} hits, avg dG: {row[3]:.2f}, strongest: {row[4]:.2f}")
        else:
            print(f"  {row[0]}: {row[1]} tx, 0 hits")

    conn.close()


if __name__ == "__main__":
    run_pipeline()
