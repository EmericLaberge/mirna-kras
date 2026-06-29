#!/usr/bin/env python3
"""
Smoke test for the LOCAL rimap pipeline inside the Docker image.

Queries 3 (transcript, miRNA) pairs via MiRScanS + MC-Flashfold, writes the
results to a new DuckDB at /tmp/work/rimap_results.duckdb. Does NOT touch
any production DB.

Runs inside the container (paths assume /app/repo and /tmp/work are mounted).
"""
from __future__ import annotations

import os
import sys
import time

import duckdb

# Repo scripts on sys.path
sys.path.insert(0, "/app/repo/scripts")
import rimap_pipeline_local as rp  # noqa: E402

DB_PATH = os.environ.get("DB_PATH", "/tmp/work/rimap_results.duckdb")
GENCODE_FASTA = os.environ.get("GENCODE_FASTA", "/tmp/work/gencode.v46.transcripts.fa.gz")
PROD_DB_PATH = "/tmp/prod.duckdb"

# A few hand-picked (transcript_id, mirna_name) pairs from KRAS-201
PAIRS = [
    ("ENST00000256078", "let-7a-5p"),
    ("ENST00000256078", "miR-21-5p"),
    ("ENST00000256078", "miR-34a-5p"),
]


def main() -> int:
    print(f"DB_PATH         = {DB_PATH}")
    print(f"GENCODE_FASTA   = {GENCODE_FASTA}")
    print(f"JAVA_CLASSPATH  = {os.environ.get('JAVA_CLASSPATH', '<unset>')}")
    print(f"MCFLASHFOLD_PATH= {os.environ.get('MCFLASHFOLD_PATH', '<unset>')}")
    print()

    print("=== Build new DB ===")
    rp.init_database()
    print(f"  created: {os.path.getsize(DB_PATH):,} bytes")

    con = duckdb.connect(DB_PATH)
    rp.load_transcripts(con)
    rp.load_mirnas(con)
    n_t = con.execute("SELECT COUNT(*) FROM transcripts").fetchone()[0]
    n_m = con.execute("SELECT COUNT(*) FROM mirnas").fetchone()[0]
    print(f"  loaded: {n_t} transcripts, {n_m} mirnas")
    mirna_map = rp.get_mirna_name_from_fasta()
    print(f"  mirna_map: {len(mirna_map)} entries")
    con.close()

    print()
    print("=== Run 3 LOCAL MiRScanS queries ===")
    con = duckdb.connect(DB_PATH)
    for enst, mirna_name in PAIRS:
        row = con.execute(
            "SELECT mirna_id FROM mirnas WHERE name = ? LIMIT 1", [mirna_name]
        ).fetchone()
        if not row:
            print(f"  SKIP {mirna_name}: not in mirnas table")
            continue
        mirna_id = row[0]
        if mirna_id not in mirna_map:
            print(f"  SKIP {mirna_name}: {mirna_id} not in mirna_map")
            continue
        _, mirna_seq = mirna_map[mirna_id]
        target_seq = rp.get_transcript_sequence(enst, GENCODE_FASTA)
        if not target_seq:
            print(f"  SKIP {enst}: not in GENCODE archive")
            continue
        print(f"  {enst} x {mirna_name}: ", end="", flush=True)
        t0 = time.time()
        rows = rp.run_mirscan("KRAS-201", enst, mirna_id, mirna_name, mirna_seq, target_seq) or []
        dt = time.time() - t0
        print(f"{len(rows)} sites ({dt:.2f}s)")
        for r in rows:
            con.execute(
                "INSERT OR IGNORE INTO results VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                [
                    r["transcript_id"], r["mirna_id"],
                    r.get("binding_dG"), r.get("Kd"),
                    r["position_start"], r["position_end"], r["region"],
                    r["seed_type"], r.get("phastcons"), r.get("bridge"),
                    r["validated"], r.get("interaction_string"),
                ],
            )
    con.close()

    print()
    print("=== New DB state ===")
    con = duckdb.connect(DB_PATH, read_only=True)
    for tbl in ("transcripts", "mirnas", "results"):
        n = con.execute(f"SELECT COUNT(*) FROM {tbl}").fetchone()[0]
        print(f"  {tbl:12s} {n:>8,} rows")
    print()
    print("=== Sample results rows ===")
    for row in con.execute(
        "SELECT transcript_id, mirna_id, position_start, binding_dG, Kd, region, seed_type "
        "FROM results LIMIT 5"
    ).fetchall():
        print(" ", row)
    con.close()

    # If production DB was mounted, compare counts
    if os.path.exists(PROD_DB_PATH):
        print()
        print(f"=== Compare to production DB ({PROD_DB_PATH}) ===")
        new_con = duckdb.connect(DB_PATH, read_only=True)
        prod_con = duckdb.connect(PROD_DB_PATH, read_only=True)
        for enst, mirna_name in PAIRS:
            row = new_con.execute(
                "SELECT mirna_id FROM mirnas WHERE name = ? LIMIT 1", [mirna_name]
            ).fetchone()
            if not row:
                continue
            mirna_id = row[0]
            new_n = new_con.execute(
                "SELECT COUNT(*) FROM results WHERE transcript_id=? AND mirna_id=?",
                [enst, mirna_id],
            ).fetchone()[0]
            prod_n = prod_con.execute(
                "SELECT COUNT(*) FROM results WHERE transcript_id=? AND mirna_id=?",
                [enst, mirna_id],
            ).fetchone()[0]
            ok = "YES" if new_n == prod_n else "NO"
            print(f"  {enst} x {mirna_name}:  new={new_n:>3d}  prod={prod_n:>3d}  {ok}")
        new_con.close()
        prod_con.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())
