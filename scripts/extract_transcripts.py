#!/usr/bin/env python3
"""Extract specific transcript sequences from GENCODE v46 transcript file"""

import os
import sys
from pathlib import Path

# Correct ENST IDs from GENCODE v46 (verified from FTP file)
TRANSCRIPTS = {
    # KRAS pathway
    "kras_pathway": [
        "ENST00000256078",  # KRAS-201
        "ENST00000275493",  # EGFR-201
        "ENST00000318493",  # MET-201
        "ENST00000259523",  # MYC-201
        "ENST00000251849",  # RAF1-201
        "ENST00000331710",  # TBK1-201
        "ENST00000315204",  # STK33-201
        "ENST00000316629",  # CDK1-201
        "ENST00000227507",  # CCND1-201
        "ENST00000366790",  # PARP1-201
        "ENST00000377687",  # KLF5-201
    ],
    # Tumor suppressors
    "tumor_suppressors": [
        "ENST00000269305",  # TP53-201
        "ENST00000267163",  # RB1-201
        "ENST00000371953",  # PTEN-201
    ],
    # MAPK pathway
    "mapk_pathway": [
        "ENST00000307102",  # MAP2K1-201
        "ENST00000262948",  # MAP2K2-201
        "ENST00000263025",  # MAPK3-201
        "ENST00000215832",  # MAPK1-201
        "ENST00000263967",  # PIK3CA-201
        "ENST00000349310",  # AKT1-201
        "ENST00000311278",  # AKT2-201
        "ENST00000263826",  # AKT3-201
        "ENST00000361445",  # MTOR-201
        "ENST00000226574",  # NFKB1-201
        "ENST00000308639",  # RELA-201
    ],
}

FASTA_GZ = Path(__file__).parent / "gencode.v46.transcripts.fa.gz"


def extract_transcripts():
    """Extract each transcript to its category directory"""
    if not FASTA_GZ.exists():
        print(f"Error: {FASTA_GZ} not found", file=sys.stderr)
        sys.exit(1)

    import gzip

    base_dir = Path(__file__).parent
    current_transcript = None
    transcript_seq = []

    with gzip.open(FASTA_GZ, 'rt') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                # Save previous transcript if it matches
                if current_transcript and transcript_seq:
                    for category, transcripts in TRANSCRIPTS.items():
                        if current_transcript in transcripts:
                            output_dir = base_dir / category
                            output_dir.mkdir(exist_ok=True)
                            safe_id = current_transcript.split('.')[0]
                            output_file = output_dir / f"{safe_id}.fasta"
                            with open(output_file, 'w') as out:
                                out.write(f">{current_transcript}\n")
                                out.write(''.join(transcript_seq) + '\n')
                            print(f"Saved: {output_file}")
                            break

                # Parse header to get transcript ID
                header = line[1:]  # remove '>'
                parts = header.split('|')
                current_transcript = parts[0]  # ENST00000256078.10
                # Remove version for matching
                current_transcript_base = current_transcript.split('.')[0]
                transcript_seq = []
            else:
                transcript_seq.append(line)

        # Don't forget last transcript
        if current_transcript and transcript_seq:
            for category, transcripts in TRANSCRIPTS.items():
                if current_transcript_base in transcripts:
                    output_dir = base_dir / category
                    output_dir.mkdir(exist_ok=True)
                    safe_id = current_transcript_base
                    output_file = output_dir / f"{safe_id}.fasta"
                    with open(output_file, 'w') as out:
                        out.write(f">{current_transcript}\n")
                        out.write(''.join(transcript_seq) + '\n')
                    print(f"Saved: {output_file}")
                    break


if __name__ == "__main__":
    extract_transcripts()
