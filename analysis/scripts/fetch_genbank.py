#!/usr/bin/env python3

import argparse
import os
from pathlib import Path
import time

from Bio import Entrez


def parse_args():
    parser = argparse.ArgumentParser(description="Fetch GenBank FASTA by accession list")
    parser.add_argument("--email", required=False, help="Entrez email, or set env NCBI_EMAIL")
    parser.add_argument("--api_key", required=False, help="NCBI API key, or set env NCBI_API_KEY")
    parser.add_argument("--acc", required=True, help="Text file with one accession per line")
    parser.add_argument("--out_fasta", required=True)
    return parser.parse_args()


def load_accessions(path):
    with open(path, encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip() and not line.startswith("#")]


def main():
    args = parse_args()

    email = args.email or os.getenv("NCBI_EMAIL")
    if not email:
        raise SystemExit("Set --email or env NCBI_EMAIL to comply with NCBI usage policy")
    Entrez.email = email

    api_key = args.api_key or os.getenv("NCBI_API_KEY")
    if api_key:
        Entrez.api_key = api_key

    accessions = load_accessions(args.acc)
    output_path = Path(args.out_fasta)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as out_handle:
        for accession in accessions:
            handle = Entrez.efetch(
                db="nucleotide",
                id=accession,
                rettype="fasta",
                retmode="text",
            )
            out_handle.write(handle.read())
            handle.close()
            time.sleep(0.34)

    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
