#!/usr/bin/env python3

import argparse
import os
from pathlib import Path

from Bio import Entrez


def parse_args():
    parser = argparse.ArgumentParser(description="Fetch single GenBank FASTA by accession")
    parser.add_argument("--email", required=False, help="Entrez email, or set env NCBI_EMAIL")
    parser.add_argument("--api_key", required=False, help="NCBI API key, or set env NCBI_API_KEY")
    parser.add_argument("--acc", required=True)
    parser.add_argument("--out_fasta", required=True)
    return parser.parse_args()


def main():
    args = parse_args()

    email = args.email or os.getenv("NCBI_EMAIL")
    if not email:
        raise SystemExit("Set --email or env NCBI_EMAIL")
    Entrez.email = email

    api_key = args.api_key or os.getenv("NCBI_API_KEY")
    if api_key:
        Entrez.api_key = api_key

    output_path = Path(args.out_fasta)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    handle = Entrez.efetch(db="nucleotide", id=args.acc, rettype="fasta", retmode="text")
    with output_path.open("w", encoding="utf-8") as out_handle:
        out_handle.write(handle.read())
    handle.close()
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
