#!/usr/bin/env python3
"""Fetch NCBI sequence reports (chromosome name, length, GenBank/RefSeq accessions)
for genome assemblies, one TSV per assembly in OUTDIR. Light HTTP only, so it is fine
on a login node. Assemblies already fetched are skipped.

Usage: scripts/fetch_sequence_reports.py OUTDIR ACCESSION [ACCESSION ...]
"""
import json
import os
import sys
import time
import urllib.parse
import urllib.request

API = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{}/sequence_reports"
COLS = ["chr_name", "genbank_accession", "refseq_accession", "length", "role",
        "assigned_molecule_location_type", "gc_percent"]


def fetch(acc):
    rows, token = [], None
    while True:
        q = {"page_size": 1000}
        if token:
            q["page_token"] = token
        url = API.format(acc) + "?" + urllib.parse.urlencode(q)
        for attempt in range(5):
            try:
                with urllib.request.urlopen(url, timeout=60) as r:
                    d = json.load(r)
                break
            except Exception:
                if attempt == 4:
                    raise
                time.sleep(2 ** attempt)
        rows += d.get("reports", [])
        token = d.get("next_page_token")
        if not token:
            return rows
        time.sleep(0.4)   # stay under NCBI's 3 requests/s without an API key


def main():
    out, accs = sys.argv[1], sys.argv[2:]
    os.makedirs(out, exist_ok=True)
    for acc in accs:
        f = os.path.join(out, acc + ".tsv")
        if os.path.exists(f):
            continue
        rows = fetch(acc)
        with open(f + ".part", "w") as fh:
            fh.write("\t".join(COLS) + "\n")
            for r in rows:
                fh.write("\t".join(str(r.get(c, "NA")) for c in COLS) + "\n")
        os.replace(f + ".part", f)
        print(acc, len(rows), "sequences", flush=True)
        time.sleep(0.4)


if __name__ == "__main__":
    main()
