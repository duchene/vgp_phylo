#!/usr/bin/awk -f
# Emit one TSV row per FASTA locus: locus_id, clade, ntips, nspecies, ncopies_max, nsites
# Accession = leading GCA_/GCF_ token including its version suffix (e.g. GCA_000001405.15).
# Sequence names look like <accession>.<chr>|copy<N>, and <chr> may itself contain dots.

function emit() {
    if (f == "") return
    ns = 0
    mx = 0
    for (k in sp) { ns++; if (sp[k] > mx) mx = sp[k] }
    base = f; sub(/.*\//, "", base); sub(/\.fa$/, "", base)
    clade = f; sub(/\/[^\/]*$/, "", clade); sub(/.*\//, "", clade)
    printf "%s\t%s\t%d\t%d\t%d\t%d\n", base, clade, n, ns, mx, w
}

FNR == 1 { emit(); f = FILENAME; n = 0; w = 0; delete sp }

/^>/ {
    n++
    id = substr($0, 2)
    if (match(id, /^GC[AF]_[0-9]+\.[0-9]+/)) acc = substr(id, 1, RLENGTH)
    else { acc = id; sub(/\|.*/, "", acc) }
    sp[acc]++
    next
}

{ if (n == 1) w += length($0) }

END { emit() }
