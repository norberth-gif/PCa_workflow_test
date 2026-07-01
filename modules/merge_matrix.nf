process MERGE_MATRIX {
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path count_files

    output:
    path "raw_counts_matrix.txt"

    script:
    """
    python3 - <<'PY'
    from pathlib import Path

    count_files = sorted(Path(".").glob("*.counts"))
    if not count_files:
        raise SystemExit("No .counts files found")

    sample_names = [path.name.removesuffix(".counts") for path in count_files]
    genes = []
    counts_by_sample = {}

    for i, path in enumerate(count_files):
        sample_counts = {}
        with path.open() as handle:
            for line_no, line in enumerate(handle, start=1):
                fields = line.rstrip("\\n").split("\\t")
                if len(fields) != 2:
                    raise SystemExit(f"{path}:{line_no} has {len(fields)} columns, expected 2")

                gene_id, count = fields
                if gene_id.startswith("__") or gene_id.startswith("N_"):
                    continue

                sample_counts[gene_id] = count
                if i == 0:
                    genes.append(gene_id)

        counts_by_sample[path.name.removesuffix(".counts")] = sample_counts

    with open("raw_counts_matrix.txt", "w") as out:
        out.write("GeneID\\t" + "\\t".join(sample_names) + "\\n")
        for gene_id in genes:
            row = [gene_id]
            for sample in sample_names:
                row.append(counts_by_sample[sample].get(gene_id, "0"))
            out.write("\\t".join(row) + "\\n")
    PY
    """
}
