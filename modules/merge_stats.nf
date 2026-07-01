process MERGE_STATS {
    publishDir "${params.outdir}/stats", mode: 'copy'

    input:
    path indiv_stats

    output:
    path "alignment_stats.csv"

    script:
    """
    echo "sample,percent" > alignment_stats.csv
    cat ${indiv_stats} >> alignment_stats.csv
    """
}
