process EXTRACT_STAT {
    tag "$star_log"

    input:
    path star_log

    output:
    path "*_stat.txt"

    script:
    """
    sample_id=\$(basename ${star_log} .Log.final.out)
    pct=\$(grep "Uniquely mapped reads %" ${star_log} | awk -F'|' '{print \$2}' | sed 's/%//g' | xargs)
    echo "\$sample_id,\$pct" > \${sample_id}_stat.txt
    """
}

