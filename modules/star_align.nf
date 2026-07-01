process STAR_ALIGN {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)
    path star_index

    output:
    path "${sample_id}.counts", emit: counts
    path "${sample_id}.Log.final.out", emit: log

    script:
    """
    case "${params.counts.strandedness}" in
        no) count_col=2 ;;
        yes) count_col=3 ;;
        reverse) count_col=4 ;;
        *)
            echo "Unsupported strandedness: ${params.counts.strandedness}" >&2
            exit 1
            ;;
    esac

    STAR --genomeDir ${star_index} \\
         --readFilesIn ${reads[0]} ${reads[1]} \\
         --readFilesCommand zcat \\
         --outSAMtype None \\
         --quantMode GeneCounts \\
         --outFileNamePrefix ${sample_id}. \\
         --runThreadN ${task.cpus}

    awk -v col="\$count_col" 'BEGIN { OFS="\\t" } !/^N_/ { print \$1, \$col }' \\
        ${sample_id}.ReadsPerGene.out.tab > ${sample_id}.counts
    """
}
