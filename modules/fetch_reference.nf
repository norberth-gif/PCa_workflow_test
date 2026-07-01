process FETCH_REFERENCE {
    tag "Ensembl_${params.reference.release}"
    publishDir "${params.outdir}/reference", mode: 'copy'

    output:
    path "${params.reference.fasta_file}", emit: fasta
    path "${params.reference.gtf_file}", emit: gtf

    script:
    """
    wget -O ${params.reference.fasta_file}.gz ${params.reference.fasta_url}
    gunzip -f ${params.reference.fasta_file}.gz

    wget -O ${params.reference.gtf_file}.gz ${params.reference.gtf_url}
    gunzip -f ${params.reference.gtf_file}.gz
    """
}
