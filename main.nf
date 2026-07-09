nextflow.enable.dsl=2

include { FETCH_REFERENCE }         from './modules/fetch_reference'
include { FASTQC }                  from './modules/fastqc'
include { MERGE_MATRIX }            from './modules/merge_matrix'
include { PREPARE_REF_GENOME }      from './modules/prepare_ref_genome'
include { RECREATE_FIG2A }          from './modules/recreate_fig2a' 
include { STAR_ALIGN }              from './modules/star_align'
include { EXTRACT_STAT }            from './modules/extract_stat'
include { MERGE_STATS }             from './modules/merge_stats'
include { FUNCTIONAL_ENRICHMENT }   from './modules/functional_enrichment'
include { MULTIQC }                 from './modules/multiqc' 
include { VOLCANO_PLOT }            from './modules/volcano_plot'

workflow {
    // first get the pig reference files
    ref_fetch_ch = FETCH_REFERENCE()

    ch_fasta = ref_fetch_ch.fasta.first()
    ch_gtf   = ref_fetch_ch.gtf.first()

    // make STAR index from the fasta and gtf
    ch_star_idx = PREPARE_REF_GENOME(ch_fasta, ch_gtf).index.first()

    // all paired fastq files, names should match params.reads
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

    // quick read quality check
    fastqc_ch = FASTQC(reads_ch)

    // align reads and then count them
    align_ch = STAR_ALIGN(reads_ch, ch_star_idx)

    // get the uniquely mapped % from the STAR log files
    single_stats_ch = EXTRACT_STAT(align_ch.log)
    stats_ch = MERGE_STATS(single_stats_ch.collect())

    counts_ch = align_ch.counts
    matrix_ch = MERGE_MATRIX(counts_ch.collect())

    // this only gets the fastqc files at the moment
    MULTIQC(fastqc_ch.collect())

    // recreate figure 2 parts
    viz_ch = RECREATE_FIG2A(matrix_ch, stats_ch, file(params.sra_csv))
    FUNCTIONAL_ENRICHMENT(viz_ch.clusters)
    
    // volcano plots for all DESeq2 contrasts
    VOLCANO_PLOT(
        matrix_ch,
        file(params.sra_csv)
    )
}


