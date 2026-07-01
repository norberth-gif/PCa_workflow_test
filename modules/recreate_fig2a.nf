process RECREATE_FIG2A {
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path counts_matrix
    path stats_file
    path sra_csv

    output:
    path "Figure_2A_Recreated.pdf", emit: plotA
    path "recreated_clusters.csv", emit: clusters

    script:
    def levels = params.experimental_design.condition_levels.collect { "\"${it}\"" }.join(', ')
    """
    ulimit -s unlimited || true

    cat > recreate_fig2a.R <<'RSCRIPT'
    library(DESeq2); library(DEGreport); library(ggplot2)
    options(expressions = 500000)

    # load data, this is needed for DESeq and plotting
    counts <- read.table("${counts_matrix}", header=TRUE, sep="\\t", row.names=1, check.names=FALSE)
    sra    <- read.csv("${sra_csv}", stringsAsFactors=FALSE)
    stats  <- read.csv("${stats_file}", stringsAsFactors=FALSE)

    # make the metadata table from the SRA file
    meta <- data.frame(
        sample = sra\$${params.experimental_design.sample_column},
        condition = sapply(
            strsplit(as.character(sra\$${params.experimental_design.treatment_column}), "${params.experimental_design.treatment_delimiter}"),
            function(x) trimws(tail(x, 1))
        ),
        row.names = sra\$${params.experimental_design.sample_column}
    )

    # keep only samples with enough unique mapping
    keep_ids <- stats\$sample[stats\$percent > ${params.analysis.min_unique_mapping_percent}]
    common <- intersect(intersect(colnames(counts), rownames(meta)), keep_ids)

    meta_sub <- meta[common, , drop=FALSE]
    target_levels <- c(${levels})
    meta_sub\$condition <- factor(meta_sub\$condition, levels = target_levels)

    # DESeq part
    dds <- DESeqDataSetFromMatrix(countData = round(counts[,common]), colData = meta_sub, design = ~ condition)
    dds <- DESeq(dds, test="LRT", reduced=~1)
    res <- results(dds)

    sig_res <- res[which(res\$padj < ${params.analysis.padj_cutoff}), ]
    sig_res_ordered <- sig_res[order(sig_res\$pvalue), ]
    final_genes <- head(rownames(sig_res_ordered), ${params.analysis.top_gene_count})

    vsd_mat   <- assay(vst(dds, blind=FALSE))
    clean_mat <- vsd_mat[final_genes, ]

    message(paste("Clustering with", length(final_genes), "significant genes."))
    pdf(NULL)
    clusters <- degPatterns(clean_mat,
                            metadata = meta_sub,
                            time = "condition",
                            col = NULL,
                            consensusCluster = FALSE,
                            minc = ${params.analysis.min_cluster_size},
                            summarize = "merge",
                            plot = TRUE)
    dev.off()

    write.csv(clusters\$normalized, "recreated_clusters.csv", row.names=FALSE)

    all_clusters <- unique(clusters\$normalized\$cluster)

    pdf("Figure_2A_Recreated.pdf", width=8, height=6)

    for (i in all_clusters) {
        cluster_df <- clusters\$normalized[clusters\$normalized\$cluster == i, ]

        p <- ggplot(cluster_df, aes(x = condition, y = value, color = condition)) +
             geom_boxplot(size = 1, outlier.shape = 16) +
             scale_color_manual(values = c("CTR"="#666666", "PFA"="#008000", "OO"="#FF0000", "PFA_OO"="#FFA500")) +
             theme_bw() +
             labs(y = "Scaled expression", x = "",
                  title = paste("Cluster", i, ":", length(unique(cluster_df\$genes)), "genes")) +
             theme(legend.position = "none")

        print(p)
    }

    dev.off()
    RSCRIPT

    Rscript recreate_fig2a.R
    """
}
