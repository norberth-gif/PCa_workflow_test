process FUNCTIONAL_ENRICHMENT {
    tag "Barplot q < ${params.enrichment.qvalue_cutoff}"
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path cluster_csv

    output:
    path "Enrichment_Barplot.pdf", emit: plotB

    script:
    """
    #!/usr/bin/env Rscript
    library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2)

    write_message_pdf <- function(message) {
        pdf("Enrichment_Barplot.pdf", width=11, height=8)
        plot.new()
        text(0.5, 0.5, message)
        dev.off()
    }

    dat <- read.csv("${cluster_csv}")
    ensembl_genes <- unique(dat[dat\$cluster == ${params.enrichment.cluster_id}, "genes"])
    ensembl_genes <- ensembl_genes[!is.na(ensembl_genes)]
    ensembl_genes <- sub("\\\\..*", "", ensembl_genes)

    if(length(ensembl_genes) == 0) {
        write_message_pdf("No genes found for cluster ${params.enrichment.cluster_id}.")
        quit(save = "no", status = 0)
    }

    gene2ensembl_url <- "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz"
    gene2ensembl_file <- "gene2ensembl.gz"

    download_ok <- FALSE
    for(attempt in seq_len(3)) {
        download_ok <- tryCatch({
            download.file(gene2ensembl_url, gene2ensembl_file, mode = "wb", quiet = TRUE)
            TRUE
        }, error = function(e) {
            message("NCBI gene2ensembl download attempt ", attempt, " failed: ", conditionMessage(e))
            FALSE
        })
        if(download_ok) break
        Sys.sleep(10 * attempt)
    }

    if(!download_ok) {
        stop("Could not download NCBI gene2ensembl mapping after 3 attempts.")
    }

    gene2ensembl <- read.delim(gzfile(gene2ensembl_file), stringsAsFactors = FALSE, check.names = FALSE)
    tax_col <- names(gene2ensembl)[1]
    pig_map <- gene2ensembl[
        gene2ensembl[[tax_col]] == 9606 &
        gene2ensembl[["Ensembl_gene_identifier"]] %in% ensembl_genes,
        c("GeneID", "Ensembl_gene_identifier")
    ]
    pig_map <- unique(pig_map[!is.na(pig_map\$GeneID), ])

    write.table(
        pig_map,
        file = "cluster_${params.enrichment.cluster_id}_ensembl_to_entrez.tsv",
        sep = "\\t",
        quote = FALSE,
        row.names = FALSE
    )

    if(nrow(pig_map) == 0) {
        write_message_pdf(
            paste(
                "No Entrez ID mapping found for cluster ${params.enrichment.cluster_id}.",
                "Mapping source: NCBI gene2ensembl, taxid 9606.",
                sep = "\\n"
            )
        )
        quit(save = "no", status = 0)
    }

    message("Mapped ", length(unique(pig_map\$GeneID)), " of ", length(ensembl_genes), " cluster genes to Entrez IDs.")

    ego <- enrichGO(
        gene          = as.character(unique(pig_map\$GeneID)),
        OrgDb         = org.Hs.eg.db,
        keyType       = "ENTREZID",
        ont           = "${params.enrichment.ontology}",
        pAdjustMethod = "${params.enrichment.p_adjust_method}",
        pvalueCutoff  = ${params.enrichment.pvalue_cutoff},
        qvalueCutoff  = ${params.enrichment.qvalue_cutoff}
    )

    pdf("Enrichment_Barplot.pdf", width=11, height=8)

    res_df <- as.data.frame(ego)
    write.table(
        res_df,
        file = "Enrichment_plot.tsv",
        sep = "\\t",
        quote = FALSE,
        row.names = FALSE
    )

    if(nrow(res_df) == 0) {
        plot.new()
        text(0.5, 0.5, "No enriched GO terms found for cluster ${params.enrichment.cluster_id}.")
        dev.off()
        quit(save = "no", status = 0)
    }

    res_df\$logQ <- -log10(res_df\$qvalue)

    res_df <- res_df[order(res_df\$logQ, decreasing = TRUE), ]
    top_df <- head(res_df, ${params.enrichment.top_terms})

    p_bar <- ggplot(top_df,
        aes(x = reorder(Description, logQ),
            y = logQ,
            fill = logQ)) +
        geom_bar(stat="identity") +
        coord_flip() +
        scale_fill_gradient(low="lightgray", high="green",
                            name="-log10(q-value)") +
        theme_minimal() +
        labs(x="GO term",
             y="-log10(q-value)",
             title="Top Enriched Functions")

    print(p_bar)

    dev.off()
    """
}
