process VOLCANO_PLOT {

    tag "Generating volcano plots"

    publishDir "${params.outdir}/volcano", mode: 'copy'

    input:
    path counts_matrix
    path sra_csv

    output:
    path "*.pdf", emit: volcano_plots
    path "*.tsv", emit: volcano_results

    script:

    def conditions = params.experimental_design.condition_levels
    def contrast_list = []

    for (i in 0..<conditions.size()) {
        for (j in (i+1)..<conditions.size()) {
            contrast_list << "\"${conditions[i]}::${conditions[j]}\""
        }
    }

    def contrasts = contrast_list.join(", ")

    """
    cat > volcano_plot.R <<'RSCRIPT'

    library(DESeq2)
    library(ggplot2)


    counts <- read.table(
        "${counts_matrix}",
        header=TRUE,
        sep="\\t",
        row.names=1,
        check.names=FALSE
    )


    sra <- read.csv(
        "${sra_csv}",
        stringsAsFactors=FALSE
    )


    meta <- data.frame(
        sample = sra\$${params.experimental_design.sample_column},
        condition = sra\$${params.experimental_design.treatment_column},
        row.names = sra\$${params.experimental_design.sample_column}
    )


    common <- intersect(
        colnames(counts),
        rownames(meta)
    )


    counts <- counts[,common]
    meta <- meta[common,,drop=FALSE]


    meta\$condition <- factor(
        meta\$condition,
        levels=c(
            ${params.experimental_design.condition_levels.collect { "\"${it}\"" }.join(",")}
        )
    )


    dds <- DESeqDataSetFromMatrix(
        countData=round(counts),
        colData=meta,
        design=~condition
    )


    dds <- DESeq(dds)


    contrasts <- c(
        ${contrasts}
    )


    for (contrast_pair in contrasts) {


        split_pair <- strsplit(
            contrast_pair,
            "::"
        )[[1]]


        condition_A <- split_pair[1]
        condition_B <- split_pair[2]


        message(
            "Running comparison: ",
            condition_A,
            " vs ",
            condition_B
        )


        res <- results(
            dds,
            contrast=c(
                "condition",
                condition_A,
                condition_B
            )
        )


        res_df <- as.data.frame(res)

        res_df\$GeneID <- rownames(res_df)


        safe_name <- paste0(
            gsub("[^A-Za-z0-9]+","_",condition_A),
            "_vs_",
            gsub("[^A-Za-z0-9]+","_",condition_B)
        )


        write.table(
            res_df,
            file=paste0(
                safe_name,
                "_DESeq2_results.tsv"
            ),
            sep="\\t",
            quote=FALSE,
            row.names=FALSE
        )


        plot_df <- res_df[
            !is.na(res_df\$padj),
            ]
        

        plot_df\$Significance <- "Not significant"


        plot_df\$Significance[
            plot_df\$padj < ${params.analysis.volcano_padj} &
            plot_df\$log2FoldChange >= ${params.analysis.volcano_log2fc}
        ] <- "Upregulated"


        plot_df\$Significance[
            plot_df\$padj < ${params.analysis.volcano_padj} &
            plot_df\$log2FoldChange <= -${params.analysis.volcano_log2fc}
        ] <- "Downregulated"



        p <- ggplot(
            plot_df,
            aes(
                x=log2FoldChange,
                y=-log10(padj),
                colour=Significance
            )
        ) +

        geom_point(
            alpha=0.6,
            size=1.8
        ) +

        scale_colour_manual(
            values=c(
                "Upregulated"="red",
                "Downregulated"="blue",
                "Not significant"="grey70"
            )
        ) +

        geom_vline(
            xintercept=c(
                -${params.analysis.volcano_log2fc},
                 ${params.analysis.volcano_log2fc}
            ),
            linetype="dashed"
        ) +

        geom_hline(
            yintercept=-log10(${params.analysis.volcano_padj}),
            linetype="dashed"
        ) +

        theme_bw() +

        labs(
            title=paste(
                condition_A,
                "vs",
                condition_B
            ),
            x="log2 Fold Change",
            y="-log10 adjusted p-value"
        )


        pdf(
            paste0(
                safe_name,
                "_Volcano_plot.pdf"
            ),
            width=8,
            height=6
        )

        print(p)

        dev.off()

    }

    RSCRIPT


    Rscript volcano_plot.R
    """
}
