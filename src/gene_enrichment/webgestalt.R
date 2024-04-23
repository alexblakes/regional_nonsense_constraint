# Run ORA with WebGestalt.

library(WebGestaltR)

run_ora <- function(interestGenes, referenceGenes, database) {
    df <- WebGestaltR(
        enrichMethod="ORA",
        organism="hsapiens",
        enrichDatabase=database,
        interestGeneFile=interestGenes,
        interestGeneType="ensembl_gene_id",
        referenceGeneFile=referenceGenes,
        referenceGeneType="ensembl_gene_id",
        minNum=10,
        maxNum=500,
        sigMethod="fdr",
        fdrThr=0.05,
        topThr=10,
        isOutput=FALSE,
        setCoverNum=10,
    )
}

run_wsc <- function(df) {

    ids <- c(strsplit(df$userId, ';'))
    names(ids) <- c(df$description)

    weightedSetCover(
        ids, 
        costs=1/(-log(df$pValue)), 
        topN=10
        )

}

main <- function(ref, refName, db, dbName, interestGenes, regionName) {

    df <- run_ora(interestGenes, ref, db)

    if (is.null(df)) {
        return(NA)
    }

    topSets = run_wsc(df)$topSets

    df <- df[df$description %in% topSets, ]
    df <- df[order(df$enrichmentRatio), ]

    outFile = paste0("data/statistics/ora_", regionName, "_vs_", refName, "_", dbName, ".tsv")

    write.table(
        df, 
        file=outFile, 
        quote=FALSE, 
        sep="\t", 
        row.names=FALSE,
    )
}

referenceFiles = c(
    rep(c("data/final/gene_list_all.txt"), 3),
    rep(c("data/final/gene_list_gnomad_constrained.txt"), 3)
)
referenceNames = c(
    rep(c("all"), 3), 
    rep(c("gnomad"), 3)
)
dbs = rep(
    c(
        "geneontology_Biological_Process_noRedundant",
        "geneontology_Molecular_Function_noRedundant",
        "phenotype_Human_Phenotype_Ontology"
    ),
    2
)
dbNames = rep(c("bp","mf","hpo"), 2)

args = commandArgs(trailingOnly = TRUE)

mapply(main, referenceFiles, referenceNames, dbs, dbNames, args[1], args[2])
