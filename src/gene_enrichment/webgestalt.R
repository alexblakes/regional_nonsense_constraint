# Run ORA with WebGestalt.

library(WebGestaltR)
library(dplyr)

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

    # Get nested list of gene IDs
    ids <- c(strsplit(df$userId, ';'))
    names(ids) <- c(df$description)

    # Run weighted set cover algorithm
    weightedSetCover(
        ids, 
        costs=1/(-log(df$pValue)), 
        topN=10
        )
}

main <- function(ref, refName, db, dbName, interestGenes, regionName) {
    # Run as script.

    df <- run_ora(interestGenes, ref, db)

    # Exit the function cleanly if no significant results found.
    if (is.null(df)) {
        return(NULL)
    }

    # Get top sets
    topSets = run_wsc(df)$topSets

    # Tidy the data
    df <- df[df$description %in% topSets, ]
    df <- df[order(df$enrichmentRatio), ]
    df$region <- regionName
    df$ref <- refName
    df$db <- dbName

    return(df)
}

# Constants
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

# Command line arguments
args = commandArgs(trailingOnly = TRUE)

#  Concatenate all tables
df <- bind_rows(mapply(main, referenceFiles, referenceNames, dbs, dbNames, args[1], args[2]))

# Write to output
write.table(
    df, 
    file=paste0("data/statistics/ora_", args[2], ".tsv"), 
    quote=FALSE, 
    sep="\t", 
    row.names=FALSE
)