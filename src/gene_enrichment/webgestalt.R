# Run ORA with WebGestalt.

library(WebGestaltR)

# enrichDatabase options: 
#   - phenotype_Human_Phenotype_Ontology
#   - geneontology_Biological_Process_noRedundant
#   - geneontology_Molecular_Function_noRedundant

run_ora <- function(interestGenes, referenceGenes) {

    WebGestaltR(
        enrichMethod="ORA",
        organism="hsapiens",
        enrichDatabase="phenotype_Human_Phenotype_Ontology",
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

wsc <- function(df) {

    ids <- c(strsplit(df$userId, ';'))
    names(ids) <- c(df$description)

    weightedSetCover(
        ids, 
        costs=1/(-log(df$pValue)), 
        topN=10
        )

}

interestGenes = "data/final/gene_list_distal_constrained.txt"
referenceGenes = "data/final/gene_list_all.txt"

df <- run_ora(interestGenes, referenceGenes)

topSets = wsc(df)$topSets
df <- df[df$description %in% topSets, ]

df <- df[order(enrichmentRatio), ]

write.table(
    df, 
    file="data/statistics/ora_test.tsv", 
    quote=FALSE, 
    sep="\t", 
    row.names=FALSE
)