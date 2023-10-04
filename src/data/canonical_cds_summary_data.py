    # Save intermediate outputs
    ## Canonical transcript list
    pd.Series(cds.transcript_id.unique()).to_csv(
        "../outputs/canonical_transcripts.txt", sep="\t", index=False, header=False
    )
    ## ENSG, ENST, and HGNC symbols
    cds[["gene_id", "transcript_id", "gene_name"]].drop_duplicates().to_csv(
        "../outputs/gene_ids.tsv", sep="\t", index=False
    )