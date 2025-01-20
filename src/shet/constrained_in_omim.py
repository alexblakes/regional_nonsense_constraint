"""Find OMIM morbid genes with regional nonsense constraint."""



import pandas as pd

import src

FILE_OMIM = "data/interim/genemap2_simple.tsv"
FILE_CONSTRAINT = "data/interim/shet_gnomad_regional_constraint.tsv"
FILE_OUT = "data/interim/shet_gnomad_regional_omim_ad.tsv"

logger = src.logger


def read_constraint(path=FILE_CONSTRAINT):
    df = pd.read_csv(path, sep="\t")
    logger.info(f"Constraint entries: {len(df)}")
    logger.info(f"Constraint unique ENSGs: {df.ensg.nunique()}")

    return df

def read_omim(path=FILE_OMIM):
    df = pd.read_csv(path, sep="\t")

    logger.info(f"OMIM data shape: {df.shape}")
    logger.info(f"Unique ENSGs in OMIM: {df.ensg.nunique()}")
    logger.info(
        f"Unique by ENSG and mode of inheritance: {len(df.drop_duplicates(['ensg','inheritance']))}"
    )

    return df


def find_unique_ad_genes(df):
    ad_genes = df.loc[df.inheritance == "Autosomal dominant"].ensg.drop_duplicates()

    logger.info(f"Unique AD genes: {len(ad_genes)}")

    return ad_genes

def write_out(df, path=FILE_OUT):
    df.to_csv(path, sep="\t", index=False)
    return df

def main():
    """Run as script."""
    constraint = read_constraint(FILE_CONSTRAINT)
    morbid_ad_genes = read_omim(FILE_OMIM).pipe(find_unique_ad_genes)

    df = constraint.merge(morbid_ad_genes, how="inner", validate="one_to_one")
    logger.info(f"ENSGs after merge: {len(df)}")

    write_out(df)

    return df

if __name__ == "__main__":
    logger = src.setup_logger(src.log_file(__file__))
    main()
