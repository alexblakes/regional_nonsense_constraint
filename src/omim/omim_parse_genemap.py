"""
Parse the genemap2.txt file and restructure it into tidy format.

In particular, the "phenotype" annotations within genemap2 are nested, and require 
manipulation to extract phenotype and inheritance annotations for each entry.

Notes:
- Entries without a phenotype annotation are dropped. These are unlikely to be 
  interesting to us.
- Where an entry / gene is associated with more than one phenotype, the data has been 
  reformated to long (tidy) format. There is one phenotype in each row. Entries / genes 
  may be duplicated.
- Only the phenotype and inheritance data has been extracted from the phenotype string. 
  Other data within the string is lost.
- Many phenotypes do not have an inheritance annotation. This is even true for some 
  phenotypes where the inheritance is explicitly given in the phenotype name (see below 
  for example). This is a limitation of the OMIM data which this script does not 
  address.
- Non-disease, susceptibility, and provisional phenotypes are still included in the end 
  output. 
"""

import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

import src

_LOGFILE = f"data/logs/{Path(__file__).stem}.log"
_FILE_IN = "data/raw/genemap2.txt"
_FILE_OUT = "data/interim/genemap2_parsed.tsv"

logger = logging.getLogger(__name__)


def read_gm(path):
    """Read genemap2.txt file into memory"""

    logger.info("Reading genemap2.txt.")

    gm = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "chr",
            "start",
            "end",
            "cyto_loc",
            "calc_cyto_loc",
            "mim",
            "symbol",
            "gene_name",
            "approved_symbol",
            "entrez",
            "ensg",
            "comment",
            "phenotype_string",
            "mouse",
        ],
    )

    logger.info(f"Entries in genemap2.txt: {len(gm)}")

    return gm


def split_phenotypes(gm):
    """Split phenotype data within genemap2"""

    logger.info("Splitting phenotype data.")

    # Phenotype data is nested, separated by ";"
    gm["phenotype_string"] = gm.phenotype_string.str.split(";")
    gm = gm.explode("phenotype_string", ignore_index=True)

    logger.info(f"Entries after exploding nested phenotypes: {len(gm)}")

    # Drop entries with no associated phenotype
    # These will rarely be interesting
    gm = gm.dropna(subset="phenotype_string")

    logger.info(f"Entries after dropping NaNs in phenotype column: {len(gm)}")

    return gm


def match_re(string, _re):
    """Find phenotypes strings matching a regular expression"""
    string = string.strip()
    return re.match(_re, string)


def parse_long_phenotypes(gm):
    """Parse the phenotype and inheritance data within the phenotype_string column.

    Phenotypes have either long or short string entries. These are treated differently.

    The regular expressions to parse long and short strings are taken from OMIM's
    GeneMap2.txt Parser, available at
    https://github.com/OMIM-org/genemap2-parser/tree/master
    """

    logger.info("Parsing long phenotypes.")

    re_long = r"^(.*),\s(\d{6})\s\((\d)\)(|, (.*))$"

    # Find long phenotype strings matching the regular expression
    pheno_long = gm.phenotype_string.apply(match_re, _re=re_long).dropna().to_frame()

    logger.info(f"Number of long phenotype entries: {len(pheno_long)}")

    # Parse the phenotype strings to extract the phenotype and inheritance annotation
    pheno_long["phenotype"] = pheno_long.phenotype_string.apply(lambda x: x.group(1))
    pheno_long["inheritance"] = pheno_long.phenotype_string.apply(lambda x: x.group(5))

    # Some entries have multiple inheritance annotations.
    # "Explode" on inheritance, so the data is in tidy format.
    pheno_long["inheritance"] = pheno_long.inheritance.str.split(", ")
    pheno_long = pheno_long.explode("inheritance")

    logger.info(
        f"Number of long phenotypes after exploding on inheritance: {len(pheno_long)}"
    )

    # Replace None with np.nan
    pheno_long = pheno_long.fillna(np.nan)

    # Drop the now unnecessary "phenotype_string" column
    pheno_long = pheno_long.drop("phenotype_string", axis=1)

    # Print summary statistics
    logger.info(
        f"Number of phenotypes lacking an inheritance annotation: {pheno_long.inheritance.isna().sum()}"
    )

    return pheno_long


def parse_short_phenotypes(gm, pheno_long):
    """Parse the phenotype and inheritance data within the phenotype_string column.
    Phenotypes have either long or short string entries.
    These are treated differently.
    The regular expressions to parse long and short strings are taken from
    OMIM's GeneMap2.txt Parser, available at
    https://github.com/OMIM-org/genemap2-parser/tree/master
    """

    logger.info("Parsing short phenotypes.")

    re_short = r"^(.*)\((\d)\)(|, (.*))$"

    # Find short phenotype strings matching the regular expression.
    # Entries matching the long expression need to be excluded.
    pheno_short = (
        gm.loc[~gm.index.isin(pheno_long.index), "phenotype_string"]
        .apply(match_re, _re=re_short)
        .to_frame()
    )

    logger.info(f"Number of short phenotype entries: {len(pheno_short)}")

    # Parse the phenotype strings to extract the phenotype and inheritance annotation
    pheno_short["phenotype"] = pheno_short.phenotype_string.apply(lambda x: x.group(1))
    pheno_short["inheritance"] = pheno_short.phenotype_string.apply(
        lambda x: x.group(3)
    )

    # Drop the unnecessary "phenotype_string" column
    pheno_short = pheno_short.drop("phenotype_string", axis=1)

    # Replace empty inheritance annotations with NA
    pheno_short = pheno_short.replace("", np.nan)

    logger.info(
        f"Number of short entries lacking an inheritance annotation {pheno_short.inheritance.isna().sum()}"
    )

    #! Entries with short phenotype strings (in pheno_short, above) are always missing
    #! an inheritance annotation. This is even true where the inheritance is explicitly
    #! given in the phenotype (e.g. "Deafness, autosomal recessive" is a phenotype, but
    #! has no associated inheritance annotation. Search MIM 603324 online to see this
    #! example.) This is a striking limitation of the OMIM data.

    return pheno_short


def concat_phenotypes(pheno_long, pheno_short):
    """Concatenate the long and short phenotype data"""

    logger.info("Combining long and short phenotype data.")

    pheno = pd.concat([pheno_long, pheno_short])

    # Print summary statistics
    logger.info(
        f"Number of non-disease phenotypes {pheno.phenotype.str.startswith('[').sum()}"
    )
    logger.info(
        f"Number of susceptibility phenotypes {pheno.phenotype.str.startswith('{').sum()}"
    )
    logger.info(
        f"Number of provisional phenotypes {pheno.phenotype.str.startswith('?').sum()}"
    )
    logger.info(
        f"Value counts of inheritance modes:\n{pheno.inheritance.value_counts()}"
    )

    return pheno


def merge_and_tidy_genemap_phenotypes(gm, pheno):
    """Merge the phenotype annotations with the genemap data.
    Save to output.
    """

    gm = gm.drop("phenotype_string", axis=1)
    gm = gm.merge(pheno, left_index=True, right_index=True)

    logger.info(f"Entries after merging genemap and phenotype data: {len(gm)}")

    return gm


def main():
    """Run as script."""

    # Parse genemap2.txt file
    gm = read_gm(_FILE_IN).pipe(split_phenotypes)
    pheno_long = parse_long_phenotypes(gm)
    pheno_short = parse_short_phenotypes(gm, pheno_long)
    pheno = concat_phenotypes(pheno_long, pheno_short)
    gm = merge_and_tidy_genemap_phenotypes(gm, pheno)

    # Write to output
    gm.to_csv(_FILE_OUT, sep="\t", index=False)


if __name__ == "__main__":
    logger = src.setup_logger(_LOGFILE)
    main()