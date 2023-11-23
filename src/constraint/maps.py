"""Get the proportion of singletons for synonymous variant contexts."""

# Imports
from pathlib import Path

import pandas as pd

from src import setup_logger
from src import constants as C

# Module constants


# Logging
logger = setup_logger(Path(__file__).stem)

# Functions


def main():
    """Run as script."""

    return df  #! Testing


if __name__ == "__main__":
    main()


# ps.to_csv("../outputs/proportion_singletons_synonymous_by_context.tsv", sep="\t", index=False)


# # ### Proportion of singletons by csq


# m1 = df["ac"] != 0
# m2 = df["variant_type"] != "CpG"
# m3 = df["csq"] == "nonsense"


# # #### Synonymous, missense, and nonsense


# # CpG only
# ps_csq_cpg = get_ps(df[m1 & ~m2].groupby("csq"))

# # Exclude CpG
# ps_csq_no_cpg = get_ps(df[m1 & m2].groupby("csq"))


# # #### Nonsense, by NMD region


# # CpG only
# ps_region_cpg = get_ps(df[m1 & m3 & ~m2].groupby("nmd"))

# # Exclude CpG
# ps_region_no_cpg = get_ps(df[m1 & m2 & m3].groupby("nmd"))


# # ### Combine CSQ and Region results


# # Reformat column names
# ps_region_cpg = ps_region_cpg.rename(columns={"nmd":"csq"})
# ps_region_no_cpg = ps_region_no_cpg.rename(columns={"nmd":"csq"})

# for a, b, c in zip([ps_csq_cpg, ps_csq_no_cpg], [ps_region_cpg, ps_region_no_cpg], ["cpg","no_cpg"]):
#     ps = pd.concat([a,b])
#     ps.to_csv(f"../outputs/proportion_singletons_by_csq_and_region_{c}.tsv", sep="\t", index=False)
#     print(f"{c}\n{ps}\n\n")
