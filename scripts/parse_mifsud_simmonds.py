#! /usr/bin/env python3

import sys
import pandas as pd

if __name__ == "__main__":
    dfm = pd.read_excel(sys.argv[1], skiprows=2)
    dfs = pd.read_excel(sys.argv[2])
    # get the accessions without version numbers from the Mifsud to compare with the Simmonds list
    accm = set(
        dfm.nucl_sequence_accession.str.replace("[.].*$", "", regex=True).unique()
    ).difference(["novel", "hou_et_al"])
    accs = set(acc for acc in dfs.Accession_number if not acc.startswith("SRR"))
    print("\n".join(accm.union(accs)))
