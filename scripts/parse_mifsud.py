#! /usr/bin/env python3

import sys
import pandas as pd

if __name__ == "__main__":
    df = pd.read_excel(sys.argv[1], skiprows=2)
    print("\n".join(df.nucl_sequence_accession.unique()))
