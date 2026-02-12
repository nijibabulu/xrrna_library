#! /usr/bin/env python3

import itertools
import re
import pandas as pd

file = "/Users/rpz/Documents/Projects/xrRNA/work/library_mifsud_simmonds_EDIT_inserts_only.cdhit.fa.clstr"


content = open(file).read()
lines = content.splitlines()

clusters = [list(g) for k, g in itertools.groupby(lines, lambda x: x.startswith(">")) if not k]
cluster_groups = [c for c in clusters if len(c) > 1]
cluster_group_accs = [
    [re.sub(r"_\d+_\d+.*$", "", s.split(">")[1]) for s in c] for c in cluster_groups
]

pd.DataFrame(cluster_group_accs).drop_duplicates().to_csv(
    "work/mifsud_simmonds_duplicate_seqs.csv", index=False, header=False, na_rep=""
)
