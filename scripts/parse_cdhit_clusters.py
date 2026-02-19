#! /usr/bin/env python3

import itertools
import re
import sys
from typing import Annotated
from pathlib import Path
import pandas as pd
import typer

# file = "/Users/rpz/Documents/Projects/xrRNA/work/library_mifsud_simmonds_EDIT_inserts_only.cdhit.fa.clstr"


def main(
    file: Annotated[Path, typer.Argument(exists=True, dir_okay=False)],
    output_csv: Annotated[Path, typer.Argument(exists=False, dir_okay=False)],
    output_redundant_accs: Annotated[Path, typer.Argument(exists=False, dir_okay=False)],
) -> None:
    content = open(file).read()
    lines = content.splitlines()

    clusters = [list(g) for k, g in itertools.groupby(lines, lambda x: x.startswith(">")) if not k]
    cluster_groups = [c for c in clusters if len(c) > 1]
    cluster_group_accs = [
        [re.sub(r"\.\.\..*$", "", s.split(">")[1]) for s in c] for c in cluster_groups
    ]
    cluster_group_lengths = [
        [s.split("nt,")[0].split()[1].strip() for s in c] for c in cluster_groups
    ]

    # construct a list of lists of alternating accessions and lengths, sorted by length in descending order within each cluster group
    sorted_seqs_and_accs = [
        list(
            itertools.chain.from_iterable(sorted(zip(accs, lens), key=lambda p: p[1], reverse=True))
        )
        for accs, lens in zip(cluster_group_accs, cluster_group_lengths)
    ]
    pd.DataFrame(sorted_seqs_and_accs).drop_duplicates().to_csv(
        output_csv, index=False, header=False, na_rep=""
    )

    with open(output_redundant_accs, "w") as f:
        for cluster in sorted_seqs_and_accs:
            # output all but the first accession
            accs = [acc for acc in cluster[2::2]]
            f.write("\n".join(accs) + "\n")


if __name__ == "__main__":
    typer.run(main)
