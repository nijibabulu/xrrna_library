#! /usr/bin/env python3

from itertools import zip_longest
from pathlib import Path

import pandas as pd
import click
from Bio import SeqIO


def get_taxon(rec):
    return dict(
        zip_longest(
            ["virus", "genus", "realm", "phylum", "class", "order", "family"],
            rec.annotations["taxonomy"],
            fillvalue="Unannotated",
        )
    )


@click.command()
@click.argument("input_dir", type=click.Path(exists=True, path_type=Path))
@click.argument("output_dir", type=click.Path(exists=True, path_type=Path))
def main(input_dir, output_dir):
    included_accessions = {"_".join(f.stem.split("_")[:2]) for f in output_dir.glob("*.fa")}
    pd.DataFrame(
        [
            get_taxon(rec)
            for file in input_dir.glob("*.gbk")
            for rec in SeqIO.parse(file, "genbank")
            if file.stem in included_accessions
        ]
    ).to_csv(output_dir / "virus_stats.csv", index=False)


if __name__ == "__main__":
    main()
