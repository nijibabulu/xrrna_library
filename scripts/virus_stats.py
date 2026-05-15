#! /usr/bin/env python3

from pathlib import Path

import pandas as pd
import click
from Bio import SeqIO


def get_taxon(rec, accession):
    ranks = ["domain", "realm", "kingdom", "phylum", "class", "order", "family", "genus"]
    taxonomy = rec.annotations.get("taxonomy", [])
    return (
        {rank: taxonomy[i] if i < len(taxonomy) else "Unannotated" for i, rank in enumerate(ranks)}
        | {"species": rec.annotations.get("organism", "Unannotated")}
        | {"accession": accession}
    )


@click.command()
@click.argument("input_dir", type=click.Path(exists=True, path_type=Path))
@click.argument("output_dir", type=click.Path(exists=True, path_type=Path))
def main(input_dir, output_dir):
    included_accessions = {
        ".".join("_".join(f.stem.split("_")[:2]).split(".")[:2]) for f in output_dir.glob("*.fa")
    }
    pd.DataFrame(
        [
            get_taxon(rec, file.stem)
            for file in input_dir.glob("*.gbk")
            for rec in SeqIO.parse(file, "genbank")
            if file.stem in included_accessions
        ]
    ).to_csv(output_dir / "virus_stats.csv", index=False)


if __name__ == "__main__":
    main()
