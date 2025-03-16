#! /usr/bin/env python3

from pathlib import Path
import click
from Bio import SeqIO
import pandas as pd


@click.command()
@click.argument("library_fa", type=click.Path(exists=True, path_type=Path))
@click.argument("sequence_dir", type=click.Path(exists=True, file_okay=False, path_type=Path))
@click.argument("output_tbl", type=click.Path(dir_okay=False))
def main(library_fa, sequence_dir, output_tbl):
    recs = list(SeqIO.parse(library_fa, "fasta"))
    tbl = (
        pd.DataFrame(
            [
                {
                    "id": rec.id,
                    "description": rec.description.split()[1],
                    "length": len(rec.seq),
                }
                for rec in recs
            ]
        )
        .groupby("description")
        .agg(length=("length", "sum"), n=("id", "size"))
        .reset_index()
    )
    tbl[["genbank_id", "version", "genome_range", "feature_type"]] = tbl.description.str.split(
        ".", expand=True
    )
    tbl[["start", "end"]] = tbl.genome_range.str.split("-", expand=True)
    tbl["start"] = tbl.start.str.replace("[^0-9]", "", regex=True).astype(int)
    tbl["end"] = tbl.end.str.replace("[^0-9]", "", regex=True).astype(int)
    tbl["genome_length"] = tbl.end - tbl.start
    tbl["genbank_id"] = tbl.genbank_id + "." + tbl.version
    tbl.drop(columns=["description", "genome_range", "version"], inplace=True)
    gbks = {f.stem: SeqIO.read(f, "genbank") for f in sequence_dir.glob("*.gbk")}
    tbl["organism"] = tbl.genbank_id.map(lambda x: gbks[x].annotations["organism"])
    tbl.groupby(["genbank_id", "organism", "feature_type"]).agg(
        length=("length", "sum"),
        genome_length=("length", "sum"),
        n=("n", "sum"),  
    ).reset_index().to_csv(output_tbl, index=False)


if __name__ == "__main__":
    main()
