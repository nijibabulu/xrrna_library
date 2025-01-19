from pathlib import Path
import click
from Bio import SeqIO, SeqFeature


@click.command()
@click.argument(
    "input_dir", type=click.Path(exists=True, file_okay=False, path_type=Path)
)
@click.argument(
    "output_dir", type=click.Path(exists=True, file_okay=False, path_type=Path)
)
def main(input_dir: Path, output_dir: Path):
    for file in input_dir.iterdir():
        if file.suffix == ".gbk":
            rec = SeqIO.read(file, "genbank")
            utr3_feats = [feat for feat in rec.features if feat.type == "3'UTR"]
            for feat in utr3_feats:
                utr_rec = feat.extract(rec)
                range_str = f"{feat.location.start}-{feat.location.end}"
                utr_rec.description += f" 3'UTR {range_str}"
                SeqIO.write(
                    utr_rec,
                    output_dir / f"{file.stem}_{range_str}.fa",
                    "fasta",
                )


if __name__ == "__main__":
    main()
