from pathlib import Path
import click
from Bio import SeqIO, SeqRecord, SeqFeature


def write_utr3(rec: SeqRecord, feat: SeqFeature, out_stem: str, output_dir: Path):
    utr_rec = feat.extract(rec)
    range_str = f"{feat.location.start}-{feat.location.end}"
    utr_rec.description += f" 3'UTR {range_str}"
    out_path = output_dir / f"{out_stem}_{range_str}.fa"
    SeqIO.write(utr_rec, out_path, "fasta")


@click.command()
@click.argument("input_dir", type=click.Path(exists=True, file_okay=False, path_type=Path))
@click.argument("output_dir", type=click.Path(exists=True, file_okay=False, path_type=Path))
def main(input_dir: Path, output_dir: Path):
    _ = [
        write_utr3(rec, feat, file.stem, output_dir)
        for file in input_dir.glob("*.gbk")
        for rec in SeqIO.parse(file, "genbank")
        for feat in rec.features
        if feat.type == "3'UTR"
    ]


if __name__ == "__main__":
    main()
