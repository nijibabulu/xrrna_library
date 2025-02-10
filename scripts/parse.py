import itertools
from pathlib import Path
from typing import List
import click
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation


def write_feat(rec: SeqRecord, feat: SeqFeature, out_stem: str, output_dir: Path):
    utr_rec = feat.extract(rec)
    feat_type_slug = feat.type.replace("'", "").replace(" ", "_")
    range_str = f"{feat.location.start}-{feat.location.end}.{feat_type_slug}"
    utr_rec.description += f" {feat.type} {range_str}"
    out_path = output_dir / f"{out_stem}_{range_str}.fa"
    SeqIO.write(utr_rec, out_path, "fasta")


def create_intergenic_feature(cds1: SeqFeature, cds2: SeqFeature):
    start = cds1.location.end + 1
    end = cds2.location.start - 1
    strand = cds1.location.strand  # should always be +
    return SeqFeature(
        location=SimpleLocation(start, end, strand=strand),
        type="intergenic",
        qualifiers={"gene": ["intergenic"]},
    )


def extract_features(feats: List[SeqFeature], min_intergenic=50):
    cdss = [feat for feat in feats if feat.type == "CDS"]
    igs = [
        create_intergenic_feature(cds1, cds2)
        for cds1, cds2 in itertools.pairwise(cdss)
        if cds2.location.start - cds1.location.end >= min_intergenic
    ]
    utrs = [feat for feat in feats if feat.type == "3'UTR"]
    return igs + utrs


@click.command()
@click.option("--min_intergenic", type=int, default=50)
@click.argument("input_dir", type=click.Path(exists=True, file_okay=False, path_type=Path))
@click.argument("output_dir", type=click.Path(exists=True, file_okay=False, path_type=Path))
def main(input_dir: Path, output_dir: Path, min_intergenic: int):
    _ = [
        write_feat(rec, feat, file.stem, output_dir)
        for file in input_dir.glob("*.gbk")
        for rec in SeqIO.parse(file, "genbank")
        for feat in extract_features(rec.features, min_intergenic)
    ]


if __name__ == "__main__":
    main()
