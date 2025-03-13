from pathlib import Path
from typing import List
import click
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation


def write_feat(rec: SeqRecord, feat: SeqFeature, out_stem: str, output_dir: Path):
    feat_rec = feat.extract(rec)
    feat_type_slug = feat.type.replace("'", "").replace(" ", "_")
    range_str = f"{feat.location.start}-{feat.location.end}.{feat_type_slug}"
    feat_rec.description += f" {feat.type} {range_str}"
    out_path = output_dir / f"{out_stem}_{range_str}.fa"
    SeqIO.write(feat_rec, out_path, "fasta")


def create_intergenic_feature(start: SeqFeature, end: SeqFeature, strand):
    return SeqFeature(
        location=SimpleLocation(start, end, strand=strand),
        type="intergenic",
        qualifiers={"gene": ["intergenic"]},
    )


def find_intergenic_features(cdss: list[SeqFeature], min_intergenic=50):
    if not cdss:
        return []

    cdss.sort(key=lambda x: x.location.start)
    igs = []

    cur_end = cdss[0].location.end
    for cds in cdss[1:]:
        if cds.location.start - 1 > cur_end + min_intergenic + 1:
            igs.append(create_intergenic_feature(cur_end + 1, cds.location.start - 1, strand=cds.location.strand))
        cur_end = max(cur_end, cds.location.end)

    return igs


def extract_features(feats: List[SeqFeature], min_intergenic=50):
    cdss = [feat for feat in feats if feat.type == "CDS"]
    igs = find_intergenic_features(cdss, min_intergenic)
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
