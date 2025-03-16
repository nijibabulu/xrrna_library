from pathlib import Path
import click
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation


def write_feat(rec: SeqRecord, feat: SeqFeature, out_stem: str, output_dir: Path):
    feat_rec = feat.extract(rec)
    feat_type_slug = feat.type.replace("'", "").replace(" ", "_")
    range_str = f"{feat.location.start}-{feat.location.end}.{feat_type_slug}"
    feat_rec.description += f" {feat.type} {range_str}"
    out_path = output_dir / f"{out_stem}.{range_str}.fa"
    SeqIO.write(feat_rec, out_path, "fasta")


def create_intergenic_feature(
    start: SeqFeature, end: SeqFeature, strand, ftype: str = "intergenic"
):
    return SeqFeature(
        location=SimpleLocation(start, end, strand=strand),
        type=ftype,
        qualifiers={"gene": ["intergenic"]},
    )


def find_intergenic_features(cdss: list[SeqFeature], seq_length, buffer_size=50):
    if not cdss:
        return []

    cdss.sort(key=lambda x: x.location.start)
    igs = []

    cur_end = cdss[0].location.end
    for cds in cdss[1:]:
        if cds.location.start > cur_end:
            intergenic_length = cds.location.start - cur_end - 1
            ftype = "intergenic" if intergenic_length > 0 else "junction"
            igs.append(
                create_intergenic_feature(
                    max(cur_end - buffer_size + 1, 0),
                    min(cds.location.start + buffer_size, seq_length),
                    strand=cds.location.strand,
                    ftype=ftype,
                )
            )
        cur_end = max(cur_end, cds.location.end)

    return igs


def format_feat(feat: SeqFeature):
    return f"{feat.type} {feat.location.start}-{feat.location.end} (len {len(feat)})"


def create_artificial_utrs(cdss: list[SeqFeature], seq_length, buffer_size=50):
    cds_end = max((cds.location.end for cds in cdss), default=0)
    if len(cdss) and cds_end < seq_length:
        res = [
            SeqFeature(
                location=SimpleLocation(
                    max(cds_end - buffer_size, 0), seq_length, strand=cdss[-1].location.strand
                ),
                type="inferred_3'UTR",
                qualifiers={"gene": ["3'UTR"]},
            )
        ]
        return res
    else:
        return []


def extract_features(rec: SeqRecord, buffer_size=50):
    feats = rec.features
    cdss = [feat for feat in feats if feat.type == "CDS"]

    # skip anything with negative strand CDS
    if any(feat.location.strand == -1 for feat in cdss):
        return []
    igs = find_intergenic_features(cdss, len(rec), buffer_size)
    utrs = [feat for feat in feats if feat.type == "3'UTR"]
    if len(utrs) == 0:
        utrs = create_artificial_utrs(cdss, len(rec), buffer_size)

    return igs + utrs


@click.command()
@click.option(
    "--buffer-size",
    type=int,
    default=50,
    help="Size on either end of every feature to include in the extracted sequence",
)
@click.option(
    "--full-seq-list",
    type=click.Path(exists=True, file_okay=True, path_type=Path),
    help="List of sequences for which no parsing is needed. The full sequences will be written to the output directory",
)
@click.option(
    "--max-feature-length", type=int, default=2000, help="Maximum length of a feature to extract"
)
@click.argument("input_dir", type=click.Path(exists=True, file_okay=False, path_type=Path))
@click.argument("output_dir", type=click.Path(exists=True, file_okay=False, path_type=Path))
def main(
    input_dir: Path,
    output_dir: Path,
    buffer_size: int,
    max_feature_length: int,
    full_seq_list: Path | None,
):
    full_seq_list = {id for id in full_seq_list.read_text().splitlines()} if full_seq_list else {}
    _ = [
        write_feat(rec, feat, file.stem, output_dir)
        for file in input_dir.glob("*.gbk")
        for rec in SeqIO.parse(file, "genbank")
        for feat in extract_features(rec, buffer_size)
        if file.stem not in full_seq_list
        and rec.annotations.get("molecule_type") not in ("DNA", "ds-RNA")
        and len(feat) <= max_feature_length
    ]

    _ = [
        SeqIO.write(rec, output_dir / f"{file.stem}.1-{len(rec) + 1}.complete.fa", "fasta")
        for file in input_dir.glob("*.gbk")
        for rec in SeqIO.parse(file, "genbank")
        if file.stem in full_seq_list
    ]


if __name__ == "__main__":
    main()
