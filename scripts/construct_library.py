#! /usr/bin/env python3

from pathlib import Path
import click

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def write_tiled_seqs(record: SeqRecord, f, oligo_length, tiling_length, seq5, seq3, file):
    seq = record.seq
    seq5 = Seq(seq5)
    seq3 = Seq(seq3)
    insert_length = oligo_length - len(seq5) - len(seq3)
    for i in range(0, len(seq) - oligo_length + tiling_length, tiling_length):
        insert = seq[i : i + insert_length]
        tiled_seq = seq5 + insert + seq3
        f.write(f">{record.id}_{i}_{i + insert_length} {file.stem}\n{tiled_seq}\n")


@click.command()
@click.option("--oligo-length", type=int, default=250, show_default=True)
@click.option("--tiling-length", type=int, default=80, show_default=True)
@click.option("--seq5", default="GATCTAATACGACTCACTATAGGATTAATATAAT", show_default=True)
@click.option("--seq3", default="AAAGAAACAACAACAACAAC", show_default=True)
@click.argument("input_dir", type=click.Path(file_okay=False, exists=True, path_type=Path))
@click.argument("output_fasta", type=click.Path(dir_okay=False))
def main(input_dir, output_fasta, oligo_length, tiling_length, seq5, seq3):
    with open(output_fasta, "w") as f:
        _ = [
            write_tiled_seqs(record, f, oligo_length, tiling_length, seq5, seq3, file)
            for file in input_dir.glob("*.fa")
            for record in SeqIO.parse(file, "fasta")
        ]


if __name__ == "__main__":
    main()
