#! /usr/bin/env python3

import collections
from pathlib import Path
import random
import click

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


AMBIG_SEQS = {
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "S": ["G", "C"],
    "W": ["A", "T"],
    "K": ["G", "T"],
    "M": ["A", "C"],
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"],
    "N": ["A", "C", "G", "T"],
}


def process_insert(seq: str, max_ambig: int):
    counter = collections.Counter(seq)
    ambig_count = sum(counter[c] for c in AMBIG_SEQS)
    if ambig_count <= max_ambig:
        return "".join(random.choice(AMBIG_SEQS.get(c, [c])) for c in seq)
    else:
        return None


def write_tiled_seqs(
    record: SeqRecord, f, oligo_length, tiling_length, seq5, seq3, file, max_ambig
):
    seq = record.seq
    seq5 = Seq(seq5)
    seq3 = Seq(seq3)
    insert_length = oligo_length - len(seq5) - len(seq3)
    for i in range(0, len(seq) - oligo_length + tiling_length, tiling_length):
        insert = process_insert(seq[i : i + insert_length], max_ambig)
        if insert is not None:
            tiled_seq = seq5 + insert + seq3
            f.write(f">{record.id}_{i}_{i + insert_length} {file.stem}\n{tiled_seq}\n")


@click.command()
@click.option("--oligo-length", type=int, default=250, show_default=True)
@click.option("--tiling-length", type=int, default=80, show_default=True)
@click.option("--seq5", default="GATCTAATACGACTCACTATAGGATTAATATAAT", show_default=True)
@click.option("--seq3", default="AAAGAAACAACAACAACAAC", show_default=True)
@click.option("--max-ambig", type=int, default=5, show_default=True)
@click.option("--random-seed", type=int, default=42, show_default=True)
@click.argument("input_dir", type=click.Path(file_okay=False, exists=True, path_type=Path))
@click.argument("output_fasta", type=click.Path(dir_okay=False))
def main(input_dir, output_fasta, oligo_length, tiling_length, seq5, seq3, max_ambig, random_seed):
    random.seed(random_seed)
    with open(output_fasta, "w") as f:
        _ = [
            write_tiled_seqs(record, f, oligo_length, tiling_length, seq5, seq3, file, max_ambig)
            for file in input_dir.glob("*.fa")
            for record in SeqIO.parse(file, "fasta")
        ]


if __name__ == "__main__":
    main()
