from pathlib import Path
import ncbi_acc_download.core
import click
from ncbi_acc_download.errors import (
    DownloadError,
    InvalidIdError,
)


@click.command()
@click.option("--molecule", default="nucleotide", help="Molecule type to download")
@click.option("--format", default="genbank", help="Format to download")
@click.option("--api-key", help="NCBI API key", default="none")
@click.argument("accession_list", type=click.Path(exists=True, path_type=Path))
@click.argument("output_dir", type=click.Path(exists=True, path_type=Path))
def main(accession_list, output_dir, molecule, format, api_key):
    config = ncbi_acc_download.core.Config(molecule=molecule, format=format, api_key=api_key)
    with open(accession_list) as f, open(output_dir / "err.txt", "w") as log:
        for line in f:
            accession = line.strip()
            out_path = output_dir / f"{accession}"
            try:
                ncbi_acc_download.core.download_to_file(accession, config, out_path)
            except InvalidIdError as _:
                print(f"Invalid ID: {accession}")
                log.write(f"{accession}\tinvalid ID")
            except DownloadError as _:
                print(f"Download error: {accession}")
                log.write(f"{accession}\tdownload error")


if __name__ == "__main__":
    main()
