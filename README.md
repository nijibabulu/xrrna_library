# xrRNA proj

## Introduction

The initial goal is downloading a number of viral sequences from NCBI and creating a strategy to tile non-coding parts of viral genomes and using those to create sequences to be ordered as RNA oligos from Twist. 

Our initial strategy for filtering the NCBI database was

* Use only RefSeq genomes (this is probably taking out too much data but it's a start)
* Use forward strand viral genomes (this criterion will likely persist)
* Require that they have been marked as complete assemblies
* Require that the nucleotide sequencing is also complete
* Require the host be a eukaryote.

The filtering is detailed [here](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&GenomicMoltype_s=ssRNA(%2B)&HostLineage_ss=Eukaryota%20(eucaryotes),%20taxid:2759&SourceDB_s=RefSeq&GenomeCompleteness_s=complete&Completeness_s=complete).

Some feature genomes will be totally covered. Those are:

- Zika virus - NC_012532 https://www.ncbi.nlm.nih.gov/nuccore/NC_012532.1
- SARS-CoV-2 - NC_045512 https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2
- Poliovirus - NC_002058 https://www.ncbi.nlm.nih.gov/nuccore/NC_002058.3 (Not in the above accession list)
- Norovirus - NC_001959 https://www.ncbi.nlm.nih.gov/nuccore/NC_001959.2 (Not in the above accession list)
- Sindbis virus - NC_001547 https://www.ncbi.nlm.nih.gov/nuccore/NC_001547.1 (Not in the above accession list)
- Chikungunya virus - NC_004162 https://www.ncbi.nlm.nih.gov/nuccore/NC_004162.2 (Not in the above accession list)

These IDs are in the file `data/tiled_genome_ids.txt`.

The download script requires a set of accession such as those in `data/sequences_v2.acc` (included in this repo). From there the sequences can be downloaded and parsed into feature-level sequences with these scripts:
```bash
python scripts/download.py --api-key=[API_KEY] data/sequences_v2.acc work/sequences_v2/
python scripts/parse.py --full-seq-list data/tiled_genome_ids.txt work/sequences_v2/ work/utrs_igs_fulls_v4
``` 

to download the sequences. *Note:* you have to get your own API key from NCBI to use the download script. This is done by logging into your NCBI account and going to `Account Settings` and going to the section **API Key Management**.

`parse.py` will construct 3' UTRs and infer 3' UTRs in cases in which none are found.  The size of these sequences can be explored in `plot.R`.

In addition, the statistics about virus genome coverage is obtained with `virus_stats.py`. We plot the coverage with `coverage.R`.

Finally we construct the library with `construct_library.py`.

```bash
python scripts/construct_library.py work/utrs_igs_fulls_v4/ work/library_v4.fa
python scripts/construct_library.py --max-ambig 0 work/utrs_igs_fulls_v4/ work/library_v5.fa
python scripts/library_table.py work/library_v2.fa work/sequences_v2/ work/library_v2.csv
```

### Check for duplicates


Make a library with inserts only and construct an index of it
```bash
python scripts/construct_library.py --seq5='' --seq3='' --oligo-length=196 work/utrs_igs_fulls_v4/ work/library_v2_inserts_only.fa 
/opt/homebrew/Cellar/bowtie2/2.5.4/bin/bowtie2-build work/library_v2_inserts_only.fa work/index/library_v2_inserts_only
bowtie2 -a  -f -x work/index/library_v2_inserts_only -U work/library_v2_inserts_only.fa -S work/library_v2_inserts_only.sam
```


Flavi-only Library:
```bash
uv run scripts/download.py --api-key=$API_KEY data/20260116_flaviAccessionList_v1.acc work/sequences_20260116_flaviAccessionList_v1/
uv run python scripts/parse.py work/sequences_20260116_flaviAccessionList_v1/ work/utrs_igs_20260116_flavi_v1/
uv run python scripts/construct_library.py work/utrs_igs_20260116_flavi_v1/ work/library_20260116_flavi_v1.fa
```

## Mifsud et al and Simmonds et al data

The Mifsud et al data can be downloaded from [https://doi.org/10.1038/s41586-024-07899-8](https://doi.org/10.1038/s41586-024-07899-8) Supplementary Table 1. And the Simmonds et all data can be downloaded from [https://doi.org/10.1038/s41564-025-02134-0](https://doi.org/10.1038/s41564-025-02134-0). 

Jeanine modified the Mifsud et al data as follows:

> I removed the TOMB outgroup sequences, filled in the accessions for the novel flavi sequences, and removed the Hou et al. sequences. There were 12 "hou et al." sequences, all from a metagenomic sampling from water and soil sources. I will check for more information and if they have an actual associated accession, but I haven't seen anything yet. It's probably ok to eliminate these without worrying we lose information. 

The sequences were retrieved as follows:

```bash
uv run scripts/parse_mifsud_simmonds.py data/Mifsud_2024_Supplementary_table_1_sequence_metadata.xlsx data/Simmonds_2025_Supplementary_Table_5.xlsx  > data/Mifsud_Simmonds.accs
```

Additional accessions from `data/additional_accs.txt` were also included.

```bash
grep GenBank data/additional_accs.txt | awk '{print $2}' > data/additional_accs.acc
cat data/Mifsud_Simmonds.accs data/additional_accs.acc > data/Mifsud_Simmonds_etc.accs
```

Noted that there was a sequence called "novel" included there.

Download the sequences and make the library 

```bash
uv run scripts/download.py --api-key=$API_KEY data/Mifsud_Simmonds_etc.accs work/sequences_mifsud_simmonds_EDIT/
uv run python scripts/parse.py --exclude-igs work/sequences_mifsud_simmonds_EDIT/ work/utrs_mifsud_simmonds_EDIT/
```

In order to discover the reason for exclusion of some sequences from teh final UTRs above, we will compare the accessions in `data/Mifsud_Simmonds_etc.accs` with the accessions in `work/utrs_mifsud_simmonds_EDIT/` and see which ones are missing.

```bash
ls work/utrs_mifsud_simmonds_EDIT/ | awk -F . '{print $1}' | sort -u > work/utrs_mifsud_simmonds_EDIT.accs
ls work/sequences_mifsud_simmonds_EDIT/ | awk -F . '{print $1}' | sort -u > work/sequences_mifsud_simmonds_EDIT.accs
comm work/sequences_mifsud_simmonds_EDIT.accs work/utrs_mifsud_simmonds_EDIT.accs
```

```bash
uv run python scripts/construct_library.py work/utrs_mifsud_simmonds_EDIT/ work/library_mifsud_simmonds_EDIT_tile80.fa
uv run python scripts/construct_library.py --tiling-length 40 work/utrs_mifsud_simmonds_EDIT/ work/library_mifsud_simmonds_EDIT_tile40.fa
uv run python scripts/construct_library.py --tiling-length 20 work/utrs_mifsud_simmonds_EDIT/ work/library_mifsud_simmonds_EDIT_tile20.fa
```

Now, take a look at an oligo-less library to see if there are any duplicates. 

```bash
uv run python scripts/construct_library.py --seq5='' --seq3='' --oligo-length=196 work/utrs_mifsud_simmonds_EDIT/ work/library_mifsud_simmonds_EDIT_inserts_only.fa 
cd-hit-est -i work/library_mifsud_simmonds_EDIT_inserts_only.fa -c 0.99 -n 10 -o work/library_mifsud_simmonds_EDIT_inserts_only.cdhit.fa 
```

With the ad-hoc script parse_cdhit_clusters.py, we can parse the clusters and outptu the duplicated accessions in `work/mifsud_simmonds_duplicate_seqs.csv`. 


Do the same with the full UTR sequences:

```bash
cat work/utrs_mifsud_simmonds_EDIT/*.fa >  work/utrs_mifsud_simmonds_EDIT.fa
cd-hit-est -i work/utrs_mifsud_simmonds_EDIT.fa -c 0.99 -n 10 -o work/utrs_mifsud_simmonds_EDIT.cdhit.fa
uv run scripts/parse_cdhit_clusters.py work/utrs_mifsud_simmonds_EDIT.cdhit.fa.clstr work/utrs_mifsud_simmonds_EDIT.cdhit.fa.dups.csv work/utrs_mifsud_simmonds_EDIT.cdhit.fa.for_removal.txt
```

Remove the accessions from work/utrs_mifsud_simmonds_EDIT.cdhit.fa.for_removal.txt before making the library. TODO: check if blacklisting works.

```bash
uv run python scripts/construct_library.py --blacklist work/utrs_mifsud_simmonds_EDIT.cdhit.fa.for_removal.txt --tiling-length 20 work/utrs_mifsud_simmonds_EDIT/ work/library_mifsud_simmonds_EDIT_tile20.dupsremoved.fa
```


