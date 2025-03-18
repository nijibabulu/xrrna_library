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

I downlooad this set of accessions into `data/sequences_v2.acc` (included in this repo) and then use 
```bash
python scripts/download.py --api-key=[API_KEY] data/sequences_v2.acc work/sequences_v2/
python scripts/parse.py --full-seq-list data/tiled_genome_ids.txt work/sequences_v2/ work/utrs_igs_fulls_v4
``` 

to download the sequences. *Note:* you have to get your own API key from NCBI to use the download script. This is done by logging into your NCBI account and going to `Account Settings` and going to the section **API Key Management**.

Currently, `parse.py` works only on explicitly labeled UTRs and intergenic sequences, which represents a subset of the seed sequences.  The size of these sequences is explored in `plot.R`.

In addition, the statistics about virus genome coverage is obtained with `virus_stats.py`. We plot the coverage with `coverage.R`.

Finally we construct the library with `construct_library.py`.

```bash
python scripts/construct_library.py work/utrs_igs_fulls_v4/ work/library_v3.fa
python scripts/library_table.py work/library_v2.fa work/sequences_v2/ work/library_v2.csv
```

### Check for duplicates


Make a library with inserts only and construct an index of it
```bash
python scripts/construct_library.py --seq5='' --seq3='' --oligo-length=196 work/utrs_igs_fulls_v4/ work/library_v2_inserts_only.fa 
/opt/homebrew/Cellar/bowtie2/2.5.4/bin/bowtie2-build work/library_v2_inserts_only.fa work/index/library_v2_inserts_only
bowtie2 -a  -f -x work/index/library_v2_inserts_only -U work/library_v2_inserts_only.fa -S work/library_v2_inserts_only.sam
```

