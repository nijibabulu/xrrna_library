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

I downlooad this set of accessions into `data/sequences_initial.acc` (included in this repo) and then use 
```bash
python scripts/download.py --api-key=[API_KEY] data/sequences_v1.acc work/sequences_v1/
python scripts/parse.py work/sequences_v1/ work/utrs_igs_v3
``` 

to download the sequences. *Note:* you have to get your own API key from NCBI to use the download script. This is done by logging into your NCBI account and going to `Account Settings` and going to the section **API Key Management**.

Currently, `parse.py` works only on explicitly labeled UTRs and intergenic sequences, which represents a subset of the seed sequences.  The size of these sequences is explored in `plot.R`.

In addition, the statistics about virus genome coverage is obtained with `virus_stats.py`. We plot the coverage with `coverage.R`.

Finally we construct teh library with `construct_library.py`.

```bash
python scripts/construct_library.py work/utrs_igs_v3/ work/library_v1.fa
```
