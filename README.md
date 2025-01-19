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

I downlooad this sequence into data/sequences_initial.acc and then use 
```bash
python scripts/download.py data/sequences_v1.acc work/sequences_v1/
python scripts/parse.py work/sequences_v1/ work/utrs_v1/
``` 

to download the sequences.

Initially, `parse.py` works only on explicitly labeled UTRs which represents a subset of the seed sequences. 
