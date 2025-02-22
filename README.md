# fuzzion2_patgen

This documentation is under construction.  

This code generates patterns for [Fuzzion2](https://www.github.com/stjude/fuzzion2/).  It can produce patterns for:

| sequencing type | event type | input |
| ---------- | ----- | ----- |
| RNA | fusion | sequence contig |
| RNA | fusion | genomic breakpoints |
| RNA | ITD/intragenic  | sequence contig |
| DNA | fusion | genomic breakpoints |

## Setup

Change to a directory where you want to keep the code, which will be referred to by $INSTALL_DIR below.

This command will retrieve a copy of the code and put it in a new "fuzzion2_patterns" subdirectory:
```
git clone https://github.com/stjude/fuzzion2_patterns.git
```

add the scripts directory to your PATH, and the Perl library directory to your PERL5LIB:
```
export PATH=$INSTALL_DIR/fuzzion2_patterns/src/scripts:$PATH
export PERL5LIB=$INSTALL_DIR/fuzzion2_patterns/src/perllib:$PERL5LIB
```


#### Dependencies

* Perl (version 5.10.1 or later)
* The following third-party Perl modules are required (this list likely needs work):
  * Set::IntSpan
  * LWP
  * Data::Compare
* BLAST (specifically the "blastn" executable), which must be available on your PATH


#### Installation test

To verify that the Perl code is runnable, execute the following command:
```
perl -cw `which fusion_contig_extension.pl`
```

This should return a message saying "syntax OK".  If error messages appear, please report them to us (see Contact section).  A common reason for errors is one or more third-party Perl modules missing from in your installation.

## Examples

The examples below show how to generate fuzzion2 patterns targeting different event and data types.  The .sh scripts may need to be updated to reference:

* a FASTA file for your genome, e.g. "GRCh37-lite.fa".  This must also have an accompanying .fai index file index file (i.e. generated by "samtools faidx FASTA_FILE")
* ncbiRefSeq.txt, a table from the UCSC Genome Annotation Database, e.g. for hg19 this is https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ncbiRefSeq.txt.gz.  This file details the mappings between RNAs and the genome sequence.
* hgnc_complete_set.txt, the HGNC database in tab-delimited format.  See file "Complete set new TSV" available from https://www.genenames.org/download/.  This is a database of gene symbols that is used for disambiguation.

## Example: RNA fusion pattern from contig sequence

The directory test/rna_fusion_contig contains scripts to generate fuzzion2 patterns from the example input file example.tsv.  Output files are also provided.  The scripts are:

```
step1_cicero_no_config.sh
# input: CICERO-format record describing a BCR-ABL1 fusion, see file example.tsv
# modify "-refflat ncbiRefSeq.txt" to point to your ncbiRefSeq.txt 
# modify "-hgnc hgnc_complete_set.txt" to point to your copy
# modify "-fasta GRCh37-lite.fa" to point to your genome FASTA file
step2_convert.sh
# processes intermediate file to yield fuzzion2 pattern file (example.tsv.extended.tab.fuzzion_extended_500.tab)
```

## Example: RNA fusion pattern from genomic coordinates

The directory test/rna_fusion_coordinates is similar to the previous example, but example.tsv specifies genomic breakpoints of the fusion.
```
step1_rna_no_config.sh
# modify "-refflat ncbiRefSeq.txt" to point to your local copy
# modify "-hgnc hgnc_complete_set.txt" to point to your copy
step2_convert.sh
```

## Example: RNA ITD or intragenic event pattern from contig sequence

The directory test/rna_intragenic contains an example of generating a pattern from a provided FLT3 ITD contig sequence.  This code generates a fuzzion2 pattern file directly rather than via an intermediate file.
```
step1_intragenic_no_config.sh
# modify "-refflat ncbiRefSeq.txt" to point to your local copy
# modify "-hgnc hgnc_complete_set.txt" to point to your copy
```

## Example: DNA fusion pattern from genomic coordinates

The directory test/dna_fusion_coordinates contains an example of generating a DNA pattern from provided breakpoint coordinate and strand information (see "test.tsv"). 

```
step1_genomic_no_config.sh
# modify "-refflat ncbiRefSeq.txt" to point to your local copy
# modify "-fasta GRCh37-lite.fa" to point to your local copy
step2_convert.sh
```


## Contact
Please contact Michael Edmonson <michael.edmonson@stjude.org> for assistance with this code.
