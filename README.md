# fuzzion2_patterns

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

## Contact
Please contact Michael Edmonson <michael.edmonson@stjude.org> for assistance with this code.
