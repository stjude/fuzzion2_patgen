#!/usr/bin/env perl
# given 2 genes and a fusion contig, provide extended contig
# MNE 7/2020 -
# TO DO:
#  - maybe detect genes if not provided?

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;
use File::Copy;
use Cwd qw(realpath);

use Digest::MD5 qw(md5_hex);
use Set::IntSpan;
use List::Util qw(min max sum);
use Bio::SeqIO;

use MiscUtils qw(dump_die build_argv_list get_hash_option unquote field_name_dwim unique_ordered_list);
use DelimitedFile;
use Reporter;
use ConfigUtils qw(config_or_manual);
use RefFlatFile;
use GeneSymbolMapper qw(new_gsm_lite);
use FAI qw(build_fai);
use BLASTer qw(dustmasker);
use BLATer;
use GenomeUtils qw(reverse_complement cook_chromosome_name);
use FileUtils qw(find_binary read_simple_file universal_open);
use SJPreferredIsoform;
use TdtConfig;
use Counter;
use GeneAnnotation;
use TranscriptFlanking;
use TemporaryFileWrangler;
use EUtilities;
use DuplicateTracker;
use DWIMColumns;
use GenBankFASTACache;
use GenBankCache;
use GenBankCDS;
use SortedColumnStreamer;
use ReadthroughFusionAnnotator;
use GeneListMatcher;
use HSPIndexer;
use DuplicationAlignment;

my %FLAGS;

my $VERBOSE = 0;
my $BLAT_MIN_SCORE;
# default probably fine, minScore 20 only required for some small AA matches

#my $EXTENDED_CHUNK_LENGTH = 500;
#my $EXTENDED_CHUNK_LENGTH = 750;
# now steve's preferred
my $EXTENDED_CHUNK_LENGTH = 500;
# 7/2021: back to the original spec

my $MASK_CONTIG_BEFORE_ALIGN = 1;

my $MASK_RESERVE_GENE_B = 0;
# whether to reserve a portion of the contig for a geneB match.  this
# was originally added to help with intragenic/ITD events, however
# those are now dealt with separately.  For events involving two
# different genes this may actually interfere with anchoring genes A
# and B.

my $GENEA_ANCHOR_MASK_FRACTION_FOR_B = 0.40;
# when anchoring geneA, what portion of the contig on the geneB side
# to mask before BLAT.  Want to capture as much as possible of geneA
# without introducing ambiguity for intragenic fusions.  Originally, slightly
# less than 50% of contig seemed about right.  However intragenic
# events are now handled differently and putative borders for ITDs
# and deletions now explicitly marked for masking.

my $F_GENE_A = "geneA";
my $F_GENE_B = "geneB";
my $F_CONTIG = "contig";
# default = CICERO-style
my $F_FUSION;
my $F_SAMPLE = "sample";
my $F_SOURCE = "source";

my $F_ITD_FEATURES = "itd_features";
# field only used by ITD/intragenic processing type

my ($F_CHR_A, $F_POS_A, $F_ORT_A,
    $F_CHR_B, $F_POS_B, $F_ORT_B);

my $TAG_ERROR_GENES_REVERSED = "ERROR_geneA_and_geneB_reversed_in_contig";

#my $CODING_ONLY = 1;
my $CODING_ONLY = 0;

my $MIN_INTERGENIC_OVERLAP_TO_MASK = 4;
# some small amount of random overlap possible due to microhomology,
# and aren't critical.  However longer overlaps can be problematic
# as they can interfere with fuzzion2's evidence gathering for both genes

my $DMZ_OVERLAP = 1;

my $LOCK_ACCESSION_FOR_INTRAGENIC = 1;
# for events within a transcript, don't combine all possible
# isoforms for genes A and B, just work through the pairs of
# the same transcript on both ends

my $MAX_NT_TO_CONSIDER_QUERY_FLUSH = 3;
# for ITDs when comparing the full-length gene sequence to the contig,
# if multiple segments for the full-length gene, max distance between
# end of one and start of next to consider contiguous (for ITDs,
# this is often 0 but there are examples of 1 nt, e.g. FLT3-FLT3)
#
# OBSOLETE, this is used in blat-based implementation onlya

#my $MAX_DUPLICATE_BLOCK_DISTANCE = 6;
#my $MAX_DUPLICATE_BLOCK_DISTANCE = 20;
my $MAX_DUPLICATE_BLOCK_DISTANCE;
# similar to the above: there is sometimes a bit of interstitial sequence
# between the duplicated sequences.  However sometimes the duplicated
# sequences are not adjacent, e.g.
# /home/medmonso/work/steve/2020_06_13_cicero_contig_extension/set2_2020_09_25/28
# - 2 HSPs
# - start of each query HSP identically matches contig in 2 places
# - however full hit for first HSP spans ~50% of contig, and
#   second contig match is after
# - first HSP matches both gene and contig, so the portion after the
#   duplicate region can't be considered interstitial
# - so, if distance between duplicated regions exceeds this threshold,
#   use the blast hit boundary instead
#
#
# originally used a small number like 6 nt, however arguably this should
# higher to introduce a longer gap in the final fuzzion2 output.
# e.g. this case:
# /home/medmonso/work/steve/2020_06_13_cicero_contig_extension/2020_09_25_set2/26
# gap is 17 nt, however contig is short and 2nd dup copy is at end.
# preserving the gap in the output will help fuzzion2 span the troublesome
# area, even though technically the "interstitial" sequence here is
# actually part of the canonical gene sequence.
#
# => final ruling: disable this feature, allow even very long gaps
#    as these must be spanned by fuzzion2 reads to confirm event

#my $MIN_INSERTION_DUPLICATE_MATCHES = 13;
#my $MIN_INSERTION_DUPLICATE_MATCHES = 9;
my $MIN_INSERTION_DUPLICATE_MATCHES = 7;
# for a region inserted into the contig (i.e. not present in the refseq),
# if this sequence doesn't perfectly repeat, minimum number of matches
# required to consider a portion of the sequence the duplicated region
#
# this example is only 9, maybe need a better approach here to avoid exception:
# /home/medmonso/work/steve/2020_06_13_cicero_contig_extension/2021_02_03_clingen_all/all/crash/test.tab


my $MIN_QUERY_GAP_TO_CONSIDER_ITD_OR_INSERTION = 5;
# minimum sequence length for putative ITDs/insertions.
# this may be too aggressive!
my $MIN_QUERY_GAP_TO_CONSIDER_ITD_CRITICAL = 10;
# if evidence can't be found for a repeat, consider spurious unless
# the gap is of a minimum size.

my $HIT_TRIM_END_RANGE = 10;
# trim hit boundaries of any gappy regions within this many nt of edge

#my $BLAST_MIN_PERCENT_IDENTITY = 90;
my $BLAST_MIN_PERCENT_IDENTITY = 85;
# reduce to 85%: one ETV6-RUNX1 fusion has a deletion vs. the reference
# which only matches at 89%:
# /home/medmonso/work/steve/2020_06_13_cicero_contig_extension/2021_01_04_all/blast/fail_coding/ETV6_RUNX1/one.tab
my $BLAST_MIN_WORD_SIZE = 20;
# BLAST 2.9.0 doesn't work at all with the SV below with the default settings:
# /home/medmonso/work/steve/2020_06_13_cicero_contig_extension/2021_01_04_all/blast/CICERO_paper/word_size_23.tab
# This does work if -word_size is specified and with a value of <= 23.
# Not sure what the default value is, some documentation seems to
# claim a value of 11, however I don't think this is still true
# because specifying a value of 11 does work (interestingly) but not
# specifying -word_size at all does not.

#my $BLAST_GAP_OPEN;
#my $BLAST_GAP_EXTEND;
my $BLAST_GAP_OPEN = 3;
my $BLAST_GAP_EXTEND = 2;
# optional parameters to blastn for gap open/extend penalties.
# See test_UBTF_gapopen3_gapext_2.tab: example where minimum 3/2 seems
# like a good idea to suppress an interfering gap between the query
# and the subject.  Setting a minimum helps suppress noise near edges
# which are especially problematic.  While the code can automatically
# resolve some of these problems, large-ish gaps may still cause
# problems and it seems better to suppress them during alignment

my $LARGE_ITD_ADD_FLANKING = 50;

my $CLEAN_SHARP_GENE_SYMBOLS = 1;

my $WGS_CONTIG_MAX_LEN = 1500;
my $WGS_CONTIG_BASE_FLANKING = 60;
# when generating contigs from WGS coordinates, how many bases on
# either side of breakpoint to use to generate basic contig

my $RNA_CONTIG_BASE_FLANKING = 80;
# same option for RNA

my $RNA_ORT_REQUIRED = 1;
# if orientation annotations are provided, whether or not they are mandatory.
# useful for sources where the annotation is not always provided.
#my $RNA_INTRAGENIC_MIN_DISTANCE_TO_PROCESS = 10000;

my $RNA_SUPPRESS_INTRAGENIC = 1;
# If set, suppress pattern generation entirely for intragenic events
# in RNA mode.  Now recommended: without knowing specific isoform to
# process this seems likely to create chaotic output:
#  - 9/2021 testing
#  - EGFRvIII which requires a specific isoform to get the result
#    matching the reported event (also this event has a dedicated
#    refSeq isoform)

my $RNA_INTRAGENIC_MIN_DISTANCE_TO_PROCESS = 2500;
# if RNA intragenic events are not suppressed completely, events where
# distances are shorter than this are considered unprocessable as they
# may be ITDs (these only currently work in contig mode)

my $RNA_SUPPRESS_DISTANT_ISOFORMS = 1;
# if a coordinate is near to exons in some isoforms but distant in
# others, only use the nearby isoforms
my $RNA_DISTANT_THRESHOLD = 20;
# TBD

my $F_PATTERN_ID = "pattern";
my $F_PATTERN_SEQUENCE = "sequence";
# Steve Rice 5/10/2021 email

my $PATTERN_CONDENSE_DUPLICATE_IF_FLANKING_WINDOW_IDENTICAL = 1;
my $PATTERN_CONDENSE_DUPLICATE_FLANKING_WINDOW_SIZE = 400;
# consider two patterns identical if they are identical (or nearly so)
# within specified flanking window.  Given insert size considerations,
# differences further away are unlikely to affect the results (Steve
# Rice analysis, email of 7/14/2021).

my $PATTERN_CONDENSE_MIN_QUERY_OVERLAP_FRAC = 0.99;
# how much of the query must overlap the target
my $PATTERN_CONDENSE_MIN_BLAST_FRAC_IDENTICAL = 0.99;
# fraction identical to consider a duplicate

my $PATTERN_CONDENSE_MULTI_HSP_MAX_QUERY_HIT_SPAN_SIZE_DIFFERENCE = 0.02;
# if multiple HSPs are involved in a match, don't allow their span sizes
# to be very different, i.e. skipping over non-aligned regions in one
# of the sequence that are not detailed in the HSP.
# for example, this file:
# /home/medmonso/work/steve/2020_06_13_cicero_contig_extension/2021_01_04_all/blast/RUNX1-RUNX1T1.tab
# - nearly all of one pattern aligns well with the other over 2 HSPs
# - HOWEVER, the interstitial region between pattern breakpoints
#   is quite different between the sequences.  It is therefore not
#   appropriate to consider one a duplicate of the other.
#
# potential tripwire: what if one HSP alone is a good enough match,
# and second smaller HSP is essentially just a distraction?

my $PATTERN_CONDENSE_MAX_GAP = 10;
# hard limit on gap size (largest gap, not total of all gap nt

my $PATTERN_CONDENSE_MAX_BREAKPOINT_DRIFT_FUSION = 5;
my $PATTERN_CONDENSE_MAX_BREAKPOINT_DRIFT_ITD = 0;

my $TRANSCRIPTOME_AMBIG_IGNORE_NON_CODING = 1;
my $TRANSCRIPTOME_AMBIG_IGNORE_READTHROUGH = 1;
my $AMBIG_FRACTION_WARNING = 0.50;
# TO DO: option to drop rows or write to rejection file if problematic
# ambiguity detected

my $DEEP_INTRONIC_MIN_CONTIG_SIDE_FRACTION = 0.70;
# if the given position is intronic in the context of the working
# isoform, distance from exon to consider deeply intronic.  This is
# half the contig length (i.e. the portion for gene A or B) times this
# multiplier.  So e.g. for a contig of 100 nt, 100 * 0.5 * 0.70 = 35
# nt: if the exon boundary is further than this a tag is added to
# the notes field.  Intended to mark cases where a contig may be
# too far away from an exon to be usable because it does not contain
# any coding sequence (e.g. a WGS CREST contig being used to target an
# RNA pattern).

my $RNA_ADD_GENE_ANNOTATIONS = 1;
# add/update genes and strand info if necessary
my $RNA_KEEP_GENE_ANNOTATIONS = 1;
# by default, keep any existing gene and strand annotation.
# the strand annotation in particular is critical information, if available.

my $RNA_FIX_GENE_ANNOTATIONS = 1;
# attempt to repair unrecognized/invalid gene symbols, e.g.
# Excel damage like "5-Sep" for "SEPT5"
my $RNA_UPDATE_GENE_STRAND = 0;
# while (re)annotating, also modify strand annotations.
# Typically we probably DON'T want to do this as there's no guarantee
# the fusion will be on the same sense as the gene model.

my $DEMUNGE_CHR_23 = 1;
# sometimes X appears as "chr23" in user data, likewise Y=24

my $REPORT_SJPI = 1;

my $READTHROUGH_MAX_DISTANCE = 200000;


#
#  parameters for 9/2023 ITD + DUP refactor:
#
my $ITD_AMBIGUITY_RESCUE_MAX_MICROHOMOLOGY_SIZE = 5;
# has been observed.  A few even higher, 6-7.  ?

my $ITD_MIN_SIZE = 10;
# TBD

my $ITD_RESCUE_MIN_SIZE = 7;
# see simple_itd_check()

#my $ITD_MAX_QUERY_GAP = 10;
#my $ITD_MAX_QUERY_GAP = 15;
my $ITD_MAX_QUERY_GAP = 55;
# TBD
# a bit of interstitial sequence is OK, but a lot may indicate
# intronic sequence in the contig, which implies an incompatible
# isoform, for example if the contig was based on an isoform
# with exons that are not present in the isoform being processed.
# e.g. FGFR1, NM_01354370.2 is missing some exons compare to
# NM_023110.3, so a contig based on the latter will not align
# well to the former.
#
# OTOH, this can sometimes include a chunk of non-refSeq sequence such
# as intronic sequence that forms a novel splice site (see
# test_large_subject_gap_JZ_wants_fusion_pattern_without_intronic.tab
# which has a 37-nt intronic inclusion).  These may be unique to some
# samples, in which case a fusion-type pattern (see below) may work
# better.  There may be no one-size-fits-all solution here.

my $ITD_LARGE_DUP_GENERATE_TWO_PATTERN_STYLES = 1;
# for large/medium dups, generate two pattern styles,
# traditional ITD and fusion

my $ITD_MIN_INTERSTITIAL_FOR_ALTERNATE_PATTERN = 8;
# for the above option, require a minimum interstitial sequence
# length to trigger generation of an alternate, fusion-style pattern.
# - if too low (or 0), we are essentially eliminating the overlap
#   requirement required by ITD-style patterns (sometimes very important
#   to help eliminate false positives)
# - if too high, an ITD-only pattern may not work for other samples
#   where the interstitial sequence is not present
# => exact threshold TBD!

my $ITD_MIN_BREAKPOINT_DISTANCE_FOR_ALTERNATE_PATTERN = 200;
# alternate patterns also require some significant distance
# between the breakpoints in the reference.  Otherwise small
# duplications can easily yield near-canonical fusion patterns
# which will produce false positive "strong-" fuzzion2 hits.
# example: test_KMT2A_small_dup_no_fusion.tab, 54 nt overlap
# of small exon12->13 duplication.

my $F_INTERNAL_ALT_PATTERN = "__alternate_pattern_generated__";
# internal tracking tag used in the above case

my $ITD_LARGE_DUP_TREAT_AS_FUSION = 1;
# if a large duplication is detected, generate a fusion pattern rather
# than an ITD pattern.  In these cases large query gaps and the
# resulting interstitial sequence is ignored, because it may be
# intronic sequence which will never align to the RNA record.
#
# The latter is an open question; might these ever be something
# important, e.g. a bit of retained intron which has some function
# during protein coding?
#
# In a small experiment, removing the interstitial sequence caught
# a few more samples, whereas the ITD-style pattern caught only
# the original sample.  The source sample had nowhere near as much
# supporting data than with the specialized pattern however.
# A hybrid pattern (interstitial preserved, but fusion-style brackets)
# may also be worth considering (see 10/27/2023 jz meeting notes)

my $ITD_BRACKET_BOTH_COPIES = 1;

my $ITD_BRACKET_ONE_COPY_IF_DUP_MINIMUM = 20;
# if we would normally bracket both copies, allow an exception to only
# bracket one copy if the duplication meets a minimum length.  This
# helps find evidence for long duplications in short-read sequencing,
# at the cost of requiring more user interpretation of the results.
# => some connection with what fuzzion2 considers a strong+ hit?  16 nt?

my $ITD_BRACKET_ONE_COPY_RESCUE = 1;
#my $ITD_BRACKET_ONE_COPY_RESCUE_MINIMUM = 60;
my $ITD_BRACKET_ONE_COPY_RESCUE_MINIMUM = 30;
# EXPERIMENTAL additional rescue:
# if the final total interstitial sequence is still long,
# fall back to using BLAST boundaries.  Can rescue e.g.
# /research/rgs01/home/clusterHome/medmonso/work/steve/2020_06_13_cicero_contig_extension/2021_01_04_all/blast/ITD_refactor/final_set/v2_fixed_annot/FLT3-FLT3-02/src.tab

my $ITD_PURE_INSERTION_MAX_SUBJECT_GAP = 3;
my $ITD_PURE_DELETION_MAX_QUERY_GAP = 3;

my $ITD_MIN_QUERY_ALIGN_FRACTION = 0.75;
#my $ITD_MIN_QUERY_ALIGN_FRACTION = 0.70;
#my $ITD_MIN_QUERY_ALIGN_FRACTION = 0.55;
# some lower alignment examples involve 5' UTR, so if that is
# used for alignment this becomes less of an issue.

# contig is expected to generally be a good match to RNA refSeq.
# If it isn't, that could indicate:
#  - inappropriate isoform being used for comparison (often works
#    fine with a different isoform, so this is OK)
#  - large amount of intronic sequence in contig, e.g.
#    test_FLT3_4_hsp_with_intronic_sequence.tab
#
# Leave some wiggle room in case query sequence contains a novel
# insertion, or there's a little bit of retained intron, etc.

my $ITD_MERGE_PATTERN_OUTPUT = 0;
# if there are multiple output rows for the same input row (e.g.
# for different transcripts) which result in the same pattern,
# report only one row with column values merged into comma-delimited lists.

my $ITD_COMPLEX_MIN_LENGTH = 9;
# if a possible complex sub detected, minimum query length affected

my $ITD_INCLUDE_UTR5_IN_MAPPING = 1;

#
#  preferences for fusion patterns:
#
#my $FUSION_MAX_INTERSTITIAL_LENGTH = 10;
# as with ITDs, long interstitial sequences can suppress legitimate
# matches in short-read sequencing due to overlap requirements.
# It's probably safer to trim these a bit more aggressively than
# intragenic events (especially ITDs).
my $FUSION_MAX_INTERSTITIAL_LENGTH = 0;
# 11/14/2024:
# disable this feature, for two reasons:
# 1. next version of fz2 will be able to skip over interstitial regions,
#    we can leave the full sequences in place to demark them.
# 2. current code doesn't modify/trim the sequence, it just shifts
#    the breakpoint so that the interstitial is no longer than X nt.
#    However this is a bad idea as the pattern can then simply fish out
#    unspliced RNA from total RNA data, leading to many false positives.
#

my $FUSION_INTERSTITIAL_TRIM_POLICY = "B";
# B: preserve boundary on B side, trim A side
# future: maybe consider different policies, e.g. center?

my $FUSION_PREFER_BEST_CONTIG_MATCHES = 1;
#my $FUSION_PREFER_BEST_CONTIG_MATCHES_BEST_RELATIVE_CUTOFF = 0.90;
my $FUSION_PREFER_BEST_CONTIG_MATCHES_BEST_RELATIVE_CUTOFF = 0.95;
# if the contig with the best alignment with the refseqs has X bases,
# minimum result to include is X * the above.  In other words tolerate
# a bit of dropoff in alignment coverage vs. the best match, but not much.
# - tighten from 90 to 95% to catch ETV6-RUNX1 SJBALL032392_D1,
#   where bogus "interstitial" result has 91% identity


my $F_GENOME = "genome";
# genome version associated with genomic coordinates given
my $F_PROCESSING_TYPE = "pattern_generation_method";
my $F_GENE_PAIR_SUMMARY = "gene_pair_summary";
my $F_GENE_PAIR_SUMMARY_EXTENDED = "gene_pair_summary_extended";
# including genomic
my $F_PATHOGENICITY_SOMATIC = "pathogenicity_somatic";

my $F_GENE_PAIR = "gene_pair";
# for fz2 grouping option

my @H_FZ2_CLONE = (
		   $F_GENOME,
		   qw(
		      genea_symbol
		      genea_acc
		      genea_feature
		      genea_chr
		      genea_pos
		      geneb_symbol
		      geneb_acc
		      geneb_feature
		      geneb_chr
		      geneb_pos
		    ),
		   $F_GENE_PAIR_SUMMARY,
		   $F_GENE_PAIR_SUMMARY_EXTENDED,
		   $F_PROCESSING_TYPE,
		   $F_PATHOGENICITY_SOMATIC,
		  );
# fields copied verbatim from fz2 intermediate file format into final
# pattern file.

my @H_FZ2_CORE = (
		  $F_PATTERN_ID,
		  $F_PATTERN_SEQUENCE,
		  # first 2 fields are fixed, required by fuzzion2
		  $F_GENE_PAIR,
		  $F_SOURCE,
		  # file patterns were merged from
		  $F_SAMPLE,
		  # TO DO: should this be a clone field?
		  # or maybe keep as is as there is some
		  # do-what-I-mean column procesing here
		  @H_FZ2_CLONE,
		 );
# core fz2 pattern file columns, other columns from earlier intermediate
# files are not passed through.

my @FUZZION_SRC_FILES;
my @FUZZION_SOURCE_APPEND_FIELDS;

my @SRC_FILES;
my @PATTERN_FILES;

use constant TAG_NO_CONTIG => "no_contig_given";
use constant TAG_DEEP_INTRONIC => "deep_intronic";
use constant TAG_ERROR_NO_INPUT_RESULTS => "ERROR_no_results_for_input_row";
use constant TAG_ERROR_NO_SEQUENCES_FOR_GENE => "ERROR_no_sequences_for_gene";
use constant TAG_ERROR_BLAST_HIT_MINUS => "ERROR_BLAST_hit_to_minus_strand";
use constant TAG_MISSING_GENE_ANNOTATION => "missing_gene_annotation";
use constant TAG_ADDED_GENE_ANNOTATION => "added_gene_annotation";
use constant TAG_REPLACED_GENE_ANNOTATION => "replaced_gene_annotation";
use constant TAG_UNRECOGNIZED_GENE_ANNOTATION => "unrecognized_gene_annotation";

use constant TAG_ADDED_STRAND => "strand_added";
use constant TAG_INCONSISTENT_STRAND_REPLACED => "inconsistent_gene_strand_replaced";
use constant TAG_INCONSISTENT_STRAND_DETECTED => "inconsistent_gene_strand_detected";
# detected a possible problem but left as-is; annotation might be legitimate,
# e.g. entisense event.

use constant BREAKPOINT_FUSION_LEFT => "]";
use constant BREAKPOINT_FUSION_RIGHT => "[";
use constant BREAKPOINT_ITD_LEFT => "}";
use constant BREAKPOINT_ITD_RIGHT => "{";

use constant VALUE_NA => "NA";

my %BREAKPOINT_CHAR = map {$_, 1} (
				   BREAKPOINT_FUSION_LEFT(),
				   BREAKPOINT_FUSION_RIGHT(),
				   BREAKPOINT_ITD_LEFT(),
				   BREAKPOINT_ITD_RIGHT()
				  );

my @CLINICAL_SRC_FILE_LISTS;
my @BLACKLIST_GENOMIC_FILES;
my @SAMPLE_SRC_FILES;

my @clopts = (
	      "-file=s",
	      "-genome=s",
	      "-blat-host=s",
	      "-blat-port=i",

	      "-refflat=s",
	      "-fasta=s",
	      "-blat-min-score=i" => \$BLAT_MIN_SCORE,
	      "-save-tempfiles",
	      "-verbose=i" => \$VERBOSE,

	      "-generate-test-file=s",
	      "-extended-length=i" => \$EXTENDED_CHUNK_LENGTH,

	      "-convert-fuzzion=s",
	      "-convert-fuzzion-src=s" => \@FUZZION_SRC_FILES,
	      # convert output for Steve Rice's use (fuzzion2)

	      "-source-append=s" => \@FUZZION_SOURCE_APPEND_FIELDS,

	      "-merge-pattern-ids-orig=s",
	      "-merge-pattern-ids-annotated=s",
	      # merge pattern ID annotations back to original file

	      "-f-gene-a=s" => \$F_GENE_A,
	      "-f-gene-b=s" => \$F_GENE_B,
	      "-f-fusion=s" => \$F_FUSION,
	      "-f-contig=s" => \$F_CONTIG,

	      "-f-chr-a=s" => \$F_CHR_A,
	      "-f-chr-b=s" => \$F_CHR_B,
	      "-f-pos-a=s" => \$F_POS_A,
	      "-f-pos-b=s" => \$F_POS_B,

	      "-repair-gene-swaps=s",

	      "-restrict-nm-a=s",
	      "-restrict-nm-b=s",
	      # debug: restrict processing to specific transcript for A/B

	      "-restrict-nm-sjpi",
	      "-report-sjpi=i" => \$REPORT_SJPI,
	      "-sjpi=s",

	      "-coding-only=i" => \$CODING_ONLY,
	      "-dmz-overlap=i" => \$DMZ_OVERLAP,

	      "-mask-contig-before-align=i" => \$MASK_CONTIG_BEFORE_ALIGN,

	      "-max-duplicate-block-distance=i" => \$MAX_DUPLICATE_BLOCK_DISTANCE,

	      "-generate-contigs=s",

	      "-min-blast-percent-identity=f" => \$BLAST_MIN_PERCENT_IDENTITY,

	      "-collate-cicero=s",
	      "-out=s",

	      "-mask-reserve-gene-b=i" => \$MASK_RESERVE_GENE_B,

	      "-hgnc=s",
	      # temporary hack during rdev deployment problems

	      "-diagnose-matches=s",
	      "-patterns=s",

	      "-genomic",
	      "-crest",
	      "-cicero",
	      "-no-chr-pos",
	      "-annotate-genes-as-breakpoints",

	      "-rna",
	      "-full-gene-pair=s",
	      "-full-gene-junctions=s",

#	      "-keep-gene-annotation",
	      "-rna-keep-gene-annotations=i" => \$RNA_KEEP_GENE_ANNOTATIONS,

	      "-rna-suppress-intragenic=i" => \$RNA_SUPPRESS_INTRAGENIC,
	      "-rna-ort-required=i" => \$RNA_ORT_REQUIRED,
	      "-rna-intragenic-min-distance=i" => \$RNA_INTRAGENIC_MIN_DISTANCE_TO_PROCESS,
	      "-rna-no-ort",
	      "-rna-add-gene-annotations=i" => \$RNA_ADD_GENE_ANNOTATIONS,
	      "-rna-fix-gene-annotations=i" => \$RNA_FIX_GENE_ANNOTATIONS,
	      "-rna-update-gene-strand=i" => \$RNA_UPDATE_GENE_STRAND,

	      "-blacklist=s",
	      # symbol-only list should be deprecated?
	      "-no-blacklist",

	      "-blacklist-genomic=s" => \@BLACKLIST_GENOMIC_FILES,

	      "-rna-flanking=i" => \$RNA_CONTIG_BASE_FLANKING,

	      "-annotate-sources",
	      "-annotate-pattern-exons=s",
	      # hsck
	      "-generate-blacklist",
	      # another hack

	      "-find-clinical-ids=s" => \@CLINICAL_SRC_FILE_LISTS,

	      "-condense=s",
	      # condense a merged pattern set
	      "-condense-no-blast",
	      "-condense-min-overlap-frac=f" => \$PATTERN_CONDENSE_MIN_QUERY_OVERLAP_FRAC,
	      "-condense-min-blast-frac=f" => \$PATTERN_CONDENSE_MIN_BLAST_FRAC_IDENTICAL,

	      "-condense-step2=s",
	      "-condense-preset=f",
	      "-condense-pattern-annot-glob=s",
	      "-condense-no-pattern-files",


	      "-summarize-samples=s" => \@SAMPLE_SRC_FILES,

	      "-find-gene-intervals=s",

	      "-generate-breakpoint-pairs=s",
	      "-step=i",

	      "-extract-svs=s",
	      "-intervals=s",

	      "-transcriptome-annotate=s",
	      "-transcriptome-fasta=s",
	      "-transcriptome-ambig-ignore-non-coding=i" => \$TRANSCRIPTOME_AMBIG_IGNORE_NON_CODING,
	      "-transcriptome-ambig-ignore-readthrough=i" => \$TRANSCRIPTOME_AMBIG_IGNORE_READTHROUGH,

	      "-readthrough-annotate=s",

	      "-readthrough-source-summary=s",

	      "-no-pos",
	      "-no-chr",

	      "-parse-cosmic=s",

	      "-extract-cicero=s",
	      "-steve-index=s",
	      # extract source CICERO calls for SVs using paths
	      # from Steve Rice's index files

	      "-supplemental-fasta=s",
	      "-add-supplemental-gene=s",
	      "-accession=s",

	      "-pattern2fasta=s",
	      "-id=s",

	      "-extract-pattern-src=s",
	      "-extract-patterns=s",
	      "-src=s",
	      "-substr",

	      "-coding-distance-counts=s",

	      "-find-missing-genes=s",

	      "-full-length-filter=s",

	      "-debug",
	      "-isofox=s",

	      "-gapopen=i" => \$BLAST_GAP_OPEN,
	      "-gapextend=i" => \$BLAST_GAP_EXTEND,

	      "-insertion2fuzzion=s",
	      "-itd2fuzzion=s",
	      # various intragenic events, not just ITDs
	      "-del2fuzzion=s",
	      # genomic deletion -> RNA pattern

	      "-itd-max-query-gap=i" => \$ITD_MAX_QUERY_GAP,
	      "-handle-large-dup-as-fusion=i" => \$ITD_LARGE_DUP_TREAT_AS_FUSION,
	      "-itd-merge-pattern-output=i" => \$ITD_MERGE_PATTERN_OUTPUT,
	      "-itd-min-query-align-fraction=f" => \$ITD_MIN_QUERY_ALIGN_FRACTION,
	      "-gene=s",

	      "-null-sample",

	      "-pattern-summary=s",

	      "-append-orig=s",
	      "-append-new=s",
	      # crude method to append a set of UNIQUELY-IDENTIFIED
	      # patterns into an existing file

	      "-add-preferred-key=s",
	      "-generate-transcript-info",

	      "-patch-breakpoints=s",

	      "-soft-clip-to-pattern",
	      "-sam=s",

	      "-merge-patterns",
	      "-f-patterns=s" => \@PATTERN_FILES,
	      # fz2.tab pattern files
	      "-f-src=s" => \@SRC_FILES,
	      # source files for fz2.tab pattern files
	      "-header-policy=s",
	      # "core" or "first"
	      "-patch=s",
	      # optionally patch fields with given key,value pairs
	      "-blank-undef",
	      "-append-source",
	      "-prepend-secondary",

	      #
	      # annotate exon info from fz2 hits file:
	      #
	      "-hit2exons=s",
	      "-read=s",
	      # filter to a specific read
	      "-sort-by-span-score=i",

	      "-itd-bracket-both-copies=i" => \$ITD_BRACKET_BOTH_COPIES,
	      "-itd-bracket-one-copy-if-dup-minimum=i" => \$ITD_BRACKET_ONE_COPY_IF_DUP_MINIMUM,

	      "-rc",

	      "-repeat-expansion-fz2",
	      # generate fuzzion2 patterns for regions flanking
	      # a repeat expansion site
	      "-repeat-expansion-ref",
	      # generate an excerpt of a reference sequence reflecting
	      # repeat expansion.  Used for custom mapping tests.
	      # seems iffy: mapper may try too hard to align everything here
	      "-repeat-expansion-whole",
	      # copy the whole genome file, adjusting one chromosome
	      # to reflect a repeat expansion.  Used for custom mapping tests.

	      "-chr=s",
	      "-start=i",
	      "-end=i",
	      "-repeat=s",
	      "-count=i",
	      "-count-list=s",
	      "-flank=i",

	      "-fusion-max-interstitial-length=i" => \$FUSION_MAX_INTERSTITIAL_LENGTH,
	      "-fusion-prefer-best-contig-matches=i" => \$FUSION_PREFER_BEST_CONTIG_MATCHES,

	      "-find-interstitial-coding",
	      "-pattern=s",
	      # optional filter

	      "-patch-features",
	      "-no-genbank-cache-check",
	      "-fuzzy-codon=i",
	      "-restrict-source=s",

	      "-patch-gene-pair-summary",
	      "-patch-gene-pair",

	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  ) || die "error parsing command-line options";

if ($FLAGS{cicero}) {
  $F_CHR_A = 'chrA';
  $F_POS_A = 'posA';
  $F_GENE_A = 'geneA';
  $F_ORT_A = 'ortA';

  $F_CHR_B = 'chrB';
  $F_POS_B = 'posB';
  $F_GENE_B = 'geneB';
  $F_ORT_B = 'ortB';
}

#
#  these utility options skip some standard checks, e.g. for required binaries:
#
if ($FLAGS{"gencode2refflat"}) {
  build_gencode2refflat();
  exit(0);
} elsif ($FLAGS{"extract-cicero"}) {
  extract_cicero();
  exit(0);
} elsif ($FLAGS{"add-supplemental-gene"}) {
  add_supplemental_gene();
  exit(0);
} elsif ($FLAGS{"pattern2fasta"}) {
  pattern2fasta();
  exit(0);
} elsif ($FLAGS{"extract-pattern-src"}) {
  extract_pattern_sources();
  exit(0);
} elsif ($FLAGS{"extract-patterns"}) {
  extract_patterns();
  exit(0);
} elsif ($FLAGS{"coding-distance-counts"}) {
  coding_distance_counts();
  exit(0);
} elsif ($FLAGS{"find-missing-genes"}) {
  find_missing_genes();
  exit(0);
} elsif ($FLAGS{"full-length-filter"}) {
  full_length_filter();
  exit(0);
} elsif ($FLAGS{"isofox"}) {
  # convert priesley "isofox" report
  isofox_convert();
  exit(0);
} elsif ($FLAGS{"insertion2fuzzion"}) {
  # treat contigs as containing insertion
  insertion2fuzzion();
  exit(0);
} elsif ($FLAGS{"itd2fuzzion"}) {
  # treat contigs as containing ITD
  itd2fuzzion();
  exit(0);
} elsif ($FLAGS{"pattern-summary"}) {
  pattern_summary();
  exit(0);
} elsif ($FLAGS{"append-orig"}) {
  append_patterns();
  exit(0);
} elsif ($FLAGS{"add-preferred-key"}) {
  add_preferred_key();
  exit(0);
} elsif ($FLAGS{"patch-breakpoints"}) {
  patch_breakpoints();
  exit(0);
} elsif ($FLAGS{"soft-clip-to-pattern"}) {
  soft_clip_to_pattern();
  exit(0);
} elsif ($FLAGS{"hit2exons"}) {
  hit2exons();
  exit(0);
} elsif ($FLAGS{"merge-patterns"}) {
  merge_pattern_files();
  exit(0);
} elsif ($FLAGS{"del2fuzzion"}) {
  # genomic deletion -> fuzzion2 pattern
  del2fuzzion();
  exit(0);
} elsif ($FLAGS{"repeat-expansion-fz2"}) {
  # generate WGS patterns to target repeat expansion
  repeat_expansion_fz2();
  exit(0);
} elsif ($FLAGS{"repeat-expansion-ref"}) {
  repeat_expansion_ref();
  exit(0);
} elsif ($FLAGS{"repeat-expansion-whole"}) {
  repeat_expansion_whole();
  exit(0);
} elsif ($FLAGS{"find-interstitial-coding"}) {
  find_interstitial_coding();
  exit(0);
} elsif ($FLAGS{"patch-gene-pair-summary"}) {
  patch_gene_pair_summary();
  exit(0);
} elsif ($FLAGS{"patch-gene-pair"}) {
  patch_gene_pair();
  exit(0);
}


my $F_EXTENDED_CHUNK = sprintf 'extended_%d', $EXTENDED_CHUNK_LENGTH;
my $F_NOTES = "fz2_notes";
# don't use "notes" as this conflicts with a header used in CICERO output

my @NEW_HEADERS = (
		   qw(
		       input_row_number
		       input_row_md5

		    ),
		   $F_GENOME,
		   qw(
		       genea_symbol
		       genea_acc
		       genea_acc_preferred
		       genea_pos_adj
		       genea_coding_distance
		       genea_feature

		       geneb_symbol
		       geneb_acc
		       geneb_acc_preferred
		       geneb_pos_adj
		       geneb_coding_distance
		       geneb_feature
		    ),
		   $F_GENE_PAIR_SUMMARY,
		   $F_GENE_PAIR_SUMMARY_EXTENDED,
		   $F_PROCESSING_TYPE,
		   qw(
		       genea_full
		       geneb_full
		       genea_contig_start
		       genea_contig_end
		       geneb_contig_start
		       geneb_contig_end

		       extended_exon
		    ),
		   $F_EXTENDED_CHUNK,
		   qw(
		       extended_max
		    ),
		   $F_NOTES
		  );
# TO DO:
# also include "sample" and other fields required by pattern conversion step?


printf STDERR "configuration:\n";
printf STDERR "  restricting to coding transcripts only?: %s\n", $CODING_ONLY ? "y" : "n";
printf STDERR "  restricting to preferred isoforms only?: %s\n", $FLAGS{"restrict-nm-sjpi"} ? "y" : "n";

#find_binary("blat", "-die" => 1);
find_binary("blastn", "-die" => 1);
# not needed for some functions below

if (my $info = $FLAGS{"generate-test-file"}) {
  generate_test_file($info);
  exit(0);
} elsif ($FLAGS{"convert-fuzzion"} or @FUZZION_SRC_FILES) {
  convert_fuzzion();
  exit(0);
} elsif ($FLAGS{"merge-pattern-ids-orig"}) {
  merge_fuzzion_ids();
  exit(0);
} elsif ($FLAGS{"repair-gene-swaps"}) {
  repair_gene_swaps();
  exit(0);
} elsif ($FLAGS{"generate-contigs"}) {
  generate_contigs();
  exit(0);
} elsif ($FLAGS{"collate-cicero"}) {
  collate_cicero_calls();
  exit(0);
} elsif ($FLAGS{"diagnose-matches"}) {
  diagnose_pattern_matches();
  exit(0);
} elsif ($FLAGS{"genomic"}) {
  # generate WGS contigs from coordinates/orientation only
  generate_genomic_contigs();
  exit(0);
} elsif ($FLAGS{"rna"}) {
  # generate RNA contigs from coordinates only
  generate_rna_contigs();
  exit(0);
} elsif ($FLAGS{"full-gene-pair"}) {
  generate_full_gene_pairs();
  exit(0);
} elsif ($FLAGS{"full-gene-junctions"}) {
  generate_full_gene_junctions();
  exit(0);
} elsif ($FLAGS{"annotate-sources"}) {
  annotate_sources();
  exit(0);
} elsif ($FLAGS{"annotate-pattern-exons"}) {
  annotate_pattern_exons();
  exit(0);
} elsif ($FLAGS{"condense"}) {
  condense_patterns();
  exit(0);
} elsif ($FLAGS{"condense-step2"}) {
  condense_from_report($FLAGS{"condense-step2"});
  exit(0);
} elsif (@CLINICAL_SRC_FILE_LISTS) {
  extract_clinical_ids();
  exit(0);
} elsif ($FLAGS{"generate-blacklist"}) {
  generate_blacklist();
  exit(0);
} elsif ($FLAGS{"find-gene-intervals"}) {
  find_gene_intervals();
  exit(0);
} elsif ($FLAGS{"generate-breakpoint-pairs"}) {
  generate_breakpoint_pairs();
  exit(0);
} elsif (@SAMPLE_SRC_FILES) {
  summarize_samples();
  exit(0);
} elsif ($FLAGS{"extract-svs"}) {
  extract_svs();
  exit(0);
} elsif ($FLAGS{"transcriptome-annotate"}) {
  transcriptome_annotate();
  exit(0);
} elsif ($FLAGS{"readthrough-annotate"}) {
  readthrough_annotate();
  exit(0);
} elsif ($FLAGS{"readthrough-source-summary"}) {
  readthrough_source_summary();
  exit(0);
} elsif ($FLAGS{"parse-cosmic"}) {
  parse_cosmic();
  exit(0);
} elsif ($FLAGS{"patch-features"}) {
  patch_features();
  exit(0);
}

my $SJPI;
if ($REPORT_SJPI or $FLAGS{"restrict-nm-sjpi"}) {
  init_sjpi();
}

my $rf = get_refflat();

my $FAI = get_fai();
# genome

my $gsm = init_gsm($rf);

my $ENABLE_SUPPLEMENTAL_REFSEQ;
# a library of sequences not in refFlat set, e.g. IGH
# or genes that are difficult to map in GRCh37, e.g. DUX4

my $ssl;
printf STDERR "  supplemental sequence library: ";
if (my $ncfa = $FLAGS{"supplemental-fasta"}) {
  printf STDERR "%s\n", $ncfa;
  $ssl = get_supplemental_cache();
  $ENABLE_SUPPLEMENTAL_REFSEQ = 1;
} else {
  printf STDERR "disabled\n";
}

my $infile = $FLAGS{file} || die "-file";

my $df = new DelimitedFile("-file" => $infile,
			   "-headers" => 1,
			     );
my $outfile = basename($infile) . ".extended.tab";

my $rpt = $df->get_reporter(
			    "-file" => $outfile,
			    "-extra" => \@NEW_HEADERS,
			    "-auto_qc" => 1,
			   );

my %rf2genomic;
my $input_row_number = 0;

my $blacklist = get_blacklist();

my $fai = get_fai();
my $tf = new TranscriptFlanking(
				"-fai" => $fai
			       );

my %suppressed;
while (my $row = $df->get_hash()) {
  my ($gene_a, $gene_b) = get_gene_a_b($row);
  $input_row_number++;

  $row->{input_row_md5} = md5_hex(map {$row->{$_}} @{$df->headers_raw});
  # if mux.pl was used for parallelization row numbers won't be unique,
  # however MD5 should be.
  # - MD5 is of the raw input row, before any tweaks (e.g. gene annotation)
  # TO DO: call populate_md5sum()

  my @gto_a;
  my @gto_b;
  if ($FLAGS{"no-chr"}) {
    @gto_a = @gto_b = ("-no-chr" => 1);
  } else {
    my $chr_a = clean_chr($row->{$F_CHR_A} || dump_die($row, "no data for $F_CHR_A. Chromosomes are recommended due to mapping ambiguity for some genes, specify -no-chr if not available"));
    my $chr_b = clean_chr($row->{$F_CHR_B} || die);
    # needed as hint for refFlat to find the correct record for gene
    # feature annotation, e.g. CRLF2 is mapped to both chr X and Y.
    @gto_a = ("-chr" => $chr_a);
    @gto_b = ("-chr" => $chr_b);
  }

  if ($ENABLE_SUPPLEMENTAL_REFSEQ) {
    push @gto_a, ("-ssl" => $ssl);
    push @gto_b, ("-ssl" => $ssl);
  }

  my $is_intragenic = $gene_a eq $gene_b;

  dump_die($row, sprintf("processing row %d: %s/%s", $input_row_number, $gene_a, $gene_b), 1) if $VERBOSE;

  my $can_process = 1;

  my @notes_general;
  my ($gene_a_transcript2rf, $gene_b_transcript2rf, %transcript2gene);
  if ($gene_a and $gene_b) {
    $gene_a_transcript2rf = get_transcripts($rf, $gsm, $gene_a, \@notes_general, \%transcript2gene, @gto_a);
    $gene_b_transcript2rf = get_transcripts($rf, $gsm, $gene_b, \@notes_general, \%transcript2gene, @gto_b);
    # mapping of transcripts to refFlat record entries

    transcript_filter($gene_a_transcript2rf, $gene_a, $FLAGS{"restrict-nm-a"});
    transcript_filter($gene_b_transcript2rf, $gene_b, $FLAGS{"restrict-nm-b"});

    unless (keys %{$gene_a_transcript2rf}) {
      $can_process = 0;
      push @notes_general, "ERROR_no_gene_A_transcripts";
    }

    if ($blacklist and is_blacklisted($blacklist, $row)) {
      $can_process = 0;
      log_error("blacklisted", \%suppressed, \@notes_general);
    }

  } else {
    $can_process = 0;
    push @notes_general, "ERROR_both_gene_A_and_B_required";
  }

  unless ($can_process) {
    # processing not possible, skip
    my %r = %{$row};
    $r{$F_GENOME} = standardize_genome($FLAGS{genome} || die "-genome");
    $r{input_row_number} = $input_row_number;
    $r{genea_acc} = "";
    $r{genea_acc_preferred} = "";
    $r{genea_symbol} = "";
    $r{genea_pos_adj} = "";
    $r{genea_coding_distance} = "";
    $r{genea_feature} = "";
    $r{geneb_acc} = "";
    $r{geneb_acc_preferred} = "";
    $r{geneb_symbol} = "";
    $r{geneb_pos_adj} = "";
    $r{geneb_coding_distance} = "";
    $r{geneb_feature} = "";
    $r{genea_full} = "";
    $r{geneb_full} = "";
    $r{extended_exon} = "";
    $r{$F_EXTENDED_CHUNK} = "";
    $r{extended_max} = "";
    $r{$F_NOTES} = join ",", @notes_general;
    $r{genea_contig_start} = "";
    $r{genea_contig_end} = "";
    $r{geneb_contig_start} = "";
    $r{geneb_contig_end} = "";
    $r{gene_pair_summary} = "";
    $r{gene_pair_summary_extended} = "";
    $r{$F_PROCESSING_TYPE} = "fusion_contig";
    $rpt->end_row(\%r);
    next;
  }
  my %saw;
  my @out_rows;

  my $no_pos_mode = $FLAGS{"no-pos"};

  my ($pos_a_raw, $pos_b_raw);
  if ($no_pos_mode) {
    $pos_a_raw = $pos_b_raw = VALUE_NA;
  } else {
    $pos_a_raw = $row->{$F_POS_A} || die;
    $pos_b_raw = $row->{$F_POS_B} || die;
  }

  foreach my $acc_a (sort keys %{$gene_a_transcript2rf}) {
    my $rf_a = $gene_a_transcript2rf->{$acc_a} || die;
    my $is_suppl_a = is_supplemental($rf_a);

    my $pos_a_adj;
    if ($no_pos_mode) {
      $pos_a_adj = VALUE_NA;
    } elsif ($is_suppl_a) {
      # not supported
      $pos_a_adj = $pos_a_raw;
    } else {
      $pos_a_adj = $tf->get_upstream_downstream(
						"-row" => $rf_a,
						"-direction" => "upstream",
						"-pos" => $pos_a_raw,
						"-pos-only" => 1
					       );
    }

    my @notes_a = @notes_general;
    my $contig_raw = $row->{$F_CONTIG} || dump_die($row, "no field $F_CONTIG");
    # put inside this loop because it may be expanded by large ITD detection,
    # and this process is specific to each isoform
    my $qc = $contig_raw;
    $qc =~ s/[ACGTN]//ig;
    dump_die($row, "ERROR: contig contains non-ACGTN characters $qc") if $qc;
    # sanity check

    my $cd_a = "";
    if ($pos_a_raw ne VALUE_NA and
	$pos_a_adj ne VALUE_NA and
	$pos_a_raw != $pos_a_adj) {
      $cd_a = abs($pos_a_raw - $pos_a_adj);
    }
    $row->{genea_coding_distance} = $cd_a;
    deep_intronic_check($row, $cd_a, "A", \$contig_raw, \@notes_a);

    #    dump_die($row, "hey now $pos_a_raw $pos_a_adj") if $pos_a_raw != $pos_a_adj;
    my $blast = get_blast();
    my ($genomic_a, $idx2exon_a);
    if ($is_suppl_a) {
      $genomic_a = get_supplemental_sequence($rf_a, $contig_raw, \@notes_a, "A");
    } else {
      ($genomic_a, $idx2exon_a) =
	get_genomic_sequence($acc_a, $rf_a, \%rf2genomic);
    }

    my ($i_start_acc, $i_end_acc, $i_start_contig, $i_end_contig);
    my $upstream_a;
    # entire upstream sequence of geneA
    my $upstream_a_exon;
    # upstream sequence of geneA for current exon only
    my $upstream_a_chunk;
    # upstream sequence limited to a maximum length

#    my $blat = get_blat();
#    $blat->verbose(1) if $VERBOSE;

    my @duplicated_blocks;
    my $is_intragenic_del;

    if ($is_intragenic) {
      #
      #  events within a gene require special handling for ITDs, etc.
      #
      printf STDERR "intragenic BLAST search for %s:\n", $acc_a if $VERBOSE;

      die "intragenic handling of supplemental sequences not yet supported" if $is_suppl_a;


      my @hsp = blast2hsps(
			    "-blast" => $blast,
			    "-query" => {
					 $acc_a => $genomic_a
					},
			    "-database" => { "contig_raw" => $contig_raw }
			   );

      if (my $new_contig = detect_large_itd(
					    "-hsp" => \@hsp,
					    "-contig" => $contig_raw,
					    "-gene" => $genomic_a,
					    "-notes" => \@notes_a
					   )) {
	push @notes_a, "expanded_contig_based_on_presumed_large_ITD";
	$contig_raw = $new_contig;

	@hsp = blast2hsps(
			     "-blast" => $blast,
			     "-query" => {
					  $acc_a => $genomic_a
					 },
			     "-database" => { "contig_raw" => $contig_raw }
			    );
      }

      #
      #  in some cases ITDs are represented within a single HSP by
      #  a gap in the query (gene) sequence vs. the hit (contig) sequence:
      #
      if (my $blocks = find_duplicate_blocks_via_query_gap(\@hsp, $contig_raw)) {
	push @duplicated_blocks, @{$blocks};
	push @notes_a, "ITD_via_gapped_HSP";
      }


      #
      #  in some cases ITDs can be identified by overlapping regions
      #  of the query sequence in multiple HSPs:
      #
      if (@hsp > 1) {
	my @hsp_q_overlap;
	# find HSPs where portions of the query overlap
	for (my $i = 0; $i < @hsp; $i++) {
	  # hacky, replace w/Algorithm::Combinatorics::combinations?
	  my $ok = 0;
	  for (my $j = 0; $j < @hsp; $j++) {
	    next if $i == $j;

#	    my $s1 = new Set::IntSpan join "-", $hsp[$i]->range("query");
#	    my $s2 = new Set::IntSpan join "-", $hsp[$j]->range("query");
	    my $s1 = new Set::IntSpan join "-", get_trimmed_range($hsp[$i], "query", 1);
	    my $s2 = new Set::IntSpan join "-", get_trimmed_range($hsp[$j], "query", 1);
	    # FIX ME: this whole block is probably hideously slow

#	    printf STDERR "idx %d %s\n", $i, $s1;
#	    printf STDERR "idx %d %s\n", $j, $s2;

 	    my $intersect = intersect $s1 $s2;

	    if ($intersect->size()) {
	      # keep HSPs whose query ranges overlap
	      $ok = 1;
	      last;
	    }
	  }
#	  printf STDERR "idx status = %d\n", $ok;
	  push @hsp_q_overlap, $hsp[$i] if $ok;
	}

	if (@hsp_q_overlap) {
#	  @hsp_q_overlap = sort {$a->start("query") <=>
#				   $b->start("query")} @hsp_q_overlap;
	  # sort blocks in query order
	  @hsp_q_overlap = sort {$a->start("hit") <=>
				   $b->start("hit")} @hsp_q_overlap;
	  # sort by HIT order so that duplicate blocks are (hopefully)
	  # ordered as well

	  printf STDERR "HSP q ranges: %s\n", join ",", map {join "-", $_->range("query")} @hsp_q_overlap;
	  printf STDERR "HSP h ranges: %s\n", join ",", map {join "-", $_->range("hit")} @hsp_q_overlap;

	  if (@hsp_q_overlap > 2) {
	    @hsp_q_overlap = prune_query_overlaps(\@hsp_q_overlap);
	  }

	  my $q1 = new Set::IntSpan join "-", $hsp_q_overlap[0]->range("query");
	  my $q2 = new Set::IntSpan join "-", $hsp_q_overlap[1]->range("query");
	  my $intersect = intersect $q1 $q2;
	  my @contig_blocks;

	  if ($intersect->size()) {
	    # a portion of the QUERY (full gene sequence) maps to
	    # multiple regions of the TARGET (contig), identify
	    # affected contig regions
	    my $i_start = $intersect->min;
	    my $i_end = $intersect->max;
	    my $i_length = $intersect->size;
	    printf STDERR "q intersect: %d-%d, %d\n", $i_start, $i_end, $i_length;

	    foreach my $hsp (@hsp_q_overlap) {
	      printf STDERR "raw hsp range: %d-%d\n", $hsp->range("query");
	    }

	    my @contig_ranges;

	    foreach my $hsp (@hsp_q_overlap) {
	      my ($q_start, $q_end) = $hsp->range("query");
	      my ($h_start, $h_end) = $hsp->range("hit");
	      push @contig_ranges, [ $h_start, $h_end ];

	      my $offset = $i_start - $q_start;

#	      my $bs = $h_start + $offset;
#	      my $be = $bs + ($i_length - 1);
	      # each of these is a duplicated block in the target/contig
	      # FLAWED:
	      # - start offset may not properly account for gaps in alignment!
	      # - need function to convert a base number(s) in the QUERY
	      #   to the corresponding base in the HIT
	      my ($bs, $be) = query_to_hit($hsp, $i_start, $i_end);

	      printf STDERR "duplicated block %d-%d\n", $bs, $be;
	      push @duplicated_blocks, [ $bs, $be ];
	    }

#	    my @qs = map {$_->start()} @hsp_q_overlap;
	    my @qs;
	    foreach my $h (@hsp_q_overlap) {
	      my $tr = get_trimmed_range($h, "query");
	      push @qs, $tr->[0];
	    }

	    my @ds = map {$_->[0]} @duplicated_blocks;

	    if ($qs[0] <= $qs[1] and $ds[0] < $ds[1]) {
	      # query blocks are ordered the same as the duplicate
	      # blocks.  This is expected, i.e. the RHS of the first
	      # query-ordered block matches the first duplicated
	      # region in the contig, and the LHS of the second
	      # query-ordered block matches the second duplicated
	      # region in the contig.
	    } else {
	      push @notes_a, sprintf "unexpected_block_order=q:%d/%d;dup:%d/%d",
		@qs, @ds;
	      # sometimes this is benign
	    }

	    die "duplicate blocks out of order" unless $duplicated_blocks[0]->[0] <= $duplicated_blocks[1]->[0];
	    # ordering sanity check

	    #
	    #  check for a large gap between the end of the first duplicated
	    #  block and the start of the second.  Sometimes the duplicated
	    #  regions are NOT adjacent!  A small amount of interstitial
	    #  sequence is fine, but beyond that we should just use the
	    #  BLAST alignment boundaries.
	    #
	    # example: /home/medmonso/work/steve/2020_06_13_cicero_contig_extension/set2_2020_09_25/28/
	    #

	    my $gap = ($duplicated_blocks[1]->[0] - $duplicated_blocks[0]->[1]) - 1;
	    printf STDERR "gap=%d\n", $gap;
	    push @notes_a, sprintf "duplicate_block_gap=%d", $gap if $gap > 0;
	    if (defined($MAX_DUPLICATE_BLOCK_DISTANCE) and
		$gap > $MAX_DUPLICATE_BLOCK_DISTANCE) {
	      push @notes_a, sprintf "raw_duplicate_blocks=%s", join ";", map {join "-", @{$_}} @duplicated_blocks;
	      # record originally-detected duplicate regions

	      my $new_end = $contig_ranges[0]->[1];
	      printf STDERR "revised dup block 1 end: %d-%d\n", @{$contig_ranges[0]};

	      $duplicated_blocks[0]->[1] = $new_end;
	      # TO DO: special handling for overlaps in this case?
	      # Or is standard overlap check below good enough?
	    }

	    resolve_block_overlap(\@duplicated_blocks, \@hsp_q_overlap);
	  }

#	  @duplicated_blocks = sort {$a->[0] <=> $b->[0]} @duplicated_blocks;
	  # ensure ordered by position in contig

	  #
	  #  look for overlaps directly in the contig sequence:
	  #
	  my $h1 = new Set::IntSpan join "-", $hsp[0]->range("hit");
	  my $h2 = new Set::IntSpan join "-", $hsp[1]->range("hit");
	  $intersect = intersect $h1 $h2;
	  if ($intersect->size()) {
	    # duplicated portion of the hit (contig)
	    printf STDERR "contig duplicated block %d-%d\n", $intersect->min(), $intersect->max();
	    push @contig_blocks, [ $intersect->min(), $intersect->max() ];
	  }

	  if (@duplicated_blocks and @contig_blocks) {
	    printf STDERR "NOTE: found duplicated blocks via 2 methods\n";
	    # ever a problem?
	  }

	  if (@duplicated_blocks) {
	    die "not exactly 2 dup blocks" unless @duplicated_blocks == 2;
	    push @notes_a, "ITD_via_multiple_HSP";
	    # TO DO: better way to to handle this?
	    # flag ENTIRE overrlap region?
	  } elsif (@contig_blocks) {
	    # 2 different portions of the query/gene sequence match
	    # 1 portion of the contig.  NOT an ITD.

#	    ($contig_gap_start, $contig_gap_end) = @{$contig_blocks[0]};
	    # should we even treat this kind of thing differently
	    # than normal processing??  artifact?
	    #
	    # NEED GOOD EXAMPLE OF THIS TYPE IF RELEVANT

#	    die "$contig_gap_start $contig_gap_end";
	  } else {
	    # no overlaps detected in either query or contig, OK
	  }
	} else {
	  #
	  #  multiple HSPs, but no regions of the query overlap.
	  #  Possible intragenic deletion.
	  #
	  $is_intragenic_del = 1;
	  push @notes_a, "intragenic_deletion";
	  # assume that's the case for now

	  my @sort_by_hit = sort {$a->start("hit") <=>
				    $b->start("hit")} @hsp;
	  # TO DO:
	  # update this to use the TRIMMED hit start

	  printf STDERR "multiple HSPs but no query overlap:\n";
	  my $prev;
	  my $query_ordered = 1;
	  foreach my $hsp (@sort_by_hit) {
	    printf STDERR "  q:%d-%d hit:%d-%d\n", $hsp->range("query"), $hsp->range("hit");
	    if ($prev and $prev->start("query") > $hsp->start("query")) {
	      $query_ordered = 0;
	    }
	    $prev = $hsp;

#	    push @duplicated_blocks, [ $hsp->range("hit") ];
	    push @duplicated_blocks, get_trimmed_range($hsp, "hit");

	  }

	  unless ($query_ordered) {
	    push @notes_a, sprintf 'unexpected_hsp_order=q:%s;hit:%s',
	      join("/", map {$_->start("query")} @sort_by_hit),
		join("/", map {$_->start("hit")} @sort_by_hit);
	  }
	  # not expected, problem???

	  resolve_block_overlap(\@duplicated_blocks, \@sort_by_hit);
	}
      }
    }

    if (@duplicated_blocks and not($is_intragenic_del)) {
      push @notes_a, sprintf "duplicate_blocks=%s", join ";", map {join "-", @{$_}} @duplicated_blocks;
    }

    #
    #  main processing, primarily for intergenic events.
    #  anchor gene A genomic sequence to contig:
    #
    my $contig_masked = $contig_raw;
    if ($MASK_CONTIG_BEFORE_ALIGN) {
      my $mask_len;
      if (@duplicated_blocks) {
	$mask_len = length($contig_raw) - $duplicated_blocks[0]->[1];
	# ITD/duplicate mode: assign first duplicate block to geneA
      } elsif ($MASK_RESERVE_GENE_B) {
	$mask_len = int(length($contig_raw) * $GENEA_ANCHOR_MASK_FRACTION_FOR_B);
      }
      substr($contig_masked, - $mask_len, $mask_len, "N" x $mask_len) if $mask_len;
    }

    if ($VERBOSE) {
      printf STDERR "   raw: %s\n", $contig_raw;
      printf STDERR "masked: %s\n", $contig_masked;
      # intra-gene events with duplicated regions can cause problems
      # with BLAT results, where a single BLAT HSP may contain multiple
      # sub-alignments that don't seem to be exposed in the API.
      # In this situation e.g. FLT3/FLT3 with a duplicated region,
      # the result matches to the entire contig even though this
      # incorporates different sub-matches.  Not sure if BLAST would
      # work better, or if the API would do the same thing.
      # In any case, we really only care about aligning the leftmost
      # side of the contig, so masking the geneB side will hopefully
      # work around this problem.
      printf STDERR "blast A (%s/%s) contig %s (%d) vs full-length %s\n", $gene_a, $acc_a, $contig_masked, length($contig_raw), $genomic_a;
    }

    my $a_start = "";
    my $a_end = "";
    my $a_start_raw = "";
    my $a_end_raw = "";

    if (my $hsps_a = run_blast(
			     "-blast" => $blast,
			     "-query" => {
					  (join "_", "A", $gene_a, $acc_a) => $genomic_a
					 },
			     "-database" => { "contig_masked_for_a_align" => $contig_masked },
			     "-notes" => \@notes_a,
			     "-end" => "A",
			      "-multi-ok" => 1,
			       "-is-suppl" => $is_suppl_a
			    )) {
#      ($a_start, $a_end) = $hsp_a->range("hit");
      ($a_start, $a_end) = hsps_range($hsps_a, "hit");

      my $range_a = sprintf '%d-%d', $a_start, $a_end;
      printf STDERR "geneA: %s\n", $range_a if $VERBOSE;
      push @notes_a, sprintf 'geneA_contig_match=%s', $range_a;
      # for geneA we are interested in the genomically left-most (5')
      # part of the match in the contig, which identifies the point
      # where we can extend the sequence further into geneA.
      # - ideally hit will cover first base of the contig sequence;
      #   how to deal if it doesn't??
      # - extension is not required to be base-perfect however, we just
      #   want some more sequence to match against

      ($i_start_acc, $i_end_acc, $i_start_contig, $i_end_contig) =
	get_hsps_indexes($hsps_a);
      # TO DO: record note/exception for geneA if $i_start_contig > 0

      push @notes_a, sprintf "geneA_contig_start_offset=%d", $i_start_contig
	unless $i_start_contig == 0;

      unless ($is_intragenic) {
	#
	#  for intergenic events only, also perform blat of RAW contig
	#  vs geneA, which will be used to find overlaps with geneB,
	#  e.g. HMBOX1_ALK
	#
	my $hsps_a_raw = run_blast(
				   "-blast" => $blast,
				   "-query" => {
						(join "_", "A", $gene_a, $acc_a) => $genomic_a
					       },
				   "-database" => { "contig_raw" => $contig_raw },
				   "-notes" => \@notes_a,
				   "-end" => "A",
				   "-multi-ok" => 1,
				  );
	# this should work since if we got this far we know
	# masked version ran OK
	($a_start_raw, $a_end_raw) = hsps_range($hsps_a_raw, "hit");
      }

      #
      #  find extended sequence for geneA:
      #
#      my $strand_a = $hsp_a->strand();
      my $strand_a = hsps_strand($hsps_a);

      if ($strand_a == 1) {
	if ($VERBOSE) {
	  printf STDERR "%s\n", substr($contig_raw, $i_start_contig, ($i_end_contig - $i_start_contig) + 1);
	  printf STDERR "%s\n", substr($genomic_a, $i_start_acc, ($i_end_acc - $i_start_acc) + 1);
	}

	# to extend sequence for a + hit, look left
	$upstream_a = substr($genomic_a, 0, $i_start_acc);
	# all upstream sequence for geneA

	$upstream_a_chunk = substr($upstream_a, - $EXTENDED_CHUNK_LENGTH);
	# ...limited to a maximum length

	if ($is_suppl_a) {
	  # not supported
	  $upstream_a_exon = $upstream_a;
	} else {
	  $upstream_a_exon = find_exonic($genomic_a, $idx2exon_a, $i_start_acc, -1);
	  # ...limited to exon boundary
	}

	if ($VERBOSE) {
	  printf STDERR "             geneA: %s\n", $genomic_a;
	  printf STDERR "          upstream: %s\n", $upstream_a;
	  printf STDERR "    upstream_chunk: %s\n", $upstream_a_chunk;
	  printf STDERR "upstream_exon_only: %s\n", $upstream_a_exon;
	  #	die join ",", $i_start_contig, $i_end_contig, $i_start_acc, $i_end_acc;
	}
      } else {
	die "geneA genomic match not on + strand";
	# this shouldn't happen as - mapped genes have their sequence
	# reverse-complemented before blat
      }
    } else {
      push @notes_a, "geneA_fail_to_align_contig_to_gene";

      my $bl = get_blast();
      $bl->strand("minus");

      if (my $hsps_a = run_blast(
				 "-blast" => $bl,
				 "-query" => {
					      (join "_", "A", $gene_a, $acc_a) => $genomic_a
					     },
				 "-database" => { "contig_masked_for_a_align" => $contig_masked },
				 "-notes" => [],
				 "-end" => "A",
				 "-multi-ok" => 1,
				 "-is-suppl" => $is_suppl_a
				)) {
	push @notes_a, sprintf "%s=A", TAG_ERROR_BLAST_HIT_MINUS;
      }


    }

    if (scalar keys %{$gene_b_transcript2rf}) {
      foreach my $acc_b (sort keys %{$gene_b_transcript2rf}) {
	#
	#  iterate through geneB accessions and anchor
	#

	next if ($is_intragenic and
		 $LOCK_ACCESSION_FOR_INTRAGENIC and
		 $acc_a ne $acc_b);

	next if $saw{$acc_a}{$acc_b};
	$saw{$acc_a}{$acc_b} = 1;

	printf STDERR "processing pair %s %s\n", $acc_a, $acc_b;

	my $rf_b = $gene_b_transcript2rf->{$acc_b} || die;
	my $is_suppl_b = is_supplemental($rf_b);

	my $pos_b_adj;
	if ($no_pos_mode) {
	  $pos_b_adj = VALUE_NA;
	} elsif ($is_suppl_b) {
	  # not supported
	  $pos_b_adj = $pos_b_raw;
	} else {
	  $pos_b_adj = $tf->get_upstream_downstream(
						    "-row" => $rf_b,
						    "-direction" => "downstream",
						    "-pos" => $pos_b_raw,
						    "-pos-only" => 1
						   );
	}

	my @notes = @notes_a;

	my $cd_b = "";
	if ($pos_b_raw ne VALUE_NA and $pos_b_adj ne VALUE_NA and
	    $pos_b_raw != $pos_b_adj) {
	  $cd_b = abs($pos_b_raw - $pos_b_adj);
	}
	$row->{geneb_coding_distance} = $cd_b;
	deep_intronic_check($row, $cd_b, "B", \$contig_raw, \@notes);

	  #	dump_die($row, "hey now $acc_b $pos_b_raw $pos_b_adj") if $pos_b_raw != $pos_b_adj and $acc_b =~ /NM_/;

	my ($genomic_b, $idx2exon_b);
	if ($is_suppl_b) {
	  $genomic_b = get_supplemental_sequence($rf_b, $contig_raw, \@notes, "B");
	} else {
	  ($genomic_b, $idx2exon_b) =
	    get_genomic_sequence($acc_b, $rf_b, \%rf2genomic);
	}

	if ($VERBOSE) {
	  printf STDERR "full-length A %s: %s\n", $acc_a, $genomic_a;
	  printf STDERR "full-length B %s: %s\n", $acc_b, $genomic_b;
	}

	my $downstream_b;
	my $downstream_b_chunk;
	my $downstream_b_exon;

	my ($b_start_raw, $b_end_raw, $gene_overlap);

	#
	#  anchor geneB:
	#
	if (not($is_intragenic)) {
	  if (my $hsps_b_raw = run_blast(
				 "-blast" => $blast,
				 "-query" => {
					      (join "_", "B", $gene_b, $acc_b) => $genomic_b
					     },
				 "-database" => { "contig_raw" => $contig_raw },
				 "-notes" => \@notes,
				 "-end" => "B",
				 "-is-suppl" => $is_suppl_b,
				"-multi-ok" => 1
				      )) {
#	    ($b_start_raw, $b_end_raw) = $hsp_b_raw->range("hit");
	    ($b_start_raw, $b_end_raw) = hsps_range($hsps_b_raw, "hit");

	    printf STDERR "raw hits: A/%s=", $acc_a;
	    if ($a_start_raw) {
	      printf STDERR "%d-%d", $a_start_raw, $a_end_raw;
	    } else {
	      print STDERR "null";
	    }
	    printf STDERR " B/%s=", $acc_b;
	    if ($b_start_raw) {
	      printf STDERR "%d-%d", $b_start_raw, $b_end_raw;
	    } else {
	      print STDERR "null";
	    }
	    print STDERR "\n";

	    if ($a_start_raw and $b_start_raw and
		$a_end_raw > $b_start_raw) {
	      my $overlap = ($a_end_raw - $b_start_raw) + 1;
	      if ($overlap >= $MIN_INTERGENIC_OVERLAP_TO_MASK) {
		# a possibly better method would be to decide which
		# of the gene A and B alignments are better, however
		# current design anchors A first.  Since extension is the
		# main goal rather than perfect breakpoints this is probably
		# good enough.
		$gene_overlap = $overlap;
		printf STDERR "OVERLAP!! %s\n", lc(substr($contig_raw, $b_start_raw - 1, $overlap));
	      }
	    }

	  }
	}

	$contig_masked = $contig_raw;
	if ($a_start and $a_end and $MASK_CONTIG_BEFORE_ALIGN) {
	  # mask out geneA match before aligning to geneB.
	  # TO DO: mask out entire left side of contig up to end of blat match?

#	  my $len = $a_end - $a_start + 1;
#	  substr($contig_masked, $a_start - 1, $len, "N" x $len);
	  #
	  # this masks only the matching portion.  Trouble with that is
	  # the upstream portion may still spuriously match esp. if
	  # both genes are in the same gene family, etc.
	  # /home/medmonso/work/steve/2020_06_13_cicero_contig_extension/2021_02_03_clingen_all/coding_only/crash7/test.tab

	  my $len = $a_end;
	  substr($contig_masked, 0, $len, "N" x $len);
	  # mask from start to end of A match (safer)

	  my $mask_frac = $len / length($contig_raw);
	  if ($mask_frac > 0.90) {
	    push @notes, sprintf "geneA_match_contig_coverage_percent=%d", $mask_frac * 100;
	  }
	}

	printf STDERR "blast full-length %s vs B contig %s\n", $genomic_b, $contig_masked if $VERBOSE;
	my $b_start = "";
	my $b_end = "";

	printf STDERR "blast B (%s/%s) contig %s (%d) vs full-length raw %s\n", $gene_b, $acc_b, $contig_masked, length($contig_masked), $genomic_b if $VERBOSE;

	if (my $hsps_b = run_blast(
				  "-blast" => $blast,
				  "-query" => {
					       (join "_", "B", $gene_b, $acc_b) => $genomic_b
					      },
				  "-database" => { "contig_masked_for_b_align" => $contig_masked },
				  "-notes" => \@notes,
				  "-end" => "B",
				   "-is-suppl" => $is_suppl_b,
				  "-multi-ok" => 1,
				 )) {
	  # TO DO: any point in masking contig for geneA hit?

#	  ($b_start, $b_end) = $hsp_b->range("hit");
	  ($b_start, $b_end) = hsps_range($hsps_b, "hit");
	  my $range_b = sprintf '%d-%d', $b_start, $b_end;
	  printf STDERR "geneB: %s\n", $range_b if $VERBOSE;
	  push @notes, sprintf 'geneB_contig_match=%s', $range_b;

	  push @notes, $TAG_ERROR_GENES_REVERSED
	    if ($a_start and $a_start > $b_start);
	  # SERIOUS ERROR, and fatal situation for fuzzion conversion.
	  # use "fusion_gene" field instead of A/B?  Are genes
	  # better ordered in this field?
	  # => there are isoform-dependent cases, so maybe not critical

	  ($i_start_acc, $i_end_acc, $i_start_contig, $i_end_contig) =
	    get_hsps_indexes($hsps_b);
	  # TO DO: record note/exception for geneB if $i_start_contig < len

	  my $r_offset = length($contig_raw) - ($i_end_contig + 1);
	  push @notes, sprintf "geneB_contig_end_offset=%d", $r_offset
	    unless $r_offset == 0;
	  my $strand_b = hsps_strand($hsps_b);
	  if ($strand_b == 1) {
	    # for + hits to geneB, we are interested in the genomically
	    # right-most (5') part of the match in the contig, which
	    # identifies the point where we can extend the sequence
	    # further into geneB.
	    if ($VERBOSE) {
	      printf STDERR "rawctg %s\n", $contig_raw;
	      printf STDERR "contig %s\n", substr($contig_raw, $i_start_contig, ($i_end_contig - $i_start_contig) + 1);
	      printf STDERR " geneB %s\n", substr($genomic_b, $i_start_acc, ($i_end_acc - $i_start_acc) + 1);
	    }

	    $downstream_b = substr($genomic_b, $i_end_acc + 1);
	    # all downstream sequence for B

	    $downstream_b_chunk = substr($downstream_b, 0, $EXTENDED_CHUNK_LENGTH);

	    if ($is_suppl_b) {
	      # not supported
	      $downstream_b_exon = $downstream_b;
	    } else {
	      $downstream_b_exon = find_exonic($genomic_b, $idx2exon_b, $i_end_acc, 1);
	    }

	    if ($VERBOSE) {
	      printf STDERR "           geneB: %s\n", $genomic_b;
	      printf STDERR "      downstream: %s\n", $downstream_b;
	      printf STDERR "downstream_chunk: %s\n", $downstream_b_chunk;
	      printf STDERR "        downExon: %s\n", $downstream_b_exon;
	    }
	  } else {
	    die "geneB genomic match not on + strand";
	  }
	} else {
	  push @notes, "geneB_fail_to_align_geneA_masked_contig_to_gene";

	  my $bl = get_blast();
	  $bl->strand("minus");
	  if (my $hsps_b = run_blast(
				  "-blast" => $bl,
				  "-query" => {
					       (join "_", "B", $gene_b, $acc_b) => $genomic_b
					      },
				  "-database" => { "contig_masked_for_b_align" => $contig_masked },
				  "-notes" => [],
				  "-end" => "B",
				   "-is-suppl" => $is_suppl_b,
				  "-multi-ok" => 1,
				 )) {
	    push @notes, sprintf "%s=B", TAG_ERROR_BLAST_HIT_MINUS;
	    # contig in wrong orientation?  Antisense?
	    # CREST-sourced data?
	  }
	}

	my $extended_max = "";
	if ($upstream_a and $downstream_b) {
	  # success
	  $extended_max = lc($upstream_a) . uc($contig_raw) . lc($downstream_b);
	  #	printf STDERR "full_length: %s\n", $extended_max;
	} else {
	  push @notes, "missing_A_upstream" unless $upstream_a;
	  push @notes, "missing_B_downstream" unless $downstream_b;
	}

	my $extended_chunk = "";
	if ($upstream_a_chunk and $downstream_b_chunk) {
	  # success
	  $extended_chunk = lc($upstream_a_chunk) . uc($contig_raw) . lc($downstream_b_chunk);
	  #	printf STDERR "full_length: %s\n", $extended_max;
	} else {
	  push @notes, "missing_A_upstream_chunk" unless $upstream_a_chunk;
	  push @notes, "missing_B_downstream_chunk" unless $downstream_b_chunk;
	}

	my $extended_exon = "";
	if ($upstream_a_exon and $downstream_b_exon) {
	  $extended_exon = lc($upstream_a_exon) . uc($contig_raw) . lc($downstream_b_exon);
	} else {
	  push @notes, "missing_A_exon_upstream" unless $upstream_a_exon;
	  push @notes, "missing_B_exon_downstream" unless $downstream_b_exon;
	}

	#
	#  report results:
	#
	my %r = %{$row};
	$r{genea_acc} = $acc_a;
	$r{genea_acc_preferred} = get_sjpi($SJPI, $acc_a);
	$r{genea_symbol} = $transcript2gene{$acc_a} || die;
	$r{genea_pos_adj} = $pos_a_adj;

	if ($no_pos_mode) {
	  $r{genea_feature} = $r{geneb_feature} = VALUE_NA;
	} else {
	  $r{genea_feature} = $is_suppl_a ? VALUE_NA : get_feature_tag($rf, $rf_a, $pos_a_adj);

	  $r{geneb_feature} = $is_suppl_b ? VALUE_NA : get_feature_tag($rf, $rf_b, $pos_b_adj);
	}

	$r{geneb_acc} = $acc_b;
	$r{geneb_acc_preferred} = get_sjpi($SJPI, $acc_b);
	$r{geneb_symbol} = $transcript2gene{$acc_b} || die;
	$r{geneb_pos_adj} = $pos_b_adj;

	$r{genea_full} = $genomic_a;
	$r{geneb_full} = $genomic_b;

	if ($DMZ_OVERLAP and $gene_overlap) {
#	  print STDERR "hey now $gene_overlap $a_start $a_end $b_start $b_end\n$a_start_raw $a_end_raw $b_start_raw $b_end_raw\n";

	  $a_end = $b_start_raw - 1 if $b_start_raw <= $a_end;
	  # -1 to set boundary before first base of B match

	  $b_start = $a_end_raw + 1 if $a_end_raw >= $b_start;
	  # +1 to set boundary after first base of A match

#	  die "adj A end=$a_end  B start=$b_start\n";
	}

	$r{genea_contig_start} = $a_start;

	if (@duplicated_blocks) {
	  # e.g. ITD
	  $r{genea_contig_end} = $duplicated_blocks[0]->[1];
	  # end of first duplicated block
	  $r{geneb_contig_start} = $duplicated_blocks[$#duplicated_blocks]->[0];
	  # start of last duplicated block
	} else {
	  $r{genea_contig_end} = $a_end;
	  $r{geneb_contig_start} = $b_start;
	}

	#
	#  find max distance between each endpoint and the center,
	#  can be helpful for finding buggy calls esp. for ITDs
	#
	my $center = int(length($contig_raw) / 2);
	my @dist;
	foreach my $f (qw(genea_contig_end geneb_contig_start)) {
	  my $v = $r{$f} || next;
	  push @dist, abs($center - $v);
	}
#	die join ",", $center, $r{genea_contig_end}, $r{geneb_contig_start};
	if (@dist) {
	  my $max_dist = max(@dist);
	  push @notes, sprintf 'max_distance_from_center=%.2f', $max_dist * 100 / $center;
	}

	if ($r{genea_contig_end} and $r{geneb_contig_start} and
	    $r{genea_contig_end} > $r{geneb_contig_start}) {
	  my $not_fatal = 0;
	  dump_die(\%r, '"that\'s unpossible": genea_contig_end > geneb_contig_start', $not_fatal);
	  printf STDERR "WARNING: MAKE THIS A FATAL ERROR\n" if $not_fatal;
	}

	$r{geneb_contig_end} = $b_end;

	$r{extended_exon} = $extended_exon;
	$r{$F_EXTENDED_CHUNK} = $extended_chunk;
	$r{extended_max} = $extended_max;
	$r{$F_NOTES} = join ",", remove_list_duplicates(\@notes);
	$r{input_row_number} = $input_row_number;
#	$rpt->end_row(\%r);
	push @out_rows, \%r;
	# will be multiple output rows per input row depending on
	# combinations

      }				# $acc_b loop
    } else {
      # no transcripts available for geneB, e.g. IGH
      my @notes = @notes_a;
      push @notes, "no_gene_B_transcripts";

      my %r = %{$row};
      $r{input_row_number} = $input_row_number;
      $r{genea_acc} = $acc_a;
      $r{genea_symbol} = $transcript2gene{$acc_a} || die;
      $r{genea_feature} = $is_suppl_a ? VALUE_NA : get_feature_tag($rf, $rf_a, $pos_a_adj);
      $r{genea_pos_adj} = $pos_a_adj;
      $r{genea_coding_distance} = "";
      $r{geneb_acc} = "";
      $r{geneb_acc_preferred} = "";
      $r{geneb_symbol} = "";
      $r{geneb_feature} = "";
      $r{geneb_pos_adj} = $pos_b_raw;
      $r{geneb_coding_distance} = "";
      $r{genea_full} = $genomic_a;
      $r{geneb_full} = "";
      $r{extended_exon} = "";
      $r{$F_EXTENDED_CHUNK} = "";
      $r{extended_max} = "";
      $r{$F_NOTES} = join ",", @notes;
      $r{genea_contig_start} = "";
      $r{genea_contig_end} = "";
      $r{geneb_contig_start} = "";
      $r{geneb_contig_end} = "";
#      $rpt->end_row(\%r);
      push @out_rows, \%r;
    }
  }				# $acc_a loop

  #
  #  QC: does at least one output row have successful mapping?
  #      If not, could indicate a problem/bug
  #
  my $all_bad = 1;
  my %fatal_reasons;

  my %fatal_search = (
		      TAG_ERROR_NO_SEQUENCES_FOR_GENE() => "no_gene_sequences",
		      TAG_ERROR_BLAST_HIT_MINUS() => "possible_antisense_contig",
		      TAG_DEEP_INTRONIC() => "deep_intronic",
		     );

  foreach my $r (@out_rows) {
    if ($r->{genea_contig_start} and $r->{geneb_contig_start}) {
      # OK result for at least one record
      $all_bad = 0;
    } else {
      my $notes = $r->{$F_NOTES} || die;

#      printf STDERR "search notes %s\n", $notes;

      my $fatal_reason;
      foreach my $tag (keys %fatal_search) {
#	printf STDERR "  check %s %d\n", $tag, index($notes, $tag);
	if (index($notes, $tag) >= 0) {
	  $fatal_reason = $fatal_search{$tag};
	  last;
	}
      }

      $fatal_reason = "intergenic" if not($fatal_reason) and $r->{genea_feature} eq "intergenic" or $r->{geneb_feature} eq "intergenic";

      $fatal_reason = "other" unless $fatal_reason;
      $fatal_reasons{$fatal_reason}++;
    }
  }
  if ($all_bad) {
    foreach my $r (@out_rows) {

      my @sorted = sort {$fatal_reasons{$b} <=> $fatal_reasons{$a}} keys %fatal_reasons;
      my $details;
      if (@sorted == 1) {
	$details = join ";", map {$_ . ":all"} @sorted;
      } else {
	$details = join ";", map {$_ . ":" . $fatal_reasons{$_}} @sorted;
      }

      $r->{$F_NOTES} .= sprintf ",%s=%s", TAG_ERROR_NO_INPUT_RESULTS, $details;
    }
  }

  #
  #  TO DO: possible to optimize results by looking for results
  #  where isoforms match contig w/fewer HSPs/gaps, etc.?
  #
  #  - if some pairs match with just 1 HSP, are these superior
  #    to pairs matching with 2+ HSPs?
  #  - if some pairs match both sides, can we just drop rows where
  #    no matches?
  #
#  printf STDERR "write %d rows\n", scalar @out_rows;
  my %hsp_counts;
  foreach my $r (@out_rows) {
    my $notes = $r->{$F_NOTES};
#    print STDERR "$notes\n";
    my %end2count;
    while ($notes =~ /raw_hsp_count_(\w)=(\d+)/g) {
      my ($end, $count) = ($1, $2);
      push @{$hsp_counts{$end}{$count}}, $r;
    }
  }

  foreach my $end (qw(A B)) {
    if (scalar keys %{$hsp_counts{$end}} > 1) {
      printf STDERR "end %s has counts %s\n", $end, join ",", sort {$a <=> $b} keys %{$hsp_counts{$end}};
#      dump_die($row, "NOTE: some isoforms appear to match better than others, prune results?", 1);
      # TO DO: option to prune here?
    }
  }

  #
  #  write final output rows for this input record:
  #
  foreach my $r (@out_rows) {
    $r->{$F_GENOME} = standardize_genome($FLAGS{genome} || die "-genome");
    $r->{genea_chr} = $r->{$F_CHR_A};
    $r->{genea_pos} = $r->{$F_POS_A};
    $r->{geneb_chr} = $r->{$F_CHR_B};
    $r->{geneb_pos} = $r->{$F_POS_B};
    # hack
    $r->{$F_PROCESSING_TYPE} = "fusion_contig";

    $r->{genea_acc_preferred} = get_sjpi($SJPI, $r->{genea_acc});
    $r->{geneb_acc_preferred} = get_sjpi($SJPI, $r->{geneb_acc});
    # TO DO: sharable sub call for this

    generate_pair_summary($r);
    $rpt->end_row($r);
  }

}  # input row
$rpt->finish();

sub get_transcripts {
  my ($rf, $gsm, $genes_raw, $notes_ref, $transcript2gene, %options) = @_;
  confess "transcript2gene not defined" unless $transcript2gene;

  my $ssl = $options{"-ssl"};

  my $chr_mode;
  my $chr;
  if ($chr = $options{"-chr"}) {
    # ensure mappings match a specific chrom
    $chr = clean_chr($chr);
    $chr_mode = 1;
  } elsif ($options{"-no-chr"}) {
    # no requirement, e.g. for full-gene-length pairs
    $chr_mode = 0;
  } else {
    die "specify -chr CHROM or -no-chr";
  }

  my $results;
  $genes_raw = unquote($genes_raw);
  my @genes = split /,/, $genes_raw;
  if ($CLEAN_SHARP_GENE_SYMBOLS) {
    # remove SJ refFlat "sharp" gene symbol modifications
    foreach (@genes) {
      my $orig = $_;
      s/_loc[A-Z]$//;
      printf STDERR "changed %s => %s\n", $orig, $_ if $orig ne $_;
    }
  }

  my %results;

  foreach my $gene (@genes) {

    my $found_gene;

    if (my $gene_rf = $gsm->find($gene)) {
      # symbol may be different in refFlat

      my $has_coding;
      my $has_noncoding;

      my $rows = $rf->find_by_gene($gene_rf) || die;
      foreach my $row (@{$rows}) {
	my $acc = $row->{name};
	my $usable = 1;
	if ($CODING_ONLY and not($acc =~ /^NM_/)) {
	  # skip non-coding transcript, e.g. NR_
	  push @{$notes_ref}, sprintf "skip_noncoding=%s", $acc;
	  $usable = 0;
	}

	if ($acc =~ /^NM_/) {
	  $has_coding++;
	} else {
	  $has_noncoding++;
	}

	my $chr_rf = clean_chr($row->{chrom} || die);
	if ($chr_mode) {
	  $usable = 0 unless $chr eq $chr_rf;
	  # e.g. for CRLF2, NM_022148 has entries on chrX and chrY
	}

	if ($usable) {
	  printf STDERR "transcript for %s => %s\n", $gene, $acc;
	  $results{$acc} = $row;
	  # only keep one row per transcript: multiple mapping locations
	  # aren't relevant for our purposes
	  $transcript2gene->{$acc} = $gene;
	  $found_gene = 1;
	}
      }
#      die "gene $gene has only non-coding" if $has_noncoding and not($has_coding);
    }

    if (not($found_gene) and $ENABLE_SUPPLEMENTAL_REFSEQ) {
      # no gene model found, look in alternate FASTA file
      die "-ssl" unless $ssl;
      if ($ssl->has_gene($gene)) {
	# contains an entry for the desired gene
	my $seqs = $ssl->get_sequences_for_gene($gene);
	foreach my $seq_ref (@{$seqs}) {
	  my $acc = $seq_ref->{accession} || die;
	  $transcript2gene->{$acc} = $gene;
	  $results{$acc} = $seq_ref;
	}
	$found_gene = 1;
      }
    }

    unless ($found_gene) {
      push @{$notes_ref}, sprintf "%s=%s", TAG_ERROR_NO_SEQUENCES_FOR_GENE, $gene;
    }
  }
  return \%results;
}

sub get_genomic_sequence {
  my ($acc, $rf, $rf2genomic) = @_;

  unless ($rf2genomic{$acc}) {
    printf STDERR "building transcript sequence for %s, strand=%s:\n", $acc, $rf->{strand};
    my $seq = "";
    my $chrom = $rf->{chrom} || die;
    my @intervals = map { [ @{$_}{qw(start end)} ] } @{$rf->{exons} || die};
    my $exon_entry = 1;
    # technically not exon # as needs to be reversed for - strand entries
    my @idx2exon;

    for (my $i = 0; $i < @intervals; $i++) {
      if ($i > 0) {
	die "unexpected ordering" unless
	  $intervals[$i]->[0] > $intervals[$i - 1]->[0];
      }
      my ($start, $end) = @{$intervals[$i]};
      die "end < start" if $end < $start;
      my $chunk = $FAI->get_chunk(
				  "-id" => $chrom,
				  "-start" => $start,
				  "-end" => $end
				 );
      printf STDERR "%s:%d-%d (%d): %s  rc=%s\n", $chrom, $start, $end, ($end - $start + 1), $chunk, reverse_complement($chunk) if $VERBOSE;
      my $start_i = length($seq);
      $seq .= $chunk;
      my $end_i = length($seq);

      for (my $j=$start_i; $j < $end_i; $j++) {
	$idx2exon[$j] = $exon_entry;
      }
      $exon_entry++;
    }

    if ($VERBOSE > 1) {
      print STDERR "c. debug $acc\n";
      my $len = length $seq;
      for (my $i = 0; $i < $len; $i++) {
	printf STDERR "  %d: %s %s\n", $i, substr($seq, $i, 1), $idx2exon[$i];
      }
    }

    my $strand = $rf->{strand} || die;
    if ($strand eq "+") {
      # no need to do anything
    } elsif ($strand eq "-") {
      $seq = reverse_complement($seq);
      @idx2exon = reverse @idx2exon;
    } else {
      die "unhandled strand $strand";
    }

    $rf2genomic{$acc} = [ $seq, \@idx2exon ];
  }

  confess "returns array" unless wantarray;
  return @{$rf2genomic{$acc}};
}

sub get_blat {
  my $blat = new BLATer();
  $blat->null_mode(1) if $FLAGS{"no-blat"};
  $blat->minScore($BLAT_MIN_SCORE);
  if ($FLAGS{"save-tempfiles"}) {
    printf STDERR "saving tempfiles\n";
    $blat->tfw->auto_unlink(0);
  }
  return $blat;
}

sub get_blast {
  my $blast = new BLASTer();
  $blast->strand("plus");
  # contig and transcript sequence have been normalized to +
  $blast->perc_identity($BLAST_MIN_PERCENT_IDENTITY);
  $blast->word_size($BLAST_MIN_WORD_SIZE);

  $blast->dust("no");
  # in some cases a low-complexity region in the contig is enough to
  # prevent BLAST alignment from working, disabling dust filter works
  # around this problem, e.g.:
  #
  # /home/medmonso/work/steve/2020_06_13_cicero_contig_extension/2021_01_04_all/blast/fail_coding
  # blastn -query query_A_KMT2A_NM_005933_d9fab58c91185f60a828a81b73b843e9.fa -subject db_contig_raw_8e2fe069fb62b4c598fb66bc83e85679.fa

#  $blast->ungapped(1);
  # don't see how to extract gap positions from API,
  # many BLAST output formats don't seem to report these blocks
  # - ungapped alignments can cause problems with duplicate region detection

  $blast->blast_2_sequences(1);
  # just one sequence vs. another, don't need to index database

  $blast->output_format("xml");
  # provides visibility for gapped alignments

  $blast->gapopen($BLAST_GAP_OPEN) if $BLAST_GAP_OPEN;
  $blast->gapextend($BLAST_GAP_EXTEND) if $BLAST_GAP_EXTEND;

  if ($FLAGS{"save-tempfiles"}) {
    printf STDERR "saving tempfiles\n";
    $blast->tfw->auto_unlink(0);
  }
  return $blast;
}


sub run_blat {
  my %options = @_;
  my $blat = $options{"-blat"} || die "-blat";
  my $notes = $options{"-notes"} || die "-notes";
  my $end = $options{"-end"} || die "-end";

  my $parser = $blat->blat(%options);
  my $result = $parser->next_result;
  # one object per query sequence (only one query seq)
  my @hsp;
  my @hits;
  if ($result) {
    while (my $hit = $result->next_hit()) {
      # hits from this query to a database sequence
      # (Bio::Search::Hit::HitI)
      push @hits, $hit;

#      die $hit->ambiguous_aln;
      # no help for FLT3/FLT3

      while (my $hsp = $hit->next_hsp) {
	# Bio::Search::HSP::PSLHSP
	printf STDERR "    score:%s strand:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s query_span:%d ref_span:%d total_span=%d query_string=%s hit_string=%s\n",
	  $hsp->score,
	    $hsp->strand,
	      $hsp->range("query"),
		$hsp->range("hit"),
		  $hsp->num_identical(),
		    $hsp->frac_identical("query"),
		      $hsp->length("query"),
			$hsp->length("hit"),
			  $hsp->length("total"),
			    $hsp->query_string(),
			      $hsp->hit_string() if $VERBOSE;

	#	    push @hsp, $hsp if hsp_filter($hsp, \%hsp_problems);
	push @hsp, $hsp;
	# TO DO: filter??

#      die $hsp->seq_inds("hit", "identical");
	# undef for this implementation??

#	die $hsp->gaps("hit");

	if ($VERBOSE) {
	  printf STDERR "gap blocks:\n";
	  foreach my $k (qw(query hit)) {
	    printf STDERR "  %s\n", $k;
	    foreach my $gb (@{$hsp->gap_blocks($k)}) {
	      printf STDERR "    %s\n", join ",", @{$gb};
	    }
	  }

	  printf STDERR  "    qs:%d qe:%d hs:%d he:%d\n",
	    $hsp->start("query"),
	      $hsp->end("query"),
		$hsp->start("hit"),
		  $hsp->end("hit");

	  my $aln = $hsp->get_aln();
	  printf STDERR "alignment sequences: %d\n", $aln->num_sequences;
	  foreach my $seq ($aln->each_seq) {
	    # Bio::LocatableSeq
	    printf STDERR "  %s: %d-%d\n", $seq->id, $seq->start, $seq->end;
	  }
	}

#	die scalar $hsp->seq_inds('query', 'identical', 'collapse' => 1);
	# implemented??

      }
    }
  }

  my $retval;
  if (@hsp) {
    my @passed;
    my @neg;
    foreach my $hsp (@hsp) {
      if ($hsp->strand() == -1) {
	push @neg, $hsp;
      } else {
	push @passed, $hsp;
      }
    }
    push @{$notes}, sprintf "%s: removed %d minus strand hit(s)", $end, scalar @neg if @neg;
    # e.g. FGFR3/TACC3 contig starts with a small ?inversion? mapped to -
    # which is a small subset of the larger hit (interestingly though
    # this is the only match to the first ~41 nt of FGFR3 end)
    @hsp = @passed;

    #
    #  TO DO: require some minimum score?
    #
    if (@hsp > 1) {
      push @{$notes}, sprintf 'hsp_count=%d', scalar @hsp;

      if ($VERBOSE) {
	foreach my $hsp (@hsp) {
	  printf STDERR "  multi_hsp:  score:%s strand:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s query_span:%d ref_span:%d total_span=%d query_string=%s hit_string=%s gap_blocks_q=%s gap_blocks_h=%s\n",
	    $hsp->score,
	      $hsp->strand,
		$hsp->range("query"),
		  $hsp->range("hit"),
		    $hsp->num_identical(),
		      $hsp->frac_identical("query"),
			$hsp->length("query"),
			  $hsp->length("hit"),
			    $hsp->length("total"),
			      $hsp->query_string(),
				$hsp->hit_string(),
				  join(",", map {@{$_}} @{$hsp->gap_blocks("query")}),
				    join(",", @{$hsp->gap_blocks("hit")->[0]});
	}
      }

      my @scores = map {$_->score} @hsp;
      foreach (my $i = 1; $i < @scores; $i++) {
	# verify blat scores reported in descending order
	printf STDERR "ERROR: blat scores not in descending order!\n"
	  if $scores[$i] > $scores[$i - 1];
	# FIX ME:
	# this could be indicative of a broader problem e.g. for ITDs,
	# i.e. multiple HSPs rather than a single HSP containing gap blocks
      }
      @hsp = $hsp[0];
      # use best hit only, e.g.
      #
      # score:4.43 strand:1 q:2779-3030 subj:1-252 num_identical:252 frac_identical_query:1 query_span:252 ref_span:252 total_span=504
      # score:0.60 strand:1 q:1034-1067 subj:523-556 num_identical:34 frac_identical_query:1 query_span:34 ref_span:34 total_span=68
    }

    if (@hsp == 1) {
      $retval = $hsp[0];
    } else {
      push @{$notes}, sprintf "ERROR_multiple_blat_hits_%s=%d", $end, scalar @hsp;
    }
  }
  return $retval;
}

sub run_blast {
  my %options = @_;
  my $blast = $options{"-blast"} || die "-blat";
  my $notes = $options{"-notes"} || die "-notes";
  my $end = $options{"-end"} || die "-end";
  my $multi_ok = $options{"-multi-ok"};
  my $is_suppl = $options{"-is-suppl"};

  printf STDERR "start BLAST for %s\n", $end if $VERBOSE;

  my $parser = $blast->blast(%options);
  my $result = $parser->next_result;
  # one object per query sequence (only one query seq)
  my @hsp;
  my @hits;
  if ($result) {
    while (my $hit = $result->next_hit()) {
      # hits from this query to a database sequence
      # (Bio::Search::Hit::HitI)
      push @hits, $hit;

#      die $hit->ambiguous_aln;
      # no help for FLT3/FLT3

      while (my $hsp = $hit->next_hsp) {
	# Bio::Search::HSP::GenericHSP
	printf STDERR "    score:%s strand:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s query_span:%d ref_span:%d total_span=%d query_string=%s hit_string=%s\n",
	  $hsp->score,
	    $hsp->strand,
	      $hsp->range("query"),
		$hsp->range("hit"),
		  $hsp->num_identical(),
		    $hsp->frac_identical("query"),
		      $hsp->length("query"),
			$hsp->length("hit"),
			  $hsp->length("total"),
			    $hsp->query_string(),
			      $hsp->hit_string() if $VERBOSE;

	#	    push @hsp, $hsp if hsp_filter($hsp, \%hsp_problems);
	push @hsp, $hsp;
	# TO DO: filter??

#      die $hsp->seq_inds("hit", "identical");
	# undef for this implementation??

#	die $hsp->gaps("hit");

	if ($VERBOSE) {
	  if (0) {
	    # only works with blat, not blast
	    printf STDERR "gap blocks:\n";
	    foreach my $k (qw(query hit)) {
	      printf STDERR "  %s\n", $k;
	      foreach my $gb (@{$hsp->gap_blocks($k)}) {
		printf STDERR "    %s\n", join ",", @{$gb};
	      }
	    }
	  }

	  printf STDERR  "    qs:%d qe:%d hs:%d he:%d\n",
	    $hsp->start("query"),
	      $hsp->end("query"),
		$hsp->start("hit"),
		  $hsp->end("hit");

	  my $aln = $hsp->get_aln();
	  printf STDERR "alignment sequences: %d\n", $aln->num_sequences;
	  foreach my $seq ($aln->each_seq) {
	    # Bio::LocatableSeq
	    printf STDERR "  %s: %d-%d\n", $seq->id, $seq->start, $seq->end;
	  }
	}

#	die scalar $hsp->seq_inds('query', 'identical', 'collapse' => 1);
	# implemented??

      }
    }
  }

  my $retval;
  if (@hsp) {

    if (0) {
      # disabled: BLASTer.pm by default is filtered to + strand hits only.
      # also in some modes (e.g. strand detection) we want to keep everything.
      my @passed;
      my @neg;
      foreach my $hsp (@hsp) {
	if ($hsp->strand() == -1) {
	  push @neg, $hsp;
	} else {
	  push @passed, $hsp;
	}
      }
      push @{$notes}, sprintf "%s: removed %d minus strand hit(s)", $end, scalar @neg if @neg;
      # e.g. FGFR3/TACC3 contig starts with a small ?inversion? mapped to -
      # which is a small subset of the larger hit (interestingly though
      # this is the only match to the first ~41 nt of FGFR3 end)
      @hsp = @passed;
    }

    if ($is_suppl and @hsp > 1) {
      # hack
      printf STDERR "supplemental sequence: using best hit only\n";
      @hsp = $hsp[0];
    }


    push @{$notes}, sprintf 'raw_hsp_count_%s=%d', $end, scalar @hsp;
    #
    #  TO DO: require some minimum score?
    #
    if (@hsp > 1) {
      push @{$notes}, sprintf 'hsp_count_%s=%d', $end, scalar @hsp;

      if ($VERBOSE) {
	foreach my $hsp (@hsp) {
#	  printf STDERR "  multi_hsp:  score:%s strand:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s query_span:%d ref_span:%d total_span=%d query_string=%s hit_string=%s gap_blocks_q=%s gap_blocks_h=%s\n",
	  printf STDERR "  multi_hsp:  score:%s strand_query:%s strand_hit:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s query_span:%d ref_span:%d total_span=%d query_string=%s hit_string=%s\n",
	    $hsp->score,
	      $hsp->strand("query"),
	      $hsp->strand("hit"),
		$hsp->range("query"),
		  $hsp->range("hit"),
		    $hsp->num_identical(),
		      $hsp->frac_identical("query"),
			$hsp->length("query"),
			  $hsp->length("hit"),
			    $hsp->length("total"),
			      $hsp->query_string(),
				$hsp->hit_string();
#				  join(",", map {@{$_}} @{$hsp->gap_blocks("query")}),
#				    join(",", @{$hsp->gap_blocks("hit")->[0]});
	}
      }

      unless ($multi_ok) {
	my @scores = map {$_->score} @hsp;
	foreach (my $i = 1; $i < @scores; $i++) {
	  # verify blat scores reported in descending order
	  printf STDERR "ERROR: blast scores not in descending order!\n"
	    if $scores[$i] > $scores[$i - 1];
	  # FIX ME:
	  # this could be indicative of a broader problem e.g. for ITDs,
	  # i.e. multiple HSPs rather than a single HSP containing gap blocks
	}
	@hsp = $hsp[0];
	# use best hit only, e.g.
	#
	# score:4.43 strand:1 q:2779-3030 subj:1-252 num_identical:252 frac_identical_query:1 query_span:252 ref_span:252 total_span=504
	# score:0.60 strand:1 q:1034-1067 subj:523-556 num_identical:34 frac_identical_query:1 query_span:34 ref_span:34 total_span=68
      }
    }

    if ($multi_ok) {
      $retval = \@hsp if @hsp;
    } elsif (@hsp == 1) {
      $retval = $hsp[0];
    } elsif (@hsp) {
      push @{$notes}, sprintf "ERROR_multiple_blast_hits_%s=%d", $end, scalar @hsp;
    }
  }
  return $retval;
}


sub get_hsp_indexes {
  my ($hsp) = @_;
  my ($i_start_acc, $i_end_acc) = $hsp->range("query");
  my ($i_start_contig, $i_end_contig) = $hsp->range("hit");
  foreach ($i_start_contig,
	   $i_end_contig,
	   $i_start_acc,
	   $i_end_acc) {
    $_--;
  }
  # TO DO: FIX ME, are these interbase coordinates???  (UCSC)

  return ($i_start_acc, $i_end_acc, $i_start_contig, $i_end_contig);
}

sub get_hsps_indexes {
  my ($hsps) = @_;
  my ($i_start_acc, $i_end_acc) = hsps_range($hsps, "query");
  my ($i_start_contig, $i_end_contig) = hsps_range($hsps, "hit");

#  printf STDERR "indexes: q=%d-%d hit=%d-%d\n", $i_start_acc, $i_end_acc, $i_start_contig, $i_end_contig;

  foreach ($i_start_contig,
	   $i_end_contig,
	   $i_start_acc,
	   $i_end_acc) {
    $_--;
  }
  # TO DO: FIX ME, are these interbase coordinates???  (UCSC)

  return ($i_start_acc, $i_end_acc, $i_start_contig, $i_end_contig);
}

sub generate_test_file {
  # generate single-SV test file for development
  my ($info) = @_;
  my @f = split /,/, $info;
  die "requires genea,geneb,tag,contig" unless @f == 4;
  my ($gene_a, $gene_b, $tag, $contig) = @f;
  my $f_out = sprintf '%s_%s_%s.tab', $tag, $gene_a, $gene_b;
  die "$f_out exists" if -e $f_out;

  my $rpt = new Reporter(
			 "-file" => $f_out,
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   geneA
					   geneB
					   contig
					)
				      ],
			 "-auto_qc" => 1,
			);
  my %r;
  $r{geneA} = $gene_a;
  $r{geneB} = $gene_b;
  $r{contig} = $contig;
  $rpt->end_row(\%r);
  $rpt->finish();
}

sub find_exonic {
  # find the index of the boundary of the current exon in the
  # specified direction
  my ($genomic, $idx2exon, $idx_start, $direction) = @_;

  my $code_start = $idx2exon->[$idx_start] || die;
  my $max = scalar @{$idx2exon};
  my $idx = $idx_start + $direction;
  # exclude start base as this refers to the alignment edge

  while (1) {
    my $i_next = $idx + $direction;
    if ($i_next < 0 or $i_next >= $max) {
      # next index is out of bounds, stop
      last;
    } else {
      my $code = $idx2exon->[$i_next] || die;
      printf STDERR "value at %d is %d\n", $i_next, $code if $VERBOSE;

      if ($code != $code_start) {
	# next base is a different exon, stop
#	die "exon change $code $code_start, final $idx";
	last;
      } else {
	$idx = $i_next;
	# OK, continue
      }
    }
  }

  my $result;
  # get sequence to exon edge, excluding
  if ($direction == -1) {
    # upstream of target
    $result = substr($genomic, $idx, $idx_start - $idx);
  } elsif ($direction == 1) {
    # downstream of target
    die unless $idx > $idx_start;
    printf STDERR "geneB down %d-%d\n", $idx_start + 1, $idx if $VERBOSE;
    $result = substr($genomic, $idx_start + 1, $idx - $idx_start);
    # TEST ME
  } else {
    die "unhandled direction";
  }
  return $result;
}

sub convert_fuzzion {
  my @infiles;

  if (my $set = $FLAGS{"convert-fuzzion"}) {
    push @infiles, split /,/, $set;
  }
  foreach my $file (@FUZZION_SRC_FILES) {
    my $src_files = read_simple_file($file);
    push @infiles, map {basename($_) . ".extended.tab"} @{$src_files};
  }
  # sources cn be mixed and matched

  foreach (@infiles) {
#    print STDERR "$_\n";
    die "where is $_" unless -s $_;
  }

  my $cicero_mode = $FLAGS{cicero};
  my %blacklist;
  my %blacklist_removed;
  unless ($FLAGS{"no-blacklist"}) {
    my $f_bl = $FLAGS{"blacklist"} || die "specify -blacklist FILE | -no-blacklist\n";
    my $set = read_simple_file($f_bl);
    foreach my $fusion (@{$set}) {
      my @genes = split /_/, $fusion;
      die "fusion $fusion doesn't have exactly 2 genes" unless @genes == 2;
      $blacklist{$genes[0]}{$genes[1]} = 1;
      $blacklist{$genes[1]}{$genes[0]} = 1;
      # also reciprocal
    }
  }


  my @fields = (
		$F_EXTENDED_CHUNK,

#		"extended_exon",
#		"extended_max"
		# not sure if these sets are even used anymore
	       );

  my $REPORT_SAMPLE = 1;
  my $REPORT_SOURCE = 1;

  # output is only a subset of the inputs, as these will be
  # in varying formats depending on origin and conversion method
  my $found_extra_source;

  foreach my $field (@fields) {
    printf STDERR "processing type %s\n", $field;

    my $base = @infiles > 1 ? "merged" : basename($infiles[0]);
    my $outfile = sprintf "%s.fuzzion_%s.tab", $base, $field;

#    my $wf = new WorkingFile($outfile);
#    my $fh = $wf->output_filehandle();

    my %pattern_seq;
    my %pattern_annot;

    my $rpt_patterns = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => \@H_FZ2_CORE,
			 "-auto_qc" => 1,
			);

    # while (my $row = $df->next("-ref" => 1)) {  # headerless

    my %saw;
    my %saw2;
    my %saw2_pid;
    my $duplicates_contig = 0;
    my %counter;
    my %errors;

    my $INCLUDE_MD5 = 0;
    # debug

    foreach my $infile (@infiles) {
      printf STDERR "  file: %s\n", $infile;
      my $infile_tag = file_to_tag($infile);

      my $df = new DelimitedFile("-file" => $infile,
				 "-headers" => 1,
				);

      my $rows_in = [];
      while (my $row = $df->get_hash()) {
	push @{$rows_in}, $row;
      }
      # FIX ME: this uses a lot of RAM.

      my $of_pattern = sprintf '%s.%s.pattern.tab', basename($infile), $field;
      my $rpt = $df->get_reporter(
				  "-file" => $of_pattern,
				  "-extra" => [
					       qw(
						   pattern
						)
					      ]
				 );
      # Annotate intermediate file with final pattern IDs.
      # There is a separate file for each pattern length style being generated.

      if ($FUSION_PREFER_BEST_CONTIG_MATCHES) {
	$rows_in = filter_by_contig_identity(
					     "-rows" => $rows_in,
					     "-df" => $df,
					     "-rpt-suppressed" => $rpt
					    );
      }
#      die scalar @{$rows_in};


      foreach my $row (@{$rows_in}) {
	dump_die($row, "new record", 1) if $VERBOSE;
	dump_die($row, "ERROR: no field $field") unless exists $row->{$field};

	my $track_annot = sub {
	  my ($pid) = @_;

	  if ($REPORT_SAMPLE) {
	    my $sample = $row->{$F_SAMPLE};
	    unless (defined $sample) {
	      if ($FLAGS{"null-sample"}) {
		$sample = "";
	      } else {
		dump_die($row, "no field $F_SAMPLE in $infile");
	      }
	    }
	    $pattern_annot{$pid}{$F_SAMPLE}{$sample} = 1 if $sample;
	  }

	  if ($REPORT_SOURCE) {
	    #
	    #  annotate source file pattern came from:
	    #
	    my $tag = $infile_tag;

	    if (@FUZZION_SOURCE_APPEND_FIELDS) {
	      # additional columns to include, e.g. project
	      my @extra;
	      foreach my $f (@FUZZION_SOURCE_APPEND_FIELDS) {
		if (my $v = $row->{$f}) {
		  $found_extra_source = 1;
		  push @extra, $v;
		}
	      }
	      if (@extra) {
		die "multiple supplemental info fields" if @extra > 1;
		$tag .= "/" . $extra[0];
	      }
	    }

	    $pattern_annot{$pid}{$F_SOURCE}{$tag} = 1;
	  }

	  cicero_clone_patch($row) if $cicero_mode;

	  foreach my $f (@H_FZ2_CLONE) {
	    dump_die($row, "can't clone input field $f; consider -cicero option?") unless exists $row->{$f};
	    my $v = $row->{$f};
	    if (length $v) {
	      $pattern_annot{$pid}{$f}{$v} = 1;
	    }
	  }

	};

	my $contig = $row->{$field};
	# skip records where no contig provided
	unless ($contig) {
	  $row->{pattern} = "";
	  $rpt->end_row($row);
	  next;
	}

	my $contig_uc = uc($contig);
	# contig sequence is reported in uppercase, rest of extended sequence
	# is in lowercase.  Uppercase to standardize, sometimes there are
	# differences in reported contig boundaries even though the
	# fusion event is identical, e.g.:
	#
	# PML   RARA    catctactgccgaggatgttccaagccgctgtgctgctcgtgcgcgctccttgacagcagccacagtgagctcaagtgcgacatcagcgcagagatccagcagcgacaggaggagctggacgccatgacgcaggcgctgcaggagcaggatagtgcctttggcgcggttcacgcgcagatgcacgcggccgtcggccagctgggccgcgcgcgtgccgagaccgaggagctgatccgcgagcgcgtgcgccaggtggtagctcacgtgcgggctcaggagcgcgagctgctggaggctgtggacgcgcggtaccagcgcgactacgaggagatggccagtcggctgggccgcctggatgctgtgctgcagcgcatccgcacgggcagcgcgctggtgcagaggatgaagtgctacgcctcggaccaggaggtgctggacatgcacggtttcctgcgccaggcgctctgccgcctgcgccaggaggaGCCCCAGAGCCTGCAAGCTGCCGTGCGCACCGATGGCTTCGACGAGTTCAAGGTGCGCCTGCAGGACCTCAGCTCTTGCATCACCCAGGGGAAAGCCATTGAGACCCAGAGCAGCAGTTCTGAAGAGATAGTGCCCAGCCCTCCCTCGCCACCCCCTCTACCCCGCATCTACAAGCCTTGCTTTGTCTGTCAGGACAAGTCCTCAGGCTACCACTATGGGGTCAGCGCCTGTGAGGGCTGCAAGGGCTTCTTCCGCCGCAGCATCCAGAAGAACATGGTGTACACGTGTCACCGGGACAAGAACTGCATCATCAACAAGGTGACCCGGAACCGCTGCCAGtactgccgactgcagaagtgctttgaagtgggcatgtccaaggagt
        # PML     RARA    catctactgccgaggatgttccaagccgctgtgctgctcgtgcgcgctccttgacagcagccacagtgagctcaagtgcgacatcagcgcagagatccagcagcgacaggaggagctggacgccatgacgcaggcgctgcaggagcaggatagtgcctttggcgcggttcacgcgcagatgcacgcggccgtcggccagctgggccgcgcgcgtgccgagaccgaggagctgatccgcgagcgcgtgcgccaggtggtagctcacgtgcgggctcaggagcgcgagctgctggaggctgtggacgcgcggtaccagcgcgactacgaggagatggccagtcggctgggccgCCTGGATGCTGTGCTGCAGCGCATCCGCACGGGCAGCGCGCTGGTGCAGAGGATGAAGTGCTACGCCTCGGACCAGGAGGTGCTGGACATGCACGGTTTCCTGCGCCAGGCGCTCTGCCGCCTGCGCCAGGAGGAGCCCCAGAGCCTGCAAGCTGCCGTGCGCACCGATGGCTTCGACGAGTTCAAGGTGCGCCTGCAGGACCTCAGCTCTTGCATCACCCAGGGGAAAGCCATTGAGACCCAGAGCAGCAGTTCTGAAGAGATAGTGCCCAGCCCTCCCTCGCCACCCCCTCTACCCCGCATCTACAAGCCTTGCTTTGTCTGTCAGGACAAGTCCTCAGGCTACCACTATGGGGTCAGCGCCTGTGAGGGCTGCAAGGGCTTCTTCCGCCGCAGCatccagaagaacatggtgtacacgtgtcaccgggacaagaactgcatcatcaacaaggtgacccggaaccgctgccagtactgccgactgcagaagtgctttgaagtgggcatgtccaaggagt

	my $md5 = md5_hex($contig);

	my $genea = $row->{genea_symbol} || dump_die($row, "no genea_symbol");
	my $geneb = $row->{geneb_symbol} || dump_die($row, "no geneb_symbol");
	# transcript pair-specific symbols.
	# - the upside of using these is that they remove a lot of
	#   potential junk in the raw entries, resolving lists to individual
	#   genes for each transcript, and removing quotes and refFlat
	#   "sharp" suffixes (locA, etc)
	# - the downside is that use of discrete symbols can introduce
	#   duplicates where source is valid list of multiple genes,
	#   e.g. gene families.  While these are filtered below, the
	#   pattern IDs will imply a single gene when in fact they may
	#   apply to multiple symbols.

	if ($blacklist{$genea}{$geneb}) {
	  my $key = join "_", $genea, $geneb;
	  $blacklist_removed{$key}++;
	  $row->{pattern} = "suppressed_by_blacklist";
	  $rpt->end_row($row);
	  next;
	}

	dump_die($row, "where is $field") unless exists $row->{$field};

	my $notes = $row->{$F_NOTES};
	dump_die($row, "missing notes field") unless defined $notes;

	my $intragenic_mode = $genea eq $geneb;
	my $brk_a = $intragenic_mode ? "}" : "]";
	my $brk_b = $intragenic_mode ? "{" : "[";
	# for intragenic events always use {} brackets so that
	# spanning the breakpoint is required

	my ($a_c_start,
	    $a_c_end,
	    $b_c_start,
	    $b_c_end) = @{$row}{qw(
				    genea_contig_start
				    genea_contig_end
				    geneb_contig_start
				    geneb_contig_end
				 )};
	die "genea/b contig start/end not defined" unless
	  defined($a_c_start) and defined($a_c_end) and
	    defined($b_c_start) and defined($b_c_end);

	my $gap = $b_c_start - $a_c_end - 1;
	if ($gap < 0) {
	  printf STDERR "ERROR: negative gap in contig between genes A (%s) and B (%s), are genes reversed?\n", $genea, $geneb;
	  $errors{genes_reversed_in_contig}++;

	  $row->{pattern} = "";
	  $rpt->end_row($row);
	  next;
	}

	if ($FUSION_MAX_INTERSTITIAL_LENGTH and
	    $gap > $FUSION_MAX_INTERSTITIAL_LENGTH) {
	  # reduce overly-long interstitial region
	  if ($FUSION_INTERSTITIAL_TRIM_POLICY eq "B") {
	    # preserve B breakpoint site
	    $a_c_end = $b_c_start - $FUSION_MAX_INTERSTITIAL_LENGTH - 1;
	    $gap = $b_c_start - $a_c_end - 1;
	  } else {
	    die "unimplemented interstitial trim policy $FUSION_INTERSTITIAL_TRIM_POLICY";
	  }
	}

	#
	#  mark breakpoint site:
	#
	# mixed upper/lowercase sequence (w/contig UC) required
	# for this operation:
	if ($gap) {
	  # contain contains interstitial sequence between genes A and B
#	  $contig =~ s/^([a-z]+[A-Z]{$a_c_end})([A-Z]{$gap})/$1\]$2\[/ || die;
	  $contig =~ s/^([a-z]+[A-Z]{$a_c_end})([A-Z]{$gap})/$1${brk_a}$2${brk_b}/ || die;
	} else {
	  # flush border between genes A and B
#	  $contig =~ s/^([a-z]+[A-Z]{$a_c_end})/$1\]\[/ || die;
	  dump_die($row, "flush border", 1);
	  $contig =~ s/^([a-z]+[A-Z]{$a_c_end})/$1${brk_a}${brk_b}/ || die "can't mark border at $a_c_end in $contig";
	}
	my $final = uc($contig);
	# final pattern sequence with brackets

	#
	#  duplicate checking:
	#
	if (my $pid = $saw{$genea}{$geneb}{$final}) {
	  # typically happens if a different combination of isoforms
	  # for this gene pair produces the same pattern
	  $duplicates_contig++;
	  $row->{pattern} = $pid;
	  &$track_annot($pid);
	  $rpt->end_row($row);
	  next;
	}

	if (my $old = $saw2{$final}) {
	  # this can happen if two different gene pairs yield the same pattern,
	  # e.g. MZT2A-BTBD2 and MZT2B-BTBD2.
	  # Much more likely to happen if we use transcript-pair-specific
	  # gene symbols rather than the raw symbol, which may contain lists.
	  printf STDERR "NOTE, duplicate contig:%s\n", join "\t",
	    $genea, $old->{genea_symbol},
	      $geneb, $old->{geneb_symbol},
	      $contig_uc, $old->{$field};
	  my $pid = $saw2_pid{$final} || die "can't get old PID";
	  # pattern ID may refer to a different gene pair...
	  $duplicates_contig++;
	  &$track_annot($pid);
	  # ...however we won't lose any annotations
	  $row->{pattern} = $pid;
	  $rpt->end_row($row);
#	  dump_die($row, "test me: after tracking change to full pattern");
	  # should be OK, just check this
	  next;
	}
	$saw2{$final} = $row;


	#
	#  assign new pattern ID:
	#
	my $count = ++$counter{$genea}{$geneb};
	my $contig_id = sprintf '%s-%s-%02d', $genea, $geneb, $count;
	dump_die($row, "new contig $contig_id", 1) if $VERBOSE;
	printf STDERR "assign %s\n", join ",", $genea, $geneb, $contig_uc;

	#
	#  add to duplicate tracking:
	#
	#	$saw{$genea}{$geneb}{$contig_uc} = $contig_id;
	$saw{$genea}{$geneb}{$final} = $contig_id;
	$saw2_pid{$final} = $contig_id;

	if ($VERBOSE) {
	  printf STDERR "debug %s\n", join ",", $a_c_start, $a_c_end, $b_c_start, $b_c_end, $gap;
	  printf STDERR "before: %s\n", $contig;
	}

	if (1) {
	  $pattern_seq{$contig_id} = $final;
	  &$track_annot($contig_id);
	} else {
	  # original format: write pattern immediately
	  # my @out;
	  # push @out, $contig_id;
	  # push @out, uc($contig);

	  # #	printf STDERR "new_record %s\n", join "\t", $genea, $geneb, $md5;

	  # if ($INCLUDE_MD5) {
	  #   push @out, $md5;
	  #   push @out, md5_hex(uc($contig));
	  # }
	  # printf $fh "%s\n", join "\t", @out;
	}

	$row->{pattern} = $contig_id;
	$rpt->end_row($row);
	# record pattern ID in source file
      }

      $rpt->finish();
      # each intermediate file for this pattern length
    }

    #
    #  write patterns and annotations:
    #
    foreach my $pid (sort keys %pattern_seq) {
      my %r;
      $r{pattern} = $pid;
      $r{sequence} = $pattern_seq{$pid};

      if ($REPORT_SOURCE) {
	$r{source} = join ",", sort keys %{$pattern_annot{$pid}{source}};
      }
      if ($REPORT_SAMPLE) {
	$r{$F_SAMPLE} = join ",", sort keys %{$pattern_annot{$pid}{$F_SAMPLE}};
      }

      foreach my $f (@H_FZ2_CLONE) {
	$r{$f} = join ",", sort keys %{$pattern_annot{$pid}{$f}};
      }

      $rpt_patterns->end_row(\%r);
    }

#    $wf->finish();
    $rpt_patterns->finish();
    # main output file for this pattern length

    printf STDERR "skipped %d duplicate source contigs\n", $duplicates_contig;

    if (%errors) {
      printf STDERR "ERRORS generating type %s:\n", $field;
      foreach my $k (sort keys %errors) {
	printf STDERR "  %s: %d\n", $k, $errors{$k};
      }
    }

  }  # file type


  if (%blacklist_removed) {
    printf STDERR "removed blacklisted fusions: %s\n", join ", ", sort keys %blacklist_removed;
  }

  if (@FUZZION_SOURCE_APPEND_FIELDS) {
    die "-source-append fields, but not found" unless $found_extra_source;
  }


}

sub get_gene_a_b {
  my ($row) = @_;
  my ($gene_a, $gene_b);
  if ($F_FUSION) {
    my $fusion = $row->{$F_FUSION} || dump_die($row, "no field $F_FUSION");
    my @f = split /_/, $fusion;
    die "$fusion not 2 genes" unless @f == 2;
    ($gene_a, $gene_b) = @f;
  } else {
    foreach my $f ($F_GENE_A, $F_GENE_B) {
      dump_die($row, "no field $f") unless exists $row->{$f};
    }
    $gene_a = $row->{$F_GENE_A};
    $gene_b = $row->{$F_GENE_B};
  }
  return ($gene_a, $gene_b);
}

sub repair_gene_swaps {
  # from an extended output file, generate a new input file repairing
  # swapped genes identified in previous run.  This is suitable for
  # re-running with cleaner output that will work with Fuzzion2 conversion.
  my $f_in = $FLAGS{"repair-gene-swaps"} || die;
  my $f_out = $f_in . ".gene_swap_repair.tab";

  my $df = new DelimitedFile("-file" => $f_in,
			     "-headers" => 1,
			     );
  my @h;
  foreach my $h (@{$df->headers_raw}) {
    last if $h eq $NEW_HEADERS[0];
    # write in original file format
    push @h, $h;
  }

  my $rpt = new Reporter(
			 "-file" => $f_out,
			 "-delimiter" => "\t",
			 "-labels" => \@h,
			 "-auto_qc" => 1,
			);
  my $repaired = 0;
  my $unrepairable = 0;
  while (my $row = $df->get_hash()) {
    my @notes = split /,/, $row->{$F_NOTES};
    if (grep {$_ eq $TAG_ERROR_GENES_REVERSED} @notes) {
      my $gene_a = $row->{$F_GENE_A} || dump_die($row, "no $F_GENE_A");
      my $gene_b = $row->{$F_GENE_B} || dump_die($row, "no $F_GENE_B");
      if ($gene_a eq $gene_b) {
	$unrepairable++;
      } else {
	$row->{$F_GENE_A} = $gene_b;
	$row->{$F_GENE_B} = $gene_a;
	$repaired++;
      }
    }
    $rpt->end_row($row);
  }
  $rpt->finish();

  printf STDERR "repaired:%d unrepairable:%d\n", $repaired, $unrepairable;
}

sub transcript_filter {
  my ($transcript2nm, $gene, $nm) = @_;

  if (not($nm) and $FLAGS{"restrict-nm-sjpi"}) {
    # want to filter by preferred isoform
    if (my $pref = $SJPI->get_preferred_isoform($gene)) {
      if (exists $transcript2nm->{$pref}) {
	# OK: preferred isoform available in our transcript set
	$nm = $pref;
      } else {
	die "isoform filter failure for $gene, can't filter to $pref, no data\n";
      }
    }
  }

  if ($nm) {
    my @delete = grep {$_ ne $nm} keys %{$transcript2nm};
    delete @{$transcript2nm}{@delete};
  }
}

sub merge_fuzzion_ids {
  my $f_orig = $FLAGS{"merge-pattern-ids-orig"} || die;
  my $f_annot = $FLAGS{"merge-pattern-ids-annotated"} || die;

  my $df = new DelimitedFile("-file" => $f_orig,
			     "-headers" => 1,
			     );
  my $headers = $df->headers_raw();
  # use all headers in original file as a key

  $df = new DelimitedFile("-file" => $f_annot,
			  "-headers" => 1,
			 );
  # gather pattern IDs from intermediate file, where multiple rows
  # appear for different isoform pairs and there may be multiple pattern IDs
  my %patterns;
  while (my $row = $df->get_hash()) {
    die unless exists $row->{pattern};
    if (my $pid = $row->{pattern}) {
      my $key = get_md5_key($row, $headers);
      $patterns{$key}{$pid} = 1;
    }
  }

  #
  #  annotate original file with all matching pattern IDs for each row
  #
  $df = new DelimitedFile("-file" => $f_orig,
			  "-headers" => 1,
			 );
  my $outfile = basename($f_orig) . ".patterns.tab";
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => [
					   qw(
					       patterns
					    )
					  ],
  			      "-auto_qc" => 1,
			     );

  while (my $row = $df->get_hash()) {
    my $key = get_md5_key($row, $headers);

    my $v = "";
    # not every row is guaranteed to have a pattern
    if (my $patterns = $patterns{$key}) {
      $v = join ",", sort keys %{$patterns};
    }
    $row->{patterns} = $v;
    $rpt->end_row($row);
  }
  $rpt->finish();

}

sub get_md5_key {
  my ($row, $headers) = @_;
  return md5_hex(@{$row}{@{$headers}});
}

sub find_duplicate_blocks_via_query_gap {
  # what positions in the HIT (contig) represent gaps in the QUERY (NM_)?
  # e.g.
  # /home/medmonso/work/steve/2020_06_13_cicero_contig_extension/2021_01_04_all/blast/gapped.tab
  my ($hsps, $contig) = @_;
  $contig = uc($contig);

  my @duplicate_hit_blocks;

  foreach my $hsp (@{$hsps}) {
    if (my $ranges = find_hit_ranges_corresponding_to_query_gaps($hsp)) {
      foreach my $range (@{$ranges}) {
	my ($start_pos, $end_pos) = @{$range};
	push @duplicate_hit_blocks, [ $start_pos, $end_pos ];
	# create a copy as this may be adjusted later

	my $start_idx = $start_pos - 1;
	my $end_idx = $end_pos - 1;
	my $len = ($end_pos - $start_pos) + 1;

	my $seq = substr($contig, $start_idx, $len);
	# putative duplicated sequence
	print STDERR "chunk $seq\n";

	my @other_idx;
	my $ptr = 0;
	while (1) {
	  my $idx = index($contig, $seq, $ptr);
	  last if $idx == -1;
	  push @other_idx, $idx unless $idx == $start_idx;
	  $ptr = $idx + $len;
	}

	my $found;

	if (@other_idx) {
	  #
	  #  direct match of insert sequence elsewhere in contig
	  #
	  die "multiple alt indexes" if @other_idx > 1;
	  # happens?
	  foreach my $alt_idx (@other_idx) {
	    if ($alt_idx > $start_idx) {
	      my $distance = $alt_idx - $end_idx;
	      die if $distance <= 0;
	      die "distance $distance" if $distance > 1;
	      # behavior TBD if not flush
	      push @duplicate_hit_blocks, [ $alt_idx + 1, $alt_idx + $len ];
	      $found = 1;
	      # 1-based position of additional duplicate
	    } else {
	      # TBD
	      die "putative duplication before blast gap";
	    }
	  }
	}

	unless ($found) {
	  # check the block immediately following: how many bases
	  # can we go before a mismatch?
	  my $next = substr($contig, $end_idx + 1, $len);

	  my $matches = 0;
	  for (my $i = 0; $i < $len; $i++) {
	    if (substr($seq, $i, 1) eq substr($next, $i, 1)) {
	      $matches++;
	    } else {
	      last;
	    }
	  }

	  my $mismatch_len = $len - $matches;
#	  if ($mismatch_len <= $MAX_INSERTION_DUPLICATE_MISMATCHES) {
	  if ($matches >= $MIN_INSERTION_DUPLICATE_MATCHES) {
	    $duplicate_hit_blocks[0]->[1] -= $mismatch_len;

	    my $block_start = $end_pos + 1;
	    my $block_end = $block_start + $matches - 1;

	    push @duplicate_hit_blocks, [ $block_start, $block_end ];
	    $found = 1;
	  }
	}

	if (not($found)) {
	  if ($len < $MIN_QUERY_GAP_TO_CONSIDER_ITD_CRITICAL) {
	    @duplicate_hit_blocks = ();
	    # consider spurious result, ignore
	    # e.g. /home/medmonso/work/steve/2020_06_13_cicero_contig_extension/2020_12_23_high_confidence_itd_and_del/14
	    # - small gap but sequence doesn't repeat, likely irrelevant
	  } else {
	    # need to study: code might need more work?
	    confess "can't find downstream matches for insert/repeat sequence";
	  }
	}

      }
    }
  }

  return @duplicate_hit_blocks ? \@duplicate_hit_blocks : undef;
}

sub find_hit_ranges_corresponding_to_query_gaps {
  # what positions in the HIT (contig) represent gaps in the QUERY (NM_)?
  my ($hsp) = @_;
  my @hit_gaps;

  if ($hsp->gaps("query")) {
    my $qs = $hsp->query_string();
    my $hs = $hsp->hit_string();

    my $hit_pos = ($hsp->range("hit"))[0];
    my $qlen = length $qs;

    my $in_gap;
    my $hit_gap;

    my $ai;
    # alignment index

    for ($ai = 0; $ai < $qlen; $ai++) {
      if (substr($qs, $ai, 1) eq "-") {
	unless ($hit_gap) {
	  # new gap
	  $hit_gap = [$hit_pos, $hit_pos];
	  push @hit_gaps, $hit_gap;
	}
	$hit_gap->[1] = $hit_pos;
	# continuation
      } else {
	$hit_gap = undef;
      }

      if (substr($hs, $ai, 1) eq "-") {
#	die "gap in hit current hit_pos = $hit_pos";
      } else {
	# aligned base
	$hit_pos++;
      }
    }
  }

  if (@hit_gaps) {
    my @filtered;
    foreach my $gap (@hit_gaps) {
      my $size = ($gap->[1] - $gap->[0]) + 1;
      push @filtered, $gap if $size >= $MIN_QUERY_GAP_TO_CONSIDER_ITD_OR_INSERTION;
    }
    @hit_gaps = @filtered;
  }

  return @hit_gaps ? \@hit_gaps : undef;
}

sub query_to_hit {
  # for a given query base number, what are corresponding hit base numbers?
  # (compensating for gaps in the alignment)
  my ($hsp, @q_positions) = @_;

  my $q_aln = $hsp->query_string();
  my $h_aln = $hsp->hit_string();

  my $q_pos = $hsp->start("query");
  my $h_pos = $hsp->start("hit");

  my $ai = 0;
  my $alen = length $q_aln;

  my %q2h;

  for (my $ai = 0; $ai < $alen; $ai++) {
    my $q_gap = substr($q_aln, $ai, 1) eq "-";
    my $h_gap = substr($h_aln, $ai, 1) eq "-";

    $q2h{$q_pos} = $h_pos unless $q_gap;

    $q_pos++ unless $q_gap;
    $h_pos++ unless $h_gap;
  }
  return map {$q2h{$_}} @q_positions;
}

sub generate_contigs {
  # find the genomic interval
  # find gene model(s) where an exon overlaps, or die
  # generate contig net of event
  # compensate for strand?
  # warn/die if too close to exon edge
  my $f_in = $FLAGS{"generate-contigs"} || die "-generate-contigs";
  my $fai = get_fai();
  my $rf = get_refflat();
  my $min_flank = 40;
  # TO DO: set global variable
  my $f_out = basename($f_in) . ".contig.tab";

  my $df = new DelimitedFile("-file" => $f_in,
			     "-headers" => 1,
			     );
  my $rpt = $df->get_reporter(
			      "-file" => $f_out,
			      "-extra" => [
					   $F_GENE_A,
					   "chrA",
					   "posA",
					   $F_GENE_B,
					   "chrB",
					   "posB",
					   $F_CONTIG,
					   "contig_debug",
					  ],
  			      "-auto_qc" => 1,
			     );

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
    #
    #  find genomic interval:
    #
    my $spec = $row->{"event_spec"} || die;
    my @f = split /::/, $spec;
    die unless @f == 2;

    my %chr;
    my ($start, $end);
    foreach my $side (@f) {
      # chr3:41266090(+)
      $side =~ /^(\w+):(\d+)\(([\+\-])\)$/ || die "can't match spec regexp in $side";
      my ($chr, $pos, $strand) = ($1, $2, $3);
      $chr{$chr} = 1;

      if (defined $start) {
	$end = $pos;
      } else {
	$start = $pos;
      }
    }
    die "only single chrom interval currently supported" unless scalar keys %chr == 1;
    my ($chr) = keys %chr;
    printf STDERR "target interval: %s:%d-%d\n", $chr, $start, $end;

    #
    #  find gene model exons containing the interval:
    #
    my $gene = $row->{gene} || die "-gene";
    my $rows = $rf->find_by_gene($gene) || die "no hits for $gene";
    # TO DO: add GSM

    my %ex_intervals;

    my %prev_intervals;
    my %next_intervals;

    # TO DO: filter to SJPI?

    foreach my $r (@{$rows}) {
      die "only + currently supported" unless $r->{strand} eq "+";

      my $exons = $r->{exons} || die;

      for (my $ei = 0; $ei < @{$exons}; $ei++) {
	my $ex = $exons->[$ei];
	printf STDERR "exon: %d-%d\n", @{$ex}{qw(start end)};

	if ($start >= $ex->{start} and $end <= $ex->{end}) {
	  # hit
	  my $key = join "-", @{$ex}{qw(start end)};
	  $ex_intervals{$key} = 1;

	  if ($ei > 0) {
	    $key = join "-", @{$exons->[$ei - 1]}{qw(start end)};
	    print STDERR "prev interval $key\n";
	    $prev_intervals{$key} = 1;
	  }

	  if (($ei + 1) < @{$exons}) {
	    $key = join "-", @{$exons->[$ei + 1]}{qw(start end)};
	    print STDERR "next interval $key\n";
	    $next_intervals{$key} = 1;
	  }

	  last;
	}
      }
    }

    die "no exon intervals" unless %ex_intervals;
    printf STDERR "ERROR: multiple exon intervals: %s\n", join ",", keys %ex_intervals if keys %ex_intervals > 1;

    my ($ex_interval) = sort keys %ex_intervals;
    my ($ex_start, $ex_end) = split /\-/, $ex_interval;

    #
    #  generate contig:
    #
    my $event_type = $row->{event_type} || die;
    my $contig;
    my $contig_debug;

    #
    #  find sequence before and after the event, and extend if necessary:
    #
    my ($before, $after);
    if ($event_type eq "del" or $event_type eq "itd") {
      $before = $fai->get_chunk(
				   "-id" => $chr,
				   "-start" => $ex_start,
				   "-end" => $start - 1
				  );

      $after = $fai->get_chunk(
				   "-id" => $chr,
				   "-start" => $end + 1,
				   "-end" => $ex_end
				  );

      printf STDERR "flanking before:%d after:%d\n", length($before), length($after);
      printf STDERR "before: %s\nafter: %s\n", $before, $after;

      if (length($before) < $min_flank) {
	#
	#  extend leading sequence into the previous exon
	#
	die "multi prev" if keys %prev_intervals > 1;
	my ($prev_start, $prev_end) = split /-/, (keys %prev_intervals)[0];
	my $needed = $min_flank - length($before);

	my $start = ($prev_end - $needed) + 1;
	die if $start < $prev_start;

	my $extra = $fai->get_chunk(
				    "-id" => $chr,
				    "-start" => $start,
				    "-end" => $prev_end
				   );

	$before = $extra . $before;
      }

      if (length($after) < $min_flank) {
	#
	#  extend trailing sequence into the next exon
	#
	die "multi next" if keys %next_intervals > 1;
	my ($next_start, $next_end) = split /-/, (keys %next_intervals)[0];
	my $needed = $min_flank - length($after);
	my $end = $next_start + $needed - 1;
	die if $end > $next_end;

	my $extra = $fai->get_chunk(
				    "-id" => $chr,
				    "-start" => $next_start,
				    "-end" => $end);
	$after .= $extra;
      }
    } else {
      die;
    }

    if ($event_type eq "del") {
      $contig = join "", $before, $after;
      $contig_debug = $before . "/" . $after;
    } elsif ($event_type eq "itd") {
      my $dup = $fai->get_chunk(
				"-id" => $chr,
				"-start" => $start,
				"-end" => $end
			       );

      $contig = join "", $before, $dup, $dup, $after;
      $contig_debug = join "/", $before, $dup, $dup, $after;
    } else {
      die "unhandled event type $event_type";
    }

    $row->{$F_GENE_A} = $gene;
    $row->{$F_GENE_B} = $gene;
    $row->{chrA} = $row->{chrB} = $chr;
    $row->{posA} = $start;
    $row->{posB} = $end;
    $row->{$F_CONTIG} = $contig;
    $row->{contig_debug} = $contig_debug;

    $rpt->end_row($row);
  }
  $rpt->finish();

}

sub get_refflat {
  my ($return_filename) = @_;

  my $f_refflat = config_or_manual(
				 "-config-type" => "genome",
				 "-config-name" => $FLAGS{genome},
				 "-parameter" => "REFSEQ_REFFLAT",
				 "-manual" => $FLAGS{refflat}
				);

  if ($return_filename) {
    return $f_refflat;
  } else {
    my $rf = new RefFlatFile();
    $rf->canonical_references_only(1);
    $rf->strip_sharp_annotations(1);
    $rf->parse_file(
		    "-refflat" => $f_refflat,
		    "-type" => "refgene",
		   );
    return $rf;
  }
}


sub resolve_block_overlap {
  my ($dbs, $hsps) = @_;

  #
  #  check to see if duplicate blocks overlap, which can happen
  #  due to indel-y alignment, particularly at edges
  #
  my $d1 = new Set::IntSpan(join "-", @{$dbs->[0]});
  my $d2 = new Set::IntSpan(join "-", @{$dbs->[1]});
  my $di = intersect $d1 $d2;
  if ($di->size()) {
    printf STDERR "conflict between dup blocks %d-%d and %d-%d\n", @{$dbs->[0]}, @{$dbs->[1]};

    die "wheeze" unless $dbs->[0]->[1] >= $dbs->[1]->[0];
    # verifya that block 1 end overlaps block 2 start

    die "gurgle" unless $dbs->[0]->[0] <= $dbs->[1]->[0];
    # verify block 1 starts at or before block2

    #
    #  resolve conflicts by prioritizing the block with the
    #  cleanest alignment.  Hacky as the issue is more the
    #  alignment in the contested area, but may be good enough.
    #
    my $fiq_1 = $hsps->[0]->frac_identical("query");
    my $fiq_2 = $hsps->[1]->frac_identical("query");

    if ($fiq_1 >= $fiq_2) {
      # first block equal or better quality than the second
      $dbs->[1]->[0] = $dbs->[0]->[1] + 1;
      # make second block start after first ends
      printf STDERR "revise second dup block to %d-%d\n", @{$dbs->[1]};
    } else {
      # second block is better
      $dbs->[0]->[1] = $dbs->[1]->[0] - 1;
      # back off first block to before second block starts
      printf STDERR "revise first dup block to %d-%d\n", @{$dbs->[0]};
    }
  }
}


sub get_trimmed_range {
  # trim HSP range start/ends by excluding gaps and mismatches within 10 nt
  my ($hsp, $type, $return_array) = @_;

  die unless $type and ($type eq "hit" or $type eq "query");

  my $type_hit = $type eq "hit";

  my $q_aln = $hsp->query_string();
  my $h_aln = $hsp->hit_string();

  printf STDERR "aligns:\n%s\n%s\n", $q_aln, $h_aln;

  #
  #  skip any gappy portion at start of hit region:
  #
  my $final_start = $hsp->start($type);
  my $pos = $final_start;
  my $max = $pos + $HIT_TRIM_END_RANGE;
  my $ai = 0;
  my $bad;
  while ($pos < $max) {
    my $q_base = substr($q_aln, $ai, 1);
    my $h_base = substr($h_aln, $ai, 1);
#    $final_start = $pos + 1 if $q_base eq '-' or $h_base eq '-' or $q_base ne $h_base;
    if ($q_base eq '-' or $h_base eq '-' or $q_base ne $h_base) {
      # alignment disagrees at ths position
      $bad = 1;
    } elsif ($bad) {
      # alignment agrees at this position, and last comparison was bad
      $final_start = $pos;
      $bad = 0;
    }
    $pos++ unless ($type_hit ? $h_base : $q_base) eq '-';
    $ai++;
  }

  #
  #  skip any gappy portion at end of hit region:
  #
  my $final_end = $hsp->end($type);
  $pos = $final_end;
  my $min = $pos - $HIT_TRIM_END_RANGE;
  $ai = length($q_aln) - 1;
  $bad = 0;
  while ($pos > $min) {
    my $q_base = substr($q_aln, $ai, 1);
    my $h_base = substr($h_aln, $ai, 1);
#    $final_end = $pos - 1 if $q_base eq '-' or $h_base eq '-' or $q_base ne $h_base;
    if ($q_base eq '-' or $h_base eq '-' or $q_base ne $h_base) {
      # alignment disagrees at ths position
      $bad = 1;
    } elsif ($bad) {
      # alignment agrees at this position, and last comparison was bad
      $final_end = $pos;
      $bad = 0;
    }

    $pos-- unless ($type_hit ? $h_base : $q_base) eq '-';
    $ai--;
  }

  printf STDERR "raw %s range:%d-%d  trimmed:%d-%d\n",
    $type, $hsp->range($type), $final_start, $final_end;

#  die "FIX ME: end trimmed, need QC for new code" unless $final_end == $hsp->end($type);

  return $return_array ? ($final_start, $final_end) :
    [ $final_start, $final_end ];
  # require explicit type rather than relying on wantarray()
}

sub remove_list_duplicates {
  my ($list) = @_;
  my %saw;
  my @out;
  foreach my $entry (@{$list}) {
    next if $saw{$entry}++;
    push @out, $entry;
  }
  return @out;
}

sub hsps_range {
  # return union of hsp range() for specified type
  my ($hsps, $type) = @_;

  my (@starts, @ends);
#  printf STDERR "ranges for $type:\n";
  foreach my $hsp (@{$hsps}) {
    my ($start, $end) = $hsp->range($type);
#    printf STDERR "  %d-%d\n", $start, $end;
    push @starts, $start;
    push @ends, $end;
  }
  return(min(@starts), max(@ends));
}

sub hsps_strand {
  # return union of hsp range() for specified type
  my ($hsps) = @_;
#  my %strands = map {$_->strand(), 1} @{$hsps};
  # FAIL: strand() without an argument refers to the query sequence
  # by default, which at least with BLAST always seems to be +,
  # making this operation a NOP.
  my %strands = map {$_->strand("hit"), 1} @{$hsps};

  my @strands = keys %strands;
  die "multi strand" unless @strands == 1;
  return $strands[0];
}

sub detect_large_itd {
  my %options = @_;
  my $hsps = $options{"-hsp"} || die "-hsp";
  my $contig = $options{"-contig"} || die "-contig";
  my $gene = $options{"-gene"} || die "-gene";
  my $notes = $options{"-notes"} || die "-notes";
  my $new_contig;

  if (@{$hsps} > 1) {
    my @hsp_by_hit = sort {$a->start("hit") <=>
			     $b->start("hit")} @{$hsps};
    # TO DO: if > 2 HSPs, filter to just the out of order ones?
    # might fix some problems without having to rewrite below code
    # only prune if some entries are ordered and some are not
    #
    # /home/medmonso/work/steve/2020_06_13_cicero_contig_extension/2021_02_03_clingen_all/coding_only/crash2/test.tab

    if ($VERBOSE) {
      printf STDERR "HSP ranges before:\n";
      foreach my $hsp (@hsp_by_hit) {
	printf STDERR "  q:%d-%d hit:%d-%d\n",
	  $hsp->range("query"), $hsp->range("hit");
      }
    }

    if (@hsp_by_hit > 2) {
      # attempt pruning: if some HSPs are ordered but some are not,
      # just keep the out-of-order ones
      push @{$notes}, "large_itd_check_more_than_2_hsp";
      my @prune;
      my $end = scalar(@hsp_by_hit) - 1;

      for (my $i = 0; $i < $end; $i++) {
	if ($hsp_by_hit[$i]->start("query") >
	    $hsp_by_hit[$i + 1]->start("query")) {
	  # blocks are out of order, keep only these two as they
	  # likely capture the critical boundary
	  @hsp_by_hit = ($hsp_by_hit[$i], $hsp_by_hit[$i + 1]);
	  last;
	}
      }
    }

    if ($VERBOSE) {
      printf STDERR "HSP ranges after:\n";
      foreach my $hsp (@hsp_by_hit) {
	printf STDERR "  q:%d-%d hit:%d-%d\n",
	  $hsp->range("query"), $hsp->range("hit");
      }
    }

    if (@hsp_by_hit > 2) {
      push @{$notes}, "cannot_evaluate_large_ITD_more_than_2_hsps";
      return "";
    }


    if ($hsp_by_hit[0]->start("query") > $hsp_by_hit[$#hsp_by_hit]->start("query")) {
      # matching blocks in gene are out of order vs. the contig.
      # this can be evidence of an ITD that's larger than the contig:
      # - first part of the contig matches the END of a duplicated block
      # - second part of the contig matches the START of the same duplicated
      #   block, which is EARLIER in the gene

      # - for each HSP, build TRIMMED start/end intervals for query/hit

      my (@qs, @qe, @hs, @he);

      my $last_rh;
      foreach my $hsp (@hsp_by_hit) {
	my $rq = get_trimmed_range($hsp, "query");
	my $rh = get_trimmed_range($hsp, "hit");
	printf STDERR "trimmed q=%d-%d h=%d-%d\n", @{$rq}, @{$rh};
	if ($last_rh and $last_rh->[1] >= $rh->[0]) {
	  # this HSP's match to the hit overlaps with the previous HSP.
	  # nudge start of this HSP match ahead to resolve.
	  # If we don't do this we'll double-count one or more bases,
	  # breaking ITD calculation.
	  # example: /home/medmonso/work/steve/2020_06_13_cicero_contig_extension/2021_01_04_all/blast/big_ITD_FGFR1.tab, NM_001174063
	  while ($last_rh->[1] >= $rh->[0]) {
	    $rh->[0]++;
	    $rq->[0]++;
	    # hack
	  }
	  printf STDERR "  adjusted: q=%d-%d h=%d-%d\n", @{$rq}, @{$rh};
	}

	push @qs, $rq->[0];
	push @qe, $rq->[1];
	push @hs, $rh->[0];
	push @he, $rh->[1];

	$last_rh = $rh;

      }
      # FIX ME: detect/resolve OVERLAP first!!


#      my @hs = map { $_->start("hit") } @{$hsps};
#      my @he = map { $_->end("hit") } @{$hsps};
      my $hs = min(@hs);
      my $he = max(@he);

      my $gap_start = $hs - 1;
      my $gap_end = length($contig) - $he;

      # QC: check that mapping to gene covers both start and end.
      # However, possible some non-matching sequence refers to a
      # different isoform??
      my $known_mismatch;
      if ($gap_start or $gap_end) {
	push @{$notes}, sprintf "large_itd_gap_at_contig_end=%d", $gap_end if $gap_end;
	push @{$notes}, sprintf "large_itd_gap_at_contig_start=%d", $gap_start if $gap_start;
#	die sprintf "hit range doesn't end at last base, len=%d he=%d gap=%d", length($contig), $he, $gap if $gap;
	$known_mismatch = 1;
      }

      my $qs = min(@qs);
      my $qe = max(@qe);

#      die "$qs $qe $hs $he";

#      my $end_first = $hsp_by_hit[0]->end("hit");
      my $end_first = $he[0];
#      die join ",", $end_first, $he[0];

#      my $start_second = $hsp_by_hit[1]->start("hit");
      my $start_second = $hs[1];

      my ($interstitial_start, $interstitial_end);
      if (($start_second - $end_first) > 1) {
	# check if a portion of the contig doesn't match the putative repeat
	# (might or might not happen)
	$interstitial_start = $end_first + 1;
	$interstitial_end = $start_second - 1;
      }

      my $dup_len = ($qe - $qs) + 1;
      my $dup_seq = substr($gene, $qs - 1, $dup_len);

#      printf STDERR "gene=%s\n", $gene;
#      printf STDERR "dup=%s\n", $dup_seq;

      my $upstream = get_flanking("-sequence" => \$gene,
				  "-direction" => -1,
				  "-edge" => $qs);
#      printf STDERR "upstream=%s\n", $upstream;
      my $downstream = get_flanking("-sequence" => \$gene,
				  "-direction" => 1,
				  "-edge" => $qe);
#      printf STDERR "downstream=%s\n", $downstream;
#      printf STDERR "joined=%s/%s\n", $dup_seq, $downstream;

      $new_contig = $upstream . $dup_seq;

      if ($interstitial_start) {
	$new_contig .= bases2chunk(\$contig, $interstitial_start, $interstitial_end);
      }
      $new_contig .= $dup_seq . $downstream;
      # second copy of duplicated sequence plus downstream flanking

      unless (index($new_contig, $contig) > 0) {
	my $msg = "large_itd_fail_to_find_old_contig_in_new_contig";
	push @{$notes}, $msg;
	# not all matches are perfect
      }
    }
  }

  printf STDERR "large ITD detected: %s\n", $new_contig ? "yup" : "nope";

  return $new_contig;
}

sub get_fai {
  my $f_fasta = config_or_manual(
				 "-config-type" => "genome",
				 "-config-name" => $FLAGS{genome},
				 "-parameter" => "FASTA",
				 "-manual" => $FLAGS{fasta}
				);

  printf STDERR "FASTA: %s\n", $f_fasta if $VERBOSE;
  return new FAI("-fasta" => $f_fasta,
		 "-chunk_cache_chromosomes" => 1,
		);
}

sub get_flanking {
  my (%options) = @_;
  my $seq_ref = $options{"-sequence"} || die;
  my $dir = $options{"-direction"} || die;
  my $edge = $options{"-edge"} || die;

  my $seq_len = length $$seq_ref;
  my ($start_base, $end_base);
  if ($dir == -1) {
    $end_base = $edge - 1;
    $start_base = ($end_base - $LARGE_ITD_ADD_FLANKING) + 1;
  } elsif ($dir == 1) {
    $start_base = $edge + 1;
    $end_base = $start_base + ($LARGE_ITD_ADD_FLANKING - 1);
  } else {
    die;
  }
  $start_base = 1 if $start_base < 1;
  $end_base = $seq_len if $end_base > $seq_len;

  return bases2chunk($seq_ref, $start_base, $end_base);
}

sub bases2chunk {
  my ($seq_ref, $start_base, $end_base) = @_;
  my $chunk_len = ($end_base - $start_base) + 1;
  return substr($$seq_ref, $start_base - 1, $chunk_len);

}

sub blast2hsps {
  my (%options) = @_;
  my $blast = $options{"-blast"} || die "-blast";
  my $parser = $blast->blast(%options);
  my $result = $parser->next_result();
  # one result per query sequence (there's only one query sequence)

  my @hsp;
  if ($result) {
    my $hit_count = 0;
    while (my $hit = $result->next_hit()) {
      # hits from this query to a database sequence
      # (Bio::Search::Hit::HitI)
      die "multiple hits not expected" if ++$hit_count > 1;
      # shouldn't happen

      while (my $hsp = $hit->next_hsp) {
	# Bio::Search::HSP::PSLHSP
	printf STDERR "    strand:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s query_span:%d ref_span:%d total_span=%d query_string=%s hit_string=%s\n",
	  $hsp->strand,
	    $hsp->range("query"),
	      $hsp->range("hit"),
		$hsp->num_identical(),
		  $hsp->frac_identical("query"),
		    $hsp->length("query"),
		      $hsp->length("hit"),
			$hsp->length("total"),
			  $hsp->query_string(),
			    $hsp->hit_string();
	push @hsp, $hsp;
      }
    }
  }

  return @hsp;
  # normally I think returning arrays rather than refs is icky,
  # but these are small lists and it helps w/retrofitting
}

sub collate_cicero_calls {
  # gather fusion calls from CICERO reports (which may have slightly
  # different header layouts) then (A) filter to usable events only
  # and (B) write out a subset of columns.
  my $files = read_simple_file($FLAGS{"collate-cicero"} || die);
  my $f_out = $FLAGS{out} || "cicero_collated.tab";

  #
  #  TO DO: prune duplicates based on geneA/geneB/contig sequence??
  #

  my $f_ort = "sv_ort";
  my $f_gene_a = "geneA";
  my $f_gene_b = "geneB";
  my $f_feature_a = "featureA";
  my $f_feature_b = "featureB";
  my $f_contig = "contig";

  my $suppress_duplicates = 1;

  printf STDERR "suppress duplicates by chr/pos/contig?: %s\n", $suppress_duplicates ? "yes" : "no";
  my @dup_key_fields = (
			"ChrA",
			"PosA",
			"ChrB",
			"PosB",
			$f_contig
		       );

  my @f_copy = (
		"Sample(s)",
		"ChrA",
		"PosA",
		$f_gene_a,
		"ChrB",
		"PosB",
		$f_gene_b,
		$f_feature_a,
		$f_feature_b,

		$f_contig
	       );

  my $rpt = new Reporter(
			 "-file" => $f_out,
			 "-delimiter" => "\t",
			 "-labels" => \@f_copy,
			);
  my $rpt_coding = new Reporter(
			 "-file" => $f_out . ".coding",
			 "-delimiter" => "\t",
			 "-labels" => \@f_copy,
			);

  my %filtered;

  my %features_usable = (
			 "intergenic" => 0,
			 "coding" => 1,

			 "intron" => 1,
			 "5utr" => 1,
			 "3utr" => 1,
			 # these might or might not be usable depending
			 # on proximity to an exon
			);

  my %saw;

  my $c = new Counter($files, "-mod" => 25);
  my $usable = 0;
  foreach my $infile (@{$files}) {
    my $df = new DelimitedFile("-file" => $infile,
			       "-headers" => 1,
			      );

    my $headers = $df->headers();
    foreach my $f (@f_copy) {
      dump_die($headers, "can't find field $f in $infile") unless defined $headers->{$f};
    }

    # while (my $row = $df->next("-ref" => 1)) {  # headerless
    while (my $row = $df->get_hash()) {
      # filter to usable events only and reocrd exceptions
      my $ort = $row->{$f_ort} || die;
      my $reject;
      my $both_coding;
      if ($ort eq ">") {
	my $feat_a = $row->{$f_feature_a} || die;
	my $feat_b = $row->{$f_feature_b} || die;

	$both_coding = 1;
	foreach my $f ($f_feature_a, $f_feature_b) {
	  my $feat = $row->{$f} || die;
	  $both_coding = 0 unless $feat eq "coding";
	  my $usable = $features_usable{$feat};
	  dump_die($row, "undefined behavior for $f $feat") unless defined $usable;
	  $reject = "unusable_feature" unless $usable;
	}
      } else {
	$reject = "unusable_ort";
      }

      if ($suppress_duplicates) {
	my $key = join ".", @{$row}{@dup_key_fields};
	$reject = "duplicate" if $saw{$key};
	$saw{$key} = 1;
      }

      if ($reject) {
	$filtered{$reject}++;
      } else {
	$rpt->end_row($row);
	$rpt_coding->end_row($row) if $both_coding;
	$usable++;
      }
    }
    $c->next();
  }
  $rpt->finish();
  $rpt_coding->finish();

  my $rejected = sum values %filtered;
  printf STDERR "usable:%d rejected:%d\n", $usable, $rejected;
  foreach my $k (sort keys %filtered) {
    printf STDERR "  %s: %d\n", $k, $filtered{$k};
  }
}

sub diagnose_pattern_matches {
  my $f_report = $FLAGS{"diagnose-matches"} || die;
  my $f_patterns = $FLAGS{"patterns"} || die "-patterns";

  open(IN, $f_patterns) || die;
  my %patterns;
  while (<IN>) {
    chomp;
    my @f = split /\t/, $_;
    die unless @f == 2;
    my ($id, $pattern) = @f;
    die if $patterns{$id};
    $patterns{$id} = $pattern;
  }

  open(IN, $f_report) || die;
  my %read_count;
  my $pattern_id;
  my %pattern2seq;
  while (<IN>) {
    chomp;
    my @f = split /\t/, $_;
    die unless @f == 2;
    my ($id, $sequence) = @f;
    next if $sequence eq "read pairs";
    # not sure what this is
    $sequence =~ s/^\s+//;
    die $sequence if $sequence =~ /\s/;
    if ($id =~ /^pattern (\S+)/) {
      # new pattern reference
      $pattern_id = $1;
      %read_count = ();
    } else {
      # matching read
      my $fasta_id = join "_", $id, ++$read_count{$id};
      die unless $pattern_id;
      push @{$pattern2seq{$pattern_id}{$id}}, [ $fasta_id, $sequence ];
    }
  }

  # for each observed pattern, write FASTA files for pattern and matching reads
  foreach my $pid (sort keys %pattern2seq) {
    #
    #  write FASTA file for pattern, indicating bracketed portion:
    #
    my $pattern_seq = $patterns{$pid} || die "can't find pattern for $pid";
    printf STDERR "pattern for %s is %s\n", $pid, $pattern_seq;

    my $i_a_end = get_edge_idx($pattern_seq, "A");
    my $i_b_start = get_edge_idx($pattern_seq, "B");

    my ($dmz_start, $dmz_end);
    # 1-based base numbers of sequence in bracketed/DMZ area

    if ($i_b_start == $i_a_end + 1) {
      die "flush $pattern_seq";
    } else {
      # gap/DMZ present
      $dmz_start = $i_a_end + 1;
      $dmz_end = $i_b_start - 1;
      # NEEDS TESTING
    }

    print STDERR "pattern $pattern_seq idx=$i_a_end $i_b_start  dmz1 = $dmz_start $dmz_end\n";

    my $pattern_nt = $pattern_seq;
    $pattern_nt =~ s/[\[\]\{\}]//g;
    die $pattern_nt if $pattern_nt !~ /[ACGT]/i;
#    die join "\n", $pattern_seq, $pattern_nt;

    my $f_pattern = sprintf "%s_pattern.fa", $pid;
    # put pattern ID first for autocomplete purposes

    open(OUT, ">" . $f_pattern) || die;
    printf OUT ">%s", $pid;
    if (defined $dmz_start) {
      printf OUT " bracketed=%d-%d", $dmz_start, $dmz_end;
    }
    print OUT "\n";
    printf OUT "%s\n", $pattern_nt;
    close OUT;

    #
    #  write FASTA file for each pair of supporting reads:
    #
    my $pair_num = 0;
    foreach my $id (sort keys %{$pattern2seq{$pid}}) {
      my $f_seqs = sprintf "%s_seqs_%d.fa", $pid, ++$pair_num;
      open(OUT, ">" . $f_seqs) || die;
      my $seqs = $pattern2seq{$pid}{$id} || die;
      foreach my $ref (@{$seqs}) {
	my ($id, $seq) = @{$ref};
	printf OUT ">%s\n%s\n", $id, $seq;
      }
      close OUT;
    }

  }

}

sub get_edge_idx {
  my ($pattern_seq, $end) = @_;
  my @search;
  if ($end eq "A") {
    @search = ( "]", "}" );
  } elsif ($end eq "B") {
    @search = ( "[", "{" );
  } else {
    die;
  }

  my $result;
  foreach my $string (@search) {
    my $i = index($pattern_seq, $string);
    if ($i >= 0) {
      $result = $i;
      last;
    }
  }
  die "can't find pattern edge" unless defined $result and $result >= 0;
  return $result;
}

sub generate_genomic_contigs {
  # from e.g. CREST WGS fusion calls
  if ($FLAGS{crest}) {
    $F_CHR_A = 'chrA';
    # in raw files, header line is commented, so may be "#chrA"
    $F_POS_A = 'posA';
    $F_ORT_A = 'ortA';
    $F_CHR_B = 'chrB';
    $F_POS_B = 'posB';
    $F_ORT_B = 'ortB';
  } else {
    die "need WGS fusion specification header setup, only -crest implemented";
  }

  my $infile = $FLAGS{file} || die "-file";
  my $genome = $FLAGS{genome} || die "-genome";

  my $df = new DelimitedFile(
			     "-file" => $infile,
			     "-headers" => 1,
			    );

  field_name_dwim(
		  "-query" => [
			       \$F_CHR_A,
			       \$F_POS_A,
			       \$F_ORT_A,
			       \$F_CHR_B,
			       \$F_POS_B,
			       \$F_ORT_B
			      ],
		  "-actual" => $df->headers_raw(),
		 );

  my $outfile = basename($infile) . ".extended.tab";

  my @extra;

  my %h = map {$_, 1} @{$df->headers_raw};

  my $f_gene_a = "genea_symbol";
  my $f_gene_b = "geneb_symbol";

  push @extra, $f_gene_a unless $h{$f_gene_a};
  push @extra, $f_gene_b unless $h{$f_gene_b};
  push @extra, $F_SAMPLE unless $h{$F_SAMPLE};
  # needed for fuzzion2 pattern generation

  my @f_null_clone = qw(
		genea_acc
		genea_feature

		geneb_acc
		geneb_feature

		gene_pair_summary
		gene_pair_summary_extended
		      );

  # TO DO: modify the below to use @NEW_HEADERS
  push @extra, (
		"contig",
		# generate similated cicero contig (small)

		"genea_contig_start",
		"genea_contig_end",
		"geneb_contig_start",
		"geneb_contig_end",

		"extended_exon",
		# this doesn't apply in this mode, clone of the next field
		$F_EXTENDED_CHUNK,
		# preferred extended flanking size
		"extended_max",
		# since this mode is not based on a gene model, just
		# hardcode a larger number
		$F_GENOME,
		$F_PROCESSING_TYPE,
		@f_null_clone,
		$F_NOTES
	       );

  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => \@extra,
  			      "-auto_qc" => 1,
			     );

  my $fai = get_fai();

  my $max_len_i = $WGS_CONTIG_MAX_LEN - 1;

  my $f_refflat = get_refflat(1);
  my $ga = new GeneAnnotation(
			      "-style" => "refgene_flatfile",
			      "-refgene_flatfile" => $f_refflat,
			      "-ignore_non_coding" => 0
			     );

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
    my @notes;
    my $chr_a = $row->{$F_CHR_A} || dump_die($row, "no $F_CHR_A");
    my $pos_a = $row->{$F_POS_A} || die;
    my $ort_a = $row->{$F_ORT_A} || die;

    my $chr_b = $row->{$F_CHR_B} || die;
    my $pos_b = $row->{$F_POS_B} || die;
    my $ort_b = $row->{$F_ORT_B} || die;

    #
    #  annotate genes, unless already present:
    #
    annotate_genes("-chr" => $chr_a,
		   "-pos" => $pos_a,
		   "-field" => $f_gene_a,
		   "-ga" => $ga,
		   "-row" => $row
		  );

    annotate_genes("-chr" => $chr_b,
		   "-pos" => $pos_b,
		   "-field" => $f_gene_b,
		   "-ga" => $ga,
		   "-row" => $row
		  );

    if ($FLAGS{"annotate-genes-as-breakpoints"}) {
      # for source data that consists of breakpoints only,
      # may be easier to review if the gene names are set to
      # genomic breakpoint details.
      $row->{$f_gene_a} = join "_", $chr_a, $pos_a, $ort_a;
      $row->{$f_gene_b} = join "_", $chr_b, $pos_b, $ort_b;
      # may need work as strand may interfere with pattern separator.
      # maybe spell out "plus" and "minus"?
    }

    #
    # get flanking sequences:
    #
    my $upstream_max;
    my $downstream_max;
    if ($ort_a eq "+") {
      $upstream_max = $fai->get_chunk(
				      "-id" => $chr_a,
				      "-start" => $pos_a - $max_len_i,
				      "-end" => $pos_a
				     );
    } elsif ($ort_a eq "-") {
      $upstream_max = $fai->get_chunk(
				      "-id" => $chr_a,
				      "-start" => $pos_a,
				      "-end" => $pos_a + $max_len_i
				     );
#      print STDERR "upstream A genomic:$upstream_max  rc:";
      $upstream_max = reverse_complement($upstream_max);
#      print STDERR "$upstream_max\n";
    } else {
      die "invalid A ort";
    }

    if ($ort_b eq "+") {
      $downstream_max = $fai->get_chunk(
				      "-id" => $chr_b,
				      "-start" => $pos_b,
				      "-end" => $pos_b + $max_len_i
				     );
#      die "$chr_b $pos_b $downstream_max";
    } elsif ($ort_b eq "-") {
      $downstream_max = $fai->get_chunk(
				      "-id" => $chr_b,
				      "-start" => $pos_b - $max_len_i,
				      "-end" => $pos_b
				     );
      $downstream_max = reverse_complement($downstream_max);
    } else {
      die "invalid B ort";
    }

    dump_die($row, "TEST ME: both strands are -")
      if $ort_a eq "-" and $ort_b eq "-";

    $upstream_max = lc($upstream_max);
    $upstream_max =~ s/(.{$WGS_CONTIG_BASE_FLANKING})$/uc($1)/e;

    $downstream_max = lc($downstream_max);
    $downstream_max =~ s/^(.{$WGS_CONTIG_BASE_FLANKING})/uc($1)/e;

    $row->{genea_contig_start} = 1;
    $row->{genea_contig_end} = $WGS_CONTIG_BASE_FLANKING;
    $row->{geneb_contig_start} = $WGS_CONTIG_BASE_FLANKING + 1;
    $row->{geneb_contig_end} = $WGS_CONTIG_BASE_FLANKING * 2;

    my $contig_base = join "",
      substr($upstream_max, - $WGS_CONTIG_BASE_FLANKING),
	substr($downstream_max, 0, $WGS_CONTIG_BASE_FLANKING);

    my $contig_max = join "", $upstream_max, $downstream_max;
    my $contig_extended = join "",
      substr($upstream_max, - $EXTENDED_CHUNK_LENGTH),
	substr($downstream_max, 0, $EXTENDED_CHUNK_LENGTH);

    $row->{contig} = $contig_base;
    $row->{$F_EXTENDED_CHUNK} = $contig_extended;
    $row->{extended_exon} = $contig_extended;
    # hack
    $row->{extended_max} = $contig_max;

#    printf STDERR "as:%d ae:%d bs:%d be:%d contig:%s up:%s down:%s\n", @{$row}{qw(genea_contig_start genea_contig_end geneb_contig_start geneb_contig_end contig)}, $upstream_max, $downstream_max;

    if ($chr_a eq $chr_b) {
      # for events on the same chrom, report distance between
      # breakpoints.  nearby events may imply problematic situations
      # like ITDs or very small deletions.
      my $distance = abs($pos_b - $pos_a);
      push @notes, sprintf "distance=%d", $distance;
    }

    push @notes, TAG_NO_CONTIG;

    $row->{$F_NOTES} = join ",", @notes;
    $row->{$F_SAMPLE} = "" unless $row->{$F_SAMPLE};
    $row->{$F_GENOME} = standardize_genome($genome);
    $row->{$F_PROCESSING_TYPE} = "genomic";

    foreach my $f (@f_null_clone) {
      $row->{$f} = "";
    }

    $rpt->end_row($row);

    if (0) {
      # for BLAST debugging
      print STDERR "debug\n";

      open(TMP, ">contig_orig.txt") || die;
      printf TMP ">contig_orig\n%s\n", $row->{contig_crest} || die;

      open(TMP, ">contig_new.txt") || die;
#      printf TMP ">contig_new\n%s\n", $contig_extended;
      printf TMP ">contig_new\n%s\n", $contig_max;

      die "X";
    }


  }

  $rpt->finish();
}

sub annotate_genes {
  my (%options) = @_;
  my $chr = $options{"-chr"} || die;
  my $pos = $options{"-pos"} || die;
  my $field = $options{"-field"} || die;
  my $ga = $options{"-ga"} || die;
  my $row = $options{"-row"} || die;

  unless ($row->{$field}) {
    # respect any existing gene annotation in this field
    $ga->find(
	      "-reference" => $chr,
	      "-start" => $pos,
	      "-end" => $pos
	     );
    my $genes_genomic = $ga->results_genes_genomic_order();
    my $v;
    if (@{$genes_genomic}) {
      $v = join ",", @{$genes_genomic};
    } else {
      $v = join ".", $chr, $pos;
      # hack
    }
    $row->{$field} = $v;
  }
}

sub generate_rna_contigs {
  #
  # from e.g. ProteinPaint RNA coordinates
  #
  my $infile = $FLAGS{file} || die "-file";
  init_sjpi();
  my $df = new DelimitedFile(
			     "-file" => $infile,
			     "-headers" => 1,
			    );
  my $outfile = basename($infile) . ".extended.tab";

#  my $keep_gene_annotation = $FLAGS{"keep-gene-annotation"};

  printf STDERR "  intragenic events suppressed: %s\n", $RNA_SUPPRESS_INTRAGENIC ? "y": "n";
  unless ($RNA_SUPPRESS_INTRAGENIC) {
    # only applies when events not suppressed completely
    printf STDERR "  for intragenic events, minimum breakpoint distance required to process: %d\n", $RNA_INTRAGENIC_MIN_DISTANCE_TO_PROCESS;
  }
  printf STDERR "  add/update gene symbol/strand: %s\n", $RNA_ADD_GENE_ANNOTATIONS ? "y" : "n";
  printf STDERR "    preserve existing info: %s\n", $RNA_ADD_GENE_ANNOTATIONS ? "y" : "n";

  field_preset_setup();

  die "chr/pos/ort fields not defined; suggest using -cicero, -crest, etc." unless ($F_CHR_A and $F_POS_A and $F_CHR_B and $F_POS_B);

  my $genome = standardize_genome($FLAGS{genome} || die "-genome");

  my @dwim_check = (
		    \$F_CHR_A,
		    \$F_POS_A,

		    \$F_CHR_B,
		    \$F_POS_B,
		   );

  push @dwim_check, (\$F_ORT_B, \$F_ORT_A) unless $FLAGS{"rna-no-ort"};
  # for orientation fields, even if the field is not mandatory we want
  # it to be specified because some sourcess are spotty,
  # e.g. ProteinPaint export mostly has the annotation, but not
  # always.  Need to be able to distinguish between a column that's
  # not specified properly (e.g. a typo in specification) and one that
  # sometimes simply doesn't contain data.

  if ($F_GENE_A and $F_GENE_B) {
    push @dwim_check, \$F_GENE_A, \$F_GENE_B;
  } else {
    if ($RNA_KEEP_GENE_ANNOTATIONS) {
      die "gene fields not defined";
    } else {
      # will be generating gene annotations, use placeholder field names
      $F_GENE_A = "__geneA";
      $F_GENE_B = "__geneB";
    }
  }

  field_name_dwim(
		  "-query" => \@dwim_check,
		  "-actual" => $df->headers_raw(),
		 );

  my $blacklist = get_blacklist();

  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => [
					   "contig",
					   @NEW_HEADERS
					  ],
  			      "-auto_qc" => 1,
			     );

  my $rf = get_refflat();

  my $gsm = init_gsm($rf);

#  my $need_gene_annotation = $keep_gene_annotation ? $RNA_ADD_GENE_ANNOTATIONS : 1;
  my $ga;
  if ($RNA_ADD_GENE_ANNOTATIONS) {
      my $f_refflat = get_refflat(1);
      $ga = new GeneAnnotation(
			       "-style" => "refgene_flatfile",
			       # will this work w/ncbiRefSeq?
			       "-refgene_flatfile" => $f_refflat,
			       "-ignore_non_coding" => 0
			      );
      # TO DO: share raw row data with $rf above?
  }


  #
  #  init reference sequence:
  #
  my $fai = get_fai();
  my $tf = new TranscriptFlanking(
				  "-fai" => $fai
				 );

  my $input_row_number = 0;
  my %suppressed;
  while (my $row = $df->get_hash()) {
    dump_die($row, "new row", 1) if $VERBOSE;

    $input_row_number++;
    my $input_row_md5 = md5_hex(map {$row->{$_}} @{$df->headers_raw});

    my @notes_general;

    add_gene_annotation($row, $ga, $gsm, \@notes_general)
      if $RNA_ADD_GENE_ANNOTATIONS;

    my $chr_a = $row->{$F_CHR_A} || die;
    my $chr_b = $row->{$F_CHR_B} || die;
    my $pos_a_raw = $row->{$F_POS_A} || die;
    my $pos_b_raw = $row->{$F_POS_B} || die;

    my $gene_a = $row->{$F_GENE_A};
    my $gene_b = $row->{$F_GENE_B};
    my $is_intragenic_simple = $gene_a eq $gene_b;
    # many intragenic records have a single gene symbol on each end.
    # however there are some compound cases the following which will
    # have to be caught later below:
    #   Project: pcgp
    #   chrA: chr8
    #   posA: 2821204
    #   geneA: LOC105377785,CSMD1

    my $ort_a = $row->{$F_ORT_A};
    my $ort_b = $row->{$F_ORT_B};
    if (not($FLAGS{"rna-no-ort"}) and $RNA_ORT_REQUIRED) {
      dump_die($row, "no $F_ORT_A") unless $ort_a;
      dump_die($row, "no $F_ORT_B") unless $ort_b;
    }

    my $can_process = 1;

    my %transcript2gene;

    #
    #  initialize transcripts for gene A (required):
    #
    my ($gene_a_transcript2rf, $gene_b_transcript2rf);

    my $intergenic_b;

    if ($gene_a) {
      $gene_a_transcript2rf = get_transcripts($rf, $gsm, $gene_a, \@notes_general, \%transcript2gene, "-chr" => $chr_a);
      transcript_filter($gene_a_transcript2rf, $gene_a, $FLAGS{"restrict-nm-a"});
      unless (keys %{$gene_a_transcript2rf}) {
	$can_process = 0;
	log_error("missing_transcripts_A", \%suppressed, \@notes_general);
      }
    }

    if ($gene_a and $gene_b) {
      $gene_b_transcript2rf = get_transcripts($rf, $gsm, $gene_b, \@notes_general, \%transcript2gene, "-chr" => $chr_b);
      # mapping of transcripts to refFlat record entries

      transcript_filter($gene_b_transcript2rf, $gene_b, $FLAGS{"restrict-nm-b"});

      unless (keys %{$gene_b_transcript2rf}) {
	$can_process = 0;
	log_error("missing_transcripts_B", \%suppressed, \@notes_general);
	# potentially double-counting
      }

      if ($blacklist and is_blacklisted($blacklist, $row)) {
	$can_process = 0;
	log_error("blacklisted", \%suppressed, \@notes_general);
      }
    } elsif ($gene_a and not($gene_b)) {
      # geneB site is intergenic
      $intergenic_b = 1;
    } else {
      $can_process = 0;

      if (not($gene_a) and not($gene_b)) {
	log_error(TAG_MISSING_GENE_ANNOTATION . "=both", \%suppressed, \@notes_general);
      } elsif (not($gene_a)) {
	log_error(TAG_MISSING_GENE_ANNOTATION . "=A", \%suppressed, \@notes_general);
      } elsif (not($gene_b)) {
	log_error(TAG_MISSING_GENE_ANNOTATION . "=B", \%suppressed, \@notes_general);
      } else {
	die "unpossible";
      }
      # TO DO: possible to salvage via fresh gene annotation?
    }

    if ($can_process and $is_intragenic_simple) {
      if ($RNA_SUPPRESS_INTRAGENIC) {
	# not allowed at all, safest route with the limited information
	# we have available
	$can_process = 0;
	log_error("intragenic_event_suppressed", \%suppressed, \@notes_general);
      } else {
	# allowed, but with a minimum distance requirement (hacky)
	my $bp_distance = abs($pos_a_raw - $pos_b_raw);
	if ($bp_distance >= $RNA_INTRAGENIC_MIN_DISTANCE_TO_PROCESS) {
	  # OK: assume a deletion, e.g. exon skip
	} else {
	  $can_process = 0;
	  log_error("intragenic_event_suppressed_distance", \%suppressed, \@notes_general);
	}
      }
    }

    unless ($can_process) {
      # processing not possible, skip
      my %r = %{$row};
      $r{contig} = "";
      $r{input_row_number} = $input_row_number;
      $r{$F_NOTES} = join ",", @notes_general;
      foreach my $f (@NEW_HEADERS) {
	$r{$f} = "" unless defined $r{$f};
      }
      $rpt->end_row(\%r);
      next;
    }

    #
    #  preprocess A accessions to remove unprocessable records:
    #
    my @acc_a_usable;
    foreach my $acc_a (sort keys %{$gene_a_transcript2rf}) {
      my $rf_a = $gene_a_transcript2rf->{$acc_a} || die;

      my @problem;
      unless (strand_consistency_check($rf_a, $ort_a)) {
	# geneA strand is expected to be consistent with gene model
	log_error("gene_A_transcript_antisense", \%suppressed, \@problem);
      }

      my ($a_up, $pos_a) = $tf->get_upstream_downstream(
						     "-row" => $rf_a,
						     "-direction" => "upstream",
						     "-pos" => $pos_a_raw,
						     "-extended" => 1
						    );
      # upstream sequence for geneA from start of gene to this base
      # be sure to use the RAW value every time as this iterates
      unless ($a_up) {
	# don't mark as an error, as this can happen
	push @problem, "A_position_outside_of_transcript";
      }

      if (@problem) {
	# some genomic positions may be outside of the gene model for
	# some isoforms
	#
	# Sample  Project ChrA    PosA    geneA   ChrB    PosB    geneB   svtype
	# SJOS040208_D1-PAPXGT    pantarget       chr1    1130200 TTLL10  chr1    2545095MMEL1

	# here for TTLL10, the position is downstream of 3' UTR for NM_153254

	my %r = %{$row};
	$r{input_row_number} = $input_row_number;
	$r{genea_acc} = $acc_a;
	$r{genea_symbol} = $transcript2gene{$acc_a} || die;
	$r{$F_NOTES} = join ",", @notes_general, @problem;
	foreach my $f (@NEW_HEADERS, "contig") {
	  $r{$f} = "" unless defined $r{$f};
	}
	$rpt->end_row(\%r);
	next;
      } else {
	push @acc_a_usable, $acc_a;
      }
    }

    if ($RNA_SUPPRESS_DISTANT_ISOFORMS) {
      #
      # optionally filter out distant isoforms for position A:
      #
      my @close;
      my @distant;

      foreach my $acc_a (@acc_a_usable) {
	my $rf_a = $gene_a_transcript2rf->{$acc_a} || die;

	my @problem;

	my ($a_up, $pos_a) = $tf->get_upstream_downstream(
							  "-row" => $rf_a,
							  "-direction" => "upstream",
							  "-pos" => $pos_a_raw,
							  "-extended" => 1
							 );
	die unless $a_up;
	# should have been detected previously

	my $dist = abs($pos_a_raw - $pos_a);
	if ($dist >= $RNA_DISTANT_THRESHOLD) {
	  push @distant, $acc_a;
	} else {
	  push @close, $acc_a;
	}
      }

      if (@close and @distant) {
	printf STDERR "%d: suppress distant %s in favor of close %s\n",
	  $pos_a_raw, join(",", @distant), join(",", @close);
	@acc_a_usable = @close;
      }
    }

    foreach my $acc_a (@acc_a_usable) {
      my $rf_a = $gene_a_transcript2rf->{$acc_a} || die;

      my @problem;

      my ($a_up, $pos_a) = $tf->get_upstream_downstream(
						     "-row" => $rf_a,
						     "-direction" => "upstream",
						     "-pos" => $pos_a_raw,
						     "-extended" => 1
						    );
      die unless $a_up;
      # should have been detected previously

      my $cd_a = set_coding_distance($row, "a", $pos_a_raw, $pos_a);

      my @iterate;
      if ($intergenic_b) {
	# gene B is intergenic: use placeholder for accession
	@iterate = "__intergenic__";
      } else {
	# gene B is in a gene model
	@iterate = sort keys %{$gene_b_transcript2rf};
      }

      my @rows_b_ok;

      foreach my $acc_b (@iterate) {
	my @problem;
	my @notes_b;

	printf STDERR "B ref: %s\n", $intergenic_b ? "intergenic" : $acc_b;

	my ($rf_b, $b_down, $pos_b);

	if ($intergenic_b) {
	  #
	  #  B site is intergenic, use genomic sequence.
	  #
	  $pos_b = $pos_b_raw;

	  $b_down = $fai->get_flanking(
				       "-direction" => "downstream",
				       "-chr" => $chr_b,
				       "-pos" => $pos_b,
				       "-ort" => $ort_b,
				       "-length" => $WGS_CONTIG_MAX_LEN
				      );

	} else {
	  #
	  #  B site in a gene model
	  #
	  my $is_intragenic = ($transcript2gene{$acc_a} || die) eq ($transcript2gene{$acc_b} || die);
	  # now that we have accession pair, check again for intergenic
	  # (can happen if geneA/B are a list)
	  if ($is_intragenic) {
	    next if $acc_a ne $acc_b;
	    # in intragenic mode, work within each single transcript rather
	    # than combinations

	    if ($RNA_SUPPRESS_INTRAGENIC) {
	      # not allowed at all, safest route with the limited information
	      # we have available
	      log_error("intragenic_event_suppressed", \%suppressed, \@problem);
	    } else {
	      # allowed, but with a minimum distance requirement (hacky)
	      my $bp_distance = abs($pos_a_raw - $pos_b_raw);
	      if ($bp_distance >= $RNA_INTRAGENIC_MIN_DISTANCE_TO_PROCESS) {
		# OK: assume a deletion, e.g. exon skip
	      } else {
		log_error("intragenic_event_suppressed_distance", \%suppressed, \@problem);
	      }
	    }
	  }

	  $rf_b = $gene_b_transcript2rf->{$acc_b} || die;

	  if (strand_consistency_check($rf_b, $ort_b)) {
	    #
	    # fusion orientation annotation matches gene model
	    #
	    ($b_down, $pos_b) = $tf->get_upstream_downstream(
							     "-row" => $rf_b,
							     "-direction" => "downstream",
							     "-pos" => $pos_b_raw,
							     "-extended" => 1
							    );
	    # downstream sequence for geneB from this base to end of gene

	    if ($b_down) {
	      my $distance = abs($pos_b_raw - $pos_b);
	      # see "geneb_coding_distance" in output
	      if ($distance > 6) {
		# For deeply intronic sites we should simply use the
		# genomic sequence on the specified strand.
		#
		# UNFINISHED
		dump_die($row, "FIX ME: adjusted distance to exon is $distance, use genomic instead?", 1);
	      }
	    } else {
	      push @problem, "B_position_outside_of_transcript";
	    }
	  } else {
	    #
	    #  fusion orientation does NOT match the gene model.
	    #  Use flanking genomic sequence for the specified strand.
	    #
	    $pos_b = $pos_b_raw;

	    $b_down = $fai->get_flanking(
					 "-direction" => "downstream",
					 "-chr" => $chr_b,
					 "-pos" => $pos_b,
					 "-ort" => $ort_b,
					 "-length" => $WGS_CONTIG_MAX_LEN
					);
	    push @notes_b, "gene_B_transcript_strand_mismatch";
	    # informational only, not an error
	  }
	}

	my $generate_out_row = sub {
	  # use closure to avoid highly repetitive code in exceptions below
	  my %r = %{$row};
	  foreach my $f (@NEW_HEADERS, "contig") {
#	    $r{$f} = "" unless defined $r{$f};
	    $r{$f} = "";
	    # wipe all columns because $row can be written to
	    # multiple times! should really not do that.
	  }

	  $r{input_row_number} = $input_row_number;
	  $r{input_row_md5} = $input_row_md5;
	  $r{genea_acc} = $acc_a;
	  $r{genea_symbol} = $transcript2gene{$acc_a} || die;
	  $r{genea_pos} = $r{genea_pos_adj} = $pos_a;
	  $r{genea_feature} = get_feature_tag($rf, $rf_a, $pos_a);
	  $r{geneb_acc} = $acc_b;
	  $r{geneb_symbol} = $transcript2gene{$acc_b} || die;
	  $r{geneb_pos} = $r{geneb_pos_adj} = $pos_b;
	  $r{geneb_feature} = get_feature_tag($rf, $rf_b, $pos_b);
	  # these annotations are populatable from here down as A/B are set up
	  # TO DO: maybe move this further up to share everywhere?
	  $r{$F_NOTES} = join ",", @notes_general, @problem, @notes_b;
	  $r{genea_chr} = $chr_a;
	  $r{geneb_chr} = $chr_b;

	  return \%r;
	};

	if (@problem) {
	  my $or = &$generate_out_row();
	  $rpt->end_row($or);
	  next;
	}

	my $cd_b = set_coding_distance($row, "b", $pos_b_raw, $pos_b);

	printf STDERR "A up for %s: %s\n", $acc_a, $a_up;
	printf STDERR "B down for %s: %s\n", $acc_b, $b_down;

	my $upstream_max = lc($a_up);
	my $flank_a = $RNA_CONTIG_BASE_FLANKING;
	if ($flank_a >= length($upstream_max)) {
	  $flank_a = length($upstream_max) - 1;
	  # in the extended contigs we need to leave at least 1 nt of
	  # lowercased sequence beyond the uppercased contig sequence
	  # (this is for fuzzion2 pattern generation).
	  push @problem, "insufficient_flanking_sequence_A" unless $flank_a > 0;
	}

	my $downstream_max = lc($b_down);
	my $flank_b = $RNA_CONTIG_BASE_FLANKING;
#	printf STDERR "dlen=%d\n", length($downstream_max);
	if ($flank_b >= length($downstream_max)) {
	  $flank_b = length($downstream_max) - 1;
	  push @problem, "insufficient_flanking_sequence_B" unless $flank_b > 0;
	}

	if (@problem) {
	  my $or = &$generate_out_row();
	  $rpt->end_row($or);
	  next;
	}

#	$upstream_max =~ s/(.{$RNA_CONTIG_BASE_FLANKING})$/uc($1)/e;
	$upstream_max =~ s/(.{$flank_a})$/uc($1)/e;
#	$downstream_max =~ s/^(.{$RNA_CONTIG_BASE_FLANKING})/uc($1)/e;
	$downstream_max =~ s/^(.{$flank_b})/uc($1)/e;
	# TO DO: share code w/WGS version above

	$row->{genea_contig_start} = 1;
	$row->{genea_contig_end} = $flank_a;
	$row->{geneb_contig_start} = $row->{genea_contig_end} + 1;
#	$row->{geneb_contig_end} = $RNA_CONTIG_BASE_FLANKING * 2;
#	$row->{geneb_contig_end} = $row->{geneb_contig_start} + ($RNA_CONTIG_BASE_FLANKING - 1);
	$row->{geneb_contig_end} = $row->{geneb_contig_start} + ($flank_b - 1);

	my $contig_base = join "",
	  substr($upstream_max, - $flank_a),
#	    substr($downstream_max, 0, $RNA_CONTIG_BASE_FLANKING);
	    substr($downstream_max, 0, $flank_b);

	my $contig_max = join "", $upstream_max, $downstream_max;
	my $contig_extended = join "",
	  substr($upstream_max, - $EXTENDED_CHUNK_LENGTH),
	    substr($downstream_max, 0, $EXTENDED_CHUNK_LENGTH);

	$row->{contig} = $contig_base;
	$row->{$F_EXTENDED_CHUNK} = $contig_extended;
	$row->{extended_exon} = $contig_extended;
	# hack
	$row->{extended_max} = $contig_max;

	$row->{genea_full} = $row->{geneb_full} = "";
	# not implemented

	$row->{genea_acc} = $acc_a;
	$row->{genea_acc_preferred} = get_sjpi($SJPI, $acc_a);
	$row->{genea_symbol} = $transcript2gene{$acc_a} || die;
	$row->{genea_pos} = $row->{genea_pos_adj} = $pos_a;
	$row->{genea_feature} = get_feature_tag($rf, $rf_a, $pos_a);
	$row->{genea_chr} = $chr_a;

	if ($intergenic_b) {
	  $row->{geneb_symbol} = "intergenic";
	  # maybe some kind of option here?
	  # downstream conversion does require a string for each
	  # side of the breakpoint.

	  foreach my $f (qw(
			     geneb_acc
			     geneb_acc_preferred
			     geneb_feature
			  )) {
	    $row->{$f} = "";
	  }
	} else {
	  $row->{geneb_acc} = $acc_b;
	  $row->{geneb_acc_preferred} = get_sjpi($SJPI, $acc_b);
	  $row->{geneb_symbol} = $transcript2gene{$acc_b} || die;
	  $row->{geneb_feature} = get_feature_tag($rf, $rf_b, $pos_b);
	}
	$row->{geneb_pos} = $row->{geneb_pos_adj} = $pos_b;
	$row->{geneb_chr} = $chr_b;

	$row->{input_row_number} = $input_row_number;
	$row->{input_row_md5} = $input_row_md5;
	$row->{$F_NOTES} = join ",", @notes_general, @notes_b;
	$row->{$F_GENOME} = $genome;

	generate_pair_summary($row);

	$row->{$F_PROCESSING_TYPE} = "rna_genomic_coordinates";

	my $r_copy = { %{$row} };

	push @rows_b_ok, $r_copy;
	# queue output as we may need to postprocess

      }  # $acc_b


      if (@rows_b_ok > 1 and $RNA_SUPPRESS_DISTANT_ISOFORMS) {
	# check for distant isoforms for gnee B:
	# if coordinates are close to exons for some isoforms
	# but not others, only use the proximal ones
	my @close;
	my @distant;
	foreach my $r (@rows_b_ok) {
	  my $dist = $r->{geneb_coding_distance};
	  if ($dist >= $RNA_DISTANT_THRESHOLD) {
	    push @distant, $r;
	  } else {
	    push @close, $r;
	  }
	}
	if (@close and @distant) {
	  @rows_b_ok = @close;
	}
      }

      foreach my $r (@rows_b_ok) {
	# write final set
	$rpt->end_row($r);
      }

    }  # $acc_a
  }
  $rpt->finish();

  if (%suppressed) {
    printf STDERR "suppressed records:\n";
    foreach my $reason (sort keys %suppressed) {
      printf STDERR "  %s: %d\n", $reason, $suppressed{$reason};
    }
  } else {
    print STDERR "no suppressed records\n";
  }

}

sub prune_query_overlaps {
  my ($hsp_q_set) = @_;
  # /home/medmonso/work/steve/2020_06_13_cicero_contig_extension/2021_01_04_all/crash3/crash.tab

  my $i_end = scalar @{$hsp_q_set} - 1;
  my %idx2intersect;
  for (my $i = 0; $i < $i_end; $i++) {
    my $q1 = new Set::IntSpan join "-", $hsp_q_set->[$i]->range("query");
    my $q2 = new Set::IntSpan join "-", $hsp_q_set->[$i + 1]->range("query");
    my $intersect = intersect $q1 $q2;
    my $isize = $intersect->size();
    printf STDERR "idx %d size %d\n", $i, $isize;
    $idx2intersect{$i} = $isize;
  }

  my ($best_i) = sort {$idx2intersect{$b} <=> $idx2intersect{$a}} keys %idx2intersect;

  my @results = @{$hsp_q_set}[$best_i, $best_i + 1];
  printf STDERR "filtered to q:%d-%d and q:%d-%d\n",
    $results[0]->range("query"),
      $results[1]->range("query");

  return @results;
}

sub log_error {
  my ($code, $suppressed, $notes) = @_;
  $suppressed->{$code}++;
  push @{$notes}, "ERROR:" . $code;
}

sub get_blacklist {
  my $result;
  if ($FLAGS{"no-blacklist"}) {
    # nothing to do
  } elsif (@BLACKLIST_GENOMIC_FILES) {
    my %blacklist;

    foreach my $f_blacklist (@BLACKLIST_GENOMIC_FILES) {
      my $df = new DelimitedFile("-file" => $f_blacklist,
				 "-headers" => 1,
				);
      # must use same header labels as main file being processed (ugh)

      confess "chrA/B fields not defined" unless defined $F_CHR_A and $F_CHR_B;
      confess "posA/B fields not defined" unless defined $F_POS_A and $F_POS_B;

      while (my $row = $df->get_hash()) {
	my $chrA = clean_chr($row->{$F_CHR_A} || dump_die($row, "no $F_CHR_A in $f_blacklist"));
	my $posA = $row->{$F_POS_A};
	my $chrB = clean_chr($row->{$F_CHR_B} || die);
	my $posB = $row->{$F_POS_B};

	$blacklist{$chrA}{$posA}{$chrB}{$posB} = 1;
      }
    }

    $result = \%blacklist;
  } else {
    die "specify -blacklist-genomic FILE [...] or -no-blacklist\n";
  }

  return $result;
  # undef if not being used
}

sub clean_chr {
  my ($chr) = @_;
  $chr =~ s/^chr//i;

  if ($DEMUNGE_CHR_23) {
    $chr = "X" if $chr eq "23";
    # don't use "== 23" in case X already
    $chr = "Y" if $chr eq "24";
  }

  return $chr;
}

sub is_blacklisted {
  my ($blacklist, $row) = @_;

  my $chrA = clean_chr($row->{$F_CHR_A} || dump_die($row, "no $F_CHR_A"));
  my $posA = $row->{$F_POS_A} || dump_die($row, "no $F_POS_A");
  my $chrB = clean_chr($row->{$F_CHR_B} || dump_die($row, "no $F_CHR_B"));
  my $posB = $row->{$F_POS_B} || dump_die($row, "no $F_POS_B");

  my $blacklisted = 0;
  if ($blacklist->{$chrA}{$posA}{$chrB}{$posB}) {
    $blacklisted = 1;
  } elsif ($blacklist->{$chrB}{$posB}{$chrA}{$posA}) {
    $blacklisted = 1;
  }
  return $blacklisted;
}

sub annotate_sources {
  # hack, make generic if used more

  my %f_patterns;

#  $f_patterns{"t3_cicero_good_jneary_20210311.txt.extended.tab.extended_750.pattern.tab"} = [qw(project_name subproject_name disease_sample_name)];
  $f_patterns{"t3_cicero_good_jneary_20210506.txt.extended.tab.extended_750.pattern.tab"} = [qw(project_name subproject_name sample)];
#  $f_patterns{"pp_sv.txt.extended.tab.extended_750.pattern.tab"} = [
#								    qw(Project Sample)
#								  ];
  $f_patterns{"pp_svfusion_2021_05_10.txt.extended.tab.extended_750.pattern.tab"} = [
								    qw(Project sample)
								  ];
  $f_patterns{"summary_benchmark184_CiceroV6_forSteve2_filtered.txt.extended.tab.extended_750.pattern.tab"} = [qw(sample)];

#  $f_patterns{"rrnaseq_clingen_20200924-all.header_fix.tsv.jz_edit.tab.extended.tab.extended_750.pattern.tab"} = [qw(project_name subproject_name disease_sample_name manual_review_call)];
  $f_patterns{"rrnaseq_clingen_20200924-all.header_fix.tsv.jz_edit.tab.extended.tab.extended_750.pattern.tab"} = [qw(project_name subproject_name sample manual_review_call)];

  #
  #  load historic JZ "verdict" calls from previous reviews:
  #
  open(VTMP, "JZVerdict.txt") || die;
  my %verdict;
  while (<VTMP>) {
    chomp;
    s/\r$//;
    # DOS
    my @f = split /\t/, $_;
    die unless @f == 3;
    my ($pair, $source, $verdict) = @f;
    my $key = join "_", $pair, $source;
    $verdict{$key} = $verdict;
  }

  #
  #  temporary hack: retrieve pathogenicity from a later file revision:
  #
  my $df = new DelimitedFile("-file" => "t3_cicero_good_jneary_20210506.txt",
			     "-headers" => 1,
			    );
  my %path;
  while (my $row = $df->get_hash()) {

    my $sample = $row->{sample} || die;
    my $chrA = $row->{chrA} || die;
    my $posA = $row->{posA} || die;
    my $chrB = $row->{chrB} || die;
    my $posB = $row->{posB} || die;

    my $key = join "_", $sample, $chrA, $posA, $chrB, $posB;
    die if $key =~ /_chr\d+/;
    # make sure chroms don't use "chr" prefix

    my $path = $row->{pathogenicity_somatic};
    if ($path) {
      printf STDERR "got %s => %s\n", $key, $path;
      $path{$key} = $path;
    }
  }

  #
  # load patterns to look for (gene pairing portion only):
  #
  my %patterns;
  my %sources;
  foreach my $f_pattern (keys %f_patterns) {
    my $source = file_to_source($f_pattern);
    $sources{$source} = 1;

    my $fields = $f_patterns{$f_pattern};

    my $df = new DelimitedFile("-file" => $f_pattern,
			       "-headers" => 1,
			      );
    while (my $row = $df->get_hash()) {
      my $pattern = $row->{pattern} || next;
      $pattern =~ s/\-\d+$//;

      my @info;
      foreach my $f (@{$fields}) {
	my $v = $row->{$f};
	dump_die($row, "no value for $f in $f_pattern") unless defined $v;
	push @info, $v;
      }

      if ($row->{subproject_name}) {
	# hack, only for t3
	my $chrA = $row->{chrA} || die;
	my $chrB = $row->{chrB} || die;
	my $posA = $row->{posA} || die;
	my $posB = $row->{posB} || die;
	my $sample = $row->{sample} || die;

	foreach ($chrA, $chrB) {
	  dump_die($row, "aha") if s/^chr//;
	}

	my $key = join "_", $sample, $chrA, $posA, $chrB, $posB;
	printf STDERR "lookup %s\n", $key;

	if (my $path = $path{$key}) {
	  push @info, $path;
	}
      }

      my $tag = join "/", @info;

      $patterns{$pattern}{$source}{$tag} = 1;
    }
  }

  my @sources = sort keys %sources;

#  my $infile = "GTEx_summary_St_Jude_RNA_patterns_2021-05-04.txt";
#  my $infile = "clingen_and_GTEx_readpairsum_St_Jude_RNA_patterns_2021-05-04.txt";
  my $infile = "clingen_and_GTEx_readpairsum_St_Jude_RNA_patterns_2021-06-21.txt";

  $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );
  my $outfile = basename($infile) . ".report.tab";
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => [
					   "Verdict",
					   @sources
					  ],
  			      "-auto_qc" => 1,
			     );

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
    my $pair = $row->{"gene-pair"} || die;
    my $source = $row->{"source"} || die;

    my $vkey = join "_", $pair, $source;
    my $verdict = $verdict{$vkey} || "";
    printf STDERR "verdict data for %s: %s\n", $vkey, $verdict ? "y" : "n";
    $row->{Verdict} = $verdict;

    foreach my $source (@sources) {
      my $value = "";
      if (my $info = $patterns{$pair}{$source}) {
	$value = join ", ", sort keys %{$info};
      }
      $row->{$source} = $value;
    }
    $rpt->end_row($row);
  }
  $rpt->finish();
}

sub get_feature_tag {
  my ($rf, $rf_row, $pos, $cache) = @_;

  my $result;
  $result = $cache->{$rf_row}{$pos} if $cache;
  unless ($result) {
    my ($feature, $feature_number, $strand) =
      $rf->get_annotation_for_position(
				       "-row" => $rf_row,
				       "-base" => $pos,
				       "-intergenic-ok" => 1,
				       "-extended" => 1,
				      );
    $feature = "exon" if $feature eq "coding";

    if ($feature eq "intergenic") {
      $result = $feature;
    } else {
      $result = join "_", $feature, $feature_number;
    }
    $cache->{$rf_row}{$pos} = $result if $cache;
  }
  return $result;
}

sub annotate_pattern_exons {
  # short-term hack
  my $f_patterns = "merged.fuzzion_extended_750.tab";
  my $pattern_root = $FLAGS{"annotate-pattern-exons"} || die;
  open(IN, $f_patterns) || die;

  my @p_ordered;
  my %wanted;

  while (<IN>) {
    chomp;
    my @f = split /\t/, $_;
    die unless @f == 2;
    my $pid = $f[0];
    if (index($pid, $pattern_root) == 0) {
      push @p_ordered, $pid;
      $wanted{$pid} = 1;
    }
  }

  my @src_files = glob("*750.pattern.tab");

  my %p2f_a;
  my %p2f_b;
  # pattern to feature

  foreach my $f_src (@src_files) {
    printf STDERR "%s...\n", $f_src;
    my $df = new DelimitedFile("-file" => $f_src,
			     "-headers" => 1,
			     );
    while (my $row = $df->get_hash()) {
      if (my $pid = $row->{pattern}) {
	if ($wanted{$pid}) {
	  my $f_a = $row->{genea_feature} || die;
	  my $f_b = $row->{geneb_feature} || die;
	  $p2f_a{$pid}{$f_a} = $f_a;
	  $p2f_b{$pid}{$f_b} = $f_b;
	}
      }
    }
  }

  my $f_out = sprintf "pattern2exon_%s.txt", $pattern_root;
  my $rpt = new Reporter(
			 "-file" => $f_out,
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   pattern
					   featureA
					   featureB
					)
				      ],
			 "-auto_qc" => 1,
			);

  foreach my $pid (@p_ordered) {
    my %r;
    $r{pattern} = $pid;
    add_feature_info(\%r, $pid, "featureA", \%p2f_a);
    add_feature_info(\%r, $pid, "featureB", \%p2f_b);
    $rpt->end_row(\%r);
  }
  $rpt->finish();

}

sub add_feature_info {
  my ($r, $pid, $label, $hash) = @_;
  if (my $set = $hash->{$pid}) {
    $r->{$label} = join ", ", sort keys %{$set};
  }
}

sub extract_clinical_ids {
  my @infiles;
  foreach my $f (@CLINICAL_SRC_FILE_LISTS) {
    my $set = read_simple_file($f);
    push @infiles, @{$set};
  }

  my @sample_fields = qw(
			  sample
			  disease_sample_name
			  Sample
		       );

  my @project_fields = qw(
			project_name
		       );
  my @subproject_fields = qw(
			   subproject_name
			  );

  #
  #  build clinical sample indexes from Steve's clinical paths:
  #
  my @samples = glob("/research/rgs01/project_space/zhanggrp/ClinicalPilot/common/cicero-post/SJ*");
  my %sample2project;
  foreach my $sample_dir (@samples) {
    # these directories already named after disease sample name
    my $sample = basename($sample_dir);
    $sample2project{$sample}{ClinicalPilot}{ClinicalPilot} = 1;
  }

  foreach my $sp_dir (glob("/research/rgs01/reference/restricted/ClinicalGenomics/clingen_curated/Clinical/*")) {
    my $sp = basename($sp_dir);
    if ($sp =~ /^\d+$/) {
      # project is named for the year only
      foreach my $sample_dir (glob($sp_dir . "/SJ*/SJ*")) {
	# subdirectories including disease codes
	my $sample = basename($sample_dir);
	$sample2project{$sample}{"Clinical"}{$sp} = 1;
      }
    }
  }

  my %samples;
  my %missing;
  my %saw;

  my $SKIP_TARGET = 1;
  # SJ IDs for TARGET data in ProteinPaint export
  printf STDERR "skip SJ TARGET IDs?: %s\n", $SKIP_TARGET ? "y" : "n";

  foreach my $fn (@infiles) {
    die "where is $fn" unless -s $fn;
    my $df = new DelimitedFile("-file" => $fn,
			       "-headers" => 1,
			     );

    my $f_sample = find_header($df, \@sample_fields);
    my $f_project = find_header($df, \@project_fields);
    my $f_subproject = find_header($df, \@subproject_fields);

    die if $f_project and not($f_subproject);

    if ($f_sample) {
      while (my $row = $df->get_hash()) {
	my $sample = $row->{$f_sample} || die;
	next unless $sample =~ /^SJ/;
	# for mixed sources like e.g. ProteinPaint export


	if ($sample =~ s/\-[A-Z]+$//) {
	  next if $SKIP_TARGET;
	}

	$samples{$sample}{$fn} = 1;

	if ($f_project and $f_subproject) {
	  my $project = $row->{$f_project} || die;
	  my $subp = $row->{$f_subproject} || die;

	  $sample2project{$sample}{$project}{$subp} = 1;
	  # merge with directory-scan annotations
	  # TO DO: maybe tweak report to directory contents only??
	}
      }
    } else {
      printf STDERR "ERROR: can't find sample field for file %s\n", $fn;
    }
  }

  my $rpt = new Reporter(
			 "-file" => "sample_ids.txt",
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   sample
					   project
					   subproject
					)
				      ],
			 "-auto_qc" => 1,
			);

  foreach my $sample (sort keys %samples) {
    if ($sample2project{$sample}) {
      my @projects = sort keys %{$sample2project{$sample}};
      die unless @projects;
      die "ambiguous project for $sample" if @projects > 1;

      my ($project) = @projects;
      my @subp = sort keys %{$sample2project{$sample}{$project}};

      printf STDERR "WARNING: ambiguous projects for %s: %s\n", $sample, join ",", @subp if @subp > 1;

      foreach my $subp (@subp) {
	my %r;
	$r{sample} = $sample;
	$r{project} = $project;
	$r{subproject} = $subp;
	$rpt->end_row(\%r);
      }
      $saw{$sample} = 1;
    } else {
      printf STDERR "WARNING: missing sample annotation for %s, files=%s\n", $sample, join ",", sort keys %{$samples{$sample}};
      foreach my $f (keys %{$samples{$sample}}) {
	$missing{$f}++;
      }
    }
  }
  $rpt->finish();

  if (%missing) {
    printf STDERR "files missing sample-to-project associations:\n";
    foreach my $f (sort keys %missing) {
      printf STDERR "  %s: %d\n", $f, $missing{$f};
    }
  }

}


sub find_header {
  my ($df, $fields) = @_;
  my $result;
  foreach my $f (@{$fields}) {
    if (exists $df->headers->{$f}) {
      $result = $f;
      last;
    }
  }
  return $result;
}

sub strand_consistency_check {
  # does the strand reported for the event match the strand of the gene model?
  my ($rf, $event_strand) = @_;
  my $result;
  $event_strand = "" unless defined $event_strand;

  if ($event_strand eq "+" or $event_strand eq "-") {
    # event strand provided and valid
    my $rf_strand = $rf->{strand} || die;
    die "invalid strand value" unless $rf_strand eq "+" or $rf_strand eq "-";
    $result = $event_strand eq $rf_strand;
  } elsif ($event_strand) {
    # bad value
    die "invalid strand value $event_strand";
  } elsif ($RNA_ORT_REQUIRED and not $FLAGS{"rna-no-ort"}) {
    die "missing strand value in user event";
  } else {
    # orientation not required, so skip check
    $result = 1;
  }
  return $result;
}

sub condense_patterns {
  # tighten up ( https://www.youtube.com/watch?v=Wro3bqi4Eb8 )
  my $infile = $FLAGS{condense} || die;
  my $tfw = new TemporaryFileWrangler();

  if (my $preset = $FLAGS{"condense-preset"}) {
    die if $preset < 0.01 or $preset > 1;
    $PATTERN_CONDENSE_MIN_QUERY_OVERLAP_FRAC = $preset;
      $PATTERN_CONDENSE_MIN_BLAST_FRAC_IDENTICAL = $preset;
#      $PATTERN_CONDENSE_MAX_GAP = 20;
#      $PATTERN_CONDENSE_MAX_BREAKPOINT_DRIFT_FUSION = 10;
  }

  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );
  my $outfile_rpt = basename($infile) . ".condense_report.tab";
  my $rpt_detail = $df->get_reporter(
			      "-file" => $outfile_rpt,
			      "-extra" => [
					   qw(
					       keep
					       subsumed_by
					       reason
					    )
					  ],
  			      "-auto_qc" => 1,
			     );

  my $glm = new GeneListMatcher("-allow_reciprocal" => 0);

  my %pid2seq;
  my %by_pair;
  my %seq2id;
  my %id2row;
  my %id2pair;

  #
  #  load sequence-based pattern blacklist:
  #
  my %pattern_blacklist;
  my ($blacklist_fa, $temp_fa);

  if ($FLAGS{debug}) {
    $blacklist_fa = "blacklist_check.fa";
    # for BLAST check of blacklist
    $temp_fa = "pattern_check.fa";
    # for query pattern(s)
  } else {
    $blacklist_fa = $tfw->get_tempfile("-append" => "blacklist");
    $temp_fa = $tfw->get_tempfile("-append" => "pattern_blast");
  }

  my %blacklist_pid2seq;
  if (my $fn = $FLAGS{"blacklist"}) {
    my $df = new DelimitedFile("-file" => $fn,
			       "-headers" => 1,
			      );
    my $wf = new WorkingFile($blacklist_fa);
    my $fh = $wf->output_filehandle();
    while (my $row = $df->get_hash()) {
      my $id = $row->{$F_PATTERN_ID} || die "no column $F_PATTERN_ID";
      my $seq = $row->{$F_PATTERN_SEQUENCE} || die;
      if ($PATTERN_CONDENSE_DUPLICATE_IF_FLANKING_WINDOW_IDENTICAL) {
	$seq = duplicate_window_trim($seq);
      }
#      $blacklisted_pattern_sequences{uc($seq)} = $id;
      die "duplicate ID $id found in blacklist $fn" if $blacklist_pid2seq{$id};

      $blacklist_pid2seq{$id} = $seq;
      write_fasta_seq($fh, $id, $seq);
    }
    $wf->finish();

  } elsif ($FLAGS{"no-blacklist"}) {
  } else {
    die "specify -blacklist FILE | -no-blacklist";
  }

  #
  #  load patterns:
  #

  my @unmergeable;

  while (my $row = $df->get_hash()) {
    my $id = $row->{$F_PATTERN_ID} || die "no column $F_PATTERN_ID";
    my $seq = $row->{$F_PATTERN_SEQUENCE} || die;
    if (($row->{fz2_hint} || "") eq "no_condense") {
      # hint from fusion2breakpoints.pl: pattern should not be condensed
      # because it contains generated interstitial sequence which already
      # contains some ambiguity
      push @unmergeable, $row;
      next;
    }

    if ($PATTERN_CONDENSE_DUPLICATE_IF_FLANKING_WINDOW_IDENTICAL) {
#      printf STDERR "before %s, len=%d: %s\n", $id, length($seq), $seq;
      $seq = duplicate_window_trim($seq);
#      printf STDERR "  after len=%d: %s\n", length($seq), $seq;
    }

#    if (0 and my $bid = $blacklisted_pattern_sequences{uc($seq)}) {
#      $pattern_blacklist{$id} = $bid;
#    }

    die "FATAL ERROR: duplicate record for pattern ID $id" if $id2row{$id};
    $id2row{$id} = $row;
    my $pair_key;

    if ($id =~ /^([\w\-\.]+)\-\d+$/) {
      # v1 pattern IDs, sometimes confounded, e.g. FHAD1-TMEM51-AS1-01
      # - period match is for CAND1.11
      #
     $pair_key = $1;
     my @genes = split /\-/, $pair_key;
     # may fail, e.g. pattern ID KCNJ14-H1-0,
     # "H1-0" is a valid gene symbol  :/
     die unless @genes >= 2;
#     printf STDERR "pair:%s  find %s\n", $pair_key, join ",", @genes;
     if (my $hits = $glm->find(\@genes)) {
       die unless @{$hits} == 1;
       my $hit = $hits->[0];
       my $match_type = $hit->{match_type_gene} || die;
       if ($match_type ne "exact") {
	 my $hit_genes = $hit->{genes} || die;
	 printf STDERR "for %s, using existing pair %s\n",
	   $pair_key,
	     join ",", @{$hit_genes};
	 my $pair_key_raw = $pair_key;
	 $pair_key = join "-", @{$hit_genes};
	 die "broken key" unless $by_pair{$pair_key};
       }
     } else {
#       printf STDERR "add glm set %s\n", join ",", @genes;
       $glm->add_set("-genes" => \@genes);
     }
    } else {
      die "can't parse pattern ID $id";
      # TO DO: use placeholder "other" key in this case?
    }
    die "ID collision $id" if $pid2seq{$id};
    $pid2seq{$id} = $seq;
    $by_pair{$pair_key}{$id} = 1;
    $seq2id{$seq}{$id} = 1;
    $id2pair{$id} = $pair_key;
  }

  my $dt_perfect = new DuplicateTracker();
  my $dt_blast = new DuplicateTracker();

  # identify pure duplicate patterns regardless of gene pairings,
  # e.g. for gene family patterns that produce the same output
  foreach my $seq (keys %seq2id) {
    my @all_ids = keys %{$seq2id{$seq}};
    if (@all_ids > 1) {
      printf STDERR "duplicates %s\n", join ", ", sort @all_ids;
#      die "implement me: perfect duplicate by sequence, change later steps to not act on patterns suppressed here";
      # this can happen if flanking window trimming has been performed,
      # or possibly with input data?
    }

    #
    #  if match spans multiple pair sets, merge them.
    #  e.g. for MELK-C18orf8-01 and MELK-RMC1-01, these are the
    #  same gene with different symbols.
    #
    my %pairs = map {$id2pair{$_} => 1} @all_ids;
    my @pairs = sort keys %pairs;
    if (@pairs > 1) {
      my @ordered = sort {scalar keys %{$by_pair{$b}} <=> scalar keys %{$by_pair{$a}}} @pairs;

      my ($largest, @others) = @ordered;
      printf STDERR "collating %s into %s\n", join(",", @others), $largest;

      foreach my $pair (@others) {
	foreach my $id (keys %{$by_pair{$pair}}) {
	  printf STDERR "move %s from %s to %s\n", $id, $pair, $largest;
	  $by_pair{$largest}{$id} = 1;
	}
	printf STDERR "deleting obsolete pair $pair\n";
	delete $by_pair{$pair};
      }
    }
  }

  my $blast = get_blast();

  #
  # identify duplicates within gene pairings:
  #
  foreach my $pair_key (sort keys %by_pair) {
    printf STDERR "processing pair %s\n", $pair_key;
    my @ids_by_len = sort {length($pid2seq{$a}) <=> length($pid2seq{$b})} keys %{$by_pair{$pair_key}};

    #
    # determine pattern type:
    #
    my $has_fusion = 0;
    my $has_itd = 0;
    foreach my $id (@ids_by_len) {
      my $seq = $pid2seq{$id} || die;
      if ($seq =~ /\[/) {
	$has_fusion = 1;
      } elsif ($seq =~ /\{/) {
	$has_itd = 1;
      } else {
	die "no breakpoint/ITD brackets found in pattern $id, $seq";
      }
    }

    #
    # substring duplicates:
    #
    my @ids = @ids_by_len;
    if (@ids > 1) {
      for (my $i = 0; $i < @ids; $i++) {
	# compare shorter or same-length patterns...
	# (if flanking window trimming is in effect many patterns may
	# have the same length)
	my $i_seq = $pid2seq{$ids[$i]};
	my $i_len = length $i_seq;

	for (my $j = @ids - 1; $j >= 0; $j--) {
	  # ...to longer patterns.  Don't stop at > $i because
	  # there may be some sequences of the same length, and
	  # qualifying comparisons may require a particular order
	  my $j_seq = $pid2seq{$ids[$j]};
	  my $j_len = length $j_seq;
	  next if $i == $j;
	  last if $j_len < $i_len;

#	  printf STDERR "i=%d j=%d\n", $i, $j;

	  my $idx = index($j_seq, $i_seq);
	  if ($idx != -1) {
	    # shorter pattern is a perfect subset of the larger pattern
	    $dt_perfect->mark_duplicate("-from" => $ids[$i],
					"-to" => $ids[$j]);
	    last;
	  }
	}
      }
    }


    if ($FLAGS{blacklist}) {
      #
      #  also BLAST pair set vs. blacklist for near-perfect matches:
      #
      @ids = grep {not(
		       $dt_perfect->is_duplicate($_)
		      )
		 } @ids_by_len;

      printf STDERR "after perfect/blast filter: %s\n", join ",", @ids;
      open(TMPQ, ">" . $temp_fa) || die;
      foreach my $id (@ids) {
	write_fasta_seq(*TMPQ, $id, $pid2seq{$id});
      }
      close TMPQ;

      my $min_q_frac_identical = 0.995;
      # typical pattern is +/- 400 nt so this will tolerate
      # a small number of mismatches only
      my $max_distance = $has_itd ? 0 : 4;
      # for ITDs, don't tolerate any shift in the breakpoint marker.
      #
      # for fusions, hack for MPP6-ATAD1: small differences in
      # breakpoint bracket locations, can happen depending on method,
      # e.g. a protein-based source might detect a small microhomology
      # at the breakpoint and bracked this off, whereas a nucleotide
      # contig approach might be seamless.  Long story short, don't
      # want to be overly strict about this as long as alignment is
      # otherwise very close.

      my $parser = $blast->blast(
				 "-query" => $temp_fa,
				 "-database" => $blacklist_fa
				);
      while (my $result = $parser->next_result()) {
	# one result per query seq
	my $q_pid = $result->query_description();
	my $q_len = $result->query_length();
	while (my $hit = $result->next_hit()) {
	  # hit to a database sequence
	  my $h_pid = $hit->name();
	  while (my $hsp = $hit->next_hsp) {
	    # find a SINGLE HSP that must:
	    # - spans nearly all of the query with very high identity
	    # - have breakpoints in the identical location
	    printf STDERR "    score:%s strand:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s query_span:%d ref_span:%d total_span=%d query_string=%s hit_string=%s\n",
	      $hsp->score,
		$hsp->strand,
		  $hsp->range("query"),
		    $hsp->range("hit"),
		      $hsp->num_identical(),
			$hsp->frac_identical("query"),
			  $hsp->length("query"),
			    $hsp->length("hit"),
			      $hsp->length("total"),
				$hsp->query_string(),
				  $hsp->hit_string();

	    my $min_q_overlap = int($q_len * 0.995);

	    if ($hsp->length("query") >= $min_q_overlap) {
	      my @ni_q = get_n_indexes($hsp->query_string());
	      my @ni_h = get_n_indexes($hsp->hit_string());
	      my $same_bp = 1;

	      my $blacklistable;

	      foreach my $idx (0, 1) {
		my $i1 = $ni_q[$idx];
		my $i2 = $ni_h[$idx];
		if ($i1 != -1 and $i2 != -1) {
		  my $distance = abs($i1 - $i2);
#		  $same_bp = 0 unless $distance == 0;
		  $same_bp = 0 unless $distance <= $max_distance;
		} else {
		  $same_bp = 0;
		}
	      }

	      if ($same_bp) {
		my $fiq = $hsp->frac_identical("query");
		printf STDERR  "CHECK ME duplicate match $q_pid $h_pid, frac_identical query=%f\n", $fiq if $fiq != 1;
		if ($fiq >= $min_q_frac_identical) {
		  $blacklistable = 1;
		} else {
		  #		  die sprintf "CHECK ME REJECT duplicate match $q_pid $h_pid, frac_identical query=%f", $fiq;
		}
	      }

	      $pattern_blacklist{$q_pid} = $h_pid if $blacklistable;

	    }
	  }
	}
      }
    }

    @ids = grep {not(
		     $dt_perfect->is_duplicate($_) or
		     $pattern_blacklist{$_}
		    )
		   } @ids_by_len;
    printf STDERR "after perfect+blacklist filter: %s\n", join ",", @ids;

    my $use_blast = 1;
    $use_blast = 0 if $FLAGS{"condense-no-blast"};

    if (@ids > 1 and $use_blast) {
      #
      # BLAST set vs. itself to find additional duplicates
      #

      # write patterns as FASTA:
      open(FATMP, ">" . $temp_fa) || die;
      foreach my $id (@ids) {
	my $seq = $pid2seq{$id};

#	$seq =~ s/\W+//g;
	# strip out pattern bracketing
	$seq =~ s/\W/N/g;
	# revision: instead, mark with Ns: if the breakpoints are in the
	# same place they'll match, otherwise will count against mismatches

	# TO DO: use different characters for fusion breakpoints
	# and ITD breakpoints, so they will count as mismatches?

	printf FATMP ">%s\n%s\n", $id, $seq;
      }
      close FATMP;

      if (1) {
	printf STDERR "temp fasta: %s\n", $temp_fa;
	copy($temp_fa, "test.fa");
      }

      my $hits = self_blast(
			    "-blast" => $blast,
			    "-fasta" => $temp_fa,
			    "-itd" => $has_itd,
			    "-pid2seq" => \%pid2seq
			 );

      for (my $i = 0; $i < @ids; $i++) {
	# compare shorter patterns...
	my $pid_short = $ids[$i];
#	die if $perfect_dup{$pid_short};
	die if $dt_perfect->is_duplicate($pid_short);
	die if $dt_blast->is_duplicate($pid_short);

	my $i_seq = $pid2seq{$ids[$i]};
	my $i_len = length $i_seq;

	for (my $j = @ids - 1; $j >= 0; $j--) {
	  # ...to same-length or longer patterns
	  my $j_seq = $pid2seq{$ids[$j]};
	  my $j_len = length $j_seq;
	  next if $i == $j;
	  last if $j_len < $i_len;

	  my $pid_long = $ids[$j];
	  if ($hits->{$pid_short}{$pid_long}) {
	    print STDERR "$pid_short vs $pid_long\n";
	    $dt_blast->mark_duplicate("-from" => $pid_short,
					"-to" => $pid_long);

	    last;
	  }
	}
      }
    }  # blast key set vs. itself



  }				# pair



  #
  #  patterns blacklisted by sequence, if additional merging
  #  has occurred, make sure the final target is blacklisted too:
  #
  foreach my $pair_key (sort keys %by_pair) {
    my @ids = sort keys %{$by_pair{$pair_key}};
    foreach my $pid (@ids) {
      if ($pattern_blacklist{$pid}) {
	# pattern has been blacklisted by sequence...
	my @alt;
	foreach my $dt ($dt_perfect, $dt_blast) {
	  if ($dt->is_duplicate($pid)) {
	    push @alt, $dt->get_final_id($pid);
	    # ...but also marked a duplicate by another method...
	  }
	}

	foreach my $pid_alt (@alt) {
	  unless ($pattern_blacklist{$pid_alt}) {
	    # ...so, those equivalent PIDs need to be blacklisted too
	    # target needs to be blacklisted too
	    $pattern_blacklist{$pid_alt} = $pid;
	    printf STDERR "pattern %s blacklisted by sequence, also blacklisting final pattern %s\n", $pid, $pid_alt;
#	    die sprintf "FIX ME: pid %s, pattern blacklist=%s but also another blacklist type in $pair_key! alt=%s", $pid, $pattern_blacklist{$pid}, join ",", @alt;
	  }
	}
      }
    }
  }


  #
  # write output:
  #
  foreach my $pair_key (sort keys %by_pair) {
    my @ids = sort keys %{$by_pair{$pair_key}};

    foreach my $pid (@ids) {
      my $row = $id2row{$pid} || die;

      my $keep = 1;
      my $subsumed;

      my @reasons;
      if (my $bid = $pattern_blacklist{$pid}) {
	$subsumed = "pattern_blacklist";
	$keep = 0;
	push @reasons, sprintf "pattern_blacklist=%s", $bid;
      } elsif ($dt_perfect->is_duplicate($pid)) {
#	$subsumed = $perfect_dup{$pid}) {
	$subsumed = $dt_perfect->get_final_id($pid);
	$keep = 0;
	push @reasons, 'perfect_dup';
      } elsif ($dt_blast->is_duplicate($pid)) {
	$subsumed = $dt_blast->get_final_id($pid);
	$keep = 0;
	push @reasons, 'blast_dup';
      } else {
	$subsumed = "";
      }

      $row->{keep} = $keep;
      $row->{reason} = join ",", @reasons;
      $row->{subsumed_by} = $subsumed;

      $rpt_detail->end_row($row);
    }
  }

  foreach my $row (@unmergeable) {
    $row->{keep} = 1;
    $row->{reason} = "marked_as_do_not_condense";
    $row->{subsumed_by} = "";
    $rpt_detail->end_row($row);
  }

  $rpt_detail->finish();

  condense_from_report($outfile_rpt);
  # step 2

}

sub self_blast {
  # blast patterns against themselves, saving usable hits
  my (%options) = @_;
  my $blast = $options{"-blast"} || die "-blast";
  my $fasta = $options{"-fasta"} || die "-fasta";
#  my $require_identical_breakpoint = $options{"-require-identical-breakpoint"};
  my $pid2seq = $options{"-pid2seq"} || die "-pid2seq";
  my $itd_mode = $options{"-itd"};
  die "-itd" unless defined $itd_mode;

  my $max_breakpoint_drift = $itd_mode ?
    $PATTERN_CONDENSE_MAX_BREAKPOINT_DRIFT_ITD :
      $PATTERN_CONDENSE_MAX_BREAKPOINT_DRIFT_FUSION;
  die "-max-breakpoint-drift" unless defined $max_breakpoint_drift;
  # - for ITDs require very strict breakpoint positioning as matching
  #   in that mode requires read overlap the site.
  # - for fusions be a little more relaxed, esp. as fusions represented
  #   both in CICERO/contig mode in RNA coordinate-only mode may
  #   be reported slightly differently but they are essentially
  #   identical

  my $parser = $blast->blast(
			     "-query" => $fasta,
			     "-database" => $fasta
			    );

  my %usable;

  while (my $result = $parser->next_result()) {
    # one result per query sequence
    # Bio::Search::Result::GenericResult
    my $q_pid = $result->query_description();
    my $q_len = $result->query_length();
    printf STDERR "query %s (%d)\n", $q_pid, $q_len;

    while (my $hit = $result->next_hit()) {
      # hits from this query to a database sequence
      # Bio::Search::Hit::Hit
      # Bio::Search::Hit::GenericHit
      my $h_pid = $hit->name();
      next if $h_pid eq $q_pid;

      my @hsp;
      while (my $hsp = $hit->next_hsp) {
	push @hsp, $hsp;
      }

      printf STDERR "  hit %s\n", $h_pid;
      my $num_identical = 0;
      my $hsp_span_qlen = 0;

      my $si_all = new Set::IntSpan();
      my @hsp_filtered;
      foreach my $hsp (@hsp) {
	my $hsp_qlen = $hsp->length("query");
	my $q_range = join '-', $hsp->range("query");

	printf STDERR "    qlen:%d ident:%d range:%s\n", $hsp_qlen, $hsp->num_identical(), $q_range;

	if ($si_all->size()) {
	  my $si_this = new Set::IntSpan $q_range;
	  if ($si_this lt $si_all) {
#	    die sprintf "hey now range %s subset of %s\n", $q_range, $si_all->run_list();
	    # HSP is a pure subset of previously-parsed range, ignore
	    next;
	  } else {
	    my $intersect = intersect $si_this $si_all;
	    if ($intersect->size()) {
	      next;
	      # hack: ignore intersecting HSPs.
	      # for ITDs, these are often inferior/lesser matches.
	      # if there is a smarter way to deal with these,
	      # it will likely add much more complexity.
	    }
	  }
	}
	push @hsp_filtered, $hsp;

	$si_all = union $si_all $q_range;
#	printf STDERR "current runlist: %s\n", join ",", $si_all;

	$hsp_span_qlen += $hsp_qlen;
	$num_identical += $hsp->num_identical();
      }

      my $q_overlap_frac = $hsp_span_qlen / $q_len;
      die "overlap fraction $q_overlap_frac is > 1, e.g. ITD parse error?"
	if $q_overlap_frac > 1;
      # unpossible

      my $non_identical = 0;
      if (@hsp_filtered == 1) {
	my $hsp = $hsp_filtered[0];
	$non_identical = $hsp->length("query") - $hsp->num_identical();
      }

      unless ($q_overlap_frac >= $PATTERN_CONDENSE_MIN_QUERY_OVERLAP_FRAC
#	      or $non_identical <= $PATTERN_CONDENSE_QUERY_MAX_MISMATCHES
	     ) {
	# - minimum percentage of the query sequence that overlaps the
	#   subject sequence.
	# - also allow a single mismatch through (e.g. for condensing
	#   of short patterns)
	printf STDERR "      reject: query overlap %% (%.3f)\n", $q_overlap_frac;
	next;
      }

      my $max_gap_size = get_max_gap_size(\@hsp);
      if ($max_gap_size > $PATTERN_CONDENSE_MAX_GAP) {
	# hard limit on largest observed gap size
	# (single, not total count gap bases)
	printf STDERR "      reject: max gap size hard limit (%d)\n", $max_gap_size;
	next;
      }

      my $si_query = new Set::IntSpan map {join "-", $_->range("query")} @hsp_filtered;
      my $si_hit = new Set::IntSpan map {join "-", $_->range("hit")} @hsp_filtered;

#      my $frac_identical_query = $num_identical / $q_len;
      # doesn't account for gaps

      my $max_span = max($si_query->size(), $si_hit->size());
      my $frac_identical_query = $num_identical / $max_span;
#      die "hey now multi HSP identical query $frac_identical_query $num_identical $max_span" if @hsp_filtered > 1;

      unless ($frac_identical_query >= $PATTERN_CONDENSE_MIN_BLAST_FRAC_IDENTICAL
#	      or $non_identical <= $PATTERN_CONDENSE_QUERY_MAX_MISMATCHES
	     ) {
	# minimum BLAST alignment identity
	printf STDERR "      reject: query identical %% (%.3f)\n", $frac_identical_query;
	next;
      }

      # compare total query and hit spans, rejecting
      # those with a large difference.  For example if a small pattern
      # mostly overlaps with a larger pattern using multiple HSPs,
      # this implies a large gap in the hit sequence which is not
      # appropriate to consider a duplicate.
      # example:
      # /home/medmonso/work/steve/2020_06_13_cicero_contig_extension/2021_01_04_all/blast/RUNX1-RUNX1T1.tab
      # patterns are largely identical upstream and downstream of
      # the breakpoint marker, which seems like a good duplicate
      # candidate.  However one pattern has a lot more interstitial
      # sequence than the other, which disqualifies it from being
      # considered a true duplicate.
      #
      # Also applies to single HSPs, where the identity is good but
      # a large-ish region at the start/end of the sequence isn't aligned,
      # these can be relevant exon changes between isoforms.
      #
      # This requirement may let us experiment with lowering the BLAST
      # identity percent, for sources containing noiser contigs.
      #

      my $cover_query_len = $si_query->cover()->size();
      my $cover_hit_len = $si_hit->cover()->size();

      my $diff_frac = abs(1 - $cover_query_len / $cover_hit_len);

      if ($diff_frac > $PATTERN_CONDENSE_MULTI_HSP_MAX_QUERY_HIT_SPAN_SIZE_DIFFERENCE) {
	printf STDERR "      reject: %s HSP, spans are too far apart (%.3f)\n",
	  (scalar @hsp_filtered == 1 ? "single" : "multi"),
	    $diff_frac;
	next;
      }

      #
      #  compare position of breakpoints in alignment:
      #
      my $bp_max_distance = 0;
      my $bp_count_found = 0;
      my %hsp_found;

      foreach my $hsp (@hsp_filtered) {
	# each HSP may contain 0, 1 or 2 breakpoints.
	# if there are multiple HSPs, there is no guarantee both
	# breakpoints will appear in one of them.
	# See also rescue procedure below.
	my @ni_q = get_n_indexes($hsp->query_string());
	my @ni_h = get_n_indexes($hsp->hit_string());

	foreach my $idx (0, 1) {
	  my $i1 = $ni_q[$idx];
	  my $i2 = $ni_h[$idx];
	  if ($i1 != -1 and $i2 != -1) {
	    # require breakpoint be observed in both alignments
	    my $distance = abs($i1 - $i2);
	    printf STDERR "    idx %d %d %d, dist=%d\n", $idx, $i1, $i2, $distance;
	    $bp_max_distance = $distance if $distance > $bp_max_distance;
	    $bp_count_found++;
	    $hsp_found{$hsp} = 1;
	  }
	}
      }

      if ($bp_count_found < 2 and @hsp_filtered > 1) {
	my @hsp_q = sort {$a->start("query") <=>
			    $b->start("query")} @hsp_filtered;
	# ensure HSPs are sorted by query order so we can examine the gaps
	# between them

	my $last_i = @hsp_q - 2;
	for (my $i = 0; $i <= $last_i; $i++) {
	  my %found_bp;

	  foreach my $ref (
			   [ "query", $q_pid ],
			   [ "hit", $h_pid ]
			  ) {
	    my ($type, $id) = @{$ref};
	    my $sequence = $pid2seq->{$id} || die;
	    my $end_this = $hsp_q[$i]->end($type);
	    my $start_next = $hsp_q[$i + 1]->start($type);
	    my $search_start = $end_this + 1;
	    my $search_end = $start_next - 1;
	    my $len = ($search_end - $search_start) + 1;
	    if ($len <= $max_breakpoint_drift) {
	      # only allow a small gap between HSP regions
	      for (my $j = $search_start; $j <= $search_end; $j++) {
		my $char = substr($sequence, $j - 1, 1);
		if ($BREAKPOINT_CHAR{$char}) {
		  # found breakpoint character between HSPs
		  $found_bp{$type} = 1;
		} elsif ($char !~ /[acgtn]/i) {
		  die "WTF: non-nucleotide, non-breakpoint character found: $char";
		}
	      }
	    }
	  }

	  if (scalar keys %found_bp == 2) {
	    # need test case!  This file is actually a bad example
	    # because interstitial sequence is much longer in one
	    # pattern vs. the other:
	    # /home/medmonso/work/steve/2020_06_13_cicero_contig_extension/2021_01_04_all/blast/RUNX1-RUNX1T1.tab
	    die "BREAKPOINT HSP RESCUE, CHECK ME";
	  }
	}
      }

      if (scalar keys %hsp_found > 1) {
	# plausible breakpoint evidence found, but involving multiple HSPs.
	# This implies some significant difference between the sequences,
	# possibly an alternative transcript sequence which doesn't
	# align directly.  These should not be considered equivalent.
	printf STDERR "plausible breakpoint match for %s => %s, but involving multiple HSPs, so disqualified\n", $q_pid, $h_pid;
	next;
      }

      if ($bp_count_found == 2) {
	printf STDERR "max distance for %s vs %s: %d\n", $q_pid, $h_pid, $bp_max_distance;
	if ($bp_max_distance > $max_breakpoint_drift) {
	  printf STDERR "      reject: %s breakpoints too far away (%d)\n",
	    ($itd_mode ? "ITD" : "fusion"),
	      $bp_max_distance;
	  next;
	}
      } elsif (0 and @hsp_filtered > 1) {
	# exempt this situation from check because Ns sometimes are
	# exempted from HSP edges, e.g.  RUNX1-RUNX1T1-133
	# vs. RUNX1-RUNX1T1-140.
	#
	# FAIL: there are counter-examples, e.g. CBFB-MYH11-52 and
	# CBFB-MYH11-62 that do NOT look good
	printf STDERR "skipping multi-HSP N distance check disqualification between %s and %s\n", $q_pid, $h_pid;
      } else {
	print STDERR "      reject: breakpoints not found in alignment\n";
	next;
      }

      printf STDERR "CHECK THIS: usable match found involving multiple HSPs!" if @hsp_filtered > 1;
      # DANGER?  what if a short query sequence mostly maps to a larger
      # subject sequence, but with a large gap?  If there are multiple
      # HSPs, compare total span lengths and die if a large difference?

      printf STDERR "    usable!\n";
      $usable{$q_pid}{$h_pid} = 1;

    }
  }
  return \%usable;
}

sub self_blast_single_hsp {
  # blast patterns against themselves, saving usable hits
  # OBSOLETE: single HSP version
  my (%options) = @_;
  my $blast = $options{"-blast"} || die "-blast";
  my $fasta = $options{"-fasta"} || die "-fasta";
#  my $require_identical_breakpoint = $options{"-require-identical-breakpoint"};
  my $max_breakpoint_drift = $options{"-max-breakpoint-drift"};
  die "-max-breakpoint-drift" unless defined $max_breakpoint_drift;

  my $parser = $blast->blast(
			     "-query" => $fasta,
			     "-database" => $fasta
			    );

  my %usable;

  while (my $result = $parser->next_result()) {
    # one result per query sequence
    # Bio::Search::Result::GenericResult
    my $q_pid = $result->query_description();
    my $q_len = $result->query_length();
    printf STDERR "query %s (%d)\n", $q_pid, $q_len;

    while (my $hit = $result->next_hit()) {
      # hits from this query to a database sequence
      # Bio::Search::Hit::Hit
      # Bio::Search::Hit::GenericHit
      my $h_pid = $hit->name();
      next if $h_pid eq $q_pid;

      while (my $hsp = $hit->next_hsp) {
	# Bio::Search::HSP::PSLHSP
	# to be considered a duplicate, a single HSP must meet QC requirements:

	my $non_identical = $hsp->length("query") - $hsp->num_identical();
	my $q_overlap_frac = $hsp->length("query") / $q_len;
	printf STDERR "  hit:%s qlen:%d HSP_qlen:%d frac:%f ident:%d gap:%d non:%d\n", $h_pid, $q_len, $hsp->length("query"), $q_overlap_frac, $hsp->num_identical(), $hsp->gaps(), $non_identical;

#	next unless $hsp->length("query") == $q_len;
	unless ($q_overlap_frac >= $PATTERN_CONDENSE_MIN_QUERY_OVERLAP_FRAC) {
	  # minimum percentage of the query sequence that overlaps the
	  # subject sequence.
	  print STDERR "    reject: query overlap %\n";
	  next;
	}

	unless ($hsp->frac_identical("query") >= $PATTERN_CONDENSE_MIN_BLAST_FRAC_IDENTICAL) {
	  # minimum BLAST alignment identity
	  print STDERR "    reject: query identical %\n";
	  next;
	}

	my $max_gap_size = get_max_gap_size($hsp);
	if ($max_gap_size > $PATTERN_CONDENSE_MAX_GAP) {
	  # hard limit on largest observed gap size
	  # (individual gap size, not total gap bases)
	  print STDERR "    reject: gap hard limit\n";
	  next;
	}

	# examine breakpoint positions and reject if sites are too far away:
	my @ni_q = get_n_indexes($hsp->query_string());
	my @ni_h = get_n_indexes($hsp->hit_string());
	my $max_distance = 0;
	foreach my $idx (0, 1) {
	  my $i1 = $ni_q[$idx];
	  my $i2 = $ni_h[$idx];
	  my $distance = abs($i1 - $i2);
	  printf STDERR "idx %d %d %d, dist=%d\n", $idx, $i1, $i2, $distance;
	  $max_distance = $distance if $distance > $max_distance;
	}
	printf STDERR "max distance for %s vs %s: %d\n", $q_pid, $h_pid, $max_distance;
	if ($max_distance > $max_breakpoint_drift) {
	  print STDERR "    reject: breakpoints too far away\n";
	  next;
	}

# 	if ($require_identical_breakpoint) {
# #	  printf STDERR "qs: %s\n", $hsp->query_string();
# #	  printf STDERR "hs: %s\n", $hsp->hit_string();
# 	  my @ni_q = get_n_indexes($hsp->query_string());
# 	  my @ni_h = get_n_indexes($hsp->hit_string());
# #	  printf STDERR "Ns: q=%d,%d, h=%d,%d\n", @ni_q, @ni_h;

# 	  if ($ni_q[0] == $ni_h[0] and $ni_q[1] == $ni_h[1]) {
# 	    # breakpoint details identical, usable
# 	  } else {
# 	    # different breakpoint site, skip
# 	    next;
# 	  }
# 	}

	printf STDERR "    usable. %s strand:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s query_span:%d ref_span:%d total_span=%d\n",
	  $h_pid,
	    $hsp->strand,
	      $hsp->range("query"),
		$hsp->range("hit"),
		  $hsp->num_identical(),
		    $hsp->frac_identical("query"),
		      $hsp->length("query"),
			$hsp->length("hit"),
			  $hsp->length("total");

	$usable{$q_pid}{$h_pid} = 1;
      }
    }
  }
  return \%usable;
}

sub generate_blacklist {
  # hack for JZ blacklist based on early pattern set.
  # gene pairs only rather than pattern IDs.
  # excel sheet does not contain annotations for one source,
  # rrnaseq_clingen_20200924-all.header_fix.tsv.jz_edit.tab.extended.tab.extended_750.pattern.tab.

  die "chr/pos/ort fields not defined; suggest using -cicero, -crest, etc." unless ($F_CHR_A and $F_POS_A and $F_CHR_B and $F_POS_B);

  my %src;
  my %already_blacklisted;
  my %pair2sample;

  foreach my $f_src (qw(
			 t3_cicero_good_jneary_20210311.txt.extended.tab.extended_750.pattern.tab
			 pp_sv.txt.extended.tab.extended_750.pattern.tab
			 rrnaseq_clingen_20200924-all.header_fix.tsv.jz_edit.tab.extended.tab.extended_750.pattern.tab
		      )) {
    # map gene pairs in sources to blacklist coordinates
    my $df = new DelimitedFile("-file" => $f_src,
			     "-headers" => 1,
			     );
    my $f_sample;
    my $h = $df->headers();
    foreach my $f (qw(
		       disease_sample_name
		       sample
		       Sample
		    )) {
      if (exists $h->{$f}) {
	$f_sample = $f;
	last;
      }
    }
    dump_die($df->headers, "no sample field found") unless $f_sample;

    while (my $row = $df->get_hash()) {
      die unless exists $row->{pattern};
      my $pair = $row->{pattern} || next;
      my $pair_raw = $pair;
      my $already_blacklisted;
      if ($pair =~ s/\-\d+$//) {
	# OK
      } elsif ($pair eq "suppressed_by_blacklist" or
	       $pair_raw eq "1") {
	# - for suppressed entries, try to save anyway, OK if there
	#   are duplicates across lists
	# - value of 1: bug?  Not sure why this happens
	$already_blacklisted = 1 if $pair eq "suppressed_by_blacklist";
	my $sym_a = $row->{"genea_symbol"} || die;
	my $sym_b = $row->{"geneb_symbol"} || die;
	$pair = join "-", $sym_a, $sym_b;
	# hack
      } else {
	dump_die($row, "can't strip suffix from $pair_raw");
      }

      my $chrA = clean_chr($row->{chrA} || die);
      my $chrB = clean_chr($row->{chrB} || die);
      my $posA = $row->{posA} || die;
      my $posB = $row->{posB} || die;

      my $bl = join ".", $chrA, $posA, $chrB, $posB;

      $src{$pair}{$bl} = 1;
      $already_blacklisted{$pair}{$bl} = 1 if $already_blacklisted;

      my $sample = $row->{$f_sample} || die;
      $pair2sample{$pair}{$sample} = 1;

    }
  }

  # iterate through JZ's blacklist spreadsheet and export matches:

  my $f_in = "GTEx_blacklist.txt";
  my $df = new DelimitedFile("-file" => $f_in,
			     "-headers" => 1,
			     );

  my %all_bl;

  my $f_pair = "gene-pair";

  my %all_pairs;

  while (my $row = $df->get_hash()) {
    my $pair = $row->{$f_pair} || die;
    $all_pairs{$pair} = 1;

    if (my $set = $src{$pair}) {
      foreach my $bl (keys %{$set}) {
	$all_bl{$bl} = $pair;
      }
    }
  }

  #
  #  write results:
  #
  my $f_out = $f_in . ".out.txt";
  my $rpt = new Reporter(
			 "-file" => $f_out,
			 "-delimiter" => "\t",
			 "-labels" => [
				       (
					$f_pair,
#					"previously_blacklisted",
					$F_CHR_A,
					$F_POS_A,
					$F_CHR_B,
					$F_POS_B,
					"samples",
					)
				      ],
			 "-auto_qc" => 1,
			);

  foreach my $bl (sort keys %all_bl) {
    my %r;
    my $pair = $all_bl{$bl};
    $all_pairs{$pair} = 2;
    # results found for this pair
    $r{$f_pair} = $all_bl{$bl};
    @r{$F_CHR_A, $F_POS_A, $F_CHR_B, $F_POS_B} = split /\./, $bl;
    $r{"previously_blacklisted"} = $already_blacklisted{$pair}{$bl} ? "y" : "n";
    $r{samples} = join ",", sort keys %{$pair2sample{$pair}};
    $rpt->end_row(\%r);
  }
  $rpt->finish();

  # QC: did we find source data for all pairs in spreadsheet?
  foreach my $pair (sort keys %all_pairs) {
    if ($all_pairs{$pair} == 1) {
      printf STDERR "ERROR: no results for %s\n", $pair;
    }
  }

}

sub generate_full_gene_pairs {
  my $infile = $FLAGS{"full-gene-pair"} || die "-full-gene-pair";

  my $rf = get_refflat();
  my $gsm = init_gsm($rf);

  #
  #  first pass -- find all refSeqs needed:
  #
  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			    );
  my %needed_rf;
  while (my $row = $df->get_hash()) {
    my ($gene_a, $gene_b) = get_gene_a_b($row);
    $input_row_number++;

    my $can_process = 1;
    my %transcript2gene;
    my @notes_general;

    if ($gene_a and $gene_b) {
      #
      # TO DO: centralize/share this code?  check usability and return iterator
      # of transcript combinations?
      #
      my $gene_a_transcript2rf = get_transcripts($rf, $gsm, $gene_a, \@notes_general, \%transcript2gene, "-no-chr" => 1);
      my $gene_b_transcript2rf = get_transcripts($rf, $gsm, $gene_b, \@notes_general, \%transcript2gene, "-no-chr" => 1);
      # mapping of transcripts to refFlat record entries

      transcript_filter($gene_a_transcript2rf, $gene_a, $FLAGS{"restrict-nm-a"});
      transcript_filter($gene_b_transcript2rf, $gene_b, $FLAGS{"restrict-nm-b"});

      foreach my $set ($gene_a_transcript2rf, $gene_b_transcript2rf) {
	foreach my $acc (sort keys %{$set}) {
	  $needed_rf{$acc} = 1;
	}
      }
    }
  }

  #
  #  fetch GenBank records for needed refseqs:
  #
  my $eu = new EUtilities();
  my $result_files = $eu->fetch_genbank("-ids" => [sort keys %needed_rf]);

  my %acc2cds;
  foreach my $file (@{$result_files}) {
    my $fh = universal_open($file);
    my $stream = Bio::SeqIO->new("-fh" => $fh,
				 -format => 'GenBank');
    while (my $record = $stream->next_seq()) {
      # Bio::SeqIO::genbank
      my $info = parse_genbank($record);
      $acc2cds{$info->{accession}} = $info->{sequence};
    }
  }
  # DRAWBACK here is incorporating genomic-based annotations e.g. mappability


  #
  #  second pass: main processing
  #
  $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			    );

  my $outfile = basename($infile) . ".extended.tab";

  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => [
					   $F_PATTERN_ID,
					   $F_PATTERN_SEQUENCE,
					   "input_row_number",
					   "genea_acc",
					   "geneb_acc",
					  ],
			      "-extra-prepend" => 1,
			      "-auto_qc" => 1,
			     );

  my $input_row_number = 0;
  my %pattern_counter;

  while (my $row = $df->get_hash()) {
    my ($gene_a, $gene_b) = get_gene_a_b($row);
    $input_row_number++;

    my %transcript2gene;
    my @notes_general;

    my $can_process = 1;

    if ($gene_a and $gene_b) {
      #
      # TO DO: centralize/share this code?  check usability and return iterator
      # of transcript combinations?
      #
      my $gene_a_transcript2rf = get_transcripts($rf, $gsm, $gene_a, \@notes_general, \%transcript2gene, "-no-chr" => 1);

      my $gene_b_transcript2rf = get_transcripts($rf, $gsm, $gene_b, \@notes_general, \%transcript2gene, "-no-chr" => 1);
      # mapping of transcripts to refFlat record entries

      transcript_filter($gene_a_transcript2rf, $gene_a, $FLAGS{"restrict-nm-a"});
      transcript_filter($gene_b_transcript2rf, $gene_b, $FLAGS{"restrict-nm-b"});
      if ($gene_a eq $gene_b) {
	$can_process = 0;
	log_error("only_fusions_supported", \%suppressed, \@notes_general);
      }

      unless (keys %{$gene_a_transcript2rf}) {
	$can_process = 0;
	log_error("missing_transcripts_A", \%suppressed, \@notes_general);
      }

      unless (keys %{$gene_b_transcript2rf}) {
	$can_process = 0;
	log_error("missing_transcripts_B", \%suppressed, \@notes_general);
	# potentially double-counting
      }

      if ($can_process) {
	# both with refSeqs?  Maybe just use GenBank records?

	foreach my $acc_a (sort keys %{$gene_a_transcript2rf}) {
	  my $rf_a = $gene_a_transcript2rf->{$acc_a} || die;

	  foreach my $acc_b (sort keys %{$gene_b_transcript2rf}) {
	    my $rf_b = $gene_b_transcript2rf->{$acc_b} || die;

	    my $upstream = $acc2cds{$acc_a} || die;
	    my $downstream = $acc2cds{$acc_b} || die;
	    my $pattern = sprintf '%s][%s', $upstream, $downstream;

	    my $counter = ++$pattern_counter{$gene_a}{$gene_b};
	    my $pid = sprintf '%s-%s-%02d', $gene_a, $gene_b, $counter;

	    $row->{$F_PATTERN_ID} = $pid;
	    $row->{$F_PATTERN_SEQUENCE} = $pattern;
	    # HACK: generating pattern IDs directly like this will clash
	    # with standard pattern generation / merging procedure.
	    # OTOH these full-length patterns are fundamentally different.

	    $row->{input_row_number} = $input_row_number;
	    $row->{genea_acc} = $acc_a;
	    $row->{geneb_acc} = $acc_b;

	    $rpt->end_row($row);
	  }
	}
      } else {
	die "can't process";
      }

    }
  }
  $rpt->finish();
}

sub init_gsm {
  # initialize refFlat gene symbol disambiguation
  my ($rf) = @_;
  my $gsm = new_gsm_lite("-hgnc" => $FLAGS{"hgnc"});
  my $genes = $rf->get_genes();
  foreach my $g (@{$genes}) {
    $gsm->add_gene("-gene" => $g);
  }
  return $gsm;
}


sub parse_genbank {
  # Bio::SeqIO::genbank
  my ($record) = @_;

  my $cds_count = 0;
  my ($cds_start, $cds_end);

  my %r;
  $r{accession} = $record->accession_number();
  $r{version} = $record->version();

  my $origin_seq = $record->seq() || die "can't get sequence";
  # POS documentation for Bio::SeqIO::genbank says Origin is available
  # via annotation() / Bio::Annotation::Collection but this is not
  # necessarily the case.  Code looks like it *sometimes* adds an
  # Annotation field, but in my tests it is actually stored in RichSeq
  # "seq" method instead

  my %genes;
  foreach my $feature ($record->get_SeqFeatures()) {
    # Bio::SeqFeature::Generic

    my $ptag = $feature->primary_tag();
    if ($ptag eq "CDS") {
      $cds_count++;
      $cds_start = $feature->start();
      $cds_end = $feature->end();
    }
  }
  die "multiple CDS" if $cds_count > 1;
  # will need to return multiple records if this happens
  die "no CDS defined" unless $cds_start and $cds_end;

  my $cds = substr($origin_seq, $cds_start - 1, ($cds_end - $cds_start) + 1);
  # TO DO: additional iterations?
  # - maybe include 5' UTR for geneB?

  $r{sequence} = $cds;

  return \%r;
}

sub get_n_indexes {
  my ($align) = @_;

  my $i1 = index($align, "N");
  my $i2 = -1;
  if ($i1 != -1) {
    $i2 = index($align, "N", $i1 + 1);
  }
  return ($i1, $i2);
}

sub summarize_samples {

  my @h_sample = qw(
		     sample
		     disease_sample_name
		     Sample
		  );

  my %disease_codes;

  foreach my $listfile (@SAMPLE_SRC_FILES) {
    my $set = read_simple_file($listfile);
    foreach my $file (@{$set}) {
      my $df = new DelimitedFile("-file" => $file,
				 "-headers" => 1,
				);
      my $h_sample;
      foreach my $f (@h_sample) {
	if (defined $df->headers->{$f}) {
	  $h_sample = $f;
	  last;
	}
      }

      if ($h_sample) {
	while (my $row = $df->get_hash()) {
	  my $sample = $row->{$h_sample} || dump_die($row, "no $h_sample");
	  if ($sample =~ /^SJ([A-Z]+)\d+/) {
	    $disease_codes{$file}{$1}++;
	  }
	}
      } else {
	printf STDERR "no sample header found in %s\n", $file;
      }
    }
  }

  foreach my $file (sort keys %disease_codes) {
    printf "%s:\n", $file;
    my $set = $disease_codes{$file};

    foreach my $code (sort {$set->{$b} <=> $set->{$a}} keys %{$set}) {
      printf "  %10s: %d\n", $code, $set->{$code};
    }
  }
}

sub generate_full_gene_junctions {
  my $infile = $FLAGS{"full-gene-junctions"} || die "-full-gene-junctions";

  my $rf = get_refflat();
  my $gsm = init_gsm($rf);

  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			    );

  my $outfile = basename($infile) . ".spec.tab";
  # output needs to be run through RNA pattern generation code

  die "chr/pos/ort fields not defined; suggest using -cicero, -crest, etc." unless ($F_CHR_A and $F_POS_A and $F_CHR_B and $F_POS_B);

  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => [
					   "genea_acc",
					   $F_CHR_A,
					   $F_POS_A,

					   "geneb_acc",
					   $F_CHR_B,
					   $F_POS_B,
					  ],
			      "-auto_qc" => 1,
			     );

  my %pattern_counter;

  while (my $row = $df->get_hash()) {
    my ($gene_a, $gene_b) = get_gene_a_b($row);

    my %transcript2gene;
    my @notes_general;

    my $can_process = 1;

    if ($gene_a and $gene_b) {
      #
      # TO DO: centralize/share this code?  check usability and return iterator
      # of transcript combinations?
      #
      my $gene_a_transcript2rf = get_transcripts($rf, $gsm, $gene_a, \@notes_general, \%transcript2gene, "-no-chr" => 1);

      my $gene_b_transcript2rf = get_transcripts($rf, $gsm, $gene_b, \@notes_general, \%transcript2gene, "-no-chr" => 1);
      # mapping of transcripts to refFlat record entries

      transcript_filter($gene_a_transcript2rf, $gene_a, $FLAGS{"restrict-nm-a"});
      transcript_filter($gene_b_transcript2rf, $gene_b, $FLAGS{"restrict-nm-b"});
      if ($gene_a eq $gene_b) {
	$can_process = 0;
	log_error("only_fusions_supported", \%suppressed, \@notes_general);
      }

      unless (keys %{$gene_a_transcript2rf}) {
	$can_process = 0;
	log_error("missing_transcripts_A", \%suppressed, \@notes_general);
      }

      unless (keys %{$gene_b_transcript2rf}) {
	$can_process = 0;
	log_error("missing_transcripts_B", \%suppressed, \@notes_general);
	# potentially double-counting
      }

      my $STEP = 3;
      my $SKIP_START_END = 25;

      my %saw;

      if ($can_process) {
	# both with refSeqs?  Maybe just use GenBank records?

	foreach my $acc_a (sort keys %{$gene_a_transcript2rf}) {
	  my $rf_a = $gene_a_transcript2rf->{$acc_a} || die;
	  my $a_chr = $rf_a->{chrom} || die;
	  my $a_digest = generate_breakpoints($rf_a, $STEP, $SKIP_START_END);

	  foreach my $acc_b (sort keys %{$gene_b_transcript2rf}) {
	    my $rf_b = $gene_b_transcript2rf->{$acc_b} || die;
	    my $b_chr = $rf_b->{chrom} || die;
	    my $b_digest = generate_breakpoints($rf_b, $STEP, $SKIP_START_END);

	    for (my $i = 0; $i < @{$a_digest}; $i++) {
	      my $a_pos = $a_digest->[$i];
	      for (my $j = 0; $j < @{$b_digest}; $j++) {
		my $b_pos = $b_digest->[$j];

		if ($saw{$a_chr}{$a_pos}{$b_chr}{$b_pos}) {
#		  die "dup";
		} else {
		  my %r;
		  $r{$F_GENE_A} = $gene_a;
		  $r{$F_GENE_B} = $gene_b;
		  $r{genea_acc} = $acc_a;
		  $r{geneb_acc} = $acc_b;
		  $r{$F_CHR_A} = $a_chr;
		  $r{$F_CHR_B} = $b_chr;
		  $r{$F_POS_A} = $a_pos;
		  $r{$F_POS_B} = $b_pos;

		  $rpt->end_row(\%r);

		  $saw{$a_chr}{$a_pos}{$b_chr}{$b_pos} = 1;
		}
	      }
	    }
	  }
	}
      } else {
	die "can't process";
      }

    }
  }
  $rpt->finish();
}

sub generate_breakpoints {
  my ($rf, $step, $skip_start_end) = @_;

  my $cds_start = $rf->{cdsStart} || die;
  my $cds_end = $rf->{cdsEnd} || die;
  my $verbose = 0;

  die if $cds_end < $cds_start;
  # in refFlat, "start" and "end" are genomic rather than transcript oriented
  # ...or should be

  printf STDERR "CDS %s raw:%d-%d ", $rf->{name}, $cds_start, $cds_end if $verbose;

  $cds_start += $skip_start_end;
  $cds_end -= $skip_start_end;
  printf STDERR " adj:%d-%d\n", $cds_start, $cds_end if $verbose;

  my $exons = $rf->{exons} || die;
  my @positions;
  foreach my $exon (@{$exons}) {
    my $start = $exon->{start} || die;
    my $end = $exon->{end} || die;
    printf STDERR "exon %d-%d\n", $start, $end if $verbose;
    if ($start >= $cds_start and $end <= $cds_end) {
      # in relevant coding interval
      for (my $i = $start; $i <= $end; $i += $step) {
	print STDERR "  $i\n" if $verbose;
	push @positions, $i;
      }
    }
  }

  return \@positions;

}


sub find_gene_intervals {
  my $infile = $FLAGS{"find-gene-intervals"} || die "-find-gene-intervals";

  my $rf = get_refflat();
  my $gsm = init_gsm($rf);

  my $JOIN_TAG = ">";

  $df = new DelimitedFile(
			  "-file" => $infile,
			  "-headers" => 1,
			 );

  die "chr/pos/ort fields not defined; suggest using -cicero, -crest, etc." unless ($F_CHR_A and $F_POS_A and $F_CHR_B and $F_POS_B);

  my %intervals;
  my %interval2feature;
  my %intervals2sites;

  my %genes_a;
  my %genes_b;

  while (my $row = $df->get_hash()) {
    my ($gene_a, $gene_b) = get_gene_a_b($row);

    my %transcript2gene;
    my @notes_general;

    my $can_process = 1;

    my $query_chr_a = clean_chr($row->{$F_CHR_A} || die);
    my $query_chr_b = clean_chr($row->{$F_CHR_B} || die);
    my $query_pos_a = clean_chr($row->{$F_POS_A} || die);
    my $query_pos_b = clean_chr($row->{$F_POS_B} || die);

    $genes_a{$gene_a} = 1;
    $genes_b{$gene_b} = 1;

    if ($gene_a and $gene_b) {
      #
      # TO DO: centralize/share this code?  check usability and return iterator
      # of transcript combinations?
      #
      my $gene_a_transcript2rf = get_transcripts($rf, $gsm, $gene_a, \@notes_general, \%transcript2gene, "-no-chr" => 1);

      my $gene_b_transcript2rf = get_transcripts($rf, $gsm, $gene_b, \@notes_general, \%transcript2gene, "-no-chr" => 1);
      # mapping of transcripts to refFlat record entries

      transcript_filter($gene_a_transcript2rf, $gene_a, $FLAGS{"restrict-nm-a"});
      transcript_filter($gene_b_transcript2rf, $gene_b, $FLAGS{"restrict-nm-b"});
      if ($gene_a eq $gene_b) {
	$can_process = 0;
	log_error("only_fusions_supported", \%suppressed, \@notes_general);
      }

      unless (keys %{$gene_a_transcript2rf}) {
	$can_process = 0;
	log_error("missing_transcripts_A", \%suppressed, \@notes_general);
      }

      unless (keys %{$gene_b_transcript2rf}) {
	$can_process = 0;
	log_error("missing_transcripts_B", \%suppressed, \@notes_general);
	# potentially double-counting
      }

      if ($can_process) {
	foreach my $acc_a (sort keys %{$gene_a_transcript2rf}) {
	  my $rf_a = $gene_a_transcript2rf->{$acc_a} || die;
	  my $a_chr = clean_chr($rf_a->{chrom} || die);
	  next unless $a_chr eq $query_chr_a;

	  # find A exon
	  my ($interval_a, $feature_a) = find_interval_and_feature($rf, $rf_a, $query_pos_a);
	  if ($interval_a) {
	    printf STDERR "OK exon for A $query_pos_a $acc_a\n";
	  } else {
	    printf STDERR "can't find exon for A $query_pos_a $acc_a\n";
	    next;
	  }

	  my $interval_a_full = sprintf '%s:%d-%d', $a_chr, @{$interval_a};
	  $interval2feature{$interval_a_full}{$feature_a} = 1;

	  foreach my $acc_b (sort keys %{$gene_b_transcript2rf}) {
	    my $rf_b = $gene_b_transcript2rf->{$acc_b} || die;
	    my $b_chr = clean_chr($rf_b->{chrom} || die);
	    next unless $b_chr eq $query_chr_b;

	    # find B exon interval
	    my ($interval_b, $feature_b) = find_interval_and_feature($rf, $rf_b, $query_pos_b);

	    if ($interval_b) {
	      printf STDERR "OK exon for B $query_pos_b $acc_b\n";
	    } else {
	      printf STDERR "can't find exon for B $query_pos_b $acc_b\n";
	      next;
	    }

	    my $interval_b_full = sprintf '%s:%d-%d', $b_chr, @{$interval_b};
	    $interval2feature{$interval_b_full}{$feature_b} = 1;

	    my $joined = join $JOIN_TAG, $interval_a_full, $interval_b_full;
	    $intervals{$joined}++;

	    my $key = join ".", $query_chr_a, $query_pos_a, $query_chr_b, $query_pos_b;

	    $intervals2sites{$joined}{$key}++;
	    # track breakpoints for diversity measurement

	  }
	}
      } else {
	die "can't process";
      }
    }
  }

  die "multiple genes for A" if scalar keys %genes_a > 1;
  die "multiple genes for B" if scalar keys %genes_b > 1;

  my ($gene_a) = keys %genes_a;
  my ($gene_b) = keys %genes_b;

  my $outfile = basename($infile) . ".intervals.tab";
  my $rpt = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   range
					   count
					   unique_junctions

					   geneA
					   chrA
					   startA
					   endA
					   featureA

					   geneB
					   chrB
					   startB
					   endB
					   featureB
					)
				      ],
			 "-auto_qc" => 1,
			);
  foreach my $k (sort {$intervals{$b} <=> $intervals{$a}} keys %intervals) {
    my %r;
    $r{range} = $k;
    $r{count} = $intervals{$k};
    $r{unique_junctions} = scalar keys %{$intervals2sites{$k}};

    my @ranges = split /$JOIN_TAG/, $k;
    my @info;
    my @feat;
    foreach my $range (@ranges) {
      my ($chr, $interval) = split /:/, $range;
      my ($start, $end) = split /-/, $interval;
      push @info, $chr, $start, $end;

      my $feature_info = join ",", sort keys %{$interval2feature{$range}};
      push @feat, $feature_info;
    }
    @r{qw(chrA startA endA chrB startB endB)} = @info;
    @r{qw(featureA featureB)} = @feat;
    $r{geneA} = $gene_a;
    $r{geneB} = $gene_b;
    $rpt->end_row(\%r);
  }
  $rpt->finish();


}

sub find_interval_and_feature {
  my ($rff, $rf, $pos) = @_;
  my $exons = $rf->{exons} || die;
  my $interval;
  for (my $i = 0; $i < @{$exons}; $i++) {
    my ($start, $end) = @{$exons->[$i]}{qw(start end)};
    if ($pos >= $start and $pos <= $end) {
      $interval = [ $start, $end ];
      last;
    }

    unless ($interval) {
      my $i_next = $i + 1;
      if ($i_next < @{$exons}) {
	my $i_start = $exons->[$i]->{end} + 1;
	my $i_end = $exons->[$i_next]->{start} - 1;
	if ($pos >= $i_start and $pos <= $i_end) {
	  # intronic
	  $interval = [ $i_start, $i_end ];
	  last;
	}
      }
    }
  }

  my $feature = get_feature_tag($rff, $rf, $pos);

  return ($interval, $feature);
}

sub generate_breakpoint_pairs {
  my $infile = $FLAGS{"generate-breakpoint-pairs"} || die;
  my $outfile = basename($infile) . ".spec.tab";
  my $step = $FLAGS{step} || die "-step";

  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );
  my $rpt = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => [
				       "geneA",
				       "chrA",
				       "posA",

				       "geneB",
				       "chrB",
				       "posB",
				      ],
			 "-auto_qc" => 1,
			);

  while (my $row = $df->get_hash()) {
    my $chr_a = $row->{chrA} || die;
    my $chr_b = $row->{chrB} || die;

    my $start_a = $row->{startA} || die;
    my $start_b = $row->{startB} || die;
    my $end_a = $row->{endA} || die;
    my $end_b = $row->{endB} || die;

    my %r = %{$row};
    $r{chrA} = $chr_a;
    $r{chrB} = $chr_b;
    for (my $ap = $start_a; $ap <= $end_a; $ap += $step) {
      $r{posA} = $ap;
      for (my $bp = $start_b; $bp <= $end_b; $bp += $step) {
	$r{posB} = $bp;
	$rpt->end_row(\%r);
      }
    }
  }
  $rpt->finish();
}

sub extract_svs {
  # extract SV records matching a set of intervals
  my $src_files = read_simple_file($FLAGS{"extract-svs"} || die);
  my $intervals = $FLAGS{intervals} || die;

  my $df = new DelimitedFile(
			     "-file" => $intervals,
			     "-headers" => 1,
			    );

  my %wanted_chr_a;
  while (my $row = $df->get_hash()) {
    my $chrA = $row->{chrA} = clean_chr($row->{chrA} || die);
    $row->{chrB} = clean_chr($row->{chrB} || die);
    push @{$wanted_chr_a{$chrA}}, $row;
  }

  foreach my $src_file (@{$src_files}) {
    my $outfile = basename($src_file) . ".excerpt.tab";
    my $df = new DelimitedFile("-file" => $src_file,
			       "-headers" => 1,
			      );
    my $rpt = $df->get_reporter(
				"-file" => $outfile,
				"-auto_qc" => 1,
			       );

    while (my $row = $df->get_hash()) {
      my $chrA = clean_chr($row->{chrA} || die);
      if (my $intervals = $wanted_chr_a{$chrA}) {
	my $posA = $row->{posA} || die;
	my $chrB = clean_chr($row->{chrB} || die);
	my $posB = $row->{posB} || die;
	foreach my $interval_row (@{$intervals}) {
	  if ($interval_row->{chrB} eq $chrB and
	      $posA >= $interval_row->{startA} and
	      $posA <= $interval_row->{endA} and
	      $posB >= $interval_row->{startB} and
	      $posB <= $interval_row->{endB}) {
	    # source row within a target interval
	    $rpt->end_row($row);
	    last;
	  }
	}
      }
    }

    $rpt->finish();
  }

}

sub condense_from_report {
  my ($report) = @_;

  my $glob = $FLAGS{"condense-pattern-annot-glob"} || "*.pattern.tab";
  my @pattern_annot_files = glob($glob);

  unless (@pattern_annot_files) {
    die "no .pattern.tab files" unless $FLAGS{"condense-no-pattern-files"};
  }

  my $df = new DelimitedFile("-file" => $report,
			     "-headers" => 1,
			     );
  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  my %non_annot = map {$_, 1} qw(
				pattern
				sequence
				keep
				subsumed_by
				reason
			       );

  my %f_annot = map {$_, 1} grep {!$non_annot{$_}} keys %{$df->headers};

  my %keep;
  my %subsumed;
  my %subsumed_by;
  my %keep2old;
  my %pid_all;

  while (my $row = $df->get_hash()) {
    my $pid = $row->{pattern} || die;
    $pid_all{$pid} = 1;
    my $keep = $row->{keep};
    if ($keep) {
      $keep{$pid} = 1;
    } else {
      my $pid_survivor = $row->{subsumed_by} || die;
      $subsumed{$pid} = $pid_survivor;
      $subsumed_by{$pid_survivor} = 1;
      push @{$keep2old{$pid_survivor}}, $row;
    }
  }

  foreach my $pid (keys %subsumed_by) {
    if (my $pid_survivor = $subsumed{$pid}) {
      # a pattern is subsumed by another, but the survivor is itself
      # subsumed.  e.g. some patterns are subsumed by perfect matches,
      # but the survivor can then be subsumed by a blast match
      printf STDERR "subsumed_by $pid itself subsumed by $pid_survivor\n";
      my $old_set = $keep2old{$pid} || die;
      push @{$keep2old{$pid_survivor}}, @{$old_set};

      die "fail" if $subsumed{$pid_survivor};
      # not sure if current code is good enough: does this ever happen?
      # if still fails, recurse?
    }
  }

  foreach my $pid (sort keys %pid_all) {
    if ($keep{$pid}) {
      # OK
    } elsif ($subsumed{$pid}) {
      # OK
    } else {
      die "no tracking info for $pid";
      # sanity check
    }
  }

  my $pattern_file = $report;
  $pattern_file =~ s/\.condense_report.tab$// || die;
  # source patterns to condense

  $df = new DelimitedFile(
			  "-file" => $pattern_file,
			  "-headers" => 1,
			 );
  my $outfile = $pattern_file . ".condensed.tab";
  my $rpt = $df->get_reporter("-file" => $outfile);

  while (my $row = $df->get_hash()) {
    # filter pattern file, keeping only surviving IDs
    my $pid = $row->{pattern} || die;
    if ($keep{$pid}) {
      if (my $obs_rows = $keep2old{$pid}) {
	# merge annotations from subsumed rows
	printf STDERR "merge %d rows to %s\n", scalar(@{$obs_rows}), $pid;
	merge_annots($row, $obs_rows, \%f_annot);
      }
      $rpt->end_row($row);
    }
  }
  $rpt->finish();

  # update .pattern.tab annotations of source files:
  foreach my $pf (@pattern_annot_files) {
    my $df = new DelimitedFile(
			    "-file" => $pf,
			    "-headers" => 1,
			   );
    my $outfile = $pf . ".condensed_pattern.tab";
    # use a suffix other than .pattern.tab to distinguish from source files
    my $rpt = $df->get_reporter("-file" => $outfile);
    while (my $row = $df->get_hash()) {
      if (my $pid = $row->{pattern}) {
	if ($keep{$pid}) {
#	  die "keep $pid";
	} elsif (my $pid_survivor = $subsumed{$pid}) {
	  my $pid_survivor = resolve_subsume($pid, \%subsumed);
	  $row->{pattern} = $pid_survivor;
	} elsif ($pid eq "suppressed_by_blacklist") {
	  # ignore
	} else {
	  # encountered unhandled ID, ???
	  # unrelated file?
	  die "no record of pattern $pid";
	}
      }
      $rpt->end_row($row);
    }
    $rpt->finish();
  }

}

sub merge_annots {
  my ($row, $obs_rows, $f_annot) = @_;
  foreach my $f (keys %{$f_annot}) {
    my %v;
    foreach my $r ($row, @{$obs_rows}) {
      my $v = $r->{$f};
      if (defined $v) {
	foreach (split /,/, $v) {
	  $v{$_} = 1;
	}
      }
    }

    my $new = join ",", sort keys %v;
    printf STDERR "%s: old=%s new=%s\n", $f, $row->{$f}, $new;
    $row->{$f} = $new;
  }
}

sub transcriptome_annotate {
  #
  # annotate regions in patterns:
  #   - that match other genes in the transcriptome
  #   - that are low complexity (via dustmasker)
  #

  printf STDERR "ignore matches to:\n";
  printf STDERR "  non-coding sequences?: %s\n", $TRANSCRIPTOME_AMBIG_IGNORE_NON_CODING ? "y" : "n";
  printf STDERR "  read-through sequences?: %s\n", $TRANSCRIPTOME_AMBIG_IGNORE_READTHROUGH ? "y" : "n";

  my $f_patterns = $FLAGS{"transcriptome-annotate"} || die;

  my $blast = get_blast();
  $blast->blast_2_sequences(0);
  # subject is transcriptome (indexed database)

  $blast->dust("yes");
  # re-enable BLAST's default dust filtering so it can run faster and
  # not trigger alignments to low-complexity regions (those will be
  # annotated separately).

  my $df = new DelimitedFile("-file" => $f_patterns,
			     "-headers" => 1,
			    );

  my @ambig_headers = qw(
			  ambig_A_dust
			  ambig_A_blast
			  ambig_A_blast_gene_count
			  ambig_A_blast_genes
			  ambig_A_blast_detail
			  ambig_A_coverage
			  ambig_A_notes

			  ambig_B_dust
			  ambig_B_blast
			  ambig_B_blast_gene_count
			  ambig_B_blast_genes
			  ambig_B_blast_detail
			  ambig_B_coverage
			  ambig_B_notes

			  ambig_any

			  ambig_homologous_gene_pair
			  ambig_warning
		       );

  my $h_needed = detect_needed_headers($df, \@ambig_headers);

  my $outfile = basename($f_patterns) . ".ambig.tab";
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => $h_needed,
  			      "-auto_qc" => 1,
			     );

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
    my $geneA = $row->{genea_symbol} || die;
    my $geneB = $row->{geneb_symbol} || die;
    my $pattern = $row->{sequence} || die;

    $pattern =~ /^(.*)[\]\}].*[\[\{](.*)$/ || die;

    my ($genea_seq, $geneb_seq) = ($1, $2);
    my $idx_b = index($pattern, $geneb_seq);

    my %regions_genes;
    my %regions_dust;
    my @notes;

    annotate_ambiguity(
		       "-row" => $row,
		       "-blast" => $blast,
		       "-end" => "A",
		       "-gene" => $geneA,
		       "-sequence" => $genea_seq,
		       "-offset" => 0);

    annotate_ambiguity(
		       "-row" => $row,
		       "-blast" => $blast,
		       "-end" => "B",
		       "-gene" => $geneB,
		       "-sequence" => $geneb_seq,
		       "-offset" => $idx_b);

    # add summary ambiguity annotation:
    my @intervals = grep {$_} @{$row}{qw(
					  ambig_A_dust
					  ambig_A_blast
					  ambig_B_dust
					  ambig_B_blast
				       )};
    my $all = "";
    if (@intervals) {
      printf STDERR "merged ci %s\n", join " / ", @intervals;
      my $si = new Set::IntSpan(@intervals);
      $all = $si->run_list();
    }
    $row->{ambig_any} = $all;

    my $warning = "";
    foreach my $f (qw(ambig_A_coverage ambig_B_coverage)) {
      my $ambig = $row->{$f};
      $warning = 1 if $ambig and $ambig >= $AMBIG_FRACTION_WARNING;
    }
    $row->{ambig_warning} = $warning;

    #
    #  cursory check for potentially homologous gene pairing, which
    #  may indicate fundamentally unreliable patterns
    #
    my $homologous = "";
    if ($geneA ne $geneB and length($geneA) == length($geneB)) {
      # hack, what if it's a larger number, e.g. SYM9 vs SYM10?

      my $len = length($geneA);

      my $mismatches = 0;
      for (my $i = 0; $i < $len; $i++) {
	$mismatches++ unless substr($geneA, $i, 1) eq substr($geneB, $i, 1);
      }
      $homologous = join ",", $geneA, $geneB if $mismatches == 1;
    }
    $row->{ambig_homologous_gene_pair} = $homologous;

    $rpt->end_row($row);
  } # $row
  $rpt->finish();
}

sub annotate_ambiguity {
  # find ambiguous features in sub-pattern
  my %options = @_;
  my $blast = $options{"-blast"} || die;
  my $row = $options{"-row"} || die;
  my $end = $options{"-end"} || die;
  my $pattern_gene = $options{"-gene"} || die;
  my $sequence = $options{"-sequence"} || die;
  my $index_offset = $options{"-offset"};
  die unless defined $index_offset;
  printf STDERR "annotating ambig for %s\n", $end;

  my $f_transcriptome = $FLAGS{"transcriptome-fasta"} || die;

  my %regions_dust;
  my %regions_genes;
  my @notes;

  #
  #  track low-complexity regions w/dustmasker:
  #
  if (my $hits = dustmasker("-sequence" => $sequence,
			    "-as-string" => 1,
			    "-offset" => $index_offset
			   )) {
    # TO DO: require minimum length?  Sometimes these are tiny
    foreach my $interval (@{$hits}) {
      $regions_dust{$interval} = 1;
    }
  }

  #
  #  BLAST vs. transcriptome:
  #
  my $parser = $blast->blast(
			     "-query" => {
					  $pattern_gene => $sequence
					 },
			     "-database" => $f_transcriptome,
			    );
  my $result = $parser->next_result();
  # one object per query sequence (only one query seq)
  my $VERBOSE = 1;

  if ($result) {
    my %hits;
    my %map;

    while (my $hit = $result->next_hit()) {
      # hits from this query to a database sequence
      # (Bio::Search::Hit::HitI)
      my $hit_name = $hit->name();

      next if $TRANSCRIPTOME_AMBIG_IGNORE_NON_CODING and $hit_name =~ /[NX]R_\d+/;

      my $hit_desc = $hit->description();
      my $hit_gene;
      if ($hit_desc =~ /\/gene=(\S+)/) {
	# FASTA built by refgene2protein.pl
	$hit_gene = $1;
      } else {
	my @paren_genes = $hit_desc =~ /\((\S+)\),/g;
	# NCBI refSeq defline
	# require trailing comma as this always seems to be present.

	# in refSeq deflines, gene symbol is typically in parens.
	# however sometimes other strings appear too, e.g.:
	# ref|NM_018702.4| / Homo sapiens adenosine deaminase RNA specific B2 (inactive) (ADARB2), mRNA.
	# Homo sapiens sphingolipid transporter 3 (putative) (SPNS3), transcript variant 2, mRNA.
	# Homo sapiens DnaJ heat shock protein family (Hsp40) member C14 (DNAJC14), transcript variant 1, mRNA.

	# can't require 2+ capital letters to avoid Hsp40, because
	# this fails to match e.g. C17orf113.

	if (@paren_genes > 1) {
	  # Homo sapiens 3'(2'), 5'-bisphosphate nucleotidase 2 (BPNT2), mRNA.
	  @paren_genes = grep {/^[A-Z]/} @paren_genes;
	}

	dump_die($row, sprintf "gene sym lookup failure in %s: %s", $hit_desc, join ",", @paren_genes) unless @paren_genes == 1;
	$hit_gene = $paren_genes[0];
      }

      if ($TRANSCRIPTOME_AMBIG_IGNORE_READTHROUGH and
	  $hit_gene =~ /\-/ and
	  $hit_desc =~ /readthrough/) {
	# e.g. NM_001199973.2
	# RPL36A-HNRNPH2,Homo sapiens RPL36A-HNRNPH2 readthrough (RPL36A-HNRNPH2), transcript variant 1, mRNA
	next;
      }


      if ($hit_desc =~ /\/map=(\S+)/) {
	$map{$1}++;
      } else {
	# some refSeq entries don't contain this information
	$map{no_map_info}++;
	# not sure how to deal with these, just ignore?
	# doesn't seem to happen that often, and when it does there
	# are always other cyto entries to bolster ambiguity case.
#	dump_die($row, "no map info in $hit_desc");
      }


      #	  next if $hit_gene eq $gene;
      # save all genes, even the expected one, for GSM purposes later

      printf STDERR "hit to %s / %s\n", $hit_name, $hit->description;

      #      die $hit->ambiguous_aln;
      # no help for FLT3/FLT3

      while (my $hsp = $hit->next_hsp) {
	# Bio::Search::HSP::PSLHSP
	printf STDERR "    score:%s strand:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s query_span:%d ref_span:%d total_span=%d query_string=%s hit_string=%s\n",
	  $hsp->score,
	    $hsp->strand,
	      $hsp->range("query"),
		$hsp->range("hit"),
		  $hsp->num_identical(),
		    $hsp->frac_identical("query"),
		      $hsp->length("query"),
			$hsp->length("hit"),
			  $hsp->length("total"),
			    $hsp->query_string(),
			      $hsp->hit_string() if $VERBOSE;

	my $qs = $hsp->start("query") + $index_offset;
	my $qe = $hsp->end("query") + $index_offset;

	my $interval = join "-", $qs, $qe;
	$hits{$hit_gene}{$interval} = 1;
	# TO DO: add minimum hit requirements here

      }				# $hsp
    }				# $hit

    #
    #  prune gene matches:
    #
    my $found_main;
    if ($hits{$pattern_gene}) {
      # don't report hit to expected gene
      $found_main = 1;
      delete $hits{$pattern_gene};
    } else {
      my $gsm = new_gsm_lite("-hgnc" => $FLAGS{"hgnc"});
      foreach my $gene (keys %hits) {
	$gsm->add_gene("-gene" => $gene);
      }
      if (my $gene_rf = $gsm->find($pattern_gene)) {
	# found a match to the pattern gene using a different symbol
	push @notes, sprintf 'gene_equiv=%s', $gene_rf;
	delete $hits{$gene_rf};
	$found_main = 1;
      }
    }

    # TO DO: prune GSM hit

    if (not($found_main) and scalar keys %hits == 1) {
      # GSM may fail if genenames.org database isn't up to date.
      # however, if there is only one gene symbol observed in the
      # results then by definition there is no ambiguity.
      my @found = keys %hits;
      push @notes, sprintf "single_gene_hit_but_no_symbol_match=%s", @found;
#    } elsif (%hits and scalar keys %map == 1) {
#      push @notes, sprintf "all_map_to_same_region=%s", join "/", (keys %map), sort keys %hits if scalar keys %map == 1;
    } elsif (%hits) {
      # merge final/pruned hits into tracker
      if (%map) {
	push @notes, sprintf 'cyto_matches=%s', join "/", sort keys %map;
	push @notes, sprintf "all_map_to_same_region=%s", join "/", (keys %map), sort keys %hits if scalar keys %map == 1;
      }
      foreach my $gene (keys %hits) {
	foreach my $interval (keys %{$hits{$gene}}) {
	  $regions_genes{$gene}{$interval} = 1;
	}
      }
    }
  }				# $result

  my $f_dust = sprintf 'ambig_%s_dust', $end;
  my $f_blast = sprintf 'ambig_%s_blast', $end;
  my $f_blast_genes = sprintf 'ambig_%s_blast_genes', $end;
  my $f_blast_gene_count = sprintf 'ambig_%s_blast_gene_count', $end;
  my $f_blast_detail = sprintf 'ambig_%s_blast_detail', $end;
  my $f_notes = sprintf 'ambig_%s_notes', $end;
  my $f_coverage = sprintf 'ambig_%s_coverage', $end;

  $row->{$f_notes} = join ",", @notes;

  $row->{$f_dust} = build_run_list(\%regions_dust);
  # dust intervals in ordered run-list format

  my %intervals_blast;
  foreach my $gene (sort keys %regions_genes) {
    my @intervals = sort keys %{$regions_genes{$gene}};
    foreach my $interval (@intervals) {
      $intervals_blast{$interval}{$gene} = 1;
#      push @other_hit, sprintf "%s=%s", $interval, $gene;
    }
  }

  $row->{$f_blast} = build_run_list(\%intervals_blast);
  # summary of BLAST intervals with intervals combined

  $row->{$f_blast_genes} = join ",", sort keys %regions_genes;
  $row->{$f_blast_gene_count} = scalar keys %regions_genes;

  # detailed BLAST hits including gene
  my @bd;
  foreach my $interval (sort keys %intervals_blast) {
    my @syms = sort keys %{$intervals_blast{$interval}};
    push @bd, sprintf '%s=%s', $interval, join "/", @syms;
  }
  $row->{$f_blast_detail} = join ",", @bd;

  # total percentage flagged by both BLAST and dust:
  my $covg = "";
  if (%intervals_blast or %regions_dust) {
    my $si = new Set::IntSpan(keys %intervals_blast, keys %regions_dust);
    my $covered = cardinality $si;
    my $plen = length $sequence;
    $covg = sprintf "%.2f", $covered / $plen;
  }
  $row->{$f_coverage} = $covg;
}

sub build_run_list {
  my ($hashref) = @_;
  my $rl = "";
  if (%{$hashref}) {
    my $si = new Set::IntSpan(keys %{$hashref});
    $rl = $si->run_list;
  }
  return $rl;
}

sub duplicate_window_trim {
  my ($seq) = @_;

  my $i1 = get_marker_idx(\$seq, "]", "}");

  # trim excess sequence before the first breakpoint marker:
  my $excess = $i1 - $PATTERN_CONDENSE_DUPLICATE_FLANKING_WINDOW_SIZE;
  $seq = substr($seq, $excess) if $excess > 0;

  # excess after second marker:
  my $i2 = get_marker_idx(\$seq, "[", "{");
  $seq = substr($seq, 0,
		$i2 + 1 + $PATTERN_CONDENSE_DUPLICATE_FLANKING_WINDOW_SIZE);
  return $seq;
}

sub get_marker_idx {
  my ($seq_ref, @markers) = @_;
  my $i;
  foreach my $marker (@markers) {
    $i = index($$seq_ref, $marker);
    last if $i >= 0;
  }
  die "breakpoint marker not found in $$seq_ref" if not(defined $i) or $i == -1;
  return $i;
}

sub get_max_gap_size {
  my ($hsps) = @_;

  my $max_gap_size = 0;

  foreach my $hsp (@{$hsps}) {
    foreach my $aln ($hsp->query_string(),
		     $hsp->hit_string()) {
      my @hits = $aln =~ /\-+/g;
      foreach my $gap (@hits) {
	my $len = length $gap;
	$max_gap_size = $len if $len > $max_gap_size;
      }
    }
  }
  return $max_gap_size;
}

sub parse_cosmic {
  # identify genomic breakpoints in COSMIC genome-based intervals
  my $f_in = $FLAGS{"parse-cosmic"} || die;
  my $f_out = basename($f_in) . ".breakpoint.tab";

  die "columns not defined, use e.g. -cicero to specify format\n" unless $F_CHR_A;


  my $df = new DelimitedFile("-file" => $f_in,
			     "-headers" => 1,
			     );
  my $rpt = $df->get_reporter(
			      "-file" => $f_out,
			      "-extra" => [
					   $F_SAMPLE,

					   $F_GENE_A,
					   $F_CHR_A,
					   $F_POS_A,
					   $F_ORT_A,

					   $F_GENE_B,
					   $F_CHR_B,
					   $F_POS_B,
					   $F_ORT_B,
					  ],
  			      "-auto_qc" => 1,
			     );

  while (my $row = $df->get_hash()) {
#    dump_die($row, "start", 1);
    my $five_strand = $row->{"5'_STRAND"} || next;
    # not all rows have these annotations
    my $three_strand = $row->{"3'_STRAND"} || die;

    $row->{$F_ORT_A} = $five_strand;
    $row->{$F_ORT_B} = $three_strand;

    $row->{sample} = $row->{SAMPLE_NAME} || die;
    $row->{$F_CHR_A} = $row->{"5'_CHROMOSOME"} || die;
    $row->{$F_CHR_B} = $row->{"3'_CHROMOSOME"} || die;

    $row->{$F_GENE_A} = $row->{"5'_GENE_NAME"} || die;
    $row->{$F_GENE_B} = $row->{"3'_GENE_NAME"} || die;

    foreach my $f ($F_GENE_A, $F_GENE_B) {
      $row->{$f} =~ s/_ENST.*//;
      # remove transcript ID from gene name
      die $row->{$f} if $row->{$f} =~ /\W/;
    }

    my $five_start_from = $row->{"5'_GENOME_START_FROM"} || die;
    my $five_start_to = $row->{"5'_GENOME_START_TO"} || die;
    my $five_stop_from = $row->{"5'_GENOME_STOP_FROM"} || die;
    my $five_stop_to = $row->{"5'_GENOME_STOP_TO"} || die;
#    dump_die($row, "5' start diff") if $five_start_from != $five_start_to;
#    dump_die($row, "5' stop diff") if $five_stop_from != $five_stop_to;
    # from the COSMIC documentation:
    #
    # [17:Q] 5'_GENOME_START_FROM - The genomic coordinate of the start (+ strand)/breakpoint (- strand) of the 5â fusion gene as described in the Translocation Name.
    #
    # [18:R] 5'_GENOME_START_TO - The range of genomic coordinates of the start (+ strand)/breakpoint (- strand) of the 5â fusion gene if it is an unknown base position.
    #
    #  - this code only processes the "from" coordinates,
    #    i.e. those matching the transcript description.
    #  - however, for records where the "to" is different, i.e. an
    #    unknown position, does this indicate the record is unreliable?

    # see also FUSION_TYPE field:
    # - only use "Observed mRNA"?
    # - what about "Inferred Breakpoint"?

    my $three_start_from = $row->{"3'_GENOME_START_FROM"} || die;
    my $three_start_to = $row->{"3'_GENOME_START_TO"} || die;
    my $three_stop_from = $row->{"3'_GENOME_STOP_FROM"} || die;
    my $three_stop_to = $row->{"3'_GENOME_STOP_TO"} || die;
#    dump_die($row, "3' start diff") if $three_start_from != $three_start_to;
#    dump_die($row, "3' stop diff") if $three_stop_from != $three_stop_to;

    my ($pos_a, $pos_b);
    # breakpoints

    #
    #  find breakpoint for geneA (end of 5' of fusion):
    #
    if ($five_strand eq "+") {
      # for + genes, breakpoint is at stop of 5' genomic interval.
      # e.g. BCR-ABL1:
      #
      # 5'_CHROMOSOME: 22
      # 5'_GENE_ID: 120293
      # 5'_GENE_NAME: BCR
      # 5'_GENOME_START_FROM: 23522397
      # 5'_GENOME_START_TO: 23522397
      # 5'_GENOME_STOP_FROM: 23654023     --> breakpoint
      # 5'_GENOME_STOP_TO: 23654023
      # 5'_LAST_OBSERVED_EXON: 19
      # 5'_STRAND: +
      $pos_a = $five_stop_to;
    } elsif ($five_strand eq "-") {
      # for - genes, gene model begins at genomic 3' end and
      # continues through breakpoint at genomic 5' end, e.g.:
      #
      # 5'_CHROMOSOME: 8
      # 5'_GENE_ID: 816
      # 5'_GENE_NAME: RGS22
      # 5'_GENOME_START_FROM: 100994165   ----> exon 20
      # 5'_GENOME_START_TO: 100994165
      # 5'_GENOME_STOP_FROM: 101118344    ----> upstream of TSS
      # 5'_GENOME_STOP_TO: 101118344
      # 5'_LAST_OBSERVED_EXON: 22
      # 5'_STRAND: -
      $pos_a = $five_start_from;
    } else {
      die "unhandled strand $five_strand";
    }

    #
    #  find breakpoint for geneB (start of 3' of fusion):
    #
    if ($three_strand eq "+") {
      # 3'_CHROMOSOME: 1
      # 3'_FIRST_OBSERVED_EXON: 24
      # 3'_GENE_ID: 15541
      # 3'_GENE_NAME: SYCP1_ENST00000369518
      # 3'_GENOME_START_FROM: 115486960   ----> exon 24
      # 3'_GENOME_START_TO: 115486960
      # 3'_GENOME_STOP_FROM: 115537988    ----> exon 31
      # 3'_GENOME_STOP_TO: 115537988
      # 3'_STRAND: +
      #
      # BCR-ABL1 ABL1:
      # 3'_CHROMOSOME: 9
      # 3'_FIRST_OBSERVED_EXON: 1
      # 3'_GENE_ID: 85234
      # 3'_GENE_NAME: ABL1_ENST00000318560
      # 3'_GENOME_START_FROM: 133729451
      # 3'_GENOME_START_TO: 133729451
      # 3'_GENOME_STOP_FROM: 133763062
      # 3'_GENOME_STOP_TO: 133763062
      # 3'_STRAND: +
      $pos_b = $three_start_from;
    } elsif ($three_strand eq "-") {
      # 3'_CHROMOSOME: 3
      # 3'_FIRST_OBSERVED_EXON: 2
      # 3'_GENE_ID: 164142
      # 3'_GENE_NAME: ETV5
      # 3'_GENOME_START_FROM: 185764097   ----> past exon 13 (end)
      # 3'_GENOME_START_TO: 185764097
      # 3'_GENOME_STOP_FROM: 185823731    ----> exon 2
      # 3'_GENOME_STOP_TO: 185823731
      # 3'_STRAND: -
      $pos_b = $three_stop_from;
    } else {
      die "unhandled strand $three_strand";
    }

    $row->{$F_POS_A} = $pos_a;
    $row->{$F_POS_B} = $pos_b;

#   dump_die($row, "$pos_a $pos_b");
    $rpt->end_row($row);
  }

  $rpt->finish();
}

sub generate_pair_summary {
  my ($r) = @_;
  my $summary_a = join "/", @{$r}{qw(genea_symbol genea_acc genea_feature)};
  my $summary_b = join "/", @{$r}{qw(geneb_symbol geneb_acc geneb_feature)};
  $r->{$F_GENE_PAIR_SUMMARY} = join "-", $summary_a, $summary_b;

  my @check = $F_GENOME;
  push @check, (
		"genea_chr",
		"genea_pos",
		"geneb_chr",
		"geneb_pos"
	       ) unless $FLAGS{"no-chr-pos"};

  foreach my $f (@check) {
    dump_die($r, "no field $f") unless $r->{$f};
  }

  my $ext_a = join "/", @{$r}{qw(genea_symbol genea_acc genea_feature),
				$F_GENOME, "genea_chr", "genea_pos"};
  my $ext_b = join "/", @{$r}{qw(geneb_symbol geneb_acc geneb_feature),
				$F_GENOME, "geneb_chr", "geneb_pos"};
  $r->{$F_GENE_PAIR_SUMMARY_EXTENDED} = join "-", $ext_a, $ext_b;
}

sub extract_cicero {
  #
  #  extract CICERO calls of interest
  #
  my $f_report = $FLAGS{"extract-cicero"} || die;
  my $f_idx = $FLAGS{"steve-index"} || die "-steve-index";
  my $df = new DelimitedFile("-file" => $f_report,
			     "-headers" => 1,
			     );

  # samples and fusions to search for:
  my %sample2fusion;
  while (my $row = $df->get_hash()) {
    my $sample = $row->{Sample} || die;
#    my $fusion = $row->{Fusion} || die;
    my $fusion = unquote($row->{Fusion} || die);
    $sample2fusion{$sample}{$fusion} = 1;
  }

  my %prefix2path = (
		     "Clinical" => "/research/rgs01/reference/restricted/ClinicalGenomics/clingen_curated/Clinical",

		     "CuratedPilot" => "/research/rgs01/reference/restricted/ClinicalGenomics/clingen_curated/ClinicalPilot/ClinicalPilot",

		     "ClinicalPilot" => "/research/rgs01/project_space/zhanggrp/ClinicalPilot/common/cicero-post"
		    );
  # translation of prefixes in index file to true locations

  # load index of file locations:
  my %sample2file;
  open(IDX, $f_idx) || die;
  while (<IDX>) {
    chomp;
    my @f = split /\t/, $_;
    die unless @f == 5;
    my ($sample, $file_raw, $file_reviewed) = @f[0,3,4];

    my $file = $file_reviewed;
    $file = $file_raw if $file_reviewed eq "none";
    # sometimes there isn't a reviewed file

    $file =~ s/^(\w+)// || die "can't find prefix in $file";
    my $prefix = $1;

    my $root = $prefix2path{$prefix} || die "can't resolve path code for $prefix in $file";
    $file = $root . $file;
    die "where is $file" unless -s $file;

    die "duplicate" if $sample2file{$sample};
    $sample2file{$sample} = $file;
  }

  my $f_out = basename($f_report) . ".extracted.tab";

  my $f_gene_a = "gene_a";
  my $f_gene_b = "gene_b";
  # use some "do what I mean" code to convert the header we want
  # (CICERO) into what we actually see in the file?

  my $dwim = get_column_dwim();
  $F_FUSION = "fusion_gene";
  $dwim->add_aliases($F_FUSION);

  my $F_REVIEW = "manual_review_status";
  $dwim->add_aliases($F_REVIEW);

  $dwim->add_aliases($F_SAMPLE, "case_sample");

  my @f_out = (
	       $F_REVIEW,
	       $F_SAMPLE,
	       $F_GENE_A,
	       $F_CHR_A,
	       $F_POS_A,
	       $F_ORT_A,
	       $F_GENE_B,
	       $F_CHR_B,
	       $F_POS_B,
	       $F_ORT_B,
	       $F_CONTIG,
	       $F_FUSION,
	       "source_file",
	       "notes",
	      );

  my $rpt = new Reporter(
			 "-file" => $f_out,
			 "-delimiter" => "\t",
			 "-labels" => \@f_out,
			 "-auto_qc" => 1,
			);

  my $errors = 0;

  foreach my $sample (sort keys %sample2fusion) {
    my %wanted = map {$_, 1} keys %{$sample2fusion{$sample}};
    my $report = $sample2file{$sample} || die "no file for $sample";

    my $df = new DelimitedFile("-file" => $report,
			       "-headers" => 1,
			      );
    # ugh: these files can have different layouts, so can't simply
    # pass through all rows

    $dwim->init_field_map($df->headers);

    my %found;

    my $review_complain;

    while (my $row = $df->get_hash()) {
      my $gene_a = $dwim->get_data($row, $F_GENE_A);
      my $gene_b = $dwim->get_data($row, $F_GENE_B);
      my $tag = join "-", $gene_a, $gene_b;

      if ($wanted{$tag}) {
	$row->{source_file} = $report;
	$dwim->populate_fields($row);
	my @notes;

	my $fusion = $row->{$F_FUSION} || die;
	# QC check: is fusion_gene
	# tab_dump /research/rgs01/reference/restricted/ClinicalGenomics/clingen_curated/Clinical/2019/SJBALL031555/SJBALL031555_D1_G1/multi_target/manual-review/SJBALL031555_D1_G1_somatic_sv.txt -match CEBPE
	# CEBPE is actually ACIN1, which is what fusion_gene indicates, ???
	my $gene_a = $row->{$F_GENE_A} || die;
	my $gene_b = $row->{$F_GENE_B} || die;
	my $ok;

	foreach my $gene_raw ($gene_a, $gene_b) {
	  my @genes = split /,/, $gene_raw;
	  foreach my $gene (@genes) {
	    if ($fusion =~ /$gene/) {
	      $ok = 1;
	      last;
	    }
	  }
	}

	unless ($ok) {
	  push @notes, sprintf "WARNING_genes_disagree_with_fusion=A:%s,B:%s,fusion:%s", $gene_a, $gene_b, $fusion;
	}

	$row->{notes} = join ",", @notes;

	unless ($row->{$F_REVIEW}) {
	  $row->{$F_REVIEW} = "n/a";
	  unless ($review_complain) {
	    # not always available
	    $review_complain = 1;
	    dump_die($row, "can't find manual review call in $report", 1);
	  }
	}

	$rpt->end_row($row);
	$found{$tag} = 1;
	# don't simply delete from %wanted at this point in case
	# there are multiple entries in the file for this event
      }
    }

    delete @wanted{keys %found};

    if (%wanted) {
      dump_die(\%wanted, "couldn't find matches", 1);
      $errors++;
    }
  }

  # TO DO: standardize headers??
  # something that can handle ALL of the many variations?

  if ($errors) {
    printf STDERR "%d errors, not writing output\n", $errors;
  } else {
    $rpt->finish();
  }
}


sub get_column_dwim {
  #
  #  alternate column names for CICERO headers
  #
  die "chr/pos/ort fields not defined; suggest using -cicero, -crest, etc." unless ($F_CHR_A and $F_POS_A and $F_CHR_B and $F_POS_B);
  my $dwim = new DWIMColumns();
  $dwim->add_aliases($F_GENE_A, "gene_a");
  $dwim->add_aliases($F_GENE_B, "gene_b");
  $dwim->add_aliases($F_CHR_A, "chr_a");
  $dwim->add_aliases($F_CHR_B, "chr_b");
  $dwim->add_aliases($F_POS_A, "pos_a");
  $dwim->add_aliases($F_POS_B, "pos_b");
  $dwim->add_aliases($F_ORT_A, "ori_a");
  $dwim->add_aliases($F_ORT_B, "ori_b");
  $dwim->add_aliases($F_CONTIG);
  return $dwim;
}

sub get_supplemental_cache {
  return new GenBankFASTACache(
			       "-file" => $FLAGS{"supplemental-fasta"} || die,
			       "-hgnc" => get_hgnc()
			      );
}

sub add_supplemental_gene {
  my $cache = get_supplemental_cache();
  my $gene = $FLAGS{"add-supplemental-gene"} || die;
  my $acc = $FLAGS{accession};

  if ($acc) {
    # already know specific accession
    $cache->add_accession($acc, $gene);
  } else {
    # find a representative sequence for the gene via HGNC
    $cache->add_gene($FLAGS{"add-supplemental-gene"} || die);
  }
}

sub get_hgnc {
  my $f_hgnc = $FLAGS{"hgnc"};

  unless ($f_hgnc) {
    my $config_species = TdtConfig::readConfig('species', 'Homo_sapiens') || die "can't find species config";
    $f_hgnc = $config_species->{HGNC} || die;
  }
  die unless $f_hgnc and -s $f_hgnc;
  return $f_hgnc;
}

sub is_supplemental {
  my ($row) = @_;
  my $is_supplemental = 0;
  if ($ENABLE_SUPPLEMENTAL_REFSEQ) {
    if ($row->{cdsStart}) {
      # refFlat-style record
      $is_supplemental = 0;
    } elsif ($row->{seq_ref}) {
      # from GenBankFASTACache
      $is_supplemental = 1;
    }
  }
  return $is_supplemental;
}

sub get_supplemental_sequence {
  my ($sr, $contig_raw, $notes_ref, $end) = @_;
  my $genomic = ${$sr->{seq_ref}};

  my $acc = $sr->{accession} || die;
  my $need_check = $acc =~ /^NG_/;
  # TO DO: other options like force check

  if ($need_check) {
    # as a hint, just use roughly the portion of the contig corresponding
    # to the gene we're trying to anchor:
    my $half = int(length($contig_raw) / 2);
    my $contig;
    if ($end eq "A") {
      $contig = substr($contig_raw, 0, $half);
    } elsif ($end eq "B") {
      $contig = substr($contig_raw, $half);
    } else {
      die;
    }

    my $blast = get_blast();
    $blast->strand("");
    # remove strand filter, we want results for both strands

    my $hsps = run_blast(
			 "-blast" => $blast,
			 "-query" => {
				      "contig" => $contig
				     },
			 "-database" => {
					 $acc => $genomic
					},
			 "-notes" => [],
			 "-end" => "stand_check",
			 "-multi-ok" => 1,
			);
    my $count_plus = 0;
    my $count_minus = 0;
    if ($hsps) {
      foreach my $hsp (@{$hsps}) {
	my $strand = $hsp->strand("hit");
	if ($strand == 1) {
	  $count_plus++;
	} elsif ($strand == -1) {
	  $count_minus++;
	} else {
	  die;
	}
      }
    }
    printf STDERR "strand check for %s, plus=%d minus=%d\n", $acc, $count_plus, $count_minus;

    my $genomic_rc = reverse_complement($genomic);

    if ($count_plus and $count_minus) {
      # deal with this if/when we see it
      my $best_hsp_strand = $hsps->[0]->strand("hit");
      if ($best_hsp_strand == -1) {
	$genomic = $genomic_rc;
	push @{$notes_ref}, sprintf "reverse_complemented_supplemental_sequence_ambiguous=%s", $acc;
	# slightly different tag
      }
    } elsif ($count_plus) {
      # already on +, nothing to do
    } elsif ($count_minus) {
      # BLAST hits only to -, so use reverse-complemented sequence.
      # e.g. IGK uses a NG_000833, a genomic sequence which does not
      # capture the strand the gene is on
      $genomic = $genomic_rc;
      push @{$notes_ref}, sprintf "reverse_complemented_supplemental_sequence=%s", $acc;
    }
  }

  return $genomic;
}

sub deep_intronic_check {
  my ($row, $exon_distance, $end, $contig_ref, $notes_ref) = @_;
  if ($exon_distance) {
    my $contig_half_len = length($$contig_ref) / 2;
    my $cutoff = $contig_half_len * $DEEP_INTRONIC_MIN_CONTIG_SIDE_FRACTION;
    if ($exon_distance > $cutoff) {
      push @{$notes_ref}, sprintf '%s=%s', TAG_DEEP_INTRONIC, $end;
    }
  }
}

sub pattern2fasta {
  # TO DO: add policy for brackets: replace by Ns or just remove?
  my $f_in = $FLAGS{"pattern2fasta"} || die;
  my $id_wanted = $FLAGS{id};

  my $f_out = $FLAGS{out};

  unless ($f_out) {
    if ($id_wanted) {
      $f_out = $id_wanted . ".fa";
    } else {
      $f_out = basename($f_in) . ".fa";
    }
  }

  my $wf = new WorkingFile($f_out);
  my $fh = $wf->output_filehandle();

  my $df = new DelimitedFile("-file" => $f_in,
			     "-headers" => 1,
			     );
  while (my $row = $df->get_hash()) {
    my $id = $row->{pattern} || die;
    next if $id_wanted and $id ne $id_wanted;

    my $seq = $row->{sequence} || die;
    write_fasta_seq($fh, $id, $seq);
  }
  $wf->finish();
}

sub extract_pattern_sources {
  # extract source rows a pattern ID was generated from
  my $plist = read_simple_file($FLAGS{"check-bad-patterns"} || die);
  my %bad_pid = map {$_ => 1} @{$plist};

  my $df = new DelimitedFile("-file" => ($FLAGS{src} || die "-src"),
			     # *.pattern.tab
			     "-headers" => 1,
			     );
  my $outfile = "src.tab";
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			     );


  while (my $row = $df->get_hash()) {
    dump_die($row, "no pattern field") unless exists $row->{pattern};
    my $pid = $row->{pattern} || next;
    $rpt->end_row($row) if $bad_pid{$pid};
  }
  $rpt->finish();
}

sub extract_patterns {
  # extract pattern sequences for a pattern list
  my $thing = $FLAGS{"extract-patterns"} || die;
  my $plist;
  if (-s $thing) {
    $plist = read_simple_file($FLAGS{"extract-patterns"} || die);
  } else {
    $plist = [ split /,/, $thing ];
  }
  my %wanted = map {$_ => 1} @{$plist};
  my $substr_mode = $FLAGS{substr};
  my %found;

  my $df = new DelimitedFile(
			     "-file" => ($FLAGS{src} || die "-src"),
			     # pattern file (e.g. merged)
			     "-headers" => 1,
			     );

  my $outfile = $FLAGS{out} || "patterns.tab";

  my $rpt;
  if (1) {
    # keep all columns
    $rpt = $df->get_reporter(
			      "-file" => $outfile,
  			      "-auto_qc" => 1,
			     );
  } else {
    # minimal
    $rpt = new Reporter(
			   "-file" => $outfile,
			   "-delimiter" => "\t",
			   "-labels" => [
					 qw(
					     pattern
					     sequence
					  )
					],
			   "-auto_qc" => 1,
			  );
  }

  while (my $row = $df->get_hash()) {
    dump_die($row, "no pattern field") unless exists $row->{pattern};
    my $pid = $row->{pattern} || next;
    my $match;
    if ($substr_mode) {
      foreach my $p (keys %wanted) {
	if (index($pid, $p) >= 0) {
	  $match = 1;
	  $found{$p} = 1;
	  last;
	}
      }
    } else {
      if ($wanted{$pid}) {
	$match = 1;
	$found{$pid} = 1;
      }
    }
    $rpt->end_row($row) if $match;
  }
  $rpt->finish();

  foreach my $pid (sort keys %wanted) {
    unless ($found{$pid}) {
      printf STDERR "WARNING: can't find pattern for %s\n", $pid;
    }
  }
}

sub add_gene_annotation {
  my ($row, $ga, $gsm, $notes_ref) = @_;

#  my $keep_existing = $FLAGS{"keep-gene-annotation"};

  foreach my $end (qw(A B)) {
    my ($chr, $pos, $gene);
    my ($f_gene, $f_strand);
    if ($end eq "A") {
      $chr = $row->{$F_CHR_A} || die;
      $pos = $row->{$F_POS_A} || die;
      $f_gene = $F_GENE_A || die;
      $f_strand = $F_ORT_A || die;
    } elsif ($end eq "B") {
      $chr = $row->{$F_CHR_B} || die;
      $pos = $row->{$F_POS_B} || die;
      $f_gene = $F_GENE_B;
      $f_strand = $F_ORT_B || die;
    } else {
      die;
    }

    $gene = $row->{$f_gene};

    my $annotate_gene;
    my $annotate_strand;
    if ($RNA_KEEP_GENE_ANNOTATIONS) {
      if ($gene) {
	# retain existing annotation...
	if ($RNA_FIX_GENE_ANNOTATIONS) {
	  # ...unless it doesn't appear to be valid
	  $annotate_gene = 2 unless $gsm->find($gene);
	}
      } else {
	$annotate_gene = 1;
      }
      $annotate_strand = 1 unless $row->{$f_strand};
    } else {
      $annotate_gene = 1;
      $annotate_strand = 1;
    }


    if ($annotate_gene or $annotate_strand) {
      $ga->find(
		"-reference" => $chr,
		"-start" => $pos,
		"-end" => $pos
	       );
      my $genes = $ga->results_genes_genomic_order();
      if (@{$genes}) {

	if ($annotate_gene) {
	  #
	  #  add/update gene symbols
	  #
	  $row->{$f_gene} = join ",", @{$genes};
	  add_tag($notes_ref, TAG_ADDED_GENE_ANNOTATION, $end);
	  if ($annotate_gene == 2) {
	    add_tag($notes_ref, TAG_REPLACED_GENE_ANNOTATION, sprintf "%s>%s", $gene, join "/", @{$genes});
	  }
	}

	#
	#  add/update strand
	#
	if ($annotate_strand) {
	  my $strand_old = $row->{$f_strand} || "";
	  my $strand_new = $ga->results_strand() || "";
	  # some sites are ambiguous, e.g. chr1.23851720 hits both E2F2
	  # (-) and LOC101928163 (+).  In these cases, better to use a
	  # blank value and look it up on a gene-by-gene basis.

	  my $update_strand = 1;
	  # add/update by default

	  my $tag;

	  if ($strand_old and $strand_old ne $strand_new) {
	    # inconsisent strand annotation.
	    #
	    # in ProteinPaint export, there is evidence some of the strand
	    # values are bogus, e.g. ClinicalPilot SJOS001_M DLG2/GRM5,
	    # raw CICERO file shows ort - but export shows strand + (???)
	    # /clinical/cgs01/clingen/prod/tartan/index/data/ClinicalPilot/ClinicalPilot/SJOS001_M/TRANSCRIPTOME/cicero-post/SJOS001_M/final_fusions.counts
	    #
	    # OTOH, could be an intentional annotation.
	    if ($RNA_UPDATE_GENE_STRAND) {
	      $tag = TAG_INCONSISTENT_STRAND_REPLACED;
	    } else {
	      # mark inconsistency but don't change.
	      $tag = TAG_INCONSISTENT_STRAND_DETECTED;
	      $update_strand = 0;
	    }
	  } elsif (not($strand_old)) {
	    $tag = TAG_ADDED_STRAND;
	  }

	  $row->{$f_strand} = $strand_new if $update_strand;
	  add_tag($notes_ref, $tag, $end) if $tag;

	}

      } elsif ($annotate_gene == 2) {
	add_tag($notes_ref, TAG_UNRECOGNIZED_GENE_ANNOTATION, $gene);
      }
    }

  }  # $end

}

sub add_tag {
  my ($notes_ref, $tag, $value) = @_;
  push @{$notes_ref}, sprintf '%s=%s', $tag, $value;
}

sub set_coding_distance {
  my ($row, $end, $pos_raw, $pos_adj) = @_;
  my $cd = "";
  if ($pos_raw ne VALUE_NA and
      $pos_adj ne VALUE_NA) {
    $cd = abs($pos_raw - $pos_adj);
  }
  die unless $end eq "a" or $end eq "b";
  my $key = sprintf 'gene%s_coding_distance', $end;
  $row->{$key} = $cd;
  return $cd;
}

sub coding_distance_counts {
  my $f = $FLAGS{"coding-distance-counts"};

  my $df = new DelimitedFileHP(
			     "-file" => $f,
			    );
  $df->prepare_query("-fields" => [
				   qw(
				       input_row_md5
				       genea_coding_distance
				       geneb_coding_distance
				       extended_500
				    )
				  ]
		    );
  # SortedColumStreamer won't work because we may see the same
  # event later from a different source

  my %all;
  my $count = 0;
  while ($df->next_row()) {
    my ($md5, $cd_a, $cd_b, $ext) = @{$df->get_query};
    my $row = [ $cd_a, $cd_b, ($ext ? 1 : 0) ];
#    printf STDERR "save %s: %s\n", $md5, join " ", @{$row};
    push @{$all{$md5}}, $row;
    if (0 and ++$count >= 20000) {
      print STDERR "debug, stopping parse\n";
      last;
    }
  }

  my $count_usable = 0;
  my $count_unusable = 0;

  my @thresholds = (0, 10, 100, 1000, 10000, 100000);
  # how many source records are usable and have coding distances
  # <= each of these thresholds?

  my %tracker;
  foreach my $pid (keys %all) {
    my $set = $all{$pid};
    my $is_usable = 0;
#    printf STDERR "load %s: %d rows\n", $pid, scalar @{$set};

    my %hit;

    foreach my $row (@{$set}) {
      die unless @{$row} == 3;
#      printf STDERR "load row %s: %s\n", $pid, join " ", @{$row};
      my ($cd_a, $cd_b, $usable) = @{$row};
      if ($usable) {
	$is_usable = 1;
	foreach my $t (@thresholds) {
	  if ($cd_a <= $t and $cd_b <= $t) {
	    $hit{$t} = 1;
	  }
	}
      }
    }

    foreach my $t (keys %hit) {
      $tracker{$t}++;
    }

    if ($is_usable) {
      $count_usable++;
    } else {
      $count_unusable++;
    }
  }

  printf STDERR "total: %d\n", scalar keys %all;
  printf STDERR "usable: %d\n", $count_usable;
  printf STDERR "usable at threshold:\n";
  foreach my $t (@thresholds) {
    printf STDERR "%7d: %d\n", $t, $tracker{$t} || 0;
  }

}

sub find_missing_genes {
  my @files = glob($FLAGS{"find-missing-genes"} || die);
  die "no files" unless @files;

  my $pattern = sprintf '%s=([^\,]+)', TAG_ERROR_NO_SEQUENCES_FOR_GENE;
  study $pattern;

  my %missing;
  foreach my $file (@files) {
    open(IN, $file) || die;
    while (<IN>) {
      foreach my $gene (/$pattern/g) {
	$missing{$gene}++;
      }
    }
  }

  printf STDERR "missing gene references:\n";
  foreach my $g (sort keys %missing) {
    printf STDERR "  %s: %d\n", $g, $missing{$g};
  }

}

sub resolve_subsume {
  my ($pid, $subsumed) = @_;

  my $level = 1;

  while ($subsumed->{$pid}) {
    printf STDERR "level %d: %s => %s\n", $level, $pid, $subsumed->{$pid};
    $pid = $subsumed->{$pid};
    $level++;
  }

  return $pid;
}

sub full_length_filter {
  # using results from a (filtered) .pattern.tab file, filter
  # full-length fusion output from CICERO sv_inframe.pl
  my $f_patterns = $FLAGS{"full-length-filter"};
  my $df = new DelimitedFile("-file" => $f_patterns,
			     "-headers" => 1,
			     );

  #
  #  find accession pairings used in .pattern.tab file:
  #
  my %wanted_pairs;
  while (my $row = $df->get_hash()) {
    my $accA = $row->{genea_acc} || die;
    my $accB = $row->{geneb_acc} || die;

    foreach ($accA, $accB) {
      s/\.\d+$//;
      # sv_inframe uses refFlat which is unversioned
    }

    $wanted_pairs{$accA}{$accB} = 1;
  }

  #
  #  filter each full_length.tab file to those pairings only:
  #
  my @flf = glob("*.full_length.tab");
  die "no full_length files" unless @flf;

  foreach my $flf (@flf) {
    my $df = new DelimitedFile("-file" => $flf,
			       "-headers" => 1,
			      );
    my $outfile = basename($flf) . ".filtered.tab";
    my $rpt = $df->get_reporter(
				"-file" => $outfile,
				"-auto_qc" => 0,
			       );

    my $usable = 0;
    my $skipped = 0;

    while (my $row = $df->get_hash()) {
      my $accA = $row->{gene_a_nm} || die;
      my $accB = $row->{gene_b_nm} || die;
      if ($wanted_pairs{$accA}{$accB}) {
	$rpt->end_row($row);
	$usable++;
      } else {
	$skipped++;
      }
    }
    $rpt->finish();
    printf STDERR "%s: usable=%d skipped=%d\n", $flf, $usable, $skipped;
  }
}

sub write_fasta_seq {
  my ($fh, $id, $seq) = @_;
  $seq =~ tr/[]{}/NNNN/;
  printf $fh ">%s /length=%d\n%s\n", $id, length($seq), uc($seq);
}

sub get_sjpi {
  my ($sjpi, $acc) = @_;
  return $REPORT_SJPI ? ($SJPI->is_preferred_nm($acc) || 0) : "";
}

sub isofox_convert {
  my $f_in = $FLAGS{isofox} || die;
  my $df = new DelimitedFile("-file" => $f_in,
			     "-headers" => 1,
			     );
  my $f_main = basename($f_in) . ".pattern.tab";
  my $f_fz2 = basename($f_in) . ".fz2.tab";

  my $rpt_main = $df->get_reporter(
			      "-file" => $f_main,
			      "-extra" => [
					   $F_PATTERN_ID
					  ],
  			      "-auto_qc" => 1,
			     );
  # source file with pattern ID appended (needed for condense step)

  my $rpt_patterns = new Reporter(
			 "-file" => $f_fz2,
			 "-delimiter" => "\t",
			 "-labels" => [
				       $F_PATTERN_ID,
				       $F_PATTERN_SEQUENCE
				      ],
			 "-auto_qc" => 1,
			);
  # fuzzion2 patterns file

  my %pcount;

  while (my $row = $df->get_hash()) {
    my $gene = $row->{GeneName} || die;
    my $pattern_number = ++$pcount{$gene};
    my $kmer = $row->{Kmer} || die;
    my $sc = $row->{SoftClipBases} || die;

    my $idx = index($kmer, $sc) || die;
    # expected to be on right side of seq
    die unless index($kmer, $sc, $idx + 1) == -1;
    # unique

    my $up = substr($kmer, 0, $idx);
    my $down = substr($kmer, $idx);

    my $pattern = join '][', $up, $down;
    my $pid = sprintf '%s-softclip-%04d', $gene, $pattern_number;

    my %r;
    $r{$F_PATTERN_ID} = $pid;
    $r{$F_PATTERN_SEQUENCE} = $pattern;
    $rpt_patterns->end_row(\%r);

    $row->{$F_PATTERN_ID} = $pid;
    $rpt_main->end_row($row);
  }

  $rpt_main->finish();
  $rpt_patterns->finish();

}

sub readthrough_annotate {
  #
  # annotate possible readthrough fusions
  #
  my $f_fz2 = $FLAGS{"readthrough-annotate"} || die;
  # fz2 patterns file

  my @rfa_options;
  if (my $genome = $FLAGS{genome}) {
    push @rfa_options, "-genome" => $genome;
  } else {
    my $host = $FLAGS{"blat-host"} || die "-genome | -blat-host HOST -blat-port PORT";
    my $port = $FLAGS{"blat-port"} || die "-genome | -blat-host HOST -blat-port PORT";
    push @rfa_options, (
			"-blat-host" => $host,
			"-blat-port" => $port,
		       );
  }
  # for blat client, replace this w/host/port

  my $rf = get_refflat();

  my $rfa = new ReadthroughFusionAnnotator(
					   "-rf" => $rf
					  );

  my $df = new DelimitedFile("-file" => $f_fz2,
			     "-headers" => 1,
			    );

  my $outfile = basename($f_fz2) . ".readthrough.tab";


  my @h_new = qw(
		  readthrough_fusion
		  readthrough_fusion_info
		  readthrough_interval
	       );
  my $h_needed = detect_needed_headers($df, \@h_new);

  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => $h_needed,
  			      "-auto_qc" => 1,
			     );

  while (my $r = $df->get_hash()) {
    my $pattern = $r->{pattern} || die;
    my $geneA_raw = pattern_gene($pattern, $r->{"genea_symbol"}) || die;
    my $geneB_raw = pattern_gene($pattern, $r->{"geneb_symbol"}) || die;
    # ideally would extract from pattern name itself, however
    # "-" characters may themselves be used in pattern gene names

    my $sequence = $r->{sequence} || die;
    $sequence =~ tr/[]{}//d;
    # fields are always present in fz2 pattern files

    my $possible_readthrough = $rfa->is_readthrough_fusion(
							   "-geneA" => $geneA_raw,
							   "-geneB" => $geneB_raw,
							   "-contig" => $sequence,
							   @rfa_options
							  );
    $r->{readthrough_fusion} = $possible_readthrough;
    $r->{readthrough_interval} = $possible_readthrough ? $rfa->rt_interval() : "";
    $r->{readthrough_fusion_info} = join ",", @{$rfa->rt_info()};

    $rpt->end_row($r);
  }

  $rpt->finish();

}

sub readthrough_source_summary {
  # summarize readthrough calls by data source
  my %source2rows;

  my $infile = $FLAGS{"readthrough-source-summary"} || die;

  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );



  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
    my @sources = split /,/, ($row->{source} || die);
    foreach my $src (@sources) {
      push @{$source2rows{$src}}, $row;
    }
  }

  my $rpt = new Reporter(
			 "-file" => basename($infile) . ".source.tab",
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   source
					   gene_pairs
					   patterns
					   rt_simple
					   rt_2l2
					   rt_total
					   rt_gene_pairs
					   rt_gene_pairs_fraction
					   rt_gene_pairs_detail
					)
				      ],
			 "-auto_qc" => 1,
			);


  foreach my $src (sort keys %source2rows) {
    my $rows = $source2rows{$src} || die;

    my @pids = map {$_->{pattern}} @{$rows};
    my %gene_pairs;
    foreach my $pid (@pids) {
      my $pair = pattern2pair($pid);
      $gene_pairs{$pair} = 1;
    }

    my @rt_rows;
    my $rt_simple = 0;
    my $rt_2l2 = 0;
    my %rt_gene_pairs;
    foreach my $r (@{$rows}) {
      if (my $rt = $r->{readthrough_fusion}) {
	push @rt_rows, $r;
	my $pair = pattern2pair($r->{pattern}) || die;
	$rt_gene_pairs{$pair} = 1;
	if ($rt eq "1") {
	  $rt_simple++;
	} elsif ($rt eq "second_to_last_exon_to_second_exon") {
	  $rt_2l2++;
	} else {
	  die "unhandled rt code $rt";
	}
      }
    }

    my %r;
    $r{source} = $src;
    $r{patterns} = scalar @{$rows};
    $r{gene_pairs} = scalar keys %gene_pairs;
    $r{rt_simple} = $rt_simple;
    $r{rt_2l2} = $rt_2l2;
    $r{rt_total} = $rt_simple + $rt_2l2;

    my $rt_pair_count = scalar keys %rt_gene_pairs;
    $r{rt_gene_pairs} = $rt_pair_count;
    $r{rt_gene_pairs_fraction} = sprintf "%.03f", $rt_pair_count / scalar keys %gene_pairs;
    $r{rt_gene_pairs_detail} = join ", ", sort keys %rt_gene_pairs;
    $rpt->end_row(\%r);
  }

  $rpt->finish();
}

sub pattern2pair {
  my ($pid) = @_;
  my $pair = $pid;
  $pair =~ s/\-\d+$// || die;
  return $pair;
}

sub insertion2fuzzion {
  # generate fz2 patterns for contigs containing insertion
  #
  # inputs:
  #  - contig sequence containing insertion or ITD vs. reference
  #  - refSeq accession to use for comparison (STRICT matching required)
  # output: fz2 pattern file
  die "TO DO: update code to work for insertions only, NOT ITD!";
  my $f_in = $FLAGS{insertion2fuzzion} || die;
  find_binary("blastn", "-die" => 1);

  my $max_subject_overlap = 35;
  # TBD

  if (not($BLAST_GAP_OPEN and $BLAST_GAP_EXTEND)) {
    # disabled for now: a little messiness is OK, e.g. UBTF:
    # if we suppress gaps this may set ITD boundaries in normal
    # part of sequence!
    # - effect on left side of contig vs gene model:
    #   - defaults: query matches gene to 1342 w/short insertion in query [OK?]
    #   - with penalty 1000: matches to only 1308 (BAD)
    # - HOWEVER: without penalty, SUBJECT REGION HIT OVERLAPS (BAD)
    $BLAST_GAP_OPEN = $BLAST_GAP_EXTEND = 1000;
    printf STDERR "setting BLAST gap open/extend penalty to %d\n", $BLAST_GAP_OPEN;
  }

  my $f_gene = "gene";
  my $f_refseq = "refseq";
  my $f_contig = "contig";

  my $force_gene = $FLAGS{gene};
  my $force_acc = $FLAGS{accession};
  # manual values for entire file

  my @f_check = $f_contig;
  push @f_check, $f_gene unless $force_gene;
  push @f_check, $f_refseq unless $force_acc;

  #
  #  pass 1: identify referenced refSeqs
  #
  my $df = new DelimitedFile("-file" => $f_in,
			     "-headers" => 1,
			     );
  my %needed_rf;
  while (my $row = $df->get_hash()) {
    foreach my $f (@f_check) {
      dump_die($row, "no data for $f") unless $row->{$f};
    }
    my $acc = $force_acc || $row->{$f_refseq} || die;
    $needed_rf{$acc} = 1;
    # ugh: complications if given accession is versioned!!
  }

  #
  #  fetch GenBank records for needed refseqs:
  #
  my $gbc = new GenBankCache();
  my $result_files = $gbc->get([sort keys %needed_rf]);

  my %acc2cds;

  foreach my $file (@{$result_files}) {
    my $gbc = new GenBankCDS("-file" => $file);
    my $info = $gbc->get_info();
    my $acc = $info->{accession} || die;
    my $cds = $info->{mrna_cds} || die;
    $acc2cds{$acc} = $cds;
  }

  my $f_out = $FLAGS{out} || basename($f_in) . ".fz2.tab";
  my $f_failed = basename($f_in) . ".fz2_failed.tab";

  my $rpt = new Reporter(
			 "-file" => $f_out,
			 "-delimiter" => "\t",
			 "-labels" => [
					$F_PATTERN_ID,
					$F_PATTERN_SEQUENCE,
				        $F_SAMPLE,
					qw(
					    genea_symbol
					    genea_acc
					    genea_feature
					    geneb_symbol
					    geneb_acc
					    geneb_feature
					 ),
				       "detected_repeat",
				       "exception",
				       $f_contig,
				       "gene_sequence"
				       ]
			);
  # fz2 output
  # TO DO: standard variable(s) for pattern file header fields?

  my $f_sample = guess_sample_column($df);


  #
  #  second pass: process
  #
  $df = new DelimitedFile(
			  "-file" => $f_in,
			  "-headers" => 1,
			 );
  my $rpt_failed = $df->get_reporter("-file" => $f_failed,
				     "-extra" => [ "exception" ]
				    );

  my $min_query_edge_anchor = 0.05;
  # expect start and end of contig mapping to be anchored with this
  # percentage tolerance.  e.g. 0.05 means the contig mapping must
  # begin within 5% of the start of the contig and end within 5% of
  # the end of the ocntig.  This ensures a good mapping between the
  # contig and the reference sequence before and after the
  # ITD/insertion.

  my %pattern_counter;
  my %saw_contig;

  while (my $row = $df->get_hash()) {
    my $contig = $row->{$f_contig} || die;
    printf STDERR "new row: %s\n", $contig;
    my $acc = $force_acc || $row->{$f_refseq} || die;
    my $gene = $force_gene || $row->{$f_gene} || die;
    my $rs_seq = $acc2cds{$acc} || die "can't get cds for $acc";

    my $blast = get_blast();
    my $parser = $blast->blast(
			       "-query" => {
					    "contig" => $contig,
					   },
			       "-database" => {
					       $acc => $rs_seq
					      }
			      );
    my $result = $parser->next_result();
    # one object per query sequence (only one query seq)

    my @hits;
    while (my $hit = $result->next_hit()) {
      push @hits, $hit;
    }
    dump_die($row, "no blast hits") unless @hits;
    die "not exactly one hit" unless @hits == 1;

    my ($hsp_left, $hsp_right);
    my $contig_len = length $contig;
    my $hsp_count = 0;
    my ($hsp_left_pos, $hsp_right_pos);
    while (my $hsp = $hits[0]->next_hsp) {
      $hsp_count++;
      printf STDERR "    score:%s strand:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s query_span:%d ref_span:%d total_span=%d query_string=%s hit_string=%s\n",
	  $hsp->score,
	    $hsp->strand,
	      $hsp->range("query"),
		$hsp->range("hit"),
		  $hsp->num_identical(),
		    $hsp->frac_identical("query"),
		      $hsp->length("query"),
			$hsp->length("hit"),
			  $hsp->length("total"),
			    $hsp->query_string(),
			      $hsp->hit_string() if $VERBOSE;

      my $qs = $hsp->start("query");
      my $qe = $hsp->end("query");

      if (not(defined $hsp_left_pos) or $qs < $hsp_left_pos) {
	$hsp_left = $hsp;
	$hsp_left_pos = $qs;
#	print STDERR "  => left\n";
      }

      if (not(defined $hsp_right_pos) or $qe > $hsp_right_pos) {
	$hsp_right = $hsp;
	$hsp_right_pos = $qe;
#	print STDERR "  => right\n";
      }
    }
    die "no left HSP" unless $hsp_left;
    die "no right HSP" unless $hsp_right;

    my $exception = "";

    #
    #  QC checks:
    #
    if ($saw_contig{$contig}) {
      $exception = "duplicate_contig";
    } elsif ($hsp_count == 1) {
      # method requires 2 sub-alignments
      $exception = "only_one_alignment";
    } else {
      #
      # ensure anchoring at contig start and end:
      #
      my $frac_start_distance = $hsp_left->start("query") / $contig_len;

      if ($frac_start_distance > $min_query_edge_anchor) {
	printf STDERR "ERROR: contig start not mapped within $min_query_edge_anchor: $frac_start_distance\n";
	$exception = "anchor_failed_left";
      }
      my $frac_end_distance = 1 - $hsp_right->end("query") / $contig_len;
      if ($frac_end_distance > $min_query_edge_anchor) {
	printf STDERR "ERROR: contig end not mapped within $min_query_edge_anchor: $frac_end_distance\n";
	$exception = "anchor_failed_right";
      }
    }

    if ($exception) {
      $row->{exception} = $exception;
      $rpt_failed->end_row($row);
      next;
    }
    $saw_contig{$contig} = 1;

    hsp_overlap_check($hsp_left, $hsp_right, "query", 1);
    # consider overlap a fatal error
    my $subject_overlap = hsp_overlap_check($hsp_left, $hsp_right, "subject", 0);
    # a small overlap due to microhomology is OK, just needs to be offset.
    # however a large overlap e.g. UBTF w/gap penalties disabled is bad!
#    die "large subject overlap of $subject_overlap nt" if $subject_overlap > $max_subject_overlap;

    my $left_contig_end = $hsp_left->end("query");
    my $left_refseq_end = $hsp_left->end("subject");

    my $right_contig_start = $hsp_right->start("query");
    my $right_refseq_start = $hsp_right->start("subject");

    if ($subject_overlap) {
      # if a very short sequence (likely microhomology) it might be
      # desirable to move breakpoint away from BOTH sites so both ends
      # of the ambiguous region will remain within ITD brackets.
      #
      # HOWEVER, in some cases the overlap is much larger and related
      # to the ITD, in which case there may not be enough alignment
      # available to do that.  Compromise: subtract the overlap from
      # the longer alignment.

#      printf STDERR "subject overlap of %d, adjusting both sides\n", $subject_overlap;

      # subtract the overlap from the longer alignment:
      if ($hsp_left->length("query") > $hsp_right->length("query")) {
	$left_contig_end -= $subject_overlap;
	$left_refseq_end -= $subject_overlap;
      } else {
	$right_contig_start += $subject_overlap;
	$right_refseq_start += $subject_overlap;
      }
    }

    printf STDERR "final: contig breaks:%d,%d refseq breaks:%d,%d\n",
      $left_contig_end, $right_contig_start,
	$left_refseq_end, $right_refseq_start;

    my $left_up_full = substr($rs_seq, 0, $left_refseq_end);
    my $right_down_full = substr($rs_seq, $right_refseq_start - 1);

    printf STDERR "up: %s\n", $left_up_full;
    printf STDERR "down: %s\n", $right_down_full;

    my $interstitial = substr($contig,
			      $left_contig_end,
			      # 1st base after left portion ends
			      (($right_contig_start - $left_contig_end) - 1)
			     );
    printf STDERR "interstitial: %s\n", $interstitial;

    my $left_up = substr($left_up_full, - $EXTENDED_CHUNK_LENGTH);
    my $right_down = substr($right_down_full, 0, $EXTENDED_CHUNK_LENGTH);

    my $fz2_pattern = sprintf '%s%s%s%s%s',
      $left_up, BREAKPOINT_ITD_LEFT,
	$interstitial,
	  BREAKPOINT_ITD_RIGHT, $right_down;

    #
    #  TO DO:
    #  attempt to detect REPEATED sequence in the ITD:
    #  - blast interstitial/ITD portion of pattern vs. the source contig
    #    - using source contig will detect repeats in the contig,
    #      however might be better to use gene model to detect
    #      REFERENCE repeats
    #  - are there multiple HSPs?
    #  - does the sequence of the shortest HSP appear in a longer HSP?
    #    If so, that's the repeat
    #

    # NOTE that interstitial sequence may not quite contain the entire repeat,
    # due to blast alignment
    my $extra = 10;
    my $check = $interstitial;
#    my $check = join "",
#      substr($left_up, - $extra),
#	$interstitial,
#	  substr($right_down, 0, $extra);

    $parser = $blast->blast(
			       "-query" => {
					    "interstitial" => $check,
					   },
			       "-database" => {
					       "contig" => $contig
#					       $acc => $rs_seq
					      }
			      );

    $result = $parser->next_result() || die;
    my $hit = $result->next_hit() || die "no hit";
    # just 1

    my @hsps;
    while (my $hsp = $hit->next_hsp()) {
      push @hsps, $hsp;
    }
    @hsps = sort {$a->length("query") <=> $b->length("query")} @hsps;
    # sort by query length chunk size, shortest first

    my $hsp2chunk = sub {
      my ($hsp) = @_;
      my $qs = $hsp->start("query");
      my $qe = $hsp->end("query");
#      my $chunk = substr($interstitial,
      return substr($check,
		    $qs - 1,
		    (($qe - $qs) + 1));
    };

    foreach my $hsp (@hsps) {
      printf STDERR "    score:%s strand:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s query_span:%d ref_span:%d total_span=%d query_string=%s hit_string=%s\n",
	  $hsp->score,
	    $hsp->strand,
	      $hsp->range("query"),
		$hsp->range("hit"),
		  $hsp->num_identical(),
		    $hsp->frac_identical("query"),
		      $hsp->length("query"),
			$hsp->length("hit"),
			  $hsp->length("total"),
			    $hsp->query_string(),
			      $hsp->hit_string() if $VERBOSE;
    }

    #
    #  detect longest repeated subsequence:
    #
    printf STDERR "chunk count: %d\n", scalar @hsps;
    my $repeat_seq = "";
    if (@hsps > 1) {
      for (my $i = 0; $i < @hsps; $i++) {
	# don't just use the first/shortest sequence because longer
	# entries may also repeat in longer-still ones
	my $shortest = &$hsp2chunk($hsps[$i]);
	print STDERR "check idx $i $shortest\n";
	for (my $j = $i + 1; $j < @hsps; $j++) {
	  if (index(&$hsp2chunk($hsps[$j]), $shortest) != -1) {
	    # found shorter sequence within longer one, promising

	    my $overlap = hsp_overlap_check($hsps[$i], $hsps[$j], "subject");
	    printf STDERR "overlap: %d\n", $overlap;
	    unless ($overlap) {
	      # don't count matches which overlap the query sequence
	      $repeat_seq = $shortest;
	      printf STDERR "idx $i, assign repeat sequence: %s\n", $repeat_seq;
	      # yay
	    }
	  }
	}
      }
    }

    my %r;
    my $pid = sprintf '%s-%s-%02d', $gene, $gene, ++$pattern_counter{$gene};

    $r{$F_PATTERN_ID} = $pid;
    $r{$F_PATTERN_SEQUENCE} = $fz2_pattern;
    $r{genea_symbol} = $r{geneb_symbol} = $gene;
    $r{genea_feature} = $r{geneb_feature} = "";
    $r{genea_acc} = $r{geneb_acc} = $acc;
    $r{detected_repeat} = $repeat_seq;
    $r{$f_contig} = $contig;
    $r{gene_sequence} = $rs_seq;
    $r{exception} = "";
    $r{$F_SAMPLE} = $f_sample ? $row->{$f_sample} : "";

    $rpt->end_row(\%r);
  }

  $rpt->finish();
  $rpt_failed->finish();
}

sub pattern_gene {
  # hack to filter a gene annotation list to names associated with
  # pattern.  example: for pattern ANKHD1-ARHGAP26-01 and gene
  # annotation ANKHD1,ANKHD1-EIF4EBP3, ANKHD1 should be used because
  # this appears in the pattern name but ANKHD1-EIF4EBP3 doesn't.
  #
  # HOWEVER this can be stymied by disambiguity updates to gene
  # symbols in pattern names  :/
  my ($pattern, $gene_raw) = @_;
  die unless $gene_raw;

  my $result = $gene_raw;
  if ($gene_raw =~ /,/) {
    my @genes = split /,/, $gene_raw;
    my @usable = grep {index($pattern, $_) != -1} @genes;
    if (@usable == 1) {
      # just a single result found
      $result = $usable[0];
    } elsif (@usable > 1) {
      # multiple possibilities, e.g.:
      #
      #  - gene list = C7orf55,C7orf55-LUC7L2
      #  - pattern = C7orf55-LUC7L2-IQSEC1-01
      #
      # in this both genes are found in pattern name, however
      # C7orf55-LUC7L2 is the correct one.
      my @sorted = sort {length($b) <=> length($a)} @usable;
      $result = $sorted[0];
      # so, consider the longest match the correct one
    } elsif (@usable == 0) {
      # check to see if gene name in pattern has been updated
      # from source gene name (i.e. for harmonization).
      # example: SPMIP7-IKZF1-01, gene annotation is for C7orf72,
      # an earlier symbol for SPMIP7
      my $p = $pattern;
      $p =~ s/\-\d+$//;
      my @pg = split /\-/, $p;
      # genes in pattern ID

      my $gsm = new_gsm_lite();
      foreach my $g (@pg) {
	$gsm->add_gene("-gene" => $g);
      }

      my %usable;
      foreach my $g (@genes) {
	if (my $hit = $gsm->find($g)) {
	  $usable{$hit} = 1;
	  # use the gene symbol from the pattern name rather than
	  # the gene annotation, as it may be more up to date
	}
      }
      @usable = sort keys %usable;
      confess "$gene_raw disambiguation failure" unless @usable == 1;
    } else {
      die "can't disambiguate $gene_raw in $pattern, usable=" . scalar @usable;
    }
    $result = $usable[0];
  }
  return $result;
}

sub field_preset_setup {
  if ($FLAGS{crest}) {
    $F_CHR_A = 'chrA';
    # in raw files, header line is commented, so may be "#chrA"
    $F_POS_A = 'posA';
    $F_ORT_A = 'ortA';
    $F_CHR_B = 'chrB';
    $F_POS_B = 'posB';
    $F_ORT_B = 'ortB';
  }
}


sub init_sjpi {
  unless ($SJPI) {
    my $f_sjpi = $FLAGS{sjpi};
    unless ($f_sjpi) {
      my $genome = $FLAGS{genome} || die "-genome";
      my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
      $f_sjpi = $config_genome->{GENE_TRANSCRIPT_MATRIX} || die;
    }
    $SJPI = new SJPreferredIsoform(
				   "-file" => $f_sjpi,
				   "-auto_gsm" => 1
				   );
  }
}

sub hsp_overlap_check {
  my ($hsp_left, $hsp_right, $type, $is_fatal) = @_;
  my $si_l = hsp2si($hsp_left, $type);
  my $si_r = hsp2si($hsp_right, $type);
  my $intersect = intersect $si_l $si_r;
  my $intersect_size = $intersect->size();
  if ($intersect_size and $is_fatal) {
    printf STDERR "overlap size=%d\n", $intersect->size();
    confess sprintf "ERROR: BLAST HSPs overlap in %s: %s\n", $type, run_list $intersect;
    # a good reason to keep gap penalty in effect
  }
  return $intersect_size;
}

sub hsp2si {
  my ($hsp, $type) = @_;
  return new Set::IntSpan(sprintf '%d-%d',
			  $hsp->start($type),
			  $hsp->end($type));
}

sub pattern_summary {
  my $f_patterns = $FLAGS{"pattern-summary"} || die;

  my $df = new DelimitedFile("-file" => $f_patterns,
			     "-headers" => 1,
			     );
  my $outfile = basename($f_patterns) . ".summary.tab";

  my %info;
  while (my $row = $df->get_hash()) {
    # do NOT use gene, feature, etc. fields as these may be condensed
    # from multiple pattern sources during duplicate removal.  Instead use
    # gene_pair_summary which preserves all the pairing info from the
    # original patterns.
    foreach my $summary (split /,/, $row->{gene_pair_summary} || die) {
      my ($gene_a, $gene_b, $exon_pair) = summary_to_pair_info($summary);

      if ($exon_pair) {
	dump_die($row, "commas in genes") if $gene_a =~ /,/ or $gene_b =~ /,/;
	my $key = join "|", $gene_a, $gene_b;
	$info{$key}{exon_pairs}{$exon_pair}++;
	$info{$key}{pattern_count}++;
	$info{$key}{geneA} = $gene_a;
	$info{$key}{geneB} = $gene_b;
      }
    }
  }

  my $rpt = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   geneA
					   geneB
					   pattern_count
					   exon_pairs
					)
				      ],
			 "-auto_qc" => 1,
			);

  foreach my $key (sort keys %info) {
    my %r = %{$info{$key}};
    my $exon_pairs = $r{exon_pairs};
    my $pair_summary = "";
    if ($exon_pairs) {
      $pair_summary = join ", ", map {$_ . "=" . $exon_pairs->{$_}} sort keys %{$exon_pairs};
    }
    $r{exon_pairs} = $pair_summary;
    $rpt->end_row(\%r);
  }
  $rpt->finish();
}

sub summary_to_pair_info {
  my ($summary) = @_;
  my @pair = split /\-/, $summary;
  if (@pair == 2) {
    # OK
  } elsif ($summary =~ /^(.*(exon|utr)_\d+)\-(.*)$/) {
    # unfortunately the gene names sometimes also contain dashes  :/
    # ABTB2/NM_145804.3/exon_1-GRM7-AS3/NR_110123.1/3utr_3
    @pair = ($1, $3);
#    die join "\n", $summary, @pair;
  } else {
    die "can't parse $summary";
  }
  my @genes;
  my @exons;
  foreach my $side (@pair) {
    my @f = split /\//, $side;
    die unless @f == 3;
    push @genes, $f[0];

    my $exno = $side =~ /(exon|utr)_(\d+)$/ ? $2 : "";
#    printf STDERR "can't get exno for $side\n" unless $exno;
    push @exons, $exno;
  }
  die unless @genes == 2;
  die "can't find pair in $summary" unless @exons == 2;
  my $ex_pair = ($exons[0] and $exons[1]) ? join("-", @exons) : "";
  return (@genes, $ex_pair);
}

sub append_patterns {
  # crude hack to append new/unique patterns to an existing set
  my $f_orig = $FLAGS{"append-orig"} || die;
  my $f_new = $FLAGS{"append-new"} || die;

  my $f_out = basename($f_orig) . ".append.tab";

  my $df = new DelimitedFile("-file" => $f_orig,
			     "-headers" => 1,
			     );
  my $rpt = $df->get_reporter(
			      "-file" => $f_out,
			     );

  my %saw;
  my @h_orig = @{$df->headers_raw};
  while (my $row = $df->get_hash()) {
    # copy old patterns first and track IDs
    my $pid = $row->{pattern};
    $saw{$pid} = 1;
    $rpt->end_row($row);
  }

  $df = new DelimitedFile("-file" => $f_new,
			  "-headers" => 1,
			 );
  my $source = file_to_source($f_new);
  while (my $row = $df->get_hash()) {
    # copy old patterns first and track IDs
    my $pid = $row->{pattern};
    die "duplicate pattern $pid" if $saw{$pid};
    $saw{$pid} = 1;
    $row->{$F_SOURCE} = $source;

    foreach my $h (@h_orig) {
      # add blanks to any additional columns.
      # TO DO: maybe use "n/a"?
      $row->{$h} = "" unless exists $row->{$h};
    }

    $rpt->end_row($row);
  }

  $rpt->finish();



}

sub guess_sample_column {
  my ($df) = @_;
  my @guess = grep {/sample/i} @{$df->headers_raw};
  my $result;
  $result = $guess[0] if @guess;
  return $result;
}

sub add_preferred_key {
  # pattern file postprocessing: create a grouping/sorting code
  # consisting of the most-preferred isoform pair contained in the pattern.
  # report for both A and B: gene (already present), transcript, feature.
  # Also include a single-column version for simple grouping.
  my $f_in = $FLAGS{"add-preferred-key"} || die;
  my $generate_transcript_info = $FLAGS{"generate-transcript-info"};
  my $dwim_checked;

  if ($generate_transcript_info) {
    # assumes CICERO format
    $F_GENE_A = "geneA";
    $F_GENE_B = "geneB";
    # these additional fields only required for generate_transcript_info()
    # call below:
    $F_CHR_A = "chrA";
    $F_CHR_B = "chrB";
    $F_POS_A = "posA";
    $F_POS_B = "posB";
  } else {
    # fz2 pattern annotations
    $F_GENE_A = "genea_symbol";
    $F_GENE_B = "geneb_symbol";
  }

  init_sjpi();

  my $df = new DelimitedFile("-file" => $f_in,
			     "-headers" => 1,
			     );
  my $outfile = basename($f_in) . ".sjpi.tab";

  my @f_new;
  push @f_new, "gene_pair_summary" if $generate_transcript_info;
  push @f_new, qw(
		  gte_A_gene
		  gte_A_transcript
		  gte_A_exon
		  gte_B_gene
		  gte_B_transcript
		  gte_B_exon
		  gte_rank
		  gte_key
	       );

  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => \@f_new,
  			      "-auto_qc" => 1,
			     );
  # can't necessarily rely on existing gene A/B symbols, as these
  # may contain a list, e.g. GCOM1,MYZAP (coding and non-coding).
  # Instead, select the most-preferred entry from "gene_pair_summary".

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {

    if ($generate_transcript_info) {
      unless ($dwim_checked) {
	# munge header name capitalization, e.g. for formatting changes
	# between cicero and cicero-post output files, which I'm sure
	# were made for a very important and well-thought-out reason
	field_name_dwim(
			"-query" => [
				     \$F_CHR_A,
				     \$F_POS_A,
				     \$F_GENE_A,
				     \$F_CHR_B,
				     \$F_POS_B,
				     \$F_GENE_B,
				    ],
			"-actual" => $df->headers_raw(),
		       );
	$dwim_checked = 1;
      }
      generate_transcript_info($row);
    }
    dump_die($row, "no gene_pair_summary field") unless exists $row->{gene_pair_summary};
    my @set = split /,/, $row->{gene_pair_summary};

    my $a_gene = "";
    my $a_transcript = "";
    my $a_feature = "";
    my $b_gene = "";
    my $b_transcript = "";
    my $b_feature = "";
    my $gte_rank = "";

    if (@set) {
      # - if set has both NM_ and NR_, choose NM_
      # - bucket each entry by rank
      # - get best rank
      # - if multiple records, sort alphanumerically
      # - QUESTION: prefer by feature annotation??
      if (grep {/NM_/} @set and grep {/NR_/} @set) {
	# exclude non-coding-only records if we have both NM and NR
	@set = grep {/NM_/} @set;
      }

      # need to consider rank for BOTH sides simultaneously.  consider:
      # A rank     B rank
      # ------     ------
      # 2          4
      # 3          1
      #
      # ranking A first ends up with an inferior pair (2+4) vs.
      # evaluating together (3+1)
      #
      # - always prefer rank 1 if available in a pair, so that
      #   at least one side gets the known-best isoform?
      # - different methods depending on ranks available?:
      #   - ranks for both isoforms known
      #   - ranks for just one isoform known
      #   - no ranks known

      my @have_pair;
      my @have_single;
      my @have_none;

      foreach my $pair (@set) {
#	printf STDERR "pair: %s\n", $pair;
	my $sides = split_pair($pair);
	my $rank_a = get_summary_rank($sides->[0]);
	my $rank_b = get_summary_rank($sides->[1]);

#	printf STDERR "%s: %s %s\n", $pair, map {$_ || "NA"} $rank_a, $rank_b;
	if ($rank_a and $rank_b) {
	  push @have_pair, [ $rank_a + $rank_b, $pair ];
	} elsif ($rank_a or $rank_b) {
	  push @have_single, [ ($rank_a || 0) + ($rank_b || 0), $pair ];
	} else {
	  push @have_none, [ 0, $pair ];
	}
      }

      my $set;
      my $score_base;
      if (@have_pair) {
	# best
	$set = \@have_pair;
	$score_base = 0;
      } elsif (@have_single) {
	# 2nd best
	$set = \@have_single;
	$score_base = 100;
      } elsif (@have_none) {
	# worst
	$set = \@have_none;
	$score_base = 200;
      } else {
	die;
      }
      die unless defined $score_base;
      my $best;
      ($best, $gte_rank) = pick_ranked($set, $score_base);
      my $sides = split_pair($best);
      ($a_gene, $a_transcript, $a_feature) = split /\//, $sides->[0];
      ($b_gene, $b_transcript, $b_feature) = split /\//, $sides->[1];
    } else {
      # can happen: use gene symbol from main annotation + "unknown"?
      $a_gene = $row->{$F_GENE_A} || die;
      $b_gene = $row->{$F_GENE_B} || die;
      $gte_rank = 300;
      # no isoform data
    }

    $row->{gte_A_gene} = $a_gene;
    $row->{gte_A_transcript} = $a_transcript;
    $row->{gte_A_exon} = $a_feature;
    $row->{gte_B_gene} = $b_gene;
    $row->{gte_B_transcript} = $b_transcript;
    $row->{gte_B_exon} = $b_feature;
    $row->{gte_key} = join ",", $a_gene, $a_transcript, $a_feature,
      $b_gene, $b_transcript, $b_feature;
    $row->{gte_rank} = $gte_rank;

    $rpt->end_row($row);
  }

  $rpt->finish();
}

sub get_summary_rank {
  my ($side) = @_;
  my @f = split(/\//, $side, 3);
  # specify 3rd parameter to split() for partially-empty entries,
  # e.g. arriba proximal events:
  # AC115283.1(19261)|AC099540.1(107270)//
  confess "error parsing sides in $side" unless @f == 3;
  return $SJPI->get_nm_rank($f[1]);
}

sub pick_ranked {
  my ($set, $score_base) = @_;
  my ($best_score) = sort {$a <=> $b} map {$_->[0]} @{$set};

  my @set_best = grep {$_->[0] == $best_score} @{$set};
  die unless @set_best;
  my $result;
  if (@set_best == 1) {
    # one winner
    $result = $set_best[0]->[1];
  } elsif (@set_best > 1) {
    # tie
    my @pairs = map {$_->[1]} @set_best;
    my $set_collapsed = collapse_features(\@pairs);

    if (@{$set_collapsed} == 1) {
      # multiple records, but the only difference was in feature annotations
      $result = $set_collapsed->[0];
    } else {
      # stochastic, just return first sorted entry
      my @sorted = sort {$a cmp $b} @{$set_collapsed};
      $result = $sorted[0];
    }
  } else {
    die;
  }
  die unless $result;
  die unless wantarray;
  return ($result, $best_score + $score_base);
}

sub collapse_features {
  my ($set_in) = @_;
  # collapse records by gene/transcript pairs only, e.g.
  #   ABL1/NM_007313.2/exon_1-BCR/NM_004327.4/exon_14
  #   ABL1/NM_007313.2/exon_1-BCR/NM_004327.4/exon_15
  # becomes
  #   ABL1/NM_007313.2/exon_1-BCR/NM_004327.4/exon_14|exon15
  my %features;
  foreach my $pair (@{$set_in}) {
    my $sides = split_pair($pair);
    my @side_a = split /\//, $sides->[0], 3;
    my @side_b = split /\//, $sides->[1], 3;
    # explicit 3: sometimes only the gene symbol is populated

    my $key_pair = join "/", @side_a[0,1], @side_b[0,1];
    $features{$key_pair}{A}{$side_a[2]} = 1;
    $features{$key_pair}{B}{$side_b[2]} = 1;
  }

  my @out;
  foreach my $pk (sort keys %features) {
    my ($a_gene, $a_nm, $b_gene, $b_nm) = split /\//, $pk;
    my @a_f = sort keys %{$features{$pk}{A}};
    my @b_f = sort keys %{$features{$pk}{B}};

    my $pair = sprintf '%s-%s',
      join("/", $a_gene, $a_nm, join "|", @a_f ),
	join("/", $b_gene, $b_nm, join "|", @b_f );
    push @out, $pair;
  }

  return \@out;
}

sub split_pair {
  my ($pair) = @_;
  my @sides = split /\-/, $pair;
  unless (@sides == 2) {
    # ANKHD1-EIF4EBP3/NM_020690.6/exon_1-ARHGAP26/NM_001135608.3/exon_5
#    $pair =~ /^(\S+\/\S+\/\S+)\-(\S+\/\S+\/\S+)$/ || confess "failed match $pair";
    $pair =~ /^(\S+\/\S*\/\S*)\-(\S+\/\S*\/\S*)$/ || confess "failed match $pair";
    # not guaranteed to be > 0 length:
    # RNU6-1_locC|RNU6-7_locC|RNU6-8_locC|RNU6-2_locC|RNU6-9_locC//-SNORA53/NR_003015/5utr_1
    @sides = ($1, $2);
  }
  return \@sides;

}

sub generate_transcript_info {
  # generate "gene_pair_summary" field for CICERO input
  my ($row) = @_;

  my $gene_a = $row->{$F_GENE_A} || die;
  my $gene_b = $row->{$F_GENE_B} || die;
  my $chr_a = $row->{$F_CHR_A} || die;
  my $chr_b = $row->{$F_CHR_B} || die;
  my $pos_a = $row->{$F_POS_A} || die;
  my $pos_b = $row->{$F_POS_B} || die;
  # TO DO: ADJUST positions for feature-tracking purposes?

  $rf = get_refflat() unless $rf;
  # global variable, hack

  foreach ($gene_a, $gene_b) {
    s/,/\|/g;
    # commas are a reserved character used by gene_pair_summary field.
    # reformat e.g. Arriba proximal entries like this:
    #   AC115283.1(19261),AC099540.1(107270)
  }

  my $rfs_a = $rf->find_by_gene($gene_a);
  my $rfs_b = $rf->find_by_gene($gene_b);
  # gene symbol ambiguity is automatically handled by default.
  # for Arriba-sourced events, these may sometimes be a list of
  # nearby/proximal genes which are outside of the gene model, e.g.:
  #
  # PNPT1(17964),EFEMP1(154093)
  #
  # don't attempt to salvage since they (a) reference more than one gene
  # and (b) are always outside the gene model.

#  printf STDERR "processing %s/%s\n", $gene_a, $gene_b;

  conditional_noncoding_filter($rfs_a);
  conditional_noncoding_filter($rfs_b);
  # to help limit enormous lists of combinations, if both NM and NR
  # records are present, use only NM

  # TO DO: if preferred isoforms are present in the list, use only those?
  # or keep everything in case we change our mind later?

  my %pairs;
  my %cache;

  if ($rfs_a and $rfs_b) {
    #
    # records for both genes in pair available
    #
    foreach my $rf_a (@{$rfs_a}) {
      next unless rf_overlap($rf_a, $chr_a, $pos_a);

      my $transcript_a = $rf_a->{name} || die;
      my $feature_a = get_feature_tag($rf, $rf_a, $pos_a, \%cache);

      foreach my $rf_b (@{$rfs_b}) {
	# build pair info
	next unless rf_overlap($rf_b, $chr_b, $pos_b);
	my $feature_b = get_feature_tag($rf, $rf_b, $pos_b, \%cache);
	my $transcript_b = $rf_b->{name} || die;

	my $pair = sprintf '%s/%s/%s-%s/%s/%s',
	  $gene_a, $transcript_a, $feature_a,
	    $gene_b, $transcript_b, $feature_b;
	$pairs{$pair} = 1;
      }
    }
  } elsif ($rfs_a) {
    #
    # records for gene A only:
    #
    my $feature_b = "";
    my $transcript_b = "";

    foreach my $rf_a (@{$rfs_a}) {
      next unless rf_overlap($rf_a, $chr_a, $pos_a);
      my $transcript_a = $rf_a->{name} || die;
      my $feature_a = get_feature_tag($rf, $rf_a, $pos_a, \%cache);
      my $pair = sprintf '%s/%s/%s-%s/%s/%s',
	$gene_a, $transcript_a, $feature_a,
	  $gene_b, $transcript_b, $feature_b;
      $pairs{$pair} = 1;
    }
  } elsif ($rfs_b) {
    #
    # records for gene B only:
    #
    my $transcript_a = "";
    my $feature_a = "";
    foreach my $rf_b (@{$rfs_b}) {
      # build pair info
      next unless rf_overlap($rf_b, $chr_b, $pos_b);
      my $feature_b = get_feature_tag($rf, $rf_b, $pos_b, \%cache);
      my $transcript_b = $rf_b->{name} || die;

      my $pair = sprintf '%s/%s/%s-%s/%s/%s',
	$gene_a, $transcript_a, $feature_a,
	  $gene_b, $transcript_b, $feature_b;
      $pairs{$pair} = 1;
    }
  }

  $row->{gene_pair_summary} = join ",", sort keys %pairs;
}

sub rf_overlap {
  my ($rfr, $chr_query, $pos_query) = @_;
  $chr_query = cook_chromosome_name($chr_query);
  my $chr_rf = cook_chromosome_name($rfr->{chrom} || die);

  my $overlap;

  if ($chr_query eq $chr_rf) {
    my $txStart = ($rfr->{txStart} || die) + 1;
    my $txEnd = $rfr->{txEnd} || die;
    die unless $txStart < $txEnd;
    # should be genomic 5'-to-3' regardless of strand
    $overlap = 1 if $pos_query >= $txStart and $pos_query <= $txEnd;
  }
  return $overlap;
}

sub conditional_noncoding_filter {
  my ($set_ref) = @_;

  my @nm;
  my @other;
  # NR or other accession types

  foreach my $rfr (@{$set_ref}) {
    my $acc = $rfr->{name} || die;
    if ($acc =~ /^NM_/) {
      push @nm, $rfr;
    } else {
      push @other, $rfr;
    }
  }

#  printf STDERR "raw:%d NM:%d other:%d\n", scalar(@{$set_ref}), scalar(@nm), scalar @other;

  if (@nm and @other) {
#    printf STDERR "keep %d NM_, discard %d other\n", scalar(@nm), scalar @other;
    @{$set_ref} = @nm;
  }

}


sub patch_breakpoints {
  # hack to patch genomic breakpoint info from intermediate files
  # back into curated pattern file
  my $f_in = $FLAGS{"patch-breakpoints"} || die;
  # file to patch

  my @pattern_src_files = qw(infiles.txt infiles_rna.txt);
  # symlinks
  my $src_dir = dirname(realpath($pattern_src_files[0]));

  #
  #  load pre-breakpoint annotation extended pattern sequences:
  #
  my %index;
 LOAD:
  foreach my $listfile (@pattern_src_files) {
    my $files = read_simple_file($listfile);
    foreach my $src_raw (@{$files}) {
      my $f_pattern = sprintf "%s/%s.extended.tab.extended_500.pattern.tab",
	$src_dir, basename($src_raw);
      # intermediate output file with extended sequence
      die "where is $f_pattern" unless -s $f_pattern;

      if (0 and %index) {
	printf STDERR "DEBUG, quitting load\n";
	last LOAD;
      }
      printf STDERR "loading %s\n", $f_pattern;

      my $df = new DelimitedFile(
				 "-file" => $f_pattern,
				 "-headers" => 1,
				);
      while (my $row = $df->get_hash()) {
	die unless exists $row->{extended_500};
	my $seq = $row->{extended_500};
	if ($seq) {
	  # not every row will have one
	  $seq = uc($seq);
	  $row->{src_file} = $f_pattern;
	  push @{$index{$seq}}, $row;
	}
      }
    }
  }

  #
  #  map final patterns back to intermediate row(s)
  #
  my $df = new DelimitedFile("-file" => $f_in,
			     "-headers" => 1,
			     );
  my $outfile = basename($f_in) . ".breakpoints.tab";
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => [
					   qw(
					       breakpoints_raw
					       breakpoints_adj
					       breakpoints_src
					    )
					  ],
  			      "-auto_qc" => 1,
			     );

  # while (my $row = $df->next("-ref" => 1)) {  # headerless

  while (my $row = $df->get_hash()) {
    my $pattern = $row->{sequence} || die;
    # fz2 pattern sequence, required

    $pattern =~ tr/[]{}//d || die;
    # remove (required) fz2 breakpoint markers
    $pattern = uc($pattern);
    # just to be sure

    my %bp_raw;
    my %bp_adj;
    my %src;
    if (my $hits = $index{$pattern}) {
      foreach my $h (@{$hits}) {
	my $chrA = $h->{chrA} || dump_die($row, "no chrA");
	my $chrB = $h->{chrB} || die;
	my $posA_raw = $h->{posA} || die;
	my $posB_raw = $h->{posB} || die;
	my $posA_adj = $h->{genea_pos_adj} || die;
	my $posB_adj = $h->{geneb_pos_adj} || die;

	my $breakpoint_raw = sprintf '%s:%d-%s:%d',
	  $chrA, $posA_raw,
	    $chrB, $posB_raw;
	# TO DO: add genome too?
	$bp_raw{$breakpoint_raw} = 1;

	my $breakpoint_adj = sprintf '%s:%d-%s:%d',
	  $chrA, $posA_adj,
	    $chrB, $posB_adj;
	$bp_adj{$breakpoint_adj} = 1;

	my $f_src = $h->{src_file} || die;
	$src{$f_src} = 1;
      }
    }

    $row->{breakpoints_raw} = join ",", sort keys %bp_raw;
    # this might be cleaner if adjusting for minor isoforms clutters the output

    $row->{breakpoints_adj} = join ",", sort keys %bp_adj;
    # OTOH, for RNA-based data this could potentially be cleaner

    $row->{breakpoints_src} = join ",", sort keys %src;
    $rpt->end_row($row);
  }

  $rpt->finish();
}

sub soft_clip_to_pattern {
  # generate fz2 pattern based on soft clipping in a SAM read
  # TO DO:
  #  - other entrypoint possibilities
  #  - deal with strand: this example assumes both aligned and clipped are +
  my $f_sam = $FLAGS{sam} || die "-sam";
  open(IN, $f_sam) || die;
  my $sam = <IN>;
  chomp $sam;
  my @f = split /\t/, $sam;
  my ($chr, $pos, $cigar, $sequence) = @f[2,3,5,9];

  my @cigars = $cigar =~ /(\d+\S)/g;

  my $cigar_softclip = $cigars[0];
  $cigar_softclip =~ /^(\d+)S$/ || die "first CIGAR entry in $cigar not softclip";
  my $sc_length = $1;
  # TO DO: also work with examples where soft-clip is at the end

  my $cigar_match = $cigars[1];
  $cigar_match =~ /^(\d+)[\=M]$/ || die "second CIGAR entry in $cigar not match";
  my $match_length = $1;

  my $regexp = sprintf '^(.{%d})(.{%d})', $sc_length, $match_length;
  $sequence =~ /$regexp/ || die;
  my ($upstream, $downstream) = ($1, $2);

  printf STDERR "up raw: %s\n", $upstream;
  printf STDERR "down raw: %s\n", $downstream;

  $upstream = substr($upstream, - $PATTERN_CONDENSE_DUPLICATE_FLANKING_WINDOW_SIZE);
  $downstream = substr($downstream, 0, $PATTERN_CONDENSE_DUPLICATE_FLANKING_WINDOW_SIZE);

  printf STDERR "up trimmed: %s\n", $upstream;
  printf STDERR "down trimmed: %s\n", $downstream;

  die "build fz2 pattern here";
}

sub itd2fuzzion {
  # generate fz2 patterns for contigs containing ITD vs. reference mRNA
  #
  # inputs:
  #  - contig sequence containing ITD vs. reference
  #  - refSeq accession to use for comparison (STRICT matching required)
  # output: fz2 pattern file
  my $f_in = $FLAGS{itd2fuzzion} || die;
  find_binary("blastn", "-die" => 1);

  my $f_gene = "gene";
  my $f_refseq = "refseq";
  my $f_contig = "contig";
  # make configurable?
  my $restrict_nm = $FLAGS{"restrict-nm-a"};

  my $cicero_mode = $FLAGS{cicero};
  # CICERO ITD call data, breakpoint genes are annotated, but
  # reference accessions are not provided
  my $no_chr_pos = $FLAGS{"no-chr-pos"};

  my ($rf, $f_refflat);
  if ($cicero_mode) {
    $rf = get_refflat() || die;
  }

  printf STDERR "configuration:\n";
  printf STDERR "  BLAST gap open/extend: %s/%s\n",
    (defined $BLAST_GAP_OPEN ? $BLAST_GAP_OPEN : "[default]"),
      (defined $BLAST_GAP_EXTEND ? $BLAST_GAP_EXTEND : "[default]");
  printf STDERR "  treat large dup as fusion: %s\n", yes_no($ITD_LARGE_DUP_TREAT_AS_FUSION);
  printf STDERR "  use 5' UTR in alignment? %s\n", yes_no($ITD_INCLUDE_UTR5_IN_MAPPING);
  printf STDERR "  merge output rows by pattern?: %s", yes_no($ITD_MERGE_PATTERN_OUTPUT);
  printf STDERR " (type %d)", $ITD_MERGE_PATTERN_OUTPUT if $ITD_MERGE_PATTERN_OUTPUT;
  print STDERR "\n";

  print STDERR "\n";

  my $force_gene = $FLAGS{gene};
  my $force_acc = $FLAGS{accession};
  # these apply to entire file, so only appropriate for single-record
  # or single-gene input

  my @f_check = $f_contig;
  push @f_check, $f_gene unless $force_gene or $cicero_mode;
  push @f_check, $f_refseq unless $force_acc or $cicero_mode;

  #
  #  pass 1: identify referenced refSeqs
  #
  my @in_rows;
  # load input into memory as it requires 2 passes, and in cicero
  # mode some annotations will be generated and added

  my $df = new DelimitedFile("-file" => $f_in,
			     "-headers" => 1,
			     );
  my %needed_rf;
  my %exceptions;

  my $f_cicero_exception = "__cicero_exception___";
  my $f_cicero_accessions = "__cicero_accessions__";
  my $f_cicero_gene = "__cicero_gene__";

  while (my $row = $df->get_hash()) {
    push @in_rows, $row;
    foreach my $f (@f_check) {
      dump_die($row, "no data for field $f") unless $row->{$f};
    }

    if ($cicero_mode) {
      my $geneA = get_gene($row, "A");
      if (not($geneA) or $geneA eq VALUE_NA) {
	$row->{$f_cicero_exception} = "geneA_not_valid";
	next;
      }

      my $geneB = get_gene($row, "B");
      if (not($geneB) or $geneB eq VALUE_NA) {
	$row->{$f_cicero_exception} = "geneB_not_valid";
	next;
      }

      if ($geneA eq $geneB) {
	$row->{$f_cicero_gene} = $geneA;
	my $acc = $FLAGS{accession};
	my @accs;
	if ($acc) {
	  # manually-specified, e.g. EGFRvIII, UBTF
	  @accs = $acc;
	} else {
	  if (my $hits = $rf->find_by_gene($geneA)) {
	    @accs = grep {/^NM_/} map {$_->{name}} @{$hits};
	    die "no usable accessions" unless @accs;
	  } else {
	    die "no refFlat matches for $geneA";
	  }
	}
	$row->{$f_cicero_accessions} = \@accs;
	foreach my $acc (@accs) {
	  $needed_rf{$acc} = 1;
	}
      } else {
	# row not usable, record exception for later
	$row->{$f_cicero_exception} = "multiple_genes_cannot_evaluate_duplication";
	next;
      }
    } else {
      my $acc = $force_acc || $row->{$f_refseq} || die;
      $needed_rf{$acc} = 1;
      # ugh: complications if given accession is versioned!!
    }
  }

  #
  #  fetch GenBank records for needed refseqs:
  #
  my $gb_cache = new GenBankCache();
  my $result_files = $gb_cache->get([sort keys %needed_rf]);

  my %acc2cds;
  my %acc2versioned;

  foreach my $file (@{$result_files}) {
    my $gbc = new GenBankCDS("-file" => $file);
    my $info = $gbc->get_info();
    my $cds;
    if ($ITD_INCLUDE_UTR5_IN_MAPPING) {
      # some duplications include an interstitial 5' UTR,
      # e.g. test_query_identity_56.tab
      $cds = $info->{mrna_cds_with_utr5} || die;
    } else {
      $cds = $info->{mrna_cds} || die;
      die "this is not compatible with feature annotation!!";
    }

    my $acc_v = $info->{accession_versioned} || die;
    my $acc_unv = $info->{accession} || die;
    foreach my $acc ($acc_v, $acc_unv) {
      # accession numbers in the input file (if provided),
      # may or may not be versioned
      $acc2cds{$acc} = $cds;
      $acc2versioned{$acc} = $acc_v;
    }
  }

  my $f_out = $FLAGS{out} || basename($f_in) . ".fz2.tab";
#  (my $f_failed = $f_out) =~ s/\.fz2\.tab$/.fz2_failed.tab/ || die;
  # can't assume suffix if -out is manually specified
  my $f_failed = $f_out . ".failed.tab";
  # so, just append

  my @h_shared = (
		  @H_FZ2_CORE,
		  qw(
		      blast_hit_count
		      blast_hsp_count
		      blast_query_coverage_fraction
		   ),
		  $F_ITD_FEATURES,
		  qw(
		      exception
		   )
		  # these fields specific to this pipeline
		 );
  # shared between main output and failure report

  my $rpt = new Reporter(
			 "-file" => $f_out,
			 "-delimiter" => "\t",
			 "-labels" => [
				       @h_shared,
				       "detected_repeat",
				       $f_contig,
#				       "gene_sequence"
				       ]
			);
  # fz2 output
  # TO DO: standard variable(s) for pattern file header fields?

  my $f_sample = guess_sample_column($df);

  #
  #  pass 2: process
  #
  my $rpt_failed = $df->get_reporter("-file" => $f_failed,
				     "-extra" => [
						  @h_shared,
						  "any_result_ok",
						 ]
				    );

  my $genome = standardize_genome($FLAGS{genome} || die "-genome");

  my $rc_mode = $FLAGS{rc};

  my %pattern_counter;
  foreach my $row_in (@in_rows) {
    my $contig = $row_in->{$f_contig} || die;
    printf STDERR "new row: %s\n", $contig;
    $contig = reverse_complement($contig) if $rc_mode;

    my $copy_row = sub {
      my %r = %{$row_in};
      # local copy of row for modifications and additions.  don't modify
      # original row hashref because output is queued rather than written
      # directly.  Also can't simple make the copy local to the row,
      # because there will be a new output row for each isoform checked.
      $r{blast_hit_count} = 0;
      $r{blast_hsp_count} = 0;
      $r{blast_query_coverage_fraction} = "";
      $r{genea_symbol} = $r{geneb_symbol} = "";
      $r{genea_acc} = $r{geneb_acc} = "";
      $r{genea_feature} = $r{geneb_feature} = "";
      # blank/null values for exceptions
      if ($no_chr_pos) {
	# not available in input
	$r{genea_chr} = "";
	$r{geneb_chr} = "";
	$r{genea_pos} = "";
	$r{geneb_pos} = "";
      } elsif ($cicero_mode) {
	$r{genea_chr} = $r{chrA} || die "no chrA";
	$r{geneb_chr} = $r{chrB} || die "no chrB";
	$r{genea_pos} = $r{posA} || die "no posA";
	$r{geneb_pos} = $r{posB} || die "no posB";
      } else {
	die "chr/pos not yet implemented for non-cicero input; see also -no-chr-pos";
      }

      $r{$F_PATHOGENICITY_SOMATIC} = "" unless exists $r{$F_PATHOGENICITY_SOMATIC};
      # only present in some input files

      $r{$F_GENOME} = $genome;
      $r{$F_SOURCE} = file_to_tag($f_in);
      $r{$f_contig} = $contig;
      $r{$F_SAMPLE} = $f_sample ? $r{$f_sample} : "";
      # populatable annotations: ADD MORE

      return \%r;
    };

    my @out_ok;
    my @out_failed;
    # all records for the input row, passed or failed

    my @accs;
    if ($cicero_mode) {
      if (my $exception = $row_in->{$f_cicero_exception}) {
	my $r = &$copy_row();
	$r->{exception} = $exception;
	$r->{$F_ITD_FEATURES} = "";
	push @out_failed, $r;
      } else {
	my $accs = $row_in->{$f_cicero_accessions} || dump_die($row_in, "no field $f_cicero_accessions");
	push @accs, @{$accs};
      }
    } else {
      my $acc = $force_acc || $row_in->{$f_refseq} || die;
      push @accs, $acc;
    }

    foreach my $acc (@accs) {

      next if $restrict_nm and $acc ne $restrict_nm;
      # typically debugging

      my $r = &$copy_row();

      unless ($acc =~ /\.\d+$/) {
	$acc = $acc2versioned{$acc} || die "can't get version for $acc";
      }

      my $gbc = new GenBankCDS("-file" => $gb_cache->get($acc, "-single" => 1));

      my $gene;
      if ($force_gene) {
	$gene = $force_gene;
      } else {
	my $f = $cicero_mode ? $f_cicero_gene : $f_gene;
	$gene = $r->{$f} || dump_die($r, "no $f");
      }

      $r->{genea_symbol} = $r->{geneb_symbol} = $gene;
      $r->{genea_acc} = $r->{geneb_acc} = $acc;
      $r->{genea_feature} = $r->{geneb_feature} = "";

#      dump_die($r, "new acc $acc", 1);

      my $rs_seq = $acc2cds{$acc} || die "can't get cds for $acc";

      my $blast = get_blast();
#      die $blast->gapextend;
      my $parser = $blast->blast(
				 "-query" => {
					      "contig" => $contig,
					     },
				 "-database" => {
						 $acc => $rs_seq
						}
				);
      my $result = $parser->next_result();
      # one object per query sequence (only one query seq)

      my @hits;
      while (my $hit = $result->next_hit()) {
	push @hits, $hit;
      }
      $r->{blast_hit_count} = scalar @hits;
      $r->{blast_hsp_count} = 0;

      unless (@hits) {
	$r->{exception} = "no_BLAST_matches";
	$r->{$F_ITD_FEATURES} = "";
	push @out_failed, $r;
	next;
      }

      dump_die($r, "no blast hits for contig") unless @hits;
      die "not exactly one hit" unless @hits == 1;

      #
      #  gather raw HSPs and calculate query coverage percentage:
      #
      my @hsps_raw;
      my @spans_q;
      my @si_subject;
      my @spans_s;
      while (my $hsp = $hits[0]->next_hsp) {
	push @hsps_raw, $hsp;
      }
      @hsps_raw = sort {$a->start("query") <=> $b->start("query")} @hsps_raw;
      # sort HSPs by query start, so references below are
      # consecutive
      # TO DO:
      # 1. sort by query END rather than start?
      # 2. detect very large ITD, where right hsp hit position is
      #    less than left hsp hit position

      $r->{blast_hsp_count} = scalar @hsps_raw;

      my $has_gaps;
      # TO DO: maybe refine by methods that choose subset of HSPs?

      foreach my $hsp (@hsps_raw) {
	printf STDERR "    score:%s strand:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s query_span:%d ref_span:%d total_span=%d query_string=%s hit_string=%s\n",
	  $hsp->score,
	    $hsp->strand,
	      $hsp->range("query"),
		$hsp->range("hit"),
		  $hsp->num_identical(),
		    $hsp->frac_identical("query"),
		      $hsp->length("query"),
			$hsp->length("hit"),
			  $hsp->length("total"),
			    $hsp->query_string(),
			      $hsp->hit_string() if $VERBOSE;
	push @spans_q, sprintf "%d-%d", $hsp->range("query");
	my $span_subject = sprintf "%d-%d", $hsp->range("subject");
	push @spans_s, $span_subject;
	push @si_subject, new Set::IntSpan($span_subject);

	my $hi = new HSPIndexer("-hsp" => $hsp);
	$has_gaps = 1 if $hi->get_gap_info("query");
	$has_gaps = 1 if $hi->get_gap_info("subject");
      }

      my $si = new Set::IntSpan(@spans_q);
      my $query_aligned_bases = size $si;
      my $query_aligned_fraction = $query_aligned_bases / length($contig);
      $r->{blast_query_coverage_fraction} = sprintf "%.3f", $query_aligned_fraction;
      $r->{gene_sequence} = $rs_seq;

      my @itd_features;

      push @itd_features, "has_hsp_gaps" if $has_gaps;

      if ($query_aligned_fraction < $ITD_MIN_QUERY_ALIGN_FRACTION) {
	push @itd_features, "low_query_alignment_identity";
	# just report this as a tag rather than rejecting,
	# some cases may be salvagable anyway
      }

      my @sub_opts = (
		      "-hsps" => \@hsps_raw,
		      "-query-sequence" => $contig,
		      "-reference-sequence" => $rs_seq,
		      "-row" => $r,
		      "-features" => \@itd_features,
		      "-acc" => $acc,
		      "-gbc" => $gbc,
		     );

      my $ok;

      if (simple_itd_check(@sub_opts)) {
	$ok = 1;
      } elsif (large_duplication_check(@sub_opts)) {
	# found implied large ITD, i.e. where HSPs are ordered by
	# query position, but subject sequences are out of order
	$ok = 1;
      } elsif (pure_insertion_check(@sub_opts)) {
	$ok = 1;
      } elsif (pure_deletion_check(@sub_opts)) {
	$ok = 1;
      } elsif (complex_indel_check(@sub_opts)) {
	$ok = 1;
      } elsif (0 and subject_overlap_check(@sub_opts,
				     "-spans-s" => \@spans_s,
				     "-si-subject" => \@si_subject,
				    )) {
	# disabled, these cases should now all handled by better code above
	die "TEST ME: subject overlap check found";
	# other code and/or prerequisites may handle all these cases now
	$ok = 1;
      }

      if (not($ok) and @hsps_raw == 1) {
	#
	#  report additional info if just one alignment.
	#  Maybe something we missed?  A small ITD, or some kind
	#  of indel not captured?:
	#
	my $hi = new HSPIndexer("-hsp" => $hsps_raw[0]);
	push @itd_features, sprintf "query_gaps=%d",
	  $hi->get_gap_info("query");
	push @itd_features, sprintf "subject_gaps=%d",
	  $hi->get_gap_info("subject");
      }

      $r->{$F_ITD_FEATURES} = join ",", @itd_features;
      $r->{detected_repeat} = "" unless $r->{detected_repeat};

      if ($ok) {
	dump_die($r, "ERROR: no genea_feature") unless $r->{genea_feature};
	dump_die($r, "ERROR: no geneb_feature") unless $r->{geneb_feature};
	# QC to ensure populated by different methods

	$r->{exception} = "";
	$r->{$F_PATTERN_ID} = generate_pid($gene, \%pattern_counter);
	push @out_ok, $r;

	if ($ITD_LARGE_DUP_GENERATE_TWO_PATTERN_STYLES and
	    $r->{$F_INTERNAL_ALT_PATTERN}) {
	  # alternate output row with different pattern is present
	  my $r_alt = $r->{$F_INTERNAL_ALT_PATTERN};
	  $r_alt->{detected_repeat} = "" unless $r_alt->{detected_repeat};
	  $r_alt->{exception} = "";
	  $r_alt->{$F_PATTERN_ID} = generate_pid($gene, \%pattern_counter);
	  push @out_ok, $r_alt;
	}
      } else {
	$r->{exception} = "no_event_detected";
	push @out_failed, $r;
      }
    } # accession

    #
    #  write results to either main or failed file.  for failures, set
    #  flag indicating whether any results could be generated for the
    #  input (e.g. another isoform)
    #

    if ($ITD_MERGE_PATTERN_OUTPUT) {
      # bucket rows by pattern:
      my %pattern2row;
      foreach my $r (@out_ok) {
	my $p = $r->{sequence} || die;
	push @{$pattern2row{$p}}, $r;
      }
      my @out_new;

      if ($ITD_MERGE_PATTERN_OUTPUT == 1) {
	# merge annotation columns, however this will mangle
	# paired information like accession/transcript/exon

	my @f_merge = @h_shared;
	# add more?

	# merge annotations in each set:
	foreach my $p (sort keys %pattern2row) {
	  my $set = $pattern2row{$p};
	  my %r = %{$set->[0]};

	  foreach my $f (@f_merge) {
	    $r{$f} = join ",", @{unique_ordered_list($set, "-key" => $f)};
	    # condense annotations
	  }
	  push @out_new, \%r;
	}

      } elsif ($ITD_MERGE_PATTERN_OUTPUT == 2) {
	# just keep one record per pattern sequence.
	# crude, but preserves separation of transcript/exon annotations
	foreach my $p (sort keys %pattern2row) {
	  my $set = $pattern2row{$p};
	  push @out_new, $set->[0];
	}
      } else {
	die "must be 1 or 2";
      }
      @out_ok = @out_new;
    }

    foreach my $r (@out_ok) {
      if (1) {
	# strip "chr" prefixes if present, so if we happen to
	# later merge patterns from sources where some use the prefix
	# and some don't, we don't get output like "chr13,13"
	foreach my $f (qw(genea_chr geneb_chr)) {
	  die unless exists $r->{$f};
	  $r->{$f} =~ s/^chr//i;
	}
      }
      generate_pair_summary($r);
      $rpt->end_row($r);
    }

    my $any_ok = @out_ok ? 1 : 0;
    foreach my $r (@out_failed) {
      generate_pair_summary($r);
      $r->{$F_PROCESSING_TYPE} = "";
      $r->{$F_PATTERN_ID} = "";
      $r->{$F_PATTERN_SEQUENCE} = "";

      $r->{any_result_ok} = $any_ok;
      # to help distinguish expected failures, e.g. for incompatible
      # isoforms, from other problems that prevent any result from
      # being generated
      $rpt_failed->end_row($r);
    }

  } # input row


  $rpt->finish();
  $rpt_failed->finish();

}

sub query_overlap_rescue {
  my ($hi_left, $hi_right, $features, $gbc) = @_;
  confess "need features" unless $features;
  # gbc: may ultimately need left and right

  my $needs_rescue;
  my $rescue_type;
  my @rescue_tags;

  my ($l_q_start, $l_q_end) = $hi_left->range("query");
  my $l_q_si = new Set::IntSpan(sprintf "%d-%d", $l_q_start, $l_q_end);
  my ($r_q_start, $r_q_end) = $hi_right->range("query");
  my $r_q_si = new Set::IntSpan(sprintf "%d-%d", $r_q_start, $r_q_end);

  # query overlap length
  my $intersect = intersect $l_q_si $r_q_si;
  if ($intersect->size()) {
    #
    # query regions overlap, rescue needed
    #
    $needs_rescue = 1;

    my ($qo_start, $qo_end) = si2range($intersect);
    printf STDERR "query overlap: %d-%d\n", $qo_start, $qo_end;

    my $olen = ($qo_end - $qo_start) + 1;
    push @rescue_tags, sprintf "query_overlap_length=%d", $olen;

    if ($l_q_end == $qo_end and
	$r_q_start == $qo_start) {
      # the contested region is at the end of the left query
      # and the start of the right query.
      #
      # compare alignment quality between the two: if one is
      # perfect and the other isn't, trim the imperfect one.
      # Since the left and right HSPs represent the anchoring
      # to the reference sequence, this will maximize the
      # upstream/downstream perfect alignment.
      my $id_left = $hi_left->get_query_range_identity_score($qo_start, $qo_end);
      my $id_right = $hi_right->get_query_range_identity_score($qo_start, $qo_end);

      if ($id_left == 1 and $id_right < 1) {
	# left alignment perfect, right is not,
	# trim right query alignment.
	$hi_right->trim_query_region($qo_start, $qo_end);
#	printf STDERR "right HSP after trim: qs=%d hs=%d\n", $hi_right->start("query"), $hi_right->start("subject");
	$rescue_type = "left";
      } elsif ($id_right == 1 and $id_left < 1) {
	# right alignment perfect, left is not,
	# trim left query alignment.
	$hi_left->trim_query_region($qo_start, $qo_end);
	# trim the rightmost side of the left alignment
	$rescue_type = "right";
      } elsif ($id_right == 1 and
	       $id_left == 1) {
	# both sides have perfect alignment
	if ($gbc) {
	  my $edge_trimmed = microhomology_feature_edge_check(
							      $hi_left,
							      $hi_right,
							      $qo_start,
							      $qo_end,
							      $gbc);
	  if ($edge_trimmed) {
	    # microhomology can be resolved by using known exon boundaries
	    $rescue_type = $edge_trimmed;
	  }

	}

	# since all genes are normalized to + strand, preserve the
	# left (upstream) alignment.  Any better handling?
	unless ($rescue_type) {
	  # unless microhomology
	  $hi_right->trim_query_region($qo_start, $qo_end);
	  if ($olen <= $ITD_AMBIGUITY_RESCUE_MAX_MICROHOMOLOGY_SIZE) {
	    # small, perfect alignment on both sides: likely microhomology
	    $rescue_type = "microhomology_left";
	  } else {
	    $rescue_type = "tie_left";
	    # arbitary
	  }
	}
      } else {
	# neither alignment is clean in the overlap region.  in this case:
	# - find the portion of the overlap that cleanly aligns
	#   in the left HSP (no mismatches, gaps, etc.)
	# - trim left HSP after this site
	# - trim right HSP before this site

	my $last_ok = $hi_left->get_last_clean_alignment_base($qo_start, $qo_end);
	$hi_left->trim_query_region($last_ok + 1, $qo_end);
	# trim left HSP after this site
	$hi_right->trim_query_region($qo_start, $last_ok);
	# trim right HSP up to and including this site

	$rescue_type = "split";
	push @{$features}, sprintf "query_overlap_split_cleanup_at_query=%d", $last_ok;

#	hsp_debug($hi_left, $hi_right);
      }
    }

    # is the contested region the end of the left hsp and the
    # start of the right hsp?  If so, is one alignment perfect
    # and the other imperfect?  If so, trim the imperfect one.
  }

  if ($needs_rescue) {
    if ($rescue_type) {
      push @rescue_tags, sprintf "query_overlap_rescue=%s", $rescue_type;
    }
  }
  push @{$features}, @rescue_tags;

  return $needs_rescue;
}

sub si2range {
  my ($intersect) = @_;
  my $rl = run_list $intersect;
  my ($start, $end);
  if ($rl =~ /^(\d+)\-(\d+)$/) {
    $start = $1;
    $end = $2;
  } elsif ($rl =~ /^(\d+)/) {
    $start = $end = $1;
  } else {
    die "intersection $rl is not a simple range";
  }
  return ($start, $end);
}

sub large_duplication_check {
  #
  #  handles cases where the subject alignment for the right HSP
  #  starts before the left HSP (very large duplication), or overlaps it
  #  (medium duplication where alignment is captured by the HSPs).
  #
  my (%options) = @_;
  my $hsps = $options{"-hsps"} || die;
  my $rs_seq = $options{"-reference-sequence"} || die;
  my $query_seq = $options{"-query-sequence"} || die;
  my $features = $options{"-features"} || die;
  my $row = $options{"-row"} || die;
  my $acc = $options{"-acc"} || die;
  my $gbc = $options{"-gbc"} || die;

#  my $gb_cache = new GenBankCache();
#  my $gbc = new GenBankCDS("-file" => $gb_cache->get($acc, "-single" => 1));

  my $i_last = scalar @{$hsps} - 1;
  my @possible_i;

  my $found_large_itd = 0;

  for (my $i = 0; $i < $i_last; $i++) {
    if ($hsps->[$i]->end("hit") > $hsps->[$i + 1]->start("hit")) {
      # left HSP ends at a point in the subject sequence *after*
      # the right HSP begins.  This implies a large repeat of an
      # earlier part of the subject sequence.
      push @possible_i, $i;
    }
  }

  if (@possible_i) {
    if (@possible_i == 1) {
      # one candidate pair
      my ($i) = @possible_i;
      my $hi_left = new HSPIndexer("-hsp" => $hsps->[$i]);
      my $hi_right = new HSPIndexer("-hsp" => $hsps->[$i + 1]);

      query_overlap_rescue($hi_left, $hi_right, $features, $gbc);

      # TO DO:
      # additional trimming if large gaps still present?

      my $q_gap = get_query_gap($hi_left, $hi_right);
      push @{$features}, sprintf "query_gap=%d", $q_gap if $q_gap;

      my $ok = 1;

      if ($q_gap > $ITD_MAX_QUERY_GAP) {
	# *might* be this type of event, but too much of a gap
	# in the query sequence vs. the reference, which might
	# indicate the contig is from a different isoform that
	# is not compatible with the current one (i.e. including
	# an exon that this isoform doesn't)
	# set some kind of status flag in output anyway?
	# maybe some isoforms would work and others not?
	#
	# OTOH, might also be a fragment of intronic sequence which
	# won't align to *any* isoform.

	unless ($ITD_LARGE_DUP_TREAT_AS_FUSION or
		$ITD_LARGE_DUP_GENERATE_TWO_PATTERN_STYLES) {
	  push @{$features}, sprintf "possible_large_dup_but_query_gap_too_long";
	  # this isoform isn't a good match, but maybe another one is
	  $ok = 0;
	}
      }

      if ($ok) {
	#
	# any gap in query alignment at breakpoint is within tolerances,
	# consider plausible large duplication
	#
	query_overlap_rescue($hi_left, $hi_right, $features, $gbc);

	my $left_subject_end = $hi_left->end("subject");
	# reference sequence from this point and upstream
	my $right_subject_start;
	# reference sequence from this point and downstream

	my $q_i_start = $hi_left->end("query") + 1;
	# interstitial start: base after left HSP query align ends.
	# the first/reference copy of the repeat should end before
	# this position.
	#
	# TO DO: alternate behavior if we want to place BOTH copies
	# of the ITD sequence within the brackets.  This may be
	# the preferred behavior to avoid false-positive "strong"
	# matches based on an overlap with a single copy only.

	my $q_i_end;
	# interstitial end
	# both base numbers, not indices

	my $implied_itd_size = -1;

	# does the left HSP's subject end site appear within the right HSP?
	# If so the ITD is captured within the alignment and should
	# be reflected in the pattern (see test_ITD_66_nt.tab):
	my $left_subj_start = $hi_left->start("subject");
#	my $left_subj_end = $hi_left->end("subject");
	my $right_subj_start = $hi_right->start("subject");
	my $right_subj_end = $hi_right->end("subject");
	my $is_large_dup;

	itd_feature_annotate($row, $gbc, $hi_left, $hi_right);
	# annotate breakpoints

	printf STDERR "trimmed:\n";
	printf STDERR "   left: q=%d-%d h=%d-%d\n",
	  $hi_left->start("query"), $hi_left->end("query"),
	    $hi_left->start("subject"), $hi_left->end("subject");
	printf STDERR "  right: q=%d-%d h=%d-%d\n",
	  $hi_right->start("query"), $hi_right->end("query"),
	    $hi_right->start("subject"), $hi_right->end("subject");


	my $allow_fusion_pattern;
	my $copies_wanted = $ITD_BRACKET_BOTH_COPIES ? 2 : 1;
	# how many copies of the dup sequence should be placed
	# within the pattern brackets.  Two safer than just one.
	#
	# TO DO: may need reset/adjustment later for long duplications,
	# as it may not be reasonable to require both copies
	my $blast_rescue_mode;

	my $da = new DuplicationAlignment(
					  "-hi_left" => $hi_left,
					  "-hi_right" => $hi_right,
					  "-subject_sequence" => $rs_seq,
					  "-bracket_copy_count" => $copies_wanted
					 );
	if ($da->has_duplication()) {
	  #
	  # duplicate is represented in HSP pair via subject overlap
	  #

	  # examples:
	  # (A) subject region of HSP fully contained within subject region
	  #  of left HSP:
	  #   test_medium_dup_right_hsp_subset_of_left.tab
	  #  /research/rgs01/home/clusterHome/medmonso/work/steve/2020_06_13_cicero_contig_extension/2021_01_04_all/blast/ITD_refactor/debug_PIK3R1/test.tab
	  #
	  # (B) first copy of repeat ends within second HSP (partial overlap)
	  # example: test_medium_dup_partial_overlap.tab

	  my $dup_sequence = $da->dup_sequence();
	  $row->{detected_repeat} = $dup_sequence;
	  $implied_itd_size = length $dup_sequence;

	  if ($copies_wanted == 2 and
	      $ITD_BRACKET_ONE_COPY_IF_DUP_MINIMUM and
	      $implied_itd_size >= $ITD_BRACKET_ONE_COPY_IF_DUP_MINIMUM) {

	    printf STDERR "before: q=%d-%d subj=%d-%d\n",
	      $da->query_bracket_start_base(),
	      $da->query_bracket_end_base(),
	      $da->subject_upstream_base(),
	      $da->subject_downstream_base();

	    $copies_wanted = 1;

	    $da = new DuplicationAlignment(
					   "-hi_left" => $hi_left,
					   "-hi_right" => $hi_right,
					   "-subject_sequence" => $rs_seq,
					   "-bracket_copy_count" => $copies_wanted,
					  );
	    # recalculate

	    printf STDERR "after: q=%d-%d subj=%d-%d\n",
	      $da->query_bracket_start_base(),
	      $da->query_bracket_end_base(),
	      $da->subject_upstream_base(),
	      $da->subject_downstream_base();


	    my $interstitial_len = ($da->query_bracket_end_base() - $da->query_bracket_start_base()) + 1;

	    if ($ITD_BRACKET_ONE_COPY_RESCUE and
		$interstitial_len >= $ITD_BRACKET_ONE_COPY_RESCUE_MINIMUM) {
	      # additional rescue to help with patterns where
	      # interstitial sequence is too long for short-read sequencing
	      $da->blast_rescue();
	      $blast_rescue_mode = 1;
	    }
	  }

	  $q_i_start = $da->query_bracket_start_base();
	  $q_i_end = $da->query_bracket_end_base();
	  # bracketed bases in query sequence

	  $left_subject_end = $da->subject_upstream_base();
	  $right_subject_start = $da->subject_downstream_base();
	  # upstream and downstream reference sequence

	  my $ptype = "medium_subject_dup_bracket_";
	  if ($copies_wanted == 1) {
	    $ptype .= "1_copy";
	    $ptype .= "_blast_rescue" if $blast_rescue_mode;
	  } elsif ($copies_wanted == 2) {
	    $ptype .= "2_copies";
	  } else {
	    die;
	  }
	  $row->{$F_PROCESSING_TYPE} = $ptype;

	  $allow_fusion_pattern = 0;
	  # if we allow this, resulting patterns are quite vulnerable
	  # to false positives especially for "strong-" type read
	  # pairs which don't overlap the site, so you're essentially
	  # just reporting matches to canonical sequence
	  #
	  # TO DO: need some kind of rescue here for long dup sequences.
	  # if we insist on both copies this may be too long for the
	  # pattern to work in short-read data.
	} else {
	  #
	  #  large duplication not completely captured in the
	  #  alignment.  In this situation a secondary, fusion-style
	  #  pattern may also be useful (e.g. for a large duplication
	  #  where a much earlier exon repeats)
	  #
	  $is_large_dup = 1;
	  $row->{$F_PROCESSING_TYPE} = "large_dup";
	  $right_subject_start = $hi_right->start("subject");
	  # here the right alignment represents the start of
	  # the duplicated sequence
	  $q_i_end = $hi_right->start("query") - 1;
	  # base number before start of repeat alignment.
	  # for large ITDs we won't attempt to contain the entire
	  # interstitial sequence within the pattern, because
	  # it can easily be longer than most short reads are
	  # able to capture.
	  $implied_itd_size = $left_subject_end - $right_subject_start;
	  $allow_fusion_pattern = 1;
	}

	die "allow_fusion_pattern not defined" unless defined $allow_fusion_pattern;

	my $ref_up = substr($rs_seq, 0, $left_subject_end);
	# reference sequence from start to end of left HSP

	my $ref_down = substr($rs_seq, $right_subject_start - 1);
	# reference sequence from start of right HSP to end

	trim_flanking(\$ref_up, \$ref_down);

	unless ($ref_up) {
	  # no upstream sequence in final pattern, e.g.
	  # SJHM031438_D1 / IKZF1 / NM_001291840.1
	  # this is not a legal fz2 format, so don't process.
	  foreach my $f ($F_PROCESSING_TYPE, "detected_repeat") {
	    delete $row->{$f};
	  }
	  push @{$features}, "error_no_upstream_sequence";
	  # can't generate a legal fz2 pattern because no upstream
	  # sequence available (may be just for a particular isoform)
	  return 0;
	}

	my $interstitial = "";
	# actually the *bracketed* sequence:
	# - if a medium dup, one or two copies of the duplicated sequence
	#   (depending on policy) plus any interstitial sequence
	# - if a large dup, the truly interstitial/unmappable sequence
	if ($q_i_start <= $q_i_end) {
	  $interstitial = substr($query_seq,
				 $q_i_start - 1,
				 ($q_i_end - $q_i_start) + 1);
	}

	my $fz2_pattern_fusion = sprintf '%s%s%s%s',
	  $ref_up, BREAKPOINT_FUSION_LEFT,
	      BREAKPOINT_FUSION_RIGHT, $ref_down;
	# fusion-style breakpoint, ignoring any interstitial sequence

	my $fz2_pattern_itd = sprintf '%s%s%s%s%s',
	  $ref_up, BREAKPOINT_ITD_LEFT,
	    $interstitial,
	      BREAKPOINT_ITD_RIGHT, $ref_down;
	# ITD-style breakpoint, retaining any interstitial sequence

	push @{$features}, sprintf "implied_itd_length=%d", $implied_itd_size;

	if ($ITD_LARGE_DUP_GENERATE_TWO_PATTERN_STYLES) {
	  my $ptype = $row->{$F_PROCESSING_TYPE} || die;
	  $row->{$F_PATTERN_SEQUENCE} = $fz2_pattern_itd;
	  $row->{$F_PROCESSING_TYPE} .= "_ITD_style";
	  # for primary result, use ITD-style pattern

	  if ($allow_fusion_pattern and
	      (length($interstitial) >= $ITD_MIN_INTERSTITIAL_FOR_ALTERNATE_PATTERN)
	     ) {
	    my $dubious_junction = dubious_junction_check($gbc, $rs_seq, $ref_up, $ref_down);
	    # only calculate this when needed as it may cause problems
	    # for smaller ITDs, where this scenario doesn't apply anyway
	    if (not($dubious_junction)) {
	      # all QC passed.
	      # create alternate output row using a supplemental,
	      # fusion-style pattern.  Only do this when there is true
	      # interstitial sequence: for these cases an ITD-style
	      # pattern will match the source sample, yet that sequence
	      # may also prevent matching for similar-but-not-identical
	      # events in other samples.  A minimum interstitial
	      # sequence length is required, if the interstitial
	      # sequence is very short we should use an ITD-style
	      # pattern only.
	      my %alt = %{$row};
	      $alt{$F_PATTERN_SEQUENCE} = $fz2_pattern_fusion;
	      my @feat = @{$features};
	      push @feat, sprintf "ignored_interstitial=%s", $interstitial if $interstitial;
	      $alt{$F_ITD_FEATURES} = join ",", @feat;
	      $alt{$F_PROCESSING_TYPE} = $ptype . "_fusion_style";

	      $row->{$F_INTERNAL_ALT_PATTERN} = \%alt;
	    }
	  }

	} elsif ($ITD_LARGE_DUP_TREAT_AS_FUSION) {
	  # just a fusion-style pattern
	  die "unhandled" unless $allow_fusion_pattern;
	  $row->{$F_PATTERN_SEQUENCE} = $fz2_pattern_fusion;
	  push @{$features}, sprintf "ignored_interstitial=%s", $interstitial if $interstitial;
	} else {
	  # just an ITD-style pattern
	  $row->{$F_PATTERN_SEQUENCE} = $fz2_pattern_itd;
	}

	$found_large_itd = 1;
      }
    } else {
      push @{$features}, "multiple_candidates_for_large_dup";
      # maybe later?
    }
  }

  return $found_large_itd;
}

sub get_query_gap {
  my ($hi_left, $hi_right) = @_;
  # HSP or HSPIndexer
  my $left_q_end = $hi_left->end("query");
  my $right_q_start = $hi_right->start("query");
  return ($right_q_start - $left_q_end) - 1;
}

sub trim_flanking {
  my ($upstream, $downstream) = @_;
#  printf STDERR "up before: %s\ndown before: %s\n", $$upstream, $$downstream;
  if (1) {
    $$upstream = substr($$upstream, - $EXTENDED_CHUNK_LENGTH);
    $$downstream = substr($$downstream, 0, $EXTENDED_CHUNK_LENGTH);
  }
#  printf STDERR "up after: %s\ndown after: %s\n", $$upstream, $$downstream;
}

sub generate_pid {
  my ($gene, $pattern_counter) = @_;
  return sprintf '%s-%s-%02d', $gene, $gene, ++$pattern_counter->{$gene};
}

sub simple_itd_check {
  # check for a simple ITD which appears within a single HSP
  # as a deletion from the reference
  my (%options) = @_;
  my $hsps = $options{"-hsps"} || die;
  my $rs_seq = $options{"-reference-sequence"} || die;
  my $query_seq = $options{"-query-sequence"} || die;
  my $features = $options{"-features"} || die;
  my $row = $options{"-row"} || die;

  #
  #  look for candidate HSPs:
  #
  my @candidate_hi;
  foreach my $hsp (@{$hsps}) {
    my $hi = new HSPIndexer(
			    "-hsp" => $hsp,
			    "-min_gap_length" => $ITD_MIN_SIZE
			   );
    push @candidate_hi, $hi if @{$hi->hit_gaps()};
  }

  my $found_simple_itd;

  if (@candidate_hi) {
    if (@candidate_hi == 1) {
      my $hi = $candidate_hi[0];
      my $gaps = $hi->hit_gaps();
      die "not exactly one gap" unless @{$gaps} == 1;
      my $gap = $gaps->[0];

      my $gbc = $options{"-gbc"} || die;

      my $gb_up = $gap->{upstream_hit_base} || die;
      my $gb_down = $gap->{downstream_hit_base} || die;

      $row->{genea_feature} = $gbc->get_feature_for_base_number($gb_up);
      $row->{geneb_feature} = $gbc->get_feature_for_base_number($gb_down);

      my $ref_up = substr($rs_seq, 0, $gb_up);
      my $ref_down = substr($rs_seq, $gb_down - 1);
      my $interstitial = substr(
			       $query_seq,
			       $gap->{query_start} - 1,
			       $gap->{gap_length}
			       );

      my $found_repeat;
      if ($ITD_BRACKET_BOTH_COPIES) {
	# original approach was simply to put brackets around the
	# inserted sequence, which is too prone to false positives,
	# especially if the flanking sequence is similar.  Really need to
	# have BOTH copies of the dup within the pattern brackets.
	my $has_up = $ref_up =~ /${interstitial}$/i;
	my $has_down = $ref_down =~ /^${interstitial}/i;
	if ($has_up and $has_down) {
	  die "test me: both up and downstream copies";
	} elsif ($has_up) {
	  die "implement me: upstream copy";
	} elsif ($has_down) {
	  $ref_down =~ s/^${interstitial}//i || die "can't trim";
	  $row->{detected_repeat} = $interstitial;
	  $interstitial .= $interstitial;
	  $found_repeat = 1;
	} else {
	  # no perfect neighboring copy
	  my $identical_down = check_identical_down($interstitial, $ref_down);

	  # TO DO:
	  # implement upstream check, when we come across a good example

	  printf STDERR "identical down: %d\n", $identical_down;

	  if ($identical_down >= $ITD_RESCUE_MIN_SIZE) {
	    # expand the interstitial sequence to include the
	    # repeated bases.  Assume any additional bases are
	    # truly interstitial as part of the insertion/dup event.
	    $ref_down =~ s/^(.{$identical_down})// || die;
	    $interstitial .= $1;
	    # move the additional repeated bases from the downstream
	    # sequence into the interstitial sequence
	    push @{$features}, sprintf "rescue_itd_site_downstream=%d", $identical_down;
	    $found_repeat = 1;
	  } else {
	    printf STDERR "WARNING: no second copy detected for interstitial sequence $interstitial\n";
	    # - if event is a simple insertion, this outcome is OK.
	    # - however if there is some kind of unhandled repeat it should
	    #   definitely be addressed.  record a different processing
	    #   type for these cases for future review/fixes
	  }
	}
      }

      die unless $ref_up and $ref_down and $interstitial;
      trim_flanking(\$ref_up, \$ref_down);

      my $fz2_pattern = sprintf '%s%s%s%s%s',
	$ref_up, BREAKPOINT_ITD_LEFT,
	  $interstitial,
	    BREAKPOINT_ITD_RIGHT, $ref_down;

      $row->{$F_PATTERN_SEQUENCE} = $fz2_pattern;
      $found_simple_itd = 1;
      $row->{$F_PROCESSING_TYPE} = $found_repeat ?
	"single_alignment_itd" : "single_alignment_insertion";
    } else {
      die "simple ITD: multiple HSPs with possible candidates!";
      # example?
    }
  }

  return $found_simple_itd;
}


sub pure_insertion_check {
  # event has a portion of the query sequence which is not aligned to
  # the reference at all.  HSPs have no overlap to the reference,
  # with near-abuting coverage.
  my (%options) = @_;
  my $hsps = $options{"-hsps"} || die;
  my $rs_seq = $options{"-reference-sequence"} || die;
  my $query_seq = $options{"-query-sequence"} || die;
  my $features = $options{"-features"} || die;
  my $row = $options{"-row"} || die;
  my $found_pure_insertion = 0;

  if (@{$hsps} >= 2) {
    my @candidate_i;
    my $final = scalar @{$hsps} - 2;
    #
    #  find HSP pairs where subject alignments abut, or nearly so
    #
    for (my $i = 0; $i <= $final; $i++) {
      my $left_ref_end = $hsps->[0]->end("subject");
      my $right_ref_start = $hsps->[$i + 1]->start("subject");

      if ($right_ref_start >= $left_ref_end) {
	# subject intervals are also ordered left to right
	my $subject_gap = $right_ref_start - $left_ref_end - 1;
	die "fix me, negative subject gap" if $subject_gap < 0;
	# maybe hsp overlap?
	if ($subject_gap <= $ITD_PURE_INSERTION_MAX_SUBJECT_GAP) {
	  push @candidate_i, $i;
	}
      }
    }

    if (@candidate_i) {
      die "multiple simple insertion candidate HSP pairs" if @candidate_i > 1;
      my ($i) = @candidate_i;
      my $hi_left = new HSPIndexer("-hsp" => $hsps->[$i]);
      my $hi_right = new HSPIndexer("-hsp" => $hsps->[$i + 1]);

      if (query_overlap_rescue($hi_left, $hi_right, $features)) {
	confess "TEST ME: query overlap for possible simple insertion";
      }

#      hsp_debug($hi_left, $hi_right);
      my $gbc = $options{"-gbc"} || die;
      itd_feature_annotate($row, $gbc, $hi_left, $hi_right);

      my $ref_up = substr($rs_seq,
			  0,
			  $hi_left->end("subject"));
      my $ref_down = substr($rs_seq,
			    $hi_right->start("subject") - 1);

      my $int_start = $hi_left->end("query") + 1;
      my $int_end = $hi_right->start("query") - 1;
      # start/end bases of interstitial sequence

      my $interstitial = substr($query_seq,
				$int_start - 1,
				($int_end - $int_start) + 1);
      trim_flanking(\$ref_up, \$ref_down);

      my $fz2_pattern = sprintf '%s%s%s%s%s',
	$ref_up, BREAKPOINT_ITD_LEFT,
	  $interstitial,
	    BREAKPOINT_ITD_RIGHT, $ref_down;

      $row->{$F_PATTERN_SEQUENCE} = $fz2_pattern;

      $row->{$F_PROCESSING_TYPE} = "insertion_multi_hsp";

      $found_pure_insertion = 1;
    }
  }

  return $found_pure_insertion;
}

sub pure_deletion_check {
  # event has 2 HSPs where query coverage is contiguous, but
  # skips a portion of the subject sequence (i.e., deletion).
  my (%options) = @_;
  my $hsps = $options{"-hsps"} || die;
  my $rs_seq = $options{"-reference-sequence"} || die;
  my $query_seq = $options{"-query-sequence"} || die;
  my $features = $options{"-features"} || die;
  my $row = $options{"-row"} || die;
  my $gbc = $options{"-gbc"} || die;

  my $found_pure_deletion;

  if (@{$hsps} >= 2) {
    my @candidate_i;
    my $final = scalar @{$hsps} - 2;
    #
    #  find HSP pairs where subject alignments abut, or nearly so
    #
    for (my $i = 0; $i <= $final; $i++) {
      my $left_ref_end = $hsps->[0]->end("subject");
      my $right_ref_start = $hsps->[$i + 1]->start("subject");

      my $left_query_end = $hsps->[0]->end("query");
      my $right_query_start = $hsps->[$i + 1]->start("query");

      if ($right_ref_start >= $left_ref_end) {
	# subject intervals are also ordered left to right
	my $query_gap = $right_query_start - $left_query_end - 1;
#	die "fix me, negative query gap" if $query_gap < 0;
	# don't think this is important at this stage, candidates
	# will be trimmed later
	if ($query_gap <= $ITD_PURE_DELETION_MAX_QUERY_GAP) {
	  # query coverage expected to be contiguous or nearly so
	  push @candidate_i, $i;
	}
      }
    }

    if (@candidate_i) {
      die "multiple simple deletion candidate HSP pairs" if @candidate_i > 1;
      my ($i) = @candidate_i;
      my $hi_left = new HSPIndexer("-hsp" => $hsps->[$i]);
      my $hi_right = new HSPIndexer("-hsp" => $hsps->[$i + 1]);

      query_overlap_rescue($hi_left, $hi_right, $features);

      hsp_debug($hi_left, $hi_right);

      #$row->{genea_feature} = $gbc->get_feature_for_base_number($left_subj_end);
      # $row->{geneb_feature} = $gbc->get_feature_for_base_number($right_subj_start);

      itd_feature_annotate($row, $gbc, $hi_left, $hi_right);

      my $ref_up = substr($rs_seq,
			  0,
			  $hi_left->end("subject"));
      my $ref_down = substr($rs_seq,
			    $hi_right->start("subject") - 1);

      my $int_start = $hi_left->end("query") + 1;
      my $int_end = $hi_right->start("query") - 1;
      # start/end bases of interstitial sequence

      my $interstitial = substr($query_seq,
				$int_start - 1,
				($int_end - $int_start) + 1);
      trim_flanking(\$ref_up, \$ref_down);

      my $fz2_pattern = sprintf '%s%s%s%s%s',
	$ref_up, BREAKPOINT_ITD_LEFT,
	  $interstitial,
	    BREAKPOINT_ITD_RIGHT, $ref_down;

      $row->{$F_PATTERN_SEQUENCE} = $fz2_pattern;

      $row->{$F_PROCESSING_TYPE} = "deletion_multi_hsp";

      $found_pure_deletion = 1;
    }
  }

  return $found_pure_deletion;
}

sub subject_overlap_check {
  my (%options) = @_;
  my $hsps = $options{"-hsps"} || die;
  my $rs_seq = $options{"-reference-sequence"} || die;
  my $query_seq = $options{"-query-sequence"} || die;
  my $features = $options{"-features"} || die;
  my $row = $options{"-row"} || die;
  my $si_subject = $options{"-si-subject"} || die;
  my $spans_s = $options{"-spans-s"} || die;

  my $found_subject_overlap = 0;

  #
  #  find cases where multiple HSPs match some subset of the reference
  #  (i.e. the ITD sequence):
  #
  my %intersect;
  for (my $i = 0; $i < @{$hsps}; $i++) {
    for (my $j = 0; $j < @{$hsps}; $j++) {
      if ($i != $j) {
	my $si_i = $si_subject->[$i];
	my $si_j = $si_subject->[$j];
	my $intersect = intersect $si_i $si_j;
	if (size $intersect) {
	  my $key = join "-", sort {$a <=> $b} ($i, $j);
	  $intersect{$key} = [ size $intersect, $intersect ];
	  printf STDERR "intersect %s %s (%s, len=%d, detail=%s)\n", $spans_s->[$i], $spans_s->[$j], $key, size $intersect, run_list $intersect;
	  # TO DO: ignore de minimis overlaps for this purpose?
	}
      }
    }
  }

  my ($hi_left, $hi_right);
  my ($subject_overlap_start, $subject_overlap_end);
  my $subject_overlap_chunk = "";

  if (%intersect) {
    #
    #  duplicated region in subject detected
    #
    my @best = sort {$intersect{$b}->[0] <=> $intersect{$a}->[0]}
      keys %intersect;

    if (@best > 1) {
      # seem dubious
      push @{$features}, sprintf "multiple_candidates_for_subject_overlap";
      # TO DO: add to warnings list instead??
    } elsif (@best == 1) {
      #
      # single candidate.
      # if we sort and pick the best one we will force dubious cases
      # like test_FLT3_4_hsp_with_intronic_sequence.tab,
      # where the contig contains a continuous run of exonic and
      # intronic sequence and therefore seems dubious.
      #
      push @{$features}, "blast_subject_overlap";
      my $key = $best[0];
      my ($i_first, $i_second) = split /\-/, $key;
      die "HSPs not consecutive" unless $i_second - $i_first == 1;
      my $hsp_left = $hsps->[$i_first];
      my $hsp_right = $hsps->[$i_second];

      $hi_left = new HSPIndexer("-hsp" => $hsp_left);
      $hi_right = new HSPIndexer("-hsp" => $hsp_right);
      # translates between query and hit base numbers and string.
      # for both of these operations we need a position translation
      # table that accounts for gaps in the alignment.
      #
      # These should be used rather than the HSP in downstream
      # code, because the alignments may undergo trimming

      if (1) {
	my $int = $intersect{$key}->[1];
	printf STDERR "debug raw subject intersection: %s\n",
	  run_list $int;
      }

      my $int;
      # subject intersection

      if (query_overlap_rescue($hi_left, $hi_right, $features)) {
	$int = $hi_left->get_intersection("subject", $hi_right);
	# recalculate subject intersection after rescue adjustment
      } else {
	# no query overlap found, so we can use already-calculated
	# subject intersection
	$int = $intersect{$key}->[1];
      }

      # after any query overlap trimming, verify HSPs are
      # ordered not just for the query sequence but for the
      # subject sequence:
      if ($hi_right->start("subject") < $hi_left->start("subject")) {
	# sanity check
	push @{$features}, "subject_overlap_but_HSPs_out_of_order";
	#	    confess sprintf "FATAL ERROR: subject regions in query-ordered HSPs are out of order (query left subject start=%d, query right subject start=%d. possible large ITD?",
	#	      $hi_left->start("subject"),
	#		$hi_right->start("subject");
      } else {
	my $left_q_end = $hi_left->end("query");
	my $right_q_start = $hi_right->start("query");

	my $q_gap = ($right_q_start - $left_q_end) - 1;
	if ($q_gap > $ITD_MAX_QUERY_GAP) {
	  # disqualified: contig not compatible with this isoform?
	  push @{$features}, sprintf "query_overlap_but_long_query_gap=%d", $q_gap;
	} else {
	  # OK, any query gap within tolerances
	  my $rl = run_list $int;
	  die "overlap $rl not a simple range" unless $rl =~ /^\d+\-\d+$/;
	  # may need updates if there are gaps
	  ($subject_overlap_start, $subject_overlap_end) = split /\-/, $rl;
	  # TO DO:
	  # if there are more than 2 HSPs, how to decide which
	  # copy of the duplication to keep?  What if the second
	  # copy has more contiguous alignment to the reference
	  # afterward?  Need an example
	  $subject_overlap_chunk = substr($rs_seq, $subject_overlap_start - 1,
					  ($subject_overlap_end - $subject_overlap_start) + 1);
	  printf STDERR "subject overlap chunk: %s\n", $subject_overlap_chunk;
	  $found_subject_overlap = 1;
	}
      }
    }
  }

  if ($found_subject_overlap) {

    my $q_interstitial_start_base = $hi_left->get_query_base_for_hit_base($subject_overlap_end) + 1;
    # left HSP: get the query string base number immediately after
    # the end of the first copy of the ITD sequence (which is part
    # of the reference sequence).  This marks the start of the
    # interstitial sequence in the pattern, i.e. the beginning of
    # the ITD proper.

    #      die substr($contig, $q_interstitial_start_base - 1);


    my $q_interstitial_end_base = $hi_right->get_query_base_for_hit_base($subject_overlap_end);
    # right HSP: the final base of the overlap range marks the end of
    # the interstitial sequence, i.e. the end of the second copy of
    # the ITD sequence.

    my $q_interstitial = substr($query_seq,
				$q_interstitial_start_base - 1,
				($q_interstitial_end_base - $q_interstitial_start_base) + 1);
    # chunk of the query sequence containing the duplicated sequence
    # and possibly some additional interstitial sequence

    # for upstream and downstream sequences, use the reference sequence
    # rather than the query sequence, as some query sequences may contain
    # intronic sequence.  This code assumes pure RNA junctions without
    # e.g. intron retention.

    my $reference_upstream = substr($rs_seq,
				    0,
				    $subject_overlap_end);
    my $reference_downstream = substr($rs_seq, $subject_overlap_end);

    die unless $rs_seq eq $reference_upstream . $reference_downstream;
    # QC/sanity

    trim_flanking(\$reference_upstream, \$reference_downstream);

    my $fz2_pattern = sprintf '%s%s%s%s%s',
      $reference_upstream, BREAKPOINT_ITD_LEFT,
	# from the reference, not the query
	$q_interstitial,
	  # from the query, 2nd copy of duplicate sequence plus
	  # any miscellaneous interstitial sequence
	  BREAKPOINT_ITD_RIGHT, $reference_downstream;
    # from the reference, not the query

    if ($fz2_pattern =~ s/($subject_overlap_chunk)}($subject_overlap_chunk){$/\}$1\{$2/) {
      # only possible for perfect duplications
      push @{$features}, "rescued_second_repeat_at_right_end";
    }

    printf STDERR "reference: %s\n", $rs_seq;
    printf STDERR "query: %s\n", $query_seq;
    printf STDERR "subject overlap: %s\n", $subject_overlap_chunk;
    printf STDERR "interstitial: %s\n", $q_interstitial;
    printf STDERR "ref upstream: %s\n", $reference_upstream;
    printf STDERR "ref downstream: %s\n", $reference_downstream;
    printf STDERR "fz2 pattern: %s\n", $fz2_pattern;

    $row->{$F_PATTERN_SEQUENCE} = $fz2_pattern;
    $row->{detected_repeat} = $subject_overlap_chunk;
  }

  return $found_subject_overlap;
}

sub yes_no {
  my ($v) = @_;
  return $v ? "yes" : "no";
}

sub complex_indel_check {
  my (%options) = @_;
  my $hsps = $options{"-hsps"} || die;
  my $rs_seq = $options{"-reference-sequence"} || die;
  my $query_seq = $options{"-query-sequence"} || die;
  my $features = $options{"-features"} || die;
  my $row = $options{"-row"} || die;
  my $gbc = $options{"-gbc"} || die;

  my $found_complex = 0;

  if (@{$hsps} == 1) {
    # find start and end of perfect portion of alignment
    my $hi = new HSPIndexer("-hsp" => $hsps->[0]);
    my ($start_q, $start_h) = find_perfect_edge($hi, "start");
    my ($end_q, $end_h) = find_perfect_edge($hi, "end");
    # these are the upstream and downstream edges, rather than
    # the edges of the interstitial region

    if ($start_q < $end_q) {
      # gap region identified
      my $ref_up = substr($rs_seq,
			  0,
			  $start_h);
      my $ref_down = substr($rs_seq, $end_h - 1);

      my $q_interstitial = substr($query_seq,
				$start_q - 1,
				($end_q - $start_q)
			       );

      if (length($q_interstitial) >= $ITD_COMPLEX_MIN_LENGTH) {
	# qualifies
	$found_complex = 1;

	$row->{genea_feature} = $gbc->get_feature_for_base_number($start_h);
	$row->{geneb_feature} = $gbc->get_feature_for_base_number($end_h);

	trim_flanking(\$ref_up, \$ref_down);
	my $fz2_pattern = sprintf '%s%s%s%s%s',
	  $ref_up, BREAKPOINT_ITD_LEFT,
	    # from the reference, not the query
	    $q_interstitial,
	      BREAKPOINT_ITD_RIGHT, $ref_down;
	# from the reference, not the query

	my $i_s_start = $start_h + 1;
        my $i_s_end = $end_h - 1;
	my $s_interstitial = substr($rs_seq,
				  $i_s_start - 1,
				  ($i_s_end - $i_s_start) + 1);

#	my $s_interstitial = substr($rs_seq,
#			    $start_h + 1,

	push @{$features}, sprintf "complex_sub=%s>%s", $s_interstitial, $q_interstitial;
	$row->{$F_PATTERN_SEQUENCE} = $fz2_pattern;
	$row->{$F_PROCESSING_TYPE} = "complex_sub";
      }
    }
  }

  return $found_complex;
}

sub find_perfect_edge {
  my ($hi, $type) = @_;
  my ($q_pos, $h_pos);
  my ($last_q_pos, $last_h_pos);
  my $direction;
  if ($type eq "start") {
    $q_pos = $hi->start("query");
    $h_pos = $hi->start("subject");
    $direction = 1;
  } elsif ($type eq "end") {
    $q_pos = $hi->end("query");
    $h_pos = $hi->end("subject");
    $direction = -1;
  } else {
    die;
  }

  my $consider_n_mismatch = 0;

  my $qb2hb = $hi->query_base_to_hit_base();
  my $qb2s = $hi->query_base_to_sequence();
  my $hb2s = $hi->hit_base_to_sequence();

  my ($q_edge, $h_edge);

  while (1) {
    my $stop;
    $h_pos = $qb2hb->{$q_pos};

    if ($h_pos) {
      # query base is aligned to subject base
      if ($last_h_pos) {
	# check for gap in subject sequence between query bases
	my $diff = abs($last_h_pos - $h_pos);
	if ($diff > 1) {
#	  print STDERR "VERIFY: gap in subject at q $q_pos, diff=$diff\n";
	  # TO DO: verify this is correct when we encounter an example
	  $stop = 1;
	}
      }

      my $q_base = $qb2s->{$q_pos} || die;
      my $h_base = $hb2s->{$h_pos} || die;
      if ($q_base ne $h_base) {
	# query and hit base don't match
#	die "mismatch for q $q_pos h $h_pos: $q_base $h_base";
	if ($consider_n_mismatch) {
	  $stop = 1;
	} else {
	  $stop = 1 unless uc($q_base) eq "N" or uc($h_base) eq "N";
	}
      }
    } else {
      # gap between query and subject
      $stop = 1;
    }

    if ($stop) {
      $q_edge = $last_q_pos;
      $h_edge = $last_h_pos;
      last;
    }

    $last_q_pos = $q_pos;
    $last_h_pos = $h_pos;
    # alignment is perfect here, save

    $q_pos += $direction;
    # next base to check
  }

  return ($q_edge, $h_edge);
}


sub hit2exons {
  # BLAST reads in a fz2 hit file vs. gene models.
  # what exons are hit?  do these match the source pattern?
  my $f_hit = $FLAGS{hit2exons} || die;
  my $wanted_read_id = $FLAGS{read};

  my $MIN_UNMAPPED_QUERY_REGION_TO_REPORT = 10;
  # setting?  maybe enough to blat?

  find_binary("blastn", "-die" => 1);

  my $blast = get_blast();
  $blast->strand("");
  # input reads may be on + or -

  my $sort_by_span_score = $FLAGS{"sort-by-span-score"};;
  $sort_by_span_score = 1 unless defined $sort_by_span_score;
  # depends: for a large file will chew up a lot of RAM

  printf STDERR "configuration:\n";
  printf STDERR "  sort output by span_score?: %s\n", yes_no($sort_by_span_score);

  my $gb_cache = new GenBankCache();

  my $max_features = 5;
  my @h_feat;
  for (my $i=1; $i <= $max_features; $i++) {
    push @h_feat, "f_" . $i;
  }

  my $rpt = new Reporter(
			 "-file" => basename($f_hit) . ".cds.tab",
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   pattern
					   read
					),
				       "note",
				       "summary",
				       "exon_span_score",
				       @h_feat
				      ],
			 "-auto_qc" => 1,
			);


  # scan results file, find reecords matching read pair ID

  open(IN, $f_hit) || die;
  my $h = <IN>;
  chomp $h;
  my @headers = split /\t/, $h;
  # a pattern file is not a standard tab-delimited file so
  # can't use DelimitedFile.pm

  my @rows_out;

  while (<IN>) {
    chomp;
    my @f = split /\t/, $_;
    my $pid;
    if ($f[0] =~ /^read\-pairs/) {
      # final line in file
      last;
    } elsif ($f[0] =~ /pattern (\S+)/) {
      $pid = $1;
    } else {
      die "can't identify record type for " . $f[0];
    }
    die sprintf "record has %d cols, %s", scalar(@f), $_ unless @f == @headers;
    # first line of each set refers to the pattern and uses the
    # same number of columns
    die sprintf "%s is not a pattern line: %s", $f[0], $_ unless $f[0] =~ "^pattern";
    my %r;
    @r{@headers} = @f;
    # pattern definition line

    for (my $i = 0; $i < 2; $i++) {
      my $l_read = <IN> || die;
      chomp $l_read;
      my @f = split /\t/, $l_read;
      die "not a read line" unless $f[0] =~ /^read (\S+)/;

      my $read_id = $1;
      # TO DO: record more details, etc.
      if ($wanted_read_id) {
	# limit to a specific pair
	next unless $read_id eq $wanted_read_id;
      }

      # TO DO:
      # maybe options to filter by fz2 hit type,
      # count of reads matching multiple features, etc.


      my $read_sequence = $f[1];
      # NOTE: this is the aligned read sequence, which may be
      # reverse-complemented from the original if mapped to -
      $read_sequence =~ s/^\s+//;
      $read_sequence =~ s/\s+$//;
      die if $read_sequence =~ /\s/;

#      dump_die(\%r, "new hit", 1);

      # to do: filter to desired read ID?

      #
      # get accession from hit pattern annotation:
      #
      my $acc_a = $r{genea_acc} || die;
      my $acc_b = $r{geneb_acc} || die;
      die "different accessions" unless $acc_a eq $acc_b;
      if ($acc_a =~ /,/) {
	printf STDERR "multiple accessions in acc_a ($acc_a), using first\n";
	my @set = split /,/, $acc_a;
	$acc_a = $set[0];
      }

      #
      # get GenBank info for refseq:
      #
      my $gbc = new GenBankCDS("-file" => $gb_cache->get($acc_a,
							 "-single" => 1));
      my $info = $gbc->get_info();
      my $ref_seq = $info->{mrna_cds_with_utr5} || die;


      # blast (either strand!) of given read
      # report annotations in read base position??

      my $parser = $blast->blast(
				 "-query" => {
					      "read" => $read_sequence,
					     },
				 "-database" => {
						 $acc_a => $ref_seq
						}
				);
      my $result = $parser->next_result();
      # one object per query sequence (only one query seq)

      my @hits;
      while (my $hit = $result->next_hit()) {
	push @hits, $hit;
      }
      die "no BLAST hits" unless @hits == 1;

      my @hsps;
      while (my $hsp = $hits[0]->next_hsp) {
	push @hsps, $hsp;
      }
      @hsps = sort {$a->start("query") <=> $b->start("query")} @hsps;
      die "no HSPs!" unless @hsps;

      my %feature2query;
      my @f_ordered;

      my $unmapped_query_region_counter = 0;

      for (my $hit_idx = 0; $hit_idx < @hsps; $hit_idx++) {
	my $hsp = $hsps[$hit_idx];

	if ($hit_idx > 0) {
	  # check for gap in query sequence between HSPs
	  my $hsp_prev = $hsps[$hit_idx - 1];

	  my $q_gap_start = $hsp_prev->end("query") + 1;
	  my $q_gap_end = $hsp->start("query") - 1;
	  my $q_gap_size = ($q_gap_end - $q_gap_start) + 1;
	  if ($q_gap_size >= $MIN_UNMAPPED_QUERY_REGION_TO_REPORT) {
	    my $feature = sprintf "unknown_%d", ++$unmapped_query_region_counter;
	    push @f_ordered, $feature;
	    for (my $qb = $q_gap_start; $qb <= $q_gap_end; $qb++) {
	      $feature2query{$feature}{$qb} = 1;
	    }
	  }
	}

	my $hi = new HSPIndexer("-hsp" => $hsp);
	my $qstr = $hsp->strand("query");
	my $hstr = $hsp->strand("hit");

	my ($qs, $qe) = $hsp->range("query");
	my ($hs, $he) = $hsp->range("hit");
	printf STDERR "strand q:%s h:%s  range q:%d-%d h:%d-%d\n",
	  $qstr, $hstr,
	      $qs, $qe, $hs, $he;

	die sprintf "query strand %s", $qstr unless $qstr == 1;
	die sprintf "hit strand %s", $hstr unless $hstr == 1;
	# I *think* the read sequence should always be on the same
	# sequence as the gene model, because fz2 reverse-complements
	# as necessary to match the pattern, which is always in the
	# gene model's orientation rather than genome orientation

	for (my $bn = $hs; $bn < $he; $bn++) {
	  my $feature = $gbc->get_feature_for_base_number($bn);
	  push @f_ordered, $feature;
	  my $qb = $hi->get_query_base_for_hit_base($bn);

	  if (my $qb = $hi->get_query_base_for_hit_base($bn)) {
	    $feature2query{$feature}{$qb} = 1;
	  }
	}
      }

      my $features_ordered = unique_ordered_list(\@f_ordered);

      die "more than $max_features features" if @{$features_ordered} > $max_features;


      my %r;
      $r{pattern} = $pid;
#      $r{read} = $read_id;
      $r{read} = excel_blat_hyperlink($read_sequence, $read_id);
      $r{accession} = $acc_a;

      my $summary = join "-", @{$features_ordered};
      $summary =~ s/exon_/ex/g;
      $summary =~ s/unknown_\d+/unknown/g;
      $r{summary} = $summary;
      # summary of events covered by this read

      # does the read show evidence of a large duplication,
      # i.e. where earlier exons are suddenly repeated?:
      my $exons_out_of_order;
      my $last_exon;
      foreach (@{$features_ordered}) {
	if (/exon_(\d+)/) {
	  my $this_exon = $1;
	  $exons_out_of_order = 1 if $last_exon and $this_exon < $last_exon;
	  $last_exon = $this_exon;
	}
      }
      my @notes;
      push @notes, "exons_out_of_order" if $exons_out_of_order;
      $r{note} = join ",", @notes;

      my %exon_lengths;
      for (my $fi = 0; $fi < $max_features; $fi++) {
	my $feature = "";
	my $interval = "";
	my $sequence = "";
	my $blat = "";

	if ($fi < @{$features_ordered}) {
	  $feature = $features_ordered->[$fi];
	  my $span = new Set::IntSpan(keys %{$feature2query{$feature}});
	  my ($qs, $qe) = si2range($span);
	  $interval = run_list $span;

	  $sequence = substr($read_sequence,
			     $qs - 1,
			     ($qe - $qs) + 1);
	  die "??? $read_sequence" unless $sequence =~ /\w/;

	  my $flen = length($sequence);
	  $exon_lengths{$feature} = $flen if $feature =~ /exon/;

	# https://genome.ucsc.edu/cgi-bin/hgBlat?org=Human&db=hg38&userSeq=AGTGCATCCATGAACTCTGGGGTTCTTCTGGTTCGG&type=Blat%27s%20guess

	  $blat = excel_blat_hyperlink($sequence,
				       sprintf "%s (%d)", $feature, $flen
				      );
	}

	my $h = "f_" . ($fi + 1);
	$r{$h} = $blat;
#	printf STDERR "%s: %s: %d-%d, %s\n", $read_id, $feature, $qs, $qe, $query_chunk;
      }

      my $span_score;
      my @span_lengths = sort {$a <=> $b} values %exon_lengths;
      if (@span_lengths == 2) {
	$span_score = $span_lengths[0] / $span_lengths[1];
	# the closer the two span sizes are, the higher the score
      } elsif (@span_lengths > 2) {
#	die sprintf "fix me, more than 2 spans found!: %s, %s", join(",", keys %exon_lengths), join ",", @span_lengths;
	$span_score = $span_lengths[$#span_lengths - 1] / $span_lengths[$#span_lengths];
	# seems odd, look into this case more
      } else {
	$span_score = 0;
      }
      $r{exon_span_score} = $span_score == 0 ? 0 : sprintf "%.2f", $span_score;

      if ($sort_by_span_score) {
	push @rows_out, \%r;
      } else {
	$rpt->end_row(\%r);
      }

      # TO DO:
      # also track portions of query read that don't match refseq,
      # e.g. intronic
    }
  }

  if ($sort_by_span_score) {
    @rows_out = sort {$b->{exon_span_score} <=> $a->{exon_span_score}} @rows_out;
    foreach my $r (@rows_out) {
      $rpt->end_row($r);
    }
  }

  $rpt->finish();
}

sub excel_blat_hyperlink {
  my ($sequence, $label) = @_;
  my %vars;
  $vars{org} = "Human";
  $vars{db} = "hg38";
  $vars{userSeq} = $sequence;
  $vars{type} = "Blat%27s%20guess";
  # hack, could encode ourselves
  my $url = sprintf 'https://genome.ucsc.edu/cgi-bin/hgBlat?%s', join "&", map {join "=", $_, $vars{$_}} sort keys %vars;
  my $blat = sprintf '=HYPERLINK("%s","%s")', $url, $label;
  # converts to native Excel link
  return $blat;
}

sub hsp_debug {
  my ($hi_left, $hi_right, $tag) = @_;
  $tag = "alignment" unless $tag;
  printf STDERR "%s:\n", $tag;
  printf STDERR "   left: q=%d-%d h=%d-%d\n",
    $hi_left->start("query"), $hi_left->end("query"),
      $hi_left->start("subject"), $hi_left->end("subject");
  printf STDERR "  right: q=%d-%d h=%d-%d\n",
    $hi_right->start("query"), $hi_right->end("query"),
      $hi_right->start("subject"), $hi_right->end("subject");
}

sub get_gene {
  my ($row, $end) = @_;
  my $f;
  if ($end eq "A") {
    $f = $F_GENE_A;
  } elsif ($end eq "B") {
    $f = $F_GENE_B;
  } else {
    die "end must be A or B";
  }

  my $gene;
  if (exists $row->{$f}) {
    $gene = $row->{$f};
  } elsif ($gene = $FLAGS{gene}) {
    # command-line override, e.g. UBTF
  } else {
    dump_die($row, "no field $f");
  }
  return $gene;
}

sub merge_pattern_files {
  # merge fz2 pattern files and:
  #  - unique-ify pattern IDs
  #  - harmonize output format by merging to core/shared columns only
  #  - write pattern.tab files?
  # output can then be run through "condense" process

  my $header_policy = $FLAGS{"header-policy"} || die "-header-policy [core|first|clone]";

  my @infiles = @PATTERN_FILES;
  foreach my $f_list (@SRC_FILES) {
    # listfiles
    my $files = read_simple_file($f_list) || die;
    foreach my $f_src (@{$files}) {
      my $f_out = basename($f_src) . ".fz2.tab";
      if (-s $f_out) {
	push @infiles, $f_out;
      } else {
	die "can't find outfile $f_out for infile $f_src";
      }
    }
  }

  my $append_source = $FLAGS{"append-source"};

  my $h_out;
  if ($header_policy eq "core") {
    $h_out = \@H_FZ2_CORE;
    # strip output to core/shared columns only
  } elsif ($header_policy eq "first") {
    # follow the format of the first file specified
    my $df = new DelimitedFile("-file" => $infiles[0],
			       "-headers" => 1,
			      );
    $h_out = $df->headers_raw();
  } elsif ($header_policy eq "clone") {
    my $df = new DelimitedFile(
			       "-file" => ($FLAGS{patterns} || die "-patterns"),
			       "-headers" => 1,
			      );
    $h_out = $df->headers_raw();
  } else {
    die "-header-policy must be core|first|clone";
  }
  die unless $h_out;

  my $f_out = $FLAGS{out} || "merged.tab";

  my %pid;
  # unique pattern IDs
  my %prefix2num;
  # observed numbers for a given pattern prefix

  my @out_main;
  my @out_prepend;
  # output rows

  my %patch = split /,/, $FLAGS{patch} || "";
  my $blank_undef = $FLAGS{"blank-undef"};
  my $prepend_secondary = $FLAGS{"prepend-secondary"};

  my $file_count = 0;
  foreach my $f_in (@infiles) {
    my $df = new DelimitedFile("-file" => $f_in,
			       "-headers" => 1,
			      );
    my $is_secondary = ++$file_count > 1;

    my $f_pattern = basename($f_in) . ".pattern.tab";
    my $rpt_pattern = $df->get_reporter(
				       "-file" => $f_pattern
				      );
    # input files, patched with updated/unique-ified pattern ID

    while (my $row = $df->get_hash()) {
      my $pid_in = $row->{pattern} || die;

      my $pid_out;
      if ($pid{$pid_in}) {
	my ($prefix, $num) = get_prefix($pid_in);
	my @set = sort {$b <=> $a} @{$prefix2num{$prefix}};
	$pid_out = sprintf '%s-%02d', $prefix, $set[0] + 1;
      } else {
	# unique already
	$pid_out = $pid_in;
      }

      my ($prefix, $num) = get_prefix($pid_out);
      push @{$prefix2num{$prefix}}, $num;
      # track observed pattern numbers

      $pid{$pid_out} = 1;
      # unique pattern ID list

      $row->{pattern} = $pid_out;

      $rpt_pattern->end_row($row);
      # copy of source file with pattern ID updated

      if (%patch) {
	# optionally patch in missing column data
	foreach my $h (keys %patch) {
	  $row->{$h} = $patch{$h} unless exists $row->{$h};
	}
      }

      if ($header_policy eq "first" or $blank_undef) {
	# output follows the format of the first file.
	# populate any columns not specified in other files with blanks.
	foreach my $h (@{$h_out}) {
	  $row->{$h} = "" unless exists $row->{$h};
	}
      }

      if ($append_source) {
	my @sources = split /,/, $row->{$F_SOURCE} || "";
	my %src = map {$_, 1} @sources;
	my $src = file_to_source($f_in);
	unless ($src{$src}) {
	  $row->{$F_SOURCE} = join ",", @sources, $src;
	}
      }

      if ($is_secondary and $prepend_secondary) {
	push @out_prepend, $row;
      } else {
	push @out_main, $row;
      }
    }
    $rpt_pattern->finish();
  }

  #
  #  write merged file:
  #
  my $rpt = new Reporter(
			 "-file" => $f_out,
			 "-delimiter" => "\t",
			 "-labels" => $h_out,
			 "-auto_qc" => 1,
			);

  foreach my $r (@out_prepend, @out_main) {
    $rpt->end_row($r);
  }

  $rpt->finish();
  # merged file


}


sub get_prefix {
  my ($pid) = @_;
  my $prefix = $pid;
  $prefix =~ s/\-(\d+)$// || die;
  return ($prefix, int $1);
}

sub cicero_clone_patch {
  # clone fz2 core annotation fields from CICERO-format input
  my ($row) = @_;
  my %map = (
	     "genea_symbol" => [ "geneA" ],
	     "genea_chr" => [ "chrA" ],
	     "genea_pos" => [ "genea_pos_adj", "posA" ],
	     # use adjusted position if field is present,
	     # otherwise use raw CICERO position

	     "geneb_symbol" => [ "geneB" ],
	     "geneb_chr" => [ "chrB" ],
	     "geneb_pos" => [ "geneb_pos_adj", "posB" ],
	    );
  foreach my $f_out (keys %map) {
    unless (exists $row->{$f_out}) {
      my $set_in = $map{$f_out} || die;
      my $v;
      foreach my $f_in (@{$set_in}) {
	if (exists $row->{$f_in}) {
	  $v = $row->{$f_in};
	  last;
	}
      }
      dump_die($row, "no mappable value for $f_out") unless defined $v;
      $row->{$f_out} = $v;
    }
  }

  $row->{pathogenicity_somatic} = "" unless exists $row->{pathogenicity_somatic};
  # doesn't exist for some sources, so not a fatal error

}

sub itd_feature_annotate {
  my ($row, $gbc, $hi_left, $hi_right) = @_;
  my $left_subj_end = $hi_left->end("subject");
  my $right_subj_start = $hi_right->start("subject");
  $row->{genea_feature} = $gbc->get_feature_for_base_number($left_subj_end);
  $row->{geneb_feature} = $gbc->get_feature_for_base_number($right_subj_start);
}

sub file_to_tag {
  my ($infile) = @_;
  my $infile_tag = basename($infile);
  $infile_tag =~ s/\..*//;
  return $infile_tag;
}

sub standardize_genome {
  my ($g) = @_;
  my $result;
  if ($g eq "GRCh37-lite" or $g eq "hg19") {
    $result = 37;
  } elsif ($g eq "GRCh38") {
    $result = 38;
  } else {
    die "can't standardize genome string $g";
  }
  return $result;
}

sub del2fuzzion {
  # create RNA pattern from genomic deletion specification
  # TO DO:
  # - support partial deletion of exons
  # - support "any_result_ok" in failed file output a-la ITD code
  my $f_in = $FLAGS{del2fuzzion} || die "-del2fuzzion";
  my $fai = get_fai();
  my $rf = get_refflat();

  my $f_out = $FLAGS{out} || basename($f_in) . ".fz2.tab";
  my $f_failed = basename($f_in) . ".fz2_failed.tab";

  my $rpt = new Reporter(
			 "-file" => $f_out,
			 "-delimiter" => "\t",
			 "-labels" => \@H_FZ2_CORE,
			);

  my $source = file_to_source($f_in);

  my $df = new DelimitedFile("-file" => $f_in,
			     "-headers" => 1,
			     );

  my $verbose = 0;
  my %pattern_counter;

  my $rpt_failed = $df->get_reporter("-file" => $f_failed,
				     "-extra" => [
						  "accession",
						  "exception"
						 ]
				    );


  while (my $row = $df->get_hash()) {
    my $gene = $row->{gene} || die "no gene";
    my $chr = $row->{chr} || die "no chr";
    my $del_start = $row->{del_start} || die "no del_start";
    my $del_end = $row->{del_end} || die "no del_end";
    my $si_del = new Set::IntSpan sprintf '%d-%d', $del_start, $del_end;

    my $hits = $rf->find_by_gene($gene);
    die "no refflat entries for $gene" unless $hits and @{$hits};
    foreach my $rfr (@{$hits}) {
      my $acc = $rfr->{name} || die;
      die "refflat accession $acc not versioned, use e.g. ncbiRefSeq table insted" unless $acc =~ /\.\d+$/;
      if (compare_chrom($chr, $rfr->{chrom})) {
	my $exons = $rfr->{exons};

	my $rl_exons = join ",", map {join "-", @{$_}{qw(start end)}} @{$exons};
	my @exons = split /,/, $rl_exons;
	printf STDERR "in: %s\n", $rl_exons if $verbose;

	my $si_exons = new Set::IntSpan $rl_exons;

	my $intersect = intersect $si_del $si_exons;
	if ($intersect->size()) {
	  # gene model intersects target deletion
	  my $exons_out = diff $si_exons $si_del;
	  # exonic ranges ex-deletion

	  printf STDERR "intersect: %s\n", $intersect->run_list() if $verbose;
	  my $rl_out = run_list $exons_out;
	  printf STDERR "non-intersect: %s\n", $rl_out if $verbose;

	  my %exons_all = map {$_, 1} @exons;

	  my @ex_out = split /,/, $rl_out;
	  my %exons_out = map {$_, 1} @ex_out;

	  foreach my $ex (keys %exons_out) {
	    die sprintf "output exon %s not canonical, only whole-exon deletions currently supported" unless $exons_all{$ex};
	    # TO DO
	  }

	  my @upstream;
	  my @downstream;

	  my $in_downstream;

	  for my $exon (@exons) {
	    if ($exons_out{$exon}) {
	      if ($in_downstream) {
		push @downstream, $exon;
	      } else {
		push @upstream, $exon;
	      }
	    } else {
	      $in_downstream = 1;
	    }
	  }

	  if ($verbose) {
	    printf STDERR "up: %s\n", join ",", @upstream;
	    printf STDERR "down: %s\n", join ",", @downstream;
	  }

	  my $exception;
	  $exception = "no_upstream_sequence" unless @upstream;
	  $exception = "no_downstream_sequence" unless @downstream;

	  if ($exception) {
	    my %r = %{$row};
	    $r{accession} = $acc;
	    $r{exception} = $exception;
	    $rpt_failed->end_row(\%r);

	    next;
	  }

	  my $strand = $rfr->{strand} || die;
	  my $is_minus;
	  if ($strand eq "+") {
	  } elsif ($strand eq "-") {
	    $is_minus = 1;
	  } else {
	    die "that's unpossible!";
	    # ralph wiggum voice
	  }

	  my $upstream = "";
	  my $breakpoint_a;
	  foreach my $exon (@upstream) {
	    my ($start, $end) = split /\-/, $exon;
	    $breakpoint_a = $end;

	    my $chunk = $fai->get_chunk(
					"-id" => $chr,
					"-start" => $start,
					"-end" => $end
				       );
	    printf STDERR "chunk %s\n", $chunk if $verbose;
	    $upstream .= $chunk;
	  }

	  my $downstream = "";
	  my $breakpoint_b;
	  foreach my $exon (@downstream) {
	    my ($start, $end) = split /\-/, $exon;
	    $breakpoint_b = $start unless $breakpoint_b;

	    my $chunk = $fai->get_chunk(
					"-id" => $chr,
					"-start" => $start,
					"-end" => $end
				       );
	    printf STDERR "chunk %s\n", $chunk if $verbose;
	    $downstream .= $chunk;
	  }

	  if ($is_minus) {
	    # adjust sequence for gene model on - strand
	    my $delim = "_";
	    my $contig = join($delim, $upstream, $downstream);
#	    print STDERR "before: $contig\n";
	    $contig = reverse_complement($contig);
#	    print STDERR "after: $contig\n";
	    ($upstream, $downstream) = split($delim, $contig);

	    swap_variables(\$breakpoint_a, \$breakpoint_b);
	  }

	  trim_flanking(\$upstream, \$downstream);

	  my $pattern = sprintf '%s%s%s%s',
	    $upstream,
	      BREAKPOINT_ITD_LEFT,
		BREAKPOINT_ITD_RIGHT,
		  $downstream;

	  printf STDERR "pattern: %s\n", $pattern if $verbose;

	  my %r;

	  my $counter = ++$pattern_counter{$gene};
	  my $pid = sprintf '%s-%s-%02d', $gene, $gene, $counter;
	  $r{$F_PATTERN_ID} = $pid;
	  $r{$F_PATTERN_SEQUENCE} = $pattern;
	  $r{$F_SOURCE} = $source;
	  $r{$F_SAMPLE} = $row->{sample} || die "no sample field in input";
	  $r{$F_GENOME} = standardize_genome($FLAGS{genome} || die "-genome");
	  $r{genea_symbol} = $r{geneb_symbol} = $gene;
	  $r{genea_acc} = $r{geneb_acc} = $acc;
	  $r{genea_chr} = $r{geneb_chr} = $chr;

	  $r{genea_pos} = $breakpoint_a;
	  $r{geneb_pos} = $breakpoint_b;
	  $r{genea_feature} = get_feature_tag($rf, $rfr, $breakpoint_a);
	  $r{geneb_feature} = get_feature_tag($rf, $rfr, $breakpoint_b);

	  $r{pattern_generation_method} = "genomic_deletion";
	  $r{$F_PATHOGENICITY_SOMATIC} = $row->{$F_PATHOGENICITY_SOMATIC} || "";

	  generate_pair_summary(\%r);

	  $rpt->end_row(\%r);

	}
      } else {
	die "TEST: mapping to different chrom";
      }

    }  # $rfr
  }  # $row

  $rpt->finish();
  $rpt_failed->finish();

}


sub compare_chrom {
  my ($chr_a, $chr_b) = @_;
  foreach ($chr_a, $chr_b) {
    s/^chr//i;
  }
  return $chr_a eq $chr_b;
}

sub file_to_source {
  my ($file) = @_;
  my $source = basename($file);
  $source =~ s/\..*// || die;
  return $source;
}

sub swap_variables {
  my ($ref1, $ref2) = @_;
  printf STDERR "before: %s %s\n", $$ref1, $$ref2;
  my $tmp = $$ref1;
  $$ref1 = $$ref2;
  $$ref2 = $tmp;
  printf STDERR "after: %s %s\n", $$ref1, $$ref2;
}

sub dubious_junction_check {
  # check for reasons to disqualify an ITD junction from having
  # a fusion-style junction equivalent
  #   - possible canonical junction
  #   - very close junction produced by small duplication

  my ($gbc, $rs_seq, $ref_up, $ref_down) = @_;

  my $idx = index($rs_seq, $ref_up);
  die if $idx == -1;
  my $up_bp = $idx + length($ref_up);
  # upstream breakpoint in reference sequence
  my $up_feature = $gbc->get_feature_for_base_number($up_bp);

  $idx = index($rs_seq, $ref_down);
  die if $idx == -1;
  my $down_bp = $idx + 1;
  # downstream breakpoint in reference sequence

#  die join "\n", $rs_seq, $ref_up, $ref_down, $up_bp, $down_bp;

  my $down_feature = $gbc->get_feature_for_base_number($down_bp);

#  for (my $i = $up_bp - 10; $i <= $down_bp + 10; $i++) {
#    printf STDERR "%d: %s\n", $i, $gbc->get_feature_for_base_number($i);
#  }

  my $dubious = 0;

  my $up_exon = $up_feature =~ /_(\d+)$/ ? $1 : die $up_feature;
  my $down_exon = $down_feature =~ /_(\d+)$/ ? $1 : die $down_feature;

  my $bp_distance = abs($down_bp - $up_bp);

  printf STDERR "features: up=%s down=%s, bp distance=%d\n",
    $up_feature, $down_feature, $bp_distance;

  $dubious = 1 if $up_exon == $down_exon;
  # due to microhomology during alignment, etc. the sites
  # may not be perfect, so may appear nearby in the same exon
  $dubious = 2 if $up_exon == $down_exon - 1;
  # classic canonical junction breakpoint

  printf STDERR "hey now: possible canonical junction detected: $up_exon $down_exon\n" if $dubious == 2;
#      die "hey now: earlier junction up=$up_exon down=$down_exon" if $down_exon < $up_exon;

  $dubious = 3 if $bp_distance < $ITD_MIN_BREAKPOINT_DISTANCE_FOR_ALTERNATE_PATTERN;

  return $dubious;
}

sub check_identical_down {
  # count of exactly matching bases from start of interstitial
  # sequence to downstream sequence
  my ($interstitial, $ref_down) = @_;
  $interstitial = uc($interstitial);
  $ref_down = uc($ref_down);
  my $count_identical = 0;
  my $len = length($interstitial);
  for (my $i = 0; $i < $len; $i++) {
    if (substr($interstitial, $i, 1) eq
	substr($ref_down, $i, 1)) {
      $count_identical++;
    } else {
      # stop at first mismatch
      last;
    }
  }
  return $count_identical;
}

sub repeat_expansion_fz2 {
  my $fai = get_fai();
  my $chr = $FLAGS{chr} || die;
  my $start = $FLAGS{start} || die;
  my $end = $FLAGS{end} || die;
  my $repeat = $FLAGS{repeat} || die;
  my $count = $FLAGS{count} || die "-count";
  # number of extra repeats to generate in pattern

  my $flank = length($repeat) * $count;
  printf STDERR "flanking sequence length: %d\n", $flank;
  # this will influence conversion of long-read data to
  # synthetic paired data when running fz2, as fz2 patterns
  # must be longer than the input reads.

  my $f_out = $FLAGS{out} || sprintf "repeat_expansion_%s_%d_%d_%s_%d.fz2.tab",
    $chr, $start, $end, $repeat, $count;

  my $rpt = new Reporter(
			 "-file" => $f_out,
			 "-delimiter" => "\t",
			 "-labels" => [
				       $F_PATTERN_ID,
				       $F_PATTERN_SEQUENCE,
				       "note",
				      ],
			 "-auto_qc" => 1,
			);


  my $chunk = $fai->get_chunk(
			      "-id" => $chr,
			      "-start" => $start,
			      "-end" => $end
			     );
  die "specified repeat $repeat doesn't perfectly match target $chunk"
    unless $chunk =~ /^($repeat)+$/i;

  my $counter = 0;
  my $get_id = sub {
    return sprintf "expand-%s-%02d", $repeat, ++$counter;
  };

  #
  #  pattern 1: upstream
  #  - left side: upstream reference sequence, including reference repeat
  #  - right side: generated repeats
  #
  #  TO DO: maybe put reference repeats with generated ones??

  my $pos_right = $end;
  my $pos_left = $end - ($flank - 1);

  my $up = $fai->get_chunk(
			   "-id" => $chr,
			   "-start" => $pos_left,
			   "-end" => $pos_right
			  );
  my $pattern = uc($up .
    BREAKPOINT_FUSION_LEFT .
      BREAKPOINT_FUSION_RIGHT .
	$repeat x $count);

  my %r;
  $r{$F_PATTERN_ID} = &$get_id();
  $r{$F_PATTERN_SEQUENCE} = $pattern;
  $r{note} = "upstream";
  $rpt->end_row(\%r);

  #
  #  pattern 1: downstream
  #  - left side: generated repeats
  #  - right side: downstream reference sequence, including reference repeat
  #
  $pos_left = $start;
  $pos_right = $pos_left + ($flank - 1);
  my $down = $fai->get_chunk(
			   "-id" => $chr,
			   "-start" => $pos_left,
			   "-end" => $pos_right
			  );
  $pattern = uc(
		$repeat x $count .
		BREAKPOINT_FUSION_LEFT .
		BREAKPOINT_FUSION_RIGHT .
		$down
	       );

  %r = ();
  $r{$F_PATTERN_ID} = &$get_id();
  $r{$F_PATTERN_SEQUENCE} = $pattern;
  $r{note} = "downstream";

  $rpt->end_row(\%r);

  $rpt->finish();
}

sub repeat_expansion_ref {
  my $fai = get_fai();
  my $chr = $FLAGS{chr} || die;
  my $start = $FLAGS{start} || die;
  my $end = $FLAGS{end} || die;
  my $repeat = $FLAGS{repeat} || die;
  my $flank = $FLAGS{flank} || die "-flank";

  my @counts;
  # number of extra repeats to generate in pattern
  if (my $count = $FLAGS{count}) {
    push @counts, $count;
  } elsif (my $counts = $FLAGS{"count-list"}) {
    push @counts, split /,/, $counts;
  }
  die "specify -count COUNT | -counts LIST" unless @counts;

  my $f_out = $FLAGS{out} || sprintf "repeat_expansion_%s_%d_%d_%d_%s_%s.fa",
    $chr, $start, $end, $flank, $repeat, join("_", @counts);

  my $wf = new WorkingFile($f_out);
  my $fh = $wf->output_filehandle();

  my $chunk = $fai->get_chunk(
			      "-id" => $chr,
			      "-start" => $start,
			      "-end" => $end
			     );
  die "specified repeat $repeat doesn't perfectly match target $chunk"
    unless $chunk =~ /^($repeat)+$/i;

  foreach my $count (@counts) {
    my $pos_right = $end;
    my $pos_left = $end - ($flank - 1);
    # upstream reference sequence, including the reference repeat region

    my $up = $fai->get_chunk(
			     "-id" => $chr,
			     "-start" => $pos_left,
			     "-end" => $pos_right
			    );

    $pos_left = $end + 1;
    $pos_right = $pos_left + ($flank - 1);
    # downstream reference sequence, after the axreference repeat region

    my $down = $fai->get_chunk(
			       "-id" => $chr,
			       "-start" => $pos_left,
			       "-end" => $pos_right
			      );

    printf $fh ">%s_%d_%d_%s_%d\n%s%s%s\n",
      $chr, $start, $end, $repeat, $count,
	$up,
	  $repeat x $count,
	    $down;
  }

  $wf->finish();

}

sub repeat_expansion_whole {
  # get reference list
  # if ref matches, adjust, otherwise copy
  my $fai = get_fai();
  my $chr_wanted = $FLAGS{chr} || die;
  my $start = $FLAGS{start} || die;
  my $end = $FLAGS{end} || die;
  my $repeat = $FLAGS{repeat} || die;
  my $count = $FLAGS{count} || die "-count";

  my $f_out = $FLAGS{out} || sprintf "repeat_expansion_%s_%d_%d_%s_%s.fa",
    $chr_wanted, $start, $end, $repeat, $count;

  #
  #  verify repeat sequence in reference:
  #
  my $chunk = $fai->get_chunk(
			      "-id" => $chr_wanted,
			      "-start" => $start,
			      "-end" => $end
			     );
  die "specified repeat $repeat doesn't perfectly match target $chunk"
    unless $chunk =~ /^($repeat)+$/i;

  my $wf = new WorkingFile($f_out);
  my $fh = $wf->output_filehandle();
  my $chrs = $fai->rnm->db_names_ordered();
  my $rnm = new ReferenceNameMapper();
  $rnm->add_name($chr_wanted);

  my $found;
  foreach my $chr (@{$chrs}) {
    if ($rnm->find_name($chr)) {
      #
      # adjust this sequence to reflect repaet expansion
      #
      $found = 1;

      my $seq_up = $fai->get_chunk(
				   "-id" => $chr_wanted,
				   "-start" => 1,
				   "-end" => $end
				  );
      # upstream: start of reference through repeat region

      my $seq_repeat = $repeat x $count;

      my $irow = $fai->get_irow("-id" => $chr_wanted);
      my $length = $irow->{sequence_length} || die;

      my $seq_down = $fai->get_chunk(
				   "-id" => $chr_wanted,
				   "-start" => $end + 1,
				   "-end" => $length
				  );
      printf $fh ">%s\n%s%s%s\n", $chr, lc($seq_up), uc($seq_repeat), lc($seq_down);
    } else {
      #
      # other reference sequence: pass through
      #
      my $seq = $fai->get_sequence("-id" => $chr) || die;
      printf $fh ">%s\n%s\n", $chr, $$seq;
      # TO DO: preserve fasta ID line
    }
  }
  $wf->finish();

  build_fai($f_out);
  # index fasta
}

sub detect_feature_edge {
  # does the given query overlap feature occur at the edge of a feature,
  # e.g. the boundary between one codon and the next?
  my ($hi, $qo_start, $qo_end, $gbc) = @_;

  my $hb_left = $hi->get_hit_base_for_query_base($qo_start);
  my $hb_right = $hi->get_hit_base_for_query_base($qo_end) || die "can't get hit base for query base $qo_end";
  # left and right boundaries of query overlap in RefSeq base space

  my $feature_left = $gbc->get_feature_for_base_number($hb_left);
  my $feature_right = $gbc->get_feature_for_base_number($hb_right);

  my $verbose = 1;

  if ($verbose) {
    printf STDERR "features left: q=%d h=%d f=%s ",
      $qo_start, $hb_left, $feature_left;

    printf STDERR "features right: q=%d h=%d f=%s\n",
      $qo_end, $hb_right, $feature_right;
  }

  my $is_edge = 0;

  if ($feature_left eq $feature_right and $feature_left !~ /utr/i) {
    # span covers a single feature.  Check exon boundaries only, not UTR
    my ($hb_up, $hb_down) = sort {$a <=> $b} ($hb_left, $hb_right);
    # hit bases might be on opposite strand from query
    die "hey now, strand" if $hb_up != $hb_left;
    # test: expected this might happen, find an example to confirm

    my $f_up = $gbc->get_feature_for_base_number($hb_up - 1);
    my $f_down = $gbc->get_feature_for_base_number($hb_down + 1);
    printf STDERR "  neighbor features: %s, %s\n", $f_up, $f_down if $verbose;

    $is_edge = 1 if $f_up ne $feature_left;
    $is_edge = 1 if $f_down ne $feature_left;
  }

  #  die "X $qo_start $feature_left   $qo_end $feature_right";

  return $is_edge;
}

sub microhomology_feature_edge_check {
  my ($hi_left, $hi_right, $qo_start, $qo_end, $gbc) = @_;
  # try various configurations of homomologous bases:
  #   - assign all to A side
  #   - assign all to B side
  #   - if more than one base, assign combinations to A and B
  # does one of these combinations yield a perfect exon edge mapping?

  my $trimmed = 0;

  if (detect_feature_edge($hi_left, $qo_start, $qo_end, $gbc)) {
    # microhomology perfectly matches exon edge on A side,
    # so keep that match and trim overlap from B side
    # /research/rgs01/home/clusterHome/medmonso/work/steve/2020_06_13_cicero_contig_extension/2021_01_04_all/blast/ITD_refactor/test_microhomology_A.tab
    if (0) {
      printf STDERR "before:\n";
      $hi_left->print_summary("left");
      $hi_right->print_summary("right");

      printf STDERR "  trim from R: %d-%d\n", $qo_start, $qo_end;
      $hi_right->trim_query_region($qo_start, $qo_end);

      printf STDERR "trimmed:\n";
      $hi_left->print_summary("left");
      $hi_right->print_summary("right");
    }
    $trimmed = "microhomology_exon_edge_A";
  } elsif (detect_feature_edge($hi_right, $qo_start, $qo_end, $gbc)) {
    # microhomology perfectly matches exon edge on B side.
    # UNTESTED, need example to verify
    if (1) {
      printf STDERR "before:\n";
      $hi_left->print_summary("left");
      $hi_right->print_summary("right");

      printf STDERR "  trim from L: %d-%d\n", $qo_start, $qo_end;
      $hi_left->trim_query_region($qo_start, $qo_end);

      printf STDERR "trimmed:\n";
      $hi_left->print_summary("left");
      $hi_right->print_summary("right");
    }
    die "untested";

    $trimmed = "microhomology_exon_edge_B";
  } elsif ($qo_start != $qo_end) {
    # multi-base microhomology
    my $len = ($qo_end - $qo_start) + 1;
    my @range = ($qo_start .. $qo_end);
    for (my $assign_left = 1; $assign_left < $len; $assign_left++) {
#      printf STDERR "index: $assign_left\n";
      my @left;
      my @right;
      for (my $i = 0; $i < @range; $i++) {
	if ($i < $assign_left) {
	  push @left, $range[$i];
	} else {
	  push @right, $range[$i];
	}
      }
      my $is_l_edge = detect_feature_edge($hi_left,
					  min(@left), max(@left),
					  $gbc);
      # check edge status for L portion
      printf STDERR "$len $assign_left L=%s R=%s L_edge:%d\n",
	join(",", @left), join(",", @right), $is_l_edge;
      if ($is_l_edge) {
	printf STDERR "before:\n";
	$hi_left->print_summary("left");
	$hi_right->print_summary("right");

	printf STDERR "  trim from L: %s\n", join ",", @right;
	printf STDERR "  trim from R: %s\n", join ",", @left;
	$hi_left->trim_query_region($right[0], $right[$#right]);
	$hi_right->trim_query_region($left[0], $left[$#left]);

	printf STDERR "trimmed:\n";
	$hi_left->print_summary("left");
	$hi_right->print_summary("right");

	$trimmed = "microhomology_exon_edge_split";
	last;
	# have to break out of loop, will fail if we continue because
	# trimming has already taken place
      }
    }
  }

  return $trimmed;
}

sub filter_by_contig_identity {
  # bucket rows by input_row_number and screen out results that have
  # significantly less contig alignment than others, as these may
  # reflect results for isoforms that sub-optionally match the contig.
  # Don't want "interstitial" sequence in patterns to actually refer
  # to a different reference isoform.
  my %options = @_;
  my $rows_in = $options{"-rows"} || die "-rows";
  my $rpt_suppressed = $options{"-rpt-suppressed"} || die;
#  my $f_key = "input_row_number";
  # FAIL: invalid if data run through mux.pl
  my $f_key = "input_row_md5";
  # use a newer algorithm?

  my %by_line;
  my %keys;
  my @keys;
  foreach my $r (@{$rows_in}) {
    unless (exists $r->{$f_key}) {
      populate_md5sum($r, $options{"-df"} || die "-df");
    }
    my $key = $r->{$f_key} || confess "no field $f_key";
    push @{$by_line{$key}}, $r;
    push @keys, $key unless $keys{$key};
    $keys{$key} = 1;
  }

  my $f_identity_frac = "__contig_identity_frac";

  my @r_out;

  foreach my $key (@keys) {
    my $set = $by_line{$key};

    foreach my $r (@{$set}) {
      # foreach row in set, compute fraction identity between contig
      # and isoform pair
      my $si = new Set::IntSpan();

      $si = si_merge_range($si, $r, "genea_contig_start", "genea_contig_end");
      $si = si_merge_range($si, $r, "geneb_contig_start", "geneb_contig_end");
      my $overlap_size = size $si;
      my $contig = $r->{contig} || die;
      $r->{$f_identity_frac} = $overlap_size / length($contig);
#      printf STDERR "identity: %f\n", $r->{$f_identity_frac};
    }

    # find best identity fraction
    my %fracs = map {$_, 1} map {$_->{$f_identity_frac}} @{$set};
    my ($best_frac) = sort {$b <=> $a} keys %fracs;
    my $cutoff = $best_frac * $FUSION_PREFER_BEST_CONTIG_MATCHES_BEST_RELATIVE_CUTOFF;

    # filter rows based on cutoff:
    foreach my $r (@{$set}) {
      if ($r->{$f_identity_frac} >= $cutoff) {
	# identity check passes, use for pattern generation
	push @r_out, $r;
      } else {
	# reject row
	$r->{$F_PATTERN_ID} = "";
	# TO DO: also some kind of "note" field?
	$rpt_suppressed->end_row($r);
      }
    }
  }

  return \@r_out;
}

sub si_merge_range {
  my ($si, $r, $f1, $f2) = @_;
  foreach ($f1, $f2) {
    dump_die($r, "where is $_") unless exists $r->{$_};
  }
  my $start = $r->{$f1};
  my $end = $r->{$f2};
  if ($start and $end) {
    if ($start <= $end) {
      $si = union $si join "-", $start, $end;
    } else {
      printf STDERR "WARNING: suspicious start %d, end %d\n", $start, $end;
      # buggy data??
      $si = union $si join "-", $end, $start;
    }
  }
  return $si;
}

sub find_interstitial_coding {
  find_binary("blastn", "-die" => 1);
  my $f_in = $FLAGS{patterns} || die "-patterns";
  my $min_interstitial_length = 20;
  # needs to be long enough to successfully blast
  my $min_hsp_identity = 0.94;
  # identity required for each HSP to be usable
  my $min_coding_identity_hq = 0.95;
  my $min_coding_identity_lq = 0.75;
  my $genome = $FLAGS{genome} || die "-genome";
  # needs to be specified for get_refflat

  my $df = new DelimitedFile(
			     "-file" => $f_in,
			     "-headers" => 1,
			    );
  my $f_out = basename($f_in) . ".coding_interstitial.tab";
  my $rpt_bad = $df->get_reporter("-file" => $f_out,
				  "-extra" => [
					       "ic_type",
					       "ic_fraction",
					       # ic = interstitial coding
					      ]
				 );
  # only bad records included

  my $gbc = new GenBankCache();
  my $rf = get_refflat();
  my $gsm = init_gsm($rf);

  my $tfw = new TemporaryFileWrangler();
  my $fa_query = $tfw->get_tempfile("-append" => ".query.fa");
  my $fa_db = $tfw->get_tempfile("-append" => ".db.fa");
  my $blast = get_blast();

  my $pattern_restrict = $FLAGS{pattern};

  my @pr;
  my %all_accs;
  while (my $row = $df->get_hash()) {
    my $interstitial = get_interstitial($row->{sequence});
    my $ilen = length($interstitial);
    next unless $interstitial and $ilen >= $min_interstitial_length;
    # checkable
    next if $row->{fz2_hint} eq "no_condense";
    # RK patterns, generated interstitial, never modify

    push @pr, $row;

    my $gene_a_list = $row->{genea_symbol} || die;
    my $gene_b_list = $row->{geneb_symbol} || die;
    foreach my $gene_list_raw ($gene_a_list, $gene_b_list) {
      foreach my $gene (split /,/, $gene_list_raw) {
	my $gene_rf = $gsm->find($gene);
	if ($gene_rf) {
	  my $rfs = $rf->find_by_gene($gene_rf);
	  foreach my $rf (@{$rfs}) {
	    my $acc = $rf->{name} || die;
	    next unless $acc =~ /NM_/;
	    # coding only
	    $all_accs{$acc} = 1;
	  }
	} else {
	  print STDERR "WARNING, can't find rfs for: $gene\n";
	}
      }
    }
  }

  $gbc->get([sort keys %all_accs]);
  # ensure all local copies of genbank records are available

  my $c = new Counter(\@pr);
  foreach my $row (@pr) {
    my $gene_a_list = $row->{genea_symbol} || die;
    my $gene_b_list = $row->{geneb_symbol} || die;
    # pattern file format
    next if $pattern_restrict and $row->{pattern} ne $pattern_restrict;

    $c->next($row->{pattern});

    my $interstitial = get_interstitial($row->{sequence});
    my $ilen = length($interstitial);
    if ($interstitial and $ilen >= $min_interstitial_length) {
      open(QTMP, ">" . $fa_query) || die;
      printf QTMP ">interstitial\n%s\n", $interstitial;
      close QTMP;
      # temp blast query file

      my $has_bad_perfect;
      my $has_bad_near_perfect;

      my @all_genes;
      foreach my $gene_list_raw ($gene_a_list, $gene_b_list) {
	foreach my $gene (split /,/, $gene_list_raw) {
	  push @all_genes, $gene;
	}
      }

      my $si_query = new Set::IntSpan();

      foreach my $gene_raw (@all_genes) {
	my $gene_rf = $gsm->find($gene_raw) || die;
	my $rfs = $rf->find_by_gene($gene_rf);
	# ignore if blank, already warned above
	foreach my $rf (@{$rfs}) {
	  my $acc = $rf->{name} || die;
	  next unless $acc =~ /NM_/;
	  # coding only
	  my $result_files = $gbc->get($acc);
	  # should be cached already so OK to query one at a time
	  die unless @{$result_files} == 1;
	  my $gbc = new GenBankCDS("-file" => $result_files->[0]);
	  my $info = $gbc->get_info();
	  my $cds = $info->{mrna_cds} || die "no mrna_cds for $acc in " . $result_files->[0];
	  # get coding sequence for transcript

	  open(QTMP, ">" . $fa_db) || die;
	  printf QTMP ">%s\n%s\n", $acc, $cds;
	  close QTMP;
	  # temp db file

	  if (0) {
	    print STDERR "DEBUG: saving query/db files\n";
	    copy($fa_query, sprintf "debug_query.fa") || die;
	    copy($fa_db, sprintf "debug_db_%s.fa", $acc) || die;
	  }

	  my @hsp = blast2hsps(
			       "-blast" => $blast,
			       "-query" => $fa_query,
			       "-database" => $fa_db
			      );

	  foreach my $hsp (@hsp) {
	    my $ni = $hsp->num_identical();
	    next unless $hsp->frac_identical("query") >= $min_hsp_identity;
#	    printf STDERR "debug q range: %d-%d\n", $hsp->range("query");
	    $si_query = union $si_query join("-", $hsp->range("query"));
	    # add query range covered by this HSP to tracker.
	    # may include more than one gene, e.g. CACUL1-RFWD2-01
	  }
	}
      }

      my $q_coverage = size $si_query;
      my $match_frac = $q_coverage / $ilen;
      # coverage of interstitial sequence vs. either gene

      my $ic_type = "";
      if ($match_frac == 1) {
	$ic_type = "perfect";
      } elsif ($match_frac >= $min_coding_identity_hq) {
	$ic_type = "HQ";
      } elsif ($match_frac >= $min_coding_identity_lq) {
	$ic_type = "LQ";
      }
      if ($ic_type) {
	$row->{ic_type} = $ic_type;
	$row->{ic_fraction} = $match_frac;
	$rpt_bad->end_row($row);
      }
    }
  }
  $rpt_bad->finish();

}

sub get_interstitial {
  my ($seq) = @_;
  my $interstitial = "";
#  if ($seq =~ /[\]\}](\w+)[\[\{]/) {
  if ($seq =~ /[\]\}](\w*)[\[\{]/) {
    $interstitial = $1;
  } else {
    die "no pattern breakpoint found";
  }
  return $interstitial;
}

sub populate_md5sum {
  my ($row, $df) = @_;
  $row->{input_row_md5} = md5_hex(map {$row->{$_}} @{$df->headers_raw});
}

sub patch_features {
  my $f_patterns = $FLAGS{"patterns"} || die "-patterns";
  my $f_out = basename($f_patterns) . ".patch_features.tab";
  my $restrict_source = $FLAGS{"restrict-source"};

  #
  # pass 1: get transcript accessions for cache
  #
  my $df = new DelimitedFile("-file" => $f_patterns,
			     "-headers" => 1,
			    );
  my %all_accs;
  while (my $row = $df->get_hash()) {
    my $pattern = $row->{sequence} || die;
    next if $restrict_source and $row->{source} ne $restrict_source;
    foreach my $side (qw(a b)) {
      my $f_acc = sprintf 'gene%s_acc', $side;
      my $acc_list = $row->{$f_acc} || dump_die($row, "no value for $f_acc");
      my @accs = split /,/, $acc_list;
      foreach my $acc (@accs) {
	$all_accs{$acc} = 1;
      }
    }
  }

  my $gbc = new GenBankCache();

  $gbc->get([sort keys %all_accs]) unless $FLAGS{"no-genbank-cache-check"};

  # pass 2: process
  $df = new DelimitedFile("-file" => $f_patterns,
			  "-headers" => 1,
			 );
  my $rpt = $df->get_reporter(
			      "-file" => $f_out,
  			      "-auto_qc" => 1,
			     );

  my $chunk_size = 80;
  # PROSC-C22orf43-01: requires 77 to not have chunk be duplicated

#  my $chunk_size = 5000;
  # TTN requires entire upstream pattern to anchor uniquely

  while (my $row = $df->get_hash()) {
    my $can_process = 1;

    $can_process = 0 if $restrict_source and $row->{source} ne $restrict_source;

    if ($can_process) {
      my $pattern = $row->{sequence} || die;
      foreach my $side (qw(a b)) {
	my $f_feature = sprintf 'gene%s_feature', $side;
	my $f_acc = sprintf 'gene%s_acc', $side;
	unless ($row->{$f_feature}) {
	  # needs to be populated
	  my $chunk;
	  my $chunk_all;
	  if ($side eq "a") {
	    if ($pattern =~ /(.{$chunk_size})[\]\}]/) {
	      $chunk = $1;
	    }
	    $pattern =~ /^(.+)[\]\}]/ || die;
	    # full side of the pattern, maybe be longer or shorter
	    # than the target chunk size.
	    # shorter example: AKAP8-SUGP2-01
	    $chunk_all = $1;
	  } elsif ($side eq "b") {
	    if ($pattern =~ /[\[\{](.{$chunk_size})/) {
	      $chunk = $1;
	    }
	    $pattern =~ /[\[\{](.*)$/ || die;
	    $chunk_all = $1;
	  } else {
	    die;
	  }

	  my $acc_list = $row->{$f_acc} || die;
	  my @accs = split /,/, $acc_list;
	  my %features;
	  foreach my $acc (@accs) {
	    my $type = $side eq "a" ? "downstream" : "upstream";
	    my $feature = get_feature_for_chunk($gbc, $acc, $chunk, $chunk_all, $type);
	    $features{$feature} = 1;
	  }

	  $row->{$f_feature} = join ",", sort keys %features;
	}
      }
    }

    $rpt->end_row($row);
  }

  $rpt->finish();
}


sub get_feature_for_chunk {
  # TO DO: "fuzzy" or "voting" version that searches annots for entire chunk
  my ($gb_cache, $acc, $chunk_small, $chunk_all, $edge) = @_;
  my $f_gb = $gb_cache->get($acc)->[0] || die;
  my $gbc = new GenBankCDS("-file" => $f_gb);
  my $info = $gbc->get_info();
  my $mrna_full = $info->{mrna_full} || die "no mrna_cds for $acc";

  my $idx;
  my $chunk;
  foreach my $try (grep {$_} $chunk_small, $chunk_all) {
    # small chunk might not be defined (smaller-than-chunk-size
    # sequence at start or end of pattern)
    my $i = index($mrna_full, $try);
    die "can't find $try in $acc $mrna_full" if $i == -1;
    my $i2 = index($mrna_full, $try, $i + 1);
    if ($i2 == -1) {
      # if chunk sequence is unique in mrna sequence
      $chunk = $try;
      $idx = $i;
      last;
    }
  }
  die unless $chunk;

  my $fuzzy_codon = $FLAGS{"fuzzy-codon"};
  die "-fuzzy-codon" unless defined $fuzzy_codon;

  my $result;
  if ($fuzzy_codon) {
    #
    # consider the breakpoint base plus 2 additional nt beyond
    # the breakpoint (upstream on 5' end, downstream on 3'),
    # and use the most frequently-occurring annotation.
    #
    # Helpful for small breakpoint imperfections due to microhomology
    # or split codons, e.g. RK pattern A2ML1-RBM39-01
    #
    my ($base_number_start, $base_number_end);
    if ($edge eq "downstream") {
      # chunk is before breakpoint, so get annotations from last base
      # number
      $base_number_end = $idx + length($chunk);
      $base_number_start = $base_number_end - 2;
    } else {
      # chunk immediately follows breakpoint, so use first base
      $base_number_start = $idx + 1;
      $base_number_end = $base_number_start + 2;
    }

    my %annot;
    for (my $i = $base_number_start; $i <= $base_number_end; $i++) {
      my $feature = $gbc->get_feature_for_base_number($i);
      $annot{$feature}++;
    }

    ($result) = sort {$annot{$b} <=> $annot{$a}} keys %annot;
#    dump_die(\%annot, $result) if scalar keys %annot > 1;
  } else {
    #
    # annotate on the breakpoint base only:
    #
    my $base_number;
    if ($edge eq "downstream") {
      # chunk is before breakpoint, so get annotations from last base
      # number
      $base_number = $idx + length($chunk);
    } else {
      # chunk immediately follows breakpoint, so use first base
      $base_number = $idx + 1;
    }
    printf STDERR "acc %s side %s: base number %d mrna: %s\n", $acc, $edge, $base_number, $mrna_full;
    $result = $gbc->get_feature_for_base_number($base_number);
  }
  return $result;
}

sub detect_needed_headers {
  my ($df, $headers) = @_;
  my %h_existing = map {$_, 1} @{$df->headers_raw};
  my @h_needed = grep {!$h_existing{$_}} @{$headers};
  # only add headers if they are not already present, e.g.
  # when refreshing incomplete annotations
  return \@h_needed;
}

sub patch_gene_pair_summary {
  # patch records missing gene_pair_summary field.
  # requires all (non-merged) summary annotations.
  my $f_patterns = $FLAGS{patterns} || die "-patterns";
  my $f_out = basename($f_patterns) . ".gps.tab";

  my %rescue_acc2symbol;
  $rescue_acc2symbol{"NM_001291398.2"} = "FAM21A";
  $rescue_acc2symbol{"NM_001367397.1"} = "FAM21C";
  $rescue_acc2symbol{"NM_001023567.5"} = "GOLGA8B";
  $rescue_acc2symbol{"NM_001368072.2"} = "GOLGA8A";

  # HACK:
  # hardcode lookups between genes and accessions,
  # because (A) there are only currently a few exceptions and
  # (B) gene symbols in the GenBank records may be newer than
  # the symbols we're processing, e.g. FAM21A --> WASHC2A
  #

  my $df = new DelimitedFile("-file" => $f_patterns,
			     "-headers" => 1,
			     );
  my $rpt = $df->get_reporter(
			      "-file" => $f_out,
  			      "-auto_qc" => 1,
			     );

  my $f_gps = "gene_pair_summary";

  my @pair_elements = qw(
		     symbol
		     acc
		     feature
		   );
  while (my $row = $df->get_hash()) {
    unless ($row->{$f_gps}) {
      # patch needed
      my %pair_side;
      foreach my $side (qw(a b)) {
	my %elements;
	foreach my $element (@pair_elements) {
	  my $f = sprintf 'gene%s_%s', $side, $element;
	  my $v = $row->{$f} || dump_die($row, "no $f");
	  $elements{$element} = $v;
	}
	my $key = join "/", map {$elements{$_}} @pair_elements;
	if ($key =~ /,/) {
	  # merged record containing lists
	  printf STDERR "list detected for %s: %s\n", $row->{pattern}, $key;
	  my %has_list;
	  foreach my $e (@pair_elements) {
	    $has_list{$e} = $elements{$e} =~ /,/ ? 1 : 0;
	  }

	  if ($has_list{symbol} and $has_list{acc} and !$has_list{feature}) {
	    # multiple gene symbols/accessions, but feature is the same
	    my $feature = $elements{feature};
	    # fixed
	    my %acc = map {$_, 1} split /,/, $elements{acc};
	    my %symbol = map {$_, 1} split /,/, $elements{symbol};

	    foreach my $acc (sort keys %acc) {
	      my $sym = $rescue_acc2symbol{$acc} || die "no gene symbol for $acc in " . $row->{pattern};
	      die unless $symbol{$sym};
	      delete $symbol{$sym};
	      $key = join "/", $sym, $acc, $feature;
	      push @{$pair_side{$side}}, $key;
	    }
	    die "leftover data" if %symbol;
	  } else {
	    die "unhandled disambiguation case";
	  }
	} else {
	  # single record
	  push @{$pair_side{$side}}, $key;
	}
      }

      my @combinations;

      foreach my $pair_a (@{$pair_side{"a"}}) {
	foreach my $pair_b (@{$pair_side{"b"}}) {
	  # at this point, lists are expected to have been
	  # split into separate entries
	  my $pair = join "-", $pair_a, $pair_b;
	  if ($pair =~ /,/) {
	    printf STDERR "ERROR: %s: summary %s contains lists!\n", $row->{pattern}, $pair;
	    die "unhandled";
	  }
	  push @combinations, $pair;
	}
      }

      my $summary = join ",", @combinations;
      $row->{$f_gps} = $summary;
    }

    $rpt->end_row($row);
  }
  $rpt->finish();
}

sub patch_gene_pair {
  my $f_patterns = $FLAGS{patterns} || die "-patterns";
  my $f_out = basename($f_patterns) . ".gene_pair.tab";

  my @h_new = ( $F_GENE_PAIR );
  my $df = new DelimitedFile("-file" => $f_patterns,
			     "-headers" => 1,
			     );
  my $h_needed = detect_needed_headers($df, \@h_new);

  my $rpt = $df->get_reporter(
			      "-file" => $f_out,
			      "-extra" => $h_needed,
  			      "-auto_qc" => 1,
			     );

  while (my $row = $df->get_hash()) {
    my $pid = $row->{pattern};
    my $pair = pid2pair($pid);
    $row->{$F_GENE_PAIR} = $pair;
    $rpt->end_row($row);
  }
  $rpt->finish();

}

sub pid2pair {
  my ($pid) = @_;
  my $gene_pair = $pid;
  $gene_pair =~ s/\-\d+$// || die;
  return $gene_pair;
}
