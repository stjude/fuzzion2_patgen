#!/bin/env perl
# frame-checking for fuzzion2 annotated patterns
# MNE 2/2023
#
# TO DO:
#   - maybe blast all frames without any score filtering, save best hit,
#     and call in-frame if best hit for both is the same frame?
#   - different approach?: convert fz2 to CICERO format, and run
#     that frame-checking?
#
# - preserve best blast score even if minimum not met?  (to detect if
#   minimum is too high)
#
# - maybe try a meshed merge of raw HSPs?  If one HSP has a clean match
#   to a region, that should be used over the other HSP
#   - align A HSP starting at left side
#   e.g. CHSY3-ACSS3-01, left side is clean match to A, then long ratty match.


use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;

use Bio::Tools::CodonTable;
use Set::IntSpan;

use FileUtils qw(universal_open);
use MiscUtils qw(dump_die);
use DelimitedFile;
use Reporter;
use GenBankCache;
use TdtConfig;
use SJPreferredIsoform;
use BLASTer;
use HSPTweaker;
use GenBankCDS qw(tolerant_protein_compare);

use constant FRAME_CALL_IN_FRAME => "in_frame";
use constant FRAME_CALL_OUT_OF_FRAME => "out_of_frame";
use constant FRAME_CALL_UNKNOWN => "unknown";

use constant FRAME_EXCEPTION_INTRAGENIC => "intragenic_events_not_supported";
# ITD or other intragenic event
use constant FRAME_EXCEPTION_POSSIBLE_A_FIRST_EXON => "possible_first_coding_exon";
use constant FRAME_EXCEPTION_POSSIBLE_B_LAST_EXON => "possible_last_coding_exon";
# gene B portion of fusion is annotated as in the last exon, so there
# may not be enough downstream sequence for BLAST alignment to work
use constant FRAME_EXCEPTION_NON_CODING_ACCESSION => "noncoding_accession";
# annotated only with non-coding refSeq accession (NR_, NG_)
use constant FRAME_EXCEPTION_POSSIBLE_NON_CODING_UTR_5 => "possible_utr5";
use constant FRAME_EXCEPTION_POSSIBLE_NON_CODING_UTR_3 => "possible_utr3";
# 5' or 3' UTR annotation
use constant FRAME_EXCEPTION_POSSIBLE_NON_CODING_OTHER => "possible_noncoding_feature";
use constant FRAME_EXCEPTION_PROTEIN_TRANSLATION_EXCEPTION => "protein_translation_exception";
use constant FRAME_EXCEPTION_MISSING_ACCESSIONS => "missing_transcript_accessions";
# e.g. NM_002537.3, standard protein translation doesn't work

my $VERBOSE;
my $PREFER_SJPI = 1;

my $TRIM_HSP = 1;
my $TRIM_HSP_MIN_RUN_LENGTH = 10;
#my $TRIM_HSP_MIN_RUN_LENGTH = 20;
# PAX5-DACH1-08: has initial match of 16 AA, but ratty run thereafter
#   => while increasing to 20 fixes this, this case arguably
#      better addressed by newer alignment edge rescue process

my $BLAST_TRACKING_MIN_SCORE = 25;
# minimum blast score to track a blast hit at all.  Set at a level
# to remove common, uninformative/junk alignments.  Hits above
# this level will still be kept and reported, to aid with tuning.

#my $BLAST_TRACKING_MIN_IDENTITIES = 5;
my $BLAST_TRACKING_MIN_IDENTITIES = 8;
# The minimum fz2 "strong" overlap match is only 15 nt / 5 AA.
# Need to test how well BLAST actually works in that case though.
# Also, code looks like it will be mostly used for fz2 patterns
# rather than hits, so a longer minimum overlap seems reasonable.
# Accepting very small matches may have unintended consequences
# as well, e.g. ASXL1-HM13-01 selects a 7-AA match over much longer ones.
# Legit?

#my $BLAST_TRACKING_MIN_IDENTITY_FRACTION = 0.95;
my $BLAST_TRACKING_MIN_IDENTITY_FRACTION = 0.90;
#my $BLAST_TRACKING_MIN_IDENTITY_FRACTION = 0.95;
# AFF1-KMT2A-36: long, strong match but with raggesd start, 88%.
#  => this one is addressed by HSP trimming.
# KMT2A-MLLT10-31: 89%, large insertion vs. reference.
#  => looks correct.  HSP trimming does NOT help here.
#     Counting gaps as mismatches fixes this, can happen due exon
#     differences far from the breakpoint (fz2 patterns considered
#     identical if sequence within 400 nt of breakpoint is the same).
#
# PAX5-DACH1-08: only 89% even counting gaps
#   => increasing trim min. run length to 20 fixes, or alignment edge rescue
# EWSR1-ERG-48: even w/trimming, ~91%

#my $BLAST_MIN_SCORE_WARN_IF_DISCARD = 1000;
my $BLAST_MIN_SCORE_WARN_IF_DISCARD = 350;
# debug: flag cases where high-scoring HSPs are discarded.
# may find more edge cases where HSP needs to be either trimmed,
# or identity fraction lowered.
# EWSR1-ERG-48: false negative w/score 377

my $BLAST_MIN_SCORE_FRAME = 35;
# minimum blast score to consider a hit for frame-checking.
# complications:
#  - what if one side of the pattern is very short?
#  - what if pattern has a small indel and there's a sync problem?
#    => presumably alignment should be OK *near the breakpoint*
#       which is our primary interest
#
# ABL1-BCR-01: score 86
# EPS15-KMT2A-01: score 50
# NOTCH2-AGL-01_lq.tab: score 38 (2nd HSP, 2 mismatches, but looks correct)
# COL12A1-COL5A2-02: score 35 (4th HSP, perfect 24, good site)

my $BLAST_MAX_A_B_GAP_AA = 5;
# max gap between end of blast hit for A and start of hit for B
# to consider adjacent.  Allows for some interstitial sequence.

#my $BLAST_PENALTY_GAP_OPEN = 32767;
#my $BLAST_PENALTY_GAP_EXTEND = 32767;
# helps clean up marginal alignments at edges
my $BLAST_PENALTY_GAP_OPEN = 0;
my $BLAST_PENALTY_GAP_EXTEND = 0;
# however: alignment trimming should deal with this pretty effectively,
# and preventing gaps may itself over-trim alignments just based on
# small indels, which shouldn't be a deal-breaker

my $EXPAND_PROTEIN_INTO_UTR = 1;
# whether to generate synthetic codons from 5' UTR upstream of CDS
# and 3' UTR downstream.
# - enables frame-checking for 5' UTR events esp. in gene B
# - makes it easier to anchor events where breakpoint is e.g. near 3' UTR

my $ALIGN_RESCUE_EDGE_MIN_LENGTH = 40;
my $ALIGN_RESCUE_EDGE_MIN_IDENTITY = 0.95;
# noisy HSP rescue: check the fraction of the alignment at the
# breakpoint edge, if this meets minimum length and identity
# requirements, accept the alighment.  Can rescue cases related to
# transcript differences distant from the +/- 400 nt identity used by fz2

my $MASK_FRAMES = 1;
# when aligning a contig to a reference gene sequence (A or B),
# whether to mask the portion of the contig referring to the other gene.
# Recommended as it can avoid ratty alignment issues for similar genes/regions

my $MULTI_HSP_RESCUE = 1;
my $MULTI_HSP_RESCUE_MIN_LENGTH = 25;
# AFF1-KMT2A-21 = 44
# MLLT10-KMT2A-05 = 26

my %FLAGS;
my @clopts = (
	      "-file=s",
	      "-genome=s",
	      "-sjpi=s",
	      "-verbose" => \$VERBOSE,
	      "-prefer-sjpi=i" => \$PREFER_SJPI,
	      "-blast-tracking-min-score=i" => \$BLAST_TRACKING_MIN_SCORE,
	      "-blast-tracking-min-identity=f" => \$BLAST_TRACKING_MIN_IDENTITY_FRACTION,

	      "-blast-min-score-frame=i" => \$BLAST_MIN_SCORE_FRAME,

	      "-build-genbank-cache",
	      "-trim-hsp=i" => \$TRIM_HSP,

	      "-trim-patterns-to-flanking=i",

	      "-expand-protein-into-utr=i" => \$EXPAND_PROTEIN_INTO_UTR,

	      "-mask-frames=i" => \$MASK_FRAMES,
	      "-multi-hsp-rescue=i" => \$MULTI_HSP_RESCUE,
	      "-multi-hsp-rescue-min-length=i" => \$MULTI_HSP_RESCUE_MIN_LENGTH
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

if ($FLAGS{"trim-patterns-to-flanking"}) {
  trim_patterns_to_flanking();
  exit(0);
}

die "minimum BLAST scores not compatible" if $BLAST_MIN_SCORE_FRAME < $BLAST_TRACKING_MIN_SCORE;

printf STDERR "configuration:\n";
printf STDERR "  filter to SJ-preferred isoform?: %s\n", $PREFER_SJPI ? "y" : "n";
printf STDERR "  generating synthetic codons for UTR sequence: %s\n", $EXPAND_PROTEIN_INTO_UTR ? "y" : "n";

my $f_in = $FLAGS{file} || die "-file";

my $sjpi;
if ($PREFER_SJPI) {
  my $f_sjpi = $FLAGS{sjpi};
  unless ($f_sjpi) {
    my $genome = $FLAGS{genome} || die "-genome";
    my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
    $f_sjpi = $config_genome->{GENE_TRANSCRIPT_MATRIX} || die;
  }
  $sjpi = new SJPreferredIsoform(
				 "-file" => $f_sjpi,
				 "-auto_gsm" => 1
				);
}

my $gbc = new GenBankCache();

if ($FLAGS{"build-genbank-cache"}) {
  build_genbank_cache();
  exit(0);
}


my $df = new DelimitedFile("-file" => $f_in,
			   "-headers" => 1,
			  );

my $f_out = basename($f_in) . ".frame.tab";

my $F_CALL = "frame_call";
my $F_EXCEPTION = "frame_exception";
my @h_new = (
	     $F_CALL,
	     $F_EXCEPTION,
	     qw(
		 frame_nm_A
		 frame_nm_B
		 frame_np_A
		 frame_np_B
		 frame_offset_inframe
		 frame_offset_A
		 frame_offset_B
		 frame_breakpoint_gap
		 frame_blast_span_A
		 frame_blast_span_B
		 frame_blast_scores_A
		 frame_blast_scores_B
		 frame_align_coverage
		 frame_notes
	      ));

my $rpt = $df->get_reporter(
			    "-file" => $f_out,
			    "-extra" => \@h_new,
			    "-auto_qc" => 1,
			   );

while (my $row = $df->get_hash()) {
  if (is_itd($row)) {
    # not supported
    foreach (@h_new) {
      $row->{$_} = "";
    }
    $row->{$F_CALL} = FRAME_CALL_UNKNOWN;
    $row->{$F_EXCEPTION} = FRAME_EXCEPTION_INTRAGENIC;
    $rpt->end_row($row);
    next;
  }

#  printf STDERR "processing: %s\n", $row->{pattern};
  my $pattern = $row->{sequence} || die;
  $pattern =~ /^(.*)[\]\}](.*)[\[\{](.*)$/ || dump_die($row, "can't find pattern start/end sites");
#  $pattern =~ /^(.*)[\[\{].*[\]\}](.*)$/ || die;
  my ($nt_a, $nt_interstitial, $nt_b) = ($1, $2, $3);

  my $nt_contiguous = join "", $nt_a, $nt_interstitial, $nt_b;
  # complete pattern sequence without breakpoint markers

  my $pi_end_a = int(length($nt_a) / 3);
  my $pi_start_b = int((length($nt_a) + length($nt_interstitial)) / 3);
  # protein indices for end of gene A / start of gene B
  # NOTE: this should have a bit of buffer due to frame offset (0, 1, 2)

  my $frames = get_protein_frames($nt_contiguous, $pi_end_a, $pi_start_b);
  # need to check alignment w/proteins to confirm ORF works for both
  # HOWEVER: what if pattern construction has bugs?
  # ironically fz2 alignment could show correct aligment, don't want
  # to toss patterns w/o review

  # TO DO:
  # select best PAIR if pairing info available!

  my $nm_a = get_transcript($row, "a");
  my $nm_b = get_transcript($row, "b");

  if (not($nm_a and $nm_b)) {
    foreach (@h_new) {
      $row->{$_} = "";
    }
    $row->{$F_CALL} = FRAME_CALL_UNKNOWN;
    $row->{$F_EXCEPTION} = FRAME_EXCEPTION_MISSING_ACCESSIONS;
    $rpt->end_row($row);
    next;
  } elsif ($nm_a =~ /^NR_/ or $nm_b =~ /^NR_/) {
    foreach (@h_new) {
      $row->{$_} = "";
    }
    $row->{$F_CALL} = FRAME_CALL_UNKNOWN;
    $row->{$F_EXCEPTION} = FRAME_EXCEPTION_NON_CODING_ACCESSION;
    $row->{"frame_notes"} = "not_coding_accession";
    $rpt->end_row($row);
    next;
  }

  my @notes;

  my ($np_a, $protein_a) = get_protein($nm_a);
  my ($np_b, $protein_b) = get_protein($nm_b);
  dump_die($row, "no protein found for $nm_a") unless $np_a;
  dump_die($row, "no protein found for $nm_b") unless $np_b;

  if (!$protein_a or!$protein_b) {
    my $acc_exception;
    if (!$protein_a) {
      $acc_exception = $nm_a;
    } elsif (!$protein_b) {
      $acc_exception = $nm_b;
    } else {
      die;
    }

    foreach (@h_new) {
      $row->{$_} = "";
    }
    $row->{$F_CALL} = FRAME_CALL_UNKNOWN;
    $row->{$F_EXCEPTION} = FRAME_EXCEPTION_PROTEIN_TRANSLATION_EXCEPTION;
    $row->{"frame_notes"} = sprintf "protein_translation_exception=%s", $acc_exception;
    $rpt->end_row($row);
    next;
  }

  unless (pair_sanity_check($row, $nm_a, $nm_b)) {
    # if source is a fz2 pattern file, verify that the selected transcripts
    # are also represented in the pairing list.  Otherwise is it possible
    # that we might choose two transcripts that weren't actually paired
    # when building the pattern?
    #
    # HOWEVER: pair details are not guaranteed in source data!
    # e.g. IKZF1-ZEB2-03 contains a merged accession from RK data
    # which is preferred, but not represented in the pair data.
    #
    # tough call: if we have pair data, maybe we should use it.
    # however as in this example, additional annotations may contain
    # a more-preferred isoform.
    push @notes, sprintf "chosen_transcripts_not_reported_in_pairs=%s/%s", $nm_a, $nm_b;
    # just make a warning
  }

  my %passed_pair;
  my %lq_pair;
  my %passed_side;

  my @blast_scores_a;
  my @blast_scores_b;
  my @blast_span_a;
  my @blast_span_b;

  #
  #  gather info for blast hits meeting the basic minimum.  These are
  #  not necessarily good enough for a frame-checking call, and may be
  #  further pruned later, but preserving the raw info is useful for
  #  debugging:
  #
  foreach my $fkey (sort keys %{$frames}) {
    my ($frame_for_a_align, $frame_for_b_align);
    my $a_suffix = "";
    my $b_suffix = "";

    if ($MASK_FRAMES) {
      $frame_for_a_align = $frames->{$fkey}{protein_mask_b} || die;
      # when aligning gene A, mask portion of contig containing B
      $a_suffix = "_mask_B";
      $frame_for_b_align = $frames->{$fkey}{protein_mask_a} || die;
      # when aligning gene B, mask portion of contig containing A
      $b_suffix = "_mask_A";
    } else {
      # use the raw frame for both alignments; can lead to problems
      $frame_for_a_align = $frame_for_b_align = $frames->{$fkey}->{protein};
    }
    # TO DO: masking for A and B portions?
    # or: just compare indices of hits to make sure matches are to A and B?
    my $hit_a = blast_frame("A", $fkey . $a_suffix, $frame_for_a_align,
			    "A_" . $np_a, $protein_a, \@notes);
    my $score_a = $hit_a ? $hit_a->score() : "";
    $passed_side{"a"}{$fkey} = $hit_a if $hit_a and
      $score_a >= $BLAST_MIN_SCORE_FRAME;
    my $hit_b = blast_frame("B", $fkey . $b_suffix, $frame_for_b_align,
			    "B_" . $np_b, $protein_b, \@notes);
    my $score_b = $hit_b ? $hit_b->score() : "";
    $passed_side{"b"}{$fkey} = $hit_b if $hit_b and
      $score_b >= $BLAST_MIN_SCORE_FRAME;

    push @blast_scores_a, $score_a;
    push @blast_scores_b, $score_b;
    push @blast_span_a, $hit_a ? join("-", $hit_a->range("query")) : "";
    push @blast_span_b, $hit_b ? join("-", $hit_b->range("query")) : "";

    printf STDERR "check %s: A=%s B=%s\n", $fkey, $score_a, $score_b;

    if ($hit_a and $hit_b) {
      printf STDERR "hits: A=%d,%d-%d B=%d,%d-%d\n",
	$hit_a->score(), $hit_a->range("query"),
	  $hit_b->score(), $hit_b->range("query");

      if ($score_a >= $BLAST_MIN_SCORE_FRAME and
	  $score_b >= $BLAST_MIN_SCORE_FRAME) {
	# this frame matches both genes with strong-enough hits
	$passed_pair{$fkey} = 1;
      } else {
	# low-quality pair: possible thresholding issue?
	$lq_pair{$fkey} = 1;
      }
      # TO DO: additional requirements/QC
    }
  }

  prune_hits(\%passed_side);
  # prune hits to only the closest pairing

  # what about lq_pair?

  my $frame_call;
  my $in_frame_offset = "";
  my $frame_exception = "";
  my $frame_breakpoint_gap = "";
  my $frame_align_coverage = "";

  if ($passed_side{"a"} and $passed_side{"b"}) {
    # high-quality matches to both sides found

    my @f_a = keys %{$passed_side{a}};
    my @f_b = keys %{$passed_side{b}};
    die unless @f_a and @f_a == 1;
    die unless @f_b and @f_b == 1;
    my ($f_a) = @f_a;
    my ($f_b) = @f_b;

    if ($f_a eq $f_b) {
      # in frame
      $frame_call = FRAME_CALL_IN_FRAME;
      $in_frame_offset = $frames->{$f_a}{offset};
    } else {
      # different frames
      $frame_call = FRAME_CALL_OUT_OF_FRAME;
    }

    my $hit_a = $passed_side{a}{$f_a} || die;
    my $hit_b = $passed_side{b}{$f_b} || die;

    #
    # calculate fraction of query covered by alignments:
    #
    my $span_a = join "-", $hit_a->range("query");
    my $span_b = join "-", $hit_b->range("query");
    my $span_all = new Set::IntSpan($span_a, $span_b);
    my $span_size = $span_all->size();

    my $frame_length = length $frames->{frame_offset_0}{protein};
    $frame_align_coverage = sprintf "%.2f", $span_size / $frame_length;

    #
    # calculate gap between gene A/B alignments:
    #
    my $a_end = $hit_a->end("query");
    my $b_start = $hit_b->start("query");
    my $gap = ($b_start - $a_end) - 1;
    $frame_breakpoint_gap = ($b_start - $a_end) - 1;
    my $qc_abut;
    if ($gap <= 0) {
      $qc_abut = "perfect";
    } elsif ($gap <= $BLAST_MAX_A_B_GAP_AA) {
      $qc_abut = "ok";
    } else {
      $qc_abut = "distant";
    }
    push @notes, sprintf "qc_blast_proximity=%s", $qc_abut || die;
  } else {
    # can't establish a good match for both sides
    $frame_call = FRAME_CALL_UNKNOWN;

    foreach my $side (qw(a b)) {
      unless ($passed_side{$side}) {
	# check sides without good evidence for annotations
	# showing exclusively non-coding features like UTR
	my $f_annot = sprintf 'gene%s_feature', $side;
	if (exists $row->{$f_annot}) {
	  # fuzzion2 pattern annotation
	  my @annot = split /,/, $row->{$f_annot};
	  my $usable;
	  my %other_feature;
	  foreach my $feature (@annot) {
	    if ($feature =~ /^exon_/) {
	      $usable = 1;
	    } else {
	      $feature =~ s/_\d+$//;
	      # strip feature #
	      $other_feature{$feature} = 1;
	    }
	  }
	  unless ($usable) {
	    printf STDERR "non-exonic features: %s\n", join ",", sort keys %other_feature;
	    if ($other_feature{"5utr"}) {
	      $frame_exception = FRAME_EXCEPTION_POSSIBLE_NON_CODING_UTR_5;
	    } elsif ($other_feature{"3utr"}) {
	      $frame_exception = FRAME_EXCEPTION_POSSIBLE_NON_CODING_UTR_3;
	    } else {
	      $frame_exception = FRAME_EXCEPTION_POSSIBLE_NON_CODING_OTHER;
#	      dump_die(\%other_feature, "need exception code for non-coding");
	    }
	  }
	}
      }
    }

    if (not($frame_exception) and not($passed_side{"a"})) {
      $frame_exception = FRAME_EXCEPTION_POSSIBLE_A_FIRST_EXON
	if cds_start_check($row, $nm_a);
    }

    if (not($frame_exception) and not($passed_side{"b"})) {
      # check if the B side begins in the last exon.  If it does, this
      # could mean there's not enough downstream sequence for blast to
      # anchor.  possible example: EWSR1-ERG-33 TO DO: similar check for
      # gene A if breakpoint is in first exon
      #
      # TO DO: revise logic for 3' UTR?
      $frame_exception = FRAME_EXCEPTION_POSSIBLE_B_LAST_EXON
	if cds_end_check($row, $nm_b);
    }

    push @notes, "LQ_pair_only" if %lq_pair;
  }
  die unless defined $frame_call;

  my @passed_a = sort keys %{$passed_side{a}};
  my $foa = join ",", map {$frames->{$_}{offset}} @passed_a;

  my @passed_b = sort keys %{$passed_side{b}};
  my $fob = join ",", map {$frames->{$_}{offset}} @passed_b;


  push @notes, "A_matches_multiple_frames" if @passed_a > 1;
  push @notes, "B_matches_multiple_frames" if @passed_b > 1;
  # TO DO: warn if in-frame frame doesn't have best scores

  $row->{$F_CALL} = $frame_call;
  $frame_exception = "unknown" if not($frame_exception) and
    $frame_call eq FRAME_CALL_UNKNOWN;
  $row->{$F_EXCEPTION} = $frame_exception;
  $row->{frame_offset_inframe} = $in_frame_offset;
  $row->{frame_breakpoint_gap} = $frame_breakpoint_gap;
  $row->{frame_align_coverage} = $frame_align_coverage;
  $row->{frame_offset_A} = $foa;
  $row->{frame_offset_B} = $fob;
  $row->{frame_nm_A} = $nm_a;
  $row->{frame_nm_B} = $nm_b;
  $row->{frame_np_A} = $np_a;
  $row->{frame_np_B} = $np_b;
  $row->{frame_blast_scores_A} = join ",", @blast_scores_a;
  $row->{frame_blast_scores_B} = join ",", @blast_scores_b;
  $row->{frame_blast_span_A} = join ",", @blast_span_a;
  $row->{frame_blast_span_B} = join ",", @blast_span_b;
  $row->{frame_notes} = join ",", @notes;

  $rpt->end_row($row);
}
$rpt->finish();

sub blast_frame {
  my ($side, $frame_id, $frame, $np_id, $protein, $notes) = @_;
  my $blast = get_blast();

  my $hsp = run_blast(
		      "-blast" => $blast,
		      "-query" => {
				   $frame_id => $frame
				  },
		      "-database" => { $np_id => $protein },
		      "-notes" => $notes,
		      "-key" => join("_", $frame_id, $np_id),
		      "-side" => $side
		     );
  # hack: with this invocation, only the best HSP is returned
  return $hsp;
}

sub get_protein {
  # get local copy of GenBank record and parse out protein annotation
  my ($nm) = @_;
  my $f_gb = $gbc->get($nm, "-single" => 1);
  my $fh = universal_open($f_gb);
  my $stream = Bio::SeqIO->new("-fh" => $fh,
			       -format => 'GenBank');
  my $record = $stream->next_seq();
  # Bio::SeqIO::genbank
  # cache files contain only one record per
  return genbank2protein($record);
}

sub get_exon_info {
  # get local copy of GenBank record and parse out exon info
  my ($nm) = @_;
  my $f_gb = $gbc->get($nm, "-single" => 1);
  my $fh = universal_open($f_gb);
  my $stream = Bio::SeqIO->new("-fh" => $fh,
			       -format => 'GenBank');
  my $record = $stream->next_seq();
  # Bio::SeqIO::genbank
  # cache files contain only one record each

  # parse:
  my $exon_count = 0;
  my ($exon_start, $exon_end);
  my ($cds_start, $cds_end);
  my ($cds_start_exon, $cds_end_exon);
  foreach my $feature ($record->get_SeqFeatures()) {
    my $ptag = $feature->primary_tag();
    if ($ptag eq "exon") {
      $exon_count++;
      $exon_start = $feature->start;
      $exon_end = $feature->end;
      # TO DO: also record/return exon containing CDS end
#      printf STDERR "exon:%d-%d CDSstart: %d\n", $exon_start, $exon_end, $cds_start || "-1";
      if ($cds_end and not($cds_end_exon) and
	  $cds_end >= $exon_start and $cds_end <= $exon_end) {
	# CDS ends in this exon
	$cds_end_exon = $exon_count;
      }
    } elsif ($ptag eq "CDS") {
      $cds_start = $feature->start();
      $cds_end = $feature->end();
      if ($cds_start >= $exon_start) {
	$cds_start_exon = $exon_count;
      } else {
	die "unexpected CDS location";
      }
    }
  }

  die unless wantarray;
  return ($exon_count, $cds_start_exon, $cds_end_exon);
}

sub get_transcript {
  my ($row, $which) = @_;
  my $f = sprintf "gene%s_acc", $which;
  my $result;
  if (my $acc_list = $row->{$f}) {
    my @set = split /,/, $acc_list;
    if (@set > 1 and $PREFER_SJPI) {
      my @pref;
      foreach my $nm (@set) {
	push @pref, $nm if $sjpi->is_preferred_nm($nm);
      }
      if (@pref and @pref == 1) {
	@set = @pref;
      } else {
	# preferred isoform not available, use highest-ranking instead
	my %rank2nm;
	foreach my $nm (@set) {
	  my $rank = $sjpi->get_nm_rank($nm);
	  unless ($rank) {
	    if ($nm =~ /^NM_/) {
	      # coding refseq
	      $rank = 1000;
	    } else {
	      $rank = 2000;
	      # NR_ / NG_ etc.
	    }
	  }
	  push @{$rank2nm{$rank}}, $nm;
	}
	my @ranks = sort {$a <=> $b} keys %rank2nm;
	my $best_rank = $ranks[0];
	@set = $rank2nm{$best_rank}->[0];
      }
    }
    $result = $set[0];
  }

  return $result;
}


sub get_protein_frames {
  # translate a nucleotide contig into 3 possible protein coding frames
  my ($contig, $pi_end_a, $pi_start_b) = @_;

  $pi_end_a--;
  $pi_start_b++;
  # add a 1-AA buffer to avoid over-masking and frame offset shifts (0,1,2)

  my $ct = new Bio::Tools::CodonTable();

  my @nt = split //, $contig;

  my %frames;
  foreach my $offset (0 .. 2) {
    my $string = substr($contig, $offset);
    my @codons = $string =~ /(...)/g;
    my $id = sprintf "frame_offset_%d", $offset;
    my $sequence = join "", map {$ct->translate($_)} @codons;

    my $sequence_mask_a = $sequence;
    for (my $i = 0; $i < $pi_end_a; $i++) {
      substr($sequence_mask_a, $i, 1) = "X";
    }
    my $sequence_mask_b = $sequence;
    my $len = length($sequence);
    for (my $i = $pi_start_b; $i < $len; $i++) {
      substr($sequence_mask_b, $i, 1) = "X";
    }


    if ($VERBOSE) {
      printf STDERR "frame %s: %s\n", $id, $sequence;
      printf STDERR "  mask_A: %s\n", $sequence_mask_a;
      printf STDERR "  mask_B: %s\n", $sequence_mask_b;
    }
    $frames{$id} = {
		    "protein" => $sequence,
		    "nucleotide" => $string,
		    "offset" => $offset,
		    "protein_mask_a" => $sequence_mask_a,
		    "protein_mask_b" => $sequence_mask_b,
		   };
  }
  return \%frames;
}


sub genbank2protein {
  # Bio::SeqIO::genbank
  #
  # TO DO: return NP_ identifier too
  my ($record) = @_;
  my $cds_count = 0;

  my ($protein_id, $protein, $protein_expanded);

  my $nucleotide_seq = $record->seq() || die;

  my ($cds_start, $cds_end);

  foreach my $feature ($record->get_SeqFeatures()) {
    my $ptag = $feature->primary_tag();
    if ($ptag eq "CDS") {
      $cds_count++;
      $cds_start = $feature->start();
      $cds_end = $feature->end();
      foreach my $tag ($feature->get_all_tags()) {
	if ($tag eq "translation") {
	  my @values = $feature->get_tag_values($tag);
	  die unless @values == 1;
	  $protein = $values[0];
	} elsif ($tag eq "protein_id") {
	  my @values = $feature->get_tag_values($tag);
	  die unless @values == 1;
	  $protein_id = $values[0];
	}
      }
    } elsif ($ptag eq "gene") {
      my $start = $feature->start();
      my $end = $feature->end();
      die "gene doesn't cover entire nt sequence" unless $start == 1 and $end == length($nucleotide_seq);
    }
  }
  die "multiple CDS" if $cds_count > 1;
  # hopefully unpossible

  die unless $cds_start;

  #
  #  manually translate CDS and compare w/canonical protein:
  #
  my $cds_string = substr($nucleotide_seq,
			  $cds_start - 1,
			  ($cds_end - $cds_start) + 1);
  my $cds_aa_recalc = nt2protein($cds_string, 1);
  if (not($cds_aa_recalc)) {
    # likely exception from standard protein translation
    return ($protein_id, undef);
  } elsif ($cds_aa_recalc eq $protein) {
    # perfect, typical case
  } elsif (tolerant_protein_compare($protein, $cds_aa_recalc)) {
    # within small tolerances for expected exceptions:
    # e.g. NM_020831:
    # /transl_except=(pos:326..328,aa:Met)
    #
    # a bonus here is that this new translation using standard genetic
    # code will also better match the sequences used during frame-checking,
    # which also use only standard genetic code
  } else {
    die join "\n", "recomputed CDS doesn't match:", $cds_aa_recalc, $protein;
  }

  my $verbose = 0;

  printf STDERR "cds start: %d\n", $cds_start if $verbose;
  my $avail_utr5_bases = $cds_start - 1;
  if (0) {
    print STDERR "DEBUG\n";
    $avail_utr5_bases = 3;
  }

  printf STDERR "avail utr5 bases: %d\n", $avail_utr5_bases if $verbose;
  my $max_synthetic_codons = int($avail_utr5_bases / 3);
  printf STDERR "max synthetic codons: %d\n", $max_synthetic_codons if $verbose;
  my $fake_cds_start = $cds_start - ($max_synthetic_codons * 3);
  printf STDERR "fake CDS start: %d\n", $fake_cds_start if $verbose;
  my $extended_nt = substr($nucleotide_seq, $fake_cds_start - 1);
  printf STDERR "  nucleotide: %s\n", $extended_nt if $verbose;

  my $ct = new Bio::Tools::CodonTable();
  my @codons = $extended_nt =~ /(...)/g;
  my $extended_aa = join "", map {$ct->translate($_)} @codons;
  printf STDERR "  AA: %s\n", $extended_aa if $verbose;

 confess sprintf "can't find canonical protein in extended protein for %s:\n%s\n%s\n", $record->accession, $protein, $extended_aa unless index($extended_aa, $cds_aa_recalc) >= 0;
  # use the recomputed sequence for lookup as it uses the standard
  # genetic code (RefSeq AA may rarely contain some exceptions)

  return ($protein_id, $EXPAND_PROTEIN_INTO_UTR ? $extended_aa : $protein);
}

sub get_blast {
  # adapted from fusion_contig_extension.pl
  my $blast = new BLASTer();
  $blast->blast_binary("blastp");

  $blast->gapopen($BLAST_PENALTY_GAP_OPEN) if $BLAST_PENALTY_GAP_OPEN;
  $blast->gapextend($BLAST_PENALTY_GAP_EXTEND) if $BLAST_PENALTY_GAP_EXTEND;

  $blast->strand("plus");
  # contig and transcript sequence have been normalized to +

#  $blast->perc_identity($BLAST_MIN_PERCENT_IDENTITY);
#  $blast->word_size($BLAST_MIN_WORD_SIZE);
#  $blast->gapopen($BLAST_GAP_OPEN) if $BLAST_GAP_OPEN;
#  $blast->gapextend($BLAST_GAP_EXTEND) if $BLAST_GAP_EXTEND;
# TO DO: revisit

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
  #
  # HOWEVER: not sure if this is necessary for proteins

  $blast->blast_2_sequences(1);
  # just one sequence vs. another, don't need to index database

  $blast->output_format("xml");
  # provides visibility for gapped alignments

  if ($FLAGS{"save-tempfiles"}) {
    printf STDERR "saving tempfiles\n";
    $blast->tfw->auto_unlink(0);
  }
  return $blast;
}


sub run_blast {
  my %options = @_;
  my $blast = $options{"-blast"} || die "-blat";
  my $multi_ok = $options{"-multi-ok"};
  my $notes = $options{"-notes"} || die "-notes";
  my $key = $options{"-key"} || die;
  my $side = $options{"-side"} || die;

  printf STDERR "start BLAST\n", if $VERBOSE;

  my $parser = $blast->blast(%options);
  my $result = $parser->next_result;
  # one object per query sequence (only one query seq)
  my @hsp;
  if ($result) {
    while (my $hit = $result->next_hit()) {
      # hits from this query to a database sequence
      # (Bio::Search::Hit::HitI)

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

	my $ok = 1;
	my $ht;
	if ($TRIM_HSP) {
	  $ht = new HSPTweaker();
	  $ht->min_run_length($TRIM_HSP_MIN_RUN_LENGTH);
	  $ht->count_gaps_as_matches(1);
	  # fz2 considers patterns identical within +/- 400 nt identical,
	  # some sequences are longer though, and duplicate-checking
	  # may bucket patterns together that have exon differences
	  # further away from the breakpoint.

	  $ht->parse($hsp);
	  if ($ht->trimmed()) {
	    $hsp = $ht;
	    # hack: only implements subset of features
	  } elsif (my $shrink_frac = $ht->trim_rejected()) {
	    # trimming rejected for excessive shrinkage,
	    # add debug message for tuning
	    push @{$notes}, sprintf "trim_rejected=%s/%.2f", $key, $shrink_frac unless @hsp;
	    # don't report if we already have a passed HSP, e.g. TCF7-SPI1-16
	  }
	}


	$ok = 0 unless $hsp->score() >= $BLAST_TRACKING_MIN_SCORE;
	$ok = 0 unless $hsp->num_identical() >= $BLAST_TRACKING_MIN_IDENTITIES;
	$ok = 0 unless $hsp->frac_identical("query") >= $BLAST_TRACKING_MIN_IDENTITY_FRACTION;

	printf STDERR "HSP debug: score=%d num_identical=%d frac_ident_query=%f ok=%d passed_count:%d\n", $hsp->score, $hsp->num_identical, $hsp->frac_identical("query"), $ok, scalar @hsp if $VERBOSE;

	if (not($ok) and $hsp->score() >= $BLAST_MIN_SCORE_WARN_IF_DISCARD) {
	  my $rescued;
	  if ($TRIM_HSP) {
	    # rescue if alignment on breakpoint edge is very clean
	    my $chunk_identity = $ht->check_edge_identity(
							  "-side" => $side,
							  "-length" => $ALIGN_RESCUE_EDGE_MIN_LENGTH,
							 );
	    if ($chunk_identity >= $ALIGN_RESCUE_EDGE_MIN_IDENTITY) {
	      $rescued = 1;
	    }
	  }

	  if ($rescued) {
	    $ok = 1;
	    push @{$notes}, sprintf "rescued_rejected_high_quality_HSP=%s/%d/%.4f", $key, $hsp->score(), $hsp->frac_identical("query");
	  } else {
	    push @{$notes}, sprintf "rejected_high_quality_HSP=%s/%d/%.4f", $key, $hsp->score(), $hsp->frac_identical("query") unless @hsp;
	    # don't add this note if we already have a passed HSP,
	    # e.g. NF1-ACE-01
	  }
	}

	push @hsp, $hsp if $ok;

#      die $hsp->seq_inds("hit", "identical");
	# undef for this implementation??

#	die $hsp->gaps("hit");

	if ($VERBOSE) {
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

  if (@hsp > 1 and not($multi_ok) and $MULTI_HSP_RESCUE) {
    my @pass;
    foreach my $hsp (@hsp) {
      my ($start, $end) = $hsp->range("query");
      my $len = ($end - $start) + 1;
      push @pass, $hsp if $len >= $MULTI_HSP_RESCUE_MIN_LENGTH;
    }

    if (@pass > 1) {
      my @sorted;
      if ($side eq "A") {
	# when anchoring to gene A, choose the rightmost-aligned HSP
	# in the hopes of being closer to B breakpoint,
	# e.g. MYB-QKI-08
	@sorted = sort {$b->end("query") <=> $a->end("query")} @pass;
	# sort by alignment end (proximity to gene B)
      } elsif ($side eq "B") {
	# when anchoring to gene B, choose the leftmost-aligned HSP
	# in the hopes of being closer to A breakpoint,
	# e.g. ARHGAP23-ERCC1-01
	@sorted = sort {$a->start("query") <=> $b->start("query")} @pass;
	# sort by alignment start (proximity to gene A)
      } else {
	die;
      }

      push @{$notes}, sprintf "multi_HSP_rescue=%s/%s/%d/%s",
	$side,
	  $key,
	    scalar(@hsp),
	      join("|", map {join ";", $_->score, join "-", $_->range("query")} @sorted);
      @hsp = $sorted[0];

    }
  }

  push @{$notes}, sprintf "multi_passed_HSP=%s/%s/%d/%s",
    $side,
      $key,
	scalar @hsp,
	  join("|", map {join ";", $_->score, join "-", $_->range("query")} @hsp),
	if @hsp > 1;

  my $retval;
  if (@hsp) {
    #
    #  TO DO: require some minimum score?
    #
    if (@hsp > 1) {
      if ($VERBOSE) {
	foreach my $hsp (@hsp) {
#	  printf STDERR "  multi_hsp:  score:%s strand:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s query_span:%d ref_span:%d total_span=%d query_string=%s hit_string=%s gap_blocks_q=%s gap_blocks_h=%s\n",
	  printf STDERR "     score:%s strand_query:%s strand_hit:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s query_span:%d ref_span:%d total_span=%d query_string=%s hit_string=%s\n",
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
      # should be unpossible due to filtering
      die "error, multiple blast hits";
    }
  }
  return $retval;
}


sub build_genbank_cache {
  my $df = new DelimitedFile("-file" => $f_in,
			     "-headers" => 1,
			    );

  my %all_nm;
  while (my $row = $df->get_hash()) {
    foreach my $end (qw(a b)) {
      my $nm = get_transcript($row, $end);
      $all_nm{$nm} = 1 if $nm;
    }
  }
  $gbc->get([ sort keys %all_nm ]);
}

sub is_itd {
  my ($row) = @_;
  my $gene_a = $row->{genea_symbol} || die;
  my $gene_b = $row->{geneb_symbol} || die;
  return $gene_a eq $gene_b;
}

sub trim_patterns_to_flanking {
  my $trim_to = $FLAGS{"trim-patterns-to-flanking"} || die;
  my $f_in = $FLAGS{file} || die;

  my $df = new DelimitedFile("-file" => $f_in,
			     "-headers" => 1,
			     );
  my $outfile = basename($f_in) . sprintf(".trim_%d.tab", $trim_to);
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
  			      "-auto_qc" => 1,
			     );

  my $verbose = 0;

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
    printf STDERR "start processing %s\n", $row->{pattern} if $verbose;
    my $seq = $row->{sequence} || die;
    printf STDERR "  seq start: %s\n", $seq if $verbose;
    my $idx_start = find_index($seq, "]", "}");
    $idx_start -= $trim_to;
    $seq = substr($seq, $idx_start) if $idx_start > 0;
    printf STDERR "  trim 1: %s\n", $seq if $verbose;

    my $idx_end = find_index($seq, "[", "{");
    $seq = substr($seq, 0, $idx_end + $trim_to + 1);
    printf STDERR "  trim 2: %s\n", $seq if $verbose;

    $row->{sequence} = $seq;
    $rpt->end_row($row);
  }
  $rpt->finish();
}

sub find_index {
  my ($seq, @chars) = @_;
  my $idx;
  foreach my $char (@chars) {
    $idx = index($seq, $char);
    last if $idx >= 0;
  }
  return $idx;
}

sub pair_sanity_check {
  my ($row, $nm_a, $nm_b) = @_;

  my $ok;

  if (my $pair_info = $row->{gene_pair_summary}) {
    # only present in fz2 pattern files
    my @pairs = split /,/, $pair_info;
    foreach my $pair (@pairs) {

      my @sides = parse_pair($pair);

      my ($side_a, $side_b) = @sides;
      my ($pair_nm_a) = (split /\//, $side_a)[1];
      my ($pair_nm_b) = (split /\//, $side_b)[1];

      if ($nm_a eq $pair_nm_a and $nm_b eq $pair_nm_b) {
	$ok = 1;
	last;
      }
    }
  } else {
    # can't check
    $ok = 1;
  }

  return $ok;

}

sub cds_end_check {
  # does the event occur near CDS end exon?
  my ($row, $nm) = @_;
  # get count of exons for this NM

  my $near_cds_end = 0;
  if (my $pair_info = $row->{gene_pair_summary}) {
    # only process if pair annotation available (fz2 input)
    my @pairs = split /,/, $pair_info;
    my %exons;
    foreach my $pair (@pairs) {
      my @pair = parse_pair($pair);
      my @i = split /\//, $pair[1];
      die unless @i == 3;
      my ($b_gene, $b_acc, $b_feature) = @i;
      if ($b_acc eq $nm and $b_feature =~ /^exon_(\d+)/) {
	$exons{$1} = 1;
      }
    }

    if (scalar(keys %exons) == 1) {
      # B end annotated with exactly one exon
      my ($b_exon) = keys %exons;
      # get count of exons in this NM_:
      my ($exon_count, $cds_start_exon, $cds_end_exon) = get_exon_info($nm);
      $near_cds_end = 1 if $b_exon >= $cds_end_exon;
    }
  }
  return $near_cds_end;
}

sub cds_start_check {
  # does the event occur in an exon <= CDS start?
  my ($row, $nm) = @_;
  # get count of exons for this NM

  my $near_cds_start;
  if (my $pair_info = $row->{gene_pair_summary}) {
    # only process if pair annotation available (fz2 input)
    my @pairs = split /,/, $pair_info;
    my %exons;
    foreach my $pair (@pairs) {
      my @pair = parse_pair($pair);
      my @i = split /\//, $pair[0];
      die unless @i == 3;
      my ($a_gene, $a_acc, $a_feature) = @i;
      if ($a_acc eq $nm and $a_feature =~ /^exon_(\d+)/) {
	$exons{$1} = 1;
      }
    }

    if (scalar(keys %exons) == 1) {
      # A start annotated with exactly one exon
      my ($a_exon) = keys %exons;
      # get count of exons in this NM_:
      my ($exon_count, $cds_start_exon) = get_exon_info($nm);
      if ($a_exon <= $cds_start_exon) {
	# event is annotated in an exon <= CDS start exon,
	# so might not be enough upstream sequence for blast to anchor
	$near_cds_start = 1;
      }
    }
  }
  return $near_cds_start;
}

sub parse_pair {
  my ($pair) = @_;
  my @f = split /\-/, $pair;
  if (@f != 2) {
    # parsiing will fail if gene symbol also contains a "-", e.g.:
    # ANKHD1-EIF4EBP3/NM_020690.6/exon_1-ARHGAP26/NM_001135608.3/exon_5
    if ($pair =~ /^(\S+\/\S+\/\S+)\-(.*)$/) {
      @f = ($1, $2);
    } else {
      die "broken pair format $pair";
    }
  }
  return @f;
}

sub nt2protein {
  my ($nt, $strip_end_stop) = @_;
  my $ct = new Bio::Tools::CodonTable();
  my @codons = $nt =~ /(...)/g;
  my $aa = join "", map {$ct->translate($_)} @codons;
  my $result = $aa;
  if ($strip_end_stop) {
    if ($result =~ s/\*$//) {
      # trailing stop codon, expected
    } else {
      # e.g. NM_002537.3:
      #   /note="protein translation is dependent on
      #   polyamine-induced +1 ribosomal frameshift; isoform 1 is
      #   encoded by transcript variant 1; ODC-Az 2"
      $result = undef;
    }
  }
  return $result;
}

sub prune_hits {
  #
  #  prune hits: save only the pairing where the end of match to A
  #  is closest to the start of match to B.
  #
  my ($passed_side) = @_;

  if ($passed_side->{"a"} and $passed_side->{"b"}) {
    # found hits to both sides
    my @f_a = sort keys %{$passed_side->{a}};
    my @f_b = sort keys %{$passed_side->{b}};
    my @combinations;

    foreach my $f_a (@f_a) {
      my $hit_a = $passed_side->{a}{$f_a} || die;
      my $a_end = $hit_a->end("query");
      foreach my $f_b (@f_b) {
	my $hit_b = $passed_side->{b}{$f_b} || die;
	my $b_start = $hit_b->start("query");
	my $distance = abs($b_start - $a_end);

	printf STDERR "%s %s: %s %s %s\n", $f_a, $f_b, $a_end, $b_start, $distance;
	push @combinations, [ $distance, $f_a, $f_b ];
      }
    }

    if (@combinations > 1) {
      # pruning needed
      my @best = sort {$a->[0] <=> $b->[0]} @combinations;
      # TO DO: what if a tie?

      my ($final_a, $final_b) = @{$best[0]}[1,2];

      foreach (["a", $final_a],
	       ["b", $final_b]) {
	my ($side, $wanted) = @{$_};
	foreach my $k (keys %{$passed_side->{$side}}) {
	  if ($k ne $wanted) {
	    # prune
	    delete $passed_side->{$side}{$k};
	    printf STDERR "prune $side $k\n";
	    # TO DO: reporting
	  }
	}
      }
    }
  }
}
