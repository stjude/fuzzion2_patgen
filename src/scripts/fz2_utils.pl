#!/bin/env perl
# fuzzion2-related utilities

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
use List::Util qw(min max sum);

use MiscUtils qw(dump_die build_argv_list get_hash_option average split_list shell_cmd unquote module_dependencies);
use FileUtils qw(read_simple_file find_binary read_simple_file);
use DelimitedFile;
use Reporter;
use ConfigUtils qw(config_or_manual);
use SJPreferredIsoform;
use TdtConfig;
use RefFlatFile;
use Fuzzion2Fuzzall;
use Fuzzion2Hits;
use ReferenceNameMapper;
use TemporaryFileWrangler;
use BLASTer;
use HSPIndexer;
use BAMCoverage;
use BAMUtils qw(get_bam_rnm);
use BAMSet;
use Counter;
use GeneSymbolMapper qw(new_gsm_lite);
use GenomeUtils qw(reverse_complement);

use constant SAM_I_QNAME => 0;
use constant SAM_I_FLAGS => 1;
use constant SAM_I_SEQ => 9;
use constant SAM_I_CIGAR => 5;

use constant SAM_FLAG_RC => 0x10;
use constant SAM_FLAG_SECONDARY_ALIGNMENT => 0x100;
use constant SAM_FLAG_SUPPLEMENTARY_ALIGNMENT => 0x800;
use constant SAM_FLAG_SEGMENT_FIRST => 0x40;
use constant SAM_FLAG_SEGMENT_LAST => 0x80;
# hack

# options for breakpoint BLAST scanning:
my $BP_SCAN_MIN_SOFT_CLIP_LENGTH_TO_CONSIDER = 10;
my $BP_SCAN_MIN_PERFECT_FLANK_TO_MATCH = 10;
# +/- breakpoint site
#my $BP_SCAN_BLAST_FLANK = 30;
my $BP_SCAN_BLAST_FLANK = 40;
# how many flanking nt to use around breakpoint in BLAST query
my $BP_SCAN_BLAST_MIN_IDENTITY = 0.95;
# minimum alignment identity fraction required to use a hit

my $BP_SCAN_MAX_BLAST_DB_SIZE = 250;
# only BLAST a maximum number of reads at a time, to avoid the
# possibility that BLAST will caps match counte (see blastn sefaults
# for "max_target_seqs" and "num_alignemnts")

my $SELF_CHECK_MIN_IDENTITY = 0.95;
# identity of the HSP alignment
my $SELF_CHECK_MIN_OVERLAP_FRAC_HQ = 0.95;
# minimum overlap between the HSP and one of the two patterns
# (one may be longer than the other)

my $SELF_CHECK_MIN_OVERLAP_FRAC_LQ = 0.75;
# lower-quality overlap requirement, for e.g. fusion clutter cases

my $SELF_CHECK_AMBIG_BP_START = "R";
my $SELF_CHECK_AMBIG_BP_END = "Y";
my $SELF_CHECK_MIN_FLANK = 200;


my @MATCH_STRINGS;
my @FA_FIELDS;

my %FLAGS;
my @clopts = (
	      "-file=s",
	      "-files=s",

	      "-describe-breakpoint-intron-sizes",
	      "-genome=s",
	      "-refflat=s",

	      "-recurrent-samples",

	      "-add-no-hit-patterns=s",
	      # to a fuzzall file, append entries for pattern IDs
	      # that are not already included, setting "IDs" to 0

	      "-fuzzall-strong-sample=s",
	      "-min-strong=i",
	      "-only-min-strong=i",
	      # filter ID and sample list to only those showing given
	      # number of strong+ hits

	      "-fuzzall-compare-source-sample=s",
	      "-patterns=s",
	      "-hits=s",
	      "-samples=s",
	      "-restrict-pid=s",

	      # start breakpoint-scan options...
	      "-breakpoint-scan",
	      "-breakpoint-scan-fuzzall=s",
	      "-breakpoint-scan-flank=i" => \$BP_SCAN_BLAST_FLANK,
	      "-bams=s",
	      "-bam=s",
	      "-pattern=s",
	      "-force",
	      "-export-fasta",
	      # write out FASTA for query and db
	      # ...end breakpoint-scan options


	      "-link-patterns",
	      "-patterns-from=s",
	      "-patterns-to=s",

	      "-tag-sjpi=s",

	      "-find-softclip-fusion",
	      # helper for generating fz2 patterns from fusion-supporting
	      # soft-clipped reads
	      "-chr=s",
	      "-pos=i",
	      "-gene=s",
	      "-min-sc-length=i",


	      "-fuzzall-add-sequence=s",
	      # add pattern sequence to a fuzzall report

	      "-count",
	      # various counts of patterns, for manuscript
	      "-count-by-source",

	      "-pattern2fasta=s",
	      "-match=s" => \@MATCH_STRINGS,
	      # restrict pattern to those where IDs match string(s)
	      "-pids=s",
	      # restrict to a file of pattern IDs
	      "-report=s" => \@FA_FIELDS,
	      # include annotation fields in FASTA header line

	      "-fuzzall-filter-to-pattern-list=s",

	      "-find-ambiguous-pattern-symbols",

	      "-pattern-filter-to-pattern-list=s",
	      # excerpt

	      "-pattern-add-interstitial-length=s",
	      # add interstitial sequence size annotation to a pattern file

	      "-pattern-cap-fusion-interstitial-length=i",
	      # patch
	      "-out=s",

	      "-hit-read-compare",
	      "-from=s",
	      "-to=s",

	      "-summary",
	      "-ignore-errors=i",
	      "-group-by-genes=i",
	      # whether to use "-group" option with gene pair in
	      # fuzzion2html / fuzzum

	      "-fuzzall-header-fix=s",
	      # unique-ify headers

	      "-snpshot",
	      # generate screenshots of breakpoints
	      "-tag=s",
	      "-manuscript-colors=i",

	      "-extract-softclip",
	      # extract soft-clipped reads at breakpoint sites
	      "-format=s",
	      "-min-soft-clip-length=i" => \$BP_SCAN_MIN_SOFT_CLIP_LENGTH_TO_CONSIDER,

	      "-bambino",

	      "-source-traceback",
	      # list of patterns to search for source record
	      "-sources=s",
	      # list of source files

	      "-verbose",

	      "-traceback-compare=s",
	      # check new patterns derived from traceback files with
	      # an existing pattern set
	      # input = list of link patterns
	      "-patterns-orig=s",
	      "-patterns-new=s",
	      "-fuzzall-orig=s",
	      "-fuzzall-new=s",
	      "-restrict=s",

	      "-false-negative-check=s",
	      "-fuzzall-pair=s",

	      "-hit2tsv=s",

	      "-add-preferred-summary=s",

	      "-patch=s",
	      "-show-module-dependencies",

	      "-fuzzum-pair-post",
	      # add additional annotations to fuzzum by gene pair
	      # e.g. SJCOGALL010876_D3_by_gene.txt

	      "-generate-minimal-cicero",

	      "-patch-breakpoints-traceback",
	      # patch genomic breakpoint info from source records,
	      # via -source-traceback output
	      "-patch-breakpoints-by-features",
	      # patch genomic breakpoint info based on annotated
	      # transcripts and exons (hack, primarily for Kriwacki
	      # and other data lacking this information in source)


	      "-fuzzum-merge",
	      # merge fuzzum results by sample and pattern/pair
	      # into an existing report containing those keys
	      "-glob=s",
	      # to match set of fuzzum file to process
	      "-target=s",
	      # file to merge into
	      "-target-sample=s",
	      "-target-id=s",
	      "-header-prefix=s",
	      # prefix to add to fuzzum columns merged in
	      "-glob-fuzzum-pattern=s",
	      # pattern-level fuzzum results

	      "-debug",

	      "-find-sj-cohort-bams",

	      "-clean-hits",
	      "-only-category=s",

	      "-sv6=s",
	      "-sample=s",
	      "-genea=s",
	      "-geneb=s",

	      "-pattern-self-identity",
	      "-index-start=i",
	      "-index-end=i",
	      "-index-out=s",
	      "-pattern-self-identity-post=s",
	      # remove duplicates

	      "-clean-pattern-categories",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  ) || die "error parsing command-line parameters";

if ($FLAGS{"describe-breakpoint-intron-sizes"}) {
  describe_breakpoint_intron_sizes();
} elsif ($FLAGS{"show-module-dependencies"}) {
  module_dependencies();
  exit(0);
} elsif ($FLAGS{"recurrent-samples"}) {
  report_recurrent_samples();
} elsif ($FLAGS{"add-no-hit-patterns"}) {
  add_no_hit_patterns();
} elsif ($FLAGS{"fuzzall-strong-sample"}) {
  fuzzall_strong_sample();
} elsif ($FLAGS{"fuzzall-compare-source-sample"}) {
  fuzzall_compare_source_sample();
} elsif ($FLAGS{"breakpoint-scan"}) {
  breakpoint_scan();
} elsif ($FLAGS{"find-softclip-fusion"}) {
  find_softclip_fusion();
} elsif ($FLAGS{"breakpoint-scan-fuzzall"}) {
  breakpoint_scan_fuzzall();
} elsif ($FLAGS{"link-patterns"}) {
  link_patterns();
} elsif ($FLAGS{"tag-sjpi"}) {
  tag_sjpi();
} elsif ($FLAGS{"fuzzall-add-sequence"}) {
  fuzzall_add_sequence();
} elsif ($FLAGS{"count"}) {
  count_patterns();
} elsif ($FLAGS{"count-by-source"}) {
  count_by_source();
} elsif ($FLAGS{"pattern2fasta"}) {
  pattern2fasta();
} elsif ($FLAGS{"fuzzall-filter-to-pattern-list"}) {
  fuzzall_filter_to_patterns();
} elsif ($FLAGS{"find-ambiguous-pattern-symbols"}) {
  find_ambiguous_pattern_symbols();
} elsif ($FLAGS{"pattern-filter-to-pattern-list"}) {
  filter_pattern_set();
} elsif ($FLAGS{"pattern-add-interstitial-length"}) {
  pattern_add_interstitial_length();
} elsif ($FLAGS{"pattern-cap-fusion-interstitial-length"}) {
  pattern_cap_fusion_interstitial_length();
} elsif ($FLAGS{"hit-read-compare"}) {
  hit_read_compare();
} elsif ($FLAGS{"summary"}) {
  build_summary();
} elsif ($FLAGS{"fuzzall-header-fix"}) {
  fuzzall_header_fix();
} elsif ($FLAGS{"snpshot"}) {
  generate_snpshots();
} elsif ($FLAGS{"extract-softclip"}) {
  extract_soft_clipped();
} elsif ($FLAGS{"bambino"}) {
  generate_bambino_viewer_scripts();
} elsif ($FLAGS{"source-traceback"}) {
  source_traceback();
} elsif ($FLAGS{"traceback-compare"}) {
  traceback_compare();
} elsif ($FLAGS{"false-negative-check"}) {
  false_negative_check();
} elsif ($FLAGS{"hit2tsv"}) {
  hit2tsv();
} elsif ($FLAGS{"add-preferred-summary"}) {
  # add filtered version of gene_pair_summary restricting to
  # results for preferred pair only
  add_preferred_summary();
} elsif ($FLAGS{"patch"}) {
  patch_patterns();
} elsif ($FLAGS{"fuzzum-pair-post"}) {
  fuzzum_pair_post();
} elsif ($FLAGS{"generate-minimal-cicero"}) {
  # skeleton/template file for contig-based pattern generation
  generate_minimal_cicero();
} elsif ($FLAGS{"patch-breakpoints-traceback"}) {
  patch_breakpoints_traceback();
} elsif ($FLAGS{"patch-breakpoints-by-features"}) {
  patch_breakpoints_by_features();
} elsif ($FLAGS{"fuzzum-merge"}) {
  fuzzum_merge();
} elsif ($FLAGS{"find-sj-cohort-bams"}) {
  find_sj_cohort_bams();
} elsif ($FLAGS{"clean-hits"}) {
  clean_hits();
} elsif ($FLAGS{"pattern-self-identity"}) {
  pattern_self_identity();
} elsif ($FLAGS{"pattern-self-identity-post"}) {
  pattern_self_identity_post();
} elsif ($FLAGS{"clean-pattern-categories"}) {
  clean_pattern_categories();
} else {
  die "unspecified command-line option";
}
exit(0);

sub describe_breakpoint_intron_sizes {
  # report sizes of introns downstream of geneA and upstream of geneB
  # for JZ/Karol 6/2022
  my $f_patterns = $FLAGS{file} || die "-file";
  my $genome = $FLAGS{genome} || die "-genome";
  my $f_refflat = $FLAGS{refflat} || die "-refflat";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $f_sjpi = $config_genome->{GENE_TRANSCRIPT_MATRIX} || die;
  my $sjpi = new SJPreferredIsoform("-file" => $f_sjpi,
				    "-auto_gsm" => 1);

  my $rff = new RefFlatFile();
  $rff->canonical_references_only(1);
  $rff->strip_sharp_annotations(1);
  $rff->parse_file(
		  "-refflat" => $f_refflat,
		  "-type" => "refgene",
		  "-generate-introns" => 1,
		 );

  my $df = new DelimitedFile("-file" => $f_patterns,
			     "-headers" => 1,
			     );
  my $outfile = basename($f_patterns) . ".intron_info.tab";

  my $rpt = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   pattern

					   geneA
					   accA
					   featureA
					   intron_name_A
					   intron_length_A
					   intron_region_A


					   geneB
					   accB
					   featureB
					   intron_name_B
					   intron_length_B
					   intron_region_B
					)
				      ],
			 "-auto_qc" => 1,
			);

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  my $found;
  while (my $row = $df->get_hash()) {
    my @summaries = split /,/, $row->{gene_pair_summary} || die;
    foreach my $summary (@summaries) {
      my @f = split /\-/, $summary;
      die unless @f == 2;
      my ($from, $to) = @f;
      @f = split /\//, $from;
      die unless @f == 3;
      my ($from_gene, $from_acc, $from_feature) = @f;

      @f = split /\//, $to;
      die unless @f == 3;
      my ($to_gene, $to_acc, $to_feature) = @f;

      if ($sjpi->is_preferred_nm($from_acc) and
	  $sjpi->is_preferred_nm($to_acc)) {
	# usable event:
	# pattern derived from SJ preferred isoforms on both sides
	$found = 1;

	my $from_wanted_feature = get_target_intron($from_feature, 0);
	# for gene A, the downstream intron # should be the same
	# as the breakpoint exon

	my $to_wanted_feature = get_target_intron($to_feature, -1);
	# for gene B, the upstream intron # should be one less

	my $i_from = get_intron_info($rff, $from_acc, $from_wanted_feature) || die;
	my $i_to = get_intron_info($rff, $to_acc, $to_wanted_feature) || die;

	my %r;
	$r{pattern} = $row->{pattern} || die;
	$r{geneA} = $from_gene;
	$r{accA} = $from_acc;
	$r{featureA} = $from_feature;
	$r{intron_name_A} = $from_wanted_feature;
	populate_intron_info(\%r, "A", $i_from);


	$r{geneB} = $to_gene;
	$r{accB} = $to_acc;
	$r{featureB} = $to_feature;
	$r{intron_name_B} = $to_wanted_feature;
	populate_intron_info(\%r, "B", $i_to);

	$rpt->end_row(\%r);
      }
      die "no matches found" unless $found;
    }
  }
  $rpt->finish();

}


sub get_target_intron {
  my ($exon_feature, $wanted_intron_offset) = @_;

  $exon_feature =~ /exon_(\d+)$/ || die "$exon_feature not exon";
  my $ex_no = $1;

  my $wanted_feature = sprintf 'intron_%d', $ex_no + $wanted_intron_offset;
  return $wanted_feature;
}

sub get_intron_info {
  my ($rff, $acc, $wanted_feature) = @_;

  my $rf_set = $rff->find_by_accession($acc) || die "can't find refflat entry for $acc";

  my $rf = $rf_set->[0];
  # there may be multiple mappings, but a single one will do

  my $result;

  foreach my $intron (@{$rf->{introns}}) {
    # iterate through intron regions, looking for the target feature
    my $pos = $intron->{start};

    my ($feature, $feature_number, $strand) =
      $rff->get_annotation_for_position(
				       "-row" => $rf,
				       "-base" => $pos,
				       "-extended" => 1,
				      );

    my $tag = join "_", $feature, $feature_number;
    if ($tag eq $wanted_feature) {
      my %result = %{$intron};
      $result{chrom} = $rf->{chrom} || die;
      $result = \%result;
    }
  }

  return $result;
}

sub populate_intron_info {
  my ($row, $tag, $ii) = @_;

  my ($chrom, $start, $end) = @{$ii}{qw(chrom start end)};
  die unless $chrom and $start and $end;

  my $tag_l = sprintf 'intron_length_%s', $tag;
  $row->{$tag_l} = ($end - $start) + 1;

  my $tag_r = sprintf 'intron_region_%s', $tag;
  $row->{$tag_r} = sprintf '%s:%d-%d', $chrom, $start, $end;
}

sub report_recurrent_samples {
  # identify what fraction of the records each unique sample ID appears in
  my $f_in = $FLAGS{file} || die;
  my $df = new DelimitedFile("-file" => $f_in,
			     "-headers" => 1,
			     );

  # find all IDs mentioned in report and count:
  my %id_count;
  my $row_count = 0;
  while (my $row = $df->get_hash()) {
    $row_count++;
    my @ids = get_id_list($row);
    foreach (@ids) {
      $id_count{$_}++;
    }
  }

  my $f_out = $FLAGS{out} || "sample_summary.tab";
  my $rpt = new Reporter(
			 "-file" => $f_out,
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   id
					   record_count
					   set_fraction
					)
				      ],
			 "-auto_qc" => 1,
			);

  foreach my $id (sort keys %id_count) {
    my $count = $id_count{$id};
    my $frac = $count / $row_count;
    my %r;
    $r{id} = $id;
    $r{record_count} = $count;
    $r{set_fraction} = $frac;
    $rpt->end_row(\%r);
  }
  $rpt->finish();

}

sub get_id_list {
  my ($row) = @_;
  my @set = split /,\s*/, $row->{"ID list"} || die;
  my @ids;
  foreach my $record (@set) {
    $record =~ /^(\S+)\(\d+\/\d+\)$/ || die;
    push @ids, $1;
  }
  return @ids;
}

sub add_no_hit_patterns {
  # - add:
  #   - a blank result record for patterns not hit in a given set
  #     (i.e. completely missed in fz2 run)
  #  - for each pattern, interstitial sequence length (helps w/debugging)
  my ($f_in, $f_patterns) = @_;

  unless ($f_in) {
    $f_in = $FLAGS{file} || die "-file";
    # fuzzall file
    $f_patterns = $FLAGS{"add-no-hit-patterns"} || die;
    # patterns to merge in, if not already present
  }

  my $df = new DelimitedFile("-file" => $f_patterns,
			     "-headers" => 1,
			     );
  #
  # load all patterns in source set:
  #
  my %all_pid;
  while (my $row = $df->get_hash()) {
    my $pid = $row->{pattern} || die;
    $all_pid{$pid} = $row;
  }

  my $pid2ilen = get_pattern2ilen($f_patterns);

  my $fa = new Fuzzion2Fuzzall("-file" => $f_in);
  my $f_out = basename($f_in) . ".add_no_hit.tab";

  my $outfile = basename($f_in) . ".add_no_hit.tab";
  my $rpt = $fa->df->get_reporter(
				  "-file" => $f_out,
				  "-auto_qc" => 1,
				  "-extra" => [
					       qw(
						   pattern_interstitial_bases
						)
					      ]
			     );

  my $example_row;
  my %need_pid = map {$_, 1} keys %all_pid;

  my $f_pid = $fa->f_pattern() || die;

  #
  # copy source file, recording observed pattern IDs:
  #
  while (my $row = $fa->next()) {
    my $pid = $row->{$f_pid} || dump_die($row, "no field $f_pid");
    delete $need_pid{$pid};
    $example_row = $row;
    $row->{pattern_interstitial_bases} = $pid2ilen->{$pid};
    $rpt->end_row($row);
  }

  #
  #  append records for missing patterns:
  #
  foreach my $pid (sort keys %all_pid) {
    if ($need_pid{$pid}) {
      my $row = $all_pid{$pid};
      my %r = %{$row};
      foreach my $f (keys %{$example_row}) {
	unless (defined $r{$f}) {
	  $r{$f} = "";
#	  print STDERR "$f\n";
	}
      }
      $r{$f_pid} = $pid;
      $r{IDs} = 0;
      $r{pattern_interstitial_bases} = $pid2ilen->{$pid};
      $rpt->end_row(\%r);
    }
  }
  $rpt->finish();
}


sub fuzzall_strong_sample {
  # filter a fuzzall report:
  #  - filter ID count to only reads with a minimum number of strong reads
  #  - optionally filter matching ID details by same minimum count
  #
  # TO DO: add option to have a strong- policy as well
  my ($f_in) = @_;
  unless ($f_in) {
    $f_in = $FLAGS{"fuzzall-strong-sample"} || die;
  }
  my $min_strong = $FLAGS{"min-strong"} || die "-min-strong";
  my $only_strong = $FLAGS{"only-min-strong"};
  die "-only-min-strong [0|1]" unless defined $only_strong;

  my $f_out = basename($f_in) . ".min_strong.tab";

  my $df = new DelimitedFile("-file" => $f_in,
			     "-headers" => 1,
			     );
  my $rpt = $df->get_reporter(
			      "-file" => $f_out,
  			      "-auto_qc" => 1,
			     );

  my $f_ids = "IDs";
  my $f_id_list = "ID list";

  my @rows_out;
  while (my $row = $df->get_hash()) {
    my $id_list = $row->{$f_id_list} || die;
    my @list = split /,\s+/, $id_list;
    my %by_strong;
    foreach my $entry (@list) {
      $entry =~ /\(\d+\/(\d+)\)/ || die;
      my $strong = $1;
      push @{$by_strong{$strong}}, $entry;
    }

    my @counts_for_ids;
    # counts we want to include in the "IDs" total
    my @counts_for_sample_detail;

    my @counts_sorted = sort {$b <=> $a} keys %by_strong;
    foreach my $count (@counts_sorted) {
      push @counts_for_ids, $count if $count >= $min_strong;

      push @counts_for_sample_detail, $count if $only_strong ? $count >= $min_strong : 1;
    }
#    printf STDERR "filtered: %s\n", join ",", @counts_for_sample_detail if @counts_for_sample_detail;

    my @all = map {@{$by_strong{$_}}} @counts_for_ids;
    my @detail = map {@{$by_strong{$_}}} @counts_for_sample_detail;

    $row->{$f_ids} = scalar @all;
    # ID count is filtered to samples showing minimum strong+ sample count
    $row->{$f_id_list} = join ", ", @detail;
    # sample detail list sorted by strong count

    push @rows_out, $row;
  }

  # write output sorted by ID count:
  foreach my $row (sort {$b->{IDs} <=> $a->{IDs}} @rows_out) {
    $rpt->end_row($row);
  }

  $rpt->finish();
}


sub fuzzall_compare_source_sample {
  #
  # check fuzzall output for evidence of matches in source samples
  #
  my $f_in = $FLAGS{"fuzzall-compare-source-sample"} || die;
  # fuzzall report, grouped by pattern, not pair
  my $f_patterns = $FLAGS{patterns} || die "-patterns FILE";

  my $min_strong_standard = $FLAGS{"min-strong"};
  die "specify -min-strong VALUE" unless defined $min_strong_standard;
  # supporting reads to consider to have solid evidence
  my $min_strong_low = 1;
  # still report low-support cases, but tag specially

  my $pid_restrict;
  if (my $f_pid_restrict = $FLAGS{"restrict-pid"}) {
    $pid_restrict = read_simple_file($f_pid_restrict, "-hash1" => 1);
  }

  my $f_samples = $FLAGS{samples};
  my $f_hits = $FLAGS{hits};
  my $f_bams = $FLAGS{bams};

  my $pid2ilen = get_pattern2ilen($f_patterns);
  # index of pattern name -> interstitial bases count

  #
  #  initialize analyzed samples set:
  #
  my $all_samples;
  if ($f_samples) {
    # pre-parsed to sample ID
    $all_samples = read_simple_file($f_samples, "-hash1" => 1);
  } elsif ($f_bams) {
    # list of available BAMs.
    # possibly sub-optimal because the list actually used in the RUN
    # may be different!  For this reason a list of hit files may be
    # preferable.
    $all_samples = {};
    my $list = read_simple_file($f_bams);
    foreach my $f (@{$list}) {
      my $sample = basename($f);
      $sample =~ s/\..*// || die;
      $all_samples->{$sample} = 1;
    }
  } elsif ($f_hits) {
    # list of fuzzion2 hits files
    $all_samples = {};
    my $list = read_simple_file($f_hits);
    foreach my $f (@{$list}) {
      my $sample = basename($f);
      if ($sample =~ s/_by_pattern\.txt$//) {
	# OK: Steve Rice standard format
      } elsif ($sample =~ s/\.hits$//) {
	# OK: MNE wrapper format
      } else {
	die "can't identify hits file suffix in $sample";
      }
      $all_samples->{$sample} = 1;
    }
  } else {
    die "specify -hits LISTFILE|-samples LISTFILE|-bams LISTFILE";
  }

  my %all_samples_by_base;
  foreach my $sample (keys %{$all_samples}) {
    my $base = sj_sample_base($sample);
    push @{$all_samples_by_base{$base}}, $sample;
  }

  my $fa = new Fuzzion2Fuzzall("-file" => $f_in);
  my $f_pattern = $fa->f_pattern || die;

  my %found_direct;
  # direct sample->pid count
  my %found_pair;
  # 1-to-many bucket of counts by sample->pair
  my %found_source_sample2pair2pid;
  # sample -> PID counts, for lookup of alternate pattern for same
  # source sample

  #
  # first pass: collate info
  #
  while (my $row = $fa->next()) {
    my $pid = $row->{$f_pattern} || die;
    next if $pid_restrict and !$pid_restrict->{$pid};

    my $gene_pair = pid2pair($pid);

    # rewrite:
    # - convert source samples to cooked format FIRST
    # - while parsing each ID, handle/track source sample patterns and
    #   non-source-sample patterns
    # - get rid of @samples_found, won't need

    my $sample_raw = $row->{sample} || "";
    # sometimes no sample provided, e.g. some Kriwacki records
    my @source_samples = split /,/, $sample_raw;
    my %source_samples;
    foreach my $sample_src (@source_samples) {
      my $sample_mapped = sample_xref($sample_src, $all_samples, \%all_samples_by_base);
      # translate source sample name to the run sample list
      $source_samples{$sample_mapped} = 1;
    }

    foreach my $r (@{$row->{IDs_parsed}}) {
      # if a minimum match is found between a pattern and a sample,
      # track the sample + gene pair hit.  Some results gravitate to
      # other patterns for the pair.
      my $csp = $r->{count_strong_plus};
      # this might actually be strong rather than strong+

      if ($csp) {
	my $sample_mapped = sample_xref($r->{id}, $all_samples, \%all_samples_by_base);
	$found_pair{$sample_mapped}{$gene_pair} += $csp;
	# OK to simply combine counts since by default fz2 will assign
	# reads to the "best" patterm only
	printf STDERR "save pair: %s %s %s %s\n", $sample_mapped, $pid, $gene_pair, $csp;

	if ($source_samples{$sample_mapped}) {
	  # this result is for a source sample for the pattern
	  #	  die $sample_mapped;
	  $found_direct{$sample_mapped}{$pid} = $csp;
	  printf STDERR "save direct: %s %s %s %s\n", $sample_mapped, $pid, $gene_pair, $csp;

	  $found_source_sample2pair2pid{$sample_mapped}{$gene_pair}{$pid} = $csp;
	}
      }
    }
  }

  #
  #  second pass: report results
  #
  $fa = new Fuzzion2Fuzzall("-file" => $f_in);
  my $f_out = basename($f_in) . ".source_check.tab";
  my $rpt = $fa->df()->get_reporter(
			      "-file" => $f_out,
			      "-extra" => [
					   qw(
					       in_sample_list
					       sample_link
					       pattern_interstitial_bases
					       found_type
					       found_count
					       found_direct
					       found_alt_pattern_for_src
					       found_pair
					    )
					  ],
  			      "-auto_qc" => 1,
			     );

  while (my $row = $fa->next()) {
    my $pid = $row->{$f_pattern} || die;
    next if $pid_restrict and !$pid_restrict->{$pid};
    my $gene_pair = pid2pair($pid);
    my @samples = sort split /,/, $row->{sample};
    # source sample, may be blank

    my @found_direct;
    my @found_pair;

    my @count_direct;
    my @count_pair;

    my @found_alt_pattern_for_src;
    my @count_alt_pattern_for_src;

    foreach my $sample_in (@samples) {
      my $sample_mapped = sample_xref($sample_in, $all_samples, \%all_samples_by_base);

      if (my $count = $found_direct{$sample_mapped}{$pid}) {
	# direct hit for this sample + pid
	push @found_direct, sprintf "%s=%d", $sample_mapped, $count;
	push @count_direct, $count;
      }

      if (my $alt_src_patterns = $found_source_sample2pair2pid{$sample_mapped}{$gene_pair}) {
	my %alt = %{$alt_src_patterns};
	delete $alt{$pid};
	if (%alt) {
	  # alternative patterns generated from the same source sample
	  # were found in the results.
	  my $best_pid;
	  my $best_count = 0;
	  foreach my $p (keys %alt) {
	    my $count = $alt{$p};
	    if ($count > $best_count) {
	      $best_count = $count;
	      $best_pid = $p;
	    }
	  }

	  push @found_alt_pattern_for_src, sprintf "%s/%s=%d", $sample_mapped, $best_pid, $best_count;
	  push @count_alt_pattern_for_src, $best_count;

#	  die "to do: check for hit to same sample but different PID (alt isoform), $pid $sample_mapped max=" . $best_count . " " . join ",", sort keys %{$alt_src_patterns};
	}
      }

      if (my $count = $found_pair{$sample_mapped}{$gene_pair}) {
	# indirect hit for this sample + gene pair
	push @found_pair, sprintf "%s=%d", $sample_mapped, $count;
	push @count_pair, $count;
      }
    }

    my ($best_count_direct) = sort {$b <=> $a} @count_direct;
    my ($best_count_pair) = sort {$b <=> $a} @count_pair;
    my ($best_count_alt_pattern_for_src) = sort {$b <=> $a} @count_alt_pattern_for_src;
    # may be multiple source samples, use best result
    $best_count_direct = 0 unless defined $best_count_direct;
    $best_count_pair = 0 unless defined $best_count_pair;
    $best_count_alt_pattern_for_src = 0 unless defined $best_count_alt_pattern_for_src;

    my %found_type;
    $found_type{direct} = $best_count_direct;
    $found_type{gene_pair} = $best_count_pair;
    $found_type{alt_pattern_for_src} = $best_count_alt_pattern_for_src;
    my ($found_type) = sort {$found_type{$b} <=> $found_type{$a}} keys %found_type;

    my $found_count = $found_type{$found_type};

    if ($found_count) {
      $found_type .= "_low" if $found_count < $min_strong_standard;
    } else{
      $found_type = "";
    }

    $row->{found_type} = $found_type;
    $row->{found_count} = $found_count;
    $row->{found_direct} = join ",", @found_direct;
    $row->{found_pair} = join ",", @found_pair;
    $row->{pattern_interstitial_bases} = $pid2ilen->{$pid};

    #
    #  check if source sample(s) present in sample run list:
    #
    my $found = 0;
    my $not_found = 0;
    my @sample_link;
    foreach my $sample_in (@samples) {
      my $sample_mapped = sample_xref($sample_in, $all_samples, \%all_samples_by_base);
      push @sample_link, sprintf "%s=%s", $sample_in, $sample_mapped if $sample_in ne $sample_mapped;
      if ($all_samples->{$sample_mapped}) {
	$found = 1;
      } else {
	$not_found = 1;
      }
    }
    $row->{sample_link} = join ",", @sample_link;

    my $in_sample_list = "";
    if ($found and $not_found) {
      $in_sample_list = "some";
    } elsif ($found) {
      $in_sample_list = "y";
    } elsif ($not_found) {
      $in_sample_list = "n";
    } elsif (@samples == 0) {
      # no usable sample data
      $in_sample_list = "n/a";
    } else {
      die "unhandled";
    }
    $row->{in_sample_list} = $in_sample_list;

    # $row->{found_alt_pattern_for_src} = $best_count_alt_pattern_for_src;
    $row->{found_alt_pattern_for_src} = join ",", @found_alt_pattern_for_src;

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

sub pid2pattern_number {
  my ($pid) = @_;
  my $gene_pair = $pid;
  $gene_pair =~ s/\-(\d+)$// || die;
  return int($1);
}

sub sample_xref {
  my ($sample_in, $all_samples, $base2samples) = @_;
  my $result;
  if ($all_samples->{$sample_in}) {
    # direct mapping exists, OK
    $result = $sample_in;
  } else {
    my $base = sj_sample_base($sample_in);
    if (my $set = $base2samples->{$base}) {
      my @sorted = sort {$a cmp $b} @{$set};
      # or maybe in reverse order to prefer e.g. D2 over D1?

      if ($sample_in =~ /_/) {
	my @hits = grep {index($_, $sample_in) != -1} @{$set};
	if (@hits) {
	  # prefer named version if present,
	  # e.g. SJST032197_D1_G1 -> SJST032197_D1
	  $result = $hits[0];
	}
      }

      $result = $sorted[0] unless $result;
      # HACK, what if one sample is better than the other??
    }
  }
  $result = $sample_in unless $result;

  return $result;
}

sub sj_sample_base {
  my ($sample) = @_;
  my $base = $sample;
  $base =~ s/_.*//;
  # obviously won't work for non-SJ samples
  return $base;
}


sub get_pattern2ilen {
  #
  #  generate map of pattern name -> interstitial bases count
  #
  my ($f_patterns) = @_;
  die unless $f_patterns;
  my $df = new DelimitedFile("-file" => $f_patterns,
			     "-headers" => 1,
			    );

  my %pattern2ilen;
  while (my $row = $df->get_hash()) {
    my $pid = $row->{pattern} || die;
    my $seq = $row->{sequence} || die;

    foreach my $style (
		       [ "]", "[" ],
		       [ "}", "{" ]
		      ) {
      my ($c1, $c2) = @{$style};
      my $i1 = index($seq, $c1);
      my $i2 = index($seq, $c2);
      if ($i1 != -1 and $i2 != -1) {
	die unless $i2 > $i1;
	my $len = ($i2 - $i1) - 1;
	$pattern2ilen{$pid} = $len;
      }
    }
  }
  return \%pattern2ilen;
}

sub breakpoint_scan {
  # single job
  # call a different/reusable sub to do the actual work once
  # all required inputs are determined
  my $f_patterns = $FLAGS{patterns} || die "-patterns";
  my $pid = $FLAGS{pattern} || die "-pattern";
  my $bam = $FLAGS{bam} || die "-bam";

  my $patterns = load_pattern_file($f_patterns);

  my $sample = bam_to_sample($bam);

  breakpoint_scan_single(
			 "-bam" => $bam,
			 "-pid" => $pid,
			 "-pattern-row" => $patterns->{$pid},
			 "-sample" => $sample
			);
}

sub breakpoint_scan_fuzzall {
  # run on a set of BAMs listed in a fuzzall report
  my $f_fuzzall = $FLAGS{"breakpoint-scan-fuzzall"};
  my $f_patterns = $FLAGS{patterns} || die "-patterns";
  my $pid = $FLAGS{pattern} || die "-pattern";
  my $f_bams = $FLAGS{bams} || die "-bams";
  my $export_fasta = $FLAGS{"export-fasta"};

  printf STDERR "configuration:\n";
  printf STDERR "  pattern excerpt flank length: %d\n", $BP_SCAN_BLAST_FLANK;
  printf STDERR "  export soft-clipped reads: %s\n",
    $export_fasta ? "y" : "n";
  print STDERR "\n";

  my $f_out = sprintf '%s.fuzzall.bam_scan.tab', $pid;
  my $rpt = get_bam_scan_reporter($pid, $f_out);
  my $bs = new BAMSet("-file" => $f_bams);

  my $fa = new Fuzzion2Fuzzall("-file" => $f_fuzzall);
  my $f_pattern = $fa->f_pattern() || die;
  my $hit_row;
  while (my $row = $fa->next()) {
    if ($row->{$f_pattern} eq $pid) {
      $hit_row = $row;
      last;
    }
  }
  die "can't find fuzzall entry for $pid" unless $hit_row;

  my $patterns = load_pattern_file($f_patterns);
  my $pr = $patterns->{$pid} || die "can't find pattern info for $pid";

  my $bam_not_found = 0;
  my $id_set = $hit_row->{IDs_parsed};
  my $c = new Counter($id_set);
  foreach my $id_ref (@{$id_set}) {
    my $sample = $id_ref->{id} || die;
    $c->next($sample);
    if (my $set = $bs->find($sample)) {
      foreach my $bam (@{$set}) {
	# might be more than one
	breakpoint_scan_single(
			       "-bam" => $bam,
			       "-rpt" => $rpt,
			       "-pid" => $pid,
			       "-pattern-row" => $pr,
			       "-sample" => $sample
			      );
      }
    } else {
      $bam_not_found++;
    }
  }
  printf STDERR "records with no BAMS available: %d\n", $bam_not_found;

  $rpt->finish();
}


sub breakpoint_scan_single {
  # maybe expand inputs to allow a larger/continuous report?
  # - target a whole gene if specific breakpoint not known
  #
  # IDEAS/TO DO:
  # - java Bambino command line for breakpoints A/B?
  #
  my %options = @_;

  my $pid = $options{"-pid"} || die;
  my $pr = $options{"-pattern-row"} || die;
  my $f_bam = $options{"-bam"} || die;
  my $sample = $options{"-sample"} || die "-sample";
  my $rpt_passthrough = $options{"-rpt"};
  my $cache = $FLAGS{force} ? 0 : 1;

  #
  #  init reference name disambiguation:
  #
  find_binary("samtools", "-die" => 1);
  find_binary("blastn", "-die" => 1);

  my $rnm = get_bam_rnm($f_bam);
  my $bn = basename($f_bam);
  $bn =~ s/\.bam$//i || die;
  my $f_out = sprintf "%s.%s.%s.bam_scan.tab",
    $pid,
    $bn,
    md5_hex(realpath $f_bam);
  # add MD5 element in case there are multiple BAM files with
  # the same basename for whatever reason

  if (not($cache and -s $f_out)) {
    #
    #  analysis needed
    #
    my $rpt = get_bam_scan_reporter($pid, $f_out);

    my %r = %{$pr};
    $r{sample} = $sample;

    #
    # is the sample being processed one of the source samples for this pattern?
    #
    my $source_samples = $pr->{sample} || die "no source samples";
    my @samples = split /,/, $source_samples;
    my %s = map {$_, 1} @samples;
    my $is_source_sample = 0;
    if ($s{$sample}) {
      # perfect match
      $is_source_sample = 1;
    } else {
      # also scan for string match, some BAM filenames may have
      # additional text in basenames.  Need an example to confirm.
      foreach my $s (@samples) {
	if (index($sample, $s) == 0 or index($s, $sample) == 0) {
	  die "TEST ME: sample $sample hits string $s";
	}
      }
    }
    $r{is_source_sample} = $is_source_sample;

    foreach my $end (qw(a b)) {
      #
      #  get chrom name in BAM and breakpoint:
      #
      my $chr_raw = $pr->{sprintf "gene%s_chr", $end} || die;
      my $chr_bam = $rnm->find_name($chr_raw) || die;

      my $pos_raw = $pr->{sprintf "gene%s_pos", $end} || die;
      my @pos = split /,/, $pos_raw;
      my $pos;
      if (@pos > 1) {
	$pos = int(average(\@pos));
      } else {
	$pos = $pos_raw;
      }

      #
      #  get coverage using standard method:
      #
      my $bc = new BAMCoverage(
			       "-bam" => $f_bam,
			      );
      $bc->allow_optical_pcr_duplicates(1);
      $bc->include_soft_clips(1);

      my $coverage_official = $bc->get_coverage(
						"-chr" => $chr_bam,
						"-pos" => $pos,
					       );
      # TO DO: maybe also check +/- 1 to see if higher?
      $r{sprintf "gene%s_coverage", $end} = $coverage_official;

      #
      #  get SAM reads overlapping this region:
      #
      my $cmd = sprintf "samtools view %s %s:%d-%d|", $f_bam, $chr_bam, $pos, $pos;
      #    printf STDERR "cmd: %s\n", $cmd;
      open(SAMTMP, $cmd) || die;

      my $count_reads_sc = 0;
      # count of reads meeting minimum soft-clipping requirement
      my %sc;
      my %sc_len;

      while (<SAMTMP>) {
	chomp;
	my @f = split /\t/, $_;
	my $cigar = $f[SAM_I_CIGAR];
	my @softs = $cigar =~ /(\d+)S/g;
	my $soft_clipped;
	foreach my $soft (@softs) {
	  $soft_clipped = 1 if $soft >= $BP_SCAN_MIN_SOFT_CLIP_LENGTH_TO_CONSIDER;
	}
	if ($soft_clipped) {
	  $count_reads_sc++;

	  my $read_id = $f[SAM_I_QNAME];
	  my $seq = $f[SAM_I_SEQ];

	  if ($sc{$read_id}) {
	    $read_id .= ".1";
	    # hack: ever breaks?
	  }
	  die "duplicate ID $read_id" if $sc{$read_id};
	  $sc{$read_id} = $seq;

	  $sc_len{$read_id} = max(@softs);
	}
      }

      my $supporting_count = 0;
      my $supporting_min_overlap = "";
      my %sc_leftover = %sc;

      if (%sc) {
	#
	# compare pattern with soft-clipped reads:
	#
	my $hits_match = pattern_blast_check(
					     "-pid" => $pid,
					     "-pattern" => $pr->{sequence},
					     "-db" => \%sc,
					     "-end" => $end,
					     "-sample" => $sample,
					     "-bam" => $f_bam
					    );
	$supporting_count = scalar keys %{$hits_match};

	my %detail;
	foreach my $id (sort keys %{$hits_match}) {
	  delete $sc_leftover{$id};
	  my $length = $hits_match->{$id};
#	  printf STDERR "ID %s overlap %d\n", $id, $length;
	  $detail{$length}++;
	}

	my @lengths;
	foreach my $length (sort {$b <=> $a} keys %detail) {
	  push @lengths, sprintf "%d=%d", $length, $detail{$length};
	}

	$supporting_min_overlap = sprintf "%s", join ",", @lengths;
      }
      $r{sprintf "gene%s_coverage_supporting", $end} = $supporting_count;
      $r{sprintf "gene%s_coverage_supporting_shortest_side_overlap", $end} = $supporting_min_overlap;
      $r{sprintf "gene%s_coverage_softclip", $end} = $count_reads_sc;

      #
      #  optionally export leftover reads with soft-clips:
      #  mine for other events?
      #
      if (0) {
	printf STDERR "leftovers: %d\n", scalar keys %sc_leftover;
	foreach my $id (sort {$sc_len{$b} <=> $sc_len{$a}} keys %sc_leftover) {
	  printf STDERR ">%s soft_clip_length=%d\n", $id, $sc_len{$id};
	  printf STDERR "%s\n", $sc{$id} || die;
	}
      }
    }

    my $coverage = $r{genea_coverage} + $r{geneb_coverage};
    my $supporting = $r{genea_coverage_supporting} + $r{geneb_coverage_supporting};
    $r{supporting_frequency_blended} = sprintf "%.4f", $supporting / $coverage;

    $r{bam} = $f_bam;


    $rpt->end_row(\%r);
    $rpt->finish();
  }

  if ($rpt_passthrough) {
    # also copy results to passthrough report
    my $df = new DelimitedFile("-file" => $f_out,
			       "-headers" => 1,
			      );
    while (my $row = $df->get_hash()) {
      $rpt_passthrough->end_row($row);
    }
  }
}


sub pattern_blast_check {
  my %options = @_;
  my $pid = $options{"-pid"} || die;
  my $pseq = $options{"-pattern"} || die;
  my $lib = $options{"-db"} || die;
  my $end = $options{"-end"} || die;

  my $fl = $BP_SCAN_BLAST_FLANK;
  my $pattern_sample;
  if ($pseq =~ /(\w{$fl})[\[\{\]\}]{2}(\w{$fl})/) {
    # pattern breakpoint doesn't contain interstitial sequence
    $pattern_sample = $1 . $2;
  } else {
    # pattern has intertital sequence.  This might be very long!
    # but the fixed size of the excerpt window "should" help target the
    # true event.
    $pseq =~ /(\w{$fl})[\]\}](.*)/ || die;
    my ($upstream, $down_all) = ($1, $2);
#    die join "\n", $pseq, $upstream, $down_all;
    $down_all =~ s/[\[\{]// || die "can't strip other bracket";
    # remove second bracket before taking sequence excerpt,
    # might fall within sampling window
    my $downstream = substr($down_all, 0, $fl);
    $pattern_sample = $upstream . $downstream;
  }
  die unless $pattern_sample;

  # stripped of breakpoint marker
  die "error" unless length($pattern_sample) == $fl * 2;

  my %hits_match;

  my $tfw = new TemporaryFileWrangler();
  my $f_query = $tfw->get_tempfile("-append" => ".query.fa");
  my $f_db = $tfw->get_tempfile("-append" => ".db.fa");

  # write query (pattern sequence):
  open(FATMP, ">" . $f_query) || die;
  printf FATMP ">%s\n%s\n", $pid, $pattern_sample;
  close FATMP;

  # - BLAST
  # - be VERY conservative about mismatches on right side, what if another
  #   site?
  # - how much overlap/identity across the breakpoint is required?
  #   - how many nt of overlap?
  #   - almost perfect identity?
  #   - maybe: perfect identity +/- 10 nt across breakpoint?

  my $blast = new BLASTer();
  $blast->blast_2_sequences(1);
  $blast->output_format("xml");
  # provides visibility for gapped alignments

  my $chunks = split_list([sort keys %{$lib}], $BP_SCAN_MAX_BLAST_DB_SIZE);

  if ($FLAGS{"export-fasta"}) {
    # export query and soft-clipped reads in FASTA format.
    # ** NOTE ** caution should be used for larger sets as BLAST may
    # cap maximum sequences matched/reported!
    copy($f_query, sprintf "query_excerpt_%s.fa", $pid) || die;
    my $sample = $options{"-sample"} || die "-sample";
    my $bam = $options{"-bam"} || die "-bam";

    my $f_lib = sprintf "lib_%s_%s_%s_%s.fa",
      $pid, $sample, $end, md5_hex($bam);
    # include a MD5 of the BAM path as it's possible there
    # may be more than one BAM file per sample
    open(FATMP, ">" . $f_lib) || die;
    foreach my $id (sort keys %{$lib}) {
      printf FATMP ">%s\n%s\n", $id, $lib->{$id};
    }
    close FATMP;
  }

  foreach my $chunk (@{$chunks}) {
    #
    # write FASTA database (soft-clipped reads) for this chunk of reads:
    #
    unlink $f_db;
    open(FATMP, ">" . $f_db) || die;
    foreach my $id (@{$chunk}) {
      printf FATMP ">%s\n%s\n", $id, $lib->{$id};
    }
    close FATMP;

    my $parser = $blast->blast(
			       "-query" => $f_query,
			       "-database" => $f_db
			      );
    my $result = $parser->next_result;
    # one object per query sequence (only one query seq)

    my $verbose = $FLAGS{verbose};

    if ($result) {
      while (my $hit = $result->next_hit()) {
	# hits from this query to a database sequence
	# (Bio::Search::Hit::HitI)
	my $hit_name = $hit->name();

	my $hsp = $hit->next_hsp();
	# can be more than one HSP per hit, however for this application
	# we're only interested in single-HSP hits

	printf STDERR "    score:%s name:%s strand:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s query_span:%d ref_span:%d total_span=%d query_string=%s hit_string=%s\n",
	  $hsp->score,
	  $hit->name,
	  $hsp->strand,
	  $hsp->range("query"),
	  $hsp->range("hit"),
	  $hsp->num_identical(),
	  $hsp->frac_identical("query"),
	  $hsp->length("query"),
	  $hsp->length("hit"),
	  $hsp->length("total"),
	  $hsp->query_string(),
	  $hsp->hit_string() if $verbose;

	my $hi = new HSPIndexer("-hsp" => $hsp);

	my $q_start = ($BP_SCAN_BLAST_FLANK - $BP_SCAN_MIN_PERFECT_FLANK_TO_MATCH) + 1;
	my $q_end = $BP_SCAN_BLAST_FLANK + $BP_SCAN_MIN_PERFECT_FLANK_TO_MATCH;

	my $identity = $hi->get_query_range_identity_fraction($q_start, $q_end);
	if ($identity >= $BP_SCAN_BLAST_MIN_IDENTITY) {
	  # minimum overlap and identity across breakpoint site is met.
	  #
	  # Next, get the minimum overlap across the breakpoint,
	  # considering both the A and B sides:
	  # move to the left until we run out of mapping
	  my $overlap_a = blast_map_match_count($hi, $BP_SCAN_BLAST_FLANK, -1);
	  # start at end of A side of breakpoint and count backwards
	  # into A
	  my $overlap_b = blast_map_match_count($hi, $BP_SCAN_BLAST_FLANK + 1, 1);
	  # check from start of B forwards
	  my $min_coverage = min($overlap_a, $overlap_b);
	  $hits_match{$hit_name} = $min_coverage;
	  # track the LEAST overlap per HSP, so we know the other
	  # side is at least as good (and maybe better)
	}
      }
    }
  }

  return \%hits_match;
}

sub blast_map_match_count {
  my ($hi, $q_pos, $direction) = @_;
  my $count = 0;
  while (1) {
    my $h_pos = $hi->get_hit_base_for_query_base($q_pos);
#    printf STDERR "query %d matches hit %d\n", $q_pos, $h_pos || -1;
    last unless $h_pos;
    $count++;
    $q_pos += $direction;
  }
  return $count;
}


sub load_pattern_file {
  my ($f_patterns) = @_;
  my $df = new DelimitedFile("-file" => $f_patterns,
			     "-headers" => 1,
			    );
  my %patterns;
  while (my $row = $df->get_hash()) {
    my $pid = $row->{pattern} || die;
    $patterns{$pid} = $row;
  }
  return \%patterns;
}

sub get_bam_scan_reporter {
  my ($pid, $f_out) = @_;
  confess "no output file specified" unless $f_out;

  my @labels = qw(
		   pattern
		   pattern_generation_method
		   sample
		   is_source_sample
		   bam

		   genea_chr
		   genea_pos
		   genea_coverage
		   genea_coverage_softclip
		   genea_coverage_supporting
		   genea_coverage_supporting_shortest_side_overlap

		   geneb_chr
		   geneb_pos
		   geneb_coverage
		   geneb_coverage_softclip
		   geneb_coverage_supporting
		   geneb_coverage_supporting_shortest_side_overlap

		   supporting_frequency_blended
		);

  my $rpt = new Reporter(
			 "-file" => $f_out,
			 "-delimiter" => "\t",
			 "-labels" => \@labels,
			 "-auto_qc" => 1,
			);
  return $rpt;
}

sub bam_to_sample {
  my ($bam) = @_;
  my $sample = basename($bam);
  $sample =~ s/\.bam$//i || die;
  return $sample;
}

sub link_patterns {
  my $f_from = $FLAGS{"patterns-from"} || die;
  my $f_to = $FLAGS{"patterns-to"} || die;

  my $df = new DelimitedFile("-file" => $f_from,
			     "-headers" => 1,
			    );

  my $search_edge_strip = 25;
  # hack: idea is shorter excerpt should be findable in longer
  # sequence even if breakpoint site is a little different due to
  # microhomology

  my %search;
  while (my $row = $df->get_hash()) {
    my $pid = $row->{pattern} || die;
    my $seq = $row->{sequence} || die;

    #    printf STDERR "before: %s\n", $seq;
    $seq = uc($seq);
    $seq = substr($seq, $search_edge_strip);
    $seq = substr($seq, 0, length($seq) - $search_edge_strip);
    $seq =~ s/[\[\]\{\}]//g;
#    printf STDERR "after: %s\n", $seq;
    $search{$pid} = $seq;
  }

  $df = new DelimitedFile("-file" => $f_to,
			  "-headers" => 1,
			 );

  my @pids = sort keys %search;
  while (my $row = $df->get_hash()) {
    my $pid = $row->{pattern} || die;
    my $seq = $row->{sequence} || die;
    $seq = uc($seq);
    $seq =~ s/[\[\]\{\}]//g;

    # possibly optimize by gene?  will be very slow otherwise

    foreach my $pid_search (@pids) {
      my $search = $search{$pid_search};
      my $idx = index($seq, $search);
      if ($idx != -1) {
	printf STDERR "%s hits %s\n", $pid_search, $pid;
      }
    }
  }

  # QC:
  # - do all source patterns match targets?
  # - any orphaned targets?
  # - do any source patterns match > 1 target?

  die "X";
}

sub tag_sjpi {
  # add field indicating for each pattern whether one of the
  # source isoform(s) is a SJ-preferred one.  Helps with reviews,
  # e.g. for PAX5 patterns exon 5 to exon 2 is expected (SJ preferred
  # isoform) but a minor isoform shows exon 4 to exon 2 for the
  # same event.  4-to-2 is not wrong because source data are the same,
  # but filtering to SJ-preferred view can help focus on
  # cases unrelated to this issue.
  my $f_in = $FLAGS{"tag-sjpi"} || die;
  my $df = new DelimitedFile("-file" => $f_in,
			     "-headers" => 1,
			     );
  my $outfile = basename($f_in) . ".sjpi.tab";
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => [
					   qw(
					       genea_acc_preferred
					       geneb_acc_preferred
					    )
					  ],
  			      "-auto_qc" => 1,
			     );

  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $f_sjpi = $config_genome->{GENE_TRANSCRIPT_MATRIX} || die;
  my $sjpi = new SJPreferredIsoform("-file" => $f_sjpi,
				    "-auto_gsm" => 1);


  while (my $row = $df->get_hash()) {

    foreach my $side (qw(a b)) {
      my $key_accs = sprintf 'gene%s_acc', $side;
      my $key_out = sprintf 'gene%s_acc_preferred', $side;
      my @accs = split /,/, $row->{$key_accs};
      die unless @accs;
      my $found_sjpi = 0;
      foreach my $acc (@accs) {
	$found_sjpi = 1 if $sjpi->is_preferred_nm($acc);
      }
      $row->{$key_out} = $found_sjpi;
    }
    $rpt->end_row($row);
  }
  $rpt->finish();
}

sub find_softclip_fusion {
  # 3. write each to a faux-cicero file with minimum annotations
  # 4. run pattern generation pipeline
  # 5. find patterns with expected annotations
  # 6. re-run on source sample: how well does pattern work?
  #    does flanking sequence need to be cleaned, etc.?
  find_binary("samtools", "-die" => 1);

  my $f_bam = $FLAGS{bam} || die "-bam";
  my $rnm = get_bam_rnm($f_bam);
  my $min_sc = $FLAGS{"min-sc-length"} || die "-min-sc-length";

  #
  #  find interval search for soft-clipped reads.
  #  TO DO: other features such as an exon of interest, etc.
  #
  my $chr = $rnm->find_name($FLAGS{chr} || die "-chr");
  my $pos = $FLAGS{pos} || die "-pos";
  my $gene = $FLAGS{gene} || die "-gene";

  #
  #  - query reads in target interval
  #  - find those with a minimum supporting soft clip length
  #  - prepare minimalist fz2 pattern generation input file
  #

  my $sample = basename($f_bam);
  $sample =~ s/\.bam$//i;
  my $f_sc = sprintf "soft_clipped_%s_%s_%d.tab", $sample, $chr, $pos;

  my $need_sc = 1;
  $need_sc = 0 if -s $f_sc;
  $need_sc = 1 if $FLAGS{force};

  if ($need_sc) {
    #
    #  query supporting reads:
    #
    my $cmd = sprintf "samtools view %s %s:%d-%d|", $f_bam, $chr, $pos, $pos;
    #    printf STDERR "cmd: %s\n", $cmd;
    open(SAMTMP, $cmd) || die;

    my $rpt = new Reporter(
			   "-file" => $f_sc,
			   "-delimiter" => "\t",
			   "-labels" => [
					 qw(
					     sample
					     geneA
					     geneB
					     read_id
					     cigar
					     contig
					  )
					],
			   "-auto_qc" => 1,
			  );


    while (<SAMTMP>) {
      chomp;
      my @f = split /\t/, $_;
      my $cigar = $f[SAM_I_CIGAR];
      my @softs = $cigar =~ /(\d+)S/g;
      my $usable_sc =0 ;
      foreach my $soft (@softs) {
	$usable_sc = $soft if $soft >= $min_sc and $soft > $usable_sc;
      }
      if ($usable_sc) {
	my %r;
	$r{sample} = $sample;
	$r{geneA} = $r{geneB} = $gene;
	$r{read_id} = $f[SAM_I_QNAME];
	$r{cigar} = $cigar;
	$r{contig} = $f[SAM_I_SEQ];
	$rpt->end_row(\%r);
      }
    }
    $rpt->finish();
  }

  # stop here, difficult to automate fz2 pattern generation call as
  # various other factors may be in play
}

sub fuzzall_add_sequence {
  my $f_fuzzall = $FLAGS{"fuzzall-add-sequence"} || die;
  my $fa = new Fuzzion2Fuzzall("-file" => $f_fuzzall);
  my $f_pattern = $fa->f_pattern() || die;
  my $hit_row;
  # fuzzall row with desired pattern

  my $f_patterns = $FLAGS{"patterns"} || die;
  my $df = new DelimitedFile("-file" => $f_patterns,
			     "-headers" => 1,
			    );
  # get sequence for patterns:
  my %patterns;
  while (my $row = $df->get_hash()) {
    my $pid = $row->{pattern} || die;
    my $sequence = $row->{sequence} || die;
    die "duplicate $pid" if $patterns{$pid};
    $patterns{$pid} = $sequence;
  }

  $df = new DelimitedFile("-file" => $f_fuzzall,
			     "-headers" => 1,
			    );
  # hack: avoid using Fuzzion2Fuzzall iterator for now as it
  # will die parsing (Excel-)truncated sample detail field
  my $f_out = basename($f_fuzzall) . ".sequence.tab";
  my $rpt = $df->get_reporter(
			      "-file" => $f_out,
			      "-extra" => [
					   qw(
					       sequence
					    )
					  ],
  			      "-auto_qc" => 1,
			     );

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
    my $pid = $row->{$f_pattern} || die;
    $row->{sequence} = $patterns{$pid} || die "can't get sequence for $pid";
    $rpt->end_row($row);
  }
  $rpt->finish();

  # TO DO: additional QC, e.g. check that all PIDs were observed?
}


sub count_patterns {
  my $f_patterns = $FLAGS{"patterns"} || die;

  my $df = new DelimitedFile("-file" => $f_patterns,
			     "-headers" => 1,
			    );
  # TO DO: option to work with multiple files, i.e.  late-stage
  # non-merged pattern files for (A) fusions and (B) ITDs

  my $SOURCE_RK = "RK_MasterFOdb_ForCleanup_100421_SDG";
  # Richard Kriwacki patterns

  #
  #  track counts gene pairs and patterns by source:
  #
  my %all;
  my %rk_exclusive;
  my %sj_exclusive;
  my %shared;
  while (my $row = $df->get_hash()) {
    my $pid = $row->{pattern} || die;
    my $gene_pair = $pid;
    $gene_pair =~ s/\-\d+$// || die;
    my $sources = $row->{source} || die;
    my %sources = map {$_, 1} split /,/, $sources;

    my $has_rk = $sources{$SOURCE_RK} || 0;
    delete $sources{$SOURCE_RK};
    my $has_sj = scalar keys %sources;

    my $shared = ($has_rk and $has_sj);

    my @track;
    push @track, \%all;

    if ($has_rk and !$has_sj) {
      # exclusive to RK set
      push @track, \%rk_exclusive;
    } elsif (!$has_rk and $has_sj) {
      # exclusive to SJ sources
      push @track, \%sj_exclusive;
    } elsif ($shared) {
      # shared between SJ and RK
      push @track, \%shared;
    } else {
      die;
    }

    foreach my $tracker (@track) {
      $tracker->{patterns}{$pid} = 1;
      $tracker->{gene_pairs}{$gene_pair} = 1;
    }
  }

  #
  #  report combinations of counts:
  #
  my $f_out = basename($f_patterns) . ".counts.tab";
  my $rpt = new Reporter(
			 "-file" => $f_out,
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   type
					   count_gene_pairs
					   count_patterns
					)
				      ],
			 "-auto_qc" => 1,
			);

  foreach my $group (
		     [ "all", \%all ],
		     [ "SJ_exclusive", \%sj_exclusive ],
		     [ "RK_exclusive", \%rk_exclusive ],
		     [ "shared", \%shared ],
		     [ "SJ_plus_shared", \%sj_exclusive, \%shared ],
		     [ "RK_plus_shared", \%rk_exclusive, \%shared ],
		    ) {
    my ($label, @trackers) = @{$group};
    my %r;
    $r{type} = $label;
    $r{count_gene_pairs} = build_counts("gene_pairs", @trackers);
    $r{count_patterns} = build_counts("patterns", @trackers);
    $rpt->end_row(\%r);
  }
  $rpt->finish();

}

sub build_counts {
  my ($key, @trackers) = @_;

  my %v;
  foreach my $tracker (@trackers) {
    # create union of values in trackers
    my $set = $tracker->{$key} || die;
    foreach my $v (keys %{$set}) {
      $v{$v} = 1;
    }
  }
  return scalar keys %v;
}

sub pattern2fasta {
  my $f_in = $FLAGS{"pattern2fasta"} || die;
  my $df = new DelimitedFile("-file" => $f_in,
			     "-headers" => 1,
			    );
  my $f_pids = $FLAGS{pids};
  my $f_out = $FLAGS{out};
  unless ($f_out) {
    $f_out = sprintf "%s.fa", basename($f_pids || $f_in);
  }
  my $wf = new WorkingFile($f_out);
  my $fh = $wf->output_filehandle();

  my $wanted;
  if (my $f_pids = $FLAGS{pids}) {
    $wanted = read_simple_file($f_pids, "-hash1" => 1);
  }

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
    my $pid = $row->{pattern} || die;
    my $seq = $row->{sequence} || die;
    $seq =~ tr/{}[]/NNNN/;

    my $usable;
    if (@MATCH_STRINGS) {
      foreach my $s (@MATCH_STRINGS) {
	$usable = 1 if $pid =~ /$s/;
      }
    } elsif ($wanted) {
      $usable = 1 if $wanted->{$pid};
    } else {
      $usable = 1;
    }

    if ($usable) {
      my $id_line = ">" . $pid;
      $id_line .= sprintf " /length=%d", length($seq);

      foreach my $f (@FA_FIELDS) {
	# report additional annotation fields in header line
	die "no field $f" unless exists $row->{$f};
	$id_line .= sprintf " /%s=%s", $f, $row->{$f};
      }
      printf $fh "%s\n%s\n", $id_line, $seq;
    }
  }

  $wf->finish();
}

sub count_by_source {
  # count patterns by sources: TCGA, SJ clingen, ProteinPaint.
  # Since a pattern may be associated with multiple sources,
  # total counts are not expected to to match the raw pattern counts
  my @files = split /,/, ($FLAGS{file} || die);
  # you don't see this

  my %source2category;

  foreach my $src (
		   "pp_svfusion_2021_09_10/pcgp",
		   "pp_svfusion_2021_09_10/target",
		   "pp_svfusion_2021_09_10/Clinical Pilot",
		   "pp_svfusion_2021_09_10/pedccl",
		   "pp_svfusion_2021_09_10/cosmic",
		   "pp_svfusion_2021_09_10/scmc",
		  ) {
    $source2category{$src} = "ProteinPaint";
  }
  foreach my $src (
		   "t3_cicero_good_jneary_20210506",
		   "rrnaseq_clingen_20200924-all",
		   "summary_benchmark184_CiceroV6_forSteve2_filtered",
		   "jenn_spec",

		   # NOTE:
		   "CICERO_only_pattern_JZcomments-final",
		   "CICERO_table_S5_ex_closs",
		   "CCSS_2022_01_05",
		   "bcd-21-0160_supplemental_tables_supp1-32_tab_24",
		   "EGFRvIII_bp_38",
		   # ...not 100% sure about these attributions
		  ) {
    $source2category{$src} = "SJ_clinical";
  }

  $source2category{"RK_MasterFOdb_ForCleanup_100421_SDG"} = "TCGA";

  my %category2pid;
  my %category2pair;

  foreach my $file (@files) {
    my $df = new DelimitedFile(
			       "-file" => $file,
			       "-headers" => 1,
			      );
    while (my $row = $df->get_hash()) {
      my $pid = $row->{pattern} || die;
      my $gene_pair = $pid;
      $gene_pair =~ s/\-\d+$// || die;
      my @sources = split /,/, $row->{source} || die;
      foreach my $source (@sources) {
	my $category = $source2category{$source} || die "no category for source $source";
	$category2pid{$category}{$pid} = 1;
	$category2pair{$category}{$gene_pair} = 1;
      }
    }
  }

  my $rpt = new Reporter(
			 "-file" => "pattern_counts_by_category.tsv",
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   category
					   patterns
					   gene_pairs
					)
				      ],
			 "-auto_qc" => 1,
			);
  foreach my $category (sort keys %category2pid) {
    my %r;
    $r{category} = $category;
    $r{patterns} = scalar keys %{$category2pid{$category}};
    $r{gene_pairs} = scalar keys %{$category2pair{$category}};
    $rpt->end_row(\%r);
  }
  $rpt->finish();
}


sub fuzzall_filter_to_patterns {
  my $f_fuzzall = $FLAGS{"fuzzall-filter-to-pattern-list"} || die;
  my $f_patterns = $FLAGS{pattern} || die "-pattern";

  my $fa = new Fuzzion2Fuzzall("-file" => $f_fuzzall);
  my $f_pattern = $fa->f_pattern() || die;

  # target pattern IDs:
  my $df = new DelimitedFile(
			     "-file" => $f_patterns,
			     "-headers" => 1,
			    );
  my %pid;
  while (my $row = $df->get_hash()) {
    $pid{$row->{pattern} || die} = 1;
  }

  # filter:
  my $f_out = basename($f_fuzzall) . ".pattern_filter.tab";
  my $rpt = $fa->df->get_reporter(
				  "-file" => $f_out,
				  "-auto_qc" => 1,
				 );

  while (my $row = $fa->next()) {
    my $pid = $row->{$f_pattern} || die;
    next unless $pid{$pid};
    $rpt->end_row($row);
  }
  $rpt->finish();
}

sub find_ambiguous_pattern_symbols {
  die "FIX ME: rewrite to use genea/geneb columns for symbols.  will be more reliable than split which can't always cleanly parse by - character";
  die "also be sure to update fz2_group field";
  die "move to fusion_contig_extension.pl to share grouping code";

  my $f_patterns = $FLAGS{patterns} || die "-patterns";
  my $df = new DelimitedFile("-file" => $f_patterns,
			     "-headers" => 1,
			    );

  my $gsm = new_gsm_lite();
  my $genes_approved = $gsm->get_hgnc_approved();
  foreach my $g (@{$genes_approved}) {
    $gsm->add_gene("-gene" => $g);
  }

  my %new2old;
  my %pair_counter;
  while (my $row = $df->get_hash()) {
    my $pid = $row->{pattern} || die;
    my $pair = pid2pair($pid);
    my $pattern_num = pid2pattern_number($pid);
    $pair_counter{$pair} = $pattern_num if $pattern_num > ($pair_counter{$pair} || 0);
    # track highest observed pattern number

    my @genes = split /\-/, $pair;
    if (@genes == 2) {
      die join ",", $pair, @genes unless @genes == 2;
      my ($gene_a, $gene_b) = @genes;

      my @genes_new;
      foreach my $g (@genes) {
	my $new = $gsm->find($g);
	push @genes_new, $new || $g;
      }
      my $pair_new = join "-", @genes_new;

      push @{$new2old{$pair_new}{$pair}}, $pid;
      # bucket each old pair ID by the latest symbol pair
    } else {
      printf STDERR "skipping pattern %s, can't parse to 2 genes\n", $pair;
    }
  }

  my $f_out = "pattern_gene_ambiguity.tsv";
  my $rpt = new Reporter(
			 "-file" => $f_out,
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   pair_new
					   pairs_in
					   patterns_affected
					   patterns_affected_count
					)
				      ],
			 "-auto_qc" => 1,
			);

  my $total_involved_gene_pairs = 0;
  my $total_involved_patterns = 0;
  my $total_affected_patterns = 0;
  my %needs_rename;

  foreach my $pair_new (sort keys %new2old) {
    my @old = sort keys %{$new2old{$pair_new}};
    if (@old > 1) {
      my %r;
      $r{pair_new} = $pair_new;
      $r{pairs_in} = join ", ", map {$_ . "=" . scalar @{$new2old{$pair_new}{$_}}} @old;

      $total_involved_gene_pairs += scalar @old;
      my @affected = grep {$_ ne $pair_new} @old;

      foreach my $pair_affected (@affected) {
	foreach my $pid (@{$new2old{$pair_new}{$pair_affected}}) {
	  $needs_rename{$pid} = $pair_new;
	}
      }

      $r{patterns_affected} = join ", ", map {$_ . "=" . scalar @{$new2old{$pair_new}{$_}}} @affected;

      my $involved_count = sum(map {scalar @{$new2old{$pair_new}{$_}}} @old);
      $total_involved_patterns += $involved_count;

      my $affected_count = sum(map {scalar @{$new2old{$pair_new}{$_}}} @affected);
      $r{patterns_affected_count} = $affected_count;
      $total_affected_patterns += $affected_count;

      $rpt->end_row(\%r);
    }
  }
  $rpt->finish();

  #
  #  summary info:
  #
  printf STDERR "involved gene pairs: %d\n", $total_involved_gene_pairs;
  printf STDERR "total patterns involved:%d affected:%d\n",
    $total_involved_patterns, $total_affected_patterns;

  #
  #  create updated pattern file with renamed pairs:
  #
  $df = new DelimitedFile(
			  "-file" => $f_patterns,
			  "-headers" => 1,
			 );
  my $outfile = basename($f_patterns) . ".ambig_fix.tab";
  $rpt = $df->get_reporter(
			   "-file" => $outfile,
			   "-auto_qc" => 1,
			  );

  my %uniq;
  my %renamed;
  while (my $row = $df->get_hash()) {
    my $pid = $row->{pattern} || die;
    if (my $pair_to = $needs_rename{$pid}) {
      my $pair_from = pid2pair($pid);
      my $counter = ++$pair_counter{$pair_to};
      my $pid_new = sprintf '%s-%02d', $pair_to, $counter;
      $row->{pattern} = $pid_new;
      $renamed{$pid} = $pid_new;
      $pid = $pid_new;
    }
    die if $uniq{$pid};
    $uniq{$pid} = 1;
    # sanity check
    $rpt->end_row($row);
  }
  $rpt->finish();

  #
  #  write log of changes:
  #
  if (%renamed) {
    my $f_out = basename($f_patterns) . ".rename_log.tab";
    my $rpt = new Reporter(
			 "-file" => $f_out,
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   pattern_old
					   pattern_new
					)
				      ],
			 "-auto_qc" => 1,
			  );
    foreach my $from (sort keys %renamed) {
      my %r;
      $r{pattern_old} = $from;
      $r{pattern_new} = $renamed{$from};
      $rpt->end_row(\%r);
    }
    $rpt->finish();
  }


}

sub filter_pattern_set {
  # extract a subset of patterns
  my $f_in = $FLAGS{patterns} || die "-patterns FILE";
  my $f_filter = $FLAGS{"pattern-filter-to-pattern-list"} || die;
  my $wanted;
  my $match;
  if (-f $f_filter) {
    $wanted = read_simple_file($f_filter, "-hash1" => 1);
  } else {
    $match = $f_filter;
  }

  my $df = new DelimitedFile("-file" => $f_in,
			     "-headers" => 1,
			     );
  my $outfile = $FLAGS{out} || basename($f_in) . ".filtered.tab";
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
  			      "-auto_qc" => 1,
			     );

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
    my $pid = $row->{pattern} || die;
    my $usable;
    if ($match) {
      $usable = 1 if $pid =~ /$match/;
    } else {
      $usable = 1 if $wanted->{$pid};
    }
    $rpt->end_row($row) if $usable;
  }

  $rpt->finish();
}

sub pattern_add_interstitial_length {
  my $f_patterns = $FLAGS{"pattern-add-interstitial-length"} || die;
  my $pid2ilen = get_pattern2ilen($f_patterns);
  my $df = new DelimitedFile("-file" => $f_patterns,
			     "-headers" => 1,
			     );
  my $outfile = basename($f_patterns) . ".ilen.tab";
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => [
					   qw(
					       pattern_interstitial_bases
					    )
					  ],
  			      "-auto_qc" => 1,
			     );

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
    my $pid = $row->{pattern} || die;
    $row->{pattern_interstitial_bases} = $pid2ilen->{$pid};
    $rpt->end_row($row);
  }
  $rpt->finish();
}

sub pattern_cap_fusion_interstitial_length {
  my $f_patterns = $FLAGS{"patterns"} || die "-patterns";
  my $f_out = $FLAGS{"out"} || die "-out FILE";
  my $i_max = $FLAGS{"pattern-cap-fusion-interstitial-length"} || die;
  my $pid2ilen = get_pattern2ilen($f_patterns);

  my $df = new DelimitedFile("-file" => $f_patterns,
			     "-headers" => 1,
			     );
  my $rpt = $df->get_reporter(
			      "-file" => $f_out,
			      "-extra" => [
					   qw(
					       pattern_source
					    )
					  ],
  			      "-auto_qc" => 1,
			     );

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
    my $gene_a = $row->{genea_symbol} || die;
    my $gene_b = $row->{geneb_symbol} || die;

    my $usable = 1;
    $usable = 0 if $gene_a eq $gene_b;
    #    $usable = 0 if $row->{fz2_hint} eq "no_condense";
    #
    # experimental Kriwacki patterns; whether or not to include?:
    # - pros:
    #   - original patterns may not work in short-read data
    #   - reversible: we can just drop these when
    #     fz2 2.0 no-interstitial-penalty code is ready
    # - cons:
    #   - including these will add ~2k patterns to the patch (a lot)
    #   - a change to the data used for the published paper
    #     (but again, this patch will be ultimately reversible
    #      with fz2 2.0)
    #
    # decision: process them too

    if ($usable) {
      # only for fusions, ITDs have different interstitial limits
      my $pid = $row->{pattern} || die;
      my $ilen = $pid2ilen->{$pid};
      if ($ilen > $i_max) {
#	printf STDERR "before: %s\n", $row->{sequence};
	my $seq_raw = $row->{sequence};
	my $seq = $seq_raw;
	$seq =~ s/\]// || die;
	# remove old interstitial start
	$seq =~ s/(.{$i_max})\[/\]$1\[/ || die "can't set interstitial site";

	# sanity check:
	my $check1 = $seq_raw;
	my $check2 = $seq;
	foreach ($check1, $check2) {
	  s/[\[\]]//g;
	}
	die unless $check1 eq $check2;

	$row->{sequence} = $seq;
	$row->{pattern_source} = $pid;
	# preserve original pattern ID, as PID in "pattern" field
	# will be updated during the merge process
	$rpt->end_row($row);
      }
    }
  }

  $rpt->finish();
}

sub hit_read_compare {
  my $f_from = $FLAGS{from} || die "-from";
  # source hits, typically a small/debug run to a limited pattern set
  my $f_to = $FLAGS{to} || die "-to";
  # target hits, results for a larger pattern set including the pattern
  # in -from

  my $from_read2pattern = get_read2pattern($f_from);
  my $to_read2pattern = get_read2pattern($f_to);

  my $f_out = $FLAGS{out} || "read_compare.tab";
  my $rpt = new Reporter(
			 "-file" => $f_out,
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   read
					   from_pattern
					   to_pattern
					)
				      ],
			 "-auto_qc" => 1,
			);

  foreach my $read (sort keys %{$from_read2pattern}) {
    my %r;
    $r{read} = $read;
    $r{from_pattern} = $from_read2pattern->{$read} || die;
    $r{to_pattern} = $to_read2pattern->{$read} || "";
    # ???
    $rpt->end_row(\%r);
  }
  $rpt->finish;

}

sub get_read2pattern {
  # for raw .hits files, not tsv
  my ($file) = @_;
  # quick hacky parser
  open(TMPR, $file) || die;
  my $h = <TMPR>;

  my %read2pattern;
  while (1) {
    my $l_pattern = <TMPR>;
    last if $l_pattern =~ /^read\-pairs /;
    # end of file summary
    $l_pattern =~ /^pattern (\S+)/ || die "can't parse pattern ID in $l_pattern";
    my $pid = $1;

    for (my $i = 0; $i <= 1; $i++) {
      my $l_read = <TMPR>;
      die unless $l_read =~ /^read (\S+)/;
      my $read = $1;
      $read2pattern{$read} = $pid;
    }
  }
  close TMPR;
  return \%read2pattern;
}

sub build_summary {
  # standard/basic postprocessing of fz2 hits files
  my $f_patterns = $FLAGS{patterns} || die "-patterns";
  my $f_out = $FLAGS{out} || "summary.tsv";
  my $ignore_log_errors = $FLAGS{"ignore-errors"};
  die "specify -ignore-errors [0|1]" unless defined $ignore_log_errors;
  my $group_by_genes = $FLAGS{"group-by-genes"};
  die "specify -group-by-genes [0|1]" unless defined $group_by_genes;
  my $mod = 50;

  foreach my $b (qw(
		     fz2_wrapper.pl
		     cluster_log_cleanup.pl
		     fuzzum
		  )) {
    find_binary($b, "-die" => 1);
  }

  #
  #  clean up cluster log files for successfully-completed jobs:
  #
  my $cmd = "cluster_log_cleanup.pl";
  $cmd .= " -ignore-errors" if $ignore_log_errors;
  # some failures expected in e.g. GTEx
  system($cmd);

  #
  #  fuzzum:
  #
  my @f_hits = glob("*.hits");
  die "no .hits files found" unless @f_hits;
  my $c = new Counter(\@f_hits, "-mod" => $mod);
  foreach my $f (@f_hits) {
    my $fo = $f . ".fuzzum";
    unless (-s $fo) {
      my $f_tmp = $fo . ".tmp";

      my @opts = sprintf "-id=%s", $f;

      push @opts, "-group=genea_symbol,geneb_symbol" if $group_by_genes;
      # if using this option, when running fuzzall:
      # - the value for the first specified column will appear in the
      #   "fuzzall v1.4.0" column
      # - additional columns will be appended to the results
      # - the counts will be aggregated by the specified columns,
      #   so rather than a row for each pattern ID there will be a single
      #   row representing all patterns with the same gene pair

#      my $cmd = sprintf "cat %s | fuzzum -id=%s > %s",
#	$f, $f, $f_tmp;
      my $cmd = sprintf "cat %s | fuzzum %s > %s",
	$f,
	join(" ", @opts),
	$f_tmp;

      system $cmd;
      die "$cmd failed with $?" if $?;
      rename($f_tmp, $fo) || die;
      # only on successful completion
    }
    $c->next($f);
  }

  #
  #  base "fuzzall" report:
  #
  unless (-s $f_out) {
    my $f_tmp = $f_out . ".tmp";
    my $cmd = sprintf 'fuzzall *.fuzzum > %s', $f_tmp;
    shell_cmd("-cmd" => $cmd);
    rename($f_tmp, $f_out) || die;
  }

  #
  # unique-ify headers to simplify tab-delimited to unique hash keys:
  #
  my $f_out_hfix = $f_out . ".header_fix.tab";
  unless (-s $f_out_hfix) {
    fuzzall_header_fix($f_out);
  }

  #
  # filter "IDs" info to only samples showing minimum strong reads
  #
  my $f_out_strong = $f_out_hfix . ".min_strong.tab";
  unless (-s $f_out_strong) {
    fuzzall_strong_sample($f_out_hfix);
  }

  my $f_out_no_hit = $f_out_strong . ".add_no_hit.tab";
  unless (-s $f_out_no_hit) {
    add_no_hit_patterns($f_out_strong, $f_patterns);
  }

  my $f_out_xlsx = $f_out_no_hit . ".xlsx";
  unless (-s $f_out_xlsx) {
    my $cmd = sprintf 'excel_utils.pl -tab2xlsx %s', $f_out_no_hit;
    shell_cmd("-cmd" => $cmd);
  }

  die $f_out_xlsx;

}


sub fuzzall_header_fix {
  my ($f_fuzzall) = @_;
  unless ($f_fuzzall) {
    $f_fuzzall = $FLAGS{"fuzzall-header-fix"} || die;
  }
  my $f_out = basename($f_fuzzall) . ".header_fix.tab";

  my $wf = new WorkingFile($f_out);
  my $fh = $wf->output_filehandle();

  open(IN, $f_fuzzall) || die "can't open $f_fuzzall: $!";
  my $header = <IN>;
  chomp $header;
  my @h = split /\t/, $header;

  my %rename = map {$_, 1} qw(
			       min
			       median
			       mean
			       max
			    );

  my $key;
  for (my $i = 0; $i < @h; $i++) {
    if ($rename{$h[$i]}) {
      $h[$i] = $key . "_" . $h[$i];
    } else {
      $key = $h[$i];
    }
  }
  printf $fh "%s\n", join "\t", @h;
  # write reformatted header

  # copy rest of file:
  while (<IN>) {
    print $fh $_;
  }

  $wf->finish();
}

sub generate_snpshots {
  my $f_patterns = $FLAGS{patterns} || die "-patterns";
  my $pid = $FLAGS{pattern} || die "-pattern";
  my $patterns = load_pattern_file($f_patterns);
  my $pr = $patterns->{$pid} || die "can't find pattern info for $pid";
  my $manuscript_colors = $FLAGS{"manuscript-colors"};
  die "-manuscript-colors [0|1]" unless defined $manuscript_colors;
  my $bams = read_simple_file($FLAGS{bams} || die "-bams");
  my $tag = $FLAGS{tag} || die "-tag (e.g. DNA/RNA)";

  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $f_fa = $config_genome->{FASTA} || die;

  foreach my $side (qw(a b)) {
    my $base = "gene" . $side;
    my $chr = $pr->{$base . "_chr"} || die;
    my $pos = $pr->{$base . "_pos"} || die;
    if ($pos =~ /,/) {
      # multiple
      my @pos = split /,/, $pos;
      $pos = int(average(@pos));
    }

    foreach my $bam (@{$bams}) {
      my $sample = sj_sample_base(basename($bam));

      my $f_out = sprintf "%s.png",
	join "_",
	$pid, $sample, $tag, $side;

      my $needed = -s $f_out ? 0 : 1;
      $needed = 1 if $FLAGS{force};

      unless ($needed) {
	printf STDERR "skipping %s, exists\n", $f_out;
	next;
      }

      my $cmd = sprintf 'java org.stjude.compbio.snpshot.SNPShot -bam %s -name %s -center %d -fasta %s -out %s',
	$bam,
	$chr,
	$pos,
	$f_fa,
	$f_out;

      if ($manuscript_colors) {
	$cmd .= " -color-reference green";
	$cmd .= " -color-sequence blue";
	$cmd .= " -color-background 250,240,230";
	$cmd .= " -color-border 160,32,240";
	$cmd .= " -color-ruler black";
	$cmd .= " -color-label black";
      }

      system $cmd;
      die "$cmd exited with $?" if $?;
    }
  }
}

sub extract_soft_clipped {
  find_binary("samtools", "-die" => 1);

  my $f_patterns = $FLAGS{patterns} || die "-patterns";
  my $pid = $FLAGS{pattern} || die "-pattern";
  my $format = $FLAGS{format} || die "-format";
  my $patterns = load_pattern_file($f_patterns);
  my $pr = $patterns->{$pid} || die "can't find pattern info for $pid";
  my $bams = read_simple_file($FLAGS{bams} || die "-bams");
  my $bam_query_buffer = 15;
  # flanking query buffer
  my $is_fasta;
  if ($format eq "fasta") {
    $is_fasta = 1;
  } else {
    die "-format only currently supports \"fasta\"";
  }


  foreach my $f_bam (@{$bams}) {
    my $rnm = get_bam_rnm($f_bam);
    my $sample = sj_sample_base(basename($f_bam));

    my %saw;
    my %ids;

    my $out_base = join "_", $pid, $sample;
    my $f_out;
    if ($is_fasta) {
      $f_out = $out_base . ".fa";
    } else {
      die;
    }
    my $wf = new WorkingFile($f_out);
    my $fh = $wf->output_filehandle();

    foreach my $side (qw(a b)) {
      my $base = "gene" . $side;
      my $chr_raw = $pr->{$base . "_chr"} || die;
      my $chr_bam = $rnm->find_name($chr_raw) || die;
      my $pos = $pr->{$base . "_pos"} || die;
      if ($pos =~ /,/) {
	# multiple
	my @pos = split /,/, $pos;
	$pos = int(average(@pos));
      }

      my $cmd = sprintf "samtools view %s %s:%d-%d|",
	$f_bam,	$chr_bam,
	$pos - $bam_query_buffer,
	$pos + $bam_query_buffer;
      # extract a region flanking the target site rather than
      # precisely overlapping it, as I'm not sure how samtools factors
      # soft-clipped regions into its "overlap" calculation
      open(SAMTMP, $cmd) || die;

      while (<SAMTMP>) {
	chomp;
	my @f = split /\t/, $_;
	if (my $max_sc = has_usable_soft_clip(\@f)) {
	  my $flags = $f[SAM_I_FLAGS];
	  die "unhandled: supplementary alignment" if $flags & SAM_FLAG_SUPPLEMENTARY_ALIGNMENT;
	  die "unhandled: secondary alignment" if $flags & SAM_FLAG_SECONDARY_ALIGNMENT;
	  my %info;
	  #	  $info{max_softclip_length} = $max_sc;
	  # might be misleading, a duplicate encountered later might have
	  # a different value
	  if ($flags & SAM_FLAG_SEGMENT_FIRST) {
	    $info{template} = "first";
	  } elsif ($flags & SAM_FLAG_SEGMENT_LAST) {
	    $info{template} = "last";
	  }

	  my @things;
	  foreach my $k (sort keys %info) {
	    push @things, sprintf "/%s=%s", $k, $info{$k};
	  }

	  if ($is_fasta) {
	    my $id = $f[SAM_I_QNAME];
	    $id .= "." . ($flags & SAM_FLAG_RC ? "R" : "F");
	    my $seq = $f[SAM_I_SEQ];

	    unless ($saw{$id}{$seq}) {
	      # the same read may appear mapped to both breakpoints,
	      # ignore if so

	      die "duplicate $id" if $ids{$id};
	      $ids{$id} = 1;
	      # however the identifier should be unique (reflecting strand)
	      # at this point

	      my $rc = reverse_complement($seq);
	      die "ugh, FIX ME" if $saw{$id}{$rc};
	      # same read mapped in a different orientation?

	      printf $fh ">%s %s\n%s\n",
		$id,
		join(" ", @things),
		$seq;
	      $saw{$id}{$seq} = 1;
	    }
	  } else {
	    die;
	  }
	}
      }

    }  # $side

    $wf->finish();
  }  # $bam
}

sub has_usable_soft_clip {
  my ($sam) = @_;
  my $cigar = $sam->[SAM_I_CIGAR];
  my @softs = $cigar =~ /(\d+)S/g;
  my $soft_clipped;
  foreach my $soft (@softs) {
    $soft_clipped = 1 if $soft >= $BP_SCAN_MIN_SOFT_CLIP_LENGTH_TO_CONSIDER;
  }

  if ($soft_clipped) {
    return max(@softs);
  } else {
    return 0;
  }
}

sub generate_bambino_viewer_scripts {
  my $f_patterns = $FLAGS{patterns} || die "-patterns";
  my $pid = $FLAGS{pattern} || die "-pattern";
  my $patterns = load_pattern_file($f_patterns);
  my $pr = $patterns->{$pid} || die "can't find pattern info for $pid";
  my $bams = read_simple_file($FLAGS{bams} || die "-bams");
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $f_fa = $config_genome->{FASTA} || die;

  foreach my $f_bam (@{$bams}) {
    my $sample = sj_sample_base(basename($f_bam));

    foreach my $side (qw(a b)) {
      my $out_base = join "_", "view", $pid, $sample, $side;
      my $f_out = $out_base . ".sh";
      my $wf = new WorkingFile($f_out);
      my $fh = $wf->output_filehandle();

      my $base = "gene" . $side;
      my $chr_raw = $pr->{$base . "_chr"} || die;
      my $pos = $pr->{$base . "_pos"} || die;
      if ($pos =~ /,/) {
	# multiple
	my @pos = split /,/, $pos;
	$pos = int(average(@pos));
      }

      my $cmd = sprintf "java Ace2.AceViewer -no-dbsnp -bam %s -name %s -center %s -fasta %s \$*\n",
	$f_bam,	$chr_raw, $pos, $f_fa;

      printf $fh "%s\n", $cmd;
      $wf->finish();
    }  # $side
  }  # $bam
}

sub source_traceback {
  # given a pattern file, a list of pattern IDs, and a list of
  # source files, find source record for given pattern
  my $f_pattern_list = $FLAGS{"pids"};
  my $f_patterns = $FLAGS{"patterns"} || die "-patterns";
  my $f_sources = $FLAGS{sources} || die "-sources";
  my $f_out = $FLAGS{out} || "traceback_source.tab";
  find_binary("blastn", "-die" => 1);

  my $max_good_results_per_pattern_and_source = 5;
  # stop searching if we get some good quality results for some samples.
  # otherwise for events with lots of samples we could report a lot of
  # redundant results.

  my $patterns = load_pattern_file($f_patterns);
  # TO DO: option for pre-parsed hash

  my $check_patterns;
  if (my $f_pattern_list = $FLAGS{pids}) {
    # check only a specific set of patterns
    $check_patterns = read_simple_file($f_pattern_list);
  } else {
    # check all patterns
    $check_patterns = [ sort keys %{$patterns} ];
  }

  my $sources = read_simple_file($f_sources);
  my %tag2row_sample;
  # index of source tag -> sample -> geneA -> geneB -> rows

  my $min_excerpt_length = 40;

  my $tfw = new TemporaryFileWrangler();
  my $f_query = $tfw->get_tempfile("-append" => ".query.fa");
  my $f_db = $tfw->get_tempfile("-append" => ".db.fa");

  my %tag2rpt;

  foreach my $f_source (@{$sources}) {
    my $tag = source2tag($f_source);

    my $df = new DelimitedFile("-file" => $f_source,
			       "-headers" => 1,
			      );
    my $f_out = basename($f_source) . ".traceback.tab";
    my $rpt = $df->get_reporter(
				"-file" => $f_out,
				"-extra" => [
					   qw(
					       link_pattern
					       link_method
					       link_rows
					    )
					  ],
  			      "-auto_qc" => 1,
			       );
    $tag2rpt{$tag} = $rpt;

    while (my $row = $df->get_hash()) {
      my $sample = $row->{sample} || die "no sample";
      dump_die($row, "unknown format") unless exists $row->{geneA};
      my $geneA = $row->{geneA};
      my $geneB = $row->{geneB};
      # may need some flexibility here
      die unless defined $geneA or defined $geneB;
      # sometimes blank
      if ($geneA and $geneB) {
	my @genes_a = parse_gene_list($geneA);
	my @genes_b = parse_gene_list($geneB);
	# might be a list  :/
	foreach my $ga (@genes_a) {
	  foreach my $gb (@genes_b) {
	    push @{$tag2row_sample{$tag}{$sample}{$ga}{$gb}}, $row;
	  }
	}
      }
    }
  }

  my $c = new Counter($check_patterns);
  foreach my $pid (@{$check_patterns}) {
    $c->next($pid);
    my $pr = $patterns->{$pid} || die;

    my %source_tags;
    foreach my $st (split /,/, $pr->{source} || die) {
      $st =~ s/\/.*//;
      # strip off project suffix, e.g.
      # pp_svfusion_2021_09_10/pcgp => pp_svfusion_2021_09_10
      $source_tags{$st} = 1;
    }
    my @source_tags = sort keys %source_tags;

    my $samples_raw = $pr->{sample};
    unless ($samples_raw) {
      printf STDERR "ERROR: no sample annotation for %s, sources=%s\n", $pid, join ",", @source_tags;
      next;
    }
    my @samples = split /,/, $pr->{sample} || dump_die($pr, "no sample");
    # event might be found in more than one sample

    my $genea_list = $pr->{"genea_symbol"} || dump_die($pr, "no genea_symbol");
    my $geneb_list = $pr->{"geneb_symbol"} || die;

    my @geneA = split /,/, $genea_list;
    my @geneB = split /,/, $geneb_list;
    # standardized pattern file column names for gene symbols
    # due to merging, this may be a list.  The pattern ID will
    # contain only one symbol, however this may also be harmonized
    # to the latest HUGO symbol, which may not match the source data.
    # These symbols should match the source files however.

    my $found_somewhere;
    # in some combination of source(s) and sample(s)

    foreach my $source_tag (@source_tags) {
      # there may be more than one source file per pattern,
      # e.g. same event observed in multiple sets

      unless ($tag2row_sample{$source_tag}) {
	# fix when happens: might be another source outside the list
	printf STDERR "ERROR: no source for $source_tag\n";
	next;
      }

      # if more than one source/sample, don't know which sample
      # appears in this source

      # TO DO: cap search at max # of good results?
      # e.g. BCR-ABL1-44 has a huge number of sample IDs in ProteinPaint
      my $count_good = 0;

      foreach my $sample (@samples) {
	next unless $tag2row_sample{$source_tag}{$sample};
	# might not be in this source

	my @rows_src;
	foreach my $gene_a (@geneA) {
	  foreach my $gene_b (@geneB) {
	    if (my $set = $tag2row_sample{$source_tag}{$sample}{$gene_a}{$gene_b}) {
	      printf STDERR "found %s-%s\n", $gene_a, $gene_b;
	      push @rows_src, @{$set};
	    }
	  }
	}

	unless (@rows_src) {
	  printf STDERR "ERROR: can't find $pid @samples $source_tag $sample $genea_list $geneb_list\n";
	  next;
	  # investigate
	}

	printf STDERR "source rows for %s / %s / %s: %d.\n", $pid, $sample, $source_tag, scalar @rows_src;

	$found_somewhere = 1;
	my $found_type;
	if (@rows_src == 1) {
	  $found_type = "single_record";
	  $count_good++;
	} elsif (@rows_src > 1) {
	  #
	  #  ambiguous source records.
	  #  Resolve by finding best match to interstitial sequence,
	  #  which should perfectly match the source contig.
	  #

	  #
	  #  find excerpt of pattern to align with source contig:
	  #
	  my $seq = $pr->{sequence} || die;
	  my $interstitial;
#	  if ($seq =~ /\](\w+)\[/) {
	  if ($seq =~ /[\]\}](\w+)[\[\{]/) {
	    $interstitial = $1;
#	  } else {
#	    die "can't find interstitial sequence in $seq";
	  }
	  my $junction_excerpt;
	  # excerpt of pattern sequence around junction to use to
	  # align with contig
	  if ($interstitial and
	      length($interstitial) >= $min_excerpt_length) {
	    $junction_excerpt = $interstitial;
	    # might be problematic if sequence is very long
	  } else {
	    $interstitial = "" unless $interstitial;
	    my $ilen = length($interstitial);
	    my $half_needed = int(($min_excerpt_length - $ilen) / 2);
	    # or maybe use ceil()
	    my $up;
	    if ($seq =~ /(\w{$half_needed})[\]\}]/) {
	      # sufficient upstream sequence
	      $up = $1;
	    } elsif ($seq =~ /(\w+)[\]\}]/) {
	      # insufficient, e.g. BLOC1S2-COL5A1-01
	      # GGAAAC][CCAAGAAAGGCTACCAGAAGACGGTTCTGGA
	      $up = $1;
	    } else {
	      die;
	    }

	    my $down;
	    if ($seq =~ /[\[\{](\w{$half_needed})/) {
	      $down= $1;
	    } elsif ($seq =~ /[\[\{](\w+)/) {
	      $down = $1;
	    } else {
	      die;
	    }
	    $junction_excerpt = join "", $up, $interstitial, $down;
#	    die "no/insufficient interstitial in $seq $ilen $half_needed $up $down $junction_excerpt";
	  }
	  die unless $junction_excerpt;

	  my $junction_excerpt_length = length($junction_excerpt);
	  # TO DO: also work with flush or nearly flush patterns,
	  # maybe by sampling +/= bracket boundaries, accounting
	  # for any interstitial length.  Danger though is that
	  # we may extract a sample that's longer than the contig
	  # and so potentially be rejected below.

	  my %id2row;
	  my $id_counter = 0;
	  my %db;
	  foreach my $r (@rows_src) {
	    die "no contig field" unless exists $r->{contig};
	    if (my $contig = $r->{contig}) {
	      # may not always be present
	      my $id = "contig_" . ++$id_counter;
	      $db{$id} = $r->{contig} || dump_die($r, "no contig info!");
	      $id2row{$id} = $r;
	    }
	  }
	  die "no contig info available for $pid" unless %db;

	  write_fasta($f_query, { "query" => $junction_excerpt });
	  write_fasta($f_db, \%db);

	  my $blast = new BLASTer();
	  $blast->blast_2_sequences(1);
	  $blast->output_format("xml");

	  my $parser = $blast->blast(
				     "-query" => $f_query,
				     "-database" => $f_db
				    );
	  my $result = $parser->next_result;
	  # one object per query sequence (only one query seq)

	  my $found_blast;
	  my $found_blast_alt;
	  my @all_q_identity_frac;
	  my $hit_count = 0;
	  if ($result) {
	    while (my $hit = $result->next_hit()) {
	      # hits from this query to a database sequence
	      # (Bio::Search::Hit::HitI)
	      my $hit_name = $hit->name();
	      $hit_count++;

	      my $hsp = $hit->next_hsp();
	      # can be more than one HSP per hit, however for this
	      # application we're only interested in single-HSP hits
	      my $verbose = 0;

	      printf STDERR "    score:%s name:%s strand:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s query_span:%d ref_span:%d total_span=%d query_string=%s hit_string=%s\n",
		$hsp->score,
		$hit->name,
		$hsp->strand,
		$hsp->range("query"),
		$hsp->range("hit"),
		$hsp->num_identical(),
		$hsp->frac_identical("query"),
		$hsp->length("query"),
		$hsp->length("hit"),
		$hsp->length("total"),
		$hsp->query_string(),
		$hsp->hit_string() if $verbose;

	      my $qfi = $hsp->frac_identical("query");

	      push @all_q_identity_frac, $qfi;
#	      printf STDERR "Debug: ni=%d exlen=%d\n", $hsp->num_identical(), $junction_excerpt_length;
	      if ($hsp->num_identical() == $junction_excerpt_length) {
		# perfect match between interstitial and a contig
		$found_blast = $hit_name;
		last;
	      } elsif ($qfi == 1 and $hsp->num_identical() >= ($junction_excerpt_length * 0.95)) {
		# perfect alignment, but fewer bases
		$found_blast_alt = $hit_name;
	      }
	    }
	  } else {
	    die "ERROR: no blast matches??";
	  }

	  if ($found_blast) {
	    my $match_row = $id2row{$found_blast} || die;
	    @rows_src = ( $match_row );
	    $found_type = "blast_perfect";
	    printf STDERR "   ...rescued via blast\n";
	    $count_good++;
	  } elsif ($found_blast_alt) {
	    # HOWEVER: what if a better match found in a different source?
	    my $match_row = $id2row{$found_blast_alt} || die;
	    @rows_src = ( $match_row );
	    printf STDERR "   ...rescued via blast (secondary)\n";
	    $found_type = "blast_shorter";
	  } else {
	    printf STDERR "    ...can't find identical match! (hcount:%d qident: %s)\n", $hit_count, join ",", sort {$b <=> $a} @all_q_identity_frac;
	    # this might happen if very similar patterns are merged
	    # via duplicate checking process
	    # TO DO: consider an error if all sources/samples searched
	    # and still not found?
	  }
	}  # multiple source rows

	#
	#  write source record for this source/sample:
	#
	my $rpt = $tag2rpt{$source_tag};
	foreach my $r (@rows_src) {
	  $r->{link_pattern} = $pid;
	  $r->{link_method} = $found_type || "ambiguous";
	  $r->{link_rows} = scalar @rows_src;
	  $rpt->end_row($r);
	}

	if ($count_good >= $max_good_results_per_pattern_and_source) {
	  printf STDERR "stopping sample search in %s after %d good matches\n", $source_tag, $count_good;
	  last;
	}
      }  # $sample
    }  # $source_tag

    printf STDERR "ERROR: source for %s not found!\n", $pid unless $found_somewhere;
  } # $pid

  foreach my $rpt (values %tag2rpt) {
    $rpt->finish();
  }
}

sub source2tag {
  my ($source) = @_;
  my $tag = basename($source);
  $tag =~ s/\..*//;
  return $tag;
}

sub write_fasta {
  my ($f_out, $db) = @_;
  unlink $f_out;
  my $wf = new WorkingFile($f_out);
  my $fh = $wf->output_filehandle();
  foreach my $id (sort keys %{$db}) {
    printf $fh ">%s\n%s\n", $id, $db->{$id};
  }
  $wf->finish();

}

sub traceback_compare {
  find_binary("blastn", "-die" => 1);
  my $f_check = $FLAGS{"traceback-compare"} || die;
  my $f_patterns_new = $FLAGS{"patterns-new"} || die "-patterns-new";
  # newly generated patterns
  my $f_patterns_orig = $FLAGS{"patterns-orig"} || die "-patterns-orig";
  # current production set
  my $f_fuzzall_orig = $FLAGS{"fuzzall-orig"} || die "-fuzzall-orig";
  my $f_fuzzall_new = $FLAGS{"fuzzall-new"} || die "-fuzzall-new";

  my $f_db = $f_patterns_orig . ".fa";
  die "don't have FASTA version of original patterns file: $f_db" unless -s $f_db;
  my $link_to_check = read_simple_file($f_check);
  my @f_new_patterns = glob("*.extended_500.pattern.tab");
  # traceback files containing link pattern ID (i.e. long
  # interstitial fusion) plus newly-generated pattern IDs for
  # that target
  die "?" unless @f_new_patterns;

  my $fa_orig_db = load_fuzzall($f_fuzzall_orig);
  my $fa_new_db = load_fuzzall($f_fuzzall_new);

  my $patterns_new = load_pattern_file($f_patterns_new);
  my $patterns_orig = load_pattern_file($f_patterns_orig);

  # load patterns db FASTA:
  open(TMP, $f_db) || die;
  my $id;
  my %db;
  while (<TMP>) {
    if (/^>(\S+)/) {
      $id = $1;
    } else {
      chomp;
      $db{$id} .= $_;
    }
  }

  # TO DO: summary info about % translated, etc.

  my $flank_nt = 30;
  my $max_hits_to_search = 3;
  my $verbose = $FLAGS{verbose};
  my $max_mismatches_tolerated = 4 + 2;
  # 4 mismatches are possible just due to N location differences,
  # allow another 2 as well
  my $max_allowed_single_bp_shift = 12;
  # the furthest distance a single breakpoint is allowed to shift
  my $exclude_RK = 1;
  my $source_RK = "RK_MasterFOdb_ForCleanup_100421_SDG";

  my $tfw = new TemporaryFileWrangler();
  my $f_query = $tfw->get_tempfile("-append" => ".query.fa");

  my %link2new;
  foreach my $f (@f_new_patterns) {
    my $df = new DelimitedFile("-file" => $f,
			     "-headers" => 1,
			     );
    while (my $row = $df->get_hash()) {
      my $pid = $row->{pattern} || next;
      # not every input row will have a pattern generated
      my $pid_link = $row->{link_pattern} || die;
      $link2new{$pid_link}{$pid} = 1;
    }
  }

  my $f_out = sprintf "%s.traceback_compare.tab", basename($f_check);

  my $rpt = new Reporter(
			 "-file" => $f_out,
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   pattern
					   pair
					   sample
					   fresh_pattern_matches_db
					   interstitial_length
					   interstitial_seq
					   count_strong_orig
					   count_strong_new_pattern
					)
				      ],
			 "-auto_qc" => 1,
			);

  my $blast = new BLASTer();
  $blast->blast_2_sequences(1);
  $blast->output_format("xml");

  printf STDERR "patterns to check: %d\n", scalar @{$link_to_check};

  if ($exclude_RK) {
    my @passed;
    foreach my $pid_link (@{$link_to_check}) {
      my $orig = $patterns_orig->{$pid_link} || die;
      push @passed, $pid_link unless $orig->{source} eq $source_RK;
      # sources can be a list, but don't exclude if pattern appears
      # in another source as well, just want to get rid of RK-exclusive
    }
    $link_to_check = \@passed;
    printf STDERR "   ...after RK filtering: %d\n", scalar @{$link_to_check};
  }

  if (my $restrict = $FLAGS{restrict}) {
    my @match = grep {/$restrict/} @{$link_to_check};
    die "restrict failed" unless @match;
    $link_to_check = \@match;
  }

  my $c = new Counter($link_to_check);
  foreach my $pid_link (@{$link_to_check}) {
    $c->next($pid_link);
    my @max_bp_shifts;
    my @count_strong_new;
    if ($link2new{$pid_link}) {
      my @set_new = sort keys %{$link2new{$pid_link}};
      foreach my $pid_new (@set_new) {
	my $pseq = $patterns_new->{$pid_new}->{sequence};
	$pseq =~ tr/[]{}/NNNN/;
	my $i_first = index($pseq, "N");
	my $i_last = rindex($pseq, "N");
	my $ilen_new = ($i_last - $i_first) - 1;

	my $excerpt = substr($pseq,
			     $i_first - $flank_nt,
			     ($i_last - $i_first) + ($flank_nt * 2));
	my $excerpt_length = length($excerpt);
	# blast vs. fasta of master set
	# check top X hits
	# require high overlap/identity percent
	die if $i_first == -1 or $i_last == -1;

	my $min_identical_nt = $excerpt_length - $max_mismatches_tolerated;

	write_fasta($f_query, { "query" => $excerpt });
	if (1) {
	  copy($f_query, "query.fa") || die;
	  print STDERR "DEBUG: copied query for $pid_link check to query.fa\n";
	}

	my $parser = $blast->blast(
				   "-query" => $f_query,
				   "-database" => $f_db
			      );
	my $result = $parser->next_result;
	# one object per query sequence (only one query seq)

	my @hits;
	if ($result) {
	  my $hit_count = 0;
	  while (my $hit = $result->next_hit()) {
	    push @hits, $hit;
	    last if @hits >= $max_hits_to_search;
	  }
	}

	foreach my $hit (@hits) {
	  my $hsp = $hit->next_hsp();
	  # just use the first one, we're only interested in excellent matches
	printf STDERR "    score:%s name:%s strand:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s query_span:%d ref_span:%d total_span=%d query_string=%s hit_string=%s\n",
	  $hsp->score,
	  $hit->name,
	  $hsp->strand,
	  $hsp->range("query"),
	  $hsp->range("hit"),
	  $hsp->num_identical(),
	  $hsp->frac_identical("query"),
	  $hsp->length("query"),
	  $hsp->length("hit"),
	  $hsp->length("total"),
	  $hsp->query_string(),
	  $hsp->hit_string() if $verbose;

	  next unless $hsp->num_identical() >= $min_identical_nt;
	  #  ensure HSP is high enough identity

	  #
	  # translate the Ns in the subject sequence to query base numbers:
	  #

	  my $hit_id = $hit->name();
	  my $hit_seq = $db{$hit_id} || die;
	  # get the hit sequence

	  my $hi = new HSPIndexer("-hsp" => $hsp);
	  my $hit_i_start = index($hit_seq, "N");
	  my $hit_i_end = rindex($hit_seq, "N");
	  die if $hit_i_start == -1 or $hit_i_end == -1;
	  my $hit_n_start = $hit_i_start + 1;
	  my $hit_n_end = $hit_i_end + 1;
	  # find base numbers of Ns in hit sequence

	  my $hit_n_start_q = h2q($hi, $hit_n_start);
	  my $hit_n_end_q = h2q($hi, $hit_n_end);

	  if ($hit_n_start_q and $hit_n_end_q) {
	    # require that both breakpoint sites (N) be present
	    # in the aligned portion of the hit sequence
	    my $q_n_start = index($excerpt, "N") + 1;
	    my $q_n_end = rindex($excerpt, "N") + 1;

	    my $bp_diff_start = abs($q_n_start - $hit_n_start_q);
	    my $bp_diff_end = abs($q_n_end - $hit_n_end_q);

	    printf STDERR "summary: %s\n", join ",", $q_n_start, $q_n_end, $hit_n_start_q, $hit_n_end_q, $bp_diff_start, $bp_diff_end;

	    my $max_bp_shift = max($bp_diff_start, $bp_diff_end);
	    printf STDERR "  max bp shift: %d\n", $max_bp_shift;
	    push @max_bp_shifts, $max_bp_shift;
	  }
	}  # $hit

	#
	#  get fz2 matches in new sample:
	#
	my $orig_info = $patterns_orig->{$pid_link} || die;
	if (my $hit_new = $fa_new_db->{$pid_new}) {
	  push @count_strong_new, split /,/, pid2fa_strong($hit_new, $orig_info->{sample});
	  # might be a list
	} else {
	  printf STDERR "WARNING: no hits for %s / %s\n", $orig_info->{sample}, $pid_new;
	  # BAM missing?
	}
      }
    } else {
      printf STDERR "WARNING: no new PIDs for %s\n", $pid_link;
    }

    my $matches_db = 0;
    if (@max_bp_shifts) {
      my $smallest = min(@max_bp_shifts);
      $matches_db = 1 if $smallest <= $max_allowed_single_bp_shift;
    }

    my %r;
    $r{pattern} = $pid_link;
    $r{pair} = pid2pair($pid_link);
    $r{fresh_pattern_matches_db} = $matches_db;

    my $orig = $patterns_orig->{$pid_link} || die;
    $r{sample} = $orig->{sample};
    my $interstitial = get_interstitial($orig->{sequence});
    $r{interstitial_seq} = $interstitial;
    $r{interstitial_length} = length($interstitial);
    # to do: detect coding, etc.?

    $r{count_strong_orig} = pid2fa_strong($fa_orig_db->{pid2pair($pid_link)},
					  $orig->{sample}
					 );
    $r{count_strong_new_pattern} = join ",", sort {$b <=> $a} @count_strong_new;

#    dump_die(\%r);

    $rpt->end_row(\%r);
  }

  $rpt->finish();

  # foreach link pid (or set of interest!)
  # get new pids
  # foreach new pid
  # get pattern sequence
  # blast vs. master db: hit or not?
  # track by pair as well?

}

sub h2q {
  my ($hi, $base) = @_;
  my $max_tries = 4;
  # - gap of 2 possible if 2 consecutive Ns mapped as gap in hit sequence
  # - tolerate a little bit more in case something else I haven't thought of
  my $lookup;
  for (my $try = 0; $try < $max_tries; $try++) {
    $lookup = $hi->get_query_base_for_hit_base($base);
    if ($lookup) {
      last;
    } else {
#      print STDERR "no index for hit base $base\n";
      $base--;
    }
  }
  return $lookup;
}

sub get_interstitial {
  my ($seq) = @_;
  my $interstitial = "";
#  if ($seq =~ /[\]\}](\w+)[\[\{]/) {
  if ($seq =~ /[\]\}](\w*)[\[\{]/) {
    $interstitial = $1;
  } else {
    die "can't find interstitial sequence";
  }
}

sub load_fuzzall {
  my ($f_fuzzall) = @_;
  my $fa = new Fuzzion2Fuzzall("-file" => $f_fuzzall);
  my $f_pattern = $fa->f_pattern() || die;
  my %fa;
  while (my $row = $fa->next()) {
    my $pkey = $row->{$f_pattern} || die;
    die if $fa{$pkey};
    $fa{$pkey} = $row;
  }
  return \%fa;
}

sub pid2fa_strong {
  my ($fa_orig, $sample_list) = @_;

  my %samples = map {$_, 1} split /,/, $sample_list;

  my @counts;
  if ($fa_orig) {
    # might not be any results for this pair in original run
    foreach my $id_ref (@{$fa_orig->{IDs_parsed}}) {
      my $id = $id_ref->{id};
      $id =~ s/\..*//;
      # hack: match SJBALL031931_D1.RNA-Seq.hits to SJBALL031931_D1
      if ($samples{$id}) {
	# this is a wanted sample
	push @counts, $id_ref->{"count_strong_plus"} || 0;
      }
    }
  }
  push @counts, 0 unless @counts;
  return join ",", @counts;
  # re-stringify
}


sub false_negative_check {
  # one-off, 11/14/2024
  my $f_fn = $FLAGS{"false-negative-check"} || die;
  my $f_fa_pair = $FLAGS{"fuzzall-pair"} || die "-fuzzall-pair";
  my $f_patterns = $FLAGS{patterns} || die "-patterns";

  my $patterns = load_pattern_file($f_patterns);
  my %pair2patterns;
  foreach my $pid (keys %{$patterns}) {
    my $pair = pid2pair($pid);
    push @{$pair2patterns{$pair}}, $pid;
  }
  my $fa_pair = load_fuzzall($f_fa_pair);

  my $df = new DelimitedFile("-file" => $f_fn,
			     "-headers" => 1,
			     );
  my $outfile = basename($f_fn) . ".fn_report.tab";
  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => [
					   qw(
					       patterns_for_pair
					       has_xiaotu_source
					       max_interstitial_length
					       fz2_distinct
					       fz2_strong_plus
					    )
					  ],
  			      "-auto_qc" => 1,
			     );
  # do we have a pattern matching this sample and pair?
  # do we have results for SPECIFIC pattern?
  # do we have results for PAIR?  e.g. alternate pattern good?

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
    my $pair = $row->{xiaotu_paper_fusion} || die;
    my $sample = $row->{sample} || die;

    my @patterns_for_pair;
    my $has_xiaotu_source = 0;
    if (my $set = $pair2patterns{$pair}) {
      foreach my $pid (@{$set}) {
	my $p = $patterns->{$pid} || die;
	push @patterns_for_pair, $p;
	my %sources = map {$_, 1} split /,/, $p->{source};
	if ($sources{"Fussion2FalseNegativeReport"}) {
	  $has_xiaotu_source = 1;
	}
      }
    }

    my @ilen = (0);
    foreach my $p (@patterns_for_pair) {
      my $seq = $p->{sequence} || die;
      my $i = get_interstitial($seq);
      push @ilen, length($i);
    }
    $row->{max_interstitial_length} = max(@ilen);

    my $fz2_distinct = "";
    my $fz2_strong_plus = "";

    if (my $pi = $fa_pair->{$pair}) {
      foreach my $idr (@{$pi->{IDs_parsed}}) {
	if ($idr->{id} eq $sample) {
	  $fz2_distinct = $idr->{count_distinct};
	  $fz2_strong_plus = $idr->{count_strong_plus};
	}
      }
    }
    $row->{has_xiaotu_source} = $has_xiaotu_source;
    $row->{fz2_distinct} = $fz2_distinct;
    $row->{fz2_strong_plus} = $fz2_strong_plus;
    $row->{patterns_for_pair} = scalar @patterns_for_pair;
    $rpt->end_row($row);
  }

  $rpt->finish();
}

sub hit2tsv {
  my $f_hit = $FLAGS{"hit2tsv"} || die;
  my $f_out = $FLAGS{out} || basename($f_hit) . ".tsv";

  my @hit_columns = qw(
			read
			alignment
			read_length
			matching_bases
			identity_percent
			strand
			left_overlap
			right_overlap
		     );

  my @headers = ("pattern", @hit_columns);

  my $rpt = new Reporter(
			 "-file" => $f_out,
			 "-delimiter" => "\t",
			 "-labels" => \@headers,
			 "-auto_qc" => 1,
			);


  open(IN, $f_hit) || die;
  my $pattern;
  while (<IN>) {
    if (/^fuzzion2/) {
      # header
    } elsif (/pattern (\S+)/) {
      $pattern = $1;
    } elsif (/read (\S+)/) {
      my $read = $1;
      die unless $pattern;
      chomp;
      my @f = split /\t/, $_;
      die "format?" unless @f == @hit_columns;

      my %r;
      @r{@hit_columns} = @f;
      $r{read} = $read;
      $r{pattern} = $pattern;
      $rpt->end_row(\%r);
    } elsif (/^read\-pairs/) {
      # final line
    } else {
      die "format? $_";
    }
  }
  $rpt->finish();

}

sub add_preferred_summary {
  my $f_in = $FLAGS{"add-preferred-summary"} || die;
  my $f_out = basename($f_in) . ".preferred.tab";
  my $sjpi = get_sjpi();

  my $df = new DelimitedFile("-file" => $f_in,
			     "-headers" => 1,
			     );
  my $rpt = $df->get_reporter(
			      "-file" => $f_out,
			      "-extra" => [
					   qw(
					       gene_pair_summary_preferred
					    )
					  ],
  			      "-auto_qc" => 1,
			     );

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
    my $gene_a = $row->{genea_symbol} || die;
    my $gene_b = $row->{geneb_symbol} || die;
    my $nm_a = $sjpi->get_preferred_isoform($gene_a);
    my $nm_b = $sjpi->get_preferred_isoform($gene_b);
    my $pairs = $row->{gene_pair_summary} || "";
    # might not be present, e.g. RK
    my $filtered;
    foreach my $pair (split /,/, $pairs) {
      my ($side_a, $side_b) = split /\-/, $pair;
      my @a = split /\//, $side_a;
      my @b = split /\//, $side_b;

      my $this_nm_a = $a[1];
      my $this_nm_b = $b[1];
      foreach ($this_nm_a, $this_nm_b) {
	s/\.\d+$//;
	# convert to unversioned
      }
      if ($this_nm_a eq $nm_a and $this_nm_b eq $nm_b) {
#	die join "\n", $pairs, $filtered if $filtered;
	$filtered = $pair;
      }
    }
    $row->{gene_pair_summary_preferred} = $filtered || "";
    $rpt->end_row($row);
  }

  $rpt->finish();

}

sub get_sjpi {
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
  my $f_sjpi = $config_genome->{GENE_TRANSCRIPT_MATRIX} || die;
  my $sjpi = new SJPreferredIsoform("-file" => $f_sjpi,
				    "-auto_gsm" => 1);
  return $sjpi;
}

sub patch_patterns {
  my $f_patterns = $FLAGS{patterns} || die "-patterns";
  my $f_patch = $FLAGS{patch} || die;
  my $f_out = $FLAGS{out} || basename($f_patterns) . ".patched.tab";

  my @f_annot_merge = qw(
			  source
			  sample
			  pathogenicity_somatic
		       );

  my $patches = read_simple_file($f_patch, "-as-hash" => 1);

  my @new_columns = map {$_->{column}} grep {$_->{operation} eq "add_column"} @{$patches};
  my %new_columns = map {$_, 1} @new_columns;

  my $df = new DelimitedFile("-file" => $f_patterns,
			     "-headers" => 1,
			    );

  my $rpt = $df->get_reporter(
			      "-file" => $f_out,
			      "-extra" => \@new_columns,
  			      "-auto_qc" => 1,
			     );

  my $patterns = read_simple_file($f_patterns, "-as-hash" => 1);

  # index by pattern ID:
  my %pid;
  foreach my $r (@{$patterns}) {
    my $pid = $r->{pattern} || die;
    die if $pid{$pid};
    $pid{$pid} = $r;
  }

  my %delete;

  # apply patches:
  foreach my $patch (@{$patches}) {
    my $op = $patch->{operation} || die;
    if ($op eq "add_column") {
      # already handled
    } elsif ($op eq "edit_pattern") {
      my $pid = $patch->{pattern} || die;
      my $f = $patch->{column} || die;
      my $v = $patch->{value};
      die unless defined $v;
      my $r = $pid{$pid} || die;
      die "unknown column" unless exists $r->{$f} or $new_columns{$f};
      $r->{$f} = $v;
    } elsif ($op eq "shift_breakpoints") {
      my @shifts = split /,/, $patch->{value};
      die unless @shifts == 2;
      my $pid = $patch->{pattern} || die;
      my $r = $pid{$pid} || die;
      my $sequence = $r->{sequence} || die;

      foreach my $shift (
			 [ $shifts[0], "]", "}" ],
			 [ $shifts[1], "[", "{" ],
			) {
	my ($amount, @brackets) = @{$shift};

	if ($amount > 0) {
	  # shift breakpoint to the right
	  my $pattern = sprintf '([%s])', join "", map {'\\' . $_} @brackets;
	  $pattern .= sprintf '(%s)', "." x $amount;
	  $sequence =~ s/$pattern/$2$1/ || die;
	} elsif ($amount < 0) {
	  # shift breakpoint to the left
	  my $pattern = sprintf '(%s)', "." x abs($amount);
	  $pattern .= sprintf '([%s])', join "", map {'\\' . $_} @brackets;
	  $sequence =~ s/$pattern/$2$1/ || die;
	}
      }
      $r->{sequence} = $sequence;
      # save changes
    } elsif ($op eq "delete_pattern") {
      my $pid = $patch->{pattern} || die;
      $delete{$pid} = 1;
    } elsif ($op eq "merge_to") {
      # merge one or more patterns into the target pattern,
      # and delete
      my $pid_target = $patch->{pattern} || die;
      my @pid_merge = split /,/, $patch->{value} || die;
      my $r_target = $pid{$pid_target} || die;

      foreach my $f_merge (@f_annot_merge) {
	my %v_target = map {$_, 1} split /,/, $r_target->{$f_merge};
	my %v_target_orig = %v_target;
	# original annotations for target pattern

	foreach my $pid_merge (@pid_merge) {
	  # patterns to merge in
	  die "merge pids identical $pid_merge" if $pid_merge eq $pid_target;
	  my $r_merge = $pid{$pid_merge} || die;
	  foreach my $v (split /,/, $r_merge->{$f_merge}) {
	    $v_target{$v} = 1;
	  }

	  $delete{$pid_merge} = 1;
	}

	printf STDERR "UPDATE %s:\n%s\n%s\n",
	  $f_merge,
	  $r_target->{$f_merge},
	  join ",", sort keys %v_target;

	$r_target->{$f_merge} = join ",", sort keys %v_target;
	# record merged results
      }
    } else {
      die "unhandled operation: $op";
    }
  }

  # if new columns were added, make sure not undef:
  foreach my $r (@{$patterns}) {
    foreach my $f (@new_columns) {
      $r->{$f} = "" unless defined $r->{$f};
    }
  }

  # write final output:
  my %delete_needed = %delete;
  foreach my $r (@{$patterns}) {
    my $pid = $r->{pattern} || die;
    my $usable = 1;
    if ($delete{$pid}) {
      $usable = 0;
      delete $delete_needed{$pid}
    }

    $rpt->end_row($r) if $usable;
  }

  if (%delete_needed) {
    die sprintf "ERROR: can't find requested patterns to delete: %s", join ",", sort keys %delete_needed;
  }

  
  $rpt->finish();
}

sub fuzzum_pair_post {
  # postprocess fuzzall results
  my $files;
  if (my $one = $FLAGS{file}) {
    $files = [ $one ];
  } elsif (my $list = $FLAGS{files}) {
    $files = read_simple_file($list);
  } elsif (my $glob = $FLAGS{glob}) {
    $files = [ glob($glob) ];
    die "$glob matches no files" unless @{$files};
  } else {
    die "specify -file FUZZUM | -files LISTFILE | -glob PATTERN";
  }
  die "no files" unless $files and @{$files};
  my $patterns = load_pattern_file($FLAGS{patterns} || die "-patterns");
  my ($k_test) = (keys %{$patterns});
  die "patterns file missing gte_rank annotation" unless $patterns->{$k_test}{gte_rank};

  my @h_new = qw(
		  exon_key
	       );

  my $format_oncoprint = 1;
  # https://docs.google.com/document/d/1klDZ0MHVkQTW2-lCu_AvpRE4_FcbhdB-yI17wNdPaOM/edit?tab=t.0

  my $f_out_gene_a;
  my $f_out_acc_a;
  my $f_out_chr_a;
  my $f_out_pos_a;
  my $f_out_strand_a;

  my $f_out_gene_b;
  my $f_out_acc_b;
  my $f_out_chr_b;
  my $f_out_pos_b;
  my $f_out_strand_b;

  if ($format_oncoprint) {
    $f_out_gene_a = "gene_a";
    $f_out_acc_a = "isoform_a";
    $f_out_chr_a = "chr_a";
    $f_out_pos_a = "position_a";
    $f_out_strand_a = "strand_a";

    $f_out_gene_b = "gene_b";
    $f_out_acc_b = "isoform_b";
    $f_out_chr_b = "chr_b";
    $f_out_pos_b = "position_b";
    $f_out_strand_b = "strand_b";
  } else {
    die "implement other formats, e.g. cicero";
  }

  push @h_new, (
		$f_out_gene_a,
		$f_out_acc_a,
		$f_out_chr_a,
		$f_out_pos_a,
		$f_out_strand_a,

		$f_out_gene_b,
		$f_out_acc_b,
		$f_out_chr_b,
		$f_out_pos_b,
		$f_out_strand_b,
	       );

  foreach my $f_fuzzum (@{$files}) {
    my $f_hits = $FLAGS{hits};
    if ($f_hits) {
      die "-hits not compatible with -files" if @{$files} > 1;
    } else {
      my $h = $f_fuzzum;
#      if ($h =~ s/_by_gene.txt/_sorted.txt/) {
#	$f_hits = $h if -s $h;
#      }
      if ($h =~ s/\.\w+\.txt$//) {
	# .by_gene.txt, .fuzzum.txt
	$f_hits = $h if -s $h;
      }
    }
    die "no hits file" unless $f_hits;
    die "invalid hits file $f_hits" unless -s $f_hits;

    my %pair2patterns;
    #
    #  parse hits file, track hits by pair and pattern ID:
    #
    my $df = new DelimitedFile("-file" => $f_hits,
			       "-headers" => 1,
			      );
    my $f_id = $df->headers_raw()->[0];
    while (my $row = $df->get_hash()) {
      if ($row->{$f_id} =~ /^pattern (\S+)/) {
	my $pid = $1;
	my $pair = $row->{gene_pair} || die;
	# should be populated
	$pair2patterns{$pair}{$pid}++;
	# track counts of hits by pattern
      }
    }

    my $f_out = $FLAGS{out} || basename($f_fuzzum . ".exons.tsv");

    $df = new DelimitedFile("-file" => $f_fuzzum,
			    "-headers" => 1,
			   );
    my $rpt = $df->get_reporter(
				"-file" => $f_out,
				"-extra" => \@h_new,
				"-auto_qc" => 1,
			       );

    while (my $row = $df->get_hash()) {
      my $pair = $row->{"pattern group"} || die;
      # gene pair

      my $set = $pair2patterns{$pair} || die;
      # pattern->count bucket

      #    my @sorted = sort {$set->{$b} <=> $set->{$a}} keys %{$set};
      my @sorted = sort {$set->{$b} <=> $set->{$a} ||
			   $patterns->{$a}{gte_rank} <=> $patterns->{$b}{gte_rank}
			 } keys %{$set};

      # multi-stage sort: first by number of reads hitting a pattern,
      # and secondarily by preferred transcript ranking.  In other
      # words, if multiple patterns have the same number of supporting
      # reads, for transcript annotation purposes choose the one with
      # the highest SJ preferred isoform ranking.

      # 1. patt
      # TO DO: if multiple patterns share the same count,
      # further tiebreaking by SJPI?
#      printf STDERR "count: %d\n", scalar @sorted;
      if (@sorted >= 5) {
	#      die join ",", map {$_, $set->{$_}, $patterns->{$_}->{gte_rank}, " "} @sorted;
      }

      my $best = $sorted[0];

      my $pattern = $patterns->{$best} || die;

      dump_die($pattern, "no gte_A_exon") unless $pattern->{gte_A_exon};
      dump_die($pattern, "no gte_B_exon") unless $pattern->{gte_B_exon};

      my $key = sprintf '%s.%s-%s',
	$pair,
	($pattern->{gte_A_exon} || die),
	($pattern->{gte_B_exon} || die);
      my $gene_a = $pattern->{gte_A_gene};
      my $acc_a = $pattern->{gte_A_transcript};
      my $gene_b = $pattern->{gte_B_gene};;
      my $acc_b = $pattern->{gte_B_transcript};

      $row->{exon_key} = $key;

      $row->{$f_out_gene_a} = $gene_a;
      $row->{$f_out_acc_a} = $acc_a;
      $row->{$f_out_chr_a} = $pattern->{genea_chr};
      $row->{$f_out_pos_a} = $pattern->{genea_pos};
      $row->{$f_out_strand_a} = $pattern->{genea_strand};

      $row->{$f_out_gene_b} = $gene_b;
      $row->{$f_out_acc_b} = $acc_b;
      $row->{$f_out_chr_b} = $pattern->{geneb_chr};
      $row->{$f_out_pos_b} = $pattern->{geneb_pos};
      $row->{$f_out_strand_b} = $pattern->{geneb_strand};

      $rpt->end_row($row);
    }

    $rpt->finish();
    # iterate through fuzzum file, create report
    # foreach pair, find most commonly-hit pattern
    # find PREFERRED exon annotation for this pattern
  }
}

sub generate_minimal_cicero {
  # generate a minimal/skeletal cicero-format file for use with
  # fz2 pattern generation
  my $f_out = $FLAGS{file} || die "-file OUTFILE";

  my @f_cicero = qw(
		     sample
		     geneA
		     chrA
		     posA
		     ortA
		     geneB
		     chrB
		     posB
		     ortB
		     contig
		     pathogenicity_somatic
		  );

  my $rpt = new Reporter(
			 "-file" => $f_out,
			 "-delimiter" => "\t",
			 "-labels" => \@f_cicero,
			 "-auto_qc" => 1,
			);
  my %r;
  # TO DO: tie() to remember previous values?
  $| = 1;

  if (my $sv6 = $FLAGS{sv6}) {
    $r{sample} = $FLAGS{sample} || die "-sample";
    $r{geneA} = $FLAGS{genea} || die "-genea";
    $r{geneB} = $FLAGS{geneb} || die "-geneb";
    my @f = split /\./, $sv6;
    die unless @f == 6;
    @r{qw(
	   chrA
	   posA
	   ortA
	   chrB
	   posB
	   ortB
	)} = @f;
    $r{contig} = "";
    $r{pathogenicity_somatic} = "";
  } else {
    foreach my $f (@f_cicero) {
      printf "%s: ", $f;
      my $v = <STDIN>;
      chomp $v;
      $v = "" unless defined $v;
      $r{$f} = $v;
    }
  }
  $rpt->end_row(\%r);
  $rpt->finish();
}

sub patch_breakpoints_traceback {
  my @f_traceback = glob("*.traceback.tab");
  die "no .traceback.tab files" unless @f_traceback;
  my $f_patterns = $FLAGS{patterns} || die "-patterns";

  # parse all traceback files
  # bucket rows by source and match type
  #  - just save chrA/B arrayrefs rather than full row contents?
  # foreach pattern row:
  #  - gather all hits by source and match type
  #  - choose the best source/match type


  #
  #  track tracebacks by pid and link method:
  #
  my %pid;
  foreach my $f_traceback (@f_traceback) {
    my $df = new DelimitedFile("-file" => $f_traceback,
			     "-headers" => 1,
			     );
    while (my $row = $df->get_hash()) {
      my $pid = $row->{link_pattern} || die;
      my $link_method = $row->{link_method} || die;

      my $chrA = $row->{chrA} || die;
      my $posA = $row->{posA} || die;
      my $chrB = $row->{chrB} || die;
      my $posB = $row->{posB} || die;

      my $ortA = "";
      my $ortB = "";
      if (exists $row->{ortA}) {
	# CICERO-style source fields
	$ortA = $row->{ortA} || "";
	$ortB = $row->{ortB} || "";
	# value is sometimes not populated e.g. in ProteinPaint db
      }

#      printf STDERR "track %s\n", join " ", $pid, $link_method, $chrA, $posA, $ortA, $chrB, $posB, $ortB;

      push @{$pid{$pid}{$link_method}}, [ $chrA, $posA, $ortA, $chrB, $posB, $ortB ];
    }
  }

  my $df = new DelimitedFile("-file" => $f_patterns,
			     "-headers" => 1,
			     );
  my $outfile = basename($f_patterns) . ".breakpoints.tab";

  my @h_new = qw(
		  genea_chr
		  genea_pos
		  genea_strand
		  geneb_chr
		  geneb_pos
		  geneb_strand
	       );
  push @h_new, "link_method" if $FLAGS{debug};

  my $rpt = $df->get_reporter(
			      "-file" => $outfile,
			      "-extra" => \@h_new,
  			      "-auto_qc" => 1,
			     );
  # use fuzzion2 pattern file column names for these annotations,
  # rather than CICERO/source format column names

  # while (my $row = $df->next("-ref" => 1)) {  # headerless

  my $rank = 0;
  my %method2rank;
  $method2rank{"single_record"} = ++$rank;
  $method2rank{"blast_perfect"} = ++$rank;
  $method2rank{"blast_shorter"} = ++$rank;
  $method2rank{"ambiguous"} = ++$rank;

  my $total = 0;
  my $has_value = 0;
  my $missing = 0;

  while (my $row = $df->get_hash()) {
    my $pid = $row->{pattern} || die;
    my $chrA = "";
    my $chrB = "";
    my $posA = "";
    my $posB = "";
    my $ortA = "";
    my $ortB = "";
    my $link_method = "missing";
    if (my $set = $pid{$pid}) {
      foreach my $method (keys %{$set}) {
	die "no rank for $method" unless $method2rank{$method};
      }
      ($link_method) = sort {$method2rank{$a} <=> $method2rank{$b}} keys %{$set};
      my $records = $set->{$link_method} || die;
      ($chrA, $posA, $ortA, $chrB, $posB, $ortB) = @{$records->[0]};

      if ($link_method eq "ambiguous") {
	my %unique = map {join(".", @{$_}), 1} @{$records};
	if (scalar keys %unique == 1) {
	  $link_method = "rescue_ambiguous_records_but_unique_positions";
	} else {
	  $link_method = sprintf "ambiguous_non_unique_positions=%d", scalar keys %unique;
	}
      }
    }
    $row->{genea_chr} = $chrA;
    $row->{genea_pos} = $posA;
    $row->{genea_strand} = $ortA;
    $row->{geneb_chr} = $chrB;
    $row->{geneb_pos} = $posB;
    $row->{geneb_strand} = $ortB;
    $row->{link_method} = $link_method;
    $rpt->end_row($row);

    $total++;
    if ($chrA and $posA and $chrB and $posB) {
      $has_value++;
    } else {
      $missing++;
    }
  }

  $rpt->finish();

  printf STDERR "total:%d populated:%d blank:%d\n", $total, $has_value, $missing;


}


sub patch_breakpoints_by_features {
  # last-gap breakpoint annotation based on transcript/exon
  my $f_patterns = $FLAGS{patterns} || die "-patterns";
  my $f_refflat = $FLAGS{refflat} || die "-refflat";
  my $f_out = basename($f_patterns) . ".breakpoints.tab";

  my $rff = new RefFlatFile();
  $rff->canonical_references_only(1);
  $rff->strip_sharp_annotations(1);
  $rff->parse_file(
		  "-refflat" => $f_refflat,
		  "-type" => "refgene",
		 );

  my $df = new DelimitedFile("-file" => $f_patterns,
			     "-headers" => 1,
			     );
  my $rpt = $df->get_reporter(
			      "-file" => $f_out,
  			      "-auto_qc" => 1,
			     );

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
    if ($row->{genea_chr} and $row->{genea_pos} and
	$row->{geneb_chr} and $row->{geneb_pos} and
	$row->{genea_strand} and $row->{geneb_strand}
       ) {
      # nothing to do
    } else {
      my $pair = $row->{gte_key} || die;
      my @f = split /,/, $pair, -1;
      # for e.g. ""MAP3K8,NM_005204.4,exon_8,intergenic,,"
      die "GTE field count mismatch: expected 6, value=$pair" unless @f == 6;
      my ($geneA, $accA, $featureA,
	  $geneB, $accB, $featureB) = @f;
      set_rf_breakpoint($row, "A", $rff, $accA, $featureA);
      set_rf_breakpoint($row, "B", $rff, $accB, $featureB);
      # if missing, set genomic breakpoint by accession feature site

      set_rf_strand($row, "A", $rff, $accA);
      set_rf_strand($row, "B", $rff, $accB);
      # if missing, set strand by refflat mapping
    }
    $rpt->end_row($row);
  }

  $rpt->finish();
}

sub set_rf_breakpoint {
  my ($row, $side, $rff, $acc, $feature) = @_;
  #  my $f_chr = "chr" . $side;
  #  my $f_pos = "pos" . $side;
  my $f_chr = sprintf "gene%s_chr", lc($side);
  my $f_pos = sprintf "gene%s_pos", lc($side);

  my $process = 1;
  $process = 0 if $row->{$f_chr} and $row->{$f_pos};
  # to to: option to clobber?

  my $set;
  my $rf;
  my $wanted_exon;
  if ($process and $acc) {
    if ($set = $rff->find_by_accession($acc, "-warn-mismatch" => 1)) {
      printf STDERR "WARNING: multiple refseq mappings for %s\n", $acc
	if @{$set} > 1;
      $rf = $set->[0];
      if ($feature =~ /_(\d+)/) {
	$wanted_exon = $1;
      } else {
	printf STDERR "WARNING: can't find feature number in %s\n", $feature;
      }
    } else {
      printf STDERR "WARNING: can't find refflat entry for %s\n", $acc;
    }
  }

  if ($process and $acc and $set and $wanted_exon) {
    # processable
    my $chr = "";
    my $pos = "";
    my $found;
    my $strand = $rf->{strand} || dump_die($rf, "invalid strand");
    die "invalid strand" unless $strand eq "+" or $strand eq "-";
    $chr = $rf->{chrom} || die;
    my $exon_count = 0;

    my @exons = @{$rf->{exons}};
    @exons = reverse(@exons) if $strand eq "-";
    # in a refflat file, exons are always ordered genomically
    # regardless of transcript strand

    foreach my $exon (@exons) {
      $exon_count++;
      if (0) {
	printf STDERR "debug %s, strand %s exon %d: %s:%d-%d\n", $acc, $strand, $exon_count, $chr, $exon->{start}, $exon->{end};
      }
      if ($exon_count == $wanted_exon) {
	die sprintf "sanity fail for %s: start=%d end=%d", $acc, $exon->{start}, $exon->{end} unless $exon->{start} <= $exon->{end};
	# sanity: assume always in genomic rather than transcript context
	# NM_001017915.3: same base # (?)

	if ($side eq "A") {
	  # upstream: want last base of exon
	  if ($strand eq "+") {
	    $pos = $exon->{end} || die;
	  } elsif ($strand eq "-") {
	    $pos = $exon->{start} || die;
	  } else {
	    die;
	  }
	} elsif ($side eq "B") {
	  # downtream: want first base of exon
	  if ($strand eq "+") {
	    $pos = $exon->{start} || die;
	  } elsif ($strand eq "-") {
	    $pos = $exon->{end} || die;
#	    printf STDERR "CHECK: $acc $strand $wanted_exon $chr $pos\n";
	  } else {
	    die "X";
	  }
	} else {
	  die;
	}
#	dump_die($exon, "debug", 1);
      }
    }

    $row->{$f_chr} = $chr;
    $row->{$f_pos} = $pos;
  }
}

sub set_rf_strand {
  my ($row, $side, $rff, $acc) = @_;
  #  my $f_chr = "chr" . $side;
  #  my $f_pos = "pos" . $side;
  my $f_strand = sprintf "gene%s_strand", lc($side);

  my $process = 1;
  $process = 0 if $row->{$f_strand};
  # to to: option to clobber?

  my $set;
  my $rf;
  if ($process and $acc) {
    if ($set = $rff->find_by_accession($acc, "-warn-mismatch" => 1)) {
      my %strands = map {$_->{strand}, 1} @{$set};
      if (scalar keys %strands > 1) {
	printf STDERR "WARNING: multiple strand mappings for %s\n", $acc;
      } else {
	my ($strand) = keys %strands;
	$row->{$f_strand} = $strand;
	printf STDERR "set strand for %s to %s\n", $acc, $strand;
      }
    } else {
      printf STDERR "WARNING: can't find refflat entry for %s\n", $acc;
    }
  }

}

sub fuzzum_merge {
  my $f_target = $FLAGS{target} || die "-target";
  my $target_sample = $FLAGS{"target-sample"} || die "-target-sample";
  my $target_id = $FLAGS{"target-id"} || die "-target-id";
  my $glob = $FLAGS{glob} || die "-glob";
  my @f_fuzzum = glob($glob);
  die "no files matching $glob" unless @f_fuzzum;
  my $prefix = $FLAGS{"header-prefix"} || die "-header-prefix STRING";

  my %fuzzum;
  #
  #  bucket fuzzum info by sample and tracking ID:
  #
  my ($h_fuzzum_sample, $h_fuzzum_id);
  my @h_add;
  foreach my $f (@f_fuzzum) {
    my $df = new DelimitedFile("-file" => $f,
			     "-headers" => 1,
			      );
    unless ($h_fuzzum_sample) {
      my $h = $df->headers_raw;
      $h_fuzzum_sample = $h->[0];
      $h_fuzzum_id = $h->[$#$h];
      # in fuzzum file, columns for sample and pattern/group can be
      # inferred from header line.
      for (my $i = 1; $i < $#$h; $i++) {
	push @h_add, $h->[$i];
      }
    }
    while (my $row = $df->get_hash()) {
      my $sample = $row->{$h_fuzzum_sample} || die;
      my $id = $row->{$h_fuzzum_id} || die;
      die if $fuzzum{$sample}{$id};
      $fuzzum{$sample}{$id} = $row;
    }
  }

  my %fuzzum_patterns;
  my $report_best_pattern = $FLAGS{"glob-fuzzum-pattern"};
  if ($report_best_pattern) {
    push @h_add, "best_pattern";
    push @h_add, "best_pattern_count";
    my @f_fuzzum_pattern = glob($report_best_pattern);
    foreach my $f (@f_fuzzum_pattern) {
      my $df = new DelimitedFile("-file" => $f,
				 "-headers" => 1,
				);
      while (my $row = $df->get_hash()) {
	my $sample = $row->{$h_fuzzum_sample} || die;
	my $pid = $row->{pattern} || die;
	my $pair = $row->{gene_pair} || die;
	my $strong_plus = $row->{"strong+"};
	die unless defined $strong_plus;
	$fuzzum_patterns{$sample}{$pair}{$pid} += $strong_plus;
      }
    }
  }

  #
  #  merge:
  #

  my $suffix = sprintf ".%s.tsv", $prefix;
  # TO DO:
  # - check for leading period
  # - convert any whitespace to _
#  my $suffix = ".merged.tsv";
  my $f_out = basename($f_target) . $suffix;
  my $df = new DelimitedFile("-file" => $f_target,
			     "-headers" => 1,
			    );
  my @h_prefixed = map {$prefix . "_" . $_} @h_add;
  my $rpt = $df->get_reporter(
			      "-file" => $f_out,
			      "-extra" => \@h_prefixed,
  			      "-auto_qc" => 1,
			     );

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
    my $sample = $row->{$target_sample} || die;
    my $id = $row->{$target_id};
    if (my $fr = $fuzzum{$sample}{$id}) {
      foreach my $h (@h_add) {
	my $h_prefixed = $prefix . "_" . $h;
	$row->{$h_prefixed} = $fr->{$h};
      }

      if ($report_best_pattern) {
	if (my $set = $fuzzum_patterns{$sample}{$id}) {
	  my ($best) = sort {$set->{$b} <=> $set->{$a}} keys %{$set};
	  $row->{$prefix . "_best_pattern"} = $best;
	  $row->{$prefix . "_best_pattern_count"} = $set->{$best};
	} else {
#	  print STDERR "ERROR: no data for $sample $id\n";
	  die "ERROR: no data for $sample $id\n";
	  # should always have this
	}
      }
    } else {
      foreach my $f (@h_prefixed) {
	$row->{$f} = "";
      }
    }
    $rpt->end_row($row);
  }
  $rpt->finish();

}

sub parse_gene_list {
  my ($list) = @_;
  my @genes = split /,/, unquote($list);
  foreach (@genes) {
    my $orig = $_;
    if (s/_loc[A-Z]$//) {
      # refflat "sharp" name modifications
#      die "$orig => $_";
    }
  }
  return @genes;
}

sub find_sj_cohort_bams {
  my $run_dir = "/research_jude/rgs01_jude/groups/zhanggrp/projects/MethodDevelopment/common/fuzzion2/PaperTests/v1.4.0_St_Jude_RNA_patterns_2024-10-11a/STJUDE/";
  # /research_jude/rgs01_jude/groups/zhanggrp/projects/MethodDevelopment/common/fuzzion2/PaperTests/v1.4.0_St_Jude_RNA_patterns_2024-10-11a/STJUDE/*by_gene.txt|sed 's/.*\///'|sed 's/_by_gene.txt//' > ~/steve_sj_sample_list.txt
  my $f_bam_cache = $ENV{HOME} . "/restricted_bam_cache.txt";

  # OTHER QUESTIONS: CSTN (childhood solid tumor network?) vs. other

  #
  #  index available BAMs:
  #
  die unless -s $f_bam_cache;
  my $files = read_simple_file($f_bam_cache);
  my %bams;
  foreach my $bam (@{$files}) {
    next unless $bam =~ /TRANSCRIPTOME/;
    next unless $bam =~ /GRCh38/;
    next if $bam =~ /\/backup\//;
    # exclude "backup" dirs, e.g.
    # /research/rgs01/reference/restricted/ClinicalGenomics/GRCh38/ClinicalPilot/FINAL/TRANSCRIPTOME/v2.0.0/backup/SJMLL001_D.RNA-Seq.bam
    my $sample = basename($bam);
    $sample =~ s/\..*// || die;
    push @{$bams{$sample}}, $bam;
  }

  #
  #  get desired sample list:
  #
  my %samples;
  if (my $f_samples = $FLAGS{samples}) {
    my $set = read_simple_file($f_samples, "-hash1" => 1);
    %samples = %{$set};
  } else {
    die unless -d $run_dir;
    my @files = glob($run_dir . "/*by_gene.txt");
    die unless @files;
    foreach my $file (@files) {
      my $sample = basename($file);
      $sample =~ s/_by_gene\.txt$// || die;
      die "duplicate sample $sample" if $samples{$sample};
      $samples{$sample} = 1;
    }
  }

  #
  #  xref samples, deal with duplicates:
  #
  my $f_cohort_bams = "SJ_cohort_bams.txt";
  my $wf = new WorkingFile($f_cohort_bams);
  my $fh = $wf->output_filehandle();
  foreach my $sample (sort keys %samples) {
    my $set = $bams{$sample} || die "no BAMs for $sample";
    my $final;
    if (@{$set} == 1) {
      $final = $set->[0];
    } else {
      # algorithm here?  larger/newer file?
      # ignore certain projects like ClinicalPilot?
      my @set = sort {-s $b <=> -s $a } @{$set};
      my $s1 = -s $set[0];
      my $s2 = -s $set[1];
      my $frac = $s1 / $s2;
      printf STDERR "largest is %.2f vs second: %s\n", $frac, join " ", @set;
      $final = $set[0];
    }

    die unless $final;
    printf $fh "%s\n", $final;

  }
  $wf->finish();

}

sub clean_hits {
  my $patterns = load_pattern_file($FLAGS{patterns} || die "-patterns");

  my @files;
  if (my $glob = $FLAGS{glob}) {
    @files = glob($glob);
  } elsif (my $lf = $FLAGS{files}) {
    my $files = read_simple_file($lf);
    @files = @{$files};
  } else {
    die "specify -glob GLOB|-files LISTFILE";
  }
  die "no files" unless @files;

  my $f_category = "category";
  # review

  my $only_category = $FLAGS{"only-category"};

  my %bad_categories = map {$_, 1} (
				    "circRNA",
				    "GTExRecurrent",
				    "ReadThrough",

				    "Problematic Patterns",
				    # old label, retired in set11
				   );
  foreach my $f_hit (@files) {
    my %removed;
    my $f_out = basename($f_hit) . ".filtered";
    # to do: maybe keep as .hits for compatibility w/other expectations?
    my $fh = new Fuzzion2Hits(
			      "-hits" => $f_hit,
			      "-out" => $f_out
			     );
    my $count_usable = 0;
    my $count_unusable = 0;
    while (my $info = $fh->next_row()) {
      my $usable = 1;
      my $pid = $fh->get_pid;
      my $category;
      if (my $p = $patterns->{$pid}) {
	dump_die($p, "no $f_category") unless exists $p->{$f_category};
	$category = $p->{$f_category};
	$usable = 0 if $bad_categories{$category};

	if ($only_category) {
	  $usable = $category eq $only_category;
	}
      } else {
	print STDERR "WARNING: can't find pattern for $pid\n";
      }
      if ($usable) {
	$count_usable++;
	$fh->write_row();
      } else {
	$count_unusable++;
	$removed{$category}++;
      }
    }
    $fh->finish;

    printf STDERR "%s: removed %s\n",
      basename($f_hit),
      join(", ", map {$_ . ":" . $removed{$_}} sort keys %removed);
  }
}

sub get_patterns {
  my $f_patterns = $FLAGS{patterns} || die "-patterns";
  return load_pattern_file($f_patterns);
}

sub pattern_self_identity {
  # compare pattern set vs. itself, looking for strong matches to
  # patterns for other gene pairs, e.g. ATXN1-NUTM*
  find_binary("blastn", "-die" => 1);

  printf STDERR "NOTE: results will need to be postprocessed to remove duplicates, e.g. A->B, B->A results\n";

  my $patterns = get_patterns();
  my $tfw = new TemporaryFileWrangler();
  my ($f_query, $f_db);
  my $verbose = $FLAGS{verbose};
  if ($verbose) {
    print STDERR "DEBUG .fa files\n";
    $f_query = "query.fa";
    $f_db = "db.fa";
  } else {
    $f_query = $tfw->get_tempfile("-append" => ".query.fa");
    $f_db = $tfw->get_tempfile("-append" => ".db.fa");
  }

  my $blast = new BLASTer();
  $blast->blast_2_sequences(1);
  $blast->output_format("xml");

  my @pid_sorted = sort keys %{$patterns};

  my $index_start = $FLAGS{"index-start"};
  my $index_end = $FLAGS{"index-end"};
  my $index_out = $FLAGS{"index-out"};

  my $f_out = $index_out ||
    basename($FLAGS{patterns}) . ".self_dupcheck.tsv";
  my $rpt = new Reporter(
			 "-file" => $f_out,
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   pattern
					   gene_pair
					   other_pattern
					   other_gene_pair
					   match_type
					   pattern_overlap
					   blast_overlap
					   blast_identity
					)
				      ],
			 "-auto_qc" => 1,
			);

  my @pid_process;
  if (defined $index_start) {
    # slice, e.g. from mux.pl -index-mode
    $rpt->write_headers();
    # force header write for mux rejoin, otherwise size 0 files
    # for sets w/no data
    @pid_process = @pid_sorted[$index_start .. $index_end];
  } else {
    # slow
    @pid_process = @pid_sorted;
  }

  my $c = new Counter(\@pid_process);
  foreach my $pid_this (@pid_process) {
    $c->next($pid_this);
    my $p_this = $patterns->{$pid_this} || die;
    my $pair_this = $p_this->{gene_pair} || die;
    my $seq_this = mask_pattern($p_this->{sequence});
    my %db;
    foreach my $pid_other (@pid_sorted) {
      # build set of all patterns for OTHER gene pairs
      my $p_other = $patterns->{$pid_other} || die;
      unless ($p_other->{gene_pair} eq $pair_this) {
	my $seq_other = mask_pattern($patterns->{$pid_other}{sequence});
	$db{$pid_other} = $seq_other;
      }
    }
    die "ERROR: no other patterns found??" unless %db;
    # shouldn't happen
    unlink($f_query, $f_db);
    write_fasta($f_query, { $pid_this => $seq_this });
    write_fasta($f_db, \%db);

    my $parser = $blast->blast(
			       "-query" => $f_query,
			       "-database" => $f_db
			      );
    my $result = $parser->next_result;
    # one object per query sequence (only one query seq)

    if ($result) {
      while (my $hit = $result->next_hit()) {
	my $hit_name = $hit->name();
	my $p_other = $patterns->{$hit_name} || die;

	my $hsp = $hit->next_hsp();
	# can be more than one HSP per hit, however for this application
	# we're only interested in single-HSP hits

	printf STDERR "    score:%s name:%s strand:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s frac_identical_hit:%s query_span:%d ref_span:%d total_span=%d query_string=%s hit_string=%s\n",
	  $hsp->score,
	  $hit->name,
	  $hsp->strand,
	  $hsp->range("query"),
	  $hsp->range("hit"),
	  $hsp->num_identical(),
	  $hsp->frac_identical("query"),
	  $hsp->frac_identical("hit"),
	  $hsp->length("query"),
	  $hsp->length("hit"),
	  $hsp->length("total"),
	  $hsp->query_string(),
	  $hsp->hit_string() if $verbose;

	next if $hsp->frac_identical("query") < $SELF_CHECK_MIN_IDENTITY;
	next if $hsp->frac_identical("hit") < $SELF_CHECK_MIN_IDENTITY;
	# minimum HSP identity checks

	my $ok_flank_up = 0;
	my $ok_flank_down = 0;

	breakpoint_flank_check($hsp->query_string(), $p_this, \$ok_flank_up, \$ok_flank_down);
	breakpoint_flank_check($hsp->hit_string(), $p_other, \$ok_flank_up, \$ok_flank_down);
	next unless $ok_flank_up and $ok_flank_down;
	# any reason this might still be problatic??

	my $found_full_length;
	my $found_overlap;

	#
	#  check for high identity match that covers almost all
	#  of one of the patterns:
	#
	my $num_identical = $hsp->num_identical();
	my $frac_q = $num_identical / length($seq_this);
	my $frac_s = $num_identical / length($patterns->{$hit_name}{sequence});
	my $frac_best = max($frac_q, $frac_s);
	# pattern lengths vary, so use the highest identity fraction
	# vs. either sequence
	if ($frac_best >= $SELF_CHECK_MIN_OVERLAP_FRAC_HQ) {
	  # high overlap for the entire pattern
	  $found_full_length = 1;
	} elsif ($frac_best >= $SELF_CHECK_MIN_OVERLAP_FRAC_LQ) {
	  # overlap over less of the pattern length, howeever
	  # this may still be significant, especially for fusion
	  # clutter cases
	  $found_overlap = 1;
	}

	if ($found_full_length or $found_overlap) {
	  my $this_gene_pair = $p_this->{gene_pair} || die;
	  my $other_gene_pair = $p_other->{gene_pair} || die;

	  my %r;
	  $r{pattern} = $pid_this;
	  $r{gene_pair} = $this_gene_pair;
	  $r{other_pattern} = $hit_name;

	  $r{other_gene_pair} = $other_gene_pair;

	  $r{pattern_overlap} = sprintf "%.3f", $frac_best;
	  $r{blast_overlap} = $hsp->length("total");
	  $r{blast_identity} = sprintf "%.3f", $hsp->frac_identical("query");
	  # simplified identity just for the blast overlap

	  if ($found_full_length) {
	    # high identity over entire pattern length
	    $r{match_type} = "full_length";
	  } elsif ($found_overlap) {
	    # high identity over a significant sub-region
	    my $match_type = "overlap";

	    if (index($this_gene_pair, $other_gene_pair) != -1 or
		index($other_gene_pair, $this_gene_pair) != -1) {
	      # e.g. ASB3-PHF6-01 vs. GPR75-ASB3-PHF6-01
	      $match_type = "fusion_clutter";
	    }
	    $r{match_type} = $match_type;
	  } else {
	    die "unimplemented";
	  }

	  $rpt->end_row(\%r);
	}

      }
    }
  }

  $rpt->finish();

}

sub mask_pattern {
  my ($seq) = @_;
  $seq =~ tr/}{][/RYRY/;
  return $seq;
}

sub breakpoint_flank_check {
  my ($sequence, $p, $ref_ok_up, $ref_ok_down) = @_;
  my $i1 = index($sequence, "R");
  my $i2 = index($sequence, "Y");
  my $ok = 0;
  if ($i1 != -1 and $i2 != -1) {
    my $b1 = $i1 + 1;
    my $b2 = $i2 + 1;
    my $flank_up = $b1;
    my $flank_down = length($sequence) - $b2;

    $$ref_ok_up = 1 if $flank_up >= $SELF_CHECK_MIN_FLANK;
    $$ref_ok_down = 1 if $flank_down >= $SELF_CHECK_MIN_FLANK;

    # rescue: consider a side OK if overlap covers most of that
    # side of the pattern, e.g. ANKHD1-MIR3655-01, B side only 83 nt
    my $seq = $p->{sequence};
    my $pi1 = index($seq, "]");
    $pi1 = index($seq, "}") if $pi1 == -1;
    die if $pi1 == -1;

    my $pi2 = index($seq, "[");
    $pi2 = index($seq, "{") if $pi2 == -1;
    die if $pi2 == -1;

    my $pb1 = $pi1 + 1;
    my $pb2 = $pi2 + 1;

    my $p_flank_up = $pb1;
    my $p_flank_down = length($seq) - $pb2;

    if ($flank_up / $p_flank_up >= $SELF_CHECK_MIN_OVERLAP_FRAC_HQ) {
      $$ref_ok_up = 1;
    }

    if ($flank_down / $p_flank_down >= $SELF_CHECK_MIN_OVERLAP_FRAC_HQ) {
      $$ref_ok_down = 1;
    }
  }
  return $ok;
}

sub pattern_self_identity_post {
  my $f_in = $FLAGS{"pattern-self-identity-post"} || die;

  my $df = new DelimitedFile("-file" => $f_in,
			     "-headers" => 1,
			     );
  my $f_out = basename($f_in) . ".post.tab";
  my $rpt = $df->get_reporter(
			      "-file" => $f_out,
  			      "-auto_qc" => 1,
			     );

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  my %saw;
  while (my $row = $df->get_hash()) {
    my $key = join "/", sort @{$row}{qw(pattern other_pattern)};
    $rpt->end_row($row) unless $saw{$key};
    $saw{$key} = 1;
  }

  $rpt->finish();
}

sub clean_pattern_categories {
  my $f_patterns = $FLAGS{"patterns"} || die "-patterns";

  my %categories_wanted = (
			   "GermlineSV" => 0,
			   "circRNA" => 0,
			   "ITD" => 1,
			   "Admix" => 1,
			   "NoMatch" => 1,
			   "RareGermline" => 0,
			   "RareLowConfidence" => 0,
			   "ReadThrough" => 0,
			   "GTExRecurrent" => 0,
			   "RareSomatic" => 1,
			   "NoReadPairDetected" => 1,
			   "LowConfidenceEx" => 0,
			   "ReadThroughEx" => 0,
			   "Intragenic deletion" => 1,
			   "BelowThreshold" => 1,
			   "RecurrentFusion" => 1,
			   "LowConfidence" => 0,
			  );
  # per JZ email 10/7/2025

  my $df = new DelimitedFile("-file" => $f_patterns,
			     "-headers" => 1,
			     );
  my $f_out = basename($f_patterns) . ".category_clean.tsv";
  my $rpt = $df->get_reporter(
			      "-file" => $f_out,
  			      "-auto_qc" => 1,
			     );

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
  while (my $row = $df->get_hash()) {
    die unless defined $row->{category};
    my $category = $row->{category};
    my $usable;
    if ($category) {
      $usable = $categories_wanted{$category};
      die "ERROR: no info for category $category" unless defined $usable;
      die unless $usable == 0 or $usable == 1;
    } else {
      # new patterns won't have a category yet
      $usable = 1;
    }
    die unless defined $usable;
    $rpt->end_row($row) if $usable;
  }

  $rpt->finish();

}

