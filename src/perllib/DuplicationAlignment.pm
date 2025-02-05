package DuplicationAlignment;
# tool for handling a duplication mapped within two ordered HSPs
# MNE 12/2023

use strict;
use Carp qw(confess);

use Exporter;

use Set::IntSpan;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);

@DuplicationAlignment::ISA = qw(Configurable Exporter);
@DuplicationAlignment::EXPORT_OK = qw();

use MethodMaker qw(
	hi_left
        hi_right

has_duplication

subject_sequence
bracket_copy_count

subject_upstream_base
subject_downstream_base
query_bracket_start_base
query_bracket_end_base

dup_sequence
		  );
# hi*: HSPIndexer.pm, with query overlaps already resolved

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  $self->setup();
  return $self;
}

sub setup {
  my ($self, %options) = @_;
  my $hi_left = $self->hi_left() || die;
  my $hi_right = $self->hi_right() || die;
  my $subject_sequence = $self->subject_sequence() || die;
  my $bcc = $self->bracket_copy_count() || die;

  my $si_s_left = $hi_left->range("subject", 1);
  my $si_s_right = $hi_right->range("subject", 1);
  my $intersect = intersect $si_s_left $si_s_right;
  my $has_duplication = 0;
  if ($intersect->size()) {
    # subject sequence duplication detected

    # TO DO: maybe some minimum?
    $has_duplication = 1;

    my ($so_start, $so_end) = si2range($intersect);
    # interval in subject sequence duplicated in both HSPs
    printf STDERR "subject overlap: %d-%d\n", $so_start, $so_end;

    my $dup_q_left_start = $hi_left->get_query_base_for_hit_base($so_start);
    my $dup_q_left_end = $hi_left->get_query_base_for_hit_base($so_end);
    # interval of overlap in query sequence for left HSP

    my $dup_q_right_start = $hi_right->get_query_base_for_hit_base($so_start);
    my $dup_q_right_end = $hi_right->get_query_base_for_hit_base($so_end);
    # interval of overlap in query sequence for right HSP

    my $dup_sequence = substr($subject_sequence,
			      $so_start - 1,
			      ($so_end - $so_start) + 1);
    # duplicated portion of the subject sequence

    # TO DO:
    # - interstitial sequence (defined as the query sequence between copies,
    #   if any, *not* the sequence ultimately bracketed in the fz2 pattern)

    my $subject_upstream_base;
    # this base number and upstream are the reference sequence
    # (before the bracketed region)

    my $subject_downstream_base;
    # this base number and downstream are the reference sequence
    # (after the bracketed region)

    my $query_bracket_start;
    my $query_bracket_end;
    # bracketed sequence base numbers in query sequence

    $subject_downstream_base = $so_end + 1;
    # first base after end of duplication

    if ($bcc == 1) {
      #
      #  target only a single copy of the duplication.
      #  Vulnerable to false positives.
      #
      $subject_upstream_base = $so_end;
      # in this mode the first copy is considered part of the reference.
      # so, left upstream is the last base number of the first copy.

      $query_bracket_start = $hi_left->get_query_base_for_hit_base($so_end) + 1;
      # base in query immediately after the end of the first copy of dup.
      # So, this will include any interstitial sequence between copies,
      # then the second copy.
      $query_bracket_end = $hi_right->get_query_base_for_hit_base($so_end);
      # last base of the second copy in query sequence.
    } elsif ($bcc == 2) {
      #
      #  include both copies of the duplication in the bracketed region.
      #
      $subject_upstream_base = $so_start - 1;
      # first base before start of dup

      $query_bracket_start = $hi_left->get_query_base_for_hit_base($so_start);
      # first base of first copy
      $query_bracket_end = $hi_right->get_query_base_for_hit_base($so_end);
      # last base of second copy
    } else {
      die "unhandled bracketed copy count $bcc";
    }

    $self->dup_sequence($dup_sequence);
    $self->subject_upstream_base($subject_upstream_base);
    $self->subject_downstream_base($subject_downstream_base);
    $self->query_bracket_start_base($query_bracket_start);
    $self->query_bracket_end_base($query_bracket_end);
  }

  $self->has_duplication($has_duplication);

#  die "X";
}

sub blast_rescue {
    # EXPERIMENTAL: FALSE POSITIVE RISK?
    # however can help with false negatives like
    # /research/rgs01/home/clusterHome/medmonso/work/steve/2020_06_13_cicero_contig_extension/2021_01_04_all/blast/ITD_refactor/final_set/v2_fixed_annot/FLT3-FLT3-02/src.tab
    my ($self) = @_;

    my $hi_left = $self->hi_left() || die;
    my $hi_right = $self->hi_right() || die;

    my $subject_upstream_base = $hi_left->end("subject");
    # this base number and upstream are the reference sequence,
    # before the bracketed region starts.
    # i.e. right edge of the upstream contig match

    my $query_bracket_start = $hi_left->end("query") + 1;
    # start base number of bracketed sequence in query

    my $subject_downstream_base = $hi_right->start("subject");
    # this base number and downstream are the reference sequence,
    # after the bracketed region ends.
    # i.e. left edge of downstream contig match

    my $query_bracket_end = $hi_right->start("query") - 1;
    # end base number of bracketed sequence in query

    # OK if query bracket end < start, will be treated
    # zero interstitial bases later

    $self->subject_upstream_base($subject_upstream_base);
    $self->query_bracket_start_base($query_bracket_start);
    $self->subject_downstream_base($subject_downstream_base);
    $self->query_bracket_end_base($query_bracket_end);
}


sub si2range {
  # STATIC
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

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
