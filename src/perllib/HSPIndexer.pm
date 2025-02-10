package HSPIndexer;
# for BLAST HSPs, translates between query and subject base numbers
# and string indexes (i.e. compensating for alignment gaps)
# MNE 9/2023

use strict;
use Exporter;
use Carp q(confess);

use Set::IntSpan;
use List::Util qw(min max);

use Configurable;
use MiscUtils qw(get_hash_option dump_die);

@HSPIndexer::ISA = qw(Configurable Exporter);
@HSPIndexer::EXPORT_OK = qw();

use constant GAP_CHARACTER => "-";

use MethodMaker qw(
		    hsp
		    query_base_to_hit_base
		    hit_base_to_query_base
		    hit_base_to_query_index
		    query_base_to_sequence
		    hit_base_to_sequence
		    hit_gaps
		    min_gap_length
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  $self->setup();
  return $self;
}

sub setup {
  my ($self) = @_;
  my $hsp = $self->hsp || die "-hsp";

  my $q_aln = $hsp->query_string();
  my $h_aln = $hsp->hit_string();

  my $q_pos = $hsp->start("query");
  my $h_pos = $hsp->start("hit");

  my $ai = 0;
  my $alen = length $q_aln;

  my %q2h;
  my %h2q;
  my %h2qi;

  my %q2nt;
  # query base number to aligned nucleotide
  my %h2nt;
  # hit base number to aligned nucleotide

  my @hit_gaps_raw;
  my $hit_gap_ref;
  my ($last_q_pos, $last_h_pos);

  for (my $ai = 0; $ai < $alen; $ai++) {
    my $q_nt = substr($q_aln, $ai, 1);
    my $h_nt = substr($h_aln, $ai, 1);

    my $q_gap = $q_nt eq GAP_CHARACTER;
    my $h_gap = $h_nt eq GAP_CHARACTER;

    unless ($q_gap or $h_gap) {
      # if a gap is present in either the query or the subject,
      # don't record base position translations
      $q2h{$q_pos} = $h_pos;

      $h2q{$h_pos} = $q_pos;
      $h2qi{$h_pos} = $ai;
      $h2nt{$h_pos} = $h_nt;

      $last_q_pos = $q_pos;
      $last_h_pos = $h_pos;
      # for gap tracking
    }

    #
    # gap tracking:
    #
    if ($h_gap) {
      # gap in hit/subject sequence
      unless ($hit_gap_ref) {
	# create new object
	die "query is also a gap" if $q_gap;
	# happens?
	$hit_gap_ref = {};
	$hit_gap_ref->{query_start} = $q_pos;
	$hit_gap_ref->{upstream_query_base} = $last_q_pos;
	$hit_gap_ref->{upstream_hit_base} = $last_h_pos;
	push @hit_gaps_raw, $hit_gap_ref;
      }
      $hit_gap_ref->{query_end} = $q_pos;
      # keep updating as gap extends
    } else {
      # no hit gap at this position
      if ($hit_gap_ref) {
	$hit_gap_ref->{downstream_query_base} = $q_pos;
	$hit_gap_ref->{downstream_hit_base} = $h_pos;
      }
      $hit_gap_ref = undef;
    }


    $q2nt{$q_pos} = $q_nt unless $q_gap;

    $q_pos++ unless $q_gap;
    $h_pos++ unless $h_gap;
  }

  my @hit_gaps;
  # final/QC passed
  my $min_gap_len = $self->min_gap_length();
  foreach my $g (@hit_gaps_raw) {
    my $gap_len = ($g->{query_end} - $g->{query_start}) + 1;
    next if $min_gap_len and $gap_len < $min_gap_len;
    $g->{gap_length} = $gap_len;
    push @hit_gaps, $g;
  }

  $self->query_base_to_hit_base(\%q2h);
  $self->hit_base_to_query_base(\%h2q);
  $self->hit_base_to_query_index(\%h2qi);
  $self->query_base_to_sequence(\%q2nt);
  $self->hit_base_to_sequence(\%h2nt);
  $self->hit_gaps(\@hit_gaps);
}

sub get_query_index_for_hit_base {
  my ($self, $hit_base) = @_;
  return $self->hit_base_to_query_index()->{$hit_base};
}

sub get_query_base_for_hit_base {
  my ($self, $hit_base) = @_;
  return $self->hit_base_to_query_base()->{$hit_base};
}

sub get_hit_base_for_query_base {
  my ($self, $query_base) = @_;
  return $self->query_base_to_hit_base()->{$query_base};
}

sub range {
  # implement similarly to HSP.
  # SLOW: PRECALCULATE?
  my ($self, $type, $return_si) = @_;
  my $hash;
  if ($type eq "query") {
    $hash = $self->query_base_to_hit_base();
  } elsif ($type eq "subject") {
    $hash = $self->hit_base_to_query_base();
  } else {
    confess "type must be either query or subject";
  }
  my @k = keys %{$hash};
  my $min = min(@k);
  my $max = max(@k);

  if ($return_si) {
    return new Set::IntSpan(sprintf "%d-%d", $min, $max);
  } else {
    return ($min, $max);
  }
}

sub end {
  # implement similarly to HSP
  # SLOW: PRECALCULATE?
  my ($self, $type) = @_;
  my $hash;
  if ($type eq "query") {
    $hash = $self->query_base_to_hit_base();
  } elsif ($type eq "subject") {
    $hash = $self->hit_base_to_query_base();
  } else {
    confess "type must be either query or subject";
  }
  my @k = keys %{$hash};
  return max @k;
}

sub start {
  # implement similarly to HSP
  # SLOW: PRECALCULATE?
  my ($self, $type) = @_;
  my $hash;
  if ($type eq "query") {
    $hash = $self->query_base_to_hit_base();
  } elsif ($type eq "subject") {
    $hash = $self->hit_base_to_query_base();
  } else {
    confess "type must be either query or subject";
  }
  my @k = keys %{$hash};
  return min @k;
}

sub get_query_range_identity_fraction {
  # identity fraction between query and hit, for the specified query range.
  # aligned bases only (i.e. match/mismatch), doesn't address gaps.
  my ($self, $query_start, $query_end) = @_;

  my @query_bases = ($query_start .. $query_end);

  my $qb2s = $self->query_base_to_sequence();
  my $hb2s = $self->hit_base_to_sequence();
  my $qb2hb = $self->query_base_to_hit_base();

  my $count_match = 0;
  my $count_mismatch = 0;
  for (my $qb = $query_start; $qb <= $query_end; $qb++) {
#    my $qs = $qb2s->{$qb} || die "no sequence for query base $qb";
    my $qs = $qb2s->{$qb};
    if ($qs) {
      if (my $hit_base = $qb2hb->{$qb}) {
	my $hs = $hb2s->{$hit_base} || die;
	if ($qs eq $hs) {
	  $count_match++;
	} else {
	  $count_mismatch++;
	}
      } else {
	# subject has a gap for this query position
	$count_mismatch++;
      }
    } else {
      # query base is outside of result, consider mismatch
      $count_mismatch++;
    }
  }
  return $count_match / ($count_match + $count_mismatch);
}

sub get_query_range_identity_score {
  # evaluate quality of alignment over a given query interval.
  # consider alignment matches and gaps.
  my ($self, $query_start, $query_end) = @_;

  my @query_bases = ($query_start .. $query_end);

  my $qb2s = $self->query_base_to_sequence();
  my $hb2s = $self->hit_base_to_sequence();
  my $qb2hb = $self->query_base_to_hit_base();

  my $count_match = 0;
  my $count_mismatch = 0;
  my @hit_bases;
  for (my $qb = $query_start; $qb <= $query_end; $qb++) {
    my $qs = $qb2s->{$qb} || die "no sequence for query base $qb";
    if (my $hit_base = $qb2hb->{$qb}) {
      push @hit_bases, $hit_base;
      my $hs = $hb2s->{$hit_base} || die;
      if ($qs eq $hs) {
	$count_match++;
      } else {
	$count_mismatch++;
      }
    } else {
      # subject has a gap for this query position
      $count_mismatch++;
    }
  }

  my $query_span = ($query_end - $query_start) + 1;
  # span of the region of interest, not the full HSP

  my $query_site = ($query_start - $self->start("query")) /
    (($self->end("query") - $self->start("query")) + 1);
  # maybe 2 metrics?:
  # - distance from overlap start to HSP start
  # - distance from overlap end to HSP end

  my $flanking_gap = 0;
  if ($query_site >= 0.8) {
    # overlap site near the end of query, so this is the left HSP.
    # check immediately upstream for a gap,
    # test_FLT3_query_overlap_rescue_in_large_ITD.tab
    if ($self->has_subject_gap($query_start - 1, $query_start)) {
      $flanking_gap = 1;
    }
  } elsif ($query_site <= 0.2) {
    # overlap start near start of query, so this is the right HSP.
    # check immediately downtream for a gap, e.g.
    # test_right_hsp_has_flanking_gap.tab
    if ($self->has_subject_gap($query_end, $query_end + 1)) {
      $flanking_gap = 1;
    }
  } elsif (0) {
    # don't consider this fatal as can happen for dubious input
    die "can't identify site, frac=$query_site q_overlap_start=$query_start qe=" . $self->end("query");
  }

  #
  #  find gaps in alignment:
  #
  @hit_bases = sort {$a <=> $b} @hit_bases;
  my $subject_span = ($hit_bases[$#hit_bases] - $hit_bases[0]) + 1;

  my $gap_count = abs($subject_span - $query_span);

  my $count_ok = $count_match;
  $count_ok -= ($gap_count / 2);
  # penalize gaps by treating each as a half-base mismatch
  $count_ok -= ($flanking_gap / 2);
  # penalize flanking gaps as a half-base mismatch

  return $count_ok / $query_span;
}

sub get_last_clean_alignment_base {
  # iterate through a query region and return the last base
  # before any flaws in the alignment (mismatch, gap).
  my ($self, $query_start, $query_end) = @_;

  my $qb2s = $self->query_base_to_sequence();
  my $hb2s = $self->hit_base_to_sequence();
  my $qb2hb = $self->query_base_to_hit_base();

  my $count_match = 0;
  my $count_mismatch = 0;
  my @hit_bases;
  my $last_ok;
  for (my $qb = $query_start; $qb <= $query_end; $qb++) {
    my $qs = $qb2s->{$qb} || die "no sequence for query base $qb";
    if (my $hit_base = $qb2hb->{$qb}) {
      push @hit_bases, $hit_base;
      my $hs = $hb2s->{$hit_base} || die;
      if ($qs eq $hs) {
	# query base is aligned and agrees
	$last_ok = $qb;
      } else {
	# query/subject mismatch
	last;
      }
    } else {
      # subject has a gap for this query position
      last;
    }
  }
  return $last_ok;
}


sub trim_query_region {
  # discard a sub-region of the alignment.
  my ($self, $q_start, $q_end) = @_;

  my $is_start = $q_start == $self->start("query");
  my $is_end = $q_end == $self->end("query");

  my $verbose = 0;

  if ($is_start or $is_end) {
    my $qb2hb = $self->query_base_to_hit_base() || die;
    my $qb2s = $self->query_base_to_sequence();
    my $hb2qb = $self->hit_base_to_query_base();
    my $hb2qi = $self->hit_base_to_query_index();
    my $hb2s = $self->hit_base_to_sequence();

    for (my $qb = $q_start; $qb <= $q_end; $qb++) {
      printf STDERR "trim q:%d", $qb if $verbose;
      if (my $hb = $qb2hb->{$qb}) {
	# only if not a gap at this query pos
	printf STDERR " h:%d", $hb if $verbose;
	delete $hb2qb->{$hb};
	delete $hb2qi->{$hb};
	delete $hb2s->{$hb};
      }
      delete $qb2hb->{$qb};
      delete $qb2s->{$qb};
      printf STDERR "\n" if $verbose;
    }

    if (0) {
      # dump remaining indexes
      foreach (sort {$b <=> $a} keys %{$hb2qb}) {
	printf STDERR "hb2qb hit:%d query:%d\n", $_, $hb2qb->{$_};
      }
      foreach (sort {$b <=> $a} keys %{$qb2hb}) {
	printf STDERR "qb2hb query:%d hit:%d\n", $_, $qb2hb->{$_};
      }
    }

  } else {
    confess sprintf "trim interval is not beginning or end of alignment: requested=%d-%d query end=%d", $q_start, $q_end, $self->end("query");
  }
}

sub get_intersection {
  my ($self, $type, $hi_other) = @_;
  my $si_this = $self->range($type, 1);
  my $si_other = $hi_other->range($type, 1);
  return intersect $si_this $si_other;
}

sub has_subject_gap {
  # check for a gap vs. the subject alignment between the given
  # query bases
  my ($self, $q_start, $q_end) = @_;
  die "must be adjacented ordered bases" unless $q_start < $q_end and $q_start == $q_end - 1;
  my $qb2hb = $self->query_base_to_hit_base();
  my $hit_start = $qb2hb->{$q_start} || die;
  my $gap;
  if (my $hit_end = $qb2hb->{$q_end}) {
    $gap = ($hit_end - $hit_start) - 1;
    # earlier gap
  } else {
    # direct gap: no mapping to query end position
    $gap = 1;
  }
  return $gap;
}

sub get_gap_info {
  my ($self, $type) = @_;

  my $hash;
  if ($type eq "query") {
    $hash = $self->query_base_to_hit_base();
  } elsif ($type eq "subject") {
    $hash = $self->hit_base_to_query_base();
  } else {
    confess "type must be either query or subject";
  }

  my @pos = sort {$a <=> $b} keys %{$hash};
  my $min = $pos[0];
  my $max = $pos[$#pos];
  my $length = ($max - $min) + 1;

  my $mapped_count = scalar keys %{$hash};
  # this includes any alignment, even mismatches
  my $gap_count = $length - $mapped_count;

  return $gap_count;
  # TO DO: option to return as fraction, etc.
}

sub print_summary {
  my ($self, $label, $indent) = @_;
  $indent = 2 unless defined $indent;

  printf STDERR "%s%s: q=%d-%d h=%d-%d\n",
    " " x $indent,
    $label,
    $self->start("query"),
    $self->end("query"),
    $self->start("subject"), $self->end("subject");
}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
