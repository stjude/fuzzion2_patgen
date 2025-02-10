package HSPTweaker;
# tweak Bio::Search::HSP::GenericHSP alignments:
#  - trim ratty edges
#  - option to not consider gaps as mismatches (relevant for alignment
#    of fz2 patterns to genes in regions distant from the breakpoint)
#
# MNE 2/2023
#
# TO DO: adapt for edge cases like fz2 minimum "strong" match, only 5 AA

use strict;
use Exporter;
use Carp qw(confess);

use Configurable;
use MiscUtils qw(get_hash_option dump_die);

@HSPTweaker::ISA = qw(Configurable Exporter);
@HSPTweaker::EXPORT_OK = qw();

my $VERBOSE = 0;

use constant BLAST_GAP_CHAR => "-";
use constant MATCH => 1;
use constant MISMATCH => 0;

#my $DEFAULT_HSP_MAX_QUERY_SHRINK_FRACTION = 0.62;
#my $DEFAULT_HSP_MAX_QUERY_SHRINK_FRACTION = 0.60;
my $DEFAULT_HSP_MAX_QUERY_SHRINK_FRACTION = 0.33;
# CHSY3-ACSS3-01 = 41%, legit
# SMARCD1-KRT3-01 = 52%, legit
# IKZF1-ZEB2-01 = 57%, legit
# EML4-NTRK3-09 = 61%, legit?
#  => these higher tolerances can be eliminated by frame gene trimming.
#
# ZNF669-ZNF124-04 = ~85%, way too much, creates another problem

use MethodMaker qw(
		    min_run_length

                    hsp
		    trimmed
		    trim_rejected

                    count_gaps_as_matches
		    max_query_shrink_fraction

		    query_start
		    query_end
		    num_identical
		    frac_identical_query
                    compare_string
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->max_query_shrink_fraction($DEFAULT_HSP_MAX_QUERY_SHRINK_FRACTION);
  # reject trimming that reduces the query span by more than this
  # fraction.  e.g. ZNF669-ZNF124-04 has a poor alignment to one frame
  # that code would otherwise trim by ~85%, confounding results for
  # other (correct) frame alignments.
  $self->configure(%options);
  return $self;
}

sub parse {
  my ($self, $hsp) = @_;
  $self->hsp($hsp);
  my $count_gaps_as_matches = $self->count_gaps_as_matches();

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

  my $qs = $hsp->query_string();
  my $qs_raw = $qs;
  my $hs = $hsp->hit_string();
  my $hs_raw = $hs;

  my $compare = "";
  my $len = length($qs);
  for (my $i = 0; $i < $len; $i++) {
    my $q = substr($qs, $i, 1);
    my $h = substr($hs, $i, 1);
    my $add;
    if ($q eq BLAST_GAP_CHAR or $h eq BLAST_GAP_CHAR) {
      $add = BLAST_GAP_CHAR;
    } elsif ($q eq $h) {
      $add = MATCH;
    } else {
      $add = MISMATCH;
    }
    $compare .= $add;
  }

  my $min_run_length = $self->min_run_length() || die "min_run_length not set";

  my $wanted = MATCH x $min_run_length;

  my $idx_start_block = index($compare, $wanted);
  my $idx_end_block = rindex($compare, $wanted);
  my $trimmed = 0;

  my $query_start = $hsp->start("query");
  my $query_end = $hsp->end("query");
  my $num_identical = $hsp->num_identical();
  my $frac_identical_query = $hsp->frac_identical("query");
  # starting values from raw HSP

  if ($idx_end_block != -1) {
    # trim sequence at end
    my $i = $idx_end_block + $min_run_length;
    if ($i < $len) {
      $qs = substr($qs, 0, $i);
      $hs = substr($hs, 0, $i);
      $compare = substr($compare, 0, $i);

      my $stop = $i - 1;

      for (my $j = $len - 1; $j > $stop; $j--) {
	my $char = substr($qs_raw, $j, 1);
#	printf STDERR "%d: %s\n", $j, $char;
	$query_end-- unless $char eq BLAST_GAP_CHAR;
      }
    }
    $trimmed = 1;
  }

  if ($idx_start_block > 0) {
    # trim sequence at start
    $qs = substr($qs, $idx_start_block);
    $hs = substr($hs, $idx_start_block);
    $compare = substr($compare, $idx_start_block);

    for (my $i = 0; $i < $idx_start_block; $i++) {
      my $char = substr($qs_raw, $i, 1);
#      printf STDERR "%d: %s\n", $i, $char;
      $query_start++ unless $char eq BLAST_GAP_CHAR;
    }
    $trimmed = 1;
  }

  if ($trimmed or $count_gaps_as_matches) {
    # recalculate identities
    my @c = split //, $compare;
    if ($count_gaps_as_matches) {
      $num_identical = grep {$_ eq MATCH or $_ eq BLAST_GAP_CHAR} @c;
    } else {
      $num_identical = grep {$_ eq MATCH} @c;
    }
    $frac_identical_query = $num_identical / scalar @c;
  }

  my $trim_rejected = 0;
  if ($trimmed) {
    my $rqs = $hsp->start("query");
    my $rqe = $hsp->end("query");
    my $raw_size = ($rqe - $rqs) + 1;
    my $cooked_size = ($query_end - $query_start) + 1;

    my $shrink_fraction = 1 - ($cooked_size / $raw_size);

    if ($shrink_fraction >= $self->max_query_shrink_fraction) {
      printf STDERR "trimming rejected: raw_size:%d cooked:%d shrunk:%.2f%%\n",
	$raw_size, $cooked_size, $shrink_fraction * 100;
      $trimmed = 0;
      $trim_rejected = $shrink_fraction;
    }

  }

  $self->compare_string($compare);
  $self->trimmed($trimmed);
  $self->trim_rejected($trim_rejected);
  $self->query_start($query_start);
  $self->query_end($query_end);
  $self->num_identical($num_identical);
  $self->frac_identical_query($frac_identical_query);
  # save results
}

sub score {
  my ($self) = @_;
  return $self->hsp()->score();
}

sub frac_identical {
  my ($self, $type) = @_;
  if ($type eq "query") {
    return $self->frac_identical_query();
  } else {
    confess "only query implemented so far";
  }
}

sub range {
  my ($self, $type) = @_;
  my ($start, $end);
  if ($type eq "query") {
    $start = $self->query_start();
    $end = $self->query_end();
  } else {
    confess "only query implemented so far";
  }
  return($start, $end);
}

sub end {
  my ($self, $type) = @_;
  my $v;
  if ($type eq "query") {
    $v = $self->query_end();
  } else {
    confess "only query implemented so far";
  }
  return $v;
}

sub start {
  my ($self, $type) = @_;
  my $v;
  if ($type eq "query") {
    $v = $self->query_start();
  } else {
    confess "only query implemented so far";
  }
  return $v;

}

sub check_edge_identity {
  my ($self, %options) = @_;
  my $side = $options{"-side"} || die;
  my $chunk_size = $options{"-length"} || die;

  my $compare_string = $self->compare_string() || die;

  my $frac_identical = 0;

  if (length($compare_string) >= $chunk_size) {
    # if enough enough sequence to meet minimum
    my $check;
    if ($side eq "A") {
      # for gene A, fusion breakpoint site is on the right edge of alignment
      $check = substr($compare_string, - $chunk_size);
    } elsif ($side eq "B") {
      # for gene B, breakpoint is on left side
      $check = substr($compare_string, 0, $chunk_size);
    } else {
      die;
    }

    my @c = split //, $check;
    my $num_identical;
    if ($self->count_gaps_as_matches) {
      $num_identical = grep {$_ eq MATCH or $_ eq BLAST_GAP_CHAR} @c;
    } else {
      $num_identical = grep {$_ eq MATCH} @c;
    }
    # TO DO: centralize/share this block
    $frac_identical = $num_identical / $chunk_size;
  }

  return $frac_identical;

}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
